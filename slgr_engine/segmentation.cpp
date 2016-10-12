//
// Copyright 2012-2015 Massachusetts Institute of Technology. All rights
// reserved. Title to and copyright in the speaker, language, gender recognition
// software, any derivatives and any associated documentation is retained by
// MIT, subject to a non-exclusive royalty-free license to the U.S. Government
// to use the software for Government purposes as defined in Federal 
// Acquisition Regulation. Distribution to or use by third parties without 
// the prior written authorization of MIT is expressly prohibited.
//
// MIT MAKES NO REPRESENTATIONS OR WARRANTIES, EXPRESS OR IMPLIED, WITH RESPECT
// TO THE FUNCTIONALITY OR USE OF THE SOFTWARE. ALL IMPLIED WARRANTIES, 
// INCLUDING BUT NOT LIMITED TO, WARRANTIES OF MERCHANTABILITY OR FITNESS FOR 
// A PARTICULAR PURPOSE ARE EXPRESSLY EXCLUDED. FURTHER, MIT MAKES NO WARRANTY 
// THAT THE USE OF THE SOFTWARE COMPONENTS OR DOCUMENTATION WILL NOT INFRINGE 
// ANY PATENTS, COPYRIGHTS, TRADEMARKS, TRADE SECRETS OR OTHER RIGHTS OF ANY 
// PARTY.
//
// MIT DISCLAIMS ANY AND ALL LIABILITY FOR SPECIAL, INCIDENTAL, OR 
// CONSEQUENTIAL DAMAGES (INCLUDING LOSS OF PROFIT) ARISING OUT OF THE 
// INSTALLATION, USE, OPERATION OR SUPPORT OF THE SOFTWARE OR ASSOCIATED 
// DOCUMENTATION, EVEN IF MIT HAS BEEN APPRISED OF THE POSSIBILITY OF SUCH 
// DAMAGES.

//
// Segmentation of a stream of feature vectors
//
// Complete rewrite of LL tool 
// BC, 12/2013
//
// Loosely based on old seg tool.  Reference -- Zhou/Hansen paper
//

#include "stdafx.h"
#include "speech_tools.h"

class compute_delta_BIC {
   vec<float> mean1, covar1, mean2, covar2, mean, covar;
   float n1, n2, n;
   int vec_dim;
public:
   float compute (Features &f, int start, int stop, int point, T2_BIC_params &params, vec<float> &stats);
   compute_delta_BIC(int vec_dim);
   ~compute_delta_BIC(void);
};

class compute_T2_distance {
   vec<float> mean1, mean2;
   vec<float> covar;
   int vec_dim, curr_start, curr_stop;
public:
   compute_T2_distance(int vec_dim);
   ~compute_T2_distance(void);
   int find_max(Features &f, int start, int start_search, int stop);
   float compute(Features &f, int start, int stop, int point);
};

int  cholesky_factor (vec<float> &A);
void Lsolve (vec<float> &A, vec<float> &b, vec<float> &x);
int matrix_inverse (vec<float> &A, float lambda);
void Usolve (vec<float> &A, vec<float> &b, vec<float> &x);

const int Interval::UNK = -1;
const int Interval::NONSPEECH = 0;
const int Interval::SPEECH = 1;
const int Segmentation::UNK = -1;
const int Segmentation::NONSPEECH = 0;
const int Segmentation::SPEECH = 1;
const int Segmentation::MALE = 0;  // obviously
const int Segmentation::FEMALE = 1;
const int Segmentation::NB = 0;
const int Segmentation::WB = 1;

Interval::Interval (int start_in, int end_in, int label_in, float score_in) {
   start = start_in;
   end = end_in;
   if (end<start)
      throw ST_exception("Start < end in segment construction.");
   label = label_in;
   score = score_in;
}

Segmentation::Segmentation (Features &feats_in) : intervals(1, Interval(0, feats_in.num_vec()-1)), feats(feats_in), bic_stats(6) {
   // Default segmentation is all frames
   speech = UNK;
   gender = UNK;
   bw = UNK;
   speaker = UNK;
}

Segmentation::Segmentation (Features &feats_in, int start, int end) : intervals(1, Interval(start, end)), feats(feats_in), bic_stats(6) {
   // Default segmentation is all frames
   speech = UNK;
   gender = UNK;
   bw = UNK;
   speaker = UNK;
}

vec<float> &Segmentation::get_bic_stats (void)
{
   return bic_stats;
}

void Segmentation::intervals_from_SAD (void) {
   if (!feats.SAD_available())
      throw ST_exception(string("Segmentation::intervals_from_SAD -- no SAD available."));

   intervals.clear();

   int seg_start, seg_end;
   bool in_seg = false;
   int num_vectors = feats.num_vec();

	for (int i=0; i<num_vectors; i++) {
      if (!in_seg && feats.SAD(i)==1) {
         seg_start = i;
         in_seg = true;
      } else if (in_seg && feats.SAD(i)==0) {
         seg_end = i;
         intervals.push_back(Interval(seg_start, seg_end));
         in_seg = false;
      }
   }

   if (in_seg) {
      seg_end = num_vectors-1;
      intervals.push_back(Interval(seg_start, seg_end));
   }

}

T2_BIC_params::T2_BIC_params(void) {
   init_search_interval = 100;
   max_seg_size = 500;
   min_seg_size = 25;
   search_inc = 50;
   bic_weight = 1.0;
   bic_thresh = 0.0;
   verbose = false;
}

void Segmentation::print(void)
{
   list<Interval>::iterator it;
   for (it=intervals.begin(); it!=intervals.end(); it++)
      printf("[%d,%d] ", it->start, it->end);
   printf("\n");
}

void Segmentation::save_marks (string &filename, float win_inc_ms, bool audacity_format) 
{
   float seg_start, seg_end;
   list<Interval>::iterator it;
   FILE *outfile = fopen(filename.c_str(), "w");
   if (outfile==NULL)
      throw ST_exception(string("Segmentation -- cannot save marks file: ") + filename);
   for (it=intervals.begin(); it!=intervals.end(); it++) {
      seg_start = 0.001*it->start*win_inc_ms;
      seg_end = 0.001*it->end*win_inc_ms;
      if (audacity_format)
         fprintf(outfile, "%.2f %.2f int\n", seg_start, seg_end);
      else
         fprintf(outfile, "int %.2f %.2f\n", seg_start, seg_end-seg_start);
   }
   fclose(outfile);
}

void Segmentation::segment_T2_BIC_one_interval (T2_BIC_params &params) {
   // One interval case
   int start = intervals.front().start;
   int start_search = start+params.min_seg_size;  // don't create intervals smaller than min_seg_size
   int end_search = start_search+params.init_search_interval;
   int imax_t2, new_cp;
   float dbic;
   compute_T2_distance t2_dist(feats.num_outfeat());
   compute_delta_BIC delta_bic(feats.num_outfeat());
   vector<int> cp;

   // Find change points
   while (end_search <= intervals.back().end) {
      // Look for maximum T^2 statistic in interval
      imax_t2 = t2_dist.find_max(feats, start, start_search, end_search);
      if (params.verbose)
         printf("start = %d, end_search = %d, point = %d\n", start, end_search, imax_t2);

      // If we're going to split and make too small of interval going forward, 
      // then stop looking for change points
      if ((intervals.back().end-imax_t2) < params.min_seg_size)
         break;

      // Now evaluate BIC at this point
      dbic = delta_bic.compute(feats, start, end_search, imax_t2, params, bic_stats);
      if (params.verbose) 
         printf("\tdbic = %f\n", dbic);

      if (dbic > params.bic_thresh) {
         if (params.verbose)
            printf("\tfound change point\n");
         bic_stats.data[5] += 1;
         cp.push_back(imax_t2);
         start = imax_t2+1;
         start_search = start+params.min_seg_size;
         end_search = start_search+params.init_search_interval;
      } else {
         end_search += params.search_inc;
         if ((end_search-start+1) > params.max_seg_size) { 
            if (params.verbose) 
               printf("\tMAX interval reached -- forcing change point\n");
            // force a change point if it doesn't create small segments
            // otherwise, the computation for change points will be too large
            int new_cp = start + params.max_seg_size;
            if ((intervals.back().end-new_cp)<params.min_seg_size)
               break;
            start = new_cp+1;
            start_search = start+params.min_seg_size;
            end_search = start_search+params.init_search_interval;
            cp.push_back(new_cp);
         }
      }
   }

   // Now chop up the interval according to the change points
   int i = 0, curr_pt, last_end, lbl;
   vector<int>::iterator it;
   lbl = intervals.front().label;
   for (it=cp.begin(); it!=cp.end(); it++) {
      curr_pt = *it;
      last_end = intervals.back().end;
      intervals.back().end = curr_pt;
      intervals.push_back(Interval(curr_pt+1, last_end));
      intervals.back().label = lbl; // propagate same label
   }
   return;

}

void Segmentation::segment_T2_BIC_multiple_intervals (T2_BIC_params &params, bool speech_only) 
{
   //
   // More than one segment case -- reduce to the one segment case and use recursion
   //
   list<Interval>::iterator it = intervals.begin();
   while (it!=intervals.end()) {
      // if speech only and not speech, skip
      if (speech_only && it->label!=it->SPEECH) {
         it++;
         continue;
      }

      // if it's too small, skip
      if ((it->end-it->start+1) < params.min_seg_size) {
         it++;
         continue;
      }

      // Reduce to one segment case
      Segmentation seg_one(feats, it->start, it->end);
      seg_one.intervals.front().label = it->label;
      seg_one.segment_T2_BIC_one_interval(params);

      // Update stats
      bic_stats.add(seg_one.get_bic_stats());

      if (seg_one.intervals.size()>1) { // Insert the chopped up interval into the vector of intervals
         intervals.insert(it, seg_one.intervals.begin(), seg_one.intervals.end());
         it = intervals.erase(it);
      } else { // No change
         it++;
      }
   }

}

void Segmentation::segment_T2_BIC (T2_BIC_params &params, bool speech_only) {

   int num_intervals = intervals.size();
   bic_stats.constant(0.0);

   segment_T2_BIC_multiple_intervals(params, speech_only);

   // Final bic stats
   if (bic_stats.data[0] > 1) {
      float n = bic_stats.data[0];
      bic_stats.scale(1.0/bic_stats.data[0]);
      bic_stats.data[0] = n;
      bic_stats.data[2] -= bic_stats.data[1]*bic_stats.data[1];
      bic_stats.data[4] -= bic_stats.data[3]*bic_stats.data[3];
      bic_stats.data[2] = sqrt(1.0e-6+bic_stats.data[2]);
      bic_stats.data[4] = sqrt(1.0e-6+bic_stats.data[4]);
   }

   if (params.verbose) {
      printf("bic_stats[0], number of candidate change points = %f\n", bic_stats.data[0]);
      printf("bic_stats[1], mean of log likelihood terms = %f\n", bic_stats.data[1]);
      printf("bic_stats[2], standard deviation of log likelihood terms = %f\n", bic_stats.data[2]);
      printf("bic_stats[3], mean of penalty term = %f\n", bic_stats.data[3]);
      printf("bic_stats[4], standard deviation of penalty term = %f\n", bic_stats.data[4]);
      printf("bic_stats[5], percentage of delta BIC above the threshold = %f\n", bic_stats.data[5]);
   }

}

compute_T2_distance::compute_T2_distance(int vec_dim_in) : mean1(vec_dim_in), mean2(vec_dim_in), 
   covar(vec_dim_in*vec_dim_in)
{
   vec_dim = vec_dim_in;
   curr_start = -1;
   curr_stop = -1;

} // constructor

compute_T2_distance::~compute_T2_distance(void)
{
   // Nothing to do
}

float compute_T2_distance::compute(Features &f, int start, int stop, int point)
{
   // Compute Hotelling T^2 statistic 
   int i, j, k;
   float sum, c, n, n1, n2;

   vec<REAL_FEAT> v(f.num_outfeat()), mean1(f.num_outfeat()), mean2(f.num_outfeat());

   n = stop-start+1;

   // If start or stop has changed, we need to recompute the covariance matrix for the segment
   if (start!=curr_start || stop!=curr_stop) {
      curr_start = start;
      curr_stop = stop;

      mean1.constant(0.0);
      covar.constant(0.0);

      // Loop and accumulate
      for (k=start; k<=stop; ++k) {
         v = f.vector(k);
         mean1.add(v);
         covar.outer_add(v);
      }

      // Now scale by number of instances and compute covariance
      c = 1.0/n;
      mean1.scale(c);
      covar.scale(c);
      covar.outer_sub(mean1);

      // Matrix inverse
      float lambda = 1.0e-6;
      int inv_success;
      for (k=0; k<4; k++) {
         // Start small with regularization
         // 1.0e-6 should be good if features are N(0,1)
         // Diagonal of covariance should be close to 1.0 as long as features are not all zero.
         inv_success = matrix_inverse(covar, 1.0e-6);
         if (inv_success)
            break;
         lambda *= 10;
      }
      if (!inv_success)
         throw ST_exception("Unexpected error in segmentation: ill-conditioned matrix in T^2 statistic");
   }

   // Compute means
   mean1.constant(0.0);
   for (k=start; k<=point; k++) {
      v = f.vector(k);
      mean1.add(v);
   }
   n1 = (float)(point-start+1);
   c = ((float) 1.0)/n1;
   mean1.scale(c);

   mean2.constant(0.0);
   for (k=point+1; k<=stop; k++) {
      v = f.vector(k);
      mean2.add(v);
   }
   n2 = (float)(stop-point);
   c = ((float) 1.0)/n2;
   mean2.scale(c);

   //  dist = (n1*n2/n)*(m1-m2)'*inv_covar*(m1-m2)
   float dist = 0;
   mean1.subtract(mean2);
   n = stop-start+1;
   for (i=0; i<vec_dim; i++) {
      sum = 0;
      for (j=0; j<vec_dim; j++)
         sum += covar.data[j+i*mean1.len]*mean1.data[j];
      dist += sum*mean1.data[i];
   }
   dist *= (n1*n2)/n;

   return dist;

}

int compute_T2_distance::find_max(Features &f, int start, int start_search, int stop)
{
   int i;
   int imx = start_search;
   float val, mx = this->compute(f, start, stop, start_search);

   for (i=start_search+1; i<stop; i++) {  // doesn't make sense when i=stop
      val = this->compute(f, start, stop, i);
      if (val > mx) {
         mx = val;
         imx = i;
      }
   }
   return imx;
}

compute_delta_BIC::compute_delta_BIC (int vec_dim_in) : mean1(vec_dim_in), mean2(vec_dim_in), mean(vec_dim_in), 
   covar1(vec_dim_in*vec_dim_in), covar2(vec_dim_in*vec_dim_in), covar(vec_dim_in*vec_dim_in)
{
   vec_dim = vec_dim_in;
}

compute_delta_BIC::~compute_delta_BIC (void)
{
   // nothing to do
}

float compute_delta_BIC::compute (Features &f, int start, int stop, int point, T2_BIC_params &params, vec<float> &stats) 
{
   int i;
   vec<float> v(f.num_outfeat());

   // Compute stats for first interval
   mean1.constant(0.0);
   covar1.constant(0.0);
   for (i=start; i<=point; i++) {
      v = f.vector(i);
      mean1.add(v);
      covar1.outer_add(v);
   }
   n1 = (float)(point-start+1);

   // Compute stats for second interval
   mean2.constant(0.0);
   covar2.constant(0.0);
   for (i=point+1; i<=stop; i++) {
      v = f.vector(i);
      mean2.add(v);
      covar2.outer_add(v);
   }
   n2 = (float)(stop-point);

   // Combined stats
   mean = mean1;
   mean.add(mean2);
   covar = covar1;
   covar.add(covar2);
   n = n1+n2;

   // Now scale and compute final means and covariances
   float c = 1.0/n1;
   mean1.scale(c);
   covar1.scale(c);
   covar1.outer_sub(mean1);
   c = 1.0/n2;
   mean2.scale(c);
   covar2.scale(c);
   covar2.outer_sub(mean2);
   c = 1.0/n;
   mean.scale(c);
   covar.scale(c);
   covar.outer_sub(mean);

   // Compute log of determinants
   float l1 = log_det_matrix(covar1, 1.0e-6);
   float l2 = log_det_matrix(covar2, 1.0e-6);
   float l12 = log_det_matrix(covar, 1.0e-6);
   
   if (params.verbose)
      printf("l1 = %f , l2 = %f , l12 = %f\n", l1, l2, l12);
      
   float delta_BIC = 0.5*n*l12-0.5*n1*l1-0.5*n2*l2;
   float vec_dim = (float) f.num_outfeat();
   if (params.verbose)
      printf("\tdbic before penalty = %f\n", delta_BIC);

   float BIC_penalty = 0.5*params.bic_weight*(vec_dim+0.5*(vec_dim*(vec_dim+1)))*log(n);
   if (params.verbose)
      printf("\tBIC penalty = %f\n", BIC_penalty);

   stats.data[0] += 1;
   stats.data[1] += delta_BIC;
   stats.data[2] += delta_BIC*delta_BIC;
   stats.data[3] += BIC_penalty;
   stats.data[4] += BIC_penalty*BIC_penalty;

   delta_BIC -= BIC_penalty;
   return delta_BIC;
}

int cholesky_factor (vec<float> &A) {
   int     i, j, k, n;
   float  *A_piv, *A_row, ip, sum, tmp;

   n = (int) sqrt((float) A.len);

   for (k=0; k<n; k++) {
      // diagonal element
      sum = A.data[k*n+k];
      A_piv = &A.data[k*n];
      for (j=0; j<k; j++) {
         tmp = *A_piv++;
         sum -= tmp*tmp;
      }

      if (sum <= 0.0) {
         // printf("Warning: Cholesky factor -- matrix not positive definite.  Diagonal element = %lf.\n", sum);
         return 0;
      }
      A.data[k*n+k] = sqrt(sum);

      // set values of column k
      for (i=k+1; i<n; i++) {
         sum = A.data[i*n+k];
         A_piv = &A.data[k*n];
         A_row = &A.data[i*n];
         ip = 0;
         for (j=0; j<k; j++)
            ip += A_row[j]*A_piv[j];
         sum -= ip;
         A.data[j*n+i] = A.data[i*n+j] = sum/A.data[k*n+k];
      }

   } // for

   return 1;

} // cholesky_factor

void Lsolve (vec<float> &A, vec<float> &b, vec<float> &x) {
   int i, j, n;

   n = x.len;
   for (i=0; i<n; i++) 
      x.data[i] = b.data[i];

   for (i=0; i<n; i++) {
      x.data[i] = x.data[i]/A.data[i*n+i];
      for (j=i+1; j<n; j++)
         x.data[j] -= A.data[i*n+j]*x.data[i];
   }

} // Lsolve 

void Usolve (vec<float> &A, vec<float> &b, vec<float> &x)
{
   int i, j, n;
   
   n = x.len;
   for (i=0; i<n; i++) 
      x.data[i] = b.data[i];

   for (i=n-1; i>=0; i--) {
      x.data[i] = x.data[i]/A.data[i*n+i];
      for (j=0; j<i; j++)
         x.data[j] -= A.data[i*n+j]*x.data[i];
   }

} // Usolve 

float log_det_matrix (vec<float> &A, float lambda) {
    // Determinant of a non-negative definite matrix with regularization
   int n = (int) sqrt((float) A.len), ok;
   vec<float> A_temp(A.len);
   float log_det = 0.0;

   A_temp = A;
   A_temp.add_diag(lambda);
   ok = cholesky_factor(A_temp);
   if (ok) {
      for (int i=0; i<n; i++)
         log_det += 2*log(A_temp.data[i*n+i]);
   } else {
      // If Cholesky doesn't work,
      // return the log det of lambda*eye(n)
      log_det = n*log(lambda);
   }
   
   return log_det;
}

int matrix_inverse (vec<float> &A, float lambda) {
   // Matrix inversion of a non-negative definite matrix with regularization
   int n = (int) sqrt((float) A.len);
   vec<float> A_temp(A.len), b(n), x(n);
   int i, j, ok;

   A_temp = A;
   A_temp.add_diag(lambda);
   ok = cholesky_factor(A_temp);
   if (!ok)
      return 0;
   b.constant(0.0);

   // Invert by solving L'*L*x = b using b as columns of the identity matrix
   for (i=0; i<n; i++) {
      b.data[i] = 1.0;
      Lsolve(A_temp, b, x);
      Usolve(A_temp, x, x);
      for (j=0; j<n; j++)
         A.data[j+i*n] = x.data[j];
      b.data[i] = 0.0;
   }
   return 1;

}
