//
// Copyright 2016 MIT Lincoln Laboratory, Massachusetts Institute of Technology
//
// Licensed under the Apache License, Version 2.0 (the "License"); you may not use these files except in compliance with
// the License.
//
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on
// an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the
// specific language governing permissions and limitations under the License.
//

//
// Clustering of intervals from a segmentation
// A segmentation is a list of intervals
// An interval is a contiguous set of feature vectors
//
// Complete rewrite of LL tool 
// Loosely based on old speaker clustering tool 
// (uses container classes and C++ style coding -- basically a data structures exercise)
// 
// References -- Zhou/Hansen paper, Reynolds paper
//
// BC, 2/2014

#include "stdafx.h"
#include <stdlib.h>
#include "speech_tools.h"
#include <sstream>

Cluster_BIC_params::Cluster_BIC_params (void)
{
   bic_thresh = 0.0;
   bic_weight = 1.0;
   min_num_clusters = 1;
   max_num_clusters = -1;
   verbose = false;
}

Segmentation_BIC_info::Segmentation_BIC_info (Segmentation &seg) : 
   s1(seg.feats.num_outfeat()), s2(seg.feats.num_outfeat()*seg.feats.num_outfeat())
{
   // This class contains info for computing delta BIC quickly:
   // stats (0,1,2) and the scaled log of the determinant of the covariance 
   int i;
   dim = seg.feats.num_outfeat();
   vec<float> v(dim);
   n = 0;
   log_det = 0;
   s1.constant(0.0);
   s2.constant(0.0);

   // Iterate over intervals in segmentation
   // Most of the time only one interval, but make this general
   list<Interval>::iterator it;
   for (it=seg.intervals.begin(); it!=seg.intervals.end(); it++) {
      n += (float) (it->end-it->start+1);
      for (i=it->start; i<=it->end; i++) {
         v = seg.feats.vector(i);
         s1.add(v);
         s2.outer_add(v);
      }
   }

   // Compute log_det
   compute_log_det();

}

void Segmentation_BIC_info::compute_log_det (void)
{
   vec<float> temp_mat(dim*dim), temp_vec(dim);
   temp_vec = s1;
   temp_vec.scale((float)(1.0/n));
   temp_mat = s2;
   temp_mat.scale((float)(1.0/n));
   temp_mat.outer_sub(temp_vec);
   log_det = (float) (0.5*n*log_det_matrix(temp_mat, (float) 1.0e-6));
}

Cluster_segmentation::Cluster_segmentation (Segmentation &seg_init) : dist((int)(seg_init.intervals.size()*seg_init.intervals.size())) 
{
   list<Interval>::iterator it_int;

   // Create list of segments.  One interval per segment.
   for (it_int=seg_init.intervals.begin(); it_int!=seg_init.intervals.end(); it_int++) {
      seg_list.push_back(Segmentation(seg_init.feats, it_int->start, it_int->end)); // default copy const. should be ok
      Segmentation &seg = seg_list.back();
      seg.bw = seg_init.bw;
      seg.gender = seg_init.gender;
      seg.speaker = seg_init.speaker;
      seg.speech = seg_init.speech;
   }
}

void Cluster_segmentation::bic_init (Cluster_BIC_params &params_in)
{
   list<Segmentation>::iterator it_seg_i, it_seg_j;
   list<Segmentation_BIC_info>::iterator it_bi_i, it_bi_j;
   params = params_in;

   // Initial BIC info stats for each segment
   for (it_seg_i=seg_list.begin(); it_seg_i!=seg_list.end(); it_seg_i++)
      bic_info.push_back(Segmentation_BIC_info(*it_seg_i));

   // Compute distance matrix
   int i, j, n;
   n = (int) seg_list.size();
   for (it_seg_i=seg_list.begin(), it_bi_i=bic_info.begin(), i=0; it_seg_i!=seg_list.end(); it_seg_i++, it_bi_i++, i++) {
      for (it_seg_j=seg_list.begin(), it_bi_j=bic_info.begin(), j=0; it_seg_j!=seg_list.end(); it_seg_j++, it_bi_j++, j++) {
         if (j>i)
            break;
         dist.data[i+j*n] = delta_bic(*it_bi_i, *it_bi_j); 
         dist.data[j+i*n] = dist.data[i+j*n]; // redundant -- should just store upper triangular portion
      }
   }
}

static bool compare_intervals (Interval &int1, Interval int2)
{
   if (int1.start<int2.start)
      return true;
   return false;
}

void Cluster_segmentation::bic_cluster (void)
{
   int num_clusters = seg_list.size();
   if (num_clusters<params.min_num_clusters)
      return;

   int i_min, j_min;
   float dbic_min;
   int local_max_clusters = params.max_num_clusters;
   if (params.max_num_clusters==-1)
      local_max_clusters = num_clusters+1;

   // Compute bic for current clusters
   find_min_dist(i_min, j_min, dbic_min);

   // Main loop -- combine one at a time until the minimum is above the threshold
   while ((num_clusters > params.min_num_clusters) && (num_clusters>1) && ((dbic_min<params.bic_thresh) || (num_clusters>local_max_clusters))) {
      // Merge closest segments
      merge(i_min, j_min);

      // Find closest segments to merge
      find_min_dist(i_min, j_min, dbic_min);
      num_clusters = seg_list.size();

      if (params.verbose) {
         printf("BIC clustering -- minimum distance is %f\n", dbic_min);
         printf("BIC clustering -- number of clusters = %d\n", num_clusters);
      }

   }

   // Sort segments
   list<Segmentation>::iterator it_seg;
   for (it_seg=seg_list.begin(); it_seg!=seg_list.end(); it_seg++)
      it_seg->intervals.sort(compare_intervals);

   // Combine adjacent segments
   for (it_seg=seg_list.begin(); it_seg!=seg_list.end(); it_seg++) {
      list<Interval>::iterator it_int, it_last;
      it_last = it_seg->intervals.begin();
      it_int = it_last;
      it_int++;
      while (it_int!=it_seg->intervals.end()) {
         if ((it_int->start-it_last->end)==1) {
            it_last->end = it_int->end;
            it_int = it_seg->intervals.erase(it_int);
         } else {
            it_last = it_int;
            it_int++;
         }
      }
   }

}

float Cluster_segmentation::delta_bic (Segmentation_BIC_info &b1, Segmentation_BIC_info &b2)
{
   vec<float> mean(b1.s1.len), covar(b1.s2.len);
   float n, log_det, penalty, dim;

   n = b1.n + b2.n;

   mean = b1.s1;
   mean.add(b2.s1);
   mean.scale(1.0/n);

   covar = b1.s2;
   covar.add(b2.s2);
   covar.scale(1.0/n);
   covar.outer_sub(mean);

   log_det = (float) (0.5*n*log_det_matrix(covar, (float) 1.0e-6));
   log_det -= b1.log_det+b2.log_det;

   dim = seg_list.front().feats.num_outfeat();
   penalty = 0.5*params.bic_weight*(dim+0.5*dim*(dim+1))*log(n);
   log_det -= penalty;

   return log_det;

}

void Cluster_segmentation::find_min_dist (int &i_min, int &j_min, float &dist_min)
{
   int i, j;
   int num_seg = (int) seg_list.size();
   float val;

   // Assume at least two segments
   i_min = 1;
   j_min = 0;
   dist_min = dist.data[num_seg];

   for (i=2; i<num_seg; i++) {
      for (j=0; j<i; j++) { // avoid the diagonal
         val = dist.data[i*num_seg+j];
         if (val < dist_min) {
            i_min = i;
            j_min = j;
            dist_min = val;
         }
      }
   }

   // swap so i_min < j_min
   i = i_min;
   i_min = j_min;
   j_min = i;

}

void Cluster_segmentation::merge (int i1, int j1)
{
   if (params.verbose)
      printf("merging: %d %d\n", i1, j1);

   // merge segments i1 and j1
   list<Segmentation>::iterator it1 = seg_list.begin();
   std::advance(it1, i1);
   list<Segmentation>::iterator it2 = seg_list.begin();
   std::advance(it2, j1);

   // Merge list of intervals
   list<Interval>::iterator iti = it1->intervals.end();
   it1->intervals.insert(iti, it2->intervals.begin(), it2->intervals.end());
   seg_list.erase(it2);

   // Merge bic_info
   list<Segmentation_BIC_info>::iterator it_bi_i, it_bi_j;
   it_bi_i = bic_info.begin();
   std::advance(it_bi_i, i1);
   it_bi_j = bic_info.begin();
   std::advance(it_bi_j, j1);
   it_bi_i->n += it_bi_j->n;
   it_bi_i->s1.add(it_bi_j->s1);
   it_bi_i->s2.add(it_bi_j->s2);
   it_bi_i->compute_log_det();
   bic_info.erase(it_bi_j);

   // Remove row and column j1
   int num_seg = (int) seg_list.size()+1;
   int i, j, i_new;
   for (i=0; i<num_seg; i++) {
      if (i==j1)
         continue;
      if (i<j1)
         i_new = i;
      else
         i_new = i-1;
      for (j=0; j<j1; j++)
         dist.data[i_new*(num_seg-1)+j] = dist.data[i*num_seg+j];
      for (j=j1+1; j<num_seg; j++)
         dist.data[i_new*(num_seg-1)+j-1] = dist.data[i*num_seg+j];
   }

   // Recompute row i1 of distance matrix
   num_seg--;
   list<Segmentation_BIC_info>::iterator it_bi = bic_info.begin();
   for (j=0; j<num_seg; j++, it_bi++) {
      dist.data[i1*num_seg+j] = delta_bic(*it_bi_i, *it_bi);
      dist.data[j*num_seg+i1] = dist.data[i1*num_seg+j];
   }

}

void Cluster_segmentation::save_marks (string &filename, float win_inc_ms, bool audacity_format)
{
   FILE *outfile = fopen(filename.c_str(), "w");
   if (outfile==NULL)
      throw ST_exception("Cluster_segmentation::save_marks -- unable to open marks file.");

   list<Segmentation>::iterator it_seg;
   float seg_start, seg_end;
   int cnum;
   for (cnum=0, it_seg=seg_list.begin(); it_seg!=seg_list.end(); it_seg++, cnum++) {
      it_seg->intervals.sort(compare_intervals);
      list<Interval>::iterator it;
      for (it=it_seg->intervals.begin(); it!=it_seg->intervals.end(); it++) {
         seg_start = 0.001*it->start*win_inc_ms;
         seg_end = 0.001*it->end*win_inc_ms;
         if (audacity_format)
            fprintf(outfile, "%.2f %.2f c%d\n", seg_start, seg_end, cnum);
         else
            fprintf(outfile, "c%d %.2f %.2f\n", cnum, seg_start, seg_end-seg_start);
      }
   }
   fclose(outfile);

}

string Cluster_segmentation::to_string (float win_inc_ms)
{
   stringstream ss;

   list<Segmentation>::iterator it_seg;
   float seg_start, seg_end;
   int cnum;
   for (cnum=0, it_seg=seg_list.begin(); it_seg!=seg_list.end(); it_seg++, cnum++) {
      it_seg->intervals.sort(compare_intervals);
      list<Interval>::iterator it;
      for (it=it_seg->intervals.begin(); it!=it_seg->intervals.end(); it++) {
         seg_start = floor(it->start*win_inc_ms);
         seg_end = floor(it->end*win_inc_ms);
         ss << fixed << "<segment start=\"" << ((int) seg_start)  << "\" end=\"" << ((int) seg_end) << "\">\n";
         ss  << "<model name=\"s" << cnum << "\">1</model>\n";
         ss << "</segment>\n";
      }
   }
   return ss.str();

}


