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
// GMM class
//

//
// Standard headers
//
#include "stdafx.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <sys/stat.h>

// Bill Campbell, 8/2008
// BC, Modified for IPDF, 12/2009
// Speech tools, 10/28/2010
// 10/28/10, Refactored from LLSpeech, compute_state_prob.cc, vec_utils.cc
// 02/13, Added GMM scoring refactored from LLSpeech, gmm_score.cc, gmm_score_subs.cc

#include "speech_tools.h"

//
// Class member functions
// 
GMM_model::GMM_model (void) {
	num_fea = 0;
	num_mix = 0;
	mean = 0;
	inv_cov = 0;
	log_weight = 0;
	log_det = 0;
	ladd_init();  
	loaded = false;
}

GMM_model::GMM_model (const GMM_model &in) { // Next two methods are here so I can use GMM models with STL
	if (in.mean==0) { // empty model input
		num_fea = 0;
		mean = 0;
		inv_cov = 0;
		log_weight = 0;
		log_det = 0;
		ladd_init();
		loaded = false;
		return;
	}

	// Copy weights and log det
	int i;
	num_mix = in.num_mix;
	log_weight = new float[num_mix];
	log_det = new float[num_mix];
	for (i=0; i<num_mix; i++) {
		log_weight[i] = in.log_weight[i];
		log_det[i] = in.log_det[i];
	}

	// Copy mean and covariance
	num_fea = in.num_fea;
	mean = new float[num_fea*num_mix];
	inv_cov = new float[num_fea*num_mix];
	for (i=0; i<(num_fea*num_mix); i++) {
		mean[i] = in.mean[i];
		inv_cov[i] = in.inv_cov[i];
	}
	ladd_init();  
	loaded = true;
}

GMM_model& GMM_model::operator= (const GMM_model &in) {
   int i;

   if (in.mean == 0) {
      mean = 0;
      inv_cov = 0;
      log_weight = 0;
      log_det = 0;
      loaded = false;
   } else {
      num_mix = in.num_mix;
      num_fea = in.num_fea;

      // Initialize and copy over data
      log_weight = new float[num_mix];
      log_det = new float[num_mix];
      for (i=0; i<num_mix; i++) {
         log_weight[i] = in.log_weight[i];
         log_det[i] = in.log_det[i];
      }
      mean = new float[num_fea*num_mix];
      inv_cov = new float[num_fea*num_mix];
      for (i=0; i<(num_fea*num_mix); i++) {
         mean[i] = in.mean[i];
         inv_cov[i] = in.inv_cov[i];
      }
      loaded = true;
   }

	ladd_init();  
	return *this;
}

GMM_model::~GMM_model (void) {
	delete[] mean;
	delete[] inv_cov;
	delete[] log_weight;
	delete[] log_det;
}

void GMM_model::compute_state_posteriors(float *x, float *state_prob)
{
	 if (!loaded)
		  throw ST_exception("GMM model not loaded.");
	 float denom;

    compute_state_prob_fast(x, state_prob,(float *) &denom);
	 for (int i=0; i<num_mix; i++) {
		  if (state_prob[i]==-1.0e6) {
				state_prob[i] = 0;
		  } else {
				state_prob[i] = (float) exp((double) (state_prob[i]-denom));
		  }
	 }
}

bool GMM_model::is_loaded () {
	return GMM_model::loaded;
}

void GMM_model::map_adapt_mean (REAL_EXP *sum, REAL_EXP *ec, float relevance_factor, REAL_EXP *mean_store) {
//
// Usage: sum, ec            expected counts and sum statistics (input)
//        relevance_factor 
//        mean               output adapted mean, if NULL replace GMM model mean with this mean
//                           if non-NULL, leave GMM alone and return adapted mean in "mean"
//
	 int i, j;
	 float alpha;

	 for (i=0; i<num_mix; i++) {
		  alpha = (float) (1/(ec[i]+relevance_factor));
		  for (j=i*num_fea; j<(i+1)*num_fea; j++) 
				mean_store[j] = relevance_factor*alpha*mean[j]+alpha*sum[j];
	 }
} 

struct element {
	 int index; 
	 float value;
};

int compare_elements (const void *ele1, const void *ele2)
{
	 if (((struct element *) ele1)->value < ((struct element *) ele2)->value) 
		  return -1;
	 else if (((struct element *) ele1)->value==((struct element *) ele2)->value)
		  return 0;
	 else
		  return 1;
}

float GMM_model::score (Features &f, vec<float> &frame_scores, vec<short> &frame_labels, int topM)
//
//  Usage: f              input features
//         frame_scores   per frame scores (optional)
//         frame_labels   labels of each frame (optional)
//
{
	 if (!loaded)
		 throw ST_exception("GMM_model::score -- Need to load GMM model first.");
	 if (f.num_outfeat()!=num_fea)
		 throw ST_exception("GMM_model::score -- GMM feature size does not match input features.");
	 if (f.num_vec()==0)
		 throw ST_exception("GMM_model::score -- No input features.");

    // Set up frame score save
    bool save_frame_scores = false;
    bool save_frame_labels = false;
    if (frame_scores.len>0) {
       if (frame_scores.len<f.num_vec()) 
          throw ST_exception("GMM_model::score -- frame_scores dimension is too small.");
       save_frame_scores = true;
    }

    if (frame_labels.len>0) {
       if (frame_labels.len != (topM*f.num_vec())) 
          throw ST_exception("GMM_model::score -- frame labels length dimension incorrect.");
       save_frame_labels = true;
    }

	 vec<float> x_feat(num_fea);
	 vec<float> state_prob(num_mix);
	 float prob;
    struct element *top = new struct element[num_mix];

	 //
	 // main loop
	 //
	 int i, j, k;
	 float s, total_score;
	 total_score = 0;

	 for (i=0; i<f.num_vec(); i++) {
		 x_feat = f.vector(i);

		 // Compute state probabilities
		 compute_state_prob_fast(x_feat.data, state_prob.data, (float *) &prob);

		 // Compute score 
		 s = -10000;
		 for (j=0; j<num_mix; j++)
			 s = linc_r(s, state_prob.data[j]);
		 total_score += s;
       
       if (save_frame_scores)
          frame_scores.data[i] = s;

       if (save_frame_labels) {
          // Sort scores 
          for (j=0; j<num_mix; j++) {
             top[j].index = j;
             top[j].value = state_prob.data[j];
          }
          qsort(top, num_mix, sizeof(struct element), compare_elements);

          // Find top M
          for (j=0, k=num_mix-1; j<topM; j++, k--)
             frame_labels.data[i*topM+j] = top[k].index;
       }

	 } // for

    delete[] top;
    total_score /= (float) f.num_vec();

    return total_score;

}

void GMM_model::score_models (Features &f, vector<GMM_model> &models, vec<float> &scores, float &ubm_score, 
   vec<float> &frame_scores, int topM, bool use_shortfall, float sf_delta)
//
//  Usage: f          input features
//         models     vector of GMM models to score
//         topM       number of top scoring mixtures to use in scoring, usually set to 5 and not easily accessible in LLSpeech
//         scores     output scores
//         ubm_score  UBM score
//  Assumption: GMM UBM is in the scoring object
//
{
	 if (!loaded)
		 throw ST_exception("GMM_model::score_models -- Need to load GMM UBM model first.");
	 if (f.num_outfeat()!=num_fea)
		 throw ST_exception("GMM_model::score_models -- GMM feature size does not match input features.");
	 if (f.num_vec()==0)
		 throw ST_exception("GMM_model::score_models -- No input features.");
	 if (scores.len != models.size())
		 throw ST_exception("GMM_model::score_models -- Output score vector dimension must be the same as the number of models.");
	 if (topM > num_mix)
		 topM = num_mix;

	 vector<GMM_model>::iterator it;
	 for (it=models.begin(); it<models.end(); it++) {
		 if (it->num_fea!=num_fea || it->num_mix!=num_mix)
			 throw ST_exception("GMM_model::score_models -- GMM model feature size or mixture size does not match UBM.");
	 }

    // Set up frame score save
    bool save_frame_scores = false;
    if (frame_scores.len>0) {
       if (frame_scores.len<(models.size()*f.num_vec())) 
          throw ST_exception("GMM_model::score_models -- frame_scores dimension is too small.");
       save_frame_scores = true;
    }

	 vec<float> x_feat(num_fea);
	 vec<float> state_prob(num_mix);
	 float prob;
	 struct element tmp;
    struct element *top = new struct element[num_mix];

	 //
	 // main loop
	 //
	 int i, j, k, l;
	 float s, sp, min;
	 ubm_score = 0;
	 for (i=0; i<models.size(); i++)
		 scores.data[i] = 0;

	 for (i=0; i<f.num_vec(); i++) {
		 x_feat = f.vector(i);

		 // Compute state probabilities
		 compute_state_prob_fast(x_feat.data, state_prob.data, (float *) &prob);

		 // Initialize topM probabilities with sorted topM elements from state_prob
		 // Increasing order -- top[0].value is smallest
		 for (j=0; j<num_mix; j++) {
			 top[j].index = j;
			 top[j].value = state_prob.data[j];
		 }
		 qsort(top, topM, sizeof(struct element), compare_elements);

		 // Funky algorithm in LLSpeech to maintain a topM list; refactored from select_index in vec_utils
		 // Not sorted, but top[0].value is the smallest
		 for (j=topM; j<num_mix; j++) {
			 if (top[j].value > top[0].value) {
				 top[0].index = j;
				 top[0].value = top[j].value;
				 for (k=1;;) {
					 l = k<<1;
					 if (l > topM) 
						 break;
					 if (l!=topM && (top[l-1].value > top[l].value)) 
						 l++;
					 if (top[k-1].value <= top[l-1].value)
						 break;
					 tmp.index = top[l-1].index;
					 tmp.value = top[l-1].value;
					 top[l-1].index = top[k-1].index;
					 top[l-1].value = top[k-1].value;
					 top[k-1].index = tmp.index;
					 top[k-1].value = tmp.value;
					 k = l;
				 }
			 }
		 }

		 // Compute score for UBM
		 s = -10000;
		 for (j=0; j<topM; j++)
			 s = linc_r(s, top[j].value);
		 ubm_score += s;

		 // Score tgt models using transcript from UBM
		 for (it=models.begin(), j=0; it<models.end(); it++, j++) {
			 s = -10000; 
			 min = (float) 1e30;
			 for (k=0; k<topM; k++) {
				 if (use_shortfall) 
					 compute_state_prob_shortfall(x_feat.data, *it, top[k].index, sp, min, sf_delta);
				 else
					 compute_state_prob(x_feat.data, *it, top[k].index, sp);
				 s = linc_r(s, sp);
			 }
          if (save_frame_scores)
             frame_scores.data[i*models.size()+j] = s;
			 scores.data[j] += s;
		 }

	 } // for

    delete[] top;

	 // Normalize scores by number of frames
	 for (i=0; i<scores.len; i++)
		 scores.data[i] /= (float) f.num_vec();
	 ubm_score /= (float) f.num_vec();

}

void GMM_model::suff_stats_ns (Features &f, vec<REAL_EXP> &sum, vec<REAL_EXP> &ec)
  //
  //  Usage: filename   input file name, raw floats
  //         sum        1st order suff stats
  //         ec         counts
  //
{
	 if (!loaded)
		 throw ST_exception("GMM_model::suff_stats_ns -- Need to load GMM model first.");

	 if (f.num_outfeat()!=num_fea)
		 throw ST_exception("GMM_model::suff_stats_ns -- GMM features size does not match input features.");
	 if (f.num_vec()==0)
		 throw ST_exception("GMM_model::suff_stats_ns -- No input features.");

	 if (sum.len!=(num_mix*num_fea) || ec.len!=num_mix) 
		 throw ST_exception("GMM_model::suff_stats_ns -- Output vectors are incorrectly dimensioned.");

	 int i, j, k;
	 vec<float> x_feat(num_fea);
	 vec<float> state_prob(num_mix);
	 for (i=0; i<(num_fea*num_mix); i++)
		  sum.data[i] = 0;
	 for (i=0; i<num_mix; i++)
		  ec.data[i] = 0;

	 //
	 // main loop
	 //
	 for (i=0; i<f.num_vec(); i++) {
		 x_feat = f.vector(i);

		 // Compute and store posteriors
		 compute_state_posteriors(x_feat.data, state_prob.data);

		 // Compute statistics
		 for (k=0; k<num_mix; k++) {
			 if (state_prob.data[k]!=0) {
				 ec.data[k] += state_prob.data[k];
				 for (j=0; j<num_fea; j++)
					 sum.data[k*num_fea+j] += state_prob.data[k]*x_feat.data[j];
			 }
		 }
	 } // for

}

#define LOOP_UNROLL_train
#define BLOCK_train 7
#define SHORTFALL_train 10
#define LOG2PI 1.837877
#define LN_MAXFLOAT 88.72283905

void GMM_model::compute_state_prob_fast (float *z, float *state_prob, float *prob)
{

	float *s_ptr, constant = 0;
	float *x, tmp, *mu, *ivar, *log_w, *log_d;
	int i, j;
	float tmp_prob;

#ifdef SHORTFALL_train
    float min_dist, delta = SHORTFALL_train;
#endif
    constant = (float) (num_fea*LOG2PI);
	 s_ptr = state_prob;
    mu = mean;
    ivar = inv_cov;
	 log_w = log_weight;
    log_d = log_det;

#ifdef SHORTFALL_train
    min_dist = (float) 1e30;
#endif
	 tmp_prob = -100000.0;
    for (i=0; i<num_mix; i++) {
		 x = z;
		 *s_ptr = *log_d + constant - 2*(*log_w++);

#ifdef LOOP_UNROLL_train
		 for (j=num_fea-1; j >= (BLOCK_train-1); j-=BLOCK_train) {
			// groups of BLOCK_train 
#if (BLOCK_train == 3)
			 tmp = (*x++ - *mu++); *s_ptr += tmp*tmp*(*ivar++);
			 tmp = (*x++ - *mu++); *s_ptr += tmp*tmp*(*ivar++);
			 tmp = (*x++ - *mu++); *s_ptr += tmp*tmp*(*ivar++);
#elif (BLOCK_train == 4)
			 tmp = (*x++ - *mu++); *s_ptr += tmp*tmp*(*ivar++);
			 tmp = (*x++ - *mu++); *s_ptr += tmp*tmp*(*ivar++);
			 tmp = (*x++ - *mu++); *s_ptr += tmp*tmp*(*ivar++);
			 tmp = (*x++ - *mu++); *s_ptr += tmp*tmp*(*ivar++);
#elif (BLOCK_train == 5)
			 tmp = (*x++ - *mu++); *s_ptr += tmp*tmp*(*ivar++);
			 tmp = (*x++ - *mu++); *s_ptr += tmp*tmp*(*ivar++);
			 tmp = (*x++ - *mu++); *s_ptr += tmp*tmp*(*ivar++);
			 tmp = (*x++ - *mu++); *s_ptr += tmp*tmp*(*ivar++);
			 tmp = (*x++ - *mu++); *s_ptr += tmp*tmp*(*ivar++);
#elif (BLOCK_train == 6)
			 tmp = (*x++ - *mu++); *s_ptr += tmp*tmp*(*ivar++);
			 tmp = (*x++ - *mu++); *s_ptr += tmp*tmp*(*ivar++);
			 tmp = (*x++ - *mu++); *s_ptr += tmp*tmp*(*ivar++);
			 tmp = (*x++ - *mu++); *s_ptr += tmp*tmp*(*ivar++);
			 tmp = (*x++ - *mu++); *s_ptr += tmp*tmp*(*ivar++);
			 tmp = (*x++ - *mu++); *s_ptr += tmp*tmp*(*ivar++);
#elif (BLOCK_train == 7)
			 tmp = (*x++ - *mu++); *s_ptr += tmp*tmp*(*ivar++);
			 tmp = (*x++ - *mu++); *s_ptr += tmp*tmp*(*ivar++);
			 tmp = (*x++ - *mu++); *s_ptr += tmp*tmp*(*ivar++);
			 tmp = (*x++ - *mu++); *s_ptr += tmp*tmp*(*ivar++);
			 tmp = (*x++ - *mu++); *s_ptr += tmp*tmp*(*ivar++);
			 tmp = (*x++ - *mu++); *s_ptr += tmp*tmp*(*ivar++);
			 tmp = (*x++ - *mu++); *s_ptr += tmp*tmp*(*ivar++);
#else
			 printf("Bad loop unrolling blocking factor. (%d) Range 3-7 currently supported\n",BLOCK_train);
			 exit(0);
#endif 

#ifdef SHORTFALL_train
			 if (*s_ptr > (min_dist+delta)) {
				 j -= (BLOCK_train-1);     
				 mu += j;
				 ivar += j;
				 *s_ptr = -1e6;
				 goto NEXT_COMP;
			 }
#endif /* ifdef SHORTFALL_train */
		 }  // for j

		 for (;j >= 0;--j)	{  // remainder
			 tmp = (*x++ - *mu++); 
			 *s_ptr += tmp*tmp*(*ivar++);
		 }
#else
		 for (j=num_fea+1; --j; ) { // dec then compare 
			 tmp = (*x++ - *mu++);
			 *s_ptr += tmp*tmp*(*ivar++);

#ifdef SHORTFALL_train
		 if (*s_ptr > (min_dist+delta)) {
			 --j;
			 mu += j;
			 ivar += j;
			 *s_ptr = -1e6;
			 goto NEXT_COMP;
		 }
#endif /* ifdef SHORTFALL_train */
		 }
#endif /* ifdef LOOP_UNROLL_train */

#ifdef SHORTFALL_train
		 if (*s_ptr < min_dist) 
			 min_dist = *s_ptr;
#endif /* ifdef SHORTFALL_train */
		 *s_ptr *= -0.5;
		 tmp_prob = linc_r(tmp_prob, *s_ptr);

NEXT_COMP:
		 ++log_d;
		 ++s_ptr;
	
	 } // loop i

	 *prob = tmp_prob;

} // compute_state_prob_fast

float GMM_model::linc_r (float x,float y)
{
  int linc_i;
  float linc_x, linc_y, linc_b, linc_d;

  linc_x = x;
  linc_y = y;
  if (linc_x > linc_y){
    linc_b = linc_x;
    linc_d = linc_x-linc_y;
  } else {
    linc_b = linc_y;
    linc_d = linc_y-linc_x;
  }
 
  if(linc_d > ladd.ladd_limit) 
	  x=linc_b;
  else {
    linc_i = (int)(linc_d*ladd.ladd_fac);
    x = linc_b+ladd.ladd_tbl[linc_i]+
		 (linc_d+ladd.ladd_idel[linc_i])*ladd.ladd_dtbl[linc_i];
  }
  return x;
}

void GMM_model::ladd_init() {
  float delta;
  int i;

  delta = -(float)RLA/(float)NLA;
  ladd.ladd_fac = (float)NLA/(float)RLA;
  for(i=0;i<=NLA+1;i++) 
	  ladd.ladd_idel[i]=i*delta;
  for(i=0;i<=NLA+1;i++) 
    ladd.ladd_tbl[i] = (float)
      log((double)1.0+exp((double)ladd.ladd_idel[i]));

  for (i=0; i<=NLA; i++) { // bias reduction
    ladd.ladd_dtbl[i]= (float) (
      (log((double)1.+exp((double)(i+.5)*delta))
	 - (ladd.ladd_tbl[i]+ladd.ladd_tbl[i+1])/2.)/4. );
  }
  ladd.ladd_tbl[0]+=ladd.ladd_dtbl[0];
  for(i=0;i<=NLA;i++) {
    ladd.ladd_tbl[i]+=ladd.ladd_dtbl[i];
    ladd.ladd_tbl[i+1]+=ladd.ladd_dtbl[i];
  }
  ladd.ladd_tbl[NLA+1]+=ladd.ladd_dtbl[NLA];
  
  for(i=0;i<=NLA;i++)
    ladd.ladd_dtbl[i]=
      (ladd.ladd_tbl[i+1]-ladd.ladd_tbl[i])*ladd.ladd_fac;
  ladd.ladd_limit = RLA;

}

void GMM_model::compute_state_prob_shortfall (float *x, GMM_model &model, int i, float &log_state_prob, float &min, float sf_delta)
// shortfall is a little bit disturbing -- it's dependent on the ordering of the topM values which in turn is dependent
// on the ordering algorithm -- the output of select_index is not sorted by mixture probability in LLSpeech
{
	int j;
	float val;
	float *mu, *ivar;
	
	mu = model.mean + i*model.num_fea;
	ivar = model.inv_cov + i*model.num_fea;
	val = (float) (log_det[i] + model.num_fea*LOG2PI-2*log_weight[i]);

	for (j=0; j<model.num_fea; j++) {
		val += (x[j]-mu[j])*(x[j]-mu[j])*ivar[j];
		if (val > min+sf_delta) 
			break;
	}
	min = (val < min) ? val : min;
	log_state_prob = (float) ((j < model.num_fea) ? -1e6 : -0.5*val);  // if it doesn't make it all the way through the loop, use the default value

} // compute_state_prob_shortfall

void GMM_model::compute_state_prob (float *x, GMM_model &model, int i, float &log_state_prob)
{
	int j;
	float val;
	float *mu, *ivar;
	
	mu = model.mean + i*model.num_fea;
	ivar = model.inv_cov + i*model.num_fea;
	val = (float) (log_det[i] + model.num_fea*LOG2PI-2*log_weight[i]);
	for (j=0; j<model.num_fea; j++)
		val += (x[j]-mu[j])*(x[j]-mu[j])*ivar[j];
	log_state_prob = (float) (-0.5*val);

} // compute_state_prob
  
void GMM_model::load (string model_file_name) {

	FILE *infile = fopen(model_file_name.c_str(), "rb");
	if (infile == NULL) 
		throw ST_exception(string("GMM_model::load_model -- unable to load GMM, ")+model_file_name);

	// Read header -- old style IO, ugly code -- probably need to update
	char nist_mark[11];
	char *d = fgets(nist_mark, 8, infile);
	if (d==NULL) {
		fclose(infile);
		throw ST_exception("GMM_model::load_model -- unable to read header for GMM model.");
	}
	int num_read = (int) strlen(nist_mark);
	if (num_read!=7 || strcmp(nist_mark,"NIST_1A")!=0) {
		fclose(infile);
		throw ST_exception("GMM_model::load_model -- unable to read header for GMM model.");
	}
	int num_bytes;
	if (fscanf(infile, " %d ", &num_bytes)!=1) {
		fclose(infile);
		throw ST_exception("GMM_model::load_model -- unable to read header for GMM model.");
	}
	rewind(infile);
	if (num_bytes<=0) {
		fclose(infile);
		throw ST_exception("GMM_model::load_model -- unable to read header for GMM model.");
	}
	char *hdr = new char[num_bytes+1];
	if (fread(hdr, sizeof(char), num_bytes, infile)!=num_bytes) {
		fclose(infile);
		delete[] hdr;
		throw ST_exception("GMM_model::load_model -- unable to read header for GMM model.");
	}
	hdr[num_bytes] = '\0';

	// Parse header 
	char *loc;
	loc = strstr(hdr, "sample_size -i");
	if (loc==NULL) {
		fclose(infile);
		delete[] hdr;
		throw ST_exception("GMM_model::load_model -- couldn't find sample size in header.");
	}
	if (sscanf(loc, "sample_size -i %d ", &num_fea)!=1) {
		fclose(infile);
		delete[] hdr;
		throw ST_exception("GMM_model::load_model -- couldn't find sample size in header.");
	}
	loc = strstr(hdr, "mixture_order -i");
	if (loc==NULL) {
		fclose(infile);
		delete[] hdr;
		throw ST_exception("GMM_model::load_model -- couldn't find mixture order in header.");
	}
	if (sscanf(loc, "mixture_order -i %d ", &num_mix)!=1) {
		fclose(infile);
		delete[] hdr;
		throw ST_exception("GMM_model::load_model -- couldn't find mixture order in header.");
	}
	if (num_fea<=0 || num_mix<=0) {
		fclose(infile);
		delete[] hdr;
		throw ST_exception("GMM_model::load_model -- header parameters don't make sense.");
	}
	delete[] hdr;

	// log mixture weights
	log_weight = new float[num_mix];
	if (fread(log_weight, sizeof(float), num_mix, infile)!=num_mix) {
		fclose(infile);
		delete[] log_weight;
		log_weight = 0;
		throw ST_exception("GMM_model::load_model -- unexpected end of file.");
	}

	// means
	mean = new float[num_fea*num_mix];
	if (fread(mean, sizeof(float), num_mix*num_fea, infile)!=(num_mix*num_fea)) {
		fclose(infile);
		delete[] log_weight;
		delete[] mean;
		log_weight = 0;
		mean = 0;
		throw ST_exception("GMM_model::load_model -- unexpected end of file.");
	}

	// log det
	log_det = new float[num_mix];
	if (fread(log_det, sizeof(float), num_mix, infile)!=num_mix) {
		fclose(infile);
		delete[] log_weight;
		delete[] mean;
		delete[] log_det;
		log_weight = 0;
		mean = 0;
		log_det = 0;
		throw ST_exception("GMM_model::load_model -- unexpected end of file.");
	}

	// inverse covariance
	inv_cov = new float[num_fea*num_mix];
	if (fread(inv_cov, sizeof(float), num_mix*num_fea, infile)!=(num_mix*num_fea)) {
		fclose(infile);
		delete[] log_weight;
		delete[] mean;
		delete[] inv_cov;
		delete[] log_det;
		inv_cov = 0;
		log_weight = 0;
		mean = 0;
		log_det = 0;
		throw ST_exception("GMM_model::load_model -- unexpected end of file.");
	}

	fclose(infile);
	loaded = true;

}

vec<float> GMM_model::ubm_mean (void) {
	vec<float> ubm_out(num_mix*num_fea);
	int i;
	for (i=0; i<(num_mix*num_fea); i++)
		ubm_out.data[i] = mean[i];
	return ubm_out;
}

vec<float> GMM_model::ubm_icov (void) {
	vec<float> ubm_out(num_mix*num_fea);
	int i;
	for (i=0; i<(num_mix*num_fea); i++)
		ubm_out.data[i] = inv_cov[i];
	return ubm_out;
}
