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
// Speech Tools header file
//
#pragma once

#include <stdint.h>
#include <list>
#include <map>
#include <set>
#include <string>
#include <vector>

using namespace std;

#include "st_exception.h"

// Types for features and fft
// Can be changed with a recompile
typedef float REAL_FEAT;
typedef float REAL_FFT;
typedef double REAL_EXP;
typedef float REAL_MODEL;

//
// Vector template
//
#include "vec.hpp"
typedef map<string, vec<REAL_MODEL> > modkey;
typedef map<string, vec<double> > zparamkey;

//
// Misc functions
//
unsigned file_len (const char *infile);

//
// General interval class
//
class Interval {
public:
   static const int SPEECH, NONSPEECH, UNK;
   int start, end;
   int label;
   float score;
   Interval (int start_in, int end_in, int label=Interval::UNK, float score=0.0);
};

class Time_Interval {
public:
   float start, end;
   int label;
   float score;
   Time_Interval (float start_in, float end_in, int label, float score=0.0);
};

//
// Signal class
//
class GMM_model;
class Signal {
 	int	num_samples; 
 	float	sampling_frequency; 
   bool  resampler_loaded;
   map<string, vec<double> > resampler_filters;
public:
 	float	*value; 
 	Signal (const Signal &x); 
   Signal (int sampling_frequency);
	Signal ();
 	Signal (const short *data, int num_samples, int sampling_frequency); 
 	Signal (const double *data, int num_samples, int sampling_frequency); 
 	~Signal(); 
 	void	comp_ifir_filter (const float *h1, int n, int filter_length_1, const float *h2, int filter_length_2); 
 	void  fir_filter (const float *h, int filter_length); 
 	void	fir_n_filter (const float *h, int filter_length, int n); 
   void  fir_updn (const double *h, const int filter_len, int L, int M);
   void  fir_zphase (const double *h, const int filter_len);
   vec<float> get_f0 (float min_f0, float max_f0, float window_dur, float &frame_step);
   list<Time_Interval> gmmsad (Signal &x, string &feat_config, vector<GMM_model> &gmm_models, int label_window, int gmmsad_keep,
      float min_gap_dur, float min_seg_dur);
   void  init_resampler (string filter_directory);
 	void	load_raw_float (string filename, int sampling_frequency);
   void  load_raw_short (string filename, int sampling_frequency);
	void  load_sph (string filename, int channel_num);
	void  load_pcm_wav (string filename, bool sum_channels=false);
	int   len();
   void  normalize (void);
 	void  preemphasis (float pre_coeff); 
 	void	remove_mean (void); 
   void  resample_8k (void);
   void  resample_16k (void);
	float sampling_freq (void);
   void  save_pcm_wav (string filename, bool scale);
 	void  save_raw_short (string filename, bool clip, bool scale);  // clip = clip at 32767/-32768
	                                                                // scale = scale so max abs() is no greater than 32767
 	void  save_raw_float (string filename); 
private:
	Signal& operator= (const Signal& src);
}; 

class Features {
protected:
	REAL_FEAT *feat;
	REAL_FEAT *energy;

	int  num_features;
	int  num_base_features;
	int  num_vectors;

	int  delta_offset;
	int  accel_offset;
	int  sdc_offset;
	int  sdc_k;

	bool features_available;
	bool energy_available;
	bool delta_available;
	bool accel_available;
	bool sdc_available;

	bool sad_available;
	short *sad;
	int sad_len;
	string outfeat;
	virtual void config (Signal &x, string &config);

public:
   Features ();
	~Features();

	// Specific to inherited classes
	virtual void process (Signal &x, string &config_str);

	// General purpose manipulation/access
	void load_raw (string filename, int num_feat);  // load from raw floats
	int num_base_feat (void);
	int num_outfeat (void);
	int num_total_feat (void);
	int num_vec (void);
	void save_raw (string filename); // save to raw floats
	void set_outfeat (string outfeat); // options are : "all" or arbitrary combinations of "feads"
	vec<REAL_FEAT> vector(int i);  // Get the ith feature vector

	// General purpose transformation
	void accel (int accel_spread);
	void delta (int delta_spread);
	void delta2point (int delta_spread);
	void feat_norm();  // Window is the entire feature set
	void rasta();
	void sdc (int sdc_p, int sdc_k_in);

	// SAD related, sloppy = ignore truncation errors in the marks
   void and_user_marks (const char *user_marks, int blocksz_ms, float in_win_inc_ms);
	void apply_SAD (void);
	void load_SAD_label (string filename, bool sloppy);    // 0/1 label file
   void load_SAD_marks (list<Time_Interval> &seg_sad, float win_inc);
	void load_SAD_marks (string filename, float win_inc);  
   string load_SAD_marks_string (string filename, float win_inc);
   int  SAD (int i); // Get ith SAD label (0/1)
	vec<uint8_t> SAD_labels(void);
   bool SAD_available(void);
	void save_SAD_label (string filename);
   void save_SAD_marks (string filename, float in_win_inc_ms);
	void xtalk (double abs_min_energy, double thresh, int med_len=1);
   bool no_speech(void);
	// TODO:
	// fNAP
private:
	Features& operator= (const Features& src); 
	Features(const Features &in);
};

class MFCC_Features : public Features {
	// General config params
	float alpha;
	float win_len_ms, win_inc_ms;
	float fb_low, fb_hi;
	float dither;
	float sampling_freq;
	bool fb_only;

	// FB params
	float **filt;
	int *filt_beg;
	int *filt_end;
	float *fc;
	int num_filt; // total number of filters
	int num_filt_bl; // actual number of filters after bandlimiting

	// Frame params
	int n_fft, n2_fft;
	int win_len, win_inc;
	vec<REAL_FFT> *preemp_filt;

	// Cep params
	float *icos_matrix;
	int num_cep;
	bool keep_c0;

	// Internal methods
	int fb_index_low, fb_index_hi;
	vec<REAL_FFT> calc_fb(vec<REAL_FFT> &in);
	vec<REAL_FEAT> icostrans (vec<REAL_FFT> &x);
	void config (Signal &x, string &config);
	void init_filtbank(float fmax, const int num_dft, bool linear=false, int num_filt=-1);
	void init_icostrans(void);
	int filtbank_get_filt_num (float freq);
public:
   void process(Signal &x, string &config_str);
   float get_win_inc_ms(void);
	MFCC_Features ();
	~MFCC_Features ();
   float duration(void);
};

class T2_BIC_params {
public:
   int init_search_interval;
   int search_inc;
   int search_offset;
   int max_seg_size;
   int min_seg_size;
   float bic_weight;
   float bic_thresh;
   bool verbose;
   T2_BIC_params(void);
};

class Cluster_BIC_params {
public:
   float bic_weight;
   float bic_thresh;
   int min_num_clusters;
   int max_num_clusters;
   bool verbose;
   Cluster_BIC_params (void);
};

class Segmentation {
   vec<float> bic_stats;
   void segment_T2_BIC_one_interval (T2_BIC_params &params);
   void segment_T2_BIC_multiple_intervals (T2_BIC_params &params, bool speech_only);
public:
   static const int SPEECH, NONSPEECH;
   static const int MALE, FEMALE, NB, WB, UNK;
   int speech;
   int gender; 
   int bw;
   int speaker;
   list<Interval> intervals;
   Features &feats;
   Segmentation(Features &feats_in);
   Segmentation (Features &feats_in, int start, int end);
   vec<float> &get_bic_stats(void);
   void intervals_from_SAD (void);
   void print(void);
   void save_marks (string &filename, float win_inc_ms, bool audacity_format=false);
   void segment_T2_BIC (T2_BIC_params &params, bool speech_only);
private:
   Segmentation& operator= (const Segmentation& src);
   // Segmentation(const Segmentation &in);
};

//
// Classes for clustering a segmentation
//
class Segmentation_BIC_info {
public:
   float n;  // number of vectors in segmentation
   float dim; // dimension of vectors
   float log_det; // log determinant of covariance
   vec<float> s1;  // first order stats
   vec<float> s2;  // second order stats
   void compute_log_det (void);
   Segmentation_BIC_info(Segmentation &seg);
};

class Cluster_segmentation {
   vec<float> dist;
   list<Segmentation> seg_list;
   list<Segmentation_BIC_info> bic_info;
   float delta_bic (Segmentation_BIC_info &b1, Segmentation_BIC_info &b2);
   void find_min_dist (int &i, int &j, float &dist_min);
   void merge (int i1, int j1);
   Cluster_BIC_params params;
public:
   Cluster_segmentation(Segmentation &seg_init);
   void bic_cluster (void);
   void bic_init (Cluster_BIC_params &params_in);
   void save_marks (string &filename, float win_inc_ms, bool audacity_format=false);
   string to_string (float win_inc_ms);
};
float log_det_matrix (vec<float> &A, float lambda);

//
// General complex FFT, radix-2, DIF
// Real FFT
//
class FFT_DIF {
   int log2n;
   int n;
   int *bit_reverse;
   REAL_FFT *sintbl;
   REAL_FFT *costbl;
public:
   FFT_DIF (int log2n);
   ~FFT_DIF();
   void fft (vec<REAL_FFT> &in, vec<REAL_FFT> &out_r);
 private:
	FFT_DIF& operator= (const FFT_DIF& src);
	FFT_DIF(const FFT_DIF &in);
};

class RFFT {
	int n, n2;
	int log2n;
   REAL_FFT *costbl; 
   REAL_FFT *sintbl;  
	FFT_DIF f;
public:
	RFFT (int n);
	~RFFT();
	void rfft(vec<REAL_FFT> &in, vec<REAL_FFT> &out_r, vec<REAL_FFT> &out_i);
private:
	RFFT& operator= (const RFFT& src);
	RFFT(const RFFT &in);
};

class mag_fft {
	int n_input;
	int n_fft;
	int n2_fft;
	int log2n;
	RFFT *r;
	vec<REAL_FFT> *in, *out_r, *out_i;
public:
	mag_fft(int n);
	~mag_fft();
	vec<REAL_FFT> calc(vec<REAL_FEAT> &x);
	int len();
private:
	mag_fft& operator= (const mag_fft& src);
	mag_fft(const mag_fft &in);
};

//
// Class for frame oriented processing
//
class Frame {
	vec<REAL_FEAT>  x;
	vec<REAL_FEAT>  hw;
	int frame_inc;
	bool frame_loaded;
	int32_t s1, s2;  // seeds for dither function, should be 32-bits
public:
	Frame (int win_len, int win_inc);
	~Frame ();
	void get_frame (Signal &x, int i);
	int get_num_frames (Signal &x);
	void load (float *x);
	void dither(float scale);
	void rm_dc();
	void hamming_window();
	REAL_FEAT energy();
	vec<REAL_FEAT>& vector();
private:
	Frame& operator= (const Frame& src); 
	Frame(const Frame &in); 
};

//
// GMM class stuff
//
#define RLA 18     // ladd_init, difference range
#define NLA 2500   // ladd_init, table size 
struct ladd_struct {
  float ladd_limit;
  float ladd_fac;
  float ladd_tbl[NLA+3];
  float ladd_dtbl[NLA+3];
  float ladd_idel[NLA+3];
};

class GMM_model {
	ladd_struct ladd;
	float linc_r (float x,float y);
	void ladd_init();
	void compute_state_prob (float *x, GMM_model &model, int i, float &log_state_prob);
	void compute_state_prob_fast (float *z, float *state_prob, float *prob);
	void compute_state_prob_shortfall (float *x, GMM_model &model, int i, float &log_state_prob, float &min, float sf_delta);
	bool loaded;
public:
	float *mean, *inv_cov, *log_weight, *log_det;
	void load (string model_file_name);
	int num_fea, num_mix;
	GMM_model (void);
	GMM_model(const GMM_model &in);
	~GMM_model (void);
	GMM_model& operator= (const GMM_model& src); 
	bool is_loaded(void);
	void compute_state_posteriors(float *x, float *state_prob);
	void map_adapt_mean (REAL_EXP *sum, REAL_EXP *ec, float relevance_factor, REAL_EXP *mean_store);
   float score (Features &f, vec<float> &frame_scores, vec<short> &frame_labels, int topM=5);
	void score_models (Features &f, vector<GMM_model> &models, vec<float> &scores, float &ubm_score, 
                      vec<float> &frame_scores, int topM=5, bool use_shortfall=true, float sf_delta=10.0);
	void suff_stats_ns (Features &f, vec<REAL_EXP> &sum, vec<REAL_EXP> &ec);
	vec<float> ubm_mean (void);
	vec<float> ubm_icov (void);
private:
};

// 
// IPDF expansion
//
typedef map<string, vec<REAL_EXP> > veckey;
class IPDF_expansion {
	int max_corank;      // Assume max corank the same across all keys for now
	GMM_model exp_gmm;
	int exp_dim;
	veckey nap_proj_map;
	void nap (vec<REAL_EXP> &x_exp, string key, int corank);
	void scale_fixed_to_vm (vec<REAL_EXP> &x_exp, vec<REAL_EXP> &ec);
public:
  IPDF_expansion (void);
  ~IPDF_expansion (void);
  vec<REAL_EXP> expansion (Features &f, float rf);
  vec<REAL_EXP> expansion_with_nap (Features &f, float rf, string nap_key, int corank);
  int dim (void);
  void load_gmm_model (string model_file_name);
  void load_nap_projection (string projection_file_name, string key);
  int num_fea (void);
  int num_mix (void);
  IPDF_expansion& operator= (const IPDF_expansion &src);
  IPDF_expansion (const IPDF_expansion &src);
};
typedef map<string, IPDF_expansion > ipdfkey;

//
// IP_score -- Basic Z-, ZT-norm scoring for inner products
// 
class IP_score {
	int exp_dim;
	modkey tnorm_map;
	modkey znorm_map;
	zparamkey zparams_mean_map;
	zparamkey zparams_std_map;
	void load_models (string filename, string key, modkey &model_map);
public:
	IP_score (IPDF_expansion &x);
	~IP_score (void);
	void fillin_nan (double *scores, int m1, int m2, double *ipdf_vec, bool same);
	int num_tnorm_models (string key);
	void load_tnorm_models (string file, string key);
	void load_znorm_models (string file, string key);
	void load_znorm_params (string file, string key);
	double sztnorm_score (vec<REAL_EXP> &x1_exp, vec<REAL_EXP> &x2_exp, string key);
	void sztnorm_score_matrix (double *scores, int m1, int m2, double *sztnorm_params);
	void znorm_params (REAL_EXP *x1_data, int x1_len, string key, double &m_z, double &s_z);
	void znorm_params_matrix (double *x1_data, int x1_len, int num_vec, string key, double *param_store);
	void ztnorm_params (REAL_EXP *x2_data, int x2_len, string key, double &m_z, double &s_z, double &m_zt, double &s_zt);
	void ztnorm_params_matrix (double *x1_data, int x1_len, int num_vec, string key, double *param_store);
	double ztnorm_score (vec<REAL_EXP> &x1_exp, vec<REAL_EXP> &x2_exp, string key);
	void ztnorm_score_matrix (double *scores, int m1, int m2, double *ztnorm_params);
private:
	IP_score& operator= (const IP_score& src);
 	IP_score (const IP_score &src); 
};

