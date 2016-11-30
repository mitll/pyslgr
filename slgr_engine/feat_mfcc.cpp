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

// MFCCs

// Refactored from LLSpeech 2.2.10 code & BC Speech Tools:
// filtbank.cc, mfb_to_cep.cc, mfb_to_cep_subs.cc, mfb_utils.cc, pcm_to_feat.cc, pcm_to_mfb_subs.cc, vec_utils.cc

// Written by BC 10/18/10
#include "stdafx.h"
#include <stdio.h>
#include <math.h>
#include "speech_tools.h"
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
using boost::property_tree::ptree;
using boost::property_tree::read_json;

#define PI 3.141592653589793

typedef float* flt_ptr;
static float warp (float freq, float alpha, float sfreq);

//
// Main processing routine
//
void MFCC_Features::process (Signal &x, string &config_str) 
{

	// MFCC configuration
	config(x, config_str);
	
	// Set up storage
	Frame f(win_len, win_inc);
	mag_fft m(win_len);
	vec<REAL_FEAT> x_frame(win_len);
	vec<REAL_FFT> x_mag(n2_fft+1);
	vec<REAL_FFT> x_fb(num_filt_bl);
	vec<REAL_FEAT> x_cep(num_features);
	REAL_FEAT *feat_ptr;

	// Allocate space for features
	num_vectors = f.get_num_frames(x);
	feat_ptr = feat = new REAL_FEAT[num_vectors*num_features];
	energy = new REAL_FEAT[num_vectors];

	// Loop across frames
	int i, j;
	for (i=0; i<num_vectors; i++) {

		// Frame based pre-processing
		f.get_frame(x, i);
		f.dither(dither);
		f.rm_dc();
		f.hamming_window();
		energy[i] = f.energy();

		// FFT magnitudes -- FFT not quite the same as LLSpeech
		x_frame = f.vector();
		x_mag = m.calc(x_frame);

		// Pre-emphasis
		x_mag.scale(*preemp_filt);

		// Calculate filterbanks with bandlimiting
		x_fb = calc_fb(x_mag);

		if (fb_only) {
			for (j=0; j<num_features; j++, feat_ptr++)
				*feat_ptr = x_fb.data[j];
		} else {
			// Cepstral coefficients
			x_cep = icostrans(x_fb);
			for (j=0; j<num_features; j++, feat_ptr++)
				*feat_ptr = x_cep.data[j];
		}

	}

	features_available = true;
	energy_available = true;

}

MFCC_Features::MFCC_Features () 
{
	filt_beg = 0;
	filt_end = 0;
	filt = 0;
	fc = 0;
	preemp_filt = 0;
	icos_matrix = 0;
}

MFCC_Features::~MFCC_Features () 
{
	delete[] filt_beg;
	delete[] filt_end;
	delete[] fc;
	if (filt!=0) {
		for (int i=0; i<num_filt; i++)
			delete[] (filt[i]);
	}
	delete[] filt;
	delete preemp_filt;
	delete[] icos_matrix;
}

void MFCC_Features::config (Signal &x, string &config) 
{
   bool linear = false;
   int tgt_num_filt = -1;

   istringstream is(config);
   ptree pt;

   try {
      read_json(is, pt);
      alpha = pt.get<float>("alpha");
      dither = pt.get<float>("dither");
      fb_low = pt.get<float>("fb_low");
      fb_hi = pt.get<float>("fb_hi");
		keep_c0 = pt.get<bool>("keep_c0");
      linear = pt.get<bool>("linear");
      num_cep = pt.get<int>("num_cep");
      win_inc_ms = pt.get<float>("win_inc_ms");
      win_len_ms = pt.get<float>("win_len_ms");
		if (pt.find("fb_only")==pt.not_found())
			fb_only = false;
		else
			fb_only = pt.get<bool>("fb_only");
      if (pt.find("tgt_num_filt")!=pt.not_found())
         tgt_num_filt = pt.get<int>("tgt_num_filt");
   } catch (exception &e) {
      throw ST_exception(string("MFCC_Features::config:  error in config string, ")+e.what());
   }

	// Sanity check
	float f_max = x.sampling_freq()/2.0f;
	if ((fb_low<0) || (fb_hi<fb_low) || (fb_hi>f_max))
      throw ST_exception("MFCC_Features::process -- filterbank parameters don't make sense.");
	sampling_freq = x.sampling_freq();

	// Find window parameters
	win_inc = (int) ((x.sampling_freq()*win_inc_ms)/1000);
	win_len = (int) ((x.sampling_freq()*win_len_ms)/1000);
	if ((win_len % 2) != 0)
		win_len++;
	if (win_inc<=0 || win_len<=0)
		throw ST_exception("MFCC_Features::process -- Window increment or length is <= 0.");

	// Initialize filterbank
	for (n_fft=1; n_fft<win_len; n_fft*= 2);
	n2_fft = n_fft/2;
	init_filtbank(f_max, n2_fft+1, linear, tgt_num_filt);
	if (fb_low==0 && fb_hi==f_max) { // Use full band
		fb_index_low = 0;
		fb_index_hi = num_filt-1;
	} else {
		fb_index_low = filtbank_get_filt_num(fb_low);
		fb_index_hi = filtbank_get_filt_num(fb_hi);
	}
   num_filt_bl = fb_index_hi-fb_index_low+1;

	// Initialize cepstral processing
	if ((num_cep <= 0) || (num_cep>(num_filt_bl-1)))
		throw ST_exception("Number of cepstral coefficients larger than number of filters.");

	num_features = 0;
	if (fb_only) {
		num_features = num_filt_bl;
	} else {
		if (keep_c0)
			num_features = 1;
		num_features += num_cep;
	}
	num_base_features = num_features;
	init_icostrans();

	// Initialize pre-emphasis
	preemp_filt = new vec<REAL_FFT> (n2_fft+1);
	float finc = (float) (f_max/(n2_fft+1.0));
	float f = 0;
	for (int j=0; j<=n2_fft; f += finc, j++)
		preemp_filt->data[j] = (REAL_FEAT) (1.0 + f*f/2.5e5);

}

vec<REAL_FFT> MFCC_Features::calc_fb (vec<REAL_FFT> &mag)
{
	int i, j, k, l;
	vec<REAL_FFT> x_fb(num_filt_bl);
	REAL_FFT val;

	for (i=fb_index_low, l=0; i<=fb_index_hi; i++, l++) {
		val = 0;
		for (j=filt_beg[i], k=0; j<filt_end[i]; j++, k++)
			val += mag.data[j]*filt[i][k];
		val = (REAL_FFT) (10.0*log10(val + 1e-20)-40.0);
		x_fb.data[l] = val;
	}
	return x_fb;

}

int MFCC_Features::filtbank_get_filt_num (float freq) 
{
  int i;
  for (i=1; i<=num_filt && fc[i]<=freq; i++);
  i--;
  if(i < 1) 
	  return 0;
  if (i==num_filt) 
	  return num_filt-1;
  if ((fc[i+1]-freq) < (freq-fc[i])) 
	  i++; // move to closest center freq
  return i-1;

}

vec<REAL_FEAT> MFCC_Features::icostrans (vec<REAL_FFT> &x)
{
	vec<REAL_FEAT> out(num_features);
	int i, j, k;
	REAL_FFT val;

	// Matrix vector multiply -- should use BLAS
	i = keep_c0 ? 0 : 1;
	for (k=0; k<num_features; k++, i++) {
		val = 0;
		for (j=0; j<num_filt_bl; j++)
			val += icos_matrix[i*num_filt_bl+j]*x.data[j];
		out.data[k] = (REAL_FEAT) val;
	}

	return out;

} 

void MFCC_Features::init_icostrans () 
{
	float tmp, W;
	int i, j;

	tmp = (float) (1.0/num_filt_bl);
	W = (float) (PI*tmp);
	icos_matrix = new REAL_FFT[(num_cep+1)*num_filt_bl];

	for (i=0; i < (num_cep+1); i++) {
		for (j=0; j < num_filt_bl; j++) {
			icos_matrix[i*num_filt_bl+j] = (float) (cos((double) i*(j+0.5)*W)*tmp);
		}
	}

}

void MFCC_Features::init_filtbank (const float fmax, const int num_dft, bool linear, int tgt_num_filt)
{
  float f, *Fptr, area;
  float finc, delta_f;
  int i, j, nf, init_len, final_len;

  finc = fmax/(num_dft-1);

  if (linear) {
     if (tgt_num_filt<=0)
        delta_f = 100;
     else
        delta_f = fmax/(tgt_num_filt+1);
     for (nf=0, f=-delta_f; f<fmax; nf++)
        f += delta_f;
     fc = new float[nf];
     for (nf=0, f=-delta_f; f<fmax; nf++)
        fc[nf] = warp(f+=delta_f, alpha, fmax);
     nf -= 2;
     if ((fc[nf]+fc[nf+1])/2 > fmax)
        nf--;
     init_len = (int) ((fc[2]-fc[0])/finc + 1);
  } else {
     // Determine length of fc
     for (nf=11, f=1000; f<fmax; nf++)
        f *= 1.1f;
     fc = new float[nf];

     // Center frequencies
     for (i=0; i<=10; i++) // 100 Hz
        fc[i] = (float) warp((float) i*100, alpha, fmax);
     for (nf=11, f=1000; f<fmax; nf++) { // 10% spacing 
        f *= 1.1f;
        fc[nf] = warp(f, alpha, fmax);
     }
     nf -= 2;
     if ((fc[nf]+fc[nf+1])/2 > fmax) 
        nf--;

     // Determine maximum initial length
     init_len = 0;
     for (i = 1; i < nf + 1; i++) {
        int x = (int) ((fc[i+1]-fc[i-1])/finc + 1);
        init_len = x > init_len ? x : init_len;
     }
  }

  filt_beg = new int[nf];
  filt_end = new int[nf];
  filt = new flt_ptr[nf];
  float *tmp_filt = new float[init_len];

  for (i=1; i<=nf; ++i) {
	  area = 0.0;
	  Fptr = tmp_filt;
	  filt_beg[i-1] = 0;
	  filt_end[i-1] = 0;
	  for (j=0, f=0.0; (f<fmax) && (j<num_dft); f+=finc, j++) {
		  if ((f <= fc[i-1]) || (f >= fc[i+1])) 
			  continue;
		  if (f < fc[i]) {
			  if (filt_beg[i-1] == 0) 
				  filt_beg[i-1] = j;
			  *Fptr = (f-fc[i-1])/(fc[i]-fc[i-1]);
			  area += (*Fptr++);
		  } else {
			  *Fptr = (fc[i+1]-f)/(fc[i+1]-fc[i]);
			  filt_end[i-1] = j+1;
			  area += (*Fptr++);
		  } 
	  } 

	  //  Single DFT sample in filter passband, happens with very small windows (e.g., 5ms)
	  if ((filt_end[i-1]==0) && (filt_beg[i-1]!=0)) 
		  filt_end[i-1] = filt_beg[i-1]+1;
	  else if ((filt_end[i-1]!=0) && (filt_beg[i-1]==0)) 
		  filt_beg[i-1] = filt_end[i-1]-1;
	  
	  final_len = filt_end[i-1] - filt_beg[i-1];
	  filt[i-1] = new float[final_len]; 
	  for (j=0; j<final_len; j++)
		  filt[i-1][j] = tmp_filt[j]/area;

  } // next i

  delete[] tmp_filt;
  num_filt = nf;

} // init_filterbank

float MFCC_Features::duration(void)
{
   return (float) (((float) num_vectors)*(0.001*win_inc_ms));
}

float MFCC_Features::get_win_inc_ms (void)
{
   return win_inc_ms;
}

//
// Local functions
//
static float warp (float freq, float alpha, float sfreq)
{
  float ret;
  float fr = (float) (freq * PI / sfreq);
  ret = (float) (fr + 2.0f * atan((double)(1.0 - alpha) * sin((double)fr) / (1.0 - (1.0 - alpha) * cos((double)fr))));
  ret *= (float) (sfreq/PI);
  return ret;
}


