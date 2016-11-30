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
// Speech Tools (ST) library
// 
// General feature functions
//

// 1/03, wmc, LL version
// 09/10, wmc, LL version 2, deltas, acceleration, feat norm, rasta
// 12/12, wmc, added xtalk refactored from xtalk.cc, xtalk_subs.cc
// 01/13, wmc, added delta2point refactored from delta_feat_subs.cc
// 01/13, wmc, added SDC refactored from sdc.cc, sdc_subs.cc

//
// Global includes.
//
#include "stdafx.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

//
// Local includes.
//
#include "speech_tools.h"
using namespace std;

//
// Local prototypes.
//
static REAL_FEAT dcep_constant (int cep_history_length);

//
// Constants
//
using namespace std;

//
// Features methods
//
Features::Features (void)  // Construct empty features object
{
	feat = 0;
	energy = 0;

   num_features = 0;
	num_base_features = 0;
	num_vectors = 0;

	delta_offset = 0;
	accel_offset = 0;
	sdc_offset = 0;
	sdc_k = 0;
	
	features_available = false;
	energy_available = false;
	delta_available = false;
	accel_available = false;
	sdc_available = false;
	
	sad_available = false;
	sad = 0;
	sad_len = 0;

	outfeat = "all";

}

Features::~Features (void)
{
	delete[] feat;
	delete[] energy;
	delete[] sad;

} // Features::~Features

void Features::accel (int accel_spread) 
{
	// Very similar to delta, but kept as a separate routine in case we want to change it later
	if (accel_spread < 1)
		throw ST_exception("Features::accel -- Accel spread must be >= 1.");
	if (!features_available || num_vectors==0)
		throw ST_exception("Features::accel -- No features available for accel computation.");
	if (!delta_available)
		throw ST_exception("Features::accel -- Delta features must be computed first.");
	if (accel_available)
		return;

	// Initialize, allocate space, copy old features
	int i, j, k;
	int new_num_feat = num_features+num_base_features;
	REAL_FEAT dcep_c = dcep_constant(accel_spread);
	REAL_FEAT *feat_new = new REAL_FEAT[new_num_feat*num_vectors];
	for (i=0; i<num_vectors; i++) {
		for (j=0; j<num_features; j++) 
			feat_new[i*new_num_feat+j] = feat[i*num_features+j];
	}
	delete[] feat;
	feat = feat_new;

	// Use first vector and last vector for boundary conditions
	REAL_FEAT *feat_ptr = feat+delta_offset;
	REAL_FEAT *delta_ptr = feat+num_features;
	for (i=0; i<num_vectors; i++) {
		for (j=0; j<num_base_features; j++) {  
			delta_ptr[j] = 0;
			for (k=-accel_spread; k<=accel_spread; k++) {
				if (i+k < 0)
					delta_ptr[j] += k*feat_new[delta_offset+j];
				else if (i+k >= num_vectors) 
					delta_ptr[j] += k*feat_new[(num_vectors-1)*new_num_feat+delta_offset+j];
				else
					delta_ptr[j] += k*feat_ptr[k*new_num_feat+j];
			}
			delta_ptr[j] *= dcep_c;
		}
		feat_ptr += new_num_feat;
		delta_ptr += new_num_feat;
	}

	// Update state
	accel_available = true;
	accel_offset = num_features;
	num_features = new_num_feat;

}

void Features::apply_SAD ()
{
	if (!features_available)
		throw ST_exception("Features::apply_SAD -- no features available.");
	if (!sad_available)
		throw ST_exception("Features::apply_SAD -- no SAD available.");
	if (sad_len != num_vectors) 
		throw ST_exception("Features::apply_SAD -- SAD number of frames does not match number of feature vectors.");
	int i, j, k, num_speech_frames = 0;
	for (i=0; i<num_vectors; i++)
		num_speech_frames += (int) sad[i];
	if (num_speech_frames == 0) 
		throw ST_exception("Features::apply_SAD -- No speech found.");

	// Apply SAD to feature vectors
	REAL_FEAT *feat_new = new REAL_FEAT[num_speech_frames*num_features];
	for (i=0, j=0; i<num_vectors; i++) {
		if (sad[i] == 1) {
			for (k=0; k<num_features; k++)
				feat_new[j*num_features+k] = feat[i*num_features+k];
			j++;
		}
	}
	num_vectors = num_speech_frames;
	delete[] feat;
	feat = feat_new;
	
	// Optionally apply SAD to energy if it's available
	if (energy_available) {
		REAL_FEAT *energy_new = new REAL_FEAT[num_speech_frames];
		for (i=0, j=0; i<num_vectors; i++) {
			if (sad[i] == 1) {
				energy_new[j] = energy[i];
				j++;
			}
		}
		delete[] energy;
		energy = energy_new;
	}

	// Change SAD to all 1's
	short *sad_new = new short[num_vectors];
	for (i=0; i<num_vectors; i++)
		sad_new[i] = 1;
	delete[] sad;
	sad = sad_new;
}

void Features::config (Signal &x, string &config) 
{
	// Doesn't do anything for the base class
}

void Features::delta (int delta_spread)
{
	if (delta_spread < 1)
		throw ST_exception("Features::delta -- Delta spread must be >= 1.");
	if (!features_available || num_vectors==0)
		throw ST_exception("Features::delta -- no features available for delta computation.");
	if (delta_available)
		return;

	// Initialize, allocate space, copy old features
	int i, j, k;
	int new_num_feat = num_features+num_base_features;
	REAL_FEAT dcep_c = dcep_constant(delta_spread);
	REAL_FEAT *feat_new = new REAL_FEAT[new_num_feat*num_vectors];
	for (i=0; i<num_vectors; i++) {
		for (j=0; j<num_features; j++) 
			feat_new[i*new_num_feat+j] = feat[i*num_features+j];
	}
	delete[] feat;
	feat = feat_new;

	// Use first vector and last vector for boundary conditions
	REAL_FEAT *feat_ptr = feat;
	REAL_FEAT *delta_ptr = feat+num_features;
	for (i=0; i<num_vectors; i++) {
		for (j=0; j<num_base_features; j++) {  // Base features are always 0 .. (num_base_feat-1)
			delta_ptr[j] = 0;
			for (k=-delta_spread; k<=delta_spread; k++) {
				if (i+k < 0)
					delta_ptr[j] += k*feat_new[j];
				else if (i+k >= num_vectors)
					delta_ptr[j] += k*feat_new[(num_vectors-1)*new_num_feat+j];
				else 
					delta_ptr[j] += k*feat_ptr[k*new_num_feat+j];
			}
			delta_ptr[j] *= dcep_c;
		}
		feat_ptr += new_num_feat;
		delta_ptr += new_num_feat;
	}

	// Update state
	delta_available = true;
	delta_offset = num_features;
	num_features = new_num_feat;

} // delta


void Features::delta2point (int delta_spread) 
{
	if (delta_spread < 1)
		throw ST_exception("Features::delta2point -- Delta parameter, delta_spread, must be >= 1.");
	if (!features_available || num_vectors==0)
		throw ST_exception("Features::delta2point -- no features available for sdc computation.");
	if (delta_available)
		return;

	// Initialize, allocate space, copy old features
	int i, j;
	int new_num_feat = num_features+num_base_features;

	REAL_FEAT *feat_new = new REAL_FEAT[new_num_feat*num_vectors];
	for (i=0; i<num_vectors; i++) {
		for (j=0; j<num_features; j++) 
			feat_new[i*new_num_feat+j] = feat[i*num_features+j];
	}
	delete[] feat;
	feat = feat_new;

	// Two point deltas, dc(t) = c(t+d)-c(t-d)
	REAL_FEAT *feat_ptr = feat;
	REAL_FEAT *delta_ptr = feat+num_features;
	for (i=0; i<num_vectors; i++) {
		for (j=0; j<num_base_features; j++) {  // Base features are always 0 .. (num_base_feat-1)
			if (i-delta_spread < 0)
				delta_ptr[j] = -feat_new[j]; // fill to left with first frame
			else 
				delta_ptr[j] = -feat_ptr[-delta_spread*new_num_feat+j];
			if (i+delta_spread >= num_vectors)
				delta_ptr[j] += feat_new[(num_vectors-1)*new_num_feat+j]; // fill to right with last frame
			else 
				delta_ptr[j] += feat_ptr[delta_spread*new_num_feat+j];

		}
		feat_ptr += new_num_feat;
		delta_ptr += new_num_feat;
	}

	// Update state
	delta_available = true;
	delta_offset = num_features;
	num_features = new_num_feat;

}

void Features::feat_norm ()
{
	static float vfloor = (float) 1.0e-6;

	if (!features_available || num_vectors==0)
		throw ST_exception("Features::feat_norm -- features not available.");

	double *mean = new double[num_features];
	double *std = new double[num_features];
	int i, j;

	for (i=0; i<num_features; i++) {
		mean[i] = 0;
		std[i] = 0;
	}

	for (i=0; i<num_vectors; i++) {
		for (j=0; j<num_features; j++) {
			mean[j] += feat[i*num_features+j];
			std[j] += feat[i*num_features+j]*feat[i*num_features+j];
		}
	}

	for (i=0; i<num_features; i++) {
		mean[i] /= num_vectors;
		std[i] /= num_vectors;
		std[i] = sqrt(std[i] - mean[i]*mean[i]);
		if (std[i] < vfloor)
			std[i] = vfloor;
		std[i] = 1/std[i];
	}

	for (i=0; i<num_vectors; i++) {
		for (j=0; j<num_features; j++) 
			feat[i*num_features+j] = std[j]*(feat[i*num_features+j]-mean[j]);
	}

	delete[] mean;
	delete[] std;

}

bool Features::no_speech (void)
{
   int i;
   bool speech = false;

	if (!sad_available)
		throw ST_exception("Features::no_speech -- no SAD available.");

   for (i=0; i<num_vectors; i++)  {
      if (sad[i]==1) {
         speech = true;
         break;
      }
   }

   return !speech;

}

vec<REAL_FEAT> Features::vector (int index) 
{
	if ((index < 0) || (index>num_vectors))
		throw ST_exception("Features::feat_vec -- index out of range.");

	if (outfeat.compare("all")==0) {
		vec<REAL_FEAT> out(&feat[index*num_features],num_features);
		return out;
	} else {
		vec<REAL_FEAT> out(0);
		for (int i=0; i<outfeat.length(); i++)  {
			switch (outfeat[i]) {
		   case 'f' : 
				if (!features_available)
					throw ST_exception("Features::feat_vec -- no features available.");
				out.append(&feat[index*num_features], num_base_features);
				break;
			case 'd' :
				if (!delta_available)
					throw ST_exception("Features::feat_vec -- no delta features available.");
				out.append(&feat[index*num_features+delta_offset], num_base_features);
				break;
			case 'a' :
				if (!accel_available)
					throw ST_exception("Features::feat_vec -- no accel features available.");
				out.append(&feat[index*num_features+accel_offset], num_base_features);
				break;
			case 'e' : 
				if (!energy_available)
					throw ST_exception("Features::feat_vec -- energy not available.");
				out.append(&energy[index], 1);
				break;
			case 's' :
				if (!sdc_available)
					throw ST_exception("Features::feat_vec -- no sdc features available.");
				out.append(&feat[index*num_features+sdc_offset], sdc_k*num_base_features);
				break;
			default:
				throw ST_exception("Features::feat_vec -- unknown feature type in outfeat.");
			}
		}
		return out;
	}

}

void Features::load_raw (string filename, int num_base_feat) 
{
	FILE  *infile;

	if (num_base_feat <= 0) 
		throw ST_exception("Features::load_raw -- Number of base features is <= 0.");

   // Find number of feature vectors
	int num_floats = file_len(filename.c_str());
	num_floats /= (int) sizeof(REAL_FEAT);

   if (num_floats==0)
      throw ST_exception("Features::load_raw -- Empty features file.");
   if ((num_floats % num_base_feat) != 0)
      throw ST_exception("Features::load_raw -- Corrupted feature file?");

	// Allow overwriting
	if (features_available) {
		if (feat!=0)
			delete[] feat;
		feat = 0;
		if (energy!=0)
			delete[] energy;
		energy = 0;

		num_features = 0;
		num_base_features = 0;
		num_vectors = 0;

		delta_offset = 0;
		accel_offset = 0;
	
		features_available = false;
		energy_available = false;
		delta_available = false;
		accel_available = false;
	
		if (sad!=0)
			delete[] sad;
		sad_available = false;
		sad = 0;
		sad_len = 0;
		outfeat = "all";
	}

   // Allocate space
   Features::num_base_features = num_base_feat;
	Features::num_features = num_base_feat;
   Features::num_vectors = num_floats/num_base_feat;
   Features::feat = new REAL_FEAT[num_features*num_vectors];

   // Open and read file.
   infile = fopen(filename.c_str(), "rb");
   if (infile==NULL)
      throw ST_exception("Features::load_raw -- Unable to open features file.");
   int num_read = fread(feat, sizeof(REAL_FEAT), num_features*num_vectors, infile);
   fclose(infile);

	features_available = true;

}  // Features::load_raw

void Features::load_SAD_label (string filename, bool sloppy) 
{ 
	// 0/1 label file
	
	FILE *infile = fopen(filename.c_str(), "r");
	if (infile==NULL)
		throw ST_exception("Features::load_SAD_label -- Unable to open SAD file.");
	if (!features_available) 
		throw ST_exception("Features::load_SAD_label -- Must compute features before loading SAD labels.");
	if (sad_available) { // Overwrite old SAD
		delete[] sad;
		sad_available = false;
		sad_len = 0;
	}
	sad = new short[num_vectors];

	int val, num_read = 0;

	while (fscanf(infile, " %d ", &val)==1) {
		num_read++;
		if (num_read > num_vectors) {
			if (sloppy)
				break;
			else
				throw ST_exception("Features::load_SAD_label -- Label file length does not match number of features.");
		} 
		if (! (val == 0 || val==1)) 
			throw ST_exception("Features::load_SAD_label -- Unexpected value in label file (should be 0/1).");
		sad[num_read-1] = val;
	}
	fclose(infile);

	if (num_read < num_vectors) {
		if (sloppy) {
			for (int i=num_read; i<num_vectors; i++)
				sad[i] = 0;
		} else
			throw ST_exception("Features::load_SAD_label -- Label file length does not match number of features.");
	}
	sad_available = true;
	sad_len = num_vectors;


} // load_SAD_label

void Features::and_user_marks (const char *user_marks, int blocksz_ms, float in_win_inc_ms)
{
	if (!sad_available) 
		throw ST_exception("Features::and_user_marks -- no SAD to and with.");

   if (user_marks==NULL || strlen(user_marks)==0)
      return;

   int i, j;

   string user_marks_str(user_marks);
   if (user_marks_str.substr(0,9).compare("markfile:")==0) {
      string mark_file = user_marks_str.substr(9,string::npos);
      short *sad_save = new short[num_vectors];
      for (i=0; i<num_vectors; i++)
         sad_save[i] = sad[i];
      try {
         load_SAD_marks(mark_file, in_win_inc_ms);
      } catch (ST_exception &e) {
         delete[] sad_save;
         throw ST_exception(string("Features::and_user_marks -- ") + e.get_error_message());
      } catch (...) {
         delete[] sad_save;
         throw ST_exception("Features::and_user_marks -- unexpected error.");
      }
      for (i=0; i<num_vectors; i++) {
         if (sad_save[i]==1 && sad[i]==1) 
            sad[i] = 1;
         else
            sad[i] = 0;
      }
      delete[] sad_save;
      return;
   }  

   float user_st, user_en;
   int sad_st, sad_en=0;
   string::iterator it;

   for (i=0, j=0, it=user_marks_str.begin(); it!=user_marks_str.end(); it++, i++) {
      user_st = i*((float) blocksz_ms)/1000.0;
      user_en = (i+1)*((float) blocksz_ms)/1000.0;
      sad_st = floor(user_st/(0.001*in_win_inc_ms));
      sad_en = floor(user_en/(0.001*in_win_inc_ms));
      if (sad_en>=num_vectors)
         sad_en = num_vectors-1;
      if ((*it)!='S' && (*it)!='s') {
         for (j=sad_st; j<=sad_en; j++)
            sad[j] = 0;
      }
   }

   // If we haven't reached the end, zero out the rest
   for (j=sad_en; j<num_vectors; j++) 
      sad[j] = 0;

}

void Features::load_SAD_marks (list<Time_Interval> &seg_sad, float win_inc)
{
   if (!features_available) 
		throw ST_exception("Features::load_SAD_marks -- must compute features before loading SAD labels.");

   // Overwrite old SAD
	if (sad_available) { 
		delete[] sad;
		sad_available = false;
		sad_len = 0;
	}
	sad = new short[num_vectors];
	int i;
	for (i=0; i<num_vectors; i++)
		sad[i] = 0;

	double start, dur;
	int start_frame, end_frame;

   list<Time_Interval>::iterator it;

   for (it=seg_sad.begin(); it!=seg_sad.end(); it++) {
		start_frame = (int) (it->start/(win_inc*1e-3)+0.5);
      end_frame = ((int) (it->end/(win_inc*1e-3)+0.5))-1;
		if (start_frame < 0) 
			throw ST_exception(string("Features::load_SAD_marks -- unexpected negative time in marks."));
		if (start_frame >= num_vectors)
			continue;
		if (end_frame >= num_vectors)
			end_frame = num_vectors-1;
		for (i=start_frame; i<=end_frame; i++)
			sad[i] = 1;
   }
	sad_available = true;
	sad_len = num_vectors;
}

void Features::load_SAD_marks (string filename, float win_inc) 
{ 
	// tag start dur
	FILE *infile = fopen(filename.c_str(), "r");
	if (infile==NULL)
		throw ST_exception(string("Features::load_SAD_marks -- unable to open SAD file ")+filename);
	if (!features_available) 
		throw ST_exception("Features::load_SAD_marks -- must compute features before loading SAD labels.");
	if (sad_available) { // Overwrite old SAD
		delete[] sad;
		sad_available = false;
		sad_len = 0;
	}
	sad = new short[num_vectors];
	int i;
	for (i=0; i<num_vectors; i++)
		sad[i] = 0;

	double start, dur;
	int start_frame, end_frame;

	while (fscanf(infile, " %*s %lf %lf ", &start, &dur)==2) {
		start_frame = (int) (start/(win_inc*1e-3)+0.5);
      end_frame = start_frame + ((int) (dur/(win_inc*1e-3)+0.5)) - 1;
		if (start_frame < 0) 
			throw ST_exception(string("Features::load_SAD_marks -- mark file does not match number of features ")+filename);
		if (start_frame >= num_vectors)
			continue;
		if (end_frame >= num_vectors)
			end_frame = num_vectors-1;
		for (i=start_frame; i<=end_frame; i++)
			sad[i] = 1;
	}
	fclose(infile);
	sad_available = true;
	sad_len = num_vectors;

} // load_SAD_marks

string Features::load_SAD_marks_string (string filename, float win_inc) 
{ 
	// tag start dur
	FILE *infile = fopen(filename.c_str(), "r");
	if (infile==NULL)
		throw ST_exception(string("Features::load_SAD_marks_string -- unable to open SAD file ")+filename);
	if (!features_available) 
		throw ST_exception("Features::load_SAD_marks -- must compute features before loading SAD labels.");

   string sad(num_vectors, 'n');  // fill with non-speech

	double start, dur;
	int i, start_frame, end_frame;

	while (fscanf(infile, " %*s %lf %lf ", &start, &dur)==2) {
		start_frame = (int) (start/(win_inc*1e-3)+0.5);
      end_frame = start_frame + ((int) (dur/(win_inc*1e-3)+0.5)) - 1;
		if (start_frame < 0) 
			throw ST_exception(string("Features::load_SAD_marks -- mark file does not match number of features ")+filename);
		if (start_frame >= num_vectors)
			continue;
		if (end_frame >= num_vectors)
			end_frame = num_vectors-1;
		for (i=start_frame; i<=end_frame; i++)
			sad[i] = 's';
	}
	fclose(infile);

   return sad;

} // load_SAD_marks_string

int Features::num_base_feat (void)
{
	return num_base_features;
}

int Features::num_outfeat (void)
{
	if (outfeat.compare("all")==0) 
		return num_features;

	char c;
	int total = 0;
	for (int i=0; i<outfeat.length(); i++) {
		c = outfeat[i]; 
		switch (c) {
		   case 'f':
		   case 'a':
		   case 'd': 
				total += num_base_features;
			   break;
		   case 's':
				total += sdc_k*num_base_features;
				break;
		   case 'e' : 
				total += 1;
				break;
		}
	}
	return total;
}

int Features::num_total_feat (void)
{
	int total_feat = num_features;

	if (sdc_available)
		total_feat += sdc_k*num_base_features;
	return total_feat;

}

int Features::num_vec (void)
{
	return num_vectors;
}

void Features::process (Signal &x, string &config_str) 
{
	// Doesn't do anything for the base class
}

void Features::rasta ()
{
	//
	// Apply the filter : (0.2 + 0.1 z^(-1) - 0.1 z^(-2) - 0.2 z^(-3))/(1-0.98 z^(-1)) 
	//
	// Note: LL filter is slightly different, pole is at 0.98 instead of 0.94

	static float coeff[5] = {0.2, 0.1, 0, -0.1, -0.2};

	if (!features_available || num_vectors==0)
		throw ST_exception("Features::rasta -- features not available.");

	// Initialize, allocate space, copy old features
	int i, j, k;
	REAL_FEAT *feat_new = new REAL_FEAT[num_features*num_vectors];
	REAL_FEAT *y = new REAL_FEAT[num_features];
	REAL_FEAT *outvec = new REAL_FEAT[num_features];
	for (i=0; i<num_features; i++)  
		y[i] = 0;

	// Filter
	REAL_FEAT *in_ptr = feat;
	REAL_FEAT *out_ptr = feat_new;
	for (i=0; i<=(num_vectors+1); i++) {
		for (j=0; j<num_features; j++) {
			outvec[j] = 0.98*y[j];
			for (k=0; k<5; k++) {
				if (i-k < 0) 
					outvec[j] += coeff[k]*feat[j]; 	// Use first vector for boundary conditions
				else if (i-k >= num_vectors)
					outvec[j] += coeff[k]*feat[(num_vectors-1)*num_features+j];
				else 
					outvec[j] += coeff[k]*in_ptr[-k*num_features+j];
			}
			y[j] = outvec[j];
		}
		in_ptr += num_features;
		if (i>1) { // Group delay of numerator is 2 samples
			for (j=0; j<num_features; j++)
				out_ptr[j] = outvec[j];
			out_ptr += num_features;
		}
	}

	// Clean up
	delete[] feat;
	feat = feat_new;
	delete[] y;
	delete[] outvec;

}

int Features::SAD (int i)
{
	if (!sad_available) 
		throw ST_exception("Features::SAD -- no SAD has been computed.");
   if (i>=num_vectors || i<0)
		throw ST_exception("Features::SAD -- index out of range.");
   return sad[i];
}

vec<uint8_t> Features::SAD_labels(void)
{
	vec<uint8_t> out(num_vectors);
	if (!sad_available) 
		throw ST_exception("Features::SAD -- no SAD has been computed.");
	int i;
	for (i=0; i<num_vectors; i++)
		out.data[i] = sad[i];
	return out;
}

bool Features::SAD_available (void)
{
   return sad_available;
}

void Features::save_raw (string filename)
{
	FILE 	*outfile;

	if (!features_available) 
		throw ST_exception("Features::save_raw -- no features to save.");
	outfile = fopen(filename.c_str(), "wb");
	if (outfile == NULL)
		throw ST_exception("Unable to open output file.\n");
	for (int i=0; i<num_vectors; i++) {
		vec<REAL_FEAT> fvec = vector(i);
		fwrite (fvec.data, sizeof(REAL_FEAT), fvec.len, outfile);
	}
	fclose(outfile);

} // Features::save_raw

void Features::save_SAD_label (string filename)
{
	FILE 	*outfile;

	if (!sad_available) 
		throw ST_exception("Features::save_sad_labels -- no SAD to save.");
	outfile = fopen(filename.c_str(), "wb");
	if (outfile == NULL)
		throw ST_exception("Unable to open output file.\n");

	for (int i=0; i<num_vectors; i++)
		fprintf(outfile, "%d\n", (int) sad[i]);

	fclose(outfile);

} // Features::save_sad_labels

void Features::save_SAD_marks (string filename, float in_win_inc_ms)
{
	FILE 	*outfile;

	if (!sad_available) 
		throw ST_exception("Features::save_sad_marks -- no SAD to save.");
	outfile = fopen(filename.c_str(), "wb");
	if (outfile == NULL)
		throw ST_exception("Unable to open output file.\n");

   float seg_start, seg_end;
   bool in_seg = false;

	for (int i=0; i<num_vectors; i++) {
      if (!in_seg && sad[i]==1) {
         seg_start = 0.001*i*in_win_inc_ms;
         in_seg = true;
      } else if (in_seg && sad[i]==0) {
         seg_end = 0.001*i*in_win_inc_ms;
         fprintf(outfile, "speech %.2f %.2f\n", seg_start, seg_end-seg_start);
         in_seg = false;
      }
   }

   if (in_seg) {
      seg_end = 0.001*(num_vectors-1)*in_win_inc_ms;
      fprintf(outfile, "speech %.2f %.2f\n", seg_start, seg_end-seg_start);
   }

   fclose(outfile);

} // Features::save_sad_labels

void Features::sdc (int sdc_p, int sdc_k_in) { 

	if (sdc_p < 1)
		throw ST_exception("Features::sdc -- Shift parameter, P,  must be >= 1.");
	if (sdc_k_in < 1)
		throw ST_exception("Features::sdc -- Stacking parameter, k, must be >= 1.");
	if (!delta_available)
		throw ST_exception("Features::sdc -- Delta coefficients must be computed before sdc can be computed.");

	// Initialize, allocate space, copy old features
	sdc_k = sdc_k_in;
	int i, j, k;
	int new_num_feat = num_features+sdc_k*num_base_features;
	REAL_FEAT *feat_new = new REAL_FEAT[new_num_feat*num_vectors];
	for (i=0; i<num_vectors; i++) {
		for (j=0; j<num_features; j++) 
			feat_new[i*new_num_feat+j] = feat[i*num_features+j];
	}
	delete[] feat;
	feat = feat_new;

	// Could also have used k=ceil(-sdc_k/2), < ceil(sdc_k/2) 
	// for both cases if there was a standard fn in C++
	int k1, k2, l;
	if ((sdc_k % 2) == 1) {
		k1 = -(sdc_k-1)/2;
		k2 = (sdc_k+1)/2;
	} else {
		k1 = -sdc_k/2;
		k2 = sdc_k/2;
	}

	// Stack sdc_k vectors
	// Warning: These are different than LLSpeech when sdc_k is even
	// LLSpeech has a variable frame delay between features and sdc features in this case
	REAL_FEAT *feat_ptr = feat;
	REAL_FEAT *sdc_ptr = feat+num_features;
	for (i=0; i<num_vectors; i++) {
		for (k=k1; k<k2; k++) {  
			j = i+k*sdc_p;
			if (j<0)
				j = 0;
			if (j>=num_vectors)
				j = num_vectors-1;
			for (l=0; l<num_base_features; l++)
				sdc_ptr[l] = feat[j*new_num_feat+delta_offset+l];
			sdc_ptr += num_base_features;
		}

		feat_ptr += new_num_feat;
		sdc_ptr = feat_ptr+num_features;
	}

	// Update features state
	sdc_available = true;
	sdc_offset = num_features;
	num_features = new_num_feat;

}

void Features::set_outfeat (string str)
{
	for (int i=0; i<str.length(); i++) {
		if (str.substr(i,1).find("feads")!=string::npos)
			throw ST_exception("Features:: set_outfeat -- unknown output feature type.");
	}
	outfeat = str;
}


void Features::xtalk (double abs_min_energy, double thresh, int med_len)
{
	int i, j, k, l;

	// xtalk refactored from LLSpeech
	// xtalk is an online algorithm, so it uses circular buffers and running estimates
	// which makes it more difficult to follow than a batch algorithm

	// Sanity checks
	if (!energy_available) 
		throw ST_exception("Features::xtalk -- no energy calculated.");
	if (med_len < 1)
		throw ST_exception("Features::xtalk -- median filter length must be positive.");

	// Works ok with "pathological" cases -- num_vectors = 1, 2
	// xtalk in LLSpeech "fails" for 1 vector inputs -- produces 2 outputs
	sad = new short[num_vectors];

	// First part refactored from xtalk_create and xtalk_reset
	// For xt_logic_off thresh=max=1
	// for xt_logic_on defaults are 15 and 30
	int counter = 0;
	int counter_thresh = 1;
	int max_counter = 1;

	// This appears to be buggy in the original code -- num_back is never set
   const int look_back = 10;
	const int num_back = 1;
	int gapc = 0;
	int adapt_interval = 500;

	float min_step = (float) 0.002166062;
	float max_step = (float) 0.002166062;
	float decay_coeff = (float) 0.9;
	float max_jump = (float) 100.0;
	float min_e, max_e, decay_e;

	// Calculated constants
	if (!(med_len % 2)) med_len++; // med_len always odd
	int delay_idx = -num_back;
	int med_idx = delay_idx - (med_len-1)/2;
	int med_half = (med_len-1)/2;;
	
	// Local memory -- small buffers, use the stack, buffers are circular
	float delay[look_back];   // array of frame energies
	int ch_active[num_back];  // unsmoothed SAD
	int *active = new int[med_len];      // median filtered SAD

	// Refactored from xtalk_execute
	bool first_time = true;
	int pnum_frames = 0;
	int ibuf, tot;
	float norm_e;

	for (int num_frames=0; med_idx<(num_vectors-1); num_frames++) {
		ibuf = num_frames % num_back;
		if (num_frames < num_vectors) {
			if (energy[num_frames] < abs_min_energy) {
				ch_active[ibuf] = 0;
			} else {
				if (first_time) {  // Don't initialize till we get a frame above abs min
					min_e = energy[num_frames];
					max_e = energy[num_frames];
					decay_e = energy[num_frames];
					counter = 0;
					for (i=0; i<look_back; i++)
						delay[i] = max_e;
					ch_active[ibuf] = 0;
				   first_time = false;
				}
				decay_e = decay_coeff*decay_e + (1.0-decay_coeff)*energy[num_frames];
				if (decay_e < min_e)
					min_e = decay_e;
				if (energy[num_frames] > max_e)
					max_e = energy[num_frames];

				if (num_frames >= adapt_interval) {
					if (num_frames >= look_back)
						i = (num_frames-look_back) % look_back;
					else
						i = (look_back-num_frames) % look_back;
					if (fabs(min_e-delay[i]) > max_jump) {
						// Sequential counter if (filling gap > 2*look_back) then probably at a new noise level
						if ((num_frames-pnum_frames) == 1) {
							gapc++;
						} else {
							gapc = 0;
						}
						pnum_frames = num_frames;
						if (gapc < (2*look_back))
							min_e = delay[i];
					}
				}		
				delay[num_frames % look_back] = min_e;

				// Normalize energy to noise floor
				norm_e = energy[num_frames] - min_e;
  
				// Slowly raise min and drop max
				min_e += min_step;
				max_e -= max_step;
  
				// Update counters
				if (norm_e > thresh) {
					counter++;
					if (counter > max_counter) {
						counter = max_counter;
					} 
				} else {
					counter--;
					if (counter<0)
						counter = 0;
				}
				if (counter == counter_thresh) {
					// non-speech <--> speech, flag current frame as active and backtrack num_back frames
					for (i=0; i<num_back; i++)
						ch_active[i] = 1;
				}
				else if (counter > counter_thresh) {
					// speech <--> speech, flag current frame as active
					ch_active[ibuf] = 1;
				} else {
					// non-speech <--> non-speech, flag current frame as inactive
					ch_active[ibuf] = 0;
				}
			}
		} // num_frames < num_vectors

		// Median filter -- filters ch_active -> active
		med_idx++;
		delay_idx++;
		if (delay_idx >= 0) { 
			i = delay_idx % num_back;
			k = delay_idx % med_len;

			// Store first pass channel activity labels
			active[k] = ch_active[i];

			// Median smoothing of channel activities 
			if (med_idx >= 0) {
				l = med_idx % med_len;
			   if (med_idx >= med_half) {
					for (j=0, tot=0; j<med_len; j++)
						tot += active[j];
					if (tot > med_half) 
						active[l] = 1;
					else 
						active[l] = 0;
				}
				sad[med_idx] = active[l];
			}
		}

	} // for 

   delete[] active;
	sad_len = num_vectors;
	sad_available = true;

}

//
// Local functions
//
static REAL_FEAT dcep_constant (int delta_spread)
{
int i;
REAL_FEAT dcep_c;

	for (dcep_c=0, i=-delta_spread; i<=delta_spread; i++)
		dcep_c += i*i;
	dcep_c = 1/dcep_c;
	return(dcep_c);

} // dcep_constant

