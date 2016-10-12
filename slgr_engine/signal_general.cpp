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
// Speech Tools (ST) library
// 
// General functions for signal processing
//

// 10/02, wmc, v1
// 08/10, wmc, v2 SCA
// 02/13, wmc, added load_pcm_wav from bhd_wav
// 11/13, wmc, added zero phase FIR, resampler

// TODO: Look into indexing into arrays using size_t (be careful):
//       better portability and handling longer duration signals.

//
// Global includes.
//
#include "stdafx.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include "speech_tools.h"
#include <string.h>
#include <cstdint>

using namespace std;

static int alaw2linear(uint8_t alawbyte);
static int ulaw2linear(uint8_t ulawbyte);

unsigned file_len (const char *infile)
{
	FILE *fp = fopen(infile, "rb");
	if (fp==NULL) 
		throw ST_exception(string("Unable to open file ")+infile);
	if (fseek(fp, 0, SEEK_END)) {
		fclose(fp);
		throw ST_exception(string("Unable to open file ")+infile);
	}
	off_t r = ftell(fp);
	if (r == ((off_t)(-1))) {
		fclose(fp);
		throw ST_exception(string("Unable to open file ")+infile);
	}
	fclose(fp);
	return ((unsigned) r);
}

//
// Constructors/destructors
// 
Signal::Signal (const Signal &x)
{
	num_samples = x.num_samples;
	sampling_frequency = x.sampling_frequency;
   if (x.value != 0) {
	   value = new float[x.num_samples];
      int i;
	   for (i=0; i<x.num_samples; i++)
		   value[i] = x.value[i];
   }
   resampler_loaded = false; 

} // Copy constructor

Signal::Signal (int sampling_frequency)
{
   Signal::sampling_frequency = (float) sampling_frequency;
   Signal::num_samples = 0;
   value = 0;
   resampler_loaded = false;  
}

Signal::Signal ()
{
	num_samples = 0;
	sampling_frequency = 0;
	value = 0;
   resampler_loaded = false;  
}

Signal::Signal (const short *data, int num_samples, int sampling_frequency)
{
int i;

Signal::sampling_frequency = (float) sampling_frequency;
	Signal::num_samples = num_samples;
	value = new float[num_samples];
	for (i=0; i<num_samples; i++)
		value[i] = data[i];
   resampler_loaded = false;  

} // Signal::Signal

Signal::Signal (const double *data, int num_samples, int sampling_frequency)
{
int i;

Signal::sampling_frequency = (float) sampling_frequency;
	Signal::num_samples = num_samples;
	value = new float[num_samples];
	for (i=0; i<num_samples; i++)
		value[i] = (float) data[i];
   resampler_loaded = false;  

} // Signal::Signal

Signal::~Signal (void)
{
	delete[] value;
} // Signal::~Signal

//
// Signal class member functions.
//
void Signal::comp_ifir_filter (const float *h1, int n, int filter_length_1, const float *h2, int filter_length_2)
//
// Implement a filter of the form 1-H_1(z^n)H_2(z).
//
{
int	i, grp_delay;
float	*new_value;

	grp_delay = ((filter_length_1-1)*n+filter_length_2-1);
	if ((grp_delay % 2) != 0)
		throw ST_exception("IFIR filter is not odd length.");
	grp_delay /= 2;

	Signal x2 = *this;
	x2.fir_n_filter(h1,filter_length_1,n);
	x2.fir_filter(h2,filter_length_2);

	new_value = new float[x2.num_samples];

	for (i=0; i<grp_delay; i++)
	  new_value[i] = -x2.value[i];

	for (i=0; i<num_samples; i++)
		new_value[i+grp_delay] = value[i]-x2.value[i+grp_delay];

	for (i=num_samples+grp_delay; i<x2.num_samples; i++)
		new_value[i] = -x2.value[i];

	delete[] value;
	value = new_value;
	num_samples = x2.num_samples;

} // comp_ifir_filter

void Signal::fir_filter (const float *h, int filter_length)
//
// FIR filters an input and creates an output.  
//
{
int  i, j, k, new_num_samples;
float	*new_value, sum;

	if (num_samples == 0) 
		throw ST_exception("No signal to filter.");

	new_num_samples = num_samples+filter_length-1;
	new_value = new float[new_num_samples];

	for (i=0; i<num_samples; i++) {
		sum = 0;
		for (k=0, j=i; k<filter_length && j>=0; k++, j--)
			sum += h[k]*value[j];
		new_value[i] = sum;
	}

	for (i=num_samples; i<new_num_samples; i++) {
		sum = 0;
		for (k=0, j=i; k<filter_length && j>=0; k++, j--)
			if (j<num_samples)
				sum += h[k]*value[j];
		new_value[i] = sum;
	}

	num_samples = new_num_samples;
	delete[] value;
	value = new_value;

} // fir_filter

void Signal::fir_n_filter (const float *h, int filter_length, int n)
//
// FIR filters an input by H(z^n) and creates an output.  
//
{
int	filter_length_n, i, j, k, new_num_samples;
float	*new_value, sum;

	if (num_samples==0)
		throw ST_exception("No signal to filter.");

   filter_length_n = (filter_length-1)*n+1;
	new_num_samples = num_samples+filter_length_n-1;
	new_value = new float[new_num_samples];

	for (i=0; i<num_samples; i++) {
		sum = 0;
		for (k=0, j=i; k<filter_length && j>=0; k++, j-=n)
			sum += h[k]*value[j];
		new_value[i] = sum;
	}

	for (i=num_samples; i<new_num_samples; i++) {
		sum = 0;
		for (k=0, j=i; k<filter_length && j>=0; k++, j-=n)
			if (j<num_samples)
				sum += h[k]*value[j];
		new_value[i] = sum;
	}

	num_samples = new_num_samples;
	delete[] value;
	value = new_value;

} // fir_n_filter

void Signal::fir_updn (const double *h, const int filter_len, int L, int M)
{
   if ((filter_len%2)!=1)
      throw ST_exception("Filter for zero phase must be odd length.");

   int i1;
   bool symmetric = true;
   int middle = (filter_len-1)/2;
   for (i1=1; (middle-i1)>=0; i1++) {
      if (h[middle-i1]!=h[middle+i1]) {
         symmetric = false;
         break;
      }
   }
   if (!symmetric)
      throw ST_exception("Filter for zero phase must be symmetric.");

   // Set up for simpler indexing
   h += middle;

   // Allocate space for output
   int y_end = (L*((int64_t)(num_samples-1)))/M;  // this can get big -- use a 64 bit unsigned int
   float *y = new float[y_end+1];

   // Main loop
   // Somewhat tricky in terms of indexing:
   // convert intermediate calculations to 64 bit to avoid overflow for long files
   int64_t top, i, y_index, h_start, h_end;
   double y_new;
   for (y_index=0; y_index<=y_end; y_index++) {
      y_new = 0;
      top = M*y_index-middle;  // need 64 bit since y_index*M can overflow 32 bits
      h_start = (int64_t) ceil(((double) top)/L);
      h_end = (M*y_index+middle)/L;
      if (h_start < 0) // value[i]=0 for i<0
         h_start = 0;
      if (h_end>=num_samples) // value[i]=0 for i>=num_samples
         h_end = num_samples-1;
      for (i=h_start; i<=h_end; i++)
         y_new += value[i]*h[M*y_index-L*i];
      y[y_index] = (float) y_new;
   }

   // Result
   delete[] value;
   value = y;
   num_samples = y_end+1;

}

void Signal::fir_zphase (const double *h, const int filter_len)
{
	if (num_samples==0)
		throw ST_exception("No signal to filter.");

   if ((filter_len%2)!=1)
      throw ST_exception("Filter for zero phase must be odd length.");

   int i;
   bool symmetric = true;
   int middle = (filter_len-1)/2;
   for (i=1; (middle-i)>=0; i++) {
      if (h[middle-i]!=h[middle+i]) {
         symmetric = false;
         break;
      }
   }
   if (!symmetric)
      throw ST_exception("Filter for zero phase must be symmetric.");

   float *new_value = new float[num_samples];  // no edge effects, zero phase
   int i1, i2, j, k;  // might use size_t here, but careful, checking for negative numbers below

   for (i=0; i<num_samples; i++) {
      new_value[i] = 0.0;
      i1 = i-middle;
      if (i1<0)
         i1 = 0;
      i2 = i+middle;
      if (i2>=num_samples)
         i2 = num_samples-1;
      for (j=middle, k=i; k>=i1 && j<filter_len; k--, j++)
         new_value[i] += (float) (h[j]*value[k]);
      for (j=middle+1, k=i+1; k<=i2 && j<filter_len; k++, j++)
         new_value[i] += (float) (h[j]*value[k]);
   }

   delete[] value;
   value = new_value;

}

void Signal::init_resampler (string filter_dir)
{
   if (resampler_loaded)
      return;

   int len;
   const char *filter_names[8] = {"h_dn2","h_dn3","h1_dn4","h2_dn4","h1_dn6","h2_dn6","interp_11k_8k","low_pass_11k_8k"};

   for (int i=0; i<8; i++) {
      string fn = filter_dir + "/" + filter_names[i] + ".dat";
      len = file_len(fn.c_str());
      len = len/sizeof(double);
      if (len==0)
         throw ST_exception(string("error initializing resampler -- empty filter file ") + fn);
      string filter_name = filter_names[i];
      resampler_filters.insert(make_pair(filter_name, vec<double>(len)));
      FILE *infile = fopen(fn.c_str(), "rb");
      if (infile==NULL)
         throw ST_exception(string("error initializing resampler -- unable to open filter file ") + fn);

      // Interesting "feature" : need empty constructor vec<double>() in order for this to work
      // operator[] for map<,> needs to be able to insert new keys if they aren't available
      if (fread(resampler_filters[filter_name].data, sizeof(double), len, infile)!=len) {
         fclose(infile);
         throw ST_exception(string("error initializing resampler -- unable to read in filter ") + fn);
      }
      fclose(infile);
   }

   resampler_loaded = true;

}

void Signal::load_raw_short (string filename, int sampling_frequency)
{
	FILE  *infile;
	short value_short;

	if (sampling_frequency <= 0)
      throw ST_exception("Signal::load_raw_short -- Sampling frequency must be >0.");

   infile = fopen(filename.c_str(), "rb");
   if (infile == NULL)
      throw ST_exception("Signal::load_raw_short -- Unable to open input file.");

   if (num_samples != 0) // Delete old data
		delete[] value;

	num_samples = file_len(filename.c_str());
	num_samples /= 2;
   if (num_samples <= 0)
      ST_exception("Signal::load_raw_short -- Empty input file.");

   value = new float[num_samples];
	size_t num_read;
   for (int i=0; i<num_samples; i++) {
      num_read = fread(&value_short,sizeof(short),1,infile);
      value[i] = (float) value_short;
   }
	fclose(infile);

	Signal::sampling_frequency = (float) sampling_frequency;

} // load_raw_short

void Signal::load_sph (string filename, int channel_num)
{
	FILE  *infile;

	unsigned st_size = file_len(filename.c_str());
	// stat(filename.c_str(), &buf);
   infile = fopen(filename.c_str(), "rb");
   if (infile == NULL)
      throw ST_exception("Signal::load_pcm_sph -- Unable to open input file.");

	// Read header -- old style IO, ugly code -- probably need to update
	char nist_mark[11];
	char *d = fgets(nist_mark, 8, infile);
	if (d==NULL) {
		fclose(infile);
		throw ST_exception("Signal::load_pcm_sph -- unable to read header for sphere file.");
	}
	int num_read = (int) strlen(nist_mark);
	if (num_read!=7 || strcmp(nist_mark,"NIST_1A")!=0) {
		fclose(infile);
		throw ST_exception("Signal::load_pcm_sph -- unable to read header for sphere file.");
	}
	int num_bytes;
	if (fscanf(infile, " %d ", &num_bytes)!=1) {
		fclose(infile);
		throw ST_exception("Signal::load_pcm_sph -- unable to read header for sphere file.");
	}
	rewind(infile);
	if (num_bytes<=0) {
		fclose(infile);
		throw ST_exception("Signal::load_pcm_sph -- unable to read header for sphere file.");
	}
	char *hdr = new char[num_bytes+1];
	if (fread(hdr, sizeof(char), num_bytes, infile)!=num_bytes) {
		fclose(infile);
		delete[] hdr;
		throw ST_exception("Signal::load_pcm_sph -- unable to read header for sphere file.");
	}
	hdr[num_bytes] = '\0';

	// Parse header 
	char *loc;
	loc = strstr(hdr, "sample_rate -i");
	if (loc==NULL) {
		fclose(infile);
		delete[] hdr;
		throw ST_exception("Signal::load_sph -- couldn't find sample rate in header.");
	}
	if (sscanf(loc, "sample_rate -i %f ", &sampling_frequency)!=1) {
		fclose(infile);
		delete[] hdr;
		throw ST_exception("Signal::load_sph -- couldn't find sample rate in header.");
	}
	char format[5]={0};
	loc = strstr(hdr, "sample_coding -s");
	if (loc==NULL) {
		fclose(infile);
		delete[] hdr;
		throw ST_exception("Signal::load_sph -- couldn't find sample coding in header.");
	}
	int num;
	if (sscanf(loc, "sample_coding -s%d %4s", &num, format)!=2) {
		fclose(infile);
		delete[] hdr;
		throw ST_exception("Signal::load_sph -- sample coding malformed in header.");
	}
	if (strcmp(format,"pcm")!=0 && strcmp(format,"ulaw")!=0 && strcmp(format,"alaw")!=0) {
		fclose(infile);
		delete[] hdr;
		throw ST_exception("Signal::load_sph -- sample coding not supported -- must be ulaw, alaw, or pcm.");
	}
	loc = strstr(hdr, "sample_n_bytes -i");
	if (loc==NULL) {
		fclose(infile);
		delete[] hdr;
		throw ST_exception("Signal::load_sph -- couldn't find sample number of bytes in header.");
	}
	int sample_n_bytes;
	if (sscanf(loc, "sample_n_bytes -i%d", &sample_n_bytes)!=1) {
		fclose(infile);
		delete[] hdr;
		throw ST_exception("Signal::load_sph -- sample number of bytes malformed in header.");
	}
	loc = strstr(hdr, "channel_count -i");
	if (loc==NULL) {
		fclose(infile);
		delete[] hdr;
		throw ST_exception("Signal::load_sph -- couldn't find channel count in header.");
	}
	int num_channels;
	if (sscanf(loc, "channel_count -i%d", &num_channels)!=1) {
		fclose(infile);
		delete[] hdr;
		throw ST_exception("Signal::load_sph -- channel_count malformed in header.");
	}
	delete[] hdr;

	printf("format: %s\n", format);
	printf("sample_n_bytes = %d\n", sample_n_bytes);
	printf("num_channels = %d\n", num_channels);

	// Sanity checks
	if (sample_n_bytes<1 || sample_n_bytes>2) {
		throw ST_exception("Signal::load_sph -- sample number of bytes out of range; must be 1 or 2.");
	}
	bool is_ulaw = false;
	if (strcmp(format,"ulaw")==0) {
		if (sample_n_bytes!=1) {
			fclose(infile);
			throw ST_exception("Signal::load_pcm_sph -- sample_n_bytes must be 1 for ulaw.");
		}
		is_ulaw = true;
	}
	bool is_alaw = false;
	if (strcmp(format,"alaw")==0) {
		if (sample_n_bytes!=1) {
			fclose(infile);
			throw ST_exception("Signal::load_pcm_sph -- sample_n_bytes must be 1 for alaw.");
		}
		is_alaw = true;
	}
	if (strcmp(format,"pcm")==0) {
		if (sample_n_bytes!=2) {
			fclose(infile);
			throw ST_exception("Signal::load_pcm_sph -- sample_n_bytes must be 2 for pcm.");
		}
	}
	if (sampling_frequency <=0) {
		fclose(infile);
		throw ST_exception("Signal::load_pcm_sph -- sampling frequency <= 0?");
	}
   if (num_samples != 0) // Delete old data
		delete[] value;
	value = 0;

	num_samples = (st_size-num_bytes)/(sample_n_bytes*num_channels);
   if (num_samples <= 0) {
		fclose(infile);
      throw ST_exception("Signal::load_pcm_sph -- Empty input file.");
	}

	short value_short, value_short2;
	uint8_t value_byte;
   value = new float[num_samples];
	size_t num_in;
   for (int i=0; i<num_samples; i++) {
		if (is_ulaw) {
			num_in = fread(&value_byte,sizeof(uint8_t),1,infile);
			value_short = ulaw2linear(value_byte);
		} else if (is_alaw) {
			num_in = fread(&value_byte,sizeof(uint8_t),1,infile);
			value_short = alaw2linear(value_byte);
		} else
			num_in = fread(&value_short,sizeof(short),1,infile);

		if (num_channels==2) {
			if (is_ulaw) {
				num_in = fread(&value_byte,sizeof(uint8_t),1,infile);
				value_short2 = ulaw2linear(value_byte);
			} else if (is_alaw) {
				num_in = fread(&value_byte,sizeof(uint8_t),1,infile);
				value_short2 = alaw2linear(value_byte);
			} else
				num_in = fread(&value_short2,sizeof(short),1,infile);
			if (channel_num==1)
				value_short = value_short2;
		}
		value[i] = (float) value_short;
   }
	fclose(infile);

} // load_pcm_sph

static const int alawshort[] = 
{-5504, -5248, -6016, -5760, -4480, -4224, -4992, -4736,
 -7552, -7296, -8064, -7808, -6528, -6272, -7040, -6784,
 -2752, -2624, -3008, -2880, -2240, -2112, -2496, -2368,
 -3776, -3648, -4032, -3904, -3264, -3136, -3520, -3392,
 -22016, -20992, -24064, -23040, -17920, -16896, -19968, -18944,
 -30208, -29184, -32256, -31232, -26112, -25088, -28160, -27136,
 -11008, -10496, -12032, -11520, -8960, -8448, -9984, -9472,
 -15104, -14592, -16128, -15616, -13056, -12544, -14080, -13568,
 -344, -328, -376, -360, -280, -264, -312, -296,
 -472, -456, -504, -488, -408, -392, -440, -424,
 -88, -72, -120, -104, -24, -8, -56, -40,
 -216, -200, -248, -232, -152, -136, -184, -168,
 -1376, -1312, -1504, -1440, -1120, -1056, -1248, -1184,
 -1888, -1824, -2016, -1952, -1632, -1568, -1760, -1696,
 -688, -656, -752, -720, -560, -528, -624, -592,
 -944, -912, -1008, -976, -816, -784, -880, -848,
 5504, 5248, 6016, 5760, 4480, 4224, 4992, 4736,
 7552, 7296, 8064, 7808, 6528, 6272, 7040, 6784,
 2752, 2624, 3008, 2880, 2240, 2112, 2496, 2368,
 3776, 3648, 4032, 3904, 3264, 3136, 3520, 3392,
 22016, 20992, 24064, 23040, 17920, 16896, 19968, 18944,
 30208, 29184, 32256, 31232, 26112, 25088, 28160, 27136,
 11008, 10496, 12032, 11520, 8960, 8448, 9984, 9472,
 15104, 14592, 16128, 15616, 13056, 12544, 14080, 13568,
 344, 328, 376, 360, 280, 264, 312, 296,
 472, 456, 504, 488, 408, 392, 440, 424,
 88, 72, 120, 104, 24, 8, 56, 40,
 216, 200, 248, 232, 152, 136, 184, 168,
 1376, 1312, 1504, 1440, 1120, 1056, 1248, 1184,
 1888, 1824, 2016, 1952, 1632, 1568, 1760, 1696,
 688, 656, 752, 720, 560, 528, 624, 592,
 944, 912, 1008, 976, 816, 784, 880, 848};

static int alaw2linear(uint8_t alawbyte) {
	return alawshort[alawbyte];
}

// This routine converts from ulaw to 16 bit linear.
//
// Craig Reese: IDA/Supercomputing Research Center
// 29 September 1989
//
// References:
// 1) CCITT Recommendation G.711  (very difficult to follow)
// 2) MIL-STD-188-113,"Interoperability and Performance Standards
//     for Analog-to_Digital Conversion Techniques,"
//     17 February 1987
//
// Input: 8 bit ulaw sample
// Output: signed 16 bit linear sample
//
static int ulaw2linear(uint8_t ulawbyte) {
  static int exp_lut[8] = {0,132,396,924,1980,4092,8316,16764};
  int sign, exponent, mantissa, sample;

  ulawbyte = ~ulawbyte;
  sign = (ulawbyte & 0x80);
  exponent = (ulawbyte >> 4) & 0x07;
  mantissa = ulawbyte & 0x0F;
  sample = exp_lut[exponent] + (mantissa << (exponent + 3));
  if (sample==0) // LLSpeech compatibility
	  sample = 4;
  if (sign != 0) sample = -sample;

  return(sample);
}

void Signal::load_pcm_wav (string filename, bool sum_channels)
{
	FILE  *infile;
	char temp[5];
	unsigned int spacer;

	// TODO: Code is not very portable -- depends on the size of int, char, etc.
	unsigned st_size = file_len(filename.c_str());
   infile = fopen(filename.c_str(), "rb");
   if (infile == NULL)
      throw ST_exception("Signal::load_pcm_wav -- Unable to open input file.");

	// Read RIFF
	if (fread(temp, sizeof(unsigned char), 4, infile)!=4) {
		fclose(infile);
      throw ST_exception("Signal::load_pcm_wav -- unexpected end of file.");
	}

	temp[4] = '\0';
	if (strcmp(temp,"RIFF")!=0) {
		fclose(infile);
      throw ST_exception("Signal::load_pcm_wav -- input file is not in RIFF format.");
	}
  
	// Read spacer
	if (fread(&spacer, sizeof(int), 1, infile)!=1) {
		fclose(infile);
      throw ST_exception("Signal::load_pcm_wav -- unexpected end of file.");
	}

	// Read WAVE
	if (fread(temp, sizeof(unsigned char), 4, infile)!=4) {
		fclose(infile);
      throw ST_exception("Signal::load_pcm_wav -- unexpected end of file.");
	}
	temp[4] = '\0';
	if (strcmp(temp,"WAVE")!=0) {
		fclose(infile);
      throw ST_exception("Signal::load_pcm_wav -- input file is not a WAVE sound file.");
	}

	// Read format chunk
	if (fread(temp, sizeof(unsigned char), 4, infile)!=4) {
		fclose(infile);
      throw ST_exception("Signal::load_pcm_wav -- unexpected end of file.");
	}
	if (strcmp(temp,"fmt ")!=0) {	
		fclose(infile);
      throw ST_exception("Signal::load_pcm_wav -- fmt not found.");
	}

	// Read spacer
	if (fread(&spacer, sizeof(int), 1, infile)!=1) {
		fclose(infile);
      throw ST_exception("Signal::load_pcm_wav -- unexpected end of file.");
	}

	// Get format data
	unsigned short cat, num_channels, bits_per_sample, block_align;
	unsigned int sampling_rate, avg_bytes_per_sec;
	if (fread(&cat, sizeof(unsigned short),1,infile)!=1 || fread(&num_channels, sizeof(unsigned short),1,infile)!=1 ||
		 fread(&sampling_rate, sizeof(unsigned int),1,infile)!=1 || fread(&avg_bytes_per_sec, sizeof(unsigned int),1,infile)!=1 ||
		 fread(&block_align, sizeof(unsigned short),1,infile)!=1) {
		fclose(infile);
		throw ST_exception("Signal::load_pcm_wav -- unexpected end of file.");
	}
	if (cat!=1) {
		fclose(infile);
		throw ST_exception("Signal::load_pcm_wav -- wav file is not PCM.");
	}
	if (num_channels!=1 && !sum_channels) {
		fclose(infile);
		throw ST_exception("Signal::load_pcm_wav -- only single channel data allowed.");
	}
	if (fread(&bits_per_sample, sizeof(unsigned short), 1, infile)!=1) {
		fclose(infile);
		throw ST_exception("Signal::load_pcm_wav -- unexpected end of file.");
	}
	if ((num_channels*sampling_rate*2 != avg_bytes_per_sec) || (bits_per_sample != 16)) {
		fclose(infile);
		throw ST_exception("Signal::load_pcm_wav -- only 16 bit linear PCM is supported.");
	}
  
	// Read data chunk
	int i, j;
	unsigned int nsize;
	if (fread(temp, sizeof(unsigned char), 4, infile)!=4) {
		fclose(infile);
		throw ST_exception("Signal::load_pcm_wav -- unexpected end of file.");
	}
	while (strcmp(temp,"data")!=0) {	
		for (i=0; i<3; i++)
			temp[i] = temp[i+1];
		if (fread(&temp[3], sizeof(unsigned char), 1, infile)!=1) {
			fclose(infile);
			throw ST_exception("Signal::load_pcm_wav -- unexpected EOF.\n");
		}
	}
	if (fread(&nsize, sizeof(unsigned int), 1, infile)!=1) { 
		fclose(infile);
		throw ST_exception("Signal::load_pcm_wav -- unexpected end of file.");
	}
   // nsize = num_channels*num_samples*(bits_per_sample/2)
   nsize /= num_channels*2;

	// Sanity checks, memory allocation, read in signal data
	sampling_frequency = (float) sampling_rate;
	if (sampling_frequency <=0) {
		fclose(infile);
		throw ST_exception("Signal::load_pcm_wav -- sampling frequency <= 0?");
	}
   if (num_samples != 0) // Delete old data
		delete[] value;
	value = 0;
	num_samples = nsize;
   if (num_samples <= 0) {
		fclose(infile);
		throw ST_exception("Signal::load_pcm_wav -- Empty input file.");
	}
   value = new float[num_samples];

	short sample;
	for (i=0; i<((int) nsize); i++) {
      value[i] = 0;
      for (j=0; j<num_channels; j++) {
		   if (fread(&sample, sizeof(short), 1, infile)!=1) {
			   fclose(infile);
			   throw ST_exception("Signal::load_pcm_wav -- corrupt wav file?\n");
		   }
         value[i] += (float) sample;
      }
	}
	fclose(infile);
}

int Signal::len ()
{
	return num_samples;
}

void Signal::normalize (void)
{
	float max = 0;
   int i;

	if (num_samples==0)
		return;

	for (i=0; i<num_samples; i++) {
		if (fabs(value[i]) > max) 
			max = fabs(value[i]);
	}
	if (max == 0) {
		max = 1;
	}

   max = (float) (32767.0/max);

   for (i=0; i<num_samples; i++)
		value[i] *= max;

} // Signal::normalize

void Signal::preemphasis (float pre_coeff)
{
float		prev_sample=0, prev_sample_temp;
int		i;

	if (num_samples == 0) 
		throw ST_exception("No signal to preemphasize.");

	for (i=0; i < num_samples; i++) {
		prev_sample_temp = value[i];
		value[i] -= pre_coeff*prev_sample;
		prev_sample = prev_sample_temp;
	}

} // Pre-emphasis

void Signal::remove_mean (void)
{
float		mean;
int		i;

	if (num_samples == 0) 
		throw ST_exception("No signal.");

	for (i=0, mean=0; i < num_samples; i++)
		mean += value[i];

	mean /= num_samples;

	for (i=0; i < num_samples; i++)
		value[i] -= mean;

}  // remove_mean 

void Signal::resample_8k (void)
{
   if (sampling_frequency==8000.0)
      return;

   if (!resampler_loaded)
      throw ST_exception("must load resampler using resampler_init");

   // Resample to 8k : filters are hard coded
   if (sampling_frequency==11025.0 || sampling_frequency==22050.0 || sampling_frequency==44100.0) {
      if (sampling_frequency==22050.0) {
         fir_updn(resampler_filters["h_dn2"].data, resampler_filters["h_dn2"].len, 1, 2);
      } else if (sampling_frequency==44100.0) {
         fir_updn(resampler_filters["h1_dn4"].data, resampler_filters["h1_dn4"].len, 1, 2);
         fir_updn(resampler_filters["h2_dn4"].data, resampler_filters["h2_dn4"].len, 1, 2);
      }
      fir_zphase(resampler_filters["low_pass_11k_8k"].data, resampler_filters["low_pass_11k_8k"].len);
      fir_updn(resampler_filters["interp_11k_8k"].data, resampler_filters["interp_11k_8k"].len, 320, 441);
   } else if (sampling_frequency==16000.0 || sampling_frequency==32000.0 || sampling_frequency==48000.0) {
      if (sampling_frequency==16000.0) {
         fir_updn(resampler_filters["h_dn2"].data, resampler_filters["h_dn2"].len, 1, 2);
      } else if (sampling_frequency==32000.0) {
         fir_updn(resampler_filters["h1_dn4"].data, resampler_filters["h1_dn4"].len, 1, 2);
         fir_updn(resampler_filters["h2_dn4"].data, resampler_filters["h2_dn4"].len, 1, 2);
      } else if (sampling_frequency==48000.0) {
         fir_updn(resampler_filters["h1_dn6"].data, resampler_filters["h1_dn6"].len, 1, 2);
         fir_updn(resampler_filters["h2_dn6"].data, resampler_filters["h2_dn6"].len, 1, 3);
      }
   } else {
      throw ST_exception("unsupported input sampling rate for resample to 8 kHz");
   } 
   sampling_frequency = 8000;

}

void Signal::resample_16k (void)
{
   if (sampling_frequency==16000.0)
      return;

   if (!resampler_loaded)
      throw ST_exception("must load resampler using resampler_init");

   // Resample to 16k : filters are hard coded
   if (sampling_frequency==22050.0 || sampling_frequency==44100.0) {
      if (sampling_frequency==44100.0)
         fir_updn(resampler_filters["h_dn2"].data, resampler_filters["h_dn2"].len, 1, 2);
      fir_zphase(resampler_filters["low_pass_11k_8k"].data, resampler_filters["low_pass_11k_8k"].len);
      fir_updn(resampler_filters["interp_11k_8k"].data, resampler_filters["interp_11k_8k"].len, 320, 441);
   } else if (sampling_frequency==16000.0 || sampling_frequency==32000.0 || sampling_frequency==48000.0) {
      if (sampling_frequency==32000.0) {
         fir_updn(resampler_filters["h_dn2"].data, resampler_filters["h_dn2"].len, 1, 2);
      } else if (sampling_frequency==48000.0) {
         fir_updn(resampler_filters["h_dn3"].data, resampler_filters["h_dn3"].len, 1, 3);
      } 
   } else {
      throw ST_exception("unsupported input sampling rate for resample to 16 kHz");
   } 
   sampling_frequency = 16000;
}

float Signal::sampling_freq (void)
{
	return sampling_frequency;
}

void Signal::save_pcm_wav (string filename, bool scale)
{
	FILE  *outfile;
	char temp[5];
	unsigned int totbytes;

	// TODO: Code is not very portable -- depends on the size of int, char, etc.
   outfile = fopen(filename.c_str(), "wb");
   if (outfile == NULL)
      throw ST_exception("Signal::save_pcm_wav -- Unable to open output file.");
	
	// Write RIFF
   strcpy(temp, "RIFF");
	if (fwrite(temp, sizeof(unsigned char), 4, outfile)!=4) {
		fclose(outfile);
      throw ST_exception("Signal::save_pcm_wav -- unexpected end of file.");
	}

	// Write total length of chunk
   totbytes = 2*num_samples + 36;  // assume 2 bytes per sample
   // 36 is the length of the preamble before the data
	if (fwrite(&totbytes, sizeof(int), 1, outfile)!=1) {
		fclose(outfile);
      throw ST_exception("Signal::save_pcm_wav -- unexpected end of file.");
	}
   totbytes -= 36;

	// Write WAVE
   strcpy(temp, "WAVE");
	if (fwrite(temp, sizeof(unsigned char), 4, outfile)!=4) {
		fclose(outfile);
      throw ST_exception("Signal::save_pcm_wav -- unexpected end of file.");
	}

	// write format chunk
   strcpy(temp, "fmt ");
	if (fwrite(temp, sizeof(unsigned char), 4, outfile)!=4) {
		fclose(outfile);
      throw ST_exception("Signal::save_pcm_wav -- unexpected end of file.");
	}

	// Write format
   unsigned int numbits = 16;  // for PCM
	if (fwrite(&numbits, sizeof(int), 1, outfile)!=1) {
		fclose(outfile);
      throw ST_exception("Signal::save_pcm_wav -- unexpected end of file.");
	}

	// Write format data
	unsigned short cat, num_channels, bits_per_sample, block_align;
	unsigned int sampling_rate, avg_bytes_per_sec;
   cat = 1;  // PCM
   num_channels = 1;
   sampling_rate = (unsigned int) sampling_frequency;
   bits_per_sample = 16;
   avg_bytes_per_sec = 2*sampling_rate;
   block_align = 2;
   bits_per_sample = 16;

	if (fwrite(&cat, sizeof(unsigned short),1,outfile)!=1 || fwrite(&num_channels, sizeof(unsigned short),1,outfile)!=1 ||
		 fwrite(&sampling_rate, sizeof(unsigned int),1,outfile)!=1 || fwrite(&avg_bytes_per_sec, sizeof(unsigned int),1,outfile)!=1 ||
		 fwrite(&block_align, sizeof(unsigned short),1,outfile)!=1 || fwrite(&bits_per_sample, sizeof(unsigned short), 1, outfile)!=1) {
		fclose(outfile);
		throw ST_exception("Signal::save_pcm_wav -- unexpected end of file.");
	}

	// Write data chunk
   strcpy(temp, "data");
	if (fwrite(temp, sizeof(unsigned char), 4, outfile)!=4 || fwrite(&totbytes, sizeof(unsigned int), 1, outfile)!=1) {
		fclose(outfile);
		throw ST_exception("Signal::save_pcm_wav -- unexpected end of file.");
	}

   // Determine max
   float max = 0.0;
   int i;
   if (scale) {
      for (i=0; i<num_samples; i++) {
		   if (fabs(value[i]) > max) 
			   max = fabs(value[i]);
	   }
	   if (max == 0) {
		   max = 1;
	   }
      max = (float) (32767.0/max);
   } else {
      max = 1.0;
   }

	short sample;
	for (i=0; i<num_samples; i++) {
      sample = (short)(max*value[i]);
		if (fwrite(&sample, sizeof(short), 1, outfile)!=1) {
			fclose(outfile);
			throw ST_exception("Signal::save_pcm_wav -- unexpected end of file -- disk full?\n");
		}
	}
	fclose(outfile);

}

void Signal::save_raw_float (string filename)
{
	FILE 	*outfile;
	
	outfile = fopen(filename.c_str(), "wb");
	if (outfile == NULL)
		throw ST_exception("Unable to open output file.\n");

	fwrite (value, sizeof(float), num_samples, outfile);
	fclose(outfile);

} // Signal::save

void Signal::save_raw_short (string filename, bool clip, bool scale)
{
	FILE 	*outfile;
	float max = 0, value_tmp;
	short value_short;
	int i;

	if (num_samples==0)
		throw ST_exception("Cannot save empty signal");

	if ((!clip) && (!scale)) 
		throw ST_exception("Must set either clipping or scaling when saving a signal to a short format.");

	for (i=0; i<num_samples; i++) {
		if (fabs(value[i]) > max) 
			max = fabs(value[i]);
	}
	if (max == 0) {
		max = 1;
	}

   max = (float) (32767.0/max);
	outfile = fopen(filename.c_str(), "wb");
	if (outfile == NULL)
		throw ST_exception("Unable to open output file.\n");

	for (i=0; i<num_samples; i++) {
		value_tmp = value[i];
		if (scale) {
			value_tmp *= max;
		}
		if (clip) {
			if (value_tmp > 32767) {
				value_tmp = 32767;
			}
			if (value_tmp < -32768) {
				value_tmp = -32768;
			}
		}
		// value_tmp = rintf(value_tmp);
      value_tmp = (float) ((value_tmp>0.0) ? floor(value_tmp+0.5) : ceil(value_tmp-0.5));
		value_short = (short) value_tmp; 
		fwrite (&value_short, sizeof(short), 1, outfile);
	}

	fclose(outfile);

} // Signal::save

