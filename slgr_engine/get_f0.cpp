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
// Classic get_f0 from Entropic for the Signal class
//

// wmc: 9/2014-11/2014, 1st version

//
// Standard headers
//
#include "stdafx.h"
#include "speech_tools.h"
#include "../get_f0_lib/f0.h"
// #include "../fftw_libs/fftw3.h"
#include "fftw3.h"

class mag_fftw {
   fftw_plan p;
   int n, n2;
   double *in;
   fftw_complex *out;
public:
   mag_fftw(int len);
   ~mag_fftw(void);
   int len(void);
   vec<float> calc(vec<float> &x);
};

// 
// Interface to legacy ESPS get_f0 function
//
vec<float> Signal::get_f0 (float min_f0, float max_f0, float win_dur, float &frame_step)
{
   int start_sample = 0;
   int end_sample = num_samples-1;

   // Call legacy ESPS function
   vector<float> f0_store;
   get_f0_session *session = init_get_f0();
   session->par->min_f0 = min_f0;
   session->par->max_f0 = max_f0;
   session->par->wind_dur = win_dur;
   session->par->frame_step = frame_step;
   int status = get_f0_esps(*this, session, start_sample, end_sample, f0_store);
   frame_step = session->par->frame_step;
   close_get_f0(session);

	if (status!=0)
		throw ST_exception("Error in get_f0");

   // Copy result if the process completed successfully
   int num_f0 = 0;
   if (status==0) // success == 0
      num_f0 = f0_store.size();
   vec<float> f0(num_f0);
   if (status == 0) { 
      int i;
      vector<float>::iterator it;
      num_f0 = f0_store.size();
      for (it=f0_store.begin(), i=0; it!=f0_store.end(); it++, i++) 
         f0.data[i] = *it;
   }

   return f0;

} // get_f0

mag_fftw::mag_fftw (int in_len) 
{
   n = in_len;
   n2 = in_len/2+1;  // one more element than mag_fft

   // Use fftw malloc to ensure proper memory alignment for arrays for SSD operations
   in = (double *) fftw_malloc(n*sizeof(double));
   if (in==NULL)
      throw ST_exception("Cannot initialize fftw.");
   out = (fftw_complex *) fftw_malloc(n2*sizeof(fftw_complex));
   if (out==NULL) {
      fftw_free(in);
      throw ST_exception("Cannot initialize fftw.");
   }

   // initialize plan
   p = fftw_plan_dft_r2c_1d (n, in, out, FFTW_ESTIMATE);
   if (p==NULL) {
      fftw_free(in);
      fftw_free(out);
      throw ST_exception("Cannot initialize fftw.");
   }
}

mag_fftw::~mag_fftw ()
{
   fftw_free(in);
   fftw_free(out);
   fftw_destroy_plan(p);
}

vec<float> mag_fftw::calc(vec<float> &x)
{
   int i;

	if (x.len!=n)
		throw ST_exception("mag_fftw::calc -- input length does not match FFT initialization.");

   for (i=0; i<x.len; i++)
      in[i] = x.data[i];

   fftw_execute(p);

   vec<float> y(n2);
   for (i=0; i<n2; i++) 
      y.data[i] = out[i][0]*out[i][0]+out[i][1]*out[i][1];

   return y;
}

int mag_fftw::len(void)
{
   return n;
}

