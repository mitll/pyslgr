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

