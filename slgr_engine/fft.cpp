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
// FFT routines
//
// Basic ideas from LLSpeech, but rewritten significantly for efficiency and
// clarity

// wmc, 10/20/10

#include "stdafx.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "speech_tools.h"

#define PI 3.1415926535897931159979634685

mag_fft::mag_fft (int n) 
{
   if (n<=0)  
      throw ST_exception("RFFT::RFFT -- length must be > 0.");
	n_fft = 1;
	log2n = 0;
   while (n_fft < n) {
      n_fft *= 2;
      log2n++;
   }
	n_input = n;
	n2_fft = n_fft/2;
	r = new RFFT(n_fft);
	in = new vec<REAL_FFT>(n_fft);
	out_r = new vec<REAL_FFT>(n2_fft);
	out_i = new vec<REAL_FFT>(n2_fft);

}

mag_fft::~mag_fft()
{
	delete r;
	delete in;
	delete out_r;
	delete out_i;
}

vec<REAL_FFT> mag_fft::calc (vec<REAL_FEAT> &x)
{
	if (x.len!=n_input)
		throw ST_exception("mag_fft::calc -- input length does not match FFT initialization.");

	// Real FFT calculation
	int i;
	for (i=0; i<n_input; i++)
		in->data[i] = x.data[i];
	for (i=n_input; i<n_fft; i++)
		in->data[i] = 0;
	r->rfft(*in, *out_r, *out_i);

	// Magnitude calculation
	vec<REAL_FFT> x_mag(n2_fft+1);
	x_mag.data[0] = out_r->data[0]*out_r->data[0];
	x_mag.data[n2_fft] = out_i->data[0]*out_i->data[0];
	for (i=1; i<n2_fft; i++)
		x_mag.data[i] = out_r->data[i]*out_r->data[i]+out_i->data[i]*out_i->data[i];

	return x_mag;

}

int mag_fft::len ()
{
	return n_fft;
}

RFFT::RFFT (int n) : f(n/2) {
   if (n<=0)  
      throw ST_exception("RFFT::RFFT -- length must be > 0.");
   log2n = 0;
   int n_compare = 1;
   while (n_compare < n) {
      n_compare *= 2;
      log2n++;
   }
   if (n_compare != n)  // not a power of two
      throw ST_exception("RFFT::RFFT -- length must be power of 2.");
	RFFT::n = n;
	n2 = n/2;

   // Cosine and sine table
   costbl = new REAL_FFT[n2+1];
   sintbl = new REAL_FFT[n2+1];
   for (int i=0; i<=n2; i++) {
      costbl[i] = cos((i*2*PI)/n);
      sintbl[i] = sin((i*2*PI)/n);
   }
}

RFFT::~RFFT () {
	if (costbl!=0)
		delete[] costbl;
	if (sintbl!=0)
		delete[] sintbl;
}

//
// FFT of a real vector
// About half the computation of doing a full complex FFT
// Scaling for this version agrees with Matlab
//
void RFFT::rfft(vec<REAL_FFT> &x, vec<REAL_FFT> &yr, vec<REAL_FFT> &yi)
{
	if (x.len!=n)
		throw ST_exception("RFFT::rfft -- input length does not match FFT initialization.");
	if ((yr.len!=n2) || (yi.len!=n2))
		throw ST_exception("RFFT::rfft -- output length must be half input length.");
		
	// First step is to do the FFT of x(0:2:end) + J*x(1:2:end)
	int i;
	for (i=0; i<n2; i++) {
		yr.data[i] = x.data[2*i];
		yi.data[i] = x.data[2*i+1];
	}
	f.fft(yr, yi);

	// Now combine with a butterfly type process
	// Store both the 0 and 0.5 frequency values in Y(0)
	int j1, j2;
	REAL_FFT tr1, ti1, tr2, ti2, xr, xi, c, s;
	xr = yr.data[0]-yi.data[0];
	yr.data[0] = yr.data[0]+yi.data[0];
	yi.data[0] = xr;
	for (j1=1; j1<(n/4)+1; j1++) {
		j2 = n/2-j1;
		c = costbl[j1];
		s = sintbl[j1];

		tr1 = yr.data[j2]+yr.data[j1];
		ti1 = yi.data[j1]-yi.data[j2];
		tr2 = yi.data[j1]+yi.data[j2];
		ti2 = yr.data[j2]-yr.data[j1];
		xr = c*tr2+s*ti2;
		xi = c*ti2-s*tr2;

		yr.data[j1] = 0.5*(tr1+xr);
		yi.data[j1] = 0.5*(ti1+xi);
		yr.data[j2] = 0.5*(tr1-xr);
		yi.data[j2] = 0.5*(xi-ti1);
	}

}

// This is adapted from Sid Burrus's code which is
// quite similar to the Markel/Gray code (pg. 160)
//
// This is standard radix 2 FFT class for complex data
//
FFT_DIF::FFT_DIF (int n)
{
   int i;

   if (n<=0)  
		throw ST_exception("FFT_DIF::FFT_DIF -- length must be > 0.");

   log2n = 0;
   int n_compare = 1;
   while (n_compare < n) {
      n_compare *= 2;
      log2n++;
   }
   if (n_compare != n)  // not a power of two
      throw ST_exception("FFT_DIF::FFT_DIF -- length must be power of 2.");
	FFT_DIF::n = n;
   int n2 = n/2;

   // Cosine and sine table
   costbl = new REAL_FFT[n2+1];
   sintbl = new REAL_FFT[n2+1];
   for (i=0; i<=n2; i++) {
      costbl[i] = cos((i*2*PI)/n);
      sintbl[i] = sin((i*2*PI)/n);
   }

   // i -> bit_reverse[i] is bit reverse
   bit_reverse = new int[n];
   for (i=0; i<n; i++)
      bit_reverse[i] = i;
   int j = 1, xt, k;
   for (i=1; i<=(n-1); i++) {
		if (i<j) {
			xt = bit_reverse[j-1];
			bit_reverse[j-1] = bit_reverse[i-1];
			bit_reverse[i-1] = xt;
      }
      k = n/2;
      while (k<j) {
			j = j-k;
			k = k/2;
      }
      j = j+k;
	}
}

FFT_DIF::~FFT_DIF () {
   if (costbl!=0)
      delete[] costbl;
   if (sintbl!=0)
      delete[] sintbl;
   if (bit_reverse!=0)
      delete[] bit_reverse;
}

void FFT_DIF::fft (vec<REAL_FFT> &in_real, vec<REAL_FFT> &in_imag) {

	if ((in_real.len!=n) || (in_imag.len!=n))
		throw ST_exception("FFT_DIF::fft -- FFT length does not match initialization.");

   int i, j, k, l, itable, n1;
   REAL_FFT xt, yt, c, s;
   REAL_FFT *x = in_real.data;
   REAL_FFT *y = in_imag.data;

   // Basic FFT loop
   int skip = 1;
   int n2 = n;
   for (k=0; k<log2n; k++) {
		n1 = n2;
      n2 = n2/2;
      itable = 0;
      for (j=0; j<n2; j++) {
			c = costbl[itable];
			s = sintbl[itable];
			for (i=j; i<n; i+=n1) {
				l = i+n2;
				xt = x[i]-x[l];
				yt = y[i]-y[l];
				x[i] += x[l];
				y[i] += y[l];
				x[l] = c*xt+s*yt;
				y[l] = c*yt-s*xt;
			}
			itable += skip;
      }
      skip *= 2;
   }

   // Bit reversal
   for (k=0; k<n; k++) {
		j = bit_reverse[k];
      if (j>k) {
			xt = x[j];
			x[j] = x[k];
			x[k] = xt;
			yt = y[j];
			y[j] = y[k];
			y[k] = yt;
      }
   }

}
