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

#include "stdafx.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <cstdint>

#include "speech_tools.h"

//
// Prototypes for local functions
//
static void	create_hamming_window (REAL_FEAT *window, int frame_size);
static double uni (int32_t &s1, int32_t &s2);

//
// Methods
//
Frame::Frame (int n, int inc) : x(n) , hw(n) {
	if (n<=0) 
		throw ST_exception("Frame::Frame -- frame size must be > 0.");
	if (inc<=0)
		throw ST_exception("Frame::Frame -- frame increment must be > 0.");
	frame_inc = inc;
	frame_loaded = false;
	create_hamming_window(hw.data, n);
	s1 = 55555;
	s2 = 99999;
		
}

Frame::~Frame () {
}

void Frame::dither(float scale)
{
	if (!frame_loaded) 
		throw ST_exception("Frame::dither -- frame not loaded.");

	for (int i=0; i<x.len; i++)
		x.data[i] += scale*(2*uni(s1, s2)-1);

}

REAL_FEAT Frame::energy () 
{
	if (!frame_loaded) 
		throw ST_exception("Frame::hamming_window -- frame not loaded.");

	REAL_FEAT e = 0;
	for (int i=0; i<x.len; i++)
		e += x.data[i]*x.data[i];
	e = 10*log10((double)(e+1.0e-20))-40;
	return e;

}

void Frame::get_frame (Signal &z, int frame_num)
{
	int offset = frame_num*frame_inc;

	if ((frame_num<0) || ((offset+x.len-1) >= z.len()))
		throw ST_exception("Frame::get_frame -- out of bounds access of signal."); 

	for (int i=0; i<x.len; i++)
		x.data[i] = z.value[offset+i];

	frame_loaded = true;

}

int Frame::get_num_frames (Signal &z)
{
	int n = z.len();
	n = (n-x.len)/frame_inc;
	n += 1;
	return n;
}

void Frame::hamming_window ()
{
	if (!frame_loaded) 
		throw ST_exception("Frame::hamming_window -- frame not loaded.");
	x.scale(hw);
}

void Frame::load (float *input) {
	for (int i=0; i<x.len; i++)
		x.data[i] = input[i];
	frame_loaded = true;
}

void Frame::rm_dc ()
{
	if (!frame_loaded) 
		throw ST_exception("Frame::rm_dc -- frame not loaded.");
	x.subtract(x.mean());
}

vec<REAL_FEAT>& Frame::vector ()
{
	return x;
}

//
// Local functions
//
static void create_hamming_window (REAL_FEAT *window, int frame_size) 
{
	int i;
	float tmp = 6.283185308/(frame_size-1);  // slightly less accurate version from LLSpeech

	for (i=0; i<frame_size; i++)
		window[i] = 0.54 - (float) 0.46*cos((double)i*tmp);

	// REAL_FEAT pi;
	// pi = (REAL_FEAT) atan(1.0)*4;
   // for (i=0; i<frame_size; i++)
	//	window[i] = (REAL_FEAT) (0.54-(0.46*cos(((2*pi)/(frame_size-1))*i)));

} // hamming_window

//
// Not very portable -- all ints should be 32-bits
// 
static double uni (int32_t &s1, int32_t &s2)
{
  double factor = 1.0/2147483563.0;
  int32_t k, z;

  // printf("s1 = %ld , s2 = %ld\n", s1, s2);

  k = s1/53668;
  s1 = 40014*(s1%53668)-k*12211;
  if (s1<0) 
	  s1 += 2147483563;
  k = s2/52774;
  s2 = 40692*(s2%52774)-k*3791;
  if (s2 < 0) 
	  s2 += 2147483399;
  
  z = (s1-2147483563) + s2;
  if (z < 1) z += 2147483562;

  return factor*((double) z);

}

