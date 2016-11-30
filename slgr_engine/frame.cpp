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

