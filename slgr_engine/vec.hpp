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
// Very simple vector class
// 
//
#pragma once

using namespace std;

#include <math.h>
#include <stdio.h>
#include <string>

#include "st_exception.h"

template<class T>
class vec {
public:
	T *data;
	int len;
	vec(T *in, int len);
	vec(int len);
	vec();
	~vec();
	vec& operator= (const vec<T>& old);
   vec (const vec<T> &in);
	double mean ();
	double std (double std_floor);
	void add (vec<T> &vec2);
   void add_diag (const T val);
	void append (T *in, int len);
   void constant (const T val);
	void divide (vec<T> &vec2);
	T ip (vec<T> &vec2);
	void outer_add (vec<T> &vec2);
	void outer_sub (vec<T> &vec2);
   void save_raw (string &outfile);
	void scale (vec<T> &vec2);
	void scale (T a);
	void subtract (vec<T> &vec2);
	void subtract (T v);
};

template<class T>
vec<T>::vec (T *in, int len)
{
	if (len < 0) 
		throw ST_exception("vec::vec -- length must be >= 0.");
	vec::len = len;
	vec::data = new T[len];
	for (int i=0; i<len; i++)
		data[i] = in[i];
}

template<class T>
vec<T>::vec ()
{
   len = 0;
   data = 0;
}

template<class T>
vec<T>::vec (int len)
{
	if (len < 0) 
		throw ST_exception("vec::vec -- length must be >= 0.");
	vec::len = len;
	if (len == 0) {
		vec::data = 0;
		return;
	}
	vec::data = new T[len];
	for (int i=0; i<len; i++)
		data[i] = 0;
}

template<class T>
vec<T>::~vec()
{
	if (data!=0)
		delete[] data;
}

template<class T>
vec<T>::vec (const vec &in)
{
	len = in.len;
	if (len == 0) {
		data = 0;
	} else {
		data = new T[len];
		for (int i=0; i<len; i++)
			data[i] = in.data[i];
	}
}

template<class T>
vec<T>& vec<T>::operator= (const vec& old)
{
	if (len != old.len) {
		if (data!=0)
			delete[] data;
		len = old.len;
		data = new T[len];
	}
	for (int i=0; i<len; i++)
		data[i] = old.data[i];
	return *this;
}

template<class T>
inline void vec<T>::add (vec &vec2)
{
	if (len != vec2.len)
		throw ST_exception("vec::add -- vectors must be the same length.");
	for (int i=0; i<len; i++) 
		data[i] += vec2.data[i];
}

template<class T>
void vec<T>::add_diag (const T val)
{
   int l2 = (int) sqrt((float) len);

   for (int i=0; i<l2; i++)
      data[i+i*l2] += val;
}


template<class T>
void vec<T>::append (T *in, int len) 
{
	if (len < 0) 
		throw ST_exception("vec::append -- length must be >= 0.");
	int new_len = (vec::len) + len;
	T *new_data = new T[new_len];
	int i = 0;
	for (; i<vec::len; i++)
		new_data[i] = data[i];
	for (int j=0; j<len; j++) 
		new_data[i+j] = in[j];
	if (data != 0)
		delete[] data;
	vec::len = new_len;
	vec::data = new_data;
}

template<class T>
double vec<T>::mean ()
{
	if (len==0)
		return 0.0;
	double sum = 0;
	for (int i=0; i<len; i++)
		sum += data[i];
	sum /= len;
	return sum;
}

template<class T>
void vec<T>::constant (const T val)
{
   for (int i=0; i<len; i++)
      data[i] = val;
}


template<class T>
void vec<T>::divide (vec<T> &vec2)
{
	if (len==0)
		return;
	for (int i=0; i<len; i++) 
		data[i] /= vec2.data[i];
}

template<class T>
T vec<T>::ip (vec &vec2)
{
	if (len != vec2.len)
		throw ST_exception("vec::ip -- vectors must be the same length.");
	double sum = 0;
	for (int i=0; i<len; i++)
		sum += data[i]*vec2.data[i];
	return ((T) sum);
}

template<class T>
inline void vec<T>::outer_add (vec &vec2)
{
	if (len != (vec2.len*vec2.len))
		throw ST_exception("vec::outer_add -- output vector must be the square of the input dimension.");

   int i, j, k;
	for (i=0, k=0; i<vec2.len; i++) {
      for (j=0; j<vec2.len; j++, k++) {
         data[k] += vec2.data[i]*vec2.data[j];
      }
   }
}

template<class T>
inline void vec<T>::outer_sub (vec &vec2)
{
	if (len != (vec2.len*vec2.len))
		throw ST_exception("vec::outer_add -- output vector must be the square of the input dimension.");

   int i, j, k;
	for (i=0, k=0; i<vec2.len; i++) {
      for (j=0; j<vec2.len; j++, k++) {
         data[k] -= vec2.data[i]*vec2.data[j];
      }
   }
}

template<class T>
void vec<T>::save_raw (string &outfile) 
{
   FILE *out = fopen(outfile.c_str(), "wb");
   fwrite(data, sizeof(T), len, out);
   fclose(out);
}

template<class T>
inline void vec<T>::scale (vec &vec2)
{
	if (len != vec2.len)
		throw ST_exception("vec::scale -- vectors must be the same length.");
	for (int i=0; i<len; i++) 
		data[i] *= vec2.data[i];
}

template<class T>
inline void vec<T>::scale (T a)
{
	for (int i=0; i<len; i++) 
		data[i] *= a;
}

template<class T>
double vec<T>::std (double std_floor)
{
	if (len==0)
		return 1.0;
	double sum = 0;
	double sum2 = 0;
	for (int i=0; i<len; i++) {
		sum += data[i];
		sum2 += data[i]*data[i];
	}
	sum /= (double) len;
	sum2 /= (double) len;
	sum2 = sqrt(sum2-sum*sum);
	if (sum2 < std_floor)
		sum2 = std_floor;
	return sum2;
}

template<class T>
inline void vec<T>::subtract (T v)
{
	if (len==0)
		return;
	for (int i=0; i<len; i++) 
		data[i] -= v;
}

template<class T>
inline void vec<T>::subtract (vec &vec2)
{
	if (len != vec2.len)
		throw ST_exception("vec::subtract -- vectors must be the same length.");
	for (int i=0; i<len; i++) 
		data[i] -= vec2.data[i];
}

