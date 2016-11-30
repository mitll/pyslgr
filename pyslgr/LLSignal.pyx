#
# Implements Cython class interface to SLGR C++ Signal class
#
#
# Copyright 2016 MIT Lincoln Laboratory, Massachusetts Institute of Technology
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not use these files except in compliance with
# the License.
#
# You may obtain a copy of the License at
#
#     http:#www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on
# an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the
# specific language governing permissions and limitations under the License.
#

from libcpp.string cimport string
from libcpp cimport bool
import numpy as np
import os
import sys
from Signal cimport Signal
from vec cimport vec
from cython.operator cimport dereference as deref

cdef class LLSignal:
    """Class to contain and process 1-dimensional signals--typically speech or audio.

    LLSignal() -- empty signal with zero samples.

    """
    def __cinit__(self, object data_mv=np.zeros(1,dtype=np.int16), int sampling_frequency = -1):

        cdef short [:] data_short
        cdef double [:] data_double

        if sampling_frequency == -1:
            self._sigPtr = new Signal()
        elif type(data_mv) is np.ndarray:
            if data_mv.dtype.type is np.short:
                data_short = data_mv
                self._sigPtr = new Signal(&data_short[0], <int> data_mv.size, <int> sampling_frequency)
            elif data_mv.dtype.type is np.float64:
                data_double = data_mv
                self._sigPtr = new Signal(&data_double[0], <int> data_mv.size, <int> sampling_frequency)
            else:
                raise TypeError('input type must be np.float64 or np.short')
        else:
            raise TypeError('input type must be np.float64 or np.short')
        if self._sigPtr == NULL:
            raise MemoryError()
        myModule = sys.modules['pyslgr']
        absPath = os.path.abspath(myModule.__file__)
        pkgDir = os.path.dirname(absPath)
        filterDir = os.path.join(pkgDir,'filters')
        self._sigPtr.init_resampler(filterDir)

    def __copy__(self):
        new_signal = type(self)(self[:], <int> self.sampling_frequency)
        return new_signal

    def __dealloc__(self):
        if self._sigPtr != NULL:
            # print "Destroying signal from LLSignal.pyx"
            del self._sigPtr
    
    def __getitem__(self, object index):
        if isinstance(index, int):   
            if (index < 0):
                index = index + self._sigPtr.len()
            if (index<0) or (index >= self._sigPtr.len()):
                raise IndexError('index out of bounds')
            return self._sigPtr.value[index]
        else:
            ind = index.indices(self._sigPtr.len())
            # Easier just to count than to try to figure out all the cases
            num_items = 0
            for i in xrange(ind[0], ind[1], ind[2]):
                num_items = num_items + 1
            result = np.empty([num_items], dtype=np.float)
            i1 = 0
            for i in xrange(ind[0], ind[1], ind[2]):
                result[i1] = self._sigPtr.value[i]
                i1 = i1 + 1
        return result

    cpdef int length(self):
        """
        Length of the signal in samples.

        x.length() -> int
        """
        return self._sigPtr.len()
   
    cpdef tuple get_f0 (self, float min_f0, float max_f0, float window_dur, float frame_step):
        """Find the fundamental frequency f0 ("pitch") from the signal using the Entropic algorithm.

        get_f0 (min_f0, max_f0, window_dur, frame_step) -> np.array(dtype=float)

        Parameters:
            min_f0 : float
                Minimum allowed fundamental frequency in Hz (e.g., 100)
            max_f0 : float
                Maximum allowed fundamental frequency in Hz (e.g., 650)
            window_duration : float
                Window duration in seconds (e.g., 0.010 -- 10 milliseconds)
            frame_step : float
                Increment of window position (e.g., 0.002 -- 2 milliseconds)
        """
        cdef vec[float] floatVecPtr = self._sigPtr.get_f0(min_f0, max_f0, window_dur, frame_step)
        result = np.empty([floatVecPtr.len], dtype=np.float)
        for i in range(floatVecPtr.len):
            result[i] = floatVecPtr.data[i]
        return (result, frame_step)

    cpdef load_pcm_wav(self, string filename, bool sum_channels=True):
        """Load a single-channel pcm-encoded Microsoft wav file format.

        load_pcm_wav(filename, sum_channels=True)

        Parameters:
            filename : string 
                Path of file to load
            sum_channels : boolean (default True)
                Default True -- sum channels if multiple present.  
                Otherwise an error will be thrown for multiple channels.
        """
        self._sigPtr.load_pcm_wav(filename, sum_channels)

    cpdef load_raw_short (self, string filename, int sampling_frequency):
        """Load a single-channel pcm-encoded file of short ints with no header.

        load_pcm_wav(filename, sampling_frequency)

        Parameters:
            filename : string 
                Path to file to load
            sampling_frequency : int
                Sampling frequency of the file in Hz (e.g., 8000 for 8 kHz)
        """
        self._sigPtr.load_raw_short(filename, sampling_frequency)

    cpdef load_sph (self, string filename, int channel_num):
        """Load a NIST sphere file -- channel is 0, 1.  Use 0 for single channel.

        load_sph(filename, channel_num)

        Parameters:
            filename : string 
                Path to file to load
            channel_num : int
                Number of channel, 0 or 1, to load
        """
        self._sigPtr.load_sph(filename, channel_num)


    cpdef normalize (self):
        """Normalize the amplitude of the waveform to 16-bits.

        normalize()
        """
        self._sigPtr.normalize()

    cpdef preemphasis (self, float pre_coeff):
        """Perform pre-emphasis on the waveform; i.e., filter with 1/(1-alpha*z^(-1))

        preemphasis (alpha)

        Parameters:
            alpha : float
                Pre-emphasis coefficient
        """
        self._sigPtr.preemphasis(pre_coeff)

    cpdef remove_mean (self):
        """Remove the mean of the signal.

        remove_mean()
        """
        self._sigPtr.remove_mean()

    cpdef resample_8k (self):
        """Resample the signal to an 8 kHz sampling rate.
        Note: The resample_init() method must be called before calling this method.

        resample_8k()
        """
        self._sigPtr.resample_8k()

    cpdef resample_16k (self):
        """Resample the signal to an 16 kHz sampling rate.
        Note: The resample_init() method must be called before calling this method.  Also,
        if the sample rate is below 16 kHz, no operation will be performed.

        resample_16k()
        """
        self._sigPtr.resample_16k()

    cpdef float sampling_frequency(self):
        """Return the sampling frequency of the currently loaded signal.

        sampling_frequency()
        """
        return self._sigPtr.sampling_freq()

    cpdef save_pcm_wav(self,string filename, bool scale=False):
        """Save the current signal in pcm-encoded Microsoft wav file format.

        save_pcm_wav(filename, scale=False)

        Parameters:
            filename : string 
                Path of file to save
            scale : boolean (default False)
                Scale the output to full 16-bit range when saving
        """
        self._sigPtr.save_pcm_wav(filename, scale)
   
    cpdef save_raw_short (self,string filename, bool clip, bool scale=False):
        """Save the current signal as short ints with no-header.

        save_raw_short(filename, clip, scale=False)

        Parameters:
            filename : string 
                Path of file to save
            clip : boolean
                Clip the output if it is greater than the largest 16-bit value
            scale : boolean (default False)
                Scale the output to full 16-bit range when saving
        """
        self._sigPtr.save_raw_short(filename, clip, scale)
    
    def __setitem__(self, object index, object new_value):
        if isinstance(index, int):
            if index >= self._sigPtr.len() or index < (-1*self._sigPtr.len()):
                raise IndexError('index out of bounds')
            if index < 0:
                index = index + self._sigPtr.len()
            if not isinstance(new_value, type(self._sigPtr.value[0])):
                raise TypeError('The object you are trying to assign is not the same as the rest of the array')
            self._sigPtr.value[index] = new_value
        else:
            if not isinstance(new_value[0], type(self._sigPtr.value[0])):
                raise TypeError('The object you are trying to assign is not the same as the rest of the array')
            ind = index.indices(self._sigPtr.len())
            num_iters = 0
            for i in xrange(ind[0], ind[1], ind[2]):
                self._sigPtr.value[i] = new_value[num_iters]
                num_iters = num_iters + 1
