# LLSignal.pxd
from Signal cimport Signal
from libcpp.string cimport string
from libcpp cimport bool

cdef class LLSignal:
    cdef Signal *_sigPtr
    cpdef int length(self)        
    cpdef float sampling_frequency(self)        

    cpdef load_pcm_wav(self, string filename, bool sum_channels=*)
    cpdef load_raw_short (self, string filename, int sampling_frequency)
    cpdef load_sph(self, string filename, int num_channels)
    cpdef save_pcm_wav(self,string filename, bool scale=*)    
    cpdef save_raw_short (self,string filename, bool clip, bool scale=*)

    cpdef resample_8k (self)
    cpdef resample_16k (self)

    cpdef tuple get_f0 (self, float min_f0, float max_f0, float window_dur, float frame_step)

    cpdef normalize (self)
    cpdef preemphasis (self, float pre_coeff)
    cpdef remove_mean (self)
