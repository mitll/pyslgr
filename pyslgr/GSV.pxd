# GSV.pxd

from libcpp cimport bool
from libcpp.string cimport string
from libcpp.vector cimport vector
from IPDF_expansion cimport IPDF_expansion
from MFCC_Features cimport MFCC_Features
from GMM_model cimport GMM_model
from LLSignal cimport LLSignal
from LLFeatures cimport LLFeatures
       
cdef class GSV:
    cdef IPDF_expansion *_mPtr
    cdef int _corank
    cdef bool _do_nap
    cdef float _rf
    cpdef process (self, LLFeatures f)
