# MFCCFeatures.pxd

from MFCC_Features cimport MFCC_Features
from LLFeatures cimport LLFeatures
from libcpp.string cimport string

cdef class MFCCFeatures(LLFeatures):
    cdef MFCC_Features *_mfccPtr
    cpdef float duration (self)
    cpdef float get_win_inc_ms(self)
    cpdef load_sad_marks (self, src)
    cpdef save_sad_marks (self, string filename)
