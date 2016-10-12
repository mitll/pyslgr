# GMMModel.pxd
from libcpp cimport bool
from libcpp.string cimport string
from libcpp.vector cimport vector
from MFCC_Features cimport MFCC_Features
from GMM_model cimport GMM_model
from LLSignal cimport LLSignal
       
cdef class GMMModel:
    cdef GMM_model *_mPtr
    cpdef load(self, string model_file_name)        
    cpdef bool is_loaded(self)      

cdef class GMMSAD:
    cdef string _config
    cdef float _min_gap_dur
    cdef float _min_seg_dur
    cdef int _idx_keep
    cdef int _label_window
    cdef string _label_keep
    cdef vector[GMM_model] _gmm_models
    cpdef gmmsad(self, LLSignal signal)
