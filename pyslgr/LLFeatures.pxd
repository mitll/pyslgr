# FeaturesPy.pxd
from libcpp cimport bool
from Features cimport Features
from libcpp.string cimport string
       
cdef class LLFeatures:
    cdef Features *_fPtr
    cpdef load_raw (self, string filename, int num_feat) 
    cpdef int num_base_feat (self) 
    cpdef int num_outfeat (self)
    cpdef int num_total_feat(self) 
    cpdef int num_vec(self) 
    cpdef save_raw (self, string filename) 
    cpdef set_outfeat (self, string outfeat) 
    cpdef accel (self, int accel_spread) 
    cpdef delta (self, int delta_spread) 
    cpdef delta2point (self, int delta_spread) 
    cpdef feat_norm(self) 
    cpdef rasta(self)
    cpdef sdc (self, int sdc_p, int sdc_k_in) 

    cpdef apply_sad (self)
    cpdef sad_labels (self)
    cpdef load_sad_labels (self, string filename, bool sloppy=*)
    cpdef save_sad_labels (self, string filename)

    cpdef xtalk (self, double abs_min_energy, double thresh, int med_len=*)

