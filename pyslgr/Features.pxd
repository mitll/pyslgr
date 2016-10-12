# Features.pxd
from libcpp cimport bool
from libcpp.string cimport string
from libcpp.list cimport list as list_cpp
from vec cimport vec
cimport numpy as np

cdef extern from "speech_tools.h":
    ctypedef float REAL_FEAT
    cdef cppclass Time_Interval:
        float start, end;
        int label;
        float score;
        Time_Interval (float start_in, float end_in, int label, float score) except +

    cdef cppclass Features:
        Features() except +
        # Features& operator= (const Features& src) except +

        # General purpose manipulation/access
        void load_raw (string filename, int num_feat) except + #load from raw floats
        int num_base_feat () except +
        int num_outfeat () except +
        int num_total_feat() except +
        int num_vec() except +
        void save_raw (string filename) except + #save to raw floats
        void set_outfeat (string outfeat) except + #options are : "all" or arbitrary combinations of "feads"
        vec[REAL_FEAT] vector(int i) # Get the ith feature vector

        # General purpose transformation
        void accel (int accel_spread) except +
        void delta (int delta_spread) except +
        void delta2point (int delta_spread) except +
        void feat_norm() except +  # Window is the entire feature set
        void rasta() except +
        void sdc (int sdc_p, int sdc_k_in) except +

        # SAD related, sloppy = ignore truncation errors in the marks
        void apply_SAD () except +
        bool SAD_available() except +
        void load_SAD_label(string filename, bool sloppy) except +
        void load_SAD_marks (string filename, float in_win_inc_ms) except +
        void load_SAD_marks (list_cpp[Time_Interval] &tlist, float in_win_inc_ms) except +
        void save_SAD_label(string filename) except +
        void save_SAD_marks (string filename, float in_win_inc_ms) except +
        vec[np.uint8_t] SAD_labels() except +	     
        void xtalk (double abs_min_energy, double thresh, int med_len) except +
