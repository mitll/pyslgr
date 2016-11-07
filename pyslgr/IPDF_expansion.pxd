# IPDF_expansion

from libcpp cimport bool
from libcpp.string cimport string
from vec cimport vec
from libcpp.vector cimport vector
from Features cimport Features

cdef extern from "speech_tools.h":
    cdef cppclass IPDF_expansion:
        IPDF_expansion() except +
        int dim () except +
        vec[double] expansion (Features &f, float rf) except +
        vec[double] expansion_with_nap (Features &f, float rf, string nap_key, int corank) except +
        void load_gmm_model (string model_file_name) except +
        void load_nap_projection (string projection_file_name, string key) except +
