#MyClassAux.pxd
from libcpp.string cimport string

#cimport Signal
from Signal cimport Signal
from Features cimport Features

cdef extern from "speech_tools.h":
    cdef cppclass MFCC_Features(Features):
        MFCC_Features()
        void process(Signal &x, string &config) except +
        float duration() except +
        float get_win_inc_ms() except +
