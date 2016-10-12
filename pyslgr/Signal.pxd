#Signal.pxd
from libcpp.string cimport string
from libcpp cimport bool
from libcpp.list cimport list
from libcpp.vector cimport vector
from vec cimport vec
from GMM_model cimport GMM_model

#ctypedef vec[float] float_vec

cdef extern from "speech_tools.h":
    cdef cppclass Time_Interval:
        float start, end;
        int label;
        float score;
        Time_Interval (float start_in, float end_in, int label, float score) except +

    cdef cppclass Signal:
        float *value
        Signal() except +
        Signal (const short *data, int num_samples, int sampling_frequency) except +
        Signal (const double *data, int num_samples, int sampling_frequency) except +
        int len()
        float sampling_freq()
        list[Time_Interval] gmmsad (Signal &x, string &feat_config, vector[GMM_model] &gmm_models, int label_window, int gmmsad_keep, float min_gap_dur, float min_seg_dur) except + 
        void load_pcm_wav (string filename, bool sum_channels) except +
        void save_pcm_wav (string filename, bool scale) except +
        void enhance_ll() except +
        void load_raw_short (string filename, int sampling_frequency) except +
        void load_sph(string filename, int num_channels) except +
        void save_raw_short (string filename, bool clip, bool scale) except +
        void init_resampler (string filter_directory) except +
        void resample_8k () except +
        void resample_16k () except +
        void remove_tones() except +
        vec[float] get_f0 (float min_f0, float max_f0, float window_dur, float &frame_step) except +
        void normalize () except +
        void preemphasis (float pre_coeff) except +
        void remove_mean () except +

