# GMMModels.pyx
from cython.operator cimport dereference as deref
from cython.operator cimport preincrement as inc
from libcpp.list cimport list as list_cpp
from libcpp.vector cimport vector
from libcpp cimport bool
from MFCCFeatures cimport MFCCFeatures
from vec cimport vec
from Features cimport Features
from LLFeatures cimport LLFeatures
from GMM_model cimport GMM_model
from GMMModel cimport GMMModel
from Signal cimport Time_Interval
import json
import numpy as np
from Scores import Scores
import os
ctypedef double REAL_EXP

cdef class GMMModel:
    def __cinit__(self):
        self._mPtr = new GMM_model()
        # print 'Creating new model!'
        if self._mPtr == NULL:
            raise MemoryError()

    def __dealloc__(self):
        if self._mPtr != NULL:
             # print "Destroying GMM_model from GMMModel.pyx"
             del self._mPtr   
             
    cpdef load(self, string model_file_name):
        """ Load a GMM model.
    
        load(model_file_name)

        Parameters:
            model_file_name : string
                File name of model to load
        """
        self._mPtr.load(model_file_name)
    
    cpdef bool is_loaded(self):
        return self._mPtr.is_loaded()

    def icov (self):
        """ Return the GMM diagonal of the inverse covariance matrix as a numpy array
        """
        cdef vec[float] ubm_cpp 
        ubm_cpp = self._mPtr.ubm_icov()
        ubm_icov = np.empty([ubm_cpp.len], dtype=np.float)
        for i1 in xrange(0, ubm_cpp.len):
            ubm_icov[i1] = ubm_cpp.data[i1]
        return ubm_icov

    def mean (self):
        """ Return the GMM mean vector as a numpy array
        """
        cdef vec[float] ubm_cpp 
        ubm_cpp = self._mPtr.ubm_mean()
        ubm_mean = np.empty([ubm_cpp.len], dtype=np.float)
        for i1 in xrange(0, ubm_cpp.len):
            ubm_mean[i1] = ubm_cpp.data[i1]
        return ubm_mean

    def num_fea(self):
        return self._mPtr.num_fea

    def num_mix(self):
        return self._mPtr.num_mix

    def score_models (self, LLFeatures f, list models, int topM=5, bool use_shortfall=True, float sf_delta=10.0):
        """ Load a GMM model.
    
        load(LLFeatures f, list models, int topM, bool use_shortfall, float sf_delta)

        Parameters:
            f  : LLFeatures
               Feature input
            models
               List of GMM models to process
            topM : int
               Number of Gaussians to score per frame (default: 5)
            use_shortfall : bool
               Use shortfall method when evaluating Gaussian mixture model (default: True)
            sf_delta : float
               Shortfall delta for pruning (default: 10.0)
        """
        lenModels = <int> len(models)
        cdef vector[GMM_model] models_cpp
        for i in range(lenModels):
            models_cpp.push_back(deref((<GMMModel>models[i])._mPtr))
        cdef vec[float] *scores_cpp  = new vec[float](lenModels)
        cdef vec[float] *frame_scores_cpp = new vec[float](0)
        cdef float ubm_score = -1.0
        self._mPtr.score_models(deref(f._fPtr), models_cpp, deref(scores_cpp), ubm_score, deref(frame_scores_cpp), topM, use_shortfall, sf_delta)
        
        # Package scores in Python objects
        # 1. Scores
        # print "Packaging scores in Python object"
        scoresPy = []
        cdef int len_scores = deref(scores_cpp).len
        for j in range(len_scores):
            scoresPy.append(deref(scores_cpp).data[j])
        
        # 2. Frame scores
        # print "Packaging FRAME scores in Python object"
        frameScoresPy = []
        cdef int len_frame_scores = deref(frame_scores_cpp).len
        for k in range(len_frame_scores):
            frameScoresPy.append(deref(frame_scores_cpp).data[k])
        
        allScores = Scores(frameScoresPy, scoresPy, ubm_score)
        
        del scores_cpp
        del frame_scores_cpp
        return allScores

    def suff_stats (self, LLFeatures f):
        # void suff_stats_ns (Features &f, vec[REAL_EXP] &sum, vec[REAL_EXP] &ec) except +

        cdef vec[REAL_EXP] *ec = new vec[REAL_EXP](self._mPtr.num_mix)   # expected counts
        cdef vec[REAL_EXP] *sum = new vec[REAL_EXP](self._mPtr.num_mix*f.num_outfeat())  # sum
        self._mPtr.suff_stats_ns(deref(f._fPtr), deref(sum), deref(ec))

        # Create and copy expected counts into numpy array
        ec_np = np.empty([ec.len], dtype=np.double)
        for i1 in xrange(0, ec.len):
            ec_np[i1] = ec.data[i1]

        # Create and copy sum into numpy array
        sum_np = np.empty([sum.len], dtype=np.double)
        for i1 in xrange(0, sum.len):
            sum_np[i1] = sum.data[i1]
        del ec
        del sum

        return ec_np, sum_np

cdef class GMMSAD:
    """Class to process a signal and produce SAD marks using a GMM

        GMMSAD (feat_config, gmm_models, label_keep, label_window, min_gap_dur=0.5, min_seg_dur=0.2)
        or GMMSAD (config) where config is a JSON string or dictionary. 

        Parameters:
            feat_config : string
                Configuration string from a JSON object -- see MFCCFeatures.process for more details
            model_dir : string
                Base directory where models are stored
            gmm_models : dict
                Dictionary of models for GMMSAD scoring.  Typical keys are 'speech', 'music', 'nonspeech'.
            label_keep : string
                The key of the model to keep (default: 'speech')
            label_window : int
                Window length for smoothing frame scores
            min_gap_dur : float (default 0.5 seconds)
                Minimum gap between segments -- segments are combined if the gap is smaller than this time
            min_seg_dur : float (default 0.2 seconds)
                Minimum segment duration -- ignore segments below this duration
    """

    def __cinit__(self, config, string model_dir="", dict gmm_models=None, string label_keep='speech', int label_window=50, float min_gap_dur=0.5, float min_seg_dur=0.2):
        if isinstance(config, str) or isinstance(config, unicode):
            config = json.loads(config)
        elif not isinstance(config, dict):
            raise ValueError("GMMSAD constructor: unknown input type")

        # Cython does not allow overloaded constructors, so we must check to make sure 
        if ('model_dir' not in config) and (model_dir==(<string> "")):
            raise ValueError("GMMSAD constructor: model_dir must be specified")
        if ('gmm_models' not in config) and (gmm_models is None):
            raise ValueError("GMMSAD constructor: model_dir must be specified")

        # Now fill in values
        if 'feat_config' in config:
            feat_config = config['feat_config']
        else:
            feat_config = config
        if model_dir==(<string> ""):
            model_dir = config['model_dir']
        if gmm_models is None:
            gmm_models = config['gmm_models']
        if 'label_keep' in config:
            label_keep = config['label_keep']
        if 'label_window' in config:
            label_window = config['label_window']
        if 'min_gap_dur' in config:
            min_gap_dur = config['min_gap_dur']
        if 'min_seg_dur' in config:
            min_seg_dur = config['min_seg_dur']

        if isinstance(feat_config, dict):
            feat_config_str = json.dumps(feat_config)
        elif isinstance(feat_config, str) or isinstance(feat_config, unicode):
            feat_config_str = feat_config
        else:
            raise ValueError("GMMSAD constructor: unknown input type")
    
        self._config = feat_config_str
        self._min_gap_dur = min_gap_dur
        self._min_seg_dur = min_seg_dur
        self._label_window = label_window
        if not gmm_models.has_key(label_keep):
            raise Exception('GMMSAD: label {} not found in gmm_models dictionary'.format(label_keep))

        # Populate _gmm_models with loaded GMM models
        i1 = 0
        cdef string ky
        cdef GMM_model *curr_model
        for ky in gmm_models.keys():
            if ky==label_keep:
                self._idx_keep = i1
            model_fn = gmm_models[ky]
            curr_model = new GMM_model()
            curr_model.load(os.path.join(model_dir, model_fn))
            self._gmm_models.push_back(deref(curr_model))
            i1 += 1
    
    cpdef process(self, LLSignal signal, LLFeatures f):
        """
           Process LLSignal and return a list of tuples.  Each tuple is a start time (seconds) and duration of the 
           detected class.

           gmmsad (signal)
              signal   LLSignal 
              f        LLFeatures object (not used, but available for uniformity)
        """
        out = signal._sigPtr.gmmsad(deref(signal._sigPtr), self._config, self._gmm_models, self._label_window, self._idx_keep, self._min_gap_dur, self._min_seg_dur)
        out_python = []
        cdef list_cpp[Time_Interval].iterator it = out.begin()
        out_python = []
        cdef float start_time, dur_time
        while it != out.end():
            start_time = deref(it).start
            dur_time = deref(it).end-start_time
            out_python.append((start_time, dur_time))
            inc(it)
        return out_python

