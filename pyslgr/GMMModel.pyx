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
import numpy as np
from Scores import Scores
import os


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
        # cdef Features features_cpp = deref(f._fPtr)
        # print "In GMMModel.pyx score_models"
        # print "Creating GMM_model vector"
        lenModels = <int> len(models)
        # print "lenModels : " , lenModels
        cdef vector[GMM_model] models_cpp
        # print "Populating..."
        for i in range(lenModels):
            models_cpp.push_back(deref((<GMMModel>models[i])._mPtr))
        # print "Creating scores empty vec"
        cdef vec[float] *scores_cpp  = new vec[float](lenModels)
        # print "Creating frame_scores empty vec"
        cdef vec[float] *frame_scores_cpp = new vec[float](0)
        #for k in range(frame_scores.size):
        #   deref(frame_scores_cpp).data[k] = frame_scores[k]
        cdef float ubm_score = -1.0
        # print "Invoking C++ score_models..."
        self._mPtr.score_models(deref(f._fPtr), models_cpp, deref(scores_cpp), ubm_score, deref(frame_scores_cpp), topM, use_shortfall, sf_delta)
        # print "... returning from C++ score_models"
        
        #Package scores in Python objects
        #1. Scores
        # print "Packaging scores in Python object"
        scoresPy = []
        cdef int len_scores = deref(scores_cpp).len
        for j in range(len_scores):
            scoresPy.append(deref(scores_cpp).data[j])
        
        #2. Frame scores
        # print "Packaging FRAME scores in Python object"
        frameScoresPy = []
        cdef int len_frame_scores = deref(frame_scores_cpp).len
        for k in range(len_frame_scores):
            frameScoresPy.append(deref(frame_scores_cpp).data[k])
        
        allScores = Scores(frameScoresPy, scoresPy, ubm_score)
        # print "Returning from score_models"
        
        del scores_cpp
        del frame_scores_cpp
        return allScores

cdef class GMMSAD:
    """Class to process a signal and produce SAD marks using a GMM
        GMMSAD (feat_config, gmm_models, label_keep, label_window, min_gap_dur=0.5, min_seg_dur=0.2)

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
    def __cinit__(self, string feat_config, string model_dir, dict gmm_models, string label_keep='speech', label_window=50, float min_gap_dur=0.5, 
                  float min_seg_dur=0.2):
        self._config = feat_config
        self._min_gap_dur = min_gap_dur
        self._min_seg_dur = min_seg_dur
        self._label_window = label_window
        if not gmm_models.has_key(label_keep):
            raise Exception('GMMSAD: label {} not found in gmm_models dictionary'.format(label_keep))

        # Populate _gmm_models with loaded GMM models
        print "Populating..."
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
    
    cpdef gmmsad(self, LLSignal signal):
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

