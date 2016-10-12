#
# LLFeatures.pyx
#
from Features cimport Features
from vec cimport vec
ctypedef float REAL_FEAT
import numpy as np
cimport numpy as np
from libcpp.string cimport string

cdef class LLFeatures:
    def __cinit__(self):
        if type(self) is LLFeatures:
            self._fPtr = new Features()        
        #if self._fPtr == NULL:
        #    raise MemoryError()

    def __dealloc__(self):
        if type(self) is LLFeatures:
            if self._fPtr != NULL:
                print "Destroying Features object from LLFeatures"
                del self._fPtr

    def __getitem__(self, object index):
        cdef vec[REAL_FEAT] vec_out 

        if isinstance(index, int):   
            if (index < 0):
                index = index + self._fPtr.num_vec()
            if index >= self._fPtr.num_vec():
                raise IndexError('index out of bounds')
            vec_out = self._fPtr.vector(index)
            result = np.empty([vec_out.len], dtype=np.float)
            for i in xrange(vec_out.len):
                result[i] = vec_out.data[i]
        else:
            ind = index.indices(self._fPtr.num_vec())
            # Easier just to count than to try to figure out all the cases
            num_items = 0
            for i in xrange(ind[0], ind[1], ind[2]):
                num_items = num_items + 1
            dim = self._fPtr.num_outfeat()
            result = np.empty([dim,num_items], dtype=np.float)
            i1 = 0
            for i in xrange(ind[0], ind[1], ind[2]):
                vec_out = self._fPtr.vector(i)
                for j in range(dim):
                    result[j][i1] = vec_out.data[j]
                i1 += 1
        return result

    cpdef accel (self, int accel_spread):
        """Calculate acceleration values from delta values.

        accel(accel_spread)

        Parameters:
            accel_spread : int
                Acceleration at time t is calculated using the frames t-k,...,t,...,t+k where k is 'accel_spread'
        """
        self._fPtr.accel(accel_spread)
             
    cpdef apply_sad (self):
        """Given that speech activity detection has been calculated or loaded, remove frames corresponding
        to non-speech.

        apply_sad()

        """
        self._fPtr.apply_SAD()

    cpdef delta (self, int delta_spread):
        """Calculate delta values from base features.

        delta(delta_spread)

        Parameters:
            delta_spread : int
                Delta at time t is calculated using the frames t-k,...,t,...,t+k where k is 'delta_spread'
        """
        self._fPtr.delta(delta_spread)
        
    cpdef delta2point (self, int delta_spread):
        """Calculate delta values from base features using only 2 values.

        delta2point(delta_spread)

        Parameters:
            delta_spread : int
                Delta at time t is calculated using the frames t-k and t+k where k is 'delta_spread'
        """
        self._fPtr.delta2point(delta_spread)
        
    cpdef feat_norm(self):
        """Normalize each feature individually to zero mean and unit variance across all frames.

        feat_norm()
        """
        self._fPtr.feat_norm()
        
    cpdef load_raw (self, string filename, int num_feat):
        """Load a file of raw floats into a feature store.  Assumes all features are base features.

        load_raw (filename, num_feat)

        Parameters:
            filename : string
                Name of file to load
            num_feat : int
                Number of base features (dimension of feature vector)
        """
        self._fPtr.load_raw(filename, num_feat)

    cpdef load_sad_labels (self, string filename, bool sloppy=False):
        """Load a file of SAD labels (0/1).  Assumes whitespace separation between labels.

        load_sad_labels(filename, sloppy)

        Parameters:
            filename : string
                Name of file to load
            sloppy : bool
                Complain if the number of labels doesn't match the number of feature vectors (default: True)
        """
        self._fPtr.load_SAD_label(filename, sloppy)
        
    cpdef int num_base_feat (self):
        """Return the number of base features in each vector.
        Base features do not include post-processing such as delta, sdc, or acceleration.

        num_base_feat() -> int
        """
        return self._fPtr.num_base_feat()
        
    cpdef int num_outfeat(self):
        """Return the number of output features in each vector.
        Output features are set with the set_outfeat() method.

        num_outfeat() -> int
        """
        return self._fPtr.num_outfeat()
        
    cpdef int num_total_feat(self):
        return self._fPtr.num_total_feat()
        
    cpdef int num_vec(self):
        """Return the number of feature vectors.

        num_vec() -> int
        """
        return self._fPtr.num_vec()
        
    cpdef rasta(self):
        """Apply RASTA to all features.

        rasta()
        """
        self._fPtr.rasta()

    cpdef sad_labels(self):
        """Return the current SAD labels (0/1) per frame.
        """
        cdef vec[np.uint8_t] vec_out 

        if not self._fPtr.SAD_available():
            raise ValueError("SAD unavailable.")
        num_vec = self._fPtr.num_vec()
        vec_out = self._fPtr.SAD_labels()
        result = np.empty([vec_out.len], dtype=np.uint8)
        for i in xrange(vec_out.len):
            result[i] = vec_out.data[i]
        return result
        
    cpdef save_raw (self, string filename):
        """Save features as raw floats to 'filename'.  The features
        saved are set by 'set_outfeat'.

        save_raw(filename)

        Parameters:
            filename : string
                Name of file to save features in
        """
        self._fPtr.save_raw(filename)
        
    cpdef save_sad_labels (self, string filename):
        """Save SAD labels (0/1) to file. 

        save_sad_labels(filename)

        Parameters:
            filename : string
                Name of output file
        """
        self._fPtr.save_SAD_label(filename)

    cpdef sdc (self, int sdc_p, int sdc_k_in): 
        """Shifted-delta features -- typically used for language recognition.
        Note: Uses available delta features which must be calculated before invoking 'sdc()'.

        sdc (sdc_p, sdc_k): 

        Parameters:
            sdc_p : int 
                p value -- shift between delta blocks (typical value, 3 )
            sdc_k : int
                k value -- number of delta blocks to stack (typical value, 7) 
        """
        self._fPtr.sdc(sdc_p, sdc_k_in)

    cpdef set_outfeat (self, string outfeat):
        """Set the features to output for typical operations.
        The value can be changed at any time.  Order in the parameter string
        determines the stacking order.  By default the 'outfeat' is set to 'all' when the
        feature object is created.

        set_outfeat(outfeat)

        Parameters:
            outfeat : string
                Set the output features.  If outfeat=='all', then all base features and
                calculated features are returned (except energy.  Otherwise, 'outfeat' 
                is examined a character at a time and features are stacked in that order.
                'f' base features, 'd' delta-features, 'a' acceleration features, 'e' energy, 
                's' sdc features.  E.g., set_outfeat('fd') would set the output to base features in indices 
                0, ..., num_base_feat-1 and delta features in num_base_feat+1, ..., -1.
        """
        self._fPtr.set_outfeat(outfeat)
        
    cpdef xtalk (self, double abs_min_energy, double thresh, int med_len=1):
        """xtalk energy based speech-activity detection.

        xtalk (abs_min_energy, thresh, med_len=1)

        Parameters:
            abs_min_energy : float
                Below this threshold is non-speech.  Typical values, -10 or 0.
            thresh : float
                Above this threshold triggers speech activity (the algorithm is adaptive).
            med_len: int (default 1)
                Median filter to smooth activity.  Large values imply less abrupt changes in
                speech activity.
        """
        self._fPtr.xtalk(abs_min_energy, thresh, med_len)

