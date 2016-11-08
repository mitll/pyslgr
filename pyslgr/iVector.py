#!/usr/bin/env python

import numpy as np
from GMMModel import GMMModel
import os

class iVector(object):
    """
    iVector extractor
    """
    def __init__ (self, config):
        """
        config : a dictionary of config parameters 
            'tv_matrix' : filename for total variability matrix -- raw floats
            'ubm_model' : UBM model file
        """
        # Load in the UBM model
        self._ubm = GMMModel()
        if 'model_dir' in config:
            ubm_fn = os.path.join(config['model_dir'], config['ubm_model'])
        else:
            ubm_fn = config['ubm_model']
        self._ubm.load(ubm_fn)

        # Load in the total variability matrix
        # self._U = np.load(config['tv_file'])
        # self._Ui = compute_block_matrices()

    def _compute_block_matrices (self):
        print 'compute blocks'

    def process(self, f):
        """
        f   :   LLFeatures object, input features 

        returns an ivector (factors with no scaling or transformation)
        """


        return np.array(xrange(50))
