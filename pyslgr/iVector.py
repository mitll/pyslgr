#!/usr/bin/env python

import numpy as np
from GMMModel import GMMModel

class iVector(object):
    """
    Implementation of an iVector extractor
    """
    def __init__ (self, config):
       """
       config : a dictionary of config parameters 
          'tv_matrix' : filename for total variability matrix -- should be a numpy array
          'ubm_model' : UBM model file
       """

       # Load in the total variability matrix
       self._U = np.load(config['tv_file'])
       self._Ui = compute_block_matrices()

    def _compute_block_matrices (self):
        print 'wow'

    def score(self, x):
        """
        x   :   a single vector as a numpy array for scoring
        """

    def train(self, x):
       """
       x    :   a numpy array where each row is a vector representing the target speaker

       returns a model representing the target speaker
       """

       
