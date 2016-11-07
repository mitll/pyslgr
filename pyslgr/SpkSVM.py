#!/usr/bin/env python

import numpy as np
from sklearn import SVC

class SpkSVM(object):
    """
    Implementation of a speaker verification system with vector inputs and SVM training
    """
    def __init__ (self, config):
       """
       config : a dictionary of config parameters 
          'bkg'  : filename of background speakers -- numpy format array of vectors
          'c'    : regularization parameter (typically: 1.0)
       """

       # Load in the background -- numpy format
       self._bkg = 3

       self._c = config['c']

    def score(self, x):
        """
        x   :   a single vector as a numpy array for scoring
        """

    def train(self, x):
       """
       x    :   a numpy array where each row is a vector representing the target speaker

       returns a model representing the target speaker
       """

       
