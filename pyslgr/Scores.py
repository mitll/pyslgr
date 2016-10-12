# -*- coding: utf-8 -*-
"""
Created on Tue Nov 10 12:45:25 2015
Class to store scores as numpy arrays
@author: JO21372
"""
import numpy as np

class Scores:
    
    """
    Expects frame scores and scores (f_scores and s respectively) as Python lists. 
    Argument u_score represents an ubm score as a float type.
    """
    def __init__(self, f_scores, s, u_score):
        
        self.frame_scores = np.array(f_scores,dtype=np.float_)
        self.scores = np.array(s, dtype=np.float_)
        self.ubm_score = u_score
        