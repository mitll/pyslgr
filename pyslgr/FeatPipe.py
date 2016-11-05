#!/usr/bin/env python

from copy import deepcopy

class FeatPipe(object):
    """
    Implementation of a full fatures extraction pipeline.
    """
    def __init__ (self, config, featClass, sad_fn=None):
       """
       config : a dictionary of config parameters with keys
           'accel_spread' : int
           'delta_spread' : int
           'delta2point'  : True/False
           'do_accel'     : True/False
           'do_delta'     : True/False
           'do_rasta'     : True/False
           'do_feat_norm' : True/False
           'do_sdc'       : True/False
           'outfeat'      : string to pass to set_outfeat
           'feat_config'  : dictionary to pass directly into LLFeatures object
           'sdc_params'   : a tuple to pass to sdc -- typically (3,7)
       featClass : an LLFeatures compatible class
       sad_fn    : a function sad_fn(x, f) that takes a signal x and features and returns marks
       """
       self.config = deepcopy(config)
       self.feat_class = featClass
       self.sad_fn = sad_fn

    def process (self, x):
        """
        Extract features 
            x        Input signal
            sad_fn   

        Returns a feature object
        """

        # Basic feature processing
        f = self.feat_class()
        f.process(x, self.config['mfcc_config'])

        # Rasta
        if ('do_rasta' in self.config) and self.config['do_rasta']:
            f.rasta()

        # Various delta style features
        if ('do_delta' in self.config) and self.config['do_delta']:
            if ('delta2point' in self.config) and self.config['delta2point']:
                f.delta2point(self.config['delta_spread'])
            else:
                f.delta(self.config['delta_spread'])
        if ('do_accel' in self.config) and self.config['do_accel']:
            f.accel(self.config['accel_spread'])
            feat_str += 'a'
        if ('do_sdc' in self.config) and self.config['do_sdc']:
            f.sdc(*self.config['sdc_params'])
        
        # Now apply SAD
        if self.sad_fn is not None:
            sad_marks = self.sad_fn(x, f)
            f.load_sad_marks(sad_marks)
            f.apply_sad()

        # Now finish off with feat norm
        if ('do_feat_norm' in self.config) and self.config['do_feat_norm']:
            f.feat_norm()

        # Set the output features
        if 'outfeat' in self.config:
            f.set_outfeat(self.config['outfeat'])

        return f

