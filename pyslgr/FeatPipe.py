#!/usr/bin/env python

from copy import deepcopy

class FeatPipe(object):
    """
    Implementation of a full fatures extraction pipeline.
    """
    def __init__ (self, config, featClass, sadClass):
       """
       config : a dictionary of config parameters with two main keys 'feat_config', 'sad_config'.
       config['feat_config'] has keys:
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
       config['sad_config'] is passed directly to the sadClass

       featClass : an LLFeatures compatible class
       sadClass  : a class with constructor sadClass(config) and method sadClass.process(LLSignal x, LLFeatures f)
       """
       self.feat_config = deepcopy(config['feat_config'])
       self.feat = featClass
       if sadClass is None:
           self.sad = None
       else:
           self.sad = sadClass(config['sad_config'])

    def process (self, x):
        """
        Extract features 
            x        Input signal

        Returns a feature object
        """

        # Basic feature processing
        f = self.feat(self.feat_config['feat_config'])
        f.process(x)

        # Rasta
        if ('do_rasta' in self.feat_config) and self.feat_config['do_rasta']:
            f.rasta()

        # Various delta style features
        if ('do_delta' in self.feat_config) and self.feat_config['do_delta']:
            if ('delta2point' in self.feat_config) and self.feat_config['delta2point']:
                f.delta2point(self.feat_config['delta_spread'])
            else:
                f.delta(self.feat_config['delta_spread'])
        if ('do_accel' in self.feat_config) and self.feat_config['do_accel']:
            f.accel(self.feat_config['accel_spread'])
            feat_str += 'a'
        if ('do_sdc' in self.feat_config) and self.feat_config['do_sdc']:
            f.sdc(*self.feat_config['sdc_params'])
        
        # Now apply SAD
        print 'Number of frames is : {}'.format(f.num_vec())
        if self.sad is not None:
            sad_marks = self.sad.process(x, f)  # returns None if applied directly to f
            if sad_marks is not None:
                print 'sad marks: {}'.format(sad_marks)
                f.load_sad_marks(sad_marks)
            f.apply_sad()
        print 'Number of frames is : {}'.format(f.num_vec())

        # Now finish off with feat norm
        if ('do_feat_norm' in self.feat_config) and self.feat_config['do_feat_norm']:
            f.feat_norm()

        # Set the output features
        if 'outfeat' in self.feat_config:
            f.set_outfeat(self.feat_config['outfeat'])

        return f

