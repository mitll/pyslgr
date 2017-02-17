#!/usr/bin/env python
#
# Copyright 2016 MIT Lincoln Laboratory, Massachusetts Institute of Technology
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not use these files except in compliance with
# the License.
#
# You may obtain a copy of the License at
#
#     http:#www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on
# an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the
# specific language governing permissions and limitations under the License.
#

from copy import deepcopy

class FeatPipe(object):
    """
    | Implementation of a full fatures extraction pipeline.
    |
    |
    | config : a dictionary of config parameters with two main keys 'pipe_config', 'sad_config'.
    | config['pipe_config'] has keys:
    |
    
        +------------------+------------------------------------------------------------------+
        | accel_spread     |  int                                                             |
        +------------------+------------------------------------------------------------------+
        | delta_spread     |  int                                                             |
        +------------------+------------------------------------------------------------------+
        | delta2point      |  True/False                                                      |
        +------------------+------------------------------------------------------------------+
        | do_accel         |  True/False                                                      |
        +------------------+------------------------------------------------------------------+
        | do_delta         |  True/False                                                      |
        +------------------+------------------------------------------------------------------+
        | do_rasta         |  True/False                                                      |
        +------------------+------------------------------------------------------------------+
        | do_feat_norm     |  True/False                                                      |
        +------------------+------------------------------------------------------------------+
        | do_sdc           |  True/False                                                      |
        +------------------+------------------------------------------------------------------+
        | outfeat          |  string to pass to set_outfeat                                   |
        +------------------+------------------------------------------------------------------+
        | feat_config      |  dictionary to pass directly into LLFeatures object              |
        +------------------+------------------------------------------------------------------+
        | sdc_params       |  a tuple to pass to sdc -- typically (3,7)                       |
        +------------------+------------------------------------------------------------------+
    |
    | config['sad_config'] is passed directly to the sadClass
    |
    | featClass : an LLFeatures compatible class
    | sadClass  : a class with constructor sadClass(config) and method sadClass.process(LLSignal x, LLFeatures f)
    """
    def __init__ (self, config, featClass, sadClass):
       self.pipe_config = deepcopy(config['pipe_config'])
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
        f = self.feat(self.pipe_config['feat_config'])
        f.process(x)

        # Rasta
        if ('do_rasta' in self.pipe_config) and self.pipe_config['do_rasta']:
            f.rasta()

        # Various delta style features
        if ('do_delta' in self.pipe_config) and self.pipe_config['do_delta']:
            if ('delta2point' in self.pipe_config) and self.pipe_config['delta2point']:
                f.delta2point(self.pipe_config['delta_spread'])
            else:
                f.delta(self.pipe_config['delta_spread'])
        if ('do_accel' in self.pipe_config) and self.pipe_config['do_accel']:
            f.accel(self.pipe_config['accel_spread'])
            feat_str += 'a'
        if ('do_sdc' in self.pipe_config) and self.pipe_config['do_sdc']:
            f.sdc(*self.pipe_config['sdc_params'])
        
        # Now apply SAD
        if self.sad is not None:
            sad_marks = self.sad.process(x, f)  # returns None if applied directly to f
            if sad_marks is not None:
                f.load_sad_marks(sad_marks)
            f.apply_sad()

        # Now finish off with feat norm
        if ('do_feat_norm' in self.pipe_config) and self.pipe_config['do_feat_norm']:
            f.feat_norm()

        # Set the output features
        if 'outfeat' in self.pipe_config:
            f.set_outfeat(self.pipe_config['outfeat'])

        return f

