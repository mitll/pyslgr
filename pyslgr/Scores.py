# -*- coding: utf-8 -*-
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
        
