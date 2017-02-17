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

class XtalkSAD (object):
    """
    | Perform energy-based speech activity detection using Xtalk
    
    
    | config : a dictionary of config parameters to pass to xtalk
    | 'abs_min_energy', 'thresh', 'med_len' (optional)
    """

    def __init__ (self, config):
        self._abs_min_energy = config['abs_min_energy']
        self._thresh = config['thresh']
        if 'med_len' in config:
            self._med_len = config['med_len']
        else:
            self._med_len = 1

    def process (self, x, f):
        f.xtalk(self._abs_min_energy, self._thresh, self._med_len)
        return None
