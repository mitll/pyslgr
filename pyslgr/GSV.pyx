# GSV.pyx
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

from cython.operator cimport dereference as deref
from cython.operator cimport preincrement as inc
from libcpp cimport bool
from vec cimport vec
from Features cimport Features
from LLFeatures cimport LLFeatures
from GMM_model cimport GMM_model
from IPDF_expansion cimport IPDF_expansion
from GMMModel cimport GMMModel

import json
import numpy as np
import os

cdef class GSV:
    def __cinit__(self, config):
        """
        GSV(config)
        
        config  dictionary or JSON string with config parameters

        """
        self._mPtr = new IPDF_expansion()
        if self._mPtr == NULL:
            raise MemoryError()
        
        if isinstance(config, str) or isinstance(config, unicode):
            config = json.loads(config)
        elif not isinstance(config, dict):
            raise ValueError("GSV constructor: unknown input type for config")

        if 'model_dir' in config:
            model_dir = config['model_dir']
        else:
            model_dir = ""

        if 'gmm_model' not in config:
            raise ValueError("GSV constructor: gmm_model must be specified")
        self._mPtr.load_gmm_model(os.path.join(model_dir, config['gmm_model']))

        if 'rf' not in config:
            raise ValueError("GSV constructor: rf (relevance factor) must be specified")
        self._rf = config['rf']

        if 'do_nap' not in config:
            self._do_nap = False
        else:
            self._do_nap = config['do_nap']

        if 'nap_projection' in config:
            self._mPtr.load_nap_projection(os.path.join(model_dir, config['nap_projection']), 'nap')

        if self._do_nap and ('corank' not in config):
            raise ValueError("GSV constructor: corank must be specified when using NAP")
        self._corank = config['corank']

    def __dealloc__(self):
        if self._mPtr != NULL:
             del self._mPtr

    cpdef process(self, LLFeatures f):
        """
        Process input features 'f' and produce a GSV expansion

        f    LLFeatures input

        """
        cdef vec[double] exp
        if self._do_nap:
            exp = self._mPtr.expansion_with_nap(deref(f._fPtr), self._rf, "nap", self._corank)
        else:
            exp = self._mPtr.expansion(deref(f._fPtr), self._rf)

        result = np.empty([exp.len], dtype=np.double)
        for i1 in xrange(0,exp.len):
            result[i1] = exp.data[i1]
        return result
             
