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

import numpy as np
from GMMModel import GMMModel
import os
import sys

class iVector(object):
    """
    | iVector extractor
    |
    |
    | config : a dictionary of config parameters
    |
    
        +--------------+------------------------------------------------------+
        | tv_matrix    | filename for total variability matrix -- raw floats  |
        +--------------+------------------------------------------------------+
        | ubm_model    | UBM model file                                       |
        +--------------+------------------------------------------------------+
    |
    """
    def __init__ (self, config):
        # Load in the UBM model
        self._ubm = GMMModel()
        if 'model_dir' in config:
            ubm_fn = os.path.join(config['model_dir'], config['ubm_model'])
        else:
            ubm_fn = config['ubm_model']
        self._ubm.load(ubm_fn)

        # Load in the total variability matrix
        if 'model_dir' in config:
            tv_fn = os.path.join(config['model_dir'], config['tv_matrix'])
        else:
            tv_fn = config['tv_matrix']
        self._U = np.fromfile(tv_fn, dtype=np.float32)
        num_rows = self._ubm.num_fea()*self._ubm.num_mix()
        self._subdim = self._U.size/num_rows
        self._U = self._U.reshape([num_rows,self._subdim])
        self._Ui = self._compute_block_matrices()

    def _compute_block_matrices (self):
        Ui = np.empty([self._ubm.num_mix(),self._subdim,self._subdim])
        icov = self._ubm.icov()
        for im in xrange(0, self._ubm.num_mix()):
            i1 = im*self._ubm.num_fea()
            i2 = (im+1)*self._ubm.num_fea()
            UtU = ((self._U[i1:i2,:].T)*icov[i1:i2]).T  # Multiply by diagonal on left -- better way to do this?
            UtU = self._U[i1:i2,:].T.dot(UtU)  # Need to add covariance here
            Ui[im] = UtU
        return Ui

    def process(self, f):
        """
        f   :   LLFeatures object, input features 

        returns an ivector (factors with no scaling or transformation)
        """
        
        # Get 0th and 1st order stats
        ec, F = self._ubm.suff_stats(f)
        sys.stdout.flush()
        
        # Now accumulate LHS -- normal equations
        lhs = np.zeros([self._subdim,self._subdim])
        for i1 in xrange(0,self._ubm.num_mix()):
            sys.stdout.flush()
            lhs += ec[i1]*self._Ui[i1]
        sys.stdout.flush()

        # RHS, part 1, F = F - diag(ec)*m_ubm
        ubm_mean = self._ubm.mean()
        for im in xrange(0, self._ubm.num_mix()):
            i1 = im*self._ubm.num_fea()
            i2 = (im+1)*self._ubm.num_fea()
            F[i1:i2] -= ec[im]*ubm_mean[i1:i2]

        # RHS, part 2, U'*Sigma^(-1)*(F-diag(ec)*m_ubm)
        F *= self._ubm.icov()
        F = self._U.T.dot(F)

        # Now solve and return result -- most of the time will be full rank, but be careful
        # Uses SVD and Linpack GELSD
        ivec, resid, rank, s = np.linalg.lstsq(lhs, F)

        return ivec
