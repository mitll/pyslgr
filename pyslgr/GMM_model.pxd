# GMM_model.pxd
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

from libcpp cimport bool
from libcpp.string cimport string
from vec cimport vec
from libcpp.vector cimport vector
from Features cimport Features
ctypedef double REAL_EXP

cdef extern from "speech_tools.h":
    cdef cppclass GMM_model:
        int num_fea
        int num_mix
        GMM_model() except +
        void load (string model_file_name) except +
        bool is_loaded() except +
        void score_models (Features &f, vector[GMM_model] &models, vec[float] &scores, float &ubm_score, vec[float] &frame_scores, int topM, bool use_shortfall, float sf_delta) except +
        void suff_stats_ns (Features &f, vec[REAL_EXP] &sum, vec[REAL_EXP] &ec) except +
        vec[float] ubm_mean () except +
        vec[float] ubm_icov () except +
