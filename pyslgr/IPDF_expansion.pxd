# IPDF_expansion
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

cdef extern from "speech_tools.h":
    cdef cppclass IPDF_expansion:
        IPDF_expansion() except +
        int dim () except +
        vec[double] expansion (Features &f, float rf) except +
        vec[double] expansion_with_nap (Features &f, float rf, string nap_key, int corank) except +
        void load_gmm_model (string model_file_name) except +
        void load_nap_projection (string projection_file_name, string key) except +
