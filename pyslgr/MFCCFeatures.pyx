#
# MFCCFeatures -- python interface
#
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
from libcpp.string cimport string
from libcpp.list cimport list as list_cpp
from MFCC_Features cimport MFCC_Features
from LLSignal cimport LLSignal
from LLFeatures cimport LLFeatures
from Features cimport Time_Interval
from cython.operator cimport dereference as deref
import json

cdef class MFCCFeatures(LLFeatures):

    def __cinit__(self, config):
        """
        config  String or dictionary containing configuration parameters for MFCCs.
                Parameters in the config are:
                alpha         Warping factor for bilinear method (no warping: 1.0)
                dither        0/1 - Add low level noise to the signal (typical: 1)
                fb_low        Lowest filter bank frequency in Hz (typical: 300)
                fb_hi         Highest filter bank frequency in Hz (typical: 3140)
                fb_only       0/1 - Instead of producing cepstral coefficients produce the 'raw' filter bank outputs instead
                keep_c0       0/1 - Keep the c0 cepstral coefficient; c0 represents frame energy (typical: 0)
                linear        true/false - linear or mel-warped scale for filter banks (typical: false)
                num_cep       int - number of cepstral coefficients (c1-c??) to output (typical: 7-19)
                tgt_num_filt  int - number of filters across the entire bandwidth; only applied for linear=true
                win_inc_ms    int - window increment in milliseconds (typical: 10)
                win_len_ms    int - window length in milliseconds (typical: 20-30)
        """
        if type(self) is MFCCFeatures:            
            self._mfccPtr = self._fPtr = new MFCC_Features()
        if self._mfccPtr == NULL:
            raise MemoryError()
        if isinstance(config, dict):
            self._config_str = json.dumps(config)
        elif isinstance(config, str) or isinstance(config, unicode):
            self._config_str = config
        else:
            raise ValueError("MFCCFeatures constructor: unknown input type for config")

    def __dealloc__(self):
        if type(self) is MFCCFeatures:
            if self._mfccPtr != NULL:
                # print "Destroying MFCC_Features object from MFCCFeatures.pyx"
                del self._mfccPtr

    cpdef float duration (self):
        """Return the duration that the MFCC data spans in seconds.

        duration() -> float seconds
        """
        return self._mfccPtr.duration()
        
    def process(self, LLSignal signal):
        """Process the signal to return mel-frequency cepstral coefficient (MFCC) features.

        process(signal) -> features

        Parameters:
                signal  Input signal -- instance of LLSignal class
        """

        self._mfccPtr.process(deref(signal._sigPtr), self._config_str)
       
    cpdef float get_win_inc_ms(self):
        """Return the window increment in milliseconds.

        get_win_inc_ms()
        """
        return self._mfccPtr.get_win_inc_ms()

    cpdef load_sad_marks (self, src):
        """Load SAD marks from a file or list

        load_sad_marks(src)

        Parameters:
            src : string
                Name of input file
            or 
            src : list of tuples 
                Tuples with start, duration in seconds: [(0.0,1.0),(2.0,1.5)] 
        """
        win_inc_float = self._mfccPtr.get_win_inc_ms()
        cdef list_cpp[Time_Interval] tlist
        
        if isinstance(src, str) or isinstance(src, unicode):
            self._fPtr.load_SAD_marks(<string> src, <float> win_inc_float)
        elif isinstance(src, list):
            for pr in src:
                ti = new Time_Interval(pr[0], pr[0]+pr[1], 1, 1.0)
                tlist.push_back(deref(ti))
            self._fPtr.load_SAD_marks(<list_cpp[Time_Interval] &> tlist, <float> win_inc_float)
        else:
            raise ValueError("load_sad_marks: unknown input type")

    cpdef save_sad_marks (self, string filename):
        """Save SAD marks to a file. 

        save_sad_marks(filename)

        Parameters:
            filename : string
                Name of output file
        """
        win_inc_float = self._mfccPtr.get_win_inc_ms()
        self._fPtr.save_SAD_marks(filename, int(win_inc_float))
        
    @staticmethod
    def get_sid_config():
        """Returns default speaker id configuration as a dictionary. User can 
        modify the entries in the dictionary or use this configuration as is. User 
        must call static method config_dict_to_str(config) before processing a signal with
        such configuration.
        """
        sid_config = dict()
        sid_config["feat_type"] = "mfcc"
        sid_config["alpha"] = "1"
        sid_config["dither"]="1"
        sid_config["fb_low"] = "300"
        sid_config["fb_hi"] = "3140"
        sid_config["keep_c0"] = "false"
        sid_config["linear"] = "false"
        sid_config["num_cep"] = "19"
        sid_config["win_inc_ms"] = "10"
        sid_config["win_len_ms"] = "20"
        return sid_config
    
    @staticmethod        
    def get_lid_config():
        """Returns default language id configuration as a dictionary. User can 
        modify the entries in the dictionary or use this configuration as is. User 
        must call static method config_dict_to_str(config) before processing a signal with
        such configuration.
        """
        lid_config = dict()
        lid_config["feat_type"] = "mfcc"
        lid_config["alpha"] = "1"
        lid_config["dither"] = "1"
        lid_config["fb_low"] = "300"
        lid_config["fb_hi"] = "3140"
        lid_config["keep_c0"] = "true"
        lid_config["linear"] = "false"
        lid_config["num_cep"] = "6"
        lid_config["win_inc_ms"] = "10"
        lid_config["win_len_ms"] = "20"
        return lid_config
    
    @staticmethod
    def config_dict_to_str(config_dict):
        """Converts a feature configuration dictionary to a json string for processing.
        """
        conf = str(config_dict)
        return conf.replace("'",'"')
