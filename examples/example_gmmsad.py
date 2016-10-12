#!/usr/bin/env python

# -*- coding: utf-8 -*-
"""
Created on Mon Nov 09 17:10:16 2015

@author: JO21372
"""
from pyslgr.LLSignal import LLSignal
from pyslgr.MFCCFeatures import MFCCFeatures
from pyslgr.GMMModel import GMMSAD
import json
import operator
import os
import sys

def process_gmmsad(fn, feat_config, model_dir, model_files):
    # Load signal
    x = LLSignal()
    x.load_pcm_wav(fn)

    gs = GMMSAD(feat_config, model_dir, model_files)
    sad_marks = gs.gmmsad(x)
    print 'sad marks : {}'.format(sad_marks)

    # JSON config file
    fp = open('../config/lid_config.json','r')
    c = json.load(fp)
    fp.close()

    # Feature config for LID
    feat_config = c['lid_config']['feat_config']
    feat_config['fb_only'] = 'true'
    mfcc_config = json.dumps(feat_config)
    print 'MFCC features configuration: {}'.format(mfcc_config)

    # Get features
    f = MFCCFeatures()
    f.process(x, mfcc_config)

    # Load and apply sad marks
    print 'Number of feature vectors before SAD: {}'.format(f.num_vec())
    f.load_sad_marks(sad_marks)
    f.save_sad_marks('tmp/wow.mark')
    f.apply_sad()
    print 'Number of feature vectors after SAD: {}'.format(f.num_vec())
    # f.save_sad_marks('tmp/eg_sad.mark')

    return
 
if __name__ == "__main__":
    slgr_config_fn = '../config/slgr_params.json'
    model_dir = '../model'

    # JSON config file
    fp = open(slgr_config_fn,'r')
    c = json.load(fp)
    fp.close()

    # Feature config for LID
    feat_config = json.dumps(c['sad_config']['feat_config'])

    # Models in order
    model_files = c['sad_config']['models']
    print 'Using models: {}'.format(model_files)

    # Processing
    fn1 = '../slgr_sad/test/0000001.wav'
    # fn1 = '../slgr_test/test/dlpt-ara-D4m-RZ.wav'
    print "Processing file: " , fn1
    process_gmmsad(fn1, feat_config, model_dir, model_files)
    
    print 'Done!  Successfully completed GMMSAD tests'
