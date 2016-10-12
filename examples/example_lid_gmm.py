#!/usr/bin/env python

# -*- coding: utf-8 -*-
"""
Created on Mon Nov 09 17:10:16 2015

@author: JO21372
"""
from pyslgr.LLSignal import LLSignal
from pyslgr.MFCCFeatures import MFCCFeatures
from pyslgr.GMMModel import GMMModel
import json
import operator
import os
import sys

def process(fn, feat_config_lid, model_path, model_names):
    # Load signal
    x = LLSignal()
    x.load_sph(fn, 0)

    # Calculate base features
    f = MFCCFeatures()
    print "Processing features..."
    f.process(x, feat_config_lid)

    # Post process features
    f.rasta()
    f.delta2point(1)
    f.sdc(3, 7)
    f.xtalk(-10.0, 5.0)
    f.apply_sad()
    f.feat_norm()
    f.set_outfeat("fs")
    print "Done!\n"

    # print "Saving raw file... "
    # f.save_raw("tmp/" + fn + ".feat")

    # Load GMM UBM
    print "Loading GMM UBM model..."
    g_ubm = GMMModel()
    g_ubm.load(os.path.join(model_path, "language", "sbkg.o1024-sdc-rasta-fn-fnap.mod"))
    print "...done"
    
    # Load language models
    print "Loading GMM language models ..."
    g_langs = []
    for i in xrange(0,len(model_names)):
        g_langs.append(GMMModel())
        g_langs[i].load(os.path.join(model_path, model_names[i]))

    # Score models
    print "Scoring models ..."
    allScores = g_ubm.score_models(f, g_langs, 5, True, 10.0)
    print "... done"
    modelScores = allScores.scores
    scores_out = []
    for i in xrange(0, modelScores.size, 2):
        scores_out.append((model_names[i], modelScores[i+1] - modelScores[i]))
    scores_out = sorted(scores_out, key=operator.itemgetter(1), reverse=True)
    print 'scores: {}'.format(scores_out)
    print "UBM score = ", allScores.ubm_score 
 
if __name__ == "__main__":
    slgr_config_fn = 'config/lid_config.json'
    model_dir = '../models'

    # JSON config file
    fp = open(slgr_config_fn,'r')
    c = json.load(fp)
    fp.close()

    # Feature config for LID
    mfcc_config = json.dumps(c['lid_config']['feat_config'])

    # Models in order
    model_names = sorted(c['lid_config']['models'].values())

    # Processing
    fn1 = 'signals/example.sph'
    print "Processing file: " , fn1
    process(fn1, mfcc_config, model_dir, model_names)
    
    print 'Done!  Successfully completed GMM model tests.'
