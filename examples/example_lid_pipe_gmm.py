#!/usr/bin/env python

from pyslgr.LLSignal import LLSignal
from pyslgr.MFCCFeatures import MFCCFeatures
from pyslgr.GMMModel import GMMModel
from pyslgr.GMMModel import GMMSAD
from pyslgr.FeatPipe import FeatPipe
from pyslgr.sad import XtalkSAD

import json
import operator
import os
import sys

def process(fn, featPipe, model_path, model_names):
    # Load signal
    x = LLSignal()
    x.load_sph(fn, 0)

    f = featPipe.process(x)
    print "Done!\n"

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
    mfcc_pipe_fn = "config/lid_mfcc+gmmsad_pipe.json"

    # initialize feature pipe
    with open(mfcc_pipe_fn, 'r') as fp:
        pipe_config = json.load(fp)
    fpipe = FeatPipe(pipe_config, MFCCFeatures, GMMSAD)

    # slgr json file
    with open(slgr_config_fn,'r') as fp:
        slgr_config = json.load(fp)

    # Models in order
    model_names = sorted(slgr_config['lid_config']['models'].values())

    # Processing
    fn = 'signals/example.sph'
    print "Processing file: " , fn
    process(fn, fpipe, model_dir, model_names)
    
    print 'Done!  Successfully completed GMM model tests.'
