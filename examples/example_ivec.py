#!/usr/bin/env python

import json

from pyslgr.GMMModel import GMMSAD
from pyslgr.iVector import iVector
from pyslgr.LLSignal import LLSignal
from pyslgr.MFCCFeatures import MFCCFeatures
from pyslgr.FeatPipe import FeatPipe

if __name__ == "__main__":

    # Load example signal
    fn = 'signals/example.sph'
    x = LLSignal()
    x.load_sph(fn, 0)

    # Initialize feature pipe
    mfcc_pipe_fn = "config/sid_mfcc+gmmsad_pipe.json"
    with open(mfcc_pipe_fn, 'r') as fp:
        pipe_config = json.load(fp)
    fpipe = FeatPipe(pipe_config, MFCCFeatures, GMMSAD)

    # Get MFCCFeatures
    f = fpipe.process(x)

    # Info
    print 'Number of output features: {}'.format(f.num_outfeat())
    # f.save_raw('tmp/ivec_feat.dat')

    # Load in and create GSV
    with open('config/ivec.json', 'r') as fp:
        config = json.load(fp)
    ivec = iVector(config)

    # Now compute a GSV expansion
    v = ivec.process(f)
    print 'A few elements of iVector: {}'.format(v[0:10])
    print 'iVector dimension = {}'.format(len(v))
