#!/usr/bin/env python

import json

from pyslgr.GMMModel import GMMSAD
from pyslgr.GSV import GSV
from pyslgr.LLSignal import LLSignal
from pyslgr.MFCCFeatures import MFCCFeatures
from pyslgr.FeatPipe import FeatPipe
from pyslgr.sad import XtalkSAD

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
    # f.save_raw('tmp/gsv_feat.dat')

    # Load in and create GSV
    with open('config/gsv.json', 'r') as fp:
        config = json.load(fp)
    gsv = GSV(config)

    # Now compute a GSV expansion
    v = gsv.process(f)
    print 'A few elements of GSV expansion: {}'.format(v[0:10])
    print 'GSV expansion dimension = {}'.format(len(v))
	
    print 'Done!  Successfully completed GSV tests'