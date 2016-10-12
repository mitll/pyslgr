#!/usr/bin/env python

import sys
import json

from pyslgr.LLSignal import LLSignal
from pyslgr.MFCCFeatures import MFCCFeatures
import os.path

import matplotlib.pyplot as plt

print '\nTesting MFCC class ...'

# JSON config file
fp = open('config/lid_config.json','r')
c = json.load(fp)
fp.close()

# Feature config for LID
feat_config = c['lid_config']['feat_config']
feat_config['fb_only'] = True
feat_config['linear'] = True
feat_config['win_len_ms'] = 20
feat_config['tgt_num_filt'] = 100
mfcc_config = json.dumps(feat_config)
print 'MFCC features configuration: {}'.format(mfcc_config)

# Load signal
fn = 'signals/example.sph'
x = LLSignal()
x.load_sph(fn, 0)
x.resample_8k()

# Get features
f = MFCCFeatures()
f.process(x, mfcc_config)
f.set_outfeat('f')
print 'Number of output features: {}'.format(f.num_outfeat())
delta_f = 0.5*x.sampling_frequency() / (feat_config['tgt_num_filt']+1)
fb_range = float(feat_config['fb_hi'])-float(feat_config['fb_low'])
print 'Approx num of feats: {}'.format(fb_range/delta_f)

b = f[:]
plt.imshow(b[-1::-1,0:500], aspect='auto')
plt.show()
