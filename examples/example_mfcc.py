#!/usr/bin/env python

import sys
import json

from pyslgr.LLSignal import LLSignal
from pyslgr.MFCCFeatures import MFCCFeatures
import os.path

print '\nTesting MFCC class ...'

# JSON config file
fp = open('config/lid_config.json','r')
c = json.load(fp)
fp.close()

# Feature config for LID -- either dictionary or string works
feat_config = c['lid_config']['feat_config']
mfcc_config = json.dumps(feat_config)
print 'MFCC features configuration: {}'.format(mfcc_config)
use_dict = True

# Load signal
fn = 'signals/example.sph'
x = LLSignal()
x.load_sph(fn, 0)

# Get features
if use_dict:
    f = MFCCFeatures(feat_config)
else:
    f = MFCCFeatures(mfcc_config)
f.process(x)
f.delta(2)
f.accel(2)
f.sdc(3, 7)
f.set_outfeat('fd')
print '\nNumber of base features: {}'.format(f.num_base_feat())
print 'Number of output features: {}'.format(f.num_outfeat())
print 'Number of frames is : {}'.format(f.num_vec())
print 'Frame rate: {}'.format(f.get_win_inc_ms())
print 'Duration in seconds is: {}'.format(f.duration())

# Indexing and slicing with features
print 'One vector, f[0]:\n{}'.format(f[0])
print 'Last vector, f[-1]:\n{}'.format(f[-1])
print 'A few vectors f[0:3]:\n{}'.format(f[0:3])
print 'Testing negative slices f[-3:-1]:\n{}'.format(f[-3:-1])
print 'Testing vector with stride f[0:5:2]:\n{}'.format(f[0:5:2])
print 'Empty vector: f[0:0]: {}'.format(f[0:0])

#
# Perform simple energy based SAD with xtalk and save labels
#
print 'Performing SAD ...'
if not os.path.exists('tmp'):
    os.mkdir('tmp')
f.xtalk(-10, 5, 5)
f.save_sad_labels(os.path.join('tmp', os.path.basename(fn) + '.lbl'))
f.save_sad_marks(os.path.join('tmp', os.path.basename(fn) + '.mark'))
lbl = f.sad_labels()
f.apply_sad()
print 'a few SAD labels: {}'.format(lbl[0:100])
print 'Number of frames is : {}'.format(f.num_vec())

# Save features
print '\nSaving features to tmp/feat.dat as raw floats.'
f.save_raw(os.path.join('tmp', os.path.basename(fn) + '.feat.dat'))

print 'Done! Successfully completed MFCC features tests.\n'
