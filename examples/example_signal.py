#!/usr/bin/env python
#
# Examples of various signal operations using pyslgr
#

import copy
import numpy as np
import os
import sys

from pyslgr.LLSignal import LLSignal

def get_f0 (infile):
    x = LLSignal()
    print 'Loading file: {}'.format(infile)
    x.load_sph(infile, 0)
    f0 = x.get_f0(100, 650, 0.010, 0.005)
    return f0

fn = 'signals/example.sph'
print 'Processing signal examples for file: {}\n'.format(fn)

print 'Loading a signal and finding the pitch ...'
f0 = get_f0(fn)
print 'slice of f0 array [130:135]: {}'.format(f0[0][130:135])

print "\nAttempting to load a signal that doesn't exist and generate an error ..."
fn2 = '../slgr_test/test/dummy.wav'
try:
    f0 = get_f0(fn2)
except RuntimeError as e:
    print 'Caught exception: {}'.format(e)

print '\nLoad signal values from an array ...'
x = np.array(range(100), dtype=np.short)
y = LLSignal(x, 8000)

print '\nTesting indexing ...'
print 'Full signal: {}'.format(y[0:])
print 'value for -1 and length()-1: {} {}'.format(y[-1], y[y.length()-1])
print 'values for slice[10:20]: {}'.format(y[10:20])
print 'empty slice [3:3]: {}'.format(y[3:3])
print 'empty slice [3:1]: {}'.format(y[3:1])
print 'negative slice [-3:]: {}'.format(y[-3:])
print 'clipping of indices [20:130]: {}'.format(y[20:130])
print 'Testing stride [20::3]: {}'.format(y[20::3])
print 'Testing reversal [10::-1]: {}'.format(y[10::-1])

print '\nLoading signal again ...'
x1 = LLSignal()
x1.load_sph(fn, 0)
if not os.path.exists('tmp/'):
    os.mkdir('tmp')

print '\nResampling signal ...'
print 'Initial sampling frequency: {}'.format(x1.sampling_frequency())
x1.resample_8k()
print 'Final sampling frequency: {}'.format(x1.sampling_frequency())
x1.save_pcm_wav('tmp/signal_enh_resampled.wav')
x1.save_pcm_wav('tmp/signal_enh_resampled_scaled.wav', True)

print '\nBasic signal operations:'
print 'mean before: {}'.format(np.mean(x1[0:]))
x1.remove_mean()
print 'mean after: {}'.format(np.mean(x1[0:]))
print 'max value before: {}'.format(np.max(x1[0:]))
x1.normalize()
print 'max value after: {}'.format(np.max(x1[0:]))
print 'Preemphasize signal:'
x1.preemphasis(0.97)
x1.save_pcm_wav('tmp/signal_preemphasized.wav')

print '\nCopying and changing signal values'
x2 = copy.copy(x1)
print 'x1[55]: {}'.format(x1[55])
print 'x2[55]: {}'.format(x2[55])
x2[55] = 22.0
print 'x1[55]: {}'.format(x1[55])
print 'x2[55]: {}\n'.format(x2[55])
print 'x1[55:58]: {}'.format(x1[55:58])
print 'x2[55:58]: {}'.format(x2[55:58])
x2[55:58] = np.array([2.0,3.0,4.0])
print 'x1[55:58]: {}'.format(x1[55:58])
print 'x2[55:58]: {}'.format(x2[55:58])

print 'Signal tests completed successfully!!!'
