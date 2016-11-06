#!/usr/bin/env python

class XtalkSAD (object):
    """
    Perform energy-based speech activity detection using Xtalk
    """

    def __init__ (self, config):
        """
        config : a dictionary of config parameters to pass to xtalk
           'abs_min_energy', 'thresh', 'med_len' (optional)
        """

        self._abs_min_energy = config['abs_min_energy']
        self._thresh = config['thresh']
        if 'med_len' in config:
            self._med_len = config['med_len']
        else:
            self._med_len = 1

    def process (self, x, f):
        f.xtalk(self._abs_min_energy, self._thresh, self._med_len)
        return None
