

"""
Detector module
===================
.. caution::
    This contains a work-in-progress implementation of a "detector sensor" module, to include noise sources into the exposure time calculations.

This is a module to implement some of the different noise sources that contributes to a reading on a detector.

This is meant to abstract some of the features such that they can be added to the final measurement through some mixin class of some sort.

This is still work in progress, so please do not import this package directly for use.
"""


import numpy as np

from castor_etc.utils import BaseClass

# TODO: implement gain profile based on measured gain from different SPI settings

class TempDependentNoise(BaseClass):
    _profile: dict

    def __init__(self, profile : dict):
        super().__init__()
        # TODO: add any checkers to make sure they're formatted properly!
        self._profile = profile

    def get_mean_value(self, target_temp : float):
        # TODO: add docstrings
        return np.interp(target_temp, self._profile['temps'], self._profile['means'])

    def get_median_value(self, target_temp : float):
        return np.interp(target_temp, self._profile['temps'], self._profile['medians'])

    def get_peak_value(self, target_temp : float):
        return np.interp(target_temp, self._profile['temps'], self._profile['peaks'])

    # TODO: add a load-profile from JSON class!
    ## The ideal behaviour is:
    ### profiles are stored as JSON files inside data, and "defaults" are assigned by the JSON it imports

# This is taken from the CIS120 test report
DARK_CURRENT_PROFILE = {
    'temps'     : np.array([20.0, 10.0, 0.0, -10.0, -20.0, -30.0, -40.0, -50.0, -60.0]),    # in Celsius
    'means'     : np.array([32.4, 9.50, 1.55, 0.340, 0.070, 0.030, 0.010, 0.004, 0.005]),   # in e-/s
    'medians'   : np.array([27.3, 7.84, 0.890, 0.130, 0.003, 0.007, 0.008, 0.002, 0.005]),  # in e-/s
    'peaks'     : np.array([30.7, 7.20, 1.20, 0.0580, 0.020, 0.009, 0.010, 0.001, 0.005])   # in e-/s
}

class DarkCurrentNoise(TempDependentNoise):
    def __init__(self, profile = DARK_CURRENT_PROFILE):
        super().__init__(profile)

READOUT_NOISE_PROFILE = {
    'temps'     : np.array([20.0, 10.0, 0.0, -10.0, -20.0, -30.0, -40.0, -50.0, -60.0]),    # in Celsius
    'means'     : np.array([6.50, 6.38, 6.52, 6.34, 5.86, 5.92, 6.00, 5.49, 5.06]),   # in e- RMS
    'medians'   : np.array([6.44, 6.32, 6.51, 6.29, 5.78, 5.88, 5.98, 5.41, 4.95]),  # in e- RMS
    'peaks'     : np.array([6.37, 6.31, 6.50, 6.25, 5.76, 5.97, 6.13, 5.36, 4.82])   # in e- RMS
}

class ReadoutNoise(TempDependentNoise):
    def __init__(self, profile = READOUT_NOISE_PROFILE):
        super().__init__(profile)
