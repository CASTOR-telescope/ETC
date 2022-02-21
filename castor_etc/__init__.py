"""
CASTOR Exposure Time Calculator (ETC)
=====================================

A Python package for easy analysis of CASTOR performance and modifications to CASTOR
parameters in Python scripts. See the `ETC_frontend` GitHub repository for a graphical
user interface to complement this package.

Includes:
  1. Astronomical source generation and background noise estimation
  2. Telescope imaging chain simulation (both photometry and spectroscopy)
  3. Convenience functions for converting between useful quantities (e.g., flux to
     electron/s to AB magnitude)

Copyright 2022, CASTOR Mission Team
Author: Isaac Cheng
Contact: isaac.cheng.ca@gmail.com
"""

__all__ = [
    "constants",
    "conversions",
    "data",
    "filepaths",
    "load_files",
    "parameters",
    "snr",
    "sources",
    "spectrum",
    "telescope",
]

from . import constants
from . import conversions
from . import data
from . import filepaths
from . import load_files
from . import parameters
from . import snr
from . import sources
from . import spectrum
from . import telescope
