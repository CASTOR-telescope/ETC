"""
background_values.py

Contains average sky background values and geocoronal emission line data.

Copyright 2022, CASTOR Mission Team
Author: Isaac Cheng
Contact: isaac.cheng.ca@gmail.com
"""

import astropy.units as u

# ------------------------------- BACKGROUND NOISE VALUES ------------------------------ #

# The flux of the geocoronal emission line [O II] 2471A.
# See <https://hst-docs.stsci.edu/stisihb/chapter-6-exposure-time-calculations/6-6-tabular-sky-backgrounds
GEOCORONAL_FLUX_HIGH = 3.0e-15  # erg/cm^2/s/arcsec^2
GEOCORONAL_FLUX_AVG = 1.5e-15  # erg/cm^2/s/arcsec^2
GEOCORONAL_FLUX_LOW = 7.5e-17  # erg/cm^2/s/arcsec^2
GEOCORONAL_WAVELENGTH = 2471 << u.AA  # angstroms
GEOCORONAL_LINEWIDTH = 0.023 << u.AA  # angstroms

# Average sky background in each passband (Earthshine + zodiacal light)
SKY_BACKGROUND = {"uv": 26.08, "u": 23.74, "g": 22.60}  # AB mag/arcsec^2

# -------------------------------------------------------------------------------------- #
