"""
parameters.py

Contains parameters for CASTOR.

Isaac Cheng - 2022
"""

from math import pi

import astropy.units as u

# -------------------------------- TELESCOPE PARAMETERS -------------------------------- #

# The telescope filter (i.e., passband) names
PASSBANDS = ["uv", "u", "g"]

# The passband wavelength limits for the different filters
PASSBAND_LIMITS = {
    "uv": [0.150, 0.300] << u.um,
    "u": [0.300, 0.400] << u.um,
    "g": [0.400, 0.550] << u.um,
}  # microns

# Total wavelength range spanned by the passbands
PASSBAND_TOT_LIMITS = [0.150, 0.550] << u.um  # microns

# Resolution of the passband response curves
PASSBAND_RESOLUTION = 1 << u.nm  # nanometres

# Photometric zeropoints for the different filters
PHOT_ZPTS = {"uv": 24.23, "u": 24.71, "g": 24.78}  # AB mag for 1 electron/s

# The PSF full-width at half-maximum
FWHM = 0.15 << u.arcsec  # arcsec

# The angular dimension covered by each pixel
PX_SCALE = 0.1 << u.arcsec  # arcsec (or arcsec/pixel), for Table 3-6. ? UPDATE ?

# The angular area of each pixel
PX_AREA = PX_SCALE * PX_SCALE  # arcsec^2

# Instantaneous field of view
IFOV_DIMEN = [0.44, 0.56] << u.deg  # degrees. Angular dimensions
IFOV_AREA = IFOV_DIMEN[0] * IFOV_DIMEN[1]  # degrees^2

# Number of pixels of CCD
MP = 960  # megapixels

# Aperture diameter
APER_DIAMETER = 100 << u.cm  # cm

# Aperture area
APER_AREA = pi * (0.5 * APER_DIAMETER) * (0.5 * APER_DIAMETER)  # cm^2

DARK_CURRENT = 0.01  # electrons/pixel/s

BIAS = 100  # electrons

# READ_NOISE = 2.9  # electrons
READ_NOISE = 2.0  # electrons, for Table 3-6

GAIN = 2.0  # electrons/ADU

# -------------------------------------------------------------------------------------- #
