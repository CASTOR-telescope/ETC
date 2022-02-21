"""
Contains parameters for CASTOR imaging chain.
"""

from math import pi
from os.path import join

import astropy.units as u

from .filepaths import DATAPATH

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
PASSBAND_TOT_LIMITS = [
    min(PASSBAND_LIMITS.values(), key=lambda x: x[0])[0],
    max(PASSBAND_LIMITS.values(), key=lambda x: x[1])[1],
]

# Resolution of the passband response curves
PASSBAND_RESOLUTION = 1 << u.nm  # nanometres

# Photometric zero-points for the different filters
PHOT_ZPTS = {"uv": 24.23, "u": 24.71, "g": 24.78}  # AB mag for 1 electron/s
# PHOT_ZPTS = {"uv": 24.463, "u": 24.511, "g": 24.766}  # AB mag for 1 electron/s

# Passband pivot wavelengths (CASTOR SMS values)
PASSBAND_PIVOTS = {"uv": 226 << u.nm, "u": 345 << u.nm, "g": 478 << u.nm}  # nanometres
# REVIEW: Passband pivot wavelengths, EE (my calculations)
# PASSBAND_PIVOTS = {"uv": 225 << u.nm, "u": 346 << u.nm, "g": 475 << u.nm}  # nanometres
# REVIEW: Passband pivot wavelengths, QE (my calculations)
# PASSBAND_PIVOTS = {"uv": 220 << u.nm, "u": 344 << u.nm, "g": 471 << u.nm}

# Filepaths to the passband response curves
PASSBAND_FILEPATHS = {
    band: join(DATAPATH, "passbands", f"passband_castor.{band}") for band in PASSBANDS
}

# Wavelength units in the passband response curve files
PASSBAND_FILEUNITS = {"uv": u.um, "u": u.um, "g": u.um}

# The PSF full-width at half-maximum
FWHM = 0.15 << u.arcsec  # arcsec

# The angular dimension covered by each pixel
PX_SCALE = 0.1 << u.arcsec  # arcsec (or arcsec/pixel), for Table 3-6. ? UPDATE ?

# The angular area of each pixel
PX_AREA = PX_SCALE * PX_SCALE  # arcsec^2

# Instantaneous field of view
IFOV_DIMEN = [0.44, 0.56] << u.deg  # degrees. Angular dimensions
IFOV_AREA = IFOV_DIMEN[0] * IFOV_DIMEN[1]  # degrees^2

# Number of pixels in CCD (x 1 million)
MP = 960  # megapixels

# Aperture diameter
MIRROR_DIAMETER = 100 << u.cm  # cm

# Aperture area
MIRROR_AREA = pi * (0.5 * MIRROR_DIAMETER) * (0.5 * MIRROR_DIAMETER)  # cm^2

DARK_CURRENT = 0.01  # electrons/pixel/s

BIAS = 100  # electrons

READ_NOISE = 2.0  # electrons

GAIN = 2.0  # electrons/ADU

# Wavelength threshold for red leak. Flux longward of this is considered to be red leak
REDLEAK_THRESHOLDS = {"uv": 3880 << u.AA, "u": 4730 << u.AA, "g": 5660 << u.AA}

# -------------------------------------------------------------------------------------- #
