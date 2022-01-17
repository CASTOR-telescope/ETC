"""
constants.py

A file to track filepaths and telescope parameters.

Isaac Cheng - January 2022

Common units (also see <https://pysynphot.readthedocs.io/en/latest/units.html>):
  - fnu :: erg/s/cm^2/Hz
  - flam :: erg/s/cm^2/A
  - photnu :: photons/s/cm^2/Hz
  - photlam :: photons/s/cm^2/A
  - enu :: electrons/s/cm^2/Hz
  - elam :: electrons/s/cm^2/A
"""

from os.path import dirname, abspath
from math import pi

# ------------------------- FILEPATHS (N.B. the trailing slash) ------------------------ #

# The location of the ETC repo
BASEPATH = dirname(dirname(abspath(__file__))) + "/"

# The directory to store output files. This directory must already exist and be
# discoverable by the container (e.g., via the "-v" bind mount for a local Docker build)
OUTPATH = "/arc/home/IsaacCheng/CASTOR/ETC_plots/"  # ! CHANGE ME !

# The directory containing the data files (e.g., passbands, sky background, etc.)
DATAPATH = BASEPATH + "data/"

# -------------------------------------------------------------------------------------- #

# -------------------------------- TELESCOPE PARAMETERS -------------------------------- #

# The telescope filter (i.e., passband) names
PASSBANDS = ["uv", "u", "g"]

# The passband wavelength limits for the different filters
PASSBAND_LIMITS = {
    "uv": [0.150, 0.300],
    "u": [0.300, 0.400],
    "g": [0.400, 0.550],
}  # microns

# The PSF full-width at half-maximum
FWHM = 0.15  # arcsec

# The angular dimension covered by each pixel
PX_SCALE = 0.1  # arcsec (or arcsec/pixel)

# The angular area of each pixel
PX_AREA = PX_SCALE * PX_SCALE  # arcsec^2

# Instantaneous field of view
IFOV = 0.25  # square degrees

# Number of pixels of CCD
MP = 960  # megapixels

# Aperture diameter
APER_DIAMETER = 100  # cm

# Aperture area
APER_AREA = pi * (0.5 * APER_DIAMETER) * (0.5 * APER_DIAMETER)  # cm^2

DARK_CURRENT = 0.01  # electrons/pixel/s

BIAS = 100  # electrons

READ_NOISE = 2.9  # electrons

GAIN = 2.0  # electrons/ADU

# -------------------------------------------------------------------------------------- #

# ----------------------------------- OTHER CONSTANTS ---------------------------------- #

# Speed of light in vacuum
LIGHTSPEED = 299792458  # m/s

# Planck's constant
PLANCK_H = 6.62607015e-34  # J.s

# Coefficient for the conversion from flam (erg/cm^2/s/A) to photlam (photon/cm^2/s/A).
# From <https://www.stsci.edu/~strolger/docs/UNITS.txt>
FLAM_TO_PHOTLAM = 5.03411250e7  # photons/erg/A

# Coefficient for the conversion from fnu (erg/cm^2/s/Hz) to photlam (photon/cm^2/s/A).
# From <https://hea-www.harvard.edu/~pgreen/figs/Conversions.pdf>
FNU_TO_PHOTLAM = 1.51e26  # photon.Hz/erg

# Coefficient for the conversion from fnu (erg/cm^2/s/Hz) to flam (erg/cm^2/s/A).
# Reference: <https://hea-www.harvard.edu/~pgreen/figs/Conversions.pdf>
FNU_TO_FLAM = 2.99792458e18  # Hz/A

# The flux of the geocoronal emission line [O II] 2471A.
# See <https://hst-docs.stsci.edu/wfc3ihb/chapter-9-wfc3-exposure-time-calculation/9-7-sky-background#id-9.7SkyBackground-9.7.29.7.2GeocoronalEmission,Airglow,andSHADOW>
GEOCORONAL_FLUX = 1.5e-15  # erg/cm^2/s/arcsec^2
GEOCORONAL_WAVELENGTH = 2471  # angstroms

# -------------------------------------------------------------------------------------- #
