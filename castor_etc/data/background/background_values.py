"""
background_values.py

Contains average sky background values and geocoronal emission line data.

Isaac Cheng - 2022
"""

import astropy.units as u

# ------------------------------- BACKGROUND NOISE VALUES ------------------------------ #

# The flux of the geocoronal emission line [O II] 2471A.
# See <https://hst-docs.stsci.edu/wfc3ihb/chapter-9-wfc3-exposure-time-calculation/9-7-sky-background#id-9.7SkyBackground-9.7.29.7.2GeocoronalEmission,Airglow,andSHADOW>
# and <https://www.stsci.edu/itt/APT_help20/STIS/c06_exptime7.html#696464>
GEOCORONAL_FLUX = 1.5e-15  # erg/cm^2/s/arcsec^2
GEOCORONAL_WAVELENGTH = 2471 << u.AA  # angstroms
GEOCORONAL_LINEWIDTH = 0.023 << u.AA  # angstroms

# Average sky background in each passband (Earthshine + zodiacal light)
SKY_BACKGROUND = {"uv": 26.08, "u": 23.74, "g": 22.60}  # AB mag/arcsec^2

# -------------------------------------------------------------------------------------- #
