"""
constants.py

A file to track filepaths and telescope parameters. N.B. use cgs units when possible.

Isaac Cheng - 2022
"""

import astropy.units as u

# ----------------------------------- SOME CONSTANTS ----------------------------------- #

# Speed of light in vacuum
LIGHTSPEED = 2.99792458e10 << (u.cm / u.s)  # cm/s

# Planck's constant
PLANCK_H = 6.62607015e-27 << (u.erg * u.s)  # erg.s

# Boltzmann's constant
K_B = 1.380649e-16 << (u.erg / u.K)  # erg/K

# The radius of the Sun
SUN_RADIUS = 6.96340e10 << u.cm  # cm

# Coefficient for the conversion from flam (erg/cm^2/s/A) to photlam (photon/cm^2/s/A).
# Reference: <https://hea-www.harvard.edu/~pgreen/figs/Conversions.pdf>
FLAM_TO_PHOTLAM = 5.0341165675427094e7  # photon/erg/A

# Coefficient for the conversion from fnu (erg/cm^2/s/Hz) to photlam (photon/cm^2/s/A).
# Reference: <https://hea-www.harvard.edu/~pgreen/figs/Conversions.pdf>
FNU_TO_PHOTLAM = 1.5091901796421519e26  # photon.Hz/erg

# Coefficient for the conversion from fnu (erg/cm^2/s/Hz) to flam (erg/cm^2/s/A).
# Reference: <https://hea-www.harvard.edu/~pgreen/figs/Conversions.pdf>
FNU_TO_FLAM = 2.99792458e18  # Hz/A

# Coefficient for the conversion from flam (erg/cm^2/s/A) to fnu (erg/cm^2/s/Hz).
# Reference: <https://hea-www.harvard.edu/~pgreen/figs/Conversions.pdf>
FLAM_TO_FNU = 3.3356409519815204e-19  # A/Hz

# 1 steradian to 1 arcsec^2
SR_TO_SQARCSEC = 4.254517029615221e10  # arcsec^2/sr

# 1 parsec
PC = 3.0856775815e18 << u.cm  # cm

# -------------------------------------------------------------------------------------- #
