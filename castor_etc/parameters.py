"""
Contains parameters for CASTOR imaging chain.

---

        GNU General Public License v3 (GNU GPLv3)

(c) 2022.                            (c) 2022.
Government of Canada                 Gouvernement du Canada
National Research Council            Conseil national de recherches
Ottawa, Canada, K1A 0R6              Ottawa, Canada, K1A 0R6
All rights reserved                  Tous droits réservés

NRC disclaims any warranties,        Le CNRC dénie toute garantie
expressed, implied, or               énoncée, implicite ou légale,
statutory, of any kind with          de quelque nature que ce
respect to the software,             soit, concernant le logiciel,
including without limitation         y compris sans restriction
any warranty of merchantability      toute garantie de valeur
or fitness for a particular          marchande ou de pertinence
purpose. NRC shall not be            pour un usage particulier.
liable in any event for any          Le CNRC ne pourra en aucun cas
damages, whether direct or           être tenu responsable de tout
indirect, special or general,        dommage, direct ou indirect,
consequential or incidental,         particulier ou général,
arising from the use of the          accessoire ou fortuit, résultant
software. Neither the name           de l'utilisation du logiciel. Ni
of the National Research             le nom du Conseil National de
Council of Canada nor the            Recherches du Canada ni les noms
names of its contributors may        de ses  participants ne peuvent
be used to endorse or promote        être utilisés pour approuver ou
products derived from this           promouvoir les produits dérivés
software without specific prior      de ce logiciel sans autorisation
written permission.                  préalable et particulière
                                     par écrit.

This file is part of the             Ce fichier fait partie du projet
FORECASTOR ETC project.              FORECASTOR ETC.

FORECASTOR ETC is free software:     FORECASTOR ETC est un logiciel
you can redistribute it and/or       libre ; vous pouvez le redistribuer
modify it under the terms of         ou le modifier suivant les termes de
the GNU General Public               la "GNU General Public
License as published by the          License" telle que publiée
Free Software Foundation,            par la Free Software Foundation :
either version 3 of the              soit la version 3 de cette
License, or (at your option)         licence, soit (à votre gré)
any later version.                   toute version ultérieure.

FORECASTOR ETC is distributed        FORECASTOR ETC est distribué
in the hope that it will be          dans l'espoir qu'il vous
useful, but WITHOUT ANY WARRANTY;    sera utile, mais SANS AUCUNE
without even the implied warranty    GARANTIE : sans même la garantie
of MERCHANTABILITY or FITNESS FOR    implicite de COMMERCIALISABILITÉ
A PARTICULAR PURPOSE. See the        ni d'ADÉQUATION À UN OBJECTIF
GNU General Public License for       PARTICULIER. Consultez la Licence
more details.                        Générale Publique GNU pour plus
                                     de détails.

You should have received             Vous devriez avoir reçu une
a copy of the GNU General            copie de la Licence Générale
Public License along with            Publique GNU avec FORECASTOR ETC ;
FORECASTOR ETC. If not, see          si ce n'est pas le cas, consultez :
<http://www.gnu.org/licenses/>.      <http://www.gnu.org/licenses/>.
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

# Resolution of the passband response curves
PASSBAND_RESOLUTION = 1 << u.nm  # nanometres

# Filepaths to the passband response curves
PASSBAND_FILEPATHS = {
    band: join(DATAPATH, "passbands", f"passband_castor.{band}") for band in PASSBANDS
}

# Wavelength units in the passband response curve files
PASSBAND_FILEUNITS = {"uv": u.um, "u": u.um, "g": u.um}

# Filepaths to the passband point spread functions (PSFs)
PSF_FILEPATHS = {
    band: join(DATAPATH, "psfs", f"{band}_psf_20x_supersampled.fits")
    for band in PASSBANDS
}

# The PSF oversampling factor. Each square pixel in the PSF file has a side length of
# PX_SCALE / PSF_SUPERSAMPLE_FACTOR (e.g., 0.1 arcsec / 20 = 0.005 arcsec)
PSF_SUPERSAMPLE_FACTOR = 20

# The PSF full-width at half-maximum
FWHM = 0.15 << u.arcsec  # arcsec

# The angular dimension covered by each pixel
PX_SCALE = 0.1 << u.arcsec  # arcsec (or arcsec/pixel)

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

# (Teledyne e2v CMOS detector)
# Dark current is 0.01 electrons/s/pixel at -50°C and halves for every reduction of 5-6°C.
# CASTOR operates at 180 K, implying:
# dark current = 0.5^((223.15K - 180K) / 6) * 0.01 ~ 1e-4 electrons/s/pixel (negligible)
DARK_CURRENT = 1e-4  # electrons/s/pixel

BIAS = 100  # electron

READ_NOISE = 2.0  # electron/pixel (high-gain). Read noise is 30 electrons for low-gain

GAIN = 2.0  # electron/ADU

# Wavelength threshold for red leak. Flux longward of this is considered to be red leak
REDLEAK_THRESHOLDS = {"uv": 3010 << u.AA, "u": 4160 << u.AA, "g": 5600 << u.AA}

# Extinction coefficients (i.e., R := A/E(B-V)) for the different telescope passbands
# Estimates taken from Yuan+2013, column 3 of Table 2, rows: NUV, (SDSS) u, (SDSS) g
# <https://ui.adsabs.harvard.edu/abs/2013MNRAS.430.2188Y/abstract>
EXTINCTION_COEFFS = {"uv": 7.06, "u": 4.35, "g": 3.31}

# -------------------------------------------------------------------------------------- #
