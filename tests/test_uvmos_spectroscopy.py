#         GNU General Public License v3 (GNU GPLv3)
#
# (c) 2022.                            (c) 2022.
# Government of Canada                 Gouvernement du Canada
# National Research Council            Conseil national de recherches
# Ottawa, Canada, K1A 0R6              Ottawa, Canada, K1A 0R6
# All rights reserved                  Tous droits réservés
#
# NRC disclaims any warranties,        Le CNRC dénie toute garantie
# expressed, implied, or               énoncée, implicite ou légale,
# statutory, of any kind with          de quelque nature que ce
# respect to the software,             soit, concernant le logiciel,
# including without limitation         y compris sans restriction
# any warranty of merchantability      toute garantie de valeur
# or fitness for a particular          marchande ou de pertinence
# purpose. NRC shall not be            pour un usage particulier.
# liable in any event for any          Le CNRC ne pourra en aucun cas
# damages, whether direct or           être tenu responsable de tout
# indirect, special or general,        dommage, direct ou indirect,
# consequential or incidental,         particulier ou général,
# arising from the use of the          accessoire ou fortuit, résultant
# software. Neither the name           de l'utilisation du logiciel. Ni
# of the National Research             le nom du Conseil National de
# Council of Canada nor the            Recherches du Canada ni les noms
# names of its contributors may        de ses  participants ne peuvent
# be used to endorse or promote        être utilisés pour approuver ou
# products derived from this           promouvoir les produits dérivés
# software without specific prior      de ce logiciel sans autorisation
# written permission.                  préalable et particulière
#                                      par écrit.
#
# This file is part of the             Ce fichier fait partie du projet
# FORECASTOR ETC project.              FORECASTOR ETC.
#
# FORECASTOR ETC is free software:     FORECASTOR ETC est un logiciel
# you can redistribute it and/or       libre ; vous pouvez le redistribuer
# modify it under the terms of         ou le modifier suivant les termes de
# the GNU General Public               la "GNU General Public
# License as published by the          License" telle que publiée
# Free Software Foundation,            par la Free Software Foundation :
# either version 3 of the              soit la version 3 de cette
# License, or (at your option)         licence, soit (à votre gré)
# any later version.                   toute version ultérieure.
#
# FORECASTOR ETC is distributed        FORECASTOR ETC est distribué
# in the hope that it will be          dans l'espoir qu'il vous
# useful, but WITHOUT ANY WARRANTY;    sera utile, mais SANS AUCUNE
# without even the implied warranty    GARANTIE : sans même la garantie
# of MERCHANTABILITY or FITNESS FOR    implicite de COMMERCIALISABILITÉ
# A PARTICULAR PURPOSE. See the        ni d'ADÉQUATION À UN OBJECTIF
# GNU General Public License for       PARTICULIER. Consultez la Licence
# more details.                        Générale Publique GNU pour plus
#                                      de détails.
#
# You should have received             Vous devriez avoir reçu une
# a copy of the GNU General            copie de la Licence Générale
# Public License along with            Publique GNU avec FORECASTOR ETC ;
# FORECASTOR ETC. If not, see          si ce n'est pas le cas, consultez :
# <http://www.gnu.org/licenses/>.      <http://www.gnu.org/licenses/>.


"""
test_uvmos_spectroscopy.py

This module contains an integrated test suite to the UVMOS spectroscopy calculations.
"""

import unittest

import astropy.units as u
import matplotlib.pyplot as plt

from castor_etc.background import Background
from castor_etc.sources import PointSource
from castor_etc.telescope import Telescope

from castor_etc.uvmos_spectroscopy import UVMOS_Spectroscopy

_TOL = 1e-5  # floating-point tolerance

class UVMOS_SpectroscopyTestCase(unittest.TestCase):
    """
    Integrated test suite to test different photometry calculations
    """
    def setUp(self):
        # Default Telescope parameters
        self.scope = Telescope()

        # Default background with one emission line
        self.bg = Background(mags_per_sq_arcsec={"uv": 26.08, "u": 23.74, "g": 22.60})
        from castor_etc.sources import ExtendedSource


        self.src = PointSource()
        self.src.use_pickles_spectrum("b0v")
        self.src.norm_to_AB_mag(1) # check this!

        # Create spectroscopy object

        self.spec = UVMOS_Spectroscopy(self.scope, self.src, self.bg)
        self.spec.specify_slit()

        self.spec.calc_source_CASTORSpectrum(extraction_width=1, extraction_lowerlim=4, extraction_upperlim = 7)
        self.spec.calc_background_CASTORSpectrum(extraction_width=1, extraction_lowerlim=4, extraction_upperlim = 7)

    def test_calc_snr(self):
        T_TARGET = 100  # seconds
        LAMBDA = 2000   # Angstroms

        snr = self.spec.calc_snr_from_t(T_TARGET, LAMBDA)
        
        # TODO: add assertions

    def test_calc_exposure_time(self):
        SNR_TARGET = 10
        LAMBDA = 2000

        t = self.spec.calc_t_from_snr(SNR_TARGET, LAMBDA)

if __name__ == '__main__':
    unittest.main()