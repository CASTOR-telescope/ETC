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
test_photometry.py

This module contains an integrated test suite to test different photometry calculations.

This is modified from the getting started Jupyter notebook to ensure that the examples there work as expected.
"""

import unittest

import astropy.units as u
import numpy as np

from castor_etc.background import Background
from castor_etc.telescope import Telescope
from castor_etc.photometry import Photometry

_TOL = 1e-2

class GalaxySourcePhotometryTestCase(unittest.TestCase):
    """
    Integrated test suite to test different photometry calculations
    """
    def setUp(self):
        # Default Telescope parameters
        self.scope = Telescope()

        # Default background with one emission line
        self.bg = Background()  # default Earthshine and zodiacal light
        self.bg.add_geocoronal_emission(
                flux=1e-15, wavelength=2345 * u.AA, linewidth=0.023 * u.AA)
              
        from castor_etc.sources import GalaxySource


        self.src = GalaxySource(
            r_eff=3 * u.arcsec,  # effective radius, sqrt(a * b)
            n=4,  # Sersic index
            axial_ratio=0.9,  # ratio of semiminor axis (b) to semimajor axis (a)
            rotation=135,  # CCW rotation from x-axis
        )
        self.src.use_galaxy_spectrum(gal_type="spiral")
        self.src.norm_luminosity_dist(luminosity=1.4e8, dist=450 * u.Mpc)  # Luminosity in solar luminosities
        self.src.redshift_wavelengths(0.1)

        self.phot = Photometry(self.scope, self.src, self.bg)
        self.phot.use_elliptical_aperture(
            a=6 * u.arcsec, b=4 * u.arcsec, center=[0, 0] * u.arcsec, rotation=31.41592654
        )

    def test_absolute_magnitude_calculation(self):
        # NOTE: This may make more sense in a test for the GalaxySource class, will move later
        abs_mag = self.src.get_AB_mag(self.scope)
        self.assertAlmostEqual(abs_mag['uv'], np.float64(25.69136659200104), delta=_TOL)
        self.assertAlmostEqual(abs_mag['u'], np.float64(25.06673659547406), delta=_TOL)
        self.assertAlmostEqual(abs_mag['g'], np.float64(23.609216549534842), delta=_TOL)

        self.assertAlmostEqual(self.src.get_AB_mag(), np.float64(23.39051359135356), delta=_TOL )

    def test_expectected_spectrum_plot(self):
        self.src.show_spectrum(plot=False)
        pass  # Just ensure no exceptions are raised

    def test_galaxy_source_snr_calculator(self):
        INTEGRATION_TIME = 4321  # seconds
        REDDENING = 0.01

        self.assertAlmostEqual(self.phot.calc_snr_or_t(t=INTEGRATION_TIME, reddening=REDDENING, quiet=True)['uv'], np.float64(1.9289396333476412), delta=_TOL)
        self.assertAlmostEqual(self.phot.calc_snr_or_t(t=INTEGRATION_TIME, reddening=REDDENING, quiet=True)["u"], np.float64(2.311967051585148), delta=_TOL)
        self.assertAlmostEqual(self.phot.calc_snr_or_t(t=INTEGRATION_TIME, reddening=REDDENING, quiet=True)["g"], np.float64(4.912896778606595), delta=_TOL)


    def test_galaxy_source_exposure_time_calculator(self):
        TARGET_SNR = {
            'uv': 1.9289396333476412,
            'u': 2.311967051585148,
            'g': 4.912896778606595,
        }
        REDDENING = 0.01

        # Test exposure times
        self.assertAlmostEqual(self.phot.calc_snr_or_t(snr=TARGET_SNR["uv"], reddening=REDDENING, quiet=True)['uv'], np.float64(4321), delta=_TOL)
        self.assertAlmostEqual(self.phot.calc_snr_or_t(snr=TARGET_SNR["u"], reddening=REDDENING, quiet=True)['u'], np.float64(4321), delta=_TOL)
        self.assertAlmostEqual(self.phot.calc_snr_or_t(snr=TARGET_SNR["g"], reddening=REDDENING, quiet=True)['g'], np.float64(4321), delta=_TOL)

if __name__ == '__main__':
    unittest.main()