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
test_conversions.py

Unit tests for the spectrum module.

OUTDATED AND NO LONGER RELEVANT TO CURRENT VERSION.
"""

import unittest

import astropy.units as u
import numpy as np

from .. import spectrum

# ALlowed differences for passing tests
_BB_TOL = 1e-13  # floating-point tolerance for BB spectrum


class TestSpectrum(unittest.TestCase):
    """
    Test castor_etc.spectrum module.
    """

    def test_generate_BB(self):
        """
        Tests the generate_BB() method.
        """
        from astropy.modeling.models import BlackBody

        #
        # Test ETC blackbody model with automatic wavelength generation
        #
        wavelengths, bb_spectrum = spectrum.generate_BB(
            T=10000, limits=[100, 4000] * u.AA, resolution=100 * u.AA
        )
        # Astropy blackbody model
        bb = BlackBody(10000 * u.K)
        bb_spectrum2 = bb(wavelengths * u.AA).to(
            u.erg / (u.s * u.cm ** 2 * u.AA * u.sr),
            equivalencies=u.spectral_density(wavelengths * u.AA),
        )
        bb_spectrum2 /= (wavelengths * u.AA).to(u.erg, equivalencies=u.spectral())
        bb_spectrum2 = bb_spectrum2.value
        # Tests
        self.assertTrue(np.allclose(bb_spectrum, bb_spectrum2, atol=_BB_TOL))
        relative_errs = np.abs(bb_spectrum - bb_spectrum2) / bb_spectrum2
        for relative_err in relative_errs:
            self.assertAlmostEqual(0, relative_err, delta=_BB_TOL)

        #
        # Test ETC blackbody model with input wavelengths
        #
        wavelengths_new = np.linspace(100, 4000, 100) << u.nm
        bb_spectrum_new = spectrum.generate_BB(
            T=3000 * u.Celsius, wavelengths=wavelengths_new
        )[1]
        # Astropy blackbody model
        bb_new = BlackBody(3000 * u.Celsius)
        bb_spectrum_new2 = bb_new(wavelengths_new).to(
            u.erg / (u.s * u.cm ** 2 * u.AA * u.sr),
            equivalencies=u.spectral_density(wavelengths_new),
        )
        bb_spectrum_new2 /= (wavelengths_new).to(u.erg, equivalencies=u.spectral())
        bb_spectrum_new2 = bb_spectrum_new2.value
        # Tests
        self.assertTrue(np.allclose(bb_spectrum_new, bb_spectrum_new2, atol=_BB_TOL))
        relative_errs_new = np.abs(bb_spectrum_new - bb_spectrum_new2) / bb_spectrum_new2
        for relative_err_new in relative_errs_new:
            self.assertAlmostEqual(0, relative_err_new, delta=_BB_TOL)

    # TODO: Add tests for other spectra and normalization methods


if __name__ == "__main__":
    unittest.main()
