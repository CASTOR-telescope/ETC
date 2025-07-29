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
test_energy.py

Unit tests for the energy module.

OUTDATED AND NO LONGER RELEVANT TO CURRENT VERSION.
"""
import unittest

import astropy.units as u
import numpy as np

from ..conversions import calc_photon_energy

_TOL = 1e-15  # floating-point tolerance


class TestEnergy(unittest.TestCase):
    """
    Test castor_etc.energy module.
    """

    def test_calc_photon_energy_wavelength_single(self):
        """
        Tests the calc_photon_energy method for a single wavelength.
        """
        results = calc_photon_energy(wavelength=100 * u.nm, wavelength_err=10 * u.nm)
        #
        true_energy = (100 * u.nm).to(u.erg, equivalencies=u.spectral()).value
        true_uncer = (10 / 100) * true_energy
        truths = [true_energy, true_uncer]
        #
        for result, truth in zip(results, truths):
            self.assertAlmostEqual(result, truth, delta=_TOL)

    def test_calc_photon_energy_wavelength_multi(self):
        """
        Tests the calc_photon_energy method for an array of wavelengths.
        """
        wavelengths = np.array([0.100, 2, 30]) * u.um
        wavelength_errs = np.array([1, 0.2, 30]) * u.um
        #
        results = calc_photon_energy(
            wavelength=wavelengths, wavelength_err=wavelength_errs
        )
        #
        true_energy = wavelengths.to(u.erg, equivalencies=u.spectral()).value
        true_uncer = (wavelength_errs / wavelengths).value * true_energy
        truths = [true_energy, true_uncer]
        #
        for result_arr, truth_arr in zip(results, truths):
            for result, truth in zip(result_arr, truth_arr):
                self.assertAlmostEqual(result, truth, delta=_TOL)

    def test_calc_photon_energy_frequency_single(self):
        """
        Tests the calc_photon_energy method for a single frequency.
        """
        raise NotImplementedError("Not implemented yet")

    def test_calc_photon_energy_frequency_multi(self):
        """
        Tests the calc_photon_energy method for an array of frequencies.
        """
        raise NotImplementedError("Not implemented yet")


if __name__ == "__main__":
    unittest.main()
