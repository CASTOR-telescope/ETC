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


"""
import unittest

import astropy.units as u
import numpy as np

from castor_etc.background import Background
from castor_etc.telescope import Telescope
from castor_etc.photometry import Photometry

_TOL = 1e-15  # floating-point tolerance

class PhotometryTestCase(unittest.TestCase):
    """
    Integrated test suite to test different photometry calculations
    """
    def setUp(self):
        # Default Telescope parameters
        self.scope = Telescope()

        # Default background with one emission line
        self.bg = Background()
        self.bg.add_geocoronal_emission(flux="high") # add a high flux default emission


    def test_point_source_photometry(self):
        """
        Basic integration test to ensure that the getting started Jupyter notebook works.
        """

        from castor_etc.sources import PointSource

        # Generate the black body from the example and add associated emission lines
        src = PointSource()
        src.generate_bb(8000 * u.K, redshift=0.06, limits=[900, 30000] * u.AA)
        src.norm_to_AB_mag(25)
        src.add_emission_line(
            center=2000 * u.AA, fwhm=200 * u.AA, peak=5e-19, shape="gaussian", abs_peak=False
        )
        src.add_absorption_line(
            center=5005 * u.AA, fwhm=40 * u.AA, dip=2e-19, shape="lorentzian", abs_dip=True
        )

        # Initialize the photometry class
        phot = Photometry(self.scope, src, self.bg)
        phot.use_optimal_aperture()

        # Confirm that the AB magnitude in each passband is consistent
        ## Each test is conducted 
        mags = src.get_AB_mag(self.scope)
        self.assertAlmostEqual(mags['uv'], np.float64(26.44580324104077),delta=_TOL) # check the uv component to 3rd decimal
        self.assertAlmostEqual(mags['u'], np.float64(24.964036810739323), delta=_TOL) # check the u component
        self.assertAlmostEqual(mags['g'], np.float64(24.391672332311423), delta=_TOL)
        
        TARGET_SNR = 10
        REDDENING = 0.01  # E(B-V)

        # Test exposure times
        exp_t = phot.calc_snr_or_t(snr=TARGET_SNR, reddening=REDDENING, quiet=True)
        self.assertAlmostEqual(exp_t['uv'], np.float64(1052.8014187348936), delta=_TOL)
        self.assertAlmostEqual(exp_t['u'], np.float64(248.21615426248252), delta=_TOL)
        self.assertAlmostEqual(exp_t['g'], np.float64(130.0356738505753), delta=_TOL)

        # Test SNR calculation for given time
        snr = phot.calc_snr_or_t(t=exp_t, reddening=REDDENING, quiet=True)
        self.assertAlmostEqual(snr['uv'], np.float64(10), delta=_TOL)
        self.assertAlmostEqual(snr['u'], np.float64(10), delta=_TOL)
        self.assertAlmostEqual(snr['g'], np.float64(10), delta=_TOL)


        





