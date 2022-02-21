"""
test_energy.py

Unit tests for the energy module.

Copyright 2022, CASTOR Mission Team
Author: Isaac Cheng
Contact: isaac.cheng.ca@gmail.com
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
