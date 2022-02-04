"""
test_conversions.py

Unit tests for the spectrum module.

Isaac Cheng - 2022
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
