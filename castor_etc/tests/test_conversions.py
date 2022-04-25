"""
test_conversions.py

Unit tests for the conversions module.
"""

import unittest

import astropy.units as u
import numpy as np

from .. import conversions as convert

# ALlowed differences for passing tests
_FLUX_MAG_COUNT = 1e-15  # floating-point tolerance for convert.convert_count_flux_mag()


class TestConversions(unittest.TestCase):
    """
    Test castor_etc.conversions module.
    """

    def test_convert_count_flux_mag(self):
        """
        Tests the convert_count_flux_mag() method by converting between every combination
        of counts, flux, and AB magnitudes.

        TODO: add tests!
        """
        raise NotImplementedError("Tests not implemented yet!")


if __name__ == "__main__":
    unittest.main()
