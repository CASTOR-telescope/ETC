"""
Utilities to calculate energies.
"""

import astropy.units as u

from . import constants as const


def calc_photon_energy(
    wavelength=None, frequency=None, wavelength_err=0.0, frequency_err=0.0
):
    """
    Calculates the energy of a photon in ergs given its wavelength or frequency.

    Parameters
    ----------
      wavelength, frequency :: scalar or `astropy.Quantity` or array
        The wavelength and frequency of the photon; only one of these should be provided.
        If values are scalars, wavelength is assumed to be in angstrom and frequency in
        hertz.

      wavelength_err, frequency_err :: scalar or `astropy.Quantity` or array
        The uncertainty in the wavelength and uncertainty; only one of these should be
        provided. If values are scalars, wavelength_err is assumed to be in angstrom and
        frequency_err in hertz.

    Returns
    -------
      energy, energy_err :: scalar or array
        The energy of the photon and its uncertainty, in ergs (1 erg = 10^-7 joule).
    """
    if wavelength is not None and frequency is not None:
        raise ValueError("Only one of wavelength and frequency should be provided")
    elif wavelength is not None:
        if isinstance(wavelength, u.Quantity):
            wavelength = wavelength.to(u.AA).value
        if isinstance(wavelength_err, u.Quantity):
            wavelength_err = wavelength_err.to(u.AA).value
        energy = (const.PLANCK_H * const.LIGHTSPEED / (wavelength * 1e-8)).value  # erg
        energy_err = energy * (wavelength_err / wavelength)
    elif frequency is not None:
        if isinstance(frequency, u.Quantity):
            frequency = frequency.to(u.Hz).value
        if isinstance(frequency_err, u.Quantity):
            frequency_err = frequency_err.to(u.Hz).value
        energy = (const.PLANCK_H * frequency).value  # erg
        energy_err = energy * (frequency_err / frequency)
    else:
        raise ValueError("Either wavelength or frequency must be provided")
    return energy, energy_err  # erg
