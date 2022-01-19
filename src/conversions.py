"""
conversions.py

Utilities to convert between different units.

Isaac Cheng - January 2022
"""

import astropy.units as u
import numpy as np

from . import constants as const


def convert_freq_wavelength(data, to="wavelength", output_unit=u.AA):
    """
    Converts between frequency and wavelength for light in a vacuum.

    Parameters
    ----------
      data :: scalar or `astropy.Quantity` or array
        The frequency of wavelength data to convert. If values are scalars, they are
        assumed to be in Hz for frequencies and in angstroms for wavelengths.

      to :: "wavelength" or "frequency"
        The quantity to convert the input data to.

      output_unit :: `astropy.Quantity`
        The unit of the returned, converted data.

    Returns
    -------
      converted_data :: scalar or array
        The converted data.
    """
    #
    # Check inputs
    #
    if to == "wavelength":
        unit = u.Hz  # data are frequencies
    elif to == "frequency":
        unit = u.AA  # data are wavelengths
    else:
        raise ValueError(f"to must be 'wavelength' or 'frequency'")
    #
    if isinstance(data, u.Quantity):
        data = data.to(unit)
    else:
        data = data * unit
    #
    # Convert
    #
    lightspeed = const.LIGHTSPEED * u.m / u.s
    converted_data = (lightspeed / data).to(output_unit).value
    return converted_data


def flam_to_photlam(flam, wavelength):
    """
    Converts from flam (erg/cm^2/s/A) to photlam (photon/cm^2/s/A).
    See <https://www.stsci.edu/~strolger/docs/UNITS.txt>.

    Parameters
    ----------
      flam :: array
        The flux in flam.

      wavelength :: array of scalars or `astropy.Quantity`
        The corresponding wavelengths of the flux. If values are scalars, they are assumed
        to be in angstroms.

    Returns
    -------
      photlam :: array
        The flux in photlam.
    """
    if isinstance(wavelength, u.Quantity):
        wavelength = wavelength.to(u.AA).value
    return const.FLAM_TO_PHOTLAM * flam * wavelength


def fnu_to_photlam(fnu, wavelength):
    """
    Converts from fnu (erg/cm^2/s/Hz) to photlam (photom/cm^2/s/A). See
    <https://hea-www.harvard.edu/~pgreen/figs/Conversions.pdf> ("F_\nu" to "f_\lambda").

    Parameters
    ----------
      fnu :: array
        The flux in fnu.

      wavelength :: array of scalars or `astropy.Quantity`
        The corresponding wavelengths of the flux. If values are scalars, they are assumed
        to be in angstroms.

    Returns
    -------
      photlam :: array
        The flux in photlam.
    """
    if isinstance(wavelength, u.Quantity):
        wavelength = wavelength.to(u.AA).value
    return const.FNU_TO_PHOTLAM * fnu / wavelength


def fnu_to_flam(fnu, wavelength):
    """
    Converts from fnu (erg/cm^2/s/Hz) to flam (erg/cm^2/s/A). See
    <https://hea-www.harvard.edu/~pgreen/figs/Conversions.pdf> ("F_\nu" to "F_\lambda").

    Parameters
    ----------
      fnu :: array
        The flux in fnu.

      wavelength :: array of scalars or `astropy.Quantity`
        The corresponding wavelengths of the flux. If values are scalars, they are assumed
        to be in angstroms.

    Returns
    -------
      photlam :: array
        The flux in photlam.
    """
    if isinstance(wavelength, u.Quantity):
        wavelength = wavelength.to(u.AA).value
    return const.FNU_TO_FLAM * fnu / (wavelength * wavelength)


def convert_rel_abs_mag(mag, dist, mag_err=0.0, dist_err=0.0, to="abs"):
    """
    Converts between relative and absolute magnitudes.

    Parameters
    ----------
      mag :: scalar or array
        The relative or absolute magnitudes.

      dist :: scalar or `astropy.Quantity` or array
        The distance to the source. If values are scalars, they are assumed to be in
        parsecs.

      mag_err :: scalar or array
        The uncertainty in the relative or absolute magnitudes.

      dist_err :: scalar or `astropy.Quantity` or array
        The uncertainty in the distance to the source. If values are scalars, they are
        assumed to be in parsecs.

      to :: "abs" or "rel"
        If "abs", convert input magnitudes to absolute magnitudes. If "rel", convert
        input magnitudes to relative magnitudes.

    Returns
    -------
      new_mag, rel_mag_err :: scalar or array
        The converted magnitudes and their uncertainties.
    """
    #
    # Check inputs
    #
    if isinstance(dist, u.Quantity):
        dist = dist.to(u.parsec).value
    if isinstance(dist_err, u.Quantity):
        dist_err = dist_err.to(u.parsec).value
    #
    # Convert magnitudes
    #
    if to == "abs":
        new_mag = mag - 5 * (np.log10(dist) - 1)
    elif to == "rel":
        new_mag = mag + 5 * (np.log10(dist) - 1)
    else:
        raise ValueError("to must be either 'abs' or 'rel'")
    new_mag_err = np.sqrt(mag_err ** 2 + (5 / np.log(10) * dist_err / dist) ** 2)
    #
    return new_mag, new_mag_err


def flux_to_mag(flux, flux_err=0.0, zpt=-48.60, calc_abs=False, dist=None, dist_err=0.0):
    """
    Calculates the relative or absolute magnitude of a source and its uncertainty given
    the monochromatic flux of the source.

    Parameters
    ----------
      flux, flux_err :: scalar or arrays
        The monochromatic flux and its uncertainty. The unit of the flux depends on the
        magnitude system.

      zpt :: scalar
        The zero point of the magnitude system.

      calc_abs :: bool
        If True, calculate the absolute magnitude; also requires the dist parameter to be
        provided. If False, calculate the relative magnitude.

      dist, dist_err :: scalar or `astropy.Quantity` or array
        The distance to the source and its uncertainty. If values are scalars, they are
        assumed to be in parsecs.

    Returns
    -------
      mag, mag_err :: float or array
        The relative or absolute magnitude and its uncertainty.
    """
    #
    # Calculate magnitude
    #
    rel_mag = -2.5 * np.log10(flux) + zpt
    rel_mag_err = 2.5 / np.log(10) * abs(flux_err / flux)
    #
    if calc_abs:
        if dist is None:
            raise ValueError("dist must be provided if calc_abs is True")
        abs_mag, abs_mag_err = convert_rel_abs_mag(
            rel_mag, dist, mag_err=rel_mag_err, dist_err=dist_err, to="abs"
        )
        return abs_mag, abs_mag_err
    #
    return rel_mag, rel_mag_err


def mag_to_flux(mag, mag_err=0.0, zpt=-48.60):
    """
    Converts magnitude to flux. The units of the flux will depend on the magnitude system.

    Parameters
    ----------
      mag, mag_err :: scalar or arrays
        The magnitudes and their uncertainties

      zpt :: scalar
        The zero point of the magnitude system. For example, zpt=-48.60 corresponds to the
        AB magnitude system and the flux will be in units of erg/s/cm^2/Hz. Likewise,
        zpt=-21.1 corresponds to the ST magnitude system and the flux will be in units of
        erg/s/cm^2/A.
    """
    flux = 10 ** (-0.4 * (mag - zpt))
    flux_err = 0.4 * np.log(10) * abs(flux * mag_err)
    return flux, flux_err
