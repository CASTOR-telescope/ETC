"""
conversions.py

Utilities to convert between different units.

Common units (also see <https://pysynphot.readthedocs.io/en/latest/units.html>):
  - fnu :: erg/s/cm^2/Hz
  - flam :: erg/s/cm^2/A
  - photnu :: photons/s/cm^2/Hz
  - photlam :: photons/s/cm^2/A

Isaac Cheng - 2022
"""

from json.tool import main
from multiprocessing.sharedctypes import Value
import astropy.units as u
import numpy as np

from . import constants as const
from . import parameters as params


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
    converted_data = (const.LIGHTSPEED / data).to(output_unit).value
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
    Converts from fnu (erg/cm^2/s/Hz) to photlam (photon/cm^2/s/A). See
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


def fnu_to_flam(fnu, wavelength, fnu_err=0.0, wavelength_err=0.0):
    """
    Converts from fnu (erg/cm^2/s/Hz) to flam (erg/cm^2/s/A). See
    <https://hea-www.harvard.edu/~pgreen/figs/Conversions.pdf> ("F_\nu" to "F_\lambda").

    Parameters
    ----------
      fnu :: array of scalar
        The flux in fnu.

      wavelength :: array of scalars or `astropy.Quantity` array
        The corresponding wavelengths of the flux. If values are scalars, they are assumed
        to be in angstroms.

      fnu_err :: array of scalars
        The absolute uncertainty in fnu.

      wavelength_err :: array of scalars or `astropy.Quantity` array
        The absolute uncertainty in wavelength. If values are scalars, they are assumed
        to be in angstroms.

    Returns
    -------
      flam, flam_err :: arrays of floats
        The flux in flam and its uncertainty.
    """
    #
    # Check inputs
    #
    if isinstance(wavelength, u.Quantity):
        wavelength = wavelength.to(u.AA).value
    if isinstance(wavelength_err, u.Quantity):
        wavelength_err = wavelength_err.to(u.AA).value
    #
    # Convert fnu (erg/cm^2/s/Hz) to flam (erg/cm^2/s/A)
    #
    flam = const.FNU_TO_FLAM * fnu / (wavelength * wavelength)
    # REVIEW: check error propagation is correct
    flam_err = flam * np.sqrt(
        (fnu_err / fnu) ** 2 + 4 * (wavelength_err / wavelength) ** 2
    )  # derived & simplified from partial derivative error propagation
    return flam, flam_err


def flam_to_fnu(flam, wavelength, flam_err=0.0, wavelength_err=0.0):
    """
    Converts from flam (erg/cm^2/s/A) to fnu (erg/cm^2/s/Hz). See
    <https://hea-www.harvard.edu/~pgreen/figs/Conversions.pdf> ("F_\lambda" to "F_\nu").

    Parameters
    ----------
      flam :: array of scalar
        The flux in flam.

      wavelength :: array of scalars or `astropy.Quantity` array
        The corresponding wavelengths of the flux. If values are scalars, they are assumed
        to be in angstroms.

      flam_err :: array of scalars
        The absolute uncertainty in flam.

      wavelength_err :: array of scalars or `astropy.Quantity` array
        The absolute uncertainty in wavelength. If values are scalars, they are assumed
        to be in angstroms.

    Returns
    -------
      fnu, fnu_err :: arrays of floats
        The flux in fnu and its uncertainty.
    """
    #
    # Check inputs
    #
    if isinstance(wavelength, u.Quantity):
        wavelength = wavelength.to(u.AA).value
    if isinstance(wavelength_err, u.Quantity):
        wavelength_err = wavelength_err.to(u.AA).value
    #
    # Convert fnu (erg/cm^2/s/Hz) to flam (erg/cm^2/s/A)
    #
    fnu = const.FLAM_TO_FNU * flam * (wavelength * wavelength)
    # REVIEW: check error propagation is correct
    fnu_err = fnu * np.sqrt(
        (flam_err / flam) ** 2 + 4 * (wavelength_err / wavelength) ** 2
    )  # derived & simplified from partial derivative error propagation
    return fnu, fnu_err


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
    the monochromatic flux of the source. This can also be used to convert counts to
    magnitudes if the photometric zero point is calibrated to 1 electron/s.

    Parameters
    ----------
      flux, flux_err :: scalar or arrays
        The monochromatic flux and its uncertainty. The unit of the flux depends on the
        magnitude system. For example, zpt=-48.60 corresponds to the AB magnitude system
        and the flux will be in units of erg/s/cm^2/Hz. Likewise, zpt=-21.1 corresponds to
        the ST magnitude system and the flux will be in units of erg/s/cm^2/A.

      zpt :: scalar
        The zero point of the magnitude system.

      calc_abs :: bool
        If True, calculate the absolute magnitude; also requires the dist parameter to be
        provided. If False, calculate the relative magnitude.

      dist, dist_err :: scalar or `astropy.Quantity` or array
        The distance to the source and its uncertainty. If values are scalars, they are
        assumed to be in parsecs. dist and dist_err are ignored if calc_abs is False.

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
    This can also be used to convert magnitudes to counts if the photometric zero point is
    calibrated to 1 electron/s.

    Parameters
    ----------
      mag, mag_err :: scalar or arrays
        The magnitudes and their uncertainties

      zpt :: scalar
        The zero point of the magnitude system. For example, zpt=-48.60 corresponds to the
        AB magnitude system and the flux will be in units of erg/s/cm^2/Hz. Likewise,
        zpt=-21.1 corresponds to the ST magnitude system and the flux will be in units of
        erg/s/cm^2/A.

    Returns
    -------
      flux, flux_err :: scalar or arrays
        The monochromatic flux and its uncertainty. The unit of the flux depends on the
        magnitude system. For example, zpt=-48.60 corresponds to the AB magnitude system
        and the flux will be in units of erg/s/cm^2/Hz. Likewise, zpt=-21.1 corresponds to
        the ST magnitude system and the flux will be in units of erg/s/cm^2/A.
    """
    flux = 10 ** (-0.4 * (mag - zpt))
    flux_err = 0.4 * np.log(10) * abs(flux * mag_err)
    return flux, flux_err


def convert_electron_flux_mag(
    var1,
    var1_type,
    var2_type,
    var1_err=0.0,
    passband=None,
    wavelengths=None,
    wavelengths_err=0.0,
):
    """
    Convert between electron rates (electron/s), flux in "fnu" units (erg/s/cm^2/Hz), flux
    in "flam" units (erg/s/cm^2/angstrom), and magnitudes (AB mags).

    Parameters
    ----------
      var1 :: scalar or array
        The electron rate (electron/s), flux (in "fnu" or "flam" units), or AB magnitude.

      var1_type, var2_type :: "electron", "fnu", "flam", or "mag"
        The type of the first and second variable.

      var1_err :: scalar or array
        The uncertainty in var1.

      passband :: "uv", "u", or "g"
        The passband of the variables. Required if var1_type or var2_type is "electron",
        ignored otherwise.

      wavelengths, wavelengths_err :: scalar or array
        The wavelengths and their uncertainties corresponding to the data represented by
        var1. Required if var1_type or var2_type is "flam", ignored otherwise.
        INFO: is there a way to convert from, e.g., AB mag to flam without this?

    Returns
    -------
      var2, var2_err :: scalar or array
        The converted var1 values and their uncertainties.
    """
    #
    # Check inputs
    #
    valid_inputs = ["electron", "fnu", "flam", "mag"]
    if var1_type not in valid_inputs or var2_type not in valid_inputs:
        raise ValueError(
            "var1_type and var2_type must be 'electron', 'fnu', 'flam', or 'mag'"
        )
    if var1_type == var2_type:
        raise ValueError("var1_type and var2_type must be different")
    if var1_type == "electron" or var2_type == "electron":
        if passband is None:
            raise ValueError(
                "passband must be provided if var1_type or var2_type is 'electron'"
            )
        elif (not isinstance(passband, str)) or (passband not in params.PASSBANDS):
            raise ValueError("passband must be 'uv', 'u', or 'g'")
    if var1_type == "flam" or var2_type == "flam":
        if wavelengths is None:
            raise ValueError(
                "wavelengths must be provided if var1_type or var2_type is 'flam'"
            )
    if np.any(var1_err < 0):
        raise ValueError("var1_err must be non-negative")

    #
    # Convert var1 & var1_err to desired quantities: var2 & var2_err
    #

    if var1_type == "mag":
        if var2_type == "electron":
            # AB magnitude to electron/s
            var2, var2_err = mag_to_flux(
                var1, mag_err=var1_err, zpt=params.PHOT_ZPTS[passband]
            )  # electron/s
        else:  # var2_type == "fnu" or var2_type == "flam"
            # AB magnitude to fnu
            var2, var2_err = mag_to_flux(var1, mag_err=var1_err, zpt=-48.60)  # fnu
            if var2_type == "flam":
                # AB magnitude to flam (via magnitude -> fnu -> flam)
                var2, var2_err = fnu_to_flam(
                    var2, wavelengths, fnu_err=var2_err, wavelength_err=wavelengths_err
                )  # flam

    elif var1_type == "fnu":
        if var2_type == "flam":
            # Fnu to flam
            var2, var2_err = fnu_to_flam(
                var1, wavelengths, fnu_err=var1_err, wavelength_err=wavelengths_err
            )  # flam
        else:  # var2_type == "mag" or var2_type == "electron"
            # Fnu to AB magnitude
            var2, var2_err = flux_to_mag(
                var1, flux_err=var1_err, zpt=-48.60, calc_abs=False
            )  # AB mags
            if var2_type == "electron":
                # Fnu to electron/s (via fnu -> AB mag -> electrons/s)
                var2, var2_err = mag_to_flux(
                    var2, mag_err=var2_err, zpt=params.PHOT_ZPTS[passband]
                )  # electron/s

    elif var1_type == "electron":
        # Electron/s to AB magnitude
        var2, var2_err = flux_to_mag(
            var1, flux_err=var1_err, zpt=params.PHOT_ZPTS[passband], calc_abs=False
        )  # AB mags (per the definition of the PHOT_ZPTS quantities)
        if var2_type == "fnu" or "flam":
            # Electron/s to fnu (via electron/s -> AB mag -> fnu)
            var2, var2_err = mag_to_flux(var2, mag_err=var2_err, zpt=-48.60)  # fnu
            if var2_type == "flam":
                # Electron/s to flam (via electron/s -> AB mag -> fnu -> flam)
                var2, var2_err = fnu_to_flam(
                    var2, wavelengths, fnu_err=var2_err, wavelength_err=wavelengths_err
                )  # flam

    elif var1_type == "flam":
        # Flam to fnu
        var2, var2_err = flam_to_fnu(
            var1, wavelengths, flam_err=var1_err, wavelength_err=wavelengths_err
        )  # flam
        if var2_type == "mag" or var2_type == "electron":
            # Flam to AB magnitude (via flam -> fnu -> AB mag)
            var2, var2_err = flux_to_mag(
                var2, flux_err=var2_err, zpt=-48.60, calc_abs=False
            )  # AB mags
            if var2_type == "electron":
                # Flam to electron/s (via flam -> fnu -> AB mag -> electron/s)
                var2, var2_err = mag_to_flux(
                    var2, mag_err=var2_err, zpt=params.PHOT_ZPTS[passband]
                )  # electron/s

    else:
        raise ValueError(
            "Something went wrong. This statement should be impossible to reach..."
        )

    return var2, var2_err