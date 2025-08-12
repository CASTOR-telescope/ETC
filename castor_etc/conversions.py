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
Unit conversion utilities
=========================

`castor_etc.conversions` provides utility methods to convert between different units and quantities.

Common units (also see <https://pysynphot.readthedocs.io/en/latest/units.html>):
  - fnu :: erg/s/cm^2/Hz
  - flam :: erg/s/cm^2/A
  - photnu :: photons/s/cm^2/Hz
  - photlam :: photons/s/cm^2/A

"""

import astropy.units as u
import numpy as np
from scipy.integrate import simpson

from . import constants as const


def calc_photon_energy(
    wavelength=None, frequency=None, wavelength_err=0.0, frequency_err=0.0):
    """
    Calculates the energy of a photon in ergs given its wavelength or frequency. Useful
    for converting between erg and photon units.

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


def convert_freq_wavelength(data, to="wavelength", output_unit = u.AA):
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
        raise ValueError("to must be 'wavelength' or 'frequency'")
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
    <https://hea-www.harvard.edu/~pgreen/figs/Conversions.pdf> ("F_\nu" to "f_\\lambda").

    Parameters
    ----------
      fnu :: array
        The flux in fnu.

      wavelength :: array of scalars or `astropy.Quantity` array
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
    <https://hea-www.harvard.edu/~pgreen/figs/Conversions.pdf> ("F_\nu" to "F_\\lambda").

    Parameters
    ----------
      fnu :: array of scalar
        The flux in fnu.

      wavelength :: scalar, `astropy.Quantity`, scalar array, or`astropy.Quantity` array
        The corresponding wavelengths of the flux. If values are scalars, they are assumed
        to be in angstroms. If not an array, wavelength should be the pivot wavelength of
        the passband.

      fnu_err :: scalar or array of scalars
        The absolute uncertainty in fnu. If a scalar, the same uncertainty is applied to
        all fnu values.

      wavelength_err :: scalar or `astropy.Quantity`
                        or array of scalars or `astropy.Quantity` array
        The absolute uncertainty in wavelength. If scalar(s), wavelength_err is assumed to
        be in angstroms. If not an array, wavelength_err should be the uncertainty in the
        pivot wavelength of the passband.

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
    flam_err = flam * np.sqrt(
        (fnu_err / fnu) ** 2 + 4 * (wavelength_err / wavelength) ** 2
    )  # derived & simplified from partial derivative error propagation
    return flam, flam_err


def flam_to_fnu(flam, wavelength, flam_err=0.0, wavelength_err=0.0):
    """
    Converts from flam (erg/cm^2/s/A) to fnu (erg/cm^2/s/Hz). See
    <https://hea-www.harvard.edu/~pgreen/figs/Conversions.pdf> ("F_\\lambda" to "F_\nu").

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
    new_mag_err = np.sqrt(mag_err**2 + (5 / np.log(10) * dist_err / dist) ** 2)
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
        and the flux will be in units of erg/s/cm^2/Hz. Likewise, zpt=-21.10 corresponds
        to the ST magnitude system and the flux will be in units of erg/s/cm^2/A.

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
        zpt=-21.10 corresponds to the ST magnitude system and the flux will be in units of
        erg/s/cm^2/A.

    Returns
    -------
      flux, flux_err :: scalar or arrays
        The monochromatic flux and its uncertainty. The unit of the flux depends on the
        magnitude system. For example, zpt=-48.60 corresponds to the AB magnitude system
        and the flux will be in units of erg/s/cm^2/Hz. Likewise, zpt=-21.10 corresponds
        to the ST magnitude system and the flux will be in units of erg/s/cm^2/A.
    """
    flux = 10 ** (-0.4 * (mag - zpt))
    flux_err = 0.4 * np.log(10) * abs(flux * mag_err)
    return flux, flux_err


def convert_AB_ST_mag(mag, wavelength, to="ABmag"):
    """
    Converts AB magnitude to ST magnitude and vice versa. See the derivation (up to Eq.
    (7)) in Casagrande & VandenBerg (2014):
    <https://ui.adsabs.harvard.edu/abs/2014MNRAS.444..392C/abstract>

    Parameters
    ----------
      mag :: scalar or arrays
        The magnitude(s) to convert.

      wavelength :: scalar or `astropy.Quantity`
                    or array of scalars or `astropy.Quantity` array
        If a single value, this is the pivot wavelength of the passband. If an array,
        these should be the wavelengths at which the magnitude is measured. If scalar(s),
        wavelength is assumed to be in angstrom.

      to :: "ABmag" or "STmag"
        If "ABmag", convert input magnitudes and their uncertainties to AB magnitudes. If
        "STmag", convert input magnitudes and their uncertainties to ST magnitudes.

    Returns
    -------
      converted_mag :: scalar or arrays
        The input magnitude(s) converted into the desired magnitude system.
    """
    #
    # Check inputs
    #
    if isinstance(wavelength, u.Quantity):
        try:
            wavelength = wavelength.to(u.AA).value
        except Exception:
            raise ValueError("`wavelength` must have the proper units (e.g., angstrom).")
    if np.any(wavelength <= 0):
        raise ValueError("All wavelengths must be positive")
    if to not in ["ABmag", "STmag"]:
        raise ValueError("`to` must be either 'ABmag' or 'STmag'")
    #
    # Do conversion
    #
    offset = 21.10 - 48.60 + 2.5 * np.log10(const.LIGHTSPEED.to(u.AA / u.s).value)
    if to == "ABmag":
        converted_mag = mag - 5 * np.log10(wavelength) + offset
    else:
        converted_mag = mag + 5 * np.log10(wavelength) - offset
    return converted_mag


def flam_to_AB_mag(wavelengths, flam, response):
    """
    Convert a spectrum in flam (erg/s/cm^2/A) to an AB magnitude. Follows Eq. (2) of
    Bessell & Murphy (2012)
    <https://ui.adsabs.harvard.edu/abs/2012PASP..124..140B/abstract>.

    Parameters
    ----------
      wavelengths :: array of floats or `astropy.Quantity` array
        The wavelengths over which to integrate the spectrum. If an array of floats, it is
        assumed to be in units of angstrom.

      flam :: array of floats
        The spectral flux density at the given wavelengths, in units of erg/s/cm^2/A.

      response :: array of floats
        The system response function at the given wavelengths. This represents the
        throughput of the whole optical path and should include the detector quantum
        efficiencies (i.e., converts from photons to electrons).

    Returns
    -------
      ab_mag :: float
        The AB magnitude of the spectrum.
    """
    #
    # Check inputs
    #
    if isinstance(wavelengths, u.Quantity):
        wavelengths = wavelengths.to(u.AA).value
    if np.any(~np.isfinite(wavelengths) | ~np.isfinite(flam) | ~np.isfinite(response)):
        raise ValueError("`wavelengths`, `flam`, and `response` must all be finite")
    #
    # Perform calculation (Eq. (2) of Bessell & Murphy (2012))
    #
    numer = simpson(y=flam * response * wavelengths, x=wavelengths)
    denom = simpson(y=response / wavelengths, x=wavelengths)
    return (
        -2.5 * np.log10(numer / (const.LIGHTSPEED.to(u.AA / u.s).value * denom)) - 48.60
    )  # AB magnitude


def convert_electron_flux_mag(
    var1,
    var1_type,
    var2_type,
    var1_err=0.0,
    phot_zpt=None,
    wavelengths=None,
    wavelengths_err=0.0,
):
    """
    Convert between electron rates (electron/s), flux in "fnu" units (erg/s/cm^2/Hz), flux
    in "flam" units (erg/s/cm^2/angstrom), and magnitudes (AB mags). To convert to/from
    electron/s, this function uses photometric zero points, which implicitly assumes the
    inputs are over the entire passband. For example, do NOT use this function to convert
    geocoronal emission line flux to electron/s.

    Parameters
    ----------
      var1 :: scalar or array
        The electron rate (electron/s), flux (in "fnu" or "flam" units), or AB magnitude.

      var1_type, var2_type :: "electron", "fnu", "flam", or "mag"
        The type of the first and second variable.

      var1_err :: scalar or array
        The uncertainty in var1.

      passband :: "uv", "u", or "g"
        The passband of the variables. Required if var1_type or var2_type is "electron"
        and phot_zpt is None, ignored otherwise.

      phot_zpt :: scalar
        The photometric zero point of the passband. This should be the AB magnitude
        required to produce 1 electron/s in the passband. Required if var1_type or
        var2_type is "electron" and passband is None, ignored otherwise. Will override any
        default zero points.

      wavelengths, wavelengths_err :: scalar or array or `astropy.Quantity` or
                                      `astropy.Quantity` array
        The wavelengths and their uncertainties corresponding to the data represented by
        var1. Required if var1_type or var2_type is "flam", ignored otherwise.
        Alternatively, if passband is specified and wavelengths is None, then the pivot
        wavelength of the passband will be used (the pivot wavelength is typically used if
        var1 is a quantity integrated over the whole passband). If scalar value(s),
        wavelengths and wavelengths_err are assumed to be in angstroms.

    Returns
    -------
      var2, var2_err :: scalar or array of scalars
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
        if phot_zpt is None:
            raise ValueError(
                "phot_zpts must be provided if var1_type or var2_type is 'electron'"
            )
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
                var1, mag_err=var1_err, zpt=phot_zpt
            )  # electron/s
        else:  # var2_type == "fnu" or var2_type == "flam"
            # AB magnitude to fnu
            var2, var2_err = mag_to_flux(var1, mag_err=var1_err, zpt=-48.60)  # fnu
            if var2_type == "flam":
                # AB magnitude to flam (via magnitude -> fnu -> flam)
                var2, var2_err = fnu_to_flam(
                    var2, wavelengths, fnu_err=var2_err, wavelength_err=wavelengths_err
                )  # flam
    #
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
                    var2, mag_err=var2_err, zpt=phot_zpt
                )  # electron/s
    #
    elif var1_type == "electron":
        # Electron/s to AB magnitude
        var2, var2_err = flux_to_mag(
            var1, flux_err=var1_err, zpt=phot_zpt, calc_abs=False
        )  # AB mags (per the definition of the PHOT_ZPTS quantities)
        if var2_type == "fnu" or "flam":
            # Electron/s to fnu (via electron/s -> AB mag -> fnu)
            var2, var2_err = mag_to_flux(var2, mag_err=var2_err, zpt=-48.60)  # fnu
            if var2_type == "flam":
                # Electron/s to flam (via electron/s -> AB mag -> fnu -> flam)
                var2, var2_err = fnu_to_flam(
                    var2, wavelengths, fnu_err=var2_err, wavelength_err=wavelengths_err
                )  # flam
    #
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
                    var2, mag_err=var2_err, zpt=phot_zpt
                )  # electron/s
    #
    else:
        raise ValueError(
            "Something went wrong. This statement should be impossible to reach..."
        )
    #
    return var2, var2_err
