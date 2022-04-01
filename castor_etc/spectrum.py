"""
Generate and handle spectral data and normalizations.

Includes:
  - redshifting wavelengths
  - blackbody radiation, power-law spectrum generation
  - user-input spectrum
  - Gaussian and Lorentzian emission/absorption lines
  - generic spiral and elliptical galaxy spectra
  - stellar spectra from the Pickles catalog
  - normalization functions:
    - normalize a blackbody spectrum to a star of given radius and distance
    - normalize a spectrum to some average value or AB magnitude, either within a passband
      or over the whole spectrum
    - normalize a spectrum to a given bolometric luminosity and distance
    - calculate the average value of a spectrum (erg/s/cm^2/A or AB mag) either within a
      passband or over the whole spectrum

---

        GNU General Public License v3 (GNU GPLv3)

(c) 2022.                            (c) 2022.
Government of Canada                 Gouvernement du Canada
National Research Council            Conseil national de recherches
Ottawa, Canada, K1A 0R6              Ottawa, Canada, K1A 0R6
All rights reserved                  Tous droits réservés

NRC disclaims any warranties,        Le CNRC dénie toute garantie
expressed, implied, or               énoncée, implicite ou légale,
statutory, of any kind with          de quelque nature que ce
respect to the software,             soit, concernant le logiciel,
including without limitation         y compris sans restriction
any warranty of merchantability      toute garantie de valeur
or fitness for a particular          marchande ou de pertinence
purpose. NRC shall not be            pour un usage particulier.
liable in any event for any          Le CNRC ne pourra en aucun cas
damages, whether direct or           être tenu responsable de tout
indirect, special or general,        dommage, direct ou indirect,
consequential or incidental,         particulier ou général,
arising from the use of the          accessoire ou fortuit, résultant
software. Neither the name           de l'utilisation du logiciel. Ni
of the National Research             le nom du Conseil National de
Council of Canada nor the            Recherches du Canada ni les noms
names of its contributors may        de ses  participants ne peuvent
be used to endorse or promote        être utilisés pour approuver ou
products derived from this           promouvoir les produits dérivés
software without specific prior      de ce logiciel sans autorisation
written permission.                  préalable et particulière
                                     par écrit.

This file is part of the             Ce fichier fait partie du projet
FORECASTOR ETC project.              FORECASTOR ETC.

FORECASTOR ETC is free software:     FORECASTOR ETC est un logiciel
you can redistribute it and/or       libre ; vous pouvez le redistribuer
modify it under the terms of         ou le modifier suivant les termes de
the GNU General Public               la "GNU General Public
License as published by the          License" telle que publiée
Free Software Foundation,            par la Free Software Foundation :
either version 3 of the              soit la version 3 de cette
License, or (at your option)         licence, soit (à votre gré)
any later version.                   toute version ultérieure.

FORECASTOR ETC is distributed        FORECASTOR ETC est distribué
in the hope that it will be          dans l'espoir qu'il vous
useful, but WITHOUT ANY WARRANTY;    sera utile, mais SANS AUCUNE
without even the implied warranty    GARANTIE : sans même la garantie
of MERCHANTABILITY or FITNESS FOR    implicite de COMMERCIALISABILITÉ
A PARTICULAR PURPOSE. See the        ni d'ADÉQUATION À UN OBJECTIF
GNU General Public License for       PARTICULIER. Consultez la Licence
more details.                        Générale Publique GNU pour plus
                                     de détails.

You should have received             Vous devriez avoir reçu une
a copy of the GNU General            copie de la Licence Générale
Public License along with            Publique GNU avec FORECASTOR ETC ;
FORECASTOR ETC. If not, see          si ce n'est pas le cas, consultez :
<http://www.gnu.org/licenses/>.      <http://www.gnu.org/licenses/>.
"""

import warnings
from numbers import Number
from os.path import join

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.io import fits
from scipy.integrate import simpson
from scipy.interpolate import interp1d

from . import constants as const
from . import parameters as params
from .conversions import calc_photon_energy, convert_electron_flux_mag
from .filepaths import DATAPATH
from .telescope import Telescope


def redshift_wavelengths(wavelengths, redshift):
    """
    Apply redshift correction to wavelengths.

    Parameters
    ----------
      wavelengths :: int or float or `astropy.Quantity` or array
        Wavelengths to redshift.

      redshift :: int or float
        Redshift to apply.

    Returns
    -------
        red_wavelengths :: array of floats
        Redshifted wavelengths.
    """
    return wavelengths * (1 + redshift)


class SpectrumMixin:
    """
    Mixin for generating spectra. To be used with `Source` object. Do not use directly!
    """

    def _check_existing_spectrum(self, overwrite):
        """
        Check for existing spectrum. If self.wavelengths and/or self.spectrum is not None,
        raise an error if overwrite is False otherwise overwrite the existing spectrum.
        """
        if self.wavelengths is not None or self.spectrum is not None:
            if overwrite:
                print(
                    "INFO: Overwriting existing wavelengths/spectrum "
                    + "with new wavelengths/spectrum."
                )
            else:
                raise ValueError(
                    "wavelengths/spectrum already exists! "
                    + "Use overwrite=True to overwrite wavelengths/spectrum."
                )

    def spectrum_erg_to_photon(self):
        """
        Convert the spectrum that has units of ergs in the numerator to units of photons
        in the numerator.

        Attributes
        ----------
          spectrum :: array of floats
            Spectrum that has units of photons in the numerator.

        Returns
        -------
          None
        """
        if self.spectrum is None or self.wavelengths is None:
            raise ValueError(
                "Please generate a spectrum before converting ergs to photons."
            )
        self.spectrum /= calc_photon_energy(wavelength=self.wavelengths)[0]

    def redshift_wavelengths(self, redshift):
        """
        Apply redshift correction to wavelengths. Does not affect the y-axis values of the
        spectrum.

        Parameters
        ----------
          redshift :: int or float
            Redshift to apply.

        Attributes
        ----------
          wavelengths :: `astropy.Quantity` array
            The redshifted wavelengths of the spectrum, in angstroms.

        Returns
        -------
          None
        """
        if not isinstance(redshift, Number):
            raise TypeError("redshift must be an int or float")
        self.wavelengths *= 1 + redshift

    def generate_uniform(self, wavelengths, value):
        """
        Generate a uniform spectrum equal to a constant value in some arbitrary unit.

        Parameters
        ----------
          wavelengths :: array of scalars or `astropy.Quantity` array
            The wavelengths over which to generate the uniform spectrum. If an array of
            scalars, it should be in angstrom.

          value :: int or float
            The value of the uniform spectrum, in arbitrary units (e.g., erg/s/cm^2/A).

        Attributes
        ----------
          wavelengths :: `astropy.Quantity` array
            The wavelengths of the spectrum, in angstroms.

          spectrum :: array of floats
            Uniform spectrum in arbitrary units.

        Returns
        -------
          None
        """
        #
        # Check inputs
        #
        if not isinstance(wavelengths, u.Quantity):
            wavelengths = wavelengths * u.AA
        else:
            wavelengths = wavelengths.to(u.AA)
        if not isinstance(value, Number):
            raise ValueError("`value` must be an int or float")
        self.wavelengths = wavelengths
        self.spectrum = np.full(np.shape(wavelengths), value, dtype=float)

    def generate_bb(
        self,
        T,
        redshift=0.0,
        emissivity=1.0,
        wavelengths=None,
        limits=[0.1, 1.2] << u.um,
        resolution=1 << u.nm,
        radius=1,
        dist=1 << u.kpc,
        overwrite=False,
    ):
        """
        Generate a blackbody (BB) spectrum (in erg/s/cm^2/A) using Planck's radiation law.
        The spectral radiance of the BB (erg/s/cm^2/A/sr) is normalized to a star of given
        radius and distance.

        Parameters
        ----------
          T :: int or float or `astropy.Quantity`
            Intrinsic blackbody temperature (i.e., the temperature of the BB at
            redshift=0). If int or float, the unit is assumed to be kelvin.

          redshift :: int or float
            Redshift of the blackbody.

          emissivity :: int or float
            Emissivity of the blackbody. (Technically, emissivity is unity per the
            definition of a BB).

          wavelengths :: array of floats or `astropy.Quantity` array
            The wavelengths over which to calculate the spectrum. If an array of floats,
            the unit is assumed to be in angstroms. If wavelengths is not None, the limits
            and resolution parameters are ignored. Note that the final wavelengths
            attribute will be an `astropy.Quantity` array in units of angstroms regardless
            of this input's units.

          limits :: list of 2 scalars or list of 2 `astropy.Quantity`
            List containing the lower (0th index) and upper (1st index) bounds for the BB
            spectrum's restframe wavelengths, inclusive. Limits should be > 0. If list
            elements are int or float, they are assumed to be in angstroms. This parameter
            is ignored if wavelengths is provided.

          resolution :: int or float or `astropy.Quantity`
            The wavelength resolution of the returned spectrum. If a scalar, it is assumed
            to be in units of angstroms. This parameter is ignored if wavelengths is
            provided.

          radius :: float or `astropy.Quantity`
            The radius of the source. If a scalar, it is assumed to be in units of solar
            radii.

          dist :: float or `astropy.Quantity`
            The distance to the blackbody. If a scalar, it is assumed to be in units of
            kpc.

          overwrite :: bool
            If True, overwrite any existing wavelengths/spectrum. If False, raise an error
            if wavelengths or spectrum is not None.

        Attributes
        ----------
          wavelengths :: `astropy.Quantity` array
            The redshifted wavelengths of the spectrum, in angstroms.

          spectrum :: array of floats
            BB spectrum in units of flam (erg/s/cm^2/A).

        Returns
        -------
          None
        """
        #
        # Check inputs
        #
        self._check_existing_spectrum(overwrite)
        if isinstance(T, u.Quantity):
            T = T.to(u.K, equivalencies=u.temperature()).value
        if wavelengths is None:
            limits = list(limits)
            for i, lim in enumerate(limits):
                if isinstance(lim, u.Quantity):
                    limits[i] = lim.to(u.AA).value
            if isinstance(resolution, u.Quantity):
                resolution = resolution.to(u.AA).value
            wavelengths = np.arange(limits[0], limits[1] + 0.5 * resolution, resolution)
        elif isinstance(wavelengths, u.Quantity):
            wavelengths = wavelengths.to(u.AA).value
        #
        # Generate BB spectrum with redshift
        #
        # Factor in redshift & convert wavelengths from angstrom to cm
        wavelengths = redshift_wavelengths(wavelengths, redshift) * 1e-8  # cm
        # Planck's radiation law
        lightspeed = const.LIGHTSPEED.value  # cm/s
        prefactor = (2 * const.PLANCK_H.value * lightspeed * lightspeed) / (
            wavelengths ** 5
        )
        denom = np.expm1(
            (const.PLANCK_H.value * lightspeed) / (wavelengths * const.K_B.value * T)
        )
        spectrum = prefactor / denom  # erg/s/cm^2/cm/sr
        #
        # Incorporate emissivity and convert per cm to per angstrom
        #
        spectrum *= 1e-8 * emissivity  # erg/s/cm^2/A/sr
        #
        # Convert wavelengths back to angstroms
        #
        wavelengths *= 1e8  # angstrom
        #
        # Assign to `Source` object attributes. Spectrum is in erg/s/cm^2/A
        #
        self.wavelengths = wavelengths * u.AA
        self.spectrum = NormMixin.norm_to_star(spectrum, radius=radius, dist=dist)  # flam

    def generate_power_law(self, ref_wavelength, wavelengths, exponent, overwrite=False):
        """
        Generate a spectrum with a shape following a power-law in some arbitrary unit. The
        flux is defined so that it is equal to 1 at the reference wavelength.

        The spectrum is calculated using the following formula:
        ```math
                    spectrum = (wavelengths / ref_wavelength)^exponent
        ```
        where each variable is as defined in the Parameters documentation below.

        Parameters
        ----------
          ref_wavelength :: scalar or `astropy.Quantity`
            The reference wavelength for the power-law. The spectrum at this wavelength
            will have a flux of 1. If a scalar, it should be in angstrom.

          wavelengths :: array of scalars or `astropy.Quantity` array
            The wavelengths over which to calculate the power-law spectrum. If an array of
            scalars, it should be in angstrom.

          exponent :: int or float
            The exponent for the power-law.

          overwrite :: bool
            If True, overwrite any existing wavelengths/spectrum. If False, raise an error
            if wavelengths or spectrum is not None.

        Attributes
        ----------
          wavelengths :: `astropy.Quantity` array
            The wavelengths of the spectrum, in angstroms.

          spectrum :: array of floats
            Spectrum following the power-law defined above, in arbitary units.

        Returns
        -------
          None
        """
        #
        # Check inputs
        #
        self._check_existing_spectrum(overwrite)
        if np.size(ref_wavelength) != 1:
            raise ValueError(
                "ref_wavelength must be a single scalar or `astropy.Quantity`."
            )
        if isinstance(wavelengths, u.Quantity):
            wavelengths = wavelengths.to(u.AA).value
        if isinstance(ref_wavelength, u.Quantity):
            ref_wavelength = ref_wavelength.to(u.AA).value
        #
        # Power-law
        #
        spectrum = (wavelengths / ref_wavelength) ** exponent
        self.wavelengths = wavelengths * u.AA
        self.spectrum = spectrum

    @staticmethod
    def _generate_gaussian(
        wavelengths,
        spectrum,
        center,
        fwhm,
        peak=None,
        tot_flux=None,
        add=True,
        abs_peak=True,
    ):
        """
        Add/subtract a Gaussian spectrum to/from an existing spectrum. This is useful for
        representing emission lines (i.e., by adding a Gaussian source) or absorption
        lines (i.e., by subtracting a Gaussian source). Note that the minimum/maximum
        wavelengths of the source spectrum will not change.

        The Gaussian spectrum can be represented by the following formulae (from
        <https://pysynphot.readthedocs.io/en/latest/spectrum.html#gaussian-emission>):
        ```math
                    gaussian = peak / exp[(wavelengths - center)^2 / (2 * sigma^2)]
        ```
        and
        ```math
                    sigma = fwhm / [2 * sqrt(2 * ln2)]
        ```
        and
        ```math
                    peak = tot_flux / sqrt(2 * pi * sigma^2) <-- see Gaussian integral
        ```
        where:
          - gaussian is the Gaussian spectrum's flux in some arbitrary unit
          - peak is the flux at the center of the Gaussian (i.e., the central wavelength)
          - center is the central wavelength of the Gaussian
          - wavelengths is the array of wavelengths over which to calculate the spectrum
          - fwhm is the full-width at half-maximum of the Gaussian
          - tot_flux is the total flux of the Gaussian under the curve

        Parameters
        ----------
          wavelengths :: array of floats or `astropy.Quantity` array
            The wavelengths over which to calculate the Gaussian spectrum.

          spectrum :: array of floats
            The spectrum to/from which to add/subtract the Gaussian spectrum.

          center :: scalar or `astropy.Quantity`
            The central wavelength of the Gaussian.

          fwhm :: scalar or `astropy.Quantity`
            The full-width at half-maximum of the Gaussian.

          peak :: int or float
            The peak flux of the Gaussian (i.e., the flux at the center wavelength). This
            determines the unit of the returned spectrum. Exactly one of peak or tot_flux
            must be specified.

          tot_flux :: int or float
            The total flux under the curve. This implicitly determines the unit of the
            returned spectrum. Exactly one of peak or tot_flux must be specified.

          add :: bool
            If True, add the Gaussian spectrum to the existing spectrum. If False,
            subtract the Gaussian from the existing spectrum.

          abs_peak :: bool
            If True, scale spectrum so that the peak of the emission or dip of the
            absorption line is at the given value. Otherwise, just naively add/subtract
            the given Gaussian peak.

        Returns
        -------
          sorted_wavelengths :: array of floats or `astropy.Quantity` array
            The wavelengths of the new spectrum. The shape of this array will be different
            from the input `wavelengths` array.

          sorted_spectrum :: array of floats
            The spectrum with the Gaussian added/subtracted. The shape of this array will
            be different from the input `spectrum` array.
        """

        def _gaussian(_peak, _wavelengths, _center, _sigma):
            _num_sigma = (_wavelengths - _center) / _sigma
            return _peak / np.exp(0.5 * _num_sigma * _num_sigma)

        #
        # Check inputs
        #
        if np.size(center) != 1:
            raise ValueError("center must be a single scalar or `astropy.Quantity`.")
        if np.size(fwhm) != 1:
            raise ValueError("fwhm must be a single scalar or `astropy.Quantity`.")
        if (peak is None and tot_flux is None) or (
            peak is not None and tot_flux is not None
        ):
            raise ValueError("Exactly one of peak or tot_flux must be specified.")
        if peak is not None and (
            np.size(peak) != 1 or not isinstance(peak, Number) or peak <= 0
        ):
            raise ValueError("peak must be a single int or float >= 0.")
        elif tot_flux is not None and (
            np.size(tot_flux) != 1 or not isinstance(tot_flux, Number) or tot_flux <= 0
        ):
            raise ValueError("tot_flux must be a single int or float >= 0.")
        # Convert lengths to angstrom
        if isinstance(wavelengths, u.Quantity):
            wavelengths_unit = wavelengths.unit
            wavelengths = wavelengths.to(u.AA).value
        else:
            wavelengths_unit = None
        if isinstance(center, u.Quantity):
            center = center.to(u.AA).value
        if isinstance(fwhm, u.Quantity):
            fwhm = fwhm.to(u.AA).value
        if fwhm <= 0:
            raise ValueError("fwhm must be >= 0.")
        #
        # Gaussian spectrum
        #
        sigma = fwhm / (2 * np.sqrt(2 * np.log(2)))
        if peak is None:
            peak = tot_flux / (np.sqrt(2 * np.pi) * sigma)
        spectrum_interp = interp1d(
            wavelengths, spectrum, kind="linear", bounds_error=False, fill_value=np.nan
        )
        if add:
            if abs_peak:
                # Ensure final peak is actually at desired value
                peak -= spectrum_interp(center)
                if peak < 0:
                    raise ValueError("peak of emission line below continuum.")
        else:
            if abs_peak:
                # Ensure final dip is actually at desired value
                center_val = spectrum_interp(center)
                if peak > center_val:
                    raise ValueError("dip of absorption line above continuum.")
                peak = center_val - peak
        # Ensure Gaussian is well sampled by evaluating at the wavelengths within +/- 5
        # sigma of the center. This also prevents overflow errors caused by calculations
        # at wavelengths too far from the center.
        gauss_wavelengths = center + np.arange(-5, 5.05, 0.1) * sigma
        sorted_wavelengths = np.unique(np.concatenate((wavelengths, gauss_wavelengths)))
        sorted_wavelengths = sorted_wavelengths[
            (sorted_wavelengths > 0) & (sorted_wavelengths <= wavelengths[-1])
        ]
        sorted_spectrum = spectrum_interp(sorted_wavelengths)
        in_range = (sorted_wavelengths >= gauss_wavelengths[0]) & (
            sorted_wavelengths <= gauss_wavelengths[-1]
        )
        if add:
            sorted_spectrum[in_range] += _gaussian(
                peak, sorted_wavelengths[in_range], center, sigma
            )
        else:
            sorted_spectrum[in_range] -= _gaussian(
                peak, sorted_wavelengths[in_range], center, sigma
            )
        is_good = np.isfinite(sorted_wavelengths) & np.isfinite(sorted_spectrum)
        sorted_wavelengths = sorted_wavelengths[is_good]
        sorted_spectrum = sorted_spectrum[is_good]
        is_negative_spectrum = sorted_spectrum < 0
        if np.any(is_negative_spectrum):
            sorted_spectrum[is_negative_spectrum] = 0.0
            warnings.warn("Setting negative flux in spectrum to zero.", RuntimeWarning)
        if wavelengths_unit is not None:
            sorted_wavelengths <<= wavelengths_unit  # convert to `astropy.Quantity` array
        return sorted_wavelengths, sorted_spectrum

    @staticmethod
    def _generate_lorentzian(
        wavelengths,
        spectrum,
        center,
        fwhm,
        peak=None,
        tot_flux=None,
        add=True,
        abs_peak=True,
    ):
        """
        Add/subtract a Lorentzian spectrum to/from an existing spectrum. This is useful
        for representing emission lines (i.e., by adding a Lorentzian source) or
        absorption lines (i.e., by subtracting a Lorentzian source). Note that the
        minimum/maximum wavelengths of the source spectrum will not change.

        The Lorentzian spectrum can be represented by the following formulae:
        ```math
                    lorentzian = peak / (1 + num_half_widths^2)
        ```
        and
        ```math
                    num_half_widths = (wavelengths - center) / probable_error
        ```
        and
        ```math
                    peak = tot_flux / (pi * probable_error) <-- see Cauchy distribution
        ```
        and
        ```math
                    probable_error = fwhm / 2
        ```
        where:
          - lorentzian is the Lorentzian spectrum's flux in some arbitrary unit
          - peak is the flux at the center (i.e., central wavelength) of the Lorentzian
          - center is the central wavelength of the Lorentzian
          - wavelengths is the array of wavelengths over which to calculate the spectrum
          - fwhm is the full-width at half-maximum of the Lorentzian
          - tot_flux is the total flux of the Lorentzian under the curve

        Parameters
        ----------
          wavelengths :: array of floats or `astropy.Quantity` array
            The wavelengths over which to calculate the Lorentzian spectrum.

          spectrum :: array of floats
            The spectrum to/from which to add/subtract the Lorentzian spectrum.

          center :: scalar or `astropy.Quantity`
            The central wavelength of the Lorentzian.

          fwhm :: scalar or `astropy.Quantity`
            The full-width at half-maximum of the Lorentzian.

          peak :: int or float
            The peak flux of the Lorentzian (i.e., the flux at the center wavelength). This
            determines the unit of the returned spectrum. Exactly one of peak or tot_flux
            must be specified.

          tot_flux :: int or float
            The total flux under the curve. This implicitly determines the unit of the
            returned spectrum. Exactly one of peak or tot_flux must be specified.

          add :: bool
            If True, add the Lorentzian spectrum to the existing spectrum. If False,
            subtract the Lorentzian from the existing spectrum.

          abs_peak :: bool
            If True, scale spectrum so that the peak of the emission or dip of the
            absorption line is at the given value. Otherwise, just naively add/subtract
            the given Lorentzian peak.

        Returns
        -------
          sorted_wavelengths :: array of floats or `astropy.Quantity` array
            The wavelengths of the new spectrum. The shape of this array will be different
            from the input `wavelengths` array.

          sorted_spectrum :: array of floats
            The spectrum with the Lorentzian added/subtracted. The shape of this array will
            be different from the input `spectrum` array.
        """

        def _lorentzian(_peak, _wavelengths, _center, _probable_error):
            _num_half_widths = (_wavelengths - _center) / _probable_error
            return _peak / (1 + _num_half_widths * _num_half_widths)

        #
        # Check inputs
        #
        if np.size(center) != 1:
            raise ValueError("center must be a single scalar or `astropy.Quantity`.")
        if np.size(fwhm) != 1:
            raise ValueError("fwhm must be a single scalar or `astropy.Quantity`.")
        if (peak is None and tot_flux is None) or (
            peak is not None and tot_flux is not None
        ):
            raise ValueError("Exactly one of peak or tot_flux must be specified.")
        if peak is not None and (
            np.size(peak) != 1 or not isinstance(peak, Number) or peak <= 0
        ):
            raise ValueError("peak must be a single int or float >= 0.")
        elif tot_flux is not None and (
            np.size(tot_flux) != 1 or not isinstance(tot_flux, Number) or tot_flux <= 0
        ):
            raise ValueError("tot_flux must be a single int or float >= 0.")
        # Convert lengths to angstrom
        if isinstance(wavelengths, u.Quantity):
            wavelengths_unit = wavelengths.unit
            wavelengths = wavelengths.to(u.AA).value
        else:
            wavelengths_unit = None
        if isinstance(center, u.Quantity):
            center = center.to(u.AA).value
        if isinstance(fwhm, u.Quantity):
            fwhm = fwhm.to(u.AA).value
        if fwhm <= 0:
            raise ValueError("fwhm must be >= 0.")
        #
        # Lorentzian spectrum
        #
        probable_error = 0.5 * fwhm
        if peak is None:
            peak = tot_flux / (np.pi * probable_error)
        spectrum_interp = interp1d(
            wavelengths, spectrum, kind="linear", bounds_error=False, fill_value=np.nan
        )
        if add:
            if abs_peak:
                # Ensure final peak is actually at desired value
                peak -= spectrum_interp(center)
                if peak < 0:
                    raise ValueError("peak of emission line below continuum.")
        else:
            if abs_peak:
                # Ensure final dip is actually at desired value
                center_val = spectrum_interp(center)
                if peak > center_val:
                    raise ValueError("dip of absorption line above continuum.")
                peak = center_val - peak
        # Ensure Lorentzian is well sampled by evaluating at the wavelengths within +/- 80
        # units of probable error from the center. This also prevents overflow errors
        # caused by calculations at wavelengths too far from the center.
        gauss_wavelengths = center + np.arange(-80, 80.25, 0.5) * probable_error
        sorted_wavelengths = np.unique(np.concatenate((wavelengths, gauss_wavelengths)))
        sorted_wavelengths = sorted_wavelengths[
            (sorted_wavelengths > 0) & (sorted_wavelengths <= wavelengths[-1])
        ]
        sorted_spectrum = spectrum_interp(sorted_wavelengths)
        in_range = (sorted_wavelengths >= gauss_wavelengths[0]) & (
            sorted_wavelengths <= gauss_wavelengths[-1]
        )
        if add:
            sorted_spectrum[in_range] += _lorentzian(
                peak, sorted_wavelengths[in_range], center, probable_error
            )
        else:
            sorted_spectrum[in_range] -= _lorentzian(
                peak, sorted_wavelengths[in_range], center, probable_error
            )
        is_good = np.isfinite(sorted_wavelengths) & np.isfinite(sorted_spectrum)
        sorted_wavelengths = sorted_wavelengths[is_good]
        sorted_spectrum = sorted_spectrum[is_good]
        is_negative_spectrum = sorted_spectrum < 0
        if np.any(is_negative_spectrum):
            sorted_spectrum[is_negative_spectrum] = 0.0
            warnings.warn("Setting negative flux in spectrum to zero.", RuntimeWarning)
        if wavelengths_unit is not None:
            sorted_wavelengths <<= wavelengths_unit  # convert to `astropy.Quantity` array
        return sorted_wavelengths, sorted_spectrum

    def add_emission_line(
        self, center, fwhm, peak=None, tot_flux=None, shape="gaussian", abs_peak=True
    ):
        """
        TODO: docstring

        Note that the minimum/maximum wavelengths of the source spectrum will not change.
        """
        if self.wavelengths is None or self.spectrum is None:
            raise ValueError("Please generate or load a spectrum first")
        if shape == "gaussian":
            spectrum_func = SpectrumMixin._generate_gaussian
        elif shape == "lorentzian":
            spectrum_func = SpectrumMixin._generate_lorentzian
        else:
            raise ValueError("Emission line shape must be 'gaussian' or 'lorentzian'")
        self.wavelengths, self.spectrum = spectrum_func(
            self.wavelengths,
            self.spectrum,
            center,
            fwhm,
            peak=peak,
            tot_flux=tot_flux,
            add=True,
            abs_peak=abs_peak,
        )

    def add_absorption_line(
        self, center, fwhm, dip=None, tot_flux=None, shape="gaussian", abs_dip=True
    ):
        """
        TODO: docstring

        Note that the minimum/maximum wavelengths of the source spectrum will not change.
        """
        if self.wavelengths is None or self.spectrum is None:
            raise ValueError("Please generate or load a spectrum first")
        if shape == "gaussian":
            spectrum_func = SpectrumMixin._generate_gaussian
        elif shape == "lorentzian":
            spectrum_func = SpectrumMixin._generate_lorentzian
        else:
            raise ValueError("Absorption line shape must be 'gaussian' or 'lorentzian'")
        self.wavelengths, self.spectrum = spectrum_func(
            self.wavelengths,
            self.spectrum,
            center,
            fwhm,
            peak=dip,
            tot_flux=tot_flux,
            add=False,
            abs_peak=abs_dip,
        )

    def use_custom_spectrum(self, filepath, wavelength_unit=u.AA, overwrite=False):
        """
        Use custom spectrum from an ASCII or FITS file.

        Parameters
        ----------
          filepath :: str
            The absolute path to the file containing the spectrum.
            If the file is in ASCII format, the first column should contain the
            wavelengths in `wavelength_units` and the second column containing the
            spectrum in flam (erg/s/cm^2/A); the columns should be separated by a constant
            number of spaces. Lines starting with a hash (#) will be ignored. The file
            extension must not be .fit or .fits.
            If the file is in FITS format, the first field (index 0) should contain the
            wavelengths in `wavelength_units` and the second field (index 1) should
            contain the spectrum in flam (erg/s/cm^2/A). The file extension must be .fit
            or .fits.

          wavelength_unit :: `astropy.Quantity` length unit
            The unit of the wavelengths in the file (e.g., u.AA for angstrom, u.nm for
            nanometer, u.um for micrometer, etc.)

          overwrite :: bool
            If True, overwrite any existing wavelengths/spectrum. If False, raise an error
            if wavelengths or spectrum is not None.

        Attributes
        ----------
          wavelengths :: `astropy.Quantity` array
            The wavelengths of the spectrum, in angstroms.

          spectrum :: array of floats
            Source spectrum in units of flam (erg/s/cm^2/A).

        Returns
        -------
          None
        """
        #
        # Check inputs
        #
        self._check_existing_spectrum(overwrite)
        if not isinstance(filepath, str):
            raise TypeError("filepath must be a string.")
        try:
            _ = wavelength_unit.to(u.AA)
        except Exception:
            raise TypeError("wavelength_units must be an `astropy.Quantity` length unit.")
        #
        # Load spectrum
        #
        file_ext = filepath.split(".")[-1].lower()
        try:
            if file_ext == "fit" or file_ext == "fits":
                data = fits.getdata(filepath)
                self.wavelengths = (data.field(0) * wavelength_unit).to(u.AA)
                self.spectrum = data.field(1)
            else:
                data = pd.read_csv(
                    filepath,
                    sep=" +",
                    header=None,
                    comment="#",
                    engine="python",
                )  # sep=" +" is Python regex to match a variable number of spaces
                self.wavelengths = (data[0].values * wavelength_unit).to(u.AA)
                self.spectrum = data[1].values
        except Exception:
            raise RuntimeError(
                "Could not read spectrum from file. File must be in ASCII or FITS format "
                + "and adhere to the guidelines specified in the docstring."
            )

    def use_galaxy_spectrum(self, gal_type, overwrite=False):
        """
        Use one of the predefined galaxy spectra. These non-uniformly sampled spectra are
        from Fioc & Rocca-Volmerange (1997)
        <https://ui.adsabs.harvard.edu/abs/1997A%26A...326..950F/abstract>. In particular,
        the data files were downloaded from the Gemini Observatory Control Software (OCS)
        GitHub repository:
        <https://github.com/gemini-hlsw/ocs/blob/develop/bundle/edu.gemini.itc/src/main/resources/sed/non_stellar/elliptical-galaxy.nm>
        and
        <https://github.com/gemini-hlsw/ocs/blob/develop/bundle/edu.gemini.itc/src/main/resources/sed/non_stellar/spiral-galaxy.nm>.

        Parameters
        ----------
          gal_type :: "elliptical" or "spiral"
            The galaxy morphology. The elliptical galaxy (class T=-5, -4) and spiral
            galaxy (type Sc, class T=5) spectra both run from 22-10000 nm.

          overwrite :: bool
            If True, overwrite any existing wavelengths/spectrum. If False, raise an error
            if wavelengths or spectrum is not None.

        Attributes
        ----------
          wavelengths :: `astropy.Quantity` array
            The wavelengths of the spectrum, in angstroms.

          spectrum :: array of floats
            Galaxy spectrum in arbitrary units, normalized such that it is equal to 1 at
            550 nm.

        Returns
        -------
          None
        """
        #
        # Check inputs
        #
        self._check_existing_spectrum(overwrite)
        if gal_type == "elliptical" or gal_type == "spiral":
            filepath = join(DATAPATH, "galaxy_spectra", f"{gal_type}_galaxy.txt")
        else:
            raise ValueError("Galaxy type must be 'elliptical' or 'spiral'")
        data = pd.read_csv(
            filepath,
            sep=" +",
            header=None,
            comment="#",
            engine="python",
        )  # sep=" +" is Python regex to match a variable number of spaces
        self.wavelengths = (data[0].values * u.nm).to(u.AA)
        self.spectrum = data[1].values

    def use_pickles_spectrum(self, spectral_class, overwrite=False):
        """
        Use a spectrum from the Pickles catalog
        (<https://ui.adsabs.harvard.edu/abs/1998PASP..110..863P/abstract>) containing
        spectra for numerous stellar spectral classes.

        A table containing valid `spectral_class` inputs is at the end of this docstring.

        Parameters
        ----------
          spectral_class :: str from the "Table of valid `spectral_class` inputs" below
            The spectral type of the star.

          overwrite :: bool
            If True, overwrite any existing wavelengths/spectrum. If False, raise an error
            if wavelengths or spectrum is not None.

        Attributes
        ----------
          wavelengths :: `astropy.Quantity` array
            The wavelengths of the spectrum, in angstroms.

          spectrum :: array of floats
            Stellar spectrum in arbitrary units, normalized such that it is equal to 1 at
            5556 angstrom.

        Returns
        -------
          None

        Table of valid `spectral_class` inputs
        -------------------------------------
        ```text
          spectral_class      Description (& wavelength range)
          --------------      --------------------------------
          "a0i"               A0 I     (1150-10620 A)
          "a0iii"             A0 III   (1150-10620 A)
          "a0iv"              A0 IV    (1150-10620 A)
          "a0v"               A0 V     (1150-10620 A)
          "a2i"               A2 I     (1150-10620 A)
          "a2v"               A2 V     (1150-10620 A)
          "a3iii"             A3 III   (1150-10620 A)
          "a3v"               A3 V     (1150-10620 A)
          "a47iv"             A4-7 IV  (1150-10620 A)
          "a5iii"             A5 III   (1150-10620 A)
          "a5v"               A5 V     (1150-10620 A)
          "a7iii"             A7 III   (1150-10620 A)
          "a7v"               A7 V     (1150-10620 A)
          "b0i"               B0 I     (1150-10620 A)
          "b0v"               B0 V     (1150-10620 A)
          "b12iii"            B1-2 III (1150-10620 A)
          "b1i"               B1 I     (1150-10620 A)
          "b1v"               B1 V     (1150-10620 A)
          "b2ii"              B2 II    (1150-10620 A)
          "b2iv"              B2 IV    (1150-10620 A)
          "b3i"               B3 I     (1150-10620 A)
          "b3iii"             B3 III   (1150-10620 A)
          "b3v"               B3 V     (1150-10620 A)
          "b57v"              B5-7 V   (1150-10620 A)
          "b5i"               B5 I     (1150-10620 A)
          "b5ii"              B5 II    (1150-10620 A)
          "b5iii"             B5 III   (1150-10620 A)
          "b6iv"              B6 IV    (1150-10620 A)
          "b8i"               B8 I     (1150-10620 A)
          "b8v"               B8 V     (1150-10620 A)
          "b9iii"             B9 III   (1150-10620 A)
          "b9v"               B9 V     (1150-10620 A)
          "f02iv"             F0-2 IV  (1150-10620 A)
          "f0i"               F0 I     (1150-10620 A)
          "f0ii"              F0 II    (1150-10620 A)
          "f0iii"             F0 III   (1150-10620 A)
          "f0v"               F0 V     (1150-10620 A)
          "f2ii"              F2 II    (1150-10620 A)
          "f2iii"             F2 III   (1150-10620 A)
          "f2v"               F2 V     (1150-10620 A)
          "f5i"               F5 I     (1150-10620 A)
          "f5iii"             F5 III   (1150-10620 A)
          "f5iv"              F5 IV    (1150-10620 A)
          "f5v"               F5 V     (1150-10620 A)
          "f6v"               F6 V     (1150-10620 A)
          "f8i"               F8 I     (1150-10620 A)
          "f8iv"              F8 IV    (1150-10620 A)
          "f8v"               F8 V     (1150-10620 A)
          "g0i"               G0 I     (1150-10620 A)
          "g0iii"             G0 III   (1150-10620 A)
          "g0iv"              G0 IV    (1150-10620 A)
          "g0v"               G0 V     (1150-10620 A)
          "g2i"               G2 I     (1150-10620 A)
          "g2iv"              G2 IV    (1150-10620 A)
          "g2v"               G2 V     (1150-10620 A)
          "g5i"               G5 I     (1150-10620 A)
          "g5ii"              G5 II    (1150-10620 A)
          "g5iii"             G5 III   (1150-10620 A)
          "g5iv"              G5 IV    (1150-10620 A)
          "g5v"               G5 V     (1150-10620 A)
          "g8i"               G8 I     (1150-10620 A)
          "g8iii"             G8 III   (1150-10620 A)
          "g8iv"              G8 IV    (1150-10620 A)
          "g8v"               G8 V     (1150-10620 A)
          "k01ii"             K0-1 II  (1150-10620 A)
          "k0iii"             K0 III   (1150-10620 A)
          "k0iv"              K0 IV    (1150-10620 A)
          "k0v"               K0 V     (1150-10620 A)
          "k1iii"             K1 III   (1150-10620 A)
          "k1iv"              K1 IV    (1150-10620 A)
          "k2i"               K2 I     (1150-10620 A)
          "k2iii"             K2 III   (1150-10620 A)
          "k2v"               K2 V     (1150-10620 A)
          "k34ii"             K3-4 II  (1150-10620 A)
          "k3i"               K3 I     (1150-10620 A)
          "k3iii"             K3 III   (1150-10620 A)
          "k3iv"              K3 IV    (1150-10620 A)
          "k3v"               K3 V     (1150-10620 A)
          "k4i"               K4 I     (1150-10620 A)
          "k4iii"             K4 III   (1150-10620 A)
          "k4v"               K4 V     (1150-10620 A)
          "k5iii"             K5 III   (1150-10620 A)
          "k5v"               K5 V     (1150-10620 A)
          "k7v"               K7 V     (1150-10620 A)
          "m0iii"             M0 III   (1150-10620 A)
          "m0v"               M0 V     (1150-10620 A)
          "m10iii"            M10 III  (1150-10620 A)
          "m1iii"             M1 III   (1150-10620 A)
          "m1v"               M1 V     (1150-10620 A)
          "m2i"               M2 I     (1150-10620 A)
          "m2iii"             M2 III   (1150-10620 A)
          "m2p5v"             M2.5 V   (1150-10620 A)
          "m2v"               M2 V     (1150-10620 A)
          "m3ii"              M3 II    (1150-10620 A)
          "m3iii"             M3 III   (1150-10620 A)
          "m3v"               M3 V     (1150-10620 A)
          "m4iii"             M4 III   (1150-10620 A)
          "m4v"               M4 V     (1150-10620 A)
          "m5iii"             M5 III   (1150-10620 A)
          "m5v"               M5 V     (1150-10620 A)
          "m6iii"             M6 III   (1150-10620 A)
          "m6v"               M6 V     (1150-10620 A)
          "m7iii"             M7 III   (1150-10620 A)
          "m8iii"             M8 III   (1150-10620 A)
          "m9iii"             M9 III   (1150-10620 A)
          "o5v"               O5 V     (1150-10620 A)
          "o8iii"             O8 III   (1150-10620 A)
          "o9v"               O9 V     (1150-10620 A)
          "rf6v"              metal-rich F6 V    (1150-10620 A)
          "rf8v"              metal-rich F8 V    (1150-10620 A)
          "rg0v"              metal-rich G0 V    (1150-10620 A)
          "rg5iii"            metal-rich G5 III  (1150-10620 A)
          "rg5v"              metal-rich G5 V    (1150-10620 A)
          "rk0iii"            metal-rich K0 III  (1150-10620 A)
          "rk0v"              metal-rich K0 V    (1150-10620 A)
          "rk1iii"            metal-rich K1 III  (1150-10620 A)
          "rk2iii"            metal-rich K2 III  (1150-10620 A)
          "rk3iii"            metal-rich K3 III  (1150-10620 A)
          "rk4iii"            metal-rich K4 III  (1150-10620 A)
          "rk5iii"            metal-rich K5 III  (1150-10620 A)
          "uka0i"             A0 I     (1150-25000 A)
          "uka0iii"           A0 III   (1150-25000 A)
          "uka0iv"            A0 IV    (1150-25000 A)
          "uka0v"             A0 V     (1150-25000 A)
          "uka2i"             A2 I     (1150-25000 A)
          "uka2v"             A2 V     (1150-25000 A)
          "uka3iii"           A3 III   (1150-25000 A)
          "uka3v"             A3 V     (1150-25000 A)
          "uka47iv"           A4-7 IV  (1150-25000 A)
          "uka5iii"           A5 III   (1150-25000 A)
          "uka5v"             A5 V     (1150-25000 A)
          "uka7iii"           A7 III   (1150-25000 A)
          "uka7v"             A7 V     (1150-25000 A)
          "ukb0i"             B0 I     (1150-25000 A)
          "ukb0v"             B0 V     (1150-25000 A)
          "ukb12iii"          B1-2 III (1150-25000 A)
          "ukb1i"             B1 I     (1150-25000 A)
          "ukb1v"             B1 V     (1150-25000 A)
          "ukb2ii"            B2 II    (1150-25000 A)
          "ukb2iv"            B2 IV    (1150-25000 A)
          "ukb3i"             B3 I     (1150-25000 A)
          "ukb3iii"           B3 III   (1150-25000 A)
          "ukb3v"             B3 V     (1150-25000 A)
          "ukb57v"            B5-7 V   (1150-25000 A)
          "ukb5i"             B5 I     (1150-25000 A)
          "ukb5ii"            B5 II    (1150-25000 A)
          "ukb5iii"           B5 III   (1150-25000 A)
          "ukb6iv"            B6 IV    (1150-25000 A)
          "ukb8i"             B8 I     (1150-25000 A)
          "ukb8v"             B8 V     (1150-25000 A)
          "ukb9iii"           B9 III   (1150-25000 A)
          "ukb9v"             B9 V     (1150-25000 A)
          "ukf02iv"           F0-2 IV  (1150-25000 A)
          "ukf0i"             F0 I     (1150-25000 A)
          "ukf0ii"            F0 II    (1150-25000 A)
          "ukf0iii"           F0 III   (1150-25000 A)
          "ukf0v"             F0 V     (1150-25000 A)
          "ukf2ii"            F2 II    (1150-25000 A)
          "ukf2iii"           F2 III   (1150-25000 A)
          "ukf2v"             F2 V     (1150-25000 A)
          "ukf5i"             F5 I     (1150-25000 A)
          "ukf5iii"           F5 III   (1150-25000 A)
          "ukf5iv"            F5 IV    (1150-25000 A)
          "ukf5v"             F5 V     (1150-25000 A)
          "ukf6v"             F6 V     (1150-25000 A)
          "ukf8i"             F8 I     (1150-25000 A)
          "ukf8iv"            F8 IV    (1150-25000 A)
          "ukf8v"             F8 V     (1150-25000 A)
          "ukg0i"             G0 I     (1150-25000 A)
          "ukg0iii"           G0 III   (1150-25000 A)
          "ukg0iv"            G0 IV    (1150-25000 A)
          "ukg0v"             G0 V     (1150-25000 A)
          "ukg2i"             G2 I     (1150-25000 A)
          "ukg2iv"            G2 IV    (1150-25000 A)
          "ukg2v"             G2 V     (1150-25000 A)
          "ukg5i"             G5 I     (1150-25000 A)
          "ukg5ii"            G5 II    (1150-25000 A)
          "ukg5iii"           G5 III   (1150-25000 A)
          "ukg5iv"            G5 IV    (1150-25000 A)
          "ukg5v"             G5 V     (1150-25000 A)
          "ukg8i"             G8 I     (1150-25000 A)
          "ukg8iii"           G8 III   (1150-25000 A)
          "ukg8iv"            G8 IV    (1150-25000 A)
          "ukg8v"             G8 V     (1150-25000 A)
          "ukk01ii"           K0-1 II  (1150-25000 A)
          "ukk0iii"           K0 III   (1150-25000 A)
          "ukk0iv"            K0 IV    (1150-25000 A)
          "ukk0v"             K0 V     (1150-25000 A)
          "ukk1iii"           K1 III   (1150-25000 A)
          "ukk1iv"            K1 IV    (1150-25000 A)
          "ukk2i"             K2 I     (1150-25000 A)
          "ukk2iii"           K2 III   (1150-25000 A)
          "ukk2v"             K2 V     (1150-25000 A)
          "ukk34ii"           K3-4 II  (1150-25000 A)
          "ukk3i"             K3 I     (1150-25000 A)
          "ukk3iii"           K3 III   (1150-25000 A)
          "ukk3iv"            K3 IV    (1150-25000 A)
          "ukk3v"             K3 V     (1150-25000 A)
          "ukk4i"             K4 I     (1150-25000 A)
          "ukk4iii"           K4 III   (1150-25000 A)
          "ukk4v"             K4 V     (1150-25000 A)
          "ukk5iii"           K5 III   (1150-25000 A)
          "ukk5v"             K5 V     (1150-25000 A)
          "ukk7v"             K7 V     (1150-25000 A)
          "ukm0iii"           M0 III   (1150-25000 A)
          "ukm0v"             M0 V     (1150-25000 A)
          "ukm10iii"          M10 III  (1150-25000 A)
          "ukm1iii"           M1 III   (1150-25000 A)
          "ukm1v"             M1 V     (1150-25000 A)
          "ukm2i"             M2 I     (1150-25000 A)
          "ukm2iii"           M2 III   (1150-25000 A)
          "ukm2p5v"           M2.5 V   (1150-25000 A)
          "ukm2v"             M2 V     (1150-25000 A)
          "ukm3ii"            M3 II    (1150-25000 A)
          "ukm3iii"           M3 III   (1150-25000 A)
          "ukm3v"             M3 V     (1150-25000 A)
          "ukm4iii"           M4 III   (1150-25000 A)
          "ukm4v"             M4 V     (1150-25000 A)
          "ukm5iii"           M5 III   (1150-25000 A)
          "ukm5v"             M5 V     (1150-25000 A)
          "ukm6iii"           M6 III   (1150-25000 A)
          "ukm6v"             M6 V     (1150-25000 A)
          "ukm7iii"           M7 III   (1150-25000 A)
          "ukm8iii"           M8 III   (1150-25000 A)
          "ukm9iii"           M9 III   (1150-25000 A)
          "uko5v"             O5 V     (1150-25000 A)
          "uko8iii"           O8 III   (1150-25000 A)
          "uko9v"             O9 V     (1150-25000 A)
          "ukrf6v"            metal-rich F6 V    (1150-25000 A)
          "ukrf8v"            metal-rich F8 V    (1150-25000 A)
          "ukrg0v"            metal-rich G0 V    (1150-25000 A)
          "ukrg5iii"          metal-rich G5 III  (1150-25000 A)
          "ukrg5v"            metal-rich G5 V    (1150-25000 A)
          "ukrk0iii"          metal-rich K0 III  (1150-25000 A)
          "ukrk0v"            metal-rich K0 V    (1150-25000 A)
          "ukrk1iii"          metal-rich K1 III  (1150-25000 A)
          "ukrk2iii"          metal-rich K2 III  (1150-25000 A)
          "ukrk3iii"          metal-rich K3 III  (1150-25000 A)
          "ukrk4iii"          metal-rich K4 III  (1150-25000 A)
          "ukrk5iii"          metal-rich K5 III  (1150-25000 A)
          "ukwf5v"            metal-weak F5 V    (1150-25000 A)
          "ukwf8v"            metal-weak F8 V    (1150-25000 A)
          "ukwg0v"            metal-weak G0 V    (1150-25000 A)
          "ukwg5iii"          metal-weak G5 III  (1150-25000 A)
          "ukwg5v"            metal-weak G5 V    (1150-25000 A)
          "ukwg8iii"          metal-weak G8 III  (1150-25000 A)
          "ukwk0iii"          metal-weak K0 III  (1150-25000 A)
          "ukwk1iii"          metal-weak K1 III  (1150-25000 A)
          "ukwk2iii"          metal-weak K2 III  (1150-25000 A)
          "ukwk3iii"          metal-weak K3 III  (1150-25000 A)
          "ukwk4iii"          metal-weak K4 III  (1150-25000 A)
          "wf5v"              metal-weak F5 V    (1150-10620 A)
          "wf8v"              metal-weak F8 V    (1150-10620 A)
          "wg0v"              metal-weak G0 V    (1150-10620 A)
          "wg5iii"            metal-weak G5 III  (1150-10620 A)
          "wg5v"              metal-weak G5 V    (1150-10620 A)
          "wg8iii"            metal-weak G8 III  (1150-10620 A)
          "wk0iii"            metal-weak K0 III  (1150-10620 A)
          "wk1iii"            metal-weak K1 III  (1150-10620 A)
          "wk2iii"            metal-weak K2 III  (1150-10620 A)
          "wk3iii"            metal-weak K3 III  (1150-10620 A)
          "wk4iii"            metal-weak K4 III  (1150-10620 A)
        ```
        """
        #
        # Check inputs
        #
        valid_spectral_classes = [
            "a0i",
            "a0iii",
            "a0iv",
            "a0v",
            "a2i",
            "a2v",
            "a3iii",
            "a3v",
            "a47iv",
            "a5iii",
            "a5v",
            "a7iii",
            "a7v",
            "b0i",
            "b0v",
            "b12iii",
            "b1i",
            "b1v",
            "b2ii",
            "b2iv",
            "b3i",
            "b3iii",
            "b3v",
            "b57v",
            "b5i",
            "b5ii",
            "b5iii",
            "b6iv",
            "b8i",
            "b8v",
            "b9iii",
            "b9v",
            "f02iv",
            "f0i",
            "f0ii",
            "f0iii",
            "f0v",
            "f2ii",
            "f2iii",
            "f2v",
            "f5i",
            "f5iii",
            "f5iv",
            "f5v",
            "f6v",
            "f8i",
            "f8iv",
            "f8v",
            "g0i",
            "g0iii",
            "g0iv",
            "g0v",
            "g2i",
            "g2iv",
            "g2v",
            "g5i",
            "g5ii",
            "g5iii",
            "g5iv",
            "g5v",
            "g8i",
            "g8iii",
            "g8iv",
            "g8v",
            "k01ii",
            "k0iii",
            "k0iv",
            "k0v",
            "k1iii",
            "k1iv",
            "k2i",
            "k2iii",
            "k2v",
            "k34ii",
            "k3i",
            "k3iii",
            "k3iv",
            "k3v",
            "k4i",
            "k4iii",
            "k4v",
            "k5iii",
            "k5v",
            "k7v",
            "m0iii",
            "m0v",
            "m10iii",
            "m1iii",
            "m1v",
            "m2i",
            "m2iii",
            "m2p5v",
            "m2v",
            "m3ii",
            "m3iii",
            "m3v",
            "m4iii",
            "m4v",
            "m5iii",
            "m5v",
            "m6iii",
            "m6v",
            "m7iii",
            "m8iii",
            "m9iii",
            "o5v",
            "o8iii",
            "o9v",
            "rf6v",
            "rf8v",
            "rg0v",
            "rg5iii",
            "rg5v",
            "rk0iii",
            "rk0v",
            "rk1iii",
            "rk2iii",
            "rk3iii",
            "rk4iii",
            "rk5iii",
            "uka0i",
            "uka0iii",
            "uka0iv",
            "uka0v",
            "uka2i",
            "uka2v",
            "uka3iii",
            "uka3v",
            "uka47iv",
            "uka5iii",
            "uka5v",
            "uka7iii",
            "uka7v",
            "ukb0i",
            "ukb0v",
            "ukb12iii",
            "ukb1i",
            "ukb1v",
            "ukb2ii",
            "ukb2iv",
            "ukb3i",
            "ukb3iii",
            "ukb3v",
            "ukb57v",
            "ukb5i",
            "ukb5ii",
            "ukb5iii",
            "ukb6iv",
            "ukb8i",
            "ukb8v",
            "ukb9iii",
            "ukb9v",
            "ukf02iv",
            "ukf0i",
            "ukf0ii",
            "ukf0iii",
            "ukf0v",
            "ukf2ii",
            "ukf2iii",
            "ukf2v",
            "ukf5i",
            "ukf5iii",
            "ukf5iv",
            "ukf5v",
            "ukf6v",
            "ukf8i",
            "ukf8iv",
            "ukf8v",
            "ukg0i",
            "ukg0iii",
            "ukg0iv",
            "ukg0v",
            "ukg2i",
            "ukg2iv",
            "ukg2v",
            "ukg5i",
            "ukg5ii",
            "ukg5iii",
            "ukg5iv",
            "ukg5v",
            "ukg8i",
            "ukg8iii",
            "ukg8iv",
            "ukg8v",
            "ukk01ii",
            "ukk0iii",
            "ukk0iv",
            "ukk0v",
            "ukk1iii",
            "ukk1iv",
            "ukk2i",
            "ukk2iii",
            "ukk2v",
            "ukk34ii",
            "ukk3i",
            "ukk3iii",
            "ukk3iv",
            "ukk3v",
            "ukk4i",
            "ukk4iii",
            "ukk4v",
            "ukk5iii",
            "ukk5v",
            "ukk7v",
            "ukm0iii",
            "ukm0v",
            "ukm10iii",
            "ukm1iii",
            "ukm1v",
            "ukm2i",
            "ukm2iii",
            "ukm2p5v",
            "ukm2v",
            "ukm3ii",
            "ukm3iii",
            "ukm3v",
            "ukm4iii",
            "ukm4v",
            "ukm5iii",
            "ukm5v",
            "ukm6iii",
            "ukm6v",
            "ukm7iii",
            "ukm8iii",
            "ukm9iii",
            "uko5v",
            "uko8iii",
            "uko9v",
            "ukrf6v",
            "ukrf8v",
            "ukrg0v",
            "ukrg5iii",
            "ukrg5v",
            "ukrk0iii",
            "ukrk0v",
            "ukrk1iii",
            "ukrk2iii",
            "ukrk3iii",
            "ukrk4iii",
            "ukrk5iii",
            "ukwf5v",
            "ukwf8v",
            "ukwg0v",
            "ukwg5iii",
            "ukwg5v",
            "ukwg8iii",
            "ukwk0iii",
            "ukwk1iii",
            "ukwk2iii",
            "ukwk3iii",
            "ukwk4iii",
            "wf5v",
            "wf8v",
            "wg0v",
            "wg5iii",
            "wg5v",
            "wg8iii",
            "wk0iii",
            "wk1iii",
            "wk2iii",
            "wk3iii",
            "wk4iii",
        ]
        if spectral_class not in valid_spectral_classes:
            raise ValueError(f"{spectral_class} is not a valid `spectral_class`.")
        try:
            data = pd.read_fwf(
                join(DATAPATH, "pickles_spectra", "dat", f"{spectral_class}.dat"),
                colspecs=[(0, 7), (7, 17)],
                header=None,
            )
            self.wavelengths = data[0].values * u.AA
            self.spectrum = data[1].values
        except Exception:
            raise RuntimeError(
                "Could not load the Pickles data for some reason (probably a formatting "
                + "quirk in the file). Please contact the developer with a minimal "
                + "working example."
            )

    def show_spectrum(self, plot=True):
        """
        Plot the spectrum (which should be in units of flam).

        Parameters
        ----------
          plot :: bool
            If True, plot the source weights and return None. If False, return the figure
            and axis instance associated with the plot.

        Returns
        -------
          None (if plot is True)

          fig, ax (if plot is False) :: `matplotlib.figure.Figure`, `matplotlib.axes.Axes`
            The figure and axis instance associated with the plot.
        """
        if self.spectrum is None or self.wavelengths is None:
            raise ValueError("Please generate a spectrum before plotting.")
        fig, ax = plt.subplots()
        ax.plot(self.wavelengths.to(u.AA).value, self.spectrum, "k", lw=1)
        ax.fill_between(self.wavelengths.to(u.AA).value, self.spectrum, alpha=0.5)
        ax.set_xlabel("Wavelength [\AA]")
        ax.set_ylabel(r"Flux Density [$\rm erg\, s^{-1}\, cm^{-2}\,$\AA$^{-1}$]")
        ax.set_ylim(bottom=0)
        if plot:
            plt.show()
        else:
            return fig, ax


class NormMixin:
    """
    Mixin for normalizing spectra.

    TODO: Make normalize to a specific value at a given wavelength
    """

    @staticmethod
    def norm_to_star(spectrum, radius=1, dist=1 << u.kpc):
        """
        Normalize the blackbody spectrum to a star of given radius and distance. By
        default, normalizes the flux to a star of 1 solar radius at 1 kpc. Reference:
        <https://github.com/spacetelescope/pysynphot/blob/925cdbac35a7851cee1bddaa2b47651235c44851/pysynphot/spectrum.py#L40>.

        Note that the spectrum's units should have a unit of steradian (sr) in the
        denominator, which will be multiplied out (e.g., erg/s/cm^2/A/sr).

        Parameters
        ----------
          spectrum :: array of floats
            The spectrum to normalize. It should have a unit of steradian in the
            denominator.

          radius :: float or `astropy.Quantity` length
            The radius of the source. If a scalar, it is assumed to be in units of solar
            radii.

          dist :: float or `astropy.Quantity` length
            The distance to the blackbody. If a scalar, it is assumed to be in units of
            kpc.

        Returns
        ----------
          norm_spectrum :: array of floats
            Spectrum normalized such that the unit of steradian in the denominator
            vanishes.
        """
        #
        # Check inputs
        #
        if not isinstance(radius, u.Quantity):
            radius = radius * const.SUN_RADIUS  # cm
        if not isinstance(dist, u.Quantity):
            dist = dist * u.kpc
        radius = radius.to(u.km).value
        dist = dist.to(u.km).value
        #
        # Normalize spectrum
        #
        radian = radius / dist
        return spectrum * np.pi * radian * radian  # multiply by projected solid angle

    def norm_to_value(
        self,
        value,
        value_type="mag",
        passband=None,
        telescopeObj=None,
        passband_lims=None,
        pivot_wavelength=None,
    ):
        """
        Normalize the spectrum so that it's average flux density (in units of "flam",
        erg/s/cm^2/A) or AB magnitude (based on mean flux density, erg/s/cm^2/A) is equal
        to the given value, either over the whole wavelength range or within a passband.
        The `Source` object should have its spectrum in units of flam (erg/s/cm^2/A) and
        wavelengths in angstrom.

        The higher the resolution of the spectrum, the more accurate the normalization.
        Also, if passband or passband_lims is not None, note that some wavelength/spectrum
        elements at the edges of the passband limits may not be included in the
        normalization due to limited floating point precision. Typically, the relative
        error should be less than 2%.

        To normalize over the whole spectrum, do not specify any passband, passband_lims,
        or pivot_wavelength. To normalize in a passband, either specify:
          - both passband + telescopeObj, or
          - passband_lims and/or pivot_wavelength (useful for normalizing in arbitrary
            passbands).

        Parameters
        ----------
          value :: int or float
            The value to which the average value of the spectrum should equal. This value
            should be in the same units as the spectrum.

          value_type :: "flam" or "mag"
            Flux density (erg/s/cm^2/A) or AB magnitude (mag).

          passband :: valid `Telescope` passband string (e.g., "uv", "u", "g") or None
            If not None, normalize the spectrum such that the average flux density within
            the given telescopeObj passband is equal to the specified value; also requires
            the telescopeObj parameter. Otherwise, normalize the spectrum such that the
            mean flux density is equal to the specified value (after converting to the
            correct units). Note that if the value_type is "mag" and passband is None,
            then the passband response is taken to be unity for the purposes of
            calculating a pivot wavelength. The pivot wavelength is required for this
            normalization otherwise a different normalization factor would apply at each
            wavelength.

          telescopeObj :: `Telescope` object
            The `Telescope` object containing the passband limits and pivot wavelengths of
            each passband. Requires the passband parameter. Only relevant if passband is
            not None.

          passband_lims :: 2-element `astropy.Quantity` array or None
            If not None, this gives the [lower, upper] bounds of the passband, inclusive.
            In this case, the function will normalize the spectrum such that the flux
            density within the given passband_lims is equal to the specified value. This
            cannot be specified if telescopeObj or passband is given. if both passband and
            passband_lims are None, normalize with respect to the entire spectrum.

          pivot_wavelength :: scalar, `astropy.Quantity` length, or None
            The pivot wavelength of the passband_lims. If a scalar, this is assumed to be
            in angstrom. This parameter is only relevant if value_type is "mag". This
            cannot be specified if telescopeObj or passband is given. If both passband and
            passband_lims are None, the normalization applies to the entire spectrum; in
            this case, the "passband response" is taken to be unity. Also see the
            `castor_etc.Telescope.calc_pivot_wavelength()` function or the
            `passband_pivots` attribute in the `Telescope` instance.

        Attributes
        ----------
          spectrum :: array
            Normalized spectrum in units of erg/s/cm^2/A.

        Returns
        -------
          None
        """

        def _get_norm_factor(_value, _flam, _wavelengths, _lims):
            # _tot_flam = simpson(y=_flam, x=_wavelengths, even="avg")  # erg/s/cm^2
            # _avg_flam = _tot_flam / (_lims[1] - _lims[0])  # erg/s/cm^2/A
            _avg_flam = simpson(y=_flam, x=_wavelengths, even="avg") / (
                _lims[1] - _lims[0]
            )  # erg/s/cm^2/A
            return _value / _avg_flam

        #
        # Check inputs
        #
        if self.spectrum is None or self.wavelengths is None:
            raise ValueError("Please generate a spectrum before normalizing.")
        if value_type not in ["flam", "mag"]:
            raise ValueError("value_type must be either 'flam' or 'mag'.")
        if (passband is not None or telescopeObj is not None) and (
            passband_lims is not None or pivot_wavelength is not None
        ):
            raise ValueError(
                "passband/telescopeObj cannot be simultaneously specified with "
                + "passband_lims or pivot_wavelength."
            )
        if passband is not None:
            if telescopeObj is None:
                raise ValueError(
                    "A `Telescope` object must be given if normalizing within a passband."
                )
            try:
                passband_lims = telescopeObj.passband_limits[passband]
                pivot_wavelength = telescopeObj.passband_pivots[passband]
            except Exception:
                raise AttributeError("Desired passband not found in given telescopeObj.")
        if isinstance(pivot_wavelength, u.Quantity):
            pivot_wavelength = pivot_wavelength.to(self.wavelengths.unit)
        if value_type == "mag" and pivot_wavelength is None:
            if passband_lims is None:
                in_passband = np.full(np.shape(self.wavelengths), True)
            else:
                in_passband = (self.wavelengths >= passband_lims[0]) & (
                    self.wavelengths <= passband_lims[1]
                )
            pivot_wavelength = (
                Telescope.calc_pivot_wavelength(
                    self.wavelengths.value[in_passband],
                    np.ones(np.shape(self.wavelengths[passband_lims]), dtype=float),
                )
                * self.wavelengths.unit
            )
        #
        # Calculate normalization factor
        #
        if value_type == "mag":
            # Convert desired AB magnitude to flam
            value = convert_electron_flux_mag(
                value, "mag", "flam", wavelengths=pivot_wavelength
            )[0]
        #
        if passband is None and passband_lims is None:
            wavelengths_AA = self.wavelengths.to(u.AA).value
            norm_factor = _get_norm_factor(
                value,
                self.spectrum,
                wavelengths_AA,
                [wavelengths_AA[0], wavelengths_AA[-1]],
            )
        else:
            if passband_lims is not None:
                if np.size(passband_lims) != 2 or not isinstance(
                    passband_lims, u.Quantity
                ):
                    raise ValueError(
                        "passband_lims must be a 2-element `astropy.Quantity` array."
                    )
            else:  # passband is None
                if not isinstance(passband, str):
                    raise TypeError("passband must be a CASTOR passband string or None.")
                if passband not in params.PASSBANDS:
                    raise ValueError(
                        "Invalid CASTOR passband. "
                        + f"Valid passbands are: {params.PASSBANDS}"
                    )
                passband_lims = params.PASSBAND_LIMITS[passband]
            is_in_passband = (self.wavelengths >= passband_lims[0]) & (
                self.wavelengths <= passband_lims[1]
            )
            wavelengths_AA = self.wavelengths[is_in_passband].to(u.AA).value
            norm_factor = _get_norm_factor(
                value,
                self.spectrum[is_in_passband],
                wavelengths_AA,
                passband_lims.to(u.AA).value,
            )
        self.spectrum *= norm_factor

    def norm_luminosity_dist(self, luminosity, dist):
        """
        Normalize the spectrum to a source of given (bolometric) luminosity and distance.
        The `Source` object should have its spectrum in units of flam (erg/s/cm^2/A) and
        wavelengths in angstrom. (Technically it is okay as long as the wavelengths are in
        some unit <U> and the spectrum is in units of erg/s/cm^2/<U>.)

        Parameters
        ----------
          luminosity :: scalar or `astropy.Quantity` unit of power
            The desired bolometric luminosity. If a scalar, this is assumed to be in units
            of solar luminosities. If an `astropy.Quantity`, it must be a unit of power
            (e.g., erg/s).

          dist :: float or `astropy.Quantity` length
            The distance to the source. If a scalar, it is assumed to be in units of kpc.

        Attributes
        ----------
          spectrum :: array
            Normalized spectrum in original flux density units (e.g., erg/s/cm^2/A).

        Returns
        -------
          None
        """
        #
        # Check inputs
        #
        if self.spectrum is None or self.wavelengths is None:
            raise ValueError("Please generate a spectrum before normalizing.")
        if not isinstance(luminosity, u.Quantity):
            luminosity = luminosity * const.SUN_LUMINOSITY
        if not isinstance(dist, u.Quantity):
            dist = dist * u.kpc
        luminosity = luminosity.to(u.erg / u.s).value
        dist = dist.to(u.cm).value
        #
        # Normalize spectrum (originally in erg/s/cm^2/A)
        #
        erg_s_A = 4 * np.pi * dist * dist * self.spectrum  # erg/s/A
        tot_luminosity = simpson(y=erg_s_A, x=self.wavelengths.value, even="avg")  # erg/s
        norm_factor = luminosity / tot_luminosity  # dimensionless
        self.spectrum *= norm_factor  # erg/s/cm^2/A

    def get_avg_value(self, value_type="mag", TelescopeObj=None):
        """
        Calculate the average value of the source spectrum (in either flam, erg/s/cm^2/A,
        or AB magnitude) either over the whole spectrum or through a telescope's passbands
        (in which case TelescopeObj is required).

        Parameters
        ----------
          value_type :: "flam" or "mag"
            The desired output value type. If "flam", the output is in units of
            erg/s/cm^2/A. If "mag", the output is in AB magnitudes and the pivot
            wavelength is automatically calculated assuming a flat response function
            (i.e., perfect throughput/response over all the entire spectrum).

          TelescopeObj :: `Telescope` object or None
            If provided, will calculate the average value in each of the telescope's
            passbands.

        Returns
        -------
          result :: scalar or dict of scalars
            If TelescopeObj is None, the result is a scalar equal to the average value
            (either erg/s/cm^2/A or AB magnitude) of the entire spectrum. If TelescopeObj
            is not None, the result is a dict of average values in each of the telescope's
            passbands.
        """
        if value_type not in ["flam", "mag"]:
            raise ValueError("value_type must be either 'flam' or 'mag'.")
        source_AA = self.wavelengths.to(u.AA).value
        if TelescopeObj is None:
            avg_flam = simpson(y=self.spectrum, x=source_AA, even="avg") / (
                source_AA[-1] - source_AA[0]
            )  # erg/s/cm^2/A
            if value_type == "mag":
                pivot_wavelength = (
                    Telescope.calc_pivot_wavelength(
                        source_AA,
                        np.ones(np.shape(source_AA), dtype=float),  # perfect response
                    )
                    * u.AA
                )
                result = convert_electron_flux_mag(
                    avg_flam,
                    "flam",
                    "mag",
                    wavelengths=pivot_wavelength,
                )[0]
            else:
                result = avg_flam
        else:
            result = dict.fromkeys(TelescopeObj.passbands)
            for band in TelescopeObj.passbands:
                passband_lims_AA = TelescopeObj.passband_limits[band].to(u.AA).value
                is_in_passband = (source_AA >= passband_lims_AA[0]) & (
                    source_AA <= passband_lims_AA[1]
                )
                avg_flam = (  # erg/s/cm^2/A
                    simpson(
                        y=self.spectrum[is_in_passband],
                        x=source_AA[is_in_passband],
                        even="avg",
                    )
                    / (passband_lims_AA[1] - passband_lims_AA[0])
                )
                if value_type == "mag":
                    result[band] = convert_electron_flux_mag(
                        avg_flam,
                        "flam",
                        "mag",
                        wavelengths=TelescopeObj.passband_pivots[band],
                    )[0]
                else:
                    result[band] = avg_flam
        return result
