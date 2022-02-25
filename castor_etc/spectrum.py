"""
Generate and handle spectral data and normalizations. Includes:
  - blackbody radiation, power-law spectrum, Gaussian spectrum
  - stellar spectral types (TODO)
  - user-input spectrum (TODO)
  - normalizations (LIST THEM)
"""

import warnings
from numbers import Number
from os.path import join

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.integrate import simpson
from scipy.interpolate import interp1d

from . import constants as const
from . import parameters as params
from .conversions import calc_photon_energy, convert_electron_flux_mag
from .filepaths import DATAPATH
from .telescope import Telescope

# TODO: Pickles spectra


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
        Use custom spectrum from an ASCII file.

        Parameters
        ----------
          filepath :: str
            The absolute path to the file containing the spectrum. The files should be in
            ASCII format with the first column containing the wavelengths in
            `wavelength_units` and the second column containing the spectrum in flam
            (erg/s/cm^2/A); the columns should be separated by a constant number of
            spaces. Lines starting with a hash (#) will be ignored.

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
        data = pd.read_csv(
            filepath,
            sep=" +",
            header=None,
            comment="#",
            engine="python",
        )  # sep=" +" is Python regex to match a variable number of spaces
        self.wavelengths = (data[0].values * wavelength_unit).to(u.AA)
        self.spectrum = data[1].values

    def use_galaxy_spectrum(self, type, overwrite=False):
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
          type :: "elliptical" or "spiral"
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
        if type == "elliptical" or type == "spiral":
            filepath = join(DATAPATH, "galaxy_spectra", f"{type}_galaxy.txt")
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
        value_type="flux",
        passband=None,
        passband_lims=None,
        pivot_wavelength=None,
    ):
        """
        Normalize the spectrum so that it's average flux density (in erg/s/cm^2/A) or AB
        magnitude is equal to the given value, either over the whole wavelength range or
        within a passband. The `Source` object should have its spectrum in units of flam
        (erg/s/cm^2/A) and wavelengths in angstrom.

        The higher the resolution of the spectrum, the more accurate the normalization.
        Also, if passband or passband_lims is not None, note that some wavelength/spectrum
        elements at the edges of the passband limits may not be included in the
        normalization due to limited floating point precision. Typically, the relative
        error should be less than 2%.

        Parameters
        ----------
          value :: int or float
            The value to which the average value of the spectrum should equal. This value
            should be in the same units as the spectrum.

          value_type :: "flux" or "mag"
            Flux density (erg/s/cm^2/A) or AB magnitude (mag).

          passband :: "uv", "u", "g", or None
            If not None, normalize the spectrum such that the flux within the given
            (predefined) CASTOR passband is equal to the specified value. At most one of
            passband or passband_lims can be specified; if both passband and passband_lims
            are None, normalize the spectrum such that the total flux is equal to the
            specified value. Note that if your `Telescope` object uses custom passbands,
            you should use the passband_lims argument instead (i.e., pass in
            `passband_limits` attribute to passband_lims).

          passband_lims :: 2-element `astropy.Quantity` array or None
            If not None, this gives the [lower, upper] bounds of the passband, inclusive.
            In this case, the function will normalize the spectrum such that the flux
            within the given passband_lims is equal to the specified value. At most one of
            passband or passband_lims can be specified; if both passband and passband_lims
            are None, normalize the spectrum such that the total flux is equal to the
            specified value.

          pivot_wavelength :: scalar, `astropy.Quantity` length, or None
            The pivot wavelength of the passband/passband_lims. If a scalar, this is
            assumed to be in angstrom. This parameter is only relevant (and must not be
            None) if value_type is "mag"; the only exception is if passband and
            passband_lims are both None so that the normalization applies to the
            bolometric flux---in this case, the "passband response" is taken to be unity.
            Also see the `castor_etc.Telescope.calc_pivot_wavelength()` function or the
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
            _avg_flam = simpson(y=_flam, x=_wavelengths, even="avg") / (
                _lims[1] - _lims[0]
            )  # erg/s/cm^2/A
            return _value / _avg_flam

        #
        # Check inputs
        #
        if self.spectrum is None or self.wavelengths is None:
            raise ValueError("Please generate a spectrum before normalizing.")
        if passband is not None and passband_lims is not None:
            raise ValueError("At most one of passband and passband_lims can be not None.")
        if value_type not in ["flux", "mag"]:
            raise ValueError("value_type must be either 'flux' or 'mag'.")
        if value_type == "mag" and pivot_wavelength is None:
            if passband is None and passband_lims is None:
                pivot_wavelength = (
                    Telescope.calc_pivot_wavelength(
                        self.wavelengths.value,
                        np.ones(np.shape(self.wavelengths), dtype=float),
                    )
                    * self.wavelengths.unit
                )
            else:
                raise ValueError(
                    "pivot_wavelength must be specified if "
                    + "passband/passband_lims is not None and value_type is 'mag'."
                )
        if isinstance(pivot_wavelength, u.Quantity):
            pivot_wavelength = pivot_wavelength.to(self.wavelengths.unit)
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


# ------------------------------------ OLD FUNCTIONS ----------------------------------- #


def generate_gaussian(center, wavelengths, fwhm, peak=None, tot_flux=None):
    """
    Generates a spectrum in the shape of a Gaussian. This is useful for representing
    emission lines (i.e., by adding a Gaussian source) or absorption lines (i.e., by
    subtracting a Gaussian source).

    The spectrum can be represented by the following formulae:
    ```math
                spectrum = peak / exp[(wavelengths - center)^2 / (2 * sigma^2)]
    ```
    and
    ```math
                sigma = fwhm / [2 * sqrt(2 * ln2)]
    ```
    and
    ```math
                peak = tot_flux / sqrt(2 * pi * sigma^2)
    ```
    where:
      - spectrum is the spectrum's flux in some arbitrary unit
      - peak is the flux at the center of the Gaussian (i.e., the central wavelength)
      - center is the central wavelength of the Gaussian
      - wavelengths is the array of wavelengths over which to calculate the spectrum
      - fwhm is the full-width at half-maximum of the Gaussian
      - tot_flux is the total flux of the Gaussian under the curve

    Parameters
    ----------
      center :: scalar or `astropy.Quantity`
        The central wavelength of the Gaussian. If an `astropy.Quantity` object, it must
        have the same units as the wavelengths array.

      wavelengths :: array of floats or `astropy.Quantity` array
        The wavelengths over which to calculate the Gaussian spectrum. If an
        `astropy.Quantity` array, it must have the same units as the center parameter.

      fwhm :: scalar or `astropy.Quantity`
        The full-width at half-maximum of the Gaussian. If an `astropy.Quantity` object,
        it must have the same units as the wavelengths array.

      peak :: scalar or `astropy.Quantity`
        The peak flux of the Gaussian (i.e., the flux at the center wavelength). This
        determines the unit of the returned spectrum. Exactly one of peak or tot_flux must
        be specified.

      tot_flux :: scalar or `astropy.Quantity`
        The total flux under the curve. This implicitly determines the unit of the
        returned spectrum. Exactly one of peak or tot_flux must be specified.

    Returns
    -------
      spectrum :: array of floats or `astropy.Quantity` array
        The flux of the source in the shape of a Gaussian curve.
    """
    #
    # Check inputs
    #
    if np.size(center) != 1:
        raise ValueError("center must be a single scalar or `astropy.Quantity`.")
    if isinstance(wavelengths, u.Quantity) and isinstance(center, u.Quantity):
        if wavelengths.unit != center.unit:
            raise ValueError("wavelengths and center must have the same units.")
    if np.size(fwhm) != 1:
        raise ValueError("fwhm must be a single scalar or `astropy.Quantity`.")
    if isinstance(fwhm, u.Quantity):
        if isinstance(center, u.Quantity) and center.unit != fwhm.unit:
            raise ValueError("center and fwhm must have the same units.")
        if isinstance(wavelengths, u.Quantity) and wavelengths.unit != fwhm.unit:
            raise ValueError("wavelengths and fwhm must have the same units.")
    if (peak is None and tot_flux is None) or (peak is not None and tot_flux is not None):
        raise ValueError("Exactly one of peak or tot_flux must be specified.")
    if peak is not None and np.size(peak) != 1:
        raise ValueError("peak must be a single scalar or `astropy.Quantity`.")
    elif np.size(tot_flux) != 1:
        raise ValueError("tot_flux must be a single scalar or `astropy.Quantity`.")
    #
    # Gaussian spectrum
    #
    sigma = fwhm / (2 * np.sqrt(2 * np.log(2)))
    if peak is None:
        peak = tot_flux / (np.sqrt(2 * np.pi) * sigma)
    spectrum = peak / np.exp((wavelengths - center) ** 2 / (2 * sigma ** 2))
    return spectrum
