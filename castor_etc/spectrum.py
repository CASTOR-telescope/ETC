"""
spectrum.py

Generate and handle spectral data and normalizations. Includes:
- blackbody radiation
- stellar spectral types (TODO)
- user-input spectrum (TODO)
- normalizations (TODO and LIST THEM)

Isaac Cheng - 2022
"""

import astropy.units as u
import numpy as np

from scipy.integrate import simpson

from . import constants as const
from . import parameters as params
from .energy import calc_photon_energy


def generate_BB(
    T,
    redshift=0.0,
    emissivity=1.0,
    wavelengths=None,
    limits=params.PASSBAND_TOT_LIMITS,
    resolution=1 << u.nm,
):
    """
    Generate a blackbody (BB) spectrum (photon/s/cm^2/A/sr) using Planck's radiation law.
    Also useful for continuum-normalization (e.g., see
    <https://pysynphot.readthedocs.io/en/latest/tutorials.html#tutorial-5-continuum-normalized-spectrum>).

    Parameters
    ----------
      T :: int or float or `astropy.Quantity`
        Intrinsic blackbody temperature (i.e., the temperature of the BB at redshift=0).
        If int or float, the unit is assumed to be kelvin.

      redshift :: int or float
        Redshift of the blackbody.

      emissivity :: int or float
        Emissivity of the blackbody. (Technically, emissivity is unity per the definition
        of a BB).

      wavelengths :: array of floats or `astropy.Quantity` array
        The wavelengths over which to calculate the spectrum. If an array of floats, the
        unit is assumed to be in angstroms. If wavelengths is not None, the limits and
        resolution parameters are ignored.

      limits :: list of 2 scalars or list of 2 `astropy.Quantity`
        List containing the lower (0th index) and upper (1st index) bounds for the BB
        spectrum's restframe wavelengths, inclusive. Limits should be > 0. If list
        elements are int or float, they are assumed to be in angstroms. This parameter is
        ignored if wavelengths is provided.

      resolution :: int or float or `astropy.Quantity`
        The wavelength resolution of the returned spectrum. If a scalar, it is assumed to
        be in units of angstroms. This parameter is ignored if wavelengths is provided.

    Returns
    -------
      wavelengths :: array of floats
        The redshifted wavelengths of the spectrum, in angstroms.
      spectrum :: array of floats
        Spectral radiance of BB in units of photon/s/cm^2/A/sr.
    """
    #
    # Check inputs
    #
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
    # Factor in redshift & convert angstrom to cm
    wavelengths = wavelengths * (1 + redshift) * 1e-8  # cm
    # Planck's radiation law
    lightspeed = const.LIGHTSPEED.value  # cm/s
    prefactor = (2 * const.PLANCK_H.value * lightspeed * lightspeed) / (wavelengths ** 5)
    denom = np.expm1(
        (const.PLANCK_H.value * lightspeed) / (wavelengths * const.K_B.value * T)
    )
    spectrum = prefactor / denom  # erg/s/cm^2/cm/sr
    #
    # Incorporate emissivity and convert per cm to per angstrom
    #
    spectrum *= 1e-8 * emissivity  # erg/s/cm^2/A/sr
    # spectrum *= 1e-8 * emissivity / const.SR_TO_SQARCSEC  # erg/s/cm^2/A/arcsec^2
    #
    # Convert to photon/s/cm^2/A/sr
    #
    spectrum /= calc_photon_energy(wavelength=wavelengths * u.cm)[0]  # photon/s/cm^2/A/sr
    #
    # Convert wavelengths back to angstroms
    #
    wavelengths *= 1e8  # angstrom
    return wavelengths, spectrum


def generate_power_law(ref_wavelength, wavelengths, exponent):
    """
    Generate a spectrum in some arbitrary unit following a power-law. The flux is defined
    so that it is equal to 1 at the reference wavelength.

    The spectrum is calculated using the following formula:
    ```math
                spectrum = (wavelength / ref_wavelength) ** exponent
    ```
    where each variable is as defined in the Parameters documentation below.

    Parameters
    ----------
      ref_wavelength :: scalar or `astropy.Quantity`
        The reference wavelength for the power-law. The spectrum at this wavelength will
        have a flux of 1. If an `astropy.Quantity` object, it must have the same units as
        the wavelengths array.

      wavelengths :: array of floats or `astropy.Quantity` array
        The wavelengths over which to calculate the power-law spectrum. If an
        `astropy.Quantity` array, it must have the same units as the ref_wavelength
        object.

      exponent :: int or float
        The exponent for the power-law.

    Returns
    -------
      spectrum :: array of floats
        The spectrum in the same units as the wavelengths array.
    """
    #
    # Check inputs
    #
    if np.size(ref_wavelength) != 1:
        raise ValueError("ref_wavelength must be a single scalar or `astropy.Quantity`.")
    if isinstance(wavelengths, u.Quantity) and isinstance(ref_wavelength, u.Quantity):
        if wavelengths.unit != ref_wavelength.unit:
            raise ValueError("wavelengths and ref_wavelength must have the same units.")
        else:
            ref_wavelength = ref_wavelength.value
            wavelengths = wavelengths.value
    #
    # Power-law
    #
    spectrum = (wavelengths / ref_wavelength) ** exponent
    return spectrum


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


def norm_to_star(spectrum, radius=const.SUN_RADIUS, dist=1 << u.kpc):
    """
    Normalize a blackbody spectrum to a star of given radius and distance. By default,
    normalizes the flux to a star of 1 solar radius at 1 kpc. Reference:
    <https://github.com/spacetelescope/pysynphot/blob/925cdbac35a7851cee1bddaa2b47651235c44851/pysynphot/spectrum.py#L40>.

    Parameters
    ----------
      spectrum :: array of floats
        The spectrum to be normalized. The spectrum's units should have a unit of
        steradian (sr) in the denominator, which will be multiplied out (e.g.,
        photon/s/cm^2/A/sr).

      radius :: float or `astropy.Quantity`
        The radius of the source. If a scalar, it is assumed to be in units of solar
        radii.

      dist :: float or `astropy.Quantity`
        The distance to the blackbody. If a scalar, it is assumed to be in units of kpc.

    Returns
    -------
      norm_spectrum :: array of floats
        The normalized spectrum.
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
    norm_factor = radius / dist  # radian
    norm_spectrum = np.pi * norm_factor * norm_factor * spectrum
    return norm_spectrum


def norm_to_value(spectrum, wavelengths, value):
    """
    Normalize the spectrum so that it has a total flux equal to some value (i.e., the area
    under the spectrum is value). This function can also be used to normalize the flux in
    a passband by simply passing in the passband's spectrum and wavelengths instead of the
    entire spectrum + wavelengths.

    Parameters
    ----------
      spectrum :: array of floats
        The spectrum (or portion of the spectrum) to be normalized.

      wavelegths :: array of floats
        The wavelengths corresponding to the spectrum.

      value :: int or float
        The value to which the area under the spectrum should equal. This value should be
        in the same units as the spectrum integrated over its wavelength range.

    Returns
    -------
      norm_spectrum :: array of floats
        The normalized spectrum.
    """
    tot_flux = simpson(y=spectrum, x=wavelengths, even="avg")
    norm_spectrum = value * spectrum / tot_flux
    return norm_spectrum
