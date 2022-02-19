"""
Calculations involving signal-to-noise ratio (SNR).
"""

import astropy.units as u
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d

from . import parameters as params
from .data.background.background_values import (
    GEOCORONAL_FLUX_AVG,
    GEOCORONAL_LINEWIDTH,
    GEOCORONAL_WAVELENGTH,
    SKY_BACKGROUND,
)
from .energy import calc_photon_energy
from .filepaths import DATAPATH
from .load_files import load_passbands, load_sky_background
from .conversions import mag_to_flux


def _calc_snr_from_t(
    t,
    signal,
    npix,
    totskynoise=0.0,
    readnoise=params.READ_NOISE,
    darkcurrent=params.DARK_CURRENT,
    redleak=0.0,
    nread=1,
):
    """
    Calculate the signal-to-noise ratio (SNR) reached given an integration time. Based on
    <https://hst-docs.stsci.edu/wfc3ihb/chapter-9-wfc3-exposure-time-calculation/9-6-estimating-exposure-times>.

    The equation to calculate the SNR is:
    ```math
                SNR = (Q*t) / sqrt(Q*t + N_pix*t*Poisson + N_pix*N_read*Read^2)
    ```
    where:
      - SNR is the signal-to-noise ratio
      - t is the integration time in seconds
      - Q is the total signal due to the source in electrons/s
      - N_pix is the number of pixels occupied by the source on the detector
      - Poisson = (B_Earthshine + B_zodiacal + B_geocoronal + Darkcurrent + Redleak) is
        the total Poisson noise due to the sky backgorund (Earthshine, zodiacal light, and
        geocoronal emission), dark current, and red leak in electrons/s (per pixel)
      - N_read is the number of detector readouts
      - Read is the detector read noise in electrons (per pixel)


    Parameters
    ----------
      t :: int or float
        The integration time in seconds.

      signal :: scalar or array of scalars
        The signals of the source in electron/s. These are the total electron rates over
        the whole npix (see below).

      npix :: int or float
        The number of pixels occupied by the source.

      totskynoise :: int or float
        The total background noise due to Earthshine + zodiacal light + geocoronal
        emission (if applicable) in electron/s (per pixel).

      readnoise :: int or float
        The CCD read noise in electron (per pixel).

      darkcurrent :: int or float
        The CCD dark current in electron/s (per pixel).

      redleak :: int or float
        The CCD read leak due to the source in electron/s (per pixel).

      nread :: int
        The number of CCD readouts.

    Returns
    -------
      snr :: float or array of floats
        The SNRs reach after t seconds.
    """
    #
    # Check inputs
    #
    variables = [t, signal, npix, totskynoise, readnoise, darkcurrent, redleak, nread]
    if np.any([isinstance(var, u.Quantity) for var in variables]):
        raise ValueError(
            "All inputs must be scalars. `astropy.Quantity` objects are not supported."
        )
    if not isinstance(nread, (int, np.integer)):
        raise ValueError("nread must be an integer")
    #
    # Calculate signal-to-noise ratio
    #
    signal_t = signal * t  # electron
    noise = np.sqrt(
        signal_t
        + npix
        * (t * (totskynoise + darkcurrent + redleak) + (readnoise * readnoise * nread))
    )  # electron
    snr = signal_t / noise
    return snr


def _calc_t_from_snr(
    snr,
    signal,
    npix,
    totskynoise=0.0,
    readnoise=params.READ_NOISE,
    darkcurrent=params.DARK_CURRENT,
    redleak=0.0,
    nread=1,
):
    """
    Calculate the time required to reach a given signal-to-noise ratio (SNR). Based on
    <https://hst-docs.stsci.edu/wfc3ihb/chapter-9-wfc3-exposure-time-calculation/9-6-estimating-exposure-times>.

    The equation to calculate the time required to reach a given SNR is:
    ```math
                t = {
                    SNR^2 * (Q + N_pix * Poisson)
                    + sqrt[
                        SNR^4 * (Q + N_pix * Poisson)^2
                        + 4 * Q^2 * SNR^2 * N_pix * N_read * Read^2
                    ]
                } / (2 * Q^2)
    ```
    where:
      - SNR is the desired signal-to-noise ratio
      - t is the integration time in seconds
      - Q is the total signal due to the source in electrons/s
      - N_pix is the number of pixels occupied by the source on the detector
      - Poisson = (B_Earthshine + B_zodiacal + B_geocoronal + Darkcurrent + Redleak) is
        the total Poisson noise due to the sky backgorund (Earthshine, zodiacal light, and
        geocoronal emission), dark current, and red leak in electrons/s (per pixel)
      - N_read is the number of detector readouts
      - Read is the detector read noise in electrons (per pixel)
    Note that this is simply the quadratic formula applied to the generic SNR equation.

    Parameters
    ----------
      snr :: int or float
        The target SNR.

      signal :: scalar or array of scalars
        The signals of the source in electron/s. These are the total electron rates over
        the whole npix (see below).

      npix :: int or float
        The number of pixels occupied by the source.

      totskynoise :: int or float
        The total background noise due to Earthshine + zodiacal light + geocoronal
        emission (if applicable) in electron/s (per pixel).

      readnoise :: int or float
        The CCD read noise in electron (per pixel).

      darkcurrent :: int or float
        The CCD dark current in electron/s (per pixel).

      redleak :: int or float
        The CCD read leak due to the source in electron/s (per pixel).

      nread :: int
        The number of CCD readouts.

    Returns
    -------
      t :: float or array of floats
        The times required to reach the given SNR in seconds.
    """
    variables = [snr, signal, npix, totskynoise, readnoise, darkcurrent, redleak, nread]
    if np.any([isinstance(var, u.Quantity) for var in variables]):
        raise ValueError(
            "All inputs must be scalars. `astropy.Quantity` objects are not supported."
        )
    if not isinstance(nread, (int, np.integer)):
        raise ValueError("nread must be an integer")
    #
    # Calculate useful quantities
    #
    snr_sq = snr * snr
    signal_sq = signal * signal
    poisson_noise = signal + npix * (totskynoise + darkcurrent + redleak)
    #
    # Calculate time to reach target SNR
    #
    numer1 = snr_sq * poisson_noise
    numer2 = snr_sq * snr_sq * poisson_noise * poisson_noise
    numer3 = 4 * snr_sq * signal_sq * (npix * readnoise * readnoise * nread)
    t = (numer1 + np.sqrt(numer2 + numer3)) / (2 * signal_sq)  # seconds
    return t


def _calc_red_leak(wavelengths, passband_response, source_flux, redleak_threshold):
    """
    Calculates the red leak of a source through a given passband response following the
    method described in
    <https://hst-docs.stsci.edu/wfc3ihb/chapter-6-uvis-imaging-with-wfc3/6-5-uvis-spectral-elements#id-6.5UVISSpectralElements-6.5.26.5.2FilterRedLeaks>.
    The passband response and source flux arrays should be the same shape and correspond
    to the same wavelengths.

    Parameters
    ----------
      wavelengths :: array of floats or `astropy.Quantity` array
        The wavelengths corresponding to the passband response and source flux arrays. If
        an array of floats, the units are assumed to be angstroms. Otherwise, the array
        can be converted to an `astropy.Quantity` array using, e.g.,
        `new_arr = arr * u.nm`.

      passband_response :: array of floats
        The passband response. This should convert the source flux to units of electron/s.
        This must have the same shape as the wavelengths array and correspond to to the
        wavelengths array element-wise.

      source_flux :: array of floats
        The flux of the source. This must have the same shape as the wavelengths array and
        correspond to to the wavelengths array element-wise.

      redleak_threshold :: int or float or `astropy.Quantity`
        The wavelength beyond which any flux will be considered "red leak". That is, if
        the wavelength > redleak_threshold, the flux convolved with the passband_response
        will be considered the red leak. If int or float, the units are assumed to be
        angstroms.

    Returns
    -------
      redleak :: float or array of floats
        The red leak of the source in electron/s.
    """
    #
    # Check inputs
    #
    if (np.shape(passband_response) != np.shape(wavelengths)) or (
        np.shape(source_flux) != np.shape(wavelengths)
    ):
        raise ValueError(
            "Passband response and source flux must have the same shape "
            + "as wavelengths array"
        )
    if isinstance(wavelengths, u.Quantity):
        wavelengths = wavelengths.to(u.AA).value
    if isinstance(redleak_threshold, u.Quantity):
        redleak_threshold = redleak_threshold.to(u.AA).value
    #
    # Convolve passband response with source flux
    #
    electron_rate = passband_response * source_flux  # electron/s
    #
    # Calculate red leak
    #
    is_redleak = wavelengths > redleak_threshold
    redleak = np.nansum(electron_rate[is_redleak])  # electron/s
    return redleak


def calc_snr_AB(
    source_AB_mags,
    snr=None,
    t=None,
    passbands="all",
    passband_response_files=None,
    passband_limits=None,
    sky_background_file=None,  # TODO: implement this (later)
    interp_resolution=params.PASSBAND_RESOLUTION,
    interp_kind="linear",
    zpts=params.PHOT_ZPTS,
    npix=(np.pi * (1.4 * params.FWHM / 2) ** 2 / (params.PX_SCALE) ** 2).value,
    readnoise=params.READ_NOISE,
    darkcurrent=params.DARK_CURRENT,
    redleak=None,  # TODO: implement this (once I figure out the cause of the discrepancy)
    include_geocoronal=True,
    nread=1,
):
    """
    ! DEPRECATED !
    Given some AB magnitudes in CASTOR passbands, calculate the signal-to-noise ratio
    (SNR) reached in a given time, t, OR the time required to reach a given SNR, snr.

    Parameters
    ----------
      npix :: int or float
        The number of pixels occupied by the source. The default value is for an "optimal"
        point source with an angular radius of 1.4 * FWHM.

    Returns
    -------
      t_lim :: dict of float arrays
        Dictionary containing the time in seconds required to reach the given SNR for each
        source. The dictionary keys are the passbands, and the values are arrays of floats
        in the order of the input sources.

    TODO: finish docstring (later. Will probably change a lot)
    """
    #
    # Check inputs
    #
    if snr is None and t is None:
        raise ValueError("Exactly one of snr or t must be specified")
    elif snr is not None and t is not None:
        raise ValueError("Only one of snr or t can be specified")
    #
    if redleak is not None:
        raise NotImplementedError("Red leak is not implemented yet!")
    redleak = 0.0
    if sky_background_file is not None:
        raise ValueError("Sky background file is not implemented yet!")
    #
    if isinstance(passbands, str):
        passbands = [passbands]
    if np.ndim(passbands) != 1:
        raise ValueError("passbands must be a 1D list of strings or 'all'")
    if passbands == ["all"]:
        passbands = params.PASSBANDS
    else:
        if any(band not in params.PASSBANDS for band in passbands):
            raise ValueError(f"Invalid filters. Valid filters are: {params.PASSBANDS}")
    #
    if passband_limits is None:
        passband_limits = {band: params.PASSBAND_LIMITS[band] for band in passbands}
    # Ensure passband_limits is astropy Quantity object
    if passband_limits is not None and not isinstance(passband_limits, u.Quantity):
        for limit in passband_limits.values():
            if not isinstance(limit, u.Quantity):
                raise ValueError(
                    "passband_limits must be an `astropy.Quantity` object "
                    + "(e.g., passband_limits={'uv': [150, 300] * u.nm})"
                )
    # Determine minimum and maximum wavelength covered by the chosen passbands
    passband_tot_limits = [
        min(passband_limits.values(), key=lambda x: x[0])[0],
        max(passband_limits.values(), key=lambda x: x[1])[1],
    ]
    #
    # Load passband response curves
    #
    passband_response = load_passbands(
        filters=passbands,
        limits=passband_limits,
        resolution=interp_resolution,
        filepaths=passband_response_files,
        kind=interp_kind,
    )
    #
    # Load Earthshine and zodiacal light
    #
    # sky_background_flam_persqarcsec = load_sky_background(
    #     resolution=interp_resolution,
    #     limits=passband_tot_limits,
    #     filepath=sky_background_file,
    #     kind=interp_kind,
    # )  # erg/cm^2/s/A/arcsec^2
    phot_zpts = pd.DataFrame(zpts, index=[0])  # AB mag for 1 electron/s
    sky_background = pd.DataFrame(SKY_BACKGROUND, index=[0])  # AB mag/arcsec^2
    sky_background_e_rate = (
        mag_to_flux(sky_background, zpt=phot_zpts)[0] * params.PX_AREA.value
    )  # electron/s
    #
    # Add geocoronal emission line
    #
    if include_geocoronal:
        geo_background = (
            GEOCORONAL_FLUX_AVG
            * params.PX_AREA.value
            * params.MIRROR_AREA.value
            / calc_photon_energy(
                wavelength=GEOCORONAL_WAVELENGTH.value, wavelength_err=0.0
            )[0]
        )  # photon/s/A
        for band in passbands:
            throughput = passband_response[band]["throughput"].values
            band_start = (
                (passband_response[band]["wavelength"].iloc[0] * u.um).to(u.AA).value
            )
            band_end = (
                (passband_response[band]["wavelength"].iloc[-1] * u.um).to(u.AA).value
            )
            # Add geocoronal emission line [O II] 2471A to the relevant passband
            if (GEOCORONAL_WAVELENGTH.value >= band_start) & (
                GEOCORONAL_WAVELENGTH.value <= band_end
            ):
                throughput_interp = interp1d(
                    passband_response[band]["wavelength"].values,
                    throughput,
                    kind=interp_kind,
                    bounds_error=False,
                    fill_value=np.nan,
                )
                geo_throughput = throughput_interp(GEOCORONAL_WAVELENGTH.to(u.um).value)
                # sky_background_e_rate[band] += (
                #     geo_background * geo_throughput * GEOCORONAL_LINEWIDTH.value
                # )  # electron/s
                # REVIEW: don't need to multiply by linewidth?
                sky_background_e_rate[band] += (
                    geo_background * geo_throughput
                )  # electron/s
    #
    # Generate signal from source
    #
    source_e_rate = dict.fromkeys(passbands)  # electron/s
    for band in source_e_rate:
        source_e_rate[band] = mag_to_flux(source_AB_mags, zpt=phot_zpts[band].values)[0]
    print(source_e_rate)
    # #
    # # Calculate red leak
    # #
    # redleaks = dict.fromkeys(passbands)
    # for band in redleaks:
    #     redleaks[band] = _calc_red_leak(passband_response[band]["wavelengths"], passband_response[band]["throughput"], )
    #
    # Calculate desired quantity
    #
    if snr is not None:
        # Find time required to reach target SNR
        t_lim = dict.fromkeys(passbands)  # s
        for band in t_lim:
            t_lim[band] = _calc_t_from_snr(
                snr=snr,
                signal=source_e_rate[band],
                npix=npix,
                totskynoise=sky_background_e_rate[band].values,
                readnoise=readnoise,
                darkcurrent=darkcurrent,
                redleak=redleak,
                nread=nread,
            )
        return t_lim
    else:
        # Find SNR reached in given time
        snr_reached = dict.fromkeys(passbands)  # dimensionless
        for band in snr_reached:
            snr_reached = _calc_snr_from_t(
                t=t,
                signal=source_e_rate[band],
                npix=npix,
                totskynoise=sky_background_e_rate[band].values,
                readnoise=readnoise,
                darkcurrent=darkcurrent,
                redleak=redleak,
                nread=nread,
            )
        return snr_reached
