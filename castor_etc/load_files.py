"""
load_files.py

Functions to load data from files.

Isaac Cheng - 2022
"""

import astropy.units as u
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d

from . import parameters as params
from .filepaths import DATAPATH


def load_passbands(
    filters="all", limits=None, resolution=None, filepaths=None, kind="linear"
):
    """
    Loads the passband throughput curves.

    Parameters
    ----------
      filters :: 1D list of strings or "all"
        The list of filters to load; valid filters are:
          - "uv"
          - "u"
          - "g"
        If "all" is given, all filters are loaded.

      limits :: dict of lists or None
        Dictionary containing key-value pairs of the filter and the wavelength where the
        value is a list containing the lower and upper wavelength limits for that passband
        (inclusive). If the values are lists of floats, the units are assumed to be
        micrometres. For example:
        {"uv": [0.15, 0.30], "u": [0.30, 0.40], "g": [0.40, 0.55]} or
        {"uv": [150, 300] * u.nm, "g": [400, 550] * u.nm}.
        If None, use the passband limits (from the parameters file) for each of the
        specified filters.

      resolution :: scalar or `astropy.Quantity` or None
        The linear interpolation resolution of the passband curves. If a scalar, the unit
        is assumed to be in micrometres. If None, use the native resolution of the
        passband files.

      filepaths :: 1D list of strings
        The absolute paths to the files containing the passband throughput curves. Each
        file should contain data for one passband; lines starting with a hash (#) will be
        ignored. If filepaths is None, use the default passband filepath for each filter.

        TODO: Explain format of each file (later).

      kind :: str
        The type of interpolation to use. See `scipy.interpolate.interp1d` for options.

    Returns
    -------
      passbands :: dict of `pandas.DataFrame`
        Dictionary containing the passband throughput curves. Each throughput curve is a
        pandas DataFrame and the wavelengths are in micrometres.

    Examples
    --------
    ```
    passbands = load_passbands(filters="all")  # load all filters
    print(passbands)  # view all data for all filters
    print(passbands["uv"])  # view all data for the "uv" filter
    print(passbands["u"]["wavelength"])  # view the wavelength values of the "u" filter
    print(passbands["g"]["throughput"])  # view the throughput values of the "g" filter
    print(passbands["uv"].iloc[0])  # view the 1st row's data for the "uv" filter
    print(passbands["uv"].iloc[0:3])  # view data from rows 0 to 2 for the "uv" filter
    ```
    """
    #
    # Check inputs
    #
    if isinstance(filters, str):
        filters = [filters]
    if np.ndim(filters) != 1:
        raise ValueError("filters must be a 1D list of strings or 'all'")
    if filters == ["all"]:
        filters = params.PASSBANDS
    else:
        if any(band not in params.PASSBANDS for band in filters):
            raise ValueError(f"Invalid filters. Valid filters are: {params.PASSBANDS}")
    #
    if limits is None:
        limits = {band: params.PASSBAND_LIMITS[band].to(u.um).value for band in filters}
    else:
        limits = limits.copy()
        for key in limits.keys():
            limits[key] = list(limits[key])
            for i, lim in enumerate(limits[key]):
                if isinstance(lim, u.Quantity):
                    limits[key][i] = lim.to(u.um).value
    if any(band not in limits.keys() for band in filters) or any(
        key not in filters for key in limits.keys()
    ):
        raise ValueError("filters and limits must be consistent")
    #
    if isinstance(resolution, u.Quantity):
        resolution = resolution.to(u.um).value
    #
    if filepaths is None:
        filepaths = [DATAPATH + f"passbands/passband_castor.{band}" for band in filters]
    else:
        if isinstance(filepaths, str):
            filepaths = [filepaths]
        if np.ndim(filepaths) != 1:
            raise ValueError("filepaths must be a 1D list of strings or None")
        elif len(filepaths) != len(filters):
            raise ValueError("filepaths and filters must be consistent")
    #
    # Load passband files
    #
    passbands = dict.fromkeys(filters)  # initialize empty dictionary
    for filepath, band in zip(filepaths, passbands):
        passbands[band] = pd.read_csv(
            filepath,
            sep=" +",
            header=None,
            comment="#",
            engine="python",
        )  # sep=" +" is Python regex to match a variable number of spaces
    #
    # Trim passband data to given limits
    #
    if resolution is not None:
        # Interpolate passband data to given resolution
        for band in passbands:
            interp = interp1d(
                passbands[band][0].values,
                passbands[band][1].values,
                kind=kind,
                bounds_error=False,
                fill_value=np.nan,
            )
            wavelengths = np.arange(
                limits[band][0], limits[band][1] + 0.5 * resolution, resolution
            )
            passbands[band] = pd.DataFrame(
                {"wavelength": wavelengths, "throughput": interp(wavelengths)},
            )
    else:
        # Do not interpolate passband data
        for band in passbands:
            is_good = (passbands[band][0] >= limits[band][0]).values & (
                passbands[band][0] <= limits[band][1]
            ).values
            passbands[band] = pd.DataFrame(
                {
                    "wavelength": (passbands[band][0].values)[is_good],
                    "throughput": (passbands[band][1].values)[is_good],
                }
            )
    return passbands


def load_sky_background(
    resolution=None,
    limits=None,
    filepath=DATAPATH + "background/high_sky_background.txt",
    kind="linear",
):
    """
    Loads the sky background noise.

    Parameters
    ----------
      resolution :: float or `astropy.Quantity` or None
        The linear interpolation resolution of the data. If a float, value must be in
        angstroms. If None, use the native resolution of the sky background data.

      limits :: list of 2 floats or list of 2 `astropy.Quantity` or None
        If resolution is None, limit specifies the minimum and maximum wavelengths of the
        returned background data, inclusive. If resolution is not None, limit specifies
        the start and end wavelengths of the interpolated background data, inclusive. If
        the list elements are floats, the values are assumed to be in angstroms. If this
        parameter is None, use the min and max wavelengths from the background file as the
        interpolation limits.

      filepath :: str
        The absolute path to the file containing the sky background data. Data should be
        separated by spaces and lines starting with a hash (#) will be ignored.

        TODO: Explain format of file (later)

      kind :: str
        The type of interpolation to use. See `scipy.interpolate.interp1d` for options.

    Returns
    -------
      background :: `pandas.DataFrame`
        DataFrame containing the sky background noise. Wavelengths are in angstroms and
        the background values are in erg/cm^2/s/A/arcsec^2.
    """
    #
    # Check inputs
    #
    if isinstance(resolution, u.Quantity):
        resolution = resolution.to(u.AA).value
    if limits is not None:
        for i, lim in enumerate(limits):
            if isinstance(lim, u.Quantity):
                limits[i] = lim.to(u.AA).value
    #
    # Load background data
    #
    background = pd.read_csv(
        filepath,
        sep=" ",
        header=0,
        comment="#",
    )
    #
    # Interpolate passband data to given resolution
    #
    if resolution is not None:
        if limits is None:
            limits = [
                np.nanmin(background["wavelength"]),
                np.nanmax(background["wavelength"]),
            ]
        wavelengths = np.arange(limits[0], limits[1] + 0.5 * resolution, resolution)
        interp_background = pd.DataFrame({"wavelength": wavelengths})
        for col in background.columns:
            if col == "wavelength":
                continue  # no need to interpolate the wavelength column
            interp = interp1d(
                background["wavelength"].values,
                background[col].values,
                kind=kind,
                bounds_error=False,
                fill_value=np.nan,
            )
            # Add column to DataFrame
            interp_background[col] = interp(wavelengths)
        background = interp_background
    elif limits is not None:
        background = background[
            (background["wavelength"] >= limits[0])
            & (background["wavelength"] <= limits[1])
        ]
    return background
