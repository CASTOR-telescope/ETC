"""
Utilities to characterize the telescope parameters.
"""

import warnings
from copy import deepcopy
from numbers import Number

import astropy.units as u
import numpy as np
import pandas as pd
from scipy.integrate import simpson
from scipy.interpolate import interp1d

from . import parameters as params

# TODO: add PSF


class Telescope:
    """
    Telescope class containing the parameters of the imaging system.

    TODO: add PSF (Moffat and Gaussian)
    """

    # def __init__(
    #     self,
    #     passbands=params.PASSBANDS,
    #     passband_limits=params.PASSBAND_LIMITS,
    #     passband_resolution=params.PASSBAND_RESOLUTION,
    #     passband_pivots=params.PASSBAND_PIVOTS,  # REVIEW: maybe remove to force calculate
    #     # passband_response_filepaths=params.PASSBAND_FILEPATHS,
    #     # passband_response_fileunits=params.PASSBAND_FILEUNITS,
    #     passband_response_filepaths=None,
    #     passband_response_fileunits=None,
    #     phot_zpts=params.PHOT_ZPTS,
    #     fwhm=params.FWHM,
    #     px_scale=params.PX_SCALE,
    #     ifov_dimen=params.IFOV_DIMEN,
    #     mp=params.MP,
    #     mirror_diameter=params.MIRROR_DIAMETER,
    #     dark_current=params.DARK_CURRENT,
    #     bias=params.BIAS,
    #     read_noise=params.READ_NOISE,
    #     gain=params.GAIN,
    #     redleak_thresholds=params.REDLEAK_THRESHOLDS,
    #     show_warnings=True,
    # ):
    #     """
    #     TODO: if passband curves provided, calculate phot_zpts and passband_pivots

    #     Create a Telescope instance.

    #     Parameters
    #     ----------
    #       passbands :: list of str
    #         The telescope filter (i.e., passband) names. CASTOR passbands are "uv", "u",
    #         and "g", but other passband names are allowed.

    #       passband_limits :: dict of 2-element `astropy.Quantity` arrays
    #         The lower and upper wavelength limits for each passband. Each element should
    #         be an `astropy.Quantity` array of length 2 (e.g., `{"uv": [0.150, 0.300] *
    #         u.um, "u": [0.300, 0.400] * u.um, "g": [0.400, 0.550] * u.um]}`).

    #       passband_resolution :: `astropy.Quantity`
    #         The desired linear interpolation resolution of the passband response curves.

    #       phot_zpts :: dict of float
    #         The photometric zeropoints for each passband.
    #         TODO: auto-calculate phot_zpts

    #       passband_pivots :: dict of `astropy.Quantity` or None
    #         The pivot wavelengths for each passband. If None, calculate the passband
    #         pivots based on the given passband limits and passband response files.

    #       passband_response_filepaths :: dict of str
    #         The absolute paths to the files containing the passband response curves. The
    #         files should be in ASCII format with the first column containing the
    #         wavelengths in `passband_response_fileunits` and the second column containing
    #         the passband response; the columns should be separated by a space. Lines
    #         starting with a hash (#) will be ignored. Note that the passband response
    #         values should account for quantum efficiency.

    #       passband_response_fileunits :: dict of `astropy.Quantity` lengths
    #         The units of the wavelength columns in the passband response files.

    #       redleak_thresholds :: dict of `astropy.Quantity` wavelengths
    #         The redleak thresholds for each passband. Flux longward of the threshold is
    #         considered red leak.

    #       show_warnings :: bool
    #         If True, print a warning when passband_limits, passband_pivots, or phot_zpts
    #         contain keys that are not in passbands.

    #     Attributes
    #     ----------
    #       passbands ::

    #       passband_limits ::

    #       passband_tot_limits ::

    #       passband_resolution ::

    #       passband_pivots ::

    #       phot_zpts ::

    #       fwhm ::

    #       px_scale ::

    #       px_area ::

    #       ifov_dimen ::

    #       ifov_area ::

    #       mp ::

    #       mirror_diameter ::

    #       mirror_area ::

    #       dark_current ::

    #       bias ::

    #       read_noise ::

    #       gain ::

    #       redleak_thresholds ::

    #     Returns
    #     -------
    #       `Telescope` instance

    #       TODO: finish docstring
    #     """

    #     def _check_dict(d, d_name, value_types, value_types_str):
    #         """
    #         Check that the input dictionary has valid keys for passbands. The passbands
    #         must all have key-value pairs in the input dictionary but the keys in the
    #         input dictionary do not necessarily have to be in passbands.

    #         Parameters
    #         ----------
    #           d :: dict
    #             The input dictionary.

    #           d_name :: str
    #             The name of the input dictionary (i.e., the variable name).

    #           value_types :: type or tuple of types
    #             The expected type of the values in the input dictionary. A TypeError will
    #             be raised if a key in the input dictionary is not any of the types in
    #             value_types. If value_types is `u.Quantity`, then the values will also be
    #             checked if they are length units (i.e., u.AA).

    #           value_types_str :: str
    #             The expected type of the values in the input dictionary. These will be
    #             printed in the case of a TypeError.

    #         Returns
    #         -------
    #           d :: dict
    #             The input dictionary with any key that is not in `passbands` removed.
    #         """
    #         d = d.copy()
    #         warning_printed = False  # track if warning has already been printed
    #         #
    #         if any(band not in d.keys() for band in passbands):
    #             raise KeyError(f"Not all keys in `passbands` are in `{d_name}`.")
    #         if value_types is u.Quantity:
    #             for value in d.values():
    #                 try:
    #                     _ = value.to(u.AA)
    #                 except Exception:
    #                     raise TypeError(f"`{d_name}` must be dict of {value_types_str}")
    #         else:
    #             if any(not isinstance(value, value_types) for value in d.values()):
    #                 raise TypeError(f"`{d_name}` must be dict of {value_types_str}")
    #         #
    #         for band in list(d):
    #             if band not in passbands:
    #                 if show_warnings and not warning_printed:
    #                     print(
    #                         f"WARNING: some key(s) in `{d_name}` are not in `passbands`! "
    #                         + "These will be removed."
    #                     )
    #                     warning_printed = True
    #                 del d[band]
    #         #
    #         return d

    #     #
    #     # Check inputs
    #     #
    #     if isinstance(passbands, str):
    #         passbands = [passbands]
    #     if (np.ndim(passbands) != 1) or (
    #         not all(isinstance(band, str) for band in passbands)
    #     ):
    #         raise TypeError("passbands must be a 1D list of strings")
    #     # if any(band not in params.PASSBANDS for band in passbands):
    #     #     raise ValueError(
    #     #         f"Invalid passbands. Valid passbands are: {params.PASSBANDS}"
    #     #     )
    #     try:
    #         _ = passband_resolution.to(u.AA)
    #     except Exception:
    #         raise TypeError(
    #             "passband_resolution must be an `astropy.Quantity` length "
    #             + "(e.g., u.nm or u.AA)"
    #         )
    #     passband_limits = _check_dict(
    #         passband_limits,
    #         "passband_limits",
    #         u.Quantity,
    #         "`astropy.Quantity` lengths (e.g., u.nm)",
    #     )
    #     if passband_pivots is not None:
    #         if (
    #             passband_response_filepaths is not None
    #             and passband_response_fileunits is not None
    #         ):
    #             print(
    #                 "INFO: both passband_response parameters and passband_pivots "
    #                 + "provided. Will ignore passband_pivots."
    #             )
    #             passband_pivots = None
    #         else:
    #             passband_pivots = _check_dict(
    #                 passband_pivots,
    #                 "passband_pivots",
    #                 u.Quantity,
    #                 "`astropy.Quantity` lengths (e.g., u.nm)",
    #             )
    #     else:
    #         if passband_response_filepaths is None or passband_response_fileunits is None:
    #             raise ValueError(
    #                 "Either passband_pivots must be given or both "
    #                 + "passband_response_filepaths and passband_response_fileunits "
    #                 + "must be given."
    #             )
    #     redleak_thresholds = _check_dict(
    #         redleak_thresholds,
    #         "redleak_thresholds",
    #         u.Quantity,
    #         "`astropy.Quantity` lengths (e.g., u.AA)",
    #     )
    #     if passband_response_filepaths is not None:
    #         passband_response_filepaths = _check_dict(
    #             passband_response_filepaths,
    #             "passband_response_filepaths",
    #             str,
    #             "absolute file path strings",
    #         )
    #     if passband_response_fileunits is not None:
    #         passband_response_fileunits = _check_dict(
    #             passband_response_fileunits,
    #             "passband_response_fileunits",
    #             u.Quantity,
    #             "`astropy.Quantity` lengths (e.g., u.um)",
    #         )
    #     if phot_zpts is not None:
    #         if (
    #             passband_response_filepaths is not None
    #             and passband_response_fileunits is not None
    #         ):
    #             print(
    #                 "INFO: both passband_response parameters and phot_zpts provided. "
    #                 + "Will ignore phot_zpts."
    #             )
    #             phot_zpts = None
    #         else:
    #             phot_zpts = _check_dict(
    #                 phot_zpts,
    #                 "phot_zpts",
    #                 Number,
    #                 "int or float",
    #             )
    #     else:
    #         if passband_response_filepaths is None or passband_response_fileunits is None:
    #             raise ValueError(
    #                 "Either phot_zpts must be given or both passband_response_filepaths "
    #                 + "and passband_response_fileunits must be given."
    #             )
    #     angles = [fwhm, px_scale, ifov_dimen]
    #     angles_str = ["fwhm", "px_scale", "ifov_dimen"]
    #     for angle, angle_str in zip(angles, angles_str):
    #         try:
    #             _ = angle.to(u.arcsec)
    #         except Exception:
    #             raise TypeError(
    #                 f"{angle_str} must be an `astropy.Quantity` angle (e.g., u.arcsec)"
    #             )
    #     try:
    #         _ = mirror_diameter.to(u.cm)
    #     except Exception:
    #         raise TypeError(
    #             "mirror_diameter must be an `astropy.Quantity` length (e.g., u.cm)"
    #         )
    #     scalars = [mp, dark_current, bias, read_noise, gain]
    #     scalars_str = ["mp", "dark_current", "bias", "read_noise", "gain"]
    #     for scalar, scalar_str in zip(scalars, scalars_str):
    #         if not isinstance(scalar, Number):
    #             raise TypeError(f"{scalar_str} must be an int or float")
    #     #
    #     # Assign attributes
    #     #
    #     self.passbands = passbands
    #     self.passband_limits = passband_limits
    #     self.passband_tot_limits = [
    #         min(passband_limits.values(), key=lambda x: x[0])[0],
    #         max(passband_limits.values(), key=lambda x: x[1])[1],
    #     ]  # total wavelength range spanned by the passbands
    #     self.passband_resolution = passband_resolution

    #     if phot_zpts is None or passband_pivots is None:
    #         passband_curves = Telescope.load_passbands(
    #             passband_response_filepaths,
    #             passband_limits,
    #             passband_response_fileunits,
    #             resolution=passband_resolution,
    #         )

    #         if phot_zpts is None:
    #             # TODO: calculate phot_zpts from passband_curves
    #             raise NotImplementedError("phot_zpts calculation not implemented yet")

    #         if passband_pivots is None:
    #             passband_pivots = dict.fromkeys(passbands)
    #             for band in passbands:
    #                 passband_pivots[band] = (
    #                     Telescope.calc_pivot_wavelength(
    #                         passband_curves[band]["wavelength"].values,
    #                         passband_curves[band]["response"].values,
    #                     )
    #                     * passband_response_fileunits[band]
    #                 )
    #     self.passband_pivots = passband_pivots

    #     # Photometric zeropoints for the different filters
    #     self.phot_zpts = phot_zpts

    #     # The PSF full-width at half-maximum
    #     self.fwhm = fwhm

    #     # The angular dimension covered by each pixel
    #     self.px_scale = px_scale
    #     # The angular area of each pixel
    #     self.px_area = px_scale * px_scale

    #     # Instantaneous field of view
    #     self.ifov_dimen = ifov_dimen
    #     self.ifov_area = ifov_dimen[0] * ifov_dimen[1]

    #     # Number of pixels of CCD (x 1 million)
    #     self.mp = mp

    #     self.mirror_diameter = mirror_diameter
    #     self.mirror_area = np.pi * 0.25 * mirror_diameter * mirror_diameter

    #     self.dark_current = dark_current
    #     self.bias = bias
    #     self.read_noise = read_noise
    #     self.gain = gain

    #     self.redleak_thresholds = redleak_thresholds

    def __init__(
        self,
        passbands=params.PASSBANDS,
        passband_limits=params.PASSBAND_LIMITS,
        passband_resolution=params.PASSBAND_RESOLUTION,
        passband_pivots=params.PASSBAND_PIVOTS,  # REVIEW: maybe remove to force calculate
        passband_response_filepaths=params.PASSBAND_FILEPATHS,
        passband_response_fileunits=params.PASSBAND_FILEUNITS,
        phot_zpts=params.PHOT_ZPTS,
        fwhm=params.FWHM,
        px_scale=params.PX_SCALE,
        ifov_dimen=params.IFOV_DIMEN,
        mp=params.MP,
        mirror_diameter=params.MIRROR_DIAMETER,
        dark_current=params.DARK_CURRENT,
        bias=params.BIAS,
        read_noise=params.READ_NOISE,
        gain=params.GAIN,
        redleak_thresholds=params.REDLEAK_THRESHOLDS,
        show_warnings=True,
    ):
        """
        TODO: if passband_pivots or phot_zpts provided, override auto-calculation from passband_response

        Create a Telescope instance.

        Parameters
        ----------
          passbands :: list of str
            The telescope filter (i.e., passband) names. CASTOR passbands are "uv", "u",
            and "g", but other passband names are allowed.

          passband_limits :: dict of 2-element `astropy.Quantity` arrays
            The lower and upper wavelength limits for each passband. Each element should
            be an `astropy.Quantity` array of length 2 (e.g., `{"uv": [0.150, 0.300] *
            u.um, "u": [0.300, 0.400] * u.um, "g": [0.400, 0.550] * u.um]}`).

          passband_resolution :: `astropy.Quantity`
            The desired linear interpolation resolution of the passband response curves.

          phot_zpts :: dict of float
            The photometric zeropoints for each passband.
            TODO: auto-calculate phot_zpts

          passband_pivots :: dict of `astropy.Quantity` or None
            The pivot wavelengths for each passband. If None, calculate the passband
            pivots based on the given passband limits and passband response files.

          passband_response_filepaths :: dict of str
            The absolute paths to the files containing the passband response curves. The
            files should be in ASCII format with the first column containing the
            wavelengths in `passband_response_fileunits` and the second column containing
            the passband response; the columns should be separated by a space. Lines
            starting with a hash (#) will be ignored. Note that the passband response
            values should account for quantum efficiency.

          passband_response_fileunits :: dict of `astropy.Quantity` lengths
            The units of the wavelength columns in the passband response files.

          redleak_thresholds :: dict of `astropy.Quantity` wavelengths
            The redleak thresholds for each passband. Flux longward of the threshold is
            considered red leak.

          show_warnings :: bool
            If True, print a warning when passband_limits, passband_pivots, or phot_zpts
            contain keys that are not in passbands.

        Attributes
        ----------
          passbands ::

          passband_limits ::

          passband_tot_limits ::

          passband_resolution ::

          passband_pivots ::

          passband_curves ::

          phot_zpts ::

          fwhm ::

          px_scale ::

          px_area ::

          ifov_dimen ::

          ifov_area ::

          mp ::

          mirror_diameter ::

          mirror_area ::

          dark_current ::

          bias ::

          read_noise ::

          gain ::

          redleak_thresholds ::

        Returns
        -------
          `Telescope` instance

          TODO: finish docstring
        """

        def _check_dict(d, d_name, value_types, value_types_str):
            """
            Check that the input dictionary has valid keys for passbands. The passbands
            must all have key-value pairs in the input dictionary but the keys in the
            input dictionary do not necessarily have to be in passbands.

            Parameters
            ----------
              d :: dict
                The input dictionary.

              d_name :: str
                The name of the input dictionary (i.e., the variable name).

              value_types :: type or tuple of types
                The expected type of the values in the input dictionary. A TypeError will
                be raised if a key in the input dictionary is not any of the types in
                value_types. If value_types is `u.Quantity`, then the values will also be
                checked if they are length units (i.e., u.AA).

              value_types_str :: str
                The expected type of the values in the input dictionary. These will be
                printed in the case of a TypeError.

            Returns
            -------
              d :: dict
                The input dictionary with any key that is not in `passbands` removed.
            """
            d = d.copy()
            warning_printed = False  # track if warning has already been printed
            #
            if any(band not in d.keys() for band in passbands):
                raise KeyError(f"Not all keys in `passbands` are in `{d_name}`.")
            if value_types is u.Quantity:
                for value in d.values():
                    try:
                        _ = value.to(u.AA)
                    except Exception:
                        raise TypeError(f"`{d_name}` must be dict of {value_types_str}")
            else:
                if any(not isinstance(value, value_types) for value in d.values()):
                    raise TypeError(f"`{d_name}` must be dict of {value_types_str}")
            #
            for band in list(d):
                if band not in passbands:
                    if show_warnings and not warning_printed:
                        warning_printed = True
                        warnings.warn(
                            f"Some key(s) in `{d_name}` are not in `passbands`! "
                            + "These will be removed.",
                            UserWarning,
                        )
                    del d[band]
            #
            return d

        #
        # Check inputs
        #
        if isinstance(passbands, str):
            passbands = [passbands]
        if (np.ndim(passbands) != 1) or (
            not all(isinstance(band, str) for band in passbands)
        ):
            raise TypeError("passbands must be a 1D list of strings")
        # if any(band not in params.PASSBANDS for band in passbands):
        #     raise ValueError(
        #         f"Invalid passbands. Valid passbands are: {params.PASSBANDS}"
        #     )
        try:
            _ = passband_resolution.to(u.AA)
        except Exception:
            raise TypeError(
                "passband_resolution must be an `astropy.Quantity` length "
                + "(e.g., u.nm or u.AA)"
            )
        passband_limits = _check_dict(
            passband_limits,
            "passband_limits",
            u.Quantity,
            "`astropy.Quantity` lengths (e.g., u.nm)",
        )
        passband_response_filepaths = _check_dict(
            passband_response_filepaths,
            "passband_response_filepaths",
            str,
            "absolute file path strings",
        )
        passband_response_fileunits = _check_dict(
            passband_response_fileunits,
            "passband_response_fileunits",
            u.Quantity,
            "`astropy.Quantity` lengths (e.g., u.um)",
        )
        if passband_pivots is not None:
            passband_pivots = _check_dict(
                passband_pivots,
                "passband_pivots",
                u.Quantity,
                "`astropy.Quantity` lengths (e.g., u.nm)",
            )
            print(
                "INFO: passband_pivots and passband response files both provided. "
                + "Will use user-supplied passband_pivots rather than calculating pivots."
            )
        redleak_thresholds = _check_dict(
            redleak_thresholds,
            "redleak_thresholds",
            u.Quantity,
            "`astropy.Quantity` lengths (e.g., u.AA)",
        )
        if phot_zpts is not None:
            phot_zpts = _check_dict(
                phot_zpts,
                "phot_zpts",
                Number,
                "int or float",
            )
            print(
                "INFO: phot_zpts and passband response files both provided. "
                + "Will use user-supplied phot_zpts rather than calculating zero-points."
            )
        angles = [fwhm, px_scale, ifov_dimen]
        angles_str = ["fwhm", "px_scale", "ifov_dimen"]
        for angle, angle_str in zip(angles, angles_str):
            try:
                _ = angle.to(u.arcsec)
            except Exception:
                raise TypeError(
                    f"{angle_str} must be an `astropy.Quantity` angle (e.g., u.arcsec)"
                )
        try:
            _ = mirror_diameter.to(u.cm)
        except Exception:
            raise TypeError(
                "mirror_diameter must be an `astropy.Quantity` length (e.g., u.cm)"
            )
        scalars = [mp, dark_current, bias, read_noise, gain]
        scalars_str = ["mp", "dark_current", "bias", "read_noise", "gain"]
        for scalar, scalar_str in zip(scalars, scalars_str):
            if not isinstance(scalar, Number):
                raise TypeError(f"{scalar_str} must be an int or float")
        #
        # Assign attributes
        #
        self.passbands = passbands
        self.passband_limits = passband_limits
        self.passband_tot_limits = [
            min(passband_limits.values(), key=lambda x: x[0])[0],
            max(passband_limits.values(), key=lambda x: x[1])[1],
        ]  # total wavelength range spanned by the passbands
        self.passband_resolution = passband_resolution

        self.passband_curves = Telescope.load_passbands(
            passband_response_filepaths,
            passband_limits,
            passband_response_fileunits,
            resolution=passband_resolution,
        )

        if phot_zpts is None:
            # TODO: calculate phot_zpts from passband_curves
            raise NotImplementedError("phot_zpts calculation not implemented yet")

        if passband_pivots is None:
            passband_pivots = dict.fromkeys(passbands)
            for band in passbands:
                passband_pivots[band] = (
                    Telescope.calc_pivot_wavelength(
                        self.passband_curves[band]["wavelength"].values,
                        self.passband_curves[band]["response"].values,
                    )
                    * passband_response_fileunits[band]
                )
        self.passband_pivots = passband_pivots

        # Photometric zeropoints for the different filters
        self.phot_zpts = phot_zpts

        # The PSF full-width at half-maximum
        self.fwhm = fwhm

        # The angular dimension covered by each pixel
        self.px_scale = px_scale
        # The angular area of each pixel
        self.px_area = px_scale * px_scale

        # Instantaneous field of view
        self.ifov_dimen = ifov_dimen
        self.ifov_area = ifov_dimen[0] * ifov_dimen[1]

        # Number of pixels of CCD (x 1 million)
        self.mp = mp

        self.mirror_diameter = mirror_diameter
        self.mirror_area = np.pi * 0.25 * mirror_diameter * mirror_diameter

        self.dark_current = dark_current
        self.bias = bias
        self.read_noise = read_noise
        self.gain = gain

        self.redleak_thresholds = redleak_thresholds

    def copy(self):
        """
        Convenience method for creating a deep copy of the Telescope object.

        Parameters
        ----------
          None

        Returns
        -------
          Telescope_copy :: `Telescope` object
            The deep copy of the `Telescope` object.
        """
        return deepcopy(self)

    @staticmethod
    def calc_pivot_wavelength(wavelengths, response, response_func="QE"):
        """
        Calculate the pivot wavelength of a passband using Eq. (A11) from Tokunaga & Vacca
        (2005) <https://ui.adsabs.harvard.edu/abs/2005PASP..117..421T/abstract>. This
        function uses Simpson's rule to approximate the integration.

        INFO: Equal energy response function vs. QE response function?

        Parameters
        ----------
          wavelengths :: array of floats
            The wavelengths in the passband

          response :: array of floats
            The passband response at each wavelength. This response should include the
            quantum efficiency of the detector.

          reponse_func :: "EE" or "QE"
            The convention used for the response function. "EE" denotes an the
            equal energy response function while "QE" denotes the quantum efficiency
            response function.

        Returns
        -------
          pivot_wavelength :: float
            The pivot wavelength of the passband (in the same units as the input
            wavelengths).
        """
        if response_func == "EE":
            numer = simpson(y=wavelengths * response, x=wavelengths, even="avg")
            denom = simpson(y=response / wavelengths, x=wavelengths, even="avg")
        elif response_func == "QE":
            numer = simpson(y=response, x=wavelengths, even="avg")
            denom = simpson(
                y=response / (wavelengths * wavelengths), x=wavelengths, even="avg"
            )
        return np.sqrt(numer / denom)

    @staticmethod
    def load_passbands(
        filepaths, limits, file_units, resolution=1 << u.nm, interp_kind="linear"
    ):
        """
        Load the passband response curves.

        Parameters
        ----------
          filepaths :: dict of str
            The absolute paths to the files containing the passband response curves. The
            files should be in ASCII format with the first column containing the
            wavelengths in `file_unit` and the second column containing the passband
            response; the columns should be separated by a space. Lines starting with a
            hash (#) will be ignored. Note that the passband response values should
            account for quantum efficiency.

          limits :: dict of 2-element `astropy.Quantity` length arrays
            Dictionary containing key-value pairs of the filter and the wavelength where
            the value is an `astropy.Quantity` array containing the lower and upper
            wavelength limits for that passband (inclusive). For example:
            `{"uv": [150, 300] * u.nm, "u": [300, 400] * u.nm, "g": [400, 550] * u.nm}.`

          file_units :: dict of `astropy.Quantity` lengths
            The units of the wavelength data in the passband response curve files.

          resolution :: `astropy.Quantity` or None
            The linear interpolation resolution of the passband curves; always recommended
            to set to not None. If None, use the native resolution of the passband files.
            There are two scenarios that can occur if the response curve file does not
            have data at the requested wavelength (i.e., `limits` go beyond the
            wavelengths provided in the file):
              1. if resolution is not None, the response curve value will be assigned NaN
              2. if resolution is None, the response curve will simply omit these
                 wavelengths from the final DataFrame.

          interp_kind :: str
            The type of interpolation to use when resolution is not None. See
            `scipy.interpolate.interp1d` for options.

        Returns
        -------
          curves :: dict of dicts
            Dictionary containing the passband response curves. Each response curve is a
            dictionary containing the wavelengths in an `astropy.Quantity` array (same
            units as `file_units`) and the response in a float array.

        Examples
        --------
        ```python
        print(passbands)                     # view all data for all passbands
        print(passbands["uv"])               # view all data for the UV-passband
        print(passbands["u"]["wavelength"])  # view the u-band wavelengths
        print(passbands["g"]["response"])    # view the g-band response
        ```
        """
        #
        # Check inputs
        #
        if (
            any(band not in limits.keys() for band in filepaths.keys())
            or any(band not in filepaths.keys() for band in limits.keys())
            or any(band not in file_units.keys() for band in filepaths.keys())
            or any(band not in filepaths.keys() for band in file_units.keys())
        ):
            raise ValueError(
                "filepaths, limits, and file_units must have consistent keys"
            )
        limits = limits.copy()
        for band in limits:
            limits[band] = (limits[band]).to(file_units[band]).value
        #
        if resolution is not None:
            resolutions = dict.fromkeys(filepaths)
            for band in resolutions:
                resolutions[band] = resolution.to(file_units[band]).value
        else:
            resolutions = None
        #
        # Load passband files
        #
        curves = dict.fromkeys(filepaths)  # initialize empty dictionary
        for band in curves:
            curves[band] = pd.read_csv(
                filepaths[band],
                sep=" +",
                header=None,
                comment="#",
                engine="python",
            )  # sep=" +" is Python regex to match a variable number of spaces
        #
        # Trim passband data to given limits
        #
        if resolutions is not None:
            # Interpolate passband data to given resolution
            for band in curves:
                interp = interp1d(
                    curves[band][0].values,
                    curves[band][1].values,
                    kind=interp_kind,
                    bounds_error=False,
                    fill_value=np.nan,
                )
                wavelengths = np.arange(
                    limits[band][0],
                    limits[band][1] + 0.5 * resolutions[band],
                    resolutions[band],
                )
                curves[band] = {
                    "wavelength": wavelengths * file_units[band],
                    "response": interp(wavelengths),
                }
        else:
            # Do not interpolate passband data
            for band in curves:
                is_good = (curves[band][0] >= limits[band][0]).values & (
                    curves[band][0] <= limits[band][1]
                ).values
                curves[band] = {
                    "wavelength": (curves[band][0].values)[is_good] * file_units[band],
                    "response": (curves[band][1].values)[is_good],
                }

        return curves
