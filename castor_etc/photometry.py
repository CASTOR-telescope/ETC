"""
Utilities for photometric calculations.
"""

import warnings
from copy import deepcopy
from numbers import Number

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from photutils.aperture import EllipticalAnnulus, EllipticalAperture, RectangularAperture
from scipy.integrate import simpson
from scipy.interpolate import interp1d

from .background import Background
from .conversions import calc_photon_energy, convert_electron_flux_mag, mag_to_flux
from .sources import ExtendedSource, PointSource, Source
from .telescope import Telescope

# TODO: deprecate self.npix
# TODO: convolve with PSF
# TODO: dithering? Maybe just a simple dither_factor=0.75 parameter to simulate oversampling?


class Photometry:
    """
    Photometry class.

    Attributes
    ----------
      aperture :: int or float or `astropy.Quantity`
        The aperture area. If an int or float, aperture is in pixel units.
    """

    def __init__(self, TelescopeObj, SourceObj, BackgroundObj):
        """
        Initialize class for photometry calculations.

        Parameters
        ----------
          TelescopeObj :: `castor_etc.Telescope` object
            The `castor_etc.Telescope` object containing the telescope parameters.

          SourceObj :: `castor_etc.Source` object
            The `castor_etc.Source` object contaning the target source parameters. The
            source's spectrum must be in units of photlam (photon/s/cm^2/A).

          BackgroundObj :: `castor_etc.Background` object
            The `castor_etc.Background` object containing the background parameters.
            Dictionary keys must match the TelescopeObj.passbands keys.

        Returns
        -------
          `Photometry` instance
        """
        #
        # Check inputs
        #
        if not isinstance(TelescopeObj, Telescope):
            raise TypeError("TelescopeObj must be a `castor_etc.Telescope` object")
        if not isinstance(SourceObj, Source):
            raise TypeError("SourceObj must be a `castor_etc.Source` object")
        if not isinstance(BackgroundObj, Background):
            raise TypeError("BackgroundObj must be a `castor_etc.Background` object")
        # TODO: check BackgroundObj dict keys
        #
        # Assign attributes
        #
        self.TelescopeObj = TelescopeObj
        self.SourceObj = SourceObj
        self.BackgroundObj = BackgroundObj

        self.aperture = None
        self.npix = None  # ! DEPRECATED number of pixels in aperture.

        self.aper_xs = None
        self.aper_ys = None
        self._aper_mask = None
        self._eff_npix = None
        self._aper_extent = None
        self.source_weights = None
        self.sky_background_weights = None
        self.dark_current_weights = None
        self.redleak_weights = None

        self._optimal_aperture_dimen = self._calc_optimal_aperture_dimen(
            SourceObj.angle_a, SourceObj.angle_b, print_info=False
        )  # only valid for point sources

    def copy(self):
        """
        Convenience method for creating a deep copy of the Photometry object.

        Parameters
        ----------
          None

        Returns
        -------
          Photometry_copy :: `Photometry` object
            The deep copy of the `Photometry` object.
        """
        return deepcopy(self)

    @staticmethod
    def angle_to_pix(angle, px_scale):
        """
        ! DEPRECATED ! (not used)
        Calculate the number of pixels subtended by an angle
        (e.g., u.arcsec NOT u.arcsec ** 2).

        Parameters
        ----------
          angle :: `astropy.Quantity` angle
            Linear angle (e.g., u.arcsec NOT u.arcsec ** 2).

          px_scale :: `astropy.Quantity` angle
            The angular size of a square pixel along one edge.

        Returns
        -------
          pix :: float
            The number of pixels lengths covered by the linear angle.
        """
        return (angle.to(u.arcsec) / px_scale.to(u.arcsec)).value

    def _assign_npix(self):
        """
        ! DEPRECATED
        Internal function. Calculate the number of pixels in the aperture.
        """
        self.npix = (
            self.aperture.to(u.arcsec ** 2) / self.TelescopeObj.px_area.to(u.arcsec ** 2)
        ).value

    def _create_aper_arrs(self, half_x, half_y, center):
        """
        Parameters
        ----------
          half_x, half_y :: floats
            The half-widths (in arcsec) of the aperture in the x- and y-directions,
            respectively. (e.g., semimajor/semiminor axes, half of a rectangle's
            length/width, etc.)

          center :: 2-element 1D array of floats
            The (x, y) center of the aperture relative to the center of the source in
            arcsec. Positive values means the source is displaced to the bottom/left
            relative to the aperture center.

        Attributes
        ----------
        REVIEW: use sparse arrays?

          TODO: update attribute docstring (no longer sparse)
          aper_xs :: (1 x N) sparse array of floats
            The aperture array containing the x-coordinates (in arcsec) of the aperture.
            The origin (0, 0) is defined to be at the center of the source. The element at
            the center of this array is the center of the aperture (i.e., at
            `[aper_xs.shape[0]//2, aper_xs.shape[1]//2]`). Together with `aper_ys`, this
            sparse array can easily broacast into a full (i.e., M x N) 2D array (e.g.,
            `arr = aper_xs + aper_ys` will produce a full (M x N) 2D array). The center of
            the aperture in the full 2D array is located at `[arr.shape[0]//2,
            arr.shape[1]//2]`.
            REVIEW: The last statement not always true depending on how I make the array.

          TODO: update attribute docstring (no longer sparse)
          aper_ys :: (M x 1) sparse array of floats
            The aperture array containing the y-coordinates (in arcsec) of the aperture.
            The origin (0, 0) is defined to be at the center of the source. The element at
            the center of this array is the center of the aperture (i.e., at
            `[aper_ys.shape[0]//2, aper_ys.shape[1]//2]`). Together with `aper_xs`, this
            sparse array can easily broacast into a full (i.e., M x N) 2D array (e.g.,
            `arr = aper_xs + aper_ys` will produce a full (M x N) 2D array). The center of
            the aperture in the full 2D array is located at `[arr.shape[0]//2,
            arr.shape[1]//2]`.
            REVIEW: The last statement not always true depending on how I make the array.

          _aper_extent :: 4-element 1D list of floats
            The [xmin, xmax, ymin, ymax] extent of the aperture in arcsec (for plotting
            the weight arrays).

        If any of the following have not been set (i.e., are None), then the array will be
        created with uniform weights:

          sky_background_weights :: (M x N) 2D array of floats
            The pixel weights for the sky background (incl. Earthshine, zodiacal light,
            geocoronal emission). This is currently all ones (1) (i.e., uniform noise).

          dark_current_weights :: (M x N) 2D array of floats
            The pixel weights for the dark current. This is currently all ones (1) (i.e.,
            uniform noise).

          redleak_weights :: (M x N) 2D array of floats
            The pixel weights for the redleak. This is currently all ones (1) (i.e.,
            uniform noise).
        """
        resolution = self.TelescopeObj.px_scale.to(u.arcsec).value
        half_res = 0.5 * resolution

        # Ensure center of aperture corresponds to center (i.e., shape // 2) pixel
        # xs_left = np.arange(-half_x, 0.0, resolution)
        # xs_right = np.arange(0.0, half_x + half_res, resolution)
        # # xs = np.concatenate((xs_left, xs_right)) + center[0]  # length N
        # xs = np.concatenate((xs_left, xs_right))  # length N
        # ys_left = np.arange(-half_y, 0.0, resolution)
        # ys_right = np.arange(0.0, half_y + half_res, resolution)
        # # ys = np.concatenate((ys_left, ys_right)) + center[1]  # length M
        # ys = np.concatenate((ys_left, ys_right)) # length M

        # (The following does not guarantee center of aperture gets its own pixel)
        xs = np.arange(-half_x, half_x + half_res, resolution)  # length N
        ys = np.arange(-half_y, half_y + half_res, resolution)  # length M
        # xs = np.arange(-half_x, half_x + half_res, resolution) + center[0]  # length N
        # ys = np.arange(-half_y, half_y + half_res, resolution) + center[1]  # length M

        self.aper_xs, self.aper_ys = np.meshgrid(xs, ys, sparse=False, indexing="xy")
        self._aper_extent = [
            xs[0] - half_res + center[0],
            xs[-1] + half_res + center[0],
            ys[0] - half_res + center[1],
            ys[-1] + half_res + center[1],
        ]  # +/- half_res will center tickmarks on pixels

        # For now, assume uniform background noise
        if self.sky_background_weights is None:
            self.sky_background_weights = np.ones_like(self.aper_xs)
        if self.dark_current_weights is None:
            self.dark_current_weights = np.ones_like(self.aper_xs)
        if self.redleak_weights is None:
            self.redleak_weights = np.ones_like(self.aper_xs)

    @staticmethod
    def calc_source_weights(profile, aper_xs, aper_ys, center):
        """
        Calculate the source weights for the given profile.

        Parameters
        ----------
          profile :: function
            TODO: explain call signature and purpose (weight function)

          half_x, half_y :: floats
            The half-widths of the aperture in x and y directions in arcsec.

          center :: 2-element 1D array of floats
            The (x, y) center of the aperture relative to the center of the source in
            arcsec. Positive values means the source is displaced to the bottom/left
            relative to the aperture center.

          res :: float
            The resolution of the weight map in arcsec.

        Returns
        -------
          aper_extent :: 4-element 1D list of floats

        TODO: finish docstring
        """
        source_weights = profile(aper_xs, aper_ys, center)
        return source_weights

    def show_source_weights(self, mark_source=False, source_markersize=4, plot=True):
        """
        Plot the source as seen through the aperture's FOV.

        Parameters
        ----------
          mark_source :: bool
            If True, mark the center of the aperture with a cyan dot.

          source_markersize :: int or float
            The markersize for the cyan point indicating the center of the source.

          plot :: bool
            If True, plot the source weights and return None. If False, return the figure,
            axis, and colorbar instance associated with the plot.

        Returns
        -------
          None (if plot is True)

          fig, ax, cbar (if plot is False) :: `matplotlib.figure.Figure`,
                                              `matplotlib.axes.Axes`,
                                              `matplotlib.colorbar.Colorbar`
            The figure, axis, and colorbar instance associated with the plot.
        """
        if self.source_weights is None or self._aper_extent is None:
            raise ValueError("Please select an aperture first.")

        rc = {"axes.grid": False}
        with plt.rc_context(rc):  # matplotlib v3.5.x has bug affecting grid + imshow
            fig, ax = plt.subplots()
            # (N.B. array already "xy" indexing. Do not transpose array)
            img = ax.imshow(
                self.source_weights,
                origin="lower",
                extent=self._aper_extent,
                interpolation=None,
                cmap="inferno",
            )
            ax.tick_params(color="grey", which="both")
            cbar = fig.colorbar(img)
            if mark_source:
                ax.plot(0, 0, "co", ms=source_markersize)
            cbar.set_label("Relative Weight to Center of Source")
            ax.set_xlabel("$x$ [arcsec]")
            ax.set_ylabel("$y$ [arcsec]")
            ax.set_title(
                "Effective No." + "\u00A0" + f"of Aperture Pixels: {self._eff_npix:.2f}"
            )
            if plot:
                plt.show()
            else:
                return fig, ax, cbar

    def set_background_weights(self, sky_background_weights):
        """
        Set the weights for the sky background.

        Parameters
        ----------
          sky_background_weights :: (M x N) 2D array of floats
            The pixel weights for the sky background (incl. Earthshine, zodiacal light,
            geocoronal emission). The shape must match the shape of the aperture array
            (see `source_weights` attribute).

        Attributes
        -------
          sky_background_weights :: (M x N) 2D array of floats
            The pixel weights for the sky background (incl. Earthshine, zodiacal light,
            geocoronal emission).
        """
        if self.source_weights is None:
            raise ValueError("Please select an aperture first.")
        if sky_background_weights.shape != self.source_weights.shape:
            raise ValueError(
                "`sky_background_weights` must have same shape as `source_weights`."
            )
        self.sky_background_weights = sky_background_weights

    def set_dark_current_weights(self, dark_current_weights):
        """
        Set the weights for the dark current.

        Parameters
        ----------
          dark_current_weights :: (M x N) 2D array of floats
            The pixel weights for the dark current. The shape must match the shape of the
            aperture array (see `source_weights` attribute).

        Attributes
        -------
          dark_current_weights :: (M x N) 2D array of floats
            The pixel weights for the dark current.
        """
        if self.source_weights is None:
            raise ValueError("Please select an aperture first.")
        if dark_current_weights.shape != self.source_weights.shape:
            raise ValueError(
                "`dark_current_weights` must have same shape as `source_weights`."
            )
        self.dark_current_weights = dark_current_weights

    def set_redleak_weights(self, redleak_weights):
        """
        Set the weights for the redleak.

        Parameters
        ----------
          redleak_weights :: (M x N) 2D array of floats
            The pixel weights for the redleak. The shape must match the shape of the
            aperture array (see `source_weights` attribute).

        Attributes
        -------
          redleak_weights :: (M x N) 2D array of floats
            The pixel weights for the redleak.
        """
        if self.source_weights is None:
            raise ValueError("Please select an aperture first.")
        if redleak_weights.shape != self.source_weights.shape:
            raise ValueError(
                "`redleak_weights` must have same shape as `source_weights`."
            )
        self.redleak_weights = redleak_weights

    def _calc_optimal_aperture_dimen(self, a, b, factor=1.4, print_info=True):
        """
        Calculates the "ideal/optimal" circular/elliptical aperture dimensions given the
        angular dimensions of a source. Note that if a dimension is smaller than half the
        telescope's FWHM, the half-FWHM value should be used for that dimension instead.

        Parameters
        ----------
          a, b :: `astropy.Quantity` angles
            The angle subtended by the semimajor and semiminor axes of the source,
            respectively. Note that if a dimension is smaller than half the telescope's
            FWHM, the half-FWHM should be used for that dimension instead.

          factor :: int or float
            The factor by which to scale the source's angular size along each dimension.

          print_info :: bool
            If True, print a message if an aperture dimension is based on the FWHM of the
            telescope's PSF instead of the point source's dimensions.

        Returns
        -------
          aper_dimen :: list of 2 `astropy.Quantity` angles
            The "ideal/optimal" circular/elliptical aperture dimensions along the
            semimajor (a) and semiminor (b) axes of the source, respectively.
        """
        if factor <= 0:
            raise ValueError("factor must be a positive scalar")

        half_fwhm = 0.5 * self.TelescopeObj.fwhm
        if self.SourceObj.angle_a < half_fwhm:
            if print_info:
                print(
                    "INFO: Source's semimajor axis (a) or radius subtends angle < "
                    + "1/2 of the FWHM of the telescope's PSF. Using 1/2 of FWHM instead."
                )
            a = half_fwhm
        if self.SourceObj.angle_b < half_fwhm:
            if print_info:
                print(
                    "INFO: Source's semiminor axis (b) or radius subtends angle < "
                    + "1/2 of the FWHM of the telescope's PSF. Using 1/2 of FWHM instead."
                )
            b = half_fwhm

        aper_dimen = [factor * a.to(u.arcsec), factor * b.to(u.arcsec)]
        return aper_dimen

    def use_optimal_aperture(self, factor=1.4, center=[0, 0] << u.arcsec, quiet=False):
        """
        Uses the "optimal" circular/elliptical aperture given a source's angular size and
        the telescope's FWHM. Note that this is only valid for point sources.

        Parameters
        ----------
          factor :: int or float
            The factor by which to scale the source's angular size along each dimension.

          center :: 2-element `astropy.Quantity` angles array
            The (x, y) center of the aperture relative to the center of the source.
            Positive values means the source is displaced to the bottom/left relative to
            the aperture center.

          quiet :: bool
            If False, print a message if an aperture dimension is based on the FWHM of the
            telescope's PSF instead of the point source's dimensions.

        Returns
        -------
          None
        """
        if not isinstance(self.SourceObj, PointSource):
            raise TypeError(
                "Only point sources are supported for optimal aperture. "
                + "Please manually define an aperture instead."
            )
        if factor <= 0:
            raise ValueError("factor must be a positive number")

        aper_dimen = self._calc_optimal_aperture_dimen(
            self.SourceObj.angle_a,
            self.SourceObj.angle_b,
            factor=factor,
            print_info=(not quiet),
        )
        self.aperture = np.pi * aper_dimen[0] * aper_dimen[1]
        self._assign_npix()

        #
        # Account for arbitrary source flux profiles through aperture
        #
        # Recall all internal angle aperture angles are in arcsec
        center = center.to(u.arcsec).value
        self._create_aper_arrs(aper_dimen[0].value, aper_dimen[1].value, center)
        source_weights = Photometry.calc_source_weights(
            self.SourceObj.profile, self.aper_xs, self.aper_ys, center
        )
        #
        # Create aperture
        #
        center_px = [
            0.5 * (source_weights.shape[1] - 1),  # x-coordinate in pixel units
            0.5 * (source_weights.shape[0] - 1),  # y-coordinate in pixel units
        ]
        aper_dimen_px = [
            (dimen / self.TelescopeObj.px_scale.to(u.arcsec)).value
            for dimen in aper_dimen
        ]
        aper = EllipticalAperture(
            positions=center_px, a=aper_dimen_px[0], b=aper_dimen_px[1], theta=0
        )
        # Restrict weight maps to aperture (which is an unrotated ellipse)
        aper_mask = aper.to_mask(method="exact").to_image(source_weights.shape)
        aper_mask[aper_mask <= 1e-12] = np.nan  # account for floating point errors
        self._aper_mask = aper_mask
        self.source_weights = source_weights * aper_mask
        self.sky_background_weights *= aper_mask
        self.dark_current_weights *= aper_mask
        self.redleak_weights *= aper_mask
        self._eff_npix = np.nansum(aper_mask)
        if abs(self._eff_npix - self.npix) > 0.1:
            # Discrepancy larger than 0.1 pixels
            warnings.warn(
                "Effective aperture area is off by more than 0.1 pixels... "
                + "Contact the developer with a minimal working example please. Thanks!"
            )

    def use_elliptical_aperture(self, a, b, center=[0, 0] << u.arcsec, rotation=0):
        """
        Use an elliptical aperture.

        Parameters
        ----------
          a, b :: int or float or `astropy.Quantity` angle
            The angular length of the semimajor and semiminor axes of the aperture,
            respectively. If int or float, a/b is assumed to be in pixel units. Otherwise,
            a/b must be an `astropy.Quantity` angle.

          center :: 2-element `astropy.Quantity` angles array
            The (x, y) center of the aperture relative to the center of the source.
            Positive values means the source is displaced to the bottom/left relative to
            the aperture center.

          rotation :: int or float
            The counter-clockwise rotation angle in degrees of the ellipse's semimajor
            axis from the positive x-axis. If rotation is 0, the semimajor axis is along
            the x-axis and the semiminor axis is along the y-axis.

        Returns
        -------
          None
        """
        if isinstance(a, u.Quantity):
            try:
                a = a.to(u.arcsec)
            except Exception:
                raise TypeError(
                    "a and b must be `astropy.Quantity` angles (e.g., u.arcsec, u.deg) "
                    + " or, if in pixel units, ints or floats"
                )
        else:
            a = a * self.TelescopeObj.px_scale.to(u.arcsec)
        if isinstance(b, u.Quantity):
            try:
                b = b.to(u.arcsec)
            except Exception:
                raise TypeError(
                    "a and b must be `astropy.Quantity` angles (e.g., u.arcsec, u.deg) "
                    + " or, if in pixel units, ints or floats"
                )
        else:
            b = b * self.TelescopeObj.px_scale.to(u.arcsec)
        if not isinstance(rotation, Number):
            raise TypeError("rotation must be an int or float")
        rotation = np.deg2rad(rotation)
        sin_rotate, cos_rotate = np.sin(rotation), np.cos(rotation)

        self.aperture = np.pi * a * b  # area of ellipse
        self._assign_npix()

        #
        # Account for arbitrary source flux profiles through aperture
        #
        # Recall all internal angle aperture angles are in arcsec
        center = center.to(u.arcsec).value
        # Modified rotation matrix to ensure full aperture is covered in 2D arrays
        x = (cos_rotate * a + sin_rotate * b).value
        y = (sin_rotate * a + cos_rotate * b).value
        x = x if x >= a.value else a.value
        y = y if y >= b.value else b.value
        self._create_aper_arrs(x, y, center)
        source_weights = Photometry.calc_source_weights(
            self.SourceObj.profile, self.aper_xs, self.aper_ys, center
        )
        #
        # Create aperture
        #
        center_px = [
            0.5 * (source_weights.shape[1] - 1),  # x-coordinate in pixel units
            0.5 * (source_weights.shape[0] - 1),  # y-coordinate in pixel units
        ]
        aper_dimen_px = [
            (a / self.TelescopeObj.px_scale.to(u.arcsec)).value,
            (b / self.TelescopeObj.px_scale.to(u.arcsec)).value,
        ]
        aper = EllipticalAperture(
            positions=center_px,
            a=aper_dimen_px[0],
            b=aper_dimen_px[1],
            theta=rotation,
        )
        # Restrict weight maps to aperture
        aper_mask = aper.to_mask(method="exact").to_image(source_weights.shape)
        aper_mask[aper_mask <= 1e-12] = np.nan  # account for floating point errors
        self._aper_mask = aper_mask
        self.source_weights = source_weights * aper_mask
        self.sky_background_weights *= aper_mask
        self.dark_current_weights *= aper_mask
        self.redleak_weights *= aper_mask
        self._eff_npix = np.nansum(aper_mask)
        if abs(self._eff_npix - self.npix) > 0.1:
            # Discrepancy larger than 0.1 pixels
            warnings.warn(
                "Effective aperture area is off by more than 0.1 pixels... "
                + "Contact the developer with a minimal working example please. Thanks!"
            )

        if isinstance(self.SourceObj, PointSource):
            aper_threshold = min(self._optimal_aperture_dimen)
            if a < aper_threshold or b < aper_threshold:
                warnings.warn(
                    "Chosen a/b is smaller than the 'ideal' aperture size "
                    + f"for this source! a & b should be at least {aper_threshold}.",
                    RuntimeWarning,
                )

    def use_rectangular_aperture(
        self, width, length, center=[0, 0] << u.arcsec, weight_res=1
    ):
        """
        Use a rectangular aperture.

        Parameters
        ----------
          width, length :: int or float or `astropy.Quantity` angle
            The width (the direction along the source's semimajor axis) and length (the
            direction along the source's semiminor axis) of the rectangular aperture. If
            int or float, length/width is assumed to be in pixel units.
            TODO: explain length/width above better

          center :: 2-element `astropy.Quantity` angles array
            The (x, y) center of the aperture relative to the center of the source.
            Positive values means the source is displaced to the bottom/left relative to
            the aperture center.

          weight_res :: float
            The resolution of the weight map in fraction of the telescope's pixel size.

        Returns
        -------
          None
        """
        #
        # Check inputs
        #
        if isinstance(width, u.Quantity):
            try:
                width = width.to(u.arcsec)
            except Exception:
                raise TypeError(
                    "length and width must be `astropy.Quantity` angles "
                    + "(e.g., u.arcsec, u.deg) or, if in pixel units, ints or floats"
                )
        else:
            width = width * self.TelescopeObj.px_scale.to(u.arcsec)
        if isinstance(length, u.Quantity):
            try:
                length = length.to(u.arcsec)
            except Exception:
                raise TypeError(
                    "length and width must be `astropy.Quantity` angles "
                    + "(e.g., u.arcsec, u.deg) or, if in pixel units, ints or floats"
                )
        else:
            length = length * self.TelescopeObj.px_scale.to(u.arcsec)
        #
        min_ifov_dimen = min(self.TelescopeObj.ifov_dimen)
        if (max(length, width) > max(self.TelescopeObj.ifov_dimen)) or (
            length > min_ifov_dimen and width > min_ifov_dimen
        ):
            raise ValueError(
                "aperture dimensions larger than telescope's IFOV dimensions"
            )
        #
        # Generate aperture
        #
        self.aperture = length * width
        self._assign_npix()

        #
        # Account for arbitrary source flux profiles through aperture
        #
        # Recall all internal angle aperture angles are in arcsec
        center = center.to(u.arcsec).value
        self._create_aper_arrs(
            0.5 * width.value,  # width is along x
            0.5 * length.value,  # length is along y
            center,
        )
        source_weights = Photometry.calc_source_weights(
            self.SourceObj.profile, self.aper_xs, self.aper_ys, center
        )
        #
        # Create aperture
        #
        center_px = [
            0.5 * (source_weights.shape[1] - 1),  # x-coordinate in pixel units
            0.5 * (source_weights.shape[0] - 1),  # y-coordinate in pixel units
        ]
        aper_dimen_px = [
            (width / self.TelescopeObj.px_scale.to(u.arcsec)).value,
            (length / self.TelescopeObj.px_scale.to(u.arcsec)).value,
        ]
        aper = RectangularAperture(
            positions=center_px, w=aper_dimen_px[0], h=aper_dimen_px[1], theta=0
        )
        # Restrict weight maps to aperture
        aper_mask = aper.to_mask(method="exact").to_image(source_weights.shape)
        aper_mask[aper_mask <= 1e-12] = np.nan  # account for floating point errors
        self._aper_mask = aper_mask
        self.source_weights = source_weights * aper_mask
        self.sky_background_weights *= aper_mask
        self.dark_current_weights *= aper_mask
        self.redleak_weights *= aper_mask
        self._eff_npix = np.nansum(aper_mask)
        if abs(self._eff_npix - self.npix) > 0.1:
            # Discrepancy larger than 0.1 pixels
            warnings.warn(
                "Effective aperture area is off by more than 0.1 pixels... "
                + "Contact the developer with a minimal working example please. Thanks!"
            )

        if isinstance(self.SourceObj, PointSource):
            aper_threshold = min(self._optimal_aperture_dimen) * 2
            if (length < aper_threshold) or (width < aper_threshold):
                warnings.warn(
                    "Chosen length/width is smaller than the 'ideal' aperture "
                    + "size for this source! "
                    + f"Length/width should be at least {aper_threshold}.",
                    RuntimeWarning,
                )

    def use_annular_aperture(self, a_out, b_out, a_in=0, b_in=None):
        """
        TODO: implement this? Maybe?
        """
        if b_in is None:
            b_in = b_out * (a_in / a_out)
        raise NotImplementedError("Not yet implemented")

    def _calc_redleaks(self, source_erate, mirror_area_cm_sq, include_redleak=True):
        """
        Calculate redleak of a source from its redleak fraction. The redleak fraction is
        defined to be the ratio of redleak electron/s to in-passband electron/s.

        This is a more robust way to calculate redleaks that is independent of how
        source_erate is determined.
        """
        #
        # Calculate redleak fraction (redleak electron/s to passband electron/s)
        #
        redleak_fracs = dict.fromkeys(self.TelescopeObj.passbands, 0.0)
        if include_redleak:
            source_wavelengths_AA = self.SourceObj.wavelengths.to(u.AA).value
            source_photon_s_A = (  # photon/s/A
                self.SourceObj.spectrum  # erg/s/cm^2/A
                * mirror_area_cm_sq  # cm^2
                / calc_photon_energy(wavelength=source_wavelengths_AA)[0]  # photon/erg
            )
            source_interp = interp1d(
                source_wavelengths_AA,
                source_photon_s_A,
                kind="linear",
                bounds_error=False,
                fill_value=np.nan,
            )  # photon/s/A
            for band in redleak_fracs:
                full_response_curve_wavelengths_AA = (
                    self.TelescopeObj.full_passband_curves[band]["wavelength"]
                    .to(u.AA)
                    .value
                )
                is_redleak = (
                    full_response_curve_wavelengths_AA
                    > self.TelescopeObj.redleak_thresholds[band].to(u.AA).value
                )
                is_in_passband = (
                    full_response_curve_wavelengths_AA
                    >= self.TelescopeObj.passband_limits[band][0].to(u.AA).value
                ) & (
                    full_response_curve_wavelengths_AA
                    <= self.TelescopeObj.passband_limits[band][1].to(u.AA).value
                )
                redleak_wavelengths = full_response_curve_wavelengths_AA[is_redleak]
                in_passband_wavelengths = full_response_curve_wavelengths_AA[
                    is_in_passband
                ]
                redleak_per_A = (
                    source_interp(redleak_wavelengths)
                    * self.TelescopeObj.full_passband_curves[band]["response"][is_redleak]
                )  # electron/s/A
                passband_erate_per_A = (
                    source_interp(in_passband_wavelengths)
                    * self.TelescopeObj.full_passband_curves[band]["response"][
                        is_in_passband
                    ]
                )  # electron/s/A
                isgood_redleak = np.isfinite(redleak_per_A)  # don't include NaNs
                isgood_passband = np.isfinite(passband_erate_per_A)  # don't include NaNs
                if not np.all(isgood_redleak) or not np.all(isgood_passband):
                    warnings.warn(
                        "Could not estimate redleak fraction "
                        + f"at 1 or more passband wavelengths in {band}-band",
                        RuntimeWarning,
                    )
                try:
                    redleak_per_px = simpson(  # electron/s (per px)
                        y=redleak_per_A[isgood_redleak],
                        x=redleak_wavelengths[isgood_redleak],
                        even="avg",
                    )
                    passband_erate_per_px = simpson(  # electron/s (per px)
                        y=passband_erate_per_A[isgood_passband],
                        x=in_passband_wavelengths[isgood_passband],
                        even="avg",
                    )
                    redleak_frac = redleak_per_px / passband_erate_per_px
                except Exception:
                    raise RuntimeError(
                        f"Unable to calculate redleak fraction for {band}-band! "
                        + "Please ensure there is at least 1 wavelength that is above "
                        + "the redleak threshold."
                    )
                if np.isfinite(redleak_frac):
                    redleak_fracs[band] = redleak_frac * self.redleak_weights
                else:
                    raise RuntimeError(
                        "Source redleak fraction could not be calculated "
                        + f"in {band}-band!"
                    )
        #
        # Calculate redleak from redleak fraction
        #
        redleaks = dict.fromkeys(self.TelescopeObj.passbands)
        for band in redleaks:
            redleaks[band] = redleak_fracs[band] * source_erate[band]
        #
        return redleaks

    @staticmethod
    def _calc_snr_from_t(
        t,
        signal,
        npix,
        totskynoise,
        readnoise,
        darkcurrent,
        redleak,
        nread=1,
    ):
        """
        Calculate the signal-to-noise ratio (SNR) reached given an integration time. Based on
        <https://hst-docs.stsci.edu/wfc3ihb/chapter-9-wfc3-exposure-time-calculation/9-6-estimating-exposure-times>.

        TODO: docstring outdated! Inputs should be arrays except for t/readnoise/npix!

        The equation to calculate the SNR is:
        ```math
                    SNR = (Q*t) / sqrt(Q*t + N_pix*t*Poisson + N_pix*N_read*Read^2)
        ```
        where:
        - SNR is the signal-to-noise ratio
        - t is the integration time in seconds
        - Q is the total signal due to the source in electrons/s
        - N_pix is the number of pixels occupied by the source on the detector
        FIXME: change name to Background or something, not Poisson
        - Poisson = (B_Earthshine + B_zodiacal + B_geocoronal + Darkcurrent + Redleak) is
            the total Poisson noise due to the sky backgorund (Earthshine, zodiacal light,
            and geocoronal emission), dark current, and red leak in electrons/s (per
            pixel)
        - N_read is the number of detector readouts
        - Read is the detector read noise in electrons (per pixel)


        Parameters
        ----------
        t :: int or float
            The integration time in seconds.

        signal :: scalar or array of scalars
            The signals of the source in electron/s. These are the total electron rates over
            the whole npix (see below).

        npix :: int
            The integer number of pixels occupied by the source. Only used for read noise.

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
        variables = [t, signal, totskynoise, readnoise, darkcurrent, redleak, nread]
        if np.any([isinstance(var, u.Quantity) for var in variables]):
            raise ValueError(
                "All inputs must be scalars or scalar arrays. "
                + "`astropy.Quantity` objects are not supported."
            )
        if not isinstance(npix, (int, np.integer)):
            raise ValueError("npix must be an integer")
        if not isinstance(nread, (int, np.integer)):
            raise ValueError("nread must be an integer")
        #
        # Calculate signal-to-noise ratio
        #
        signal_t = np.nansum(signal * t)  # electron
        noise = np.sqrt(
            signal_t
            + np.nansum(t * (totskynoise + darkcurrent + redleak))
            + (npix * readnoise * readnoise * nread)
        )  # electron
        snr = signal_t / noise
        # signal_t = signal * t  # electron
        # noise = np.sqrt(
        #     signal_t
        #     + (
        #         t * (totskynoise + darkcurrent + redleak)
        #         + (readnoise * readnoise * nread)
        #     )
        # )  # electron
        # snr = np.nansum(signal_t) / np.nansum(noise)
        return snr

    @staticmethod
    def _calc_t_from_snr(
        snr,
        signal,
        npix,
        totskynoise,
        readnoise,
        darkcurrent,
        redleak,
        nread=1,
    ):
        """
        Calculate the time required to reach a given signal-to-noise ratio (SNR). Based on
        <https://hst-docs.stsci.edu/wfc3ihb/chapter-9-wfc3-exposure-time-calculation/9-6-estimating-exposure-times>.

        TODO: docstring outdated! Inputs should be arrays except for snr/readnoise/npix!

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
        FIXME: change name to Background or something, not Poisson
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
            The signals of the source in electron/s. These are the total electron rates
            over the whole npix (see below).

        npix :: int
            The integer number of pixels occupied by the source. Only used for read noise.

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
        variables = [snr, signal, totskynoise, readnoise, darkcurrent, redleak, nread]
        if np.any([isinstance(var, u.Quantity) for var in variables]):
            raise ValueError(
                "All inputs must be scalars or scalar arrays. "
                + "`astropy.Quantity` objects are not supported."
            )
        if not isinstance(npix, (int, np.integer)):
            raise ValueError("npix must be an integer")
        if not isinstance(nread, (int, np.integer)):
            raise ValueError("nread must be an integer")
        # #
        # # Calculate useful quantities
        # #
        # snr_sq = snr * snr
        # signal_sq = signal * signal
        # poisson_noise = signal + npix * (totskynoise + darkcurrent + redleak)
        # #
        # # Calculate time to reach target SNR
        # #
        # numer1 = snr_sq * poisson_noise
        # numer2 = snr_sq * snr_sq * poisson_noise * poisson_noise
        # numer3 = 4 * snr_sq * signal_sq * (npix * readnoise * readnoise * nread)
        # t = (numer1 + np.sqrt(numer2 + numer3)) / (2 * signal_sq)  # seconds
        #
        # Calculate useful quantities
        #
        snr_sq = snr * snr
        tot_sig = np.nansum(signal)
        signal_sq = np.nansum(tot_sig * tot_sig)
        poisson_noise = tot_sig + np.nansum(totskynoise + darkcurrent + redleak)
        #
        # Calculate time to reach target SNR
        #
        numer1 = snr_sq * poisson_noise
        numer2 = snr_sq * snr_sq * poisson_noise * poisson_noise
        numer3 = 4 * snr_sq * signal_sq * (npix * readnoise * readnoise * nread)
        t = (numer1 + np.sqrt(numer2 + numer3)) / (2 * signal_sq)  # seconds
        return t

    def calc_snr_or_t(self, t=None, snr=None, nread=1, include_redleak=True):
        """
        Calculate the signal-to-noise ratio (SNR) reached in a given time or the
        integration time (t) required to reach a given SNR.

        The `Source` object associated with the `Photometry` instance should have its flux
        in units of flam (erg/s/cm^2/A).

        Parameters
        ----------
          t :: float
            Integration time in seconds.

          snr :: float
            Desired signal-to-noise ratio.

          nread :: int
            Number of detector readouts.

          include_redleak :: bool
            If False, the redleak contribution to the total noise is not included. May be
            useful if the source spectrum is, for example, uniform.

        Returns
        -------
          results :: dict
            If t is given, this is the snr reached after t seconds. If snr is given, this
            is the time in seconds required to reach the given snr.
        """
        if (t is None and snr is None) or (t is not None and snr is not None):
            raise ValueError("Exactly one of t or snr must be specified")
        #
        # Make some useful variables
        #
        response_curve_wavelengths_AA = dict.fromkeys(self.TelescopeObj.passbands)
        # passband_resolution_AA = self.TelescopeObj.passband_resolution.to(u.AA).value
        for band in response_curve_wavelengths_AA:
            response_curve_wavelengths_AA[band] = (
                self.TelescopeObj.passband_curves[band]["wavelength"].to(u.AA).value
            )
        #
        # Calculate sky background electron/s
        # (incl. Earthshine & zodiacal light, excl. geocoronal emission lines)
        #
        sky_background_erate = dict.fromkeys(self.TelescopeObj.passbands, 0.0)
        px_area_arcsec_sq = self.TelescopeObj.px_area.to(u.arcsec ** 2).value
        mirror_area_cm_sq = self.TelescopeObj.mirror_area.to(u.cm ** 2).value
        #
        if self.BackgroundObj.mags_per_sq_arcsec is None:
            #
            # Use passband response curves to calculate sky background electron/s
            # (will linearly interp sky background spectra to passband curve resolution)
            #
            # Earthshine contribution
            if self.BackgroundObj.earthshine_flam is not None:
                es_photon_s_A = (  # photon/s/A
                    self.BackgroundObj.earthshine_flam  # erg/s/cm^2/A
                    * mirror_area_cm_sq  # cm^2
                    / calc_photon_energy(  # photon/erg
                        wavelength=self.BackgroundObj.earthshine_wavelengths
                    )[0]
                )
                es_interp = interp1d(
                    self.BackgroundObj.earthshine_wavelengths,
                    es_photon_s_A,
                    kind="linear",
                    bounds_error=False,
                    fill_value=np.nan,
                )
                # Convolve with response curve (photon -> electron) & integrate over band
                for band in self.TelescopeObj.passbands:
                    es_erate_A = (
                        es_interp(response_curve_wavelengths_AA[band])
                        * self.TelescopeObj.passband_curves[band]["response"]
                    )  # electron/s/A
                    isgood_es = np.isfinite(es_erate_A)  # don't integrate over NaNs
                    if not np.all(isgood_es):
                        warnings.warn(
                            "Could not estimate Earthshine noise (electron/s/A) "
                            + f"at 1 or more passband wavelengths in {band}-band",
                            RuntimeWarning,
                        )
                    es_erate = simpson(
                        y=es_erate_A,
                        x=response_curve_wavelengths_AA[band][isgood_es],
                        even="avg",
                    ) / np.sum(
                        isgood_es
                    )  # ? I have to do this... Why ?
                    if np.isfinite(es_erate):
                        sky_background_erate[band] += es_erate
                    else:
                        raise RuntimeError(
                            "Earthshine noise (electron/s) could not be estimated "
                            + f"at all in {band}-band!"
                        )
            # Zodiacal light contribution
            if self.BackgroundObj.zodi_flam is not None:
                zodi_photon_s_A = (  # photon/s/A
                    self.BackgroundObj.zodi_flam  # erg/s/cm^2/A
                    * mirror_area_cm_sq  # cm^2
                    / calc_photon_energy(  # photon/erg
                        wavelength=self.BackgroundObj.zodi_wavelengths
                    )[0]
                )
                zodi_interp = interp1d(
                    self.BackgroundObj.zodi_wavelengths,
                    zodi_photon_s_A,
                    kind="linear",
                    bounds_error=False,
                    fill_value=np.nan,
                )
                # Convolve with response curve (photon -> electron) & integrate over band
                for band in self.TelescopeObj.passbands:
                    zodi_erate_A = (
                        zodi_interp(response_curve_wavelengths_AA[band])
                        * self.TelescopeObj.passband_curves[band]["response"]
                    )  # electron/s/A
                    isgood_zodi = np.isfinite(zodi_erate_A)  # don't integrate over NaNs
                    if not np.all(isgood_zodi):
                        warnings.warn(
                            "Could not estimate zodiacal light noise (electron/s/A) "
                            + f"at 1 or more passband wavelengths in {band}-band",
                            RuntimeWarning,
                        )
                    zodi_erate = simpson(
                        y=zodi_erate_A[isgood_zodi],
                        x=response_curve_wavelengths_AA[band][isgood_zodi],
                        even="avg",
                    ) / np.sum(
                        isgood_zodi
                    )  # ? I have to do this... Why ?
                    if np.isfinite(zodi_erate):
                        sky_background_erate[band] += zodi_erate
                    else:
                        raise RuntimeError(
                            "Zodiacal light noise (electron/s) could not be estimated "
                            + f"at all in {band}-band!"
                        )
        else:
            #
            # Use passband photometric zero points to calculate sky background electron/s
            #
            for band in self.TelescopeObj.passbands:
                try:
                    # Convert sky background AB mag per arcsec^2 (per pixel) to electron/s
                    sky_background_erate[band] = (
                        mag_to_flux(
                            self.BackgroundObj.mags_per_sq_arcsec[band],
                            zpt=self.TelescopeObj.phot_zpts[band],
                        )[0]
                        * px_area_arcsec_sq
                    )
                except Exception:
                    raise KeyError(
                        "No sky background magnitude (`mags_per_sq_arcsec` from "
                        + "`Background` object) or photometric zero point (`phot_zpts` "
                        + f"from `Telescope` object) for {band}-band!\n"
                        + "(The issue is likely with the `Background` object...)"
                    )
        #
        # Add geocoronal emission line contribution to sky background
        #
        for gw, gf, gl in zip(
            self.BackgroundObj.geo_wavelength,
            self.BackgroundObj.geo_flux,
            self.BackgroundObj.geo_linewidth,
        ):
            for band in self.TelescopeObj.passbands:
                # Add geocoronal emission (electron/s) to proper passband(s)
                if (gw >= response_curve_wavelengths_AA[band][0]) and (
                    gw <= response_curve_wavelengths_AA[band][-1]
                ):
                    # (Doing this in the if statement since each geocoronal emission line
                    # is likely only in 1 band. Reduces unnecessary computation)
                    geo_photon_rate = (
                        gf
                        * px_area_arcsec_sq
                        * mirror_area_cm_sq
                        / calc_photon_energy(wavelength=gw)[0]
                    )  # photon/s
                    response_interp = interp1d(
                        response_curve_wavelengths_AA[band],
                        self.TelescopeObj.passband_curves[band]["response"],
                        kind="linear",
                        bounds_error=False,
                        fill_value=np.nan,
                    )
                    geo_erate = response_interp(gw) * geo_photon_rate  # electron/s
                    if not np.isfinite(geo_erate):
                        warnings.warn(
                            "Could not estimate geocoronal emission noise contribution "
                            + f"(electron/s) in {band}-band!",
                            RuntimeWarning,
                        )
                    else:
                        sky_background_erate[band] += geo_erate
                    # (Now don't break out of loop: check other bands in case geocoronal
                    # emission line is in multiple bands)
        #
        # Sky background noise (electron/s) is present in every pixel in aperture
        #
        for band in sky_background_erate:
            sky_background_erate[band] *= self.sky_background_weights
        #
        # Calculate signal in each passband (flam -> electron/s)
        #
        if self.source_weights is None:
            raise ValueError("Please choose an aperture first.")
        #
        source_erate = dict.fromkeys(self.TelescopeObj.passbands, np.nan)

        # source_wavelengths_AA = self.SourceObj.wavelengths.to(u.AA).value
        # source_photon_s_A = (  # photon/s/A
        #     self.SourceObj.spectrum  # erg/s/cm^2/A
        #     * mirror_area_cm_sq  # cm^2
        #     / calc_photon_energy(wavelength=source_wavelengths_AA)[0]  # photon/erg
        # )
        # source_interp = interp1d(
        #     source_wavelengths_AA,
        #     source_photon_s_A,
        #     kind="linear",
        #     bounds_error=False,
        #     fill_value=np.nan,
        # )  # photon/s/A
        # # Convolve with response curve (photon -> electron) & integrate over band
        # for band in source_erate:
        #     source_erate_A = (
        #         source_interp(response_curve_wavelengths_AA[band])
        #         * self.TelescopeObj.passband_curves[band]["response"]
        #     )  # electron/s/A
        #     isgood_source = np.isfinite(source_erate_A)  # don't integrate over NaNs
        #     if not np.all(isgood_source):
        #         warnings.warn(
        #             "Could not estimate source signal (electron/s/A) "
        #             + f"at 1 or more passband wavelengths in {band}-band",
        #             RuntimeWarning,
        #         )
        #     source_erate_per_px = simpson(
        #         y=source_erate_A[isgood_source],
        #         x=response_curve_wavelengths_AA[band][isgood_source],
        #         even="avg",
        #     )
        #     print("source_erate_per_px:", band, source_erate_per_px)
        #     if np.isfinite(source_erate_per_px):
        #         source_erate[band] = source_erate_per_px * self.source_weights
        #     else:
        #         raise RuntimeError(
        #             "Signal of source (in electron/s) could not be calculated "
        #             + f"in {band}-band!"
        #         )
        #     print("(total) source_erate", band, np.nansum(source_erate[band]))
        #     print(
        #         "redleak fraction",
        #         band,
        #         np.nansum(redleaks[band]) / np.nansum(source_erate[band]),
        #     )

        for band in source_erate:
            is_in_passband = (
                self.SourceObj.wavelengths >= self.TelescopeObj.passband_limits[band][0]
            ) & (self.SourceObj.wavelengths <= self.TelescopeObj.passband_limits[band][1])
            tot_passband_flam_per_px = (
                np.nansum(self.SourceObj.spectrum[is_in_passband]) / self._eff_npix
            )
            # Convert flux to electron/s
            source_erate[band] = convert_electron_flux_mag(
                tot_passband_flam_per_px * self.source_weights,
                "flam",
                "electron",
                phot_zpt=self.TelescopeObj.phot_zpts[band],
                wavelengths=self.TelescopeObj.passband_pivots[band],
            )[0]
        #
        # Calculate redleak
        #
        redleaks = self._calc_redleaks(
            source_erate, mirror_area_cm_sq, include_redleak=include_redleak
        )
        #
        # Calculate desired results (either integration time given SNR or SNR given time)
        #
        results = dict.fromkeys(self.TelescopeObj.passbands)
        dark_current = self.TelescopeObj.dark_current * self.dark_current_weights
        if t is not None:
            for band in results:
                results[band] = Photometry._calc_snr_from_t(
                    t=t,
                    signal=source_erate[band],
                    # npix=np.ceil(self._eff_npix).astype(int),  # only used for read noise
                    npix=int(np.ceil(self._eff_npix)),  # only used for read noise
                    totskynoise=sky_background_erate[band],
                    readnoise=self.TelescopeObj.read_noise,  # constant per pixel
                    darkcurrent=dark_current,
                    redleak=redleaks[band],
                    nread=nread,
                )
        else:
            for band in results:
                results[band] = Photometry._calc_t_from_snr(
                    snr=snr,
                    signal=source_erate[band],
                    # npix=np.ceil(self._eff_npix).astype(int),  # only used for read noise
                    npix=int(np.ceil(self._eff_npix)),  # only used for read noise
                    totskynoise=sky_background_erate[band],
                    readnoise=self.TelescopeObj.read_noise,  # constant per pixel
                    darkcurrent=dark_current,
                    redleak=redleaks[band],
                    nread=nread,
                )
        return results
