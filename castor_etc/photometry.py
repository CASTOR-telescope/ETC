"""
Utilities for photometric calculations.

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
from copy import deepcopy
from numbers import Number

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from photutils.aperture import EllipticalAperture, RectangularAperture
from scipy.integrate import simpson
from scipy.interpolate import interp1d

from .background import Background
from .conversions import calc_photon_energy, mag_to_flux
from .sources import PointSource, Source
from .telescope import Telescope

# TODO: convolve with PSF (waiting for data file)
# TODO: get encircled energy as function of aperture radius (for point sources)


# The optimal aperture for a point source is a circular aperture with a radius equal
# to the factor below times half the telescope's FWHM
_OPTIMAL_APER_FACTOR = 1.4


class Photometry:
    """
    Photometry class.
    """

    def __init__(self, TelescopeObj, SourceObj, BackgroundObj):
        """
        Initialize class for photometry calculations.

        Note that the `Source` object (i.e., the `SourceObj` parameter) should have its
        spectrum in units of flam (erg/s/cm^2/A).

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

        Attributes
        ----------
          TelescopeObj :: `castor_etc.Telescope` object
            The `castor_etc.Telescope` object containing the telescope parameters.

          SourceObj :: `castor_etc.Source` object
            The `castor_etc.Source` object contaning the target source parameters.

          BackgroundObj :: `castor_etc.Background` object
            The `castor_etc.Background` object containing the background parameters.

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
        #
        # Assign attributes
        #
        self.TelescopeObj = TelescopeObj
        self.SourceObj = SourceObj
        self.BackgroundObj = BackgroundObj
        #
        # Initialize attributes that will be used in the future
        #
        # Weights for each pixel in aperture. NaNs are excluded from all calculations
        self.source_weights = None
        self.sky_background_weights = None
        self.dark_current_weights = None
        self.redleak_weights = None
        # Attributes for internal use
        self._aper_area = None  # exact area of the aperture from given aperture params
        self._aper_xs = None  # array containing x-coordinates of pixels
        self._aper_ys = None  # array containing y-coordinates of pixels
        self._aper_mask = None  # photutils aperture mask
        self._exact_npix = None  # number of pixels calculated from given aperture params
        self._eff_npix = None  # number of pixels calculated from photutils mask
        self._aper_extent = None  # the [xmin, xmax, ymin, ymax] extent of the imshow plot
        self._source_aper_overlap_arr = None  # array showing the source-aperture overlap
        self._source_aper_overlap_frac = None  # the source-aperture overlap fraction
        # Default optimal aperture dimensions for a point source
        self._optimal_aperture_radius = _OPTIMAL_APER_FACTOR * 0.5 * TelescopeObj.fwhm

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

    def _assign_exact_npix(self):
        """
        Internal function. Calculate the exact number of pixels in the aperture based on
        user parameters. Will compare this value to photutils' number of pixels.

        Parameters
        ----------
          None

        Attributes
        __________
          _exact_npix :: float
            The exact number of pixels in the aperture calculated from user parameters.

        Returns
        -------
          None
        """
        self._exact_npix = (
            self._aper_area.to(u.arcsec ** 2)
            / self.TelescopeObj.px_area.to(u.arcsec ** 2)
        ).value

    def _create_aper_arrs(self, half_x, half_y, center, overwrite=False):
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

          overwrite :: bool
            If True, allow overwriting of any existing aperture arrays.

        Attributes
        ----------
          _aper_xs :: (M x N) 2D array of floats
            The aperture array containing the x-coordinates (in arcsec) relative to the
            center of the aperture.

          _aper_ys :: (M x N) 2D array of floats
            The aperture array containing the y-coordinates (in arcsec) relative to the
            center of the aperture.

          _aper_extent :: 4-element 1D list of floats
            The [xmin, xmax, ymin, ymax] extent of the aperture in arcsec (for plotting
            the weight arrays).

          sky_background_weights :: (M x N) 2D array of floats
            The pixel weights for the sky background (incl. Earthshine, zodiacal light,
            geocoronal emission). This is currently all ones (1) (i.e., uniform noise).

          dark_current_weights :: (M x N) 2D array of floats
            The pixel weights for the dark current. This is currently all ones (1) (i.e.,
            uniform noise).

          redleak_weights :: (M x N) 2D array of floats
            The pixel weights for the redleak. This is currently all ones (1) (i.e.,
            uniform noise).

        Returns
        -------
          center_px :: 2-element 1D list of floats
            The (x, y) center index of the aperture coordinates arrays (i.e., _aper_xs and
            _aper_ys). To be very explicit, this is not (0, 0) but rather the center of
            the 2D arrays.

          source_center_px :: 2-element 1D list of floats
            The (x, y) center of the source in pixel units. These values can be negative
            or larger than the aperture coordinate array size, which means the center of
            the source is outside the aperture array. N.B. the aperture array may be
            slightly larger than the aperture (but these excess pixels are mapped to
            NaNs).
        """
        if not overwrite and (self._aper_xs is not None or self._aper_ys is not None):
            raise ValueError(
                "An aperture for this `Photometry` object already exists. "
                + "Use `overwrite=True` to allow overwriting of the aperture "
                + " and all associated weights (i.e., source, sky background, "
                + "dark current, and red leak weights will all be reset)."
            )
        px_scale_arcsec = self.TelescopeObj.px_scale.to(u.arcsec).value
        half_px_scale = 0.5 * px_scale_arcsec  # to ensure correct arange/extent

        xs = np.arange(-half_x, half_x + half_px_scale, px_scale_arcsec)  # length N
        ys = np.arange(-half_y, half_y + half_px_scale, px_scale_arcsec)  # length M

        self._aper_xs, self._aper_ys = np.meshgrid(xs, ys, sparse=False, indexing="xy")
        self._aper_extent = [
            xs[0] - half_px_scale + center[0],
            xs[-1] + half_px_scale + center[0],
            ys[0] - half_px_scale + center[1],
            ys[-1] + half_px_scale + center[1],
        ]  # we use +/- half_px_scale to center matplotlib tickmarks on pixels

        # Assume uniform background noise, dark current, and red leak
        self.sky_background_weights = np.ones_like(self._aper_xs)
        self.dark_current_weights = np.ones_like(self._aper_xs)
        self.redleak_weights = np.ones_like(self._aper_xs)

        # Find pixel that corresponds to center of the source
        center_px = [
            0.5 * (self._aper_xs.shape[1] - 1),  # x-coordinate in pixel units
            0.5 * (self._aper_ys.shape[0] - 1),  # y-coordinate in pixel units
        ]
        # Recall point-slope form of a line: (y - y0) = m * (x - x0)
        num_px_offset_x = (-center[0] + xs[0]) / px_scale_arcsec + center_px[0]
        num_px_offset_y = (-center[1] + ys[0]) / px_scale_arcsec + center_px[1]
        source_center_px_x = center_px[0] + num_px_offset_x
        source_center_px_y = center_px[1] + num_px_offset_y

        return center_px, [source_center_px_x, source_center_px_y]

    @staticmethod
    def _calc_source_weights(profile, aper_xs, aper_ys, center):
        """
        Calculate the source weights for the given profile.

        Parameters
        ----------
          profile :: function with a header of `profile(x, y, aper_center)`
            Function with 3 positional parameters that describes the flux of the
            source at each aperture pixel relative to the flux at the source center.
              - `x` and `y` are 2D arrays representing the angular distance, in arcsec,
              of each pixel from the aperture's center,
              - `aper_center` is a 2-element 1D array of floats representing the x- and y-
              coordinates of the aperture's center relative to the source's center,
              respectively.

          half_x, half_y :: floats
            The half-widths of the aperture in x- and y- directions in arcsec. (e.g.,
            semimajor/semiminor axes, half of a rectangle's length/width, etc.)

          center :: 2-element 1D array of floats
            The (x, y) center of the aperture relative to the center of the source in
            arcsec. Positive values means the source is displaced to the bottom/left
            relative to the aperture center.

        Returns
        -------
          source_weights :: (M x N) 2D array of floats
            The source weights for each pixel in the aperture. These represent the flux of
            the source at each pixel relative to the flux at the center of the source.
        """
        return profile(aper_xs, aper_ys, center)

    def _do_source_aper_overlap(self, source_center_px, px_scale_arcsec):
        """
        Internal function to find the fraction of the source within the aperture. Note
        that the `source_weights` attribute should already be created and masked with the
        `_aper_mask` attribute.

        Parameters
        ----------
          source_center_px :: 2-element 1D list of floats
            The (x, y) center of the source in pixel units. These values can be negative
            or larger than the aperture coordinate array size, which means the center of
            the source is outside the aperture array. N.B. the aperture array may be
            slightly larger than the aperture (but these excess pixels are mapped to
            NaNs).

          px_scale_arcsec :: float
            The `TelescopeObj` object's pixel scale in units of arcsec. Only passing this
            as a parameter to avoid having to calculate it multiple times.

        Attributes
        ----------
          ! (`source_weights` NOT MODIFIED ANYMORE) !
          source_weights :: (M x N) 2D array of floats
            The source weights for each pixel in the aperture. Any pixel within the
            aperture but outside of the source is now zero. Note that the extent of the
            source is determined by the source's semimajor and semiminor axes.

          _source_aper_overlap_arr :: (M x N) 2D array of floats
            The array showing the source and the aperture overlap.

          _source_aper_overlap_frac :: float
            The fraction of the source within the aperture.

        Returns
        -------
          None
        """
        #
        # Create an ellipse representing the source
        #
        source_ellipse = (
            EllipticalAperture(
                positions=source_center_px,
                a=self.SourceObj.angle_a.to(u.arcsec).value / px_scale_arcsec,
                b=self.SourceObj.angle_b.to(u.arcsec).value / px_scale_arcsec,
                theta=self.SourceObj.rotation,  # already in radians
            )
            .to_mask(method="exact")
            .to_image(self.source_weights.shape)
        )
        if source_ellipse is None:
            raise ValueError("Source is not contained within the aperture!")
        source_ellipse[source_ellipse >= 1e-12] = 1.0
        source_ellipse[source_ellipse <= 1e-12] = 0.0
        #
        # Account for fraction of source contained within aperture
        #
        # self.source_weights *= source_ellipse  # REVIEW: don't do this b/c sharp cutoff? Remember that uniform profile is uniform over whole aperture.
        self._source_aper_overlap_arr = source_ellipse * self._aper_mask
        self._source_aper_overlap_frac = np.nansum(self._source_aper_overlap_arr) / (
            self.SourceObj.area.to(u.arcsec ** 2).value
            / self.TelescopeObj.px_area.to(u.arcsec ** 2).value
        )
        if self._source_aper_overlap_frac > 1.0:
            # Sometimes point sources are too small and result in >100% overlap
            self._source_aper_overlap_frac = 1.0

        # TODO: actually integrate PSF pixel-by-pixel in the future
        # For now, just assume 95% encircled energy for point source
        if isinstance(self.SourceObj, PointSource):
            self._source_aper_overlap_frac = 0.95

    def show_source_weights(
        self, mark_source=False, source_markersize=4, norm=None, plot=True
    ):
        """
        Plot the source as seen through the photometry aperture. The pixels are colored by
        the flux of the source at each pixel relative to the flux at the center of the
        source. Coloring also includes the effects of fractional pixel weights (i.e., from
        the aperture mask, visualized using `show_aper_weights()`). These two effects
        combined give the "source weights".

        The "percentage of source within aperture" is based on the `SourceObj` object's
        angular area, which was specified when the `SourceObj` object was created. This
        percentage is is the projected area of the source contained within the aperture
        divided by the (total) projected area of the source.

        Parameters
        ----------
          mark_source :: bool
            If True, mark the center of the aperture with a cyan dot.

          source_markersize :: int or float
            The markersize for the cyan point indicating the center of the source.

          norm :: `matplotlib.colors` normalization class (e.g., `LogNorm`) or None
            The scaling and normalization to use for the colorbar. If None, then a linear
            scaling with the default (min pixel value, max pixel value) bounds are used.

          plot :: bool
            If True, plot the source weights and return None. If False, return the figure,
            axis, and colorbar instance associated with the plot.

        Returns
        -------
        If plot is True:
          None

        If plot is False:
          fig, ax :: `matplotlib.figure.Figure` and `matplotlib.axes.Axes` objects
            The figure and axis instance associated with the plot.

          cbar :: `matplotlib.colorbar.Colorbar` object
            The colorbar instance associated with the plot.
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
                norm=norm,
            )
            ax.tick_params(color="grey", which="both")
            cbar = fig.colorbar(img)
            if mark_source:
                ax.plot(0, 0, "co", ms=source_markersize)
            cbar.set_label("Flux Relative to Center of Source")
            ax.set_xlabel("$x$ [arcsec]")
            ax.set_ylabel("$y$ [arcsec]")
            # ax.set_title(
            #     "Effective No." + "\u00A0" + f"of Aperture Pixels: {self._eff_npix:.2f}"
            # )
            ax.set_title(
                "Percentage of Source Within Aperture: "
                + f"{100 * self._source_aper_overlap_frac:.2f}"
                + r"\%"
            )
            if plot:
                plt.show()
            else:
                return fig, ax, cbar

    def show_aper_weights(self, plot=True):
        """
        Plot the aperture photometry mask. The aperture mask shows the fractional overlap
        between the aperture and the pixel, ranging from (0, 1]. A pixel that is wholly
        contained in the aperture has an aperture weight equal to 1. Similarly, a pixel
        that is partially contained in the aperture has a weight between 0 and 1
        (exclusive) directly proportional to the area of the aperture that overlaps the
        pixel. A pixel that is wholly outside the aperture is assigned a weight of NaN.

        The "effective number of aperture pixels" is the sum of these pixel values
        (excluding all NaNs). This approach using fractional pixel weights gives a more
        accurate pixel count for non-rectangular aperture shapes while maintaning the
        benefit of allowing pixel-by-pixel modifications, if desired. This is to emulate
        the physical process of reading out and summing up each pixel in the aperture.

        Parameters
        ----------
          plot :: bool
            If True, plot the aperture weights and return None. If False, return the
            figure, axis, and colorbar instance associated with the plot.

        Returns
        -------
        If plot is True:
          None

        If plot is False:
          fig, ax :: `matplotlib.figure.Figure` and `matplotlib.axes.Axes` objects
            The figure and axis instance associated with the plot.

          cbar :: `matplotlib.colorbar.Colorbar` object
            The colorbar instance associated with the plot.
        """
        if self._aper_mask is None or self._aper_extent is None:
            raise ValueError("Please select an aperture first.")

        rc = {"axes.grid": False}
        with plt.rc_context(rc):  # matplotlib v3.5.x has bug affecting grid + imshow
            fig, ax = plt.subplots()
            # (N.B. array already "xy" indexing. Do not transpose array)
            img = ax.imshow(
                self._aper_mask,
                origin="lower",
                extent=self._aper_extent,
                interpolation=None,
                cmap="viridis",
            )
            ax.tick_params(color="grey", which="both")
            cbar = fig.colorbar(img)
            cbar.set_label("Relative Weight to Center of Aperture")
            ax.set_xlabel("$x$ [arcsec]")
            ax.set_ylabel("$y$ [arcsec]")
            ax.set_title(
                "Effective No." + "\u00A0" + f"of Aperture Pixels: {self._eff_npix:.2f}"
            )
            if plot:
                plt.show()
            else:
                return fig, ax, cbar

    def show_source_aper_overlap(self, plot=True):
        """
        Plot the pixel-by-pixel overlap between the source object and the photometry
        aperture (the latter includes fractional pixel weights).

        The "percentage of source within aperture" is based on the `SourceObj` object's
        angular area, which was specified when the `SourceObj` object was created. This
        percentage is is the projected area of the source contained within the aperture
        divided by the (total) projected area of the source.

        Parameters
        ----------
          plot :: bool
            If True, plot the overlap between the source object and the aperture, then
            return None. If False, return the figure, axis, and colorbar instance
            associated with the plot.

        Returns
        -------
        If plot is True:
          None

        If plot is False:
          fig, ax :: `matplotlib.figure.Figure` and `matplotlib.axes.Axes` objects
            The figure and axis instance associated with the plot.

          cbar :: `matplotlib.colorbar.Colorbar` object
            The colorbar instance associated with the plot.
        """
        rc = {"axes.grid": False}
        with plt.rc_context(rc):  # matplotlib v3.5.x has bug affecting grid + imshow
            fig, ax = plt.subplots()
            # (N.B. array already "xy" indexing. Do not transpose array)
            img = ax.imshow(
                self._source_aper_overlap_arr,
                origin="lower",
                extent=self._aper_extent,
                interpolation=None,
                cmap="copper",
            )
            ax.tick_params(color="grey", which="both")
            cbar = fig.colorbar(img)
            cbar.set_label(
                "Fraction of Source Overlapped with Aperture Pixel\n"
                + "(including fractional pixel weights from aperture)"
            )
            ax.set_xlabel("$x$ [arcsec]")
            ax.set_ylabel("$y$ [arcsec]")
            ax.set_title(
                "Percentage of Source Within Aperture: "
                + f"{100 * self._source_aper_overlap_frac:.2f}"
                + r"\%"
            )
            if plot:
                plt.show()
            else:
                return fig, ax, cbar

    def set_background_weights(self, sky_background_weights):
        """
        Set the pixel weights for the sky background. The value of these weights represent
        the amount of sky background noise at each pixel. By default, the sky background
        weights are uniform (i.e., equal to 1.0) over the aperture (except at the aperture
        edges, where there may be fractional pixel weighting as described in the
        `show_aper_weights()` docstring).

        When doing photometry calculations, these sky background weights will be
        multiplied with the calculated sky background noise per pixel (includes
        Earthshine, zodiacal light, and geocoronal emission) and summed pixel-by-pixel to
        determine the total background noise contribution; any pixels with a value of NaN
        are excluded from the summation.

        No restrictions are placed on the value of the sky background weights, but the
        user is responsible for ensuring they make sense in the relevant context. For
        example, if one pixel in the aperture has double the sky background noise while
        the rest of the pixels have ordinary levels of background noise, the user should
        set that one pixel to have a background weight of 2.0 and keep the other
        background weights at 1.0 (as opposed to setting that one pixel weight to 1.0
        while changing the other weights to 0.5). As a reminder, this is because the total
        sky background noise is based on multiplying the sky background noise per pixel by
        these background weights---the latter scenario would decrease the total sky
        background noise as opposed to slightly increasing it...

        Parameters
        ----------
          sky_background_weights :: (M x N) 2D array of floats
            The pixel weights for the sky background (incl. Earthshine, zodiacal light,
            geocoronal emission). The shape must match the shape of the aperture array
            (see the `source_weights` attribute).

        Attributes
        -------
          sky_background_weights :: (M x N) 2D array of floats
            The pixel weights for the sky background. Note that the sky background
            includes Earthshine, zodiacal light, and geocoronal emission.

        Returns
        -------
          None
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
        Set the pixel weights for the dark current. The value of these weights represent
        the amount of dark current noise at each pixel. By default, the dark current
        weights are uniform (i.e., equal to 1.0) over the aperture (except at the aperture
        edges, where there may be fractional pixel weighting as described in the
        `show_aper_weights()` docstring).

        When doing photometry calculations, these dark current weights will be multiplied
        with the dark current value per pixel (specified in the `Telescope` object) and
        summed pixel-by-pixel to determine the total dark current contribution; any pixels
        with a value of NaN are excluded from the summation.

        No restrictions are placed on the value of the dark current weights, but the user
        is responsible for ensuring they make sense in the relevant context. For example,
        if one pixel in the aperture has double the dark current while the rest of the
        pixels have the nominal dark current, the user should set that one pixel to have a
        dark current weight of 2.0 and keep the other dark current weights at 1.0 (as
        opposed to setting that one pixel weight to 1.0 while changing the other weights
        to 0.5). As a reminder, this is because the total dark current noise is based on
        multiplying the dark current per pixel by these dark current weights---the latter
        scenario would decrease the total dark current noise as opposed to slightly
        increasing it...

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
        Set the pixel weights for the red leak. The value of these weights represent the
        amount of red leak contamination in each pixel. By default, the red leak weights
        are uniform (i.e., equal to 1.0) over the aperture. Unlike all the other pixel
        weights, the red leak weights do NOT include fractional pixel weighting based on
        the aperture mask to prevent double-counting of the aperture weights during the
        red leak calculation.

        When doing photometry calculations, these red leak weights will be multiplied with
        the calculated red leak value per pixel and summed pixel-by-pixel to determine the
        total red leak noise; any pixels with a value of NaN are excluded from the
        summation.

        No restrictions are placed on the value of the red leak weights, but the user is
        responsible for ensuring they make sense in the relevant context. For example, if
        one pixel in the aperture has double the red leak flux while the rest of the
        pixels have the average red leak value, the user should set that one pixel to have
        a red leak weight of 2.0 and keep the other red leak weights at 1.0 (as opposed to
        setting that one pixel weight to 1.0 while changing the other weights to 0.5). As
        a reminder, this is because the total red leak contamination is based on
        multiplying the red leak per pixel by these red leak weights---the latter scenario
        would decrease the total red leak value as opposed to slightly increasing it...

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

    def use_optimal_aperture(
        self, factor=_OPTIMAL_APER_FACTOR, quiet=False, overwrite=False
    ):
        """
        Uses the "optimal" circular aperture calculated from the telescope's FWHM. Note
        that this is only valid for point sources.

        Note that the encircled energy for a point source with this aperture is going to
        be approximated as 95% for the photometry calculation. This will change in the
        future as we acquire the necessary PSF data.

        Parameters
        ----------
          factor :: int or float
            The factor by which to scale the telescope's FWHM.

          quiet :: bool
            If False, print a message if the point source's diameter is larger than the
            telescope's FWHM.

          overwrite :: bool
            If True, allow overwriting of any existing aperture associated with this
            `Photometry` object.

        Attributes
        ----------
          source_weights :: (M x N) 2D array of floats
            The pixel weights for the source.

          sky_background_weights :: (M x N) 2D array of floats
            The pixel weights for the sky background.

          dark_current_weights :: (M x N) 2D array of floats
            The pixel weights for the dark current.

          redleak_weights :: (M x N) 2D array of floats
            The pixel weights for the red leak.

        Returns
        -------
          None
        """
        #
        # Check inputs
        #
        if not isinstance(self.SourceObj, PointSource):
            raise TypeError(
                "Only point sources are supported for optimal aperture. "
                + "Please manually define an aperture instead."
            )
        if factor <= 0:
            raise ValueError("factor must be a positive number")
        if not quiet and self.SourceObj.angle_a > 0.5 * self.TelescopeObj.fwhm:
            warnings.warn(
                "The point source's diameter is larger than the telescope's FWHM. "
                + "Will use the telescope's FWHM for calculating the optimal aperture "
                + "and, for now, the encircled energy fraction.",
                RuntimeWarning,
            )
        #
        # Calculate the optimal aperture dimensions, aperture area, and exact # pixels
        #
        aper_radius_arcsec = factor * 0.5 * self.TelescopeObj.fwhm.to(u.arcsec)
        self._aper_area = np.pi * aper_radius_arcsec * aper_radius_arcsec
        self._assign_exact_npix()
        #
        # Create source weights with arbitrary source flux profile through aperture
        #
        # Recall all internal angle aperture angles are in arcsec
        center = [0, 0]  # arcsec
        # N.B. must round to nearest multiple of px_scale or else rounding errors will
        # affect source weights
        px_scale_arcsec = self.TelescopeObj.px_scale.to(u.arcsec).value
        center_px, source_center_px = self._create_aper_arrs(
            np.ceil(aper_radius_arcsec.value / px_scale_arcsec) * px_scale_arcsec,
            np.ceil(aper_radius_arcsec.value / px_scale_arcsec) * px_scale_arcsec,
            center,
            overwrite=overwrite,
        )
        source_weights = Photometry._calc_source_weights(
            self.SourceObj.profile, self._aper_xs, self._aper_ys, center
        )
        #
        # Create aperture
        #
        aper_radius_px = (aper_radius_arcsec / px_scale_arcsec).value
        aper = EllipticalAperture(
            positions=center_px, a=aper_radius_px, b=aper_radius_px, theta=0
        )
        # Restrict weight maps to aperture (which is an unrotated ellipse)
        aper_mask = aper.to_mask(method="exact").to_image(source_weights.shape)
        aper_mask[aper_mask <= 1e-12] = np.nan  # account for floating point errors
        self._aper_mask = aper_mask
        self.source_weights = source_weights * aper_mask
        self.sky_background_weights *= aper_mask
        self.dark_current_weights *= aper_mask
        # (Do NOT include fractional pixel weights for red leak weights to prevent
        # double-counting of fractional pixel weights in `_calc_redleaks()`)
        self.redleak_weights *= np.ma.masked_where(
            np.isfinite(aper_mask), aper_mask, copy=True
        ).filled(1.0)
        self._eff_npix = np.nansum(aper_mask)
        #
        # Find the fraction of the source within the aperture
        #
        self._do_source_aper_overlap(source_center_px, px_scale_arcsec)
        #
        # Final sanity check
        #
        if abs(self._eff_npix - self._exact_npix) > 0.1:
            # Discrepancy larger than 0.1 pixels
            warnings.warn(
                "Effective aperture area is off by more than 0.1 pixels... "
                + "Contact the developer with a minimal working example please. Thanks!",
                RuntimeWarning,
            )

    def use_elliptical_aperture(
        self, a, b, center=[0, 0] << u.arcsec, rotation=0, overwrite=False
    ):
        """
        Use an elliptical aperture.

        Parameters
        ----------
          a, b :: int or float or `astropy.Quantity` angle
            The angular length of the semimajor and semiminor axes of the aperture,
            respectively. If int or float, a/b is assumed to be in pixel units.

          center :: 2-element `astropy.Quantity` angles array
            The (x, y) center of the aperture relative to the center of the source.
            Positive values means the source is displaced to the bottom/left relative to
            the aperture center.

          rotation :: int or float
            The counter-clockwise rotation angle in degrees of the ellipse's semimajor
            axis from the positive x-axis. If rotation is 0, the semimajor axis is along
            the x-axis and the semiminor axis is along the y-axis.

          overwrite :: bool
            If True, allow overwriting of any existing aperture associated with this
            `Photometry` object.

        Attributes
        ----------
          source_weights :: (M x N) 2D array of floats
            The pixel weights for the source.

          sky_background_weights :: (M x N) 2D array of floats
            The pixel weights for the sky background.

          dark_current_weights :: (M x N) 2D array of floats
            The pixel weights for the dark current.

          redleak_weights :: (M x N) 2D array of floats
            The pixel weights for the red leak.

        Returns
        -------
          None
        """
        #
        # Check inputs and
        #
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
        #
        # Calculate exact aperture area and number of pixels in aperture
        #
        self._aper_area = np.pi * a * b  # area of ellipse
        self._assign_exact_npix()
        #
        # Create source weights with arbitrary source flux profile through aperture
        #
        # Recall all internal angle aperture angles are in arcsec
        center = center.to(u.arcsec).value
        px_scale_arcsec = self.TelescopeObj.px_scale.to(u.arcsec).value
        # Below is modified rotation matrix to ensure full aperture is covered in arrays
        # N.B. must round to nearest multiple of px_scale or else rounding errors will
        # affect source weights
        if abs(a.value - b.value) >= 1e-15:
            # Non-circular aperture
            rotation = np.deg2rad(rotation)
            abs_sin_rotate, abs_cos_rotate = abs(np.sin(rotation)), abs(np.cos(rotation))
            x = (
                np.ceil((abs_cos_rotate * a + abs_sin_rotate * b) / px_scale_arcsec)
                * px_scale_arcsec
            ).value
            y = (
                np.ceil((abs_sin_rotate * a + abs_cos_rotate * b) / px_scale_arcsec)
                * px_scale_arcsec
            ).value
            x = x if x >= a.value else a.value
            y = y if y >= b.value else b.value
        else:
            # Circular aperture
            x = a.value
            y = b.value
        #
        center_px, source_center_px = self._create_aper_arrs(
            x, y, center, overwrite=overwrite
        )
        source_weights = Photometry._calc_source_weights(
            self.SourceObj.profile, self._aper_xs, self._aper_ys, center
        )
        #
        # Create aperture
        #
        aper_dimen_px = [(a / px_scale_arcsec).value, (b / px_scale_arcsec).value]
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
        # (Do NOT include fractional pixel weights for red leak weights to prevent
        # double-counting of fractional pixel weights in `_calc_redleaks()`)
        self.redleak_weights *= np.ma.masked_where(
            np.isfinite(aper_mask), aper_mask, copy=True
        ).filled(1.0)
        self._eff_npix = np.nansum(aper_mask)
        #
        # Find the fraction of the source within the aperture
        #
        self._do_source_aper_overlap(source_center_px, px_scale_arcsec)
        #
        # Final sanity checks
        #
        if abs(self._eff_npix - self._exact_npix) > 0.1:
            # Discrepancy larger than 0.1 pixels
            warnings.warn(
                "Effective aperture area is off by more than 0.1 pixels... "
                + "Contact the developer with a minimal working example please. Thanks!",
                RuntimeWarning,
            )
        if isinstance(self.SourceObj, PointSource):
            if a < self._optimal_aperture_radius or b < self._optimal_aperture_radius:
                warnings.warn(
                    "Chosen a/b is smaller than the 'ideal' aperture size "
                    + "for this source! a & b should be at least "
                    + f"{self._optimal_aperture_radius}.",
                    UserWarning,
                )

    def use_rectangular_aperture(
        self, width, length, center=[0, 0] << u.arcsec, overwrite=False
    ):
        """
        Use a rectangular aperture.

        Parameters
        ----------
          width, length :: int or float or `astropy.Quantity` angle
            The width (along the x-axis) and length (along the y-axis) of the rectangular
            aperture. If int or float, length/width is assumed to be in pixel units.

          center :: 2-element `astropy.Quantity` angles array
            The (x, y) center of the aperture relative to the center of the source.
            Positive values means the source is displaced to the bottom/left relative to
            the aperture center.

          overwrite :: bool
            If True, allow overwriting of any existing aperture associated with this
            `Photometry` object.

        Attributes
        ----------
          source_weights :: (M x N) 2D array of floats
            The pixel weights for the source.

          sky_background_weights :: (M x N) 2D array of floats
            The pixel weights for the sky background.

          dark_current_weights :: (M x N) 2D array of floats
            The pixel weights for the dark current.

          redleak_weights :: (M x N) 2D array of floats
            The pixel weights for the red leak.

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
        # Calculate exact aperture area and number of pixels in aperture
        #
        self._aper_area = length * width
        self._assign_exact_npix()
        #
        # Create source weights with arbitrary source flux profile through aperture
        #
        # Recall all internal angle aperture angles are in arcsec
        center = center.to(u.arcsec).value
        center_px, source_center_px = self._create_aper_arrs(
            0.5 * width.value,  # width is along x
            0.5 * length.value,  # length is along y
            center,
            overwrite=overwrite,
        )
        source_weights = Photometry._calc_source_weights(
            self.SourceObj.profile, self._aper_xs, self._aper_ys, center
        )
        #
        # Create aperture
        #
        px_scale_arcsec = self.TelescopeObj.px_scale.to(u.arcsec).value
        aper_dimen_px = [
            (width / px_scale_arcsec).value,
            (length / px_scale_arcsec).value,
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
        # (Do NOT include fractional pixel weights for red leak weights to prevent
        # double-counting of fractional pixel weights in `_calc_redleaks()`)
        self.redleak_weights *= np.ma.masked_where(
            np.isfinite(aper_mask), aper_mask, copy=True
        ).filled(1.0)
        self._eff_npix = np.nansum(aper_mask)
        #
        # Find the fraction of the source within the aperture
        #
        self._do_source_aper_overlap(source_center_px, px_scale_arcsec)
        #
        # Final sanity checks
        #
        if abs(self._eff_npix - self._exact_npix) > 0.1:
            # Discrepancy larger than 0.1 pixels
            #
            # BUG in photutils (as of v1.4.0)?:
            # Try making an aperture with
            # `width=2.12 * u.arcsec, length=1.8 * u.arcsec, center=[0, 0] * u.arcsec`
            # and the `Telescope` pixel size is `0.1 * arcsec`... _eff_npix will be wrong
            # regardless of the aper array sizes (i.e., not a padding problem).
            # TODO: report this bug to photutils...
            #
            warnings.warn(
                "Effective aperture area is off by more than 0.1 pixels... "
                + "Contact the developer with a minimal working example please. Thanks!"
                + "\nNOTE: As of photutils-v1.4.0, there seems to be a bug where, in some "
                + "cases, the rectangular aperture mask will return the wrong number of "
                + "pixels regardless of the array size used to house the aperture "
                + "mask!! This requires a fix from the `photutils` team, unfortunately.",
                RuntimeWarning,
            )
            print("eff_npix, exact_npix", self._eff_npix, self._exact_npix)
        if isinstance(self.SourceObj, PointSource):
            aper_threshold = self._optimal_aperture_radius * 2
            if (length < aper_threshold) or (width < aper_threshold):
                warnings.warn(
                    "Chosen length/width is smaller than the 'ideal' aperture "
                    + "size for this source! "
                    + f"Length/width should be at least {aper_threshold}.",
                    UserWarning,
                )

    def use_annular_aperture(self, a_out, b_out, a_in=0, b_in=None):
        """
        Not implemented!

        TODO: maybe implement this if enough demand?
        """
        if b_in is None:
            b_in = b_out * (a_in / a_out)
        raise NotImplementedError("Not yet implemented")

    def calc_redleak_frac(self, quiet=False):
        """
        Calculate a source's red leak fraction. The red leak fraction is defined to be the
        ratio of the electron rate (i.e., electron/s) induced by red leak flux to the
        total electron rate induced by the entire spectrum.

        Parameters
        ----------
          quiet :: bool
            If True, suppress warnings from red leak fraction calculations.

        Returns
        -------
          redleak_fracs :: dict of floats
            Dictionary containing the red leak fraction in each passband.
        """
        #
        # Calculate red leak fraction (red leak electron/s to total electron/s)
        #
        redleak_fracs = dict.fromkeys(self.TelescopeObj.passbands, 0.0)
        #
        # Make useful source spectrum-derived quantities
        #
        source_wavelengths_AA = self.SourceObj.wavelengths.to(u.AA).value
        source_photon_s_A = (  # photon/s/A
            self.SourceObj.spectrum  # erg/s/cm^2/A
            * self.TelescopeObj.mirror_area.to(u.cm ** 2).value  # cm^2
            / calc_photon_energy(wavelength=source_wavelengths_AA)[0]  # photon/erg
        )
        source_interp = interp1d(
            source_wavelengths_AA,
            source_photon_s_A,
            kind="linear",
            bounds_error=False,
            fill_value=np.nan,
        )  # photon/s/A
        #
        # Find red leak fraction per band
        #
        for band in redleak_fracs:
            full_response_curve_wavelengths_AA = (
                self.TelescopeObj.full_passband_curves[band]["wavelength"].to(u.AA).value
            )
            is_redleak = (
                full_response_curve_wavelengths_AA
                > self.TelescopeObj.redleak_thresholds[band].to(u.AA).value
            )
            redleak_wavelengths = full_response_curve_wavelengths_AA[is_redleak]
            redleak_per_A = (
                source_interp(redleak_wavelengths)
                * self.TelescopeObj.full_passband_curves[band]["response"][is_redleak]
            )  # electron/s/A
            total_erate_per_A = (
                source_interp(full_response_curve_wavelengths_AA)
                * self.TelescopeObj.full_passband_curves[band]["response"]
            )  # electron/s/A
            isgood_redleak = np.isfinite(redleak_per_A)  # don't include NaNs
            isgood_total = np.isfinite(total_erate_per_A)  # don't include NaNs
            if not quiet and (not np.all(isgood_redleak) or not np.all(isgood_total)):
                warnings.warn(
                    "Could not estimate red leak fraction "
                    + f"at 1 or more wavelengths in {band}-band. "
                    + "This may just be caused by the source spectrum not being "
                    + f"defined at all wavelengths present in the {band}-band definition "
                    + f"file (which goes up to {round(max(redleak_wavelengths), 2)} A) "
                    + "and is typically not a reason to worry. "
                    + "This warning can be suppressed with `quiet=True`.",
                    RuntimeWarning,
                )
            try:
                redleak_per_px = simpson(  # electron/s (per px)
                    y=redleak_per_A[isgood_redleak],
                    x=redleak_wavelengths[isgood_redleak],
                    even="avg",
                )
                total_erate_per_px = simpson(  # electron/s (per px)
                    y=total_erate_per_A[isgood_total],
                    x=full_response_curve_wavelengths_AA[isgood_total],
                    even="avg",
                )
                redleak_frac = redleak_per_px / total_erate_per_px
            except Exception:
                raise RuntimeError(
                    f"Unable to calculate red leak fraction for {band}-band! "
                    + "Please ensure there is at least 1 wavelength that is above "
                    + "the red leak threshold."
                )
            if np.isfinite(redleak_frac):
                redleak_fracs[band] = redleak_frac
            else:
                raise RuntimeError(
                    "Source red leak fraction could not be calculated "
                    + f"in {band}-band!"
                )
        return redleak_fracs

    def _calc_redleaks(
        self, source_erate, mirror_area_cm_sq, include_redleak=True, quiet=False
    ):
        """
        Calculate the red leak of a source (in electron/s) from its in-passband red leak
        fraction. The in-passband red leak fraction is defined to be the ratio of red leak
        electron/s to in-passband electron/s.

        This is a robust way to calculate red leaks that is independent of how
        source_erate is determined and independent of any normalizations. It simply scales
        each pixel's electron rate (i.e., electron/s) in each passband by its respective
        in-passband red leak fraction to get the red leak electron rate.

        Parameters
        ----------
          source_erate :: dict of 2D arrays
            Dictionary containing the pixel-by-pixel electron rate (i.e., electron/s)
            induced by the source in each passband. In other words, this is the
            pixel-by-pixel signal in each passband.

          mirror_area_cm_sq :: float
            The area of the telescope's mirror in square centimeters.

          include_redleak :: bool
            If False, do not include red leak (i.e., the returned dictionary will be
            zeroes).

          quiet :: bool
            If True, suppress warnings from in-passband red leak fraction calculations.

        Returns
        -------
          redleaks :: dict of 2D arrays
            Dictionary containing the pixel-by-pixel red leak values (in electron/s).
        """
        #
        # Calculate in-passband redleak fraction (red leak electron/s to passband
        # electron/s)
        #
        redleak_fracs = dict.fromkeys(self.TelescopeObj.passbands, 0.0)
        if include_redleak:
            #
            # Make useful source spectrum-derived quantities
            #
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
            #
            # Find in-passband redleak fraction per band
            #
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
                if not quiet and (
                    not np.all(isgood_redleak) or not np.all(isgood_passband)
                ):
                    warnings.warn(
                        "Could not estimate in-passband red leak fraction "
                        + f"at 1 or more wavelengths in {band}-band. "
                        + "This may just be caused by the source spectrum not being "
                        + f"defined at all wavelengths present in the {band}-band "
                        + "definition file (which goes up to "
                        + f"{round(max(redleak_wavelengths), 2)} A) "
                        + "and is typically not a reason to worry. "
                        + "This warning can be suppressed with `quiet=True`.",
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
                        "Unable to calculate in-passband red leak fraction for "
                        + f"{band}-band! Please ensure there is at least 1 wavelength "
                        + "that is above the red leak threshold."
                    )
                if np.isfinite(redleak_frac):
                    # Important: note that redleak_weights should not have fractional
                    # pixel weighting based on the aperture mask or else the aperture
                    # weights will be double-counted when the redleak_fracs arrays are
                    # multiplied by the source_erates...
                    redleak_fracs[band] = redleak_frac * self.redleak_weights
                else:
                    raise RuntimeError(
                        "Source in-passband red leak fraction could not be calculated "
                        + f"in {band}-band!"
                    )
        #
        # Calculate red leak from in-passband red leak fraction
        #
        redleaks = dict.fromkeys(self.TelescopeObj.passbands)
        for band in redleaks:
            # (Recall source_erate includes fractional pixel weighting from aperture mask)
            redleaks[band] = redleak_fracs[band] * source_erate[band]
        #
        return redleaks

    @staticmethod
    def _calc_snr_from_t(
        t,
        signal,
        totskynoise,
        darkcurrent,
        redleak,
        readnoise,
        read_npix,
        nread=1,
    ):
        """
        Calculate the signal-to-noise ratio (SNR) reached given an integration time.

        The equation to calculate the SNR is:
        ```math
                    SNR = (Q*t) / sqrt(Q*t + N_pix*t*Other_Poisson + N_pix*N_read*Read^2)
        ```
        where:
          - SNR is the signal-to-noise ratio
          - t is the integration time in seconds
          - Q is the total signal due to the source in electrons/s
          - N_pix is the number of pixels occupied by the source on the detector
          - Other_Poisson = (B_Earthshine + B_zodi + B_geocoronal + Darkcurrent + Redleak)
          is the total Poisson noise due to the sky backgorund (Earthshine, zodiacal
          light, and geocoronal emission), dark current, and red leak in electrons/s (per
          pixel)
          - N_read is the number of detector readouts
          - Read is the detector read noise in electrons (per pixel)
        The equation above is based on the formula detailed here:
        <https://hst-docs.stsci.edu/wfc3ihb/chapter-9-wfc3-exposure-time-calculation/9-6-estimating-exposure-times>.

        In practice, the calculation takes in 2D arrays for the source signal, sky
        background, dark current, and red leak. The factor of N_pix is already included
        for these arrays since they have fractional pixel weighting. Therefore, simply
        summing up these arrays (excluding NaNs) is enough to calculate the total
        contribution from that element without further multiplying by N_pix. The only
        exception is N_read and the read noise, which are scalars and must be multiplied
        by an integer number of pixels.


        Parameters
        ----------
          t :: int or float
            The integration time in seconds.

          signal :: (M x N) 2D array of scalars
            2D array giving the signal of the source (in electron/s) for each pixel in the
            aperture. This should include fractional pixel weighting from the aperture
            mask.

          totskynoise :: (M x N) 2D array of scalars
            2D array giving the total background noise due to Earthshine + zodiacal light
            + geocoronal emission (if applicable) in electron/s for each pixel in the
            aperture. This should include fractional pixel weighting from the aperture
            mask.

          darkcurrent :: (M x N) 2D array of scalars
            2D array giving the detector dark current in electron/s for each pixel in the
            aperture. This should include fractional pixel weighting from the aperture
            mask.

          darkcurrent :: (M x N) 2D array of scalars
            2D array giving the red leak in electron/s for each pixel in the
            aperture. This should include fractional pixel weighting from the aperture
            mask.

          readnoise :: int or float
            The detector read noise in electron (per pixel).

          read_npix :: int
            The integer number of pixels in the aperture. Only used for multiplying with
            the read noise and number of detector readouts.

          nread :: int
            The number of detector readouts.

        Returns
        -------
          snr :: float
            The SNR reached after t seconds.
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
        if not isinstance(read_npix, (int, np.integer)):
            raise ValueError("read_npix must be an integer")
        if not isinstance(nread, (int, np.integer)):
            raise ValueError("nread must be an integer")
        #
        # Calculate signal-to-noise ratio
        #
        signal_t = np.nansum(signal * t)  # electron
        noise = np.sqrt(
            signal_t
            + np.nansum(t * (totskynoise + darkcurrent + redleak))
            + (read_npix * readnoise * readnoise * nread)
        )  # electron
        snr = signal_t / noise
        return snr

    @staticmethod
    def _calc_t_from_snr(
        snr,
        signal,
        totskynoise,
        darkcurrent,
        redleak,
        readnoise,
        read_npix,
        nread=1,
    ):
        """
        Calculate the time required to reach a given signal-to-noise ratio (SNR).

        The theoretical equation to calculate the time required to reach a given SNR is:
        ```math
                    t = {
                        SNR^2 * (Q + N_pix * Other_Poisson)
                        + sqrt[
                            SNR^4 * (Q + N_pix * Other_Poisson)^2
                            + 4 * Q^2 * SNR^2 * N_pix * N_read * Read^2
                        ]
                    } / (2 * Q^2)
        ```
        where:
          - SNR is the desired signal-to-noise ratio
          - t is the integration time in seconds
          - Q is the total signal due to the source in electrons/s
          - N_pix is the number of pixels occupied by the source on the detector
          - Other_Poisson = (B_Earthshine + B_zodi + B_geocoronal + Darkcurrent + Redleak)
          is the total Poisson noise due to the sky backgorund (Earthshine, zodiacal
          light, and geocoronal emission), dark current, and red leak in electrons/s (per
          pixel)
          - N_read is the number of detector readouts
          - Read is the detector read noise in electrons (per pixel)
        Note that this is simply the quadratic formula applied to the generic SNR
        equation. Also, the equation above is based on the formula detailed here:
        <https://hst-docs.stsci.edu/wfc3ihb/chapter-9-wfc3-exposure-time-calculation/9-6-estimating-exposure-times>.

        In practice, the calculation takes in 2D arrays for the source signal, sky
        background, dark current, and red leak. The factor of N_pix is already included
        for these arrays since they have fractional pixel weighting. Therefore, simply
        summing up these arrays (excluding NaNs) is enough to calculate the total
        contribution from that element without further multiplying by N_pix. The only
        exception is N_read and the read noise, which are scalars and must be multiplied
        by an integer number of pixels.

        Parameters
        ----------
          snr :: int or float
            The target SNR.

          signal :: (M x N) 2D array of scalars
            2D array giving the signal of the source (in electron/s) for each pixel in the
            aperture. This should include fractional pixel weighting from the aperture
            mask.

          totskynoise :: (M x N) 2D array of scalars
            2D array giving the total background noise due to Earthshine + zodiacal light
            + geocoronal emission (if applicable) in electron/s for each pixel in the
            aperture. This should include fractional pixel weighting from the aperture
            mask.

          darkcurrent :: (M x N) 2D array of scalars
            2D array giving the detector dark current in electron/s for each pixel in the
            aperture. This should include fractional pixel weighting from the aperture
            mask.

          darkcurrent :: (M x N) 2D array of scalars
            2D array giving the red leak in electron/s for each pixel in the
            aperture. This should include fractional pixel weighting from the aperture
            mask.

          readnoise :: int or float
            The detector read noise in electron (per pixel).

          read_npix :: int
            The integer number of pixels in the aperture. Only used for multiplying with
            the read noise and number of detector readouts.

          nread :: int
            The number of detector readouts.

        Returns
        -------
          t :: float
            The time required to reach the given SNR in seconds.
        """
        variables = [snr, signal, totskynoise, readnoise, darkcurrent, redleak, nread]
        if np.any([isinstance(var, u.Quantity) for var in variables]):
            raise ValueError(
                "All inputs must be scalars or scalar arrays. "
                + "`astropy.Quantity` objects are not supported."
            )
        if not isinstance(read_npix, (int, np.integer)):
            raise ValueError("read_npix must be an integer")
        if not isinstance(nread, (int, np.integer)):
            raise ValueError("nread must be an integer")
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
        numer3 = 4 * snr_sq * signal_sq * (read_npix * readnoise * readnoise * nread)
        t = (numer1 + np.sqrt(numer2 + numer3)) / (2 * signal_sq)  # seconds
        return t

    def calc_snr_or_t(
        self,
        t=None,
        snr=None,
        reddening=0,
        encircled_energy=None,
        npix=None,
        nread=1,
        include_redleak=True,
        quiet=False,
    ):
        """
        Calculate the signal-to-noise ratio (SNR) reached in a given time or the
        integration time (t) required to reach a given SNR.

        Again, note that the `Source` object associated with the `Photometry` instance
        should have its spectrum in units of flam (erg/s/cm^2/A).

        Parameters
        ----------
          t :: float
            Integration time in seconds.

          snr :: float
            Desired signal-to-noise ratio.

          reddening :: int or float
            The reddening (i.e., E(B-V) in AB mags) based on the pointing. This value will
            be multiplied with each of the telescope's extinction coefficients to get the
            total extinction in each passband.

          encircled_energy :: int or float or None
            The fraction of the source flux enclosed within the aperture. If None, use the
            fraction of the source'a area enclosed within the aperture, which is a value
            between (0, 1]; this value can be seen on the source weights plot (i.e., via
            the `show_source_weights()` method).

          npix :: int or float or None
            The number of pixels occupied by the source to use for the noise calculations
            (i.e., affects sky background, dark current, and read noise). If None, use the
            automatically calculated "effective pixel" count. Note that the number of
            pixels used for the read noise will always be npix rounded up to the nearest
            integer. Finally, this value does not affect the total signal of the
            source---only the non-signal-derived noise values are affected (e.g., higher
            npix means higher total sky background, dark current, and read noise but not
            higher signal nor redleak).

          nread :: int
            Number of detector readouts.

          include_redleak :: bool
            If False, the redleak contribution to the total noise is not included. May be
            useful if the source spectrum is, for example, uniform.

          quiet :: bool
            If True, suppress warning messages from any red leak calculations. Only
            matters if include_redleak is True.

        Returns
        -------
          results :: dict
            If t is given, this is the snr reached after t seconds. If snr is given, this
            is the time in seconds required to reach the given snr.
        """
        #
        # Check some inputs (npix and nread will be checked later, in other functions)
        #
        if (t is None and snr is None) or (t is not None and snr is not None):
            raise ValueError("Exactly one of t or snr must be specified")
        if not isinstance(reddening, Number):
            raise TypeError("reddening must be an int or float")
        if self.source_weights is None:
            raise ValueError("Please choose an aperture first.")
        if encircled_energy is not None:
            if (
                not isinstance(encircled_energy, Number)
                or encircled_energy > 1
                or encircled_energy <= 0
            ):
                raise ValueError("encircled_energy must be a number between (0, 1]")
        #
        # Make some useful variables
        #
        response_curve_wavelengths_AA = dict.fromkeys(self.TelescopeObj.passbands)
        for band in response_curve_wavelengths_AA:
            response_curve_wavelengths_AA[band] = (
                self.TelescopeObj.passband_curves[band]["wavelength"].to(u.AA).value
            )
        if npix is None:
            aper_weight_scale = 1.0  # don't scale aperture weights
            npix = self._eff_npix
        else:
            aper_weight_scale = npix / self._eff_npix  # scale aperture weights
        #
        # Calculate sky background electron/s
        # (incl. Earthshine & zodiacal light, excl. geocoronal emission lines)
        #
        sky_background_erate = dict.fromkeys(self.TelescopeObj.passbands, 0.0)
        px_area_arcsec_sq = self.TelescopeObj.px_area.to(u.arcsec ** 2).value
        mirror_area_cm_sq = self.TelescopeObj.mirror_area.to(u.cm ** 2).value
        # Use passband photometric zero points + sky background AB magnitudes (per sq.
        # arcsec) to calculate sky background electron/s.
        # Note that this is completely equivalent to convolving the sky background spectra
        # with passband response curves, which was the previous method (now removed
        # because this is much simpler)! Compare results if you want!
        if self.BackgroundObj.mags_per_sq_arcsec is None:
            background_mags_per_sq_arcsec = self.BackgroundObj._get_mags_per_sq_arcsec(
                self.TelescopeObj
            )
            for band in self.TelescopeObj.passbands:
                # Convert sky background AB mag per arcsec^2 to electron/s (per pixel)
                sky_background_erate[band] = (
                    mag_to_flux(
                        background_mags_per_sq_arcsec[band],
                        zpt=self.TelescopeObj.phot_zpts[band],
                    )[0]
                    * px_area_arcsec_sq
                )
        else:
            for band in self.TelescopeObj.passbands:
                try:
                    # Convert sky background AB mag per arcsec^2 to electron/s (per pixel)
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
        # REVIEW: remove requirement for linewidth? Unused here... Maybe need it for spectroscopy?
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
                        gf  # erg/cm^2/s/arcsec^2
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
            sky_background_erate[band] *= self.sky_background_weights * aper_weight_scale
        #
        # Calculate signal in each passband (flam -> electron/s)
        #
        source_erate = dict.fromkeys(self.TelescopeObj.passbands, np.nan)
        source_ab_mags = self.SourceObj.get_AB_mag(self.TelescopeObj)
        # (N.B. Unlike background noise, aper_weight_scale and npix cancel out for the
        # source electron/s calculation below. Thus, just use the original _eff_npix and
        # source_weights below; no need for npix or aper_weight_scale)
        if encircled_energy is None:
            encircled_energy = self._source_aper_overlap_frac
        for band in source_erate:
            # Account for extinction
            source_passband_mag = source_ab_mags[band] + (
                self.TelescopeObj.extinction_coeffs[band] * reddening
            )
            # Convert extinction-corrected AB mag to electron/s using passband's
            # photometric zero point and account for fraction of source within aperture
            # REVIEW: is this calculation correct???
            source_erate[band] = (
                mag_to_flux(
                    source_passband_mag,
                    zpt=self.TelescopeObj.phot_zpts[band],
                )[0]
                * (self.source_weights / self._eff_npix)  # must norm by _eff_npix
                * encircled_energy  # either user-provided or the source-aperture overlap
            )
            # NOTE: normalizing by _eff_npix and multiplying by the source-aperture
            # fractional overlap is actually equivalent to the surface brightness approach
            # (i.e., dividing by the area of the source and multiplying by the pixel area
            # instead of dividing by _eff_npix and multiplying by the source-aperture
            # fractional overlap) for extended sources! (Try it if you don't believe me.)
            # There are 2 important things to note, however. First, since the "surface
            # brightness" of our simulated source is based on the user-supplied
            # dimensions, the surface brightness approach will differ from our current
            # approach above when the area of the aperture is larger than the source. In
            # this case, I believe that our approach above makes more sense since the
            # surface brightness approach would lead to a source that increases in total
            # brightness as the aperture increases! Second, the approach above works for
            # point sources as well; further, the user can also specify the precise
            # encircled energy fraction to replace the source-aperture area overlap
            # approximation.
        # # See Eq. (2) and Eq. (9) of
        # # <https://hst-docs.stsci.edu/acsihb/chapter-9-exposure-time-calculations/9-2-determining-count-rates-from-sensitivities#id-9.2DeterminingCountRatesfromSensitivities-9.2.1>.
        # if isinstance(self.SourceObj, PointSource):
        #     for band in source_erate:
        #         # Account for extinction
        #         source_passband_mag = source_ab_mags[band] + (
        #             self.TelescopeObj.extinction_coeffs[band] * reddening
        #         )
        #         # Convert extinction-corrected AB mag to electron/s using Eq. (2) from the
        #         # link above and the passband's photometric zero point
        #         # REVIEW: is this calculation correct???
        #         source_erate[band] = (
        #             mag_to_flux(
        #                 source_passband_mag,
        #                 zpt=self.TelescopeObj.phot_zpts[band],
        #             )[0]
        #             * (self.source_weights / self._eff_npix)  # must norm by _eff_npix
        #             * self._source_aper_overlap_frac
        #         )
        # else:
        #     for band in source_erate:
        #         # Account for extinction
        #         source_passband_mag = source_ab_mags[band] + (
        #             self.TelescopeObj.extinction_coeffs[band] * reddening
        #         )
        #         # Convert extinction-corrected AB mag to electron/s using Eq. (9) from the
        #         # link above and the passband's photometric zero point
        #         # REVIEW: is this calculation correct???
        #         # ! The following calculation doesn't make sense! Try having an aperture > source size...
        #         # source_erate[band] = (
        #         #     mag_to_flux(
        #         #         source_passband_mag,
        #         #         zpt=self.TelescopeObj.phot_zpts[band],
        #         #     )[0]
        #         #     * self.source_weights
        #         #     * (
        #         #         px_area_arcsec_sq / self.SourceObj.area.to(u.arcsec ** 2).value
        #         #     )  # convert to surface brightness
        #         # )
        #
        # Calculate red leak (electron/s) at each pixel
        #
        redleaks = self._calc_redleaks(
            source_erate, mirror_area_cm_sq, include_redleak=include_redleak, quiet=quiet
        )
        #
        # Calculate desired results (either integration time given SNR or SNR given time)
        #
        results = dict.fromkeys(self.TelescopeObj.passbands)
        dark_current = (
            self.TelescopeObj.dark_current * self.dark_current_weights * aper_weight_scale
        )
        if t is not None:
            for band in results:
                results[band] = Photometry._calc_snr_from_t(
                    t=t,
                    signal=source_erate[band],
                    totskynoise=sky_background_erate[band],
                    darkcurrent=dark_current,
                    redleak=redleaks[band],
                    readnoise=self.TelescopeObj.read_noise,  # constant per pixel
                    read_npix=int(np.ceil(npix)),  # only used for readnoise and nread
                    nread=nread,
                )
        else:
            for band in results:
                results[band] = Photometry._calc_t_from_snr(
                    snr=snr,
                    signal=source_erate[band],
                    totskynoise=sky_background_erate[band],
                    darkcurrent=dark_current,
                    redleak=redleaks[band],
                    readnoise=self.TelescopeObj.read_noise,  # constant per pixel
                    read_npix=int(np.ceil(npix)),  # only used for readnoise and nread
                    nread=nread,
                )
        return results
