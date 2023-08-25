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
from scipy.interpolate import interp1d
from scipy.signal import oaconvolve

from .background import Background
from .conversions import calc_photon_energy, mag_to_flux
from .sources import CustomSource, ExtendedSource, GalaxySource, PointSource, Source
from .telescope import Telescope

# The optimal aperture for a point source is a circular aperture with a radius equal
# to the factor below times half the telescope's FWHM
_OPTIMAL_APER_FACTOR = 1.4
# The supersampling factor along each axis
_SUPERSAMPLE_FACTOR = 20


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
        if (
            isinstance(SourceObj, CustomSource)
            and SourceObj.passband not in TelescopeObj.passbands
        ):
            raise ValueError(
                "The `CustomSource` object's surface brightness profile passband "
                + f"('{SourceObj.passband}') is not a valid `TelescopeObj` "
                + f"passband ({TelescopeObj.passbands})"
            )
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
        self.source_weights = dict.fromkeys(TelescopeObj.passbands, None)
        self.source_weights["noiseless"] = None  # source weights before PSF convolution
        self.sky_background_weights = None
        self.dark_current_weights = None
        # Attributes for internal use
        self._aper_area = None  # exact area of the aperture from given aperture params
        self._aper_xs = None  # array containing x-coordinates of pixels
        self._aper_ys = None  # array containing y-coordinates of pixels
        self._aper_mask = None  # photutils aperture mask
        self._exact_npix = None  # number of pixels calculated from given aperture params
        self._eff_npix = None  # number of pixels calculated from photutils mask
        self._aper_extent = None  # the [xmin, xmax, ymin, ymax] extent of the imshow plot
        self._xdim = None  # the number of pixels along the x-direction of the aperture
        self._ydim = None  # the number of pixels along the y-direction of the aperture
        # Default optimal aperture dimensions for a point source. For error-checking
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
            self._aper_area.to(u.arcsec**2)
            / self.TelescopeObj.px_area.to(u.arcsec**2)
        ).value

    def _create_aper_arrs(
        self,
        half_x,
        half_y,
        center,
        overwrite=False,
    ):
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
            relative to the aperture center. This is because the source will always be at
            (0, 0) and the aperture center will be at this `center` value.

          overwrite :: bool
            If True, allow overwriting of any existing aperture arrays.

        Attributes
        ----------
          _xdim, _ydim :: ints
            The true, non-supersampled number of pixels along the x- and y-directions.
            This should be the dimensions of the final arrays after binning.

          _aper_xs :: (M x N) 2D array of floats
            The aperture array containing the x-coordinates (in arcsec) relative to the
            center of the aperture. Includes supersampling.

          _aper_ys :: (M x N) 2D array of floats
            The aperture array containing the y-coordinates (in arcsec) relative to the
            center of the aperture. Includes supersampling.

          _aper_extent :: 4-element 1D list of floats
            The [xmin, xmax, ymin, ymax] extent of the aperture in arcsec (for plotting
            the weight arrays). Includes supersampling.

          sky_background_weights :: (M x N) 2D array of floats
            The pixel weights for the sky background (incl. Earthshine, zodiacal light,
            geocoronal emission). This is currently all ones (1) (i.e., uniform noise).
            Includes supersampling.
            The sky background noise is given per sq. arcsec, so a 1 sq. arcsec pixel with
            a sky background weight of 1 means the pixel receives the full sky background
            noise (for that pixel), while a sky background weight of 0.4 means the pixel
            only receives 40% of the full sky background noise for that pixel.

          dark_current_weights :: (M x N) 2D array of floats
            The pixel weights for the dark current. This is currently all ones (1) (i.e.,
            uniform noise). Includes supersampling.
            A dark current weight of 1 means the pixel experiences the full dark current
            rate (per pixel), while a dark current weight of 0.4 means the dark current
            noise for that pixel is 40% of the telescope's dark current value.

        Returns
        -------
          center_px :: 2-element 1D list of floats
            The (x, y) center index of the aperture coordinates arrays (i.e., _aper_xs and
            _aper_ys). To be very explicit, this is not (0, 0) but rather the center of
            the 2D arrays.
        """
        #
        # Check inputs
        #
        if not overwrite and (self._aper_xs is not None or self._aper_ys is not None):
            raise ValueError(
                "An aperture for this `Photometry` object already exists. "
                + "Use `overwrite=True` to allow overwriting of the aperture "
                + " and all associated weights (i.e., source, sky background, "
                + "and dark current weights will all be reset)."
            )
        #
        px_scale_arcsec = self.TelescopeObj.px_scale.to(u.arcsec).value
        half_px_scale = 0.5 * px_scale_arcsec  # to ensure correct arange/extent
        #
        # Ensure full aperture is contained within array, including making sure extent is
        # possible to tile with pixels. Aperture mask will take care of pixel weighting
        #
        # The true number of pixels along the x- and y-directions
        self._xdim = np.arange(-half_x, half_x + half_px_scale, px_scale_arcsec).size
        self._ydim = np.arange(-half_y, half_y + half_px_scale, px_scale_arcsec).size
        #
        # Supersample aperture at PSF supersampling resolution for binning down to
        # (y, x) = (self._ydim, self._xdim).
        # (The binning will be handled by the `_bin_arrs_remove_nans()` method, which is
        # called in the various `use_<???>_aperture()` methods)
        #
        xs = np.linspace(
            -half_x, half_x, self._xdim * self.TelescopeObj.psf_supersample_factor
        )  # length N
        ys = np.linspace(
            -half_y, half_y, self._ydim * self.TelescopeObj.psf_supersample_factor
        )  # length M
        #
        self._aper_xs, self._aper_ys = np.meshgrid(xs, ys, sparse=False, indexing="xy")
        self._aper_extent = [
            xs[0] + center[0],
            xs[-1] + center[0],
            ys[0] + center[1],
            ys[-1] + center[1],
        ]
        #
        # Assume uniform background noise and dark current over photometry aperture. Need
        # to divide by the supersampling factor of the PSFs since each final pixel should
        # have a weight of 1
        #
        self.sky_background_weights = (
            np.ones_like(self._aper_xs) / self.TelescopeObj.psf_supersample_factor**2
        )
        self.dark_current_weights = (
            np.ones_like(self._aper_xs) / self.TelescopeObj.psf_supersample_factor**2
        )
        #
        # Find pixel that corresponds to center of the source
        #
        center_px = [
            0.5 * (self._aper_xs.shape[1] - 1),  # x-coordinate in pixel units
            0.5 * (self._aper_ys.shape[0] - 1),  # y-coordinate in pixel units
        ]

        return center_px

    def _calc_source_weights(self, center):
        """
        Calculate the source weights (before PSF convolution) for the given profile.

        Parameters
        ----------
          center :: 2-element 1D array of floats
            The (x, y) center of the aperture relative to the center of the source in
            arcsec. Positive values means the source is displaced to the bottom/left
            relative to the aperture center. This is because the source will always be at
            (0, 0) and the aperture center will be at this `center` value.

        Returns
        -------
          source_weights :: (M x N) 2D array of floats
            The source weights for each pixel in the aperture. These (currently) represent
            the flux of the source at each pixel relative to the flux at the center of the
            source. After convolution & normalization, the source weights will represent
            the fraction of the source's flux contained within the pixel.
        """
        if isinstance(self.SourceObj, PointSource):
            # Set closest pixel to center to 1 and all other pixels to 0
            # If multiple pixels are equally closest, set these to 1 / (num equally close)
            dists = np.sqrt(
                (self._aper_xs + center[0]) ** 2 + (self._aper_ys + center[1]) ** 2
            )
            min_dist = np.nanmin(dists)
            min_idxs = np.abs(dists - min_dist) < 1e-15  # float comparison
            source_weights = np.ones(self._aper_xs.shape) * min_idxs / np.nansum(min_idxs)
        else:
            source_weights = self.SourceObj.profile(self._aper_xs, self._aper_ys, center)
        return source_weights.astype(np.float64)

    def show_source_weights(
        self,
        passband,
        mark_source=False,
        source_markersize=4,
        norm=None,
        plot=True,
    ):
        """
        Plot the source as seen through the photometry aperture. The pixels are colored by
        the fraction of the source's flux contained within each pixel. Coloring also
        includes the effects of fractional pixel weights (i.e., from the aperture mask,
        visualized using `show_aper_weights()`). These two effects combined give the
        "source weights".

        Changing the source weights will affect the fraction of flux contained within the
        aperture and thus affect photometry calculations.

        Parameters
        ----------
          passband :: "noiseless" or one of the `TelescopeObj.passbands` keys
            If "noiseless", then the source weights are shown without any point spread
            function (PSF) convolution. If one of the `TelescopeObj.passbands` keys, then
            the source weights are shown after convolution with the passband's PSF.
            (Technically, these source weights are binned down from the supersampled
            version used for convolution with PSFs. That is, the true "noiseless" image is
            binned down from the PSF's oversampled resolution to produce this "noiseless"
            image.)

          mark_source :: bool
            If True, mark the center of the source with a cyan dot.

          source_markersize :: int or float
            The markersize for the cyan point indicating the center of the source.

          norm :: `matplotlib.colors` normalization class (e.g., `LogNorm`) or None
            The scaling and normalization to use for the colorbar. If None, then a linear
            scaling with the default (min pixel value, max pixel value) bounds are used.

          plot :: bool
            If True, plot the source weights and return None. If False, return the figure,
            axis, image, and colorbar instance associated with the plot.

        Returns
        -------
        If plot is True:
          None

        If plot is False:
          fig, ax :: `matplotlib.figure.Figure` and `matplotlib.axes.Axes` objects
            The figure and axis instance associated with the plot.

          img :: `matplotlib.image.AxesImage` object
            The image instance returned by `ax.imshow`.

          cbar :: `matplotlib.colorbar.Colorbar` object
            The colorbar instance associated with the plot.
        """
        try:
            source_weights = self.source_weights[passband]
        except KeyError:
            raise KeyError(f"{passband} is not 'noiseless' or a valid passband.")
        if source_weights is None or self._aper_extent is None:
            raise ValueError("Please select an aperture first.")
        #
        rc = {"axes.grid": False}
        with plt.rc_context(rc):  # matplotlib v3.5.x has bug affecting grid + imshow
            fig, ax = plt.subplots()
            # (N.B. array already "xy" indexing. Do not transpose array)
            img = ax.imshow(
                source_weights,
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
            if isinstance(self.SourceObj, CustomSource):
                cbar.set_label("Custom Surface Brightness (incl. fractional pixels)")
            else:
                cbar.set_label(
                    "Fraction of Flux Contained Within Pixel\n(incl. fractional pixels)"
                )
            ax.set_xlabel("$x$ [arcsec]")
            ax.set_ylabel("$y$ [arcsec]")
            percent = r"\%" if plt.rcParams["text.usetex"] else "%"
            if passband == "noiseless":
                ax.set_title(
                    "Fraction of Flux Contained Within Aperture:\n"
                    + f"{np.nansum(source_weights) * 100:.2f}{percent} (noiseless)"
                )
            else:
                ax.set_title(
                    "Fraction of Flux Contained Within Aperture:\n"
                    + f"{np.nansum(source_weights) * 100:.2f}{percent} ({passband}-band)"
                )
            if plot:
                plt.show()
            else:
                return fig, ax, img, cbar

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
            figure, axis, image, and colorbar instance associated with the plot.

        Returns
        -------
        If plot is True:
          None

        If plot is False:
          fig, ax :: `matplotlib.figure.Figure` and `matplotlib.axes.Axes` objects
            The figure and axis instance associated with the plot.

          img :: `matplotlib.image.AxesImage` object
            The image instance returned by `ax.imshow`.

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
            cbar.set_label("Fraction of Pixel Overlapping Aperture")
            ax.set_xlabel("$x$ [arcsec]")
            ax.set_ylabel("$y$ [arcsec]")
            ax.set_title(
                "Effective No." + "\u00A0" + f"of Aperture Pixels: {self._eff_npix:.2f}"
            )
            if plot:
                plt.show()
            else:
                return fig, ax, img, cbar

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
            The sky background noise is given per sq. arcsec, so a 1 sq. arcsec pixel with
            a sky background weight of 1 means the pixel receives the full sky background
            noise (for that pixel), while a sky background weight of 0.4 means the pixel
            only receives 40% of the full sky background noise for that pixel.

        Attributes
        -------
          sky_background_weights :: (M x N) 2D array of floats
            The pixel weights for the sky background. Note that the sky background
            includes Earthshine, zodiacal light, and geocoronal emission.

        Returns
        -------
          None
        """
        if None in self.source_weights.values():
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
            A dark current weight of 1 means the pixel experiences the full dark current
            rate (per pixel), while a dark current weight of 0.4 means the dark current
            noise for that pixel is 40% of the telescope's dark current value.

        Attributes
        -------
          dark_current_weights :: (M x N) 2D array of floats
            The pixel weights for the dark current.
        """
        if None in self.source_weights.values():
            raise ValueError("Please select an aperture first.")
        if dark_current_weights.shape != self.source_weights.shape:
            raise ValueError(
                "`dark_current_weights` must have same shape as `source_weights`."
            )
        self.dark_current_weights = dark_current_weights

    def _bin_arrs_remove_nans(self, center):
        """
        Bin arrays to self._xdim (N') and self._ydim (M'), as well as calculate the
        effective number of pixels within the aperture. This method also removes columns
        and rows of the aperture mask containing all NaNs.

        This method modifies the following attributes: `_eff_npix`, `_aper_mask`,
        `source_weights`, `sky_background_weights`, `dark_current_weights`, `_aper_xs`,
        `_aper_ys`. Also updates `_aper_extent` to reflect changed arrays.

        Parameters
        ----------
          center :: 2-element 1D array of floats
            The (x, y) center of the aperture relative to the center of the source in
            arcsec. Positive values means the source is displaced to the bottom/left
            relative to the aperture center. This is because the source will always be at
            (0, 0) and the aperture center will be at this `center` value.

        Attributes
        ----------
          _eff_npix :: float
            The effective number of pixels in the aperture.

          _aper_mask :: (M' x N') 2D array of floats
            The aperture mask, including fractional pixel weights. The aperture mask shows
            the fractional overlap between the aperture and the pixel, ranging from (0,
            1]. A pixel that is wholly contained in the aperture has an aperture weight
            equal to 1. Similarly, a pixel that is partially contained in the aperture has
            a weight between 0 and 1 (exclusive) directly proportional to the area of the
            aperture that overlaps the pixel. A pixel that is wholly outside the aperture
            has a weight of NaN. Now without rows and columns containing only NaNs.

          source_weights :: dict of (M' x N') 2D array of floats
            The pixel weights for the source representing the fraction of flux contained
            in each pixel, including fractional pixel weights. Now without rows and
            columns containing only NaNs (based on _aper_mask).

          sky_background_weights :: (M' x N') 2D array of floats
            The pixel weights for the sky background. Now without rows and columns
            containing only NaNs (based on _aper_mask).

          dark_current_weights :: (M' x N') 2D array of floats
            The pixel weights for the dark current. Now without rows and columns
            containing only NaNs (based on _aper_mask).

          _aper_xs :: (M' x N') 2D array of floats
            The aperture array containing the x-coordinates (in arcsec) relative to the
            center of the aperture. Now without rows and columns containing only NaNs
            (based on _aper_mask)

          _aper_ys :: (M' x N') 2D array of floats
            The aperture array containing the y-coordinates (in arcsec) relative to the
            center of the aperture. Now without rows and columns containing only NaNs
            (based on _aper_mask)

          _aper_extent :: 4-element 1D list of floats
            The [xmin, xmax, ymin, ymax] extent of the aperture in arcsec (for plotting
            the weight arrays). Updated to reflect the extent containing only rows and
            columns with at least 1 non-NaN value.

        Returns
        -------
          None
        """

        def _bin(arr):
            """
            Bin down 2D arrays, based on <https://stackoverflow.com/a/36102436>.
            """

            binned_arr = arr.reshape(
                self._ydim,
                self.TelescopeObj.psf_supersample_factor,
                self._xdim,
                self.TelescopeObj.psf_supersample_factor,
            )
            return np.nansum(np.nansum(binned_arr, axis=3), axis=1)

        #
        # Bin arrays
        #
        self._aper_mask = (
            _bin(self._aper_mask) / self.TelescopeObj.psf_supersample_factor**2
        )
        self._aper_mask[self._aper_mask < 1e-14] = np.nan
        for band in self.source_weights:
            self.source_weights[band] = _bin(self.source_weights[band])
        self.sky_background_weights = _bin(self.sky_background_weights)
        self.dark_current_weights = _bin(self.dark_current_weights)
        self._eff_npix = np.nansum(self._aper_mask)
        #
        # Ensure new aperture coordinate arrays reflect binned resolution
        #
        px_scale_arcsec = self.TelescopeObj.px_scale.to(u.arcsec).value
        xs = np.linspace(
            -0.5 * self._xdim * px_scale_arcsec,
            0.5 * self._xdim * px_scale_arcsec,
            self._xdim,
        )
        ys = np.linspace(
            -0.5 * self._ydim * px_scale_arcsec,
            0.5 * self._ydim * px_scale_arcsec,
            self._ydim,
        )
        self._aper_xs, self._aper_ys = np.meshgrid(xs, ys, sparse=False, indexing="xy")
        #
        # Remove columns that were originally all NaN
        #
        isgood_columns = ~(np.isnan(self._aper_mask).all(axis=0))
        self._aper_mask = self._aper_mask[:, isgood_columns]
        isgood_rows = ~(np.isnan(self._aper_mask).all(axis=1))
        self._aper_mask = self._aper_mask[isgood_rows, :]
        finite_mask = np.isfinite(self._aper_mask).astype(np.float64)
        finite_mask[np.abs(finite_mask) < 1e-14] = np.nan
        for band in self.source_weights:
            # Remove all NaN columns
            self.source_weights[band] = self.source_weights[band][:, isgood_columns]
            # Remove all NaN rows
            self.source_weights[band] = self.source_weights[band][isgood_rows, :]
            self.source_weights[band] *= finite_mask
        for arr, arr_name in zip(
            [
                self.sky_background_weights,
                self.dark_current_weights,
                self._aper_xs,
                self._aper_ys,
            ],
            [
                "sky_background_weights",
                "dark_current_weights",
                "_aper_xs",
                "_aper_ys",
            ],
        ):
            if arr is not None:
                arr = arr[:, isgood_columns]  # remove all NaN columns
                arr = arr[isgood_rows, :]  # remove all NaN rows
                if arr_name != "_aper_xs" and arr_name != "_aper_ys":
                    arr *= finite_mask
                setattr(self, arr_name, arr)
        #
        # Update _aper_extent to reflect new arrays with NaN rows & columns removed
        #
        if (
            self._aper_extent is not None
            and self._aper_xs is not None
            and self._aper_ys is not None
        ):
            first_column = self._aper_xs[:, 0][0]
            last_column = self._aper_xs[:, -1][0]
            first_row = self._aper_ys[0, :][0]
            last_row = self._aper_ys[-1, :][0]
            self._aper_extent = [
                first_column + center[0],
                last_column + center[0],
                first_row + center[1],
                last_row + center[1],
            ]  # we use +/- half_px_scale to center matplotlib tickmarks on pixels

    @staticmethod
    def _rotate_ab_to_xy(a, b, rotation, px_scale_arcsec):
        """
        Convert a (possibly rotated) source's semimajor- and semiminor-axis dimensions to
        the x- and y- pixel dimensions needed for
        `castor_etc.Photometry._create_aper_arrs()`.

        Below is modified rotation matrix to ensure the full aperture is covered in
        subsequent arrays.

        Parameters
        ----------
          a, b :: `astropy.units.Quantity` objects in units of arcsec
            The angular dimensions of the source's semimajor and semiminor axis.

          rotation :: float
            The CCW angle of the source with respect to the x-axis in radians.

          px_scale_arcsec :: float
            The pixel scale in arcsec/pixel.

        Returns
        -------
          x, y :: floats
            The x- and y- pixel dimensions needed for
            `castor_etc.Photometry._create_aper_arrs()`.

        """
        if abs(a.value - b.value) >= 1e-15:
            # Non-circular aperture
            abs_sin_rotate = abs(np.sin(rotation))
            abs_cos_rotate = abs(np.cos(rotation))
            x = (abs_cos_rotate * a + abs_sin_rotate * b).value
            y = (abs_sin_rotate * a + abs_cos_rotate * b).value
            # # Round to nearest multiple of _px_scale_arcsec
            # x = (
            #     np.ceil((abs_cos_rotate * a + abs_sin_rotate * b) / px_scale_arcsec)
            #     * px_scale_arcsec
            # ).value
            # y = (
            #     np.ceil((abs_sin_rotate * a + abs_cos_rotate * b) / px_scale_arcsec)
            #     * px_scale_arcsec
            # ).value
            x = x if x >= a.value else a.value
            y = y if y >= b.value else b.value
        else:
            # Circular aperture
            x = a.value
            y = b.value
            # # Round to nearest multiple of _px_scale_arcsec
            # x = np.ceil(a.value / px_scale_arcsec) * px_scale_arcsec
            # y = np.ceil(b.value / px_scale_arcsec) * px_scale_arcsec
        return x, y

    def use_optimal_aperture(
        self, factor=_OPTIMAL_APER_FACTOR, quiet=False, overwrite=False
    ):
        """
        Uses the "optimal" circular aperture calculated from the telescope PSF's
        full-width at half-maximum (FWHM). Note that this aperture is only valid for point
        sources.

        The default optimal aperture factor assumes the point spread function (PSF) is a
        2D Gaussian (more specifically, a 2D multivariate Normal distribution) with the
        same FWHM as the telescope's FWHM. The "true" optimal aperture differs between
        passbands since they all have different PSFs.

        We estimate the fraction of flux contained within the aperture (i.e., the
        encircled energy) by supersampling the user's aperture at the PSF's (supersampled)
        resolution. We then compare the flux within the aperture to the total flux given
        by the sum of the PSF values. Note that the noiseless image (before PSF
        convolution) as well as the PSF array itself should both sum to roughly 1. The sum
        of the `source_weights` attribute for a particular passband gives the encircled
        energy within that passband.

        Parameters
        ----------
          factor :: int or float
            The factor by which to scale the telescope's FWHM. The radius of the optimal
            aperture will be `R = factor * (FWHM/2)`.

          quiet :: bool
            If False, print a warning if the point source's diameter is larger than the
            telescope's FWHM.

          overwrite :: bool
            If True, allow overwriting of any existing aperture associated with this
            `Photometry` object.

        Attributes
        ----------
          _aper_mask :: (M x N) 2D array of floats
            The aperture mask, including fractional pixel weights. The aperture mask shows
            the fractional overlap between the aperture and the pixel, ranging from (0,
            1]. A pixel that is wholly contained in the aperture has an aperture weight
            equal to 1. Similarly, a pixel that is partially contained in the aperture has
            a weight between 0 and 1 (exclusive) directly proportional to the area of the
            aperture that overlaps the pixel. A pixel that is wholly outside the aperture
            has a weight of NaN.

          source_weights :: dict of (M x N) 2D array of floats
            The pixel weights for the source for each of the telescope's passbands,
            including fractional pixel weights. These weights represent the fraction of
            flux from the source contained in each pixel. A pixel with a source weight of
            0.1 means that 10% of the flux from the source is contained within that pixel.

          sky_background_weights :: (M x N) 2D array of floats
            The pixel weights for the sky background (incl. Earthshine, zodiacal light,
            geocoronal emission).
            The sky background noise is given per sq. arcsec, so a 1 sq. arcsec pixel with
            a sky background weight of 1 means the pixel receives the full sky background
            noise (for that pixel), while a sky background weight of 0.4 means the pixel
            only receives 40% of the full sky background noise for that pixel.

          dark_current_weights :: (M x N) 2D array of floats
            The pixel weights for the dark current.
            A dark current weight of 1 means the pixel experiences the full dark current
            rate (per pixel), while a dark current weight of 0.4 means the dark current
            noise for that pixel is 40% of the telescope's dark current value.

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
        aper_radius_arcsec = factor * 0.5 * self.TelescopeObj.fwhm.to(u.arcsec).value
        self._aper_area = (
            np.pi * aper_radius_arcsec * aper_radius_arcsec * (u.arcsec**2)
        )
        self._assign_exact_npix()
        px_scale_arcsec = self.TelescopeObj.px_scale.to(u.arcsec).value
        #
        # Source weights representing 100% of the flux from the source
        # (This is simply the sum of the PSF values since this is a point source)
        #
        sum_tot_source_weights = {
            band: np.nansum(self.TelescopeObj.psfs[band])
            for band in self.TelescopeObj.passbands
        }
        sum_tot_source_weights["noiseless"] = 1.0
        #
        # Create source weights with arbitrary source flux profile through aperture
        #
        # Recall all internal angle aperture angles are in arcsec
        center = [0, 0]  # arcsec
        center_px = self._create_aper_arrs(
            # # Round to nearest multiple of px_scale
            # np.ceil(aper_radius_arcsec / px_scale_arcsec) * px_scale_arcsec,
            # np.ceil(aper_radius_arcsec / px_scale_arcsec) * px_scale_arcsec,
            aper_radius_arcsec,
            aper_radius_arcsec,
            center,
            overwrite=overwrite,
        )
        source_weights = self._calc_source_weights(center)
        #
        # Convolve source weights with PSF and normalize so source weights represent
        # fraction of flux contained within pixel
        #
        self.source_weights["noiseless"] = (
            source_weights / sum_tot_source_weights["noiseless"]
        )
        for band in self.TelescopeObj.passbands:
            self.source_weights[band] = (
                oaconvolve(source_weights, self.TelescopeObj.psfs[band], mode="same")
                / sum_tot_source_weights[band]
            )
        #
        # Create aperture
        #
        aper_radius_px = (
            aper_radius_arcsec
            / px_scale_arcsec
            * self.TelescopeObj.psf_supersample_factor
        )
        aper = EllipticalAperture(
            positions=center_px, a=aper_radius_px, b=aper_radius_px, theta=0
        )
        #
        # Restrict weight maps to aperture (which is an unrotated ellipse)
        #
        aper_mask = aper.to_mask(method="exact").to_image(source_weights.shape)
        self._aper_mask = aper_mask
        self.sky_background_weights *= aper_mask
        self.dark_current_weights *= aper_mask
        for band in self.source_weights:
            self.source_weights[band] *= aper_mask
        self._bin_arrs_remove_nans(center)
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
        self, a, b, center=[0, 0] << u.arcsec, rotation=0, quiet=False, overwrite=False
    ):
        """
        Use an elliptical aperture.

        The fraction of the source's flux enclosed by the aperture is estimated by
        supersampling the specified aperture at the PSF's supersampled resolution and
        comparing the enclosed flux to the flux from a sufficiently large aperture that is
        centered on the source. The sum of the `source_weights` attribute for a particular
        passband gives the fraction of flux enclosed within the aperture for that
        passband.

        Parameters
        ----------
          a, b :: int or float or `astropy.Quantity` angle
            The angular length of the semimajor and semiminor axes of the aperture,
            respectively. If int or float, a/b is assumed to be in pixel units.

          center :: 2-element `astropy.Quantity` angles array
            The (x, y) center of the aperture relative to the center of the source.
            Positive values means the source is displaced to the bottom/left relative to
            the aperture center. This is because the source will always be at (0, 0) and
            the aperture center will be at this `center` value.

          rotation :: int or float
            The counter-clockwise rotation angle in degrees of the ellipse's semimajor
            axis from the positive x-axis. If rotation is 0, the semimajor axis is along
            the x-axis and the semiminor axis is along the y-axis.

          quiet :: bool
            If True, do not print a warning if the source is an `ExtendedSource` object.

          overwrite :: bool
            If True, allow overwriting of any existing aperture associated with this
            `Photometry` object.

        Attributes
        ----------
          _aper_mask :: (M x N) 2D array of floats
            The aperture mask, including fractional pixel weights. The aperture mask shows
            the fractional overlap between the aperture and the pixel, ranging from (0,
            1]. A pixel that is wholly contained in the aperture has an aperture weight
            equal to 1. Similarly, a pixel that is partially contained in the aperture has
            a weight between 0 and 1 (exclusive) directly proportional to the area of the
            aperture that overlaps the pixel. A pixel that is wholly outside the aperture
            has a weight of NaN.

          source_weights :: dict of (M x N) 2D array of floats
            The pixel weights for the source for each of the telescope's passbands,
            including fractional pixel weights. These weights represent the fraction of
            flux from the source contained in each pixel. A pixel with a source weight of
            0.1 means that 10% of the flux from the source is contained within that pixel.

          sky_background_weights :: (M x N) 2D array of floats
            The pixel weights for the sky background (incl. Earthshine, zodiacal light,
            geocoronal emission).
            The sky background noise is given per sq. arcsec, so a 1 sq. arcsec pixel with
            a sky background weight of 1 means the pixel receives the full sky background
            noise (for that pixel), while a sky background weight of 0.4 means the pixel
            only receives 40% of the full sky background noise for that pixel.

          dark_current_weights :: (M x N) 2D array of floats
            The pixel weights for the dark current.
            A dark current weight of 1 means the pixel experiences the full dark current
            rate (per pixel), while a dark current weight of 0.4 means the dark current
            noise for that pixel is 40% of the telescope's dark current value.

        Returns
        -------
          None
        """
        #
        # Check inputs
        #
        px_scale_arcsec = self.TelescopeObj.px_scale.to(u.arcsec).value
        if isinstance(a, u.Quantity):
            try:
                a = a.to(u.arcsec)
            except Exception:
                raise TypeError(
                    "a and b must be `astropy.Quantity` angles (e.g., u.arcsec, u.deg) "
                    + " or, if in pixel units, ints or floats"
                )
        else:
            a = a * px_scale_arcsec * u.arcsec
        if isinstance(b, u.Quantity):
            try:
                b = b.to(u.arcsec)
            except Exception:
                raise TypeError(
                    "a and b must be `astropy.Quantity` angles (e.g., u.arcsec, u.deg) "
                    + " or, if in pixel units, ints or floats"
                )
        else:
            b = b * px_scale_arcsec * u.arcsec
        if not isinstance(rotation, Number):
            raise TypeError("rotation must be an int or float")
        #
        # Calculate exact aperture area and number of pixels in aperture
        #
        self._aper_area = np.pi * a * b  # area of ellipse
        self._assign_exact_npix()
        rotation = np.deg2rad(rotation)
        #
        # Get sum of source weights if aperture was centred on source. Must do this before
        # creating actual aperture
        #
        need_overwrite = False  # True only if calculating sum_tot_source_weights
        if isinstance(self.SourceObj, PointSource):
            sum_tot_source_weights = {
                band: np.nansum(self.TelescopeObj.psfs[band])
                for band in self.TelescopeObj.passbands
            }
            sum_tot_source_weights["noiseless"] = 1.0
        elif isinstance(self.SourceObj, (GalaxySource, ExtendedSource)):
            if not quiet and isinstance(self.SourceObj, ExtendedSource):
                warnings.warn(
                    "The ExtendedSource calculation assumes 100% of the flux is "
                    + "contained within the size defined in the ExtendedSource (i.e., "
                    + "by angle_a and angle_b).",
                    UserWarning,
                )
            _tot_x, _tot_y = Photometry._rotate_ab_to_xy(
                self.SourceObj.angle_a.to(u.arcsec),
                self.SourceObj.angle_b.to(u.arcsec),
                self.SourceObj.rotation,
                px_scale_arcsec,
            )
            tot_center_px = self._create_aper_arrs(
                _tot_x, _tot_y, [0, 0], overwrite=overwrite
            )
            tot_source_weights = self._calc_source_weights([0, 0])
            tot_aper = EllipticalAperture(
                positions=tot_center_px,
                a=(self.SourceObj.angle_a.to(u.arcsec) / px_scale_arcsec).value
                * self.TelescopeObj.psf_supersample_factor,
                b=(self.SourceObj.angle_b.to(u.arcsec) / px_scale_arcsec).value
                * self.TelescopeObj.psf_supersample_factor,
                theta=self.SourceObj.rotation,
            )
            # Restrict weight maps to aperture
            tot_aper_mask = tot_aper.to_mask(method="exact").to_image(
                tot_source_weights.shape
            )
            tot_aper_mask[tot_aper_mask < 1e-14] = np.nan
            tot_source_weights *= tot_aper_mask
            # N.B. we do not need to convolve the tot_source_weights with the PSFs because
            # the PSFs just spread out the source flux. We are only interested in the
            # total source flux, and since this quantity is a property of the source
            # (i.e., independent of aperture and passbands), we don't need to convolve the
            # source with the PSFs. We do, however, include the sum of the PSF values
            # because they might not sum to 1, in which case a convolution with this PSF
            # will actually increase/decrease the total flux by this factor.
            _nansum_tot_source_weights = np.nansum(tot_source_weights)  # avoid redoing
            sum_tot_source_weights = {
                band: (
                    np.nansum(self.TelescopeObj.psfs[band]) * _nansum_tot_source_weights
                )
                for band in self.TelescopeObj.passbands
            }
            sum_tot_source_weights["noiseless"] = _nansum_tot_source_weights
            if isinstance(self.SourceObj, GalaxySource):
                # The total flux is twice the flux contained in the half-light radius
                for band in sum_tot_source_weights:
                    sum_tot_source_weights[band] *= 2
            need_overwrite = True
        else:  # CustomSource
            sum_tot_source_weights = {band: 1.0 for band in self.source_weights}
        #
        # Create source weights with arbitrary source flux profile through aperture
        #
        # Recall all internal aperture angles are in arcsec
        center = center.to(u.arcsec).value
        x, y = Photometry._rotate_ab_to_xy(a, b, rotation, px_scale_arcsec)
        center_px = self._create_aper_arrs(x, y, center, overwrite=need_overwrite)
        source_weights = self._calc_source_weights(center)
        #
        # Convolve source weights with PSF and normalize so source weights represent
        # fraction of flux contained within pixel
        #
        self.source_weights["noiseless"] = (
            source_weights / sum_tot_source_weights["noiseless"]
        )
        for band in self.TelescopeObj.passbands:
            self.source_weights[band] = (
                oaconvolve(source_weights, self.TelescopeObj.psfs[band], mode="same")
                / sum_tot_source_weights[band]
            )
        #
        # Create aperture
        #
        aper = EllipticalAperture(
            positions=center_px,
            a=(a / px_scale_arcsec).value * self.TelescopeObj.psf_supersample_factor,
            b=(b / px_scale_arcsec).value * self.TelescopeObj.psf_supersample_factor,
            theta=rotation,
        )
        #
        # Restrict weight maps to aperture (which is an unrotated ellipse)
        #
        aper_mask = aper.to_mask(method="exact").to_image(source_weights.shape)
        self._aper_mask = aper_mask
        self.sky_background_weights *= aper_mask
        self.dark_current_weights *= aper_mask
        for band in self.source_weights:
            self.source_weights[band] *= aper_mask
        self._bin_arrs_remove_nans(center)
        # Check for potential small normalization errors.
        # (Ensure sum of source_weights is <= 1)
        for band in self.source_weights:
            _nansum_source_weights = np.nansum(self.source_weights[band])
            if _nansum_source_weights > 1.0:
                self.source_weights[band] /= _nansum_source_weights
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
        self,
        width,
        length,
        center=[0, 0] << u.arcsec,
        quiet=False,
        overwrite=False,
    ):
        """
        Use a rectangular aperture.

        The fraction of the source's flux enclosed by the aperture is estimated by
        supersampling the specified aperture at the PSF's supersampled resolution and
        comparing the enclosed flux to the flux from a sufficiently large aperture that is
        centered on the source. The sum of the `source_weights` attribute for a particular
        passband gives the fraction of flux enclosed within the aperture for that
        passband.

        Parameters
        ----------
          width, length :: int or float or `astropy.Quantity` angle
            The width (along the x-axis) and length (along the y-axis) of the rectangular
            aperture. If int or float, length/width is assumed to be in pixel units.

          center :: 2-element `astropy.Quantity` angles array
            The (x, y) center of the aperture relative to the center of the source.
            Positive values means the source is displaced to the bottom/left relative to
            the aperture center. This is because the source will always be at (0, 0) and
            the aperture center will be at this `center` value.

          quiet :: bool
            If True, do not print a warning if the source is an `ExtendedSource` object.

          overwrite :: bool
            If True, allow overwriting of any existing aperture associated with this
            `Photometry` object.

        Attributes
        ----------
          _aper_mask :: (M x N) 2D array of floats
            The aperture mask, including fractional pixel weights. The aperture mask shows
            the fractional overlap between the aperture and the pixel, ranging from (0,
            1]. A pixel that is wholly contained in the aperture has an aperture weight
            equal to 1. Similarly, a pixel that is partially contained in the aperture has
            a weight between 0 and 1 (exclusive) directly proportional to the area of the
            aperture that overlaps the pixel. A pixel that is wholly outside the aperture
            has a weight of NaN.

          source_weights :: dict of (M x N) 2D array of floats
            The pixel weights for the source for each of the telescope's passbands,
            including fractional pixel weights. These weights represent the fraction of
            flux from the source contained in each pixel. A pixel with a source weight of
            0.1 means that 10% of the flux from the source is contained within that pixel.

          sky_background_weights :: (M x N) 2D array of floats
            The pixel weights for the sky background (incl. Earthshine, zodiacal light,
            geocoronal emission).
            The sky background noise is given per sq. arcsec, so a 1 sq. arcsec pixel with
            a sky background weight of 1 means the pixel receives the full sky background
            noise (for that pixel), while a sky background weight of 0.4 means the pixel
            only receives 40% of the full sky background noise for that pixel.

          dark_current_weights :: (M x N) 2D array of floats
            The pixel weights for the dark current.
            A dark current weight of 1 means the pixel experiences the full dark current
            rate (per pixel), while a dark current weight of 0.4 means the dark current
            noise for that pixel is 40% of the telescope's dark current value.

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
        px_scale_arcsec = self.TelescopeObj.px_scale.to(u.arcsec).value
        #
        # Get sum of source weights if aperture was centred on source. Must do this before
        # creating actual aperture
        #
        need_overwrite = False  # True only if calculating sum_tot_source_weights
        if isinstance(self.SourceObj, PointSource):
            sum_tot_source_weights = {
                band: np.nansum(self.TelescopeObj.psfs[band])
                for band in self.TelescopeObj.passbands
            }
            sum_tot_source_weights["noiseless"] = 1.0
        elif isinstance(self.SourceObj, (GalaxySource, ExtendedSource)):
            if not quiet and isinstance(self.SourceObj, ExtendedSource):
                warnings.warn(
                    "The ExtendedSource calculation assumes 100% of the flux is "
                    + "contained within the size defined in the ExtendedSource (i.e., "
                    + "by angle_a and angle_b).",
                    UserWarning,
                )
            _tot_x, _tot_y = Photometry._rotate_ab_to_xy(
                self.SourceObj.angle_a.to(u.arcsec),
                self.SourceObj.angle_b.to(u.arcsec),
                self.SourceObj.rotation,
                px_scale_arcsec,
            )
            tot_center_px = self._create_aper_arrs(
                _tot_x, _tot_y, [0, 0], overwrite=overwrite
            )
            tot_source_weights = self._calc_source_weights([0, 0])
            tot_aper = EllipticalAperture(  # sources are elliptical, not rectangular
                positions=tot_center_px,
                a=(self.SourceObj.angle_a.to(u.arcsec) / px_scale_arcsec).value
                * self.TelescopeObj.psf_supersample_factor,
                b=(self.SourceObj.angle_b.to(u.arcsec) / px_scale_arcsec).value
                * self.TelescopeObj.psf_supersample_factor,
                theta=self.SourceObj.rotation,
            )
            # Restrict weight maps to aperture
            tot_aper_mask = tot_aper.to_mask(method="exact").to_image(
                tot_source_weights.shape
            )
            tot_source_weights *= tot_aper_mask
            # N.B. we do not need to convolve the tot_source_weights with the PSFs because
            # the PSFs just spread out the source flux. We are only interested in the
            # total source flux, and since this quantity is a property of the source
            # (i.e., independent of aperture and passbands), we don't need to convolve the
            # source with the PSFs. We do, however, include the sum of the PSF values
            # because they might not sum to 1, in which case a convolution with this PSF
            # will actually increase/decrease the total flux by this factor.
            _nansum_tot_source_weights = np.nansum(tot_source_weights)  # avoid redoing
            sum_tot_source_weights = {
                band: (
                    np.nansum(self.TelescopeObj.psfs[band]) * _nansum_tot_source_weights
                )
                for band in self.TelescopeObj.passbands
            }
            sum_tot_source_weights["noiseless"] = _nansum_tot_source_weights
            if isinstance(self.SourceObj, GalaxySource):
                # The total flux is twice the flux contained in the half-light radius
                for band in sum_tot_source_weights:
                    sum_tot_source_weights[band] *= 2
            need_overwrite = True
        else:  # CustomSource
            sum_tot_source_weights = {band: 1.0 for band in self.source_weights}
        #
        # Create source weights with arbitrary source flux profile through aperture
        #
        # Recall all internal aperture angles are in arcsec
        center = center.to(u.arcsec).value
        center_px = self._create_aper_arrs(
            0.5 * width.value,  # width is along x
            0.5 * length.value,  # length is along y
            center,
            overwrite=need_overwrite,
        )
        source_weights = self._calc_source_weights(center)
        #
        # Convolve source weights with PSF and normalize so source weights represent
        # fraction of flux contained within pixel
        #
        self.source_weights["noiseless"] = (
            source_weights / sum_tot_source_weights["noiseless"]
        )
        for band in self.TelescopeObj.passbands:
            self.source_weights[band] = (
                oaconvolve(source_weights, self.TelescopeObj.psfs[band], mode="same")
                / sum_tot_source_weights[band]
            )
        #
        # Create aperture
        #
        aper = RectangularAperture(
            positions=center_px,
            w=(width / px_scale_arcsec).value * self.TelescopeObj.psf_supersample_factor,
            h=(length / px_scale_arcsec).value * self.TelescopeObj.psf_supersample_factor,
            theta=0,
        )
        #
        # Restrict weight maps to aperture (which is an unrotated ellipse)
        #
        aper_mask = aper.to_mask(method="exact").to_image(source_weights.shape)
        self._aper_mask = aper_mask
        self.sky_background_weights *= aper_mask
        self.dark_current_weights *= aper_mask
        for band in self.source_weights:
            self.source_weights[band] *= aper_mask
        self._bin_arrs_remove_nans(center)
        # Check for potential small normalization errors.
        # (Ensure sum of source_weights is <= 1)
        for band in self.source_weights:
            _nansum_source_weights = np.nansum(self.source_weights[band])
            if _nansum_source_weights > 1.0:
                self.source_weights[band] /= _nansum_source_weights
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
                + "\nNOTE: As of photutils-v1.4.0, there seems to be a bug where, in "
                + "some cases, the rectangular aperture mask will return the wrong "
                + "number of pixels regardless of the array size used to house the "
                + "aperture mask!! This requires a fix from the `photutils` team, "
                + "unfortunately.",
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

    @staticmethod
    def _calc_snr_from_t(
        t,
        signal,
        totskynoise,
        darkcurrent,
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
          - Other_Poisson = (B_Earthshine + B_zodi + B_geocoronal + Darkcurrent)
            is the total Poisson noise due to the sky backgorund (Earthshine, zodiacal
            light, and geocoronal emission), and dark current in electrons/s (per
            pixel)
          - N_read is the number of detector readouts
          - Read is the detector read noise in electrons (per pixel)
        The equation above is based on the formula detailed here:
        <https://hst-docs.stsci.edu/wfc3ihb/chapter-9-wfc3-exposure-time-calculation/9-6-estimating-exposure-times>.

        In practice, the calculation takes in 2D arrays for the source signal, sky
        background, and dark current. The factor of N_pix is already included for these
        arrays since they have fractional pixel weighting. Therefore, simply summing up
        these arrays (excluding NaNs) is enough to calculate the total contribution from
        that element without further multiplying by N_pix. The only exception is N_read
        and the read noise, which are scalars and must be multiplied by an integer number
        of pixels.

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
        variables = [t, signal, totskynoise, readnoise, darkcurrent, nread]
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
            + np.nansum(t * (totskynoise + darkcurrent))
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
          - Other_Poisson = (B_Earthshine + B_zodi + B_geocoronal + Darkcurrent)
            is the total Poisson noise due to the sky backgorund (Earthshine, zodiacal
            light, and geocoronal emission), and dark current in electrons/s (per
            pixel)
          - N_read is the number of detector readouts
          - Read is the detector read noise in electrons (per pixel)
        Note that this is simply the quadratic formula applied to the generic SNR
        equation. Also, the equation above is based on the formula detailed here:
        <https://hst-docs.stsci.edu/wfc3ihb/chapter-9-wfc3-exposure-time-calculation/9-6-estimating-exposure-times>.

        In practice, the calculation takes in 2D arrays for the source signal, sky
        background, and dark current. The factor of N_pix is already included for these
        arrays since they have fractional pixel weighting. Therefore, simply summing up
        these arrays (excluding NaNs) is enough to calculate the total contribution from
        that element without further multiplying by N_pix. The only exception is N_read
        and the read noise, which are scalars and must be multiplied by an integer number
        of pixels.

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
        variables = [snr, signal, totskynoise, readnoise, darkcurrent, nread]
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
        signal_sq = tot_sig * tot_sig
        poisson_noise = tot_sig + np.nansum(totskynoise + darkcurrent)
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
        npix=None,
        nread=1,
        quiet=False,
    ):
        """
        Calculate the signal-to-noise ratio (SNR) reached in a given time or the
        integration time (t) required to reach a given SNR.

        Again, note that the `Source` object associated with the `Photometry` instance
        should have its spectrum in units of flam (erg/s/cm^2/A) unless it is a
        `CustomSource` object.

        Note that, for `CustomSource` objects, the `reddening` parameter will be ignored
        since the `CustomSource` object's surface brightness profile should already be in
        units of electron/s (given for each pixel). Additionally, the SNR/integration time
        calculation will only be performed for the passband corresponding to the
        `CustomSource` object's surface brightness profile.

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

          npix :: int or float or None
            The number of pixels occupied by the source to use for the noise calculations
            (i.e., affects sky background, dark current, and read noise). If None, use the
            automatically calculated "effective pixel" count. Note that the number of
            pixels used for the read noise will always be npix rounded up to the nearest
            integer. Finally, this value does not affect the total signal of the
            source---only the non-signal-derived noise values are affected (e.g., higher
            npix means higher total sky background, dark current, and read noise but not
            higher signal).

          nread :: int
            Number of detector readouts.

          quiet :: bool
            If True, do not print an info messages about the fraction of flux contained
            within the aperture and do not print an info message about ignoring the
            `reddening` parameter when using a `CustomSource`.

        Returns
        -------
          results :: dict
            If `t` is given, this is the SNR reached after `t` seconds. If `snr` is given,
            this is the time in seconds required to reach the given SNR. The results are
            given for each `TelescopeObj` passband or, if using a `CustomSource`, for the
            passband defined in the `CustomSource` object.
        """
        # TODO: can further optimize case when source is a `CustomSource` (i.e., only
        # calculate quantities in available passband)
        #
        # Check some inputs (npix and nread will be checked later, in other functions)
        #
        if (t is None and snr is None) or (t is not None and snr is not None):
            raise ValueError("Exactly one of `t` or `snr` must be specified")
        if not isinstance(reddening, Number):
            raise TypeError("`reddening` must be an int or float")
        if self.source_weights is None:
            raise ValueError("Please choose an aperture first.")
        if isinstance(self.SourceObj, CustomSource) and not quiet:
            print(
                "INFO: The `reddening` parameter is ignored for `CustomSource` objects. "
                + "You can silence this message by setting `quiet=True`."
            )
        #
        # Make some useful variables
        #
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
        px_area_arcsec_sq = self.TelescopeObj.px_area.to(u.arcsec**2).value
        mirror_area_cm_sq = self.TelescopeObj.mirror_area.to(u.cm**2).value
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
        passband_limits_AA = dict.fromkeys(self.TelescopeObj.passband_limits)
        for band in passband_limits_AA:
            passband_limits_AA[band] = (
                self.TelescopeObj.passband_limits[band].to(u.AA).value
            )
        # REVIEW: remove requirement for linewidth?
        # (Unused here... Maybe need it for spectroscopy?)
        for gw, gf, gl in zip(
            self.BackgroundObj.geo_wavelength,
            self.BackgroundObj.geo_flux,
            self.BackgroundObj.geo_linewidth,
        ):
            for band in self.TelescopeObj.passbands:
                # Add geocoronal emission (electron/s) to proper passband(s)
                if (gw >= passband_limits_AA[band][0]) and (
                    gw <= passband_limits_AA[band][-1]
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
                        self.TelescopeObj.full_passband_curves[band]["wavelength"]
                        .to(u.AA)
                        .value,
                        self.TelescopeObj.full_passband_curves[band]["response"],
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
        if isinstance(self.SourceObj, CustomSource):
            source_erate = dict.fromkeys([self.SourceObj.passband], np.nan)
            # The user's surface brightness profile should already give the electron/s
            # induced by the source in the given passband for each pixel. The source
            # weights are simply the user's inputted data, linearly interpolated to the
            # TelescopeObj's pixel scale, and masked with the aperture mask.
            source_erate[self.SourceObj.passband] = self.source_weights
        else:
            source_erate = dict.fromkeys(self.TelescopeObj.passbands, np.nan)
            source_ab_mags = self.SourceObj.get_AB_mag(self.TelescopeObj)
            # (N.B. Unlike background noise, aper_weight_scale and npix cancel out for the
            # source electron/s calculation below. Thus, just use the original _eff_npix
            # and source_weights below; no need for npix or aper_weight_scale)
            #
            # Use fraction of flux contained in aperture + normalized spectrum to
            # calculate the signal in each passband
            for band in source_erate:
                encircled_energy = np.nansum(self.source_weights[band])
                if not quiet:
                    print(
                        f"INFO: Fraction of flux within aperture in {band}-band = "
                        + f"{encircled_energy}"
                    )
                # Account for extinction
                source_passband_mag = source_ab_mags[band] + (
                    self.TelescopeObj.extinction_coeffs[band] * reddening
                )
                # Convert extinction-corrected AB mag to electron/s using Eq. (2) of
                # <https://hst-docs.stsci.edu/acsihb/chapter-9-exposure-time-calculations/9-2-determining-count-rates-from-sensitivities#id-9.2DeterminingCountRatesfromSensitivities-9.2.1>
                # and the passband's photometric zero point
                erate_per_px = (
                    mag_to_flux(
                        source_passband_mag,
                        zpt=self.TelescopeObj.phot_zpts[band],
                    )[0]
                    / self._eff_npix
                )  # electron/s/pixel
                source_erate[band] = (
                    erate_per_px * self._aper_mask * encircled_energy
                )  # array containing the source-produced electron/s for each pixel
        #
        # Calculate desired results (either integration time given SNR or SNR given time)
        #
        if isinstance(self.SourceObj, CustomSource):
            results = dict.fromkeys([self.SourceObj.passband])
        else:
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
                    readnoise=self.TelescopeObj.read_noise,  # constant per pixel
                    read_npix=int(np.ceil(npix)),  # only used for readnoise and nread
                    nread=nread,
                )
        return results
