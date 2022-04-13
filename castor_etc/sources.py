"""
Simulates different astronomical sources and contains different flux profiles.

Predefined flux profiles:
  - uniform
  - ellipse (exponentially-decaying flux with tunable scale lengths along the major and
    minor axes)
  - sersic

Predefined source types:
  - point (e.g., for stars)
  - extended (e.g., for supernova remnants, galaxies)

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
from abc import ABCMeta, abstractmethod
from copy import deepcopy
from numbers import Number

import astropy.units as u
import numpy as np
from astropy.modeling.models import Sersic2D
from photutils.aperture import EllipticalAperture

from . import constants as const
from .parameters import FWHM
from .spectrum import NormMixin, SpectrumMixin


class Profiles:
    """
    Predefined flux profiles for different astronomical sources.
    """

    def __init__():
        pass

    @staticmethod
    def uniform(a=None, b=None, angle=0):
        """
        Uniform flux profile within some elliptical region or over the entire aperture
        (note that the aperture will still have fractional pixel weighting). This profile
        is typically used for a point source, hence the default values should be
        appropriate (and likely the best choices) in most cases.

        Parameters
        ----------
          a, b :: scalar or `astropy.Quantity` angles or None
            The semimajor and semiminor axes of the ellipse, respectively. If scalars,
            they are assumed to be in arsec. If a pixel partially overlaps the ellipse,
            the pixel will be weighted by its fractional overlap with the ellipse. If both
            None, generate uniform profile over the entire aperture.

          angle :: scalar
            The counter-clockwise angle to rotate the ellipse, in degrees. At an angle of
            0, a and b are aligned with the x- and y-axes, respectively.

        Returns
        -------
          profile :: function
            Function that generates a uniform flux profile.
        """
        if (a is None and b is not None) or (a is not None and b is None):
            raise ValueError("a and b must either both be None or both be not None")
        if isinstance(a, u.Quantity):
            try:
                a = a.to(u.arcsec).value
            except Exception:
                raise TypeError("a must be an angular quantity (e.g., u.arcsec)")
        if isinstance(b, u.Quantity):
            try:
                b = b.to(u.arcsec).value
            except Exception:
                raise TypeError("b must be an angular quantity (e.g., u.arcsec)")
        if (a is not None and a <= 0) or (b is not None and b <= 0):
            raise ValueError("a and b must be > 0 arcseconds")
        if not isinstance(angle, Number):
            raise TypeError("angle must be a scalar")
        angle = np.deg2rad(angle)

        def _uniform_ellipse(x, y, center):
            """
            Uniform flux profile over an ellipse (with fractional pixel weighting).

            Parameters
            ----------
              x, y :: 2D arrays of floats
                Elements are angular distance, in arcsec, from aperture center.

              center :: 2-element 1D array of floats
                The x- and y-angular offset, in arcsec, of the aperture center from the
                source center.
                Not needed for this function but required to have this as a parameter.

            Returns
            -------
              weights :: 2D array of floats
                The flux at each pixel relative to the flux at the source center.
            """
            px_scale_x = (x[0][-1] - x[0][0]) / (len(x) - 1)  # arcsec per pixel
            px_scale_y = (y[-1][-1] - y[0][0]) / (len(y) - 1)  # arcsec per pixel
            center_px = [0.5 * (x.shape[1] - 1), 0.5 * (y.shape[0] - 1)]
            aper = EllipticalAperture(
                positions=center_px, a=a / px_scale_x, b=b / px_scale_y, theta=angle
            )
            weights = aper.to_mask(method="exact").to_image(shape=x.shape)
            weights[weights < 1e-12] = np.nan
            return weights

        def _uniform_aper(x, y, center):
            """
            Uniform flux profile over the entire aperture.

            Parameters
            ----------
              x, y :: 2D arrays of floats
                Elements are angular distance, in arcsec, from aperture center.

              center :: 2-element 1D array of floats
                The x- and y-angular offset, in arcsec, of the aperture center from the
                source center.
                Not needed for this function but required to have this as a parameter.

            Returns
            -------
              weights :: 2D array of floats
                The flux at each pixel relative to the flux at the source center.
            """
            return np.ones(x.shape, dtype=float)

        if a is None and b is None:
            return _uniform_aper
        else:
            return _uniform_ellipse

    @staticmethod
    def ellipse(
        a0,
        b0,
        angle=0,
    ):
        """
        Exponentially-decaying elliptical profile with arbitrary rotation. This profile is
        designed such that the center of the ellipse is always at (0, 0) but the aperture
        center (specified in the `castor_etc.Photometry` object) relative to this ellipse
        center can be set to arbitrary coordinates.

        The weights are calculated as:
        ```math
                weights = exp[-( (source_a/a0)^2 + (source_b/b0)^2 )]
        ```
        where
        ```math
                source_a = (x + x0) * cos(angle) + (y + y0) * sin(angle)
        ```
        and
        ```math
                source_b = (y + y0) * cos(angle) - (x + x0) * sin(angle)
        ```
        - x and y are the x- and y-coordinates of the aperture pixel (in arcsec)
        relative to the aperture (not source) center,
        - x0 and y0 are the x- and y-coordinates of the aperture center relative to the
        source center, in arcsec,
        - angle is the counter-clockwise rotation of the ellipse, in degrees. At an
        angle of 0, the semimajor axis (corresponding to the same line as a0) is
        aligned with the x-axis and the semiminor axis (corresponding to the same line
        as b0) is aligned with the y-axis.


        Parameters
        ----------
          a0, b0 :: scalar or `astropy.Quantity` angles
            The "scale-length" semimajor and semiminor axes of the ellipse, respectively.
            That is, the flux at the elliptical level curve defined by a0 and b0 is a
            factor of e less than the flux at the center of the ellipse. If scalars, a0
            and b0 are assumed to be in arcsec.

          angle :: scalar
            The counter-clockwise angle to rotate the ellipse, in degrees. At an angle of
            0, a0 and b0 are aligned with the x- and y-axes, respectively.

        Returns
        -------
          profile :: function
            Function that generates an exponentially-decaying elliptical profile.
        """
        if isinstance(a0, u.Quantity):
            try:
                a0 = a0.to(u.arcsec).value
            except Exception:
                raise TypeError("a0 must be an angular quantity (e.g., u.arcsec)")
        if isinstance(b0, u.Quantity):
            try:
                b0 = b0.to(u.arcsec).value
            except Exception:
                raise TypeError("b0 must be an angular quantity (e.g., u.arcsec)")
        if not isinstance(angle, Number):
            raise TypeError("angle must be a scalar")
        if a0 <= 0 or b0 <= 0:
            raise ValueError("a0 and b0 must be > 0 arcseconds")
        angle = np.deg2rad(angle)
        sin_angle, cos_angle = np.sin(angle), np.cos(angle)

        def _ellipse(x, y, center):
            """
            Exponentially-decaying elliptical flux profile.

            Parameters
            ----------
              x, y :: 2D arrays of floats
                Elements are angular distance, in arcsec, from aperture center.

              center :: 2-element 1D array of floats
                The x- and y-angular offset, in arcsec, of the aperture center from the
                source center.

            Returns
            -------
              weights :: 2D array of floats
                The flux at each pixel relative to the flux at the source center.
            """
            source_x = ((x + center[0]) * cos_angle + (y + center[1]) * sin_angle) / a0
            source_y = ((y + center[1]) * cos_angle - (x + center[0]) * sin_angle) / b0
            weights = np.exp(-(source_x * source_x + source_y * source_y))
            return weights

        return _ellipse

    @staticmethod
    def sersic(r_eff, n=1, e=0, angle=0):
        """
        2D Sersic profile for describing how the flux of a galaxy changes as a function of
        distance from its center. This is actually just a wrapper function for astropy's
        Sersic2D model
        (<https://docs.astropy.org/en/stable/api/astropy.modeling.functional_models.Sersic2D.html>).

        The following docstring is adapted from the astropy Sersic2D docstring.

        Parameters
        ----------
          r_eff :: scalar or `astropy.Quantity` angle
            Effective (half-light) radius. If a scalar, r_eff is assumed to be in arsec.

          n :: scalar
            Sersic index of the galaxy.

          e :: scalar
            Eccentricity of the galaxy.

          angle :: scalar
            The counter-clockwise rotation angle in degrees. At an angle of 0, the
            semimajor and semiminor axes of the galaxy are aligned with the x- and y-axes,
            respectively.

        Returns
        -------
          profile :: function
            Sersic-index profile function.
        """
        #
        # Check inputs
        #
        if isinstance(r_eff, u.Quantity):
            try:
                r_eff = r_eff.to(u.arcsec).value
            except Exception:
                raise TypeError("r_eff must be an angular quantity (e.g., u.arcsec)")
        for param, param_name in zip([n, e, angle], ["n", "e", "angle"]):
            if not isinstance(param, Number):
                raise TypeError(f"{param_name} must be a scalar")
        angle = np.deg2rad(angle)

        def _sersic(x, y, center):
            """
            2D Sersic profile.

            Parameters
            ----------
              x, y :: 2D arrays of floats
                Elements are angular distance, in arcsec, from aperture center.

              center :: 2-element 1D array of floats
                The x- and y-angular offset, in arcsec, of the aperture center from the
                source center.

            Returns
            -------
              weights :: 2D array of floats
                The flux at each pixel relative to the flux at the source center.
            """
            model = Sersic2D(
                amplitude=1,  # arbitrary. We only care about relative weights
                r_eff=r_eff,
                n=n,
                x_0=-center[0],
                y_0=-center[1],
                ellip=e,
                theta=angle,
            )
            weights = model(x, y)
            central_flux = model(-center[0], -center[1])
            weights /= central_flux
            return weights

        return _sersic


class Source(SpectrumMixin, NormMixin, metaclass=ABCMeta):
    """
    Base (abstract) class for astronomical sources.
    """

    @abstractmethod
    def __init__(self, profile, init_dimensions=True):
        """
        Abstract method for all Source subclasses.

        Parameters
        ----------
          profile :: function with a header of
                     `(x, y, aper_center) -> (2D array of floats)`

            `profile` is a function that describes how the flux of the source varies as a
            function of x- & y-axis distance (in arcsec) from the center of the source
            relative to the flux at the source's center.

            The function signature for `profile` must have exactly 3 positional arguments,
            `(x, y, aper_center)`, and return an array of floats.
              - `x` and `y` are 2D arrays representing the angular distance, in arcsec,
              of each pixel from the aperture's center,
              - `aper_center` is a 2-element 1D array of floats representing the x- and y-
              coordinates of the aperture's center relative to the source's center,
              respectively.

            Given these three arguments, the `profile` function should return a 2D array
            of the same shape as `x` and `y`, where each element is the flux at that pixel
            relative to the flux at the center of the source (not aperture). It may be
            helpful to look at how some of the predefined profiles are implemented (e.g.,
            `Profiles.ellipse`).

          init_dimensions :: bool
            If True, initialize attributes related to the source's dimensions (i.e., a, b,
            dist, angle_a, angle_b, rotation) to None. If False, do not initialize these
            attributes.

        Attributes
        ----------
          profile :: function with a header of
                     `(x, y, aper_center) -> (2D array of floats)`
            The 2D profile function specifying how the flux changes as a function of pixel
            coordinates.

          wavelengths, spectrum :: None (for now)
            In the future, these will be 1D arrays characterizing the spectrum of the
            source.

        Returns
        -------
          Instance of a base `Source` object
        """
        #
        # Check profile
        #
        try:
            _test_arr = np.arange(4, dtype=float).reshape((2, 2))
            _test_result = profile(_test_arr, _test_arr, np.array([0.0, 0.0]))
            if not isinstance(_test_result, np.ndarray) or (
                np.shape(_test_result) != np.shape(_test_arr)
            ):
                raise TypeError(
                    "`profile` should return a 2D array of floats. "
                    + "See the docstring below for more details.\n"
                    + Source.__init__.__doc__
                )
        except Exception:
            raise ValueError(
                "profile may not be of the correct form. "
                + "profile must be a function with 3 positional arguments "
                + "`(x, y, aper_center)` and return a 2D array of floats with the same "
                + "shape as x and y. Also check the traceback for more details."
            )
        self.profile = profile
        #
        # Initialize some potential future attributes
        #
        self.spectrum = None
        self.wavelengths = None
        if init_dimensions:
            self.a = None
            self.b = None
            self.dist = None
            self.angle_a = None
            self.angle_b = None
            self.area = None
            self.rotation = None

    def copy(self):
        """
        Convenience method for creating a deep copy of the `Source` object.

        Parameters
        ----------
          None

        Returns
        -------
          Source_copy :: `Source` object
            The deep copy of the `Source` object.
        """
        return deepcopy(self)


class PointSource(Source):
    """
    Point source.
    """

    def __init__(self, profile=None, angle=None, radius=None, dist=None):
        """
        TODO: deprecate `angle`, `radius`, and `dist` parameters

        Generate a circular point source (e.g., a star). Can optionally give either:
          1. the angle subtended by the diameter of the source, or
          2. the physical radius (i.e., with units of length) of the source & distance to
          the source.

        If item 2, this will calculate the angle subtended by a circular source of a given
        radius at a certain distance.

        Note that the angular diameter of a point source should be no larger than the
        telescope's full-width at half-maximum (FWHM).

        Parameters
        ----------
          profile :: function with a header of `profile(x, y, aper_center)` or None
            Function with 3 positional parameters that describes the flux of the
            source at each aperture pixel relative to the flux at the source center.
              - `x` and `y` are 2D arrays representing the angular distance, in arcsec,
              of each pixel from the aperture's center,
              - `aper_center` is a 2-element 1D array of floats representing the x- and y-
              coordinates of the aperture's center relative to the source's center,
              respectively.
            If None, use a uniform flux profile with default arguments (i.e.,
            `Profiles.uniform()`). This should be the preferred option in most, if not
            all, cases.

          angle :: int or float or `astropy.Quantity` or None
            The angle subtended by the circular point source's diameter. If a scalar,
            `angle` is assumed to be in arcsec. If `angle` is given, `radius` and `dist`
            must both be None. If `angle`, `radius`, and `dist` are all None, use the
            default telescope's FWHM (0.15 arcsec).
            Internally, this will be assigned to the `angle_a` and `angle_b` attributes.

          radius :: int or float or `astropy.Quantity` or None
            The physical radius of the circular point source. If a scalar, `radius` is
            assumed to be in units of solar radii. If `radius` is given, `dist` must also
            be provided and `angle` must be None.
            Internally, this will be assigned to the `a` (semimajor axis) and `b`
            (semiminor axis) attributes as well as used to calculate the `angle_a` and
            `angle_b` attributes.

          dist :: int or float or `astropy.Quantity` or None
            The distance to the source. If a scalar, `dist` is assumed to be in units of
            kpc. If `dist` is given, `radius` must also be provided and `angle` must be
            None.
            Internally, this will be assigned to the `dist` attribute as well as used to
            calculate the `angle_a` and `angle_b` attributes.

        Attributes
        ----------
          profile :: function with a header of
                     `(x, y, aper_center) -> (2D array of floats)`
            The 2D profile function specifying how the flux changes as a function of pixel
            coordinates.

          a, b :: `astropy.Quantity` lengths or None
            The physical semimajor and semiminor axis lengths of the circular point
            source; since this source is circular, these will be equal. If radius was not
            specified, these attributes will be None.

          dist :: `astropy.Quantity` lengths or None
            The distance to the source. If dist was not specified, this attribute will be
            None.

          angle_a, angle_b :: `astropy.Quantity` angles
            The angles subtended by the circular point source's radius. Since this source
            is circular, these will be equal.

          rotation :: float
            The CCW rotation of the source relative to the x-axis, in radians. Since this
            is a circular point source, this rotation does not have any effect and is set
            to zero for simplicity.

        Returns
        -------
          `PointSource` instance
        """
        if profile is None:
            profile = Profiles.uniform()
        super().__init__(profile)
        #
        # Check inputs
        #
        if angle is not None or radius is not None or dist is not None:
            warnings.warn(
                "The `angle`, `radius`, and `dist` parameters are deprecated.",
                DeprecationWarning,
            )
        if angle is None and radius is None and dist is None:
            # angle = deepcopy(FWHM)  # must copy or else will pass by reference
            angle = 1e-12 * u.arcsec  # super small value that is unlikely to be > FWHM
        elif ((radius is None or dist is None) and angle is None) or (
            angle is not None and (radius is not None or dist is not None)
        ):
            raise ValueError(
                "`angle` and `radius`/`dist` cannot both be simultaneously specified."
            )
        #
        # Assign attributes to source
        #
        if angle is None:
            if not isinstance(radius, u.Quantity):
                radius = radius * const.SUN_RADIUS  # cm
            if not isinstance(dist, u.Quantity):
                dist = dist * u.kpc
            self.a = self.b = radius
            self.dist = dist
            if radius <= 0 or dist <= 0:
                raise ValueError("`radius` and `dist` must be greater than zero.")
            radius = radius.to(u.AU).value
            dist = dist.to(u.pc).value
            # Calculate the angle subtended by source radius (from definition of parsec)
            angle = np.arctan(radius.value / dist.value) * u.arcsec
        else:
            if not isinstance(angle, u.Quantity):
                angle = angle * u.arcsec
            if angle <= 0:
                raise ValueError("`angle` must be greater than zero.")
            angle *= 0.5  # convert to angular radius
        self.angle_a = self.angle_b = angle
        self.area = np.pi * self.angle_a * self.angle_b
        self.rotation = 0.0


class ExtendedSource(Source):
    """
    Extended source.
    """

    def __init__(
        self, profile, angle_a=None, angle_b=None, a=None, b=None, dist=None, rotation=0.0
    ):
        """
        Generate an extended source given either:
          1. the angles subtended by the source's semimajor (`angle_a`) and semiminor
          (`angle_b`) axes,
          2. the physical semimajor (`a`) and semiminor (`b`) axis lengths and distance
          (`dist`) to the source.

        Note that users should use the `GalaxySource` class to generate galaxies.

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

          angle_a, angle_b :: int or float or `astropy.Quantity`
            The angle subtended by the extended source's semimajor and semiminor axes,
            respectively. If scalars, they are assumed to be in arcsec. If given, `a`,
            `b`, and `dist` must all be None. If not given, `a`, `b`, and `dist` must all
            be provided.

          a, b :: int or float or `astropy.Quantity`
            The physical semimajor and semiminor axes of the extended source,
            respectively. If scalars, they are assumed to be in units of kpc. If given,
            `dist` must also be provided and `angle_a` & `angle_b` must both be None. If
            not given, `angle_a` and `angle_b` must both be provided. Internally, these
            values along with the distance will be used to calculate the `angle_a` and
            `angle_b` attributes.

          dist :: int or float or `astropy.Quantity`
            The distance to the source. If a scalar, `dist` is assumed to be in units of
            kpc. If given, `a` & `b` must also be provided and `angle_a` & `angle_b` must
            both be None. If not given, `angle_a` & `angle_b` must both be provided.
            Internally, these values along with the semimajor and semiminor axis lengths
            will be used to calculate the `angle_a` and `angle_b` attributes.

          rotation :: int or float
            The counter-clockwise (CCW) rotation of the source's semimajor axis relative
            to the x-axis, in degrees. At a rotation of 0 degrees, the source's semimajor
            axis is along the x-axis and the source's semiminor axis is along the y-axis.

        Attributes
        ----------
          profile :: function with a header of
                     `(x, y, aper_center) -> (2D array of floats)`
            The 2D profile function specifying how the flux changes as a function of pixel
            coordinates.

          a, b :: `astropy.Quantity` lengths or None
            The physical semimajor and semiminor axis lengths of the extended source. If
            parameters `a` and `b` were not specified, these attributes will be None.

          dist :: `astropy.Quantity` lengths or None
            The distance to the source. If the parameter `dist` was not specified, this
            attribute will be None.

          angle_a, angle_b :: `astropy.Quantity` angles
            The angles subtended by the extended source's semimajor and semiminor axes.

          rotation :: float
            The CCW rotation of the source's semimajor axis relative to the x-axis, in
            radians. At a rotation of 0 radians, the source's semimajor axis is along the
            x-axis and the source's semiminor axis is along the y-axis.

        Returns
        -------
          `ExtendedSource` instance
        """
        super().__init__(profile)
        #
        # Check inputs
        #
        if (
            (a is None or b is None or dist is None)
            and (angle_a is None or angle_b is None)
        ) or (
            (angle_a is not None or angle_b is not None)
            and (a is not None or b is not None or dist is not None)
        ):
            raise ValueError(
                "Exactly one of (angle_a, angle_b) or (a, b, dist) must be provided"
            )
        if not isinstance(rotation, Number):
            raise ValueError("The rotation angle must be an int/float.")
        #
        # Assign attributes to source
        #
        if angle_a is None:
            if not isinstance(a, u.Quantity):
                a = a * u.kpc
            self.a = a
            if not isinstance(b, u.Quantity):
                b = b * u.kpc
            self.b = b
            if not isinstance(dist, u.Quantity):
                dist = dist * u.kpc
            self.dist = dist
            a = a.to(u.AU).value
            b = b.to(u.AU).value
            dist = dist.to(u.pc).value
            if a <= 0 or b <= 0 or dist <= 0:
                raise ValueError("`a`, `b`, and `dist` must all be greater than zero.")
            # Calculate the angle subtended by source radius from definition of parsec
            self.angle_a = np.arctan(a.value / dist.value) * u.arcsec
            self.angle_b = np.arctan(a.value / dist.value) * u.arcsec
        else:
            if angle_a <= 0 or angle_b <= 0:
                raise ValueError("`angle_a` and `angle_b` must be greater than zero.")
            if isinstance(angle_a, u.Quantity):
                self.angle_a = angle_a.to(u.arcsec)
            else:
                self.angle_a = angle_a * u.arcsec
            if isinstance(angle_b, u.Quantity):
                self.angle_b = angle_b.to(u.arcsec)
            else:
                self.angle_b = angle_b * u.arcsec
        self.area = np.pi * self.angle_a * self.angle_b
        self.rotation = np.deg2rad(rotation)


class GalaxySource(Source):
    """
    A galaxy.
    """

    def __init__(
        self,
        r_eff,
        n,
        angle_a=None,
        angle_b=None,
        a=None,
        b=None,
        dist=None,
        rotation=0.0,
    ):
        """
        Generate a galaxy with a specific effective radius and Sersic index given either:
          1. the angles subtended by the source's semimajor (`angle_a`) and semiminor
          (`angle_b`) axes,
          2. the physical semimajor (`a`) and semiminor (`b`) axis lengths and distance
          (`dist`) to the source.

        Parameters
        ----------
          r_eff :: scalar or `astropy.Quantity` angle
            Effective (half-light) radius. If a scalar, r_eff is assumed to be in arsec.

          n :: scalar
            Sersic index of the galaxy.

          angle_a, angle_b :: int or float or `astropy.Quantity`
            The angle subtended by the extended source's semimajor and semiminor axes,
            respectively. If scalars, they are assumed to be in arcsec. If given, `a`,
            `b`, and `dist` must all be None. If not given, `a`, `b`, and `dist` must all
            be provided.

          a, b :: int or float or `astropy.Quantity`
            The physical semimajor and semiminor axes of the extended source,
            respectively. If scalars, they are assumed to be in units of kpc. If given,
            `dist` must also be provided and `angle_a` & `angle_b` must both be None. If
            not given, `angle_a` and `angle_b` must both be provided. Internally, these
            values along with the distance will be used to calculate the `angle_a` and
            `angle_b` attributes.

          dist :: int or float or `astropy.Quantity`
            The distance to the source. If a scalar, `dist` is assumed to be in units of
            kpc. If given, `a` & `b` must also be provided and `angle_a` & `angle_b` must
            both be None. If not given, `angle_a` & `angle_b` must both be provided.
            Internally, these values along with the semimajor and semiminor axis lengths
            will be used to calculate the `angle_a` and `angle_b` attributes.

          rotation :: int or float
            The counter-clockwise (CCW) rotation of the source's semimajor axis relative
            to the x-axis, in degrees. At a rotation of 0 degrees, the source's semimajor
            axis is along the x-axis and the source's semiminor axis is along the y-axis.

        Attributes
        ----------
          profile :: function with a header of
                     `(x, y, aper_center) -> (2D array of floats)`
            The 2D profile function specifying how the flux changes as a function of pixel
            coordinates. In this case, it is a Sersic profile.

          a, b :: `astropy.Quantity` lengths or None
            The physical semimajor and semiminor axis lengths of the extended source. If
            parameters `a` and `b` were not specified, these attributes will be None.

          dist :: `astropy.Quantity` lengths or None
            The distance to the source. If the parameter `dist` was not specified, this
            attribute will be None.

          angle_a, angle_b :: `astropy.Quantity` angles
            The angles subtended by the extended source's semimajor and semiminor axes.
            These attributes were also used to calculate the eccentricity of the galaxy's
            Sersic profile.

          rotation :: float
            The CCW rotation of the source's semimajor axis relative to the x-axis, in
            radians. At a rotation of 0 radians, the source's semimajor axis is along the
            x-axis and the source's semiminor axis is along the y-axis.

        Returns
        -------
          `GalaxySource` instance
        """
        #
        # Check inputs
        #
        if (
            (a is None or b is None or dist is None)
            and (angle_a is None or angle_b is None)
        ) or (
            (angle_a is not None or angle_b is not None)
            and (a is not None or b is not None or dist is not None)
        ):
            raise ValueError(
                "Exactly one of (angle_a, angle_b) or (a, b, dist) must be provided"
            )
        if not isinstance(rotation, Number):
            raise ValueError("The rotation angle must be an int/float.")
        #
        # Assign attributes to source
        #
        if angle_a is None:
            if not isinstance(a, u.Quantity):
                a = a * u.kpc
            self.a = a
            if not isinstance(b, u.Quantity):
                b = b * u.kpc
            self.b = b
            if not isinstance(dist, u.Quantity):
                dist = dist * u.kpc
            self.dist = dist
            a = a.to(u.AU).value
            b = b.to(u.AU).value
            dist = dist.to(u.pc).value
            if a <= 0 or b <= 0 or dist <= 0:
                raise ValueError("`a`, `b`, and `dist` must all be greater than zero.")
            # Calculate the angle subtended by source radius from definition of parsec
            self.angle_a = np.arctan(a.value / dist.value) * u.arcsec
            self.angle_b = np.arctan(a.value / dist.value) * u.arcsec
        else:
            if isinstance(angle_a, u.Quantity):
                self.angle_a = angle_a.to(u.arcsec)
            else:
                self.angle_a = angle_a * u.arcsec
            if isinstance(angle_b, u.Quantity):
                self.angle_b = angle_b.to(u.arcsec)
            else:
                self.angle_b = angle_b * u.arcsec
            if angle_a <= 0 or angle_b <= 0:
                raise ValueError("`angle_a` and `angle_b` must be greater than zero.")
        eccentricity = np.sqrt(1 - (self.angle_b.value / self.angle_a.value) ** 2)
        profile = Profiles.sersic(r_eff, n, e=eccentricity, angle=rotation)
        super().__init__(profile, init_dimensions=False)
        self.area = np.pi * self.angle_a * self.angle_b
        self.rotation = np.deg2rad(rotation)
