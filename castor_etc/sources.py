"""
Simulates different astronomical sources and contains different flux profiles.

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

from abc import ABCMeta, abstractmethod
from copy import deepcopy
from numbers import Number

import astropy.units as u
import numpy as np
from astropy.modeling.models import Sersic2D
from photutils.aperture import EllipticalAperture

from . import constants as const
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
        (note that the aperture will still have fractional pixel weighting).

        Parameters
        ----------
          a, b :: scalar or `astropy.Quantity` angles or None
            The semimajor and semiminor axes of the ellipse, respectively. If scalars,
            they are assumed to be in arsec. If a pixel partially overlaps the ellipse,
            the pixel will be weighted by its fractional overlap with the ellipse. If both
            None, generate uniform profile over the entire aperture.

          angle :: scalar
            The counter-clockwise angle to rotate the ellipse, in degrees. At an angle of
            0, a0 and b0 are aligned with the x- and y-axes, respectively.

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
            Parameters
            ----------
              x, y :: 2D arrays of floats
                Elements are angular distance from aperture center.

              center :: 2-element 1D array of floats
                The x- and y-angular offset of the aperture center from the source center.
                Not needed for this function but required to have this as a parameter.
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
            Parameters
            ----------
              x, y :: 2D arrays of floats
                Elements are angular distance from aperture center.
                Not needed for this function but required to have this as a parameter.

              center :: 2-element 1D array of floats
                The x- and y-angular offset of the aperture center from the source center.
                Not needed for this function but required to have this as a parameter.
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
        Exponentially-decaying elliptical profile with arbitrary rotation. The center of
        the ellipse is always at (0, 0), but the aperture center relative to this ellipse
        center (set in the `castor_etc.Photometry` object) can be set to arbitrary
        coordinates.

        Parametric equation of ellipse from: <https://math.stackexchange.com/q/2645689>.
        Also see this visualization: <https://www.desmos.com/calculator/xyiew6ioct>.

        The weights are calculated as:
        ```math
                    weights = [INSERT CORRECT WEIGHTS FORMULA HERE]
                    https://docs.astropy.org/en/stable/api/astropy.modeling.functional_models.Ellipse2D.html#astropy.modeling.functional_models.Ellipse2D
        ```
        TODO: explain docstring better

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

        TODO: make sure docstring is accurate
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
            Parameters
            ----------
              x, y :: 2D arrays of floats
                Elements are angular distance from aperture center.

              center :: 2-element 1D array of floats
                The x- and y-angular offset of the aperture center from the source center.
            """
            source_x = ((x + center[0]) * cos_angle + (y + center[1]) * sin_angle) / a0
            source_y = ((y + center[1]) * cos_angle - (x + center[0]) * sin_angle) / b0
            weights = np.exp(-(source_x * source_x + source_y * source_y))
            return weights

        return _ellipse

    @staticmethod
    def sersic(r_eff, n=4, e=0, angle=0):
        """
        Wrapper function for astropy's Sersic2D model
        (<https://docs.astropy.org/en/stable/api/astropy.modeling.functional_models.Sersic2D.html>).

        This function is for describing how the flux of a source (i.e., a galaxy) changes
        with distance from its center.

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
            Parameters
            ----------
              x, y :: 2D arrays of floats
                Elements are angular distance from aperture center.

              center :: 2-element 1D array of floats
                The x- and y-angular offset of the aperture center from the source center.
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
    Base class for astronomical sources.
    """

    @abstractmethod
    def __init__(self, profile):
        """
        REVIEW: ENSURE THIS DOCSTRING IS CORRECT. CURRENTLY OUTDATED. CENTER AND ANGLE?

        Parameters
        ----------
          profile :: function with header `(x, y, aper_center) -> (2D array of floats, float)`
            `profile` is a function that describes how the flux of the source varies as a
            function of x- & y-axis distance from the center of the source (units of
            arcsec) relative to the flux at the center.

            The function signature for `profile` must have exactly 3 positional arguments,
            `(x, y, aper_center)`, and return an array of floats. That is, the profile
            function should give the relative weights of the flux as it varies with x- and
            y- Cartesian coordinates (origin is at the center of the array). Note that
            the parameters `x` and `y` will be 2D arrays of floats with the same shape and
            have (implied) units of arcsec. `aper_center` is a 2-element 1D array of
            floats specifying the (x, y) center of the aperture relative to the center of
            the source, in (implied) units of arcsec.
        """
        try:
            _in_arr = np.arange(4, dtype=float).reshape((2, 2))
            _result = profile(_in_arr, _in_arr, np.array([0.0, 0.0]))
            if not isinstance(_result, np.ndarray) or (
                np.shape(_result) != np.shape(_in_arr)
            ):
                # TODO: relax the check for shape...
                raise TypeError(
                    "profile should return a 2D np.ndarray like the (future) inputs. "
                    + "See the docstring below for more details.\n"
                    + Source.__init__.__doc__
                )
        except Exception:
            raise ValueError(
                "profile may not be of the correct form. "
                + "profile must be a function with 3 positional arguments (x, y, center) "
                + "and return an array of floats with the same shape as x and y. "
                + "Also check the traceback for more details."
            )
        self.profile = profile
        self.spectrum = None
        self.wavelengths = None
        self.a = None
        self.b = None
        self.dist = None
        self.angle_a = None
        self.angle_b = None

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

    def __init__(self, profile, angle=None, radius=None, dist=None):
        """
        Generate a circular point source given either:
          1. the angle subtended by the source's diamater,
          2. radius & distance.

        If item 2, calculate the angle subtended by a circular source of a given radius at
        a certain distance.

        Parameters
        ----------
          profile :: function with signature `profile(a, b, center)`
            Function with 3 positional parameters that describes how the flux of the
            source (arbitrary units) varies as a function of semimajor (`a`) and semiminor
            (`b`) axis distance from the center of the aperture (arcsec units). The
            `center` parameter is the (x, y) center of the aperture relative to the source
            in arcsec.

          angle :: int or float or `astropy.Quantity`
            The angle subtended by the circular point source's diameter. If a scalar,
            angle is assumed to be in arcsec. If given, radius and dist must both be None.
            If not given, radius and dist must both be provided. This will be assigned to
            the `angle_a` and `angle_b` attributes.

          radius :: int or float or `astropy.Quantity`
            The physical radius of the circular point source. If a scalar, radius is
            assumed to be in units of solar radii. If given, dist must also be provided
            and angle must be None. If not given, angle must be provided. This will be
            assigned to the `a` (semimajor axis) and `b` (semiminor axis) attributes.

          dist :: int or float or `astropy.Quantity`
            The distance to the source. If a scalar, dist is assumed to be in units of
            kpc. If given, radius must also be provided and angle must be None. If not
            given, angle must be provided. This will be assigned to the `dist` attribute.

        Attributes
        ----------
          profile ::

          spectrum ::

          wavelengths ::

          a ::

          b ::

          dist ::

          angle_a ::

          angle_b ::


        Returns
        -------
          `PointSource` instance

        TODO: finish docstring
        """
        super().__init__(profile)
        #
        # Check inputs
        #
        if ((radius is None or dist is None) and angle is None) or (
            angle is not None and (radius is not None or dist is not None)
        ):
            raise ValueError(
                "Either both radius and dist must be provided "
                + "or only npix must be provided"
            )
        #
        # Calculate angular area of source
        #
        if angle is None:
            if not isinstance(radius, u.Quantity):
                radius = radius * const.SUN_RADIUS  # cm
            if not isinstance(dist, u.Quantity):
                dist = dist * u.kpc
            self.a = self.b = radius
            self.dist = dist
            radius = radius.to(u.AU).value
            dist = dist.to(u.pc).value
            # Calculate the angle subtended by source radius from definition of parsec
            angle = np.arctan(radius.value / dist.value) * u.arcsec
        else:
            if not isinstance(angle, u.Quantity):
                angle = angle * u.arcsec
            angle *= 0.5  # convert to angular radius
        self.angle_a = self.angle_b = angle


class ExtendedSource(Source):
    """
    Extended source.
    """

    def __init__(self, profile, angle_a=None, angle_b=None, a=None, b=None, dist=None):
        """
        Generate an extended source given either:
          1. the angle subtended by the source's semimajor (`angle_a`) and semiminor
             (`angle_b`) axes,
          2. semimajor (`a`) and semiminor (`b`) axes and distance (`dist`).

        Note that the orientation of the source can be specified in the `profile` function.

        Parameters
        ----------
          profile :: function with signature `profile(a, b, center)`
            Function with 3 positional parameters that describes how the flux of the
            source (arbitrary units) varies as a function of semimajor (`a`) and semiminor
            (`b`) axis distance from the center of the aperture (arcsec units). The
            `center` parameter is the (x, y) center of the aperture relative to the source
            in arcsec.

          angle_a, angle_b :: int or float or `astropy.Quantity`
            The angle subtended by the extended source's semimajor and semiminor axes,
            respectively. If scalars, they are assumed to be in arcsec. If given, `a`,
            `b`, and `dist` must all be None. If not given, `a`, `b`, and `dist` must all
            be provided.

          a, b :: int or float or `astropy.Quantity`
            The physical semimajor and semiminor axes of the extended source,
            respectively. If scalars, they are assumed to be in units of kpc. If given,
            `dist` must also be provided and `angle_a` & `angle_b` must both be None. If
            not given, `angle_a` and `angle_b` must both be provided.

          dist :: int or float or `astropy.Quantity`
            The distance to the source. If a scalar, `dist` is assumed to be in units of
            kpc. If given, `a` & `b` must also be provided and `angle_a` & `angle_b` must
            both be None. If not given, `angle_a` & `angle_b` must both be provided.

        Attributes
        ----------
          profile ::

          spectrum ::

          wavelengths ::

          a ::

          b ::

          dist ::

          angle_a ::

          angle_b ::

        Returns
        -------
          `ExtendedSource` instance

        TODO: finish docstring
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
        #
        # Calculate angular area of source
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
