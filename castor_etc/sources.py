"""
Simulates different astronomical sources.
"""

from abc import ABCMeta, abstractmethod
from copy import deepcopy

import astropy.units as u
import numpy as np

from . import constants as const
from .spectrum import NormMixin, SpectrumMixin


def _rotate_ccw(x, y, theta, origin=(0, 0)):
    """
    Rotate a point/array counter-clockwise by theta degrees about the origin. Theta starts
    at zero on the positive x-axis (right) and increases toward the positive y-axis (up).

    Parameters
    ----------
      x, y :: float or array-like
        The x- and y-coordinates of the point/array to rotate.

      theta :: float or array-like
        The angle of rotation in degrees.

      origin :: 2-tuple of floats or array-like with shape (2, shape(x)) (optional)
        The point about which to rotate. Default is (0, 0). If x and y are arrays, this
        should be a 2D array where the first index (row) is the x-coordinate origin and
        the second index (column) is the y-coordinate of the origin.

    Returns
    -------
      x_rot, y_rot :: float or array-like
        The rotated x- and y-coordinates of the point/array.
    """
    theta = np.deg2rad(theta)
    xnew = x - origin[0]
    ynew = y - origin[1]
    xnew2 = np.cos(theta) * xnew - np.sin(theta) * ynew
    ynew2 = np.sin(theta) * xnew + np.cos(theta) * ynew
    x_rot = xnew2 + origin[0]
    y_rot = ynew2 + origin[1]
    return x_rot, y_rot


def _rotate_cw(x, y, theta, origin=(0, 0)):
    """
    Rotate a point/array clockwise by theta degrees about the origin. Theta starts at zero
    on the positive x-axis (right) and increases toward the negative y-axis (down).

    Parameters
    ----------
      x, y :: float or array-like
        The x- and y-coordinates of the point/array to rotate.

      theta :: float or array-like
        The angle of rotation in degrees.

      origin :: 2-tuple of floats or array-like with shape (2, shape(x)) (optional)
        The point about which to rotate. Default is (0, 0). If x and y are arrays, this
        should be a 2D array where the first index (row) is the x-coordinate origin and
        the second index (column) is the y-coordinate of the origin.

    Returns
    -------
      x_rot, y_rot :: float or array-like
        The rotated x- and y-coordinates of the point/array.
    """
    theta = np.deg2rad(theta)
    xnew = x - origin[0]
    ynew = y - origin[1]
    xnew2 = np.cos(theta) * xnew + np.sin(theta) * ynew
    ynew2 = -np.sin(theta) * xnew + np.cos(theta) * ynew
    x_rot = xnew2 + origin[0]
    y_rot = ynew2 + origin[1]
    return x_rot, y_rot


class Profiles:
    """
    Predefined flux profiles for different astronomical sources.
    """

    def __init__():
        pass

    @staticmethod
    def uniform():
        """
        Uniform flux profile.

        Returns
        -------
          profile :: function
            Function that generates a uniform flux profile.
        """

        def _uniform(x, y, center):
            """
            Generate equal weights.
            """
            if np.shape(x) != np.shape(x):
                raise ValueError("Input arr shape mismatch!")
            return np.ones(np.shape(x), dtype=float)

        return _uniform, 0.0

    @staticmethod
    def ellipse(
        a0,
        b0,
        angle=0,
        # center=[0, 0] << u.arcsec
    ):
        """
        Exponentially-decaying elliptical profile with arbitrary rotation. Please set the
        center of the ellipse in when creating the aperture.

        Parametric equation of ellipse from: <https://math.stackexchange.com/q/2645689>.
        Also see this visualization: <https://www.desmos.com/calculator/xyiew6ioct>.

        The weights are calculated as:
                    ```math
                    weights = exp[-((x/a) + (y/b))]
                    ```
        TODO: explain docstring better

        Parameters
        ----------
          a0, b0 :: `astropy.Quantity` angles
            The "scale-length" semimajor and semiminor axes of the ellipse, respectively.
            That is, the flux at the elliptical level curve defined by a0 and b0 is a
            factor of e less than the flux at the center of the ellipse.

          angle :: scalar
            The counter-clockwise angle to rotate the ellipse, in degrees. At an angle of
            0, a0 and b0 are aligned with the x- and y-axes, respectively.

          center :: `astropy.Quantity` angle array
            The x- & y-coordinates of the ellipse center.

            ! SET THE CENTER IN PHOTOMETRY APERTURE !

        Returns
        -------
          profile :: function
            Function that generates an exponentially-decaying elliptical profile.

        TODO: make sure docstring is accurate
        """
        try:
            a0 = a0.to(u.arcsec).value
            b0 = b0.to(u.arcsec).value
            # center = center.to(u.arcsec).value
        except Exception:
            raise TypeError("a0 and b0 must be angular quantities (e.g., u.arcsec)")
        angle = np.deg2rad(angle)  # -angle because xy-indexing + imshow is flipped
        sin_angle, cos_angle = np.sin(angle), np.cos(angle)

        def _ellipse(x, y, center):
            """
            Parameters
            ----------
              x, y :: 2D arrays of floats
                Elements are angular distance from source center. Center of the arrays is
              center :: 2-element 1D array of floats
                (NO LONGER NEEDED) The x- and y-angular offset of the aperture (NOT SOURCE) center.
            """
            # aper_center_idx = [x.shape[0] // 2, x.shape[1] // 2]
            # x = x - center[0]
            # y = y - center[1]

            # source_x = (
            #     x * (a0 * cos_angle * cos_param - b0 * sin_angle * sin_param) + center[0]
            # )
            # source_y = (
            #     y * (a0 * sin_angle * cos_param + b0 * cos_angle * sin_param) + center[1]
            # )

            # ang_param = np.arctan2(y, x)
            # sin_param, cos_param = np.sin(ang_param), np.cos(ang_param)
            # # source_x and source_y give the "scale length" ellipse
            # source_x = a0 * cos_angle * cos_param - b0 * sin_angle * sin_param - center[0]
            # source_y = a0 * sin_angle * cos_param + b0 * cos_angle * sin_param - center[1]

            source_x = ((x + center[0]) * cos_angle + (y + center[1]) * sin_angle) / a0
            source_y = ((y + center[1]) * cos_angle - (x + center[0]) * sin_angle) / b0

            # source_x *= x
            # source_y += y
            # x /= a0
            # y /= b0
            # source_x = x * cos_angle * cos_param - y * sin_angle * sin_param + center[0]
            # source_y = x * sin_angle * cos_param + y * cos_angle * sin_param + center[1]
            # weights = np.exp(-((source_x / a0) ** 2 + (source_y / b0) ** 2))  # "X" pattern
            # weights = np.exp(-(abs(source_x / a0) + abs(source_y / b0)))  # "X" pattern

            # source_x, source_y = _rotate_ccw(x, y, np.rad2deg(angle), origin=center)
            # source_x, source_y = source_y, source_x  # for plotting...
            # weights = source_x
            # x = x + center[0]
            # y = y + center[1]
            # weights = np.exp(-abs(a0 / source_x) - abs(b0 / source_y))
            # x, y = _rotate_ccw(x, y, np.rad2deg(angle), origin=center)
            # weights = y

            # weights = np.exp(-abs(x / source_x) - abs(y / source_y))
            # weights = np.exp(-abs(source_x) - abs(source_y))
            weights = np.exp(-(source_x * source_x + source_y * source_y))
            # weights = np.ones_like(weights)
            # weights = np.exp(-x / a0 - y / b0) * (source_x + source_y)
            # weights = (source_x + source_y)

            # cos_angle == 1 => weights == 1
            # weights = (source_x / a0) ** 2 + (source_y / b0) ** 2
            # weights = (cos_angle * cos_param) ** 2 + (cos_angle * sin_param) ** 2

            # weights = np.exp(-(source_x / a0) - (source_y / b0))

            return weights

        return _ellipse, angle

    # TODO: Sersic profile (N.B. effective radius is an angle)


class Source(SpectrumMixin, NormMixin, metaclass=ABCMeta):
    """
    Base class for astronomical sources.
    """

    @abstractmethod
    def __init__(self, profile_tuple):
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
        profile, profile_angle = profile_tuple
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
                "profile is not of the correct form. "
                + "profile must be a function with 3 positional arguments (x, y, center) "
                + "and return an array of floats with the same shape as x and y."
            )
        self.profile = profile
        self.profile_angle = profile_angle
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

    Attributes
    ----------
      profile ::

      spectrum ::

      wavelengths ::

      a, b ::

      dist ::

      angle_a, angle_b ::

    TODO: finish docstring
    """

    def __init__(self, profile, radius=None, dist=None, angle=None):
        """
        Generate a circular point source given either:
          1. radius & distance,
          2. the (full) angle subtended by the source.

        If item 1, calculate the angle subtended by a circular source of a given radius at
        a certain distance.

        Parameters
        ----------
          profile :: function with signature `profile(flux, a, b, **kwargs)`
            Function with 2 positional parameters that describes how the flux of the
            source (arbitrary units) varies as a function of semimajor (a) and semiminor
            (b) axis distance (pixel units). It must return a single value or an array of
            values (i.e., it should return some multiple of flux).

          radius :: int or float or `astropy.Quantity`
            The physical radius of the circular point source. If a scalar, radius is
            assumed to be in units of solar radii. If given, dist must also be provided
            and angle must be None. If not given, angle must be provided. This will be
            assigned to the `a` (semimajor axis) and `b` (semiminor axis) attributes.

          dist :: int or float or `astropy.Quantity`
            The distance to the source. If a scalar, dist is assumed to be in units of
            kpc. If given, radius must also be provided and angle must be None. If not
            given, angle must be provided. This will be assigned to the `dist` attribute.

          angle :: int or float or `astropy.Quantity`
            The angle subtended by the circular point source's radius. If a scalar, angle
            is assumed to be in arcsec. If given, radius and dist must both be None. If
            not given, radius and dist must both be provided. This will be assigned to the
            `angle_a` and `angle_b` attributes.

          TODO: finish docstring

        Returns
        -------
          `PointSource` instance
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
            radius = radius.to(u.AU).value
            if isinstance(dist, u.Quantity):
                dist = dist.to(u.kpc).value
            self.a = self.b = radius
            self.dist = dist
            # Calculate the angle subtended by source radius from definition of parsec
            angle = np.arctan(radius / (dist * 1000)) * u.arcsec
        else:
            if not isinstance(angle, u.Quantity):
                angle = angle * u.arcsec
            angle *= 0.5  # convert to angular radius
        self.angle_a = self.angle_b = angle


class ExtendedSource(Source):
    """
    Extended source.

    Attributes
    ----------
      profile ::

      spectrum ::

      wavelengths ::

      a, b ::

      dist ::

      angle_a, angle_b ::

    TODO: finish docstring
    """

    def __init__(self, profile, a=None, b=None, dist=None, angle_a=None, angle_b=None):
        """
        Generate an extended source given either:
          1. semimajor (a) and semiminor (b) axes and distance,
          2. the angle subtended by the source's semimajor (angle_a) and semiminor
             (angle_b) axes.

        TODO: finish docstring
        """
        super().__init__(profile)
        raise NotImplementedError("Extended source not implemented yet")


# def circ_source(radius=10.0, dist=0.5 << u.kpc):
#     """
#     Calculate the angle subtended by a circular source of a given radius at a certain
#     distance. If the source subtends an angle less than the full width at half maximum
#     (FWHM) of the telescope's point-spread function (PSF), then the PSF's FWHM will be
#     returned instead.

#     Parameters
#     ----------
#       radius :: int or float or `astropy.Quantity`
#         The physical radius of the circular source. If a scalar, radius is assumed to be
#         in units of solar radii.

#       dist :: int or float or `astropy.Quantity`
#         The distance to the source. If a scalar, dist is assumed to be in units of kpc.

#     Returns
#     -------
#       angle :: float
#         The angle subtended by the circular source, in arcsec. If the source subtends an
#         angle less than the FWHM of the telescope's PSF, the PSF's FWHM will be returned
#         instead.
#     """
#     if not isinstance(radius, u.Quantity):
#         radius = radius * const.SUN_RADIUS  # cm
#     radius = radius.to(u.AU).value
#     if isinstance(dist, u.Quantity):
#         dist = dist.to(u.kpc).value
#     angle = np.arctan(radius / (dist * 1000))  # arcsec, from definition of parsec
#     if angle < params.FWHM.to(u.arcsec).value:
#         return params.FWHM.to(u.arcsec).value
#     return angle
