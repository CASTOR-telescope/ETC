"""
sources.py

Simulates different astronomical sources.

Isaac Cheng - 2022
"""

import astropy.units as u
import numpy as np

from . import constants as const
from . import parameters as params


def circ_source(radius=10.0, dist=0.5 << u.kpc):
    """
    Calculate the angle subtended by a circular source of a given radius at a certain
    distance. If the source subtends an angle less than the full width at half maximum
    (FWHM) of the telescope's point-spread function (PSF), then the PSF's FWHM will be
    returned instead.

    Parameters
    ----------
      radius :: int or float or `astropy.Quantity`
        The physical radius of the circular source. If a scalar, radius is assumed to be
        in units of solar radii.

      dist :: int or float or `astropy.Quantity`
        The distance to the source. If a scalar, dist is assumed to be in units of kpc.

    Returns
    -------
      angle :: float
        The angle subtended by the circular source, in arcsec. If the source subtends an
        angle less than the FWHM of the telescope's PSF, the PSF's FWHM will be returned
        instead.
    """
    if not isinstance(radius, u.Quantity):
        radius = radius * const.SUN_RADIUS  # cm
    radius = radius.to(u.AU).value
    if isinstance(dist, u.Quantity):
        dist = dist.to(u.kpc).value
    angle = np.arctan(radius / (dist * 1000))  # arcsec, from definition of parsec
    if angle < params.FWHM.to(u.arcsec).value:
        return params.FWHM.to(u.arcsec).value
    return angle
