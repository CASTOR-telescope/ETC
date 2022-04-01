"""
Utilities to characterize the sky background (Earthshine, zodiacal light, geocoronal
emission).

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

from copy import deepcopy
from numbers import Number
from os.path import join

import astropy.units as u
import numpy as np
from astropy.io import fits
from scipy.integrate import simpson

from .conversions import convert_electron_flux_mag
from .data.background.background_values import (
    GEOCORONAL_FLUX_AVG,
    GEOCORONAL_FLUX_HIGH,
    GEOCORONAL_FLUX_LOW,
    GEOCORONAL_LINEWIDTH,
    GEOCORONAL_WAVELENGTH,
)
from .filepaths import DATAPATH
from .telescope import Telescope


class Background:
    """
    Object to characterize the sky background.
    """

    def __init__(
        self,
        earthshine_file=join(DATAPATH, "background", "earthshine.fits"),
        zodi_file=join(DATAPATH, "background", "zodi.fits"),
        mags_per_sq_arcsec=None,
    ):
        """
        Create a `Background` object that characterizes the sky background. Can contain
        Earthshine, zodiacal light, and geocoronal emission lines.

        Parameters
        ----------
          earthshine_file :: str or None
            The absolute path to the file containing the Earthshine data. It must be a
            FITS file with the first field (index zero) containing the wavelengths in
            angstroms and the second field (index one) containing the Earthshine flux in
            flam (erg/cm^2/s/A). If None, omit Earthshine component. If mags_per_sq_arcsec
            provided, the data from this file will not be used.

          zodi_file :: str or None
            The absolute path to the file containing the zodiacal light data. It must be a
            FITS file with the first field (index zero) containing the wavelengths in
            angstroms and the second field (index one) containing the Earthshine flux in
            flam (erg/cm^2/s/A). If None, omit zodiacal light component. If
            mags_per_sq_arcsec provided, the data from this file will not be used.

          mags_per_sq_arcsec :: dict of floats or None
            The sky background AB magnitudes per square arcsecond (incl. Earthshine &
            zodiacal light) in each of the telescope's passbands. If provided, will use
            these values & the telescope's photometric zero-points to convert the sky
            background to electron/s. If None, will use the earthshine file and zodi file
            & passband response curves to convert the sky background spectra to
            electron/s.
            Example: `mags_per_sq_arcsec={"uv": 26.08, "u": 23.74, "g": 22.60}`

        Attributes
        ----------
          earthshine_wavelengths :: array of floats or None
            Earthshine spectrum wavelengths in angstrom.

          earthshine_flam :: array of floats or None
            Earthshine flux density in flam (erg/cm^2/s/A).

          zodi_wavelengths :: array of floats or None
            Zodiacal light spectrum wavelengths in angstrom.

          zodi_flam :: array of floats or None
            Zodiacal light flux density in flam (erg/cm^2/s/A).

          mags_per_sq_arcsec :: dict of floats or None
            The sky background AB magnitudes per square arcsecond (incl. Earthshine &
            zodiacal light) in each of the telescope's passbands.

        Returns
        -------
          background :: `Background` object
            The `Background` object containing earthshine and zodiacal light data.
            Geocoronal emission lines can be added via the `add_geocoronal_emission()`
            function.
        """
        #
        # Check inputs before assining to attributes
        #
        if earthshine_file is not None:
            earthshine = fits.getdata(earthshine_file)
            self.earthshine_wavelengths = earthshine.field(0)
            self.earthshine_flam = earthshine.field(1)
        else:
            self.earthshine_wavelengths = None
            self.earthshine_flam = None
        if zodi_file is not None:
            zodi = fits.getdata(zodi_file)
            self.zodi_wavelengths = zodi.field(0)
            self.zodi_flam = zodi.field(1)
        else:
            self.zodi_wavelengths = None
            self.zodi_flam = None
        if mags_per_sq_arcsec is not None:
            for mag in mags_per_sq_arcsec.values():
                if not isinstance(mag, Number):
                    raise ValueError(
                        "mags_per_sq_arcsec must be a dict of floats. For example, "
                        + "`mags_per_sq_arcsec={'uv': 26.08, 'u': 23.74, 'g': 22.60}`."
                    )
        self.mags_per_sq_arcsec = mags_per_sq_arcsec
        #
        # Initialize attributes for geocoronal emission lines
        #
        self.geo_flux = []
        self.geo_wavelength = []
        self.geo_linewidth = []

    def copy(self):
        """
        Convenience function for creating a deep copy of the `Background` object.

        Parameters
        ----------
          None

        Returns
        -------
          Background_copy :: `Background` object
            The deep copy of the `Background` object.
        """
        return deepcopy(self)

    def add_geocoronal_emission(
        self,
        flux="avg",
        wavelength=GEOCORONAL_WAVELENGTH,
        linewidth=GEOCORONAL_LINEWIDTH,
    ):
        """
        Add a geocoronal emission line.

        Parameters
        ----------
          flux :: float or "high" or "avg" or "low"
            The flux of the geocoronal emission line in erg/cm^2/s/arcsec^2. If "high",
            "avg", or "low", use the pre-defined geocoronal emission line values
            (3.0e-15, 1.5e-15, and 7.5e-17 erg/cm^2/s/arcsec^2, respectively).

          wavelength :: int or float or `astropy.Quantity` length
            The central wavelength of the geocoronal emission line. If an int or float, it
            is assumed to be in angstrom.

          linewidth :: int or float or `astropy.Quantity` length
            The linewidth of the geocoronal emission line. If an int or float, it is
            assumed to be in angstrom.
            This is currently not used in any calculations.

        Attributes
        ----------
          geo_flux :: list of floats
            The fluxes of the geocoronal emission lines in erg/cm^2/s/arcsec^2.

          geo_wavelength :: list of floats
            The central wavelengths of the geocoronal emission lines in angstrom.

          geo_linewidth :: list of floats
            The linewidths of the geocoronal emission lines in angstrom.
            This is currently not used in any calculations.

        Returns
        -------
          None
        """
        if flux == "high":
            flux = GEOCORONAL_FLUX_HIGH
        elif flux == "avg":
            flux = GEOCORONAL_FLUX_AVG
        elif flux == "low":
            flux = GEOCORONAL_FLUX_LOW
        elif not isinstance(flux, Number):
            raise ValueError(
                "`flux` must be a scalar or one of 'high', 'medium', or 'low'"
            )
        if isinstance(wavelength, u.Quantity):
            try:
                wavelength = wavelength.to(u.AA).value
            except Exception:
                raise TypeError("geo_wavelength must be an `astropy.Quantity` length.")
        if isinstance(linewidth, u.Quantity):
            try:
                linewidth = linewidth.to(u.AA).value
            except Exception:
                raise TypeError("geo_linewidth must be an `astropy.Quantity` length.")
        self.geo_flux.append(flux)
        self.geo_wavelength.append(wavelength)
        self.geo_linewidth.append(linewidth)  # currently not used in any calculations

    def calc_mags_per_sq_arcsec(self, TelescopeObj):
        """
        Calculates the sky background AB magnitudes per square arcsecond (not including
        geocoronal emission lines) through the telescope's passbands. This is useful if
        reusing the same sky background object for multiple Photometry/Spectroscopy
        instances with the same `Telescope` object.

        Parameters
        ----------
          TelescopeObj :: `Telescope` object
            The `castor_etc.Telescope` object to use for the sky background calculations
            (needed for passband limits and pivots).

        Attributes
        ----------
          mags_per_sq_arcsec :: dict of floats
            The total sky background AB magnitudes per square arcsecond in the telescope's
            passbands.

        Returns
        -------
          None
        """
        if not isinstance(TelescopeObj, Telescope):
            raise TypeError("TelescopeObj must be a `castor_etc.Telescope` object.")

        def _sum_flam(wavelengths, flam):
            for band in avg_sky_flam:
                in_passband = (wavelengths >= passband_limits_AA[band][0]) & (
                    wavelengths <= passband_limits_AA[band][1]
                )
                passband_range = passband_limits_AA[band][1] - passband_limits_AA[band][0]
                # Mean value of the flux density in the passband
                avg_sky_flam[band] += (
                    simpson(y=flam[in_passband], x=wavelengths[in_passband], even="avg")
                    / passband_range
                )  # erg/cm^2/s/A

        passband_limits_AA = {
            band: limits.to(u.AA).value
            for band, limits in TelescopeObj.passband_limits.items()
        }
        avg_sky_flam = dict.fromkeys(
            TelescopeObj.passband_limits, 0.0
        )  # average flam through band

        if self.earthshine_wavelengths is not None and self.earthshine_flam is not None:
            _sum_flam(self.earthshine_wavelengths, self.earthshine_flam)
        if self.zodi_wavelengths is not None and self.zodi_flam is not None:
            _sum_flam(self.zodi_wavelengths, self.zodi_flam)
        avg_sky_mags = dict.fromkeys(TelescopeObj.passband_limits, np.nan)
        for band in avg_sky_flam:
            # (No need to return uncertainty)
            avg_sky_mags[band] = convert_electron_flux_mag(
                avg_sky_flam[band],
                "flam",
                "mag",
                wavelengths=TelescopeObj.passband_pivots[band],
            )[0]

        # REVIEW:
        # At this point, although the sky background magnitudes are technically just AB
        # magnitudes (converted from erg/s/cm^2/A), these sky backgrounds only make sense
        # if they are in units of AB mags per square arcsecond (recall the sky background
        # values assumed by Dr. Patrick Côté). Thus we will just assume these values are
        # in AB mags per sq. arcsec...
        self.mags_per_sq_arcsec = avg_sky_mags
