#         GNU General Public License v3 (GNU GPLv3)
#
# (c) 2022.                            (c) 2022.
# Government of Canada                 Gouvernement du Canada
# National Research Council            Conseil national de recherches
# Ottawa, Canada, K1A 0R6              Ottawa, Canada, K1A 0R6
# All rights reserved                  Tous droits réservés
#
# NRC disclaims any warranties,        Le CNRC dénie toute garantie
# expressed, implied, or               énoncée, implicite ou légale,
# statutory, of any kind with          de quelque nature que ce
# respect to the software,             soit, concernant le logiciel,
# including without limitation         y compris sans restriction
# any warranty of merchantability      toute garantie de valeur
# or fitness for a particular          marchande ou de pertinence
# purpose. NRC shall not be            pour un usage particulier.
# liable in any event for any          Le CNRC ne pourra en aucun cas
# damages, whether direct or           être tenu responsable de tout
# indirect, special or general,        dommage, direct ou indirect,
# consequential or incidental,         particulier ou général,
# arising from the use of the          accessoire ou fortuit, résultant
# software. Neither the name           de l'utilisation du logiciel. Ni
# of the National Research             le nom du Conseil National de
# Council of Canada nor the            Recherches du Canada ni les noms
# names of its contributors may        de ses  participants ne peuvent
# be used to endorse or promote        être utilisés pour approuver ou
# products derived from this           promouvoir les produits dérivés
# software without specific prior      de ce logiciel sans autorisation
# written permission.                  préalable et particulière
#                                      par écrit.
#
# This file is part of the             Ce fichier fait partie du projet
# FORECASTOR ETC project.              FORECASTOR ETC.
#
# FORECASTOR ETC is free software:     FORECASTOR ETC est un logiciel
# you can redistribute it and/or       libre ; vous pouvez le redistribuer
# modify it under the terms of         ou le modifier suivant les termes de
# the GNU General Public               la "GNU General Public
# License as published by the          License" telle que publiée
# Free Software Foundation,            par la Free Software Foundation :
# either version 3 of the              soit la version 3 de cette
# License, or (at your option)         licence, soit (à votre gré)
# any later version.                   toute version ultérieure.
#
# FORECASTOR ETC is distributed        FORECASTOR ETC est distribué
# in the hope that it will be          dans l'espoir qu'il vous
# useful, but WITHOUT ANY WARRANTY;    sera utile, mais SANS AUCUNE
# without even the implied warranty    GARANTIE : sans même la garantie
# of MERCHANTABILITY or FITNESS FOR    implicite de COMMERCIALISABILITÉ
# A PARTICULAR PURPOSE. See the        ni d'ADÉQUATION À UN OBJECTIF
# GNU General Public License for       PARTICULIER. Consultez la Licence
# more details.                        Générale Publique GNU pour plus
#                                      de détails.
#
# You should have received             Vous devriez avoir reçu une
# a copy of the GNU General            copie de la Licence Générale
# Public License along with            Publique GNU avec FORECASTOR ETC ;
# FORECASTOR ETC. If not, see          si ce n'est pas le cas, consultez :
# <http://www.gnu.org/licenses/>.      <http://www.gnu.org/licenses/>.


import numpy as np
import astropy.units as u

from castor_etc.sources import GalaxySource, PointSource, CustomSource


class Scene:
    """
    Scene class.

    Stores multiple sources and renders them onto
    a detector image.
    """
    def __init__(self, telescope):
        """
        Initialize the Scene.

        Parameters
        ----------
        telescope : Telescope
            Telescope object used to determine the detector dimensions.
        """
        self.shape = (
            telescope.transit_ccd_dim[1],
            telescope.transit_ccd_dim[0],
        )
        # List containing all sources added to the scene
        self.sources_list = []

    def addSource(self, source, mag=None, delta_x=0, delta_y=0, source_name = "source_1"):
        """
        Add a source to the scene.

        Parameters
        ----------
        source : Source
            Source object to add.

        mag : float
            Apparent magnitude of the source.

        delta_x, delta_y : float
            Pixel offsets from the centre of the detector.
        """
         # Store the source and its properties

        if type(source_name) != str:
            raise TypeError("source_name must be a string")

        
        self.sources_list.append(
            (source, mag, int(round(delta_x)), int(round(delta_y)), source_name)
        )

    @staticmethod
    def Gaussian2D(x, y, sigma, a=1.0, x0=0.0, y0=0.0):
        term1 = (x - x0) ** 2 / (2.0 * sigma**2)
        term2 = (y - y0) ** 2 / (2.0 * sigma**2)
        return a * np.exp(-(term1 + term2))

    def displaySource(self, telescope, passband, source_index, exptime=1.0):
        """
        Generate detector image for a single source.

        Parameters
        ----------
        telescope : Telescope
            Telescope object.

        passband : str
            photometric filter.

        source_index : int
            Index of the source in sources_list.

        exptime : float
            Exposure time in seconds.

        Returns
        -------
        array
            Detector image of the chosen source in electrons.
        """

        # Retrieve chosen source
        source, mag, delta_x, delta_y, source_name = self.sources_list[source_index]

        # Exit if CustomSource passband doesn't match
        if isinstance(source, CustomSource) and source.passband != passband:
            return np.zeros(self.shape)

        y_indices, x_indices = np.indices(self.shape)

        # Determine the detector centre
        center_y = self.shape[0] // 2
        center_x = self.shape[1] // 2

        # Source postion on detector
        target_y = center_y + delta_y
        target_x = center_x + delta_x

        # Determine total flux
        if mag is None:
            if isinstance(source, CustomSource):
                total_flux = 1.0 * exptime
            else:
                raise ValueError(
                    f"Magnitude must be specified for {type(source).__name__}."
                )
        else:
            # calc flux from mag w/ zero points
            zpt = telescope.phot_zpts[passband]
            electron_rate = 10 ** (-0.4 * (mag - zpt))
            total_flux = electron_rate * exptime

        # ----- Point Source -----
        if isinstance(source, PointSource):
            # Convert FWHM to pixels and then to Sigma
            fwhm_pixels = (
                telescope.fwhm.to(u.arcsec).value
                / telescope.px_scale.to(u.arcsec).value
            )
            sigma_pixels = fwhm_pixels / (2.0 * np.sqrt(2.0 * np.log(2.0)))

            source_grid = self.Gaussian2D(
                x_indices,
                y_indices,
                sigma=sigma_pixels,
                x0=target_x,
                y0=target_y,
            )

            grid_sum = np.sum(source_grid)

            # Normalize/scale to our actual calculated electron count
            if grid_sum > 0:
                source_grid = (source_grid / grid_sum) * total_flux
            else:
                source_grid = np.zeros(self.shape)

        # ----- Galaxy Source -----
        elif isinstance(source, GalaxySource):
            # Convert detector coordinates into arcseconds
            px_scale_arcsec = telescope.px_scale.to(u.arcsec).value

            x_arcsec = (x_indices - center_x) * px_scale_arcsec
            y_arcsec = (y_indices - center_y) * px_scale_arcsec

            galaxy_center_arcsec = np.array(
                [delta_x * px_scale_arcsec, delta_y * px_scale_arcsec]
            )

            # galaxy profile
            source_grid = source.profile(
                x_arcsec,
                y_arcsec,
                center=-galaxy_center_arcsec,
            )

            # Normalize
            grid_sum = np.sum(source_grid)
            if grid_sum > 0:
                source_grid = (source_grid / grid_sum) * total_flux
            else:
                source_grid = np.zeros(self.shape)

        # ----- Custom Source -----
        elif isinstance(source, CustomSource):
            # convert detector coordinates into arcsecond
            px_scale_arcsec = telescope.px_scale.to(u.arcsec).value

            x_arcsec = (x_indices - center_x) * px_scale_arcsec
            y_arcsec = (y_indices - center_y) * px_scale_arcsec

            source_center_arcsec = np.array(
                [delta_x * px_scale_arcsec, delta_y * px_scale_arcsec]
            )

            # custom source profile
            source_grid = source.profile(
                x_arcsec,
                y_arcsec,
                center=-source_center_arcsec,
            )

            # If a mag was specified, normalize and scale.
            # Otherwise, use the counts/sec in the FITS data multiplied by exptime
            if mag is not None:
                grid_sum = np.sum(source_grid)
                if grid_sum > 0:
                    source_grid = (source_grid / grid_sum) * total_flux
                else:
                    source_grid = np.zeros(self.shape)
            else:
                source_grid = source_grid * exptime

        else:
            raise TypeError(
                "Scene supports only PointSource, GalaxySource, and CustomSource objects."
            )

        

        return source_grid

    def displayScene(self, telescope, passband, exptime=1.0):
        """
        Render all sources in the scene.

        Parameters
        ----------
        telescope : Telescope
            Telescope object.

        passband : str
            Photometric filter.

        exptime : float
            Exposure time in seconds.

        Returns
        -------
        array
            image containing all sources.
        """
        final_image = np.zeros(self.shape)

        # Add each source to the image
        for source_index in range(len(self.sources_list)):
            final_image += self.displaySource(
                telescope,
                passband,
                source_index,
                exptime=exptime,
            )

        return final_image