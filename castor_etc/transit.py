"""
Transit calculations.

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

import numpy as np

import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

from astropy.coordinates import SkyCoord
import astropy.units as u
from matplotlib.colors import LogNorm

from tqdm import tqdm

from scipy.ndimage import interpolation, gaussian_filter
from skimage.transform import downscale_local_mean
from .conversions import calc_photon_energy, mag_to_flux,flam_to_AB_mag, convert_electron_flux_mag
from scipy.interpolate import interp1d
from astropy.constants import R_sun, pc


from os.path import join

from . import psf
from astroquery.gaia import Gaia


from castor_etc.background import Background
from castor_etc.sources import PointSource
from castor_etc.telescope import Telescope
from castor_etc.spectrum import getStarData

from .filepaths import DATAPATH
from pytransit import QuadraticModel
from photutils.aperture import EllipticalAperture, aperture_photometry

# The optimal aperture for a point source is a circular aperture with a radius equal to the factor below times half the telescope's FWHM
_OPTIMAL_APER_FACTOR = 1.4

class Observation:
    """
    Observation class.
    """

    def __init__(self, TelescopeObj, SourceObj, BackgroundObj):
        """
        Initialize class for transit calculations.

        Note that the `Source` object (i.e., the `SourceObj` parameter) should be a point source.

        Parameters
        ----------
          TelescopeObj :: `castor_etc.Telescope` object
            The `castor_etc.Telescope` object containing the telescope parameters.

          SourceObj :: `castor_etc.Source` object
            The `castor_etc.Source` object contaning the target source parameters.

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
          `Transit` instance
        """

        #
        # Check inputs
        #

        if not isinstance(TelescopeObj, Telescope):
            raise TypeError("TelescopeObj must be a `castor_etc.Telescope` Object")
        
        if not isinstance(SourceObj, PointSource):
            raise TypeError("SourceObj must be a `castor_etc.PointSource` object. Other sources are not currently supported")
        
        if not isinstance(BackgroundObj, Background):
            raise TypeError("BackgroundObj must be a `castor_etc.Background` object")
        

        #
        # Assign attributes
        #
        self.TelescopeObj = TelescopeObj
        self.SourceObj = SourceObj
        self.BackgroundObj = BackgroundObj

        self.gaia = {} 

        self.t_id = '' # Target name. Used for plots

        self.gs_criteria = {'SN_pix':10, 'exptime': 0.5 * u.second, 'ccd_aperture': 4} # Criteria for guide star

        self.exptime = 60 * u.second # exposure time
        self.nstack = 10 # number of exposures to stack (t_tot = exptime * nstack)
        self.tstart = 0.0 * u.d # light curve start time [days]
        self.tend = 6.0/24. * u.d # light curve end time [days]

        self.ccd_dim = [self.TelescopeObj.ccd_dim[0],self.TelescopeObj.ccd_dim[1]] # [pxl,pxl] used for plotting (including converting ra,dec -> x_pxl, y_pxl)

        self.noise_sources = ['illum','shot','read','dark','jitter', 'sky_background'] # List of noise sources to include
        
        self.passband_name = None # 'uv', 'u', 'g'

        self.gain = 1 # self.TelescopeObj.gain = 2 but we are overiding this.
        self.saturation = 45000.0
        self.read_noise = self.TelescopeObj.read_noise
        self.dark_noise = self.TelescopeObj.dark_current
        self.sky_background_noise = None

        self.xjit = 0.0 # standard deviation of jitter sampling [pxl]
        self.yjit  = 0.0 # standard deviation of jitter sampling [pxl]
        self.xpad = 10 # padding to deal with convolution fall-off
        self.ypad = 10 # padding to deal with convolution fall-off
        self.xout = self.ccd_dim[0] # x-axis
        self.yout = self.ccd_dim[1] # y-axis
        self.noversample = 2 # CCD oversampling factor
        
        self.detector_aperture = 1. * u.m
        self.detector_pixscale = 10. * u.um
        self.illumfile = 'illum_pattern.csv' # CCD illumination pattern file, currently a 2D array of one ones

        self.source_weights = None
        self.sky_background_weights = None
        self.dark_current_weights = None
        self._aper_area = None
        self._aper_xs = None
        self._aper_ys = None
        self._aper_mask = None
        self._eff_npix = None
        self._aper_extent = None
        self._encircled_energy = None
        self._optimal_aperture_radius = _OPTIMAL_APER_FACTOR * 0.5 * TelescopeObj.fwhm

        self.ccd_aperture = (_OPTIMAL_APER_FACTOR * 0.5 * TelescopeObj.fwhm.to(u.arcsec).value)/TelescopeObj.px_scale.to(u.arcsec).value # radius [pxl] for aperture photometry

        self.quiet = False
  
    def calc_sky_background_erate(self):
        """
        Noise associated with the sky background. This calculation is copy and pasted from calc_snr_or_t function in the photometry class.

        Attributes
        ----------
            sky_background_noise :: dict of float
                Sky background noise (electron/s/pixel) associated with each passband name

        Returns
        -------
            None
        """
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
        passband_limits_AA = dict.fromkeys(self.TelescopeObj.passband_limits)
        for band in passband_limits_AA:
            passband_limits_AA[band] = (
                self.TelescopeObj.passband_limits[band].to(u.AA).value
            )
        # REVIEW: remove requirement for linewidth? Unused here... Maybe need it for spectroscopy?
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
        self.sky_background_noise = sky_background_erate

    def _photon_count(self,temp=5780., metallicity=0.0, logg=4.44, 
                        Gmag=7.0, Gmag_abs=4.635,
                        radius=1.0
    ):
        """
        Calculate photon count from model for chosen passband name.
        
        Parameters
        ----------
            temp :: float
                Effective temperature of the stellar atmosphere model
            
            metallicity :: float
                Metallicity of the stellar atmosphere model

            logg :: float
                logg value of the stellar atmosphere model

            Gmag :: float
                Gaia G magnitude of the stellar atmoshpere model

            Gmag_abs :: float
                Absolute Gaia G magnitude of the stellar atmosphere model
            
            radius :: float
                Radius of the star

            dpc :: float
                Distance to the star. Included for testing *Not currently used*

        Returns
            cnt_rate :: 

            gaia_flux ::
        """

        wl, fl = getStarData(temp, metallicity, logg, model_grid=self.SourceObj.stellar_model_grid) # Angstrom, flam

        wavelength_AA = wl.to(u.AA).value

        # Estimate distance using Gmag and estimated Gabs
        dpc = 10**( (Gmag - Gmag_abs)/5. + 1 )

        # Scale surface flux to observed flux using esimated dpc and radius
        fl *= ( ( radius * R_sun ) / ( dpc * pc ) )**2

        # Calculate Gaia bandpass fluxes
        tr_fcn = np.loadtxt(join(DATAPATH,"transit_data/instrument_data/transmission_functions/","GaiaDR2_Passbands.dat"))

        gaia_flux = []
        for j in [1,3,5]: # Select columns listing transmission values, not the errors
            ti = np.where(tr_fcn[:,j] < 99)[0]
            _int_tr = np.interp(wavelength_AA, tr_fcn[ti,0], tr_fcn[ti,j], left=0, right=0)
            gaia_flux.append( np.trapz(_int_tr * fl, wavelength_AA) )

        # Rendition of get_AB_mag method in spectrum.py; calculates the AB magnitude of the source through the telescope's passband.
        passband_wavelength = self.TelescopeObj.full_passband_curves[self.passband_name]['wavelength'].to(u.AA).value

        passband_response = self.TelescopeObj.full_passband_curves[self.passband_name]['response']

        passband_interp = interp1d(
            x=passband_wavelength,
            y=passband_response,
            kind='linear',
            bounds_error=False,
            fill_value=np.nan,
        )

        passband_response = passband_interp(wavelength_AA)
        # Do not integrate NaNs
        isgood_passband = np.isfinite(passband_response)
        isgood_spectrum = np.isfinite(fl)
        if np.any(~isgood_passband):
            if np.all(~isgood_passband):
                raise RuntimeError(
                    "Source spectrum could not be estimated"
                    + f"at any {self.passband_name}-band wavelength"
                )
            elif np.any(
                ~isgood_passband[(wavelength_AA >= passband_wavelength.min()) & (wavelength_AA <= passband_wavelength.max())]
            ):
                warnings.warn(
                    f"{self.passband_name}-band response could not be estimated" + "at some source spectrum wavelengths",
                    RuntimeWarning,
                )

            if np.any(~isgood_spectrum):
                if np.all(~isgood_spectrum):
                    raise RuntimeError("Source spectrum values are all non-finite")
                elif np.any(
                    ~isgood_spectrum[
                        (wavelength_AA >= passband_wavelength.min()) & (wavelength_AA <= passband_wavelength.max())
                    ]
                ): # only warn if there are NaNs/infs in the passband range
                    warnings.warn(
                        f"Source spectrum values are not finite at some wavelengths",
                        RuntimeError,
                    )
                    
        ab_mag = flam_to_AB_mag(
            wavelengths=wavelength_AA[isgood_passband & isgood_spectrum],
            flam=fl[isgood_passband & isgood_spectrum],
            response=passband_response[isgood_passband & isgood_spectrum]
        ) 

        #response depends on the bandpass id chosen.
        e_rate = convert_electron_flux_mag(var1=ab_mag,var1_type="mag",var2_type="electron", phot_zpt=self.TelescopeObj.phot_zpts[self.passband_name])[0]

        cnt_rate = e_rate / self.gain

        return cnt_rate, gaia_flux

    def specify_bandpass(self,passband_name = None):
        """
        Set bandpass specific parameters

        Parameters
        ---------- 
            passband_name :: str
                The bandpass id to determine bandpass specific parameters

        Attributes
        ----------
            passband_name :: str
                The bandpass id to determine bandpass specific parameters
        
        Return
        ------
            None

        """

        #
        # Check inputs
        #
        if (passband_name is None):
            raise ValueError("Bandpass must be specified, i.e., ['uv','u','g']")
        if passband_name not in self.TelescopeObj.passbands:
            raise TypeError("Specified bandpass must belong to CASTOR telescope filter, i.e., ['uv','u','g']")
        else:
            self.passband_name = passband_name


    def id_guide_stars(self,gs_criteria=None,plot_SN=False):
        """
        Parameters
        ----------
            gs_criteria :: dict
                Criteria for guide star
                    SN_pix :: int or float
                        Minimum S/N per pixel associated with each guide star
                    exptime :: int or float
                        Assume 2 Hz for guide camera, i.e., 0.5 sec
                    ccd_aperture :: int or float
                        Radius [pxl] to evaluate max S/N associated with each source, ask Dr. Woods
            
            plot_SN :: bool
                If True, then a Max. S/N per pixel vs. Gaia G [mag] graph is plotted.

        Return
        ------
            None
        """

        if gs_criteria != None:
            if isinstance(gs_criteria['exptime'],u.Quantity):
                try:
                    self.gs_criteria['exptime'] = gs_criteria['exptime'].to(u.second)
                except Exception:
                    raise TypeError(" guide start exposure time must be a `astropy.Quantity` time (e.g., sec,day) ")

        if not isinstance(plot_SN, bool):
            raise TypeError("plot_SN must be either True or False")
        

        # Generate scene simulation using guide star exposure criteria
        _nstack, _exptime, _quiet = self.nstack, self.exptime.to(u.second), self.quiet # Store previous nstack and exptime
        self.nstack, self.exptime, self.quiet = 1, self.gs_criteria['exptime'].to(u.second), True
        SN_pix, _SN_ap = self.scene_sim(return_SN_only=True) # Get guide star S/N per pixel

        # Restore nstack, exptime settings
        self.nstack, self.exptime, self.quiet = _nstack, _exptime, _quiet

        if len(self.SourceObj.gaia['ra'] > 1):
            gs_i = np.arange(1,len(self.SourceObj.gaia['ra']))

            _ti = np.where( (self.SourceObj.gaia['x'][gs_i] > self.gs_criteria['ccd_aperture']) \
                            & (self.SourceObj.gaia['x'][gs_i] < (self.ccd_dim[0] - self.gs_criteria['ccd_aperture']) ) \
                            & (self.SourceObj.gaia['y'][gs_i] > self.gs_criteria['ccd_aperture']) \
                            & (self.SourceObj.gaia['y'][gs_i] < (self.ccd_dim[1] - self.gs_criteria['ccd_aperture']) ) )[0]
            gs_i = gs_i[_ti]
            # Get maximum S/N associated with each potential guide star
            _xx, _yy = np.meshgrid( np.arange(self.ccd_dim[0]), np.arange(self.ccd_dim[1]) )
            _xx, _yy = _xx.T, _yy.T
            SN_pix_max = np.zeros(len(gs_i))
            SN_ap = np.zeros(len(gs_i))
            for i in range(len(gs_i)):
                ti = np.sqrt( (_xx - self.SourceObj.gaia['x'][gs_i[i]])**2 \
                            + (_yy - self.SourceObj.gaia['y'][gs_i[i]])**2 ) < self.gs_criteria['ccd_aperture']
                SN_pix_max[i] = np.max( SN_pix[ti] )
                SN_ap[i] = _SN_ap[gs_i[i]]
            ti = np.where( SN_pix_max > self.gs_criteria['SN_pix'] )[0]
            gs_i = gs_i[ ti ]
            SN_pix_max = SN_pix_max[ ti ]
            SN_ap = SN_ap[ ti ]

            if len(gs_i) > 0:
                # Remove guide stars that are too close together (retaining the brightest one)
                ki = []
                for i in range(len(gs_i)):
                    if (i in ki) == False:
                        _sep = np.sqrt( (self.SourceObj.gaia['x'][gs_i] - self.SourceObj.gaia['x'][gs_i[i]])**2 \
                                      + (self.SourceObj.gaia['y'][gs_i] - self.SourceObj.gaia['y'][gs_i[i]])**2 )
                        ti = np.where( _sep < self.gs_criteria['ccd_aperture'] )[0]
                        ki.append( ti[ np.argmin(self.SourceObj.gaia['Gmag'][gs_i[ti]]) ] )
                ki = np.hstack(ki)
                gs_i = gs_i[ki]
                SN_pix_max = SN_pix_max[ki]
                SN_ap = SN_ap[ki]
            else:
                self.gaia['gs_i'] = []
                self.gaia['gs_SN_pix_max'] = []
        else:
            self.gaia['gs_i'] = []
            self.gaia['gs_SN_pix_max'] = []

        if self.quiet == False:
            print('{:.0f} guide star(s) identified.'.format(len(gs_i)))

        self.gaia['gs_i'] = gs_i
        self.gaia['gs_SN_pix_max'] = SN_pix_max

        if plot_SN:
            fig = plt.figure(num=2)
            plt.clf()
            ax = fig.add_subplot(111)

            ax.scatter(self.SourceObj.gaia['Gmag'][self.gaia['gs_i']],self.gaia['gs_SN_pix_max'])

            ax.set_xlabel('Gaia G [mag]',fontsize=16)
            ax.set_ylabel('Max. S/N per pixel',fontsize=16)

            ylim = np.array(ax.get_ylim())
            ax.set_ylim(ylim)

            xlim = np.array(ax.get_xlim())
            ax.set_xlim(np.flip(xlim))

            ax.xaxis.set_minor_locator(AutoMinorLocator())
            ax.yaxis.set_minor_locator(AutoMinorLocator())
            ax.tick_params(axis='both',labelsize=16)
            ax.tick_params(axis='both',which='major',length=6)
            ax.tick_params(axis='both',which='minor',length=3)
            ax.tick_params(axis='x',which='both',direction='inout')
            ax.tick_params(axis='x',which='both',top=True,direction='in')
            ax.tick_params(axis='y',which='both',direction='inout')
            ax.tick_params(axis='y',which='both',right=True,direction='in')

            fig.tight_layout()

    def _remove_aper_mask_nan_row_col(self,center):
        """
        Remove columns and rows of the aperture mask containing all NaNs. Modifies the
        following attributes: `_aper_mask`, `source_weights`, `sky_background_weights`,
        `dark_current_weights`, `_aper_xs`, `_aper_ys`. Also updates `_aper_extent` to
        reflect changed arrays.

        Parameters
        ----------
          center :: 2-element 1D array of floats
            The (x, y) center of the aperture relative to the center of the source in
            arcsec. Positive values means the source is displaced to the bottom/left
            relative to the aperture center.

        Attributes
        ----------
          _aper_mask :: (M x N) 2D array of floats
            The aperture mask, including fractional pixel weights. The aperture mask shows
            the fractional overlap between the aperture and the pixel, ranging from (0,
            1]. A pixel that is wholly contained in the aperture has an aperture weight
            equal to 1. Similarly, a pixel that is partially contained in the aperture has
            a weight between 0 and 1 (exclusive) directly proportional to the area of the
            aperture that overlaps the pixel. A pixel that is wholly outside the aperture
            has a weight of NaN. Now without rows and columns containing only NaNs.

          source_weights :: (M x N) 2D array of floats
            The pixel weights for the source. Now without rows and columns containing only
            NaNs (based on _aper_mask). Note that changing the source weights for a point
            source will not affect the final photometry calculation. Instead, set the
            `encircled_energy` parameter in the `calc_snr_or_t()` method.

          sky_background_weights :: (M x N) 2D array of floats
            The pixel weights for the sky background. Now without rows and columns
            containing only NaNs (based on _aper_mask).

          dark_current_weights :: (M x N) 2D array of floats
            The pixel weights for the dark current. Now without rows and columns
            containing only NaNs (based on _aper_mask).

          _aper_xs :: (M x N) 2D array of floats
            The aperture array containing the x-coordinates (in arcsec) relative to the
            center of the aperture. Now without rows and columns containing only NaNs
            (based on _aper_mask).

          _aper_ys :: (M x N) 2D array of floats
            The aperture array containing the y-coordinates (in arcsec) relative to the
            center of the aperture. Now without rows and columns containing only NaNs
            (based on _aper_mask).

          _aper_extent :: 4-element 1D list of floats
            The [xmin, xmax, ymin, ymax] extent of the aperture in arcsec (for plotting
            the weight arrays). Updated to reflect the extent containing only rows and
            columns with at least 1 non-NaN value.

        Returns
        -------
          None
        """
        if self._aper_mask is not None:
            isgood_columns = ~np.isnan(self._aper_mask).all(axis=0)
            self._aper_mask = self._aper_mask[:, isgood_columns]
            isgood_rows = ~np.isnan(self._aper_mask).all(axis=1)
            self._aper_mask = self._aper_mask[isgood_rows, :]

            for arr,arr_name in zip(
                [
                    self.source_weights,
                    self.sky_background_weights,
                    self.dark_current_weights,
                    self._aper_xs,
                    self._aper_ys
                ],
                [
                    "source_weights",
                    "sky_background_weights",
                    "dark_current_weights",
                    "_aper_xs",
                    "_aper_ys"
                ],
            ):
                if arr is not None:
                    arr = arr[:, isgood_columns]
                    arr = arr[isgood_rows, :]
                    setattr(self,arr_name,arr)

            if (
                self._aper_extent is not None
                and self._aper_xs is not None
                and self._aper_ys is not None
            ):
                half_px_scale = 0.5 * self.TelescopeObj.px_scale.to(u.arcsec).value
                first_column = self._aper_xs[:,0][0]
                last_column = self._aper_xs[:,-1][0]
                first_row = self._aper_ys[0,:][0]
                last_row = self._aper_ys[-1,:][0]
                self._aper_extent = [
                    first_column - half_px_scale + center[0],
                    last_column + half_px_scale + center[0],
                    first_row - half_px_scale + center[1],
                    last_row + half_px_scale + center[1]
                ]

    def _calc_source_weights(self,center):
        """
        Calculate the source weights for the given profile.

        Parameters
        ----------
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
        telescope_standard_dev_sq = (
            self.TelescopeObj.fwhm.to(u.arcsec).value ** 2
        ) / (4*2*np.log(2))

        source_weights = np.exp(
            -((self._aper_xs + center[0])**2 + (self._aper_ys + center[1])**2) / (2 * telescope_standard_dev_sq)
        )
        
        return source_weights

    def _create_aper_arrs(self,half_x,half_y,center,overwrite=False):
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

        Returns
        -------
          center_px :: 2-element 1D list of floats
            The (x, y) center index of the aperture coordinates arrays (i.e., _aper_xs and
            _aper_ys). To be very explicit, this is not (0, 0) but rather the center of
            the 2D arrays.
        """
        px_scale_arcsec = self.TelescopeObj.px_scale.to(u.arcsec).value
        half_px_scale = 0.5*px_scale_arcsec

        xs = np.arange(-half_x,half_x+half_px_scale,px_scale_arcsec)
        ys = np.arange(-half_y,half_y+half_px_scale,px_scale_arcsec)

        self._aper_xs, self._aper_ys = np.meshgrid(xs,ys,sparse=False,indexing='xy')
        self._aper_extent = [
            xs[0] - half_px_scale + center[0],
            xs[-1] + half_px_scale + center[0],
            ys[0] - half_px_scale + center[1],
            ys[0] + half_px_scale + center[1],
        ]

        self.sky_background_weights = np.ones_like(self._aper_xs)
        self.dark_current_weights = np.ones_like(self._aper_xs)

        center_px = [
            0.5*(self._aper_xs.shape[1]-1),
            0.5*(self._aper_ys.shape[0]-1)
        ]

        return center_px

    def _use_optimal_aperture(self,factor=_OPTIMAL_APER_FACTOR,overwrite=False):
        """
        Parameters
        ----------
          factor :: int or float
            The factor by which to scale the telescope's FWHM. The radius of the optimal
            aperture will be `R = factor * (FWHM/2)`.

          quiet :: bool
            If False, print a message if the point source's diameter is larger than the
            telescope's FWHM as well as print the encircled energy fraction.

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

          source_weights :: (M x N) 2D array of floats
            The pixel weights for the source.
            Note that changing the source weights for a point source will not affect the
            final photometry calculation. Instead, set the `encircled_energy` parameter in
            the `calc_snr_or_t()` method.

          sky_background_weights :: (M x N) 2D array of floats
            The pixel weights for the sky background.

          dark_current_weights :: (M x N) 2D array of floats
            The pixel weights for the dark current.

        Returns
        -------
          None 
        """

        aper_radius_arcsec = factor * 0.5 * self.TelescopeObj.fwhm.to(u.arcsec).value
        self._aper_area = (
            np.pi * aper_radius_arcsec * aper_radius_arcsec * (u.arcsec ** 2)
        )

        center = [0,0]

        px_scale_arcsec = self.TelescopeObj.px_scale.to(u.arcsec).value
        center_px = self._create_aper_arrs(
            np.ceil(aper_radius_arcsec / px_scale_arcsec) * px_scale_arcsec,
            np.ceil(aper_radius_arcsec / px_scale_arcsec) * px_scale_arcsec,
            center,
            overwrite=overwrite
        )

        source_weights = self._calc_source_weights(center)

        aper_radius_px = aper_radius_arcsec / px_scale_arcsec
        aper = EllipticalAperture(
            positions=center_px, a = aper_radius_px, b=aper_radius_px, theta=0
        )

        aper_mask = aper.to_mask(method="exact").to_image(source_weights.shape)
        aper_mask[aper_mask <= 1e-12] = np.nan
        self._aper_mask = aper_mask
        self.source_weights = source_weights * aper_mask
        self.sky_background_weights *= aper_mask
        self.dark_current_weights *= aper_mask
        self._eff_npix = np.nansum(aper_mask)

        self._encircled_energy = 1 - 0.5**(factor*factor)

        self._remove_aper_mask_nan_row_col(center)

    def _point_source_sim(self, target_flux_fraction,scene_phot_count):
        """
        Parameters
        ----------
        target_flux_fraction :: int or float
            Used to reduce target star flux when injecting transit signal

        scene_phot_count :: int or float
            Photon count from the source

        Returns
        -------
            pixels_final.T :: array of floats
                point source simulation
            
            pixels_err.T :: array of floats
                error associated with the point source simulation 
        """

        scene_phot_count *= target_flux_fraction 

        encircled_energy = self._encircled_energy

        count_per_px = scene_phot_count / self._eff_npix # count/pix

        pixels_final = np.zeros_like(self._aper_mask)

        for istack in range(self.nstack):

            source_count = count_per_px * encircled_energy * self._aper_mask
            
            if 'illum' in self.noise_sources:
                il = np.genfromtxt(fname=join(DATAPATH, "transit_data/instrument_data",self.illumfile), delimiter=',')
                z = len(source_count) / len(il)
                il_int = interpolation.zoom(il,z)
                source_count *= il_int

            # Add additional noise

            source_count_noise = np.copy(source_count)

            if 'shot' in self.noise_sources:
                source_count_noise += np.sqrt(np.abs(source_count)) * np.random.normal(size=(source_count.shape[0], source_count.shape[1]))

            if ('read' in self.noise_sources) | ('dark' in self.noise_sources) | ('sky_background' in self.noise_sources):
                _std_dev = 0.
                if 'read' in self.noise_sources:
                    _std_dev += self.read_noise * self.exptime.to(u.second).value / self.gain
                if 'dark' in self.noise_sources:
                    _std_dev += self.dark_noise * self.exptime.to(u.second).value / self.gain
                if 'sky_background' in self.noise_sources:
                    _std_dev += self.sky_background_noise[self.passband_name] * self.exptime.to(u.second).value / self.gain

                source_count_noise += _std_dev * np.random.normal(size=(source_count_noise.shape[0], source_count_noise.shape[1]))


            source_count_noise_int = np.floor(source_count_noise)

            ti = source_count_noise_int < 0
            if len(np.where(ti)[0]) > 0:
                source_count_noise_int[ti] = 0

            pixels_final += source_count_noise_int

        pixels_err = np.sqrt(np.abs(pixels_final)/float(self.nstack))

        return pixels_final.T,pixels_err.T

    def scene_sim(self,
                  all_sources=True,
                  return_scene=False,
                  update_gaia=True,
                  quiet=None,
                  return_SN_only=False,
                  ):
        """
        Generate simulated scene using sources listed in self.SourceObj.gaia

        Parameter
        --------
            all_sources :: Boolean
                True -> Include all Gaia sources in scene
                False -> Include only the target
            return_scene :: Booleans
                True -> scene_flux = self.scene_sim()
            quiet :: Boolean
                To print 'Generating scene...'
            return_SN_only :: Boolean
                Return S/N without updating object (used for guide star identification)

        Return 
        ------
            if return_SN_only:
                SN_pix,SN_ap
            if return_scene:
                pixels_final.T,pixels_err.T
        """

        #
        # Check inputs (Need to implement)
        #
        if (self.passband_name is None):
            raise ValueError("Bandpass must be specified, i.e., ['uv','u','g']")
        
        # Need to check other parameters

        telsecope_fwhm_arcsec = self.TelescopeObj.fwhm.to(u.arcsec).value
        sigma = (telsecope_fwhm_arcsec) / (2 * np.sqrt(2 * np.log(2)))

        if all_sources:
            nsrc = len(self.SourceObj.gaia['ra'])
        else:
            nsrc = 1

        if quiet != None:
            self.quiet = quiet


        # Calculate flux using ETC for each Gaia source
        scene_phot_count = np.zeros(nsrc)
        gaia_flux = np.zeros( (nsrc,3) ) # [G, Bp, Rp] flux

        for n in range(nsrc):
            scene_phot_count[n], gaia_flux[n] = self._photon_count(
            Gmag = self.SourceObj.gaia['Gmag'][n],
            Gmag_abs = self.SourceObj.gaia['Gmag_abs'][n],
            temp = self.SourceObj.gaia['Teff'][n],
            radius = self.SourceObj.gaia['radius'][n],
            metallicity = 0.0, #currently assumes solar
            logg = self.SourceObj.gaia['logg'][n]
            )

            scene_phot_count[n] *= self.exptime.to(u.second).value

        starmodel_flux = np.copy(scene_phot_count) # [counts]

        #
        # Generate scene for single observation
        #
        if self.quiet == False:
            print('Generating scene...')
            pbar = tqdm(total=self.nstack)

        pixels_final = np.zeros( (self.xout, self.yout) ) # Final simulated CCD image
        for istack in range(self.nstack):
            # Add pointing jitter
            if 'jitter' in self.noise_sources:
                xjit, yjit = self.xjit, self.yjit
            else:
                xjit, yjit = 0.0, 0.0

            xpad = self.xpad * self.noversample
            ypad = self.ypad * self.noversample

            # Add each unconvolved source one by one
            if (istack == 0) | ( (xjit > 1.e-3) | (yjit > 1.e-3) ):
                for isource in range(nsrc):
                    # Get (x,y) coordinates with jitter
                    xcoo = self.SourceObj.gaia['x'][isource] + xjit + xpad/2
                    ycoo = self.SourceObj.gaia['y'][isource] + yjit + ypad/2
                    pixels1 = psf.gen_unconv_image(self,starmodel_flux[isource],xcoo,ycoo)
                    if (isource == 0):
                        pixels = np.copy(pixels1)
                    else:
                        pixels += pixels1
                pixels_unconvolved = np.copy(pixels)
            else:
                pixels = np.copy(pixels_unconvolved)

            # Convolve image 
            if ((xjit < 1.e-3) & (yjit < 1.e-3)):
                if istack == 0:
                    pixels_conv = gaussian_filter(pixels,sigma=sigma)
            else:
                pixels_conv = gaussian_filter(pixels,sigma=sigma)

            # Remove padding
            pshape = pixels_conv.shape
            xpad = self.xpad * self.noversample
            ypad = self.ypad * self.noversample
            pixels_conv_ras = pixels_conv[ypad:pshape[0] - ypad, xpad:pshape[1] - xpad]

            # Scale to native resolution (remove oversampling)
            pixels_conv_ras_nav = downscale_local_mean(pixels_conv_ras, (self.noversample, self.noversample)) #array of e/pix

            # Calculate S/N per pixel
            # Currently only calculated for nstack=1, which is what's used for guide star evaluation
            self.calc_sky_background_erate() # sky noise

            SN_pix = np.zeros_like(pixels_conv_ras_nav)
            ti = pixels_conv_ras_nav > 0
            SN_pix[ti] = pixels_conv_ras_nav[ti] / np.sqrt( pixels_conv_ras_nav[ti] + np.sqrt(self.nstack) * self.read_noise + (self.dark_noise + self.sky_background_noise[self.passband_name]) * self.nstack * self.exptime.to(u.second).value)

            # Calculate target's S/N within ccd_aperture
            if istack == 0:
                x0 = self.SourceObj.gaia['x'] + (self.xout - self.ccd_dim[0])/2 - 0.5
                y0 = self.SourceObj.gaia['y'] + (self.yout - self.ccd_dim[1])/2 + 0.5
                _coord = []
                for ll in range(len(x0)):
                    _coord.append([x0[ll], y0[ll]])
            
                aper_radius_px = self.ccd_aperture
                
                aper = EllipticalAperture(positions=_coord, a=aper_radius_px, b=aper_radius_px, theta=0)

                _ap_phot = aperture_photometry(np.abs(pixels_conv_ras_nav), 
                                    aper, method='center')
                _fl = _ap_phot['aperture_sum']
                n_pix = np.floor(aper.area) # Number of pixels in the aperture
                SN_ap = np.array( _fl / np.sqrt( _fl + np.sqrt(self.nstack) * self.read_noise + \
                                                (self.dark_noise + self.sky_background_noise[self.passband_name]) * self.nstack * n_pix * self.exptime.to(u.second).value ) )
            
            # Add illumination pattern
            if 'illum' in self.noise_sources:
                #add illumination pattern
                #---------------------------
                #read illumination pattern from illumfile
                il = np.genfromtxt(fname= join(DATAPATH, "transit_data/instrument_data",self.illumfile), delimiter=',')
                # difference in size between scene and illum
                z = len(pixels_conv_ras_nav) / len(il)
                # lin interpolate illum pattern to scene size
                il_int = interpolation.zoom(il,z)
                # multiply scene by this pattern
                pixels_conv_ras_nav *= il_int

            # Add additional noise sources

            pixels_conv_ras_nav_noise = np.copy(pixels_conv_ras_nav)

            if 'shot' in self.noise_sources:
                # shot-noise
                pixels_conv_ras_nav_noise += np.sqrt(np.abs(pixels_conv_ras_nav)) * np.random.normal(size=(pixels_conv_ras_nav.shape[0], pixels_conv_ras_nav.shape[1]))

            if ('read' in self.noise_sources) | ('dark' in self.noise_sources) | ('sky_background' in self.noise_sources):
                _std_dev = 0.
                if 'read' in self.noise_sources:
                    _std_dev += self.read_noise * self.exptime.to(u.second).value / self.gain
                if 'dark' in self.noise_sources:
                    _std_dev += self.dark_noise * self.exptime.to(u.second).value / self.gain
                if 'sky_background' in self.noise_sources:
                    _std_dev += self.sky_background_noise[self.passband_name] * self.exptime.to(u.second).value / self.gain

                pixels_conv_ras_nav_noise += _std_dev * np.random.normal(size=(pixels_conv_ras_nav_noise.shape[0], pixels_conv_ras_nav_noise.shape[1]))

            # Quantize Image
            pixels_conv_ras_nav_noise_int = np.floor(pixels_conv_ras_nav_noise)

            # np.floor() was returning -1 for very small, positive numbers
            # If pixels_conv_ras_nav_int contains -1, set to 0
            ti = pixels_conv_ras_nav_noise_int < 0
            if len( np.where(ti)[0] ) > 0:
                pixels_conv_ras_nav_noise_int[ti] = 0

            # Add for each frame in stack
            pixels_final += pixels_conv_ras_nav_noise_int

            if self.quiet == False:
                pbar.update()

        if return_SN_only:
            return SN_pix, SN_ap

        # Calculate flux error
        pixels_err = np.sqrt( np.abs(pixels_final) / float(self.nstack) )

        if update_gaia:
            self.gaia['scene'] = pixels_final.T
            self.gaia['scene_err'] = pixels_err.T
            self.gaia['SN_pix'] = SN_pix.T
            self.gaia['SN_ap'] = SN_ap
            self.gaia['counts'] = scene_phot_count
            self.gaia['gaia_flux'] = gaia_flux

        if self.quiet == False:
            pbar.close()

        if return_scene:
            return pixels_final.T, pixels_err.T

    def plot_fov(self,plot_grid=True,
                    plot_guide_stars=True,
                    vmin=None,vmax=None,add_scene_sim=True):
        """
        Ask James

        Parameters
        ----------
            pa :: 
            plot_guide_stars :: Boolean
            vmin ::
            vmax ::
            plot_grid :: Boolean
            add_scene_sim :: Boolean

        Return
        ------
            None
        """

        plt.close(1)
        fig = plt.figure(num=1,figsize=(5.5,5),facecolor='w')
        ax = fig.add_subplot(111)

        ax_pos = ax.get_position()
        ax_pos.x0 = 0.175
        ax_pos.x1 = 0.99
        ax_pos.y0 = 0.108
        ax_pos.y1 = 0.99
        ax.set_position(ax_pos)

        Gmag_lim = [np.min(self.SourceObj.gaia['Gmag']), np.max(self.SourceObj.gaia['Gmag'])]
        sym_lim = [0.5, 100]
        if Gmag_lim[1] == Gmag_lim[0]:
            sym_size = np.zeros_like(self.SourceObj.gaia['Gmag']) + 5
        else:
            sym_size = ((Gmag_lim[1] - self.SourceObj.gaia['Gmag']) / (Gmag_lim[1] - Gmag_lim[0] )) \
                            * (sym_lim[1] - sym_lim[0]) + sym_lim[0]

        # Plot target
        ax.scatter(self.SourceObj.gaia['x'][0],self.SourceObj.gaia['y'][0],edgecolor='r',
                        marker='o',s=80,facecolor='None')

        # Plot FoV boundary
        _x = [0, self.ccd_dim[0], self.ccd_dim[0], 0, 0]
        _y = [0, 0, self.ccd_dim[1], self.ccd_dim[1], 0]
        ax.plot(_x,_y,c='r')

        if plot_guide_stars:
            if ( ('gs_i' in self.gaia.keys()) == False ) & hasattr(self,'gs_criteria'):
                self.id_guide_stars()
            if 'gs_i' in self.gaia.keys():
                ax.scatter(self.SourceObj.gaia['x'][ self.gaia['gs_i'] ],
                           self.SourceObj.gaia['y'][ self.gaia['gs_i'] ],
                           marker='s',zorder=1,facecolor='None',edgecolor='r',lw=0.7)

        ax.set_xlabel(r'x [pxl]',fontsize=16)
        ax.set_ylabel(r'y [pxl]',fontsize=16)

        xlim = int(self.ccd_dim[0]/2) + self.xout * 0.7 * np.array([-1.0,1.0])
        ylim = int(self.ccd_dim[1]/2) + self.yout * 0.7 * np.array([-1.0,1.0])

        if add_scene_sim:
            if ('scene' in self.gaia.keys()) == False:
                self.scene_sim()

            extent = np.array([ int(self.ccd_dim[0]/2) - self.xout/2,
                                int(self.ccd_dim[0]/2) + self.xout/2,
                                int(self.ccd_dim[1]/2) - self.yout/2 - 2,
                                int(self.ccd_dim[1]/2) + self.yout/2 - 2])

            cmap = plt.get_cmap('cividis')
            _f = self.gaia['scene']-np.min(self.gaia['scene'])+1
            _f = np.flip(_f,axis=0) # Go through code to figure this out...

            ax.imshow(_f,
                    norm=LogNorm(vmin=vmin,vmax=vmax),
                    extent=extent,
                    interpolation=None,
                    cmap=cmap,zorder=-100)

        ax.set_xlim(xlim)
        ax.set_ylim(ylim)

        ax.xaxis.set_minor_locator(AutoMinorLocator())
        ax.yaxis.set_minor_locator(AutoMinorLocator())
        ax.tick_params(axis='both',labelsize=16)
        ax.tick_params(axis='both',which='major',length=6)
        ax.tick_params(axis='both',which='minor',length=3)
        ax.tick_params(axis='x',which='both',direction='inout')
        ax.tick_params(axis='x',which='both',top=True,direction='in')
        ax.tick_params(axis='y',which='both',direction='inout')
        ax.tick_params(axis='y',which='both',right=True,direction='in')

        plt.ion()
        plt.show()
        plt.pause(1.e-6)
        
    def specify_pl_model(self,RpRs, P, t0, b, aRs):
        """
        Model planet for transit simulation

        Parameters
        ----------
            RpRs :: int or float
                Planet-star radius ratio
            P :: int or float
                Orbital period <Need to read pytransit for units>
            t0 :: int of float
                Zero epoch <Need to read pytransit for units>
            b :: int or float
                Impact parameter, which is the projected distance between the planet disk and the stellar disk at mid-transit
            aRs :: int of float
               Ratio of the semi-major axis to the stellar radius

        Attributes
        ---------
            pl_model :: dict of arrays of floats
                Planet model parameters

        Return
        ------
            None
        """

        #
        # Check inputs 
        #

        if not isinstance(RpRs, list):
            raise TypeError("The radius ratio can either be a scalar, 1D vector, or a 2D array")
        if not isinstance(P, list):
            raise TypeError("Orbital period can be either be scalars or vectors")
        if not isinstance(t0, list):
            raise TypeError("Zero epoch can be either scaalrs or vectors")
        if not isinstance(b, list):
            raise TypeError("b must be an int or float")
        if not isinstance(aRs, list):
            raise TypeError("aRs must be an int or float")

        pl_model  = {
            'RpRs': RpRs,
            'P': P,
            't0': t0,
            'b': b,
            'aRs': aRs
        }
        self.pl_model = pl_model

    def specify_exposure_parameters(self,exptime=60*u.second, nstack=10, tstart=0.0*u.d, tend=6.0/24. * u.d):
        """
        Exposure parameters for planet transit simulation

        Parameters
        ----------
            exptime :: int
                Exposure time
            nstack :: int 
                Number of exposures to stack (t_tot = exptime*nstack)
            tstart :: int of float
                Light curve start time [days]
            tend :: int or float
                Light curve end time [days]

        Attributes
        ----------
            exptime :: int
                Exposure time
            nstack :: int 
                Number of exposures to stack (t_tot = exptime*nstack)
            tstart :: int of float
                Light curve start time [days]
            tend :: int or float
                Light curve end time [days]
                
        Return
        ------
            None
            
        """

        #
        # Check inputs 
        #

        if isinstance(exptime, u.Quantity):
            try:
                self.exptime = exptime.to(u.second)
            except Exception:
                raise TypeError("exptime must be a `astropy.Quantity` time (e.g., sec, day)")
        if isinstance(nstack, int):
            try:
                self.nstack = nstack
            except Exception:
                raise TypeError("nstack must be an int")
        if isinstance(tstart, u.Quantity):
            try:
                self.tstart = tstart.to(u.d)
            except Exception:
                raise TypeError("tstart must be a `astropy.Quantity` time (e.g., sec, day)")
        if isinstance(tend, u.Quantity):
            try:
                self.tend = tend.to(u.d)
            except Exception:
                raise TypeError("tend must be a `astropy.Quantity` time (e.g., sec,day)")
        
    def calc_pl_model(self,model='pytransit_QuadraticModel',t_grid=[],exp_time=-1):
        """
        Parameters
        ----------
            Planet model
            t_grid :: array
            exp_time :: 

        Return
        ------
            _pl_lc_tot ::

        """
        # To do: include additional/alternative transit models
        if len(t_grid) == 0:
            cadence = self.nstack * self.exptime.to(u.second).value
            t_grid = np.arange(self.tstart.to(u.d).value, self.tend.to(u.d).value, cadence / (24. * 3600.))
        nt = len(t_grid)

        const_G = 6.6743e-8 # [cm3 / (g s2)]

        _pl_lc_tot = np.zeros(nt)
        if hasattr(self,'pl_model'):
            _pl_lc = []
            n_pl = len(self.pl_model['RpRs'])

            if model == 'pytransit_QuadraticModel':
                
                tm = QuadraticModel()

                if exp_time < 0:
                    tot_exptime = self.nstack * self.exptime.to(u.second).value / (24. * 3600.)
                else:
                    tot_exptime = exp_time / (24. * 3600.)
                tm.set_data(time=t_grid, exptimes=tot_exptime, nsamples=10)

                for n in range(n_pl):
                    if 'u' in self.pl_model.keys():
                        ldc = self.pl_model['u']
                    else:
                        ldc = [0.6, 0.3]
                    if 'aRs' in self.pl_model.keys():
                        aRs = self.pl_model['aRs'][n]
                    elif 'rho_s' in self.pl_model.keys():
                        aRs = ( self.pl_model['rho_s'] * const_G * (self.pl_model['P'][n]*24.*3600.)**2 / (3. * np.pi) )**(1./3.)
                    _irad = np.arccos( self.pl_model['b'][n] / self.pl_model['aRs'][n] )
                    if 'e' in self.pl_model.keys():
                        _e, _w = self.pl_model['e'][n], self.pl_model['w'][n]
                    else:
                        _e, _w = 0., 0.
                    _pv = [self.pl_model['RpRs'][n], self.pl_model['t0'][n],
                            self.pl_model['P'][n], self.pl_model['aRs'][n],
                            _irad, _e, _w ]
                    _pl_lc.append( tm.evaluate_ps(_pv[0], ldc, _pv[1], _pv[2], _pv[3],
                                    _pv[4], _pv[5], _pv[6]) - 1. )

                # Add calculation of total model
                _pl_lc_tot = np.hstack(_pl_lc) # Need to update for multi-planet

        return _pl_lc_tot

    def lc_sim(self,quiet=False,return_lc=False,cadence=-1,exp_time=-1):
        """
        Light curve simulation.

        Parameters
        ---------
            quite :: Boolean
                For conditional rendering of some lines of code
            return_lc :: Boolean
            cadence ::
            exp_time ::

        Return 
            if return_lc:
                return t_grid, fl
        """
        # cadence: [s] cadence of observations. If < 0, use total exposure time (nstack * iframe * exptime)

        nsrc = len(self.SourceObj.gaia['ra'])

        if cadence < 0:
            cadence = self.nstack * self.exptime.to(u.second).value
        t_grid = np.arange(self.tstart.to(u.d).value, self.tend.to(u.d).value, cadence / (24. * 3600.))
        nt = len(t_grid)

        fl = np.zeros( (nt, nsrc) )
        fl_err = np.zeros_like(fl)

        # Calculate planet transit model
        # Flux decrease due to transit is included in _point_source_sim via target_flux_fraction parameter
        pl_lc = self.calc_pl_model(t_grid=t_grid,exp_time=exp_time) + 1.

        if quiet == False:
            print('Generating simulated lightcurve...')
            pbar = tqdm(total=nt)

        self._use_optimal_aperture()

        # Index 0 because we are concerned about the target only and the band information is already inputted in.
        scene_phot_count = self._photon_count(
            Gmag = self.SourceObj.gaia['Gmag'][0],
            Gmag_abs = self.SourceObj.gaia['Gmag_abs'][0],
            temp = self.SourceObj.gaia['Teff'][0],
            radius = self.SourceObj.gaia['radius'][0],
            metallicity = 0.0, #currently assumes solar
            logg = self.SourceObj.gaia['logg'][0]

        )[0]

        scene_phot_count *= self.exptime.to(u.second).value

        for n in range(nt):
            obs, err = self._point_source_sim(target_flux_fraction=pl_lc[n], scene_phot_count=scene_phot_count)
            i=0
            fl[n,i] = np.nansum(obs) 
            fl_err[n,i] = np.nansum(err)

            if quiet == False:
                pbar.update()
        if quiet == False:
            pbar.close()

        self.lc_t = t_grid
        self.lc_fl = fl
        self.lc_err = fl_err

        # Normalize to out-of-transit median
        # ti = np.where( np.abs(pl_lc - 1.) < 1.e-8 )[0]
        # if len(ti) > 0:
        #     ti = np.arange(len(pl_lc))
        ti = np.arange(len(pl_lc))
        for i in range(1):
            med_fl = np.median(self.lc_fl[ti,i])
            self.lc_fl[:,i] /= med_fl
            self.lc_err[:,i] /= med_fl

        if return_lc:
            return t_grid, fl

    def plot_lc(self,plot_model=True,exp_time=-1,
                    t_from_mid=True,t_unit='d'):
        """
        Rough plot. Need to update

        Parameters
        ----------
            plot_model :: Boolean
            exp_time
            t_from_mid :: Boolean
            t_unit :: str

        Return
        ------

        """
        t_offset = 0.
        if t_from_mid:
            t_offset = -(self.lc_t[-1] - self.lc_t[0]) / 2.

        t_scale = 1.
        if t_unit == 'hr':
            t_scale = 24.

        y_offset = -1.0
        y_scale = 1.e3

        plt.close(5)

        fig = plt.figure(num=5,figsize=(8,5),facecolor='w')
        ax = fig.add_subplot(111)


        ax.scatter( (self.lc_t + t_offset) * t_scale, (self.lc_fl[:,0] + y_offset) * y_scale,label='Simulated CASTOR')
        ax.errorbar( (self.lc_t + t_offset) * t_scale, (self.lc_fl[:,0] + y_offset) * y_scale, self.lc_err[:,0] * y_scale,ls='None')
        if t_from_mid:
            ax.set_xlabel(r'$t-t_0$ [' + t_unit + ']',fontsize=16)
        else:
            ax.set_xlabel('Time [d]',fontsize=16)
        ax.set_ylabel('Normalized Flux [ppt]',fontsize=16)
        xlim = np.array([self.lc_t[0], self.lc_t[-1]])
        xlim += t_offset
        xlim *= t_scale
        xlim += 0.5 * t_scale * (self.lc_t[1] - self.lc_t[0]) * np.array([-1.0, 1.0])
        ax.set_xlim(xlim)

        ax.xaxis.set_minor_locator(AutoMinorLocator())
        ax.yaxis.set_minor_locator(AutoMinorLocator())
        ax.tick_params(axis='both',labelsize=16)
        ax.tick_params(axis='both',which='major',length=6)
        ax.tick_params(axis='both',which='minor',length=3)
        ax.tick_params(axis='x',which='both',direction='inout')
        ax.tick_params(axis='x',which='both',top=True,direction='in')
        ax.tick_params(axis='y',which='both',direction='inout')
        ax.tick_params(axis='y',which='both',right=True,direction='in')

        if plot_model & hasattr(self,'pl_model'):
            _t_grid = np.linspace(self.lc_t[0],self.lc_t[-1],1000)
            pl_lc = self.calc_pl_model(t_grid=_t_grid,exp_time=exp_time) + 1.
            ax.plot( (_t_grid + t_offset) * t_scale,(pl_lc + y_offset) * y_scale,c='k',label='Transit Model')

            ax.legend(handletextpad=0.3,handlelength=2,
                                    labelspacing=0.2,facecolor='w',framealpha=0.0,
                                    fontsize=10,frameon=True,ncol=1)

        fig.tight_layout(pad=0.2)

        plt.ion()
        plt.show()
