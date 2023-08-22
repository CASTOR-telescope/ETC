"""
Grism calculations.

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

import astropy.units as u

from .conversions import calc_photon_energy, mag_to_flux, flam_to_photlam
from scipy.interpolate import interp1d
from astropy.modeling.models import Sersic2D


from os.path import join

from .filepaths import DATAPATH
import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
import spectres
from astropy.io import ascii

from castor_etc.background import Background
from castor_etc.sources import GalaxySource
from castor_etc.telescope import Telescope

class Grism(object):
    """
    Grism class.
    """
    
    def __init__(self,TelescopeObj, SourceObj, BackgroundObj):
        """
        Initialize class for Grism calculations.
        """
        
        if not isinstance(TelescopeObj, Telescope):
            raise TypeError("TelescopeObj must be a `castor_etc.Telescope` Object")
            
        if not isinstance(SourceObj, GalaxySource):
            raise TypeError("SourceObj must be a `castor_etc.GalaxySource`")
            
        if not isinstance(BackgroundObj, Background):
            raise TypeError("BackgroundObj must be a `castor_etc.Background`")
            
        # Assign attributes
        self.TelescopeObj = TelescopeObj
        self.SourceObj = SourceObj
        self.BackgroundObj = BackgroundObj

        # Initialize attributes that will be used in the future
        self.grism_channel = None
        self.grism_box = None
        self.sky_background_noise = None
        self.total_unif_noise = None
        self.grism_noise_total = None
        self.integrated_grism_box_count = None
        self.source_image = None
        self.source_seg = None
        
        self.add_bkgrd_noise_count = 0 #add noise per reads e-/pix
        self.exposure_time = None
        
        self.source_spectrum = [self.SourceObj.wavelengths.to(u.AA).value,self.SourceObj.spectrum]

    def _create_segmentation_map(self):
        """
        Internal function.
        
        """
        angle_a = self.SourceObj.angle_a #arcsec
        angle_b = self.SourceObj.angle_b #arcsec
        angle = self.SourceObj.rotation #rad
        n = self.SourceObj.sersic_index

        axial_ratio = angle_b.to(u.arcsec).value / angle_a.to(u.arcsec).value
        r_eff = np.sqrt(angle_a.to(u.arcsec).value * angle_b.to(u.arcsec).value) * u.arcsec
        e = np.sqrt(1 - (axial_ratio * axial_ratio))

        if r_eff.to(u.arcsec).value > 1.:
            raise ValueError("Effective radius of Sersic galaxy source should be less than or equal to 1 arcsec to do grism spectroscopy.")
        
        r_eff_pix = r_eff.to(u.arcsec).value / self.TelescopeObj.px_scale.value
        box_size = r_eff_pix * 4

        x,y = np.meshgrid(np.arange(box_size), np.arange(box_size))

        sersic_model = Sersic2D(amplitude=1, r_eff=r_eff_pix,
                                n=n,x_0=box_size/2,y_0=box_size/2,ellip=e,
                                theta=angle)
        
        self.source_image = sersic_model(x,y)

        self.source_seg = self.source_image / np.max(self.source_image) > 1e-2

    def disperse(self, grism_channel="u", check=True):
        """
        Grism disperser function

        source_image
            2D array of fluxes. Relative of absolute.
            Pixel scale needs to be CASTOR pixel scale (no oversampling).

        source_disperse_region
            Boolean 2D array, same size as source_image.
            Pixels with False will be masked.
            Pixels with True will be dispersed.

        source_spectrum
            Source spectrum, in flux densities.
            Source spectrum should be normalized (through spectrum.normalize_spectrum or equivalent), and have same wavelength grid as filter_transmission.

        grism_channel
            Which grism? "uv" or "u"

        """
        self.grism_channel = grism_channel
        
        self._create_segmentation_map()

        if self.source_image is None:
            source_image = np.copy(0)
        else:
            source_image = np.copy(self.source_image)

        if self.source_seg is None:
            source_disperse_region = np.copy(0)
        else:
            source_disperse_region = np.copy(self.source_seg)

        source_spectrum = np.copy(self.source_spectrum)

        #get filter transmission
        filter_transmission = ascii.read(join(DATAPATH,"passbands",f"passband_castor.{grism_channel}"))
        #get full grism throughtput
        grism_throughtput = ascii.read(join(DATAPATH,"grism_data",f"castor_grism_efficiency_.{grism_channel}.txt"))
        #get grism_dispersion
        #grism_dispersion = ascii.read(join(DATAPATH,"grism_data",f"grism_dispersion_{grism_channel}_extrap.txt"))
        grism_dispersion = ascii.read(join(DATAPATH,"grism_data",f"grism_dispersion_{grism_channel}.txt"))
        #get grism psf profile
        grism_psf_profile = ascii.read(join(DATAPATH,"grism_data","grism_approx_profile_uv.txt"))

        #Normalize image to get relative flux in each pixel (to scale spectum)
        norm = np.nansum(source_image[source_disperse_region])
        source_image_norm = source_image / norm

        #Get pixels of source that need to be dispered (i.e., pixels in non-masked region)
        nb_rows_img, nb_columns_img = source_image.shape
        r_i, c_i = np.where(source_disperse_region)
        indices = [(r,c) for r,c in zip(r_i, c_i)]

        #Create grism box that will be populated as we disperse the source sprectum onto it.
        margin_grism_profile = 10./100. #in percent
        margin_grism_dispersion = 5./100. #in percent
        #spatial profile extent
        nb_rows_grism = nb_rows_img + int(np.nanmax(grism_psf_profile['col2']) + np.abs(np.nanmin(grism_psf_profile['col2'])))
        nb_rows_grism = nb_rows_grism + int(nb_rows_grism*margin_grism_profile)*2
        #dispersion direction extent
        nb_columns_grism = int(np.nanmax(grism_dispersion['col1']))
        nb_columns_grism = nb_columns_grism + int(nb_columns_grism*margin_grism_dispersion)*2
        #create box
        self.grism_box = np.zeros((nb_rows_grism, nb_columns_grism))

        #grism offsets wrt direct imaging
        y_offset = int((nb_rows_grism - nb_rows_img) / 2) #in pixels
        #arbitrary offset in dispersion direction
        x_offset = 0  #in pixels

            #NOW USING FULL GRISM THROUGHTPUT
            #    #Multiply spectrum with filter response curve and grism efficiency curve
            #    wave, flux = source_spectrum
            #    grid_ftrans = np.interp(wave, filter_transmission['col1']*1e4, filter_transmission['col2'])
            #    flux *= grid_ftrans
            #
            #    #Use dummy grism_efficiency for now (flat 25% transmission)
            #    grism_efficiency = 0.25
            #    flux *= grism_efficiency
        wave, flux = np.copy(source_spectrum)
        #Converts ergs/cm2/s/A to photons/cm2/s/A
        spectrum_photon_density = flam_to_photlam(flux,wave) #in photons/cm2/s/A
        #Multiply spectrum (in photons/cm2/s/A) with full grism throughtput to get e-/cm2/s/A
        wavelength_key = grism_throughtput.keys()[0]
        thput_key = grism_throughtput.keys()[1]
        grid_gthput = np.interp(wave, grism_throughtput[wavelength_key]*10, 10**grism_throughtput[thput_key])
        spectrum_electron_density = spectrum_photon_density * grid_gthput   #"On detector" spectrum in e-/cm2/s/A
        flux = spectrum_electron_density

        #Resample spectrum to pixel dispersion
        #Pixel 1 corresponds to 3000 Angstrom for u channel and 1500 Angstrom for uv channel
        if grism_channel=="u":
            wave_zp = 3000
        if grism_channel=="uv":
            wave_zp = 1500
        wavelength_array = np.array([wave_zp+np.sum(grism_dispersion['col2'][:i]*10) for i in range(len(grism_dispersion['col2']))])
        #If spectrum doesn't fully overlap with wavelength_array, values of 0 are assumed.
        flux_resamp = spectres.spectres(wavelength_array, wave, flux, fill=0)

        if check:
            #show 1D spectrum as seen "on detector"
            plt.figure()
            plt.plot(source_spectrum[0],source_spectrum[1],'-k', label='Emitted Spectrum')
            plt.xlim(min(wavelength_array)*0.9, max(wavelength_array)*1.1)
            norm_y = 1.1*max(source_spectrum[1][(source_spectrum[0]>plt.gca().get_xlim()[0]) & (source_spectrum[0]<plt.gca().get_xlim()[1])])
            plt.ylim(0, norm_y)
            plt.plot(grism_throughtput[wavelength_key]*10, 10**grism_throughtput[thput_key] * norm_y, '-', color='grey', label='End-to-end grism throughtput', zorder=-1)
            plt.xlabel('Wavelength (angstroms)')
            plt.ylabel('Flux (ergs/cm2/s/A)')
            plt.legend()
            plt.figure()
            plt.plot(wavelength_array,flux_resamp,'-r', ds='steps-mid', label='"On Detector" Spectrum')
            plt.xlim(min(wavelength_array)*0.9, max(wavelength_array)*1.1)
            plt.xlabel('Wavelength (angstroms)')
            plt.ylabel('"On detector" Electron Flux (e-/cm2/s/A)')
            plt.show()

        #convert flux from electron flux densities (e-/cm2/s/A) to electron flux of e-/cm2/s (multiply by dispersion).
        #This is done here to be able to convert flux densities to counts (electrons) during the exposure.
        #If done at the exposure stage, then it complicates things as different wavelengths with different dispersions then overlap on the same pixels, which makes the conversion impracticle at this stage.
        flux_resamp *= grism_dispersion['col2']*10 #in e-/cm2/s

        #Make psf profile irradiance spectrum array. Currently assumes no PSF variation across the different grisms.
        psf_profile_pix_oversampled = grism_psf_profile['col2']
        psf_profile_pix = np.arange(np.min(psf_profile_pix_oversampled), np.max(psf_profile_pix_oversampled)+1, 1)
        psf_profile_oversampled = grism_psf_profile['col4']
        psf_profile = np.interp(psf_profile_pix, psf_profile_pix_oversampled, psf_profile_oversampled)
        #normalize to conserve total flux smeared due to psf profile
        norm = np.nansum(psf_profile)
        psf_profile_norm = psf_profile / norm
        #Smear spectrum using psf_profile in spatial direction
        spectrum_spatial = np.ones((len(psf_profile_norm),len(wavelength_array))) * flux_resamp * psf_profile_norm[:, None]
        
        #Loop over the pixel and disperse them (scale and add all 2D irradiance spectrum arrays).
        for indice in indices:

            #get spatial position (pixel indice) in the grism box of the source indice.
            y_indice = indice[0] + y_offset
            #get position of pixel1 in grism box
            x_indice = indice[1] + x_offset

            #Scale 2D array according to relative flux of each pixel in direct imaging.
            scale = source_image_norm[indice]
            spectrum_spatial_scaled = spectrum_spatial * scale

            #Add to grism box.
            y_len_spectrum_spatial = len(spectrum_spatial)
            y_start = int(y_indice - ((y_len_spectrum_spatial-1)/ 2))
            y_stop = y_start + y_len_spectrum_spatial
            x_len_spectrum_spatial = len(spectrum_spatial[0])
            x_stop = x_indice+x_len_spectrum_spatial
            #print(y_len_spectrum_spatial, y_stop, y_len_spectrum_spatial, y_indice)
            #print(x_len_spectrum_spatial, x_stop, x_len_spectrum_spatial, x_indice)
            #return grism_box, indice, spectrum_spatial

            self.grism_box[y_start:y_stop, x_indice:x_stop] += spectrum_spatial_scaled

    #        if check:
    #            plt.figure()
    #            plt.imshow(self.grism_box, aspect="auto", interpolation="none")
    #            plt.show()

        #mirror area
        mirror_diameter = 100 #cm
        mirror_area = np.pi * 0.25 * mirror_diameter**2 #cm2

        #Multiply by mirror area to get a count rate in e-/s
        self.grism_box *= mirror_area

        if check:
            plt.figure()
            plt.imshow(self.grism_box, aspect="auto", interpolation="none")
            plt.colorbar(label='"On detector" Count Rate (e-/s)')
            plt.xlabel('Pixels (Dispersion direction)')
            plt.ylabel('Pixels (Spatial direction)')
            plt.show()
            
    def expose(self, exposure_time=1000):
        """
        Function to simulate a noiseless grism observation for a given integration time.

        exposure_time
            Exposure time in seconds.

        """ 

        #self.grism_box is in e-/s
        #self.integrated_grism_box_count is then in units of e-
        self.exposure_time = exposure_time
        self.integrated_grism_box_count = self.grism_box * exposure_time

        #NEW METHOD GIVES 'self.grism_box * exposure_time' IN e- (INSTEAD OF ergs/cm2/A BEFORE).
        #NO LONGER NEEDED TO USE THE INVERSE SENSITIVITY FUNCTION.
        #integrated_grism_box_count is in electrons
    
    def _calc_sky_background_erate(self):
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
    
    def _calculate_tot_unif_noise(self, Nreads=1, Nbin=1):
        """
         
        Nreads: total number of read-outs
        Nbin: the number of detector pixels binned to one read-out pixel when on-chip binning is used.
        add_bkgrd_noise: additionnal noise, eg, HST: the background added using the post-flash option in e− pixel-1
        total_unif_noise is per pixel
        """
        self._calc_sky_background_erate() #Calculates e/sec/pix associated with the sky background 
        
        self.total_unif_noise = self.sky_background_noise[self.grism_channel] * self.exposure_time + self.TelescopeObj.dark_current * self.exposure_time + self.TelescopeObj.read_noise**2 * Nreads / Nbin + self.add_bkgrd_noise_count * Nreads
        
    def total_noise(self, Nreads=1, Nbin=1):
        """
        Function to generate the total noise of a grism observation for a given integration time.

        Nreads: total number of read-outs (int).

        Nbin: the number of detector pixels binned to one read-out pixel when on-chip binning is used (int).

        """
        
        self.grism_noise_total = np.zeros_like(self.grism_box)

        #compute the uniform background
        self._calculate_tot_unif_noise(Nreads=Nreads, Nbin=Nbin)

        #compute total background
        self.grism_noise_total = np.sqrt( self.integrated_grism_box_count + self.total_unif_noise
                                            )
        
    def show_2d_snr_per_resolution(self):
        """
        
        """

        plt.grid(False)
        plt.imshow(self.integrated_grism_box_count/self.grism_noise_total, aspect="auto",interpolation="none")
        plt.colorbar(label='SNR')
        plt.xlabel('Pixels (Dispersion direction)')
        plt.ylabel('Pixels (Spatial direction)')

    def show_1d_snr_per_resolution(self):
        """
        

        """

        box_center = int((self.integrated_grism_box_count.shape[0]-1) / 2)
        half_source_size = int((self.source_image.shape[0]-1) / 2)

        grism_1d = np.sum(self.integrated_grism_box_count[box_center-half_source_size:box_center+half_source_size+1,:],
                  axis=0)
        grism_1d_x = np.arange(0,len(grism_1d), 1)

        plt.figure()

        sum_signal_1d = np.sum(self.integrated_grism_box_count[box_center-half_source_size:box_center+half_source_size+1,:], 
                            axis=0)
        quad_error_1d = np.sqrt(np.sum(self.grism_noise_total[box_center-half_source_size:box_center+half_source_size+1,:]**2, 
                                    axis=0))
        snr_1d = sum_signal_1d / quad_error_1d

        plt.grid(False)
        plt.plot(grism_1d_x, snr_1d, '-k')
        plt.ylabel('SNR 1D ')
        plt.xlabel('Pixels (Dispersion direction)')"""
Grism calculations.

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

import astropy.units as u

from .conversions import calc_photon_energy, mag_to_flux, flam_to_photlam
from scipy.interpolate import interp1d
from astropy.modeling.models import Sersic2D


from os.path import join

from .filepaths import DATAPATH
import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
import spectres
from astropy.io import ascii

from castor_etc.background import Background
from castor_etc.sources import GalaxySource
from castor_etc.telescope import Telescope

class Grism(object):
    """
    Grism class.
    """
    
    def __init__(self,TelescopeObj, SourceObj, BackgroundObj):
        """
        Initialize class for Grism calculations.
        """
        
        if not isinstance(TelescopeObj, Telescope):
            raise TypeError("TelescopeObj must be a `castor_etc.Telescope` Object")
            
        if not isinstance(SourceObj, GalaxySource):
            raise TypeError("SourceObj must be a `castor_etc.GalaxySource`")
            
        if not isinstance(BackgroundObj, Background):
            raise TypeError("BackgroundObj must be a `castor_etc.Background`")
            
        # Assign attributes
        self.TelescopeObj = TelescopeObj
        self.SourceObj = SourceObj
        self.BackgroundObj = BackgroundObj

        # Initialize attributes that will be used in the future
        self.grism_channel = None
        self.grism_box = None
        self.sky_background_noise = None
        self.total_unif_noise = None
        self.grism_noise_total = None
        self.integrated_grism_box_count = None
        self.source_image = None
        self.source_seg = None
        
        self.add_bkgrd_noise_count = 0 #add noise per reads e-/pix
        self.exposure_time = None
        
        self.source_spectrum = [self.SourceObj.wavelengths.to(u.AA).value,self.SourceObj.spectrum]

    def _create_segmentation_map(self):
        """
        Internal function.
        
        """
        angle_a = self.SourceObj.angle_a #arcsec
        angle_b = self.SourceObj.angle_b #arcsec
        angle = self.SourceObj.rotation #rad
        n = self.SourceObj.sersic_index

        axial_ratio = angle_b.to(u.arcsec).value / angle_a.to(u.arcsec).value
        r_eff = np.sqrt(angle_a.to(u.arcsec).value * angle_b.to(u.arcsec).value) * u.arcsec
        e = np.sqrt(1 - (axial_ratio * axial_ratio))

        if r_eff.to(u.arcsec).value > 1.:
            raise ValueError("Effective radius of Sersic galaxy source should be less than or equal to 1 arcsec to do grism spectroscopy.")
        
        r_eff_pix = r_eff.to(u.arcsec).value / self.TelescopeObj.px_scale.value
        box_size = r_eff_pix * 4

        x,y = np.meshgrid(np.arange(box_size), np.arange(box_size))

        sersic_model = Sersic2D(amplitude=1, r_eff=r_eff_pix,
                                n=n,x_0=box_size/2,y_0=box_size/2,ellip=e,
                                theta=angle)
        
        self.source_image = sersic_model(x,y)

        self.source_seg = self.source_image / np.max(self.source_image) > 1e-2

    def disperse(self, grism_channel="u", check=True):
        """
        Grism disperser function

        source_image
            2D array of fluxes. Relative of absolute.
            Pixel scale needs to be CASTOR pixel scale (no oversampling).

        source_disperse_region
            Boolean 2D array, same size as source_image.
            Pixels with False will be masked.
            Pixels with True will be dispersed.

        source_spectrum
            Source spectrum, in flux densities.
            Source spectrum should be normalized (through spectrum.normalize_spectrum or equivalent), and have same wavelength grid as filter_transmission.

        grism_channel
            Which grism? "uv" or "u"

        """
        self.grism_channel = grism_channel
        
        self._create_segmentation_map()

        if self.source_image is None:
            source_image = np.copy(0)
        else:
            source_image = np.copy(self.source_image)

        if self.source_seg is None:
            source_disperse_region = np.copy(0)
        else:
            source_disperse_region = np.copy(self.source_seg)

        source_spectrum = np.copy(self.source_spectrum)

        #get filter transmission
        filter_transmission = ascii.read(join(DATAPATH,"passbands",f"passband_castor.{grism_channel}"))
        #get full grism throughtput
        grism_throughtput = ascii.read(join(DATAPATH,"grism_data",f"castor_grism_efficiency_.{grism_channel}.txt"))
        #get grism_dispersion
        #grism_dispersion = ascii.read(join(DATAPATH,"grism_data",f"grism_dispersion_{grism_channel}_extrap.txt"))
        grism_dispersion = ascii.read(join(DATAPATH,"grism_data",f"grism_dispersion_{grism_channel}.txt"))
        #get grism psf profile
        grism_psf_profile = ascii.read(join(DATAPATH,"grism_data","grism_approx_profile_uv.txt"))

        #Normalize image to get relative flux in each pixel (to scale spectum)
        norm = np.nansum(source_image[source_disperse_region])
        source_image_norm = source_image / norm

        #Get pixels of source that need to be dispered (i.e., pixels in non-masked region)
        nb_rows_img, nb_columns_img = source_image.shape
        r_i, c_i = np.where(source_disperse_region)
        indices = [(r,c) for r,c in zip(r_i, c_i)]

        #Create grism box that will be populated as we disperse the source sprectum onto it.
        margin_grism_profile = 10./100. #in percent
        margin_grism_dispersion = 5./100. #in percent
        #spatial profile extent
        nb_rows_grism = nb_rows_img + int(np.nanmax(grism_psf_profile['col2']) + np.abs(np.nanmin(grism_psf_profile['col2'])))
        nb_rows_grism = nb_rows_grism + int(nb_rows_grism*margin_grism_profile)*2
        #dispersion direction extent
        nb_columns_grism = int(np.nanmax(grism_dispersion['col1']))
        nb_columns_grism = nb_columns_grism + int(nb_columns_grism*margin_grism_dispersion)*2
        #create box
        self.grism_box = np.zeros((nb_rows_grism, nb_columns_grism))

        #grism offsets wrt direct imaging
        y_offset = int((nb_rows_grism - nb_rows_img) / 2) #in pixels
        #arbitrary offset in dispersion direction
        x_offset = 0  #in pixels

            #NOW USING FULL GRISM THROUGHTPUT
            #    #Multiply spectrum with filter response curve and grism efficiency curve
            #    wave, flux = source_spectrum
            #    grid_ftrans = np.interp(wave, filter_transmission['col1']*1e4, filter_transmission['col2'])
            #    flux *= grid_ftrans
            #
            #    #Use dummy grism_efficiency for now (flat 25% transmission)
            #    grism_efficiency = 0.25
            #    flux *= grism_efficiency
        wave, flux = np.copy(source_spectrum)
        #Converts ergs/cm2/s/A to photons/cm2/s/A
        spectrum_photon_density = flam_to_photlam(flux,wave) #in photons/cm2/s/A
        #Multiply spectrum (in photons/cm2/s/A) with full grism throughtput to get e-/cm2/s/A
        wavelength_key = grism_throughtput.keys()[0]
        thput_key = grism_throughtput.keys()[1]
        grid_gthput = np.interp(wave, grism_throughtput[wavelength_key]*10, 10**grism_throughtput[thput_key])
        spectrum_electron_density = spectrum_photon_density * grid_gthput   #"On detector" spectrum in e-/cm2/s/A
        flux = spectrum_electron_density

        #Resample spectrum to pixel dispersion
        #Pixel 1 corresponds to 3000 Angstrom for u channel and 1500 Angstrom for uv channel
        if grism_channel=="u":
            wave_zp = 3000
        if grism_channel=="uv":
            wave_zp = 1500
        wavelength_array = np.array([wave_zp+np.sum(grism_dispersion['col2'][:i]*10) for i in range(len(grism_dispersion['col2']))])
        #If spectrum doesn't fully overlap with wavelength_array, values of 0 are assumed.
        flux_resamp = spectres.spectres(wavelength_array, wave, flux, fill=0)

        if check:
            #show 1D spectrum as seen "on detector"
            plt.figure()
            plt.plot(source_spectrum[0],source_spectrum[1],'-k', label='Emitted Spectrum')
            plt.xlim(min(wavelength_array)*0.9, max(wavelength_array)*1.1)
            norm_y = 1.1*max(source_spectrum[1][(source_spectrum[0]>plt.gca().get_xlim()[0]) & (source_spectrum[0]<plt.gca().get_xlim()[1])])
            plt.ylim(0, norm_y)
            plt.plot(grism_throughtput[wavelength_key]*10, 10**grism_throughtput[thput_key] * norm_y, '-', color='grey', label='End-to-end grism throughtput', zorder=-1)
            plt.xlabel('Wavelength (angstroms)')
            plt.ylabel('Flux (ergs/cm2/s/A)')
            plt.legend()
            plt.figure()
            plt.plot(wavelength_array,flux_resamp,'-r', ds='steps-mid', label='"On Detector" Spectrum')
            plt.xlim(min(wavelength_array)*0.9, max(wavelength_array)*1.1)
            plt.xlabel('Wavelength (angstroms)')
            plt.ylabel('"On detector" Electron Flux (e-/cm2/s/A)')
            plt.show()

        #convert flux from electron flux densities (e-/cm2/s/A) to electron flux of e-/cm2/s (multiply by dispersion).
        #This is done here to be able to convert flux densities to counts (electrons) during the exposure.
        #If done at the exposure stage, then it complicates things as different wavelengths with different dispersions then overlap on the same pixels, which makes the conversion impracticle at this stage.
        flux_resamp *= grism_dispersion['col2']*10 #in e-/cm2/s

        #Make psf profile irradiance spectrum array. Currently assumes no PSF variation across the different grisms.
        psf_profile_pix_oversampled = grism_psf_profile['col2']
        psf_profile_pix = np.arange(np.min(psf_profile_pix_oversampled), np.max(psf_profile_pix_oversampled)+1, 1)
        psf_profile_oversampled = grism_psf_profile['col4']
        psf_profile = np.interp(psf_profile_pix, psf_profile_pix_oversampled, psf_profile_oversampled)
        #normalize to conserve total flux smeared due to psf profile
        norm = np.nansum(psf_profile)
        psf_profile_norm = psf_profile / norm
        #Smear spectrum using psf_profile in spatial direction
        spectrum_spatial = np.ones((len(psf_profile_norm),len(wavelength_array))) * flux_resamp * psf_profile_norm[:, None]
        
        #Loop over the pixel and disperse them (scale and add all 2D irradiance spectrum arrays).
        for indice in indices:

            #get spatial position (pixel indice) in the grism box of the source indice.
            y_indice = indice[0] + y_offset
            #get position of pixel1 in grism box
            x_indice = indice[1] + x_offset

            #Scale 2D array according to relative flux of each pixel in direct imaging.
            scale = source_image_norm[indice]
            spectrum_spatial_scaled = spectrum_spatial * scale

            #Add to grism box.
            y_len_spectrum_spatial = len(spectrum_spatial)
            y_start = int(y_indice - ((y_len_spectrum_spatial-1)/ 2))
            y_stop = y_start + y_len_spectrum_spatial
            x_len_spectrum_spatial = len(spectrum_spatial[0])
            x_stop = x_indice+x_len_spectrum_spatial
            #print(y_len_spectrum_spatial, y_stop, y_len_spectrum_spatial, y_indice)
            #print(x_len_spectrum_spatial, x_stop, x_len_spectrum_spatial, x_indice)
            #return grism_box, indice, spectrum_spatial

            self.grism_box[y_start:y_stop, x_indice:x_stop] += spectrum_spatial_scaled

    #        if check:
    #            plt.figure()
    #            plt.imshow(self.grism_box, aspect="auto", interpolation="none")
    #            plt.show()

        #mirror area
        mirror_diameter = 100 #cm
        mirror_area = np.pi * 0.25 * mirror_diameter**2 #cm2

        #Multiply by mirror area to get a count rate in e-/s
        self.grism_box *= mirror_area

        if check:
            plt.figure()
            plt.imshow(self.grism_box, aspect="auto", interpolation="none")
            plt.colorbar(label='"On detector" Count Rate (e-/s)')
            plt.xlabel('Pixels (Dispersion direction)')
            plt.ylabel('Pixels (Spatial direction)')
            plt.show()
            
    def expose(self, exposure_time=1000):
        """
        Function to simulate a noiseless grism observation for a given integration time.

        exposure_time
            Exposure time in seconds.

        """ 

        #self.grism_box is in e-/s
        #self.integrated_grism_box_count is then in units of e-
        self.exposure_time = exposure_time
        self.integrated_grism_box_count = self.grism_box * exposure_time

        #NEW METHOD GIVES 'self.grism_box * exposure_time' IN e- (INSTEAD OF ergs/cm2/A BEFORE).
        #NO LONGER NEEDED TO USE THE INVERSE SENSITIVITY FUNCTION.
        #integrated_grism_box_count is in electrons
    
    def _calc_sky_background_erate(self):
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
    
    def _calculate_tot_unif_noise(self, Nreads=1, Nbin=1):
        """
         
        Nreads: total number of read-outs
        Nbin: the number of detector pixels binned to one read-out pixel when on-chip binning is used.
        add_bkgrd_noise: additionnal noise, eg, HST: the background added using the post-flash option in e− pixel-1
        total_unif_noise is per pixel
        """
        self._calc_sky_background_erate() #Calculates e/sec/pix associated with the sky background 
        
        self.total_unif_noise = self.sky_background_noise[self.grism_channel] * self.exposure_time + self.TelescopeObj.dark_current * self.exposure_time + self.TelescopeObj.read_noise**2 * Nreads / Nbin + self.add_bkgrd_noise_count * Nreads
        
    def total_noise(self, Nreads=1, Nbin=1):
        """
        Function to generate the total noise of a grism observation for a given integration time.

        Nreads: total number of read-outs (int).

        Nbin: the number of detector pixels binned to one read-out pixel when on-chip binning is used (int).

        """
        
        self.grism_noise_total = np.zeros_like(self.grism_box)

        #compute the uniform background
        self._calculate_tot_unif_noise(Nreads=Nreads, Nbin=Nbin)

        #compute total background
        self.grism_noise_total = np.sqrt( self.integrated_grism_box_count + self.total_unif_noise
                                            )
        
    def show_2d_snr_per_resolution(self):
        """
        
        """

        plt.grid(False)
        plt.imshow(self.integrated_grism_box_count/self.grism_noise_total, aspect="auto",interpolation="none")
        plt.colorbar(label='SNR')
        plt.xlabel('Pixels (Dispersion direction)')
        plt.ylabel('Pixels (Spatial direction)')

    def show_1d_snr_per_resolution(self):
        """
        

        """

        box_center = int((self.integrated_grism_box_count.shape[0]-1) / 2)
        half_source_size = int((self.source_image.shape[0]-1) / 2)

        grism_1d = np.sum(self.integrated_grism_box_count[box_center-half_source_size:box_center+half_source_size+1,:],
                  axis=0)
        grism_1d_x = np.arange(0,len(grism_1d), 1)

        plt.figure()

        sum_signal_1d = np.sum(self.integrated_grism_box_count[box_center-half_source_size:box_center+half_source_size+1,:], 
                            axis=0)
        quad_error_1d = np.sqrt(np.sum(self.grism_noise_total[box_center-half_source_size:box_center+half_source_size+1,:]**2, 
                                    axis=0))
        snr_1d = sum_signal_1d / quad_error_1d

        plt.grid(False)
        plt.plot(grism_1d_x, snr_1d, '-k')
        plt.ylabel('SNR 1D ')
        plt.xlabel('Pixels (Dispersion direction)')