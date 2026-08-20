"""
UVMOS Spectroscopy calculations.

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

import astropy.units as u
import astropy
from regions import RectangleSkyRegion, RectanglePixelRegion, PixCoord 
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
from astropy import constants as const
from astropy.wcs import WCS
from scipy.interpolate import RectBivariateSpline, interp1d
import scipy
from scipy.optimize import curve_fit
import math
from os.path import join
from astropy.coordinates import angular_separation, Angle, SkyCoord
from astropy.stats import gaussian_fwhm_to_sigma
from scipy.ndimage import rotate
from astropy.coordinates import offset_by


from castor_etc.background import Background
from castor_etc.sources import PointSource
from castor_etc.telescope import Telescope

from castor_etc import DATAPATH

import sys

def Gaussian(x, x0, sigma, a):
    return a * np.e**(- (x-x0)**2 / (2*sigma**2)  )

def Gaussian2D(x, y, sigma, a=1, x0=0, y0=0): 
    term1 = ( (x-x0)**2 )/(2*sigma**2)  
    term2 = ( (y-y0)**2 )/(2*sigma**2)  
    return a * np.e**( - (term1 + term2 ) )

class UVMOS_Spectroscopy:
    """
    UVMOS_Spectroscopy class.
    """

    def __init__(self, TelescopeObj, SceneObj, BackgroundObj, deltaXY, **kwargs): 
        """
        Initializes the UVMOS_Spectroscopy class

        params
        ------
        TelescopeObj = The telescope object.  Must be a TelescopeObj class object.
        
        SceneObj = The scene object.  Must be a Scene class object.
        
        BackgroundObj = The background object.  Must be a BackgroundObj class object.
                    
        deltaXY = (array) The distance in pixels from the center of field to the centre of the slit.  A numpy array of coordinates, the first column is X and                   the second column if Y

        kwargs (optional params)
        ------
        thetaList = (list) A list of angles that can range from 0 to 180 degrees.  If thetaList==None, then all slits will have an orientation
                    of 0 degrees.

        case = (int) The two different dmd orientations.  Case 1 is for the DMD at a 45 degree tilt and case 2 is for the DMD at a 0 degree tilt.  
                See documentation for more details.

        """
        # Get kwargs
        thetaList = kwargs.get("theta", np.zeros(np.shape(deltaXY)[0])) # defaults to 0 degree angle
        case = kwargs.get("case", 1) # defaults to case 1 if not otherwise specified
        slit_source_name = kwargs.get("source_name", ["source_1"]) # defaults to source_1 is not otherwise specified

        # Check if case has the correct type and value
        if type(case) != int:
            raise TypeError("case must be an integer")
        if case != 1 and case !=2:
            raise ValueError("case must be 1 or 2")

        #Check TelescopeObj
        if not isinstance(TelescopeObj, Telescope):
            raise TypeError("TelescopeObj must be a `castor_etc.Telescope` object")

        #Initialize BackgroundObj
        if not isinstance(BackgroundObj, Background):
            raise TypeError("BackgroundObj must be a `castor_etc.Background` object")

        # Check if source_name has the correct type and is inputted when there are multiple sources.
        for i in slit_source_name:
            if type(i) != str:
                raise TypeError("source_name must be a string")
        if np.shape(deltaXY)[0] > 1:
            if slit_source_name == ["source_1"]:
                raise ValueError("Source names need to be inputted when you have multiple sources")

        # Check if thetaList has the correct dimensions and the correct type.
        if thetaList is not None:
            if len(thetaList)!=np.shape(deltaXY)[0]:
                raise ValueError("thetaList must have the same length as SourceList")

            for theta in thetaList:
                if type(theta)!=float and type(theta)!=int:
                    raise TypeError("thetaList must be comprised of floats or integers")
                elif theta<0 or theta>180:
                    raise ValueError("theta must be between 0 and 180 degrees")

        # Initialize inputs
        self.TelescopeObj = TelescopeObj
        self.SceneObj = SceneObj 
        self.BackgroundObj = BackgroundObj
        self.case = case
        self.deltaXY = deltaXY
        self.slit_source_name = slit_source_name

        if thetaList is not None:
            self.thetaList = thetaList 
        else:
            self.thetaList = None

        # Unpack Scene class
        self.detector_dimensions = np.shape(self.SceneObj.displayScene(self.TelescopeObj, passband = "uv"))
        self.source_array = self.SceneObj.displayScene(self.TelescopeObj, passband = "uv")
        
        #  These attributes may later be moved to the Telescope object
        self.gain = 1
        self.min_wave = (150 * u.nm).to(u.AA) # Minimum wavelength 
        self.max_wave = (300 * u.nm).to(u.AA) # Maximum wavelength 
        self.slit_width =  0.214 * u.arcsec
        self.slit_height = 1 * u.arcsec

        self.FWHM = TelescopeObj.fwhm 
        self.dispersion = 0.061 * u.nm
        self.pixel_size = TelescopeObj.px_scale #0.1 arcsec/pixel
        
        self.source_detector = None
        self.background_detector = None
        self.source_CASTORSpectrum = None
        self.background_CASTORSpectrum = None
        self.waves_CASTORSpectrum  = None
        self.waves_CASTORSpectrumBackground = None
        
        self.slit_width_pix = None
        self.slit_height_pix = None
        self.source_extracted_numpixs = None
        self.background_extracted_numpixs = None

        # The DMD attributes should be part of the parameters script.
        self.dmd_x_dimen_mm = 20.7 * u.mm
        self.dmd_y_dimen_mm = 11.7 * u.mm
        self.dmd_x_dimen_arc = 213 * u.arcsec
        self.dmd_y_dimen_arc = 121 * u.arcsec
        self.detector_pix_width = 10 * u.nm
        self.plate_scale = 0.097 * u.mm / u.arcsec

        self.detectorDMDScale = 10.8 / 10.0  # DMD pixel width / detector pixel width
        self.DMDPixelDimen = (1920, 1080) # pixels

        self.slit_width = 0.214 * u.arcsec # Specify slit dimensions
        self.slit_height = 1 * u.arcsec
        
        self.slit_width_pix = math.ceil((self.slit_width/self.pixel_size).value) # Generate pixel dimensions
        self.slit_height_pix = math.ceil((self.slit_height/self.pixel_size).value) 

        self.pixel_FWHM = self.FWHM.value * self.plate_scale.value / (self.detector_pix_width.value * 10e-3)

    
    def _specify_DMD_FOV(self):
        """
        Creates an array the size of the DMD

        returns
        -------
        dmdMask = (array) array the dimensions of the DMD
        """
        
        # case 1, DMD is on 45 degree angle, dispersion is aligned with detector 
        if self.case == 1:

            self.detectorArray = np.zeros((self.detector_dimensions[0], self.detector_dimensions[1]))
            self.centre = PixCoord(round(self.detector_dimensions[0] / 2), round(self.detector_dimensions[1] / 2))
            self.dmdArea = RectanglePixelRegion(self.centre, self.DMDPixelDimen[0], self.DMDPixelDimen[1], angle = 45.0 * u.deg)
            dmdMask = self.dmdArea.to_mask(mode = "center") 

        else: # case 2, DMD is aligned with the detector, dispersion is at a 45 degree angle

            self.detectorArray = np.zeros((self.detector_dimensions[1], self.detector_dimensions[0])) # Have the detector be oriented vertically
            self.centre = PixCoord(round(self.detector_dimensions[1] / 2), round(self.detector_dimensions[0] / 2))
            self.dmdArea = RectanglePixelRegion(self.centre, self.DMDPixelDimen[0], self.DMDPixelDimen[1], angle = 90.0 * u.deg)
            dmdMask = self.dmdArea.to_mask(mode = "center")

        return dmdMask

    def _get_slit_coordinates(self):
        '''
        Makes a catalog of sources with their ra and decs in pixel coordinates.

        returns
        -------
        sourceCatalog = (list of PixCoords) the source coordinates in pixels
        '''

        dmd = self._specify_DMD_FOV() # this defines self.centre

        sourceCatalog = []
        
        for sourceNum in range(np.shape(self.deltaXY)[0]):
            delta_x = self.deltaXY[sourceNum][0]
            delta_y = self.deltaXY[sourceNum][1]

            x_cent = self.centre.x 
            y_cent = self.centre.y

            x = x_cent + delta_x
            y = y_cent + delta_y
                
            sourceCoord = PixCoord(round(x), round(y))
            sourceCatalog = np.append(sourceCatalog, sourceCoord)       

        return sourceCatalog

    
    def _generate_slits(self):   
        """
        Creates a list of slits (RectanglePixelRegions)

        returns
        -------
        slitCatalog = (list of RectanglePixelRegions) the slit regions
        """
   
        slitCatalog = []

        slitCoords = self._get_slit_coordinates()

        if self.thetaList != None:
            for source in range(len(slitCoords)):
                if self.case == 1:
                    slitRegion = RectanglePixelRegion(center = slitCoords[source], height = self.slit_height_pix, width = self.slit_width_pix, angle = (self.thetaList[source] + 45) * u.deg)
                else:
                    slitRegion = RectanglePixelRegion(center = slitCoords[source], height = self.slit_height_pix, width = self.slit_width_pix, angle = self.thetaList[source] * u.deg)
                slitCatalog = np.append(slitCatalog, slitRegion)
        else:
            for source in range(len(slitCoords)):
                if self.case == 2:
                    slitRegion = RectanglePixelRegion(center = slitCoords[source], height = self.slit_height_pix, width = self.slit_width_pix, angle = self.thetaList * 45 * u.deg)
                else:
                    slitRegion = RectanglePixelRegion(center = slitCoords[source], height = self.slit_height_pix, width = self.slit_width_pix, angle = self.thetaList * u.deg)
            
                slitCatalog = np.append(slitCatalog, slitRegion)

        return slitCatalog
        

    def _check_slits(self):   
        """
        Checks that no slits are off the detector and that no slits are overlapping.
        """

        slitCatalog = self._generate_slits()
        dmdFOV = self._specify_DMD_FOV()

        detector_array = np.zeros((self.detector_dimensions[0], self.detector_dimensions[1]))
        dmdMask = dmdFOV.to_image((self.detector_dimensions[0], self.detector_dimensions[1]))
        dmdMask = np.where(dmdMask == 0, 1, 0)
        detector_array = detector_array + dmdMask

        for i in range(len(slitCatalog)):
            slitMask = slitCatalog[i].to_mask(mode = "center") 
            slitMask = slitMask.to_image((self.detector_dimensions[0], self.detector_dimensions[1]))
            detector_array = detector_array + slitMask
               
        results = np.where(detector_array > 1)[0]
        if len(results) > 0:
            raise ValueError("Two or more slits are overlapping or are outside of the field of view of the DMD")


    def show_slits(self):
        """
        Shows the position and orientation of the slits and sources on the detector.
        """

        self._check_slits()

        slitCatalog = self._generate_slits()
        dmdMask = self._specify_DMD_FOV() 

        fig, ax = plt.subplots(figsize = (10, 10))

        dmdMask = dmdMask.to_image((self.detector_dimensions[0], self.detector_dimensions[1]))
        
        for slit in range(len(slitCatalog)):

            # Plot the slits
            bottom_corner = [slitCatalog[slit].center.x - self.slit_width_pix/2, slitCatalog[slit].center.y - self.slit_height_pix/2]
            slit_patch = patches.Rectangle(bottom_corner, self.slit_width_pix, 
                                           self.slit_height_pix, angle = self.thetaList[slit],
                                           rotation_point = "center", edgecolor = "r", fill = False)
            ax.add_patch(slit_patch)

        self.source_array = np.where(dmdMask == 1, self.source_array, - np.max(self.source_array) / 4)

        plt.imshow(self.source_array, cmap = "gray", origin = "lower")
        plt.colorbar()
        plt.show()

    def _calc_slit_transmission(self, print_transmission_fact=False): 
        """
        Calculates the slit transmission percentage.

        params
        ------
        print_transmission_fact = (bool) if True, prints the transmission.  

        returns
        -------
        transmission = (list of floats) the transmissions for each slit
        """

        slitCatalog = self._generate_slits()
        sourceCatalog = self._get_slit_coordinates()

        transmissions = []

        for i in range(len(sourceCatalog)):
            self.size = 10

            sourceImg = self.source_array 
            source_cutout = sourceImg[sourceCatalog[i].y - 2 : sourceCatalog[i].y + self.size * 2 - 1, sourceCatalog[i].x - 2 : sourceCatalog[i].x + self.size * 2 - 1] # Gaussian without the slit

            x = np.linspace(1, np.shape(source_cutout)[0], np.shape(source_cutout)[0])
            y = np.linspace(1, np.shape(source_cutout)[1], np.shape(source_cutout)[1])

            slitMask = slitCatalog[i].to_mask(mode = "center") 
            
            slitImg = slitMask.to_image((self.detector_dimensions[0], self.detector_dimensions[1])) 
            
            img_slit = np.where(slitImg == 0, 1, sourceImg) 
            slit_cutout = img_slit[sourceCatalog[i].y - self.size : sourceCatalog[i].y + self.size + 1, sourceCatalog[i].x - self.size : sourceCatalog[i].x + self.size + 1] # The slit mask

            interp_spline = RectBivariateSpline(x, y, source_cutout, kx=1, ky=1) # Make the spline. 

        # ------- Total transmitted intensity
            vol1 = interp_spline.integral(-self.size, self.size, -self.size, self.size)
                                            
            pix_weights = np.zeros((self.size*2 + 1, self.size*2 + 1))
    
            vals = []
            for j in range(-self.size, self.size + 1):    
                for k in range(-self.size, self.size + 1):

                # Intensity transmitted to this "pixel"
                    vol2 = interp_spline.integral(j - 1, j + 1, k - 1, k + 1)
                    pix_weights[j + self.size - 2, k + self.size - 2] = pix_weights[j + self.size - 2, k + self.size - 2] + vol2/vol1 
            
            pix_dist = np.where(slit_cutout == 1, 0, pix_weights) 

            vol1 = np.sum(pix_weights)
            vol2 = np.sum(pix_dist)

            fractional_slit_transmission_2D = vol2 / vol1
            transmissions.append(fractional_slit_transmission_2D)

            if print_transmission_fact:
                print('Fractional slit transmission (2D estimation) {}'.format(fractional_slit_transmission_2D))
                print('Fractional slit loss (2D estimation) {}'.format(1 - fractional_slit_transmission_2D))

        return transmissions


    def calc_source_pix_weights(self):
        """
        Determines the source pixel weights

        returns
        -------
        pix_dist_list = (list of arrays) list of the slit arrays
        """

        slitCatalog = self._generate_slits()
        sourceCatalog = self._get_slit_coordinates()

        pix_dist_list = []

        for i in range(len(sourceCatalog)):

            sourceImg = self.source_array
            source_cutout = sourceImg[sourceCatalog[i].y - self.size - 1 : sourceCatalog[i].y + self.size, sourceCatalog[i].x - self.size -1 : sourceCatalog[i].x + self.size] 

            x = np.linspace(1, np.shape(source_cutout)[0], np.shape(source_cutout)[0])
            y = np.linspace(1, np.shape(source_cutout)[1], np.shape(source_cutout)[1])

            slitMask = slitCatalog[i].to_mask(mode = "center")
            slitImg = slitMask.to_image((self.detector_dimensions[0], self.detector_dimensions[1])) 
            img_slit = np.where(slitImg == 0, slitImg, sourceImg) 
            slit_cutout = img_slit[sourceCatalog[i].y - self.size - 1 : sourceCatalog[i].y + self.size, sourceCatalog[i].x - self.size - 1 : sourceCatalog[i].x + self.size] 

            interp_spline = RectBivariateSpline(x, y, source_cutout, kx=1, ky=1)

        # ------- Total transmitted intensity
            vol1 = interp_spline.integral(0, self.size * 2 + 1, 0, self.size * 2 + 1)
                                            
        # -------
            pix_dist = np.zeros(np.shape(slit_cutout))
            # This is the per pixel transmitted intensity
            # Then divide the per pixel transmitted intensity by the total transmitted intensity
            vals = []
            for j in range(0, self.size * 2 + 1):    
                for k in range(0, self.size * 2 + 1):

                # Intensity transmitted to this "pixel"
                    vol2 = interp_spline.integral(j - 1, j + 1, k - 1, k + 1)
                    pix_dist[j - 1, k - 1] = pix_dist[j - 1, k - 1] + vol2/vol1 # was 2
            
            pix_dist = np.where(pix_dist < 0, -1 * pix_dist, pix_dist)
            pix_dist = np.where(slit_cutout == 0, 0, pix_dist)
            pix_dist_list.append(pix_dist)
         
        return pix_dist_list
    
    def show_source_pix_weights(self):
        """
        This plots the results of the calc_source_pix_weights() function.
        """
        
        pix_dist_list = self.calc_source_pix_weights()

        for i in range(len(pix_dist_list)):
        # -----------------------------
            fig = plt.figure()
            ax = plt.subplot(111)

            plt.title('Source viewed on the detector')
            plt.xlabel('Pixel')
            plt.ylabel('Pixel')

            im = plt.imshow(pix_dist_list[i], origin='lower', cmap='viridis', extent = (-self.size, self.size+1, -self.size, self.size+1))

            plt.grid()

            cbar = plt.colorbar( )
            cbar.set_label('Normalized Intensity')
            ax.set_aspect('equal')

            plt.savefig("source_pix_weights_"+str(i)) 

    
    def show_slit_image(self, wavelength): 
        """ 
        Shows the smeared out slit.  This is accomplished by rotating the slit to be vertical, then smearing out the slit           horizontally, then rotating the smeared out slit back to it's original position.

        params
        ------
        wavelength = (float) wavelength in angstroms

        returns
        -------
        detector_list = (list of arrays) a list of the smeared out slits

        center_pix_list = (list of int) a list of the center pixels of the smeared out slits
        """

        pix_dist_list = self.calc_source_pix_weights() 

        detector_list = []
        center_pix_list = []

        for i in range(len(pix_dist_list)):

            source_key = self.slit_source_name[i]

            count = 0
            for j in range(len(self.SceneObj.sources_list)):
                if self.SceneObj.sources_list[j][-1] == source_key:
                    y = self.SceneObj.sources_list[j][0].spectrum
                    x = self.SceneObj.sources_list[j][0].wavelengths.value
                    count += 1
            if count == 0:
                raise ValueError("Inputted sources names must match inputted source name in the scene class")
                    

            # ----- Extract only the wavelengths in the spectral region
            ind = np.where( (x >= self.min_wave.value) & (x <= self.max_wave.value) )
            xx = x[ind]
            yy = y[ind] 

            # ------- Determine the fractional slit transmission
            fractional_slit_transmission = self._calc_slit_transmission()[i]

        # --------- Find detected flux
            fact = 1/(const.h*const.c).to('erg angstrom')
            T = self._getTransmission(xx) 
            yy_phot = yy*(xx*(fact.value))*self.TelescopeObj.mirror_area.value*T*fractional_slit_transmission*self.gain # Detected flux in counts / s / angstrom

            interp_source_spectrum = interp1d(xx, yy_phot, kind='cubic',  ) 

            # --------- If the inputted wavelength is near the extremes, i.e., (1500,2995) nm of the spectrum, then in order to calculate the bands, we have to shift the range while creating detector heat map.
            
            wave = (wavelength * u.nm).to(u.AA).value
            if (wave - self.size < min(xx)): 
                wave = wave + self.size + self.dispersion.to(u.AA).value
        
            if (wave + self.slit_width_pix > max(xx)):
                wave = wave - self.size - self.dispersion.to(u.AA).value

            # -------- Create the detector heat map
            if self.case == 1:
                angle = 45 - self.thetaList[i] 
            else:
                angle = 90 - self.thetaList[i]

            straight_detector = rotate(pix_dist_list[i], angle = angle, reshape = True, mode = "constant")

            slit_idxs = []

            for k in range(len(straight_detector[0, :])): # 
                slit_idx = np.where(straight_detector[k, :] > 0)[0]
                for j in slit_idx:
                    slit_idxs.append(j) 

            max_slit = 0 
            min_slit = 50 

            for j in slit_idxs:
                if max_slit < j:
                    max_slit = j
                if min_slit > j:
                    min_slit = j

            max_slit_width = max_slit - min_slit - 1

            straight_detector = straight_detector[:, min_slit:max_slit - 1]

            bands = np.arange(wave - max_slit_width, wave + max_slit_width, self.dispersion.to('angstrom').value)
            tot_xpixels = len(bands) + max_slit_width


            detector = np.zeros((np.shape(straight_detector)[0], tot_xpixels)) 

            for k in range(len(bands)-1):

                if (bands[k] < wave < bands[k+1]):
                    center_pix = k
                    center_pix_list.append(center_pix)

                x = np.linspace(bands[k], bands[k+1], 1000 )
                y = interp_source_spectrum(x) 
                y = np.where(y < 0, 0, y) 

                # Flux in photons per second
                phot_in_band = np.array(scipy.integrate.cumulative_trapezoid(y, x )[-1])

                # Create the pixel map for this band
                pix_map = straight_detector * phot_in_band

                for j in range(max_slit_width): 
                    detector[:,k+j] = detector[:,k+j] + pix_map[:,j]

            rotated_detector = rotate(detector, angle = -angle, reshape = True, mode = "constant")
            detector_list.append(rotated_detector)
            
        return detector_list, center_pix_list
        

    def _extraction(self, detector, pix_waves, extraction_width, extraction_lowerlim, extraction_upperlim):
        """
        Extracts the wavelengths, flux, and number of pixels from the smeared slit.

        params
        ------
        detector = (array of floats) 
        
        pix_waves = ()

        extraction_width = (int) width of the extraction box along the wavelength axis in pixels.  By default this is 1.
        
        extraction_lowerlim = (int) Lower limit of the extraction box along the spatial axis in pixels.
        
        extraction_upperlim = (int) Upper limit of the extraction box along the spatial axis in pixels.

        returns
        -------

        extracted_waves = (array of floats) wavelengths
        
        extracted_flux = (array of floats) the flux extracted from the slit 
        
        extracted_numpix = (int) number of pixels in slit
        
        """

        # --- Warnings
        if extraction_lowerlim < 1:
            raise TypeError('Lower limit of extraction box cannot be negative.')

        if extraction_upperlim > np.shape(detector)[1]: 
            raise TypeError('Upper limit of the extraction box cannot be greater than the slit height in pixels')

        if extraction_width > np.shape(detector)[0]:
            raise TypeError('The width of the extraction box cannot be greater than the total number of pixels along the wavelength axis')

        if not isinstance( extraction_lowerlim, int ):
            raise TypeError('The lower edge of the extraction box must be an integer')

        if not isinstance( extraction_upperlim, int ):
            raise TypeError('The upper edge of the extraction box must be an integer')

        if not isinstance( extraction_width, int ):
            raise TypeError('The width of the extraction box must be an integer')


        # ---------
        num_extractions = math.floor( np.shape(detector)[1] / extraction_width) 

        extracted_numpix = (extraction_width*( extraction_upperlim - extraction_lowerlim +1 ) ) * num_extractions

        extracted_waves = [ np.mean(pix_waves[i:i + extraction_width]) for i in range(0, len(pix_waves), extraction_width) if len(pix_waves[i:i + extraction_width])== extraction_width ]

        extracted_flux = [] # counts / extraction box

        for i in range(num_extractions): 

            x_left = extraction_width* i
            x_right = extraction_width* i + extraction_width

            extracted_flux.append(sum(detector[extraction_lowerlim - 1:extraction_upperlim,x_left:x_right].flatten()) )

        return extracted_waves, extracted_flux, extracted_numpix

    
    def _getTransmission(self, x):
        """
        Unpacks spectrographTransmission_CASTOR.txt.  Fits it with a cubic spline.

        params
        ------
        x = (array) an array of wavelengths

        returns
        -------
        the transmission values from the cubic spline
        """

        transmissionpath = "spectrographTransmission_CASTOR.txt"
        lines = open(transmissionpath).readlines() #  NOTE: Could be moved to 'data' folder  later
        trans_w = [float(line.split()[0]) * 10 for line in lines if '#' not in line.split()]
        trans_T = [float(line.split()[1]) for line in lines if '#' not in line.split()]

        func = interp1d( trans_w, trans_T, kind='cubic'  )

        return func(x)

    def showTransmission(self):
        """
        Plots the transmission as a function of wavelength
        """
        x = np.linspace( self.min_wave.value, self.max_wave.value, 1000 )

        plt.figure(figsize=(10, 5))
        plt.plot( x, self._getTransmission(x), color='tab:red' )
        plt.title("CASTOR Spectrographic Throughput")
        plt.ylabel('Efficiency')
        plt.xlabel(r"Wavelength [$\AA$]")

    def _getDispersion(self, x):
        """
        Unpacks dispersion_resolution2.dat and fits it with a cubic spline.

        params
        ------
        x = (array) array of wavelengths in angstroms

        returns
        -------
        the dispersion values from the cubic spline
        """

        dispersionpath = "dispersion_resolution2.dat"
        lines = open(dispersionpath).readlines() #  NOTE: Could be moved to 'data' folder  later
        disp_w = [float(line.split()[0]) for line in lines if '#' not in line]
        disp_v = [float(line.split()[2]) for line in lines if '#' not in line]

        func = interp1d( disp_w, disp_v, kind='cubic', fill_value='extrapolate'  )

        return func(x)

    
    def _calc_sigmaPix(self, dispersion):
        """
        Calculate the sigma value of a Gaussian fit to a pixel with a width in wavelength equal
        to the provided dispersion

        params
        ------
        dispersion = ()

        returns
        -------
        absolute value of cf_sigma = ()
        """

        # ----- Build the gate function representing a pixel
        w = np.linspace(- dispersion / 2, dispersion / 2)
        f = np.array([1] * len(w))

        w1 = np.linspace(- dispersion, - dispersion / 2)
        w2 = np.linspace(dispersion / 2, dispersion)

        f1 = np.array([0] * len(w1))
        f2 = np.array([0] * len(w2))

        gate_w = np.concatenate((w1, w, w2))
        gate_f = np.concatenate((f1, f, f2))

        # ----- Fit a gaussian to the gate function
        cf_x0, cf_sigma, cf_a = curve_fit(Gaussian, gate_w, gate_f)[0]  

        return abs(cf_sigma)

    
    def _calc_sigmaPSF(self, dispersion):
        """
        Calculate the sigma value of the Gaussian PSF on the detector

        params
        ------
        dispersion = ()

        returns
        -------
        sigma_psf_wave = ()
        """

        sigma_psf = self.FWHM / (2 * np.sqrt(2 * np.log(2)))
        sigma_psf_wave = sigma_psf * (1 / self.pixel_size) * dispersion

        return sigma_psf_wave

    
    def _calcR(self, x, disp ):
        """
        Calculate the resolving power

        Parameters
        ----------
        x = array of wavelength values
        
        disp = array of dispersion values corresponding to the wavelenghts in x

        Returns
        -------
        Array of resolving power values corresponding to x
        """

        fact = 2 * np.sqrt(2 * np.log(2))
        dL =  [ fact *  np.sqrt( self._calc_sigmaPSF(d)**2 + self._calc_sigmaPix(d)**2  )for d in disp]

        return [x[i] / dL[i] for i in range(len(x))]

    def showResolvingPower(self):

        minn = self.min_wave.to('nm').value
        maxx = self.max_wave.to('nm').value
        x = np.linspace(minn, maxx, 500)
        dispersion = self._getDispersion(x)
        R = self._calcR(x, dispersion)

        fig = plt.figure(figsize=(10,4))
        ax = plt.subplot(111)
        plt.plot(x * 10, R)

        plt.title('Resolving power as a function of wavelength for the UVMOS spectrograph')
        plt.ylabel('Resolving Power (R)')
        plt.xlabel(r'Wavelength [$\AA$]')  

    
    def calc_source_CASTORSpectrum(self, extraction_width=1, extraction_lowerlim=1 , extraction_upperlim='max' ):
        """
        Calculates the source spectrum.  Done for non-vertical slits by rotating the slit to be vertical and then smearing
        out the slit horizontally.
        
        params
        ------
        
        extraction_width = (int) width of the extraction box along the wavelength axis in pixels.  By default this is 1.
        
        extraction_lowerlim = (int) Lower limit of the extraction box along the spatial axis in pixels.
        
        extraction_upperlim = (int) Upper limit of the extraction box along the spatial axis in pixels.
        """

        self.waves_CASTORSpectrum_list = []
        self.source_CASTORSpectrum_list = []
        self.source_extracted_numpixs_list = []
        self.source_detector_list = []
        self.rot_pix_dist_list = []

        pix_dist_list = self.calc_source_pix_weights()

        for k in range(len(pix_dist_list)):

            source_key = self.slit_source_name[k]

            count = 0
            for j in range(len(self.SceneObj.sources_list)):
                if self.SceneObj.sources_list[j][-1] == source_key:
                    y = self.SceneObj.sources_list[j][0].spectrum
                    x = self.SceneObj.sources_list[j][0].wavelengths.value
                    count += 1
            if count == 0:
                raise ValueError("Inputted sources names must match inputted source name in the scene class")

            # ------ Extract only the wavelengths in the spectral region
            ind = np.where(  (x >= self.min_wave.value) & (x <= self.max_wave.value)  )
            xx = x[ind]
            yy = y[ind]

            # ------ Determine the fractional slit transmission
            fractional_slit_transmission = self._calc_slit_transmission()[k]

            # ------  Find detected flux
            fact = 1/( const.h*const.c ).to('erg angstrom') # Erg to photons conversion factor
            T = self._getTransmission( xx ) # Get transmission array
            y_phot = yy*(xx*fact.value)*self.TelescopeObj.mirror_area.value*T*fractional_slit_transmission*self.gain # Flux in photons / s / angstrom

            interp_source_spectrum = scipy.interpolate.interp1d( xx, y_phot,  kind='cubic' )

            # ------ Create the detector heat map
            if self.case == 1:
                angle = -(270 - (self.thetaList[k] - 45))
            else:
                angle = -(270 - self.thetaList[k])

            rot_pix_dist = rotate(pix_dist_list[k], angle = angle, reshape = True, mode = "constant") # Straighten slit
            rot_pix_dist = np.where(rot_pix_dist > 0.001, rot_pix_dist, 0) # Mask out noise 

            slit_idxs = []

            for i in range(len(rot_pix_dist[0, :])):
                slit_idx = np.where(rot_pix_dist[i, :] > 0)[0]
                for j in slit_idx:
                    slit_idxs.append(j) 

            max_slit = 0 
            min_slit = 10 

            for j in slit_idxs:
                if max_slit < j:
                    max_slit = j
                if min_slit > j:
                    min_slit = j

            max_slit_width = max_slit - min_slit - 1
            rot_pix_dist = rot_pix_dist[:, min_slit:max_slit - 1]

            self.rot_pix_dist_list.append(rot_pix_dist)

            bands = np.arange(min(xx), max(xx), self.dispersion.to('angstrom').value)
            
            tot_xpixels =  math.ceil(( (self.max_wave - self.min_wave ) / self.dispersion ).si.value) + max_slit_width 

            pix_waves = [ round((min(xx)-self.dispersion.to('angstrom').value/2) + (i*self.dispersion.to('angstrom').value ),2) for i  in range(tot_xpixels) ] 

            source_detector = np.zeros((np.shape(rot_pix_dist)[0], tot_xpixels))

            for i in range(len(bands)-1):

                x = np.linspace(bands[i], bands[i+1], 500 )
                y = interp_source_spectrum(x)

                # Flux in photons per second
                phot_in_band = np.array(scipy.integrate.cumulative_trapezoid(y, x )[-1])
                # phot_in_band is an integer

                # Create the pixel map for this band
                pix_map = rot_pix_dist * phot_in_band

                # Update the detctor heat map
                for j in range(max_slit_width): 
                    source_detector[:,i+j] = source_detector[:,i+j] + pix_map[:,j]

            # ------ Extract the spectrum from the detector heat map
            if extraction_upperlim == 'max':
                extraction_upperlim = self.slit_height_pix

            waves_CASTORSpectrum, source_CASTORSpectrum, source_extracted_numpixs =  self._extraction(source_detector, pix_waves, extraction_width, extraction_lowerlim, extraction_upperlim)

            self.waves_CASTORSpectrum_list.append(waves_CASTORSpectrum) # Append all outputs to a list
            self.source_CASTORSpectrum_list.append(source_CASTORSpectrum)
            self.source_extracted_numpixs_list.append(source_extracted_numpixs)
            self.source_detector_list.append(source_detector)

            

    def calc_background_CASTORSpectrum(self, extraction_width=1, extraction_lowerlim=1 , extraction_upperlim='max' ):
        """
        Calculates the background spectrum.  Assumes an evenly illuminated background.
        
        params
        ------
        
        extraction_width = (int) width of the extraction box along the wavelength axis in pixels.  By default this is 1.
        
        extraction_lowerlim = (int) Lower limit of the extraction box along the spatial axis in pixels.
        
        extraction_upperlim = (int) Upper limit of the extraction box along the spatial axis in pixels.

        """

        self.waves_CASTORSpectrumBackground_list = []
        self.background_CASTORSpectrum_list = []
        self.background_extracted_numpixs_list = []

        for k in range(len(self.slit_source_name)):

            y = self.BackgroundObj.earthshine_flam + self.BackgroundObj.zodi_flam  # Flux in erg / s / cm^2 / angstrom / square arcsecond
            x = self.BackgroundObj.earthshine_wavelengths # Wavelength in angstrom

            # ------ Extract only the wavelengths in the spectral region
            ind = np.where(  (x >= self.min_wave.value) & (x <= self.max_wave.value)  )
            xx = x[ind]
            yy = y[ind]

            # ------ Determine the fractional slit transmission
            fractional_slit_transmission = self._calc_slit_transmission()[k] 

            # ------  Find detected flux
            fact = 1/( const.h*const.c ).to('erg angstrom') # Erg to photons conversion factor
            slit_area = (self.slit_width*self.slit_height).value
            T = self._getTransmission(xx)  # Get transmission array
            y_phot = yy*(xx*fact.value)*self.TelescopeObj.mirror_area.value*T*fractional_slit_transmission*slit_area*self.gain # Flux in counts / s / angstrom

            interp_background_spectrum = scipy.interpolate.interp1d(xx, y_phot,  kind='cubic')

            # ------ Determine the source pixel weights (evenly illuminated slit)

            rot_pix = self.rot_pix_dist_list[k] 
            max_slit_width = np.shape(rot_pix)[1]
            rot_pix_dist = np.ones(np.shape(rot_pix)) * (1/(self.slit_height_pix * max_slit_width)) 

            # ------ Create the detector heat map
            tot_xpixels = math.ceil(( (self.max_wave - self.min_wave ) / self.dispersion ).si.value) + max_slit_width

            background_detector = np.zeros((np.shape(rot_pix_dist)[0], tot_xpixels ))

            pix_waves = [ round((min(xx)-self.dispersion.to('angstrom').value/2) + (i*self.dispersion.to('angstrom').value ),2) for i  in range(tot_xpixels) ]

            bands = np.arange(min(xx), max(xx), self.dispersion.to('angstrom').value)

            for i in range(len(bands)-1):

                x = np.linspace(bands[i], bands[i+1], 500 )
                y = interp_background_spectrum(x)

                # Flux in photons per second
                phot_in_band = np.array(scipy.integrate.cumulative_trapezoid(y, x )[-1])

                # Create the pixel map for this band
                pix_map = rot_pix_dist * phot_in_band

                # Update the detctor heat map
                for j in range(max_slit_width):
                    background_detector[:,i+j] = background_detector[:,i+j] + pix_map[:,j]

            # ------ Extract the spectrum from the detector heat map
            background_CASTORSpectrum = [] # in photons / pixel column

            # ------ Extract the spectrum from the detector heat map
            if extraction_upperlim == 'max':
                extraction_upperlim = self.slit_height_pix

            waves_CASTORSpectrumBackground, background_CASTORSpectrum, background_extracted_numpixs =  self._extraction(background_detector, pix_waves, extraction_width, extraction_lowerlim, extraction_upperlim)

        # Get background spectrum in counts / s / pixel
            background_CASTORSpectrum  = np.array(background_CASTORSpectrum) / ( extraction_width * ( extraction_upperlim - extraction_lowerlim +1 )  )

            self.waves_CASTORSpectrumBackground_list.append(waves_CASTORSpectrumBackground)
            self.background_CASTORSpectrum_list.append(background_CASTORSpectrum)
            self.background_extracted_numpixs_list.append(background_extracted_numpixs)


    
    def calc_snr_from_t(self, t, wave, nread=1):
        """
        Calculates the signal to noise ratio from the exposure time.

        params
        ------
        
        t = (float) time in seconds

        wave = (float) wavelength in angstroms

        Returns
        -------
        
        snr_list = (list of floats) signal to noise ratio
        """

        self.calc_source_CASTORSpectrum()
        self.calc_background_CASTORSpectrum()

        snr_list = []

        for i in range(len(self.slit_source_name)):

            read_npix = self.source_extracted_numpixs_list[i]

            signal = np.array(self.source_CASTORSpectrum_list[i]) 
            signal_t = signal * t 

            noise = signal_t + t * read_npix * (np.array(self.background_CASTORSpectrum_list[i]) + self.TelescopeObj.dark_current) + read_npix * self.TelescopeObj.read_noise**2 * nread

            S_N = np.divide(signal_t, np.sqrt(noise)  , where=np.sqrt(noise) > 0 )

            interp_SN = scipy.interpolate.interp1d(self.waves_CASTORSpectrum_list[i], S_N, kind='cubic')

            result = interp_SN(wave)

            if np.isnan(result):
                result = np.array(['NaN'])
            else:
                result = abs(result)

            snr_list.append(result.astype(float))

        return snr_list

        

    def calc_t_from_snr(self, snr, wave, nread=1):
        """
        Calculates the exposure time required to obtain a specified signal to noise ratio

        params
        ------

        snr = (float) signal to noise ratio

        wave = (float) wavelength in angstroms

        Returns
        -------

        t_list = (list of floats) exposure time in seconds
        """

        self.calc_source_CASTORSpectrum()
        self.calc_background_CASTORSpectrum()

        t_list = []

        for i in range(len(self.slit_source_name)):

            read_npix = self.source_extracted_numpixs_list[i]

            snr_sq = snr * snr

            signal = np.array(self.source_CASTORSpectrum_list[i])
            signal_sq =  signal * signal

            poisson_noise = signal + np.array(self.background_CASTORSpectrum_list[i])*read_npix + self.TelescopeObj.dark_current*read_npix

            numer1 = snr_sq * poisson_noise
            numer2 = snr_sq * snr_sq * poisson_noise * poisson_noise
            numer3 = 4 * snr_sq * signal_sq * (read_npix * self.TelescopeObj.read_noise**2 * nread)

            t =  np.divide( numer1 + np.sqrt(numer2 + numer3) , 2 * signal_sq  , where=signal_sq > 0 )

            interp_t = scipy.interpolate.interp1d(self.waves_CASTORSpectrum_list[i], t , kind='cubic') 

            result = interp_t(wave)

            if np.isnan(result):
                result = np.array(['NaN'])
            else:
                result = abs(result)

            t_list.append(result.astype(float))

        return t_list

            


