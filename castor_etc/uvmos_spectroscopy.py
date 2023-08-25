import astropy.units as u
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from astropy import constants as const
from scipy import integrate
from scipy.interpolate import interp1d
from scipy.interpolate import RectBivariateSpline
import scipy
from scipy.optimize import curve_fit
import math
from os.path import join

from castor_etc.background import Background
from castor_etc.photometry import Photometry
from castor_etc.sources import ExtendedSource, GalaxySource, PointSource, Source
from castor_etc.telescope import Telescope

from .filepaths import DATAPATH

def Gaussian(x, x0, sigma, a):
    return a * np.e**(- (x-x0)**2 / (2*sigma**2)  )

def Gaussian2D(x, y, sigma, a=1, x0=0, y0=0):
    term1 = ( (x-x0)**2 )/(2*sigma**2)  
    term2 = ( (y-y0)**2 )/(2*sigma**2)  
    return a * np.e**( - (term1 + term2 ) )

class UVMOS_Spectroscopy:
    """
    Spectroscopy class.
    """
    def __init__(self, TelescopeObj, SourceObj, BackgroundObj):
        #
        # Check inputs
        #
        

        if not isinstance(TelescopeObj, Telescope):
            raise TypeError("TelescopeObj must be a `castor_etc.Telescope` object")
       
        if not isinstance(SourceObj, PointSource):
            raise TypeError("SourceObj must be a `castor_etc.PointSource` object.  Other Sources are not currently supported.")
        
        if not isinstance(BackgroundObj, Background):
            raise TypeError("BackgroundObj must be a `castor_etc.Background` object")
        

        #
        # Assign attributes
        #
        self.TelescopeObj = TelescopeObj
        self.SourceObj = SourceObj
        self.BackgroundObj = BackgroundObj
        
        
        #
        #  These attributes may later be moved to the Telescope object
        #
        self.gain = 1
        self.min_wave = (150 * u.nm).to(u.AA) # Minimum wavelength 
        self.max_wave = (300 * u.nm).to(u.AA) # Maximum wavelength 
        self.slit_width =  0.214 * u.arcsec
        self.slit_height = 1 * u.arcsec
        self.FWHM = TelescopeObj.fwhm
        self.dispersion = 0.061 * u.nm
        self.pixel_size = TelescopeObj.px_scale
        
        
        #
        #
        #
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

    
    def specify_slit(self, slit_width=0.214*u.arcsec, slit_height = 1 * u.arcsec ):
        """
        Specify the size of the slit 
        """
        self.slit_width = slit_width
        self.slit_height = slit_height

        self.slit_width_pix = math.ceil( ( slit_width/self.pixel_size).value )
        self.slit_height_pix = math.ceil( ( slit_height/self.pixel_size).value )

    def show_slit(self):
        """
        Plot the slit transmission visualization 
        
        """
        fig = plt.figure()
        ax = plt.subplot(111)

        plt.title('Source veiwed through the slit')
        plt.xlabel('Arcseconds')
        plt.ylabel('Arcseconds')


        #-------- Plot entire light distribution 
        nt = 500
        x = np.linspace(-1,1 ,nt)
        y = np.linspace( -1,1,nt)
        XX, YY = np.meshgrid(x,y)
        
        sigma = self.FWHM.value / (2*np.sqrt(2*np.log(2)))
       
        ZZ = Gaussian2D(XX,YY,sigma)

        cax = ax.pcolor(x,y,ZZ,cmap='binary',shading='auto')

        #-------- Plot light that falls on the slit
        nt = 600
        x_onslit = np.linspace(-self.slit_width.value/2,self.slit_width.value/2 ,nt)
        y_onslit = np.linspace( -self.slit_height.value/2,self.slit_height.value/2,nt)
        XX_onslit, YY_onslit = np.meshgrid(x_onslit,y_onslit)
        ZZ_onslit = Gaussian2D(XX_onslit,YY_onslit,sigma)

        cax = ax.pcolor(x_onslit,y_onslit,ZZ_onslit,cmap='viridis',shading='auto')


        # -------- Plot the FWHM
        r = self.FWHM.value / 2
        theta = np.linspace(0,2*np.pi)
        plt.plot(r*np.cos(theta), r*np.sin(theta),color='k',linestyle='--',alpha=0.5,label='FWHM')

        # --------- Plot the slit
        lenn = 10
        # Plot left side of slit
        plt.plot( [-self.slit_width.value/2]*lenn, np.linspace(-self.slit_height.value/2 , self.slit_height.value/2, lenn), color='tab:red', linewidth=2, label='Slit')
        # Plot right side of slit
        plt.plot( [self.slit_width.value/2]*lenn, np.linspace(-self.slit_height.value/2 , self.slit_height.value/2, lenn), color='tab:red', linewidth=2 )
        # Plot upper side of slit
        plt.plot( np.linspace(-self.slit_width.value /2, self.slit_width.value /2,lenn),[self.slit_height.value/2]*lenn, color='tab:red', linewidth=2 )
        # Plot lower side of slit
        plt.plot( np.linspace(-self.slit_width.value /2, self.slit_width.value /2,lenn),[-self.slit_height.value/2]*lenn, color='tab:red', linewidth=2 )


        plt.xlim(-self.slit_width.value - 0.25,+self.slit_width.value + 0.25)
        plt.ylim(-self.slit_height.value - 0.25,+self.slit_height.value + 0.25)

        cbar = plt.colorbar(cax)
        cbar.set_label('Normalized Intensity')
        ax.set_aspect('equal')

        plt.legend(loc='upper left')
        
    def _calc_slit_transmission(self, print_transmission_fact=False):

        nt = 1000
        x = np.linspace(-10, 10, nt)
        y = np.linspace(-10, 10, nt)
        XX, YY = np.meshgrid(x, y)

        sigma = self.FWHM.value / (2 * np.sqrt(2 * np.log(2)))
        ZZ = Gaussian2D(XX, YY, sigma)

        interp_spline = RectBivariateSpline(x, y, ZZ, kx=5, ky=5)

        vol1 = interp_spline.integral(-10,10,-10,10)
        vol2 = interp_spline.integral( -self.slit_width.value/2, self.slit_width.value/2, -self.slit_height.value/2,
                                       self.slit_height.value/2)

        fractional_slit_transmission_2D = vol2 / vol1

        if print_transmission_fact:
            print('Fractional slit transmission (2D estimation) {}'.format(fractional_slit_transmission_2D))
            print('Fractional slit loss (2D estimation) {}'.format(1 - fractional_slit_transmission_2D))

        return fractional_slit_transmission_2D
    
    def calc_source_pix_weights(self):
        """
        Calculate the pixel weights for the source
        """

        nt = 1000
        x = np.linspace(-10, 10, nt)
        y = np.linspace(-10, 10, nt)
        XX, YY = np.meshgrid(x, y)

        sigma = self.FWHM.value / (2 * np.sqrt(2 * np.log(2)))
        ZZ = Gaussian2D(XX, YY, sigma)

        # ------- Mask out the light which is blocked by the slit
        ZZ[:, np.where(x > self.slit_width.value / 2)[0]] = 0
        ZZ[:, np.where(x < -self.slit_width.value / 2)[0]] = 0
        ZZ[np.where(y > self.slit_height.value / 2)[0], :] = 0
        ZZ[np.where(y < -self.slit_height.value / 2)[0], :] = 0

        interp_spline = RectBivariateSpline(x, y, ZZ, kx=5, ky=5)

        # -------
        grid_x = [(- self.slit_width_pix * self.pixel_size.value / 2) + self.pixel_size.value * i for i in
                  range(0, self.slit_width_pix + 1)]
        grid_y = [(- self.slit_height_pix * self.pixel_size.value / 2) + self.pixel_size.value * i for i in
                  range(0, self.slit_height_pix + 1)]

        # ------- Total transmitted intensity
        vol1 = interp_spline.integral( -(self.slit_width_pix * self.pixel_size.value)/2,
                                            (self.slit_width_pix * self.pixel_size.value)/2,
                                            -(self.slit_height_pix * self.pixel_size.value)/2,
                                            (self.slit_height_pix * self.pixel_size.value)/2)
        # -------
        vals = []
        for j in range(1,len(grid_y)):    
            for i in range(1,len(grid_x)):

                # Intensity transmitted to this "pixel"
                vol2 = interp_spline.integral( grid_x[i-1] , grid_x[i] , grid_y[j-1] , grid_y[j] )
                vals.append(vol2/vol1)

        pix_dist = np.reshape(vals,(self.slit_height_pix ,self.slit_width_pix ))
         
        return pix_dist
    
    def show_source_pix_weights(self):
        
        pix_dist = self.calc_source_pix_weights()
        
        # -----------------------------
        fig = plt.figure()
        ax = plt.subplot(111)

        plt.title('Source veiwed on the detector')
        plt.xlabel('Pixel')
        plt.ylabel('Pixel')

        im = plt.imshow(pix_dist,origin='lower',cmap='magma')

        ax.set_xticks(np.arange(-0.5,self.slit_width_pix,1))
        ax.set_yticks(np.arange(-0.5,self.slit_height_pix,1))

        ax.set_xticklabels(np.arange(0,self.slit_width_pix+1))
        ax.set_yticklabels(np.arange(0,self.slit_height_pix+1))

        plt.grid()

        cbar = plt.colorbar( )
        cbar.set_label('Normalized Intensity')
        ax.set_aspect('equal')
        
    def _extraction(self, detector, pix_waves, extraction_width, extraction_lowerlim, extraction_upperlim):
        
        # --- Warnings
        if extraction_lowerlim < 0:
            raise TypeError('Lower limit of extraction box cannot be negative.')
        
        if extraction_upperlim > self.slit_height_pix:
            raise TypeError('Upper limit of the extraction box cannot be greater than the slit height in pixels')
       
        if extraction_width > len(self.source_detector[0]):
            raise TypeError('The width of the extraction box cannot be greater than the total number of pixels along the wavelength axis')

        if not isinstance( extraction_lowerlim, int ):
            raise TypeError('The lower edge of the extraction box must be an integer')

        if not isinstance( extraction_upperlim, int ):
            raise TypeError('The upper edge of the extraction box must be an integer')

        if not isinstance( extraction_width, int ):
            raise TypeError('The width of the extraction box must be an integer')
            
            
        # ---------
        num_extractions = math.floor( len(detector[0]) / extraction_width)
            
        extracted_numpix = (extraction_width*( extraction_upperlim- extraction_lowerlim +1 ) ) * num_extractions
        
        extracted_waves = [ np.mean(pix_waves[i:i + extraction_width]) for i in range(0, len(pix_waves), extraction_width) if len(pix_waves[i:i + extraction_width])== extraction_width ]
        
        extracted_flux = [] # counts / extraction box 

        for i in range(num_extractions):

            x_left = extraction_width* i
            x_right = extraction_width* i + extraction_width

            extracted_flux.append(sum(detector[extraction_lowerlim - 1:extraction_upperlim,x_left:x_right].flatten()) )

        return extracted_waves, extracted_flux, extracted_numpix

    def _getTransmission(self, x):

        transmissionpath = join(DATAPATH, "UVMOS_data", "spectrographTransmission_CASTOR.txt")
        lines = open(transmissionpath).readlines() #  NOTE: Could be moved to 'data' folder  later
        trans_w = [float(line.split()[0]) * 10 for line in lines if '#' not in line.split()]
        trans_T = [float(line.split()[1]) for line in lines if '#' not in line.split()]

        func = interp1d( trans_w, trans_T, kind='cubic'  )

        return func(x)

    def _getDispersion(self, x):

        dispersionpath = join(DATAPATH, "UVMOS_data", "dispersion_resolution2.dat")
        lines = open(dispersionpath).readlines() #  NOTE: Could be moved to 'data' folder  later
        disp_w = [float(line.split()[0]) for line in lines if '#' not in line]
        disp_v = [float(line.split()[2]) for line in lines if '#' not in line]

        func = interp1d( disp_w, disp_v, kind='cubic', fill_value='extrapolate'  )

        return func(x)

    def showTransmission(self):
        x = np.linspace( self.min_wave.value, self.max_wave.value, 1000 )

        plt.figure(figsize=(10, 5))
        plt.plot( x, self._getTransmission(x), color='tab:red' )
        plt.title("CASTOR Spectrographic Throughput")
        plt.ylabel('Efficiency')
        plt.xlabel(r"Wavelength [$\AA$]")

    def _calc_sigmaPix(self, dispersion):
        """
        Calculate the sigma value of a Gaussian fit to a pixel with a width in wavelength equal
        to the provided dispersion
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
        cf_x0, cf_sigma, cf_a = curve_fit(Gaussian, gate_w, gate_f)[0]  # cf_x0 ,

        return abs(cf_sigma)

    def _calc_sigmaPSF(self, dispersion):
        """
        Calculate the sigma value of the Gaussian PSF on the detector
        """

        sigma_psf = self.FWHM / (2 * np.sqrt(2 * np.log(2)))
        sigma_psf_wave = sigma_psf * (1 / self.pixel_size) * dispersion

        return sigma_psf_wave

    def _calcR(self, x, disp ):
        """
        Calculate the resolving power

        Parameters
        ----------
        x array of wavelength values
        disp array of dispersion values corresponding to the wavelenghts in x

        Returns
        -------
        Array of resolving power values corresponding to x
        """

        fact = 2 * np.sqrt(2 * np.log(2))
        dL =  [ fact *  np.sqrt( self._calc_sigmaPSF(d)**2 + self._calc_sigmaPix(d)**2  )for d in disp]

        return [x[i] / dL[i] for i in range(len(x))]

    def showResolvingPower(self,):

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

    def calc_source_CASTORSpectrum(self, extraction_width=int(1), extraction_lowerlim=0 , extraction_upperlim='max' ):
        """
        
        params
        ------
        
        extraction_width = (int) width of the extraction box along the wavelength axis in pixels.  By default this is 1.
        
        extraction_lowerlim = (int) Lower limit of the extraction box along the spatial axis in pixels.
        
        extraction_upperlim = (int) Upper limit of the extraction box along the spatial axis in pixels.
        """
        
        y = self.SourceObj.spectrum  # Flux in erg / s / cm^2 / angstrom
        x = self.SourceObj.wavelengths.to('angstrom').value # Wavelength in angstrom

        # ------ Extract only the wavelengths in the spectral region 
        ind = np.where(  (x >= self.min_wave.value) & (x <= self.max_wave.value)  )
        xx = x[ind]
        yy = y[ind]
        
        # ------ Determine the fractional slit transmission
        fractional_slit_transmission = self._calc_slit_transmission()
     
        # ------  Find detected flux 
        fact = 1/( const.h*const.c ).to('erg angstrom') # Erg to photons conversion factor
        T = self._getTransmission( xx ) # Get transmission array
        y_phot = yy*(xx*fact.value)*self.TelescopeObj.mirror_area.value*T*fractional_slit_transmission*self.gain # Flux in photons / s / angstrom

        interp_source_spectrum = scipy.interpolate.interp1d( xx, y_phot,  kind='cubic'  )
        
        # ------ Determine the source pixel weights
        pix_dist = self.calc_source_pix_weights()
        
        # ------ Create the detector heat map       
        tot_xpixels = math.ceil(( (self.max_wave - self.min_wave ) / self.dispersion ).si.value) + self.slit_width_pix

        self.source_detector = np.zeros((self.slit_height_pix,tot_xpixels ))

        pix_waves = [ round((min(xx)-self.dispersion.to('angstrom').value/2) + (i*self.dispersion.to('angstrom').value ),2) for i  in range(tot_xpixels) ]
       
        
        bands = np.arange(min(xx), max(xx), self.dispersion.to('angstrom').value)
        
        for i in range(0,len(bands)-1):

            x = np.linspace(bands[i], bands[i+1], 500 )
            y = interp_source_spectrum(x)

            # Flux in photons per second
            phot_in_band = np.array(scipy.integrate.cumtrapz(y, x )[-1])

            # Create the pixel map for this band
            pix_map = pix_dist * phot_in_band


            # Update the detctor heat map
            for j in range(self.slit_width_pix):
                self.source_detector[:,i+j] = self.source_detector[:,i+j] + pix_map[:,j]   

        # ------ Extract the spectrum from the detector heat map
        if extraction_upperlim == 'max':
            extraction_upperlim = self.slit_height_pix
        self.waves_CASTORSpectrum, self.source_CASTORSpectrum, self.source_extracted_numpixs =  self._extraction(self.source_detector, pix_waves, extraction_width, extraction_lowerlim, extraction_upperlim)

    def calc_background_CASTORSpectrum(self, extraction_width=int(1), extraction_lowerlim=0 , extraction_upperlim='max' ):
        """
        
        params
        ------
        
        extraction_width = (int) width of the extraction box along the wavelength axis in pixels.  By default this is 1.
        
        extraction_lowerlim = (int) Lower limit of the extraction box along the spatial axis in pixels.
        
        extraction_upperlim = (int) Upper limit of the extraction box along the spatial axis in pixels.
        """
        
        y = self.BackgroundObj.earthshine_flam + self.BackgroundObj.zodi_flam  # Flux in erg / s / cm^2 / angstrom / square arcsecond
        x = self.BackgroundObj.earthshine_wavelengths # Wavelength in angstrom

        # ------ Extract only the wavelengths in the spectral region 
        ind = np.where(  (x >= self.min_wave.value) & (x <= self.max_wave.value)  )
        xx = x[ind]
        yy = y[ind]

        # ------ Determine the fractional slit transmission
        fractional_slit_transmission = self._calc_slit_transmission()

        # ------  Find detected flux 
        fact = 1/( const.h*const.c ).to('erg angstrom') # Erg to photons conversion factor
        slit_area = (self.slit_width*self.slit_height).value
        T = self._getTransmission(xx)  # Get transmission array
        y_phot = yy*(xx*fact.value)*self.TelescopeObj.mirror_area.value*T*fractional_slit_transmission*slit_area*self.gain # Flux in counts / s / angstrom
        
        interp_background_spectrum = scipy.interpolate.interp1d( xx, y_phot,  kind='cubic'  )
        
        # ------ Determine the source pixel weights (evenly illuminated slit)
        pix_dist = pix_dist = np.full( (self.slit_height_pix ,self.slit_width_pix), 1/(self.slit_width_pix*self.slit_height_pix) )
        
        # ------ Create the detector heat map   
        tot_xpixels = math.ceil(( (self.max_wave - self.min_wave ) / self.dispersion ).si.value) + self.slit_width_pix
        
        self.background_detector = np.zeros((self.slit_height_pix,tot_xpixels ))
        
        pix_waves = [ round((min(xx)-self.dispersion.to('angstrom').value/2) + (i*self.dispersion.to('angstrom').value ),2) for i  in range(tot_xpixels) ]

        bands = np.arange(min(xx), max(xx), self.dispersion.to('angstrom').value)
        
        for i in range(0,len(bands)-1):

            x = np.linspace(bands[i], bands[i+1], 500 )
            y = interp_background_spectrum(x)

            # Flux in photons per second
            phot_in_band = np.array(scipy.integrate.cumtrapz(y, x )[-1])

            # Create the pixel map for this band
            pix_map = pix_dist * phot_in_band


            # Update the detctor heat map
            for j in range(self.slit_width_pix):
                self.background_detector[:,i+j] = self.background_detector[:,i+j] + pix_map[:,j]   

        # ------ Extract the spectrum from the detector heat map     
        self.background_CASTORSpectrum = [] # in photons / pixel column
        
        # ------ Extract the spectrum from the detector heat map
        if extraction_upperlim == 'max':
            extraction_upperlim = self.slit_height_pix
        
        self.waves_CASTORSpectrumBackground, self.background_CASTORSpectrum, self.background_extracted_numpixs =  self._extraction(self.background_detector, pix_waves, extraction_width, extraction_lowerlim, extraction_upperlim)
        
        # Get background spectrum in counts / s / pixel
        self.background_CASTORSpectrum  = np.array(self.background_CASTORSpectrum) / ( extraction_width * ( extraction_upperlim - extraction_lowerlim +1 )  )

    def calc_snr_from_t(self, t, wave, nread=1):   
        
        read_npix = self.source_extracted_numpixs 
        
        signal_t = np.array(self.source_CASTORSpectrum)*t
        
        noise = signal_t + t*read_npix*(np.array(self.background_CASTORSpectrum)+self.TelescopeObj.dark_current) + read_npix*self.TelescopeObj.read_noise**2*nread 
        
        S_N = np.divide(signal_t, np.sqrt(noise)  , where=np.sqrt(noise) > 0 )
        
        interp_SN = scipy.interpolate.interp1d(self.waves_CASTORSpectrum, S_N, kind='cubic')
        
        return interp_SN(wave)
    
    def calc_t_from_snr(self, snr, wave, nread=1):
        
        read_npix = self.source_extracted_numpixs 
        
        snr_sq = snr * snr

        signal = np.array(self.source_CASTORSpectrum)
        signal_sq =  signal* signal

        poisson_noise = signal + np.array(self.background_CASTORSpectrum)*read_npix + self.TelescopeObj.dark_current*read_npix

        numer1 = snr_sq * poisson_noise
        numer2 = snr_sq * snr_sq * poisson_noise * poisson_noise
        numer3 = 4 * snr_sq * signal_sq * (read_npix * self.TelescopeObj.read_noise**2 * nread)


        t =  np.divide( numer1 + np.sqrt(numer2 + numer3) , 2 * signal_sq  , where=signal_sq > 0 )


        interp_t = scipy.interpolate.interp1d(self.waves_CASTORSpectrum, t , kind='cubic')
    
        return interp_t(wave)
