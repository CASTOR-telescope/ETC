"""
Generate and handle spectral data and normalizations.

Includes:
  - redshifting wavelengths
  - blackbody radiation, power-law spectrum generation
  - user-input spectrum
  - Gaussian and Lorentzian emission/absorption lines
  - generic spiral and elliptical galaxy spectra
  - stellar spectra from the Pickles catalog
  - stellar spectra from the 'ATLAS9' and 'BTSettl' catalogs
  - normalization functions:
    - normalize a blackbody spectrum to a star of given radius and distance
    - normalize a spectrum to some average value or AB magnitude, either within a passband
      or over the whole spectrum
    - normalize a spectrum to a given bolometric luminosity and distance
    - calculate the average value of a spectrum (erg/s/cm^2/A or AB mag) either within a
      passband or over the whole spectrum

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
from numbers import Number
from os.path import join
import os

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.io import fits
from scipy.integrate import simpson
from scipy.interpolate import interp1d

from . import constants as const
from .conversions import calc_photon_energy, flam_to_AB_mag, fnu_to_flam, mag_to_flux
from .filepaths import DATAPATH
from .telescope import Telescope

from astroquery.gaia import Gaia
from astropy.coordinates import SkyCoord
from astropy.constants import R_sun, pc
from astropy.wcs import WCS

# Change the path variable according to the build instructions.
path = '/Users/dhananjhaybansal/Applications/POET/stellar_models/'
#path = '/opt/mypackage/castor_etc/data/transit_data/stellar_models/'

def getStarData(temperature, metallicity, logg, 
                    model_grid='ATLAS9',
                    ):
    """
    Reads in star data from atlas directory, according to temperature, metallicity, and log(G) values.
    
    Parameters
    ---------- 
        temperature :: float
            Temperature of the reference star in kelvin.

        metallicity :: 'p/m00','p/m05','p/m10','p/m15','p/m20', where p is plus, m is minus.
            Metallicity of the reference star. Accepted inputs. 

        logg :: g00,g05,g10,g15,g20,g25,g30,g35,g40,g45,g50
            log(G) value for the reference star.
    
    Returns
    -------- 
        wavelengths :: `astropy.Quantity` wavelength  
            wavelength of the reference star in AA. 

        flux :: array of floats
         Flux of the reference star in units erg s^-1 cm^-2 A^-1
    """

    # Default to ATLAS9 model
    #     If Teff < 3500 K (i.e., min Teff in ATLAS9 grid, switch to BTSettl grid)
    if ( (model_grid == 'BTSettl') | (temperature < 3500) ):
        if not os.path.exists(os.path.join(path,"BTSettl_CIFIST/")):
            raise RuntimeError(
                f"The specified directory path `{path}` does not contain the stellar_models directory." + "Please check the path variable."
            )
        else:
          grid_dir = os.path.join(path,"BTSettl_CIFIST/")

        # Scan directory to inventory model grid
        # Assumes all files starting with 'lte' are model files
        # Only stored solar metallicity models; may need to update for other metallicities
        dir_flist = os.listdir(grid_dir)
        model_flist = []
        grid_param = []
        for n,fname in enumerate(dir_flist):
            if fname.startswith('lte'):
                _teff = float(fname[3:6]) * 100
                if fname[6] == '+':
                    _logg = -1 * float(fname[7:10])
                else:
                    _logg = float(fname[7:10])
                _mh = float(fname[11:14])
                model_flist.append(fname)
                grid_param.append([_teff, _logg, _mh])
        grid_param = np.vstack(grid_param)
        model_flist = np.array(model_flist)

        # Identify closest model prioritizing teff, logg, then metallicity
        mi = np.argmin(np.abs(grid_param[:,0] - temperature))
        ti = np.where(grid_param[:,0] == grid_param[mi,0])[0]
        grid_param = grid_param[ti]
        model_flist = model_flist[ti]

        mi = np.argmin(np.abs(grid_param[:,1] - logg))
        ti = np.where(grid_param[:,1] == grid_param[mi,1])[0]
        grid_param = grid_param[ti]
        model_flist = model_flist[ti]

        mi = np.argmin(np.abs(grid_param[:,2] - metallicity))
        ti = np.where(grid_param[:,2] == grid_param[mi,2])[0]
        grid_param = grid_param[ti]
        model_flist = model_flist[ti]

        # Read in selected model
        d = np.loadtxt(grid_dir + model_flist[0])
        wv = d[:,0] * u.AA
        flux = d[:,1] # [erg/s/cm2/angstrom]

    elif model_grid == 'ATLAS9':
        if not os.path.exists(os.path.join(path, "ATLAS9/ck04models/")):
            raise RuntimeError(
                f"The specified directory path `{path}` does not contain the stellar_model directory." + "Please modify the path variable. "
            )
        else:
            grid_dir = os.path.join(path, "ATLAS9/ck04models/")

        teff_grid = np.array([3500, 3750, 4000, 4250, 4500, 4750, 5000, 5250, 5500, 5750, 
                              6000, 6250, 6500, 6750, 7000, 7250, 7500, 7750, 8000, 8250, 
                              8500, 8750, 9000, 9250, 9500, 9750, 10000, 10250, 10500, 10750, 
                              11000, 11250, 11500, 11750, 12000, 12250, 12500, 12750, 13000, 14000, 
                              15000, 16000, 17000, 18000, 19000, 20000, 21000, 22000, 23000, 24000, 
                              25000, 26000, 27000, 28000, 29000, 30000, 31000, 32000, 33000, 34000, 
                              35000, 36000, 37000, 38000, 39000, 40000, 41000, 42000, 43000, 44000, 
                              45000, 46000, 47000, 48000, 49000, 50000])
        mh_grid = np.array([-0.25, -0.20, -0.15, -0.10, -0.05, 0.0, 0.02, 0.05])
        logg_grid = np.array([0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0])

        mi = [ np.argmin( np.abs(teff_grid - temperature) ),
               np.argmin( np.abs(mh_grid - metallicity) ),
               np.argmin( np.abs(logg_grid - logg) ) ]

        teff_str = '{:.0f}'.format( teff_grid[mi[0]] )
        if mh_grid[mi[1]] < 0:
            mh_str = 'm{:0>2d}'.format( int(np.abs( mh_grid[mi[1]] ) * 100) )
        else:
            mh_str = 'p{:0>2d}'.format( int(np.abs( mh_grid[mi[1]] ) * 100) )
        logg_str = 'g{:0>2d}'.format( int( logg_grid[mi[2]] * 10) )

        specfile = grid_dir + 'ck' + mh_str + '/' + 'ck' + mh_str + '_' + teff_str + '.fits'
        if os.path.isfile(specfile):
            hdul = fits.open(specfile)
            data = hdul[1].data #the first extension has table
            wv = data['WAVELENGTH'] * u.AA
            flux = data[logg_str] # [erg/s/cm^2/A]
        else:
            print('Spectrum not found: ',specfile)

    return wv, flux

def interp_EEM_table(Teff=[], Gmag=[], Bpmag=[], Rpmag=[]):
    """
    Interpolate EEM_table using input Teff

    Parameters
      ----------
          Teff :: array of floats
            Effective temperature of the identified sources

          Gmag :: array of floats
            Gaia G magnitude of the identified sources

          Bpmag :: array of floats
            Gaia BP magnitude of the identifies sources

          Rpmag :: array of floats
            Gaia RP magnitude of the idenitifies sources
      
      Return
      ------
          interp_output
            Interpolation result
    
    """

    table_fname = join(DATAPATH,"transit_data","EEM_dwarf_UBVIJHK_colors_Teff.txt")

    _logg_sun = 4.438

    # Read in EEM table
    eem = {'Teff':[], 'radius':[], 'mass':[], 'logg':[], 'Gmag_abs':[],
           'B-V':[], 'G-V':[], 'Bp-Rp':[], 'G-Rp':[], 'U-B':[], 'V-I':[]}
    eem_hdr = []
    with open(table_fname,'r') as file:
        data_block = False
        for line in file:
            # Only read in between lines starting with '$SpT'
            if line.startswith('#SpT'):
                if (data_block == False):
                    data_block = True
                else:
                    break
            if data_block:
                _line = line.rstrip().split()

                if len(eem_hdr) == 0:
                    eem_hdr = np.array(_line.copy())
                else:
                    ci = np.where(eem_hdr == 'Teff')[0][0]
                    eem['Teff'].append( float(_line[ci]) )

                    ci = np.where(eem_hdr == 'R_Rsun')[0][0]
                    if _line[ci].startswith('...'):
                        eem['radius'].append( np.nan )
                    else:
                        eem['radius'].append( float(_line[ci]) )

                    ci = np.where(eem_hdr == 'Msun')[0][0]
                    if _line[ci].startswith('...'):
                        eem['mass'].append( np.nan )
                    else:
                        eem['mass'].append( float(_line[ci]) )

                    ci = np.where(eem_hdr == 'M_G')[0][0]
                    if _line[ci].startswith('...'):
                        eem['Gmag_abs'].append( np.nan )
                    elif _line[ci].endswith(':'):
                        eem['Gmag_abs'].append( float(_line[ci][:-1]) )
                    else:
                        eem['Gmag_abs'].append( float(_line[ci]) )

                    ci = np.where(eem_hdr == 'B-V')[0][0]
                    if _line[ci].startswith('...'):
                        eem['B-V'].append( np.nan )
                    else:
                        eem['B-V'].append( float(_line[ci]) )

                    ci = np.where(eem_hdr == 'G-V')[0][0]
                    if _line[ci].startswith('...'):
                        eem['G-V'].append( np.nan )
                    else:
                        eem['G-V'].append( float(_line[ci]) )

                    ci = np.where(eem_hdr == 'Bp-Rp')[0][0]
                    if _line[ci].startswith('...'):
                        eem['Bp-Rp'].append( np.nan )
                    else:
                        eem['Bp-Rp'].append( float(_line[ci]) )

                    ci = np.where(eem_hdr == 'G-Rp')[0][0]
                    if _line[ci].startswith('...'):
                        eem['G-Rp'].append( np.nan )
                    else:
                        eem['G-Rp'].append( float(_line[ci]) )

                    ci = np.where(eem_hdr == 'U-B')[0][0]
                    if _line[ci].startswith('...'):
                        eem['U-B'].append( np.nan )
                    else:
                        eem['U-B'].append( float(_line[ci]) )

                    ci = np.where(eem_hdr == 'V-Ic')[0][0]
                    if _line[ci].startswith('...'):
                        eem['V-I'].append( np.nan )
                    else:
                        eem['V-I'].append( float(_line[ci]) )
    for _lbl in eem.keys():
        eem[_lbl] = np.array(eem[_lbl])

    # Calculate logg
    eem['logg'] = _logg_sun + np.log10( eem['mass'] )  - 2. * np.log10( eem['radius']**2 )

    # Sort by Teff
    si = np.argsort(eem['Teff'])
    for _lbl in eem.keys():
        eem[_lbl] = eem[_lbl][si]

    # If Teff not specified, use Bpmag, Rpmag to get Teff
    if len(Teff) == 0:
        si = np.argsort(eem['Bp-Rp'])
        Teff = np.interp(Bpmag - Rpmag, eem['Bp-Rp'][si], eem['Teff'][si])

    # Interpolate table using input Teff
    # Use grid edges for missing values
    interp_output = {'Teff':Teff}
    for _lbl in eem.keys():
        if _lbl != 'Teff':
            vi = np.isfinite(eem[_lbl])
            interp_output[_lbl] = np.interp(np.log10(Teff), 
                                        np.log10(eem['Teff'][vi]), eem[_lbl][vi])

    # Calculate U, I and add to table (if Gmag provided)
    if len(Gmag) > 0:
        interp_output['Rpmag'] = Gmag - interp_output['G-Rp']
        interp_output['Bpmag'] = interp_output['Bp-Rp'] + interp_output['Rpmag']

        interp_output['Umag'] = interp_output['U-B'] + interp_output['B-V'] - \
                                    interp_output['G-V'] + Gmag
        interp_output['Imag'] = Gmag - ( interp_output['G-V'] + interp_output['V-I'] )

    return interp_output


def redshift_wavelengths(wavelengths, redshift):
    """
    Apply redshift correction to wavelengths.

    Parameters
    ----------
      wavelengths :: int or float or `astropy.Quantity` or array
        Wavelengths to redshift.

      redshift :: int or float
        Redshift to apply.

    Returns
    -------
        red_wavelengths :: array of floats
        Redshifted wavelengths.
    """
    return wavelengths * (1 + redshift)


class SpectrumMixin:
    """
    Mixin for generating spectra. To be used with `Source` object. Do not use directly!
    """

    def _check_existing_spectrum(self, overwrite, quiet=False):
        """
        Check for existing spectrum and if the source is a `CustomSource` instance. If the
        object is a `CustomSource` instance, raise an error. If self.wavelengths and/or
        self.spectrum is not None, raise an error if overwrite is False; otherwise print a
        message notifying user that the method will overwrite the existing spectrum
        (unless `quiet` is True).

        Parameters
        ----------
          overwrite :: bool
            If False and self.wavelengths and/or self.spectrum is not None, raise an
            error. If True and self.wavelengths and/or self.spectrum is not None, print a
            message informing the user that the existing spectrum will be overwritten
            (unless `quiet` is True).

          quiet :: bool
            If True, do not print a message when overwriting an existing spectrum.

        Returns
        -------
          None
        """
        from .sources import CustomSource  # avoid circular import error

        if isinstance(self, CustomSource):
            raise ValueError("A `CustomSource` object does not support a spectrum.")
        if self.wavelengths is not None or self.spectrum is not None:
            if not overwrite:
                raise ValueError(
                    "wavelengths/spectrum already exists! "
                    + "Use overwrite=True to overwrite wavelengths/spectrum."
                )
            elif not quiet:
                print(
                    "INFO: Overwriting existing wavelengths/spectrum "
                    + "with new wavelengths/spectrum."
                )

    def spectrum_erg_to_photon(self):
        """
        Convert the spectrum that has units of ergs in the numerator to units of photons
        in the numerator.

        Attributes
        ----------
          spectrum :: array of floats
            Spectrum that has units of photons in the numerator.

        Returns
        -------
          None
        """
        if self.spectrum is None or self.wavelengths is None:
            raise ValueError(
                "Please generate a spectrum before converting ergs to photons."
            )
        self.spectrum /= calc_photon_energy(wavelength=self.wavelengths)[0]

    def redshift_wavelengths(self, redshift):
        """
        Apply redshift correction to wavelengths. Does not affect the y-axis values of the
        spectrum.

        Parameters
        ----------
          redshift :: int or float
            Redshift to apply.

        Attributes
        ----------
          wavelengths :: `astropy.Quantity` array
            The redshifted wavelengths of the spectrum, in angstroms.

        Returns
        -------
          None
        """
        if not isinstance(redshift, Number):
            raise TypeError("redshift must be an int or float")
        self.wavelengths *= 1 + redshift

    def generate_uniform(
        self, wavelengths, value, unit="ABmag", overwrite=False, quiet=False
    ):
        """
        Generate a uniform spectrum equal to a constant value in either flam
        (erg/s/cm^2/A), fnu (erg/s/cm^2/Hz), ABmag (AB magnitude), or STmag (ST
        magnitude). Note that the computed (and stored) spectrum will always be in units
        of flam.

        Parameters
        ----------
          wavelengths :: array of scalars or `astropy.Quantity` array
            The wavelengths over which to generate the uniform spectrum. If an array of
            scalars, it should be in angstrom.

          value :: int or float
            The value of the uniform spectrum in the specified `unit`.

          unit :: "flam" or "fnu" or "ABmag" or "STmag"
            The unit of the `value`: either flam (erg/s/cm^2/A), fnu (erg/s/cm^2/Hz),
            ABmag (AB magnitude), or STmag (ST magnitude).

          overwrite :: bool
            If True, overwrite any existing wavelengths/spectrum. If False, raise an error
            if wavelengths or spectrum is not None.

          quiet :: bool
            If True, do not print a message when overwriting an existing spectrum.

        Attributes
        ----------
          wavelengths :: `astropy.Quantity` array
            The given wavelengths of the spectrum, converted into units of angstrom.

          spectrum :: array of floats
            Spectrum, in erg/s/cm^2/A, that is uniform in the specified unit.

        Returns
        -------
          None
        """
        #
        # Check inputs
        #
        self._check_existing_spectrum(overwrite, quiet=quiet)
        if not isinstance(wavelengths, u.Quantity):
            wavelengths = wavelengths * u.AA
        else:
            wavelengths = wavelengths.to(u.AA)
        if not isinstance(value, Number):
            raise ValueError("`value` must be an int or float")
        if unit not in ["flam", "fnu", "ABmag", "STmag"]:
            raise ValueError("`unit` must be one of 'flam', 'fnu', 'ABmag', or 'STmag'")
        if unit == "flam" or unit == "fnu":
            if value <= 0:
                raise ValueError("`value` must be > 0 if `unit` is 'flam' or 'fnu'")
        #
        # Generate spectrum
        #
        spectrum = np.full(np.shape(wavelengths), value, dtype=float)
        if unit == "fnu":
            spectrum = fnu_to_flam(
                fnu=spectrum, wavelength=wavelengths, fnu_err=0.0, wavelength_err=0.0
            )[0]
        elif unit == "ABmag":
            # Convert to fnu
            spectrum = mag_to_flux(mag=spectrum, mag_err=0.0, zpt=-48.60)[0]
            # Convert fnu to flam
            spectrum = fnu_to_flam(
                fnu=spectrum, wavelength=wavelengths, fnu_err=0.0, wavelength_err=0.0
            )[0]
        elif unit == "STmag":
            # Convert directly to flam
            spectrum = mag_to_flux(mag=spectrum, mag_err=0.0, zpt=-21.10)[0]
        self.wavelengths = wavelengths
        self.spectrum = spectrum

    def generate_bb(
        self,
        T,
        redshift=0.0,
        emissivity=1.0,
        wavelengths=None,
        limits=[0.09, 1.2] << u.um,
        resolution=1 << u.nm,
        radius=1,
        dist=1 << u.kpc,
        overwrite=False,
        quiet=False,
    ):
        """
        Generate a blackbody (BB) spectrum (in erg/s/cm^2/A) using Planck's radiation law.
        The spectral radiance of the BB (erg/s/cm^2/A/sr) is normalized to a star of given
        radius and distance.

        Parameters
        ----------
          T :: int or float or `astropy.Quantity`
            Intrinsic blackbody temperature (i.e., the temperature of the BB at
            redshift=0). If int or float, the unit is assumed to be kelvin.

          redshift :: int or float
            Redshift of the blackbody.

          emissivity :: int or float
            Emissivity of the blackbody. (Technically, emissivity is unity per the
            definition of a BB).

          wavelengths :: array of floats or `astropy.Quantity` array
            The wavelengths over which to calculate the spectrum. If an array of floats,
            the unit is assumed to be in angstroms. If wavelengths is not None, the limits
            and resolution parameters are ignored. Note that the final wavelengths
            attribute will be an `astropy.Quantity` array in units of angstroms regardless
            of this input's units.

          limits :: list of 2 scalars or list of 2 `astropy.Quantity`
            List containing the lower (0th index) and upper (1st index) bounds for the BB
            spectrum's restframe wavelengths, inclusive. Limits should be > 0. If list
            elements are int or float, they are assumed to be in angstroms. This parameter
            is ignored if wavelengths is provided.

          resolution :: int or float or `astropy.Quantity`
            The wavelength resolution of the returned spectrum. If a scalar, it is assumed
            to be in units of angstroms. This parameter is ignored if wavelengths is
            provided.

          radius :: float or `astropy.Quantity`
            The radius of the source. If a scalar, it is assumed to be in units of solar
            radii.

          dist :: float or `astropy.Quantity`
            The distance to the blackbody. If a scalar, it is assumed to be in units of
            kpc.

          overwrite :: bool
            If True, overwrite any existing wavelengths/spectrum. If False, raise an error
            if wavelengths or spectrum is not None.

          quiet :: bool
            If True, do not print a message when overwriting an existing spectrum.

        Attributes
        ----------
          wavelengths :: `astropy.Quantity` array
            The redshifted wavelengths of the spectrum, in angstroms.

          spectrum :: array of floats
            BB spectrum in units of flam (erg/s/cm^2/A).

        Returns
        -------
          None
        """
        #
        # Check inputs
        #
        self._check_existing_spectrum(overwrite, quiet=quiet)
        if isinstance(T, u.Quantity):
            T = T.to(u.K, equivalencies=u.temperature()).value
        if wavelengths is None:
            limits = list(limits)
            for i, lim in enumerate(limits):
                if isinstance(lim, u.Quantity):
                    limits[i] = lim.to(u.AA).value
            if isinstance(resolution, u.Quantity):
                resolution = resolution.to(u.AA).value
            wavelengths = np.arange(limits[0], limits[1] + 0.5 * resolution, resolution)
        elif isinstance(wavelengths, u.Quantity):
            wavelengths = wavelengths.to(u.AA).value
        #
        # Generate BB spectrum with redshift
        #
        # Convert wavelengths from angstrom to cm
        wavelengths = wavelengths * 1e-8  # cm
        # Planck's radiation law
        lightspeed = const.LIGHTSPEED.value  # cm/s
        prefactor = (2 * const.PLANCK_H.value * lightspeed * lightspeed) / (
            wavelengths ** 5
        )
        denom = np.expm1(
            (const.PLANCK_H.value * lightspeed) / (wavelengths * const.K_B.value * T)
        )
        spectrum = prefactor / denom  # erg/s/cm^2/cm/sr
        #
        # Incorporate emissivity and convert per cm to per angstrom
        #
        spectrum *= 1e-8 * emissivity  # erg/s/cm^2/A/sr
        #
        # Factor in redshift and convert wavelengths back to angstroms
        #
        wavelengths = redshift_wavelengths(wavelengths, redshift) * 1e8  # angstrom
        #
        # Assign to `Source` object attributes. Spectrum is in erg/s/cm^2/A
        #
        self.wavelengths = wavelengths * u.AA
        self.spectrum = NormMixin.norm_to_star(spectrum, radius=radius, dist=dist)  # flam

    def generate_power_law(
        self, ref_wavelength, wavelengths, exponent, overwrite=False, quiet=False
    ):
        """
        Generate a spectrum with a shape following a power-law in some arbitrary unit. The
        flux is defined so that it is equal to 1 at the reference wavelength.

        The spectrum is calculated using the following formula:
        ```math
                    spectrum = (wavelengths / ref_wavelength)^exponent
        ```
        where each variable is as defined in the Parameters documentation below.

        Parameters
        ----------
          ref_wavelength :: scalar or `astropy.Quantity`
            The reference wavelength for the power-law. The spectrum at this wavelength
            will have a flux of 1. If a scalar, it should be in angstrom.

          wavelengths :: array of scalars or `astropy.Quantity` array
            The wavelengths over which to calculate the power-law spectrum. If an array of
            scalars, it should be in angstrom.

          exponent :: int or float
            The exponent for the power-law.

          overwrite :: bool
            If True, overwrite any existing wavelengths/spectrum. If False, raise an error
            if wavelengths or spectrum is not None.

          quiet :: bool
            If True, do not print a message when overwriting an existing spectrum.

        Attributes
        ----------
          wavelengths :: `astropy.Quantity` array
            The wavelengths of the spectrum, in angstroms.

          spectrum :: array of floats
            Spectrum following the power-law defined above, in arbitary units.

        Returns
        -------
          None
        """
        #
        # Check inputs
        #
        self._check_existing_spectrum(overwrite, quiet=quiet)
        if np.size(ref_wavelength) != 1:
            raise ValueError(
                "ref_wavelength must be a single scalar or `astropy.Quantity`."
            )
        if isinstance(wavelengths, u.Quantity):
            wavelengths = wavelengths.to(u.AA).value
        if isinstance(ref_wavelength, u.Quantity):
            ref_wavelength = ref_wavelength.to(u.AA).value
        #
        # Power-law
        #
        spectrum = (wavelengths / ref_wavelength) ** exponent
        self.wavelengths = wavelengths * u.AA
        self.spectrum = spectrum

    @staticmethod
    def _generate_gaussian(
        wavelengths,
        spectrum,
        center,
        fwhm,
        peak=None,
        tot_flux=None,
        add=True,
        abs_peak=True,
    ):
        """
        Add/subtract a Gaussian spectrum to/from an existing spectrum. This is useful for
        representing emission lines (i.e., by adding a Gaussian source) or absorption
        lines (i.e., by subtracting a Gaussian source). Note that the minimum/maximum
        wavelengths of the source spectrum will not change.

        The Gaussian spectrum can be represented by the following formulae (from
        <https://pysynphot.readthedocs.io/en/latest/spectrum.html#gaussian-emission>):
        ```math
                    gaussian = peak / exp[(wavelengths - center)^2 / (2 * sigma^2)]
        ```
        and
        ```math
                    sigma = fwhm / [2 * sqrt(2 * ln2)]
        ```
        and
        ```math
                    peak = tot_flux / sqrt(2 * pi * sigma^2) <-- see Gaussian integral
        ```
        where:
          - gaussian is the Gaussian spectrum's flux in some arbitrary unit
          - peak is the flux at the center of the Gaussian (i.e., the central wavelength)
          - center is the central wavelength of the Gaussian
          - wavelengths is the array of wavelengths over which to calculate the spectrum
          - fwhm is the full-width at half-maximum of the Gaussian
          - tot_flux is the total flux of the Gaussian under the curve

        Parameters
        ----------
          wavelengths :: array of floats or `astropy.Quantity` array
            The wavelengths over which to calculate the Gaussian spectrum.

          spectrum :: array of floats
            The spectrum to/from which to add/subtract the Gaussian spectrum.

          center :: scalar or `astropy.Quantity`
            The central wavelength of the Gaussian. If a scalar, it is assumed to be in
            angstrom.

          fwhm :: scalar or `astropy.Quantity`
            The full-width at half-maximum of the Gaussian. If a scalar, it is assumed to
            be in angstrom.

          peak :: int or float
            The peak flux of the Gaussian (i.e., the flux at the center wavelength).
            Exactly one of peak or tot_flux must be specified.

          tot_flux :: int or float
            The total flux under the curve. Exactly one of peak or tot_flux must be
            specified.

          add :: bool
            If True, add the Gaussian spectrum to the existing spectrum. If False,
            subtract the Gaussian from the existing spectrum.

          abs_peak :: bool
            If True, ensure that the peak of the emission line or dip of the absorption
            line is at the given value. Otherwise, just add/subtract the given Gaussian
            peak to/from the continuum.

        Returns
        -------
          sorted_wavelengths :: array of floats or `astropy.Quantity` array
            The wavelengths of the new spectrum. The shape of this array will be different
            from the input `wavelengths` array.

          sorted_spectrum :: array of floats
            The spectrum with the Gaussian added/subtracted. The shape of this array will
            be different from the input `spectrum` array.
        """

        def _gaussian(_peak, _wavelengths, _center, _sigma):
            _num_sigma = (_wavelengths - _center) / _sigma
            return _peak / np.exp(0.5 * _num_sigma * _num_sigma)

        #
        # Check inputs
        #
        if np.size(center) != 1:
            raise ValueError("center must be a single scalar or `astropy.Quantity`.")
        if np.size(fwhm) != 1:
            raise ValueError("fwhm must be a single scalar or `astropy.Quantity`.")
        if (peak is None and tot_flux is None) or (
            peak is not None and tot_flux is not None
        ):
            raise ValueError("Exactly one of peak or tot_flux must be specified.")
        if peak is not None and (
            np.size(peak) != 1 or not isinstance(peak, Number) or peak <= 0
        ):
            if add:
                raise ValueError("peak must be a single int or float >= 0.")
            else:
                raise ValueError("dip must be a single int or float >= 0.")
        elif tot_flux is not None and (
            np.size(tot_flux) != 1 or not isinstance(tot_flux, Number) or tot_flux <= 0
        ):
            raise ValueError("tot_flux must be a single int or float >= 0.")
        # Convert lengths to angstrom
        if isinstance(wavelengths, u.Quantity):
            wavelengths_unit = wavelengths.unit
            wavelengths = wavelengths.to(u.AA).value
        else:
            wavelengths_unit = None
        if isinstance(center, u.Quantity):
            center = center.to(u.AA).value
        if isinstance(fwhm, u.Quantity):
            fwhm = fwhm.to(u.AA).value
        if fwhm <= 0:
            raise ValueError("fwhm must be >= 0.")
        #
        # Gaussian spectrum
        #
        sigma = fwhm / (2 * np.sqrt(2 * np.log(2)))
        if peak is None:
            peak = tot_flux / (np.sqrt(2 * np.pi) * sigma)
        spectrum_interp = interp1d(
            wavelengths, spectrum, kind="linear", bounds_error=False, fill_value=np.nan
        )
        if add:
            if abs_peak:
                # Ensure final peak is actually at desired value
                peak -= spectrum_interp(center)
                if peak < 0:
                    raise ValueError("peak of emission line below continuum.")
        else:
            if abs_peak:
                # Ensure final dip is actually at desired value
                center_val = spectrum_interp(center)
                if peak > center_val:
                    raise ValueError("dip of absorption line above continuum.")
                peak = center_val - peak
        # Ensure Gaussian is well sampled by evaluating at the wavelengths within +/- 5
        # sigma of the center. This also prevents overflow errors caused by calculations
        # at wavelengths too far from the center.
        gauss_wavelengths = center + np.arange(-5, 5.05, 0.1) * sigma
        sorted_wavelengths = np.unique(np.concatenate((wavelengths, gauss_wavelengths)))
        sorted_wavelengths = sorted_wavelengths[
            (sorted_wavelengths >= wavelengths[0])
            & (sorted_wavelengths <= wavelengths[-1])
        ]
        sorted_spectrum = spectrum_interp(sorted_wavelengths)
        in_range = (sorted_wavelengths >= gauss_wavelengths[0]) & (
            sorted_wavelengths <= gauss_wavelengths[-1]
        )
        if add:
            sorted_spectrum[in_range] += _gaussian(
                peak, sorted_wavelengths[in_range], center, sigma
            )
        else:
            sorted_spectrum[in_range] -= _gaussian(
                peak, sorted_wavelengths[in_range], center, sigma
            )
        is_good = np.isfinite(sorted_wavelengths) & np.isfinite(sorted_spectrum)
        sorted_wavelengths = sorted_wavelengths[is_good]
        sorted_spectrum = sorted_spectrum[is_good]
        is_negative_spectrum = sorted_spectrum < 0
        if np.any(is_negative_spectrum):
            sorted_spectrum[is_negative_spectrum] = 0.0
            warnings.warn("Setting negative flux in spectrum to zero.", RuntimeWarning)
            # print(
            #     "wavelengths with negative spectrum:",
            #     sorted_wavelengths[is_negative_spectrum],
            # )
        if wavelengths_unit is not None:
            sorted_wavelengths <<= wavelengths_unit  # convert to `astropy.Quantity` array
        return sorted_wavelengths, sorted_spectrum

    @staticmethod
    def _generate_lorentzian(
        wavelengths,
        spectrum,
        center,
        fwhm,
        peak=None,
        tot_flux=None,
        add=True,
        abs_peak=True,
    ):
        """
        Add/subtract a Lorentzian spectrum to/from an existing spectrum. This is useful
        for representing emission lines (i.e., by adding a Lorentzian source) or
        absorption lines (i.e., by subtracting a Lorentzian source). Note that the
        minimum/maximum wavelengths of the source spectrum will not change.

        The Lorentzian spectrum can be represented by the following formulae:
        ```math
                    lorentzian = peak / (1 + num_half_widths^2)
        ```
        and
        ```math
                    num_half_widths = (wavelengths - center) / probable_error
        ```
        and
        ```math
                    peak = tot_flux / (pi * probable_error) <-- see Cauchy distribution
        ```
        and
        ```math
                    probable_error = fwhm / 2
        ```
        where:
          - lorentzian is the Lorentzian spectrum's flux in some arbitrary unit
          - peak is the flux at the center (i.e., central wavelength) of the Lorentzian
          - center is the central wavelength of the Lorentzian
          - wavelengths is the array of wavelengths over which to calculate the spectrum
          - fwhm is the full-width at half-maximum of the Lorentzian
          - tot_flux is the total flux of the Lorentzian under the curve

        Parameters
        ----------
          wavelengths :: array of floats or `astropy.Quantity` array
            The wavelengths over which to calculate the Lorentzian spectrum.

          spectrum :: array of floats
            The spectrum to/from which to add/subtract the Lorentzian spectrum.

          center :: scalar or `astropy.Quantity`
            The central wavelength of the Lorentzian. If a scalar, it is assumed to be in
            angstrom.

          fwhm :: scalar or `astropy.Quantity`
            The full-width at half-maximum of the Lorentzian. If a scalar, it is assumed
            to be in angstrom.

          peak :: int or float
            The peak flux of the Lorentzian (i.e., the flux at the center wavelength).
            Exactly one of peak or tot_flux must be specified.

          tot_flux :: int or float
            The total flux under the curve. Exactly one of peak or tot_flux must be
            specified.

          add :: bool
            If True, add the Lorentzian spectrum to the existing spectrum. If False,
            subtract the Lorentzian from the existing spectrum.

          abs_peak :: bool
            If True, ensure that the peak of the emission line or dip of the absorption
            line is at the given value. Otherwise, just add/subtract the given Lorentzian
            peak to/from the continuum.

        Returns
        -------
          sorted_wavelengths :: array of floats or `astropy.Quantity` array
            The wavelengths of the new spectrum. The shape of this array will be different
            from the input `wavelengths` array.

          sorted_spectrum :: array of floats
            The spectrum with the Lorentzian added/subtracted. The shape of this array
            will be different from the input `spectrum` array.
        """

        def _lorentzian(_peak, _wavelengths, _center, _probable_error):
            _num_half_widths = (_wavelengths - _center) / _probable_error
            return _peak / (1 + _num_half_widths * _num_half_widths)

        #
        # Check inputs
        #
        if np.size(center) != 1:
            raise ValueError("center must be a single scalar or `astropy.Quantity`.")
        if np.size(fwhm) != 1:
            raise ValueError("fwhm must be a single scalar or `astropy.Quantity`.")
        if (peak is None and tot_flux is None) or (
            peak is not None and tot_flux is not None
        ):
            raise ValueError("Exactly one of peak or tot_flux must be specified.")
        if peak is not None and (
            np.size(peak) != 1 or not isinstance(peak, Number) or peak <= 0
        ):
            if add:
                raise ValueError("peak must be a single int or float >= 0.")
            else:
                raise ValueError("dip must be a single int or float >= 0.")
        elif tot_flux is not None and (
            np.size(tot_flux) != 1 or not isinstance(tot_flux, Number) or tot_flux <= 0
        ):
            raise ValueError("tot_flux must be a single int or float >= 0.")
        # Convert lengths to angstrom
        if isinstance(wavelengths, u.Quantity):
            wavelengths_unit = wavelengths.unit
            wavelengths = wavelengths.to(u.AA).value
        else:
            wavelengths_unit = None
        if isinstance(center, u.Quantity):
            center = center.to(u.AA).value
        if isinstance(fwhm, u.Quantity):
            fwhm = fwhm.to(u.AA).value
        if fwhm <= 0:
            raise ValueError("fwhm must be >= 0.")
        #
        # Lorentzian spectrum
        #
        probable_error = 0.5 * fwhm
        if peak is None:
            peak = tot_flux / (np.pi * probable_error)
        spectrum_interp = interp1d(
            wavelengths, spectrum, kind="linear", bounds_error=False, fill_value=np.nan
        )
        if add:
            if abs_peak:
                # Ensure final peak is actually at desired value
                peak -= spectrum_interp(center)
                if peak < 0:
                    raise ValueError("peak of emission line below continuum.")
        else:
            if abs_peak:
                # Ensure final dip is actually at desired value
                center_val = spectrum_interp(center)
                if peak > center_val:
                    raise ValueError("dip of absorption line above continuum.")
                peak = center_val - peak
        # Ensure Lorentzian is well sampled by evaluating at the wavelengths within +/- 80
        # units of probable error from the center. This also prevents overflow errors
        # caused by calculations at wavelengths too far from the center.
        lorentz_wavelengths = center + np.arange(-80, 80.25, 0.5) * probable_error
        sorted_wavelengths = np.unique(np.concatenate((wavelengths, lorentz_wavelengths)))
        sorted_wavelengths = sorted_wavelengths[
            (sorted_wavelengths >= wavelengths[0])
            & (sorted_wavelengths <= wavelengths[-1])
        ]
        sorted_spectrum = spectrum_interp(sorted_wavelengths)
        in_range = (sorted_wavelengths >= lorentz_wavelengths[0]) & (
            sorted_wavelengths <= lorentz_wavelengths[-1]
        )
        if add:
            sorted_spectrum[in_range] += _lorentzian(
                peak, sorted_wavelengths[in_range], center, probable_error
            )
        else:
            sorted_spectrum[in_range] -= _lorentzian(
                peak, sorted_wavelengths[in_range], center, probable_error
            )
        is_good = np.isfinite(sorted_wavelengths) & np.isfinite(sorted_spectrum)
        sorted_wavelengths = sorted_wavelengths[is_good]
        sorted_spectrum = sorted_spectrum[is_good]
        is_negative_spectrum = sorted_spectrum < 0
        if np.any(is_negative_spectrum):
            sorted_spectrum[is_negative_spectrum] = 0.0
            warnings.warn("Setting negative flux in spectrum to zero.", RuntimeWarning)
            # print(
            #     "wavelengths with negative spectrum:",
            #     sorted_wavelengths[is_negative_spectrum],
            # )
        if wavelengths_unit is not None:
            sorted_wavelengths <<= wavelengths_unit  # convert to `astropy.Quantity` array
        return sorted_wavelengths, sorted_spectrum

    def add_emission_line(
        self, center, fwhm, peak=None, tot_flux=None, shape="gaussian", abs_peak=False
    ):
        """
        Add a well-sampled emission line to the spectrum. Note that the minimum/maximum
        wavelengths of the source spectrum will not change.

        For generating an emission line spectrum (i.e., not adding/subtracting an emission
        line to/from a spectrum), see the `generate_emission_line()` method instead.

        N.B. the order in which emission/absorption lines are added will affect the final
        spectrum if using the abs_peak/abs_dip flag. For instance, adding an emission line
        on top of a continuum then specifying an absorption line with an absolute dip is
        not the same as specifying an absorption line with an absolute dip then adding an
        emission line on top of the new continuum.

        Parameters
        ----------
          center :: scalar or `astropy.Quantity`
            The central wavelength of the emission line. If a scalar, it is assumed to be
            in angstrom.

          fwhm :: scalar or `astropy.Quantity`
            The full-width at half-maximum of the emission line. If a scalar, it is
            assumed to be in angstrom.

          peak :: int or float
            The peak flux of the emission line (i.e., the flux at the center wavelength).
            Exactly one of peak or tot_flux must be specified.

          tot_flux :: int or float
            The total flux under the curve. Exactly one of peak or tot_flux must be
            specified.

          shape :: "gaussian" or "lorentzian"
            The emission line profile.

          abs_peak :: bool
            If True, ensure that the peak of the emission line is at the given value.
            Otherwise, just add the given emission line to the spectrum.

        Attributes
        ----------
          wavelengths :: `astropy.Quantity` array
            The wavelengths of the spectrum including the emission line, in angstroms.
            This wavelengths array will have a different shape than the previous
            wavelengths array.

          spectrum :: array of floats
            The spectrum with the emission line added. This spectrum array will have a
            different shape than the previous spectrum array.

        Returns
        -------
          None
        """
        if self.wavelengths is None or self.spectrum is None:
            raise ValueError("Please generate or load a spectrum first")
        if shape == "gaussian":
            spectrum_func = SpectrumMixin._generate_gaussian
        elif shape == "lorentzian":
            spectrum_func = SpectrumMixin._generate_lorentzian
        else:
            raise ValueError("Emission line shape must be 'gaussian' or 'lorentzian'")
        self.wavelengths, self.spectrum = spectrum_func(
            self.wavelengths,
            self.spectrum,
            center,
            fwhm,
            peak=peak,
            tot_flux=tot_flux,
            add=True,
            abs_peak=abs_peak,
        )

    def add_absorption_line(
        self, center, fwhm, dip=None, tot_flux=None, shape="gaussian", abs_dip=False
    ):
        """
        Add a well-sampled absorption line to the spectrum. Note that the minimum/maximum
        wavelengths of the source spectrum will not change.

        N.B. the order in which emission/absorption lines are added will affect the final
        spectrum if using the abs_peak/abs_dip flag. For instance, adding an emission line
        on top of a continuum then specifying an absorption line with an absolute dip is
        not the same as specifying an absorption line with an absolute dip then adding an
        emission line on top of the new continuum.

        Parameters
        ----------
          center :: scalar or `astropy.Quantity`
            The central wavelength of the absorption line. If a scalar, it is assumed to
            be in angstrom.

          fwhm :: scalar or `astropy.Quantity`
            The full-width at half-maximum of the absorption line. If a scalar, it is
            assumed to be in angstrom.

          dip :: int or float
            The minimum flux of the absorption line (i.e., the flux at the center
            wavelength). Exactly one of dip or tot_flux must be specified.

          tot_flux :: int or float
            The total flux under (above) the curve. Exactly one of dip or tot_flux must
            be specified.

          shape :: "gaussian" or "lorentzian"
            The absorption line profile.

          abs_dip :: bool
            If True, ensure that the dip of the absorption line is at the given value.
            Otherwise, just subtract the given absorption line from the spectrum.

        Attributes
        ----------
          wavelengths :: `astropy.Quantity` array
            The wavelengths of the spectrum including the absorption line, in angstroms.
            This wavelengths array will have a different shape than the previous
            wavelengths array.

          spectrum :: array of floats
            The spectrum with the absorption line subtracted. This spectrum array will
            have a different shape than the previous spectrum array.

        Returns
        -------
          None
        """
        if self.wavelengths is None or self.spectrum is None:
            raise ValueError("Please generate or load a spectrum first")
        if shape == "gaussian":
            spectrum_func = SpectrumMixin._generate_gaussian
        elif shape == "lorentzian":
            spectrum_func = SpectrumMixin._generate_lorentzian
        else:
            raise ValueError("Absorption line shape must be 'gaussian' or 'lorentzian'")
        self.wavelengths, self.spectrum = spectrum_func(
            self.wavelengths,
            self.spectrum,
            center,
            fwhm,
            peak=dip,
            tot_flux=tot_flux,
            add=False,
            abs_peak=abs_dip,
        )

    def generate_emission_line(
        self,
        center,
        fwhm,
        peak=None,
        tot_flux=None,
        shape="gaussian",
        limits=[100, 1200] << u.nm,
        overwrite=False,
        quiet=False,
    ):
        """
        Generate a spectrum representing a single emission line. The resolution of the
        spectrum is at least 1% of the wavelength range.

        To add/subtract a spectral line to/from a spectrum, see the `add_emission_line()`
        and `add_absorption_line()` methods instead.

        Parameters
        ----------
          center :: scalar or `astropy.Quantity`
            The central wavelength of the emission line. If a scalar, it is assumed to be
            in angstrom.

          fwhm :: scalar or `astropy.Quantity`
            The full-width at half-maximum of the emission line. If a scalar, it is
            assumed to be in angstrom.

          peak :: int or float
            The peak flux of the emission line (i.e., the flux at the center wavelength).
            Exactly one of peak or tot_flux must be specified.

          tot_flux :: int or float
            The total flux under the curve. Exactly one of peak or tot_flux must be
            specified.

          shape :: "gaussian" or "lorentzian"
            The emission line profile.

          limits :: 2-element 1D array of `astropy.Quantity` or scalars
            The [min, max] wavelengths of the spectrum. If the elements are scalars, it is
            assumed to be in angstrom. The center of the emission line must be within
            these limits.

          overwrite :: bool
            If True, overwrite any existing wavelengths/spectrum. If False, raise an error
            if wavelengths or spectrum is not None.

          quiet :: bool
            If True, do not print a message when overwriting an existing spectrum.

        Attributes
        ----------
          wavelengths :: `astropy.Quantity` array
            The wavelengths of the spectrum, in angstroms.

          spectrum :: array of floats
            The emission line spectrum in units of flam (erg/s/cm^2/A).

        Returns
        -------
          None
        """
        #
        # Check inputs
        #
        self._check_existing_spectrum(overwrite, quiet=quiet)
        if isinstance(center, u.Quantity):
            center = center.to(u.AA).value
        limits_AA = []
        for limit in limits:
            if isinstance(limit, u.Quantity):
                limit = limit.to(u.AA).value
            limits_AA.append(limit)
        if np.shape(limits_AA) != (2,):
            raise ValueError(
                "limits must be a 2-element 1D array of `astropy.Quantity` or scalars"
            )
        if (center < limits_AA[0]) or (center > limits_AA[1]):
            raise ValueError("center must be within limits")
        if shape == "gaussian":
            spectrum_func = SpectrumMixin._generate_gaussian
        elif shape == "lorentzian":
            spectrum_func = SpectrumMixin._generate_lorentzian
        else:
            raise ValueError("Emission line shape must be 'gaussian' or 'lorentzian'")
        #
        # Generate spectrum
        #
        # Assume flat (zero) spectrum as baseline. Then add emission line on top of this
        resolution = (limits_AA[1] - limits_AA[0]) / 100
        wavelengths = np.arange(limits_AA[0], limits_AA[1] + 0.5 * resolution, resolution)
        self.wavelengths, self.spectrum = spectrum_func(
            wavelengths=wavelengths * u.AA,
            spectrum=np.zeros(len(wavelengths), dtype=float),
            center=center,
            fwhm=fwhm,
            peak=peak,
            tot_flux=tot_flux,
            add=True,
            abs_peak=False,
        )

    def set_spectrum(self, wavelengths, spectrum, unit, overwrite=False, quiet=False):
        """
        Set the spectrum of the source based on the input arrays. To use a spectrum from a
        file, see the `use_custom_spectrum()` method.

        The input spectrum should have units of either flam (erg/s/cm^2/A), fnu
        (erg/s/cm^2/Hz), ABmag (AB magnitude), or STmag (ST magnitude). Note that the
        computed (and stored) spectrum will always be in units of flam.

        Parameters
        ----------
          wavelengths :: array of scalars or `astropy.Quantity` array
            The wavelengths over which to generate the uniform spectrum. If an array of
            scalars, it should be in angstrom. This should be a 1D array with the same
            length as the `spectrum` array.

          specturm :: array of scalars
            The value of the spectrum at the given wavelengths, in units of `unit`. This
            should be a 1D array with the same length as the `wavelengths` array.

          unit :: "flam" or "fnu" or "ABmag" or "STmag"
            The unit of the `spectrum` array: either flam (erg/s/cm^2/A), fnu
            (erg/s/cm^2/Hz), ABmag (AB magnitude), or STmag (ST magnitude).

          overwrite :: bool
            If True, overwrite any existing wavelengths/spectrum. If False, raise an error
            if wavelengths or spectrum is not None.

          quiet :: bool
            If True, do not print a message when overwriting an existing spectrum.

        Attributes
        ----------
          wavelengths :: astropy.Quantity` array
            The given wavelengths of the spectrum, converted into units of angstrom.

          spectrum :: array of floats
            The given spectrum, converted into units of erg/s/cm^2/A.

        Returns
        -------
          None
        """
        #
        # Check inputs
        #
        self._check_existing_spectrum(overwrite, quiet=quiet)
        if not isinstance(wavelengths, u.Quantity):
            wavelengths = wavelengths * u.AA
        else:
            wavelengths = wavelengths.to(u.AA)
        if isinstance(spectrum, u.Quantity):
            raise TypeError("`spectrum` must be an array of scalars")
        if np.shape(spectrum) != np.shape(wavelengths) or np.ndim(spectrum) != 1:
            raise ValueError(
                "`wavelengths` and `spectrum` must be 1D arrays of the same shape"
            )
        if unit not in ["flam", "fnu", "ABmag", "STmag"]:
            raise ValueError("`unit` must be one of 'flam', 'fnu', 'ABmag', or 'STmag'")
        if unit == "flam" or unit == "fnu":
            if np.any(spectrum) <= 0:
                raise ValueError(
                    "All `spectrum` values must be > 0 if `unit` is 'flam' or 'fnu'"
                )
        #
        # Convert spectrum to units of flam (erg/s/cm^2/A)
        #
        if unit == "fnu":
            spectrum = fnu_to_flam(
                fnu=spectrum, wavelength=wavelengths, fnu_err=0.0, wavelength_err=0.0
            )[0]
        elif unit == "ABmag":
            # Convert to fnu
            spectrum = mag_to_flux(mag=spectrum, mag_err=0.0, zpt=-48.60)[0]
            # Convert fnu to flam
            spectrum = fnu_to_flam(
                fnu=spectrum, wavelength=wavelengths, fnu_err=0.0, wavelength_err=0.0
            )[0]
        elif unit == "STmag":
            # Convert directly to flam
            spectrum = mag_to_flux(mag=spectrum, mag_err=0.0, zpt=-21.10)[0]
        self.wavelengths = wavelengths
        self.spectrum = spectrum

    def use_custom_spectrum(
        self, filepath, wavelength_unit=u.AA, overwrite=False, quiet=False
    ):
        """
        Use custom spectrum from an ASCII or FITS file. To use a spectrum from an array,
        use the `set_spectrum()` method.

        Parameters
        ----------
          filepath :: str
            The absolute path to the file containing the spectrum.
            If the file is in ASCII format, the first column should contain the
            wavelengths in `wavelength_units` and the second column containing the
            spectrum in flam (erg/s/cm^2/A); the columns should be separated by a constant
            number of spaces. Lines starting with a hash (#) will be ignored. The file
            extension must not be .fit or .fits.
            If the file is in FITS format, the first field (index 0) should contain the
            wavelengths in `wavelength_units` and the second field (index 1) should
            contain the spectrum in flam (erg/s/cm^2/A). The file extension must be .fit
            or .fits.

          wavelength_unit :: `astropy.Quantity` length unit
            The unit of the wavelengths in the file (e.g., u.AA for angstrom, u.nm for
            nanometer, u.um for micrometer, etc.)

          overwrite :: bool
            If True, overwrite any existing wavelengths/spectrum. If False, raise an error
            if wavelengths or spectrum is not None.

          quiet :: bool
            If True, do not print a message when overwriting an existing spectrum.

        Attributes
        ----------
          wavelengths :: `astropy.Quantity` array
            The wavelengths of the spectrum, in angstroms.

          spectrum :: array of floats
            Source spectrum in units of flam (erg/s/cm^2/A).

        Returns
        -------
          None
        """
        #
        # Check inputs
        #
        self._check_existing_spectrum(overwrite, quiet=quiet)
        if not isinstance(filepath, str):
            raise TypeError("filepath must be a string.")
        try:
            _ = wavelength_unit.to(u.AA)
        except Exception:
            raise TypeError("wavelength_units must be an `astropy.Quantity` length unit.")
        #
        # Load spectrum
        #
        file_ext = filepath.split(".")[-1].lower()
        try:
            if file_ext == "fit" or file_ext == "fits":
                data = fits.getdata(filepath)
                self.wavelengths = (data.field(0)[0] * wavelength_unit).to(u.AA)
                self.spectrum = data.field(1)[0] # index 0 to access the corresponding arrays; data variable contains embedded arrays
            else:
                data = pd.read_csv(
                    filepath,
                    sep="\s+",
                    header=None,
                    comment="#",
                    engine="python",
                )  # sep=" +" is Python regex to match a variable number of spaces
                self.wavelengths = (data[0].values * wavelength_unit).to(u.AA)
                self.spectrum = data[1].values
        except Exception:
            raise RuntimeError(
                "Could not read spectrum from file. File must be in ASCII or FITS format "
                + "and adhere to the guidelines specified in the docstring."
            )

    def use_galaxy_spectrum(self, gal_type, overwrite=False, quiet=False):
        """
        Use one of the predefined galaxy spectra. These non-uniformly sampled spectra are
        from Fioc & Rocca-Volmerange (1997)
        <https://ui.adsabs.harvard.edu/abs/1997A%26A...326..950F/abstract>. In particular,
        the data files were downloaded from the Gemini Observatory Control Software (OCS)
        GitHub repository:
        <https://github.com/gemini-hlsw/ocs/blob/develop/bundle/edu.gemini.itc/src/main/resources/sed/non_stellar/elliptical-galaxy.nm>
        and
        <https://github.com/gemini-hlsw/ocs/blob/develop/bundle/edu.gemini.itc/src/main/resources/sed/non_stellar/spiral-galaxy.nm>.

        Parameters
        ----------
          gal_type :: "elliptical" or "spiral"
            The galaxy morphology. The elliptical galaxy (class T=-5, -4) and spiral
            galaxy (type Sc, class T=5) spectra both run from 22-9698 nm.

          overwrite :: bool
            If True, overwrite any existing wavelengths/spectrum. If False, raise an error
            if wavelengths or spectrum is not None.

          quiet :: bool
            If True, do not print a message when overwriting an existing spectrum.

        Attributes
        ----------
          wavelengths :: `astropy.Quantity` array
            The wavelengths of the spectrum, in angstroms.

          spectrum :: array of floats
            Galaxy spectrum in arbitrary units, normalized such that it is equal to 1 at
            550 nm.

        Returns
        -------
          None
        """
        #
        # Check inputs
        #
        self._check_existing_spectrum(overwrite, quiet=quiet)
        if gal_type == "elliptical" or gal_type == "spiral":
            filepath = join(DATAPATH, "galaxy_spectra", f"{gal_type}_galaxy.txt")
        else:
            raise ValueError("Galaxy type must be 'elliptical' or 'spiral'")
        data = pd.read_csv(
            filepath,
            sep=" +",
            header=None,
            comment="#",
            engine="python",
        )  # sep=" +" is Python regex to match a variable number of spaces
        self.wavelengths = (data[0].values * u.nm).to(u.AA)
        self.spectrum = data[1].values

    def _calc_xy(self):
        """
        Intenal function that converts (ra,dec) of the identified Gaia sources to pixel positions on the CCD. 

        Attributes
        ----------
          gaia :: dict of arrays
            Dictionary containing the interpolated parameters of the queried Gaia sources.

        Return
        ------
            None
        """

        # Get projected CCD coordinates (centered on target)
        wcs_input_dict = {
            'CTYPE1': 'RA---TAN', 
            'CUNIT1': 'deg', 
            'CDELT1': self.fov.value/self.ccd_dim[0], 
            'CRPIX1': int(self.ccd_dim[0]/2), 
            'CRVAL1': self.gaia['ra'][0], 
            'NAXIS1': self.ccd_dim[0],
            'CTYPE2': 'DEC--TAN', 
            'CUNIT2': 'deg', 
            'CDELT2': self.fov.value/self.ccd_dim[1], 
            'CRPIX2': int(self.ccd_dim[1]/2), 
            'CRVAL2': self.gaia['dec'][0], 
            'NAXIS2': self.ccd_dim[1],
            'CROTA2': self.fov_pa.value,
        }

        wcs = WCS(wcs_input_dict)

        # Convert Gaia source ra,dec -> x,y
        g_coord = SkyCoord(ra=self.gaia['ra'], dec=self.gaia['dec'], 
                            unit=(u.degree, u.degree), frame='icrs')
        _x, _y = wcs.world_to_pixel(g_coord)

        self.gaia['x'], self.gaia['y'] = _x, _y
        self.gaia['wcs'] = wcs

    def _search_gaia(self):
        """
        Internal Function. Search Gaia catalog for target and background stars.

        Parameters
        ----------
            TelescopeObj :: `castor_etc.Telescope` instance
              The `Telescope` object passed down to the _calc_xy function.

        Attributes
        ----------
          gaia :: dict of arrays
            Dictionary containing the interpolated parameters of the queried Gaia sources.

        Returns
        -------
          None

        """
        try: 
            srch_str = "SELECT *, DISTANCE(POINT({:.6f},{:.6f}), POINT(ra,dec))".format(self.ra.value,self.dec.value) + "AS ang_sep FROM gaiadr2.gaia_source " + "WHERE 1 = CONTAINS( POINT({:.6f},{:.6f}), ".format(self.ra.value,self.dec.value) + "CIRCLE(ra,dec,{:.2f}))".format(self.srch_rad.value) + "AND phot_g_mean_mag <={:.2f}".format(self.srch_Gmax) + "AND parallax IS NOT NULL ORDER BY ang_sep ASC"

            job = Gaia.launch_job(srch_str)
            results = job.get_results()

            ra, dec, Gmag, Bpmag, Rpmag = np.array(results['ra']), np.array(results['dec']), np.array(results['phot_g_mean_mag']), np.array(results['phot_bp_mean_mag']), np.array(results['phot_rp_mean_mag'])

            print('{:.0f} Gaia source(s) found'.format(len(ra)))

            # Identify target within Gaia search (brightest target within sep_max)
            sep_max = (8.0 / 3600.) #* u.deg
            c_coord = SkyCoord(ra=self.ra.value, dec=self.dec.value, unit=(u.degree,u.degree), frame='icrs')
            g_coord = SkyCoord(ra=ra,dec=dec,unit=(u.degree,u.degree), frame='icrs')
            sep = g_coord.separation(c_coord).degree
            ti = np.where(sep < sep_max)[0]

            if len(ti) > 0:
                ti = ti[ np.argmin(Gmag[ti]) ]
            else:
                ti = np.argmin(sep)
                print('No Gaia soueces found within {:.2f} arcsec of specified RA and DEC'.format(sep_max*3600.))
            
            gi = np.delete(np.arange(len(ra)),ti)

            # Place target at top of list, sort remaining sources by brightness.
            si = np.argsort(Gmag[gi])
            scalars = [ra,dec,Gmag,Bpmag,Rpmag,sep]
            for scalar in scalars:
                scalar = np.hstack([ scalar[ti], scalar[gi[si]] ])

            nmax = self.srch_nmax # I believe this way is more efficient than using self.srch_nmax multiple time <Need to check>.

            if len(ra) > nmax:
                ra = ra[:nmax]
                dec = dec[:nmax]
                Gmag = Gmag[:nmax]
                Bpmag = Bpmag[:nmax]
                Rpmag = Rpmag[:nmax]
                sep = sep[:nmax]

            self.gaia = {'ra': ra, 'dec': dec, 'x': np.zeros_like(ra), 'y': np.zeros_like(ra), 'Gmag': Gmag, 'Bpmag': Bpmag, 'Rpmag': Rpmag, 'results': None, 'wcs': None}

            self._calc_xy()

        except Exception:
            raise RuntimeError('Gaia search failed. Please revise the input target parameters.')

    def _specify_target_parameters(self,run_gaia_search = True):
        """
        Internal function. Used to specify target Gaia parameters and populate the gaia attribute of the PointSource. If temperature, Teff, and Gaia G magnitude, Gmag, of the target are specified, then Gaia catalog query will not run. 

        Parameters
        ----------
          run_gaia_search :: bool
            If Teff and Gmag of the target are specified, then this attributes is set to False and gaia query run is called off.
        
        Attributes
        ----------
          gaia :: dict of arrays
            Dictionary containing the interpolated parameters of the queried Gaia sources.

        Returns
        -------
          None

        """

        xout = self.ccd_dim[0]
        yout = self.ccd_dim[1]

        if (self.Gmag != None) & (self.Teff != None):
            run_gaia_search = False

        if run_gaia_search:
            self._search_gaia()
        else:
            self.gaia = {
                'ra': np.array([self.ra.value]), 
                'dec': np.array([self.dec.value]),
                'x': np.zeros(1) + xout/2,
                'y': np.zeros(1) + yout/2, 
                'Gmag': np.array([np.nan]),
                'Bpmag': np.array([np.nan]), 
                'Rpmag': np.array([np.nan]), 
                'Gmag_abs': np.array([np.nan]),
                'logg': np.array([np.nan]),
                'radius': np.array([np.nan]),
                'metallicity': np.array([np.nan]),
                'Teff': np.array([np.nan]), 
                'results': None, # Full results from Gaia query
                'wcs': None # wcs.world_to_pixel object used for projection transformation    
            }

            self._calc_xy()
        

        # Replace target parameters with user-specified values
        scalars = [self.Teff, self.Gmag, self.logg, self.radius, self.metallicity, self.Bpmag, self.Rpmag]
        scalars_str = ["Teff", "Gmag", "logg", "radius", "metallicity","Bpmag","Rpmag"]
        for scalar, scalar_str in zip(scalars, scalars_str):
            if scalar != None:    
              self.gaia[scalar_str][0] = scalar

        # Interpolate EEM_table using Teff or (Bpmag, Rpmag)
        if 'Teff' in self.gaia.keys():
            interp_eem = interp_EEM_table(Teff=self.gaia['Teff'], Gmag=self.gaia['Gmag'])
        else:
            interp_eem = interp_EEM_table(Gmag=self.gaia['Gmag'], Bpmag=self.gaia['Bpmag'], Rpmag=self.gaia['Rpmag'])
            
        if run_gaia_search:   
          self.gaia['Teff'] = interp_eem['Teff']
          self.gaia['logg'] = interp_eem['logg']
          self.gaia['radius'] = interp_eem['radius']
          self.gaia['metallicity'] = np.zeros_like(self.gaia['Teff'])
          self.gaia['Gmag_abs'] = interp_eem['Gmag_abs']

        # If certain values are missing, replace with interpolated EEM_table values
        for _lbl in ['Bpmag','Rpmag','logg','radius','Gmag_abs']:
            ti  = np.isnan(self.gaia[_lbl])
            self.gaia[_lbl][ti] = interp_eem[_lbl][ti]

        # Replace missing metallicities with solar value
        ti = np.isnan(self.gaia['metallicity'])
        self.gaia['metallicity'][ti] = 0.

    def use_gaia_spectrum(self, TelescopeObj, ra = None, dec = None, srch_Gmax=21., srch_nmax = 100, srch_rad = None, Teff = None, Gmag = None, logg = None, radius = None, metallicity = None, Bpmag = None, Rpmag = None, stellar_model_grid = 'ATLAS9', bkg_sources = True, fov = None, fov_pa = 0 * u.deg, overwrite=False, quiet=False):
        """
        Use a spectrum from either the ATLAS9 catalog <Ask James for the link> or the Bt-Sett1 catalog <Ask James for the link> containing flux spectra for numerous stellar model atmospheres. 

        Parameters
        ----------
          TelescopeObj :: `castor_etc.Telescope` instance
            'The `Telescope` object containing the field of view (FoV) value required to calculate the search radius for querying the Gaia database.

          ra :: `astropy.Quantity` int or float
            Right ascension of the target in degree

          dec  :: `astropy.Quantity` int or float
            Declination of the target in degree

          srch_Gmax :: int or float
            Maximum Gaia G magnitude for Gaia catalog query (applies to both the target and the guide stars)

          srch_nmax :: int
            Maximum Gaia sources to include.

          srch_rad :: `astropy.Quantity`
            Search radius, in degree, for Gaia catalog query (applies to both the target and the guide stars)

          Teff :: int or float 
            Effective temperature of the target in Kelvin

          Gmag :: int or float
            Gaia G magnitude of the target used to query the Gaia catalog 

          logg :: int or float
            log(g) value of the target

          radius :: int or float
            radius of the target <ask james about units>

          metallicity :: int or float
            metallicity of the target

          Bpmag :: int or float
            Gaia BP magnitude of the target

          Rpmag :: int or float
            Gaia RP magnitude of the target

          stellar_model_grid :: str
            Stellar model grid, 'ATLAS9' or 'BTSettl', for selecting spectrum of the source according to the interpolated stellar atmosphere model.
          
          bkg_sources :: bool
            If True, then background Gaia sources are included during the transit simulation calculation. If False, nmax is set to 1.

          fov :: int or float
            Full width FoV in degree

          fov_pa :: `astropy.Quantity`
            Field of view position angle

          overwrite :: bool
            If True, overwrite any existing wavelengths/spectrum. If False, raise an error
            if wavelengths or spectrum is not None.

          quiet :: bool
            If True, do not print a message when overwriting an existing spectrum.

        Attributes
        ----------
          wavelengths :: `astropy.Quantity` array
            The wavelengths of the spectrum, in angstroms,

          spectrum :: array of floats
            Stellar spectrum in flam units

          ra :: `astropy.Quantity` angle
            Right ascension of the target in degree

          dec  :: `astropy.Quantity` angle
            Declination of the target in degree

          srch_Gmax :: int or float
            Maximum Gaia G magnitude for Gaia catalog query (applies to both the target and the guide stars)

          srch_nmax :: int
            Maximum Gaia sources to include.

          srch_rad :: `astropy.Quantity`
            Search radius, in degree, for Gaia catalog query (applies to both the target and the guide stars)

          Teff :: int or float 
            Effective temperature of the target in Kelvin

          Gmag :: int or float
            Gaia G magnitude of the target used to query the Gaia catalog 

          logg :: int or float
            log(g) value of the target

          radius :: int or float
            radius of the target <ask james about units>

          metallicity :: int or float
            metallicity of the target

          Bpmag :: int or float
            Gaia BP magnitude of the target

          Rpmag :: int or float
            Gaia RP magnitude of the target

          stellar_model_grid :: str
            Stellar model grid, 'ATLAS9' or 'BTSettl', for selecting spectrum of the source according to the interpolated stellar atmosphere model.
          
          bkg_sources :: bool
            If True, then background Gaia sources are included during the transit simulation calculation. If False, nmax is set to 1.

          fov :: int or float
            Full width FoV in degree

          fov_pa :: `astropy.Quantity`
            Field of view position angle

          ccd_dim :: int list
            CCD dimesions adopted from TelescopeObj 

          Returns
          -------
            None

        """

        if not isinstance(stellar_model_grid, str):
            raise TypeError(
                "stellar_model_grid should be a string"
                + "(e.g., 'ATLAS9' or 'BTSettl')"
                )
        
        if not isinstance(bkg_sources, bool):
            raise TypeError(
                "bkg_sources should be either True or False"
                )
        
        if not isinstance(srch_nmax, int):
            raise TypeError("srch_nmax must an integer number.")
            
        scalars = [ra,dec,fov,fov_pa,srch_rad]
        scalars_str = ["ra","dec","fov","fov_pa","srch_rad"]
        for scalar,scalar_str in zip(scalars,scalars_str):
            if scalar != None:
                if not isinstance(scalar, u.Quantity):
                      raise TypeError(
                          f"{scalar_str} must `astropy.Quantity` angles" + "(e.g., u.arcsec, u.deg)"
                      )

        scalars = [srch_Gmax, Teff, Gmag, logg, radius, metallicity, Bpmag, Rpmag]
        scalars_str = ["srch_Gmax","Teff", "Gmag", "logg", "radius", "metallicity","Bpmag","Rpmag"]
        for scalar, scalar_str in zip(scalars, scalars_str):
            if scalar != None:    
              if not isinstance(scalar, Number):
                  raise TypeError(f"{scalar_str} must be an int or float")
              
        if fov == None:
          fov = TelescopeObj.fov.to(u.deg) 
        
        if srch_rad == None:
            srch_rad = (fov / 2.) * np.sqrt(2.)

        if bkg_sources == False:
            srch_nmax = 1 # Only the target Gaia source.

        # Intializing PointSource attributes after checking the inputs
        self.ra = ra.to(u.deg)
        self.dec = dec.to(u.deg)
        self.srch_Gmax = srch_Gmax
        self.srch_nmax = srch_nmax
        self.srch_rad = srch_rad.to(u.deg)
        self.Teff = Teff
        self.Gmag = Gmag
        self.logg = logg
        self.radius = radius
        self.metallicity = metallicity
        self.Bpmag = Bpmag
        self.Rpmag = Rpmag
        self.stellar_model_grid = stellar_model_grid
        self.bkg_sources = bkg_sources
        self.fov = fov.to(u.deg)
        self.fov_pa = fov_pa.to(u.deg)
        self.ccd_dim = TelescopeObj.ccd_dim

  
        # Performs a gaia query and interpolates potential sources' parameters, i.e., Teff, radius...
        self._specify_target_parameters()

        self._check_existing_spectrum(overwrite,quiet=quiet)

        try:
          # Zero index because we are only concerned with the target here. 
          wavelengths, spectrum = getStarData(
              temperature = self.gaia['Teff'][0],
              metallicity = self.gaia['metallicity'][0],
              logg = self.gaia['logg'][0],
              model_grid = stellar_model_grid
              )

          # Estimate distance using Gmag and estimated Gabs
          dpc = 10**( (self.gaia['Gmag'][0] - self.gaia['Gmag_abs'][0])/5. + 1 )

          # Scale surface flux to observed flux using esimated dpc and radius
          spectrum *= ( ( self.gaia['radius'][0] * R_sun ) / ( dpc * pc ) )**2 

          # CASTOR wavelength range.
          _range = np.where ( (wavelengths.to(u.AA).value <= 11000)  )

          self.wavelengths = wavelengths[_range]
          self.spectrum = spectrum[_range]

        except Exception:
            raise RuntimeError(
                f"Could not load the {stellar_model_grid} data for some reason. Please contact the developer with a minimal" + "working example."
            )

    def use_pickles_spectrum(self, spectral_class, overwrite=False, quiet=False):
        """
        Use a spectrum from the Pickles catalog
        (<https://ui.adsabs.harvard.edu/abs/1998PASP..110..863P/abstract>) containing
        spectra for numerous stellar spectral classes.

        A table containing valid `spectral_class` inputs is at the end of this docstring.

        Parameters
        ----------
          spectral_class :: str from the "Table of valid `spectral_class` inputs" below
            The spectral type of the star.

          overwrite :: bool
            If True, overwrite any existing wavelengths/spectrum. If False, raise an error
            if wavelengths or spectrum is not None.

          quiet :: bool
            If True, do not print a message when overwriting an existing spectrum.

        Attributes
        ----------
          wavelengths :: `astropy.Quantity` array
            The wavelengths of the spectrum, in angstroms.

          spectrum :: array of floats
            Stellar spectrum in arbitrary units, normalized such that it is equal to 1 at
            5556 angstrom.

        Returns
        -------
          None

        Table of valid `spectral_class` inputs
        -------------------------------------
        ```text
          spectral_class      Description (& wavelength range)
          --------------      --------------------------------
          "a0i"               A0 I     (1150-10620 A)
          "a0iii"             A0 III   (1150-10620 A)
          "a0iv"              A0 IV    (1150-10620 A)
          "a0v"               A0 V     (1150-10620 A)
          "a2i"               A2 I     (1150-10620 A)
          "a2v"               A2 V     (1150-10620 A)
          "a3iii"             A3 III   (1150-10620 A)
          "a3v"               A3 V     (1150-10620 A)
          "a47iv"             A4-7 IV  (1150-10620 A)
          "a5iii"             A5 III   (1150-10620 A)
          "a5v"               A5 V     (1150-10620 A)
          "a7iii"             A7 III   (1150-10620 A)
          "a7v"               A7 V     (1150-10620 A)
          "b0i"               B0 I     (1150-10620 A)
          "b0v"               B0 V     (1150-10620 A)
          "b12iii"            B1-2 III (1150-10620 A)
          "b1i"               B1 I     (1150-10620 A)
          "b1v"               B1 V     (1150-10620 A)
          "b2ii"              B2 II    (1150-10620 A)
          "b2iv"              B2 IV    (1150-10620 A)
          "b3i"               B3 I     (1150-10620 A)
          "b3iii"             B3 III   (1150-10620 A)
          "b3v"               B3 V     (1150-10620 A)
          "b57v"              B5-7 V   (1150-10620 A)
          "b5i"               B5 I     (1150-10620 A)
          "b5ii"              B5 II    (1150-10620 A)
          "b5iii"             B5 III   (1150-10620 A)
          "b6iv"              B6 IV    (1150-10620 A)
          "b8i"               B8 I     (1150-10620 A)
          "b8v"               B8 V     (1150-10620 A)
          "b9iii"             B9 III   (1150-10620 A)
          "b9v"               B9 V     (1150-10620 A)
          "f02iv"             F0-2 IV  (1150-10620 A)
          "f0i"               F0 I     (1150-10620 A)
          "f0ii"              F0 II    (1150-10620 A)
          "f0iii"             F0 III   (1150-10620 A)
          "f0v"               F0 V     (1150-10620 A)
          "f2ii"              F2 II    (1150-10620 A)
          "f2iii"             F2 III   (1150-10620 A)
          "f2v"               F2 V     (1150-10620 A)
          "f5i"               F5 I     (1150-10620 A)
          "f5iii"             F5 III   (1150-10620 A)
          "f5iv"              F5 IV    (1150-10620 A)
          "f5v"               F5 V     (1150-10620 A)
          "f6v"               F6 V     (1150-10620 A)
          "f8i"               F8 I     (1150-10620 A)
          "f8iv"              F8 IV    (1150-10620 A)
          "f8v"               F8 V     (1150-10620 A)
          "g0i"               G0 I     (1150-10620 A)
          "g0iii"             G0 III   (1150-10620 A)
          "g0iv"              G0 IV    (1150-10620 A)
          "g0v"               G0 V     (1150-10620 A)
          "g2i"               G2 I     (1150-10620 A)
          "g2iv"              G2 IV    (1150-10620 A)
          "g2v"               G2 V     (1150-10620 A)
          "g5i"               G5 I     (1150-10620 A)
          "g5ii"              G5 II    (1150-10620 A)
          "g5iii"             G5 III   (1150-10620 A)
          "g5iv"              G5 IV    (1150-10620 A)
          "g5v"               G5 V     (1150-10620 A)
          "g8i"               G8 I     (1150-10620 A)
          "g8iii"             G8 III   (1150-10620 A)
          "g8iv"              G8 IV    (1150-10620 A)
          "g8v"               G8 V     (1150-10620 A)
          "k01ii"             K0-1 II  (1150-10620 A)
          "k0iii"             K0 III   (1150-10620 A)
          "k0iv"              K0 IV    (1150-10620 A)
          "k0v"               K0 V     (1150-10620 A)
          "k1iii"             K1 III   (1150-10620 A)
          "k1iv"              K1 IV    (1150-10620 A)
          "k2i"               K2 I     (1150-10620 A)
          "k2iii"             K2 III   (1150-10620 A)
          "k2v"               K2 V     (1150-10620 A)
          "k34ii"             K3-4 II  (1150-10620 A)
          "k3i"               K3 I     (1150-10620 A)
          "k3iii"             K3 III   (1150-10620 A)
          "k3iv"              K3 IV    (1150-10620 A)
          "k3v"               K3 V     (1150-10620 A)
          "k4i"               K4 I     (1150-10620 A)
          "k4iii"             K4 III   (1150-10620 A)
          "k4v"               K4 V     (1150-10620 A)
          "k5iii"             K5 III   (1150-10620 A)
          "k5v"               K5 V     (1150-10620 A)
          "k7v"               K7 V     (1150-10620 A)
          "m0iii"             M0 III   (1150-10620 A)
          "m0v"               M0 V     (1150-10620 A)
          "m10iii"            M10 III  (1150-10620 A)
          "m1iii"             M1 III   (1150-10620 A)
          "m1v"               M1 V     (1150-10620 A)
          "m2i"               M2 I     (1150-10620 A)
          "m2iii"             M2 III   (1150-10620 A)
          "m2p5v"             M2.5 V   (1150-10620 A)
          "m2v"               M2 V     (1150-10620 A)
          "m3ii"              M3 II    (1150-10620 A)
          "m3iii"             M3 III   (1150-10620 A)
          "m3v"               M3 V     (1150-10620 A)
          "m4iii"             M4 III   (1150-10620 A)
          "m4v"               M4 V     (1150-10620 A)
          "m5iii"             M5 III   (1150-10620 A)
          "m5v"               M5 V     (1150-10620 A)
          "m6iii"             M6 III   (1150-10620 A)
          "m6v"               M6 V     (1150-10620 A)
          "m7iii"             M7 III   (1150-10620 A)
          "m8iii"             M8 III   (1150-10620 A)
          "m9iii"             M9 III   (1150-10620 A)
          "o5v"               O5 V     (1150-10620 A)
          "o8iii"             O8 III   (1150-10620 A)
          "o9v"               O9 V     (1150-10620 A)
          "rf6v"              metal-rich F6 V    (1150-10620 A)
          "rf8v"              metal-rich F8 V    (1150-10620 A)
          "rg0v"              metal-rich G0 V    (1150-10620 A)
          "rg5iii"            metal-rich G5 III  (1150-10620 A)
          "rg5v"              metal-rich G5 V    (1150-10620 A)
          "rk0iii"            metal-rich K0 III  (1150-10620 A)
          "rk0v"              metal-rich K0 V    (1150-10620 A)
          "rk1iii"            metal-rich K1 III  (1150-10620 A)
          "rk2iii"            metal-rich K2 III  (1150-10620 A)
          "rk3iii"            metal-rich K3 III  (1150-10620 A)
          "rk4iii"            metal-rich K4 III  (1150-10620 A)
          "rk5iii"            metal-rich K5 III  (1150-10620 A)
          "uka0i"             A0 I     (1150-25000 A)
          "uka0iii"           A0 III   (1150-25000 A)
          "uka0iv"            A0 IV    (1150-25000 A)
          "uka0v"             A0 V     (1150-25000 A)
          "uka2i"             A2 I     (1150-25000 A)
          "uka2v"             A2 V     (1150-25000 A)
          "uka3iii"           A3 III   (1150-25000 A)
          "uka3v"             A3 V     (1150-25000 A)
          "uka47iv"           A4-7 IV  (1150-25000 A)
          "uka5iii"           A5 III   (1150-25000 A)
          "uka5v"             A5 V     (1150-25000 A)
          "uka7iii"           A7 III   (1150-25000 A)
          "uka7v"             A7 V     (1150-25000 A)
          "ukb0i"             B0 I     (1150-25000 A)
          "ukb0v"             B0 V     (1150-25000 A)
          "ukb12iii"          B1-2 III (1150-25000 A)
          "ukb1i"             B1 I     (1150-25000 A)
          "ukb1v"             B1 V     (1150-25000 A)
          "ukb2ii"            B2 II    (1150-25000 A)
          "ukb2iv"            B2 IV    (1150-25000 A)
          "ukb3i"             B3 I     (1150-25000 A)
          "ukb3iii"           B3 III   (1150-25000 A)
          "ukb3v"             B3 V     (1150-25000 A)
          "ukb57v"            B5-7 V   (1150-25000 A)
          "ukb5i"             B5 I     (1150-25000 A)
          "ukb5ii"            B5 II    (1150-25000 A)
          "ukb5iii"           B5 III   (1150-25000 A)
          "ukb6iv"            B6 IV    (1150-25000 A)
          "ukb8i"             B8 I     (1150-25000 A)
          "ukb8v"             B8 V     (1150-25000 A)
          "ukb9iii"           B9 III   (1150-25000 A)
          "ukb9v"             B9 V     (1150-25000 A)
          "ukf02iv"           F0-2 IV  (1150-25000 A)
          "ukf0i"             F0 I     (1150-25000 A)
          "ukf0ii"            F0 II    (1150-25000 A)
          "ukf0iii"           F0 III   (1150-25000 A)
          "ukf0v"             F0 V     (1150-25000 A)
          "ukf2ii"            F2 II    (1150-25000 A)
          "ukf2iii"           F2 III   (1150-25000 A)
          "ukf2v"             F2 V     (1150-25000 A)
          "ukf5i"             F5 I     (1150-25000 A)
          "ukf5iii"           F5 III   (1150-25000 A)
          "ukf5iv"            F5 IV    (1150-25000 A)
          "ukf5v"             F5 V     (1150-25000 A)
          "ukf6v"             F6 V     (1150-25000 A)
          "ukf8i"             F8 I     (1150-25000 A)
          "ukf8iv"            F8 IV    (1150-25000 A)
          "ukf8v"             F8 V     (1150-25000 A)
          "ukg0i"             G0 I     (1150-25000 A)
          "ukg0iii"           G0 III   (1150-25000 A)
          "ukg0iv"            G0 IV    (1150-25000 A)
          "ukg0v"             G0 V     (1150-25000 A)
          "ukg2i"             G2 I     (1150-25000 A)
          "ukg2iv"            G2 IV    (1150-25000 A)
          "ukg2v"             G2 V     (1150-25000 A)
          "ukg5i"             G5 I     (1150-25000 A)
          "ukg5ii"            G5 II    (1150-25000 A)
          "ukg5iii"           G5 III   (1150-25000 A)
          "ukg5iv"            G5 IV    (1150-25000 A)
          "ukg5v"             G5 V     (1150-25000 A)
          "ukg8i"             G8 I     (1150-25000 A)
          "ukg8iii"           G8 III   (1150-25000 A)
          "ukg8iv"            G8 IV    (1150-25000 A)
          "ukg8v"             G8 V     (1150-25000 A)
          "ukk01ii"           K0-1 II  (1150-25000 A)
          "ukk0iii"           K0 III   (1150-25000 A)
          "ukk0iv"            K0 IV    (1150-25000 A)
          "ukk0v"             K0 V     (1150-25000 A)
          "ukk1iii"           K1 III   (1150-25000 A)
          "ukk1iv"            K1 IV    (1150-25000 A)
          "ukk2i"             K2 I     (1150-25000 A)
          "ukk2iii"           K2 III   (1150-25000 A)
          "ukk2v"             K2 V     (1150-25000 A)
          "ukk34ii"           K3-4 II  (1150-25000 A)
          "ukk3i"             K3 I     (1150-25000 A)
          "ukk3iii"           K3 III   (1150-25000 A)
          "ukk3iv"            K3 IV    (1150-25000 A)
          "ukk3v"             K3 V     (1150-25000 A)
          "ukk4i"             K4 I     (1150-25000 A)
          "ukk4iii"           K4 III   (1150-25000 A)
          "ukk4v"             K4 V     (1150-25000 A)
          "ukk5iii"           K5 III   (1150-25000 A)
          "ukk5v"             K5 V     (1150-25000 A)
          "ukk7v"             K7 V     (1150-25000 A)
          "ukm0iii"           M0 III   (1150-25000 A)
          "ukm0v"             M0 V     (1150-25000 A)
          "ukm10iii"          M10 III  (1150-25000 A)
          "ukm1iii"           M1 III   (1150-25000 A)
          "ukm1v"             M1 V     (1150-25000 A)
          "ukm2i"             M2 I     (1150-25000 A)
          "ukm2iii"           M2 III   (1150-25000 A)
          "ukm2p5v"           M2.5 V   (1150-25000 A)
          "ukm2v"             M2 V     (1150-25000 A)
          "ukm3ii"            M3 II    (1150-25000 A)
          "ukm3iii"           M3 III   (1150-25000 A)
          "ukm3v"             M3 V     (1150-25000 A)
          "ukm4iii"           M4 III   (1150-25000 A)
          "ukm4v"             M4 V     (1150-25000 A)
          "ukm5iii"           M5 III   (1150-25000 A)
          "ukm5v"             M5 V     (1150-25000 A)
          "ukm6iii"           M6 III   (1150-25000 A)
          "ukm6v"             M6 V     (1150-25000 A)
          "ukm7iii"           M7 III   (1150-25000 A)
          "ukm8iii"           M8 III   (1150-25000 A)
          "ukm9iii"           M9 III   (1150-25000 A)
          "uko5v"             O5 V     (1150-25000 A)
          "uko8iii"           O8 III   (1150-25000 A)
          "uko9v"             O9 V     (1150-25000 A)
          "ukrf6v"            metal-rich F6 V    (1150-25000 A)
          "ukrf8v"            metal-rich F8 V    (1150-25000 A)
          "ukrg0v"            metal-rich G0 V    (1150-25000 A)
          "ukrg5iii"          metal-rich G5 III  (1150-25000 A)
          "ukrg5v"            metal-rich G5 V    (1150-25000 A)
          "ukrk0iii"          metal-rich K0 III  (1150-25000 A)
          "ukrk0v"            metal-rich K0 V    (1150-25000 A)
          "ukrk1iii"          metal-rich K1 III  (1150-25000 A)
          "ukrk2iii"          metal-rich K2 III  (1150-25000 A)
          "ukrk3iii"          metal-rich K3 III  (1150-25000 A)
          "ukrk4iii"          metal-rich K4 III  (1150-25000 A)
          "ukrk5iii"          metal-rich K5 III  (1150-25000 A)
          "ukwf5v"            metal-weak F5 V    (1150-25000 A)
          "ukwf8v"            metal-weak F8 V    (1150-25000 A)
          "ukwg0v"            metal-weak G0 V    (1150-25000 A)
          "ukwg5iii"          metal-weak G5 III  (1150-25000 A)
          "ukwg5v"            metal-weak G5 V    (1150-25000 A)
          "ukwg8iii"          metal-weak G8 III  (1150-25000 A)
          "ukwk0iii"          metal-weak K0 III  (1150-25000 A)
          "ukwk1iii"          metal-weak K1 III  (1150-25000 A)
          "ukwk2iii"          metal-weak K2 III  (1150-25000 A)
          "ukwk3iii"          metal-weak K3 III  (1150-25000 A)
          "ukwk4iii"          metal-weak K4 III  (1150-25000 A)
          "wf5v"              metal-weak F5 V    (1150-10620 A)
          "wf8v"              metal-weak F8 V    (1150-10620 A)
          "wg0v"              metal-weak G0 V    (1150-10620 A)
          "wg5iii"            metal-weak G5 III  (1150-10620 A)
          "wg5v"              metal-weak G5 V    (1150-10620 A)
          "wg8iii"            metal-weak G8 III  (1150-10620 A)
          "wk0iii"            metal-weak K0 III  (1150-10620 A)
          "wk1iii"            metal-weak K1 III  (1150-10620 A)
          "wk2iii"            metal-weak K2 III  (1150-10620 A)
          "wk3iii"            metal-weak K3 III  (1150-10620 A)
          "wk4iii"            metal-weak K4 III  (1150-10620 A)
        ```
        """
        #
        # Check inputs
        #
        self._check_existing_spectrum(overwrite, quiet=quiet)
        valid_spectral_classes = [
            "a0i",
            "a0iii",
            "a0iv",
            "a0v",
            "a2i",
            "a2v",
            "a3iii",
            "a3v",
            "a47iv",
            "a5iii",
            "a5v",
            "a7iii",
            "a7v",
            "b0i",
            "b0v",
            "b12iii",
            "b1i",
            "b1v",
            "b2ii",
            "b2iv",
            "b3i",
            "b3iii",
            "b3v",
            "b57v",
            "b5i",
            "b5ii",
            "b5iii",
            "b6iv",
            "b8i",
            "b8v",
            "b9iii",
            "b9v",
            "f02iv",
            "f0i",
            "f0ii",
            "f0iii",
            "f0v",
            "f2ii",
            "f2iii",
            "f2v",
            "f5i",
            "f5iii",
            "f5iv",
            "f5v",
            "f6v",
            "f8i",
            "f8iv",
            "f8v",
            "g0i",
            "g0iii",
            "g0iv",
            "g0v",
            "g2i",
            "g2iv",
            "g2v",
            "g5i",
            "g5ii",
            "g5iii",
            "g5iv",
            "g5v",
            "g8i",
            "g8iii",
            "g8iv",
            "g8v",
            "k01ii",
            "k0iii",
            "k0iv",
            "k0v",
            "k1iii",
            "k1iv",
            "k2i",
            "k2iii",
            "k2v",
            "k34ii",
            "k3i",
            "k3iii",
            "k3iv",
            "k3v",
            "k4i",
            "k4iii",
            "k4v",
            "k5iii",
            "k5v",
            "k7v",
            "m0iii",
            "m0v",
            "m10iii",
            "m1iii",
            "m1v",
            "m2i",
            "m2iii",
            "m2p5v",
            "m2v",
            "m3ii",
            "m3iii",
            "m3v",
            "m4iii",
            "m4v",
            "m5iii",
            "m5v",
            "m6iii",
            "m6v",
            "m7iii",
            "m8iii",
            "m9iii",
            "o5v",
            "o8iii",
            "o9v",
            "rf6v",
            "rf8v",
            "rg0v",
            "rg5iii",
            "rg5v",
            "rk0iii",
            "rk0v",
            "rk1iii",
            "rk2iii",
            "rk3iii",
            "rk4iii",
            "rk5iii",
            "uka0i",
            "uka0iii",
            "uka0iv",
            "uka0v",
            "uka2i",
            "uka2v",
            "uka3iii",
            "uka3v",
            "uka47iv",
            "uka5iii",
            "uka5v",
            "uka7iii",
            "uka7v",
            "ukb0i",
            "ukb0v",
            "ukb12iii",
            "ukb1i",
            "ukb1v",
            "ukb2ii",
            "ukb2iv",
            "ukb3i",
            "ukb3iii",
            "ukb3v",
            "ukb57v",
            "ukb5i",
            "ukb5ii",
            "ukb5iii",
            "ukb6iv",
            "ukb8i",
            "ukb8v",
            "ukb9iii",
            "ukb9v",
            "ukf02iv",
            "ukf0i",
            "ukf0ii",
            "ukf0iii",
            "ukf0v",
            "ukf2ii",
            "ukf2iii",
            "ukf2v",
            "ukf5i",
            "ukf5iii",
            "ukf5iv",
            "ukf5v",
            "ukf6v",
            "ukf8i",
            "ukf8iv",
            "ukf8v",
            "ukg0i",
            "ukg0iii",
            "ukg0iv",
            "ukg0v",
            "ukg2i",
            "ukg2iv",
            "ukg2v",
            "ukg5i",
            "ukg5ii",
            "ukg5iii",
            "ukg5iv",
            "ukg5v",
            "ukg8i",
            "ukg8iii",
            "ukg8iv",
            "ukg8v",
            "ukk01ii",
            "ukk0iii",
            "ukk0iv",
            "ukk0v",
            "ukk1iii",
            "ukk1iv",
            "ukk2i",
            "ukk2iii",
            "ukk2v",
            "ukk34ii",
            "ukk3i",
            "ukk3iii",
            "ukk3iv",
            "ukk3v",
            "ukk4i",
            "ukk4iii",
            "ukk4v",
            "ukk5iii",
            "ukk5v",
            "ukk7v",
            "ukm0iii",
            "ukm0v",
            "ukm10iii",
            "ukm1iii",
            "ukm1v",
            "ukm2i",
            "ukm2iii",
            "ukm2p5v",
            "ukm2v",
            "ukm3ii",
            "ukm3iii",
            "ukm3v",
            "ukm4iii",
            "ukm4v",
            "ukm5iii",
            "ukm5v",
            "ukm6iii",
            "ukm6v",
            "ukm7iii",
            "ukm8iii",
            "ukm9iii",
            "uko5v",
            "uko8iii",
            "uko9v",
            "ukrf6v",
            "ukrf8v",
            "ukrg0v",
            "ukrg5iii",
            "ukrg5v",
            "ukrk0iii",
            "ukrk0v",
            "ukrk1iii",
            "ukrk2iii",
            "ukrk3iii",
            "ukrk4iii",
            "ukrk5iii",
            "ukwf5v",
            "ukwf8v",
            "ukwg0v",
            "ukwg5iii",
            "ukwg5v",
            "ukwg8iii",
            "ukwk0iii",
            "ukwk1iii",
            "ukwk2iii",
            "ukwk3iii",
            "ukwk4iii",
            "wf5v",
            "wf8v",
            "wg0v",
            "wg5iii",
            "wg5v",
            "wg8iii",
            "wk0iii",
            "wk1iii",
            "wk2iii",
            "wk3iii",
            "wk4iii",
        ]
        if spectral_class not in valid_spectral_classes:
            raise ValueError(f"{spectral_class} is not a valid `spectral_class`.")
        try:
            data = pd.read_fwf(
                join(DATAPATH, "pickles_spectra", "dat", f"{spectral_class}.dat"),
                colspecs=[(0, 7), (7, 17)],
                header=None,
            )
            self.wavelengths = data[0].values * u.AA
            self.spectrum = data[1].values
        except Exception:
            raise RuntimeError(
                "Could not load the Pickles data for some reason (probably a formatting "
                + "quirk in the file). Please contact the developer with a minimal "
                + "working example."
            )

    def show_spectrum(self, plot=True):
        """
        Plot the spectrum (which should be in units of flam).

        Parameters
        ----------
          plot :: bool
            If True, plot the source weights and return None. If False, return the figure
            and axis instance associated with the plot.

        Returns
        -------
          None (if plot is True)

          fig, ax (if plot is False) :: `matplotlib.figure.Figure`, `matplotlib.axes.Axes`
            The figure and axis instance associated with the plot.
        """
        if self.spectrum is None or self.wavelengths is None:
            raise ValueError("Please generate a spectrum before plotting.")
        fig, ax = plt.subplots()
        ax.plot(self.wavelengths.to(u.AA).value, self.spectrum, "k", lw=1)
        ax.fill_between(self.wavelengths.to(u.AA).value, self.spectrum, alpha=0.5)
        if plt.rcParams["text.usetex"]:
            ax.set_xlabel("Wavelength [\AA]")
            ax.set_ylabel(r"Flux Density [$\rm erg\, s^{-1}\, cm^{-2}\,$\AA$^{-1}$]")
        else:
            ax.set_xlabel(r"Wavelength [$\rm \AA$]")
            ax.set_ylabel(r"Flux Density [$\rm erg\, s^{-1}\, cm^{-2}\,\AA^{-1}$]")
        ax.set_ylim(bottom=0)
        if plot:
            plt.show()
        else:
            return fig, ax

    def calc_redleak_frac(self, TelescopeObj, quiet=False):
        """
        Calculate a source's red leak fraction. The red leak fraction is defined to be the
        ratio of the electron rate (i.e., electron/s) induced by red leak flux to the
        total electron rate induced by the entire spectrum.

        Parameters
        ----------
          TelescopeObj :: `castor_etc.Telescope` instance
            The `Telescope` object containing the passband response curves to use for the
            red leak calculation.

          quiet :: bool
            If True, suppress warnings from red leak fraction calculations.

        Returns
        -------
          redleak_fracs :: dict of floats
            Dictionary containing the red leak fraction in each passband.
        """
        from .sources import CustomSource  # avoid circular import error

        if isinstance(self, CustomSource):
            raise AttributeError("Custom sources do not have red leak fractions!")
        if self.wavelengths is None or self.spectrum is None:
            raise ValueError("Please generate or load a spectrum first")
        #
        # Calculate red leak fraction (red leak electron/s to total electron/s)
        #
        redleak_fracs = dict.fromkeys(TelescopeObj.passbands, 0.0)
        #
        # Make useful source spectrum-derived quantities
        #
        source_wavelengths_AA = self.wavelengths.to(u.AA).value
        source_photon_s_A = (  # photon/s/A
            self.spectrum  # erg/s/cm^2/A
            * TelescopeObj.mirror_area.to(u.cm ** 2).value  # cm^2
            / calc_photon_energy(wavelength=source_wavelengths_AA)[0]  # photon/erg
        )
        source_interp = interp1d(
            source_wavelengths_AA,
            source_photon_s_A,
            kind="linear",
            bounds_error=False,
            fill_value=np.nan,
        )  # photon/s/A
        #
        # Find red leak fraction per band
        #
        for band in redleak_fracs:
            full_response_curve_wavelengths_AA = (
                TelescopeObj.full_passband_curves[band]["wavelength"].to(u.AA).value
            )
            is_redleak = (
                full_response_curve_wavelengths_AA
                > TelescopeObj.redleak_thresholds[band].to(u.AA).value
            )
            redleak_wavelengths = full_response_curve_wavelengths_AA[is_redleak]
            redleak_per_A = (
                source_interp(redleak_wavelengths)
                * TelescopeObj.full_passband_curves[band]["response"][is_redleak]
            )  # electron/s/A
            total_erate_per_A = (
                source_interp(full_response_curve_wavelengths_AA)
                * TelescopeObj.full_passband_curves[band]["response"]
            )  # electron/s/A
            isgood_redleak = np.isfinite(redleak_per_A)  # don't include NaNs
            isgood_total = np.isfinite(total_erate_per_A)  # don't include NaNs
            if not quiet and (not np.all(isgood_redleak) or not np.all(isgood_total)):
                warnings.warn(
                    "Could not estimate red leak fraction "
                    + f"at 1 or more wavelengths in {band}-band. "
                    + "This may just be caused by the source spectrum not being "
                    + f"defined at all wavelengths present in the {band}-band definition "
                    + "file (which runs from "
                    + f"{round(min(full_response_curve_wavelengths_AA), 2)} A "
                    + f"to {round(max(full_response_curve_wavelengths_AA), 2)} A)."
                    + "and is typically not a reason to worry. "
                    + "This warning can be suppressed with `quiet=True`.",
                    RuntimeWarning,
                )
            try:
                redleak_per_px = simpson(  # electron/s (per px)
                    y=redleak_per_A[isgood_redleak],
                    x=redleak_wavelengths[isgood_redleak],
                    even="avg",
                )
                total_erate_per_px = simpson(  # electron/s (per px)
                    y=total_erate_per_A[isgood_total],
                    x=full_response_curve_wavelengths_AA[isgood_total],
                    even="avg",
                )
                redleak_frac = redleak_per_px / total_erate_per_px
            except Exception:
                raise RuntimeError(
                    f"Unable to calculate red leak fraction for {band}-band! "
                    + "Please ensure there is at least 1 wavelength that is above "
                    + "the red leak threshold."
                )
            if np.isfinite(redleak_frac):
                if redleak_frac > 1:
                    # Catch errors caused by very small electron/s rates
                    # (e.g., single emission line spectrum)
                    redleak_fracs[band] = 1.0
                else:
                    redleak_fracs[band] = redleak_frac
            elif not quiet:
                warnings.warn(
                    "Source red leak fraction could not be calculated "
                    + f"in {band}-band!",
                    RuntimeWarning,
                )
        return redleak_fracs


class NormMixin:
    """
    Mixin for normalizing spectra.

    TODO: Make normalization to a specific value at a given wavelength?
    """

    @staticmethod
    def norm_to_star(spectrum, radius=1, dist=1 << u.kpc):
        """
        Normalize the blackbody spectrum to a star of given radius and distance. By
        default, normalizes the flux to a star of 1 solar radius at 1 kpc. Reference:
        <https://github.com/spacetelescope/pysynphot/blob/925cdbac35a7851cee1bddaa2b47651235c44851/pysynphot/spectrum.py#L40>.

        Note that the spectrum's units should have a unit of steradian (sr) in the
        denominator, which will be multiplied out (e.g., erg/s/cm^2/A/sr).

        Parameters
        ----------
          spectrum :: array of floats
            The spectrum to normalize. It should have a unit of steradian in the
            denominator.

          radius :: float or `astropy.Quantity` length
            The radius of the source. If a scalar, it is assumed to be in units of solar
            radii.

          dist :: float or `astropy.Quantity` length
            The distance to the blackbody. If a scalar, it is assumed to be in units of
            kpc.

        Returns
        ----------
          norm_spectrum :: array of floats
            Spectrum normalized such that the unit of steradian in the denominator
            vanishes.
        """
        #
        # Check inputs
        #
        if not isinstance(radius, u.Quantity):
            radius = radius * const.SUN_RADIUS  # cm
        if not isinstance(dist, u.Quantity):
            dist = dist * u.kpc
        radius = radius.to(u.km).value
        dist = dist.to(u.km).value
        #
        # Normalize spectrum
        #
        radian = radius / dist
        return spectrum * np.pi * radian * radian  # multiply by projected solid angle

    def norm_to_AB_mag(
        self,
        ab_mag,
        passband=None,
        TelescopeObj=None,
    ):
        """
        Normalize a spectrum to a given AB magnitude in a specified passband or to a given
        bolometric AB magnitude. The spectrum should be in units of flam (erg/s/cm^2/A).
        The bolometric magnitude calculation assumes a perfect (unity) passband response.

        WARNING: if the spectrum does not vanish at the edges (e.g., a uniform spectrum),
        then the bolometric magnitude will depend on the length of the spectrum! This is
        because the area under the curve does not converge!

        Parameters
        ----------
          ab_mag :: int or float
            The desired AB magnitude.

          passband :: valid `Telescope` passband string (e.g., "uv", "u", "g") or None
            If not None, normalize the spectrum such that it has the desired AB magnitude
            in this passband; `TelescopeObj` must also be provided. If None, normalize the
            spectrum such that the spectrum's bolometric magnitude equals the given AB
            magnitude; `TelescopeObj` must not be provided.

          TelescopeObj :: `Telescope` object or None
            The `Telescope` object containing the limits and response curves of the
            different passbands. If not None, `passband` must also be provided. If None,
            `passband` must not be provided.

        Attributes
        ----------
          spectrum :: array of floats
            The renormalized spectrum in units of erg/s/cm^2/A (identical to the previous
            spectrum units).

        Returns
        -------
          None
        """
        #
        # Check inputs
        #
        if self.spectrum is None or self.wavelengths is None:
            raise ValueError("Please generate a spectrum before normalizing.")
        if (TelescopeObj is None and passband is not None) or (
            TelescopeObj is not None and passband is None
        ):
            raise ValueError(
                "Please either specify both `TelescopeObj` and `passband`"
                + "or neither of them."
            )
        #
        # Normalize
        #
        current_ab_mag = self.get_AB_mag(TelescopeObj=TelescopeObj)
        if passband is not None:
            if passband not in TelescopeObj.passbands:
                raise ValueError(
                    f"Invalid `passband`. Valid passbands are: {TelescopeObj.passbands}."
                )
            current_ab_mag = current_ab_mag[passband]
        norm_factor = 10 ** (-0.4 * (ab_mag - current_ab_mag))
        self.spectrum *= norm_factor

    def norm_luminosity_dist(self, luminosity, dist):
        """
        Normalize the spectrum to a source of given (bolometric) luminosity and distance.
        The `Source` object should have its spectrum in units of flam (erg/s/cm^2/A) and
        wavelengths in angstrom. (Technically it is okay as long as the wavelengths are in
        some unit <U> and the spectrum is in units of erg/s/cm^2/<U>.)

        Parameters
        ----------
          luminosity :: scalar or `astropy.Quantity` unit of power
            The desired bolometric luminosity. If a scalar, this is assumed to be in units
            of solar luminosities. If an `astropy.Quantity`, it must be a unit of power
            (e.g., erg/s).

          dist :: float or `astropy.Quantity` length
            The distance to the source. If a scalar, it is assumed to be in units of kpc.

        Attributes
        ----------
          spectrum :: array
            Normalized spectrum in original flux density units (e.g., erg/s/cm^2/A).

        Returns
        -------
          None
        """
        #
        # Check inputs
        #
        if self.spectrum is None or self.wavelengths is None:
            raise ValueError("Please generate a spectrum before normalizing.")
        if not isinstance(luminosity, u.Quantity):
            luminosity = luminosity * const.SUN_LUMINOSITY
        if not isinstance(dist, u.Quantity):
            dist = dist * u.kpc
        luminosity = luminosity.to(u.erg / u.s).value
        dist = dist.to(u.cm).value
        #
        # Normalize spectrum (originally in erg/s/cm^2/A)
        #
        erg_s_A = 4 * np.pi * dist * dist * self.spectrum  # erg/s/A
        tot_luminosity = simpson(
            y=erg_s_A, x=self.wavelengths.to(u.AA).value, even="avg"
        )  # erg/s
        norm_factor = luminosity / tot_luminosity  # dimensionless
        self.spectrum *= norm_factor  # erg/s/cm^2/A

    def get_AB_mag(self, TelescopeObj=None):
        """
        Calculate the AB magnitude of the source either through the telescope's passbands
        or over the whole spectrum (bolometric magnitude). Note that the source spectrum
        should be in units of erg/s/cm^2/A. The bolometric magnitude calculation assumes a
        perfect (unity) passband response.

        WARNING: if the spectrum does not vanish at the edges (e.g., a uniform
        spectrum), then the bolometric magnitude will depend on the length of the
        spectrum! This is because the area under the curve does not converge!

        Parameters
        ----------
          TelescopeObj :: `Telescope` object or None
            If provided, will calculate the AB magnitude of the source in each of the
            telescope's passbands. Otherwise will calculate the source's bolometric AB
            magnitude.

        Returns
        -------
          ab_mags :: scalar or dict of scalars
            If `TelescopeObj` is None, the result is a scalar equal to the bolometric AB
            magnitude of the source. If TelescopeObj is not None, the result is a dict of
            scalars representing the source's AB magnitude through each of the telescope's
            passbands.
        """
        if self.spectrum is None or self.wavelengths is None:
            raise ValueError(
                "Please generate a spectrum before calculating AB magnitude(s)."
            )
        wavelengths_AA = self.wavelengths.to(u.AA).value
        if TelescopeObj is None:
            #
            # Calculate bolometric AB magnitude
            #
            return flam_to_AB_mag(
                wavelengths_AA,
                self.spectrum,
                np.ones_like(wavelengths_AA, dtype=float),  # perfect "passband" response
            )
        else:
            #
            # Calculate the AB magnitude through each of the telescope's passbands
            #
            if not isinstance(TelescopeObj, Telescope):
                raise TypeError(
                    "`TelescopeObj` must be a `castor_etc.telescope.Telescope` object."
                )
            ab_mags = dict.fromkeys(TelescopeObj.passbands)
            spectrum_interp = interp1d(
                x=wavelengths_AA,
                y=self.spectrum,
                kind="linear",
                bounds_error=False,
                fill_value=np.nan,
            )
            for band in TelescopeObj.passbands:
                passband_wavelengths = (
                    TelescopeObj.full_passband_curves[band]["wavelength"].to(u.AA).value
                )
                passband_spectrum = spectrum_interp(passband_wavelengths)
                isgood_passband = np.isfinite(passband_spectrum)  # do not integrate NaNs
                if np.any(~isgood_passband):
                    if np.all(~isgood_passband):
                        raise RuntimeError(
                            "Source spectrum could not be estimated "
                            + f"at any {band}-band wavelength"
                        )
                    else:
                        warnings.warn(
                            "Source spectrum could not be estimated "
                            + f"at some {band}-band wavelengths",
                            RuntimeWarning,
                        )
                ab_mags[band] = flam_to_AB_mag(
                    passband_wavelengths[isgood_passband],
                    passband_spectrum[isgood_passband],
                    TelescopeObj.full_passband_curves[band]["response"][isgood_passband],
                )
            return ab_mags
