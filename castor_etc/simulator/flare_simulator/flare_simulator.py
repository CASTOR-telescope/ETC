
# TODO: refactor this code into a class if possible

# Function created with use of fiducial_flare package

import numpy as np
from astropy import table
from astropy import units as u
from matplotlib import pyplot as plt

from .fiducial_flare import *

t = 100 * u.s
plt.ion()


def compute_flux_wavelength(star_radius, star_eff_temp, dist_to_star, show_figure):
    """ 
    Calculates flux from input values, picks corresponding spectra for M0-M9.
    Adds corresponding correction factors for chosen MUSCLES model.
    Generate time-evolving spectra from a random series of flares.

    Parameters
    ----------
    {Star_radius} in solar radius
    {star_eff_temp} in Kelvin
    {dist_to_star} in light-years
    {print_figure} prints the full spectrum, 
    Returns
    -------
    spectra : 2D array
        Array of spectra in each tbin, where the array has dimensions
        (len(tbins)-1, len(wbins)-1). Units will match the product of
        the eqd and SiIV_quiescent units, divided by time and length.
        Returns wavelength, flux
        
    """

    #Units (W / m2) converted later
    R_sun = 6.95700e8  # Mean Radius in meters
    sigma_sb = 5.67037e-8 # Units W/m^2/k^4
    ly = 9.46073e15 #Units meters

    Calculated_flux = (((R_sun*star_radius)**2)*sigma_sb*(star_eff_temp**4))/((ly*dist_to_star)**2)

    # Depending on star temperature it will choose which flare profile has the closest match
    # Files correspond to Main Sequence M star flare spectra
    # Files are from the following link for MUSCLES spectra for M0-M5.5 'https://archive.stsci.edu/missions/hlsp/muscles/'
    # Files taken in Feb 2026

    # Partially Convective
    # M0-M1
    if star_eff_temp >= 3600:
        path_sed = 'ETC/castor_etc/data/flare_simulator_data/M1.5V_hlsp_muscles_multi_multi_gj667c_broadband_v22_adapt-const-res-sed.fits'

    # M2
    elif star_eff_temp < 3600 and star_eff_temp >= 3500:
        path_sed = 'ETC/castor_etc/data/flare_simulator_data/M2V_hlsp_muscles_multi_multi_gj176_broadband_v22_adapt-const-res-sed.fits'

    # M3
    elif star_eff_temp < 3500 and star_eff_temp >= 3400:
        path_sed = 'ETC/castor_etc/data/flare_simulator_data/M3V_hlsp_muscles_multi_multi_gj581_broadband_v22_adapt-const-res-sed.fits'

    # Fully Convective
    # Fully Convective produce stronger more energetic flares
    # M4
    elif star_eff_temp < 3400 and star_eff_temp >= 3200:

        path_sed = 'ETC/castor_etc/data/flare_simulator_data/M4_hlsp_muscles_multi_multi_gj876_broadband_v22_adapt-const-res-sed.fits'

    # M5.5
    else:
        path_sed = 'ETC/castor_etc/data/flare_simulator_data/M5.5_hlsp_muscles_multi_multi_gj551_broadband_v22_adapt-const-res-sed.fits'

    # No MUSCLES model available for M stars M6-M9, default will be M5.5


    # Spectral Energy Distribution (SED), the base quiescent SED is from the MSUCLES spectrum different for each star class
    sed = table.Table.read(path_sed, hdu=1, unit_parse_strict="silent")

    # Next part adds correction factors based on research paper "Optically Quiet, But FUV Loud:" cited below:
    # Paper used archival FUV observations of low-mass stars to test the UV predictions of literature flare models
    # Flare model used in this funtion is the MUSCLES model
    # MUSCLES -> (Measurements of the Ultraviolet Spectral Characteristics of Low-mass Exoplanetary Systems) Treasury Survey
    # Wavelength values below were taken from Table 1 in "Optically Quiet, But FUV Loud:" paper
    # Corretion factors used was from Table 4 and Table 5 for MUSCLES Model in "Optically Quiet, But FUV Loud:" paper
    # Corretion factors were split based on if the star was partially or fully convective

    #### Cited paper: ####
    # Jackman, J. A. G., Shkolnik, E. L., Loyd, R. O. P., and Richey-Yowell, T.,
    # “Optically quiet, but FUV loud: results from comparing the far-ultraviolet predictions of flare models with TESS and HST”,
    # <i>Monthly Notices of the Royal Astronomical Society</i>, vol. 533, no. 2, OUP, pp. 1894–1906, 2024.
    # doi:10.1093/mnras/stae1570.
    # 'https://ui.adsabs.harvard.edu/abs/2024MNRAS.533.1894J/abstract'

    # Paper accessed Feb 2026


    w = sed["WAVELENGTH"]
    f = sed['BOLOFLUX']

    # Changing file format
    w = np.asarray(w)
    f = np.asarray(f)

    if star_eff_temp >= 3400: #Partially Convective Correction Factor Field Age M star M0~M3

    #####  White Light Correction for partially convective ######

        # how this works: 1. condition, 2. if true, 3. if false
        f = np.where((w > 1173.65) & (w < 1198.49), 0.42 * f, f)
        f = np.where((w > 1201.71) & (w < 1206.50), 0.42 * f, f)
        f = np.where((w > 1223.00) & (w < 1273.50), 0.42 * f, f)
        f = np.where((w > 1328.20) & (w < 1354.49), 0.42 * f, f)
        f = np.where((w > 1356.71) & (w < 1357.59), 0.42 * f, f)
        f = np.where((w > 1359.51) & (w < 1362.70), 0.42 * f, f)


    ##### FUV Correction for partially convective #####

        f = np.where((w > 1173.65) & (w < 1198.49), 0.5 * f, f)
        f = np.where((w > 1201.71) & (w < 1212.16), 0.5 * f, f)
        f = np.where((w > 1219.18) & (w < 1274.04), 0.5 * f, f)
        f = np.where((w > 1329.25) & (w < 1354.49), 0.5 * f, f)
        f = np.where((w > 1356.71) & (w < 1357.59), 0.5 * f, f)
        f = np.where((w > 1359.51) & (w < 1428.90), 0.5 * f, f)


    ##### Pseudo-continuum 130 Correction for partially convective #####

        f = np.where((w > 1173.65) & (w < 1174.50), 0.68 * f, f)
        f = np.where((w > 1176.80) & (w < 1190.00), 0.68 * f, f)
        f = np.where((w > 1223.00) & (w < 1238.40), 0.68 * f, f)
        f = np.where((w > 1239.30) & (w < 1242.00), 0.68 * f, f)
        f = np.where((w > 1243.50) & (w < 1273.50), 0.68 * f, f)
        f = np.where((w > 1329.00) & (w < 1334.00), 0.68 * f, f)
        f = np.where((w > 1336.00) & (w < 1354.49), 0.68 * f, f)
        f = np.where((w > 1356.71) & (w < 1357.59), 0.68 * f, f)
        f = np.where((w > 1359.51) & (w < 1362.70), 0.68 * f, f)


     ##### Si IV Correction for partially convective #####

        f = np.where((w > 1393.76) & (w < 1402.77), 0.19 * f, f)

    ##### Si III Correction for partially convective #####

        f = np.where((w == 1206.51), 0.10 * f, f)

    ##### C III Correction for partially convective #####

        f = np.where((w == 1174.93), 0.84 * f, f)
        f = np.where((w == 1175.25), 0.84 * f, f)
        f = np.where((w == 1175.59), 0.84 * f, f)
        f = np.where((w == 1175.71), 0.84 * f, f)
        f = np.where((w == 1175.99), 0.84 * f, f)
        f = np.where((w == 1176.37), 0.84 * f, f)


    # Fully Convective Correction Factor Field Age M Star M4~M9
    else:

    #####  White Light Correction for fully convective ######

        f = np.where((w > 1173.65) & (w < 1198.49), 4.7 * f, f)
        f = np.where((w > 1201.71) & (w < 1206.50), 4.7 * f, f)
        f = np.where((w > 1223.00) & (w < 1273.50), 4.7 * f, f)
        f = np.where((w > 1328.20) & (w < 1354.49), 4.7 * f, f)
        f = np.where((w > 1356.71) & (w < 1357.59), 4.7 * f, f)
        f = np.where((w > 1359.51) & (w < 1362.70), 4.7 * f, f)


    ##### FUV Correction for fully convective #####


        f = np.where((w > 1173.65) & (w < 1198.49), 2.9 * f, f)
        f = np.where((w > 1201.71) & (w < 1212.16), 2.9 * f, f)
        f = np.where((w > 1219.18) & (w < 1274.04), 2.9 * f, f)
        f = np.where((w > 1329.25) & (w < 1354.49), 2.9 * f, f)
        f = np.where((w > 1356.71) & (w < 1357.59), 2.9 * f, f)
        f = np.where((w > 1359.51) & (w < 1428.90), 2.9 * f, f)


    ##### Pseudo-continuum 130 Correction for fully convective #####

        f = np.where((w > 1173.65) & (w < 1174.50), 4.71 * f, f)
        f = np.where((w > 1176.80) & (w < 1190.00), 4.71 * f, f)
        f = np.where((w > 1223.00) & (w < 1238.40), 4.71 * f, f)
        f = np.where((w > 1239.30) & (w < 1242.00), 4.71 * f, f)
        f = np.where((w > 1243.50) & (w < 1273.50), 4.71 * f, f)
        f = np.where((w > 1329.00) & (w < 1334.00), 4.71 * f, f)
        f = np.where((w > 1336.00) & (w < 1354.49), 4.71 * f, f)
        f = np.where((w > 1356.71) & (w < 1357.59), 4.71 * f, f)
        f = np.where((w > 1359.51) & (w < 1362.70), 4.71 * f, f)


    ##### Si IV Correction for fully convective #####

        f = np.where((w > 1393.76) & (w < 1402.77), 3.09 * f, f)

    ##### C II Correction for fully convective #####

        f = np.where((w == 1334.53), 4.70 * f, f)
        f = np.where((w == 1335.71), 4.70 * f, f)

    ##### C III Correction for fully convective #####

        f = np.where((w == 1174.93), 6.41 * f, f)
        f = np.where((w == 1175.25), 6.41 * f, f)
        f = np.where((w == 1175.59), 6.41 * f, f)
        f = np.where((w == 1175.71), 6.41 * f, f)
        f = np.where((w == 1175.99), 6.41 * f, f)
        f = np.where((w == 1176.37), 6.41 * f, f)

    ##### N V Correction for fully convective #####

        f = np.where((w == 1238.82), 8.08 * f, f)
        f = np.where((w == 1242.80), 8.08 * f, f)

    # The rest of the function was modified from fiducial_flare package, refer to top of document for more information

    # The SED has higher resolution needed.
    # In fiducial flare spectrum, lines are represented by single 200 km/s wide bins (about 1 Å in the FUV).
    # Higher resolutions will yield odd-looking output once we add flare spectra to the SED.
    # We'll rebin to 1 Å bins from 100-3000 Å and then 10 Å through the 5.5 µm limit of the SED
    # The original SED also has variable binning. The left and and right edges of the bins are given in the
    #'WAVELENGTH0' and 'WAVELENGTH1' columns.
    # The rebin function is one of the few in fiducial_flare that specifically requires input without units.

    wbins_sed = np.append(sed['WAVELENGTH0'], sed['WAVELENGTH1'][-1]) * u.AA

    wbins_fuv = np.arange(100, 3000, 1) * u.AA
    wbins_red = np.arange(3000, 5.5e4, 10) * u.AA
    wbins = np.hstack((wbins_fuv, wbins_red))

    # Now rebin, f is 'BOLOFLUX'
    Flux_quiescent_bolo = rebin(wbins.value, wbins_sed.value, f)

    # Add the proper units since rebin can't work with them.
    Flux_quiescent_bolo = Flux_quiescent_bolo * u.Unit('AA-1')


    # The flux here is normalized by the bolometric luminosity of the star (units of Å-1).
    # Calculating flux of inputted star using flux calculated at beginning.
    Flux_bolo_your_star = Calculated_flux * u.Unit('W m-2')
    Flux_quiescent = Flux_quiescent_bolo * Flux_bolo_your_star
    Flux_quiescent = Flux_quiescent.to('erg s-1 cm-2 AA-1')


    # Now that we've got the quiescent SED all prepped, let's simulate some flares to go on top of that.
    # Blackbody emission is added to simulate flux emitted at longer wavelengths.
    # Flares are scaled to the quisecent flux of the star in the Si IV 1393,1402 Å emission line doublet.

    # Flares are scaled to this value for SiIV.
    wbin_SiIV = [1390, 1410] * u.AA
    Fq_SiIV = rebin(wbin_SiIV.value, wbins.value, Flux_quiescent.value) * u.Unit('erg s-1 cm-2 AA-1')

    # This actually spits out the flux density, but what we want is the flux.
    Fq_SiIV = Fq_SiIV * np.diff(wbin_SiIV)


    # We need some time bins to simulate flares over. We will use 60 s bins covering a full day.
    # sadly arange can't handle unit input
    tbins = np.arange(0, 24*60*60, 60) * u.s


    # Simulating time series of fluxes from a random series of flares.
    np.random.seed(42)
    Flux_flare = flare_series_spectra(wbins, tbins, SiIV_quiescent=Fq_SiIV)


    # The resulting array has dimensions of (no. time bins) x (no. wavelength bins), in this case 1439x8099.
    # Each row (e.g. Flux_flare[0,:]) is a spectrum for the corresponding time bin.

    # To get the spectrum CASTOR will actually see, we need to add these spectra onto the quiescent SED.
    Flux_tot = Flux_quiescent[None,:] + Flux_flare

    # Highest peak
    imax = np.argmax(Flux_flare[:,0])

    # Plot spectra at the highest peak and at quiescence
    if show_figure == "yes" or show_figure == "Yes":
        plt.figure()

        # Quiescence
        line_quiescence, = plt.step(wbins[:-1], Flux_quiescent, where='pre', label='quiescence')

        # Highest peak
        line_peak, = plt.step(wbins[:-1], Flux_tot[imax,:], where='pre', label='max')


        plt.legend(handles=(line_quiescence, line_peak))
        plt.xlabel('Wavelength (Å)')
        plt.ylabel('Flux (erg s-1 cm-2 Å-1)')
        plt.xlim(0, 5500) #Limit is set to CASTORS wavelengths
        plt.yscale('log')

    # Skips plotting
    else: show_figure == "no" or show_figure == "No"
    # Converts to units CASTOR_UVMOS can use
    wavelength = wbins[:-1].value
    flux = Flux_tot[imax,:].value

    #UVMOS has wavelengths from 1500Å to 5500Å, this is sliced to reflect that.
    #If you have a different wavelength range comment out and return wavelength, flux
    wavelength_sliced = wavelength[1400:3101]
    flux_sliced = flux[1400:3101]

    # Wavlength and flux can be directly inputted into CASTOR_UVMOS to use
    return wavelength_sliced, flux_sliced

    #If not using UVMOS range
    #return wavelength, flux

#endregion
