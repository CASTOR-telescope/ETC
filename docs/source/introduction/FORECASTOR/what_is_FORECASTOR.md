# FORECASTOR

The first tool developed for the *Finding Optics Requirements and Exposure times for CASTOR (FORECASTOR)* project, is a dedicated pixel-based photometric exposure time calculator and associated web-based
graphical user interface written entirely in Python. The *CASTOR* team has set out to make the package as maintainable as possible while providing a simple user experience. The critical components are separated by object classes. In python, object classes are a type of 'blueprint' to keep track of all of the data about a particular object as well as any functions that we want on use with that data. The user may separately define the `Telescope`, `Background`, and `Source` objects, and modify the parameters as needed. These objects are then passed to the `Photometry` object which returns the desired signal-to-noise calculation. 


## Telescope

When using the `castor_etc` Python package, you will first want to define an instance of a `Telescope` object. This object is completely customizable to your desired characteristics. Some features you may change are pixel scale, read noise, the number and name of passbands, the extinction coefficients and more. If nothing is changed the `Telescope` object is initialized with default values according to the most up-to-date information available at the time. 


*CASTORs* default photometric passbands are UV, u and g. These may be substituted for any arbitrary throughput appropriate to a given detector, optics, and filter combination. The passband response curve, the photometric zero-point and pivot wavelength are calculated automatically. The photometric zero-point (which is defined as the AB magnitude which produces 1 {$e^{-}s^{-1}$} in a given passband) is calculated assuming a flat spectrum in AB magnitude and converges on the zero-point using either the secant or bisection method which can be chosen by the user. The pivot wavelength is computed following Eq. (A11) from Tokunaga & Vacca (2005) and uses Simpson’s rule to approximate the integration.


For each of *CASTORs* passbands, `castor_etc` also includes a default wavelength-averaged PSF sampled at 20× the telescope’s pixel scale (i.e., each pixel has a side length of 0.1/20′′ = 0.005′′). It is theoretically possible for other missions to use this package as long as the telescope detector is CCD or CMOS based or similar since *CASTOR* has such a wide range of customization.


## Background

The second object to initialize is the sky `Background`, this is characterized by three parameters, Earthshine, zodiacal light, and geocoronal emission (i.e. airglow). Users can customize by inputting their own Earthshine and zodiacal light spectra or users can describe the sky background in AB magnitude per square arcsecond in each passband. If none is supplied the default will be spectra from the Hubble Space Telescope. The default airglow value is set to represent [Oii] 2471 Å emission line, which is centred at 2471 Å with a linewidth of 0.023 Å. `castor_etc` provides three predefined flux values of "high", "avg" and "low". Users may add geocoronal emission lines of arbitrary wavelength, linewidth, and flux if the default isn’t to their needs. Finally, the sky background is considered uniform but there is options for spatially non-uniform backgrounds defined pixel-by-pixel.

## Source

The next step is to define a `Source` object for mock observations, there are three steps to follow.

1. Determining the type of the source (a point source, an
extended source like a diffuse nebula, a galaxy, etc.).
2. Describing the physical properties of the source, such as
its spectrum (including any emission/absorption lines),
redshift, distance, surface brightness profile (for extended sources or galaxies), etc.
3. (Optional) Renormalizing the source spectrum. There
are several normalization schemes available: normalize
to an AB magnitude within a passband, normalize to
a total luminosity and distance, normalize a blackbody
spectrum to a star at a given radius and distance. Note
that these normalizations can be applied at any time (e.g.,
can be before or after the addition of spectral lines).

There are currently three source type classes: `PointSource`, `ExtendedSource`, and `GalaxySource` which facilitate the creation of mock point sources, diffuse extended sources, and galaxies, respectively. Users may also upload their own data in either a FITS or ASCII file to create a CustomSource instance where the data will be automatically interpolated to the Telescope’s pixel scale.

Along with these sources, `castor_etc` provides the means to generate the following spectra: blackbody, power-law, emission line with different line shapes, and a flat spectrum as well as use of spectra from the Pickles library and spectra for spiral and elliptical galaxies from Fioc, M., & Rocca-Volmerange, B. 1997. Finally, we support adding emission and absorption lines with various line profiles to any spectra listed above.


# Photometry 


`castor_ect` is capable of obtaining photometric estimates, which is precise measurements of the brightness of celestial objects to get information such as apparent magnitude, temperature, distance, stellar pulsations, exoplanet transits. This is done with the `Photometry` class which is initialized by giving it the `Telescope`, `Source`, and `Background` objects. Photutils (Bradley et al. 2022) is used to generate apertures with fractional pixel contributions for signal-to-noise (S/N) measurements. The contribution of a pixel is directly proportional to how much of the pixel is contained within the aperture. The shape of the aperture can currently be customized for rectangular apertures, elliptical apertures, and the point sources or “optimal” apertures which is a circular aperture centred on the source that maximizes the S/N. The optimal aperture diameter is by default set to 1.4x the telescopes FWHM, this is roughly if the *Point Spread Function* (PSF) is a 2D Gaussian. This multiplicative factor can be changed by the user.


The fraction of flux enclosed within the aperture is also determined by comparing the flux within the aperture to some reference flux value. These enclosed flux fractions are used to determine the number of electrons produced per second on the detector in a given passband. 


1. Point Source: The reference flux is simply the sum of the PSF pixel values, which represents 100% of the flux from the source.
2. Galaxies: The reference flux is the flux contained within a centred elliptical aperture that is the same size as the galaxy’s half-light radius and of the same orientation. We then assume the total flux from the galaxy to be twice this reference flux. 
3. Extended Sources: The reference value representing 100% of the flux for an extended source is the signal obtained through a centred elliptical aperture of the same dimensions and orientation as the extended source. We assume 100% of the flux is contained within the source’s angular extent.

An enclosed flux fraction of 100% corresponds to the magnitude of the source that the user set. If the user normalizes a spectrum to an AB magnitude of 25, then this AB magnitude will be the AB magnitude of the source if 100% of its flux was contained within the aperture. If the user selects an aperture that only contains 50% of the flux, however, the effective AB magnitude will be dimmer and produce half as much signal as an aperture that contains 100% of the flux.
Currently the castor_etc can only simulate single objects from start to finish, users may upload custom images (e.g., a crowded field) and use the ETC’s tools to obtain photometry and spectroscopy (Glover et al. 2022) estimates. As the `castor_etc` software progress, more tools will become available for public use.
