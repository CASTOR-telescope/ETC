# Changelog

See [semantic versioning](https://semver.org/spec/v2.0.0.html) for the rationale behind
the version numbers.

## Current version (as of 2025-02-24)

- Fixed a bug ([#31](https://github.com/CASTOR-telescope/ETC/issues/31)) in the definition
  of `ellipticity`/`eccentricity` for a `GalaxySource`

## [1.3.2](https://github.com/CASTOR-telescope/ETC/tree/v1.3.2) (2024-11-10)

- Fixed the broken JupyterLab terminal on CANFAR
- Removed the `even` parameter in all `scipy.integrate.simpson()` calls, since this
  parameter is removed in the newest versions of scipy.

## [1.3.1](https://github.com/CASTOR-telescope/ETC/tree/v1.3.1) (2024-01-02)

- Added an `_encircled_energies` attribute to the `Photometry` object. This is needed for
  v1.1.1 of the [ETC frontend](https://github.com/CASTOR-telescope/ETC_frontend)
  ([#26](https://github.com/CASTOR-telescope/ETC/pull/26)).

## [1.3.0](https://github.com/CASTOR-telescope/ETC/tree/v1.3.0) (2023-11-13)

- Change `MP` and `READ_NOISE` parameters in parameters.py to 930 and 3.0, respectively.
  Also add some comments about parameters
  ([#22](https://github.com/CASTOR-telescope/ETC/pull/22))
- Introduced a new feature of Grism spectroscopy. This is very much still in development
  and currently only supports Sersic galaxy sources. Available results are: 2D dispersed
  spectrum and 1D SNR dispersed spectrum.
- Adapted transit numerical simulation from the POET mission. Using CASTOR tools, the
  light curve simulation calculation has been reduced from ~ 4 minutes to 1 second.
- Added a new feature that allows users to query GaiaDR2 database on the basis of
  (ra,dec), and threshold Gaia G magnitude, and select an appropriate normalized spectrum.
  This is achievable by selecting `Gaia` from the Predefined Spectra drop-down menu for a
  Point Source. To do transit simulation, one has to choose a Gaia Point Source.
- Added a new function `show_slit_image` in the uvmos_spectroscopy.py file to display a
  cropped view of the source on the detector.

## [1.2.1](https://github.com/CASTOR-telescope/ETC/tree/v1.2.1) (2023-09-18)

- Add ability to control the colour bar normalization when visualizing PSFs
  ([#16](https://github.com/CASTOR-telescope/ETC/pull/16))
- Update the bundled CASTOR PSFs ([#17](https://github.com/CASTOR-telescope/ETC/pull/17))

## [1.2.0](https://github.com/CASTOR-telescope/ETC/tree/v1.2.0) (2023-08-25)

- Changed AB magnitude calculation to interpolate the passband to the spectrum resolution
  ([#7](https://github.com/CASTOR-telescope/ETC/pull/7))
- Added PSF convolution to all sources using CASTOR PSFs
  ([#11](https://github.com/CASTOR-telescope/ETC/pull/11))
  - The fraction of flux enclosed within the aperture is estimated from the supersampled
    simulated images before binning down to the telescope's resolution. This works for any
    aperture with that is arbitrarily centred on the source. In this way, this change
    fixes the normalization limitation from v1.1.1, whose reference flux assumed galaxy
    sources were centred on and unrotated relative to the aperture
  - Now all photometry calculations use the fraction of flux enclosed within the aperture
    to calculate the electrons produced per second on the detector due to the source. Note
    that the reference flux, which determines what we consider to be 100% of the flux, is
    always based on the _noiseless_ image, i.e., with no PSF convolution
    - An enclosed flux fraction of 100% corresponds to the magnitude of the source that
      the user set. So if the user normalizes a spectrum to, say, an AB magnitude of 25,
      then this AB magnitude will be the AB magnitude of the source if 100% of its flux
      was contained within the aperture. If the user selects an aperture that only
      contains 50% of the flux, however, the effective AB magnitude will be dimmer
  - The photometry ETC is now 100% pixel-based, even for point sources. Changing any pixel
    in the `source_weights` attribute, which gives the fraction of the total flux
    contained within the pixel, _will_ change the photometry calculations

## [1.1.1](https://github.com/CASTOR-telescope/ETC/tree/v1.1.1) (2022-09-09)

- Jennifer Glover added UVMOS capabilities to ETC
- Tyrone Woods fixed a bug (noticed by Michael Balogh?) in `photometry.py` re: high Sérsic
  index galaxies

## [1.0.0](https://github.com/CASTOR-telescope/ETC/tree/v1.0.0) (2022-05-11)

Initial release!
