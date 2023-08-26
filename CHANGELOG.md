# Changelog

See [semantic versioning](https://semver.org/spec/v2.0.0.html) for the rationale behind
the version numbers.

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
