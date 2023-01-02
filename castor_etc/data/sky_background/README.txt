Note that although the FITS headers for `earthshine.fits` and `zodi.fits` say the spectrum
is in units of erg/s/cm^2/A, the actual units are in erg/s/cm^2/A/arcsec^2. It does not
make sense if the sky background values are not in units of surface brightness. What has
likely happened is that these sky background values are given for a specific aperture
(i.e., an aperture of 1 sq. arcsec), hence the header units are in erg/s/cm^2/A.

P.S. The data are from:
<https://github.com/spacetelescope/pysynphot/blob/master/pysynphot/data/generic/Zodi.fits>
and
<https://github.com/spacetelescope/pysynphot/blob/master/pysynphot/data/generic/earthshine.fits>.

The following description is from Table 6.4 of the STIS Instrument Handbook, v20.0
(<https://hst-docs.stsci.edu/stisihb/chapter-6-exposure-time-calculations/6-6-tabular-sky-backgrounds>):
"The high sky values [of Earthshine and zodiacal light] are defined as the earthshine at
38° from the limb and by the high zodiacal light of m_V = 22.1 arcsec^-2."

I've checked that the Earthshine and zodiacal light values in the link above closely match
those given in `earthshine.fits` and `zodi.fits` (albeit the values in Table 6.4 are less
precise and of lower resolution than the data in this folder).
