Note that although the FITS headers for `earthshine.fits` and `zodi.fits` say the spectrum
is in units of erg/s/cm^2/A, the actual units are in erg/s/cm^2/A/arcsec^2. It does not
make sense if the sky background values are not in units of surface brightness. What has
likely happened is that these sky background values are given for a specific aperture
(i.e., an aperture of 1 sq. arcsec), hence the header units are in erg/s/cm^2/A.

P.S. The data are from:
<https://github.com/spacetelescope/pysynphot/blob/master/pysynphot/data/generic/Zodi.fits>
and
<https://github.com/spacetelescope/pysynphot/blob/master/pysynphot/data/generic/earthshine.fits>.
