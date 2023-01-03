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

And the following is from Figure 6.1 of the STIS Instrument Handbook, v20.0
(<https://hst-docs.stsci.edu/stisihb/chapter-6-exposure-time-calculations/6-6-tabular-sky-backgrounds>):
"In the ETCs and in this Handbook, the choices for earthshine of 'shadow,' 'average,' and
'extremely high' correspond to 0, 50% of, and twice the 'high' values in Table 6.4. For
the zodiacal sky background, the values in Table 6.4 correspond to a high value of m_V =
22.1 arcsec^-2 from Table 6.2, while the low and average zodiacal light are scaled to m_V
= 23.3 arcsec^-2 and 22.7 arcsec^-2, respectively."

I've checked that the "high" Earthshine values in Table 6.4 of the STIS Instrument
Handbook closely match those given in `earthshine.fits` (albeit the values in Table 6.4
are less precise and of lower resolution than the data in this folder).

In contrast, the "high" zodiacal light numbers tabulated in Table 6.4 are noticeably
higher than those in `zodi.fits`. Upon thinking about it more, I don't actually think the
`zodi.fits` spectrum corresponds to a particular "high"/"average"/"low" level of zodiacal
light. Instead, I think the `zodi.fits` spectrum is renormalized to match certain
heuristics (also see <https://github.com/CASTOR-telescope/ETC/issues/6>). FYI, I
calculated a bolometric AB mag of 23.39, a V-band AB magnitude of 23.16 using the coarse
Bessel V passband curve available here:
<http://spiff.rit.edu/classes/phys440/lectures/filters/bess-v.pass>, and a V-band AB
magnitude of 23.19 using the Johnson V passband curve available here:
<https://github.com/spacetelescope/pysynphot/blob/master/pysynphot/data/generic/johnson_v.fits>.
So I guess `zodi.fits` is appproximately the "low" zodiacal light estimate (m_V = 23.3).
