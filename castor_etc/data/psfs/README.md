# CASTOR PSFs

This folder contains the point spread functions (PSFs) for each of CASTOR's passbands. The
following files are averaged 2D PSFs prepared by Madeline Marshall that are, strictly
speaking, only valid for the center of CASTOR's field of view:

- `NUVsamples_median_withJitter_X00-000d_Y00-000d_S0-001mm.fits`,
- `Usamples_median_withJitter_X00-000d_Y00-000d_S0-001mm.fits`,
- `Gsamples_median_withJitter_X00-000d_Y00-000d_S0-001mm.fits`.

These files are sampled at 10x the telescope's resolution (i.e., each pixel has a side
length of 0.1/10=0.01 arcsec).

Since the accuracy of our calculations depends on having a finely sampled PSF, we also
prepared 3 files that linearly interpolate the data above to obtain PSFs that are 20x the
telescope's resolution (i.e., each pixel has a side length of 0.005 arcsec). These files
are:

- `uv_psf_20x_sampled.fits`,
- `u_psf_20x_sampled.fits`,
- `g_psf_20x_sampled.fits`.

We use these 20x supersampled PSFs as our default PSFs in the exposure time calculator.
The process to create these 3 files is documented in the
[`make_oversampled_psfs.ipynb`](https://github.com/CASTOR-telescope/ETC_notebooks/blob/master/make_oversampled_psfs.ipynb)
notebook.
