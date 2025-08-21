:py:mod:`castor_etc.telescope`
==============================

.. py:module:: castor_etc.telescope

.. autodoc2-docstring:: castor_etc.telescope
   :allowtitles:

Module Contents
---------------

Classes
~~~~~~~

.. list-table::
   :class: autosummary longtable
   :align: left

   * - :py:obj:`Telescope <castor_etc.telescope.Telescope>`
     - .. autodoc2-docstring:: castor_etc.telescope.Telescope
          :summary:

Functions
~~~~~~~~~

.. list-table::
   :class: autosummary longtable
   :align: left

   * - :py:obj:`secant_method <castor_etc.telescope.secant_method>`
     - .. autodoc2-docstring:: castor_etc.telescope.secant_method
          :summary:
   * - :py:obj:`bisection_method <castor_etc.telescope.bisection_method>`
     - .. autodoc2-docstring:: castor_etc.telescope.bisection_method
          :summary:

Data
~~~~

.. list-table::
   :class: autosummary longtable
   :align: left

   * - :py:obj:`_ds9heat <castor_etc.telescope._ds9heat>`
     - .. autodoc2-docstring:: castor_etc.telescope._ds9heat
          :summary:
   * - :py:obj:`_ds9heat_cmap <castor_etc.telescope._ds9heat_cmap>`
     - .. autodoc2-docstring:: castor_etc.telescope._ds9heat_cmap
          :summary:

API
~~~

.. py:data:: _ds9heat
   :canonical: castor_etc.telescope._ds9heat
   :value: None

   .. autodoc2-docstring:: castor_etc.telescope._ds9heat

.. py:data:: _ds9heat_cmap
   :canonical: castor_etc.telescope._ds9heat_cmap
   :value: 'LinearSegmentedColormap(...)'

   .. autodoc2-docstring:: castor_etc.telescope._ds9heat_cmap

.. py:function:: secant_method(f, x0, x1, tol=1e-06, max_iter=100)
   :canonical: castor_etc.telescope.secant_method

   .. autodoc2-docstring:: castor_etc.telescope.secant_method

.. py:function:: bisection_method(f, x0, x1, tol=1e-06, max_iter=100)
   :canonical: castor_etc.telescope.bisection_method

   .. autodoc2-docstring:: castor_etc.telescope.bisection_method

.. py:class:: Telescope(passbands=params.PASSBANDS, passband_limits=params.PASSBAND_LIMITS, passband_response_filepaths=params.PASSBAND_FILEPATHS, passband_response_fileunits=params.PASSBAND_FILEUNITS, passband_resolution=params.PASSBAND_RESOLUTION, passband_pivots=None, phot_zpts=None, phot_zpts_kwargs={'ab_mags': {'uv': [25.5, 23.5], 'u': [25.5, 23.5], 'g': [25.5, 23.5]}, 'method': 'secant', 'tol': 0.0002, 'max_iter': 100}, psf_filepaths=params.PSF_FILEPATHS, psf_supersample_factor=params.PSF_SUPERSAMPLE_FACTOR, fwhm=params.FWHM, px_scale=params.PX_SCALE, transit_fov=params.TRANSIT_FOV, ifov_dimen=params.IFOV_DIMEN, transit_ccd_dim=params.TRANSIT_CCD_DIMENSIONS, mp=params.MP, mirror_diameter=params.MIRROR_DIAMETER, dark_current=params.DARK_CURRENT, bias=params.BIAS, read_noise=params.READ_NOISE, gain=params.GAIN, redleak_thresholds=params.REDLEAK_THRESHOLDS, extinction_coeffs=params.EXTINCTION_COEFFS, show_warnings=True)
   :canonical: castor_etc.telescope.Telescope

   .. autodoc2-docstring:: castor_etc.telescope.Telescope

   .. rubric:: Initialization

   .. autodoc2-docstring:: castor_etc.telescope.Telescope.__init__

   .. py:method:: copy()
      :canonical: castor_etc.telescope.Telescope.copy

      .. autodoc2-docstring:: castor_etc.telescope.Telescope.copy

   .. py:method:: show_psf(passband, norm=None, plot=True)
      :canonical: castor_etc.telescope.Telescope.show_psf

      .. autodoc2-docstring:: castor_etc.telescope.Telescope.show_psf

   .. py:method:: calc_pivot_wavelength(wavelengths, response, response_func='EE')
      :canonical: castor_etc.telescope.Telescope.calc_pivot_wavelength
      :staticmethod:

      .. autodoc2-docstring:: castor_etc.telescope.Telescope.calc_pivot_wavelength

   .. py:method:: load_passbands(filepaths, limits, file_units, resolution=1 << u.nm, interp_kind='linear')
      :canonical: castor_etc.telescope.Telescope.load_passbands
      :staticmethod:

      .. autodoc2-docstring:: castor_etc.telescope.Telescope.load_passbands

   .. py:method:: calc_phot_zpts(passband_curves, ab_mags, mirror_area, method='secant', tol=0.001, max_iter=100)
      :canonical: castor_etc.telescope.Telescope.calc_phot_zpts
      :staticmethod:

      .. autodoc2-docstring:: castor_etc.telescope.Telescope.calc_phot_zpts
