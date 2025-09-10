:py:mod:`castor_etc.photometry`
===============================

.. py:module:: castor_etc.photometry

.. autodoc2-docstring:: castor_etc.photometry
   :parser: myst
   :allowtitles:

Module Contents
---------------

Classes
~~~~~~~

.. list-table::
   :class: autosummary longtable
   :align: left

   * - :py:obj:`Photometry <castor_etc.photometry.Photometry>`
     - .. autodoc2-docstring:: castor_etc.photometry.Photometry
          :parser: myst
          :summary:

Data
~~~~

.. list-table::
   :class: autosummary longtable
   :align: left

   * - :py:obj:`_OPTIMAL_APER_FACTOR <castor_etc.photometry._OPTIMAL_APER_FACTOR>`
     - .. autodoc2-docstring:: castor_etc.photometry._OPTIMAL_APER_FACTOR
          :parser: myst
          :summary:
   * - :py:obj:`_SUPERSAMPLE_FACTOR <castor_etc.photometry._SUPERSAMPLE_FACTOR>`
     - .. autodoc2-docstring:: castor_etc.photometry._SUPERSAMPLE_FACTOR
          :parser: myst
          :summary:

API
~~~

.. py:data:: _OPTIMAL_APER_FACTOR
   :canonical: castor_etc.photometry._OPTIMAL_APER_FACTOR
   :value: 1.4

   .. autodoc2-docstring:: castor_etc.photometry._OPTIMAL_APER_FACTOR
      :parser: myst

.. py:data:: _SUPERSAMPLE_FACTOR
   :canonical: castor_etc.photometry._SUPERSAMPLE_FACTOR
   :value: 20

   .. autodoc2-docstring:: castor_etc.photometry._SUPERSAMPLE_FACTOR
      :parser: myst

.. py:class:: Photometry(TelescopeObj, SourceObj, BackgroundObj)
   :canonical: castor_etc.photometry.Photometry

   .. autodoc2-docstring:: castor_etc.photometry.Photometry
      :parser: myst

   .. rubric:: Initialization

   .. autodoc2-docstring:: castor_etc.photometry.Photometry.__init__
      :parser: myst

   .. py:method:: copy()
      :canonical: castor_etc.photometry.Photometry.copy

      .. autodoc2-docstring:: castor_etc.photometry.Photometry.copy
         :parser: myst

   .. py:method:: _assign_exact_npix()
      :canonical: castor_etc.photometry.Photometry._assign_exact_npix

      .. autodoc2-docstring:: castor_etc.photometry.Photometry._assign_exact_npix
         :parser: myst

   .. py:method:: _create_aper_arrs(half_x, half_y, center, overwrite=False)
      :canonical: castor_etc.photometry.Photometry._create_aper_arrs

      .. autodoc2-docstring:: castor_etc.photometry.Photometry._create_aper_arrs
         :parser: myst

   .. py:method:: _calc_source_weights(center)
      :canonical: castor_etc.photometry.Photometry._calc_source_weights

      .. autodoc2-docstring:: castor_etc.photometry.Photometry._calc_source_weights
         :parser: myst

   .. py:method:: show_source_weights(passband, mark_source=False, source_markersize=4, norm=None, plot=True)
      :canonical: castor_etc.photometry.Photometry.show_source_weights

      .. autodoc2-docstring:: castor_etc.photometry.Photometry.show_source_weights
         :parser: myst

   .. py:method:: show_aper_weights(plot=True)
      :canonical: castor_etc.photometry.Photometry.show_aper_weights

      .. autodoc2-docstring:: castor_etc.photometry.Photometry.show_aper_weights
         :parser: myst

   .. py:method:: set_background_weights(sky_background_weights)
      :canonical: castor_etc.photometry.Photometry.set_background_weights

      .. autodoc2-docstring:: castor_etc.photometry.Photometry.set_background_weights
         :parser: myst

   .. py:method:: set_dark_current_weights(dark_current_weights)
      :canonical: castor_etc.photometry.Photometry.set_dark_current_weights

      .. autodoc2-docstring:: castor_etc.photometry.Photometry.set_dark_current_weights
         :parser: myst

   .. py:method:: _bin_arrs_remove_nans(center)
      :canonical: castor_etc.photometry.Photometry._bin_arrs_remove_nans

      .. autodoc2-docstring:: castor_etc.photometry.Photometry._bin_arrs_remove_nans
         :parser: myst

   .. py:method:: _rotate_ab_to_xy(a, b, rotation, px_scale_arcsec)
      :canonical: castor_etc.photometry.Photometry._rotate_ab_to_xy
      :staticmethod:

      .. autodoc2-docstring:: castor_etc.photometry.Photometry._rotate_ab_to_xy
         :parser: myst

   .. py:method:: use_optimal_aperture(factor=_OPTIMAL_APER_FACTOR, quiet=False, overwrite=False)
      :canonical: castor_etc.photometry.Photometry.use_optimal_aperture

      .. autodoc2-docstring:: castor_etc.photometry.Photometry.use_optimal_aperture
         :parser: myst

   .. py:method:: use_elliptical_aperture(a, b, center=[0, 0] << u.arcsec, rotation=0, quiet=False, overwrite=False)
      :canonical: castor_etc.photometry.Photometry.use_elliptical_aperture

      .. autodoc2-docstring:: castor_etc.photometry.Photometry.use_elliptical_aperture
         :parser: myst

   .. py:method:: use_rectangular_aperture(width, length, center=[0, 0] << u.arcsec, quiet=False, overwrite=False)
      :canonical: castor_etc.photometry.Photometry.use_rectangular_aperture

      .. autodoc2-docstring:: castor_etc.photometry.Photometry.use_rectangular_aperture
         :parser: myst

   .. py:method:: _calc_snr_from_t(t, signal, totskynoise, darkcurrent, readnoise, read_npix, nread=1)
      :canonical: castor_etc.photometry.Photometry._calc_snr_from_t
      :staticmethod:

      .. autodoc2-docstring:: castor_etc.photometry.Photometry._calc_snr_from_t
         :parser: myst

   .. py:method:: _calc_t_from_snr(snr, signal, totskynoise, darkcurrent, readnoise, read_npix, nread=1)
      :canonical: castor_etc.photometry.Photometry._calc_t_from_snr
      :staticmethod:

      .. autodoc2-docstring:: castor_etc.photometry.Photometry._calc_t_from_snr
         :parser: myst

   .. py:method:: calc_snr_or_t(t=None, snr=None, reddening=0, npix=None, nread=1, quiet=False)
      :canonical: castor_etc.photometry.Photometry.calc_snr_or_t

      .. autodoc2-docstring:: castor_etc.photometry.Photometry.calc_snr_or_t
         :parser: myst
