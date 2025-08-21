:py:mod:`castor_etc.uvmos_spectroscopy`
=======================================

.. py:module:: castor_etc.uvmos_spectroscopy

.. autodoc2-docstring:: castor_etc.uvmos_spectroscopy
   :allowtitles:

Module Contents
---------------

Classes
~~~~~~~

.. list-table::
   :class: autosummary longtable
   :align: left

   * - :py:obj:`UVMOS_Spectroscopy <castor_etc.uvmos_spectroscopy.UVMOS_Spectroscopy>`
     - .. autodoc2-docstring:: castor_etc.uvmos_spectroscopy.UVMOS_Spectroscopy
          :summary:

Functions
~~~~~~~~~

.. list-table::
   :class: autosummary longtable
   :align: left

   * - :py:obj:`Gaussian <castor_etc.uvmos_spectroscopy.Gaussian>`
     - .. autodoc2-docstring:: castor_etc.uvmos_spectroscopy.Gaussian
          :summary:
   * - :py:obj:`Gaussian2D <castor_etc.uvmos_spectroscopy.Gaussian2D>`
     - .. autodoc2-docstring:: castor_etc.uvmos_spectroscopy.Gaussian2D
          :summary:

API
~~~

.. py:function:: Gaussian(x, x0, sigma, a)
   :canonical: castor_etc.uvmos_spectroscopy.Gaussian

   .. autodoc2-docstring:: castor_etc.uvmos_spectroscopy.Gaussian

.. py:function:: Gaussian2D(x, y, sigma, a=1, x0=0, y0=0)
   :canonical: castor_etc.uvmos_spectroscopy.Gaussian2D

   .. autodoc2-docstring:: castor_etc.uvmos_spectroscopy.Gaussian2D

.. py:class:: UVMOS_Spectroscopy(TelescopeObj, SourceObj, BackgroundObj)
   :canonical: castor_etc.uvmos_spectroscopy.UVMOS_Spectroscopy

   .. autodoc2-docstring:: castor_etc.uvmos_spectroscopy.UVMOS_Spectroscopy

   .. rubric:: Initialization

   .. autodoc2-docstring:: castor_etc.uvmos_spectroscopy.UVMOS_Spectroscopy.__init__

   .. py:method:: specify_slit(slit_width=0.214 * u.arcsec, slit_height=1 * u.arcsec)
      :canonical: castor_etc.uvmos_spectroscopy.UVMOS_Spectroscopy.specify_slit

      .. autodoc2-docstring:: castor_etc.uvmos_spectroscopy.UVMOS_Spectroscopy.specify_slit

   .. py:method:: show_slit()
      :canonical: castor_etc.uvmos_spectroscopy.UVMOS_Spectroscopy.show_slit

      .. autodoc2-docstring:: castor_etc.uvmos_spectroscopy.UVMOS_Spectroscopy.show_slit

   .. py:method:: _calc_slit_transmission(print_transmission_fact=False)
      :canonical: castor_etc.uvmos_spectroscopy.UVMOS_Spectroscopy._calc_slit_transmission

      .. autodoc2-docstring:: castor_etc.uvmos_spectroscopy.UVMOS_Spectroscopy._calc_slit_transmission

   .. py:method:: calc_source_pix_weights()
      :canonical: castor_etc.uvmos_spectroscopy.UVMOS_Spectroscopy.calc_source_pix_weights

      .. autodoc2-docstring:: castor_etc.uvmos_spectroscopy.UVMOS_Spectroscopy.calc_source_pix_weights

   .. py:method:: show_source_pix_weights()
      :canonical: castor_etc.uvmos_spectroscopy.UVMOS_Spectroscopy.show_source_pix_weights

      .. autodoc2-docstring:: castor_etc.uvmos_spectroscopy.UVMOS_Spectroscopy.show_source_pix_weights

   .. py:method:: show_slit_image(wave)
      :canonical: castor_etc.uvmos_spectroscopy.UVMOS_Spectroscopy.show_slit_image

      .. autodoc2-docstring:: castor_etc.uvmos_spectroscopy.UVMOS_Spectroscopy.show_slit_image

   .. py:method:: _extraction(detector, pix_waves, extraction_width, extraction_lowerlim, extraction_upperlim)
      :canonical: castor_etc.uvmos_spectroscopy.UVMOS_Spectroscopy._extraction

      .. autodoc2-docstring:: castor_etc.uvmos_spectroscopy.UVMOS_Spectroscopy._extraction

   .. py:method:: _getTransmission(x)
      :canonical: castor_etc.uvmos_spectroscopy.UVMOS_Spectroscopy._getTransmission

      .. autodoc2-docstring:: castor_etc.uvmos_spectroscopy.UVMOS_Spectroscopy._getTransmission

   .. py:method:: _getDispersion(x)
      :canonical: castor_etc.uvmos_spectroscopy.UVMOS_Spectroscopy._getDispersion

      .. autodoc2-docstring:: castor_etc.uvmos_spectroscopy.UVMOS_Spectroscopy._getDispersion

   .. py:method:: showTransmission()
      :canonical: castor_etc.uvmos_spectroscopy.UVMOS_Spectroscopy.showTransmission

      .. autodoc2-docstring:: castor_etc.uvmos_spectroscopy.UVMOS_Spectroscopy.showTransmission

   .. py:method:: _calc_sigmaPix(dispersion)
      :canonical: castor_etc.uvmos_spectroscopy.UVMOS_Spectroscopy._calc_sigmaPix

      .. autodoc2-docstring:: castor_etc.uvmos_spectroscopy.UVMOS_Spectroscopy._calc_sigmaPix

   .. py:method:: _calc_sigmaPSF(dispersion)
      :canonical: castor_etc.uvmos_spectroscopy.UVMOS_Spectroscopy._calc_sigmaPSF

      .. autodoc2-docstring:: castor_etc.uvmos_spectroscopy.UVMOS_Spectroscopy._calc_sigmaPSF

   .. py:method:: _calcR(x, disp)
      :canonical: castor_etc.uvmos_spectroscopy.UVMOS_Spectroscopy._calcR

      .. autodoc2-docstring:: castor_etc.uvmos_spectroscopy.UVMOS_Spectroscopy._calcR

   .. py:method:: showResolvingPower()
      :canonical: castor_etc.uvmos_spectroscopy.UVMOS_Spectroscopy.showResolvingPower

      .. autodoc2-docstring:: castor_etc.uvmos_spectroscopy.UVMOS_Spectroscopy.showResolvingPower

   .. py:method:: calc_source_CASTORSpectrum(extraction_width=int(1), extraction_lowerlim=0, extraction_upperlim='max')
      :canonical: castor_etc.uvmos_spectroscopy.UVMOS_Spectroscopy.calc_source_CASTORSpectrum

      .. autodoc2-docstring:: castor_etc.uvmos_spectroscopy.UVMOS_Spectroscopy.calc_source_CASTORSpectrum

   .. py:method:: calc_background_CASTORSpectrum(extraction_width=int(1), extraction_lowerlim=0, extraction_upperlim='max')
      :canonical: castor_etc.uvmos_spectroscopy.UVMOS_Spectroscopy.calc_background_CASTORSpectrum

      .. autodoc2-docstring:: castor_etc.uvmos_spectroscopy.UVMOS_Spectroscopy.calc_background_CASTORSpectrum

   .. py:method:: calc_snr_from_t(t, wave, nread=1)
      :canonical: castor_etc.uvmos_spectroscopy.UVMOS_Spectroscopy.calc_snr_from_t

      .. autodoc2-docstring:: castor_etc.uvmos_spectroscopy.UVMOS_Spectroscopy.calc_snr_from_t

   .. py:method:: calc_t_from_snr(snr, wave, nread=1)
      :canonical: castor_etc.uvmos_spectroscopy.UVMOS_Spectroscopy.calc_t_from_snr

      .. autodoc2-docstring:: castor_etc.uvmos_spectroscopy.UVMOS_Spectroscopy.calc_t_from_snr
