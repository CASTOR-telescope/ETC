:py:mod:`castor_etc.background`
===============================

.. py:module:: castor_etc.background

.. autodoc2-docstring:: castor_etc.background
   :parser: myst
   :allowtitles:

Module Contents
---------------

Classes
~~~~~~~

.. list-table::
   :class: autosummary longtable
   :align: left

   * - :py:obj:`Background <castor_etc.background.Background>`
     - .. autodoc2-docstring:: castor_etc.background.Background
          :parser: myst
          :summary:

API
~~~

.. py:class:: Background(earthshine_file: str | None = join(DATAPATH, 'sky_background', 'earthshine.fits'), zodi_file: str | None = join(DATAPATH, 'sky_background', 'zodi.fits'), mags_per_sq_arcsec: dict[float] | None = None)
   :canonical: castor_etc.background.Background

   .. autodoc2-docstring:: castor_etc.background.Background
      :parser: myst

   .. rubric:: Initialization

   .. autodoc2-docstring:: castor_etc.background.Background.__init__
      :parser: myst

   .. py:method:: copy()
      :canonical: castor_etc.background.Background.copy

      .. autodoc2-docstring:: castor_etc.background.Background.copy
         :parser: myst

   .. py:method:: add_geocoronal_emission(flux='avg', wavelength=GEOCORONAL_WAVELENGTH, linewidth=GEOCORONAL_LINEWIDTH)
      :canonical: castor_etc.background.Background.add_geocoronal_emission

      .. autodoc2-docstring:: castor_etc.background.Background.add_geocoronal_emission
         :parser: myst

   .. py:method:: _get_mags_per_sq_arcsec(TelescopeObj)
      :canonical: castor_etc.background.Background._get_mags_per_sq_arcsec

      .. autodoc2-docstring:: castor_etc.background.Background._get_mags_per_sq_arcsec
         :parser: myst

   .. py:method:: calc_mags_per_sq_arcsec(TelescopeObj, overwrite=False)
      :canonical: castor_etc.background.Background.calc_mags_per_sq_arcsec

      .. autodoc2-docstring:: castor_etc.background.Background.calc_mags_per_sq_arcsec
         :parser: myst
