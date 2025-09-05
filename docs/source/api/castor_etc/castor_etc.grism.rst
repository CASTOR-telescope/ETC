:py:mod:`castor_etc.grism`
==========================

.. py:module:: castor_etc.grism

.. autodoc2-docstring:: castor_etc.grism
   :allowtitles:

Module Contents
---------------

Classes
~~~~~~~

.. list-table::
   :class: autosummary longtable
   :align: left

   * - :py:obj:`Grism <castor_etc.grism.Grism>`
     - .. autodoc2-docstring:: castor_etc.grism.Grism
          :summary:

API
~~~

.. py:class:: Grism(TelescopeObj, SourceObj, BackgroundObj)
   :canonical: castor_etc.grism.Grism

   Bases: :py:obj:`object`

   .. autodoc2-docstring:: castor_etc.grism.Grism

   .. rubric:: Initialization

   .. autodoc2-docstring:: castor_etc.grism.Grism.__init__

   .. py:method:: _create_segmentation_map()
      :canonical: castor_etc.grism.Grism._create_segmentation_map

      .. autodoc2-docstring:: castor_etc.grism.Grism._create_segmentation_map

   .. py:method:: disperse(grism_channel='u', check=True)
      :canonical: castor_etc.grism.Grism.disperse

      .. autodoc2-docstring:: castor_etc.grism.Grism.disperse

   .. py:method:: expose(exposure_time=1000)
      :canonical: castor_etc.grism.Grism.expose

      .. autodoc2-docstring:: castor_etc.grism.Grism.expose

   .. py:method:: _calc_sky_background_erate()
      :canonical: castor_etc.grism.Grism._calc_sky_background_erate

      .. autodoc2-docstring:: castor_etc.grism.Grism._calc_sky_background_erate

   .. py:method:: _calculate_tot_unif_noise(Nreads=1, Nbin=1)
      :canonical: castor_etc.grism.Grism._calculate_tot_unif_noise

      .. autodoc2-docstring:: castor_etc.grism.Grism._calculate_tot_unif_noise

   .. py:method:: total_noise(Nreads=1, Nbin=1)
      :canonical: castor_etc.grism.Grism.total_noise

      .. autodoc2-docstring:: castor_etc.grism.Grism.total_noise

   .. py:method:: show_2d_snr_per_resolution()
      :canonical: castor_etc.grism.Grism.show_2d_snr_per_resolution

      .. autodoc2-docstring:: castor_etc.grism.Grism.show_2d_snr_per_resolution

   .. py:method:: show_1d_snr_per_resolution()
      :canonical: castor_etc.grism.Grism.show_1d_snr_per_resolution

      .. autodoc2-docstring:: castor_etc.grism.Grism.show_1d_snr_per_resolution
