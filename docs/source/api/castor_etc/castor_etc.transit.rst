:py:mod:`castor_etc.transit`
============================

.. py:module:: castor_etc.transit

.. autodoc2-docstring:: castor_etc.transit
   :allowtitles:

Module Contents
---------------

Classes
~~~~~~~

.. list-table::
   :class: autosummary longtable
   :align: left

   * - :py:obj:`Observation <castor_etc.transit.Observation>`
     - .. autodoc2-docstring:: castor_etc.transit.Observation
          :summary:

Functions
~~~~~~~~~

.. list-table::
   :class: autosummary longtable
   :align: left

   * - :py:obj:`addflux2pix <castor_etc.transit.addflux2pix>`
     - .. autodoc2-docstring:: castor_etc.transit.addflux2pix
          :summary:
   * - :py:obj:`gen_unconv_image <castor_etc.transit.gen_unconv_image>`
     - .. autodoc2-docstring:: castor_etc.transit.gen_unconv_image
          :summary:

Data
~~~~

.. list-table::
   :class: autosummary longtable
   :align: left

   * - :py:obj:`_OPTIMAL_APER_FACTOR <castor_etc.transit._OPTIMAL_APER_FACTOR>`
     - .. autodoc2-docstring:: castor_etc.transit._OPTIMAL_APER_FACTOR
          :summary:

API
~~~

.. py:data:: _OPTIMAL_APER_FACTOR
   :canonical: castor_etc.transit._OPTIMAL_APER_FACTOR
   :value: 1.4

   .. autodoc2-docstring:: castor_etc.transit._OPTIMAL_APER_FACTOR

.. py:function:: addflux2pix(px, py, pixels, fmod)
   :canonical: castor_etc.transit.addflux2pix

   .. autodoc2-docstring:: castor_etc.transit.addflux2pix

.. py:function:: gen_unconv_image(pars, starmodel_flux, xcoo, ycoo)
   :canonical: castor_etc.transit.gen_unconv_image

   .. autodoc2-docstring:: castor_etc.transit.gen_unconv_image

.. py:class:: Observation(TelescopeObj, SourceObj, BackgroundObj, stellar_model_dir)
   :canonical: castor_etc.transit.Observation

   .. autodoc2-docstring:: castor_etc.transit.Observation

   .. rubric:: Initialization

   .. autodoc2-docstring:: castor_etc.transit.Observation.__init__

   .. py:method:: calc_sky_background_erate()
      :canonical: castor_etc.transit.Observation.calc_sky_background_erate

      .. autodoc2-docstring:: castor_etc.transit.Observation.calc_sky_background_erate

   .. py:method:: _photon_count(temp=5780.0, metallicity=0.0, logg=4.44, Gmag=7.0, Gmag_abs=4.635, radius=1.0)
      :canonical: castor_etc.transit.Observation._photon_count

      .. autodoc2-docstring:: castor_etc.transit.Observation._photon_count

   .. py:method:: specify_bandpass(passband_name=None)
      :canonical: castor_etc.transit.Observation.specify_bandpass

      .. autodoc2-docstring:: castor_etc.transit.Observation.specify_bandpass

   .. py:method:: id_guide_stars(gs_criteria=None, plot_SN=False)
      :canonical: castor_etc.transit.Observation.id_guide_stars

      .. autodoc2-docstring:: castor_etc.transit.Observation.id_guide_stars

   .. py:method:: _remove_aper_mask_nan_row_col(center)
      :canonical: castor_etc.transit.Observation._remove_aper_mask_nan_row_col

      .. autodoc2-docstring:: castor_etc.transit.Observation._remove_aper_mask_nan_row_col

   .. py:method:: _calc_source_weights(center)
      :canonical: castor_etc.transit.Observation._calc_source_weights

      .. autodoc2-docstring:: castor_etc.transit.Observation._calc_source_weights

   .. py:method:: _create_aper_arrs(half_x, half_y, center, overwrite=False)
      :canonical: castor_etc.transit.Observation._create_aper_arrs

      .. autodoc2-docstring:: castor_etc.transit.Observation._create_aper_arrs

   .. py:method:: _use_optimal_aperture(factor=_OPTIMAL_APER_FACTOR, overwrite=False)
      :canonical: castor_etc.transit.Observation._use_optimal_aperture

      .. autodoc2-docstring:: castor_etc.transit.Observation._use_optimal_aperture

   .. py:method:: _point_source_sim(target_flux_fraction, scene_phot_count)
      :canonical: castor_etc.transit.Observation._point_source_sim

      .. autodoc2-docstring:: castor_etc.transit.Observation._point_source_sim

   .. py:method:: scene_sim(all_sources=True, return_scene=False, update_gaia=True, quiet=None, return_SN_only=False)
      :canonical: castor_etc.transit.Observation.scene_sim

      .. autodoc2-docstring:: castor_etc.transit.Observation.scene_sim

   .. py:method:: plot_fov(plot_guide_stars=True, vmin=None, vmax=None, add_scene_sim=True)
      :canonical: castor_etc.transit.Observation.plot_fov

      .. autodoc2-docstring:: castor_etc.transit.Observation.plot_fov

   .. py:method:: specify_pl_model(RpRs, P, t0, b, aRs)
      :canonical: castor_etc.transit.Observation.specify_pl_model

      .. autodoc2-docstring:: castor_etc.transit.Observation.specify_pl_model

   .. py:method:: specify_exposure_parameters(exptime=60 * u.second, nstack=10, tstart=0.0 * u.d, tend=6.0 / 24.0 * u.d)
      :canonical: castor_etc.transit.Observation.specify_exposure_parameters

      .. autodoc2-docstring:: castor_etc.transit.Observation.specify_exposure_parameters

   .. py:method:: calc_pl_model(model='pytransit_QuadraticModel', t_grid=[], exp_time=-1)
      :canonical: castor_etc.transit.Observation.calc_pl_model

      .. autodoc2-docstring:: castor_etc.transit.Observation.calc_pl_model

   .. py:method:: lc_sim(quiet=False, return_lc=False, cadence=-1, exp_time=-1)
      :canonical: castor_etc.transit.Observation.lc_sim

      .. autodoc2-docstring:: castor_etc.transit.Observation.lc_sim

   .. py:method:: plot_lc(plot_model=True, exp_time=-1, t_from_mid=True, t_unit='d')
      :canonical: castor_etc.transit.Observation.plot_lc

      .. autodoc2-docstring:: castor_etc.transit.Observation.plot_lc
