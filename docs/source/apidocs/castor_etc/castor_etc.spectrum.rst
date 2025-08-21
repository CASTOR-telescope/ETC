:py:mod:`castor_etc.spectrum`
=============================

.. py:module:: castor_etc.spectrum

.. autodoc2-docstring:: castor_etc.spectrum
   :allowtitles:

Module Contents
---------------

Classes
~~~~~~~

.. list-table::
   :class: autosummary longtable
   :align: left

   * - :py:obj:`SpectrumMixin <castor_etc.spectrum.SpectrumMixin>`
     - .. autodoc2-docstring:: castor_etc.spectrum.SpectrumMixin
          :summary:
   * - :py:obj:`NormMixin <castor_etc.spectrum.NormMixin>`
     - .. autodoc2-docstring:: castor_etc.spectrum.NormMixin
          :summary:

Functions
~~~~~~~~~

.. list-table::
   :class: autosummary longtable
   :align: left

   * - :py:obj:`getStarData <castor_etc.spectrum.getStarData>`
     - .. autodoc2-docstring:: castor_etc.spectrum.getStarData
          :summary:
   * - :py:obj:`interp_EEM_table <castor_etc.spectrum.interp_EEM_table>`
     - .. autodoc2-docstring:: castor_etc.spectrum.interp_EEM_table
          :summary:
   * - :py:obj:`redshift_wavelengths <castor_etc.spectrum.redshift_wavelengths>`
     - .. autodoc2-docstring:: castor_etc.spectrum.redshift_wavelengths
          :summary:

API
~~~

.. py:function:: getStarData(temperature, metallicity, logg, stellar_model_dir, model_grid='ATLAS9')
   :canonical: castor_etc.spectrum.getStarData

   .. autodoc2-docstring:: castor_etc.spectrum.getStarData

.. py:function:: interp_EEM_table(Teff=[], Gmag=[], Bpmag=[], Rpmag=[])
   :canonical: castor_etc.spectrum.interp_EEM_table

   .. autodoc2-docstring:: castor_etc.spectrum.interp_EEM_table

.. py:function:: redshift_wavelengths(wavelengths, redshift)
   :canonical: castor_etc.spectrum.redshift_wavelengths

   .. autodoc2-docstring:: castor_etc.spectrum.redshift_wavelengths

.. py:class:: SpectrumMixin
   :canonical: castor_etc.spectrum.SpectrumMixin

   .. autodoc2-docstring:: castor_etc.spectrum.SpectrumMixin

   .. py:method:: _check_existing_spectrum(overwrite, quiet=False)
      :canonical: castor_etc.spectrum.SpectrumMixin._check_existing_spectrum

      .. autodoc2-docstring:: castor_etc.spectrum.SpectrumMixin._check_existing_spectrum

   .. py:method:: spectrum_erg_to_photon()
      :canonical: castor_etc.spectrum.SpectrumMixin.spectrum_erg_to_photon

      .. autodoc2-docstring:: castor_etc.spectrum.SpectrumMixin.spectrum_erg_to_photon

   .. py:method:: redshift_wavelengths(redshift)
      :canonical: castor_etc.spectrum.SpectrumMixin.redshift_wavelengths

      .. autodoc2-docstring:: castor_etc.spectrum.SpectrumMixin.redshift_wavelengths

   .. py:method:: generate_uniform(wavelengths, value, unit='ABmag', overwrite=False, quiet=False)
      :canonical: castor_etc.spectrum.SpectrumMixin.generate_uniform

      .. autodoc2-docstring:: castor_etc.spectrum.SpectrumMixin.generate_uniform

   .. py:method:: generate_bb(T, redshift=0.0, emissivity=1.0, wavelengths=None, limits=[0.09, 1.2] << u.um, resolution=1 << u.nm, radius=1, dist=1 << u.kpc, overwrite=False, quiet=False)
      :canonical: castor_etc.spectrum.SpectrumMixin.generate_bb

      .. autodoc2-docstring:: castor_etc.spectrum.SpectrumMixin.generate_bb

   .. py:method:: generate_power_law(ref_wavelength, wavelengths, exponent, overwrite=False, quiet=False)
      :canonical: castor_etc.spectrum.SpectrumMixin.generate_power_law

      .. autodoc2-docstring:: castor_etc.spectrum.SpectrumMixin.generate_power_law

   .. py:method:: _generate_gaussian(wavelengths, spectrum, center, fwhm, peak=None, tot_flux=None, add=True, abs_peak=True)
      :canonical: castor_etc.spectrum.SpectrumMixin._generate_gaussian
      :staticmethod:

      .. autodoc2-docstring:: castor_etc.spectrum.SpectrumMixin._generate_gaussian

   .. py:method:: _generate_lorentzian(wavelengths, spectrum, center, fwhm, peak=None, tot_flux=None, add=True, abs_peak=True)
      :canonical: castor_etc.spectrum.SpectrumMixin._generate_lorentzian
      :staticmethod:

      .. autodoc2-docstring:: castor_etc.spectrum.SpectrumMixin._generate_lorentzian

   .. py:method:: add_emission_line(center, fwhm, peak=None, tot_flux=None, shape='gaussian', abs_peak=False)
      :canonical: castor_etc.spectrum.SpectrumMixin.add_emission_line

      .. autodoc2-docstring:: castor_etc.spectrum.SpectrumMixin.add_emission_line

   .. py:method:: add_absorption_line(center, fwhm, dip=None, tot_flux=None, shape='gaussian', abs_dip=False)
      :canonical: castor_etc.spectrum.SpectrumMixin.add_absorption_line

      .. autodoc2-docstring:: castor_etc.spectrum.SpectrumMixin.add_absorption_line

   .. py:method:: generate_emission_line(center, fwhm, peak=None, tot_flux=None, shape='gaussian', limits=[100, 1200] << u.nm, overwrite=False, quiet=False)
      :canonical: castor_etc.spectrum.SpectrumMixin.generate_emission_line

      .. autodoc2-docstring:: castor_etc.spectrum.SpectrumMixin.generate_emission_line

   .. py:method:: set_spectrum(wavelengths, spectrum, unit, overwrite=False, quiet=False)
      :canonical: castor_etc.spectrum.SpectrumMixin.set_spectrum

      .. autodoc2-docstring:: castor_etc.spectrum.SpectrumMixin.set_spectrum

   .. py:method:: use_custom_spectrum(filepath, wavelength_unit=u.AA, overwrite=False, quiet=False)
      :canonical: castor_etc.spectrum.SpectrumMixin.use_custom_spectrum

      .. autodoc2-docstring:: castor_etc.spectrum.SpectrumMixin.use_custom_spectrum

   .. py:method:: use_galaxy_spectrum(gal_type, overwrite=False, quiet=False)
      :canonical: castor_etc.spectrum.SpectrumMixin.use_galaxy_spectrum

      .. autodoc2-docstring:: castor_etc.spectrum.SpectrumMixin.use_galaxy_spectrum

   .. py:method:: _calc_xy()
      :canonical: castor_etc.spectrum.SpectrumMixin._calc_xy

      .. autodoc2-docstring:: castor_etc.spectrum.SpectrumMixin._calc_xy

   .. py:method:: _search_gaia()
      :canonical: castor_etc.spectrum.SpectrumMixin._search_gaia

      .. autodoc2-docstring:: castor_etc.spectrum.SpectrumMixin._search_gaia

   .. py:method:: _specify_target_parameters(run_gaia_search=True)
      :canonical: castor_etc.spectrum.SpectrumMixin._specify_target_parameters

      .. autodoc2-docstring:: castor_etc.spectrum.SpectrumMixin._specify_target_parameters

   .. py:method:: use_gaia_spectrum(TelescopeObj, ra=None, dec=None, srch_Gmax=21.0, srch_nmax=100, srch_rad=None, Teff=None, Gmag=None, logg=None, radius=None, metallicity=None, Bpmag=None, Rpmag=None, stellar_model_grid='ATLAS9', stellar_model_dir=None, bkg_sources=True, fov=None, fov_pa=0 * u.deg, overwrite=False, quiet=False)
      :canonical: castor_etc.spectrum.SpectrumMixin.use_gaia_spectrum

      .. autodoc2-docstring:: castor_etc.spectrum.SpectrumMixin.use_gaia_spectrum

   .. py:method:: use_pickles_spectrum(spectral_class, overwrite=False, quiet=False)
      :canonical: castor_etc.spectrum.SpectrumMixin.use_pickles_spectrum

      .. autodoc2-docstring:: castor_etc.spectrum.SpectrumMixin.use_pickles_spectrum

   .. py:method:: show_spectrum(plot=True)
      :canonical: castor_etc.spectrum.SpectrumMixin.show_spectrum

      .. autodoc2-docstring:: castor_etc.spectrum.SpectrumMixin.show_spectrum

   .. py:method:: calc_redleak_frac(TelescopeObj, quiet=False)
      :canonical: castor_etc.spectrum.SpectrumMixin.calc_redleak_frac

      .. autodoc2-docstring:: castor_etc.spectrum.SpectrumMixin.calc_redleak_frac

.. py:class:: NormMixin
   :canonical: castor_etc.spectrum.NormMixin

   .. autodoc2-docstring:: castor_etc.spectrum.NormMixin

   .. py:method:: norm_to_star(spectrum, radius=1, dist=1 << u.kpc)
      :canonical: castor_etc.spectrum.NormMixin.norm_to_star
      :staticmethod:

      .. autodoc2-docstring:: castor_etc.spectrum.NormMixin.norm_to_star

   .. py:method:: norm_to_AB_mag(ab_mag, passband=None, TelescopeObj=None)
      :canonical: castor_etc.spectrum.NormMixin.norm_to_AB_mag

      .. autodoc2-docstring:: castor_etc.spectrum.NormMixin.norm_to_AB_mag

   .. py:method:: norm_luminosity_dist(luminosity, dist)
      :canonical: castor_etc.spectrum.NormMixin.norm_luminosity_dist

      .. autodoc2-docstring:: castor_etc.spectrum.NormMixin.norm_luminosity_dist

   .. py:method:: get_AB_mag(TelescopeObj=None)
      :canonical: castor_etc.spectrum.NormMixin.get_AB_mag

      .. autodoc2-docstring:: castor_etc.spectrum.NormMixin.get_AB_mag
