:py:mod:`castor_etc.conversions`
================================

.. py:module:: castor_etc.conversions

.. autodoc2-docstring:: castor_etc.conversions
   :parser: myst
   :allowtitles:

Module Contents
---------------

Functions
~~~~~~~~~

.. list-table::
   :class: autosummary longtable
   :align: left

   * - :py:obj:`calc_photon_energy <castor_etc.conversions.calc_photon_energy>`
     - .. autodoc2-docstring:: castor_etc.conversions.calc_photon_energy
          :parser: myst
          :summary:
   * - :py:obj:`convert_freq_wavelength <castor_etc.conversions.convert_freq_wavelength>`
     - .. autodoc2-docstring:: castor_etc.conversions.convert_freq_wavelength
          :parser: myst
          :summary:
   * - :py:obj:`flam_to_photlam <castor_etc.conversions.flam_to_photlam>`
     - .. autodoc2-docstring:: castor_etc.conversions.flam_to_photlam
          :parser: myst
          :summary:
   * - :py:obj:`fnu_to_photlam <castor_etc.conversions.fnu_to_photlam>`
     - .. autodoc2-docstring:: castor_etc.conversions.fnu_to_photlam
          :parser: myst
          :summary:
   * - :py:obj:`fnu_to_flam <castor_etc.conversions.fnu_to_flam>`
     - .. autodoc2-docstring:: castor_etc.conversions.fnu_to_flam
          :parser: myst
          :summary:
   * - :py:obj:`flam_to_fnu <castor_etc.conversions.flam_to_fnu>`
     - .. autodoc2-docstring:: castor_etc.conversions.flam_to_fnu
          :parser: myst
          :summary:
   * - :py:obj:`convert_rel_abs_mag <castor_etc.conversions.convert_rel_abs_mag>`
     - .. autodoc2-docstring:: castor_etc.conversions.convert_rel_abs_mag
          :parser: myst
          :summary:
   * - :py:obj:`flux_to_mag <castor_etc.conversions.flux_to_mag>`
     - .. autodoc2-docstring:: castor_etc.conversions.flux_to_mag
          :parser: myst
          :summary:
   * - :py:obj:`mag_to_flux <castor_etc.conversions.mag_to_flux>`
     - .. autodoc2-docstring:: castor_etc.conversions.mag_to_flux
          :parser: myst
          :summary:
   * - :py:obj:`convert_AB_ST_mag <castor_etc.conversions.convert_AB_ST_mag>`
     - .. autodoc2-docstring:: castor_etc.conversions.convert_AB_ST_mag
          :parser: myst
          :summary:
   * - :py:obj:`flam_to_AB_mag <castor_etc.conversions.flam_to_AB_mag>`
     - .. autodoc2-docstring:: castor_etc.conversions.flam_to_AB_mag
          :parser: myst
          :summary:
   * - :py:obj:`convert_electron_flux_mag <castor_etc.conversions.convert_electron_flux_mag>`
     - .. autodoc2-docstring:: castor_etc.conversions.convert_electron_flux_mag
          :parser: myst
          :summary:

API
~~~

.. py:function:: calc_photon_energy(wavelength=None, frequency=None, wavelength_err=0.0, frequency_err=0.0)
   :canonical: castor_etc.conversions.calc_photon_energy

   .. autodoc2-docstring:: castor_etc.conversions.calc_photon_energy
      :parser: myst

.. py:function:: convert_freq_wavelength(data, to='wavelength', output_unit=u.AA)
   :canonical: castor_etc.conversions.convert_freq_wavelength

   .. autodoc2-docstring:: castor_etc.conversions.convert_freq_wavelength
      :parser: myst

.. py:function:: flam_to_photlam(flam, wavelength)
   :canonical: castor_etc.conversions.flam_to_photlam

   .. autodoc2-docstring:: castor_etc.conversions.flam_to_photlam
      :parser: myst

.. py:function:: fnu_to_photlam(fnu, wavelength)
   :canonical: castor_etc.conversions.fnu_to_photlam

   .. autodoc2-docstring:: castor_etc.conversions.fnu_to_photlam
      :parser: myst

.. py:function:: fnu_to_flam(fnu, wavelength, fnu_err=0.0, wavelength_err=0.0)
   :canonical: castor_etc.conversions.fnu_to_flam

   .. autodoc2-docstring:: castor_etc.conversions.fnu_to_flam
      :parser: myst

.. py:function:: flam_to_fnu(flam, wavelength, flam_err=0.0, wavelength_err=0.0)
   :canonical: castor_etc.conversions.flam_to_fnu

   .. autodoc2-docstring:: castor_etc.conversions.flam_to_fnu
      :parser: myst

.. py:function:: convert_rel_abs_mag(mag, dist, mag_err=0.0, dist_err=0.0, to='abs')
   :canonical: castor_etc.conversions.convert_rel_abs_mag

   .. autodoc2-docstring:: castor_etc.conversions.convert_rel_abs_mag
      :parser: myst

.. py:function:: flux_to_mag(flux, flux_err=0.0, zpt=-48.6, calc_abs=False, dist=None, dist_err=0.0)
   :canonical: castor_etc.conversions.flux_to_mag

   .. autodoc2-docstring:: castor_etc.conversions.flux_to_mag
      :parser: myst

.. py:function:: mag_to_flux(mag, mag_err=0.0, zpt=-48.6)
   :canonical: castor_etc.conversions.mag_to_flux

   .. autodoc2-docstring:: castor_etc.conversions.mag_to_flux
      :parser: myst

.. py:function:: convert_AB_ST_mag(mag, wavelength, to='ABmag')
   :canonical: castor_etc.conversions.convert_AB_ST_mag

   .. autodoc2-docstring:: castor_etc.conversions.convert_AB_ST_mag
      :parser: myst

.. py:function:: flam_to_AB_mag(wavelengths, flam, response)
   :canonical: castor_etc.conversions.flam_to_AB_mag

   .. autodoc2-docstring:: castor_etc.conversions.flam_to_AB_mag
      :parser: myst

.. py:function:: convert_electron_flux_mag(var1, var1_type, var2_type, var1_err=0.0, phot_zpt=None, wavelengths=None, wavelengths_err=0.0)
   :canonical: castor_etc.conversions.convert_electron_flux_mag

   .. autodoc2-docstring:: castor_etc.conversions.convert_electron_flux_mag
      :parser: myst
