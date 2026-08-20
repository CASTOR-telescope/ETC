:py:mod:`castor_etc.detector`
=============================

.. py:module:: castor_etc.detector

.. autodoc2-docstring:: castor_etc.detector
   :parser: myst
   :allowtitles:

Module Contents
---------------

Classes
~~~~~~~

.. list-table::
   :class: autosummary longtable
   :align: left

   * - :py:obj:`TempDependentNoise <castor_etc.detector.TempDependentNoise>`
     -
   * - :py:obj:`DarkCurrentNoise <castor_etc.detector.DarkCurrentNoise>`
     -
   * - :py:obj:`ReadoutNoise <castor_etc.detector.ReadoutNoise>`
     -

Data
~~~~

.. list-table::
   :class: autosummary longtable
   :align: left

   * - :py:obj:`DARK_CURRENT_PROFILE <castor_etc.detector.DARK_CURRENT_PROFILE>`
     - .. autodoc2-docstring:: castor_etc.detector.DARK_CURRENT_PROFILE
          :parser: myst
          :summary:
   * - :py:obj:`READOUT_NOISE_PROFILE <castor_etc.detector.READOUT_NOISE_PROFILE>`
     - .. autodoc2-docstring:: castor_etc.detector.READOUT_NOISE_PROFILE
          :parser: myst
          :summary:

API
~~~

.. py:class:: TempDependentNoise(profile: dict)
   :canonical: castor_etc.detector.TempDependentNoise

   Bases: :py:obj:`castor_etc.BaseClass`

   .. py:attribute:: _profile
      :canonical: castor_etc.detector.TempDependentNoise._profile
      :type: dict
      :value: None

      .. autodoc2-docstring:: castor_etc.detector.TempDependentNoise._profile
         :parser: myst

   .. py:method:: get_mean_value(target_temp: float)
      :canonical: castor_etc.detector.TempDependentNoise.get_mean_value

      .. autodoc2-docstring:: castor_etc.detector.TempDependentNoise.get_mean_value
         :parser: myst

   .. py:method:: get_median_value(target_temp: float)
      :canonical: castor_etc.detector.TempDependentNoise.get_median_value

      .. autodoc2-docstring:: castor_etc.detector.TempDependentNoise.get_median_value
         :parser: myst

   .. py:method:: get_peak_value(target_temp: float)
      :canonical: castor_etc.detector.TempDependentNoise.get_peak_value

      .. autodoc2-docstring:: castor_etc.detector.TempDependentNoise.get_peak_value
         :parser: myst

.. py:data:: DARK_CURRENT_PROFILE
   :canonical: castor_etc.detector.DARK_CURRENT_PROFILE
   :value: None

   .. autodoc2-docstring:: castor_etc.detector.DARK_CURRENT_PROFILE
      :parser: myst

.. py:class:: DarkCurrentNoise(profile=DARK_CURRENT_PROFILE)
   :canonical: castor_etc.detector.DarkCurrentNoise

   Bases: :py:obj:`castor_etc.detector.TempDependentNoise`

.. py:data:: READOUT_NOISE_PROFILE
   :canonical: castor_etc.detector.READOUT_NOISE_PROFILE
   :value: None

   .. autodoc2-docstring:: castor_etc.detector.READOUT_NOISE_PROFILE
      :parser: myst

.. py:class:: ReadoutNoise(profile=READOUT_NOISE_PROFILE)
   :canonical: castor_etc.detector.ReadoutNoise

   Bases: :py:obj:`castor_etc.detector.TempDependentNoise`
