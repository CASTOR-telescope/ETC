:py:mod:`castor_etc.tests.test_energy`
======================================

.. py:module:: castor_etc.tests.test_energy

.. autodoc2-docstring:: castor_etc.tests.test_energy
   :allowtitles:

Module Contents
---------------

Classes
~~~~~~~

.. list-table::
   :class: autosummary longtable
   :align: left

   * - :py:obj:`TestEnergy <castor_etc.tests.test_energy.TestEnergy>`
     - .. autodoc2-docstring:: castor_etc.tests.test_energy.TestEnergy
          :summary:

Data
~~~~

.. list-table::
   :class: autosummary longtable
   :align: left

   * - :py:obj:`_TOL <castor_etc.tests.test_energy._TOL>`
     - .. autodoc2-docstring:: castor_etc.tests.test_energy._TOL
          :summary:

API
~~~

.. py:data:: _TOL
   :canonical: castor_etc.tests.test_energy._TOL
   :value: 1e-15

   .. autodoc2-docstring:: castor_etc.tests.test_energy._TOL

.. py:class:: TestEnergy(methodName='runTest')
   :canonical: castor_etc.tests.test_energy.TestEnergy

   Bases: :py:obj:`unittest.TestCase`

   .. autodoc2-docstring:: castor_etc.tests.test_energy.TestEnergy

   .. rubric:: Initialization

   .. autodoc2-docstring:: castor_etc.tests.test_energy.TestEnergy.__init__

   .. py:method:: test_calc_photon_energy_wavelength_single()
      :canonical: castor_etc.tests.test_energy.TestEnergy.test_calc_photon_energy_wavelength_single

      .. autodoc2-docstring:: castor_etc.tests.test_energy.TestEnergy.test_calc_photon_energy_wavelength_single

   .. py:method:: test_calc_photon_energy_wavelength_multi()
      :canonical: castor_etc.tests.test_energy.TestEnergy.test_calc_photon_energy_wavelength_multi

      .. autodoc2-docstring:: castor_etc.tests.test_energy.TestEnergy.test_calc_photon_energy_wavelength_multi

   .. py:method:: test_calc_photon_energy_frequency_single()
      :canonical: castor_etc.tests.test_energy.TestEnergy.test_calc_photon_energy_frequency_single
      :abstractmethod:

      .. autodoc2-docstring:: castor_etc.tests.test_energy.TestEnergy.test_calc_photon_energy_frequency_single

   .. py:method:: test_calc_photon_energy_frequency_multi()
      :canonical: castor_etc.tests.test_energy.TestEnergy.test_calc_photon_energy_frequency_multi
      :abstractmethod:

      .. autodoc2-docstring:: castor_etc.tests.test_energy.TestEnergy.test_calc_photon_energy_frequency_multi
