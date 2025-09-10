:py:mod:`castor_etc.sources`
============================

.. py:module:: castor_etc.sources

.. autodoc2-docstring:: castor_etc.sources
   :parser: myst
   :allowtitles:

Module Contents
---------------

Classes
~~~~~~~

.. list-table::
   :class: autosummary longtable
   :align: left

   * - :py:obj:`Profiles <castor_etc.sources.Profiles>`
     - .. autodoc2-docstring:: castor_etc.sources.Profiles
          :parser: myst
          :summary:
   * - :py:obj:`Source <castor_etc.sources.Source>`
     - .. autodoc2-docstring:: castor_etc.sources.Source
          :parser: myst
          :summary:
   * - :py:obj:`PointSource <castor_etc.sources.PointSource>`
     - .. autodoc2-docstring:: castor_etc.sources.PointSource
          :parser: myst
          :summary:
   * - :py:obj:`ExtendedSource <castor_etc.sources.ExtendedSource>`
     - .. autodoc2-docstring:: castor_etc.sources.ExtendedSource
          :parser: myst
          :summary:
   * - :py:obj:`GalaxySource <castor_etc.sources.GalaxySource>`
     - .. autodoc2-docstring:: castor_etc.sources.GalaxySource
          :parser: myst
          :summary:
   * - :py:obj:`CustomSource <castor_etc.sources.CustomSource>`
     - .. autodoc2-docstring:: castor_etc.sources.CustomSource
          :parser: myst
          :summary:

API
~~~

.. py:class:: Profiles()
   :canonical: castor_etc.sources.Profiles

   .. autodoc2-docstring:: castor_etc.sources.Profiles
      :parser: myst

   .. rubric:: Initialization

   .. autodoc2-docstring:: castor_etc.sources.Profiles.__init__
      :parser: myst

   .. py:method:: uniform(a=None, b=None, angle=0)
      :canonical: castor_etc.sources.Profiles.uniform
      :staticmethod:

      .. autodoc2-docstring:: castor_etc.sources.Profiles.uniform
         :parser: myst

   .. py:method:: ellipse(a0, b0, angle=0)
      :canonical: castor_etc.sources.Profiles.ellipse
      :staticmethod:

      .. autodoc2-docstring:: castor_etc.sources.Profiles.ellipse
         :parser: myst

   .. py:method:: sersic(r_eff, n=1, e=0, angle=0)
      :canonical: castor_etc.sources.Profiles.sersic
      :staticmethod:

      .. autodoc2-docstring:: castor_etc.sources.Profiles.sersic
         :parser: myst

.. py:class:: Source(profile, init_dimensions=True, check_profile=True)
   :canonical: castor_etc.sources.Source

   Bases: :py:obj:`castor_etc.spectrum.SpectrumMixin`, :py:obj:`castor_etc.spectrum.NormMixin`

   .. autodoc2-docstring:: castor_etc.sources.Source
      :parser: myst

   .. rubric:: Initialization

   .. autodoc2-docstring:: castor_etc.sources.Source.__init__
      :parser: myst

   .. py:method:: copy()
      :canonical: castor_etc.sources.Source.copy

      .. autodoc2-docstring:: castor_etc.sources.Source.copy
         :parser: myst

.. py:class:: PointSource()
   :canonical: castor_etc.sources.PointSource

   Bases: :py:obj:`castor_etc.sources.Source`

   .. autodoc2-docstring:: castor_etc.sources.PointSource
      :parser: myst

   .. rubric:: Initialization

   .. autodoc2-docstring:: castor_etc.sources.PointSource.__init__
      :parser: myst

.. py:class:: ExtendedSource(angle_a=None, angle_b=None, a=None, b=None, dist=None, rotation=0.0, profile='uniform', exponential_scale_lengths=None)
   :canonical: castor_etc.sources.ExtendedSource

   Bases: :py:obj:`castor_etc.sources.Source`

   .. autodoc2-docstring:: castor_etc.sources.ExtendedSource
      :parser: myst

   .. rubric:: Initialization

   .. autodoc2-docstring:: castor_etc.sources.ExtendedSource.__init__
      :parser: myst

.. py:class:: GalaxySource(r_eff, n, axial_ratio, rotation=0.0)
   :canonical: castor_etc.sources.GalaxySource

   Bases: :py:obj:`castor_etc.sources.Source`

   .. autodoc2-docstring:: castor_etc.sources.GalaxySource
      :parser: myst

   .. rubric:: Initialization

   .. autodoc2-docstring:: castor_etc.sources.GalaxySource.__init__
      :parser: myst

.. py:class:: CustomSource(profile_filepath, passband, center=None, px_scale_unit=u.deg)
   :canonical: castor_etc.sources.CustomSource

   Bases: :py:obj:`castor_etc.sources.Source`

   .. autodoc2-docstring:: castor_etc.sources.CustomSource
      :parser: myst

   .. rubric:: Initialization

   .. autodoc2-docstring:: castor_etc.sources.CustomSource.__init__
      :parser: myst
