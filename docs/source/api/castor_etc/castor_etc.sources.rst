:py:mod:`castor_etc.sources`
============================

.. py:module:: castor_etc.sources

.. autodoc2-docstring:: castor_etc.sources
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
          :summary:
   * - :py:obj:`Source <castor_etc.sources.Source>`
     - .. autodoc2-docstring:: castor_etc.sources.Source
          :summary:
   * - :py:obj:`PointSource <castor_etc.sources.PointSource>`
     - .. autodoc2-docstring:: castor_etc.sources.PointSource
          :summary:
   * - :py:obj:`ExtendedSource <castor_etc.sources.ExtendedSource>`
     - .. autodoc2-docstring:: castor_etc.sources.ExtendedSource
          :summary:
   * - :py:obj:`GalaxySource <castor_etc.sources.GalaxySource>`
     - .. autodoc2-docstring:: castor_etc.sources.GalaxySource
          :summary:
   * - :py:obj:`CustomSource <castor_etc.sources.CustomSource>`
     - .. autodoc2-docstring:: castor_etc.sources.CustomSource
          :summary:

API
~~~

.. py:class:: Profiles()
   :canonical: castor_etc.sources.Profiles

   .. autodoc2-docstring:: castor_etc.sources.Profiles

   .. rubric:: Initialization

   .. autodoc2-docstring:: castor_etc.sources.Profiles.__init__

   .. py:method:: uniform(a=None, b=None, angle=0)
      :canonical: castor_etc.sources.Profiles.uniform
      :staticmethod:

      .. autodoc2-docstring:: castor_etc.sources.Profiles.uniform

   .. py:method:: ellipse(a0, b0, angle=0)
      :canonical: castor_etc.sources.Profiles.ellipse
      :staticmethod:

      .. autodoc2-docstring:: castor_etc.sources.Profiles.ellipse

   .. py:method:: sersic(r_eff, n=1, e=0, angle=0)
      :canonical: castor_etc.sources.Profiles.sersic
      :staticmethod:

      .. autodoc2-docstring:: castor_etc.sources.Profiles.sersic

.. py:class:: Source(profile, init_dimensions=True, check_profile=True)
   :canonical: castor_etc.sources.Source

   Bases: :py:obj:`castor_etc.spectrum.SpectrumMixin`, :py:obj:`castor_etc.spectrum.NormMixin`

   .. autodoc2-docstring:: castor_etc.sources.Source

   .. rubric:: Initialization

   .. autodoc2-docstring:: castor_etc.sources.Source.__init__

   .. py:method:: copy()
      :canonical: castor_etc.sources.Source.copy

      .. autodoc2-docstring:: castor_etc.sources.Source.copy

.. py:class:: PointSource()
   :canonical: castor_etc.sources.PointSource

   Bases: :py:obj:`castor_etc.sources.Source`

   .. autodoc2-docstring:: castor_etc.sources.PointSource

   .. rubric:: Initialization

   .. autodoc2-docstring:: castor_etc.sources.PointSource.__init__

.. py:class:: ExtendedSource(angle_a=None, angle_b=None, a=None, b=None, dist=None, rotation=0.0, profile='uniform', exponential_scale_lengths=None)
   :canonical: castor_etc.sources.ExtendedSource

   Bases: :py:obj:`castor_etc.sources.Source`

   .. autodoc2-docstring:: castor_etc.sources.ExtendedSource

   .. rubric:: Initialization

   .. autodoc2-docstring:: castor_etc.sources.ExtendedSource.__init__

.. py:class:: GalaxySource(r_eff, n, axial_ratio, rotation=0.0)
   :canonical: castor_etc.sources.GalaxySource

   Bases: :py:obj:`castor_etc.sources.Source`

   .. autodoc2-docstring:: castor_etc.sources.GalaxySource

   .. rubric:: Initialization

   .. autodoc2-docstring:: castor_etc.sources.GalaxySource.__init__

.. py:class:: CustomSource(profile_filepath, passband, center=None, px_scale_unit=u.deg)
   :canonical: castor_etc.sources.CustomSource

   Bases: :py:obj:`castor_etc.sources.Source`

   .. autodoc2-docstring:: castor_etc.sources.CustomSource

   .. rubric:: Initialization

   .. autodoc2-docstring:: castor_etc.sources.CustomSource.__init__
