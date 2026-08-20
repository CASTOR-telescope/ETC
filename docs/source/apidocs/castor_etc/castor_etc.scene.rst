:py:mod:`castor_etc.scene`
==========================

.. py:module:: castor_etc.scene

.. autodoc2-docstring:: castor_etc.scene
   :parser: myst
   :allowtitles:

Module Contents
---------------

Classes
~~~~~~~

.. list-table::
   :class: autosummary longtable
   :align: left

   * - :py:obj:`Scene <castor_etc.scene.Scene>`
     - .. autodoc2-docstring:: castor_etc.scene.Scene
          :parser: myst
          :summary:

API
~~~

.. py:class:: Scene(telescope)
   :canonical: castor_etc.scene.Scene

   .. autodoc2-docstring:: castor_etc.scene.Scene
      :parser: myst

   .. rubric:: Initialization

   .. autodoc2-docstring:: castor_etc.scene.Scene.__init__
      :parser: myst

   .. py:method:: addSource(source, mag=None, delta_x=0, delta_y=0, source_name='source_1')
      :canonical: castor_etc.scene.Scene.addSource

      .. autodoc2-docstring:: castor_etc.scene.Scene.addSource
         :parser: myst

   .. py:method:: Gaussian2D(x, y, sigma, a=1.0, x0=0.0, y0=0.0)
      :canonical: castor_etc.scene.Scene.Gaussian2D
      :staticmethod:

      .. autodoc2-docstring:: castor_etc.scene.Scene.Gaussian2D
         :parser: myst

   .. py:method:: displaySource(telescope, passband, source_index, exptime=1.0)
      :canonical: castor_etc.scene.Scene.displaySource

      .. autodoc2-docstring:: castor_etc.scene.Scene.displaySource
         :parser: myst

   .. py:method:: displayScene(telescope, passband, exptime=1.0)
      :canonical: castor_etc.scene.Scene.displayScene

      .. autodoc2-docstring:: castor_etc.scene.Scene.displayScene
         :parser: myst
