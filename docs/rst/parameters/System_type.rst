System_type
===========
The parameter is provides the program with a broad
overview of the type of system that will be simulated, and is used
by Python to initialize certain variable, and to control what variables
are asked for later.

Type
  Enumerator

Values
  star
    System in which the central object is a star

  cv
    System with a secondary star, which can occult the central object and disk depending on phase

  bh
    System with a black hole binary

  agn
    AGN

  previous
    In this case, one is starting from a previous run with python, and one want to either continue the
    run or change some parameters associated with radiation sources


File
  `python.c <https://github.com/agnwinds/python/blob/master/source/python.c>`_


Child(ren)
  * :ref:`Binary.mass_sec`

  * :ref:`Central_object.luminosity`

  * :ref:`Wind.number_of_components`

  * :ref:`Spectrum.orbit_phase`

  * :ref:`Central_object.temp`

  * :ref:`Binary.period`

  * :ref:`Atomic_data`

  * :ref:`Central_object.geometry_for_source`

  * :ref:`Wind.old_windfile`

  * :ref:`Central_object.blackbody_temp`

  * :ref:`Boundary_layer.radiation`

