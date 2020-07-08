.. _tutorial:

Tutorial
========

.. warning::

   Under construction. Contributions very welcome!

*SeismicMesh* supports the generation of both 2D and 3D meshes in
either serial or parallel. It also supports the generation of
complex mesh sizing function that are relevant to Seismology.

In order to use these sizing functions, it is assumed you have a seismic velocity model
defined on a structured grid. This seismic velocity model is passed to the *MeshSizeFunction*
class along with the domain extents ::

    import SeismicMesh

    fname = "velocity_models/vel_z6.25m_x12.5m_exact.segy"
    bbox = (-12e3, 0, 0, 67e3)

    # Construct mesh sizing object from velocity model
    ef = SeismicMesh.MeshSizeFunction(
        bbox=bbox,
        model=fname,
        other-args-go-here,...
    )



.. warning::

    All of the mesh size functions detailed below assume you pass the :mod:`bbox` and :mod:`fname` to the *MeshSizeFunction* class constructor.

Mesh size function
-------------------------------------------

Given a coordinate in :math:`R^n` where :math:`n= 2,3`, the sizing map returns the desired mesh size :mod:`h`. The mesh sizing capability provides a method to draft new meshes in a consistent and repeatable manner. The sizing map is built on a Cartesian grid, which simplifies implementation details especially in regard to distributed memory parallelism.  Furthermore, seismic velocity models are available on structured grids and thus the same grid can be used to build the sizing map on.The notion of an adequate mesh size is determined by a combination of the physics of acoustic/elastic wave propagation, the desired numerical accuracy of the solution (e.g., polynomial order,timestepping method, etc.), and the computational cost of the model. In the following sub-sections,each mesh adaptation strategy is described briefly and how to use it.


Wavelength-to-gridscale
^^^^^^^^^^^^^^^^^^^^^^^
The highest frequency of the source wavelet :math:`f_{max}` and the smallest value of the velocity model :math:`v_{min}` define the shortest scale length of the problem since the shortest spatial wavelength :math:`\lambda_{min}` is equal to the :math:`\frac{v_{min}}{f_{max}}`. For marine domains, :math:`v_{min}` is approximately 1,484 m/s, which is the speed of sound in seawater, thus the finest mesh resolution is near the water layer.

The user is able to specify the number of vertices per wavelength :math:`\alpha_{wl}` the peak source frequency :math:`f_{max}`.  This sizing heuristic also  can be used to take into account varying polynomial orders for finite elements. For instance if using quadratic P=2 elements, :math:`\alpha_{wl}` can be safely  be set to 5 to avoid excessive dispersion and dissipatation::

   import SeismicMesh
   fname = "velocity_models/vel_z6.25m_x12.5m_exact.segy"
   bbox = (-12e3, 0, 0, 67e3)

   # Construct mesh sizing object from velocity model
   ef = SeismicMesh.MeshSizeFunction(
       bbox=bbox,
       model=fname,
       freq=2, # maximum source frequency
       wl=3, # :math:`\alpha_{wl}` number of grid points per wavelength
   )



Resolving seismic velocity gradients
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


Courants-Friedrichs-Lewey (CFL) condition
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* Acoustic wave equation


* Elastic wave equation



Mesh size gradation
^^^^^^^^^^^^^^^^^^^^^^^


Mesh generation
-------------------------------------------

.. warning:
    Results can be made fully deterministic by specifying the argument `seed=0` to the generator. This ensures that all
    stochastic operations will be repeated in the same way as the random number used as the `seed` is fixed.

.. note:
    Parallelism is activated by passing the :mod:`COMM` to the *MeshSizeFunction* constructor ::

  ef = ef.build(comm=COMM)
  ef = ef.construct_lambdas(COMM)

    and also to the *MeshGenerator* constructor ::

  mshgen = SeismicMesh.MeshGenerator(ef, method="cgal")
  points, cells = mshgen.build(COMM=COMM)


2D Mesh generation
^^^^^^^^^^^^^^^^^^^^^^^


3D Mesh generation
^^^^^^^^^^^^^^^^^^^^^^^


Mesh improvement
-------------------------------------------

3D *Sliver* removal
^^^^^^^^^^^^^^^^^^^^^^^
