.. _tutorial:
Tutorial
========

*SeismicMesh* supports the generation of both 2D and 3D meshes in
either serial or parallel from seismic velocity models. It also supports the generation of
complex mesh sizing functions that are relevant to seismological applications.

Here I show how to build a 2D mesh adapted to the BP2004 benchmark model. Nearly the same code extends to 3D with only changes to the data and domain extents. For a 3D example, see the file ``example/example_3D.py``

Distributed memory parallelism can be used by importing ``mpi4py`` and declaring the following three lines near the top of the script (after other imports)::

    from mpi4py import MPI

    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()

Otherwise, the application program interface remains the same between serial and parallel execution. Parallel exeuction takes place by typing::

    mpirun -np N python your_script.py

where `N` is the number of processors (e.g., 2,3,4 etc.)

.. warning::
    Oversubscribing the problem to too many processors will surely lead to problems! Generally, keeping the number of vertices per rank between 10-30k/ rank results in optimal performance.


.. warning::

    Under construction. Contributions very welcome!

Data
--------

Data for this 2D tutorial can be downloaded for the BP2004 benchmark model::

    wget http://s3.amazonaws.com/open.source.geoscience/open_data/bpvelanal2004/vel_z6.25m_x12.5m_exact.segy.gz


For more details about this data, see the following link https://wiki.seg.org/wiki/2004_BP_velocity_estimation_benchmark_model


The original data for the 3D example can be downloaded here::

    wget https://s3.amazonaws.com/open.source.geoscience/open_data/seg_eage_models_cd/salt_and_overthrust_models.tar.gz


For more details about this 3D data, see the following link: https://wiki.seg.org/wiki/SEG/EAGE_Salt_and_Overthrust_Models
The 3D Salt model was used from that file.


Some things to know
---------------------

In order to use these sizing functions, it is assumed you have a seismic velocity model
defined on a structured grid as was mentioned in the overview section.

This seismic velocity model is passed to the *MeshSizeFunction* class along with the domain extents ::

    import SeismicMesh

    fname = "velocity_models/vel_z6.25m_x12.5m_exact.segy"
    bbox = (-12e3, 0, 0, 67e3)

    # Construct mesh sizing object from velocity model
    ef = SeismicMesh.MeshSizeFunction(
        bbox=bbox,
        model=fname,
        other-args-go-here,...
    )

* If the user uses the *MeshSizeFunction* class, the software makes the assumption the domain can be approximated by a rectangle/cube. Thus, the user specifies the domain geometry as a tuple of coordinates in meters representing the corners of the domain::

.. math::

    bbox = (z_{min}, z_{max}, x_{min}, x_{max})

* In 3D::

.. math::

    bbox = (z_{min}, z_{max}, x_{min}, x_{max}, y_{min}, y_{max})`

.. note :: The program automatically generates the rectangle/cube domain geometry used during meshing if a *MeshSizeFunction* object is passed to the generator.


.. warning::

    All of the mesh size functions detailed below assume you pass the :mod:`bbox` and :mod:`fname` to the *MeshSizeFunction* class constructor.

* If the user wants to define a mesh with an irregular boundary (other than a cube), then they will have to use the :class:`SignedDistanceFunctionGenerator` and either pass it
  (x,y,z) scatterd points that define the location of the desired boundary and/or a velocity interval range that represents the region of the domain they want meshed. See the section immeditately below for instructions on how to create and use a custom signed distance function from a seismic velocity model.

Creating your own signed distance functions
-----------------------------------------------

Mesh size function
-------------------------------------------

Given a coordinate in :math:`R^n` where :math:`n= 2,3`, the sizing map returns the desired mesh size :mod:`h`. The mesh sizing capability provides a method to draft new meshes in a consistent and repeatable manner from available seismic velocity models. The sizing map is built on a Cartesian grid, which simplifies implementation details especially in regard to distributed memory parallelism. Furthermore, seismic velocity models are available on structured grids and thus the same grid can be used to build the sizing map on.

The notion of an adequate mesh size is determined by a combination of the physics of acoustic/elastic wave propagation, the desired numerical accuracy of the solution (e.g., spatial polynomial order, timestepping method, etc.), and allowable computational cost of the model. In the following sub-sections, each mesh strategy is described briefly and how to use it by calling the :class:`MeshSizeFunction` class constructor.


Wavelength-to-gridscale
^^^^^^^^^^^^^^^^^^^^^^^
The highest frequency of the source wavelet :math:`f_{max}` and the smallest value of the velocity model :math:`v_{min}` define the shortest scale length of the problem since the shortest spatial wavelength :math:`\lambda_{min}` is equal to the :math:`\frac{v_{min}}{f_{max}}`. For marine domains, :math:`v_{min}` is approximately 1,484 m/s, which is the speed of sound in seawater, thus the finest mesh resolution is near the water layer.

The user is able to specify the number of vertices per wavelength :math:`\alpha_{wl}` the peak source frequency :math:`f_{max}`. This sizing heuristic also  can be used to take into account varying polynomial orders for finite elements. For instance if using quadratic P=2 elements, :math:`\alpha_{wl}` can be safely be set to 5 to avoid excessive dispersion and dissipatation otherwise that would occur with P=1 elements::

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

Seismic domains are known for sharp gradients in material properties, such as seismic velocity. These sharp gradients lead to reflections and refractions in propagated waves, which are critical for successful imaging. Thus, finer mesh resolution can be deployed inversely proportional to the local standard deviation of P-wave velocity. The local standard deviation of seismic P-wave velocity is calculated in a sliding window around each point on the velocity model. The user chooses the mapping relationship between the local standard deviation of the seismic velocity model and the values of the corresponding mesh size nearby it. This parameter is referred to as the :math:`grad` and is specified in meters.
For instance a :math:`grad` of 50, would imply that the largest gradient in seismic P-wave velocity is mapped to a minimum resolution of 50-m.::

    import SeismicMesh

    fname = "velocity_models/vel_z6.25m_x12.5m_exact.segy"
    bbox = (-12e3, 0, 0, 67e3)

    # Construct mesh sizing object from velocity model
    ef = SeismicMesh.MeshSizeFunction(
        bbox=bbox,
        model=fname,
        grad=50, # the desired mesh size in meters near the shaprest gradient in the domain
    )

.. note:

    The mapping of the local standard deviation of the gradient of seismic velocity is normalized to an interval of :math:`[0,1]` so that the largest gradient is assigned the mesh resolution indicated by :math`grad` and all other grad-to-mesh-sizes are associated using a linear relationship (with a slope of 1 and y-intercept of 0).



Courant-Friedrichs-Lewey (CFL) condition
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Almost all numerical wave propagators utilize explicit numerical methods in the seismic domain. The major advantage for these explicit methods is computational speed. However, it is well-known that all explicit or semi-explicit methods require that the Courant number be bounded above by the Courant-Friedrichs-Lewey (CFL) condition. Ignoring this condition will lead to a numerical instability and a useless unstable simulation. Thus, one thing we must be careful of when using the above mesh size functions is that the CFL condition is indeed bounded.

After sizing functions have been activated, a conservative maximum Courant number is enforced.

For the linear acoustic wave equation assuming isotropic mesh resolution, the CFL condition is commonly described by

.. math::

    C_{r}(x) = \frac{(\Delta t*v_p(x))}{dim*h(x)}

where :math:`h` is the diameter of the circumball that inscribes the element either calculated from :math:`f(h)` or from the actual mesh cells, :math:`dim` is the spatial dimension of the problem (2 or 3), :math:`\Delta t` is the intended simulation time step in seconds and :math:`v_p` is the local seismic P-wave velocity. The above equation can be rearranged to find the minimum mesh size possible for a given :math:`v_p` and :math:`\Delta t`, based on some user-defined value of :math:`Cr \leq 1`. If there are any violations of the CFL, they can bed edited before building the mesh so to satisfy that the maximum :math:`Cr` is less than some conservative threshold. We normally apply :math:`Cr = 0.5`, which provides a solid buffer but this can but this can be controlled by the user like the following::

    import SeismicMesh
    fname = "velocity_models/vel_z6.25m_x12.5m_exact.segy"
    bbox = (-12e3, 0, 0, 67e3)

    # Construct mesh sizing object from velocity model
    ef = SeismicMesh.MeshSizeFunction(
        bbox=bbox,
        model=fname,
        cr=0.5, # maximum bounded Courant number to be bounded in the mesh sizing function
        dt=0.001, # for the given :math:`\Delta t` of 0.001 seconds
        ...
    )

Further, the space order of the method (:math:`p`) can also be incorporated into the above formula to consider the higher spatial order that the simulation will use::

    ef = SeismicMesh.MeshSizeFunction(
        bbox=bbox,
        model=fname,
        cr=0.5, # maximum bounded Courant number :math:`Cr_{max}` in the mesh
        dt=0.001, # for the given :math:`\Delta t` of 0.001 seconds
        space_order = 2, # assume quadratic elements :math:`P=2`
        ...
    )

The above code implies that the mesh will be used in a simulation with :math:`P=2` quadratic elements, and thus will ensure the :math:`Cr_{max}` is divided by :math:`\frac{1}{space\_order}`


Mesh size gradation
^^^^^^^^^^^^^^^^^^^^^^^

In regions where there are sharp material contrasts, the variation in element size can become substantially large, especially using the aforementioned sizing strategies such as the wavelength-to-gridscale. Attempting to construct a mesh with such large spatial variations in mesh sizes would result in low-geometric quality elements that compromise the numerical stability of a model.

Thus, the final stage of the development of a mesh size function :math:`h(x)` involves ensuring a size smoothness limit, :math:`g` such that for any two points :math:`x_i`, :math:`x_j`, the local increase in size is bounded such as:

 :math:`h(\boldsymbol{x_j}) \leq h(\boldsymbol{x_i}) + \alpha_g||\boldsymbol{x_i}-\boldsymbol{x_j}||`

A smoothness criteria is necessary to produce a mesh that can simulate physical processes with a practical time step as sharp gradients in mesh resolution typically lead to highly skewed angles that result in poor numerical performance.

We adopt the method to smooth the mesh size function originally proposed by [grading]_. A smoother sizing function is congruent with a higher overall element quality but with more triangles in the mesh. Generally, setting :math:`0.2 \leq \alpha_g \leq 0.3` produces good results::

   import SeismicMesh
   fname = "velocity_models/vel_z6.25m_x12.5m_exact.segy"
   bbox = (-12e3, 0, 0, 67e3)

   # Construct mesh sizing object from velocity model
   ef = SeismicMesh.MeshSizeFunction(
       bbox=bbox,
       model=fname,
       grade=0.15, # :math:`g` cell-to-cell size rate growth bound
       ...
   )

Domain extension
^^^^^^^^^^^^^^^^^^^

.. note::

    It is assumed that the top side of the domain represents the free-surface thus no domain extension is applied there.

In seismology applications, the goal is often to model the propagation of an elastic or acoustic wave through an infinite domain. However, this is obviously not possible so the domain is approximated by a finite region of space. This can lead to undeseriable artifical reflections off the sides of the domain however. A common approach to avoid these artifical reflections is to extend the domain and enforce abosrbing boundary conditions in this extension. In terms of meshing to take this under consideration, the user has the option to specify a domain extension of variable width on all three sides of the domain like so::

   import SeismicMesh
   fname = "velocity_models/vel_z6.25m_x12.5m_exact.segy"
   bbox = (-12e3, 0, 0, 67e3)

   # Construct mesh sizing object from velocity model
   ef = SeismicMesh.MeshSizeFunction(
       bbox=bbox,
       model=fname,
       domain_extension=250, # domain will be extended by 250-m on all three sides
       ...
   )

In this domain extension region, mesh resolution can be adapted according to following three different styles.

 * ``Linear`` - extends the seismic velocities on the edges of the domain linearly into the domain extension.

 * ``Constant`` - places a constant velocity of 4000 m/s in the domain extension.

 * ``Edge`` - reflects the seismic velocity about the domain boundary so that velocity profile is symmetric w.r.t domain boudnaries.

An example of the ``edge`` style is below::

   # Construct mesh sizing object from velocity model
   ef = SeismicMesh.MeshSizeFunction(
       bbox=bbox,
       model=fname,
       domain_extension=250, # domain will be extended by 250-m on all three sides
       padstyle="edge", # velocity will be reflected about the edges of the domain
       ...
   )

.. note::

    In our experience, the ``edge`` option works the best at reducing reflections with absorbing boundary conditions.



Mesh generation
-------------------------------------------

.. warning:
    Connectivity is made approximately deterministic as each instance of mesh generation uses
    the same ``seed=0``.

2D Mesh generation
^^^^^^^^^^^^^^^^^^^^^^^


3D Mesh generation
^^^^^^^^^^^^^^^^^^^^^^^


Mesh improvement
-------------------------------------------

3D *Sliver* removal
^^^^^^^^^^^^^^^^^^^^^^^




References
______________

.. [grading] Persson, Per-Olof. "Mesh size functions for implicit geometries and PDE-based gradient limiting."
                Engineering with Computers 22.2 (2006): 95-109.

.. [firedrake] Florian Rathgeber, David A. Ham, Lawrence Mitchell, Michael Lange, Fabio Luporini, Andrew T. T. Mcrae, Gheorghe-Teodor Bercea, Graham R. Markall, and Paul H. J. Kelly. Firedrake: automating the finite element method by composing abstractions. ACM Trans. Math. Softw., 43(3):24:1–24:27, 2016. URL: http://arxiv.org/abs/1501.01809, arXiv:1501.01809, doi:10.1145/2998441.
