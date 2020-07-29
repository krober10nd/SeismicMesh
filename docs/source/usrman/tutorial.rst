.. _tutorial:

.. warning::

    Under construction. Contributions very welcome!

Tutorial
========

*SeismicMesh* supports the generation of both 2D and 3D meshes in
either serial or parallel from seismic velocity models. It also supports the generation of
complex mesh sizing functions that are relevant to seismological applications.

Here I show how to build meshes from sizing functions created with the software and explain what the options mean. Nearly the same code extends to 3D with only changes to the data and domain extents. For a 3D example, see the file ``example/example_3D.py``

Distributed memory parallelism can be used by first importing ``mpi4py`` and declaring the following three lines near the top of the script (after your other imports)::

    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()

Otherwise, the application program interface remains the same between serial and parallel execution in both 2D and 3D. Parallel exeuction takes place by typing::

    mpirun -np N python your_script.py

where `N` is the number of processors (e.g., 2,3,4 etc.)

.. warning::
    Oversubscribing the problem to too many processors will surely lead to problems and slow downs. Generally, keeping the number of vertices per rank between 10-30k/ rank results in optimal performance.



Data for examples
-------------------

Data for this 2D tutorial can be downloaded for the BP2004 benchmark model::

    wget http://s3.amazonaws.com/open.source.geoscience/open_data/bpvelanal2004/vel_z6.25m_x12.5m_exact.segy.gz

and for the irregular free surface example here::

    wget http://s3.amazonaws.com/open.source.geoscience/open_data/bpmodel94/velocity.segy.gz

For more details about these two models, see the following link::

    https://wiki.seg.org/wiki/2004_BP_velocity_estimation_benchmark_model
    https://wiki.seg.org/wiki/1994_BP_migration_from_topography

The original data for the 3D example can be downloaded here::

    wget https://s3.amazonaws.com/open.source.geoscience/open_data/seg_eage_models_cd/salt_and_overthrust_models.tar.gz


For more details about this 3D data, see the following link: https://wiki.seg.org/wiki/SEG/EAGE_Salt_and_Overthrust_Models
The 3D Salt model was used from that archive.

The following procedure was used to turn the original binary data into a format that was usable for *SeismicMesh*::

    import zipfile
    import numpy as np
    # Dimensions
    nx, ny, nz = 676, 676, 210

    # Extract Saltf@@ from SALTF.ZIP
    zipfile.ZipFile('./data/SALTF.ZIP', 'r').extract('Saltf@@', path='./data/')

    # Load data
    with open('./data/Saltf@@', 'r') as file:
       v = np.fromfile(file, dtype=np.dtype('float32').newbyteorder('>'))
       v = v.reshape(nx, ny, nz, order='F')
    # Write the v to a binary file
    file = open("EAGE_Salt.bin", "wb")
    file.write(v)
    file.close()


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

* The user specifies the domain extents as a tuple of coordinates in meters representing the corners of the domain::

.. math::

    bbox = (z_{min}, z_{max}, x_{min}, x_{max})

* In 3D::

.. math::

    bbox = (z_{min}, z_{max}, x_{min}, x_{max}, y_{min}, y_{max})`

.. note :: The program automatically generates the rectangle/cube domain geometry used during meshing if a *MeshSizeFunction* object is passed to the generator.


.. warning::

    All of the mesh size functions detailed below assume you pass the :mod:`bbox` and :mod:`fname` to the *MeshSizeFunction* class constructor.

* If the user wants to define a mesh with an irregular boundary (other than a cube), then they will have to use the :class:`SignedDistanceFunctionGenerator` and pass a velocity interval range that represents the region of the domain they want meshed. See the section immeditately below for instructions on how to create and use a custom signed distance function from a seismic velocity model.

Creating your own signed distance functions
-----------------------------------------------

Using the :class:`SignedDistanceFunctionGenerator` tool, the user can threshold their velocity model to produce a signed distance function which can then be used to mesh with. For example, in the case of an irregular free surface boundary in the 2D Foothills model, a SDF can be created by thresholding a modified velocity model. Specifically, we mesh the region of the domain with a P-wave velocity greater than 4,000 m/s::

     # Bulid a signed distance function from the seismic velocity model
     # Some pockets of velocity < 4000 exist, fill those in.
     vp2 = ef.vp.copy()
     vp2 = np.where(vp2 < 4000, 4001, vp2)
     SDF = SdfGen(
         bbox=bbox, field=vp2, min_threshold=4000.0, gridspacing=(10.0, 15.0),
     ).SDF

.. note :: We modify the velocity model so that some small pockets of area with velocity less than 4,000 m/s don't appear as holes in the mesh when we use 4,000 m/s to create the signed distance function.

The output of the call to the ``SdfGen`` gives a handle to a function object, which can be passed to the mesh generator as detailed below in the section describing how to call the mesh generator.

Images shown below are the seismic P-wave velocity model for the foothills example and the resulting SDF from exeucting the code directly above. In this second image below, the color yellow indicates the region to-be-meshed and the greenish color indicates the region outside of the domain.

.. image:: Foothills.png

.. image:: ExampleOfSDF.png

Mesh size function
-------------------------------------------

Given a coordinate in :math:`R^n` where :math:`n= 2,3`, the sizing function returns the desired mesh size :mod:`h` near to that point. The mesh sizing capability provides is a convenience class that helps draft new meshes in a consistent and repeatable manner directly from available seismic velocity models. The sizing map is built on a Cartesian grid, which simplifies implementation details especially in regard to distributed memory parallelism. Furthermore, seismic velocity models are available on structured grids and thus the same grid can be used to build the sizing map on.

.. note:
    Seismic velocity models often have constant grid spacing in each dimension. The software considers this automatically.

The notion of an adequate mesh size is determined by a combination of the physics of acoustic/elastic wave propagation, the desired numerical accuracy of the solution (e.g., spatial polynomial order, timestepping method, etc.), and allowable computational cost of the model amongst other things. In the following sub-sections, each available mesh strategy is briefly described and psuedo-code regarding how to call the :class:`MeshSizeFunction` class constructor.

.. note :: The final mesh size map is taken as the minimum of all supplied sizing functions.

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
For instance a :math:`grad` of 50 would imply that the largest gradient in seismic P-wave velocity is mapped to a minimum resolution of 50-m.::

    import SeismicMesh

    fname = "velocity_models/vel_z6.25m_x12.5m_exact.segy"
    bbox = (-12e3, 0, 0, 67e3)

    # Construct mesh sizing object from velocity model
    ef = SeismicMesh.MeshSizeFunction(
        bbox=bbox,
        model=fname,
        grad=50, # the desired mesh size in meters near the shaprest gradient in the domain
    )

.. image:: SlopeStrat3D.jpg

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

.. image:: ExGrade3D.jpg

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

 * ``Constant`` - places a constant velocity of the users selection in the domain extension.

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

.. image:: domainext.png



Mesh generation
-------------------------------------------

.. warning:
    Connectivity is made approximately deterministic as each instance of mesh generation uses
    the same ``seed=0``. The user can specify the seed if they like.

The user generates the mesh in either two ways: The user passes a mesh size function object like this::

    # Construct a mesh generator object
    mshgen = SeismicMesh.MeshGenerator(ef)

And then they call the ``build`` method specifying the number of iterations they want the generator to perform::

    # Build the mesh
    points, cells = mshgen.build(max_iter=75, axis=1)

.. note :: Generally setting max_iter to between 50 to 100 iterations works best. By default it runs 50 iterations.

.. note :: For parallel exeuction, the user can choose which axis (0, 1, or 2 (if 3D)) to decompose the domain.

Or, the second way the user specified their own mesh size function ``f(h)`` and/or ``f(d)``::

    # Example of a signed distance function for a cylinder in unit space.
    def fd(p):
        # sizing function of a cylinder in [-1., -1., -.1] x [1., 1., 1.]
        r, z = np.sqrt(p[:, 0] ** 2 + p[:, 1] ** 2), p[:, 2]
        d1, d2, d3 = r - 1.0, z - 1.0, -z - 1.0
        d4, d5 = np.sqrt(d1 ** 2 + d2 ** 2), np.sqrt(d1 ** 2 + d3 ** 2)
        d = np.maximum.reduce([d1, d2, d3])
        ix = (d1 > 0) * (d2 > 0)
        d[ix] = d4[ix]
        ix = (d1 > 0) * (d3 > 0)
        d[ix] = d5[ix]
        return d


    # Example of an uniform resolution sizing function
    def fh(p):
        # note for parallel execution this logic is required
        # since the decomposition of the sizing function passes a tuple to fh
        if type(p) == tuple:
            h = np.zeros_like(p[0]) + hmin
        else:
            h = np.zeros_like(p) + hmin
        return h


    # Construct mesh generator
    mshgen = SeismicMesh.MeshGenerator(hmin=hmin, bbox=bbox, fd=fd, fh=fh)

.. warning :: In this second way, the user *must* pass the arguments `hmin` the desired minimum resolution in the domain and the tuple of domain extents `bbox` to the mesh generator class constructor.

And then they call the ``build`` method specifying the number of iterations they want the generator to perform::

    # Build the mesh
    points, cells = mshgen.build(max_iter=50)

Mesh improvement
-------------------------------------------


3D *Sliver* removal
^^^^^^^^^^^^^^^^^^^^^^^

If the intended usage of the mesh is for numerical simulation, it is strongly encouraged to execute the sliver removal method after generating a mesh in 3D like so::

    points, cells = mshgen.build(
        points=points, mesh_improvement=True, max_iter=50, min_dh_bound=5,
    )

Note that here we pass it the points from the previous call to build and specify the flag ``mesh_improvement`` to *True*. The option ``min_dh_bound`` represents the target lower bound for the dihedral angle. By default, ``min_dh_angle`` is set to :math:`10`.  The sliver removal algorithm will attempt 50 iterations but will terminate earlier if no slivers are detected.

.. warning:: Do not set the minimum dihedral angle bound greater than 15 unless you've already succesfully ran the mesh with a lower threshold. Otherwise, the method will likely not converge.


References
______________

.. [grading] Persson, Per-Olof. "Mesh size functions for implicit geometries and PDE-based gradient limiting."
                Engineering with Computers 22.2 (2006): 95-109.

.. [firedrake] Florian Rathgeber, David A. Ham, Lawrence Mitchell, Michael Lange, Fabio Luporini, Andrew T. T. Mcrae, Gheorghe-Teodor Bercea, Graham R. Markall, and Paul H. J. Kelly. Firedrake: automating the finite element method by composing abstractions. ACM Trans. Math. Softw., 43(3):24:1–24:27, 2016. URL: http://arxiv.org/abs/1501.01809, arXiv:1501.01809, doi:10.1145/2998441.
