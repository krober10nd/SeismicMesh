Overview
========

This software aims to create end-to-end workflows (e.g., from seismic velocity models to mesh) to build quality two- and three-dimensional (2D and 3D) unstructured triangular and tetrahedral meshes for seismic domains suitable for numerical wave propagators. A focus is placed on parallel unstructured mesh generation and adaptation capabilities for numerical simulations using the finite element method. The capabilities of the software are exposed to the user via a simple application program interfaces (API) written in Python. The software can be scaled up on distributed memory clusters.

*SeismicMesh* has been used to generate meshes in 2D and 3D acoustic and elastic wave propagators written in *Firedrake's* [firedrake]_, which are used in Full Waveform Inversion, Reverse Time Migration, and Time Travel Tomography applications.


Software architecture
-------------------------------------------

The software is implemented in a mixed language environment (Python and C++). The Python language is used for the API while computationally expensive operations are performed in C++. The two languages are linked together with *pybind11* and installation is handled using *cmake*. The Computational Geometry Algorithms Library [cgal]_ is used to perform geometric operations that use floating point arithmetic to avoid numerical precision issues. Besides this, several common Python packages: *Numpy*, *Scipy*, *MeshIO*, *SegyIO*, and *MPI4py* are used. In particular, *MeshIO* enables to the program to output meshes in dozens of commonly used different formats.

A simple mesh generation workflow consists of two steps: first calling the *MeshSizeFunction* class followed by calling the *MeshGenerator* class.

*MeshSizeFunction*
^^^^^^^^^^^^^^^^^^^^^^^

This class creates an isotropic mesh size function and defines the domain through a signed distance function. The user passes a number of arguments to the class based on the intended application of the mesh.

*MeshGeneration*
^^^^^^^^^^^^^^^^^^^^^^^

This class takes as input a signed distance function and mesh sizing function OR a *MeshSizeFunction* object. It constructs a modified mesh outputting points and simplices, while preserving domain geometry and the mesh sizing distribution. It supports both serial and distributed memory parallel execution. Additionally, it has support a mesh improvement strategy [slivers]_ aimed at removing degenerate elements that are ubiqituous in 3D Delaunay mesh generation methods.


Inputs
-------------------------------------------

* The only required input file to generate a mesh in 2D is a binary file containing the velocity data on a structured grid. In 2D, the SEG-y format containing the seismic velocities of the domain is used. To store the seismic velocity model as a SEG-y file (if it isn't already in this format), the traces represent columns of the seismic velocity model. In 3D, the seismic velocity model is stored as a binary file as well but we do not use the SEG-y format. Instead in 3D, data is stored contiguously in memory in the format z,x,y following the little-endian format. The user must specify how the data will be reshaped in memory by passing the number of rows, columns in the x-direction, and columns in the y-direction.


* If the user uses the *MeshSizeFunction* class, the software makes the assumption the domain can be approximated by a rectangle/cube. Thus, the user specifies the domain geometry as a tuple of coordinates in meters::
    :math:`bbox = (z_{min}, z_{max}, x_{min}, x_{max})`

* In 3D::
    :math:`bbox = (z_{min}, z_{max}, x_{min}, x_{max}, y_{min}, y_{max})`

.. note :: The program automatically generates the rectangle/cube domain geometry used during meshing if a *MeshSizeFunction* object is passed to the generator.


*DistMesh*
-------------------------------------------

For the generation of triangular meshes in 2D and 3D, we use a modified version of the *DistMesh* algorithm [distmesh]_. The algorithm is both simple (approximately 30 lines of code) and useful as it can produce meshes in N-dimensional space, on parametric surfaces in :math:`R^3`, and can enforce point constraints. By using our approach to produce mesh size functions, the mesh generation algorithm is capable of producing high-geometric quality meshes, faithful to our target sizing fields, and that are numerically stable.

The mesh generation algorithm is iterative. It commences with an initial distribution of vertices in the domain and rearranges the vertices to create higher-geometric quality elements. The bars of the triangular mesh act as *springs* that obey a constitutive law (e.g., Hooke's Law) otherwise referred to as a *force function*. During each meshing iteration, the discrepancy between the length of the edges in the mesh connectivity and their target length produce movement in the triangles' vertices. After a sufficient number of iterations, an equilibrium-like state is approached and the movement of the vertices ceases. The equilibrium-like state corresponds to a mesh that contains mostly isotropic equilateral triangles, which is important for accurate numerical simulation. As with most mesh generators, a sequence of mesh improvement strategies are applied after mesh generation terminates to ensure the mesh will be numerically robust for simulation.

The primary computational bottleneck in the *DistMesh* algorithm is the time spent in the Delauany triangulation algorithm, which occurs each iteration of the mesh generation step as the points are incrementally move. The other steps involving the calculation of the target sizing field and signed distance function can be executed in constant time provided the user defines these function as gridded interpolants. This  fact motivated us to pursue a distributed memory parallelism approach similar to [hpc_del]_ to parallize the Delaunay algorithm, which reduces the time spent performing each meshing iteration and makes the approach feasible for large-scale 3D mesh generation.


.. References
.. ..........

.. [hpc_del] Peterka, Tom, Dmitriy Morozov, and Carolyn Phillips. "High-performance computation of distributed-memory parallel 3D Voronoi and Delaunay tessellation." SC'14: Proceedings of the International Conference for High Performance Computing, Networking, Storage and Analysis. IEEE, 2014.

.. [distmesh] P.-O. Persson, G. Strang, A Simple Mesh Generator in MATLAB.
              SIAM Review, Volume 46 (2), pp. 329-345, June 2004 (PDF)

.. [firedrake] Florian Rathgeber, David A. Ham, Lawrence Mitchell, Michael Lange, Fabio Luporini, Andrew T. T. Mcrae, Gheorghe-Teodor Bercea, Graham R. Markall, and Paul H. J. Kelly. Firedrake: automating the finite element method by composing abstractions. ACM Trans. Math. Softw., 43(3):24:1–24:27, 2016. URL: http://arxiv.org/abs/1501.01809, arXiv:1501.01809, doi:10.1145/2998441.

.. [cgal] The CGAL Project. CGAL User and Reference Manual. CGAL Editorial Board, 5.0.2 edition, 2020

.. [slivers] Tournois, Jane, Rahul Srinivasan, and Pierre Alliez. "Perturbing slivers in 3D Delaunay meshes." Proceedings of the 18th international meshing roundtable. Springer, Berlin, Heidelberg, 2009. 157-173.
