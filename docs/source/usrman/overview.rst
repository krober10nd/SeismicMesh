Overview
========


This software aims to create end-to-end workflows (e.g., from seismic velocity model to simulation ready mesh) to build quality two- and three-dimensional (2D and 3D) unstructured triangular and tetrahedral meshes for seismic domains. The main advantage of triangular meshes over cubes or rectangles is their ability to cost-effectively resolve variable material properties. These generated meshes are suitable for acoustic and elastic numerical wave propagators and a focus is placed on parallel unstructured mesh generation capabilities. A simple application program interfaces (API) written in Python enables the user to call parallel meshing algorithms that can be scaled up on distributed memory clusters.


*SeismicMesh* is currently being used to generate simplical meshes in 2D and 3D for acoustic and elastic wave propagators written using a Domain Specific Language called *Firedrake* [firedrake]_. These type of numerical simulations are used in Full Waveform Inversion, Reverse Time Migration, and Time Travel Tomography applications.


Mesh definition
-------------------------------------------
.. note ::
    Triangular and tetrahedral conformal meshes.

The domain :math:`\Omega` is partitioned into a finite set of cells :math:`\mathcal{T}_{h} = {T}` with disjoint interiors
such that

:math:`\cup_{T \in \mathcal{T}_{h}} T = \Omega`

Together, these cells form a mesh of the domain :math:`\Omega`. In our case, the cells are triangles and thus the mesh is a triangulation composed of :math:`nt` triangles and :math:`np` vertices in either two or three dimensional space. Note in 3D, the triangulation is composed of tetrahedral but we still refer to it as a triangulation. In 2D, a triangle has 3 vertices, 3 edges, and 1 face. In 3D, a tetrahedral has 4 vertices, 6 edges, and 4 facets. These cells :math:`t` are obtained by tessellating a set of vertices that lie in two or three dimensional space using the well-known and efficient Delaunay triangulation algorithms.


*High quality cells*
^^^^^^^^^^^^^^^^^^^^^^^
.. note ::
    A high quality cell has a cell quality of 1.

Any set of points can be triangulated, but the resulting triangulation will likely not be useable for numerical simulation. Thus, we strive to produce high-quality geometric meshes.

There are many definitions of mesh quality. Here I use the following formula to quantify how *well-shaped* the cells where cell quality is defined as d * circumcircle_radius / incircle_radius (where d is 2 for triangles and 3 for tetrahedra). The value is between 0 and 1, where 1 is a perfectly symmetrical simplex.


Software architecture
-------------------------------------------
.. note ::
    Python calls peformant C++ libraries like CGAL and Boost.

The software is implemented in a mixed language environment (Python and C++). The Python language is used for the API while computationally expensive operations are performed in C++. The two languages are linked together with *pybind11* and installation is carried out using *cmake*. The Computational Geometry Algorithms Library [cgal]_ is used to perform geometric operations that use floating point arithmetic to avoid numerical precision issues. Besides this, several common Python packages: *numpy*, *scipy*, *meshio*, *segyio*, and *mpi4py* are used.


Inputs
-------------------------------------------

Seismic velocity model
^^^^^^^^^^^^^^^^^^^^^^^^

A seismic velocity model is defined on an axis-aligned regular Cartesian grid in either 2D or 3D containing scalar values (typically the P-wave velocity speed at each grid point).

Currently, *SeismicMesh* can read seismic velocity models from SEG-y files or in binary format. The latter requires some more information; see the docstring.

Signed distance function
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Given a point :math:`x`, the signed distance function returns the :math:`d` distance to the boundary of the domain :math:`∂ \Omega`.

Let :math:`\Omega ⊂ D ∈ R^N` be the domain in :math:`N` dimensions. The boundary of the domain to-be-meshed is represented as the 0-level set of a continuous function:

:math:`φ(·) : D → R.`

such that:

:math:`\Omega := {x ∈ D, φ(x) < 0}`

where :math:`φ : D × R+ → R` is Lipschitz continuous and called the level set function. If we assume :math:`|∇φ(·)| = 0` on the set :math:`{x ∈ D, φ(x) = 0}`, then we have :math:`∂ \Omega = {x ∈ D, φ(x) = 0}` i.e., the boundary :math:`∂ \Omega` is the 0-level set of :math:`φ(·)`. The property that :math:`|∇φ(·)| = 0` is satisfied if :math:`φ(·)` is a signed distance function.

.. note :: 
    We provide several simple signed distance functions: such as a Rectangle, Cube, and Disk. See the geometry module. 

Mesh sizing function
^^^^^^^^^^^^^^^^^^^^^^^^^^

Given a point :math:`x`, the sizing function :math:`f(h)` returns the isotropic mesh size defined at :math:`x`. By mesh size, we specifically mean the triangular edge length nearby `x` assuming the triangles will be close to equilateral in the finished mesh.

The purpose of `get_sizing_function_from_segy` to build this function directly from the seismic velocity model provided.


*DistMesh* algorithm
-------------------------------------------

.. note ::
    This program uses a modified version of the *DistMesh* algorithm [distmesh]_ to generate simplical meshes.

For the generation of triangular meshes in 2D and 3D, we use a modified version of the *DistMesh* algorithm [distmesh]_. The algorithm is both simple and practically useful as it can produce high-geometric quality meshes in N-dimensional space. Further, by utilizing our approach to produce mesh size functions, the mesh generation algorithm is capable of generating high-quality meshes faithful to user-defined target sizing fields. A benefit of this is that mesh sizes can be built to respect numerically stability requirements a priori.

Briefly, the mesh generation algorithm is iterative and terminates after a pre-set number of iterations (e.g., 50-100). It commences with an initial distribution of vertices in the domain and iteratively relocates the vertices to create higher-geometric quality elements. The edges of the mesh act as *springs* that obey a constitutive law (e.g., Hooke's Law) otherwise referred to as a *force function*. During each meshing iteration, the discrepancy between the length of the edges in the mesh connectivity and their target length from the sizing function produce movement in the triangles' vertices.

The boundary of the domain is enforced by projecting any points that leave the domain back into it each meshing iteration. After a sufficient number of iterations, an equilibrium-like state is almost always approached and the movement of the vertices becomes relatively small. The equilibrium-like state of the mesh connectivity corresponds to a mesh that contains mostly isotropic equilateral triangles, which is critical for numerical simulation. However, as with nearly all mesh generators, a sequence of mesh improvement strategies are applied after mesh generation terminates to ensure the mesh will be robust for simulation.


Mesh adaptation
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. warning ::
    Functionality to adapt an existing mesh is a work in progress


3D *Sliver* removal
^^^^^^^^^^^^^^^^^^^^^^^^^^

3D Delaunay mesh generation algorithms form degenerate elements called *slivers*. If any *sliver* exists in a 3D mesh, the numerical solution can become unstable. Fortunately, this problem does not occur in 2D and, in 2D, a high quality mesh free of degenerate elements is easily achieved. To tackle this problem in 3D, a method similar to that of [slivers]_ was implemented. This algorithm aims at removing *slivers* while preserving the triangulation sizing distribution and domain boundary.


The *sliver* removal technique fits well within the *DistMesh* framework. For example, like the mesh generation approach, the algorithm operates iteratively. Each meshing iteration, it perturbs *only* vertices associated with *slivers* so that the circumspheres' radius of the *sliver* tetrahedral increases rapidly (i.e.., gradient ascent of the circumsphere radius) [slivers]_. The method operates on an existing mesh that ideally already has a high-mesh quality and is efficient since it uses CGAL's incremental Delaunay capabilities. The perturbation of a vertex of the *sliver* leads to a local combinational change in the nearby mesh connectivity to maintain Delaunayhood and almost always destroys the *sliver* in lieu of elements with larger dihedral angles.

.. note ::
    A *sliver* element is defined by their dihedral angle (i.e., angle between two surfaces) of which a tetrahedral has :math:`6`. Generally, if a 3D mesh has a minimum dihedral angle less than 1 degree, it will be numerically unstable. We've had success in simulating with meshes that have minimum dihedral angles of minimally around 5 degrees.


Parallelism and speed
-------------------------------------------

.. note ::
    This code uses distributed memory parallelism with the MPI4py package.

When constructing models at scale, the primary computational bottleneck in the *DistMesh* algorithm becomes the time spent in the Delauany triangulation algorithm, which occurs each iteration of the mesh generation step. The other steps involving the formation and calculation of the target sizing field and signed distance function are far less demanding. Using *mpi4py*, I implemented a simplified version of the [hpc_del]_ to parallelize the Delaunay triangulation algorithm. This approach scales well and reduces the time spent performing each meshing iteration thus making the approach feasible for large-scale 3D mesh generation applications. The domain is decomposed into axis-aligned *slices* than cut one axis of the domain. While this strategy doesn't fare well with load balancing, it simplifies the implementation and runtime communication cost associated with neighboring processor exchanges.

When possible, *SeismicMesh* uses low-level functionality from the CGAL package including the evaluation of geometric predicates, circumball calculations, polygonal intersection tests, and incremental triangulation capabilities.


References
-------------------------------------------

.. [hpc_del] Peterka, Tom, Dmitriy Morozov, and Carolyn Phillips. "High-performance computation of distributed-memory parallel 3D Voronoi and Delaunay tessellation." SC'14: Proceedings of the International Conference for High Performance Computing, Networking, Storage and Analysis. IEEE, 2014.

.. [distmesh] P.-O. Persson, G. Strang, A Simple Mesh Generator in MATLAB.
              SIAM Review, Volume 46 (2), pp. 329-345, June 2004 (PDF)

.. [firedrake] Florian Rathgeber, David A. Ham, Lawrence Mitchell, Michael Lange, Fabio Luporini, Andrew T. T. Mcrae, Gheorghe-Teodor Bercea, Graham R. Markall, and Paul H. J. Kelly. Firedrake: automating the finite element method by composing abstractions. ACM Trans. Math. Softw., 43(3):24:1–24:27, 2016. URL: http://arxiv.org/abs/1501.01809, arXiv:1501.01809, doi:10.1145/2998441.

.. [cgal] The CGAL Project. CGAL User and Reference Manual. CGAL Editorial Board, 5.0.2 edition, 2020

.. [slivers] Tournois, Jane, Rahul Srinivasan, and Pierre Alliez. "Perturbing slivers in 3D Delaunay meshes." Proceedings of the 18th international meshing roundtable. Springer, Berlin, Heidelberg, 2009. 157-173.
