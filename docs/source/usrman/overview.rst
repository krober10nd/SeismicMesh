Overview
========

This software aims to create end-to-end workflows (e.g., from seismic velocity model to mesh) to build quality two- and three-dimensional (2D and 3D) unstructured triangular and tetrahedral meshes for seismic domains. These meshes are suitable for acoustic and elastic numerical wave propagators. A focus is placed on parallel unstructured mesh generation and adaptation capabilities. The capabilities of the software are exposed to the user via a simple application program interfaces (API) written in Python and the software feature parallel algorithms that can be scaled up on distributed memory clusters.

*SeismicMesh* has been used to generate meshes in 2D and 3D for acoustic and elastic wave propagators written in *Firedrake* [firedrake]_. These type of numerical simulations are used in Full Waveform Inversion, Reverse Time Migration, and Time Travel Tomography applications.

Mesh?
-------------------------------------------
A mesh in our context is an unstructured triangulation composed of :math:`nt` triangles (otherwise referred to as entities) and :math:`np` vertices in either two or three dimensional space. These entities :math:`t` are obtained by tessellating a set of vertices that lie in two or three dimensional space using the well-known Delaunay triangulation. 

*High quality mesh?*
^^^^^^^^^^^^^^^^^^^^^^^
 We use the following formula to quantify how *equilateral* the mesh is [Bank1998]_
 
.. math::
  q_E = 4\sqrt{3}A_E\left(\sum_{i = 1}^{3}(\lambda_{E}^2)_i\right)^{-1}

where :math:`A_E` is the area of the element and :math:`\lambda_{E}_i` is the length of the :math:`i^{th}` edge of the entity. :math:`q_E = 1` corresponds to an equilateral triangular element and :math:`q_E = 0` indicates a completely degenerated element.

The consideration of what constitutes a *high quality* mesh rests on the statistical distribution of entity quality, :math:`q_E`. Generally a mesh with :math:`q_E > 0.95` and 
:math:`q_E - 3\sigma_{q_E} > 0.75` (where the over-line and :math:`\sigma` denote the mean and standard deviation respectively) is considered *high quality* and can be simulated without any changes to the mesh topology. Thus, when :math:`\overline{q_E} - 3\sigma_{q_E} > 0.75` is achieved, the mesh generation terminates and this generally occurs between 30-100 iterations in most seismic domains.

Software architecture
-------------------------------------------

The software is implemented in a mixed language environment (Python and C++). The Python language is used for the API while computationally expensive operations are performed in C++. The two languages are linked together with *pybind11* and installation is carried out using *cmake*. The Computational Geometry Algorithms Library [cgal]_ is used to perform geometric operations that use floating point arithmetic to avoid numerical precision issues. Besides this, several common Python packages: *Numpy*, *Scipy*, *MeshIO*, *SegyIO*, and *MPI4py* are used. In particular, *MeshIO* enables to the program to output meshes in dozens of commonly used different formats.

A typical workflow consists of two steps: first calling the *MeshSizeFunction* class followed by calling the *MeshGenerator* class.

*MeshSizeFunction*
^^^^^^^^^^^^^^^^^^^^^^^

This class creates an isotropic mesh size function and defines the domain through a signed distance function. The user passes a number of arguments to the class based on the intended application of the mesh. Specific arguments are detailed in the API section. 

*MeshGeneration*
^^^^^^^^^^^^^^^^^^^^^^^

This class takes as input a *MeshSizeFunction* object. However, it can also accept a user-defined signed distance function and mesh sizing function.  It constructs a mesh in either serial or parallel outputting points and simplices. Additionally, it has support a mesh improvement strategy [slivers]_ aimed at removing degenerate elements that are ubiqituous in 3D Delaunay mesh generation methods.


Inputs
-------------------------------------------

* The only required input file to generate a mesh in 2D is a binary file containing the velocity data on a structured grid. In 2D, the SEG-y format containing the seismic velocities of the domain is used. To store the seismic velocity model as a SEG-y file (if it isn't already in this format), the traces represent columns of the seismic velocity model. In 3D, the seismic velocity model is stored as a binary file but we do not use the SEG-y format. Instead in 3D, data is stored contiguously in memory in the format z,x,y following the little-endian format. For 3D, the user must specify how the data will be reshaped in memory by passing the number of rows, columns in the x-direction, and columns in the y-direction.


Signed distance function
^^^^^^^^^^^^^^^^^^^^^^^^^^

Mesh sizing function
^^^^^^^^^^^^^^^^^^^^^^^^^^


*DistMesh* algorithm
-------------------------------------------

For the generation of triangular meshes in 2D and 3D, we use a modified version of the *DistMesh* algorithm [distmesh]_. The algorithm is both simple (approximately 30 lines of code) and practically useful as it can produce high-geometric quality meshes in N-dimensional space. Further, by utilizing our approach to produce mesh size functions, the mesh generation algorithm is capable of producing high-geometric quality meshes that are faithful to our target sizing fields and that are numerically stable.

Briefly, the mesh generation algorithm is iterative. It commences with an initial distribution of vertices in the domain and iteratively relocates the vertices to create higher-geometric quality elements. The edges of the mesh act as *springs* that obey a constitutive law (e.g., Hooke's Law) otherwise referred to as a *force function*. During each meshing iteration, the discrepancy between the length of the edges in the mesh connectivity and their target length produce movement in the triangles' vertices. After a sufficient number of iterations, an equilibrium-like state is approached and the movement of the vertices becomes relatively small. The equilibrium-like state corresponds to a mesh that contains mostly isotropic equilateral triangles, which is critical for numerical simulation. As with most mesh generators, a sequence of mesh improvement strategies are applied after mesh generation terminates to ensure the mesh will be robust for simulation.


Mesh improvement
-------------------------------------------

Mesh adaptation
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. warning:
    Functionality to adapt an existing mesh is a work in progress


3D *Sliver* removal
^^^^^^^^^^^^^^^^^^^^^^^^^^

It is well-known that 3D mesh generation becomes significantly more challenging than 2D. Firstly, the computational cost is increased related to the additional spatial dimension and, for Delaunay approaches as we use here b) the presence of degenerate elements called *slivers* appears. If any sliver exists in a 3D mesh, the FEM solution can become numerically unstable and thus unusable. Fortunately, this problem does not occur in 2D but obviously needs to be dealt with in 3D. To tackle this problem in 3D, a method similar to that of [slivers ] aimed at removing low-quality elements while preserving the triangulation sizing distribution and domain structure was implemented.

The *sliver* removal technique fits within the *DistMesh* framework. Like the mesh generation approach, the algorithm operates iteratively. However, in this approach, it perturbs only vertices associated with *slivers* so that the circumsphere's radius of the *sliver* tetrahedral increases rapidly (e.g., gradient ascent of the circumsphere radius) [slivers]_. Futher, the method can operate on an existing mesh that already has a high-geometry mesh quality. The perturbation of a vertex of the *sliver* leads to a local combinational change in the nearby mesh connectivity to maintain Delaunayhood and almost always destroys the *sliver* in lieu of elements with larger dihedral angles.

Note here, we define *sliver* elements by their dihedral angle (i.e., angle between two surfaces) of which a tetrahedral has $6$. Generally, if a 3D mesh has a minimum dihedral angle less than 1 degree, it will not be numerically stable to simulate with.


Parallelism
-------------------------------------------

All algoirthms support distributed memory parallelism. When constructing models at scale, the primary computational bottleneck in the *DistMesh* algorithm becomes the time spent in the Delauany triangulation algorithm, which occurs each iteration of the mesh generation step. The other steps involving the formation and calculation of the target sizing field and signed distance function can be executed in constant time. Using *MPI4py* we implemented a simplified version of the [hpc_del]_ to parallize the Delaunay triangulation algorithm, and we later show this scales well and reduces the time spent performing each meshing iteration making the approach feasible for large-scale 3D mesh generation.


.. References
.. ..........

.. [hpc_del] Peterka, Tom, Dmitriy Morozov, and Carolyn Phillips. "High-performance computation of distributed-memory parallel 3D Voronoi and Delaunay tessellation." SC'14: Proceedings of the International Conference for High Performance Computing, Networking, Storage and Analysis. IEEE, 2014.

.. [distmesh] P.-O. Persson, G. Strang, A Simple Mesh Generator in MATLAB.
              SIAM Review, Volume 46 (2), pp. 329-345, June 2004 (PDF)

.. [firedrake] Florian Rathgeber, David A. Ham, Lawrence Mitchell, Michael Lange, Fabio Luporini, Andrew T. T. Mcrae, Gheorghe-Teodor Bercea, Graham R. Markall, and Paul H. J. Kelly. Firedrake: automating the finite element method by composing abstractions. ACM Trans. Math. Softw., 43(3):24:1–24:27, 2016. URL: http://arxiv.org/abs/1501.01809, arXiv:1501.01809, doi:10.1145/2998441.

.. [cgal] The CGAL Project. CGAL User and Reference Manual. CGAL Editorial Board, 5.0.2 edition, 2020

.. [slivers] Tournois, Jane, Rahul Srinivasan, and Pierre Alliez. "Perturbing slivers in 3D Delaunay meshes." Proceedings of the 18th international meshing roundtable. Springer, Berlin, Heidelberg, 2009. 157-173.
