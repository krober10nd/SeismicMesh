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

This class takes as input a signed distance function and mesh sizing function OR a *MeshSizeFunction* object. It constructs a mesh outputting points and simplices. It supports both serial and distributed memory parallel execution. Additionally, it has support a mesh improvement strategy [slivers]_ aimed at removing degenerate elements that are ubiqituous in 3D Delaunay mesh generation methods.


Inputs
-------------------------------------------

* The only required input file to generate a mesh is a binary file (e.g., SEG-y) containing the P-wave velocities of the domain.


* If the user uses the *MeshSizeFunction* class, the software makes the assumption the domain can be approximated by a rectangle/cube. Thus, the user specifies the domain geometry as a tuple of coordinates in meters::
    :math:`bbox = (z_{min}, z_{max}, x_{min}, x_{max})`

* In 3D::
    :math:`bbox = (z_{min}, z_{max}, x_{min}, x_{max}, y_{min}, y_{max})`

.. note :: The program automatically generates the rectangle/cube domain geometry used during meshing if a *MeshSizeFunction* object is passed to the generator.


*DistMesh*
-------------------------------------------



[distmesh]_

.. References
.. ..........

.. [distmesh] P.-O. Persson, G. Strang, A Simple Mesh Generator in MATLAB.
              SIAM Review, Volume 46 (2), pp. 329-345, June 2004 (PDF)

.. [firedrake] Florian Rathgeber, David A. Ham, Lawrence Mitchell, Michael Lange, Fabio Luporini, Andrew T. T. Mcrae, Gheorghe-Teodor Bercea, Graham R. Markall, and Paul H. J. Kelly. Firedrake: automating the finite element method by composing abstractions. ACM Trans. Math. Softw., 43(3):24:1–24:27, 2016. URL: http://arxiv.org/abs/1501.01809, arXiv:1501.01809, doi:10.1145/2998441.

.. [cgal] The CGAL Project. CGAL User and Reference Manual. CGAL Editorial Board, 5.0.2 edition, 2020

.. [slivers] Tournois, Jane, Rahul Srinivasan, and Pierre Alliez. "Perturbing slivers in 3D Delaunay meshes." Proceedings of the 18th international meshing roundtable. Springer, Berlin, Heidelberg, 2009. 157-173.
