Overview
========

This software aims to create end-to-end workflows to build quality two- and three-dimensional (2D and 3D) unstructured triangular and tetrahedral meshes for seismic domains suitable for numerical wave propagators. An open-source Python software package called *SeismicMesh* has been developed to this end. A focus is placed on parallel unstructured mesh generation and adaptation capabilities for numerical simulations using the finite element method. The capabilities of the software are exposed to the user via a simple application program interfaces (API) written in Python. The software can be scaled up on distributed memory clusters. The technology has been used to generate meshes in 2D and 3D acoustic and elastic wave propagators written in *Firedrake's* [firedrake]_ Uniformed Form Language, which are used in Full Waveform Inversion algorithms.


Software architecture
-------------------------------------------

The software is implemented in a mixed language environment (Python and C++). The Python language is used for the API while computationally expensive operations are performed in C++. The two languages are linked together with *pybind11*. The Computational Geometry Algorithms Library [_cgal] is used to perform geometric operations that use floating point arithmetic to avoid numerical precision issues. Besides this, we use several very common Python packages: *Numpy*, *Scipy*, *MeshIO*, *SegyIO*, and *MPI4py*. In particular, *MeshIO* enables to the program to output meshes in dozens of commonly used different formats.


Inputs
-------------------------------------------

The only required input to generate a mesh is a binary file (e.g., SEG-y) containing the P-wave velocities of the domain. The software makes the assumption the domain can be approximated by a rectangle/cube. Thus, the user specifies the domain geometry as a tuple of coordinates in meters: the bottom front corner and the top back corner. The program automatically generates the rectangle/cube domain geometry used during meshing. Interior structures are meshed using signed distance functions.


*DistMesh*
-------------------------------------------

[distmesh]_

.. References
.. ..........

.. [distmesh] P.-O. Persson, G. Strang, A Simple Mesh Generator in MATLAB.
              SIAM Review, Volume 46 (2), pp. 329-345, June 2004 (PDF)

.. [firedrake] P.-O. Persson, G. Strang, A Simple Mesh Generator in MATLAB.
              SIAM Review, Volume 46 (2), pp. 329-345, June 2004 (PDF)
