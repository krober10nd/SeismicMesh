.. image:: https://circleci.com/gh/krober10nd/SeismicMesh/tree/parallel.svg?style=shield
        :target: https://circleci.com/gh/krober10nd/SeismicMesh/tree/parallel 

.. image:: https://codecov.io/gh/krober10nd/SeismicMesh/branch/parallel/graph/badge.svg
  	:target: https://codecov.io/gh/krober10nd/SeismicMesh
    
.. image:: https://img.shields.io/badge/code%20style-black-000000.svg
        :target: https://github.com/ambv/black


.. image:: http://www.repostatus.org/badges/latest/active.svg
	:alt: Project Status: Active – The project has reached a stable, usable state and is being actively developed.
	:target: http://www.repostatus.org/#active


SeismicMesh_: Mesh generation for Seismology in Python
==============================================
2D/3D triangular meshing for a slab of Earth based on modifications to the DistMesh_ algorithm. SeismicMesh is distributed under the GNU-GPL_ license.

.. _SeismicMesh: https://github.com/krober10nd/SeismicMesh
.. _DistMesh: http://persson.berkeley.edu/distmesh/
.. _`GNU-GPL`: http://www.gnu.org/copyleft/gpl.html

Installation :
==============================================

This software requires, a C++ compiler (GNU or Intel) with support for std++14, CMake >=3.0, CGAL >= 5.0 and Boost > 1.4.8. Once these packages are installed, you can run: 

1.  ``git submodule update --init --recursive``

2. ``pip install .``

If you do not have administrative rights on your system, use the flag ``--user`` 

1. ``git submodule update --init --recursive`` 

2. ``pip install --user .``

Notice the file ``Requirements.txt`` which indicates all the dependencies and their respective version numbers. If installing on a cluster with a local installation of ``CGAL`` and ``Boost``, make a directory called ``build`` in the root of the project (i.e., `SeismicMesh/``) and from within the build directory type ``ccmake ..`` and then specify the absolute paths of Boost, CGAL, MPFR, and GMP. Also, make sure that your compiler is not set to the default on the Linux system (i.e., an old version of GNU), which is typically much older than what is required for this program. 

Testing:
==============================================
To run tests, install ``pytest``i.e., ``pip install pytest

1. ``cd tests/``
2. ``pytest .``

Gallery:
==============================================
.. image:: imgs/seismic_example.png

