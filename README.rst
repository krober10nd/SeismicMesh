.. image:: https://circleci.com/gh/krober10nd/SeismicMesh/tree/master.svg?style=shield
        :target: https://circleci.com/gh/krober10nd/SeismicMesh/tree/master 

.. image:: https://codecov.io/gh/krober10nd/SeismicMesh/branch/master/graph/badge.svg
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

This software requires, CMake >=3.0, CGAL >= 5.0 and Boost. Once these packages are installed, you can run: 

1.  ``git submodule update --init --recursive``

2. ``pip install -r requirements.txt``

2. ``pip install .``

If you do not have administrative rights on your system, use the flag ``--user`` 

1. ``git submodule update --init --recursive`` 

2. ``pip install --user -r requirements.txt``

3. ``pip install --user .``

Testing:
==============================================
To run tests, install ``pytest``. 

1. ``pytest tests/``

Gallery:
==============================================
.. image:: imgs/seismic_example.png

