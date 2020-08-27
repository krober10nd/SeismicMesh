.. image:: https://circleci.com/gh/krober10nd/SeismicMesh/tree/par3d.svg?style=shield
        :target: https://circleci.com/gh/krober10nd/SeismicMesh/tree/par3d 

.. image:: https://codecov.io/gh/krober10nd/SeismicMesh/branch/par3d/graph/badge.svg
  	:target: https://codecov.io/gh/krober10nd/SeismicMesh
    
.. image:: https://img.shields.io/badge/code%20style-black-000000.svg
        :target: https://github.com/ambv/black


.. image:: http://www.repostatus.org/badges/latest/active.svg
	:target: http://www.repostatus.org/#active

.. image:: https://readthedocs.org/projects/seismicmesh/badge/?version=par3d
        :target: https://seismicmesh.readthedocs.io/en/par3d/?badge=par3d
	
.. image:: https://img.shields.io/badge/License-GPLv3-blue.svg
	:target: https://www.gnu.org/licenses/gpl-3.0
	
.. image:: https://zenodo.org/badge/216707188.svg
   :target: https://zenodo.org/badge/latestdoi/216707188



SeismicMesh_: Mesh generation for Seismology in Python
=========================================================
2D/3D triangular meshing for a slab of Earth based on modifications to the DistMesh_ algorithm. SeismicMesh is distributed under GPL3.

.. _SeismicMesh: https://github.com/krober10nd/SeismicMesh
.. _DistMesh: http://persson.berkeley.edu/distmesh/
.. _`GNU-GPL`: http://www.gnu.org/copyleft/gpl.html


Installation
=====================

For installation, SeismicMesh needs CGAL (https://www.cgal.org/)::

    sudo apt install libcgal-dev

and pybind11 (https://github.com/pybind/pybind11):: 

    sudo apt-get install python-pybind11

After that, SeismicMesh can be installed from the Python Package
Index (https://pypi.org/project/SeismicMesh/), so with::

    pip install -U SeismicMesh

For more detailed information about installation and requirements see: 

`Install <https://seismicmesh.readthedocs.io/en/par3d/install.html>`_
- How to install SeismicMesh. 


Example 
===========

The user can quickly build quality 2D/3D meshes from seismic velocity models in serial/parallel like so. 
First we take care of our imports::

    import numpy as np
    import zipfile
    
    from mpi4py import MPI
    import meshio

    import SeismicMesh
    
    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()

Here we download a benchmark 3D velocity model and turn it into a Numpy array::

    # https://s3.amazonaws.com/open.source.geoscience/open_data/seg_eage_models_cd/Salt_Model_3D.tar.gz

    if rank == 0: 
        # Dimensions of model (number of grid points in z, x, and y)
        nx, ny, nz = 676, 676, 210

        path = "velocity_models/Salt_Model_3D/3-D_Salt_Model/VEL_GRIDS/"
        # Extract Saltf@@ from SALTF.ZIP
        zipfile.ZipFile(path + "SALTF.ZIP", "r").extract("Saltf@@", path=path)

        # Load data into a numpy array
        with open(path + "Saltf@@", "r") as file:
            vp = np.fromfile(file, dtype=np.dtype("float32").newbyteorder(">"))
            vp = vp.reshape(nx, ny, nz, order="F")
            vp = np.flipud(vp.transpose((2, 0, 1)))  # z, x and then y

The domain is defined (in this case) as a cube and domain extents are provided in meters::

    # Bounding box describing domain extents (corner coordinates)
    bbox = (-4200, 0, 0, 13520, 0, 13520)

A graded sizing function is created from the velocity model along with a signed distance function by passing
the velocity grid that we created above::

    ef = SeismicMesh.MeshSizeFunction(
        bbox=bbox,
        velocity_grid=vp,
        dt=0.001,
        freq=2,
        wl=5,
        grade=0.25,
        hmin=150,
        hmax=5e3,
        domain_ext=250,
        padstyle="linear_ramp",
    )

    ef = ef.build()

The user then calls the mesh generator::

    # Construct a mesh generator object
    mshgen = SeismicMesh.MeshGenerator(ef)

    # Build the mesh
    points, cells = mshgen.build(max_iter=75, axis=1)

For 3D mesh generation, we provide an implementation to bound the minimum dihedral angle::

    points, cells = mshgen.build(
        points=points, mesh_improvement=True, max_iter=50, min_dh_bound=5,
    )

Meshes can be written quickly to disk using meshio and visualized with Paraview::

    if rank == 0:
        meshio.write_points_cells(
            "EAGE_Salt.vtk", points / 1000.0, [("tetra", cells)],
        )

 
More information
==================

All other information is available at: https://seismicmesh.readthedocs.io

`Getting started <https://seismicmesh.readthedocs.io/en/par3d/overview.html>`_
- Learn the basics about the program and the application domain. 

`Tutorials <https://seismicmesh.readthedocs.io/en/par3d/tutorial.html>`_
- Tutorials that will guide you through the main features.


Gallery:
==============================================

.. image:: https://github.com/krober10nd/SeismicMesh/raw/par3d/imgs/seismic_example3.png
.. image:: https://github.com/krober10nd/SeismicMesh/raw/par3d/imgs/seismic_example.png

