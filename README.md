![image](https://circleci.com/gh/krober10nd/SeismicMesh/tree/par3d.svg?style=shield%0A%20%20%20%20%20:target:%20https://circleci.com/gh/krober10nd/SeismicMesh/tree/par3d)

![image](https://codecov.io/gh/krober10nd/SeismicMesh/branch/par3d/graph/badge.svg%0A%20:target:%20https://codecov.io/gh/krober10nd/SeismicMesh)

![image](https://img.shields.io/badge/code%20style-black-000000.svg%0A%20%20%20%20%20:target:%20https://github.com/ambv/black)

![image](http://www.repostatus.org/badges/latest/active.svg%0A%20:target:%20http://www.repostatus.org/#active)

![image](https://readthedocs.org/projects/seismicmesh/badge/?version=par3d%0A%20%20%20%20%20:target:%20https://seismicmesh.readthedocs.io/en/par3d/?badge=par3d)

![image](https://img.shields.io/badge/License-GPLv3-blue.svg%0A%20:target:%20https://www.gnu.org/licenses/gpl-3.0)

[![image](https://zenodo.org/badge/216707188.svg)](https://zenodo.org/badge/latestdoi/216707188)

[![image](https://img.shields.io/pypi/pyversions/SeismicMesh.svg?style=flat-square)](https://pypi.org/pypi/SeismicMesh)

[![image](https://img.shields.io/pypi/v/SeismicMesh.svg?style=flat-square)](https://pypi.org/project/SeismicMesh)

[![image](https://img.shields.io/pypi/dm/SeismicMesh.svg?style=flat-square)](https://pypistats.org/packages/seismicmesh)

[SeismicMesh](https://github.com/krober10nd/SeismicMesh): Mesh generation for Seismology in Python
==================================================================================================

2D/3D triangular meshing for a slab of Earth based on modifications to
the [DistMesh](http://persson.berkeley.edu/distmesh/) algorithm.
SeismicMesh is distributed under GPL3.

Installation
============

For installation, SeismicMesh needs [CGAL](https://www.cgal.org/) and
[pybind11](https://github.com/pybind/pybind11):

    sudo apt install libcgal-dev python-pybind11

After that, SeismicMesh can be installed from the Python Package Index
([pypi](https://pypi.org/project/SeismicMesh/)), so with:

    pip install -U SeismicMesh

For more detailed information about installation and requirements see:

[Install](https://seismicmesh.readthedocs.io/en/par3d/install.html) -
How to install SeismicMesh.

Example
=======

The user can quickly build quality 2D/3D meshes from seismic velocity
models in serial/parallel.

**WARNING: To run the code snippet below you must download the 2D BP2004
seismic velocity model and then you must uncompress it (e.g., gunzip).
This file can be downloaded from**
[here](http://s3.amazonaws.com/open.source.geoscience/open_data/bpvelanal2004/vel_z6.25m_x12.5m_exact.segy.gz)

![Above shows the mesh from running the code below. Note, the seismic
velocity data has been interpolated onto the vertices of the mesh using
a different program with linear
interpolation.](https://user-images.githubusercontent.com/18619644/91577721-82722c80-e91f-11ea-82f2-519687722e7b.jpg)

``` {.sourceCode .python}
from mpi4py import MPI
import meshio

import SeismicMesh

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

"""
Build a mesh of the BP2004 benchmark velocity model in serial or parallel
Takes roughly 1 minute with 2 processors and less than 1 GB of RAM.
"""

# Name of SEG-Y file containg velocity model.
fname = "vel_z6.25m_x12.5m_exact.segy"
# Read in it
vp = SeismicMesh.ReadSegy(fname)

# Bounding box describing domain extents (corner coordinates)
bbox = (-12000, 0.0, 0.0, 67000.0)

# Construct mesh sizing object from velocity model
ef = SeismicMesh.MeshSizeFunction(
    bbox=bbox,
    velocity_grid=vp,
    freq=2,
    wl=10,
    dt=0.001,
    hmin=75.0,
    grade=0.15,
    domain_ext=1e3,
    padstyle="linear_ramp",
)

# Build mesh size function
ef = ef.build()

# Construct a mesh generator object
mshgen = SeismicMesh.MeshGenerator(ef)

# Build the mesh
points, facets = mshgen.build(axis=1)

if rank == 0:
    # Write the mesh in a vtk format for visualization in ParaView
    # NOTE: SeismicMesh outputs assumes the domain is (z,x) so for visualization
    # in ParaView, we swap the axes so it appears as in the (x,z) plane.
    meshio.write_points_cells(
        "BP2004.vtk",
        points[:,[0,1]]/ 1000,
        [("triangle", facets)],
        file_format="vtk",
    )
```

**WARNING: To run the code snippet below you must download the 3D EAGE
seismic velocity model from (WARNING: File is \~500 MB)**
[here](https://s3.amazonaws.com/open.source.geoscience/open_data/seg_eage_models_cd/Salt_Model_3D.tar.gz)

**WARNING: Computationaly demanding! Running this example requires
around 8 GB of RAM due to the 3D nature of the problem and the domain
size.**

![Above shows the mesh from running the code below. Note, the seismic
velocity data has been interpolated onto the vertices of the mesh using
a different program with linear
interpolation.](https://user-images.githubusercontent.com/18619644/91485472-4be5d480-e881-11ea-9abf-75ae2fb6b2b1.jpg)

``` {.sourceCode .python}
import numpy as np
import zipfile

from mpi4py import MPI
import meshio

import SeismicMesh

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()


if rank == 0:
    # Dimensions of model (number of grid points in z, x, and y)
    nx, ny, nz = 676, 676, 210

    path = "Salt_Model_3D/3-D_Salt_Model/VEL_GRIDS/"
    # Extract Saltf@@ from SALTF.ZIP
    zipfile.ZipFile(path + "SALTF.ZIP", "r").extract("Saltf@@", path=path)

    # Load data into a numpy array
    with open(path + "Saltf@@", "r") as file:
        vp = np.fromfile(file, dtype=np.dtype("float32").newbyteorder(">"))
        vp = vp.reshape(nx, ny, nz, order="F")
        vp = np.flipud(vp.transpose((2, 0, 1)))  # z, x and then y
else:
    vp = np.zeros(shape=(1, 1, 1))
    vp[:] = 1500.0

# The domain is defined (in this case) as a cube and domain extents are provided in meters

# Bounding box describing domain extents (corner coordinates)
bbox = (-4200, 0, 0, 13520, 0, 13520)

# A graded sizing function is created from the velocity model along with a signed distance function by passing
# the velocity grid that we created above. More details for the :class:`MeshSizeFunction` can be found here
# https://seismicmesh.readthedocs.io/en/par3d/api.html#seimsicmesh-meshsizefunction

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

# The user then calls the mesh generator

# Construct a mesh generator object
mshgen = SeismicMesh.MeshGenerator(ef)

# Build the mesh
points, cells = mshgen.build(max_iter=75, axis=1)

# For 3D mesh generation, we provide an implementation to bound the minimum dihedral angle::

points, cells = mshgen.build(
    points=points, mesh_improvement=True, max_iter=50, min_dh_bound=5,
)

# Meshes can be written quickly to disk using meshio and visualized with ParaView::

if rank == 0:

    # NOTE: SeismicMesh outputs assumes the domain is (z,x,y) so for visualization
    # in ParaView, we swap the axes so it appears as in the (x,y,z) plane.
    meshio.write_points_cells(
        "EAGE_Salt.vtk", points[:,[1,2,0]]/ 1000.0, [("tetra", cells)],
    )
```

More information
================

All other information is available at:
<https://seismicmesh.readthedocs.io>

[Getting
started](https://seismicmesh.readthedocs.io/en/par3d/overview.html) -
Learn the basics about the program and the application domain.

[Tutorials](https://seismicmesh.readthedocs.io/en/par3d/tutorial.html) -
Tutorials that will guide you through the main features.
