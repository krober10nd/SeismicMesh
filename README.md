<p align="center">
  <a href="https://github.com/krober10nd/SeismicMesh"><img alt="SeismicMesh" src="https://user-images.githubusercontent.com/18619644/92964244-28f31d00-f44a-11ea-9aa0-3d8ed0a1b60e.jpg" width="40%"></a>
  <p align="center">Create high-quality 2D/3D meshes from seismic velocity models.</p>
</p>


[![CircleCI](https://img.shields.io/circleci/project/github/krober10nd/SeismicMesh/master.svg?style=flat-square)](https://circleci.com/gh/krober10nd/SeismicMesh/tree/master)
[![CodeCov](https://codecov.io/gh/krober10nd/SeismicMesh/branch/master/graph/badge.svg)](https://codecov.io/gh/krober10nd/SeismicMesh)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg?style=flat-square)](https://github.com/psf/black)
[![PyPI pyversions](https://img.shields.io/pypi/pyversions/SeismicMesh.svg?style=flat-square)](https://pypi.org/pypi/SeismicMesh/)
[![PyPi downloads](https://img.shields.io/pypi/dm/SeismicMesh.svg?style=flat-square)](https://pypistats.org/packages/seismicmesh)
[![ReadTheDocs](https://readthedocs.org/projects/seismicmesh/badge/?version=master)](https://seismicmesh.readthedocs.io/en/master/?badge=master)
[![Zenodo](https://zenodo.org/badge/216707188.svg)](https://zenodo.org/badge/latestdoi/216707188)
[![PyPi]( https://img.shields.io/pypi/v/SeismicMesh.svg?style=flat-square)](https://pypi.org/project/SeismicMesh)
[![GPL](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)


[SeismicMesh](https://github.com/krober10nd/SeismicMesh): Mesh generation for Seismology in Python
==================================================================================================

*SeismicMesh* is a tool to generate 2D/3D triangular meshes that are used for acoustic and elastic wave propagators discretized
with the finite element method. With short scripts, variable resolution meshes are built that have size transitions which reflect variations in P-wave or S-wave velocities. Seismic P-wave/S-wave data are typically provided on [regular Cartesian grids for global and regional domains](https://ds.iris.edu/spud/earthmodel).

Generating variable resolution unstructured meshes of complex and large-scale 2D/3D geophysical domains with popular mesh generation tools such as [gmsh](https://gmsh.info) or [cgal](https://doc.cgal.org/latest/Mesh_3/index.html) requires deciding how to size elements in the domain, which can become a complex task. Typically users must either create their own mesh sizing function that reflect the geophysical data in the domain or rely on analytical mesh sizing functions. However, for seismological problems with mesh size variations mostly in the interior of the domain, mesh sizing heuristics like [boundary layer/attractor adaptation](https://gmsh.info/doc/texinfo/gmsh.html) or [characteristic size calculations](https://gmsh.info/doc/texinfo/gmsh.html) may not be directly relevant. Instead, in a typical seismologial domain, variations in mesh size generally reflect internal material properties such as P-wave or S-wave velocity, which are used to cost-effectively model waves by ensuring there are sufficient number of grid-points per wavelength.

Thus in this work we provide an approach to build large scale 2D/3D mesh sizing functions with the [mesh sizing function module](https://seismicmesh.readthedocs.io/en/master/api.html#seimsicmesh-meshsizefunction). This tool maps variations in seismic velocities from a seismic velocity model to triangular mesh sizes. Importantly, the sizing module can ensure that mesh size transitions vary smoothly (e.g., are graded) and an estimate of the [Courant number](https://en.wikipedia.org/wiki/Courant%E2%80%93Friedrichs%E2%80%93Lewy_condition) can be bounded above by a constant--amongst other capabilities--which are important considerations for accurate and successful simulations with finite elements. Our sizing functions are defined on the same regular Cartesian grids as the original seismic velocity model bypassing the need for the user to create their own sizing function on a triangular mesh. Their structured nature also enables efficient performance.

*SeismicMesh* supports both 2D and 3D triangular meshing in either serial or using distributed memory parallelism relying on the [Message Passing Interface](https://computing.llnl.gov/tutorials/mpi/). In 3D mesh generation, a mesh improvement method ([sliver removal](https://hal.inria.fr/inria-00430202/document)) is used to ensure a minimum quality bound (e.g. [minimum dihedral angle](https://en.wikipedia.org/wiki/Dihedral_angle#:~:text=A%20dihedral%20angle%20is%20the,line%20as%20a%20common%20edge)) can be enforced and will lead to numerically stable simulations.

Our mesh generation approach provided in this package can be operated standalone (e.g., without the sizing function module). It is based off modifications to the [DistMesh](http://persson.berkeley.edu/distmesh/) algorithm. Thus in its most basic operation, *SeismicMesh* can mesh *any* domain that can be defined by a [signed distance function](https://en.wikipedia.org/wiki/Signed_distance_function#:~:text=In%20mathematics%20and%20its%20applications,whether%20x%20is%20in%20%CE%A9.) with mesh sizes that follow variations described by a user-defined [mesh sizing function](http://persson.berkeley.edu/pub/persson05qualmesh.pdf)

*SeismicMesh* is distributed under the GPL3 license and more details can be found in our [short paper](https://github.com/krober10nd/SeismicMesh/blob/master/paper/paper.md).

Installation
============

For installation, SeismicMesh needs [CGAL](https://www.cgal.org/) and
[pybind11](https://github.com/pybind/pybind11):

    sudo apt install libcgal-dev python3-pybind11

After that, SeismicMesh can be installed from the Python Package Index
([pypi](https://pypi.org/project/SeismicMesh/)), so with:

    pip install -U SeismicMesh

For more detailed information about installation and requirements see:

[Install](https://seismicmesh.readthedocs.io/en/master/install.html) -
How to install SeismicMesh.

Examples
========

The user can quickly build quality 2D/3D meshes from seismic velocity
models in serial/parallel.

**WARNING: To run the code snippet below you must download the 2D BP2004
seismic velocity model and then you must uncompress it (e.g., gunzip).
This file can be downloaded from**
[here](http://s3.amazonaws.com/open.source.geoscience/open_data/bpvelanal2004/vel_z6.25m_x12.5m_exact.segy.gz)

![Above shows the mesh in ParaView that results from running the code below](https://user-images.githubusercontent.com/18619644/91606181-004a2e00-e948-11ea-83a4-e5ce05c7d82f.png)

```python
from mpi4py import MPI
import meshio

from SeismicMesh import get_sizing_function_from_segy, generate_mesh, Rectangle

comm = MPI.COMM_WORLD

"""
Build a mesh of the BP2004 benchmark velocity model in serial or parallel
Takes roughly 1 minute with 2 processors and less than 1 GB of RAM.
"""

# Name of SEG-Y file containg velocity model.
fname = "vel_z6.25m_x12.5m_exact.segy"

# Bounding box describing domain extents (corner coordinates)
bbox = (-12000.0, 0.0, 0.0, 67000.0)

# Desired minimum mesh size in domain
hmin = 75.0

rectangle = Rectangle(bbox)

# Construct mesh sizing object from velocity model
ef = get_sizing_function_from_segy(
    fname,
    bbox,
    hmin=hmin,
    wl=10,
    freq=2,
    dt=0.001,
    grade=0.15,
    domain_pad=1e3,
    pad_style="edge",
)

points, cells = generate_mesh(domain=rectangle, edge_length=ef, h0=hmin)

if comm.rank == 0:
    # Write the mesh in a vtk format for visualization in ParaView
    # NOTE: SeismicMesh outputs assumes the domain is (z,x) so for visualization
    # in ParaView, we swap the axes so it appears as in the (x,z) plane.
    meshio.write_points_cells(
        "BP2004.vtk",
        points[:, [1, 0]] / 1000,
        [("triangle", cells)],
        file_format="vtk",
    )
```

**WARNING: To run the code snippet below you must download (and uncompress) the 3D EAGE
seismic velocity model from (WARNING: File is \~500 MB)**
[here](https://s3.amazonaws.com/open.source.geoscience/open_data/seg_eage_models_cd/Salt_Model_3D.tar.gz)

**WARNING: Computationaly demanding! Running this example takes around 3 minutes in serial and requires
around 2 GB of RAM due to the 3D nature of the problem and the domain
size.**

![Above shows the mesh in ParaView that results from running the code below.](https://user-images.githubusercontent.com/18619644/91606008-c5e09100-e947-11ea-97e2-58e4b2f23d2b.jpg)

<!--exdown-skip-->
```python
from mpi4py import MPI
import zipfile
import meshio

from SeismicMesh import (
    get_sizing_function_from_segy,
    generate_mesh,
    sliver_removal,
    Cube,
)

comm = MPI.COMM_WORLD

# Bounding box describing domain extents (corner coordinates)
bbox = (-4200.0, 0.0, 0.0, 13520.0, 0.0, 13520.0)

# Desired minimum mesh size in domain.
hmin = 150.0

# This file is in a big Endian binary format, so we must tell the program the shape of the velocity model.
path = "Salt_Model_3D/3-D_Salt_Model/VEL_GRIDS/"
if comm.rank == 0:
    # Extract binary file Saltf@@ from SALTF.ZIP
    zipfile.ZipFile(path + "SALTF.ZIP", "r").extract("Saltf@@", path=path)

fname = path + "Saltf@@"

# Dimensions of model (number of grid points in z, x, and y)
nx, ny, nz = 676, 676, 210

cube = Cube(bbox)

# A graded sizing function is created from the velocity model along with a signed distance function by passing
# the velocity grid that we created above.
# More details can be found here: https://seismicmesh.readthedocs.io/en/master/api.html

ef = get_sizing_function_from_segy(
    fname,
    bbox,
    hmin=hmin,
    dt=0.001,
    freq=2,
    wl=5,
    grade=0.15,
    hmax=5e3,
    domain_pad=250,
    pad_style="linear_ramp",
    nz=nz,
    nx=nx,
    ny=ny,
    byte_order="big",
    axes_order=(2, 0, 1),  # order for EAGE (x, y, z) to default order (z,x,y)
    axes_order_sort="F", # binary is packed in a FORTRAN-style
)

points, cells = generate_mesh(domain=cube, h0=hmin, edge_length=ef, max_iter=75)

# For 3D mesh generation, we provide an implementation to bound the minimum dihedral angle::
points, cells = sliver_removal(
    points=points, bbox=bbox, h0=hmin, domain=cube, edge_length=ef
)

# Meshes can be written quickly to disk using meshio and visualized with ParaView::
if comm.rank == 0:

    # NOTE: SeismicMesh outputs assumes the domain is (z,x,y) so for visualization
    # in ParaView, we swap the axes so it appears as in the (x,y,z) plane.
    meshio.write_points_cells(
        "EAGE_Salt.vtk",
        points[:, [1, 2, 0]] / 1000.0,
        [("tetra", cells)],
    )
```

**The user can still specify their own signed distance functions and sizing functions to `generate_mesh` (in serial or parallel) just like the original DistMesh algorithm. Try the codes below!**

![Above shows the mesh in ParaView that results from running the code below.](https://user-images.githubusercontent.com/18619644/93465337-05542a80-f8c1-11ea-8774-a059e215088f.png)

```python
# Mesh a unit cylinder
from mpi4py import MPI
from numpy import array, maximum, sqrt, zeros_like
import meshio

from SeismicMesh import generate_mesh, sliver_removal

comm = MPI.COMM_WORLD


hmin = 0.10
bbox = (-1.0, 1.0, -1.0, 1.0, -1.0, 1.0)


def cylinder(p):
    r, z = sqrt(p[:, 0] ** 2 + p[:, 1] ** 2), p[:, 2]
    d1, d2, d3 = r - 1.0, z - 1.0, -z - 1.0
    d4, d5 = sqrt(d1 ** 2 + d2 ** 2), sqrt(d1 ** 2 + d3 ** 2)
    d = maximum.reduce([d1, d2, d3])
    ix = (d1 > 0) * (d2 > 0)
    d[ix] = d4[ix]
    ix = (d1 > 0) * (d3 > 0)
    d[ix] = d5[ix]
    return d


def fh(p):
    # Note: for parallel execution this logic is required
    # since the decomposition of the sizing function passes a tuple to fh
    if type(p) == tuple:
        h = zeros_like(p[0]) + hmin
    else:
        h = array([hmin] * len(p))
    return h


points, cells = generate_mesh(
    bbox=bbox,
    domain=cylinder,
    h0=hmin,
    edge_length=fh,
    max_iter=100,
)

points, cells = sliver_removal(
    points=points,
    domain=cylinder,
    edge_length=fh,
    h0=hmin,
    min_dh_angle_bound=5.0,
    bbox=bbox,
)


if comm.rank == 0:
    meshio.write_points_cells(
        "Cylinder.vtk",
        points,
        [("tetra", cells)],
        file_format="vtk",
    )
```

![Disk](https://user-images.githubusercontent.com/18619644/95608173-1f51da80-0a33-11eb-90be-170beda85b5a.png)

```python
import meshio
import SeismicMesh

disk = SeismicMesh.Disk([0.0, 0.0], 1.0)
points, cells = SeismicMesh.generate_mesh(domain=disk, edge_length=0.1)
meshio.write_points_cells(
    "disk.vtk",
    points,
    [("triangle", cells)],
    file_format="vtk",
)
```

![Rect](https://user-images.githubusercontent.com/18619644/95607603-5d023380-0a32-11eb-9c6f-41fac9e00aa7.png)

```python
import meshio
import SeismicMesh

bbox = (0.0, 1.0, 0.0, 1.0, 0.0, 1.0)
cube = SeismicMesh.Cube(bbox)
points, cells = SeismicMesh.generate_mesh(domain=cube, edge_length=0.05)
points, cells = SeismicMesh.sliver_removal(points=points, domain=cube, edge_length=0.05)
meshio.write_points_cells(
    "cube.vtk",
    points,
    [("tetra", cells)],
    file_format="vtk",
)
```

![cube](https://user-images.githubusercontent.com/18619644/95621214-b3c63800-0a47-11eb-9600-a80ef2410334.png)


```python
import meshio
import SeismicMesh

bbox = (0.0, 1.0, 0.0, 1.0, 0.0, 1.0)
cube = SeismicMesh.Cube(bbox)
points, cells = SeismicMesh.generate_mesh(domain=cube, edge_length=0.05)
meshio.write_points_cells(
    "cube.vtk",
    points,
    [("tetra", cells)],
    file_format="vtk",
)
```

![torus](https://user-images.githubusercontent.com/18619644/96498978-25aa3880-1223-11eb-9738-8a4e86c44dbc.png)


```python
# mesh a torus
import numpy as np
import meshio
import SeismicMesh


bbox = (-1.0, 1.0, -1.0, 1.0, -10.0, 1.0)
hmin = 0.05

def length(x):
    return np.sum(np.abs(x) ** 2, axis=-1) ** (1.0 / 2)

def Torus(p, t=(0.5, 0.2)):
    xz = np.column_stack((p[:, 0], p[:, 2]))
    q = np.column_stack((length(xz) - t[0], p[:, 1]))
    return length(q) - t[1]

points, cells = SeismicMesh.generate_mesh(
    domain=Torus, edge_length=hmin, bbox=bbox, verbose=2, max_iter=100
)
points, cells = SeismicMesh.sliver_removal(
    points=points, domain=Torus, edge_length=hmin, bbox=bbox, verbose=2
)
meshio.write_points_cells(
    "torus.vtk",
    points,
    [("tetra", cells)],
)
```

![prism](https://user-images.githubusercontent.com/18619644/96511116-f69cc280-1234-11eb-984e-b0001b15c7b5.png)


```python
# mesh a prism
import numpy as np
import meshio

import SeismicMesh


bbox = (-1.0, 1.0, -1.0, 1.0, -1.0, 1.0)

hmin = 0.05


def length(x):
    return np.sum(np.abs(x) ** 2, axis=-1) ** (1.0 / 2)


def sdTriPrism(p, h=[0.5, 0.5]):
    q = np.abs(p)
    return np.maximum(
        q[:, 2] - h[1],
        np.maximum(q[:, 0] * 0.866025 + p[:, 1] * 0.5, -p[:, 1]) - h[0] * 0.5,
    )


points, cells = SeismicMesh.generate_mesh(
    bbox=bbox,
    domain=sdTriPrism,
    edge_length=hmin,
    verbose=2,
    max_iter=100,
)
points, cells = SeismicMesh.sliver_removal(
    bbox=bbox, points=points, domain=sdTriPrism, edge_length=hmin
)
meshio.write_points_cells(
    "prism.vtk",
    points,
    [("tetra", cells)],
    file_format="vtk",
)
```

![Union](https://user-images.githubusercontent.com/18619644/96755280-e3aaff00-13a8-11eb-9f88-95a6684e928b.png)


```python
# Compute the union of several SDFs to create more complex geometries
import meshio
import SeismicMesh

h = 0.10
rect0 = SeismicMesh.Rectangle((0.0, 1.0, 0.0, 0.5))
rect1 = SeismicMesh.Rectangle((0.0, 0.5, 0.0, 1.0))
disk0 = SeismicMesh.Disk([0.5, 0.5], 0.5)
union = SeismicMesh.Union([rect0, rect1, disk0])
points, cells = SeismicMesh.generate_mesh(domain=union, edge_length=h)
meshio.write_points_cells(
    "Lshape_wDisk.vtk",
    points,
    [("triangle", cells)],
    file_format="vtk",
)
```
![Leaf](https://user-images.githubusercontent.com/18619644/96755336-f6bdcf00-13a8-11eb-99ec-bd7e7d9cad1d.png)

```python
# Compute the intersection of several SDFs to create more complex geometries
import meshio
import SeismicMesh

h = 0.05
rect0 = SeismicMesh.Rectangle((0.0, 1.0, 0.0, 1.0))
disk0 = SeismicMesh.Disk([0.25, 0.25], 0.5)
disk1 = SeismicMesh.Disk([0.75, 0.75], 0.5)
intersection = SeismicMesh.Intersection([rect0, disk0, disk1])
points, cells = SeismicMesh.generate_mesh(domain=intersection, edge_length=h)
meshio.write_points_cells(
    "Leaf.vtk",
    points,
    [("triangle", cells)],
    file_format="vtk",
)
```

![Hole](https://user-images.githubusercontent.com/18619644/96766828-0ab9fe80-13b2-11eb-8bca-6306934008d4.png)

```python
# Compute the difference of two SDFs to create more complex geometries.
import meshio
import SeismicMesh

h = 0.05
rect0 = SeismicMesh.Rectangle((0.0, 1.0, 0.0, 1.0))
disk0 = SeismicMesh.Disk([0.5, 0.5], 0.1)
disk1 = SeismicMesh.Disk([0.75, 0.75], 0.20)
difference = SeismicMesh.Difference([rect0, disk0, disk1])
points, cells = SeismicMesh.generate_mesh(domain=difference, edge_length=h)
meshio.write_points_cells(
    "Hole.vtk",
    points,
    [("triangle", cells)],
    file_format="vtk",
)
```

![Cube_wHoles](https://user-images.githubusercontent.com/18619644/96785337-0a772e80-13c5-11eb-88fb-311b5bfdfed4.png)

```python
import meshio
import SeismicMesh

h = 0.10
cube0 = SeismicMesh.Cube((0.0, 1.0, 0.0, 1.0, 0.0, 1.0))
ball1 = SeismicMesh.Ball([0.5, 0.0, 0.5], 0.30)
ball2 = SeismicMesh.Ball([0.5, 0.5, 0.0], 0.30)
ball3 = SeismicMesh.Ball([0.0, 0.5, 0.5], 0.30)
ball4 = SeismicMesh.Ball([0.5, 0.5, 0.5], 0.45)
difference = SeismicMesh.Difference([cube0, ball1, ball2, ball3, ball4])
points, cells = SeismicMesh.generate_mesh(domain=difference, edge_length=h, verbose=1)
points, cells = SeismicMesh.sliver_removal(
    points=points, domain=difference, edge_length=h, verbose=1
)
meshio.write_points_cells(
    "Ball_wHoles.vtk",
    points,
    [("tetra", cells)],
    file_format="vtk",
)
```

How does performance and cell quality compare to `gmsh` and `cgal` mesh generators?
===================================================================================

Here we use SeismicMesh 3.0.4, [pygalmesh](https://github.com/nschloe/pygalmesh) 0.8.2, and [pygmsh](https://github.com/nschloe/pygmsh) 7.0.0 (more details in the benchmarks folder).

Some key findings:

* Mesh generation in 2D and 3D using analytical sizing functions is quickest when using `gmsh` followed by `cgal` and then `SeismicMesh`.
* However, using mesh sizing functions defined on gridded interpolants significantly slow down both `gmsh` and `cgal`. In these cases, `SeismicMesh` and `gmsh` perform similarly both outperforming `cgal`'s 3D mesh generator in terms of mesh generation time.
* `SeismicMesh` produces often comparable or higher mean cell qualities in 2D/3D than either `gmsh` or `cgal` and this may have implications on mesh improvement strategies as higher minimum quality may be realizable with some common mesh improvement strategies (e.g., NetGen)
* All methods produce 3D triangulations that have a minimum dihedral angle > 10 degrees enabling stable numerical simulation.
* Head over to the `benchmarks` folder for more detailed information on these experiments.

![Summary of the benchmarks](https://user-images.githubusercontent.com/18619644/95696741-b923ae00-0c12-11eb-9d96-e52e5e9de7ae.jpg)


* **In the figure for the panels that show cell quality, solid lines indicate the mean and dashed lines indicate the minimum cell quality in the mesh.**

* Note: it's important to point out here that a significant speed-up can be achieved for moderate to large problems using the [parallel capabilities](https://seismicmesh.readthedocs.io/en/master/tutorial.html#basics) provided in `SeismicMesh`.


Changelog
=========

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project (tries to) adhere to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Unreleased

## [3.0.6] - 2020-10-21
### Fixed
- Silence messages about pfix when verbose=0
### Added
- Added more examples on README
- New unions/intersections/differences with several SDF primivitives
- Automatic corner constraints in serial

## [3.0.5] - 2020-10-18

### Fixed
- Preserving fixed points in serial.
- Units in km-s detection warning bug.
- Docstring fixes to `generate_mesh`
- Improved mesh quality in 3D

### Added
- Automatic corner point extraction for cubes and rectangles.
- More support for reading binary files packed in a binary format.
- Check to make sure bbox is composed of all floats.

## [3.0.4] - 2020-10-12

### Added
- Improve conformity of level-set in final mesh through additional set of Newton boundary projection iterations.


More information
================

All other information is available at:
<https://seismicmesh.readthedocs.io>

[Getting
started](https://seismicmesh.readthedocs.io/en/master/overview.html) -
Learn the basics about the program and the application domain.

[Tutorials](https://seismicmesh.readthedocs.io/en/master/tutorial.html) -
Tutorials that will guide you through the main features.
