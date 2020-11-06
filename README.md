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


[SeismicMesh](https://github.com/krober10nd/SeismicMesh): Triangular Mesh generation in Python
==================================================================================================

*SeismicMesh* is a tool to generate unstructured and high geometric quality 2D/3D triangular meshes that *can be* used for acoustic and elastic wave propagators discretized
with the finite element method. *SeismicMesh* however can also be used as a general purpose simplicial mesh generation tool in 2D/3D. With short scripts, variable resolution meshes are built that have size transitions which reflect variations in P-wave or S-wave velocities from seismic velocity models. Seismic P-wave/S-wave data are typically provided on regular Cartesian grids for global and regional domains.

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


Contributing
============

All contributions are welcome!

To contribute to the software:

1. [Fork](https://docs.github.com/en/free-pro-team@latest/github/getting-started-with-github/fork-a-repo) the repository.
2. Clone the forked repository, add your contributions and push the changes to your fork.
3. Create a [Pull request](https://github.com/krober10nd/SeismicMesh/pulls)

Before creating the pull request, make sure that the tests pass by running
```
tox
```
Some things that will increase the chance that your pull request is accepted:
-  Write tests.
- Add Python docstrings that follow the [Sphinx](https://sphinx-rtd-tutorial.readthedocs.io/en/latest/docstrings.html).
- Write good commit and pull request messages.


[style]: https://sphinx-rtd-tutorial.readthedocs.io/en/latest/docstrings.html

Problems?
==========

If something isn't working as it should or you'd like to recommend a new addition/feature to the software, please let me know by starting an issue through the [issues](https://github.com/krober10nd/SeismicMesh/issues) tab. I'll try to get to it as soon as possible.


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
    axes_order_sort="F",  # binary is packed in a FORTRAN-style
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

**The user can still specify their own signed distance functions and sizing functions to `generate_mesh` (in serial or parallel) just like the original DistMesh algorithm but now with quality bounds in 3D. Try the codes below!**

<img alt="Cylinder" src="https://user-images.githubusercontent.com/18619644/97082301-0e7e9880-15df-11eb-9055-15394213d755.png" width="30%">

```python
# Mesh a cylinder
from mpi4py import MPI
import meshio

import SeismicMesh

comm = MPI.COMM_WORLD

hmin = 0.10

cylinder = SeismicMesh.Cylinder(h=1.0, r=0.5)

points, cells = SeismicMesh.generate_mesh(
    domain=cylinder,
    edge_length=hmin,
)

points, cells = SeismicMesh.sliver_removal(
    points=points,
    domain=cylinder,
    edge_length=hmin,
)

if comm.rank == 0:
    meshio.write_points_cells(
        "Cylinder.vtk",
        points,
        [("tetra", cells)],
        file_format="vtk",
    )
```

<img alt="Disk" src="https://user-images.githubusercontent.com/18619644/97063883-b9a83700-1578-11eb-9cd7-3ff0cbac20d9.png" width="30%">


```python
# mesh a disk
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

<img alt="Square" src="https://user-images.githubusercontent.com/18619644/97063852-7b127c80-1578-11eb-97d5-cfe07cc969ec.png" width="30%">

```python
# mesh a square/rectangle
import meshio
import SeismicMesh

bbox = (0.0, 1.0, 0.0, 1.0)
square = SeismicMesh.Rectangle(bbox)
points, cells = SeismicMesh.generate_mesh(domain=square, edge_length=0.05)
meshio.write_points_cells(
    "square.vtk",
    points,
    [("triangle", cells)],
    file_format="vtk",
)
```
<img alt="Cube" src="https://user-images.githubusercontent.com/18619644/97063751-e1e36600-1577-11eb-9387-613f3ae04bff.png" width="30%">

```python
# mesh a cuboid/cube
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
<img alt="Torus" src="https://user-images.githubusercontent.com/18619644/97063588-eeb38a00-1576-11eb-8cff-8e77ea4d2946.png" width="30%">


```python
# mesh a torus
import meshio
import SeismicMesh

hmin = 0.10

torus = SeismicMesh.Torus(r1=1.0, r2=0.5)
points, cells = SeismicMesh.generate_mesh(
    domain=torus,
    edge_length=hmin,
)
points, cells = SeismicMesh.sliver_removal(
    points=points, domain=torus, edge_length=hmin
)
meshio.write_points_cells(
    "torus.vtk",
    points,
    [("tetra", cells)],
)
```

<img alt="Torus" src="https://user-images.githubusercontent.com/18619644/97081705-8ac2ad00-15da-11eb-9466-a86216b8908c.png" width="30%">

```python
# mesh a prism
import meshio

import SeismicMesh

hmin = 0.05
prism = SeismicMesh.Prism(b=0.5, h=0.5)

points, cells = SeismicMesh.generate_mesh(
    domain=prism,
    edge_length=hmin,
)
points, cells = SeismicMesh.sliver_removal(
    points=points, domain=prism, edge_length=hmin
)
meshio.write_points_cells(
    "prism.vtk",
    points,
    [("tetra", cells)],
    file_format="vtk",
)
```
<img alt="Union" src="https://user-images.githubusercontent.com/18619644/97081772-045a9b00-15db-11eb-8356-7863cdf274a3.png" width="30%">

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
<img alt="Leaf" src="https://user-images.githubusercontent.com/18619644/97081808-41bf2880-15db-11eb-9333-2d1230621c01.png" width="30%">

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
<img alt="Hole" src="https://user-images.githubusercontent.com/18619644/97081829-69ae8c00-15db-11eb-815d-a8302f822337.png" width="30%">

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
<img alt="Cube wHoles" src="https://user-images.githubusercontent.com/18619644/97081862-ad08fa80-15db-11eb-94b2-801001137f1a.png" width="30%">

```python
# Compute the difference of several SDFs in 3D
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
    "Cube_wHoles.vtk",
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
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Unreleased 
- None

## [3.1.3] - 2020-11-06
### Fixed
- Cylinder radius and height are now correct.
- Torus, Prism, and Cylinder now have `dim` tag.

### Improved
- More control over the `grad` option in the mesh sizing function.

## [3.1.2] - 2020-11-04
### Improved
- Faster calculation of boundary vertices.
- More robust sliver removal in 3D.
### Fixed
- Corners are only constrained for constant resolution meshes

## [3.1.0] - 2020-10-28
### Added
- New geometric primitives--torus, wedge/prism, and cylinder.
- Updated images on README.
### Fixed
- Only constrain corners near 0-level set.
- Bug fix to 3D binary velocity reading.

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
