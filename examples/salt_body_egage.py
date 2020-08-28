import zipfile

import meshio
import numpy as np

import SeismicMesh
from SeismicMesh.geometry import (
    SignedDistanceFunctionGenerator as SdfGen,
)  # import the tool used to generate the SDF

# This example highlights how build a 3D mesh from a signed distance function built from an iso-contour in serial.

# The velocity model was downloaded from here:
# https://s3.amazonaws.com/open.source.geoscience/open_data/seg_eage_models_cd/Salt_Model_3D.tar.gz

# Dimensions of model
nx, ny, nz = 676, 676, 210

path = "velocity_models/Salt_Model_3D/3-D_Salt_Model/VEL_GRIDS/"
# Extract Saltf@@ from SALTF.ZIP
zipfile.ZipFile(path + "SALTF.ZIP", "r").extract("Saltf@@", path=path)

# Load data
with open(path + "Saltf@@", "r") as file:
    vp = np.fromfile(file, dtype=np.dtype("float32").newbyteorder(">"))
    vp = vp.reshape(nx, ny, nz, order="F")
    vp = np.flipud(vp.transpose((2, 0, 1)))  # z, x and then y


# Bounding box describing domain extents (corner coordinates)
bbox = (-4200, 0, 0, 13520, 0, 13520)

# Minimum mesh resolution
hmin = 100.0

ef = SeismicMesh.MeshSizeFunction(
    bbox=bbox,
    velocity_grid=vp,
    dt=0.001,
    freq=2,
    wl=5,
    grade=0.25,
    hmin=hmin,
    hmax=5e3,
)

# Build the sizing function
ef.build()

# Bulid a signed distance function from the geodataset
# Only capturing features with a minimum velocity greater than the speed of sound in water
# and less than the maximum velocity in the model
SDF = SdfGen(
    bbox=bbox,
    field=vp,
    min_threshold=3500,
    gridspacing=(20.0, 20.0, 20.0),
    narrow=1000,
).SDF

# Construct mesh generator
mshgen = SeismicMesh.MeshGenerator(hmin=hmin, bbox=bbox, fd=SDF, fh=ef.fh)

# Build the mesh
points, cells = mshgen.build(max_iter=50)

# Mesh improvement
points, cells = mshgen.build(points=points, mesh_improvement=True)

# Write the mesh to disk
# NOTE: SeismicMesh outputs assumes the domain is (z,x,y) so for visualization
# in ParaView, we swap the axes so it appears as in the (x,y,z) plane.
meshio.write_points_cells(
    "eage_salt_body.vtk",
    points[:, [1, 2, 0]] / 1000.0,
    [("tetra", cells)],
    file_format="vtk",
)
