import meshio

import SeismicMesh
from SeismicMesh.geometry import (
    SignedDistanceFunctionGenerator as SdfGen,
)  # import the tool used to generate the SDF

# This example highlights how build a 3D mesh from a signed distance function built from an iso-contour in serial.

# Name of SEG-Y file containg velocity model.
fname = "velocity_models/EGAGE_Salt.bin"
# Bounding box describing domain extents (corner coordinates)
bbox = (-4200, 0, 0, 13520, 0, 13520)

# Minimum mesh resolution
hmin = 100.0

ef = SeismicMesh.MeshSizeFunction(
    bbox=bbox,
    model=fname,
    nx=676,  # size of velocity model in x-direction
    ny=676,  # size of velocity model in y-direction
    nz=210,  # size of velocity model in z-direction
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
    field=ef.vp,
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
meshio.write_points_cells(
    "eage_salt_body.vtk", points / 1000.0, [("tetra", cells)], file_format="vtk"
)
