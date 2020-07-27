import meshio

import SeismicMesh
from SeismicMesh.geometry import (
    SignedDistanceFunctionGenerator as SdfGen,
)  # import the tool used to generate the SDF

# This example highlights how build a mesh from a custom
# signed distance function. It constructs a mesh of the BP2004 model
# without a water layer by thresholding the velocity model.


# Name of SEG-Y file containg velocity model.
fname = "velocity_models/vel_z6.25m_x12.5m_exact.segy"
# Bounding box describing domain extents (corner coordinates)
bbox = (-12e3, 0, 0, 67e3)
# Minimum mesh resolution
hmin = 75.0

# Construct mesh sizing object from velocity model
ef = SeismicMesh.MeshSizeFunction(
    bbox=bbox, model=fname, freq=2, wl=10, dt=0.001, hmin=hmin, grade=0.15,
)

# Build the sizing function
ef.build()
ef.plot()

# Bulid a signed distance function from the geodataset
# Only capturing features with a minimum velocity greater than the speed of sound in water
# and less than the maximum velocity in the model
SDF = SdfGen(
    bbox=bbox,
    field=ef.vp,
    min_threshold=1487,
    max_threshold=6000,
    gridspacing=(6.25, 12.5),
).SDF

# Construct mesh generator
mshgen = SeismicMesh.MeshGenerator(hmin=hmin, bbox=bbox, fd=SDF, fh=ef.fh)

# Build the mesh
points, cells = mshgen.build(max_iter=100)

# Write the mesh to disk
meshio.write_points_cells(
    "no_water_layer_2d.vtk", points / 1000.0, [("triangle", cells)], file_format="vtk"
)
