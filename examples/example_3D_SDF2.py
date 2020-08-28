import meshio
import numpy as np
from mpi4py import MPI

import SeismicMesh

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

hmin = 0.10
bbox = (-1.0, 1.0, -1.0, 1.0, -1.0, 1.0)


# signed distancem function
def fd(p):
    return np.sqrt((p ** 2).sum(1)) - 1.0


# sizing function
def fh(p):
    # note for parallel execution this logic is required
    # since the decomposition of the sizing function passes a tuple to fh
    if type(p) == tuple:
        h = np.zeros_like(p[0]) + hmin
    else:
        h = np.zeros_like(p) + hmin
    return h


# Construct mesh generator
mshgen = SeismicMesh.MeshGenerator(hmin=hmin, bbox=bbox, fd=fd, fh=fh)

# Build the mesh
points, cells = mshgen.build(max_iter=75, axis=1)

# Remove degenerate slivers
points, cells = mshgen.build(points=points, mesh_improvement=True)

if rank == 0:
    meshio.write_points_cells(
        "cylinder_3d.vtk",
        points / 1000,
        [("tetra", cells)],
        file_format="vtk",
    )
