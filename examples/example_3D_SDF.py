import meshio
import numpy as np
from mpi4py import MPI

import SeismicMesh

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

hmin = 0.10
bbox = (-1.0, 1.0, -1.0, 1.0, -1.0, 1.0)

# This examples builds a small and simple cylinder in either serial or parallel.


def fd(p):
    """sizing function of a cylinder in [-1., -1., -.1] x [1., 1., 1.]"""
    r, z = np.sqrt(p[:, 0] ** 2 + p[:, 1] ** 2), p[:, 2]
    d1, d2, d3 = r - 1.0, z - 1.0, -z - 1.0
    d4, d5 = np.sqrt(d1 ** 2 + d2 ** 2), np.sqrt(d1 ** 2 + d3 ** 2)
    d = np.maximum.reduce([d1, d2, d3])
    ix = (d1 > 0) * (d2 > 0)
    d[ix] = d4[ix]
    ix = (d1 > 0) * (d3 > 0)
    d[ix] = d5[ix]
    return d


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
points, cells = mshgen.build(max_iter=75)

# Remove degenerate slivers
points, cells = mshgen.build(points=points, mesh_improvement=True)

if rank == 0:
    meshio.write_points_cells(
        "cylinder_3d.vtk",
        points / 1000,
        [("tetra", cells)],
        file_format="vtk",
    )
