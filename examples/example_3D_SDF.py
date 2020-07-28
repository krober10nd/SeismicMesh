import numpy as np
import meshio
from mpi4py import MPI

import SeismicMesh

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

hmin = 0.10
bbox = (-1, 1, -1, 1, -1, 1)


def fd(p):
    r, z = np.sqrt(p[:, 0] ** 2 + p[:, 1] ** 2), p[:, 2]
    d1, d2, d3 = r - 1.0, z - 1.0, -z - 1.0
    return np.maximum.reduce([d1, d2, d3])


def fh(p):
    return np.array([hmin] * len(p))


# Construct mesh generator
mshgen = SeismicMesh.MeshGenerator(hmin=hmin, bbox=bbox, fd=fd, fh=fh)

# Build the mesh
points, cells = mshgen.build(max_iter=100)

if rank == 0:
    meshio.write_points_cells(
        "cylinder_3d.vtk", points / 1000, [("tetra", cells)], file_format="vtk",
    )
