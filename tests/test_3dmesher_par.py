import os

import numpy as np
from mpi4py import MPI

import SeismicMesh

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()


def test_3dpar_mesher():

    fname = os.path.join(os.path.dirname(__file__), "test3D.bin")
    wl = 10
    hmin = 50
    freq = 3
    grade = 0.35
    nz, nx, ny = 20, 10, 10
    ef = SeismicMesh.MeshSizeFunction(
        bbox=(-2e3, 0, 0, 1e3, 0, 1e3),
        nx=nx,
        ny=ny,
        nz=nz,
        endianness="big",
        grade=grade,
        freq=freq,
        wl=wl,
        model=fname,
        hmin=hmin,
    )
    # Build mesh size function (in parallel)
    ef = ef.build(comm=comm)
    # All processors get same information
    ef = comm.bcast(ef, 0)
    # Build lambda functions
    ef = ef.construct_lambdas()

    mshgen = SeismicMesh.MeshGenerator(
        ef, method="qhull"
    )  # parallel currently only works in qhull

    # Build the mesh with all axis combinations
    points, cells = mshgen.build(max_iter=30, nscreen=1, seed=0, COMM=comm, axis=0)

    if rank == 0:
        vol = SeismicMesh.geometry.simpvol(points / 1000, cells)
        assert np.abs(2 - np.sum(vol)) < 0.10  # km2

        assert np.abs(3965 - len(points)) < 10
        assert np.abs(21017 - len(cells)) < 10

        # intersections = SeismicMesh.geometry.doAnyOverlap(points, cells, dim=3)


if __name__ == "__main__":
    test_3dpar_mesher()
