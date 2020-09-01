import os

import numpy as np
import pytest
from mpi4py import MPI

import SeismicMesh

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()


@pytest.mark.parallel
def test_3dpar_mesher_adapt():
    fname = os.path.join(os.path.dirname(__file__), "test3D.bin")
    nz, nx, ny = 20, 10, 10

    # Load data
    with open(fname, "r") as file:
        vp = np.fromfile(file, dtype=np.dtype("float32").newbyteorder("<"))
        vp = vp.reshape(nx, ny, nz, order="F")
        vp = np.flipud(vp.transpose((2, 0, 1)))  # z, x and then y

    wl = 10
    hmin = 50
    freq = 4
    grade = 0.15
    nz, nx, ny = 20, 10, 10
    ef = SeismicMesh.MeshSizeFunction(
        bbox=(-2e3, 0, 0, 1e3, 0, 1e3),
        nx=nx,
        ny=ny,
        nz=nz,
        grade=grade,
        freq=freq,
        wl=wl,
        velocity_grid=vp,
        hmin=hmin,
    )
    ef = ef.build()

    mshgen = SeismicMesh.MeshGenerator(ef)

    # generate a mesh
    points, cells = mshgen.build(
        max_iter=10,
        axis=0,
    )
    points = comm.bcast(points, 0)

    # pass the points and restart with a different axis
    points, cells = mshgen.build(
        points=points,
        max_iter=10,
        axis=1,
    )

    points = comm.bcast(points, 0)

    # pass the points and restart with a different axis
    points, cells = mshgen.build(
        points=points,
        max_iter=10,
        axis=2,
    )

    if rank == 0:
        # import meshio

        # meshio.write_points_cells(
        #    "foo3D_V3.vtk",
        #    points,
        #    [("tetra", cells)],
        # )

        vol = SeismicMesh.geometry.simpvol(points / 1000, cells)
        assert np.abs(2 - np.sum(vol)) < 0.10  # km2


if __name__ == "__main__":
    test_3dpar_mesher_adapt()
