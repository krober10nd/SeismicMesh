import pytest
import os

import numpy as np
from mpi4py import MPI

import SeismicMesh

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()


@pytest.mark.parallel
def test_3dpar_mesher():

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

    points, cells = mshgen.build(
        max_iter=50,
        nscreen=1,
        seed=0,
        axis=0,
    )

    points, cells = mshgen.build(
        points=points,
        max_iter=15,
        mesh_improvement=True,
    )

    if rank == 0:
        # import meshio

        # meshio.write_points_cells(
        #    "foo3D_V3.vtk", points, [("tetra", cells)],
        # )

        vol = SeismicMesh.geometry.simpvol(points / 1000, cells)
        assert np.abs(2 - np.sum(vol)) < 0.10  # km2
        print(len(points), len(cells))
        assert np.abs(9220 - len(points)) < 1000
        assert np.abs(49156 - len(cells)) < 1000


if __name__ == "__main__":
    test_3dpar_mesher()
