import os

import numpy as np
import pytest
from mpi4py import MPI

import SeismicMesh

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()


@pytest.mark.parallel
def test_2dpar_mesher_adapt():
    fname = os.path.join(os.path.dirname(__file__), "testing.segy")
    vp = SeismicMesh.ReadSegy(fname)
    wl = 5
    hmin = 100
    grade = 0.005
    ef = SeismicMesh.MeshSizeFunction(
        bbox=(-10e3, 0, 0, 10e3), grade=grade, wl=wl, velocity_grid=vp, hmin=hmin
    )
    ef = ef.build()

    mshgen = SeismicMesh.MeshGenerator(ef)
    points, cells = mshgen.build(max_iter=20)

    points = comm.bcast(points, 0)

    # pass the points and restart with a different axis
    points, cells = mshgen.build(
        points=points, max_iter=20, axis=1, perform_checks=True
    )

    if rank == 0:
        import meshio

        meshio.write_points_cells(
            "test2d.vtk", points / 1000, [("triangle", cells)], file_format="vtk"
        )
        area = SeismicMesh.geometry.simpvol(points / 1000, cells)
        assert np.abs(100 - np.sum(area)) < 0.50  # km2


if __name__ == "__main__":
    test_2dpar_mesher_adapt()
