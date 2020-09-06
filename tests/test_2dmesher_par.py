import os

import numpy as np
import pytest
from mpi4py import MPI

import SeismicMesh

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()


@pytest.mark.parallel
def test_2dmesher_par():

    fname = os.path.join(os.path.dirname(__file__), "testing.segy")
    vp = SeismicMesh.read_segy(fname)
    wl = 5
    hmin = 100
    grade = 0.005
    ef = SeismicMesh.MeshSizeFunction(
        bbox=(-10e3, 0, 0, 10e3), grade=grade, wl=wl, velocity_grid=vp, hmin=hmin
    )
    ef = ef.build()

    # test cgal
    mshgen = SeismicMesh.MeshGenerator(ef)
    points, cells = mshgen.build(max_iter=100, seed=0, perform_checks=True)
    if rank == 0:
        # import meshio

        # meshio.write_points_cells(
        #    "test2d.vtk", points / 1000, [("triangle", cells)], file_format="vtk",
        # )
        area = SeismicMesh.geometry.simpvol(points / 1000, cells)
        # 7658 vertices and 14965
        assert np.abs(100 - np.sum(area)) < 0.50  # km2
        assert np.abs(7658 - len(points)) < 100
        assert np.abs(14965 - len(cells)) < 100


if __name__ == "__main__":
    test_2dmesher_par()
