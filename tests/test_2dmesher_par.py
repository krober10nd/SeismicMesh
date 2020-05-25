import os

from mpi4py import MPI
import numpy as np

import SeismicMesh

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()


def test_2dmesher_par():

    fname = os.path.join(os.path.dirname(__file__), "testing.segy")
    wl = 5
    hmin = 100
    grade = 0.005
    ef = SeismicMesh.MeshSizeFunction(
        bbox=(-10e3, 0, 0, 10e3), grade=grade, wl=wl, model=fname, hmin=hmin
    )
    ef = ef.build(comm=comm)
    # Build lambda functions
    ef = ef.construct_lambdas(comm)

    mshgen = SeismicMesh.MeshGenerator(ef, method="qhull")
    points, cells = mshgen.build(max_iter=100, seed=0, COMM=comm)
    if rank == 0:
        # import meshio

        # meshio.write_points_cells(
        #    "test2d.vtk", points / 1000, [("triangle", cells)], file_format="vtk",
        # )
        area = SeismicMesh.geometry.simpvol(points / 1000, cells)
        assert np.abs(100 - np.sum(area)) < 0.10  # km2
        assert np.abs(8788 - len(points)) < 20
        assert np.abs(17203 - len(cells)) < 20


if __name__ == "__main__":
    test_2dmesher_par()
