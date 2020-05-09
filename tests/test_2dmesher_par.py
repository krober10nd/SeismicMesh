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
    ef = ef.build()
    mshgen = SeismicMesh.MeshGenerator(ef, method="qhull")
    points, facets = mshgen.parallel_build(max_iter=100, seed=0, COMM=comm)
    if rank == 0:
        # should have points (7690, 2)
        # should have facets (15047, 3)
        assert np.abs(7690 - len(points)) < 5
        assert np.abs(15047 - len(facets)) < 5


if __name__ == "__main__":
    test_2dmesher_par()
