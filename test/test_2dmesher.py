import os

import numpy as np

import SeismicMesh


def test_2dmesher():

    fname = os.path.join(os.path.dirname(__file__), "testing.segy")
    wl = 5
    hmin = 100
    grade = 0.005
    ef = SeismicMesh.MeshSizeFunction(
        bbox=(-10e3, 0, 0, 10e3), grade=grade, wl=wl, segy=fname, hmin=hmin
    )
    ef = ef.build()
    mshgen = SeismicMesh.MeshGenerator(ef)
    points, facets = mshgen.build(max_iter=100)
    # should have: 7601 vertices and 14873 cells
    assert np.abs((len(points) - 7601)) < 50
    assert (len(facets) - 14873) < 50


if __name__ == "__main__":
    test_2dmesher()
