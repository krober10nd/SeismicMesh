import pytest
import os

import numpy as np

import SeismicMesh


@pytest.mark.serial
def test_2dmesher():

    fname = os.path.join(os.path.dirname(__file__), "testing.segy")
    wl = 5
    hmin = 100
    grade = 0.005
    ef = SeismicMesh.MeshSizeFunction(
        bbox=(-10e3, 0, 0, 10e3), grade=grade, wl=wl, model=fname, hmin=hmin
    )
    ef = ef.build()
    mshgen = SeismicMesh.MeshGenerator(ef, method="qhull")
    points, facets = mshgen.build(max_iter=100, seed=0)
    # should have: 7690 vertices and 15045 cells
    assert np.abs(len(points) - 7690) < 20
    assert np.abs(len(facets) - 15045) < 20

    mshgen = SeismicMesh.MeshGenerator(ef, method="cgal")
    points, facets = mshgen.build(max_iter=100, seed=0)
    # should have: 7690 vertices and 15045 cells
    assert np.abs(len(points) - 7690) < 20
    assert np.abs(len(facets) - 15045) < 20


if __name__ == "__main__":
    test_2dmesher()
