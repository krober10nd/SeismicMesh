import os

import numpy as np
import pytest

import SeismicMesh


@pytest.mark.serial
def test_2dmesher():

    fname = os.path.join(os.path.dirname(__file__), "testing.segy")
    vp = SeismicMesh.read_segy(fname)
    wl = 5
    hmin = 100
    grade = 0.005
    ef = SeismicMesh.MeshSizeFunction(
        bbox=(-10e3, 0, 0, 10e3),
        grade=grade,
        grad=50.0,
        wl=wl,
        velocity_grid=vp,
        hmin=hmin,
    )
    ef = ef.build()
    ef.plot()
    mshgen = SeismicMesh.MeshGenerator(ef)
    points, facets = mshgen.build(max_iter=100)
    ef.write_velocity_model("testing")
    ef.save_mesh_size_function_options("sizingOpts")
    ef.plot()
    # should have: 7690 vertices and 15045 cells
    assert np.abs(len(points) - 7690) < 20
    assert np.abs(len(facets) - 15045) < 20


if __name__ == "__main__":
    test_2dmesher()
