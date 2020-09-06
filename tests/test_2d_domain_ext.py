import os

import numpy as np
import pytest

import SeismicMesh


@pytest.mark.serial
def test_2d_domain_ext():

    fname = os.path.join(os.path.dirname(__file__), "testing.segy")
    vp = SeismicMesh.read_segy(fname)
    wl = 5
    freq = 10
    hmin = 100
    hmax = 1000
    grade = 0.005
    domain_ext = 1e3
    ef = SeismicMesh.MeshSizeFunction(
        bbox=(-10e3, 0, 0, 10e3),
        grade=grade,
        domain_ext=domain_ext,
        wl=wl,
        freq=freq,
        velocity_grid=vp,
        hmin=hmin,
        hmax=hmax,
    )
    ef = ef.build()
    mshgen = SeismicMesh.MeshGenerator(ef)
    points, cells = mshgen.build(max_iter=100, seed=0, nscreen=1)
    # 12284 vertices and 24144 cells
    print(points.shape, cells.shape)
    assert np.abs(len(points) - 12284) < 20
    assert np.abs(len(cells) - 24130) < 20

    ef = SeismicMesh.MeshSizeFunction(
        bbox=(-10e3, 0, 0, 10e3),
        grade=grade,
        domain_ext=500,
        pad_style="constant",
        wl=wl,
        freq=freq,
        velocity_grid=vp,
        hmin=hmin,
        hmax=hmax,
    )
    ef = ef.build()
    mshgen = SeismicMesh.MeshGenerator(ef)
    points, cells = mshgen.build(max_iter=100, seed=0, nscreen=1)
    # (10992, 2) (21569, 3)
    assert np.abs(len(points) - 10992) < 20
    assert np.abs(len(cells) - 21569) < 20

    ef = SeismicMesh.MeshSizeFunction(
        bbox=(-10e3, 0, 0, 10e3),
        grade=grade,
        domain_ext=2000,
        pad_style="linear_ramp",
        wl=wl,
        freq=freq,
        velocity_grid=vp,
        hmin=hmin,
        hmax=hmax,
    )
    ef = ef.build()
    mshgen = SeismicMesh.MeshGenerator(ef)
    points, cells = mshgen.build(max_iter=100, seed=0, nscreen=1)
    # 15260 30029

    print(len(points), len(cells))
    assert np.abs(len(points) - 15260) < 20
    assert np.abs(len(cells) - 30029) < 20


if __name__ == "__main__":
    test_2d_domain_ext()
