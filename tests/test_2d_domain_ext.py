import os
import numpy as np

import SeismicMesh


def test_2d_domain_ext():

    fname = os.path.join(os.path.dirname(__file__), "testing.segy")
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
        model=fname,
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
        padstyle="constant",
        wl=wl,
        freq=freq,
        model=fname,
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
        padstyle="linear_ramp",
        wl=wl,
        freq=freq,
        model=fname,
        hmin=hmin,
        hmax=hmax,
    )
    ef = ef.build()
    mshgen = SeismicMesh.MeshGenerator(ef)
    points, cells = mshgen.build(max_iter=100, seed=0, nscreen=1)
    # (14773, 2) (29076, 3)
    assert np.abs(len(points) - 14773) < 20
    assert np.abs(len(cells) - 29076) < 20


if __name__ == "__main__":
    test_2d_domain_ext()
