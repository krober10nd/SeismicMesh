import os

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
    points, facets = mshgen.build(max_iter=100, seed=0, nscreen=1)
    # 12284 vertices and 24144 cells
    assert len(points) == 12284
    assert len(facets) == 24144


if __name__ == "__main__":
    test_2d_domain_ext()
