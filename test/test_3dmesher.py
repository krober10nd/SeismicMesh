import os
import SeismicMesh


def test_3dmesher():

    fname = os.path.join(os.path.dirname(__file__), "test3D.bin")
    wl = 10
    hmin = 50
    freq = 5
    grade = 0.005
    nz, nx, ny = 20, 10, 10
    ef = SeismicMesh.MeshSizeFunction(
        bbox=(-2e3, 0, 0, 1e3, 0, 1e3),
        nx=nx,
        ny=ny,
        nz=nz,
        endianness="big",
        grade=grade,
        freq=freq,
        wl=wl,
        model=fname,
        hmin=hmin,
    )
    ef = ef.build()
    mshgen = SeismicMesh.MeshGenerator(ef, method="cgal")
    points, cells = mshgen.build(nscreen=1, max_iter=20, seed=0)
    # should have  17577 vertices and 111191 cells
    assert len(points) == 17577
    assert len(cells) == 111191


if __name__ == "__main__":
    test_3dmesher()
