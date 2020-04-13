import os
import SeismicMesh


def test_3dmesher():

    fname = os.path.join(os.path.dirname(__file__), "test3D.bin")
    wl = 10
    hmin = 50
    freq = 2
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
    mshgen = SeismicMesh.MeshGenerator(ef, method="qhull")
    points, cells = mshgen.build(nscreen=1, max_iter=20, seed=0)
    # 16459 vertices and 102868 cells
    assert len(points) == 16459
    assert len(cells) == 102868
    import meshio

    # Write to disk (see meshio for more details)
    meshio.write_points_cells(
        "foo3D_V3.vtk", points, [("tetra", cells)],
    )


if __name__ == "__main__":
    test_3dmesher()
