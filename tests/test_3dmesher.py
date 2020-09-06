import os

import numpy as np
import pytest

import SeismicMesh


@pytest.mark.serial
def test_3dmesher():

    fname = os.path.join(os.path.dirname(__file__), "test3D.bin")
    nz, nx, ny = 20, 10, 10

    # Load data
    with open(fname, "r") as file:
        vp = np.fromfile(file, dtype=np.dtype("float32").newbyteorder("<"))
        vp = vp.reshape(nx, ny, nz, order="F")
        vp = np.flipud(vp.transpose((2, 0, 1)))  # z, x and then y

    wl = 10
    hmin = 50
    freq = 2
    grade = 0.005
    ef = SeismicMesh.MeshSizeFunction(
        bbox=(-2e3, 0, 0, 1e3, 0, 1e3),
        grade=grade,
        grad=hmin,
        freq=freq,
        wl=wl,
        velocity_grid=vp,
        hmin=hmin,
    )
    ef = ef.build()
    ef.WriteVelocityModel("foo3d")
    mshgen = SeismicMesh.MeshGenerator(ef)
    points, cells = mshgen.build(nscreen=1, max_iter=50, seed=0)
    print(len(points), len(cells))
    # 16459 90552

    assert len(points) == 16459
    assert np.abs(len(cells) - 90552) < 600

    points, cells = mshgen.build(points=points, mesh_improvement=True)

    import meshio

    meshio.write_points_cells("food3D.vtk", points, [("tetra", cells)])
    meshio.write_points_cells(
        "foo3d.msh",
        points / 1000.0,
        [("tetra", cells)],
        file_format="gmsh22",
        binary=False,
    )


if __name__ == "__main__":
    test_3dmesher()
