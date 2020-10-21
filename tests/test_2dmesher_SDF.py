import pytest
from numpy import allclose, sum

from SeismicMesh import generate_mesh, geometry


@pytest.mark.serial
def test_2dmesher_SDF():
    """Unit disk"""
    hmin = 0.2
    bbox = (-1.0, 1.0, -1.0, 1.0)

    disk = geometry.Disk([0.0, 0.0], 1)

    def EF(p):
        d = disk.eval(p)
        return hmin - d * 0.15

    points, cells = generate_mesh(
        bbox=bbox,
        domain=disk,
        h0=hmin,
        edge_length=EF,
        max_iter=100,
    )

    # print(len(points), len(cells))
    assert allclose([len(points), len(cells)], [63, 93], atol=10)
    assert allclose(sum(geometry.simp_vol(points, cells)), 3.14, atol=hmin)

    import meshio

    meshio.write_points_cells(
        "blah.vtk",
        points[:, [1, 0]] / 1000,
        [("triangle", cells)],
        file_format="vtk",
    )


if __name__ == "__main__":
    test_2dmesher_SDF()
