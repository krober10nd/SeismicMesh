from numpy import allclose, sum, sqrt, array
import pytest

from SeismicMesh import generate_mesh, geometry


@pytest.mark.serial
def test_2dmesher_SDF():
    """Unit circle"""
    hmin = 0.2
    bbox = (-1.0, 1.0, -1.0, 1.0)

    def circle(p):
        """Signed distance to circle centered at xc, yc with radius r."""
        return sqrt(((p - array([0, 0])) ** 2).sum(-1)) - 1

    def EF(p):
        d = circle(p)
        return hmin - d * 0.15

    points, cells = generate_mesh(
        bbox=bbox,
        domain=circle,
        h0=hmin,
        cell_size=EF,
        max_iter=100,
    )

    print(len(points), len(cells))
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
