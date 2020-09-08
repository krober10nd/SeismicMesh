from numpy import allclose, sum
import pytest

from SeismicMesh import generate_mesh, geometry


@pytest.mark.serial
def test_2d_SDF():
    """Unit circle"""
    hmin = 0.2
    bbox = (-1.0, 1.0, -1.0, 1.0)

    def circle(p):
        return geometry.dcircle(p, 0.0, 0.0, 1.0)

    def EF(p):
        d = circle(p)
        return hmin - d * 0.15

    points, cells = generate_mesh(
        bbox=bbox,
        signed_distance_function=circle,
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