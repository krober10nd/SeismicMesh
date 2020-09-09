from numpy import allclose, sum, array, maximum, sqrt
import pytest

from SeismicMesh import generate_mesh, sliver_removal, geometry


@pytest.mark.serial
def test_3dmesher_SDF():
    """Unit cylinder"""

    hmin = 0.10
    bbox = (-1.0, 1.0, -1.0, 1.0, -1.0, 1.0)

    def cylinder(p):
        r, z = sqrt(p[:, 0] ** 2 + p[:, 1] ** 2), p[:, 2]
        d1, d2, d3 = r - 1.0, z - 1.0, -z - 1.0
        d4, d5 = sqrt(d1 ** 2 + d2 ** 2), sqrt(d1 ** 2 + d3 ** 2)
        d = maximum.reduce([d1, d2, d3])
        ix = (d1 > 0) * (d2 > 0)
        d[ix] = d4[ix]
        ix = (d1 > 0) * (d3 > 0)
        d[ix] = d5[ix]
        return d

    def EF(p):
        return array(
            [
                hmin,
            ]
            * len(p)
        )

    points, cells = generate_mesh(
        bbox=bbox,
        domain=cylinder,
        h0=hmin,
        cell_size=EF,
        max_iter=100,
    )

    points, cells = sliver_removal(
        points=points, domain=cylinder, cell_size=EF, h0=hmin
    )

    assert allclose(sum(geometry.simp_vol(points, cells)), 6.28, atol=hmin)

    print(len(points), len(cells))
    allclose([len(points), len(cells)], [6825, 36206], atol=100)

    import meshio

    meshio.write_points_cells(
        "blah.vtk",
        points / 1000,
        [("tetra", cells)],
        file_format="vtk",
    )


if __name__ == "__main__":
    test_3dmesher_SDF()
