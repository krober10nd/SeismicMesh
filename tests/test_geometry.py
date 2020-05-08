import numpy as np
from SeismicMesh import geometry as geo
from SeismicMesh.plots import simpplot as smp
from SeismicMesh.plots import points as pts


def test_geometry3D():

    # tetraheadrals in a unit cube with l, w, h of 2 units
    points = np.array(
        [
            [-1, -1, -1],
            [1, -1, -1],
            [1, 1, -1],
            [-1, 1, -1],
            [-1, -1, 1],
            [1, -1, 1],
            [1, 1, 1],
            [-1, 1, 1],
            [0, 0, 0],
        ],
        dtype=np.float,
    )

    cells = np.array(
        [
            [3, 4, 0, 8],
            [3, 7, 8, 2],
            [6, 8, 7, 2],
            [6, 5, 8, 2],
            [6, 5, 7, 8],
            [8, 7, 4, 5],
            [3, 7, 4, 8],
            [3, 8, 1, 2],
            [8, 5, 1, 2],
            [8, 4, 1, 5],
            [8, 4, 0, 1],
            [3, 8, 0, 1],
        ],
        dtype=np.int,
    )

    facets = geo.get_facets(cells)
    assert len(facets) == (12 * 4)

    boundary_facets = geo.get_boundary_facets(cells)
    boundary_facets = geo.remove_duplicates(boundary_facets)
    assert len(boundary_facets) == 12

    boundary_vertices = geo.get_boundary_vertices(cells, dim=3)
    assert len(boundary_vertices) == 8

    vol = geo.simpvol(points, cells)
    assert np.sum(vol) == 8.0


if __name__ == "__main__":
    test_geometry3D()
