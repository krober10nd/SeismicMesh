import os

import numpy as np
import pytest

from SeismicMesh import geometry as geo


@pytest.mark.serial
def test_geometry2():
    # two equilateral triangle laid next to each other
    points = np.array([[0, 0], [6, 0], [3, 3], [9, 3]], dtype=float)
    cells = np.array([[3, 2, 1], [2, 0, 1]], dtype=int)
    ce = geo.calc_re_ratios(points, cells, dim=2)
    assert np.allclose(ce, [0.70710678, 0.70710678])

    # get edges (shared are dup'ed)
    edges = geo.get_edges(cells)
    assert len(edges) == 6

    # get boundary edges (e.g., edge [1,2] should not be present)
    wedges = geo.get_winded_boundary_edges(cells)
    assert np.allclose(wedges, [[0, 1], [1, 3], [2, 3], [0, 2]])

    intersections = geo.do_any_overlap(points, cells, dim=2)
    assert len(intersections) == 0

    # test Laplacian2
    pdata = os.path.join(os.path.dirname(__file__), "points.txt")
    cdata = os.path.join(os.path.dirname(__file__), "cells.txt")
    p0 = np.loadtxt(pdata, delimiter=",")
    c0 = np.loadtxt(cdata, delimiter=",", dtype=int)

    geo.laplacian2(p0, c0)


if __name__ == "__main__":
    test_geometry2()
