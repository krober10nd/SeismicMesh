import numpy as np
import pytest

import SeismicMesh


def quarter_annulus(h):
    rect = SeismicMesh.Rectangle((0.0, 1.0, 0.0, 1.0))
    disk0 = SeismicMesh.Disk([0.0, 0.0], 0.25)
    diff0 = SeismicMesh.Difference([rect, disk0])

    disk1 = SeismicMesh.Disk([0.0, 0.0], 1.0)
    quarter = SeismicMesh.Intersection([diff0, disk1])

    points, cells = SeismicMesh.generate_mesh(
        domain=quarter,
        edge_length=lambda x: h + 0.1 * np.abs(disk0.eval(x)),
        h0=h,
        verbose=0,
    )
    return points, cells


@pytest.mark.serial
def test_2d_min_qual():

    hmins = np.logspace(-1.0, -2.0, num=5)
    for h in hmins:
        points, cells = quarter_annulus(h)
        quals = SeismicMesh.geometry.simp_qual(points, cells)
        assert np.amin(quals) > 0.60
