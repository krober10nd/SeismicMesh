import numpy as np
import math
import pytest

import SeismicMesh


def cylinder(h):
    cylinder = SeismicMesh.Cylinder(h=1.0, r=1.0)
    points, cells = SeismicMesh.generate_mesh(domain=cylinder, edge_length=h, verbose=1)
    points, cells = SeismicMesh.sliver_removal(
        points=points, domain=cylinder, edge_length=h, verbose=1
    )
    return points, cells


@pytest.mark.serial
def test_3d_sliver():

    hmins = np.logspace(-1.0, -1.5, num=5)
    min_dh_bound = 10 * math.pi / 180
    max_dh_bound = 180 * math.pi / 180
    for h in hmins:
        points, cells = cylinder(h)
        dh_angles = SeismicMesh.geometry.calc_dihedral_angles(points, cells)
        out_of_bounds = np.argwhere(
            (dh_angles[:, 0] < min_dh_bound) | (dh_angles[:, 0] > max_dh_bound)
        )
        ele_nums = np.floor(out_of_bounds / 6).astype("int")
        ele_nums, ix = np.unique(ele_nums, return_index=True)
        assert len(ele_nums) == 0
