import math

import numpy as np
import pytest

import SeismicMesh

min_dh_bound = 10 * math.pi / 180
max_dh_bound = 180 * math.pi / 180


def calc_dh_angles(points, cells):
    dh_angles = SeismicMesh.geometry.calc_dihedral_angles(points, cells)
    out_of_bounds = np.argwhere(
        (dh_angles[:, 0] < min_dh_bound) | (dh_angles[:, 0] > max_dh_bound)
    )
    ele_nums = np.floor(out_of_bounds / 6).astype("int")
    ele_nums, ix = np.unique(ele_nums, return_index=True)
    return ele_nums


def box_with_refinement(h):
    cube = SeismicMesh.geometry.Cube((-1.0, 1.0, -1.0, 1.0, -1.0, 1.0))

    def edge_length(x):
        return h + 0.1 * np.sqrt(x[:, 0] ** 2 + x[:, 1] ** 2 + x[:, 2] ** 2)

    points, cells = SeismicMesh.generate_mesh(
        domain=cube, h0=h, edge_length=edge_length, verbose=1
    )
    points, cells = SeismicMesh.sliver_removal(
        domain=cube,
        points=points,
        h0=h,
        edge_length=edge_length,
        verbose=1,
    )

    return points, cells


def cylinder(h):
    cylinder = SeismicMesh.Cylinder(h=1.0, r=1.0)
    points, cells = SeismicMesh.generate_mesh(domain=cylinder, edge_length=h, verbose=0)
    points, cells = SeismicMesh.sliver_removal(
        points=points, domain=cylinder, edge_length=h, verbose=0
    )
    return points, cells


def l_shape_3d(h):
    tol = h / 10

    cube0 = SeismicMesh.Cube((-1.0, 1.0, -1.0, 1.0, -1.0, tol))
    cube1 = SeismicMesh.Cube((-1.0, 1.0, 0.0, 1.0, -tol, 1.0))
    cube2 = SeismicMesh.Cube((-tol, 1.0, -1.0, 1.0, -tol, 1.0))
    domain = SeismicMesh.Union([cube0, cube1, cube2])

    points, cells = SeismicMesh.generate_mesh(domain=domain, edge_length=h, verbose=0)
    points, cells = SeismicMesh.sliver_removal(
        domain=domain, points=points, edge_length=h, verbose=0
    )
    return points, cells


@pytest.mark.serial
def test_3d_sliver_cylinder():

    hmins = np.logspace(-1.0, -1.3, num=5)
    for h in hmins:
        points, cells = cylinder(h)
        ele_nums = calc_dh_angles(points, cells)
        assert len(ele_nums) == 0


@pytest.mark.serial
def test_3d_sliver_box3d():
    hmins = np.logspace(-1.0, -1.2, num=5)
    for h in hmins:
        points, cells = l_shape_3d(h)
        ele_nums = calc_dh_angles(points, cells)
        assert len(ele_nums) == 0


@pytest.mark.serial
def test_3d_sliver_box3d_refined():
    hmins = np.logspace(-1.0, -1.3, num=5)
    for h in hmins:
        points, cells = box_with_refinement(h)
        ele_nums = calc_dh_angles(points, cells)
        assert len(ele_nums) == 0
