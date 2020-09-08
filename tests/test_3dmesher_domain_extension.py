import os

from numpy import allclose
import pytest

from SeismicMesh import (
    generate_mesh,
    sliver_removal,
    geometry,
    get_sizing_function_from_segy,
    write_velocity_model,
)


@pytest.mark.serial
@pytest.mark.parametrize(
    "style_answer",
    (
        ("linear_ramp", [1388, 6592]),
        ("edge", [1383, 6545]),
        ("constant", [1406, 6622]),
    ),
)
def test_3dmesher_domain_ext(style_answer):
    style, answer = style_answer
    fname = os.path.join(os.path.dirname(__file__), "test3D.bin")
    wl = 5
    freq = 2
    hmin = 150
    grade = 0.005
    bbox = (-2e3, 0.0, 0.0, 1e3, 0, 1e3)
    ef, bbox = get_sizing_function_from_segy(
        fname,
        bbox,
        grade=grade,
        grad=hmin,
        freq=freq,
        wl=wl,
        hmin=hmin,
        nz=20,
        nx=10,
        ny=10,
        byte_order="little",
        domain_ext=200,
        pad_style=style,
    )

    write_velocity_model(
        fname, nz=20, nx=10, ny=10, byte_order="little", ofname="testing"
    )

    def cube(p):
        return geometry.dblock(p, *bbox)

    points, cells = generate_mesh(
        bbox=bbox,
        signed_distance_function=cube,
        cell_size=ef,
        h0=hmin,
        perform_checks=False,
    )

    points, cells = sliver_removal(
        points=points, bbox=bbox, signed_distance_function=cube, h0=hmin
    )
    print(len(points), len(cells))
    allclose([len(points), len(cells)], answer, atol=100)

    import meshio

    meshio.write_points_cells("foo3D" + style + ".vtk", points, [("tetra", cells)])
