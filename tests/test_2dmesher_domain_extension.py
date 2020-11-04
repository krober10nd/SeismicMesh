import os

import pytest
from numpy import allclose

from SeismicMesh import (
    Rectangle,
    generate_mesh,
    get_sizing_function_from_segy,
    plot_sizing_function,
    write_velocity_model,
)


@pytest.mark.serial
@pytest.mark.parametrize(
    "style_answer",
    (
        ("linear_ramp", [9428, 18525]),
        ("edge", [9724, 19078]),
        ("constant", [9428, 18525]),
    ),
)
def test_2dmesher_domain_extension(style_answer):
    style, answer = style_answer
    fname = os.path.join(os.path.dirname(__file__), "testing.segy")
    bbox = (-10e3, 0.0, 0.0, 10e3)
    wl = 5
    freq = 5.0
    hmin = 100
    hmax = 10e6
    grade = 0.005
    grad = 50.0
    rectangle = Rectangle(bbox)
    ef = get_sizing_function_from_segy(
        fname,
        bbox=bbox,
        grade=grade,
        grad=grad,
        stencil_size=100,
        wl=wl,
        freq=freq,
        hmin=hmin,
        hmax=hmax,
        pad_style=style,
        domain_pad=1e3,
    )
    plot_sizing_function(ef)
    write_velocity_model(fname, bbox=bbox)

    points, cells = generate_mesh(
        rectangle,
        ef,
        h0=hmin,
        perform_checks=True,
    )
    # import meshio

    # meshio.write_points_cells(
    #    "TEST.vtk", points, [("triangle", cells)], file_format="vtk"
    # )
    print(len(points), len(cells))
    assert allclose([len(points), len(cells)], answer, atol=100)


if __name__ == "__main__":
    test_2dmesher_domain_extension()
