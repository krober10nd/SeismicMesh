import os

from numpy import allclose
import pytest

from SeismicMesh import (
    generate_mesh,
    Rectangle,
    get_sizing_function_from_segy,
    plot_sizing_function,
    write_velocity_model,
)


@pytest.mark.serial
def test_2dmesher():

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
        wl=wl,
        freq=freq,
        hmin=hmin,
        hmax=hmax,
    )
    plot_sizing_function(ef)
    write_velocity_model(fname)

    points, cells = generate_mesh(
        rectangle,
        ef,
        hmin,
        perform_checks=True,
    )
    # should have: 7690 vertices and 15045 cells
    print(len(points), len(cells))
    allclose([len(points), len(cells)], [7690, 15045], atol=100)

    # import meshio

    # meshio.write_points_cells(
    #    "blah.vtk",
    #    points[:, [1, 0]] / 1000,
    #    [("triangle", cells)],
    #    file_format="vtk",
    # )


if __name__ == "__main__":
    test_2dmesher()
