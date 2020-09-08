import os

import numpy as np
import pytest

from SeismicMesh import (
    generate_mesh,
    geometry,
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
    ef, bbox = get_sizing_function_from_segy(
        fname,
        bbox=bbox,
        grade=grade,
        grad=grad,
        wl=wl,
        freq=freq,
        hmin=hmin,
        hmax=hmax,
    )
    plot_sizing_function(ef, bbox)
    write_velocity_model(fname)

    def rectangle(p):
        return geometry.drectangle(p, *bbox)

    points, cells = generate_mesh(
        bbox=bbox,
        signed_distance_function=rectangle,
        cell_size=ef,
        h0=hmin,
        perform_checks=True,
    )
    # should have: 7690 vertices and 15045 cells
    print(len(points), len(cells))
    assert np.abs(len(points) - 7690) < 20
    assert np.abs(len(cells) - 15045) < 20

    # import meshio

    # meshio.write_points_cells(
    #    "blah.vtk",
    #    points[:, [1, 0]] / 1000,
    #    [("triangle", cells)],
    #    file_format="vtk",
    # )


if __name__ == "__main__":
    test_2dmesher()
