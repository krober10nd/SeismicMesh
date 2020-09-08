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
    wl = 5
    hmin = 100
    grade = 0.005
    ef, bbox = get_sizing_function_from_segy(
        fname,
        bbox=(-10e3, 0, 0, 10e3),
        grade=grade,
        grad=50.0,
        wl=wl,
        hmin=hmin,
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
    assert np.abs(len(points) - 7906) < 20
    assert np.abs(len(cells) - 15487) < 20


if __name__ == "__main__":
    test_2dmesher()
