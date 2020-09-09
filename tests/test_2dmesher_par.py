import os

import numpy as np
import pytest
from mpi4py import MPI

from SeismicMesh import (
    get_sizing_function_from_segy,
    generate_mesh,
    Rectangle,
    geometry,
)

comm = MPI.COMM_WORLD


@pytest.mark.parallel2
def test_2dmesher_par():

    fname = os.path.join(os.path.dirname(__file__), "testing.segy")
    bbox = (-10e3, 0.0, 0.0, 10e3)
    freq = 2
    wl = 10
    hmin = 75
    grade = 0.005

    rectangle = Rectangle(bbox)
    ef = get_sizing_function_from_segy(
        fname, bbox, hmin=hmin, wl=wl, freq=freq, grade=grade
    )

    points, cells = generate_mesh(
        rectangle,
        ef,
        hmin,
        max_iter=100,
        perform_checks=False,
    )

    if comm.rank == 0:
        import meshio

        meshio.write_points_cells(
            "test2d.vtk",
            points / 1000,
            [("triangle", cells)],
            file_format="vtk",
        )
        area = geometry.simp_vol(points / 1000, cells)
        assert np.abs(100 - np.sum(area)) < 0.50  # km2


if __name__ == "__main__":
    test_2dmesher_par()
