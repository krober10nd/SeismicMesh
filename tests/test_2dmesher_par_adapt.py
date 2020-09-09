import os

import numpy as np
import pytest
from mpi4py import MPI

from SeismicMesh import (
    get_sizing_function_from_segy,
    generate_mesh,
    geometry,
    Rectangle,
)

comm = MPI.COMM_WORLD


@pytest.mark.parallel2
def test_2dmesher_par_adapt():
    fname = os.path.join(os.path.dirname(__file__), "testing.segy")
    bbox = (-10e3, 0.0, 0.0, 10e3)
    wl = 5
    freq = 5
    hmin = 100
    grade = 0.005
    ef = get_sizing_function_from_segy(
        fname,
        bbox,
        hmin=hmin,
        wl=wl,
        freq=freq,
        grade=grade,
    )

    rectangle = Rectangle(bbox)

    points, cells = generate_mesh(
        cell_size=ef,
        domain=rectangle,
        h0=hmin,
        perform_checks=False,
    )
    points = comm.bcast(points, 0)

    # pass the points and restart with a different axis
    points, cells = generate_mesh(
        points=points,
        cell_size=ef,
        h0=hmin,
        domain=rectangle,
        perform_checks=False,
    )

    if comm.rank == 0:
        import meshio

        meshio.write_points_cells(
            "test2d.vtk", points / 1000, [("triangle", cells)], file_format="vtk"
        )
        area = geometry.simp_vol(points / 1000, cells)
        assert np.abs(100 - np.sum(area)) < 0.50  # km2


if __name__ == "__main__":
    test_2dmesher_par_adapt()
