import os

import numpy as np
import pytest
from mpi4py import MPI

from SeismicMesh import get_sizing_function_from_segy, generate_mesh, Cube, geometry

comm = MPI.COMM_WORLD


@pytest.mark.parallel3
def test_3dmesher_par_adapt():
    fname = os.path.join(os.path.dirname(__file__), "test3D.bin")
    nz, nx, ny = 20, 10, 10
    bbox = (-2e3, 0.0, 0.0, 1e3, 0.0, 1e3)
    cube = Cube(bbox)
    hmin = 50
    wl = 10
    freq = 4
    grade = 0.15
    ef = get_sizing_function_from_segy(
        fname,
        bbox,
        hmin=hmin,
        wl=wl,
        freq=freq,
        grade=grade,
        nx=nx,
        ny=ny,
        nz=nz,
        byte_order="little",
    )

    points, cells = generate_mesh(
        cube,
        ef,
        hmin,
        max_iter=10,
        perform_checks=False,
    )

    points = comm.bcast(points, 0)

    points, cells = generate_mesh(
        points=points,
        domain=cube,
        h0=hmin,
        cell_size=ef,
        axis=1,
        max_iter=10,
        perform_checks=False,
    )

    points = comm.bcast(points, 0)

    points, cells = generate_mesh(
        points=points,
        h0=hmin,
        cell_size=ef,
        domain=cube,
        axis=2,
        max_iter=10,
        perform_checks=False,
    )

    if comm.rank == 0:
        import meshio

        meshio.write_points_cells(
            "foo3D_V3.vtk",
            points,
            [("tetra", cells)],
        )

        vol = geometry.simp_vol(points / 1000, cells)
        assert np.abs(2 - np.sum(vol)) < 0.10  # km2


if __name__ == "__main__":
    test_3dmesher_par_adapt()
