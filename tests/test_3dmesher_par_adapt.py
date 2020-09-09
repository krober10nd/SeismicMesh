import os

import numpy as np
import pytest
from mpi4py import MPI

from SeismicMesh import get_sizing_function_from_segy, generate_mesh, geometry

comm = MPI.COMM_WORLD


@pytest.mark.parallel3
def test_3dmesher_par_adapt():
    fname = os.path.join(os.path.dirname(__file__), "test3D.bin")
    nz, nx, ny = 20, 10, 10
    bbox = (-2e3, 0.0, 0.0, 1e3, 0.0, 1e3)
    hmin = 50
    wl = 10
    freq = 4
    grade = 0.15
    ef, bbox = get_sizing_function_from_segy(
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

    def cube(p):
        return geometry.dblock(p, *bbox)

    points, cells = generate_mesh(
        bbox=bbox,
        h0=hmin,
        cell_size=ef,
        signed_distance_function=cube,
        max_iter=10,
        perform_checks=False,
    )

    points = comm.bcast(points, 0)

    points, cells = generate_mesh(
        points=points,
        bbox=bbox,
        h0=hmin,
        cell_size=ef,
        signed_distance_function=cube,
        axis=1,
        max_iter=10,
        perform_checks=False,
    )

    points = comm.bcast(points, 0)

    points, cells = generate_mesh(
        points=points,
        bbox=bbox,
        h0=hmin,
        cell_size=ef,
        signed_distance_function=cube,
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
