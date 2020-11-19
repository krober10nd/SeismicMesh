import os

import numpy as np
import pytest
from mpi4py import MPI

from SeismicMesh import (
    Cube,
    generate_mesh,
    geometry,
    get_sizing_function_from_segy,
    sliver_removal,
)

comm = MPI.COMM_WORLD


@pytest.mark.parallel3
def test_3dmesher_par():

    fname = os.path.join(os.path.dirname(__file__), "test3D.bin")
    nz, nx, ny = 20, 10, 10
    bbox = (-2e3, 0.0, 0.0, 1e3, 0.0, 1e3)

    hmin = 50
    wl = 10
    freq = 4
    grade = 0.15
    grad = 50.0
    cube = Cube(bbox)
    ef = get_sizing_function_from_segy(
        fname,
        bbox,
        hmin=hmin,
        grade=grade,
        grad=grad,
        wl=wl,
        freq=freq,
        nz=nz,
        nx=nx,
        ny=ny,
        byte_order="little",
        axes_order=(2, 0, 1),  # order for EAGE (x, y, z) to default order (z,x,y)
    )

    points, cells = generate_mesh(edge_length=ef, domain=cube, verbose=2)

    points, cells = sliver_removal(
        points=points,
        domain=cube,
        edge_length=ef,
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
        print(len(points), len(cells))
        assert np.abs(9220 - len(points)) < 5000
        assert np.abs(49156 - len(cells)) < 5000


if __name__ == "__main__":
    test_3dmesher_par()
