import os

import numpy as np
import pytest

from SeismicMesh import (
    generate_mesh,
    geometry,
    get_sizing_function_from_segy,
    write_velocity_model,
)


@pytest.mark.serial
def test_3dmesher():

    fname = os.path.join(os.path.dirname(__file__), "test3D.bin")

    wl = 10
    hmin = 50
    freq = 2
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
    )

    write_velocity_model(fname, nz=20, nx=10, ny=10)

    def cube(p):
        return geometry.dblock(p, *bbox)

    points, cells = generate_mesh(
        bbox=bbox, signed_distance_function=cube, cell_size=ef, h0=hmin
    )
    print(len(points), len(cells))
    # 16459 90552

    assert len(points) == 16459
    assert np.abs(len(cells) - 90552) < 600

    import meshio

    meshio.write_points_cells("food3D.vtk", points, [("tetra", cells)])


if __name__ == "__main__":
    test_3dmesher()
