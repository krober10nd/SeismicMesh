import pytest
import numpy as np

from SeismicMesh import (
    Rectangle,
    generate_mesh,
    get_sizing_function_from_segy,
    plot_sizing_function,
)


@pytest.mark.serial
def test_2dmesher_vs_water():

    fname = None
    bbox = (-10000.0, 0.0, 0.0, 10000.0)
    wl = 5
    freq = 2.0
    hmin = 10
    hmax = 10e6
    grade = 0.0
    grad = 0.0
    rectangle = Rectangle(bbox)
    nz = 200
    nx = 200
    vs = np.zeros((nz, nx))  # it will be replaced by 1500 m/s
    vs[0:150, :] = 1000  # m/s
    ef = get_sizing_function_from_segy(
        fname,
        bbox=bbox,
        grade=grade,
        grad=grad,
        wl=wl,
        freq=freq,
        hmin=hmin,
        hmax=hmax,
        velocity_data=vs,
        nz=nz,
        nx=nx,
    )
    plot_sizing_function(ef)
    assert ef.eval((-5000, 5000)) == 100
    assert ef.eval((-1, 5000)) == 150  # water layer
    ef.hmin = None  # to force h0 equal to 125 m (see below)
    points, cells = generate_mesh(
        rectangle,
        ef,
        h0=125,
        perform_checks=True,
    )
    # should have: 5616 vertices and 10955 cells
    print(len(points), len(cells))
    assert np.allclose([len(points), len(cells)], [5616, 10955], atol=100)

    # if True:
    #    import meshio

    #    meshio.write_points_cells(
    #       "blah.vtk",
    #       points[:, [1, 0]] / 1000,
    #       [("triangle", cells)],
    #       file_format="vtk",
    #    )


if __name__ == "__main__":
    test_2dmesher_vs_water()
