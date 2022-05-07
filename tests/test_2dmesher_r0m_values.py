import pytest
import numpy as np

from SeismicMesh import (
    Rectangle,
    generate_mesh,
    get_sizing_function_from_segy,
    plot_sizing_function,
)


@pytest.mark.serial
def test_2dmesher_r0m_values():

    fname = None
    bbox = (-10000.0, 0.0, 0.0, 10000.0)
    wl = 5
    freq = 2.0
    hmin = 75
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
    assert ef.eval((-5000, 5000)) == 100  # hmin* that appears in mesh sizing
    assert ef.eval((-1, 5000)) == 150     # h in the water layer
    
    # This tests if r0m will be equal to min(r0), where min(r0) = hmin*,  
    # and hmin* is the minimum value that appears in mesh sizing.
    # Note that hmin* >= hmin (hmin is given to get_sizing_function_from_segy).
    # This is the default behavior of generate_mesh, and could generate
    # elements smaller than hmin*.
    # h0 is used to generate the initial set of points:
    #   if ef.hmin = None, h0 will be defined by the user;
    #   if ef.hmin != None, h0 will be equal to ef.hmin.
    # Here, r0m should be equal to 100.
    points, cells = generate_mesh(
        rectangle,
        ef,
        h0=125, # since ef.hmin = 75 (i.e., it is not None), h0 will be equal to ef.hmin = 75 
        perform_checks=True,
        r0m_is_h0 = False, # by default, r0m will be min(r0)=himn*; in this case, r0m = 100
    )
    # should have: 15460 vertices and 30464 cells
    print(len(points), len(cells))
    assert np.allclose([len(points), len(cells)], [15460, 30464], atol=100)

    #if True:
    #   import meshio

    #   meshio.write_points_cells(
    #      "blah1.vtk",
    #      points[:, [1, 0]] / 1000,
    #      [("triangle", cells)],
    #      file_format="vtk",
    #   )
    
    # This tests if r0m will be equal to h0, i.e., the size used to
    # generate the initial set of points.
    # By forcing r0m=h0 would generate elements with sizes similar 
    # to mesh sizing data (or maybe a little bit larger than mesh sizing).
    # h0 is used to generate the initial set of points:
    #   if ef.hmin = None, h0 will be defined by the user;
    #   if ef.hmin != None, h0 will be equal to ef.hmin.
    # Here, r0m should be equal to 75.
    # Ideally, r0m should be between h0 and hmin*.
    points, cells = generate_mesh(
        rectangle,
        ef,
        h0=125, # since ef.hmin = 75 (i.e., it is not None), h0 will be equal to ef.hmin = 75 
        perform_checks=True,
        r0m_is_h0 = True, # with True, r0m will be equal to h0; in this case, r0m = 75
    )
    # should have: 10566 vertices and 20771 cells
    print(len(points), len(cells))
    assert np.allclose([len(points), len(cells)], [10566, 20771], atol=100)

    #if True:
    #   import meshio

    #   meshio.write_points_cells(
    #      "blah2.vtk",
    #      points[:, [1, 0]] / 1000,
    #      [("triangle", cells)],
    #      file_format="vtk",
    #   )
    
    # This tests if r0m will be equal to h0, i.e., the size used to
    # generate the initial set of points.
    # By forcing r0m=h0 would generate elements with sizes similar 
    # to mesh sizing data (or maybe a little bit larger than mesh sizing).
    # h0 is used to generate the initial set of points:
    #   if ef.hmin = None, h0 will be defined by the user;
    #   if ef.hmin != None, h0 will be equal to ef.hmin.
    # Here, we force hmin*=hmin=100, therefore
    # r0m should be equal to 100.
    # Ideally, r0m should be between h0 and hmin*.
    ef.hmin = 100
    points, cells = generate_mesh(
        rectangle,
        ef,
        h0=125, # since ef.hmin = 100 (i.e., it is not None), h0 will be equal to ef.hmin = 100 
        perform_checks=True,
        r0m_is_h0 = True, # with True, r0m will be equal to h0; in this case, r0m = 100
    )
    # should have: 8739 vertices and 17131 cells
    print(len(points), len(cells))
    assert np.allclose([len(points), len(cells)], [8739, 17131], atol=100)
    
    #if True:
    #   import meshio

    #   meshio.write_points_cells(
    #      "blah3.vtk",
    #      points[:, [1, 0]] / 1000,
    #      [("triangle", cells)],
    #      file_format="vtk",
    #   )

if __name__ == "__main__":
    test_2dmesher_r0m_values()
