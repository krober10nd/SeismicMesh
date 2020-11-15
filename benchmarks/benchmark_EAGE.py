# benchmark the generation and quality of EAGE mesh
import time
import zipfile
import argparse
import numpy

import meshplex

# import meshio
import pygalmesh
import pygmsh
import gmsh

from SeismicMesh import (
    get_sizing_function_from_segy,
    Cube,
    generate_mesh,
    sliver_removal,
)


from helpers import print_stats_3d

# Bounding box describing domain extents (corner coordinates)
bbox = (-4200.0, 0.0, 0.0, 13520.0, 0.0, 13520.0)


def _build_sizing(HMIN=150.0, FREQ=2):

    # This file is in a big Endian binary format, so we must tell the program the shape of the velocity model.
    path = "Salt_Model_3D/3-D_Salt_Model/VEL_GRIDS/"
    # Extract binary file Saltf@@ from SALTF.ZIP
    zipfile.ZipFile(path + "SALTF.ZIP", "r").extract("Saltf@@", path=path)

    fname = path + "Saltf@@"

    # Dimensions of model (number of grid points in z, x, and y)
    nx, ny, nz = 676, 676, 210

    ef = get_sizing_function_from_segy(
        fname,
        bbox,
        hmin=HMIN,
        dt=0.001,
        freq=FREQ,
        wl=5,
        grade=0.15,
        hmax=5e3,
        nz=nz,
        nx=nx,
        ny=ny,
        byte_order="big",
        axes_order=[2, 0, 1],  # default order z, x, y -> order for EAGE x, y, z
        axes_order_sort="F",  # binary is packed in a FORTRAN-style
    )

    return ef


# for pytest-benchmark
# ef = _build_sizing()


def test_seismic_mesh(benchmark):
    angles, quality, elapsed, num_vertices, num_cells = benchmark.pedantic(
        run_SeismicMesh, args=(ef), iterations=1, rounds=5, warmup_rounds=0
    )
    assert numpy.amin(angles / numpy.pi * 180) > 10.0


def test_gmsh(benchmark):
    angles, quality, elapsed, num_vertices, num_cells = benchmark.pedantic(
        run_gmsh, args=(ef), iterations=1, rounds=5, warmup_rounds=0
    )
    assert numpy.amin(angles / numpy.pi * 180) > 10.0


def test_cgal(benchmark):
    angles, quality, elapsed, num_vertices, num_cells = benchmark.pedantic(
        run_cgal, args=(ef), iterations=1, rounds=5, warmup_rounds=0
    )
    assert numpy.amin(angles / numpy.pi * 180) > 10.0


def run_cgal(ef, HMIN=75.0):

    print("generating a mesh with cgal...")
    t1 = time.time()
    mesh = pygalmesh.generate_mesh(
        pygalmesh.Cuboid([-4200.0, 0.0, 0.0], [0.0, 13520.0, 13520.0]),
        facet_angle=30,
        cell_radius_edge_ratio=2.0,
        cell_size=lambda x: ef.eval(x) / 1.1,
        edge_size=HMIN,
    )
    elapsed = time.time() - t1

    # mesh.write("cgal_EAGE.vtk")

    plex = meshplex.MeshTetra(mesh.points, mesh.cells[1][1])
    angles = plex.q_min_sin_dihedral_angles
    quality = plex.q_radius_ratio

    num_cells = len(mesh.cells[1][1])
    num_vertices = len(mesh.points)

    return angles, quality, elapsed, num_vertices, num_cells


def run_gmsh(ef, HMIN=75.0):
    with pygmsh.geo.Geometry() as geom:

        geom.add_box(-4200.0, 0.0, 0.0, 13520.0, 0.0, 13520.0, HMIN)

        geom.set_mesh_size_callback(
            lambda dim, tag, x, y, z: (ef.eval([x, y, z])) / 1.1
        )
        gmsh.option.setNumber("Mesh.CharacteristicLengthExtendFromBoundary", 0)

        print("generating a mesh with gmsh...")
        t1 = time.time()
        mesh = geom.generate_mesh()
        elapsed = time.time() - t1

        points = mesh.points
        cells = mesh.cells[2].data

        num_cells = len(cells)
        num_vertices = len(points)

        plex = meshplex.MeshTetra(points, cells)
        angles = plex.q_min_sin_dihedral_angles
        quality = plex.q_radius_ratio

        # mesh.write("gmsh_EAGE.vtk")

        return angles, quality, elapsed, num_vertices, num_cells


def run_SeismicMesh(ef, HMIN=75.0):

    cube = Cube(bbox)

    t1 = time.time()
    points, cells = generate_mesh(
        domain=cube,
        edge_length=ef,
        max_iter=25,
    )

    points, cells = sliver_removal(
        points=points,
        domain=cube,
        edge_length=ef,
    )
    elapsed = time.time() - t1

    # meshio.write_points_cells(
    #    "sm_EAGE.vtk",
    #    points,
    #    [("tetra", cells)],
    #    file_format="vtk",
    # )

    plex = meshplex.MeshTetra(points, cells)
    angles = plex.q_min_sin_dihedral_angles
    quality = plex.q_radius_ratio

    num_cells = len(cells)
    num_vertices = len(points)

    return angles, quality, elapsed, num_vertices, num_cells


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument(
        "--method",
        dest="method",
        type=str,
        default=None,
        required=False,
        help="Run benchmark with method=('cgal','sm', 'gmsh')",
    )

    args = parser.parse_args()

    ef = _build_sizing()
    if args.method == "cgal":
        a1, q1, t1, nv, nc = run_cgal(ef)
        print_stats_3d(a1, q1, "cgal", t1, nv, nc)
    elif args.method == "sm":
        a1, q1, t1, nv, nc = run_SeismicMesh(ef)
        print_stats_3d(a1, q1, "SeismicMesh", t1, nv, nc)
    elif args.method == "gmsh":
        a1, q1, t1, nv, nc = run_gmsh(ef)
        print_stats_3d(a1, q1, "gmsh", t1, nv, nc)
    else:
        # a1, q1, t1, nv1, nc1 = run_cgal(ef)
        a2, q2, t2, nv2, nc2 = run_SeismicMesh(ef)
        a3, q3, t3, nv3, nc3 = run_gmsh(ef)
        # print_stats_3d(a1, q1, "CGAL", t1, nv1, nc1)
        print_stats_3d(a2, q2, "SeismicMesh", t2, nv2, nc2)
        print_stats_3d(a3, q3, "gmsh", t3, nv3, nc3)
