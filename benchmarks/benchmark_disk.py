import time
import argparse
import numpy

import meshplex

# import meshio

import pygalmesh
import pygmsh
import SeismicMesh

from helpers import print_stats_2d


def test_seismic_mesh(benchmark):
    quality, elapsed, num_vertices, num_cells = benchmark.pedantic(
        run_SeismicMesh, iterations=1, rounds=5, warmup_rounds=0
    )
    assert numpy.amin(quality) > 0.10


def test_gmsh(benchmark):
    quality, elapsed, num_vertices, num_cells = benchmark.pedantic(
        run_gmsh, iterations=1, rounds=5, warmup_rounds=0
    )
    assert numpy.amin(quality) > 0.10


def test_cgal(benchmark):
    quality, elapsed, num_vertices, num_cells = benchmark.pedantic(
        run_cgal, iterations=1, rounds=5, warmup_rounds=0
    )
    assert numpy.amin(quality) > 0.10


def run_gmsh(HMIN=0.01):

    with pygmsh.geo.Geometry() as geom:
        geom.add_circle([0.0, 0.0], 1.0, mesh_size=HMIN)
        mesh = geom.generate_mesh()

        t1 = time.time()
        mesh = geom.generate_mesh()
        elapsed = time.time() - t1

    # mesh.write("gmsh_circle.vtk")

    points = mesh.points
    cells = mesh.cells[1].data

    num_cells = len(cells)
    num_vertices = len(points)

    plex = meshplex.MeshTri(points, cells)
    quality = numpy.abs(plex.cell_quality)

    return quality, elapsed, num_vertices, num_cells


def run_cgal(HMIN=0.01):

    n = 50
    points = numpy.array(
        [
            [numpy.cos(alpha), numpy.sin(alpha)]
            for alpha in numpy.linspace(0.0, 2 * numpy.pi, n, endpoint=False)
        ]
    )
    constraints = [[k, k + 1] for k in range(n - 1)] + [[n - 1, 0]]

    t1 = time.time()
    mesh = pygalmesh.generate_2d(points, constraints, edge_size=HMIN)
    elapsed = time.time() - t1

    # mesh.write("cgal_circle.vtk")

    points = mesh.points
    cells = mesh.cells[0].data

    num_cells = len(cells)
    num_vertices = len(points)

    plex = meshplex.MeshTri(points, cells)
    quality = numpy.abs(plex.cell_quality)

    return quality, elapsed, num_vertices, num_cells


def run_SeismicMesh(HMIN=0.01):

    bbox = (-1.0, 1.0, -1.0, 1.0)

    disk = SeismicMesh.geometry.Disk([0.0, 0.0], 1.0)

    def fh(x):
        return numpy.array([HMIN] * len(x))

    t1 = time.time()
    points, cells = SeismicMesh.generate_mesh(
        bbox=bbox,
        h0=HMIN,
        domain=disk,
        edge_length=fh,
        max_iter=25,
        verbose=0,
    )
    elapsed = time.time() - t1

    # meshio.write_points_cells(
    #    "SeismicMesh_circle.vtk",
    #    points,
    #    [("triangle", cells)],
    # )

    plex = meshplex.MeshTri(points, cells)
    quality = numpy.abs(plex.cell_quality)

    num_cells = len(cells)
    num_vertices = len(points)

    return quality, elapsed, num_vertices, num_cells


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

    if args.method == "cgal":
        q1, t1, nv, nc = run_cgal()
        print_stats_2d(q1, "cgal", t1, nv, nc)
    elif args.method == "sm":
        q1, t1, nv, nc = run_SeismicMesh()
        print_stats_2d(q1, "SeismicMesh", t1, nv, nc)
    elif args.method == "gmsh":
        q1, t1, nv, nc = run_gmsh()
        print_stats_2d(q1, "gmsh", t1, nv, nc)
    else:
        q1, t1, nv1, nc1 = run_cgal()
        q2, t2, nv2, nc2 = run_SeismicMesh()
        q3, t3, nv3, nc3 = run_gmsh()
        print_stats_2d(q1, "CGAL", t1, nv1, nc1)
        print_stats_2d(q2, "SeismicMesh", t2, nv2, nc2)
        print_stats_2d(q3, "gmsh", t3, nv3, nc3)
