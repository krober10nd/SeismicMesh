# benchmark the generation and quality of BP2004 mesh
import time
import argparse
import numpy

import meshplex
import meshio
import pygmsh
import gmsh

from SeismicMesh import (
    get_sizing_function_from_segy,
    Rectangle,
    generate_mesh,
    geometry,
)


from helpers import print_stats_2d

HMIN = 75.0


gmsh.option.setNumber("Mesh.Algorithm", 5)


def run_gmsh(ef):
    with pygmsh.geo.Geometry() as geom:

        geom.add_polygon(
            [
                [-12e3, 0.0],
                [-12e3, 67e3],
                [0.0, 67e3],
                [0.0, 0.0],
            ]
        )

        geom.set_mesh_size_callback(lambda dim, tag, x, y, z: (ef.eval([x, y])))

        t1 = time.time()
        mesh = geom.generate_mesh()
        elapsed = time.time() - t1

        points = mesh.points
        cells = mesh.cells[1].data

        num_cells = len(cells)
        num_vertices = len(points)

        plex = meshplex.MeshTri(points, cells)
        quality = numpy.abs(plex.cell_quality)

        mesh.write("BP2004_gmsh.vtk")
        return quality, elapsed, num_vertices, num_cells


def run_SeismicMesh(ef):

    bbox = (-12000, 0.0, 0.0, 67000.0)

    rectangle = Rectangle(bbox)

    t1 = time.time()
    points, cells = generate_mesh(
        domain=rectangle,
        cell_size=ef,
        h0=HMIN,
        max_iter=25,
        delta_t=0.3,
        nscreen=10,
    )
    points, cells = geometry.delete_boundary_entities(
        points, cells, dim=2, min_qual=0.10
    )
    elapsed = time.time() - t1

    meshio.write_points_cells(
        "BP2004_sm.vtk",
        points,
        [("triangle", cells)],
        file_format="vtk",
    )

    plex = meshplex.MeshTri(points, cells)
    quality = numpy.abs(plex.cell_quality)

    num_cells = len(cells)
    num_vertices = len(points)

    return quality, elapsed, num_vertices, num_cells


def _build_sizing():

    # Name of SEG-Y file containg velocity model.
    fname = "vel_z6.25m_x12.5m_exact.segy"

    # Bounding box describing domain extents (corner coordinates)
    bbox = (-12000, 0.0, 0.0, 67000.0)

    # Construct mesh sizing object from velocity model
    ef = get_sizing_function_from_segy(
        fname,
        bbox,
        hmin=HMIN,
        wl=10,
        freq=2,
        dt=0.001,
        grade=0.15,
    )
    return ef


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
        print("not done")
        # q1, t1, nv, nc = run_cgal()
        # print_stats_2d(q1, "cgal", t1, nv, nc)
    elif args.method == "sm":
        q1, t1, nv, nc = run_SeismicMesh(ef)
        print_stats_2d(q1, "SeismicMesh", t1, nv, nc)
    elif args.method == "gmsh":
        q1, t1, nv, nc = run_gmsh(ef)
        print_stats_2d(q1, "gmsh", t1, nv, nc)
    else:
        # q1, t1, nv1, nc1 = run_cgal()
        q2, t2, nv2, nc2 = run_SeismicMesh(ef)
        q3, t3, nv3, nc3 = run_gmsh(ef)
        # print_stats_2d(q1, "CGAL", t1, nv1, nc1)
        print_stats_2d(q2, "SeismicMesh", t2, nv2, nc2)
        print_stats_2d(q3, "gmsh", t3, nv3, nc3)