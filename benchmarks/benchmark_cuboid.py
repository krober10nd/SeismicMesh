import time
import argparse

import meshplex
import meshio

import pygalmesh
import SeismicMesh

from helpers import print_stats_3d

HMIN = 0.025
GRADE = HMIN / 3.0


def run_gmsh():

    return angles, quality, num_verts, num_cells


def run_cgal():
    class Field(pygalmesh.SizingFieldBase):
        def eval(self, x):
            return (1 - x[0]) * GRADE + HMIN

    t1 = time.time()
    mesh = pygalmesh.generate_mesh(
        pygalmesh.Cuboid([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]),
        cell_size=Field(),
        edge_size=HMIN,
    )
    elapsed = time.time() - t1

    mesh.write("cgal_cuboid.vtk")

    plex = meshplex.MeshTetra(mesh.points, mesh.cells[1][1])
    angles = plex.q_min_sin_dihedral_angles
    quality = plex.q_radius_ratio
    num_verts = len(mesh.points)
    num_cells = len(mesh.cells[1][1])

    return angles, quality, elapsed, num_verts, num_cells


def run_SeismicMesh():

    bbox = (0.0, 1.0, 0.0, 1.0, 0.0, 1.0)

    cube = SeismicMesh.Cube(bbox)

    def fh(x):
        return (1 - x[:, 0]) * GRADE + HMIN

    t1 = time.time()
    points, cells = SeismicMesh.generate_mesh(
        bbox=bbox, h0=HMIN, domain=cube, cell_size=fh, nscreen=10, max_iter=30
    )
    points, cells = SeismicMesh.sliver_removal(
        points=points,
        bbox=bbox,
        h0=HMIN,
        domain=cube,
        cell_size=fh,
        min_dh_angle_bound=10,
        nscreen=10,
    )
    elapsed = time.time() - t1

    meshio.write_points_cells(
        "SeismicMesh_cuboid.vtk",
        points,
        [("tetra", cells)],
    )

    plex = meshplex.MeshTetra(points, cells)
    angles = plex.q_min_sin_dihedral_angles
    quality = plex.q_radius_ratio

    num_verts = len(points)
    num_cells = len(cells)

    return angles, quality, elapsed, num_verts, num_cells


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
        help="Run benchmark with method=('cgal','sm')",
    )

    args = parser.parse_args()

    if args.method == "cgal":
        a1, q1, t1, nv, nc = run_cgal()
        print_stats_3d(a1, q1, "CGAL", t1, nv, nc)
    elif args.method == "sm":
        a1, q1, t1, nv, nc = run_SeismicMesh()
        print_stats_3d(a1, q1, "SeismicMesh", t1, nv, nc)
    else:
        a1, q1, t1, nv1, nc1 = run_cgal()
        a2, q2, t2, nv2, nc2 = run_SeismicMesh()
        print_stats_3d(a1, q1, "CGAL", t1, nv1, nc1)
        print_stats_3d(a2, q2, "SeismicMesh", t2, nv2, nc2)
