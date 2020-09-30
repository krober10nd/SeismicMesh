from math import sqrt
import time
import argparse
import numpy

import meshplex
import meshio

import pygalmesh
import pygmsh
import SeismicMesh

from helpers import print_stats_3d


HMIN = 0.025


def run_gmsh():
    with pygmsh.occ.Geometry() as geom:
        geom.add_ball([0.0, 0.0, 0.0], 1.0)

        geom.set_mesh_size_callback(
            lambda dim, tag, x, y, z: abs(sqrt(x ** 2 + y ** 2 + z ** 2) - 0.5) + HMIN
        )
        t1 = time.time()
        mesh = pygmsh.generate_mesh(geom)
        elapsed = time.time() - t1

    mesh.write("gmsh_sphere.vtk")

    plex = meshplex.MeshTetra(mesh.points, mesh.cells[1][1])
    angles = plex.q_min_sin_dihedral_angles
    quality = plex.q_radius_ratio

    return angles, quality, elapsed


def run_cgal():
    class Field(pygalmesh.SizingFieldBase):
        def eval(self, x):
            return abs(numpy.sqrt(numpy.dot(x, x)) - 0.5) / 5 + HMIN

    t1 = time.time()
    mesh = pygalmesh.generate_mesh(
        pygalmesh.Ball([0.0, 0.0, 0.0], 1.0),
        cell_size=Field(),
        edge_size=HMIN,
    )
    elapsed = time.time() - t1

    mesh.write("cgal_sphere.vtk")

    plex = meshplex.MeshTetra(mesh.points, mesh.cells[1][1])
    angles = plex.q_min_sin_dihedral_angles
    quality = plex.q_radius_ratio

    num_cells = len(mesh.cells[1][1])
    num_vertices = len(mesh.points)

    return angles, quality, elapsed, num_cells, num_vertices


def run_SeismicMesh():

    bbox = (-1.0, 1.0, -1.0, 1.0, -1.0, 1.0)

    def dsphere(p, xc, yc, zc, r):
        """Signed distance function for a sphere centered at xc,yc,zc with radius  r."""
        return (
            numpy.sqrt((p[:, 0] - xc) ** 2 + (p[:, 1] - yc) ** 2 + (p[:, 2] - zc) ** 2)
            - r
        )

    def sphere(x):
        return dsphere(x, 0, 0, 0, 1)

    def fh(x):
        return numpy.abs(numpy.sqrt(numpy.einsum("ij, ij->i", x, x)) - 0.5) / 5 + HMIN

    t1 = time.time()
    points, cells = SeismicMesh.generate_mesh(
        bbox=bbox, h0=HMIN, domain=sphere, cell_size=fh, nscreen=10, max_iter=35
    )
    points, cells = SeismicMesh.sliver_removal(
        points=points,
        bbox=bbox,
        h0=HMIN,
        domain=sphere,
        cell_size=fh,
        min_dh_angle_bound=10,
        nscreen=10,
    )
    elapsed = time.time() - t1

    meshio.write_points_cells(
        "SeismicMesh_sphere.vtk",
        points,
        [("tetra", cells)],
    )

    plex = meshplex.MeshTetra(points, cells)
    angles = plex.q_min_sin_dihedral_angles
    quality = plex.q_radius_ratio

    num_cells = len(cells)
    num_vertices = len(points)

    return angles, quality, elapsed, num_cells, num_vertices


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
        a1, q1, t1, nv, nc = run_cgal()
        print_stats_3d(a1, q1, "CGAL", t1, nv, nc)
    elif args.method == "sm":
        a1, q1, t1, nv, nc = run_SeismicMesh()
        print_stats_3d(a1, q1, "SeismicMesh", t1, nv, nc)
    elif args.method == "gmsh":
        a1, q1, t1, nv, nc = run_gmsh()
        print_stats_3d(a1, q1, "gmsh", t1, nv, nc)
    else:
        a1, q1, t1, nv1, nc1 = run_cgal()
        a2, q2, t2, nv2, nc2 = run_SeismicMesh()
        a3, q3, t3, nv3, nc3 = run_gmsh()
        print_stats_3d(a1, q1, "CGAL", t1, nv1, nc1)
        print_stats_3d(a2, q2, "SeismicMesh", t2, nv2, nc2)
        print_stats_3d(a3, q3, "gmsh", t3, nv3, nc3)
