import time
import argparse
import numpy
from mpi4py import MPI

import meshplex
import pygalmesh
import termplotlib as tpl
import meshio
import SeismicMesh


def print_stats(angles, quality, extra_cols=None):
    """https://github.com/nschloe/optimesh/blob/master/optimesh/helpers.py"""
    extra_cols = [] if extra_cols is None else extra_cols

    angles = angles / numpy.pi * 180
    angles_hist, angles_bin_edges = numpy.histogram(
        angles, bins=numpy.linspace(0.0, 180.0, num=73, endpoint=True)
    )

    q = quality
    q_hist, q_bin_edges = numpy.histogram(
        q, bins=numpy.linspace(0.0, 1.0, num=41, endpoint=True)
    )

    grid = tpl.subplot_grid((1, 4))
    grid[0, 0].hist(angles_hist, angles_bin_edges, bar_width=1, strip=True)
    grid[0, 1].aprint("min angle:     {:7.3f}".format(numpy.min(angles)))
    grid[0, 1].aprint("avg angle:     {:7.3f}".format(60))
    grid[0, 1].aprint("max angle:     {:7.3f}".format(numpy.max(angles)))
    grid[0, 1].aprint("std dev angle: {:7.3f}".format(numpy.std(angles)))
    grid[0, 2].hist(q_hist, q_bin_edges, bar_width=1, strip=True)
    grid[0, 3].aprint("min quality: {:5.3f}".format(numpy.min(q)))
    grid[0, 3].aprint("avg quality: {:5.3f}".format(numpy.average(q)))
    grid[0, 3].aprint("max quality: {:5.3f}".format(numpy.max(q)))
    for k, col in enumerate(extra_cols):
        grid[0, 4 + k].aprint(col)

    grid.show()


HMIN = 0.05


def run_cgal():
    # CGAL

    class Field(pygalmesh.SizingFieldBase):
        def eval(self, x):
            return HMIN

    t1 = time.time()
    mesh = pygalmesh.generate_mesh(
        pygalmesh.Cuboid([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]),
        cell_size=Field(),
        edge_size=HMIN,
    )
    print(time.time() - t1)

    mesh.write("cgal.vtk")

    plex = meshplex.MeshTetra(mesh.points, mesh.cells[1][1])
    angles = plex.q_min_sin_dihedral_angles
    quality = plex.q_radius_ratio

    print_stats(angles, quality)


def run_SeismicMesh():
    comm = MPI.COMM_WORLD

    # SeismicMesh
    bbox = (0.0, 1.0, 0.0, 1.0, 0.0, 1.0)

    cube = SeismicMesh.Cube(bbox)

    def fh(x):
        # Note: for parallel execution this logic is required
        # since the decomposition of the sizing function passes a tuple to fh
        if type(x) == tuple:
            h = numpy.zeros_like(x[0]) + HMIN
        else:
            h = numpy.array([HMIN] * len(x))
        return h

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
    print(time.time() - t1)

    if comm.rank == 0:

        plex = meshplex.MeshTetra(points, cells)
        angles = plex.q_min_sin_dihedral_angles
        quality = plex.q_radius_ratio

        print_stats(angles, quality)

        meshio.write_points_cells(
            "SeismicMesh_cuboid.vtk",
            points,
            [("tetra", cells)],
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument(
        "--method",
        dest="method",
        type=str,
        default=None,
        required=True,
        help="Run benchmark with method=('cgal','sm')",
    )

    args = parser.parse_args()

    if args.method == 'cgal':
        run_cgal()
    if args.method == 'sm':
        run_SeismicMesh()
