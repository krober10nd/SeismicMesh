import numpy

import termplotlib as tpl

__all__ = ["print_stats_3d"]


def print_stats_3d(
    angles, quality, method="not specified", elapsed=0.0, num_verts=0.0, num_cells=0.0
):
    """https://github.com/nschloe/optimesh/blob/master/optimesh/helpers.py"""

    angles = angles / numpy.pi * 180
    angles_hist, angles_bin_edges = numpy.histogram(
        angles, bins=numpy.linspace(0.0, 180.0, num=73, endpoint=True)
    )

    q = quality
    q_hist, q_bin_edges = numpy.histogram(
        q, bins=numpy.linspace(0.0, 1.0, num=41, endpoint=True)
    )

    grid = tpl.subplot_grid((2, 4))
    grid[0, 0].aprint(method)
    grid[0, 1].aprint("Mesh creation time (seconds):  {:7.2f}".format(elapsed))
    grid[0, 1].aprint(
        "Mesh creation speed (vertices/seconds):  {:7.2f}".format(num_verts / elapsed)
    )
    grid[0, 2].aprint("Number of vertices:      {:7.0f}".format(num_verts))
    grid[0, 2].aprint("Number of cells:         {:7.0f}".format(num_cells))
    grid[0, 3].aprint("")
    grid[1, 0].aprint("dihedral angles")
    grid[1, 0].hist(angles_hist, angles_bin_edges, grid=[15, 25])
    grid[1, 1].aprint("min angle:     {:7.3f}".format(numpy.min(angles)))
    grid[1, 1].aprint("avg angle:     {:7.3f}".format(numpy.mean(angles)))
    grid[1, 1].aprint("max angle:     {:7.3f}".format(numpy.max(angles)))
    grid[1, 1].aprint("std dev angle: {:7.3f}".format(numpy.std(angles)))
    grid[1, 2].aprint("mesh quality metric")
    grid[1, 2].hist(q_hist, q_bin_edges, grid=[15, 25])
    grid[1, 3].aprint("min quality: {:5.3f}".format(numpy.min(q)))
    grid[1, 3].aprint("avg quality: {:5.3f}".format(numpy.average(q)))
    grid[1, 3].aprint("max quality: {:5.3f}".format(numpy.max(q)))

    grid.show()
