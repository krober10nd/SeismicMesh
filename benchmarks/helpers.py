import math

import numpy

import termplotlib as tpl

__all__ = ["print_stats_3d", "print_stats_2d"]


def print_stats_2d(
    quality, method="not specified", elapsed=0.0, num_verts=0.0, num_cells=0.0
):
    """https://github.com/nschloe/optimesh/blob/master/optimesh/helpers.py
    NTS: terminal dimension: 286x93, font size 11, Menlo Regular, Line spacing 0.8, Char. spacing 1.044
    """

    q = quality
    q_hist, q_bin_edges = numpy.histogram(
        q, bins=numpy.linspace(0.0, 1.0, num=41, endpoint=True)
    )

    grid = tpl.subplot_grid((2, 1))
    grid[0, 0].aprint(method)
    grid[0, 0].aprint("Mesh creation time (sec.):  {:7.2f}".format(elapsed))
    grid[0, 0].aprint(
        "Mesh creation speed (vertices/sec.):  {:7.2f}".format(num_verts / elapsed)
    )
    grid[0, 0].aprint("Number of vertices:      {:7.0f}".format(num_verts))
    grid[0, 0].aprint("Number of cells:         {:7.0f}".format(num_cells))

    grid[1, 0].aprint("min quality: {:5.3f}".format(numpy.min(q)))
    grid[1, 0].aprint("avg quality: {:5.3f}".format(numpy.average(q)))
    grid[1, 0].aprint("max quality: {:5.3f}".format(numpy.max(q)))
    grid[1, 0].hist(q_hist, q_bin_edges, grid=[15, 25])

    grid.show()

    stats = numpy.array(
        [
            elapsed,
            num_verts / elapsed,
            num_verts,
            num_cells,
            numpy.min(q),
            numpy.average(q),
            numpy.max(q),
        ]
    )
    numpy.savetxt(method + "_2d.txt", stats, delimiter=",")


def print_stats_3d(
    angles, quality, method="not specified", elapsed=0.0, num_verts=0.0, num_cells=0.0
):
    """https://github.com/nschloe/optimesh/blob/master/optimesh/helpers.py
    NTS: terminal dimension: 286x93, font size 11, Menlo Regular, Line spacing 0.8, Char. spacing 1.044
    """
    # convert to dihedral angles in radians
    const = 2 * math.sqrt(2) / 3
    angles = numpy.arcsin(const * angles)
    # convert to degrees
    angles = 180 * angles / numpy.pi
    angles_hist, angles_bin_edges = numpy.histogram(
        angles, bins=numpy.linspace(0.0, 180.0, num=73, endpoint=True)
    )

    q = quality
    q_hist, q_bin_edges = numpy.histogram(
        q, bins=numpy.linspace(0.0, 1.0, num=41, endpoint=True)
    )

    grid = tpl.subplot_grid((2, 2))
    grid[0, 0].aprint(method)
    grid[0, 0].aprint("Mesh creation time (sec.):  {:7.2f}".format(elapsed))
    grid[0, 0].aprint(
        "Mesh creation speed (vertices/sec.):  {:7.2f}".format(num_verts / elapsed)
    )
    grid[0, 1].aprint("Number of vertices:      {:7.0f}".format(num_verts))
    grid[0, 1].aprint("Number of cells:         {:7.0f}".format(num_cells))

    grid[1, 0].aprint("min. dihedral angles:     {:7.3f}".format(numpy.min(angles)))
    grid[1, 0].aprint("avg. dihedral angles:     {:7.3f}".format(numpy.mean(angles)))
    grid[1, 0].aprint("max  dihedral angles:     {:7.3f}".format(numpy.max(angles)))
    grid[1, 0].aprint("std dev. dihedral angles: {:7.3f}".format(numpy.std(angles)))
    grid[1, 0].hist(angles_hist, angles_bin_edges, grid=[15, 25])

    grid[1, 1].aprint("min quality: {:5.3f}".format(numpy.min(q)))
    grid[1, 1].aprint("avg quality: {:5.3f}".format(numpy.average(q)))
    grid[1, 1].aprint("max quality: {:5.3f}".format(numpy.max(q)))
    grid[1, 1].hist(q_hist, q_bin_edges, grid=[15, 25])

    grid.show()

    stats = numpy.array(
        [
            elapsed,
            num_verts / elapsed,
            num_verts,
            num_cells,
            numpy.min(angles),
            numpy.mean(angles),
            numpy.max(angles),
            numpy.std(angles),
            numpy.min(q),
            numpy.average(q),
            numpy.max(q),
        ]
    )
    numpy.savetxt(method + ".txt", stats, delimiter=",")
