import copy

import numpy as np
from mpi4py import MPI
from scipy.interpolate import RegularGridInterpolator

from .. import geometry
from .cpp import cpputils

"""
Migration routines for moving things during parallel Delaunay
"""


def localize_sizing_function(fh, h0, bbox, dim, axis, comm):
    """
    Localize the global sizing function into local chunks
    that cover the span of the local point set.
    """
    rank = comm.Get_rank()
    size = comm.Get_size()
    for r in range(0, size):
        _bbox = copy.deepcopy(bbox)
        if rank == 0:
            # form local point set
            for i in range(dim):
                if i == axis:
                    new_lims = np.linspace(_bbox[i, 0], _bbox[i, 1], size + 1)
                    _bbox[i, :] = new_lims[r : r + 2]
            grid = np.mgrid[
                tuple(slice(min - 3 * h0, max + 3 * h0, h0) for min, max in _bbox)
            ].astype(float)
            # interpolate global --> local sizing grid
            lh = fh(tuple(grid[d] for d in range(dim)))
            # updated grid vectors
            vecs = [np.arange(min - 3 * h0, max + 3 * h0, h0) for min, max in _bbox]
            # form local interpolant
            _lfh = RegularGridInterpolator(
                vecs, lh, bounds_error=False, fill_value=None
            )
            if r == 0:
                lfh = copy.deepcopy(_lfh)
                continue
            else:
                # send local interpolant to r
                comm.send(_lfh, dest=r, tag=11)
        else:
            if rank == r:
                # recv local interpolant from rank 0
                lfh = comm.recv(source=0, tag=11)
    return lfh


def localize_points(blocks, extents, comm, dim):
    """ Distribute points to local subdomains """
    rank = comm.Get_rank()
    size = comm.Get_size()
    for local in range(1, size):
        if rank == 0:
            comm.send(blocks[local], dest=local, tag=12)
        elif rank == local:
            points = np.reshape(comm.recv(source=0, tag=12), (-1, dim))
        if rank == 0:
            points = blocks[0]
        if rank != 0:
            extents = None
        extents = comm.bcast(extents, 0)
    return points, extents


def aggregate(points, faces, comm, size, rank, dim=2):
    """
    Collect global triangulation onto rank 0
    """
    # geometry.dump_mesh(points, faces, rank)

    soff_p = np.zeros((size), dtype=int)
    soff_t = np.zeros((size), dtype=int)

    soff_p[rank] = len(points)
    soff_t[rank] = len(faces)

    off_p = np.zeros((size), dtype=int)
    off_t = np.zeros((size), dtype=int)

    comm.Reduce(soff_p, off_p, op=MPI.SUM, root=0)
    comm.Reduce(soff_t, off_t, op=MPI.SUM, root=0)

    if rank == 0:
        csum_p = np.cumsum(off_p)
        gpoints = points
        gfaces = faces
    for r in np.arange(1, size):
        if rank == r:
            comm.send(points, dest=0, tag=12)
            comm.send(faces, dest=0, tag=13)
        if rank == 0:
            tmp = np.reshape(comm.recv(source=r, tag=12), (off_p[r], dim))
            tmp2 = (
                np.reshape(comm.recv(source=r, tag=13), (off_t[r], dim + 1))
                + csum_p[r - 1]
            )
            gpoints = np.append(gpoints, tmp, axis=0)
            gfaces = np.append(gfaces, tmp2, axis=0)
    if rank == 0:
        upoints, ufaces, ix = geometry.fixmesh(gpoints, gfaces, delunused=True, dim=dim)
        return upoints, ufaces
    else:
        return True, True


def enqueue(extents, points, faces, rank, size, dim=2):
    """
    Return ranks that cell sites (vertices of triangulation) need to be sent
    """
    # determine llc (lower left corner) and urc (upper right corner)
    if rank == 0:
        le = [extents[rank + 1][0:dim]]
        re = [extents[rank + 1][dim : dim * 2]]
    elif rank == size - 1:
        le = [extents[rank - 1][0:dim]]
        re = [extents[rank - 1][dim : dim * 2]]
    else:
        le = [extents[rank - 1][0:dim], extents[rank + 1][0:dim]]
        re = [extents[rank - 1][dim : dim * 2], extents[rank + 1][dim : dim * 2]]

    # add dummy box above if rank==0 or under if rank=size-1
    if rank == size - 1:
        le = np.append(le, [-999999999] * dim)
        re = np.append(re, [-999999998] * dim)

    if rank == 0:
        le = np.insert(le, 0, [-999999999] * dim)
        re = np.insert(re, 0, [-999999998] * dim)

    vtoe, ptr = geometry.vertex_to_entities(points, faces, dim=dim)

    if dim == 2:
        exports = cpputils.where_to2(points, faces, vtoe, ptr, le, re, rank)
    elif dim == 3:
        exports = cpputils.where_to3(points, faces, vtoe, ptr, le, re, rank)

    return exports


def exchange(comm, rank, size, exports, dim=2):
    """
    Exchange data via MPI using P2P comm
    """
    NSB = int(exports[0, 0])
    NSA = int(exports[0, 1])

    tmp = []
    # send points below
    if NSB != 0 and (rank != 0):  # rank 0 can't send below
        comm.send(exports[1 : NSB + 1, 0:dim], dest=rank - 1, tag=11)

    # recv  points from above
    if rank != size - 1:
        tmp = np.append(tmp, comm.recv(source=rank + 1, tag=11))

    # send points above
    if NSA != 0 and rank != (size - 1):  # topmost rank can't send to above
        # all put the topmost rank receive from above
        comm.send(exports[NSB + 1 : NSB + 1 + NSA, 0:dim], dest=rank + 1, tag=11)

    # receive points from below
    if rank != 0:
        # all but the bottommost rank receive from below
        tmp = np.append(tmp, comm.recv(source=rank - 1, tag=11))

    new_points = np.reshape(tmp, (int(len(tmp) / dim), dim))

    return new_points
