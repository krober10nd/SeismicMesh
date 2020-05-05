import numpy as np
from mpi4py import MPI

from .cpp import cpputils

from .. import geometry

"""
Migration routines for moving points during parallel Delaunay
"""


def aggregate(points, faces, comm, size, rank):
    """
    Collect global triangulation onto rank 0
    """
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
            tmp = np.reshape(comm.recv(source=r, tag=12), (off_p[r], 2))
            tmp2 = (
                np.reshape(comm.recv(source=r, tag=13), (off_t[r], 3)) + csum_p[r - 1]
            )
            gpoints = np.append(gpoints, tmp, axis=0)
            gfaces = np.append(gfaces, tmp2, axis=0)
    if rank == 0:
        upoints, ufaces, ix = geometry.fixmesh(gpoints, gfaces, delunused=True)
        return upoints, ufaces
    else:
        return True, True


def enqueue(extents, points, faces, rank, size):
    """
    Return ranks that cell sites (vertices of triangulation) need to be sent
    """
    # determine llc (lower left corner) and urc (upper right corner)
    if rank == 0:
        le = [extents[rank + 1][0:2]]
        re = [extents[rank + 1][2:4]]
    elif rank == size - 1:
        le = [extents[rank - 1][0:2]]
        re = [extents[rank - 1][2:4]]
    else:
        le = [extents[rank - 1][0:2], extents[rank + 1][0:2]]
        re = [extents[rank - 1][2:4], extents[rank + 1][2:4]]

    # add dummy box above if rank==0 or under if rank=size-1
    if rank == size - 1:
        le = np.append(le, [-999999, -999999])
        re = np.append(re, [-999998, -999998])

    if rank == 0:
        le = np.insert(le, 0, [-999999, -999999])
        re = np.insert(re, 0, [-999998, -999998])

    vtoe, ptr = geometry.vertex_to_elements(points, faces)
    exports = cpputils.where_to2(points, faces, vtoe, ptr, le, re, rank)

    return exports


def enqueue3(extents, points, faces, rank, size):
    """
    Return ranks that cell sites (vertices of triangulation) need to be sent
    """

    # determine llc (lower left corner) and urc (upper right corner)
    if rank == 0:
        le = [extents[rank + 1][0:3]]
        re = [extents[rank + 1][3:6]]
    elif rank == size - 1:
        le = [extents[rank - 1][0:3]]
        re = [extents[rank - 1][3:6]]
    else:
        le = [extents[rank - 1][0:3], extents[rank + 1][0:3]]
        re = [extents[rank - 1][3:6], extents[rank + 1][3:6]]

    # add dummy box above if rank==0 or under if rank=size-1
    if rank == size - 1:
        le = np.append(le, [-999999, -999999, -999999])
        re = np.append(re, [-999998, -999998, -999998])

    if rank == 0:
        le = np.insert(le, 0, [-999999, -999999, -999999])
        re = np.insert(re, 0, [-999998, -999998, -999998])

    vtoe, ptr = geometry.vertex_to_elements(points, faces, dim=3)
    exports = cpputils.where_to3(points, faces, vtoe, ptr, le, re, rank)

    return exports


def exchange(comm, rank, size, exports):
    """
    Exchange data via MPI using P2P comm
    """
    NSB = int(exports[0, 0])
    NSA = int(exports[0, 1])

    tmp = []
    # send points below
    if NSB != 0:
        comm.send(exports[1 : NSB + 1, 0:2], dest=rank - 1, tag=11)

    # send points above
    if rank != size - 1:
        tmp = np.append(tmp, comm.recv(source=rank + 1, tag=11))

    # receive points from above
    if NSA != 0:
        # all put the top receive from above
        comm.send(exports[NSB + 1 : NSB + 1 + NSA, 0:2], dest=rank + 1, tag=11)

    # receive points from below
    if rank != 0:
        # all but the bottom receive from below
        tmp = np.append(tmp, comm.recv(source=rank - 1, tag=11))

    new_points = np.reshape(tmp, (int(len(tmp) / 2), 2))

    return new_points
