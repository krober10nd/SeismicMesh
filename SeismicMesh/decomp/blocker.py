import numpy as np


def blocker(points, rank, nblocks, axis=0):  # noqa: C901
    """ Decompose point coordinates into # of blocks
        Blocks are orientated parallel to axis and have two neighbors
        above and below/left and right to the block depending on axis.
    """
    num_points, dim = points.shape

    assert dim > 2 or dim < 3, "dimensions of points are incorrect"
    assert num_points // nblocks > 1, "too few points for nblocks"
    assert axis > 0 or axis < 2, " axis is incorrect"

    x_sorted = np.argsort(points[:, 0])
    y_sorted = np.argsort(points[:, 1])
    if dim > 2:
        z_sorted = np.argsort(points[:, 2])

    num_points = points.shape[0]

    step = num_points // nblocks
    blocks = []
    if dim == 2:
        # decompose x-axis
        if axis == 0:
            ixx, iyy = np.meshgrid(
                np.arange(num_points, step=step), np.arange(num_points, step=num_points)
            )
            k = 0
            ma = np.amax(ixx.shape)
            for idx, idy in zip(ixx.ravel(), iyy.ravel()):
                # if last entry and odd number of points, add remainder to last block
                # and skip last grid
                if k == (ma - 2):
                    step += num_points % nblocks
                if k == (ma - 1) and (num_points % nblocks) != 0:
                    continue
                k += 1
                common = set(x_sorted[idx : idx + step]).intersection(y_sorted)
                blocks.append(points.take(list(common), axis=0))
        # decompose y-axis
        elif axis == 1:
            ixx, iyy = np.meshgrid(
                np.arange(num_points, step=num_points), np.arange(num_points, step=step)
            )
            k = 0
            ma = np.amax(ixx.shape)
            for idx, idy in zip(ixx.ravel(), iyy.ravel()):
                # if last entry and odd number of points, add remainder to last block
                # and skip last grid
                if k == (ma - 2):
                    step += num_points % nblocks
                if k == (ma - 1) and (num_points % nblocks) != 0:
                    continue
                k += 1

                common = set(y_sorted[idy : idy + step]).intersection(x_sorted)
                blocks.append(points.take(list(common), axis=0))
    elif dim == 3:
        # decompose x axis
        if axis == 0:
            ixx, iyy, izz = np.meshgrid(
                np.arange(num_points, step=step),
                np.arange(num_points, step=num_points),
                np.arange(num_points, step=num_points),
            )
            ma = np.amax(ixx.shape)
            k = 0
            for idx, idy, idz in zip(ixx.ravel(), iyy.ravel(), izz.ravel()):
                # if last entry and odd number of points, add remainder to last block
                # and skip last grid
                if k == (ma - 2):
                    step += num_points % nblocks
                if k == (ma - 1) and (num_points % nblocks) != 0:
                    continue
                k += 1
                common = set(
                    set(x_sorted[idx : idx + step]).intersection(y_sorted)
                ).intersection(z_sorted)
                blocks.append(points.take(list(common), axis=0))
        # decompose y axis
        elif axis == 1:
            ixx, iyy, izz = np.meshgrid(
                np.arange(num_points, step=num_points),
                np.arange(num_points, step=step),
                np.arange(num_points, step=num_points),
            )
            ma = np.amax(ixx.shape)
            k = 0
            for idx, idy, idz in zip(ixx.ravel(), iyy.ravel(), izz.ravel()):
                # if last entry and odd number of points, add remainder to last block
                # and skip last grid
                if k == (ma - 2):
                    step += num_points % nblocks
                if k == (ma - 1) and (num_points % nblocks) != 0:
                    continue
                k += 1
                common = set(
                    set(y_sorted[idy : idy + step]).intersection(x_sorted)
                ).intersection(z_sorted)
                blocks.append(points.take(list(common), axis=0))
        # decompose z axis
        elif axis == 2:
            ixx, iyy, izz = np.meshgrid(
                np.arange(num_points, step=num_points),
                np.arange(num_points, step=num_points),
                np.arange(num_points, step=step),
            )
            ma = np.amax(ixx.shape)
            k = 0
            for idx, idy, idz in zip(ixx.ravel(), iyy.ravel(), izz.ravel()):
                # if last entry and odd number of points, add remainder to last block
                # and skip last grid
                if k == (ma - 2):
                    step += num_points % nblocks
                if k == (ma - 1) and (num_points % nblocks) != 0:
                    continue
                k += 1
                common = set(
                    set(z_sorted[idz : idz + step]).intersection(x_sorted)
                ).intersection(y_sorted)
                blocks.append(points.take(list(common), axis=0))

    # delete zero length partitions
    blocks = [x for x in blocks if len(x) != 0]

    block_extents = []
    for block in blocks:
        tmpm = np.amin(block, axis=0)
        tmpp = np.amax(block, axis=0)
        # min x min y max x max y
        if dim == 2:
            block_extents.append([tmpm[0], tmpm[1], tmpp[0], tmpp[1]])
        # min x min y max x max y min z max z
        elif dim == 3:
            block_extents.append([tmpm[0], tmpm[1], tmpm[2], tmpp[0], tmpp[1], tmpp[2]])

    return blocks[rank], block_extents
