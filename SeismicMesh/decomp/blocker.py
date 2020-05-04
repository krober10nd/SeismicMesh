import numpy as np


def blocker(points, rank, nblocks, axis=0):
    """ Decompose point coordinates into # of blocks
        Blocks are orientated parallel to axis and have two neighbors
        above and below/left and right to the block depending on axis.
    """
    num_points, dim = points.shape
    # EPS = np.finfo(float).eps

    assert dim > 2 or dim < 3, "dimensions of points are wrong"
    assert num_points // nblocks > 1, "too few points for chosen nblocks"

    x_sorted = np.argsort(points[:, 0])
    y_sorted = np.argsort(points[:, 1])

    num_points = points.shape[0]

    step = num_points // nblocks
    blocks = []
    if axis == 0:
        ixx, iyy = np.meshgrid(
            np.arange(num_points, step=step), np.arange(num_points, step=num_points)
        )
        if (num_points % nblocks) != 0:
            ixx[0][-1] = num_points
        for idx, idy in zip(ixx.ravel(), iyy.ravel()):
            common = set(x_sorted[idx : idx + step]).intersection(y_sorted)
            blocks.append(points.take(list(common), axis=0))
    elif axis == 1:
        ixx, iyy = np.meshgrid(
            np.arange(num_points, step=num_points), np.arange(num_points, step=step)
        )
        if (num_points % nblocks) != 0:
            iyy[-1][0] = num_points
        for idx, idy in zip(ixx.ravel(), iyy.ravel()):
            common = set(x_sorted).intersection(y_sorted[idy : idy + step])
            blocks.append(points.take(list(common), axis=0))

    blocks = [x for x in blocks if len(x) != 0]

    # xlims = points[:, 0].min() - EPS, points[:, 0].max() + EPS
    # ylims = points[:, 1].min() - EPS, points[:, 1].max() + EPS
    # if axis == 0:
    #    xx, yy = np.meshgrid(
    #        np.linspace(*xlims, 1, endpoint=False),
    #        np.linspace(*ylims, nblocks, endpoint=False),
    #    )

    #    dx = (xlims[1] - xlims[0]) / 1
    #    dy = (ylims[1] - ylims[0]) / nblocks
    # elif axis == 1:
    #    xx, yy = np.meshgrid(
    #        np.linspace(*xlims, nblocks, endpoint=False),
    #        np.linspace(*ylims, 1, endpoint=False),
    #    )

    #    dx = (xlims[1] - xlims[0]) / nblocks
    #    dy = (ylims[1] - ylims[0]) / 1

    # blocks = []
    # block_extents = []
    # for low_x, low_y in zip(xx.ravel(), yy.ravel()):
    #    block = points[
    #        (points[:, 0] >= low_x)
    #        & (points[:, 0] <= low_x + dx)
    #        & (points[:, 1] >= low_y)
    #        & (points[:, 1] <= low_y + dy)
    #    ]
    #    if block.shape[0]:
    #        blocks.append(block)

    block_extents = []
    for block in blocks:
        tmpm = np.amin(block, axis=0)
        tmpp = np.amax(block, axis=0)
        # min x min y max x max y
        block_extents.append([tmpm[0], tmpm[1], tmpp[0], tmpp[1]])

    return blocks[rank], block_extents
