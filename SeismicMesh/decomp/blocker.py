import numpy as np


def blocker(points, rank, nblocks, bbox):
    """ Decompose point coordinates into # of blocks
        Blocks are orientated parallel to x-axis and have a neighbor
        above and below +-y the block.
    """
    num_points, dim = points.shape

    assert dim > 2 or dim < 3, "dimensions of points are wrong"
    assert num_points // nblocks > 1, "too few points for chosen nblocks"

    xlims = bbox[0], bbox[2]  # points[:, 0].min(), points[:, 0].max()
    ylims = bbox[1], bbox[3]  # points[:, 1].min(), points[:, 1].max()
    xx, yy = np.meshgrid(
        np.linspace(*xlims, 1, endpoint=False),
        np.linspace(*ylims, nblocks, endpoint=False),
    )

    dx = (xlims[1] - xlims[0]) / 1
    dy = (ylims[1] - ylims[0]) / nblocks

    blocks = []
    block_extents = []
    for low_x, low_y in zip(xx.ravel(), yy.ravel()):
        block = points[
            (points[:, 0] >= low_x)
            & (points[:, 0] <= low_x + dx)
            & (points[:, 1] > low_y)
            & (points[:, 1] < low_y + dy)
        ]
        if block.shape[0]:
            blocks.append(block)

    block_extents = []
    for block in blocks:
        tmpm = np.amin(block, axis=0)
        tmpp = np.amax(block, axis=0)
        # min x min y max x max y
        block_extents.append([bbox[0], tmpm[1], bbox[3], tmpp[1]])

    # correcting block extents so bbox is respected
    block_extents[0][1] = bbox[1]  # min y
    block_extents[-1][3] = bbox[3]  # max y

    return blocks[rank], block_extents
