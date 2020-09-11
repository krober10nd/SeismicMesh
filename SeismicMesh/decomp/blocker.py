import numpy as np


def blocker(points, rank, num_blocks, axis=0):
    """Decompose point coordinates into # of blocks
    Blocks are orientated parallel to axis and have two neighbors
    above and below/left and right to the block depending on axis.
    """
    num_points, dim = points.shape
    EPS = np.finfo(float).eps

    if dim < 2 or dim > 3:
        raise ValueError("Dimensions of points are not supported")
    assert num_points // num_blocks > 1, "too few points for chosen num_blocks"

    x_lims = points[:, 0].min() - EPS, points[:, 0].max() + EPS
    y_lims = points[:, 1].min() - EPS, points[:, 1].max() + EPS
    if dim > 2:
        z_lims = points[:, 2].min() - EPS, points[:, 2].max() + EPS

    if axis == 0:
        if dim == 2:
            xx, yy = np.meshgrid(
                np.linspace(*x_lims, 1, endpoint=False),
                np.linspace(*y_lims, num_blocks, endpoint=False),
            )

            dx = (x_lims[1] - x_lims[0]) / 1
            dy = (y_lims[1] - y_lims[0]) / num_blocks
        elif dim == 3:
            xx, yy, zz = np.meshgrid(
                np.linspace(*x_lims, 1, endpoint=False),
                np.linspace(*y_lims, num_blocks, endpoint=False),
                np.linspace(*z_lims, 1, endpoint=False),
            )

            dx = (x_lims[1] - x_lims[0]) / 1
            dy = (y_lims[1] - y_lims[0]) / num_blocks
            dz = (z_lims[1] - z_lims[0]) / 1

    elif axis == 1:
        if dim == 2:
            xx, yy = np.meshgrid(
                np.linspace(*x_lims, num_blocks, endpoint=False),
                np.linspace(*y_lims, 1, endpoint=False),
            )

            dx = (x_lims[1] - x_lims[0]) / num_blocks
            dy = (y_lims[1] - y_lims[0]) / 1
        elif dim == 3:
            xx, yy, zz = np.meshgrid(
                np.linspace(*x_lims, num_blocks, endpoint=False),
                np.linspace(*y_lims, 1, endpoint=False),
                np.linspace(*z_lims, 1, endpoint=False),
            )

            dx = (x_lims[1] - x_lims[0]) / num_blocks
            dy = (y_lims[1] - y_lims[0]) / 1
            dz = (z_lims[1] - z_lims[0]) / 1
    elif axis == 2:
        xx, yy, zz = np.meshgrid(
            np.linspace(*x_lims, 1, endpoint=False),
            np.linspace(*y_lims, num_blocks, endpoint=False),
            np.linspace(*z_lims, num_blocks, endpoint=False),
        )

        dx = (x_lims[1] - x_lims[0]) / 1
        dy = (y_lims[1] - y_lims[0]) / 1
        dz = (z_lims[1] - z_lims[0]) / num_blocks

    blocks = []
    block_extents = []

    if dim == 2:
        for low_x, low_y in zip(xx.ravel(), yy.ravel()):
            block = points[
                (points[:, 0] >= low_x)
                & (points[:, 0] <= low_x + dx)
                & (points[:, 1] >= low_y)
                & (points[:, 1] <= low_y + dy)
            ]
            if block.shape[0]:
                blocks.append(block)

        block_extents = []
        for block in blocks:
            tmpm = np.amin(block, axis=0)
            tmpp = np.amax(block, axis=0)
            # min x min y max x max y
            block_extents.append([tmpm[0], tmpm[1], tmpp[0], tmpp[1]])
    elif dim == 3:
        for low_x, low_y, low_z in zip(xx.ravel(), yy.ravel(), zz.ravel()):
            block = points[
                (points[:, 0] >= low_x)
                & (points[:, 0] <= low_x + dx)
                & (points[:, 1] >= low_y)
                & (points[:, 1] <= low_y + dy)
                & (points[:, 2] >= low_z)
                & (points[:, 2] <= low_z + dz)
            ]
            if block.shape[0]:
                blocks.append(block)

        block_extents = []
        for block in blocks:
            tmpm = np.amin(block, axis=0)
            tmpp = np.amax(block, axis=0)
            # min x min y min z max x max y max z
            block_extents.append([tmpm[0], tmpm[1], tmpm[2], tmpp[0], tmpp[1], tmpp[2]])

    return blocks, block_extents
