import numpy as np


def dblock_v2(p, x1, x2, y1, y2, z1, z2, r=0.0):
    """Correct 3d sdf for a block, b is half the size of the box in x,y,z with corner rounding r
       note query points p are centered to 0
    """
    max = np.maximum
    min = np.minimum
    b = np.abs(np.array([x2 - x1, y2 - y1, z2 - z1]) / 2.0)
    q = np.abs(p) - b
    return (
        np.sqrt(np.sum(max(q, 0.0) ** 2, axis=1))
        + min(max(q[:, 0], max(q[:, 1], q[:, 2])), 0.0)
        - r
    )


def drectangle(p, x1, x2, y1, y2):
    min = np.minimum
    """Signed distance function for rectangle with corners (x1,y1), (x2,y1),
    (x1,y2), (x2,y2).
    This has an incorrect distance to the four corners but that isn't a big deal
    """
    return -min(min(min(-y1 + p[:, 1], y2 - p[:, 1]), -x1 + p[:, 0]), x2 - p[:, 0])


def dblock(p, x1, x2, y1, y2, z1, z2):
    min = np.minimum
    return -min(
        min(
            min(min(min(-z1 + p[:, 2], z2 - p[:, 2]), -y1 + p[:, 1]), y2 - p[:, 1]),
            -x1 + p[:, 0],
        ),
        x2 - p[:, 0],
    )
