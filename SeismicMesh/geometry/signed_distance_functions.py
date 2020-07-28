import numpy as np


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


def dintersect(d1, d2):
    """Signed distance to set intersection of two regions described by signed
    distance functions d1 and d2.
    Not exact the true signed distance function for the difference,
    for example around corners.
    """
    return np.maximum(d1, d2)


def ddiff(d1, d2):
    """Signed distance to set difference between two regions described by
    signed distance functions d1 and d2.
    Not exact the true signed distance function for the difference,
    for example around corners.
    """
    return np.maximum(d1, -d2)
