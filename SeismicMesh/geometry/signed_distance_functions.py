import numpy as np


class Circle:
    def __init__(self, xc, yc, r):
        self.xc = xc
        self.yc = yc
        self.r = r

    def eval(self, x):
        return _dcircle(x, self.xc, self.yc, self.r)


class Rectangle:
    def __init__(self, bbox):
        self.x1 = bbox[0]
        self.x2 = bbox[1]
        self.y1 = bbox[2]
        self.y2 = bbox[3]

    def eval(self, x):
        return drectangle(x, self.x1, self.x2, self.y1, self.y2)


class Cube:
    def __init__(self, bbox):
        self.x1 = bbox[0]
        self.x2 = bbox[1]
        self.y1 = bbox[2]
        self.y2 = bbox[3]
        self.z1 = bbox[4]
        self.z2 = bbox[5]

    def eval(self, x):
        return dblock(x, self.x1, self.x2, self.y1, self.y2, self.z1, self.z2)


def _dcircle(p, xc, yc, r):
    """Signed distance to circle centered at xc, yc with radius r."""
    return np.sqrt(((p - np.array([xc, yc])) ** 2).sum(-1)) - r


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
