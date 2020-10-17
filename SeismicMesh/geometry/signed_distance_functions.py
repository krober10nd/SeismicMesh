import numpy as np


class Disk:
    def __init__(self, x0, r):
        self.xc = x0[0]
        self.yc = x0[1]
        self.r = r
        self.x1 = x0[0] - r
        self.x2 = x0[0] + r
        self.y1 = x0[1] - r
        self.y2 = x0[1] + r

    def eval(self, x):
        return _ddisk(x, self.xc, self.yc, self.r)


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


def _ddisk(p, xc, yc, r):
    """Signed distance to disk centered at xc, yc with radius r."""
    return np.sqrt(((p - np.array([xc, yc])) ** 2).sum(-1)) - r


def drectangle(p, x1, x2, y1, y2):
    min = np.minimum
    """Signed distance function for rectangle with corners (x1,y1), (x2,y1),
    (x1,y2), (x2,y2).
    This has an incorrect distance to the four corners but that isn't a big deal
    """
    return -min(min(min(-y1 + p[:, 1], y2 - p[:, 1]), -x1 + p[:, 0]), x2 - p[:, 0])


def dblock0(p, x1, x2, y1, y2, z1, z2):
    # adapted from:
    # https://github.com/nschloe/dmsh/blob/3305c417d373d509c78491b24e77409411aa18c2/dmsh/geometry/rectangle.py#L31
    # outside dist
    # https://gamedev.stackexchange.com/a/44496
    w = x2 - x1
    h = y2 - y1
    d = z2 - z1
    cx = (x1 + x2) / 2
    cy = (y1 + y2) / 2
    cz = (z1 + z2) / 2
    dx = np.abs(p[:, 0] - cx) - w / 2
    dy = np.abs(p[:, 1] - cy) - h / 2
    dz = np.abs(p[:, 2] - cz) - d / 2
    is_inside = (dx <= 0) & (dy <= 0) & (dz <= 0)
    dx[dx < 0.0] = 0.0
    dy[dy < 0.0] = 0.0
    dz[dz < 0.0] = 0.0
    dist = np.sqrt(dx ** 2 + dy ** 2 + dz ** 2)
    # inside dist
    a = np.array(
        [
            p[is_inside, 0] - x1,
            x2 - p[is_inside, 0],
            p[is_inside, 1] - y1,
            y2 - p[is_inside, 1],
            p[is_inside, 2] - z1,
            z2 - p[is_inside, 2],
        ]
    )
    dist[is_inside] = -np.min(a, axis=0)
    return dist


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
