import numpy as np
import itertools

from .cpp.fast_geometry import drectangle_fast, dblock_fast


def _length(x):
    return np.sum(np.abs(x) ** 2, axis=-1) ** (1.0 / 2)


def corners(bbox):
    """Get the corners of a box in N-dim"""
    mins = bbox[::2]
    maxs = bbox[1::2]
    return np.array(list(itertools.product(*zip(mins, maxs))))


def _gather_corners(domains):
    corners = [d.corners for d in domains if d.corners is not None]
    if len(corners) == 0:
        return None
    else:
        return np.concatenate(corners)


def _build_rotation2(object):
    object.R = np.array(
        [
            [+np.cos(object.rotation), -np.sin(object.rotation)],
            [+np.sin(object.rotation), +np.cos(object.rotation)],
        ]
    )
    object.R_inv = object.R.T
    if object.rotation != 0.0:
        object.corners = corners(object.bbox)
        tmp = np.dot(object.R, object.corners.T).T
        object.bbox = (
            np.min(tmp[:, 0]),
            np.max(tmp[:, 0]),
            np.min(tmp[:, 1]),
            np.max(tmp[:, 1]),
        )
        object.corners = corners(object.bbox)
    return object


def _build_rotation3(object):
    object.R = np.array(
        [
            [+1, +0, +0],
            [+0, +np.cos(object.rotation), -np.sin(object.rotation)],
            [+0, +np.sin(object.rotation), +np.cos(object.rotation)],
        ]
    )
    object.R_inv = object.R.T
    if object.rotation != 0.0:
        object.corners = corners(object.bbox)
        tmp = np.dot(object.R, object.corners.T).T
        object.bbox = (
            np.min(tmp[:, 0]),
            np.max(tmp[:, 0]),
            np.min(tmp[:, 1]),
            np.max(tmp[:, 1]),
            np.min(tmp[:, 2]),
            np.max(tmp[:, 2]),
        )
        object.corners = corners(object.bbox)
    return object


class Repeat:
    def __init__(self, bbox, domain, period):
        self.bbox = bbox
        self.domain = domain
        self.corners = None
        self.period = np.array(period)
        self.parent = Cube(bbox)

    def eval(self, x):
        q = np.mod(x + 0.5 * self.period, self.period) - 0.5 * self.period
        return np.maximum(self.domain.eval(q), self.parent.eval(x))


class Union:
    def __init__(self, domains):
        geom_dim = [d.dim for d in domains]
        assert np.all(geom_dim != 2) or np.all(geom_dim != 3)
        self.dim = geom_dim[0]
        if self.dim == 2:
            self.bbox = (
                min(d.bbox[0] for d in domains),
                max(d.bbox[1] for d in domains),
                min(d.bbox[2] for d in domains),
                max(d.bbox[3] for d in domains),
            )
        elif self.dim == 3:
            self.bbox = (
                min(d.bbox[0] for d in domains),
                max(d.bbox[1] for d in domains),
                min(d.bbox[2] for d in domains),
                max(d.bbox[3] for d in domains),
                min(d.bbox[4] for d in domains),
                max(d.bbox[5] for d in domains),
            )
        self.corners = _gather_corners(domains)
        self.domains = domains

    def eval(self, x):
        return np.minimum.reduce([d.eval(x) for d in self.domains])


class Intersection:
    def __init__(self, domains):
        geom_dim = [d.dim for d in domains]
        assert np.all(geom_dim != 2) or np.all(geom_dim != 3)
        self.dim = geom_dim[0]
        if self.dim == 2:
            self.bbox = (
                min(d.bbox[0] for d in domains),
                max(d.bbox[1] for d in domains),
                min(d.bbox[2] for d in domains),
                max(d.bbox[3] for d in domains),
            )
        elif self.dim == 3:
            self.bbox = (
                min(d.bbox[0] for d in domains),
                max(d.bbox[1] for d in domains),
                min(d.bbox[2] for d in domains),
                max(d.bbox[3] for d in domains),
                min(d.bbox[4] for d in domains),
                max(d.bbox[5] for d in domains),
            )
        self.corners = _gather_corners(domains)
        self.domains = domains

    def eval(self, x):
        return np.maximum.reduce([d.eval(x) for d in self.domains])


class Difference:
    def __init__(self, domains):
        geom_dim = [d.dim for d in domains]
        assert np.all(geom_dim != 2) or np.all(geom_dim != 3)
        self.dim = geom_dim[0]
        if self.dim == 2:
            self.bbox = (
                min(d.bbox[0] for d in domains),
                max(d.bbox[1] for d in domains),
                min(d.bbox[2] for d in domains),
                max(d.bbox[3] for d in domains),
            )
        elif self.dim == 3:
            self.bbox = (
                min(d.bbox[0] for d in domains),
                max(d.bbox[1] for d in domains),
                min(d.bbox[2] for d in domains),
                max(d.bbox[3] for d in domains),
                min(d.bbox[4] for d in domains),
                max(d.bbox[5] for d in domains),
            )
        self.corners = _gather_corners(domains)
        self.domains = domains

    def eval(self, x):
        return np.maximum.reduce(
            [-d.eval(x) if n > 0 else d.eval(x) for n, d in enumerate(self.domains)]
        )


class Disk:
    def __init__(self, x0, r, rotate=0.0):
        self.dim = 2
        self.xc = x0[0]
        self.yc = x0[1]
        self.r = r
        self.bbox = (x0[0] - r, x0[0] + r, x0[1] - r, x0[1] + r)
        self.rotation = rotate
        self = _build_rotation2(self)
        self.corners = None

    def eval(self, x):
        if self.rotation != 0.0:
            x = np.dot(self.R_inv, x.T).T
        return _ddisk(x, self.xc, self.yc, self.r)


class Ball:
    def __init__(self, x0, r, rotate=0.0):
        self.dim = 3
        self.xc = x0[0]
        self.yc = x0[1]
        self.zc = x0[2]
        self.r = r
        self.bbox = (x0[0] - r, x0[0] + r, x0[1] - r, x0[1] + r, x0[2] - r, x0[2] + r)
        self.rotation = rotate
        self = _build_rotation3(self)
        self.corners = None

    def eval(self, x):
        if self.rotation != 0.0:
            x = np.dot(self.R_inv, x.T).T
        return dball(x, self.xc, self.yc, self.zc, self.r)


class Rectangle:
    def __init__(self, bbox, rotate=0.0):
        self.dim = 2
        self.corners = corners(bbox)
        self.bbox0 = bbox
        self.bbox = bbox
        self.rotation = rotate
        self = _build_rotation2(self)

    def eval(self, x):
        if self.rotation != 0.0:
            x = np.dot(self.R_inv, x.T).T
        return drectangle_fast(
            x, self.bbox0[0], self.bbox0[1], self.bbox0[2], self.bbox0[3]
        )


class Cube:
    def __init__(self, bbox, rotate=0.0):
        self.dim = 3
        self.corners = corners(bbox)
        self.bbox0 = bbox
        self.bbox = bbox
        self.rotation = rotate
        self = _build_rotation3(self)

    def eval(self, x):
        if self.rotation != 0.0:
            x = np.dot(self.R_inv, x.T).T
        return dblock_fast(
            x,
            self.bbox0[0],
            self.bbox0[1],
            self.bbox0[2],
            self.bbox0[3],
            self.bbox0[4],
            self.bbox0[5],
        )


class Torus:
    def __init__(self, r1, r2, rotate=0.0):
        """A torus with outer radius `r1` and inner radius of `r2`"""
        assert r1 > 0.0 and r2 > 0.0
        z = 2 * max(r1, r2)
        self.dim = 3
        self.bbox = (-2 * z, 2 * z, -2 * z, 2 * z, -2 * z, 2 * z)
        self.t = (r1, r2)
        self.rotation = rotate
        self = _build_rotation3(self)
        self.corners = None

    def eval(self, x):
        if self.rotation != 0.0:
            x = np.dot(self.R_inv, x.T).T
        xz = np.column_stack((x[:, 0], x[:, 2]))
        q = np.column_stack((_length(xz) - self.t[0], x[:, 1]))
        return _length(q) - self.t[1]


class Prism:
    def __init__(self, b, h, rotate=0.0):
        self.bbox = (-b, +b, -b, +b, -h, +h)
        self.h = (b, h)
        self.dim = 3
        self.rotation = rotate
        self = _build_rotation3(self)
        self.corners = None

    def eval(self, x):
        if self.rotation != 0.0:
            x = np.dot(self.R_inv, x.T).T
        q = np.abs(x)
        return np.maximum(
            q[:, 2] - self.h[1],
            np.maximum(q[:, 0] * 0.866025 + x[:, 1] * 0.5, -x[:, 1]) - self.h[0] * 0.5,
        )


class Cylinder:
    def __init__(self, h=1.0, r=0.5, rotate=0.0):
        assert h > 0.0 and r > 0.0
        h /= 2.0
        sz = max(h, r)
        self.dim = 3
        self.bbox = (-2 * sz, 2 * sz, -2 * sz, 2 * sz, -2 * sz, 2 * sz)
        self.h = (r, h)
        self.rotation = rotate
        self = _build_rotation3(self)
        self.corners = None

    def eval(self, x):
        if self.rotation != 0.0:
            x = np.dot(self.R_inv, x.T).T
        xz = np.column_stack((x[:, 0], x[:, 2]))
        lxz = np.column_stack((_length(xz), x[:, 1]))
        d = np.abs(lxz) - self.h
        return np.minimum(np.maximum(d[:, 0], d[:, 1]), 0.0) + _length(
            np.maximum(d, 0.0)
        )


def _ddisk(p, xc, yc, r):
    """Signed distance to disk centered at xc, yc with radius r."""
    return np.sqrt(((p - np.array([xc, yc])) ** 2).sum(-1)) - r


def dball(p, xc, yc, zc, r):
    """Signed distance function for a ball centered at xc,yc,zc with radius  r."""
    return np.sqrt((p[:, 0] - xc) ** 2 + (p[:, 1] - yc) ** 2 + (p[:, 2] - zc) ** 2) - r


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
