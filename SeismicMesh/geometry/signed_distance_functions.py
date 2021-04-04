import numpy as np
import matplotlib.pyplot as plt
import itertools

from _fast_geometry import drectangle_fast, dblock_fast

import random


def _generate_samples(bbox, dim, N):
    N = int(N)
    points = []
    _xrange = (bbox[0] - 0.01, bbox[1] + 0.01)
    _yrange = (bbox[2] - 0.01, bbox[3] + 0.01)
    if dim == 2:
        points.append(
            [
                (
                    random.uniform(*_xrange),
                    random.uniform(*_yrange),
                )
                for i in range(N)
            ]
        )
    elif dim == 3:
        _zrange = (bbox[4] - 0.01, bbox[5] + 0.01)
        points.append(
            [
                (
                    random.uniform(*_xrange),
                    random.uniform(*_yrange),
                    random.uniform(*_zrange),
                )
                for i in range(N)
            ]
        )
    points = np.asarray(points)
    points = points.reshape(-1, dim)
    return points


def _show(geo, filename=None, samples=10000):
    p = _generate_samples(geo.bbox, geo.dim, N=samples)
    d = geo.eval(p)
    ix = np.logical_and(d > -0.01, d < 0.01)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    if geo.dim == 2:
        im = ax.scatter(p[ix, 0], p[ix, 1], p[ix, 0] * 0.0, c=d[ix], marker=".", s=5.0)
        ax.set_xlabel("X-axis")
        ax.set_ylabel("Y-axis")
    elif geo.dim == 3:
        im = ax.scatter(p[ix, 0], p[ix, 1], p[ix, 2], c=d[ix], marker=".", s=5.0)
        ax.set_xlabel("X-axis")
        ax.set_ylabel("Y-axis")
        ax.set_zlabel("Z-axis")
    plt.title("Approximate 0-level set")
    fig.colorbar(im, ax=ax)
    im.set_clim(-0.1, 0.1)
    ax.set_aspect("auto")

    if filename is None:
        plt.show()
    else:
        plt.savefig(filename)


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


def _manipulate(object, x):
    if object.translation_vec is not None:
        x = _translate_back(object, x)
    if object.rotation != 0.0:
        x = _rotate_back(object, x)
    if object.v is not None:
        x = _scale_back(object, x)
    return x


def _configure_manipulations(object):
    object = _build_stretch(object)
    object = _build_rotation(object)
    object = _build_translation(object)
    return object


def _scale_back(object, x):
    sd = object.dim
    x = x.T
    x_shape = x.shape
    x = x.reshape(sd, -1)
    vx = np.multiply.outer(np.dot(object.v, x), object.v)
    x = vx / object.alpha + (x.T - vx)
    x = x.T.reshape(x_shape)
    x = x.T
    return x


def _translate_back(object, x):
    return x - object.translation_vec


def _rotate_back(object, x):
    return np.dot(object.R_inv, x.T).T


def _build_stretch(object):
    sd = object.dim
    if object.v is not None:
        assert (
            len(object.v) == object.dim
        ), "Length of stretch vector does not equal dimension"
    if object.v is not None:
        object.alpha = np.sqrt(np.dot(object.v, object.v))
        object.v /= object.alpha
        _corners = corners(object.bbox)
        vx = np.multiply.outer(np.dot(object.v, _corners.T), object.v)
        tmp = (vx * object.alpha + (_corners - vx)).T
        if sd == 2:
            object.bbox = (
                np.min(tmp[0]),
                np.max(tmp[0]),
                np.min(tmp[1]),
                np.max(tmp[1]),
            )
        elif sd == 3:
            object.bbox = (
                np.min(tmp[0]),
                np.max(tmp[0]),
                np.min(tmp[1]),
                np.max(tmp[1]),
                np.min(tmp[2]),
                np.max(tmp[2]),
            )
        object.corners = corners(object.bbox)
    return object


def _build_rotation(object):
    if object.dim == 2:
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
    elif object.dim == 3:
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


def _build_translation(object):
    v = object.translation_vec
    if v is not None:
        assert (
            len(v) == object.dim
        ), "Length of translation vector does not equal dimension"
    if v is not None:
        if object.dim == 2:
            object.bbox = (
                object.bbox[0] + v[0],
                object.bbox[1] + v[0],
                object.bbox[2] + v[1],
                object.bbox[3] + v[1],
            )
        elif object.dim == 3:
            object.bbox = (
                object.bbox[0] + v[0],
                object.bbox[1] + v[0],
                object.bbox[2] + v[1],
                object.bbox[3] + v[1],
                object.bbox[4] + v[2],
                object.bbox[5] + v[2],
            )
    return object


class Repeat:
    def __init__(self, bbox, domain, period):
        self.bbox = bbox
        self.domain = domain
        self.corners = None
        self.period = np.array(period)
        self.parent = Cube(bbox)
        self.dim = 3

    def eval(self, x):
        q = np.mod(x + 0.5 * self.period, self.period) - 0.5 * self.period
        return np.maximum(self.domain.eval(q), self.parent.eval(x))

    def show(self, filename=None, samples=10000):
        _show(self, filename=None, samples=samples)


def _loop_call(func, d):
    tmp = d[0]
    for i in range(0, len(d) - 1):
        tmp = func(tmp, d[i + 1])
    return tmp


class Union:
    def __init__(self, domains, smoothness=0.0):
        geom_dim = [d.dim for d in domains]
        assert np.all(geom_dim != 2) or np.all(geom_dim != 3)
        self.dim = geom_dim[0]
        self.k = smoothness
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

    def _smooth_union(self, d1, d2):
        h = np.maximum(self.k - np.abs(d1 - d2), 0.0)
        return np.minimum(d1, d2) - np.divide(h * h * 0.25, self.k)

    def eval(self, x):
        d = [d.eval(x) for d in self.domains]
        if self.k == 0.0:
            return np.minimum.reduce(d)
        else:
            return _loop_call(self._smooth_union, d)

    def show(self, filename=None, samples=10000):
        _show(self, filename=None, samples=samples)


class Intersection:
    def __init__(self, domains, smoothness=0.0):
        geom_dim = [d.dim for d in domains]
        assert np.all(geom_dim != 2) or np.all(geom_dim != 3)
        self.dim = geom_dim[0]
        self.k = smoothness
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

    def _smooth_intersection(self, d1, d2):
        h = np.maximum(self.k - np.abs(d1 - d2), 0.0)
        return np.maximum(d1, d2) + h * h * 0.25 / self.k

    def eval(self, x):
        d = [d.eval(x) for d in self.domains]
        if self.k == 0.0:
            return np.maximum.reduce(d)
        else:
            return _loop_call(self._smooth_intersection, d)

    def show(self, filename=None, samples=10000):
        _show(self, filename=None, samples=samples)


class Difference:
    def __init__(self, domains, smoothness=0.0):
        geom_dim = [d.dim for d in domains]
        assert np.all(geom_dim != 2) or np.all(geom_dim != 3)
        self.dim = geom_dim[0]
        self.k = smoothness
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

    def _smooth_difference(self, d1, d2):
        h = np.maximum(self.k - np.abs(-d1 - d2), 0.0)
        return np.maximum(-d1, d2) + np.divide(h * h * 0.25, self.k)

    def eval(self, x):
        if self.k == 0.0:
            return np.maximum.reduce(
                [-d.eval(x) if n > 0 else d.eval(x) for n, d in enumerate(self.domains)]
            )
        else:
            d = [d.eval(x) for d in self.domains]
            return _loop_call(self._smooth_difference, d[::-1])

    def show(self, filename=None, samples=10000):
        _show(self, filename=None, samples=samples)


class Disk:
    def __init__(self, x0, r, rotate=0.0, stretch=None, translate=None):
        self.dim = 2
        self.xc = x0[0]
        self.yc = x0[1]
        self.r = r
        self.bbox = (x0[0] - r, x0[0] + r, x0[1] - r, x0[1] + r)
        self.rotation = rotate
        self.v = stretch
        self.translation_vec = translate
        self = _configure_manipulations(self)
        self.corners = None

    def eval(self, x):
        x = _manipulate(self, x)
        return _ddisk(x, self.xc, self.yc, self.r)

    def show(self, filename=None, samples=10000):
        _show(self, filename=None, samples=samples)


class Ball:
    def __init__(self, x0, r, rotate=0.0, stretch=None, translate=None):
        if stretch is not None:
            assert len(stretch) == 3
        self.dim = 3
        self.xc = x0[0]
        self.yc = x0[1]
        self.zc = x0[2]
        self.r = r
        self.bbox = (x0[0] - r, x0[0] + r, x0[1] - r, x0[1] + r, x0[2] - r, x0[2] + r)
        self.rotation = rotate
        self.v = stretch
        self.translation_vec = translate
        self = _configure_manipulations(self)
        self.corners = None

    def eval(self, x):
        x = _manipulate(self, x)
        return dball(x, self.xc, self.yc, self.zc, self.r)

    def show(self, filename=None, samples=10000):
        _show(self, filename=None, samples=samples)


class Rectangle:
    def __init__(self, bbox, rotate=0.0, stretch=None, translate=None):
        self.dim = 2
        self.corners = corners(bbox)
        self.bbox0 = bbox
        self.bbox = bbox
        self.rotation = rotate
        self.v = stretch
        self.translation_vec = translate
        self = _configure_manipulations(self)

    def eval(self, x):
        x = _manipulate(self, x)
        return drectangle_fast(
            x, self.bbox0[0], self.bbox0[1], self.bbox0[2], self.bbox0[3]
        )

    def show(self, filename=None, samples=10000):
        _show(self, filename=None, samples=samples)


class Cube:
    def __init__(self, bbox, rotate=0.0, stretch=None, translate=None):
        self.dim = 3
        self.corners = corners(bbox)
        self.bbox0 = bbox
        self.bbox = bbox
        self.rotation = rotate
        self.v = stretch
        self.translation_vec = translate
        self = _configure_manipulations(self)

    def eval(self, x):
        x = _manipulate(self, x)
        return dblock_fast(
            x,
            self.bbox0[0],
            self.bbox0[1],
            self.bbox0[2],
            self.bbox0[3],
            self.bbox0[4],
            self.bbox0[5],
        )

    def show(self, filename=None, samples=10000):
        _show(self, filename=None, samples=samples)


class Torus:
    def __init__(self, r1, r2, rotate=0.0, stretch=None, translate=None):
        """A torus with outer radius `r1` and inner radius of `r2`"""
        assert r1 > 0.0 and r2 > 0.0
        z = 2 * max(r1, r2)
        self.dim = 3
        self.bbox = (-2 * z, 2 * z, -2 * z, 2 * z, -2 * z, 2 * z)
        self.t = (r1, r2)
        self.rotation = rotate
        self.v = stretch
        self.translation_vec = translate
        self = _configure_manipulations(self)
        self.corners = None

    def eval(self, x):
        x = _manipulate(self, x)
        xz = np.column_stack((x[:, 0], x[:, 2]))
        q = np.column_stack((_length(xz) - self.t[0], x[:, 1]))
        return _length(q) - self.t[1]

    def show(self, filename=None, samples=10000):
        _show(self, filename=None, samples=samples)


class Prism:
    def __init__(self, b, h, rotate=0.0, stretch=None, translate=None):
        self.bbox = (-b, +b, -b, +b, -h, +h)
        self.h = (b, h)
        self.dim = 3
        self.rotation = rotate
        self.v = stretch
        self.translation_vec = translate
        self = _configure_manipulations(self)
        self.corners = None

    def eval(self, x):
        x = _manipulate(self, x)
        q = np.abs(x)
        return np.maximum(
            q[:, 2] - self.h[1],
            np.maximum(q[:, 0] * 0.866025 + x[:, 1] * 0.5, -x[:, 1]) - self.h[0] * 0.5,
        )

    def show(self, filename=None, samples=10000):
        _show(self, filename=None, samples=samples)


class Cylinder:
    def __init__(self, h=1.0, r=0.5, rotate=0.0, stretch=None, translate=None):
        assert h > 0.0 and r > 0.0
        h /= 2.0
        sz = max(h, r)
        self.dim = 3
        self.bbox = (-2 * sz, 2 * sz, -2 * sz, 2 * sz, -2 * sz, 2 * sz)
        self.h = (r, h)
        self.rotation = rotate
        self.v = stretch
        self.translation_vec = translate
        self = _configure_manipulations(self)
        self.corners = None

    def eval(self, x):
        x = _manipulate(self, x)
        xz = np.column_stack((x[:, 0], x[:, 2]))
        lxz = np.column_stack((_length(xz), x[:, 1]))
        d = np.abs(lxz) - self.h
        return np.minimum(np.maximum(d[:, 0], d[:, 1]), 0.0) + _length(
            np.maximum(d, 0.0)
        )

    def show(self, filename=None, samples=10000):
        _show(self, filename=None, samples=samples)


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
