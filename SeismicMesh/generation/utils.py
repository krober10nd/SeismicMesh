import copy

import numpy as np
import scipy.sparse as spsparse


def odd(r):
    odds = []
    for i in range(r):
        if i % 2 != 0:
            odds.append(i)
    return odds


def create_staggered_grid(h0, dim, bbox):
    points = np.mgrid[tuple(slice(min, max + h0, h0) for min, max in bbox)].astype(
        float
    )
    odds_rows = odd(points[0].shape[0])
    odds_cols = odd(points[1].shape[0])
    points[1][odds_rows] += h0 / 2
    if dim == 3:
        points[2][odds_cols] += h0 / 2
    points = points.reshape(dim, -1).T
    return points


# def create_staggered_grid_2d(h, bounding_box):
#    """https://github.com/nschloe/dmsh/blob/master/dmsh/main.py"""
#    bounding_box = bounding_box.flatten()
#    x_step = h
#    y_step = h * np.sqrt(3) / 2
#    bb_width = bounding_box[1] - bounding_box[0]
#    bb_height = bounding_box[3] - bounding_box[2]
#    midpoint = [
#        (bounding_box[0] + bounding_box[1]) / 2,
#        (bounding_box[2] + bounding_box[3]) / 2,
#    ]
#
#    num_x_steps = int(bb_width / x_step)
#    if num_x_steps % 2 == 1:
#        num_x_steps -= 1
#    num_y_steps = int(bb_height / y_step)
#    if num_y_steps % 2 == 1:
#        num_y_steps -= 1
#
#    # Generate initial (staggered) point list from bounding box.
#    # Make sure that the midpoint is one point in the grid.
#    x2 = num_x_steps // 2
#    y2 = num_y_steps // 2
#    x, y = np.meshgrid(
#        midpoint[0] + x_step * np.arange(-x2, x2 + 1),
#        midpoint[1] + y_step * np.arange(-y2, y2 + 1),
#    )
#    # Staggered, such that the midpoint is not moved.
#    # Unconditionally move to the right, then add more points to the left.
#    offset = (y2 + 1) % 2
#    x[offset::2] += h / 2
#
#    out = np.column_stack([x.reshape(-1), y.reshape(-1)])
#
#    # add points in the staggered lines to preserve symmetry
#    n = 2 * (-(-y2 // 2))
#    extra = np.empty((n, 2))
#    extra[:, 0] = midpoint[0] - x_step * x2 - h / 2
#    extra[:, 1] = midpoint[1] + y_step * np.arange(-y2 + offset, y2 + 1, 2)
#
#    out = np.concatenate([out, extra])
#    return out


def make_init_points(bbox, rank, size, axis, h0, dim):
    """Create a structured grid in parallel of the entire domain
    Each processor owns a part of the domain.
    """
    _bbox = copy.deepcopy(bbox)

    for i in range(dim):
        if i == axis:
            new_lims = np.linspace(_bbox[i, 0], _bbox[i, 1], size + 1)
            _bbox[i, :] = new_lims[rank : rank + 2]
            if rank != 0:
                # starting point must be lasts + h0
                prev_lims = new_lims[rank - 1 : rank - 1 + 2]
                tmp = np.mgrid[slice(prev_lims[0], prev_lims[1] + h0, h0)]
                _bbox[i, 0] = tmp[-1] + h0

    points = np.mgrid[tuple(slice(min, max + h0, h0) for min, max in _bbox)].astype(
        float
    )
    points = points.reshape(dim, -1).T
    return points


def dense(Ix, J, S, shape=None, dtype=None):
    """
    Similar to MATLAB's SPARSE(I, J, S, ...), but instead returning a
    dense array.
    """

    # Advanced usage: allow J and S to be scalars.
    if np.isscalar(J):
        x = J
        J = np.empty(Ix.shape, dtype=int)
        J.fill(x)
    if np.isscalar(S):
        x = S
        S = np.empty(Ix.shape)
        S.fill(x)

    # Turn these into 1-d arrays for processing.
    S = S.flat
    II = Ix.flat
    J = J.flat
    return spsparse.coo_matrix((S, (II, J)), shape, dtype).toarray()
