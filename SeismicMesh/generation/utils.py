import copy

import numpy as np
import scipy.sparse as spsparse


def closestNumber(n, m):
    # Find the quotient
    q = int(n / m)

    # 1st possible closest number
    n1 = m * q

    # 2nd possible closest number
    if (n * m) > 0:
        n2 = m * (q + 1)
    else:
        n2 = m * (q - 1)

    # if true, then n1 is the required closest number
    if abs(n - n1) < abs(n - n2):
        return n1

    # else n2 is the required closest number
    return n2


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


def __setdiff_rows(A, B, return_index=False):
    """
    Similar to MATLAB's setdiff(A, B, 'rows'), this returns C, I
    where C are the row of A that are not in B and I satisfies
    C = A[I,:].

    Returns I if return_index is True.
    """
    A = np.require(A, requirements="C")
    B = np.require(B, requirements="C")

    assert A.ndim == 2, "array must be 2-dim'l"
    assert B.ndim == 2, "array must be 2-dim'l"
    assert A.shape[1] == B.shape[1], "arrays must have the same number of columns"
    assert A.dtype == B.dtype, "arrays must have the same data type"

    # NumPy provides setdiff1d, which operates only on one dimensional
    # arrays. To make the array one-dimensional, we interpret each row
    # as being a string of characters of the appropriate length.
    orig_dtype = A.dtype
    ncolumns = A.shape[1]
    dtype = np.dtype((np.character, orig_dtype.itemsize * ncolumns))
    C = (
        np.setdiff1d(A.view(dtype), B.view(dtype))
        .view(A.dtype)
        .reshape((-1, ncolumns), order="C")
    )
    if return_index:
        raise NotImplementedError
    else:
        return C


def unique_rows(ar):
    ar_row_view = ar.view("|S%d" % (ar.itemsize * ar.shape[1]))
    _, unique_row_indices = np.unique(ar_row_view, return_index=True)
    ar_out = ar[unique_row_indices]
    return ar_out


def dense(Ix, J, S, shape=None, dtype=None):
    """
    Similar to MATLAB's SPARSE(I, J, S, ...), but instead returning a
    dense array.

    Usage
    -----
    >>> shape = (m, n)
    >>> A = dense(I, J, S, shape, dtype)
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
