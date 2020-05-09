import numpy as np
import scipy.sparse as spsparse


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


def unique_rows(a):
    # The cleaner alternative `numpy.unique(a, axis=0)` is slow; cf.
    # <https://github.com/numpy/numpy/issues/11136>.
    b = np.ascontiguousarray(a).view(np.dtype((np.void, a.dtype.itemsize * a.shape[1])))
    a_unique, inv, cts = np.unique(b, return_inverse=True, return_counts=True)
    a_unique = a_unique.view(a.dtype).reshape(-1, a.shape[1])
    return a_unique, inv, cts


def dense(I, J, S, shape=None, dtype=None):
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
        J = np.empty(I.shape, dtype=int)
        J.fill(x)
    if np.isscalar(S):
        x = S
        S = np.empty(I.shape)
        S.fill(x)

    # Turn these into 1-d arrays for processing.
    S = S.flat
    II = I.flat
    J = J.flat
    return spsparse.coo_matrix((S, (II, J)), shape, dtype).toarray()
