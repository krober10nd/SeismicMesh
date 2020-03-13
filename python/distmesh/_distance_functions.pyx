# encoding: utf-8
"""Distance functions implemented in C."""

#-----------------------------------------------------------------------------
#  Copyright (C) 2004-2012 Per-Olof Persson
#  Copyright (C) 2012 Bradley Froehle

#  Distributed under the terms of the GNU General Public License. You should
#  have received a copy of the license along with this program. If not,
#  see <http://www.gnu.org/licenses/>.
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Cython imports
#-----------------------------------------------------------------------------

cimport numpy as np
np.import_array()

cdef extern from "src/distance_functions.c":
    double _dellipse "dellipse" (double x0, double y0, double a, double b)
    double _dellipsoid "dellipsoid" (double x0, double y0, double z0,
            double a, double b, double z)
    double _dsegment "dsegment" (double x0, double y0, double p1x, double p1y,
            double p2x, double p2y)

#-----------------------------------------------------------------------------
# Functions
#-----------------------------------------------------------------------------

def dellipse(p, axes):
    """
    d = dellipse(p, axes)

    Parameters
    ----------
    p : array, shape (np, 2)
        points
    axes : array, shape (2,)

    Returns
    -------
    d = array, shape (np, )
        distance from each point to the ellipse
    """
    cdef double a, b
    cdef np.npy_intp n, i

    a, b = axes

    cdef np.ndarray[np.double_t, ndim=2] P = \
         np.PyArray_FromObject(p, np.NPY_DOUBLE, 2, 2)
    n = P.shape[0]
    assert P.shape[1] == 2, "array should have shape (np, 2)"

    cdef np.ndarray[np.double_t, ndim=1] D = \
         np.PyArray_SimpleNew(1, &n, np.NPY_DOUBLE)

    for i in range(n):
        D[i] = _dellipse(P[i,0], P[i,1], a, b)
    return D

def dellipsoid(p, axes):
    """
    d = dellipsoid(p, axes)

    Parameters
    ----------
    p : array, shape (np, 3)
        points
    axes : array, shape (3,)

    Returns
    -------
    d = array, shape (np, )
        distance from each point to the ellipsoid
    """
    cdef double a, b, c
    cdef np.npy_intp n, i

    a, b, c = axes

    cdef np.ndarray[np.double_t, ndim=2] P = \
         np.PyArray_FromObject(p, np.NPY_DOUBLE, 2, 2)
    n = P.shape[0]
    assert P.shape[1] == 3, "array should have shape (np, 3)"

    cdef np.ndarray[np.double_t, ndim=1] D = \
         np.PyArray_SimpleNew(1, &n, np.NPY_DOUBLE)

    for i in range(n):
        D[i] = _dellipsoid(P[i,0], P[i,1], P[i,2], a, b, c)
    return D

def dsegment(p, v):
    """
    d = dsegment(p, v)

    Parameters
    ----------
    p : array, shape (np, 2)
        points
    v : array, shape (nv, 2)
        vertices of a closed array, whose edges are v[0]..v[1],
        ... v[nv-2]..v[nv-1]

    Output
    ------
    ds : array, shape (np, nv-1)
        distance from each point to each edge
    """
    cdef double p1x,p1y, p2x,p2y
    cdef np.npy_intp n, nv1, i, iv

    cdef np.ndarray[np.double_t, ndim=2] P = \
         np.PyArray_FromObject(p, np.NPY_DOUBLE, 2, 2)
    n = P.shape[0]
    assert P.shape[1] == 2, "array should have shape (np, 2)"

    cdef np.ndarray[np.double_t, ndim=2] V = \
         np.PyArray_FromObject(v, np.NPY_DOUBLE, 2, 2)
    nv1 = V.shape[0]-1
    assert V.shape[1] == 2, "array should have shape (nv, 2)"

    cdef np.npy_intp DS_dims[2]
    DS_dims[0] = n; DS_dims[1] = nv1

    cdef np.ndarray[np.double_t, ndim=2] DS = \
         np.PyArray_SimpleNew(2, DS_dims, np.NPY_DOUBLE)

    for iv in range(nv1):
        p1x = V[iv,0]
        p1y = V[iv,1]
        p2x = V[iv+1,0]
        p2y = V[iv+1,1]

        for i in range(n):
            DS[i,iv] = _dsegment(P[i,0], P[i,1], p1x, p1y, p2x, p2y)
    return DS
