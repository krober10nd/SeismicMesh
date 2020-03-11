# encoding: utf-8
"""Utilities for manipulating meshes."""

#-----------------------------------------------------------------------------
#  Copyright (C) 2004-2012 Per-Olof Persson
#  Copyright (C) 2012 Bradley Froehle

#  Distributed under the terms of the GNU General Public License. You should
#  have received a copy of the license along with this program. If not,
#  see <http://www.gnu.org/licenses/>.
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------

from __future__ import division

import itertools
import numpy as np

# Local imports.
import distmesh.mlcompat as ml

# Py3k
try: range = xrange
except: pass

__all__ = [
    'boundedges',
    'boundedgesnd',
    'uniref',
    'bndproj',
    'mkt2t',
    'assert_t2t_t2n',
    'circumcenter',
    'uniformity',
    'fixmesh',
    'simpqual',
    'simpvol',
    ]

#-----------------------------------------------------------------------------
# Functions
#-----------------------------------------------------------------------------

def boundedges(p, t):
    """Find boundary edges of a triangular mesh.

    Usage
    -----
    >>> e = boundedges(p,t)
    """
    edges = np.vstack((t[:,[0,1]],
                       t[:,[0,2]],
                       t[:,[1,2]]))
    node3 = np.hstack((t[:,2],t[:,1],t[:,0]))
    edges.sort(1)
    _, ix, jx = ml.unique_rows(edges, True, True)
    vec, _ = np.histogram(jx, np.arange(max(jx)+2))
    qx = (vec == 1).nonzero()
    e = edges[ix[qx]]
    node3 = node3[ix[qx]]

    # Orientation
    v1 = p[e[:,1]]-p[e[:,0]]
    v2 = p[node3]-p[e[:,1]]
    ix, = (v1[:,0]*v2[:,1]-v1[:,1]*v2[:,0] > 0).nonzero()
    e[ix,:] = e[ix, ::-1]

    return e

def boundedgesnd(t):
    """Find boundary edges in n-dims, but no guarantees of orientation"""
    dim = t.shape[1]-1
    vertices = set(range(dim+1))
    faces = list(itertools.combinations(vertices, dim))

    edges = np.vstack((t[:,face] for face in faces))
    edges.sort(1)
    _, ix, jx = ml.unique_rows(edges, True, True)
    vec, _ = np.histogram(jx, np.arange(max(jx)+2))
    qx = (vec == 1).nonzero()
    e = edges[ix[qx]]

    return e

def uniref(p, t, nref=1, fd=None):
    """Uniformly refine simplicial mesh.

    Usage
    -----
    >>> p, t = uniref(p, t, nref, fd)

    Parameters
    ----------
    p : array, shape (np, dim)
        Nodes
    t : array, shape (nt, dim+1)
        Triangulation
    nref : int, optional
        Number of uniform refinements
    fd : callable as fd(p), optional
        Boundary distance function
    """

    dim = p.shape[1]
    assert dim == t.shape[1]-1

    for i in range(nref):
        n = p.shape[0]
        nt = t.shape[0]

        if dim == 1:
            pmid = (p[t[:,0]]+p[t[:,1]])/2
            t1 = t[:,[0]]
            t2 = t[:,[1]]
            t12 = arange(n, n+nt)[:,None]
            t = np.vstack((np.hstack((t1,t12)),
                           np.hstack((t12,t2))))
            p = np.vstack((p, pmid))
        elif dim == 2:
            pair = np.vstack((t[:,[0,1]],
                              t[:,[0,2]],
                              t[:,[1,2]]))
            pair.sort(1)
            pair, pairj = ml.unique_rows(pair, return_inverse=True)
            pmid = (p[pair[:,0]] + p[pair[:,1]])/2
            t1 = t[:,[0]]
            t2 = t[:,[1]]
            t3 = t[:,[2]]
            t12 = pairj[0*nt:1*nt, None] + n
            t13 = pairj[1*nt:2*nt, None] + n
            t23 = pairj[2*nt:3*nt, None] + n
            t = np.vstack((np.hstack((t1,t12,t13)),
                           np.hstack((t12,t23,t13)),
                           np.hstack((t2,t23,t12)),
                           np.hstack((t3,t13,t23))))
            p = np.vstack((p,pmid))
        else:
            raise NotImplementedError

        if fd is not None:
            for i in range(5):
                bndproj(p, t, fd)

    return p, t

def bndproj(p, t, fd):
    """Project boundary points to true boundary.

    Updates p in place.

    Parameters
    ----------
    p : array, shape (np, dim)
        nodes
    t : array, shape (nt, dim+1)
        triangulation
    fd : callable, as fd(p)
        distance function
    """

    deps = np.sqrt(np.finfo(np.double).eps)*max(p.max(0)-p.min(0))

    dim = p.shape[1]

    if dim == 2:
        e = boundedges(p,t)
        e = np.unique(e.flat)

        d = fd(p[e])
        dgradx = (fd(p[e]+[deps,0])-d)/deps
        dgrady = (fd(p[e]+[0,deps])-d)/deps
        dgrad2 = dgradx**2 + dgrady**2
        dgrad2[dgrad2 == 0] = 1
        p[e] -= (d*np.vstack((dgradx, dgrady))/dgrad2).T

    elif dim == 3:
        if t.shape[1] == 3:
            tri = t
        else:
            tri = surftri(p, t)
        tri = np.unique(tri.flat)

        d = fd(p[tri])
        dgradx = (fd(p[e]+[deps,0,0])-d)/deps
        dgrady = (fd(p[e]+[0,deps,0])-d)/deps
        dgradz = (fd(p[e]+[0,0,deps])-d)/deps
        dgrad2 = dgradx**2 + dgrady**2 + dgradz**2
        dgrad2[dgrad2 == 0] = 1
        p[e] -= (d*np.vstack((dgradx, dgrady, dgradz))/dgrad2).T

    else:
        raise NotImplementedError

def mkt2t(t):
    """Compute element connectivities from element indices.

    t2t, t2n = mkt2t(t)
    """

    nt = t.shape[0]
    dim = t.shape[1]-1

    if dim == 1:
        edges = np.vstack((t[:,[1]],
                           t[:,[0]]))
    elif dim == 2:
        edges = np.vstack((t[:,[1,2]],
                           t[:,[2,0]],
                           t[:,[0,1]]))
    elif dim == 3:
        edges = np.vstack((t[:,[1,2,3]],
                           t[:,[2,3,0]],
                           t[:,[3,0,1]],
                           t[:,[0,1,2]]))
    else:
        raise NotImplementedError

    # Each row of ts is the index of the entry of t when reading the array.
    ts = np.indices(t.shape).reshape(t.ndim,-1,order='F').T

    edges.sort(1)
    _, jx = ml.unique_rows(edges, return_inverse=True)
    ix = jx.argsort()

    jx = jx[ix]
    ts = ts[ix]

    ix, = (np.diff(jx) == 0).nonzero()

    ts1 = ts[ix]
    ts2 = ts[ix+1]

    t2t = np.empty((nt, dim+1), dtype=int); t2t.fill(-1)
    t2n = np.empty((nt, dim+1), dtype=int); t2n.fill(-1)

    t2t[ts1[:,0], ts1[:,1]] = ts2[:,0]
    t2t[ts2[:,0], ts2[:,1]] = ts1[:,0]

    t2n[ts1[:,0], ts1[:,1]] = ts2[:,1]
    t2n[ts2[:,0], ts2[:,1]] = ts1[:,1]

    return t2t, t2n

def assert_t2t_t2n(t2t, t2n):
    """Raises an AssertionError if t2t/t2n are not consistent."""

    # Negative entries correspond to boundaries. Check that t2t and t2n
    # agree on where the boundaries are.
    mt2t = np.ma.masked_less(t2t, 0)
    mt2n = np.ma.masked_less(t2n, 0)
    assert (mt2t.mask == mt2n.mask).all(), "inconsistent boundaries"

    # Build arrays which satisfy t[i,j] = i, n[i,j] = j
    nt, dim1 = t2t.shape
    t, n = np.mgrid[:nt, :dim1]
    mt = np.ma.masked_array(t, mask=mt2t.mask)
    mn = np.ma.masked_array(n, mask=mt2n.mask)

    # For each non-boundary edge 'j' of element 'i', we expect that
    # t2t[t2t[i,j], t2n[i,j]] == i
    # t2n[t2t[i,j], t2n[i,j]] == j
    assert (t2t[mt2t.compressed(), mt2n.compressed()]
            == mt.compressed()).all(), "inconsistent t2t"
    assert (t2n[mt2t.compressed(), mt2n.compressed()]
            == mn.compressed()).all(), "inconsistent t2n"

def circumcenter(p, t):
    """Find circumcenter of each triangle.

    Parameters
    ----------
    p : array, shape (np, 2)
        nodes
    t : array, shape (nt, 3)
        triangulation

    Returns
    -------
    pc : array, shape (nt, 2)
        circumcenters
    r : array, shape (nt, )
        radii
    """
    nt = len(t)
    pc = np.zeros((nt, 2))
    r = np.zeros((nt,1))

    for it in range(nt):
        ct = t[it]
        dp1 = p[ct[1]]-p[ct[0]]
        dp2 = p[ct[2]]-p[ct[0]]

        mid1 = (p[ct[1]]+p[ct[0]])/2
        mid2 = (p[ct[2]]+p[ct[0]])/2

        s = np.linalg.solve(np.array([[-dp1[1], dp2[1]],
                                      [dp1[0], -dp2[0]]]),
                            -mid1+mid2)
        cpc = mid1+s[0]*np.array([-dp1[1],dp1[0]])
        cr = np.linalg.norm(p[ct[0]]-cpc)

        pc[it] = cpc
        r[it] = cr

    return pc, r

def uniformity(p, t, fh):
    """Uniformity measure: how close the element sizes in the mesh are to the
    desired mesh size function.
    """
    pc, r = circumcenter(p, t)
    sz = r / fh(pc)
    return np.std(sz)/np.mean(sz)



def simpvol(p, t):
    """Signed volumes of the simplex elements in the mesh."""
    dim = p.shape[1]
    if dim == 1:
        d01 = p[t[:,1]]-p[t[:,0]]
        return d01
    elif dim == 2:
        d01 = p[t[:,1]]-p[t[:,0]]
        d02 = p[t[:,2]]-p[t[:,0]]
        return (d01[:,0]*d02[:,1]-d01[:,1]*d02[:,0])/2
    else:
        raise NotImplementedError

def simpqual(p, t):
    """Simplex quality.

    Usage
    -----
    q = simpqual(p, t)

    Parameters
    ----------
    p : array, shape (np, dim)
        nodes
    t : array, shape (nt, dim+1)
        triangulation

    Returns
    -------
    q : array, shape (nt, )
        qualities
    """
    assert (p.ndim == 2
            and t.ndim == 2
            and p.shape[1]+1 == t.shape[1])

    dim = p.shape[1]
    if dim == 1:
        return np.ones((1, nt))

    elif dim == 2:
        length = lambda p1: np.sqrt((p1**2).sum(1))
        a = length(p[t[:,1]]-p[t[:,0]])
        b = length(p[t[:,2]]-p[t[:,0]])
        c = length(p[t[:,2]]-p[t[:,1]])
        r = 0.5*np.sqrt((b+c-a)*(c+a-b)*(a+b-c)/(a+b+c))
        R = a*b*c/np.sqrt((a+b+c)*(b+c-a)*(c+a-b)*(a+b-c))
        return 2*r/R

    else:
        raise NotImplementedError

def fixmesh(p, t, ptol=2e-13):
    """Remove duplicated/unused nodes and fix element orientation.

    Parameters
    ----------
    p : array, shape (np, dim)
    t : array, shape (nt, nf)

    Usage
    -----
    p, t = fixmesh(p, t, ptol)
    """
    snap = (p.max(0)-p.min(0)).max()*ptol
    _, ix, jx = ml.unique_rows(np.round(p/snap)*snap, True, True)

    p = p[ix]
    t = jx[t]

    flip = simpvol(p,t)<0
    t[flip, :2] = t[flip, 1::-1]

    return p, t
