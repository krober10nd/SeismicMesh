# encoding: utf-8
"""DistMesh 2D"""

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

import numpy as np
import scipy.spatial as spspatial

# Local imports
import distmesh.mlcompat as ml
import distmesh.utils as dmutils

__all__ = ['distmesh2d']

#-----------------------------------------------------------------------------
# Functions
#-----------------------------------------------------------------------------

def distmesh2d(fd, fh, h0, bbox, pfix=None, fig='gcf', dptol=.001, ttol=.1,
               Fscale=1.2, deltat=.2, geps_multiplier=.001, densityctrlfreq=30):
    """
    distmesh2d: 2-D Mesh Generator using Distance Functions.

    Usage
    -----
    >>> p, t = distmesh2d(fd, fh, h0, bbox, pfix)

    Parameters
    ----------
    fd:        Distance function d(x,y)
    fh:        Scaled edge length function h(x,y)
    h0:        Initial edge length
    bbox:      Bounding box, (xmin, ymin, xmax, ymax)
    pfix:      Fixed node positions, shape (nfix, 2)
    fig:       Figure to use for plotting, or None to disable plotting.

    Returns
    -------
    p:         Node positions (Nx2)
    t:         Triangle indices (NTx3)

    Example: (Uniform Mesh on Unit Circle)
    >>> fd = lambda p: sqrt((p**2).sum(1))-1.0
    >>> p, t = distmesh2d(fd, huniform, 2, (-1,-1,1,1))

    Example: (Rectangle with circular hole, refined at circle boundary)
    >>> fd = lambda p: ddiff(drectangle(p,-1,1,-1,1), dcircle(p,0,0,0.5))
    >>> fh = lambda p: 0.05+0.3*dcircle(p,0,0,0.5)
    >>> p, t = distmesh2d(fd, fh, 0.05, (-1,-1,1,1),
                          [(-1,-1), (-1,1), (1,-1), (1,1)])

    Example: (Polygon)
    >>> pv=[(-0.4, -0.5), (0.4, -0.2), (0.4, -0.7), (1.5, -0.4), (0.9, 0.1),
            (1.6, 0.8), (0.5, 0.5), (0.2, 1.0), (0.1, 0.4), (-0.7, 0.7),
            (-0.4, -0.5)]
    >>> fd = lambda p: dpoly(p, pv)
    >>> p, t = distmesh2d(fd, huniform, 0.1, (-1,-1, 2,1), pv)

    Example: (Ellipse)
    >>> fd = lambda p: p[:,0]**2/2**2 + p[:,1]**2/1**2 - 1
    >>> p, t = dm.distmesh2d(fd, dm.huniform, 0.2, (-2,-1, 2,1))

    Example: (Square, with size function point and line sources)
    >>> fd = lambda p: dm.drectangle(p,0,1,0,1)
    >>> fh = lambda p: np.minimum(np.minimum(
            0.01+0.3*abs(dm.dcircle(p,0,0,0)),
            0.025+0.3*abs(dm.dpoly(p,[(0.3,0.7),(0.7,0.5)]))), 0.15)
    >>> p, t = dm.distmesh2d(fd, fh, 0.01, (0,0,1,1), [(0,0),(1,0),(0,1),(1,1)])

    Example: (NACA0012 airfoil)
    >>> hlead=0.01; htrail=0.04; hmax=2; circx=2; circr=4
    >>> a=.12/.2*np.array([0.2969,-0.1260,-0.3516,0.2843,-0.1036])
    >>> a0=a[0]; a1=np.hstack((a[5:0:-1], 0.0))
    >>> fd = lambda p: dm.ddiff(
        dm.dcircle(p,circx,0,circr),
        (abs(p[:,1])-np.polyval(a1, p[:,0]))**2-a0**2*p[:,0])
    >>> fh = lambda p: np.minimum(np.minimum(
            hlead+0.3*dm.dcircle(p,0,0,0),
            htrail+0.3*dm.dcircle(p,1,0,0)), hmax)

    >>> fixx = 1.0-htrail*np.cumsum(1.3**np.arange(5))
    >>> fixy = a0*np.sqrt(fixx)+np.polyval(a1, fixx)
    >>> fix = np.vstack((
            np.array([(circx-circr,0),(circx+circr,0),
                      (circx,-circr),(circx,circr),
                      (0,0),(1,0)]),
            np.vstack((fixx, fixy)).T,
            np.vstack((fixx, -fixy)).T))
    >>> box = (circx-circr,-circr, circx+circr,circr)
    >>> h0 = min(hlead, htrail, hmax)
    >>> p, t = dm.distmesh2d(fd, fh, h0, box, fix)
    """

    if fig == 'gcf':
        import matplotlib.pyplot as plt
        fig = plt.gcf()

    geps = geps_multiplier * h0
    deps=np.sqrt(np.finfo(np.double).eps)*h0;

    # Extract bounding box
    xmin, ymin, xmax, ymax = bbox
    if pfix is not None:
        pfix = np.array(pfix, dtype='d')

    # 0. Prepare a figure.
    if fig is not None:
        from distmesh.plotting import SimplexCollection
        fig.clf()
        ax = fig.gca()
        c = SimplexCollection()
        ax.add_collection(c)
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)
        ax.set_aspect('equal')
        ax.set_axis_off()
        fig.canvas.draw()

    # 1. Create initial distribution in bounding box (equilateral triangles)
    x, y = np.mgrid[xmin:(xmax+h0):h0,
                    ymin:(ymax+h0*np.sqrt(3)/2):h0*np.sqrt(3)/2]
    x[:, 1::2] += h0/2                               # Shift even rows
    p = np.vstack((x.flat, y.flat)).T                # List of node coordinates

    # 2. Remove points outside the region, apply the rejection method
    p = p[fd(p)<geps]                                # Keep only d<0 points
    r0 = 1/fh(p)**2                                  # Probability to keep point
    p = p[np.random.random(p.shape[0])<r0/r0.max()]  # Rejection method
    if pfix is not None:
        p = ml.setdiff_rows(p, pfix)                 # Remove duplicated nodes
        pfix = ml.unique_rows(pfix); nfix = pfix.shape[0]
        p = np.vstack((pfix, p))                     # Prepend fix points
    else:
        nfix = 0
    N = p.shape[0]                                   # Number of points N

    count = 0
    pold = float('inf')                              # For first iteration

    while True:
        count += 1

        # 3. Retriangulation by the Delaunay algorithm
        dist = lambda p1, p2: np.sqrt(((p1-p2)**2).sum(1))
        if (dist(p, pold)/h0).max() > ttol:          # Any large movement?
            pold = p.copy()                          # Save current positions
            t = spspatial.Delaunay(p).vertices       # List of triangles
            pmid = p[t].sum(1)/3                     # Compute centroids
            t = t[fd(pmid) < -geps]                  # Keep interior triangles
            # 4. Describe each bar by a unique pair of nodes
            bars = np.vstack((t[:, [0,1]],
                              t[:, [1,2]],
                              t[:, [2,0]]))          # Interior bars duplicated
            bars.sort(axis=1)
            bars = ml.unique_rows(bars)              # Bars as node pairs
            # 5. Graphical output of the current mesh
            if fig is not None:
                c.set_simplices((p, t))
                fig.canvas.draw()

        # 6. Move mesh points based on bar lengths L and forces F
        barvec = p[bars[:,0]] - p[bars[:,1]]         # List of bar vectors
        L = np.sqrt((barvec**2).sum(1))              # L = Bar lengths
        hbars = fh(p[bars].sum(1)/2)
        L0 = (hbars*Fscale
              *np.sqrt((L**2).sum()/(hbars**2).sum()))  # L0 = Desired lengths

        # Density control - remove points that are too close
        if (count % densityctrlfreq) == 0 and (L0 > 2*L).any():
            ixdel = np.setdiff1d(bars[L0 > 2*L].reshape(-1), np.arange(nfix))
            p = p[np.setdiff1d(np.arange(N), ixdel)]
            N = p.shape[0]; pold = float('inf')
            continue

        F = L0-L; F[F<0] = 0                         # Bar forces (scalars)
        Fvec = F[:,None]/L[:,None].dot([[1,1]])*barvec # Bar forces (x,y components)
        Ftot = ml.dense(bars[:,[0,0,1,1]],
                        np.repeat([[0,1,0,1]], len(F), axis=0),
                        np.hstack((Fvec, -Fvec)),
                        shape=(N, 2))
        Ftot[:nfix] = 0                              # Force = 0 at fixed points
        p += deltat*Ftot                             # Update node positions

        # 7. Bring outside points back to the boundary
        d = fd(p); ix = d>0                          # Find points outside (d>0)
        if ix.any():
            dgradx = (fd(p[ix]+[deps,0])-d[ix])/deps # Numerical
            dgrady = (fd(p[ix]+[0,deps])-d[ix])/deps # gradient
            dgrad2 = dgradx**2 + dgrady**2
            p[ix] -= (d[ix]*np.vstack((dgradx, dgrady))/dgrad2).T # Project

        # 8. Termination criterion: All interior nodes move less than dptol (scaled)
        if (np.sqrt((deltat*Ftot[d<-geps]**2).sum(1))/h0).max() < dptol:
            break

    # Clean up and plot final mesh
    p, t = dmutils.fixmesh(p, t)

    if fig is not None:
        c.set_simplices((p, t))
        fig.canvas.draw()

    return p, t
