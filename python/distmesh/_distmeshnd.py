# encoding: utf-8
"""DistMesh ND"""

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
import scipy.spatial as spspatial

# Local imports
import distmesh.mlcompat as ml

__all__ = ['distmeshnd']

#-----------------------------------------------------------------------------
# Functions
#-----------------------------------------------------------------------------

def distmeshnd(fd, fh, h0, bbox, pfix=None, fig='gcf'):
    """
    distmeshnd: N-D Mesh Generator using Distance Functions.

    Usage
    -----
    >>> p, t = distmesh2d(fd, fh, h0, bbox, pfix)

    Parameters
    ----------
    fd:        Distance function d(x,y)
    fh:        Scaled edge length function h(x,y)
    h0:        Initial edge length
    bbox:      Bounding box, (xmin, ymin, zmin, ..., xmax, ymax, zmax, ...)
    pfix:      Fixed node positions, shape (nfix, dim)
    fig:       Figure to use for plotting, or None to disable plotting.

    Returns
    -------
    p:         Node positions (np, dim)
    t:         Triangle indices (nt, dim+1)

    Example: (Unit ball)
    >>> dim = 3
    >>> fd = lambda p: sqrt((p**2).sum(1))-1.0
    >>> bbox = np.vstack((-np.ones(dim), np.ones(dim)))
    >>> p, t = distmeshnd(fd, huniform, 2, bbox)
    """

    if fig == 'gcf':
        import matplotlib.pyplot as plt
        fig = plt.gcf()

    bbox = np.array(bbox).reshape(2, -1)
    dim = bbox.shape[1]

    ptol=.001; ttol=.1; L0mult=1+.4/2**(dim-1); deltat=.1; geps=1e-1*h0;
    deps=np.sqrt(np.finfo(np.double).eps)*h0;

    if pfix is not None:
        pfix = np.array(pfix, dtype='d')
        nfix = len(pfix)
    else:
        pfix = np.empty((0, dim))
        nfix = 0

    # 0. Prepare a figure.
    if fig is not None:
        if dim == 2:
            from distmesh.plotting import SimplexCollection
            fig.clf()
            ax = fig.gca()
            c = SimplexCollection()
            ax.add_collection(c)
            ax.set_xlim(bbox[:,0])
            ax.set_ylim(bbox[:,1])
            ax.set_aspect('equal')
            ax.set_axis_off()
            fig.canvas.draw()
        elif dim == 3:
            import mpl_toolkits.mplot3d
            from distmesh.plotting import axes_simpplot3d
            fig.clf()
            ax = fig.gca(projection='3d')
            ax.set_xlim(bbox[:,0])
            ax.set_ylim(bbox[:,1])
            ax.set_zlim(bbox[:,2])
            fig.canvas.draw()
        else:
            print("Plotting only supported in dimensions 2 and 3.")
            fig = None

    # 1. Create initial distribution in bounding box (equilateral triangles)
    p = np.mgrid[tuple(slice(min, max+h0, h0) for min, max in bbox.T)]
    p = p.reshape(dim, -1).T

    # 2. Remove points outside the region, apply the rejection method
    p = p[fd(p)<geps]                                # Keep only d<0 points
    r0 = fh(p)
    p = np.vstack((pfix,
                   p[np.random.rand(p.shape[0]) < r0.min()**dim / r0**dim]))
    N = p.shape[0]

    count = 0
    pold = float('inf')                              # For first iteration

    while True:

        # 3. Retriangulation by the Delaunay algorithm
        dist = lambda p1, p2: np.sqrt(((p1-p2)**2).sum(1))
        if (dist(p, pold)/h0).max() > ttol:          # Any large movement?
            pold = p.copy()                          # Save current positions
            t = spspatial.Delaunay(p).vertices       # List of triangles
            pmid = p[t].sum(1)/(dim+1)               # Compute centroids
            t = t[fd(pmid) < -geps]                  # Keep interior triangles
            # 4. Describe each bar by a unique pair of nodes
            bars = np.vstack((t[:,pair] for pair in
                              itertools.combinations(range(dim+1), 2)))
            bars.sort(axis=1)
            bars = ml.unique_rows(bars)              # Bars as node pairs
            # 5. Graphical output of the current mesh
            if fig is not None:
                if dim == 2:
                    c.set_simplices((p,t))
                    fig.canvas.draw()

                elif dim == 3:
                    if count % 5 == 0:
                        ax.cla()
                        axes_simpplot3d(ax, p, t, p[:,1] > 0)
                        ax.set_title('Retriangulation #%d' % count)
                        fig.canvas.draw()
            else:
                print('Retriangulation #%d' % count)
            count += 1

        # 6. Move mesh points based on bar lengths L and forces F
        barvec = p[bars[:,0]] - p[bars[:,1]]         # List of bar vectors
        L = np.sqrt((barvec**2).sum(1))              # L = Bar lengths
        hbars = fh(p[bars].sum(1)/2)
        L0 = (hbars*L0mult*((L**dim).sum()/(hbars**dim).sum())**(1.0/dim))
        F = L0-L; F[F<0] = 0                         # Bar forces (scalars)
        Fvec = F[:,None]/L[:,None].dot(np.ones((1,dim)))*barvec # Bar forces (x,y components)
        Ftot = ml.dense(bars[:,[0]*dim+[1]*dim],
                        np.repeat([list(range(dim))*2], len(F), axis=0),
                        np.hstack((Fvec, -Fvec)),
                        shape=(N, dim))
        Ftot[:nfix] = 0                              # Force = 0 at fixed points
        p += deltat*Ftot                             # Update node positions

        # 7. Bring outside points back to the boundary
        d = fd(p); ix = d>0                          # Find points outside (d>0)
        if ix.any():
            def deps_vec(i): a = [0]*dim; a[i] = deps; return a
            dgrads = [(fd(p[ix]+deps_vec(i))-d[ix])/deps for i in range(dim)]
            dgrad2 = sum(dgrad**2 for dgrad in dgrads)
            p[ix] -= (d[ix]*np.vstack(dgrads)/dgrad2).T # Project

        # 8. Termination criterion: All interior nodes move less than dptol (scaled)
        maxdp = deltat*np.sqrt((Ftot[d<-geps]**2).sum(1)).max()
        if maxdp < ptol*h0:
            break

    return p, t
