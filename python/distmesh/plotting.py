# encoding: utf-8
"""Plotting routines."""

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

import collections
import functools

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cbook
from matplotlib.collections import PathCollection
from matplotlib.patches import PathPatch
from matplotlib.path import Path

# Local imports.
import distmesh.utils as dmutils

# Py3k
try: range = xrange
except: pass

__all__ = ['SimplexCollection', 'simpplot']

#-----------------------------------------------------------------------------
# Classes
#-----------------------------------------------------------------------------

class SimplexCollection(PathCollection):
    """A collection of triangles."""
    def __init__(self, simplices=None, **kwargs):
        kwargs.setdefault('linewidths', 0.5)
        kwargs.setdefault('edgecolors', 'k')
        kwargs.setdefault('facecolors', (0.8, 0.9, 1.0))
        PathCollection.__init__(self, [], **kwargs)
        if simplices is not None:
            self.set_simplices(simplices)

    def set_simplices(self, simplices):
        """Usage: set_simplices((p, t))"""
        p, t = simplices
        code = [Path.MOVETO, Path.LINETO, Path.LINETO, Path.CLOSEPOLY]
        self.set_paths([Path(edge, code) for edge in p[t[:,[0,1,2,0]]]])

#-----------------------------------------------------------------------------
# 2D Plotting
#-----------------------------------------------------------------------------

def axes_simpplot2d(ax, p, t, nodes=False, annotate='', **kwargs):
    """Plot a triangulation

    Parameters
    ----------
    p : array, shape (np, 2)
        nodes
    t : array, shape (nt, 3)
        simplices
    nodes : bool, optional
        draw a marker at each node
    annotate : str, optional
        'p' : annotate nodes
        't' : annotate simplices
    **kwargs : dict
        additional arguments to pass to SimplexCollection
    """
    scalex = kwargs.pop('scalex', True)
    scaley = kwargs.pop('scaley', True)
    if not ax._hold: ax.cla()

    assert p.shape[1] == 2

    c = SimplexCollection((p, t), **kwargs)
    ax.add_collection(c)
    if nodes:
        ax.plot(p[:,0], p[:,1], '.k', markersize=16)
    if 'p' in annotate:
        for i in range(len(p)):
            ax.annotate(str(i), p[i], ha='center', va='center')
    if 't' in annotate:
        for it in range(len(t)):
            pmid = p[t[it]].mean(0)
            ax.annotate(str(it), pmid, ha='center', va='center')
    ax.set_axis_off()
    ax.set_aspect('equal')
    ax.autoscale_view(scalex=scalex, scaley=scaley)
    return c

#-----------------------------------------------------------------------------
# 3D Plotting
#-----------------------------------------------------------------------------

def _mkfaces(t):
    """All exterior faces, as a set."""
    return set(tuple(sorted(f)) for f in dmutils.boundedgesnd(t))

def _trimesh(ax, t, x, y, z, **kwargs):
    """Display a 3D triangular mesh.

    Ignores ax._hold.
    """
    from mpl_toolkits.mplot3d.art3d import pathpatch_2d_to_3d
    patches = []
    code = [Path.MOVETO, Path.LINETO, Path.LINETO, Path.CLOSEPOLY]
    for f in t:
        bdry = np.take(f, range(4), mode='wrap')
        pp = PathPatch(Path(np.column_stack((x[bdry], y[bdry])), code),
                       **kwargs)
        ax.add_patch(pp)
        pathpatch_2d_to_3d(pp, z[bdry])
        patches.append(pp)
    return patches

def axes_simpplot3d(ax, p, t, pmask=None, **kwargs):
    """Plot a surface or volume triangulation.

    Parameters
    ----------
    p : array, shape (np, 3)
    t : array, shape (nt, 3) or (nt, 4)
    pmask : callable or bool array of shape (np,)

    Additional keyword arguments
    ----------------------------
    facecolor : facecolor
    ifacecolor : facecolor for faces exposed by pmask
    """
    if not ax._hold: ax.cla()
    had_data = ax.has_data()

    facecolor = kwargs.pop('facecolor', (0.8, 0.9, 1.0))
    ifacecolor = kwargs.pop('ifacecolor', (0.9, 0.8, 1.0))
    xs, ys, zs = p.T

    ret = cbook.silent_list('mpl_toolkits.mplot3d.art3d.PathPatch3D')

    if t.shape[1] == 4:
        tri1 = _mkfaces(t)

        if pmask is not None:
            if isinstance(pmask, collections.Callable):
                pmask = pmask(p)
            t = t[pmask[t].any(1)]
            tri2 = _mkfaces(t)
            tri1 = tri1.intersection(tri2)
            tri2 = tri2.difference(tri1)
            c = _trimesh(ax, tri2, xs, ys, zs, facecolor=ifacecolor)
            ret.extend(c)
    else:
        tri1 = t
        if pmask is not None:
            if isinstance(pmask, collections.Callable):
                pmask = pmask(p)
            tri1 = t[pmask[t].any(1)]

    c = _trimesh(ax, tri1, xs, ys, zs, facecolor=facecolor)
    ret.extend(c)
    ax.auto_scale_xyz(xs, ys, zs, had_data)

    return ret

#-----------------------------------------------------------------------------
# pyplot interface
#-----------------------------------------------------------------------------

def simpplot(p, t, *args, **kwargs):
    """Plot a simplicial mesh

    Parameters
    ----------
    p : array, shape (np, dim)
        nodes
    t : array, shape (nt, dim+1)
        elements

    Additional 2D parameters
    ------------------------
    nodes : bool, optional
        draw a marker at each node
    annotate : str, optional
        'p' : annotate nodes
        't' : annotate simplices

    Additional 3D parameters
    ------------------------
    pmask : callable or bool array of shape (np,)
    """
    ax = plt.gca()
    # allow callers to override the hold state by passing hold=True|False
    washold = ax.ishold()
    hold = kwargs.pop('hold', None)
    if hold is not None:
        ax.hold(hold)
    try:
        dim = p.shape[1]
        if dim == 2:
            ret = axes_simpplot2d(ax, p, t, *args, **kwargs)
        elif dim == 3:
            import mpl_toolkits.mplot3d
            ax = plt.gca(projection='3d')
            ret = axes_simpplot3d(ax, p, t, *args, **kwargs)
        else:
            raise NotImplementedError("Unknown dimension.")
        plt.draw_if_interactive()
    finally:
        ax.hold(washold)
    return ret
