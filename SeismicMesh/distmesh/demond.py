# encoding: utf-8
"""Distmesh ND examples."""

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

import numpy as np
import matplotlib.pyplot as plt

# Local imports.
import distmesh as dm

#-----------------------------------------------------------------------------
# Demo functions
#-----------------------------------------------------------------------------

def unit_ball():
    """3-D Unit ball"""
    fd = lambda p: np.sqrt((p**2).sum(1))-1.0
    return dm.distmeshnd(fd, dm.huniform, 0.2, (-1,-1,-1, 1,1,1))

def cylinder_with_hole():
    """Cylinder with hole"""
    def fd10(p):
        r = np.sqrt(p[:,0]**2 + p[:,1]**2)
        z = p[:,2]

        d1 = r-1.0
        d2 = z-1.0
        d3 = -z-1.0
        d4 = np.sqrt(d1**2+d2**2)
        d5 = np.sqrt(d1**2+d3**2)
        d = dm.dintersect(dm.dintersect(d1, d2), d3)

        ix = (d1>0)*(d2>0)
        d[ix] = d4[ix]
        ix = (d1>0)*(d3>0)
        d[ix] = d5[ix]

        return dm.ddiff(d, dm.dsphere(p, 0,0,0, 0.5))

    def fh10(p):
        h1 = 4*np.sqrt((p**2).sum(1))-1.0
        return np.minimum(h1, 2.0)

    return dm.distmeshnd(fd10, fh10, 0.1, (-1,-1,-1, 1,1,1))

def meshdemond(pause=None):
    """Run all Distmesh 2D examples."""
    if pause is None:
        pause = lambda : None

    plt.ion()
    np.random.seed(1) # Always the same results

#    def fstats(p, t):
#        print('%d nodes, %d elements, min quality %.2f'
#              % (len(p), len(t), dm.simpqual(p,t).min()))

    fstats = lambda p, t: None

    print('3-D Unit ball')
    p, t = unit_ball()
    fstats(p, t)
    pause(); print('')

    print('Cylinder with hole')
    p, t = cylinder_with_hole()
    fstats(p, t)
    pause(); print('')

if __name__ == "__main__":
    # Py3k
    try: input = raw_input
    except: pass

    pause = lambda : input('(press enter to continue)')
    meshdemond(pause)
