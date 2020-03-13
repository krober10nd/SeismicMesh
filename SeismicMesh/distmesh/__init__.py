# encoding: utf-8
"""PyDistMesh: A Simple Mesh Generator in Python

https://github.com/bfroehle/pydistmesh
"""

#-----------------------------------------------------------------------------
#  Copyright (C) 2004-2012 Per-Olof Persson
#  Copyright (C) 2012 Bradley Froehle

#  Distributed under the terms of the GNU General Public License. You should
#  have received a copy of the license along with this program. If not,
#  see <http://www.gnu.org/licenses/>.
#-----------------------------------------------------------------------------

__version__ = '1.2'

#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------

from distmesh._distmesh2d import distmesh2d
from distmesh._distmeshnd import distmeshnd
from distmesh.utils import *
from distmesh.plotting import *
