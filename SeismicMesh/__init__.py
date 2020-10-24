# encoding: utf-8
# -----------------------------------------------------------------------------
#  Copyright (C) 2020 Keith Jared Roberts

#  Distributed under the terms of the GNU General Public License. You should
#  have received a copy of the license along with this program. If not,
#  see <http://www.gnu.org/licenses/>.

from . import decomp, geometry, migration
from .geometry import (
    Ball,
    Disk,
    Cube,
    Cylinder,
    Rectangle,
    Union,
    Intersection,
    Difference,
    Torus,
    Prism,
)
from .generation import generate_mesh, sliver_removal
from .sizing import (
    get_sizing_function_from_segy,
    write_velocity_model,
    plot_sizing_function,
)
from .__about__ import __version__


__all__ = [
    "__version__",
    "geometry",
    "Rectangle",
    "Cube",
    "Cylinder",
    "Disk",
    "Union",
    "Torus",
    "Prism",
    "Ball",
    "Intersection",
    "Difference",
    "get_sizing_function_from_segy",
    "write_velocity_model",
    "plot_sizing_function",
    "generate_mesh",
    "sliver_removal",
    "decomp",
    "migration",
]
