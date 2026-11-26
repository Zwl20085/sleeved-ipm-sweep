"""sleeved-ipm-sweep: coupled Maxwell 2D / MAPDL parametric sweeps for sleeved IPM rotors.

Licence-free modules: params, geometry, materials, analytic, sweep (with a custom
solver), post, mapdl deck writing. Solving needs Ansys Electronics Desktop (via
PyAEDT) for the EM part and Ansys Mechanical APDL for the stress part.
"""

from .config import Config
from .geometry import GeometryError, build_rotor, build_stator
from .params import TOPOLOGIES, Params, default_params
from .sweep import Sweep, grid

__version__ = "0.1.0"

__all__ = ["Config", "GeometryError", "Params", "Sweep", "TOPOLOGIES", "build_rotor", "build_stator",
           "default_params", "grid", "__version__"]
