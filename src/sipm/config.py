"""Run-time configuration: where Ansys lives and where working files go.

Nothing here is hard-coded to a machine. Ansys installations are discovered from
the standard environment variables (``AWP_ROOTxyz`` for MAPDL, ``ANSYSEM_ROOTxyz``
for Electronics Desktop); the newest is used unless a version or an explicit path
is given. The working directory defaults to ``<tempdir>/sipm_runs`` and can be set
with the ``SIPM_WORKDIR`` environment variable or per call.
"""

from __future__ import annotations

import os
import platform
import re
import tempfile
from dataclasses import dataclass, field
from pathlib import Path


def _versions(prefix: str) -> dict[str, str]:
    out = {}
    for k, v in os.environ.items():
        m = re.fullmatch(prefix + r"(\d{3})", k.upper())
        if m and v:
            out[m.group(1)] = v
    return out


def find_mapdl(version: str | None = None) -> Path | None:
    """Path of the MAPDL executable, or None if no installation is found."""
    vers = _versions("AWP_ROOT")
    if not vers:
        return None
    ver = version or max(vers)
    if ver not in vers:
        raise FileNotFoundError(f"AWP_ROOT{ver} is not set; found versions {sorted(vers)}")
    root = Path(vers[ver])
    if platform.system() == "Windows":
        exe = root / "ansys" / "bin" / "winx64" / f"ANSYS{ver}.exe"
    else:
        exe = root / "ansys" / "bin" / f"ansys{ver}"
    return exe if exe.exists() else None


def find_aedt_version(version: str | None = None) -> str | None:
    """AEDT version string such as '2024.1' (newest installed unless given)."""
    if version:
        return version
    vers = _versions("ANSYSEM_ROOT")
    if not vers:
        return None
    v = max(vers)
    return f"20{v[:2]}.{v[2]}"


def default_workdir() -> Path:
    return Path(os.environ.get("SIPM_WORKDIR") or Path(tempfile.gettempdir()) / "sipm_runs")


@dataclass
class Config:
    """Execution settings shared by the drivers and the sweep runner."""

    workdir: Path = field(default_factory=default_workdir)
    mapdl_exe: Path | None = None        # None -> discover
    mapdl_version: str | None = None     # e.g. "241"
    mapdl_cores: int = 2
    mapdl_timeout: float = 1800.0        # s per point
    aedt_version: str | None = None      # e.g. "2024.1"; None -> newest installed
    aedt_non_graphical: bool = True
    aedt_cores: int = 4
    keep_aedt_projects: bool = False     # delete .aedt files after a point is solved
    keep_mapdl_files: bool = False       # delete .rst/.db etc. after a point is solved
    # memory protection (see sipm.memory) ---------------------------------------------
    min_free_ram_gb: float = 4.0         # wait before a point while free RAM is below this
    min_free_ram_fraction: float = 0.0   # ... or below this fraction of total RAM (larger wins)
    ram_wait_timeout_s: float = 1800.0   # then stop the sweep cleanly (resumable)
    ram_poll_s: float = 10.0
    max_solver_ram_gb: float | None = None  # kill a solver tree above this RSS; point -> failed
    max_parallel: int = 1                # points solved at the same time (only 1 is supported)
    aedt_restart_every: int = 10         # release and restart AEDT every N EM points (0 = never)

    def __post_init__(self):
        self.workdir = Path(self.workdir)
        if self.max_parallel != 1:
            raise ValueError("max_parallel: only sequential sweeps (1) are supported")
        if self.min_free_ram_gb < 0 or not 0 <= self.min_free_ram_fraction < 1:
            raise ValueError("min_free_ram_gb must be >= 0 and min_free_ram_fraction in [0, 1)")
        if self.max_solver_ram_gb is not None and self.max_solver_ram_gb <= 0:
            raise ValueError("max_solver_ram_gb must be positive (or None)")
        if self.aedt_restart_every < 0:
            raise ValueError("aedt_restart_every must be >= 0")

    def resolved_mapdl(self) -> Path:
        exe = Path(self.mapdl_exe) if self.mapdl_exe else find_mapdl(self.mapdl_version)
        if exe is None or not Path(exe).exists():
            raise FileNotFoundError(
                "MAPDL executable not found: set AWP_ROOTxyz or Config.mapdl_exe")
        return Path(exe)

    def resolved_aedt_version(self) -> str:
        v = find_aedt_version(self.aedt_version)
        if v is None:
            raise FileNotFoundError("AEDT not found: set ANSYSEM_ROOTxyz or Config.aedt_version")
        return v


_MATLAB_VARS = ("KMP_STACKSIZE", "KMP_BLOCKTIME", "KMP_HANDLE_SIGNALS")


def child_environment(base: dict | None = None) -> dict:
    """Environment for Ansys child processes.

    When Python runs under MATLAB (``pyenv``), MATLAB's OpenMP settings (notably
    ``KMP_STACKSIZE=512k``) and its library folders on PATH are inherited; MAPDL's
    solver threads then overflow their stacks. These are removed here.
    """
    env = dict(os.environ if base is None else base)
    if not (env.get("MW_RESOURCE_PATH") or env.get("MATLAB_ARCH")):
        return env
    for k in _MATLAB_VARS:
        env.pop(k, None)
    sep = os.pathsep
    roots = []
    if env.get("MW_RESOURCE_PATH"):
        roots.append(env["MW_RESOURCE_PATH"])
    for p in env.get("PATH", "").split(sep):
        # <matlabroot>/bin/win64 (system() children) or <matlabroot>/extern/bin/win64 (pyenv host)
        for up in (os.path.join(p, os.pardir, os.pardir), os.path.join(p, os.pardir, os.pardir, os.pardir)):
            if p and os.path.isfile(os.path.join(up, "bin", "matlab.exe")):
                roots.append(up)
    roots = [os.path.normcase(os.path.abspath(r)) for r in roots]
    keep = []
    for p in env.get("PATH", "").split(sep):
        if p and not any(os.path.normcase(os.path.abspath(p)).startswith(r) for r in roots):
            keep.append(p)
    env["PATH"] = sep.join(keep)
    return env


def sanitize_process_environment() -> None:
    """Apply :func:`child_environment` to this process (for launchers we do not control)."""
    clean = child_environment()
    for k in list(os.environ):
        if k not in clean:
            del os.environ[k]
    os.environ.update(clean)


def ansys_available(kind: str = "any") -> bool:
    """True if the requested Ansys product can be located (licence not checked)."""
    has_mapdl = find_mapdl() is not None
    has_aedt = bool(_versions("ANSYSEM_ROOT"))
    if kind == "mapdl":
        return has_mapdl
    if kind == "aedt":
        return has_aedt
    return has_mapdl or has_aedt
