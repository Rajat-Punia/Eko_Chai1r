"""DockQ of a predicted complex against its experimental structure.

wrapper over the DockQ command-line tool
"""

import shutil
import subprocess
from pathlib import Path

from . import config as C


def dockq_available():
    return shutil.which("DockQ") is not None


def run_dockq(predicted_pdb, experimental_pdb, mapping=C.DOCKQ_MAPPING):
    """DockQ score, or None if the tool is absent, a file is missing, or
    DockQ itself fails. Failures are reported, never silently scored 0."""
    exe = shutil.which("DockQ")
    if exe is None:
        raise RuntimeError("DockQ is not on PATH; see environment.yml")
    for p in (predicted_pdb, experimental_pdb):
        if not Path(p).is_file():
            return None
    try:
        out = subprocess.check_output(
            [exe, str(predicted_pdb), str(experimental_pdb), "--mapping", mapping],
            stderr=subprocess.STDOUT, text=True)
    except subprocess.CalledProcessError:
        return None
    for line in out.splitlines():
        s = line.strip()
        if s.startswith("DockQ:"):
            try:
                return float(s.split()[1])
            except (IndexError, ValueError):
                return None
    return None


def quality(score):
    """Conventional DockQ bands, as drawn in panel D."""
    if score is None:
        return "failed"
    if score >= C.DOCKQ_MEDIUM:
        return "medium or better"
    if score >= C.DOCKQ_ACCEPTABLE:
        return "acceptable"
    return "incorrect"
