"""Open a structure file 
"""

import gzip
from contextlib import contextmanager
from pathlib import Path


@contextmanager
def open_structure(path):
    """Text handle for a .cif/.pdb file, transparently gunzipping a .gz."""
    path = Path(path)
    if path.suffix == ".gz":
        fh = gzip.open(path, "rt")
    else:
        fh = open(path, "r")
    try:
        yield fh
    finally:
        fh.close()


def resolve(path):
    """The given path, or its .gz sibling if only that exists."""
    path = Path(path)
    if path.exists():
        return path
    gz = path.with_suffix(path.suffix + ".gz")
    if gz.exists():
        return gz
    return path
