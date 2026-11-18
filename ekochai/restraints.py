"""Build Chai-1 pocket restraints from an aligned reference complex.

"""

from pathlib import Path

import numpy as np
import pandas as pd
from Bio.PDB import PDBParser, is_aa

from . import config as C
from .files import open_structure, resolve

_PARSER = PDBParser(QUIET=True)


def _ca_residues(structure, chain_id):
    out = {}
    model = next(iter(structure))
    if chain_id not in model:
        return out
    for res in model[chain_id]:
        if is_aa(res, standard=True) and "CA" in res:
            num, icode = res.get_id()[1], res.get_id()[2].strip()
            out[f"{num}{icode}"] = (num, res["CA"].coord)
    return out


def epitope_residues(reference_pdb, antigen_chain, heavy_chain, light_chain,
                     pdb_id=None, contact=C.EPITOPE_CONTACT_ANGSTROM):
    pdb_id = pdb_id or Path(reference_pdb).stem.split("_")[1]
    shift = C.NUMBERING_SHIFT.get(pdb_id.lower(), 0)

    with open_structure(resolve(reference_pdb)) as fh:
        st = _PARSER.get_structure("complex", fh)
    antigen = _ca_residues(st, antigen_chain)
    partners = [("B", _ca_residues(st, heavy_chain)),
                ("C", _ca_residues(st, light_chain))]

    hits = []
    for chain_label, ab in partners:
        for _, (resnum, ag_xyz) in antigen.items():
            idx = resnum + shift - C.E2_REFERENCE_START
            near = any(np.linalg.norm(xyz - ag_xyz) <= contact
                       for _, xyz in ab.values())
            if not near:
                continue
            aa = (C.E2_SEQUENCE[idx - 1] if 0 <= idx - 1 < len(C.E2_SEQUENCE)
                  else "X")
            hits.append((chain_label, f"{aa}{idx}"))
    return hits


def write_restraints(hits, cutoff, out_path):
    rows = [{
        "chainA": chain,
        "res_idxA": "",
        "chainB": C.CHAIN_IDS["antigen"],
        "res_idxB": label,
        "connection_type": "pocket",
        "confidence": 1.0,
        "min_distance_angstrom": 0.0,
        "max_distance_angstrom": float(cutoff),
        "comment": "antigen-antibody",
        "restraint_id": f"restraint_{i}",
    } for i, (chain, label) in enumerate(hits)]

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(out_path, index=False)
    return out_path


def build(reference_pdb, antigen_chain, heavy_chain, light_chain, cutoff,
          out_path, pdb_id=None):
    """epitope_residues + write_restraints in one call."""
    hits = epitope_residues(reference_pdb, antigen_chain, heavy_chain,
                            light_chain, pdb_id=pdb_id)
    return write_restraints(hits, cutoff, out_path), hits
