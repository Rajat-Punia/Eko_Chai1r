"""Contact-weighted ipTM score.

Chai-1 reports `per_chain_pair_iptm`, a square matrix over the input chains
[A, B, C] = [antigen, heavy, light], so the two
interfaces are weighted by how much contact they actually make.
"""

import ast

import numpy as np
from Bio.PDB import MMCIFParser

from . import config as C
from .files import open_structure, resolve

_PARSER = MMCIFParser(QUIET=True)


def chain_pair_entries(row):
    """The four antigen<->antibody entries of one model's ipTM matrix.

    `row` is a record from Chai-1's scores.csv; its `per_chain_pair_iptm`
    field is a stringified nested list whose first element is the matrix.
    """
    m = ast.literal_eval(row["per_chain_pair_iptm"])[0]
    return m[0][1], m[0][2], m[1][0], m[2][0]   # A>B, A>C, B>A, C>A


def iptm(row):
    """Unweighted ipTM: the plain mean of the four entries. NOT the published
    metric; provided because the published data carries both columns."""
    return float(np.mean(chain_pair_entries(row)))


def ca_coords(cif_path, chain_id):
    """Ca coordinates of one chain of a Chai-1 model."""
    with open_structure(resolve(cif_path)) as fh:
        st = _PARSER.get_structure("model", fh)
    model = next(st.get_models())
    xyz = [r["CA"].coord for ch in model if ch.id == chain_id
           for r in ch if "CA" in r]
    return np.vstack(xyz).astype(np.float32) if xyz else np.zeros((0, 3), np.float32)


def count_contacts(a, b, cutoff=C.WIPTM_CONTACT_ANGSTROM):
    """Number of Ca pairs within `cutoff`, one from each chain."""
    if len(a) == 0 or len(b) == 0:
        return 0
    d2 = np.sum((a[:, None, :] - b[None, :, :]) ** 2, axis=-1)
    return int(np.sum(d2 <= cutoff ** 2))


def wiptm(row, cif_path):
    """Contact-weighted ipTM -- the metric plotted in panels B and D.

        sAB = (A>B + B>A)/2            antigen-heavy, direction-averaged
        sAC = (A>C + C>A)/2            antigen-light
        wAB, wAC                       Ca-Ca pairs within 8 A at each interface
        wiPTM = (wAB*sAB + wAC*sAC) / (wAB + wAC)


    Returns (wiptm, wAB, wAC).
    """
    a_b, a_c, b_a, c_a = chain_pair_entries(row)
    s_ab, s_ac = 0.5 * (a_b + b_a), 0.5 * (a_c + c_a)
    try:
        ag = ca_coords(cif_path, C.CHAIN_IDS["antigen"])
        hv = ca_coords(cif_path, C.CHAIN_IDS["heavy"])
        lt = ca_coords(cif_path, C.CHAIN_IDS["light"])
    except Exception:
        return 0.0, 0, 0
    w_ab, w_ac = count_contacts(ag, hv), count_contacts(ag, lt)
    total = w_ab + w_ac
    if total == 0:
        return 0.0, w_ab, w_ac
    return float((w_ab * s_ab + w_ac * s_ac) / total), w_ab, w_ac


def select_best(frame, metric="wiptm"):
    return frame.loc[frame[metric].idxmax()]
