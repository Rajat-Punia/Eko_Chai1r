"""IMGT CDR-loop Ca RMSD, predicted against experimental.

The measurement is made in the ANTIGEN frame: the predicted model is
superposed on the experimental complex using the rigid E2 core only.
CDR residues are matched between the two structures by IMGT number.
"""

import numpy as np
from Bio.PDB import MMCIFParser, PDBParser, Superimposer
from Bio.Align import PairwiseAligner, substitution_matrices
from anarci import run_anarci

from . import config as C
from .files import open_structure, resolve

IMGT_CDR = {"1": (27, 38), "2": (56, 65), "3": (105, 117)}
CDR_ORDER = ["H1", "H2", "H3", "L1", "L2", "L3"]

# Restrict H3 to its central 8 residues -- the apex, which is also the part
# that contacts the epitope.
CENTRAL_N = {"H3": 8}

# E2 core refinement (fitting).
TRIM_CUTOFF, TRIM_ROUNDS, TRIM_MIN = 2.0, 5, 50

AA3 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
       "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
       "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
       "TYR": "Y", "VAL": "V", "MSE": "M"}

_ALIGNER = PairwiseAligner()
_ALIGNER.substitution_matrix = substitution_matrices.load("BLOSUM62")
_ALIGNER.open_gap_score, _ALIGNER.extend_gap_score = -11, -1
_ALIGNER.mode = "global"


def residues(chain):
    return [r for r in chain
            if r.id[0] == " " and "CA" in r and r.get_resname() in AA3]


def _seq(res_list):
    return "".join(AA3[r.get_resname()] for r in res_list)


def align_pairs(a, b):
    """Matched residue pairs from a global BLOSUM62 alignment."""
    al = _ALIGNER.align(_seq(a), _seq(b))[0]
    return [(a[i0 + k], b[j0 + k])
            for (i0, i1), (j0, _) in zip(al.aligned[0], al.aligned[1])
            for k in range(i1 - i0)]


def imgt_map(res_list):
    """residue -> (IMGT position, insertion code)."""
    _, num, ali, _ = run_anarci([("x", _seq(res_list))], scheme="imgt")
    if num[0] is None:
        raise RuntimeError("ANARCI found no antibody domain")
    out, i = {}, ali[0][0]["query_start"]
    for (pos, ins), aa in num[0][0][0]:
        if aa == "-":
            continue
        out[res_list[i]] = (pos, ins.strip())
        i += 1
    return out


def _cdr_of(pos):
    for name, (lo, hi) in IMGT_CDR.items():
        if lo <= pos <= hi:
            return name
    return None


def _rmsd(pairs):
    if not pairs:
        return float("nan"), 0
    d = np.array([p.coord - q.coord for p, q in pairs])
    return float(np.sqrt((d ** 2).sum(axis=1).mean())), len(pairs)


def superpose_on_e2_core(pred_model, exp_model, pred_e2, exp_e2):
    """Superpose the whole predicted model on the rigid E2 core.

    Returns (rms, n_core, n_aligned). Every predicted atom is moved, so all
    later RMSDs are in the antigen frame.
    """
    pairs = align_pairs(residues(pred_model[pred_e2]), residues(exp_model[exp_e2]))
    keep = pairs
    for _ in range(TRIM_ROUNDS):
        s = Superimposer()
        s.set_atoms([q["CA"] for _, q in keep], [p["CA"] for p, _ in keep])
        rot, tran = s.rotran
        d = [np.linalg.norm(np.dot(p["CA"].coord, rot) + tran - q["CA"].coord)
             for p, q in keep]
        nxt = [pr for pr, dd in zip(keep, d) if dd <= TRIM_CUTOFF]
        if len(nxt) < TRIM_MIN or len(nxt) == len(keep):
            break
        keep = nxt
    sup = Superimposer()
    sup.set_atoms([q["CA"] for _, q in keep], [p["CA"] for p, _ in keep])
    sup.apply([a for ch in pred_model for res in ch for a in res])
    return sup.rms, len(keep), len(pairs)


def cdr_rmsd(predicted_cif, experimental_pdb, exp_chains,
             pred_chains=("A", "B", "C")):
    """Per-loop Ca RMSD for one antibody.

    exp_chains  (antigen, heavy, light) chain ids of the experimental file
    pred_chains the same for the Chai-1 model; the CIF keeps three chains.

    Returns (rows, superposition) where rows is a list of dicts with
    antibody-agnostic keys cdr / rmsd / n / rmsd_framework.
    """
    p_e2, p_h, p_l = pred_chains
    x_e2, x_h, x_l = exp_chains
    with open_structure(resolve(predicted_cif)) as fh:
        pred = next(MMCIFParser(QUIET=True).get_structure("p", fh).get_models())
    with open_structure(resolve(experimental_pdb)) as fh:
        exp = next(PDBParser(QUIET=True).get_structure("e", fh).get_models())

    core_rms, n_core, n_aligned = superpose_on_e2_core(pred, exp, p_e2, x_e2)

    per_cdr, framework = {}, []
    for tag, pv, xv in (("H", residues(pred[p_h]), residues(exp[x_h])),
                        ("L", residues(pred[p_l]), residues(exp[x_l]))):
        pm, xm = imgt_map(pv), imgt_map(xv)
        inv = {v: k for k, v in xm.items()}
        for pres, key in pm.items():
            xres = inv.get(key)
            if xres is None:
                continue
            loop = _cdr_of(key[0])
            if loop:
                per_cdr.setdefault(tag + loop, []).append((pres["CA"], xres["CA"]))
            else:
                framework.append((pres["CA"], xres["CA"]))

    # imgt_map is built in sequence order, so each list is already N -> C
    for loop, n_keep in CENTRAL_N.items():
        pr = per_cdr.get(loop)
        if pr and len(pr) > n_keep:
            lo = (len(pr) - n_keep) // 2
            per_cdr[loop] = pr[lo:lo + n_keep]

    fw = Superimposer()
    fw.set_atoms([q for _, q in framework], [p for p, _ in framework])
    rot, tran = fw.rotran

    rows = []
    for loop in CDR_ORDER:
        pairs = per_cdr.get(loop, [])
        v, n = _rmsd(pairs)
        if pairs:
            d = np.array([(np.dot(p.coord, rot) + tran) - q.coord
                          for p, q in pairs])
            v_fw = float(np.sqrt((d ** 2).sum(axis=1).mean()))
        else:
            v_fw = float("nan")
        rows.append(dict(cdr=loop, rmsd=v, n=n, rmsd_framework=v_fw,
                         n_framework=n))
    sup = dict(core_rmsd=core_rms, n_core=n_core, n_aligned=n_aligned,
               framework_rmsd=fw.rms, n_framework_pairs=len(framework))
    return rows, sup
