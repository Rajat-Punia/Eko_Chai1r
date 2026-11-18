#!/usr/bin/env python3
import argparse
import subprocess
import sys
from pathlib import Path

import pandas as pd
from tqdm import tqdm

try:
    import gemmi
except ImportError:
    sys.stderr.write("ERROR: Python package 'gemmi' is required.\n"
                     "Install with: pip install gemmi\n")
    raise


def convert_cif_to_pdb_merge_AB(cif_path: Path, pdb_out: Path) -> None:
    """
    Load CIF, create a new structure where:
      - chain A stays A
      - chain H and chain L are merged into a single chain B
      - residues in L are renumbered by +300 (insertion code preserved)
    Save as PDB to pdb_out.
    """
    if not cif_path.exists():
        raise FileNotFoundError(f"CIF not found: {cif_path}")

    st = gemmi.read_structure(str(cif_path))
    if len(st) == 0:
        raise ValueError(f"No models in structure: {cif_path}")

    model = st[0]

    # Capture H and L before removing anything
    chain_H = model.find_chain("H")
    chain_L = model.find_chain("L")

    # If a chain B already exists, remove by NAME (gemmi expects str)
    try:
        model.remove_chain("B")
    except Exception:
        pass

    new_B = gemmi.Chain("B")

    def add_chain_residues_into(chain_obj, shift_l=False):
        if chain_obj is None:
            return
        for res in chain_obj:
            r = res.clone()
            if shift_l:
                sid = r.seqid
                sid.num = sid.num + 300  # keep insertion code (e.g., 27A)
                r.seqid = sid
            new_B.add_residue(r)

    add_chain_residues_into(chain_H, shift_l=False)
    add_chain_residues_into(chain_L, shift_l=True)

    # Remove original H and L by NAME (avoid TypeError)
    try:
        model.remove_chain("H")
    except Exception:
        pass
    try:
        model.remove_chain("L")
    except Exception:
        pass

    # Add merged B
    model.add_chain(new_B)

    pdb_out.parent.mkdir(parents=True, exist_ok=True)
    st.write_minimal_pdb(str(pdb_out))


def run_dockq(model_pdb: Path, exp_pdb: Path) -> float:
    cmd = ["DockQ", str(model_pdb), str(exp_pdb), "--mapping", "AB:AB"]
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
    except FileNotFoundError:
        raise RuntimeError("DockQ not found in PATH.")
    out = (proc.stdout or "") + "\n" + (proc.stderr or "")

    score = float("nan")
    for line in out.splitlines():
        if "DockQ:" in line:
            parts = line.strip().split()
            try:
                idx = parts.index("DockQ:") + 1
                score = float(parts[idx])
            except Exception:
                try:
                    score = float(parts[-1])
                except Exception:
                    pass
    return score


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--antibody", required=True)
    parser.add_argument("--base_dir", default=".")
    parser.add_argument("--csv_name", default="ranking_scores.csv")
    args = parser.parse_args()

    base = Path(args.base_dir).resolve()
    antibody = args.antibody

    ranking_csv = base / antibody / args.csv_name
    exp_pdb = base / "exp_structs_Chai1r" / f"{antibody}.pdb"
    out_csv = base / f"DockQ_{antibody}.csv"

    if not ranking_csv.exists():
        sys.stderr.write(f"ERROR: ranking CSV not found: {ranking_csv}\n")
        sys.exit(1)
    if not exp_pdb.exists():
        sys.stderr.write(f"ERROR: experimental PDB not found: {exp_pdb}\n")
        sys.exit(1)

    df = pd.read_csv(ranking_csv)
    required_cols = {"seed", "sample", "ranking_score"}
    if not required_cols.issubset(df.columns):
        raise ValueError(f"CSV must contain columns: {required_cols}. Found: {set(df.columns)}")

    dockq_scores = []
    for _, row in tqdm(df.iterrows(), total=len(df), desc=f"Processing {antibody}"):
        seed = int(row["seed"])
        sample = int(row["sample"])

        cif_path = base / antibody / f"seed-{seed}_sample-{sample}" / "model.cif"
        pdb_path = base / antibody / f"seed-{seed}_sample-{sample}" / "model.pdb"

        try:
            convert_cif_to_pdb_merge_AB(cif_path, pdb_path)
            dockq = run_dockq(pdb_path, exp_pdb)
        except Exception as e:
            sys.stderr.write(f"[WARN] seed {seed} sample {sample}: {e}\n")
            dockq = float("nan")

        dockq_scores.append(dockq)

    df_out = df.copy()
    df_out["dockq"] = dockq_scores
    df_out.to_csv(out_csv, index=False)
    print(f"Wrote: {out_csv}")


if __name__ == "__main__":
    main()

