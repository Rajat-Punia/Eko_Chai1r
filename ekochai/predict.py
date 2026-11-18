"""Run Chai-1 on one antigen-antibody pair, with or without restraints.
"""

import logging
import shutil
from pathlib import Path

import numpy as np
import pandas as pd

from . import config as C


def run(fasta_path, output_dir, seed, restraints_path=None, overwrite=True,
        **kwargs):
    """Predict one complex and write Chai-1's scores alongside the models.

    fasta_path       3-chain FASTA: antigen, Ab_heavy, Ab_light
    restraints_path  pocket restraint file, or None for the unrestrained arm
    seed             the only thing varied between replicate runs

    Inference arguments are fixed in config.CHAI_INFERENCE_KWARGS; `kwargs`
    overrides them only if you mean to depart from the published protocol.

    Writes pred.model_idx_{0..4}.cif, scores.model_idx_{0..4}.npz and a
    flattened scores.csv, and returns the path to scores.csv.
    """
    from chai_lab.chai1 import run_inference      # imported late: needs a GPU

    output_dir = Path(output_dir)
    if output_dir.exists() and overwrite:
        logging.warning("removing old output directory: %s", output_dir)
        shutil.rmtree(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    args = dict(C.CHAI_INFERENCE_KWARGS)
    args.update(kwargs)
    if restraints_path is not None:
        args["constraint_path"] = Path(restraints_path)

    run_inference(fasta_file=Path(fasta_path), output_dir=output_dir,
                  seed=seed, **args)

    rows = []
    for idx in range(args["num_diffn_samples"]):
        npz = output_dir / f"scores.model_idx_{idx}.npz"
        if not npz.exists():
            logging.warning("score file not found: %s", npz)
            continue
        with np.load(npz) as s:
            rec = {k: (s[k].item() if s[k].size == 1 else s[k].tolist())
                   for k in s.files}
        rec["model_idx"] = idx
        rows.append(rec)

    out = output_dir / "scores.csv"
    pd.DataFrame(rows).to_csv(out, index=False)
    return out
