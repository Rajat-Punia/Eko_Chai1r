# AlphaFold 3 

## How it was run

One job per antibody, three protein chains (A = E2 ectodomain, H = heavy,
L = light), **1000 model seeds × 5 samples = 5000 models per antibody**.
`hcab44_input_template.json` is the input for one antibody, with the embedded
MSAs and templates stripped. Substitute the heavy and light sequences from
`../data/sequences/` for the other six antibodies.

```bash
python run_alphafold.py --json_path=hcab44_input_template.json --output_dir=hcab44/
```

Output is one directory per model, `seed-{N}_sample-{M}/`, plus a
`ranking_scores.csv` listing every model's AF3 ranking score.

## Scoring

`compute_dockq_af3.py` walks `ranking_scores.csv`, runs DockQ on each model
against the experimental complex, and writes `DockQ_{antibody}.csv`:

```bash
python compute_dockq_af3.py --antibody hcab44
```

Those files are concatenated into `../data/summary/af3_dockq_per_model.csv`
(7 antibodies * 1000 seeds * 5 samples = 35,000 rows).



## Notes

* "AF3 (1 seed)": for each seed, the best of its 5 samples by
  ranking score; mean ± SD over the 1000 seeds.
* "AF3 (1000 seeds)": the single best model by ranking score over all 5000.

