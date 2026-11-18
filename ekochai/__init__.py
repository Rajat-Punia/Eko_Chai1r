"""Restraint-guided prediction of HCV E2-antibody complex structures using Chai-1.

Modules
-------
config      every constant that defines the method
restraints  reference complex -> Chai-1 pocket restraints
predict     thin wrapper over chai_lab.chai1.run_inference
scoring     ipTM for a two-chain antibody, and model selection
structures  Chai-1 CIF -> the two-chain PDB DockQ is run on
dockq       DockQ against the experimental complex
cdr         IMGT CDR-loop Ca RMSD
"""

__version__ = "1.0.0"
