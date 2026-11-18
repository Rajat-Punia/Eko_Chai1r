"""Convert a Chai-1 model to the two-chain form DockQ is run on.

"""

from pathlib import Path

from Bio.PDB import MMCIFParser, PDBIO

from . import config as C
from .files import open_structure, resolve

_PARSER = MMCIFParser(QUIET=True)


def cif_to_two_chain_pdb(cif_path, out_pdb):
    """Write `cif_path` as a two-chain PDB: A = antigen, B = heavy + light."""
    with open_structure(resolve(cif_path)) as fh:
        st = _PARSER.get_structure("model", fh)
    model = next(st.get_models())
    for chain in model:
        if chain.id == C.CHAIN_IDS["light"]:
            chain.id = C.CHAIN_IDS["heavy"]
            for res in chain:
                het, seqid, icode = res.id
                res.id = (het, seqid + C.LIGHT_CHAIN_RESIDUE_OFFSET, icode)
    out_pdb = Path(out_pdb)
    out_pdb.parent.mkdir(parents=True, exist_ok=True)
    io = PDBIO()
    io.set_structure(st)
    io.save(str(out_pdb))
    return out_pdb
