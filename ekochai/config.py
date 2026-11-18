"""Constants that define the method. 

"""

# --------------------------------------------------------------- antigen ----
# E2 ectodomain of HCV genotype 1a strain 1a154 (GenBank KJ187983), as modelled
# (chain A of the Chai-1 input FASTA). 301 residues. Residue NUMBERS throughout this file follow H77 polyprotein numbering.
E2_SEQUENCE = (
    "TPGAKQNIQLINTNGSWHINSTALNCNESLNTGWLAGLFYQHKFNSSGCPERLASCRRLT"
    "DFAQGWGPISYANGSGLDERPYCWHYPPRPCGIVPAKSVCGPVYCFTPSPVVVGTTDRSG"
    "APTYSWGANDTDVFVLNNTRPPLGNWFGCTWMNSTGFTKVCGAPPCVIGGVGNNTLLCPT"
    "DCFRKHPEATYSRCGSGPWITPRCMVNYPYRLWHYPCTINYTIFKVRMYVGGVEHRLEAA"
    "CNWTRGERCDLEDRDRSELSPLLLSTTQWQVLPCSFTTLPALSTGLIHLHQNIVDVQYLY"
    "G"
)


# E2_SEQUENCE[n - E2_REFERENCE_START - 1].
E2_REFERENCE_START = 404

# Two reference structures are numbered five residues out from the rest and are
NUMBERING_SHIFT = {"7mww": -5, "7mwx": -5}

# ------------------------------------------------------------- restraints ----
# An antigen residue joins the epitope if any of its Ca lies within this
# distance of any antibody Ca in the aligned reference complex. 
EPITOPE_CONTACT_ANGSTROM = 8.0

# Value written to `max_distance_angstrom` in the Chai-1 pocket restraint: how
# far the antibody chain may sit from the reference epitope. 
RESTRAINT_CUTOFFS = (8, 12, 16, 20, 24)
PUBLISHED_CUTOFF = 16

# ------------------------------------------------------------ Chai-1 run ----
# Arguments to chai_lab.chai1.run_inference, unchanged for every run.
CHAI_INFERENCE_KWARGS = dict(
    num_trunk_recycles=3,
    num_diffn_timesteps=200,
    num_diffn_samples=5,
    use_esm_embeddings=True,
    use_msa_server=True,
)
N_DIFFUSION_SAMPLES = 5        # models per run, = num_diffn_samples
SEEDS_PUBLISHED = tuple(range(1, 9))   # seeds 1-8 carry every condition

# ---------------------------------------------------------------- chains ----
# Chai-1 input FASTA order. 
CHAIN_IDS = {"antigen": "A", "heavy": "B", "light": "C"}

# Cutoff for the Ca-Ca contact counts that weight the two antibody interfaces
# in wiPTM. 
WIPTM_CONTACT_ANGSTROM = 8.0

# Before DockQ, heavy and light are merged into ONE chain so that DockQ scores
# a single antigen-Fv interface: the light chain is renamed C -> B and its
# residue numbers are offset by this amount. 
LIGHT_CHAIN_RESIDUE_OFFSET = 300
DOCKQ_MAPPING = "AB:AB"

# Conventional DockQ quality thresholds.
DOCKQ_ACCEPTABLE = 0.23
DOCKQ_MEDIUM = 0.49

# ------------------------------------------------- reference epitope bank ----
# Every reference complex that the restraint generator can build from. Each is
# "<pdb>_<antigenChain>_<heavyChain>_<lightChain>"; the aligned coordinates are
# in data/reference_structures/.
REFERENCE_EPITOPES = [
    "4mwf_C_H_L", "4web_E_H_L", "6bkb_E_H_L", "6bkc_E_H_L", "6bkd_E_H_L",
    "6meh_C_H_L", "6mei_C_H_L", "6mej_C_A_B", "6mej_C_H_L", "6mek_C_F_E",
    "6mek_C_H_L", "6urh_C_H_L", "6uyd_E_H_L", "6uyf_E_H_L", "6uyg_E_A_B",
    "6uyg_E_H_L", "6uym_E_H_L", "6wo3_E_H_L", "6wo4_E_H_L", "6wo5_E_A_B",
    "6wo5_E_H_L", "6woq_E_A_B", "6woq_E_J_K", "7jtf_E_A_B", "7jtg_E_A_B",
    "7mww_E_H_L", "7mwx_C_A_B", "7rfb_D_H_L", "7rfc_D_H_L", "7t6x_U_C_D",
    "7t6x_U_H_L", "8fsj_E_C_B", "8fsj_E_H_L", "8tgv_A_E_F", "8thz_A_E_F",
    "8txq_C_A_B", "8w0v_C_H_L", "8w0w_C_H_L", "8w0x_C_H_L", "8w0y_C_H_L",
    "hepc108_A_H_L"
]

# Which antigenic region each reference epitope belongs to. An antibody is only
# ever scanned against epitopes of its OWN region - see EPITOPE_SCAN below.
REGION_EPITOPES = {
    "FRLY": [
        "4mwf_C_H_L", "6bkb_E_H_L", "6bkc_E_H_L", "6bkd_E_H_L",
        "6meh_C_H_L", "6mei_C_H_L", "6mej_C_H_L", "6mek_C_H_L",
        "6urh_C_H_L", "6uyd_E_H_L", "6uyf_E_H_L", "6uyg_E_H_L",
        "6uym_E_H_L", "6wo3_E_H_L", "6wo4_E_H_L", "6wo5_E_H_L",
        "6woq_E_J_K", "7jtf_E_A_B", "7jtg_E_A_B", "7rfb_D_H_L",
        "7rfc_D_H_L", "7t6x_U_C_D", "8fsj_E_H_L", "8tgv_A_E_F",
        "8thz_A_E_F", "8w0v_C_H_L", "8w0w_C_H_L", "hepc108_A_H_L"
    ],
    "BSW": [
        "6mej_C_A_B", "6mek_C_F_E", "6uyg_E_A_B", "6wo5_E_A_B",
        "6woq_E_A_B", "8w0y_C_H_L"
    ],
    "BKLY": [
        "4web_E_H_L", "7mww_E_H_L", "7mwx_C_A_B", "7t6x_U_H_L",
        "8fsj_E_C_B", "8thz_A_H_L", "8w0x_C_H_L", "8txq_C_A_B"
    ],
}

# ------------------------------------------------------------- antibodies ----
# Region assignment read off the ELISA AUC panel (data/summary/elisa_auc.csv):
# each mAb is knocked out by exactly one Eko probe.
ANTIBODY_REGION = {
    "HEPC108": "FRLY", "hcab55": "FRLY", "hcab64": "FRLY",
    "hcab17": "BSW", "hcab44": "BSW",
    "hcab23": "BKLY", "hcab40": "BKLY",
}
ANTIBODIES = ["HEPC108", "hcab55", "hcab64", "hcab17", "hcab44",
              "hcab23", "hcab40"]

# Each antibody's own deposited structure, excluded from its own scan.
# hcab23 and hcab44 have no deposited structure in the reference bank.
SELF_STRUCTURE = {
    "HEPC108": "hepc108_A_H_L",
    "hcab55": "8w0v_C_H_L",
    "hcab64": "8w0w_C_H_L",
    "hcab40": "8w0x_C_H_L",
    "hcab17": "8w0y_C_H_L",
}

# ------------------------------------------------------- the epitope scan ----

#
# Two properties:
#   * REGION RESTRICTION. An antibody is scanned only against epitopes of its
#     own antigenic region. 
#   * SELF-EXCLUSION. An antibody's own deposited structure never appears in
#     its own scan.

EPITOPE_SCAN = {
    "HEPC108": [
        "4mwf_C_H_L", "6bkb_E_H_L", "6bkc_E_H_L", "6bkd_E_H_L",
        "6meh_C_H_L", "6mei_C_H_L", "6mej_C_H_L", "6mek_C_H_L",
        "6urh_C_H_L", "6uyd_E_H_L", "6uyf_E_H_L", "6uyg_E_H_L",
        "6uym_E_H_L", "6wo3_E_H_L", "6wo4_E_H_L", "6wo5_E_H_L",
        "6woq_E_J_K", "7jtf_E_A_B", "7jtg_E_A_B", "7rfb_D_H_L",
        "7rfc_D_H_L", "7t6x_U_C_D", "8fsj_E_H_L", "8tgv_A_E_F",
        "8thz_A_E_F", "8w0v_C_H_L", "8w0w_C_H_L"
    ],
    "hcab55": [
        "4mwf_C_H_L", "6bkb_E_H_L", "6bkc_E_H_L", "6bkd_E_H_L",
        "6meh_C_H_L", "6mei_C_H_L", "6mej_C_H_L", "6mek_C_H_L",
        "6urh_C_H_L", "6uyd_E_H_L", "6uyf_E_H_L", "6uyg_E_H_L",
        "6uym_E_H_L", "6wo3_E_H_L", "6wo4_E_H_L", "6wo5_E_H_L",
        "6woq_E_J_K", "7jtf_E_A_B", "7jtg_E_A_B", "7rfb_D_H_L",
        "7rfc_D_H_L", "7t6x_U_C_D", "8fsj_E_H_L", "8tgv_A_E_F",
        "8thz_A_E_F", "hepc108_A_H_L"
    ],
    "hcab64": [
        "4mwf_C_H_L", "6bkb_E_H_L", "6bkc_E_H_L", "6bkd_E_H_L",
        "6meh_C_H_L", "6mei_C_H_L", "6mej_C_H_L", "6mek_C_H_L",
        "6urh_C_H_L", "6uyd_E_H_L", "6uyf_E_H_L", "6uyg_E_H_L",
        "6uym_E_H_L", "6wo3_E_H_L", "6wo4_E_H_L", "6wo5_E_H_L",
        "6woq_E_J_K", "7jtf_E_A_B", "7jtg_E_A_B", "7rfb_D_H_L",
        "7rfc_D_H_L", "7t6x_U_C_D", "8fsj_E_H_L", "8tgv_A_E_F",
        "8thz_A_E_F", "hepc108_A_H_L"
    ],
    "hcab17": [
        "6mej_C_A_B", "6mek_C_F_E", "6uyg_E_A_B", "6wo5_E_A_B",
        "6woq_E_A_B"
    ],
    "hcab44": [
        "6mej_C_A_B", "6mek_C_F_E", "6uyg_E_A_B", "6wo5_E_A_B",
        "6woq_E_A_B", "8w0y_C_H_L"
    ],
    "hcab23": [
        "4web_E_H_L", "7mww_E_H_L", "7mwx_C_A_B", "7t6x_U_H_L",
        "8fsj_E_C_B", "8txq_C_A_B", "8w0x_C_H_L"
    ],
    "hcab40": [
        "4web_E_H_L", "7mww_E_H_L", "7mwx_C_A_B", "7t6x_U_H_L",
        "8fsj_E_C_B", "8txq_C_A_B"
    ],
}


def epitope_scan(antibody):
    """Reference epitopes scanned for one antibody (see EPITOPE_SCAN)."""
    return list(EPITOPE_SCAN[antibody])
