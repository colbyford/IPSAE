# cli.py
# Command-line interface for the ipsae package. The interface is
# backwards-compatible with the original ipsae.py script:
#
#   ipsae <path_to_pae_file> <path_to_pdb_or_cif_file> <pae_cutoff> <dist_cutoff>
#
# Derived from ipsae.py by Roland Dunbrack, Fox Chase Cancer Center.
# https://www.biorxiv.org/content/10.1101/2025.02.10.637595v2
# MIT license: script can be modified and redistributed for non-commercial and
# commercial use, as long as this information is reproduced.

import sys

from . import __version__
from .api import score_interactions

USAGE = """Usage for AF2 (PDB format):
   ipsae <path_to_pae_json_file> <path_to_pdb_file> <pae_cutoff> <dist_cutoff>
   ipsae RAF1_KSR1_scores_rank_001_alphafold2_multimer_v3_model_4_seed_003.json RAF1_KSR1_unrelaxed_rank_001_alphafold2_multimer_v3_model_4_seed_003.pdb 10 15

Usage for AF3 (mmCIF format):
   ipsae <path_to_pae_json_file> <path_to_mmcif_file> <pae_cutoff> <dist_cutoff>
   ipsae fold_aurka_tpx2_full_data_0.json  fold_aurka_tpx2_model_0.cif 10 15

Usage for Boltz (PDB or mmCIF format):
   ipsae <path_to_pae_npz_file> <path_to_mmcif_file> <pae_cutoff> <dist_cutoff>
   ipsae <path_to_pae_npz_file> <path_to_pdb_file> <pae_cutoff> <dist_cutoff>
   ipsae pae_AURKA_TPX2_model_0.npz  AURKA_TPX2_model_0.cif 10 15
   ipsae pae_AURKA_TPX2_model_0.npz  AURKA_TPX2_model_0.pdb 10 15

Options:
   --csv        also write the chain-pair scores as a CSV file
   --version    print the ipsae version and exit
   -h, --help   print this message and exit

All output files are written next to the cif/pdb file. The same interface is
available as "python -m ipsae" and through the legacy "python ipsae.py" script.
"""


def main(argv=None):
    """Console entry point; returns a process exit code."""
    if argv is None:
        argv = sys.argv[1:]

    positional = []
    write_csv = False
    for arg in argv:
        if arg in ("-h", "--help"):
            print(USAGE)
            return 0
        elif arg == "--version":
            print(f"ipsae {__version__}")
            return 0
        elif arg == "--csv":
            write_csv = True
        elif arg.startswith("--"):
            print(f"Unknown option: {arg}")
            print(USAGE)
            return 1
        else:
            positional.append(arg)

    if len(positional) < 4:
        print(USAGE)
        return 1

    pae_file, structure_file = positional[0], positional[1]
    try:
        pae_cutoff = float(positional[2])
        dist_cutoff = float(positional[3])
    except ValueError:
        print(f"PAE and distance cutoffs must be numbers; got: {positional[2]} {positional[3]}")
        return 1

    try:
        result = score_interactions(pae_file, structure_file, pae_cutoff, dist_cutoff)
    except (OSError, ValueError) as error:
        print(error)
        return 1

    result.write_outputs()
    if write_csv:
        result.to_csv()
    return 0


if __name__ == "__main__":
    sys.exit(main())
