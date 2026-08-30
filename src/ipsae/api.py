# api.py
# High-level programmatic interface for the ipsae package.
#
# Derived from ipsae.py by Roland Dunbrack, Fox Chase Cancer Center.
# https://www.biorxiv.org/content/10.1101/2025.02.10.637595v2
# MIT license: script can be modified and redistributed for non-commercial and
# commercial use, as long as this information is reproduced.

from dataclasses import dataclass, field
from typing import Dict, List, Optional

from .confidence import ConfidenceData, load_confidence
from .outputs import (
    build_chain_pair_groups,
    build_residue_records,
    write_byres_file,
    write_csv_file,
    write_pair_scores_file,
    write_pml_file,
)
from .parsers import Structure, load_structure
from .scoring import ScoreResults, compute_scores


def detect_model_type(pae_file, structure_file):
    """Detect the prediction software from the file extensions.

    Returns (model_type, file_format) where model_type is 'af2', 'af3', or
    'boltz' and file_format is 'pdb' or 'cif'. Gzip-compressed JSON PAE files
    (``.json.gz``) and AF2 pickle files (``.pkl``) are also recognized.
    """
    pae_is_json = pae_file.endswith(".json") or pae_file.endswith(".json.gz")

    if ".pdb" in structure_file and (pae_is_json or pae_file.endswith(".pkl")):
        return "af2", "pdb"
    elif ".cif" in structure_file and pae_is_json:
        return "af3", "cif"
    elif ".cif" in structure_file and pae_file.endswith(".npz"):  # Boltz1/2 in cif format
        return "boltz", "cif"
    elif ".pdb" in structure_file and pae_file.endswith(".npz"):  # Boltz1/2 in pdb format
        return "boltz", "pdb"
    raise ValueError(f"Wrong PDB or PAE file type: {structure_file} / {pae_file}")


def cutoff_string(cutoff):
    """Format a cutoff for file names and output columns (zero-padded below 10)."""
    string = str(int(cutoff))
    if cutoff < 10:
        string = "0" + string
    return string


@dataclass
class IPSAEResult:
    """Scores for one model, with programmatic access and file writers."""

    structure: Structure
    confidence: ConfidenceData
    scores: ScoreResults
    pae_file: str
    structure_file: str
    pae_cutoff: float
    dist_cutoff: float
    model_type: str
    pdb_stem: str
    path_stem: str
    pae_string: str
    dist_string: str
    _groups: Optional[List] = field(default=None, repr=False)

    def _chain_pair_groups(self):
        if self._groups is None:
            self._groups = build_chain_pair_groups(
                self.structure, self.confidence, self.scores,
                self.pae_string, self.dist_string, self.pdb_stem)
        return self._groups

    @property
    def chain_pairs(self) -> List[Dict]:
        """Chain-pair score records (dicts), one per output row (asym A->B, asym B->A, max)."""
        records = []
        for group in self._chain_pair_groups():
            for rec in group:
                records.append({key: value for key, value in rec.items() if not key.startswith("_")})
        return records

    @property
    def residues(self) -> List[Dict]:
        """By-residue score records (dicts), in the order of the _byres.txt file."""
        return build_residue_records(self.structure, self.confidence, self.scores)

    def get_score(self, chain1, chain2, metric="ipSAE", score_type="max"):
        """Return a single score for a chain pair.

        ``metric`` is one of 'ipSAE', 'ipSAE_d0chn', 'ipSAE_d0dom', 'ipTM_af',
        'ipTM_d0chn', 'pDockQ', 'pDockQ2', or 'LIS'; ``score_type`` is 'asym'
        (chain1 aligned, chain2 scored) or 'max' (maximum over both directions).
        """
        for rec in self.chain_pairs:
            if rec["Type"] != score_type:
                continue
            if score_type == "max" and {rec["Chn1"], rec["Chn2"]} != {str(chain1), str(chain2)}:
                continue
            if score_type == "asym" and (rec["Chn1"], rec["Chn2"]) != (str(chain1), str(chain2)):
                continue
            if metric not in rec:
                raise KeyError(f"Unknown metric: {metric}; choose from "
                               "ipSAE, ipSAE_d0chn, ipSAE_d0dom, ipTM_af, ipTM_d0chn, pDockQ, pDockQ2, LIS")
            return rec[metric]
        raise KeyError(f"No {score_type} scores for chain pair {chain1}/{chain2}")

    def write_outputs(self, output_stem=None, byres=True, pml=True):
        """Write the chain-pair score file (.txt), by-residue file (_byres.txt),
        and PyMOL script (.pml). Returns a dict of the written file paths."""
        stem = output_stem if output_stem is not None else self.path_stem
        groups = self._chain_pair_groups()
        paths = {"scores": stem + ".txt"}
        write_pair_scores_file(groups, paths["scores"])
        if pml:
            paths["pml"] = stem + ".pml"
            write_pml_file(groups, paths["pml"])
        if byres:
            paths["byres"] = stem + "_byres.txt"
            write_byres_file(self.structure, self.confidence, self.scores, paths["byres"])
        return paths

    def to_csv(self, csv_path=None):
        """Write the chain-pair scores to a CSV file and return its path."""
        if csv_path is None:
            csv_path = self.path_stem + ".csv"
        write_csv_file(self._chain_pair_groups(), csv_path)
        return csv_path


def score_interactions(pae_file, structure_file, pae_cutoff=10.0, dist_cutoff=10.0, model_type=None):
    """Score all pairwise chain-chain interactions of an AF2/AF3/Boltz model.

    Parameters
    ----------
    pae_file : str
        PAE file: AF2 scores ``.json``/``.json.gz``/``.pkl``, AF3
        full-data/confidences ``.json``/``.json.gz``, or Boltz ``pae_*.npz``.
    structure_file : str
        Model coordinates: ``.pdb`` (AF2/Boltz) or ``.cif`` (AF3/Boltz).
    pae_cutoff : float
        PAE cutoff (Angstroms) for the ipSAE aligned-residue pair selection.
    dist_cutoff : float
        CA-CA distance cutoff (Angstroms) for interface residue counts.
    model_type : str, optional
        'af2', 'af3', or 'boltz'; detected from the file extensions if omitted.

    Returns
    -------
    IPSAEResult
    """
    if model_type is None:
        model_type, file_format = detect_model_type(pae_file, structure_file)
    else:
        _, file_format = detect_model_type(pae_file, structure_file)

    pae_cutoff = float(pae_cutoff)
    dist_cutoff = float(dist_cutoff)
    pae_string = cutoff_string(pae_cutoff)
    dist_string = cutoff_string(dist_cutoff)

    if file_format == "pdb":
        pdb_stem = structure_file.replace(".pdb", "")
    else:
        pdb_stem = structure_file.replace(".cif", "")
    path_stem = f'{pdb_stem}_{pae_string}_{dist_string}'

    structure = load_structure(structure_file, file_format)
    confidence = load_confidence(pae_file, structure, model_type)
    scores = compute_scores(structure, confidence, pae_cutoff, dist_cutoff)

    return IPSAEResult(
        structure=structure,
        confidence=confidence,
        scores=scores,
        pae_file=pae_file,
        structure_file=structure_file,
        pae_cutoff=pae_cutoff,
        dist_cutoff=dist_cutoff,
        model_type=model_type,
        pdb_stem=pdb_stem,
        path_stem=path_stem,
        pae_string=pae_string,
        dist_string=dist_string,
    )
