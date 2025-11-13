"""
Useful functions for peptide pre-processing.
"""

from io import StringIO
import os

import numpy as np
import pandas as pd
import requests as r
from tqdm import tqdm

from Bio import SeqIO
from structuremap.processing import format_alphafold_data, annotate_accessibility, get_smooth_score


def get_complete_sequence(proteinID: str) -> str:
    """
    Retrieves the UniProt amino acid sequence for a protein.

    TODO: handle invalid protein sequence response
    """

    url = "http://www.uniprot.org/uniprot/" + proteinID + ".fasta"
    response = r.post(url)
    cData = "".join(response.text)

    Seq = StringIO(cData)
    pSeq = list(SeqIO.parse(Seq, "fasta"))
    return str(pSeq[0].seq)

def create_sequence_cache(sequence_cache_path: os.path, verbose: bool=True) -> None:
    """
    Create an empty DataFrame mapping prtein UniProt IDs to
    full amino acid sequences, only if the cache doesn't already
    exist.
    """

    if os.path.isfile(sequence_cache_path):
        if verbose: print("Protein sequence cache already exists!")
        return

    empty_df = pd.DataFrame(columns=["Protein ID", "Complete Sequence"])
    empty_df.to_csv(sequence_cache_path)
    if verbose: print("Created new protein sequence cache.")
    return

def update_sequence_cache(sequence_cache_path: os.path, uniprotIDs: list[str], verbose: bool=True) -> None:
    """
    Update a DataFrame (UniProt ID -> AA sequence mapping) with the unknown sequences
    for the given uniprotIDs.

    TODO: speed this up with multiple threads?
    """

    # Load cache of mapping from UniProtIDs -> complete amino acid sequence
    sequence_cache_df = pd.read_csv(sequence_cache_path).set_index("Unnamed: 0")
    sequence_cache_df.index.name = None
    if verbose: print(sequence_cache_df)

    # Determine unknown proteins
    #unknown_uniprotIDs = uniprotIDs[~np.isin(uniprotIDs, sequence_cache_df["Protein ID"].values)]
    unknown_uniprotIDs = np.setdiff1d(uniprotIDs, sequence_cache_df["Protein ID"].values)
    n_unknown = len(unknown_uniprotIDs)
    if n_unknown == 0:
        print("All proteins have a known sequence!")
        return

    unknown_sequences_df = pd.DataFrame({"Protein ID": unknown_uniprotIDs})
    print(f"{n_unknown} unknown sequence(s) to retrieve.")

    # Retrieve sequences for unknown proteins
    tqdm.pandas()
    unknown_sequences_df["Complete Sequence"] = unknown_sequences_df["Protein ID"].progress_apply(get_complete_sequence)
    if verbose: print(unknown_sequences_df)

    # Update cache with newly retrieved sequences
    sequence_cache_df_updated = pd.concat([sequence_cache_df, unknown_sequences_df])
    sequence_cache_df_updated.to_csv(sequence_cache_path)
    print(f"Updated {sequence_cache_path} with {n_unknown} new protein sequence(s).")
    return

def create_modifications_pattern(amino_acid: str, modifications: list[str]) -> str:
    """
    Create a regex string to match a set of amino acid modifications.
    """

    whole, mantissa = modifications[0].split(".")
    pattern = r"{}\[{}\.{}\]".format(amino_acid, whole, mantissa)

    for i in range(1, len(modifications)):
        whole, mantissa = modifications[i].split(".")
        pattern += r"|{}\[{}\.{}\]".format(amino_acid, whole, mantissa)

    return pattern

def filter_amino_acid_sequence(string: str) -> str:
    """
    Filters a string down to only IUPAC amino acid characters.
    """

    IUPACCodes = "ACDEFGHIKLMNPQRSTVWY"
    return "".join([c for c in string if c in IUPACCodes])

def extract_sites(peptides: pd.DataFrame, amino_acid: str="M", amino_acid_str: str="Met", analysis_threshold: int=20, modifications: list[str]=["649.3660", "655.3735"]) -> pd.DataFrame:
    """
    Process a standard isoTOP-ABPP dataset to extract peptide sites
    and sequences (left and right) adjacent to the amino acid of interest.
    The resulting amino acid location column (e.g. "Met Location") will be
    zero-indexed, while the resulting site column will be one-indexed.

    Assumes columns (in peptides):
     - Protein ID: the UniProt ID corresponding to the peptide's full protein
     - Peptide Sequence: the peptide sequence (no modifications)
     - Light Modified Peptide: the peptide sequence (with modification)
     - Complete Sequence: the full-length protein sequence for each corresponding peptide
    """

    # Find peptide in full protein sequence
    peptides["Peptide Location"] = pd.Series([a.find(b) for a, b in zip(peptides["Complete Sequence"], peptides["Peptide Sequence"])])
    peptides["Peptide Length"] = peptides["Peptide Sequence"].str.len()

    # Create regex pattern to identify desired modifications
    modifications_pattern = create_modifications_pattern(amino_acid, modifications)
    left_prefix_pattern = "(.*)(" + modifications_pattern + ")"

    # Extract left prefix of modified target site (for subsequent indexing)
    peptides["Left Prefix"] = peptides["Light Modified Peptide"].str.extract(left_prefix_pattern)[0]
    peptides["Left Prefix"] = peptides["Left Prefix"].map(filter_amino_acid_sequence)
    peptides["Left Prefix Length"] = peptides["Left Prefix"].str.len()

    # Find target site in full protein sequence (zero-indexed)
    peptides[f"{amino_acid_str} Location"] = peptides["Peptide Location"] + peptides["Left Prefix Length"]

    # Annotate target site (one-indexed)
    peptides["Site"] = peptides["Protein ID"] + "_M" + (peptides[f"{amino_acid_str} Location"] + 1).astype(str)

    # Extract sequences (left and right) adjacent to Met site, based on given analysis threshold
    peptides[f"Left {analysis_threshold}"] = [A[B-analysis_threshold:B] if (B - analysis_threshold >= 0) else A[0:B-1] for A, B in zip(peptides["Complete Sequence"], peptides[f"{amino_acid_str} Location"])]
    peptides[f"Right {analysis_threshold}"] = [A[B+1:B+1+analysis_threshold] for A, B in zip(peptides["Complete Sequence"], peptides[f"{amino_acid_str} Location"])]

    return peptides

def calculate_accessibilities(cif_dir: os.path, pae_dir: os.path, uniprotIDs: list[str], radii: list[int], angles: list[int]=[180], smooth: bool=True) -> pd.DataFrame:
    """
    Calculate the amino acid crowdedness (pPSE) of each residue site at various radii and angles.
    """

    assert len(uniprotIDs) > 0
    assert len(radii) > 0
    assert len(angles) > 0 

    # Format AlphaFold protein structure data into DataFrame
    alphafold_annotation = format_alphafold_data(
        directory=cif_dir, 
        protein_ids=uniprotIDs
    )

    # Calculate accessibilities at various radii and angles
    accessibilities = alphafold_annotation.copy()

    for radius in radii: 
        for angle in angles:

            exposure = annotate_accessibility(
                df=alphafold_annotation, 
                max_dist=radius, 
                max_angle=angle, 
                error_dir=pae_dir
            )

            accessibilities = accessibilities.merge(
                exposure, how="left", on=["protein_id", "AA", "position"]
            )

    if smooth:
        cols_to_smooth = [col for col in accessibilities.columns if "_pae" in col]
        accessibilities = get_smooth_score(accessibilities, cols_to_smooth, [10])

    if "nAA_24_180_pae" not in accessibilities.columns:
        exposure = annotate_accessibility(
            df=alphafold_annotation, 
            max_dist=24, 
            max_angle=180, 
            error_dir=pae_dir
        )
        accessibilities = accessibilities.merge(
            exposure, how="left", on=["protein_id", "AA", "position"]
        )
    if "nAA_24_180_pae_smooth10" not in accessibilities.columns:
        accessibilities = get_smooth_score(accessibilities, ["nAA_24_180_pae"], [10])

    accessibilities["IDR"] = np.where(accessibilities["nAA_24_180_pae_smooth10"] <= 34.27, 1, 0)

    return accessibilities
