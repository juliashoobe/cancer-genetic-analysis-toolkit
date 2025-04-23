"""Holds the core functions for the Cancer Genetic Analysis Toolkit."""

from collections.abc import Generator
from typing import Optional

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from scipy import stats

eqRecordGenerator = Generator[SeqRecord, None, None]


def parse_data(
    mutation_file: Optional[str] = None,
    expression_file: Optional[str] = None,
    fasta_file: Optional[str] = None,
) -> tuple[Optional[pd.DataFrame], Optional[pd.DataFrame], dict[str, str]]:
    """Parsing mutation, expression, and sequence data."""
    mutation_df = pd.read_csv(mutation_file) if mutation_file else None
    expression_df = pd.read_csv(expression_file) if expression_file else None

    sequences: dict[str, str] = {}
    if fasta_file:
        for record in SeqIO.parse(fasta_file, "fasta"):
            sequences[record.id] = str(record.seq)

    return mutation_df, expression_df, sequences


def mutation_detection(
    reference_seq: str, mutated_seq: str
) -> list[dict[str, str]]:
    """Compares a reference DNA sequence with a mutated target DNA sequence.

    Identifies position where mutations occur, indicating differences in the
    base pairs.
    """
    if len(reference_seq) != len(mutated_seq):
        raise ValueError("Sequences must be the same length.")

    mutations = []
    for i in range(len(reference_seq)):
        if reference_seq[i] != mutated_seq[i]:
            mutation = {
                "position": str(i + 1),  # Convert to string here
                "reference_base": reference_seq[i],
                "mutated_base": mutated_seq[i],
            }
            mutations.append(mutation)

    if not mutations:
        mutations.append({"message": "No mutations detected."})

    return mutations


def gene_expression(
    gene: str, expression_df: pd.DataFrame
) -> Optional[pd.Series]:
    """Analyzes data to retrieve expression values for specific gene.

    Returns expression data for the specified gene.
    """
    if gene in expression_df.columns:
        return expression_df[[gene]]
    else:
        print(f"Gene '{gene}' not found in expression dataset.")
        return None


def differential_expression(
    mutation_df: pd.DataFrame, expression_df: pd.DataFrame, gene: str
) -> Optional[dict[str, float]]:
    """Analyzes mutation-driven differential gene expression."""
    if gene not in expression_df.columns:
        print(f"Gene '{gene}' not found in expression dataset.")
        return None

    # Update to use 'label' for cancer vs non-cancer differentiation
    cancer_samples = mutation_df[mutation_df["label"] == "cancer"]
    non_cancer_samples = mutation_df[mutation_df["label"] == "non-cancer"]

    if len(cancer_samples) < 2 or len(non_cancer_samples) < 2:
        print(
            "Insufficient samples in one or both groups for statistical test."
        )
        return None

    cancer_expression = expression_df.loc[cancer_samples.index, gene].dropna()
    non_cancer_expression = expression_df.loc[
        non_cancer_samples.index, gene
    ].dropna()

    if len(cancer_expression) < 2 or len(non_cancer_expression) < 2:
        print("Not enough valid data after dropping NaNs")
        return None

    t_stat, p_value = stats.ttest_ind(cancer_expression, non_cancer_expression)

    if (
        pd.isna(t_stat)
        or pd.isna(p_value)
        or t_stat == float("inf")
        or p_value == float("inf")
    ):
        print("Invalid t-test result (NaN or infinite values).")
        return None

    return {"t_statistic": t_stat, "p_value": p_value}
