"""Holds the core functions for the Cancer Genetic Analysis Toolkit."""

from collections.abc import Generator
from typing import Optional

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

eqRecordGenerator = Generator[SeqRecord, None, None]


def parse_data(
    mutation_file: str, expression_file: str, fasta_file: str
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, str]]:
    """Parsing mutation, expression, and sequence data."""
    mutation_df = pd.read_csv(mutation_file)
    expression_df = pd.read_csv(expression_file)

    # Parse sequences directly within the function
    records = list(SeqIO.parse(fasta_file, "fasta"))

    sequences: dict[str, str] = {
        record.id: str(record.seq)
        for record in records
        if record.id is not None
    }

    return mutation_df, expression_df, sequences


def mutation_detection(
    reference_seq: str, mutated_seq: str
) -> list[dict[str, str]]:
    """Compares a reference DNA sequence with a mutated target DNA sequence.

    Identifies position where mutations occur, indicating differences in the
    base pairs.
    """
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
