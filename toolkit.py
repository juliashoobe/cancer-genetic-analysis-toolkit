"""Holds the core functions for the Cancer Genetic Analysis Toolkit."""

from collections.abc import Generator
from typing import Optional

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

eqRecordGenerator = Generator[SeqRecord, None, None]


def parse_sequences(fasta_file: str) -> list[SeqRecord]:
    """Returns a list of SeqRecord objects from a FASTA file."""
    return list(SeqIO.parse(fasta_file, "fasta"))


def parse_data(
    mutation_file: str, expression_file: str, fasta_file: str
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, str]]:
    """Parses mutation, expression, and sequence data from files."""
    mutation_df = pd.read_csv(mutation_file)
    expression_df = pd.read_csv(expression_file)

    records = parse_sequences(fasta_file)

    sequences: dict[str, str] = {
        record.id: str(record.seq)
        for record in records
        if record.id is not None
    }

    return mutation_df, expression_df, sequences


def mutation_detection(
    reference_seq: str, mutated_seq: str
) -> list[dict[str, object]]:
    """Compares a reference DNA sequence with a mutated target DNA sequence.

    Identifies position where mutations occur, indicating differences in the
    base pairs.
    """
    mutations = []
    for i in range(len(reference_seq)):
        if reference_seq[i] != mutated_seq[i]:
            mutation = {
                "position": i + 1,
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


mutation_file = "data/mutation_data.csv"
expression_file = "data/expression_data.csv"
fasta_file = "data/sequences.fasta"

mutation_df = pd.DataFrame(
    {
        "gene": ["BRCA1", "TP53", "EGFR"],
        "mutation_type": ["missense", "nonsense", "frameshift"],
        "position": [100, 150, 200],
    }
)

expression_df = pd.DataFrame(
    {
        "BRCA1": [5.2, 3.1, 4.5],
        "TP53": [1.2, 2.3, 2.4],
        "EGFR": [7.5, 6.2, 8.1],
    },
    index=["sample1", "sample2", "sample3"],
)

# Sample FASTA sequence (normally, this would be loaded from the FASTA file)
sequences = {
    "BRCA1_REF": "GACAGGTACAAGAAGGAGTATGCATCAATGTGGTCGTGTGGAACAAACGCCACTGGAG"
    "ACTGGGTTAACCATTCGCTCCAGCGTCA",
    "BRCA1_VAR1": "GACAGGTACAAGAAGGAGTATGCATCAATGTGGTCGTGTGGAAC"
    "AAACGCCACTGGAGACTGGGTTAACCATTCGCTCCAGCGTCA",
}

mutations = mutation_detection(sequences["BRCA1_REF"], sequences["BRCA1_VAR1"])
print("Detected Mutations:")
print(mutations)

mutation_df, expression_df, sequences = parse_data(
    mutation_file, expression_file, fasta_file
)

# Test gene expression retrieval
gene = "BRCA1"
expression_data = gene_expression(gene, expression_df)
if expression_data is not None:
    print(f"Expression data for {gene}:")
    print(expression_data)
