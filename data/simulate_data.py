"""Simulate mutation_data.csv, expression_data.csv, and sequences.fasta."""

import pandas as pd
import numpy as np
import random
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import os

np.random.seed(42)
random.seed(42)

# --- Settings ---
SAMPLES = 200  # Number of samples in the simulated data
GENES = [
    "BRCA1",
    "BRCA2",
    "TP53",
    "PIK3CA",
    "PTEN",
    "EGFR",
    "HER2",
    "AKT1",
    "CDH1",
    "BRAF",
]  # 10 most common genes associated with breast cancer

LABELS = ["cancer", "non-cancer"]
CANCER_RATES = 0.7  # 70% of samples will be cancer-related (30% non-cancer)

# Gene-specific mutation probabilities
MUTATION_PROBS = {
    "BRCA1": 0.55,
    "BRCA2": 0.25,
    "TP53": 0.40,
    "PIK3CA": 0.30,
    "PTEN": 0.25,
    "EGFR": 0.15,
    "HER2": 0.20,
    "AKT1": 0.10,
    "CDH1": 0.05,
    "BRAF": 0.05,
}


# --- 1. Mutation Data ---
# Create a dictionary to store mtuation data for each sample
mutation_data = {
    "sample_id": [f"S{i + 1}" for i in range(SAMPLES)],
    "label": np.random.choice(
        LABELS, size=SAMPLES, p=[CANCER_RATES, 1 - CANCER_RATES]
    ),  # Randomly assign cancer (70%) and non-cancer (30%) labels to samples
}

# For each gene, generate mutation data based on the label
for gene in GENES:
    probs = []
    for label in mutation_data["label"]:
        if label == "cancer":
            p = MUTATION_PROBS[gene]
        else:
            p = MUTATION_PROBS[gene] * 0.1
        probs.append(np.random.binomial(1, p))
    mutation_data[gene] = probs

mutation_df = pd.DataFrame(mutation_data)
mutation_df.to_csv("data/mutation_data.csv", index=False)


# --- 2. Gene Expression Data ---
expression_data = {
    "sample_id": mutation_df["sample_id"],
    "label": mutation_df["label"],
}

# For each gene, generate expression data based on mutation and cancer status
for gene in GENES:
    expression = []
    for idx, row in mutation_df.iterrows():
        label = row["label"]
        mutated = row[gene]

        base_mean = 8
        noise = np.random.normal(0, 1)

        if label == "cancer":
            if mutated:
                if gene in ["BRCA1", "BRCA2", "TP53", "PTEN"]:
                    shift = np.random.normal(-2, 1)
                else:
                    shift = np.random.normal(2, 1)
            else:
                shift = np.random.normal(1, 0.5)
        else:
            shift = 0

        expression.append(base_mean + shift + noise)

    expression_data[gene] = expression


expression_df = pd.DataFrame(expression_data)
expression_df.to_csv("data/expression_data.csv", index=False)


# --- 3. DNA Sequences ---
def simulate_dna_sequence(length: int = 100) -> str:
    """Generate a random DNA sequence of the given length."""
    return "".join(random.choices("ACGT", k=length))


reference_sequences = {gene: simulate_dna_sequence() for gene in GENES}
records = []

for gene, seq in reference_sequences.items():
    records.append(
        SeqRecord(Seq(seq), id=f"{gene}_REF", description="Reference")
    )

for gene, ref_seq in reference_sequences.items():
    for i in range(5):
        seq_list = list(ref_seq)
        for _ in range(4):
            idx = random.randint(0, len(seq_list) - 1)
            orig = seq_list[idx]
            options = ["A", "C", "G", "T"]
            options.remove(orig)
            seq_list[idx] = random.choice(options)
        records.append(
            SeqRecord(
                Seq("".join(seq_list)),
                id=f"{gene}_VAR{i + 1}",
                description="Mutated",
            )
        )

SeqIO.write(records, "data/sequences.fasta", "fasta")
