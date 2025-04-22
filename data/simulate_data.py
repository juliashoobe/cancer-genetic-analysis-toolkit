"""Simulate mutation_data.csv, expression_data.csv, and sequences.fasta."""

import random

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# --- Set random seeds for reproducibility ---
np.random.seed(42)
random.seed(42)

# --- Settings ---
SAMPLES = 200
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
]
LABELS = ["cancer", "non-cancer"]
CANCER_RATES = 0.7  # 70% cancer, 30% non-cancer

# Realistic mutation probabilities based on research
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
sample_ids = [f"S{i + 1}" for i in range(SAMPLES)]
labels = np.random.choice(
    LABELS, size=SAMPLES, p=[CANCER_RATES, 1 - CANCER_RATES]
)

# Create DataFrame with sample_id and label
mutation_df = pd.DataFrame({"sample_id": sample_ids, "label": labels})

# Add gene mutation columns
for gene in GENES:
    mutation_df[gene] = mutation_df["label"].apply(
        lambda label, g=gene: np.random.binomial(
            1,
            MUTATION_PROBS[g]
            if label == "cancer"
            else MUTATION_PROBS[g] * 0.1,
        )
    )


# Save mutation data
mutation_df.to_csv("data/mutation_data.csv", index=False)

# --- 2. Gene Expression Data ---
expression_data = {
    "sample_id": mutation_df["sample_id"],
    "label": mutation_df["label"],
}

for gene in GENES:
    gene_expression = []
    for _, row in mutation_df.iterrows():
        label = row["label"]
        mutated = row[gene]
        base_mean = 8
        noise = np.random.normal(0, 1)

        if label == "cancer":
            if mutated:
                # Tumor suppressors go down, oncogenes go up
                if gene in ["BRCA1", "BRCA2", "TP53", "PTEN", "CDH1"]:
                    shift = np.random.normal(-2, 1)
                else:
                    shift = np.random.normal(2, 1)
            else:
                shift = np.random.normal(1, 0.5)
        else:
            shift = 0

        expression = base_mean + shift + noise
        gene_expression.append(expression)

    expression_data[gene] = gene_expression

expression_df = pd.DataFrame(expression_data)
expression_df.to_csv("data/expression_data.csv", index=False)


# --- 3. DNA Sequences ---
def simulate_dna_sequence(length: int = 1000) -> str:
    """Creates DNA sequences with 1000 bp and common mutations."""
    return "".join(random.choices("ACGT", k=length))


reference_sequences = {gene: simulate_dna_sequence() for gene in GENES}
records = []

# Reference sequences
for gene, seq in reference_sequences.items():
    records.append(
        SeqRecord(Seq(seq), id=f"{gene}_REF", description="Reference")
    )

# Mutated sequences using real-looking common mutations
COMMON_MUTATIONS = {
    "BRCA1": [(100, "A"), (250, "G"), (430, "T")],
    "BRCA2": [(150, "C"), (520, "G"), (890, "T"), (900, "A")],
    "TP53": [(75, "G"), (300, "A"), (650, "T")],
    "PIK3CA": [(120, "C"), (400, "T")],
    "PTEN": [(60, "A"), (200, "G"), (310, "C")],
    "EGFR": [(50, "G"), (700, "A"), (750, "T")],
    "HER2": [(180, "T"), (340, "C"), (560, "G")],
    "AKT1": [(90, "T"), (150, "G")],
    "CDH1": [(25, "A"), (330, "C"), (475, "T")],
    "BRAF": [(200, "G"), (500, "A")],
}

# Generate 5 variants per gene, each with 1–5 mutations randomly chosen
for gene, ref_seq in reference_sequences.items():
    known_mutations = COMMON_MUTATIONS[gene]
    for i in range(5):
        seq_list = list(ref_seq)
        num_mutations = random.randint(1, len(known_mutations))
        selected_mutations = random.sample(known_mutations, k=num_mutations)

        for pos, new_base in selected_mutations:
            if pos < len(seq_list):
                seq_list[pos] = new_base  # Overwrite base at known position

        mutated_seq = "".join(seq_list)
        records.append(
            SeqRecord(
                Seq(mutated_seq),
                id=f"{gene}_VAR{i + 1}",
                description=f"Mutated with {num_mutations} known changes",
            )
        )

# Save to FASTA
SeqIO.write(records, "data/sequences.fasta", "fasta")
