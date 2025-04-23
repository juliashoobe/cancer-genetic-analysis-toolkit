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

# --- COMMON_MUTATIONS (Point mutations) ---
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

# --- 1. Mutation Data ---
sample_ids = [f"S{i + 1}" for i in range(SAMPLES)]
labels = np.random.choice(
    LABELS, size=SAMPLES, p=[CANCER_RATES, 1 - CANCER_RATES]
)

mutation_records = []
for sample_id, label in zip(sample_ids, labels):
    row = {"sample_id": sample_id, "label": label}
    mutation_descriptions = []

    for gene in GENES:
        if label == "cancer":
            prob = MUTATION_PROBS[gene]
        else:
            prob = MUTATION_PROBS[gene] * 0.1
        is_mutated = np.random.binomial(1, prob)
        row[gene] = bool(is_mutated)

        if is_mutated:
            known_mutations = COMMON_MUTATIONS[gene]
            k = random.randint(1, len(known_mutations))
            selected = random.sample(known_mutations, k=k)

            desc = ";".join([f"{gene}:{pos}{base}" for pos, base in selected])
            mutation_descriptions.append(desc)

    if mutation_descriptions:
        row["mutation"] = ";".join(mutation_descriptions)
    else:
        row["mutation"] = "none"

    mutation_records.append(row)


mutation_df = pd.DataFrame(mutation_records)
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

    # Add the gene expression list to the dictionary
    expression_data[gene] = gene_expression

# Convert the expression data into a DataFrame
expression_df = pd.DataFrame(expression_data)
expression_df.to_csv("data/expression_data.csv", index=False)


# --- 3. DNA Sequences (Point Mutations Only) ---
def simulate_dna_sequence(length: int = 1000) -> str:
    """Creates DNA sequences with 1000 bp and common point mutations."""
    return "".join(random.choices("ACGT", k=length))


reference_sequences = {gene: simulate_dna_sequence() for gene in GENES}
records = []

# Reference sequences (unmutated)
for gene, seq in reference_sequences.items():
    records.append(SeqRecord(Seq(seq), id=f"{gene}_REF", description=""))

# Mutated sequences with only point mutations
for gene, ref_seq in reference_sequences.items():
    known_mutations = COMMON_MUTATIONS[gene]
    for i in range(5):  # 5 mutated variants per gene
        seq_list = list(ref_seq)  # Convert to list for mutation
        num_mutations = random.randint(1, len(known_mutations))
        selected_mutations = random.sample(known_mutations, k=num_mutations)

        for pos, new_base in selected_mutations:
            if pos < len(seq_list):
                seq_list[pos] = new_base

        mutated_seq = "".join(seq_list)
        records.append(
            SeqRecord(
                Seq(mutated_seq),
                id=f"{gene}_VAR{i + 1}",
                description="",
            )
        )

# Save to FASTA file
SeqIO.write(records, "data/sequences.fasta", "fasta")
