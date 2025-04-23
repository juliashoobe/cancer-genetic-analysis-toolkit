"""Tests CLI."""

import subprocess
from pathlib import Path
from typing import Any

import pandas as pd
import pytest


def run_cli(command: list[str]) -> subprocess.CompletedProcess[str]:
    """Runs the CLI command and captures the output."""
    return subprocess.run(
        command,
        text=True,
        capture_output=True,
    )


@pytest.fixture
def mock_files(tmp_path: Path) -> tuple[Any, Any, dict[str, str]]:
    """Fixture to create mock data for testing."""
    mutation_data = {
        "Gene": ["BRCA1", "TP53"],
        "Mutation": ["Del", "Substitution"],
    }
    expression_data = {"BRCA1": [1.2, 3.4], "TP53": [2.1, 4.5]}

    mutation_df = pd.DataFrame(mutation_data)
    expression_df = pd.DataFrame(expression_data)

    sequences = {"seq1": "ATCGATAG", "seq2": "ATCGATCG"}

    mutation_file = tmp_path / "mock_mutation.csv"
    expression_file = tmp_path / "mock_expression.csv"
    fasta_file = tmp_path / "mock_sequences.fasta"

    mutation_df.to_csv(mutation_file, index=False)
    expression_df.to_csv(expression_file, index=False)

    with open(fasta_file, "w") as f:
        for seq_id, seq in sequences.items():
            f.write(f">{seq_id}\n{seq}\n")

    return mutation_df, expression_df, sequences


def test_cli_mutation_detection(
    mock_files: tuple[pd.DataFrame, pd.DataFrame, dict[str, str]],
) -> None:
    """Test the mutation detection CLI command."""
    fasta_file = "mock_sequences.fasta"

    result = run_cli(
        [
            "python",
            "toolkit.py",
            "mutation",
            "--ref-id",
            "seq1",
            "--target-id",
            "seq2",
            "--fasta",
            fasta_file,
        ]
    )

    assert result.returncode == 0


def test_cli_gene_expression(
    mock_files: tuple[pd.DataFrame, pd.DataFrame, dict[str, str]],
) -> None:
    """Test the gene expression CLI command."""
    expression_file = "mock_expression.csv"

    result = run_cli(
        [
            "python",
            "toolkit.py",
            "expression",
            "--gene",
            "BRCA1",
            "--expression-file",
            str(expression_file),
        ]
    )
    print(result.stderr)

    assert "BRCA1" in result.stdout


def test_cli_differential_expression(
    mock_files: tuple[pd.DataFrame, pd.DataFrame, dict[str, str]],
) -> None:
    """Test the differential expression CLI command."""
    mutation_file, expression_file = "mock_mutation.csv", "mock_expression.csv"

    result = run_cli(
        [
            "python",
            "toolkit.py",
            "differential",
            "--gene",
            "BRCA1",
            "--mutation-file",
            str(mutation_file),
            "--expression-file",
            str(expression_file),
        ]
    )

    assert "t-statistic" in result.stdout
