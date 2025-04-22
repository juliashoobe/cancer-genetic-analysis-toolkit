"""Tests parse_data() function."""

from collections.abc import Iterator
from typing import Any

import pandas as pd
import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from toolkit import parse_data


@pytest.fixture
def mock_files() -> tuple[pd.DataFrame, pd.DataFrame, dict[str, str]]:
    """Fixture to create mock data for testing parse_data function."""
    mutation_data = {
        "Gene": ["BRCA1", "TP53"],
        "Mutation": ["Del", "Substitution"],
    }
    expression_data = {"BRCA1": [1.2, 3.4], "TP53": [2.1, 4.5]}

    mutation_df = pd.DataFrame(mutation_data)
    expression_df = pd.DataFrame(expression_data)

    sequences = {"seq1": "ATCGATAG", "seq2": "ATCGATCG"}

    return mutation_df, expression_df, sequences


def test_parse_data(
    mock_files: tuple[pd.DataFrame, pd.DataFrame, dict[str, str]],
) -> None:
    """Test the parse_data function to ensure it correctly parses."""
    mutation_df, expression_df, sequences = mock_files
    mutation_file = "mock_mutation.csv"
    expression_file = "mock_expression.csv"
    fasta_file = "mock_sequences.fasta"

    # Mocking pd.read_csv and SeqIO.parse
    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(
            pd,
            "read_csv",
            lambda file_path: mutation_df
            if "mutation" in file_path
            else expression_df,
        )

        def mock_parse(file_path: Any, format: Any) -> Iterator[SeqRecord]:
            return iter(
                [
                    SeqRecord(Seq("ATCGATAG"), id="seq1"),
                    SeqRecord(Seq("ATCGATCG"), id="seq2"),
                ]
            )

        mp.setattr("Bio.SeqIO.parse", mock_parse)

        result = parse_data(mutation_file, expression_file, fasta_file)

        expected_result = (mutation_df, expression_df, sequences)

        assert result == expected_result
