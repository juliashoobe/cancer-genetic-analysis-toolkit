"""Tests mutation_detection() function."""

import pytest

from toolkit import mutation_detection


@pytest.fixture
def mock_sequences() -> dict[str, str]:
    """Fixture to provide mock sequences for testing."""
    return {
        "seq1": "ATCGATAG",
        "seq2": "ATCGATCG",
        "seq3": "ATCGATAG",
        "seq4": "ATCGGTTG",
    }


def test_single_mutation(mock_sequences: dict[str, str]) -> None:
    """Test mutation_detection with single mutation."""
    reference_seq: str = mock_sequences["seq1"]
    mutated_seq: str = mock_sequences["seq2"]

    expected_result: list[dict[str, str]] = [
        {"position": "7", "reference_base": "A", "mutated_base": "C"}
    ]

    result = mutation_detection(reference_seq, mutated_seq)

    assert result == expected_result


def test_no_mutations(mock_sequences: dict[str, str]) -> None:
    """Test mutation_detection with no mutations."""
    reference_seq: str = mock_sequences["seq1"]
    mutated_seq: str = mock_sequences["seq3"]

    expected_result: list[dict[str, str]] = [
        {"message": "No mutations detected."}
    ]

    result = mutation_detection(reference_seq, mutated_seq)

    assert result == expected_result


def test_multiple_mutations(mock_sequences: dict[str, str]) -> None:
    """Test mutation_detection with multiple mutations."""
    reference_seq: str = mock_sequences["seq1"]
    mutated_seq: str = mock_sequences["seq4"]

    expected_result: list[dict[str, str]] = [
        {"position": "5", "reference_base": "A", "mutated_base": "G"},
        {"position": "7", "reference_base": "A", "mutated_base": "T"},
    ]

    result = mutation_detection(reference_seq, mutated_seq)

    assert result == expected_result
