"""Tests differntial_expression() function."""

import pandas as pd
import pytest

from toolkit import differential_expression


@pytest.fixture
def mock_data() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Creates mock mutation and gene expression data for testing."""
    mutation_data = {
        "sample_id": ["sample1", "sample2", "sample3", "sample4", "sample5"],
        "label": [
            "cancer",
            "non-cancer",
            "cancer",
            "non-cancer",
            "cancer",
        ],
    }
    mutation_df = pd.DataFrame(mutation_data)

    expression_data = {
        "TP53": [2.5, 3.0, 2.8, 4.5, 5.0],
        "BRCA1": [5.1, 5.3, 5.2, 4.0, 4.8],
        "MYC": [1.2, 1.4, 1.3, 1.1, 1.0],
    }
    expression_df = pd.DataFrame(expression_data)

    return mutation_df, expression_df


def test_differential_expression_valid(
    mock_data: tuple[pd.DataFrame, pd.DataFrame],
) -> None:
    """Tests function with valid gene and comparison between groups."""
    mutation_df, expression_df = mock_data
    result = differential_expression(mutation_df, expression_df, "BRCA1")

    assert result is not None and all(
        key in result for key in ["t_statistic", "p_value"]
    )


def test_differential_expression_gene_not_found(
    mock_data: tuple[pd.DataFrame, pd.DataFrame],
) -> None:
    """Tests function when gene is not found in expression dataset."""
    mutation_df, expression_df = mock_data
    result = differential_expression(mutation_df, expression_df, "EGFR")

    assert result is None


def test_differential_expression_small_groups(
    mock_data: tuple[pd.DataFrame, pd.DataFrame],
) -> None:
    """Tests function with small sample sizes."""
    mutation_df, expression_df = mock_data
    mutation_df = mutation_df.head(2)
    expression_df = expression_df.head(2)

    result = differential_expression(mutation_df, expression_df, "TP53")

    assert result is None
