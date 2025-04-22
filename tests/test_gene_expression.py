"""Tests gene_expression() function."""

import pandas as pd
import pytest

from toolkit import gene_expression


@pytest.fixture
def mock_expression() -> pd.DataFrame:
    """Creates a mock gene expression DataFrame."""
    data = {
        "TP53": [2.5, 3.0, 2.8],
        "BRCA1": [5.1, 5.3, 5.2],
        "MYC": [1.2, 1.4, 1.3],
    }
    return pd.DataFrame(data)


def test_gene_found(mock_expression: pd.DataFrame) -> None:
    """Check that correct values returned for existing gene."""
    result = gene_expression("BRCA1", mock_expression)
    expected = pd.DataFrame({"BRCA1": [5.1, 5.3, 5.2]})
    pd.testing.assert_frame_equal(result, expected)


def test_gene_not_found(
    mock_expression: pd.DataFrame, capsys: pytest.CaptureFixture[str]
) -> None:
    """Test when the gene is not found in the DataFrame."""
    result = gene_expression("EGFR", mock_expression)
    captured = capsys.readouterr()
    assert result is None
    assert "Gene 'EGFR' not found in expression dataset." in captured.out


def test_gene_case_sensitivity(mock_expression: pd.DataFrame) -> None:
    """Ensures gene matching is case-sensitive."""
    result = gene_expression("brca1", mock_expression)
    assert result is None
