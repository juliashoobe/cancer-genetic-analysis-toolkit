"""Tests CLI."""

import contextlib
import sys
from io import StringIO
from pathlib import Path
from typing import Any

import pytest

from cli import main


@pytest.fixture
def mock_files(tmp_path: Path) -> tuple[str, str, str]:
    """Create temporary CSV files."""
    mutation_path = tmp_path / "mutation.csv"
    expression_path = tmp_path / "expression.csv"
    fasta_path = tmp_path / "sequences.fasta"

    mutation_data = """sample_id,label
sample1,cancer
sample2,non-cancer
sample3,cancer
sample4,non-cancer
sample5,cancer"""
    expression_data = """TP53,BRCA1,MYC
2.5,5.1,1.2
3.0,5.3,1.4
2.8,5.2,1.3
4.5,4.0,1.1
5.0,4.8,1.0"""
    fasta_data = """>ref
ATCGATCGATCG
>mut
ATCGATTGATCG"""

    mutation_path.write_text(mutation_data)
    expression_path.write_text(expression_data)
    fasta_path.write_text(fasta_data)

    return str(mutation_path), str(expression_path), str(fasta_path)


def run_cli(monkeypatch: Any, cli_args: Any) -> str:
    """Tests the CLI."""
    monkeypatch.setattr(sys, "argv", cli_args)
    captured_output = StringIO()
    monkeypatch.setattr(sys, "stdout", captured_output)
    monkeypatch.setattr(sys, "stderr", captured_output)
    with contextlib.suppress(SystemExit):
        main()
    return captured_output.getvalue()


def test_mutation_command(
    monkeypatch: pytest.MonkeyPatch, mock_files: tuple[str, str, str]
) -> None:
    """Tests the CLI mutation_detection command."""
    mutation_file, _, fasta_file = mock_files
    args = [
        "cli.py",
        "mutation",
        "--ref-id",
        "ref",
        "--target-id",
        "mut",
        "--fasta",
        fasta_file,
        "--mutation-file",
        mutation_file,
    ]
    output = run_cli(monkeypatch, args)
    assert "position" in output
    assert "reference_base" in output
    assert "mutated_base" in output


def test_expression_command(
    monkeypatch: pytest.MonkeyPatch, mock_files: tuple[str, str, str]
) -> None:
    """Tests the CLI gene_expression command."""
    _, expression_file, _ = mock_files
    args = [
        "cli.py",
        "expression",
        "--gene",
        "TP53",
        "--expression-file",
        expression_file,
    ]
    output = run_cli(monkeypatch, args)
    assert "2.5" in output
    assert "5.0" in output


def test_differential_expression_command(
    monkeypatch: pytest.MonkeyPatch, mock_files: tuple[str, str, str]
) -> None:
    """Tests the CLI differential_expression command."""
    mutation_file, expression_file, _ = mock_files
    args = [
        "cli.py",
        "differential",
        "--gene",
        "TP53",
        "--mutation-file",
        mutation_file,
        "--expression-file",
        expression_file,
    ]
    output = run_cli(monkeypatch, args)
    assert "t-statistic" in output
    assert "p-value" in output
