"""CLI for Cancer Genetic Analysis Toolkit."""

import argparse
import logging
import sys

from toolkit import (
    differential_expression,
    gene_expression,
    mutation_detection,
    parse_data,
)

# Set up logging for debug information
logging.basicConfig(level=logging.INFO, format="%(message)s")


def run_mutation_detection(args: argparse.Namespace) -> None:
    """Run mutation detection given reference and target sequence ID."""
    # Load the sequences and mutations from the data files
    _, _, sequences = parse_data(
        mutation_file=args.mutation_file,
        expression_file="",
        fasta_file=args.fasta,
    )

    reference_seq = sequences.get(args.ref_id)
    target_seq = sequences.get(args.target_id)

    if reference_seq is None or target_seq is None:
        print(
            "Error: One or both sequence IDs not found in FASTA file.",
            file=sys.stderr,
        )
        sys.exit(1)

    # Detect mutations between the reference and target sequences
    mutations = mutation_detection(reference_seq, target_seq)
    if mutations:
        for mutation in mutations:
            print(mutation)
    else:
        print("No mutations detected.", file=sys.stderr)


def run_gene_expression(args: argparse.Namespace) -> None:
    """Retrieve and print expression values for specified gene."""
    # Load expression data from the CSV file
    _, expression_df, _ = parse_data(
        mutation_file="", expression_file=args.expression_file, fasta_file=""
    )

    expression = gene_expression(args.gene, expression_df)
    if expression is not None:
        print(expression.to_string(index=False))
    else:
        print(
            f"No expression data found for gene: {args.gene}", file=sys.stderr
        )


def run_differential_expression(args: argparse.Namespace) -> None:
    """Perform mutation-driven differential gene expression analysis."""
    # Load mutation and expression data from the CSV files
    mutation_df, expression_df, _ = parse_data(
        mutation_file=args.mutation_file,
        expression_file=args.expression_file,
        fasta_file="",
    )

    result = differential_expression(mutation_df, expression_df, args.gene)
    if result is not None:
        print(f"t-statistic: {result['t_statistic']:.4f}")
        print(f"p-value: {result['p_value']:.4e}")
    else:
        print(
            "No results found for differential expression analysis.",
            file=sys.stderr,
        )


def main() -> None:
    """Main function to parse arguments and dispatch subcommands."""
    parser = argparse.ArgumentParser(
        description="Cancer Genetic Analysis Toolkit CLI"
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    # Mutation Detection
    mutation_parser = subparsers.add_parser(
        "mutation", help="Run mutation detection"
    )
    mutation_parser.add_argument(
        "--ref-id", required=True, help="Reference sequence ID"
    )
    mutation_parser.add_argument(
        "--target-id", required=True, help="Target sequence ID"
    )
    mutation_parser.add_argument(
        "--fasta", required=True, help="path to FASTA file"
    )
    mutation_parser.add_argument(
        "--mutation-file", required=True, help="Path to mutation CSV file"
    )
    mutation_parser.set_defaults(func=run_mutation_detection)

    # Gene Expression
    expression_parser = subparsers.add_parser(
        "expression", help="Retrieve gene expression values"
    )
    expression_parser.add_argument("--gene", required=True, help="Gene name")
    expression_parser.add_argument(
        "--expression-file", required=True, help="Path to expression CSV file"
    )
    expression_parser.set_defaults(func=run_gene_expression)

    # Differential Expression
    diff_parser = subparsers.add_parser(
        "differential", help="Run mutation-driven differential expression"
    )
    diff_parser.add_argument("--gene", required=True, help="Gene name")
    diff_parser.add_argument(
        "--mutation-file", required=True, help="Path to mutation CSV file"
    )
    diff_parser.add_argument(
        "--expression-file", required=True, help="Path to expression CSV file"
    )
    diff_parser.set_defaults(func=run_differential_expression)

    args = parser.parse_args()

    if not args.command:
        parser.print_help()
        sys.exit(1)

    # Check required arguments for each command
    if (
        (
            args.command == "mutation"
            and (
                not args.ref_id
                or not args.target_id
                or not args.fasta
                or not args.mutation_file
            )
        )
        or (
            args.command == "expression"
            and (not args.gene or not args.expression_file)
        )
        or (
            args.command == "differential"
            and (
                not args.gene
                or not args.mutation_file
                or not args.expression_file
            )
        )
    ):
        print(
            "Error: the following arguments are required: ",
            end="",
            file=sys.stderr,
        )
        if args.command == "mutation":
            print(
                "--ref-id, --target-id, --fasta, --mutation-file",
                file=sys.stderr,
            )
        elif args.command == "expression":
            print("--gene, --expression-file", file=sys.stderr)
        elif args.command == "differential":
            print(
                "--gene, --mutation-file, --expression-file", file=sys.stderr
            )
        sys.exit(1)

    # Call the function associated with the command
    args.func(args)


if __name__ == "__main__":
    main()
