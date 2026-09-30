"""Command-line entry point for the public demonstration."""

from __future__ import annotations

import argparse

from mpi4py import MPI

from .config import load_config
from .solver import run_simulation


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run the nondimensional FEniCSx hydraulic-fracture demo."
    )
    parser.add_argument(
        "--config",
        default="configs/demo_case.yaml",
        help="Path to the YAML configuration file.",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    config = load_config(args.config)
    output_dir = run_simulation(config)
    if MPI.COMM_WORLD.rank == 0:
        print(f"Simulation completed. Results: {output_dir}")


if __name__ == "__main__":
    main()
