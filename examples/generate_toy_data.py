#!/usr/bin/env python3
"""
generate_toy_data.py -- Generate a synthetic 2D symmetric matrix dataset for smoke-testing WaveQTL-2D.

Creates two synthetic targets with 5 replicates per condition (Cond_A vs Cond_B):
  - Target_Null: 36x36 shared grid (Cond_A spans 0..38, Cond_B spans 0..35 to demonstrate
    automatic shared-coordinate alignment) with no difference between conditions.
  - Target_Signal: 40x40 grid where Cond_B exhibits a synthetic off-diagonal block shift (+0.9)
    over cells [10..17, 24..31] (and its symmetric mirror).
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def _base_symmetric_pattern(dim: int) -> np.ndarray:
    idx = np.arange(dim, dtype=np.float64)
    ii, jj = np.meshgrid(idx, idx, indexing="ij")
    band = 0.25 * np.exp(-np.abs(ii - jj) / 6.0)
    checker = 0.12 * np.sin(ii / 4.5) * np.cos(jj / 4.5)
    base = 0.5 * (band + checker + (band + checker).T)
    return base


def _matrix_to_records(mat: np.ndarray, target: str, condition: str, sample: int) -> list:
    records = []
    dim = mat.shape[0]
    for i in range(dim):
        for j in range(i, dim):
            records.append(
                {
                    "target": target,
                    "condition": condition,
                    "sample": sample,
                    "cell_i": i,
                    "cell_j": j,
                    "value": round(float(mat[i, j]), 6),
                }
            )
    return records


def generate_toy_dataset(output_tsv: Path, seed: int = 42, n_samples: int = 5) -> None:
    rng = np.random.default_rng(seed)
    rows = []

    # 1. Target_Null: mismatched bounds (Cond_A: 39x39, Cond_B: 36x36 -> intersection 0..35)
    base_null_full = _base_symmetric_pattern(39)
    for cond, dim in (("Cond_A", 39), ("Cond_B", 36)):
        base = base_null_full[:dim, :dim]
        for s in range(1, n_samples + 1):
            noise = rng.normal(0.0, 0.14, size=(dim, dim))
            sym_noise = 0.5 * (noise + noise.T)
            rows.extend(_matrix_to_records(base + sym_noise, "Target_Null", cond, s))

    # 2. Target_Signal: 40x40 with localized off-diagonal block perturbation in Cond_B
    dim_sig = 40
    base_sig = _base_symmetric_pattern(dim_sig)
    delta = np.zeros((dim_sig, dim_sig), dtype=np.float64)
    delta[9:16, 25:32] = 0.90
    delta[25:32, 9:16] = 0.90

    for cond in ("Cond_A", "Cond_B"):
        mean_mat = base_sig if cond == "Cond_A" else (base_sig + delta)
        for s in range(1, n_samples + 1):
            noise = rng.normal(0.0, 0.14, size=(dim_sig, dim_sig))
            sym_noise = 0.5 * (noise + noise.T)
            rows.extend(_matrix_to_records(mean_mat + sym_noise, "Target_Signal", cond, s))

    output_tsv.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(output_tsv, sep="\t", index=False)
    print(f"Wrote synthetic dataset ({len(rows):,} rows) -> {output_tsv}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate synthetic 2D sample matrix TSV(.gz) for WaveQTL-2D.")
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("examples/data/toy_2d_matrices.tsv.gz"),
        help="Destination path for synthetic TSV(.gz) file.",
    )
    parser.add_argument("--seed", type=int, default=42, help="Random seed.")
    args = parser.parse_args()
    generate_toy_dataset(args.output, seed=args.seed)


if __name__ == "__main__":
    main()
