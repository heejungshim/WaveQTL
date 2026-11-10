"""
cli.py -- Simple command-line entrypoint for WaveQTL-2D.
"""

import argparse
from pathlib import Path
import time

from waveqtl2d.execute import execute_all_comparisons
from waveqtl2d.plot_samples import plot_all_sample_matrices
from waveqtl2d.prepare import prepare_all_inputs
from waveqtl2d.visualize import summarize_and_visualize_all


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="waveqtl-2d",
        description="Compare two groups of 2D symmetric sample matrices using 2D Haar wavelets and WaveQTL.",
    )
    parser.add_argument(
        "--input",
        type=Path,
        required=False,
        help="Path to input per-sample 2D matrix table (.tsv or .tsv.gz). Required when running 'all' or 'prepare'.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Output directory for figures, summary TSVs, and intermediate WaveQTL working files.",
    )
    parser.add_argument(
        "--step",
        choices=["all", "prepare", "execute", "plot-samples", "visualize"],
        default="all",
        help="Pipeline stage to execute (default: all).",
    )
    parser.add_argument(
        "--group-by",
        nargs="*",
        default=["target"],
        help="Zero, one, or more columns grouping rows into comparisons (default: target). Pass empty `--group-by` for a single-comparison TSV.",
    )
    parser.add_argument("--cond-col", type=str, default="condition", help="Condition column name (default: condition).")
    parser.add_argument("--sample-col", type=str, default="sample", help="Sample replicate column name (default: sample).")
    parser.add_argument("--row-col", type=str, default="cell_i", help="2D row coordinate column name (default: cell_i).")
    parser.add_argument("--col-col", type=str, default="cell_j", help="2D column coordinate column name (default: cell_j).")
    parser.add_argument("--val-col", type=str, default="value", help="Matrix cell value column name (default: value).")
    parser.add_argument("--cond1", type=str, default=None, help="Optional baseline condition label (coded 0.0).")
    parser.add_argument("--cond2", type=str, default=None, help="Optional comparison condition label (coded 1.0; effect = cond2 - cond1).")
    parser.add_argument("--min-grid-size", type=int, default=64, help="Minimum dyadic 2^J grid size for 2D DWT padding (default: 64).")
    parser.add_argument("--num-perm", type=int, default=1000, help="Number of WaveQTL permutations (default: 1000; set 0 for analytical BF only).")
    parser.add_argument("--waveqtl-bin", type=Path, default=None, help="Optional path to compiled WaveQTL C++ binary.")
    parser.add_argument(
        "--full-matrix-tsv",
        type=Path,
        default=None,
        help="Optional TSV(.gz) of unpartitioned condition matrices to display in Columns 1-5 of comparison plots (defaults to sample averages).",
    )
    parser.add_argument(
        "--full-cond-col",
        type=str,
        default=None,
        help="Optional condition column name inside --full-matrix-tsv if different from --cond-col.",
    )
    parser.add_argument("--vmax", type=float, default=None, help="Optional fixed symmetric heatmap color limit [-vmax, +vmax].")
    parser.add_argument("--dpi", type=int, default=200, help="Output PNG resolution in DPI (default: 200).")
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    out_dir = args.output_dir
    work_dir = out_dir / "waveqtl_work"
    sample_plots_dir = out_dir / "sample_plots"

    t0 = time.time()

    if args.step in ("all", "prepare"):
        if args.input is None:
            parser.error("--input is required when --step is 'all' or 'prepare'.")
        print(f"=== [1/4] Preparing 2D Haar Wavelet Inputs from {args.input} ===")
        prepare_all_inputs(
            input_tsv=args.input,
            work_dir=work_dir,
            group_by=args.group_by,
            cond_col=args.cond_col,
            sample_col=args.sample_col,
            row_col=args.row_col,
            col_col=args.col_col,
            val_col=args.val_col,
            cond1=args.cond1,
            cond2=args.cond2,
            min_grid_size=args.min_grid_size,
        )

    if args.step in ("all", "plot-samples"):
        print("\n=== [2/4] Rendering Individual Sample Matrix Grids ===")
        plot_all_sample_matrices(
            work_dir=work_dir,
            sample_plots_dir=sample_plots_dir,
            vmax=args.vmax,
            dpi=min(args.dpi, 150),
        )

    if args.step in ("all", "execute"):
        print(f"\n=== [3/4] Running WaveQTL Association Testing (num_perm={args.num_perm}) ===")
        execute_all_comparisons(
            work_dir=work_dir,
            num_perm=args.num_perm,
            waveqtl_bin=args.waveqtl_bin,
        )

    if args.step in ("all", "visualize"):
        print("\n=== [4/4] Reconstructing 2D Effects & Generating Summary Figures ===")
        summary_df, _ = summarize_and_visualize_all(
            work_dir=work_dir,
            output_dir=out_dir,
            full_matrix_tsv=args.full_matrix_tsv,
            full_cond_col=args.full_cond_col,
            row_col=args.row_col,
            col_col=args.col_col,
            val_col=args.val_col,
            vmax=args.vmax,
            dpi=args.dpi,
        )
        if not summary_df.empty:
            print("\n=== Ranked Summary ===")
            cols = ["comp_id", "log_lr", "p_perm", "max_log10_bf", "top1_domain_bf", "frac_cells_z2", "frac_cells_z3"]
            print(summary_df[cols].to_string(index=False))

    print(f"\nCompleted in {time.time() - t0:.1f}s. Outputs saved to: {out_dir}")


if __name__ == "__main__":
    main()
