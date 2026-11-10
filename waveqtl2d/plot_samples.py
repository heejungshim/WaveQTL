"""
plot_samples.py -- Render side-by-side individual sample 2D matrix grids per comparison.

Adapts grid geometry automatically to sample size:
  - Small sample count (max(N1, N2) <= 6): K rows x 2 columns (Left: cond1, Right: cond2).
  - Larger sample count (max(N1, N2) > 6): ceil(K/5) rows x 10 columns
    (Left 5 columns: cond1, Right 5 columns: cond2, separated by a dashed divider).
"""

import json
import math
from pathlib import Path
from typing import List, Optional, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def _auto_symmetric_vmax(values: np.ndarray, default: float = 0.60) -> float:
    """Compute a clean symmetric color limit from observed values."""
    clean = np.abs(values[~np.isnan(values)])
    if clean.size == 0:
        return default
    q = float(np.percentile(clean, 98.0))
    if q <= 1e-9:
        return default
    return round(max(q, 0.05), 2)


def build_raw_sample_matrix(
    sample_df: pd.DataFrame,
    start_cell: int,
    end_cell: int,
) -> Tuple[np.ndarray, float, float]:
    """Construct a symmetric (dim x dim) matrix for one sample (NaN = unobserved)."""
    dim = end_cell - start_cell + 1
    mat = np.full((dim, dim), np.nan, dtype=np.float64)
    if sample_df.empty:
        return mat, 0.0, np.nan

    rows = sample_df["cell_i"].to_numpy(dtype=np.int64)
    cols = sample_df["cell_j"].to_numpy(dtype=np.int64)
    vals = sample_df["value"].to_numpy(dtype=np.float64)

    mask = (
        (rows >= start_cell)
        & (rows <= end_cell)
        & (cols >= start_cell)
        & (cols <= end_cell)
        & (~np.isnan(vals))
    )
    gi = rows[mask] - start_cell
    gj = cols[mask] - start_cell
    v = vals[mask]
    mat[gi, gj] = v
    mat[gj, gi] = v

    n_obs = int(np.sum(~np.isnan(mat)))
    cov_frac = n_obs / (dim * dim)
    obs_mean = float(np.nanmean(mat)) if n_obs > 0 else np.nan
    return mat, cov_frac, obs_mean


def plot_single_comparison_samples(
    comp_dir: Path,
    output_png: Path,
    vmax: Optional[float] = None,
    dpi: int = 150,
) -> None:
    """Render individual sample 2D matrices for a single prepared comparison."""
    sample_file = comp_dir / "sample_matrices.tsv.gz"
    meta_file = comp_dir / "comparison_meta.json"
    if not sample_file.exists() or not meta_file.exists():
        print(f"Skipping sample plot for {comp_dir.name}: missing prepared sample files.")
        return

    with open(meta_file, "r", encoding="utf-8") as f:
        meta = json.load(f)

    df = pd.read_csv(sample_file, sep="\t")
    c1 = str(meta["cond1"])
    c2 = str(meta["cond2"])
    n1 = int(meta["n_cond1"])
    n2 = int(meta["n_cond2"])
    start_cell = int(meta["start_cell"])
    end_cell = int(meta["end_cell"])
    dim = int(meta["grid_dim"])

    plot_vmax = float(vmax) if vmax is not None else _auto_symmetric_vmax(df["value"].to_numpy(dtype=np.float64))
    plot_vmin = -plot_vmax

    cmap = plt.get_cmap("RdBu_r").copy()
    cmap.set_bad(color="#cccccc")

    max_n = max(n1, n2)
    output_png.parent.mkdir(parents=True, exist_ok=True)

    if max_n <= 6:
        fig, axes = plt.subplots(max_n, 2, figsize=(8.0, 3.6 * max_n), dpi=dpi, squeeze=False)
        fig.subplots_adjust(left=0.08, right=0.88, top=0.93, bottom=0.04, hspace=0.35, wspace=0.25)
        fig.suptitle(
            f"{meta['comp_id']} — Individual Sample 2D Matrices ({n1} {c1} vs {n2} {c2})\n"
            f"Shared Cells: [{start_cell}..{end_cell}] ({dim}x{dim}) | Grey = Unobserved",
            fontsize=12.5,
            fontweight="bold",
            y=0.98,
        )

        im_ref = None
        for r_idx in range(max_n):
            s_order = r_idx + 1
            for col_idx, (cond_label, n_cond) in enumerate(((c1, n1), (c2, n2))):
                ax = axes[r_idx, col_idx]
                if s_order > n_cond:
                    ax.axis("off")
                    continue
                ax.set_facecolor("#cccccc")
                sub = df[(df["condition"].astype(str) == cond_label) & (df["sample_order"] == s_order)]
                s_name = str(sub["sample_id"].iloc[0]) if not sub.empty else str(s_order)
                mat, cov_f, m_val = build_raw_sample_matrix(sub, start_cell, end_cell)
                im = ax.imshow(mat, origin="upper", cmap=cmap, vmin=plot_vmin, vmax=plot_vmax, interpolation="nearest")
                if im_ref is None:
                    im_ref = im
                m_str = f"{m_val:+.2f}" if not np.isnan(m_val) else "N/A"
                ax.set_title(
                    f"{cond_label} #{s_name} (cov {cov_f * 100:.0f}%, mean {m_str})",
                    fontsize=9.5,
                    fontweight="bold",
                )
                ax.set_xticks([])
                ax.set_yticks([])

        if im_ref is not None:
            cbar_ax = fig.add_axes([0.91, 0.25, 0.025, 0.50])
            cbar = fig.colorbar(im_ref, cax=cbar_ax)
            cbar.set_label("Matrix Value", fontsize=11, fontweight="bold")
    else:
        cols_per_cond = 5
        nrows = int(math.ceil(max_n / cols_per_cond))
        fig, axes = plt.subplots(
            nrows=nrows,
            ncols=cols_per_cond * 2,
            figsize=(28, max(3.2 * nrows, 6.5)),
            dpi=dpi,
            squeeze=False,
            gridspec_kw={"wspace": 0.15, "hspace": 0.35, "right": 0.92, "left": 0.04, "top": 0.86, "bottom": 0.06},
        )
        fig.suptitle(
            f"{meta['comp_id']} — Individual Sample 2D Matrices ({n1} {c1} vs {n2} {c2}) | "
            f"Shared Cells [{start_cell}..{end_cell}] ({dim}x{dim})",
            fontsize=14,
            fontweight="bold",
            y=0.96,
        )
        fig.text(0.26, 0.905, f"{c1} ({n1} Samples)", fontsize=13, fontweight="bold", ha="center", color="#1f4e79")
        fig.text(0.70, 0.905, f"{c2} ({n2} Samples)", fontsize=13, fontweight="bold", ha="center", color="#800000")

        im_ref = None
        for r in range(nrows):
            for c in range(cols_per_cond * 2):
                axes[r, c].axis("off")

        for cond_offset, cond_label, n_cond, title_color in (
            (0, c1, n1, "#1f4e79"),
            (cols_per_cond, c2, n2, "#800000"),
        ):
            for s_order in range(1, n_cond + 1):
                r = (s_order - 1) // cols_per_cond
                c = cond_offset + ((s_order - 1) % cols_per_cond)
                ax = axes[r, c]
                ax.axis("on")
                ax.set_facecolor("#cccccc")
                sub = df[(df["condition"].astype(str) == cond_label) & (df["sample_order"] == s_order)]
                s_name = str(sub["sample_id"].iloc[0]) if not sub.empty else str(s_order)
                mat, cov_f, _ = build_raw_sample_matrix(sub, start_cell, end_cell)
                im = ax.imshow(mat, origin="upper", cmap=cmap, vmin=plot_vmin, vmax=plot_vmax, interpolation="nearest")
                if im_ref is None:
                    im_ref = im
                obs_extra = ""
                if "n_obs" in sub.columns and not sub.empty:
                    obs_extra = f"N={int(sub['n_obs'].iloc[0])}, "
                ax.set_title(
                    f"Sample {s_name} ({obs_extra}cov={cov_f * 100:.0f}%)",
                    fontsize=8.5,
                    fontweight="bold",
                    color=title_color,
                )
                ax.set_xticks([])
                ax.set_yticks([])

        sep_line = plt.Line2D([0.48, 0.48], [0.06, 0.89], transform=fig.transFigure, color="#888888", linestyle="--", linewidth=1.5)
        fig.add_artist(sep_line)

        if im_ref is not None:
            cbar_ax = fig.add_axes([0.935, 0.20, 0.015, 0.60])
            cbar = fig.colorbar(im_ref, cax=cbar_ax)
            cbar.set_label("Matrix Value", fontsize=11, fontweight="bold")

    plt.savefig(output_png, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved sample matrix plot -> {output_png.name}")


def plot_all_sample_matrices(
    work_dir: Path,
    sample_plots_dir: Path,
    vmax: Optional[float] = None,
    dpi: int = 150,
) -> List[Path]:
    """Render sample matrix plots for all prepared comparisons in work_dir."""
    manifest_path = work_dir / "manifest.json"
    if manifest_path.exists():
        with open(manifest_path, "r", encoding="utf-8") as f:
            comp_dirs = [work_dir / item["comp_id"] for item in json.load(f)]
    else:
        comp_dirs = sorted(p for p in work_dir.iterdir() if p.is_dir() and (p / "sample_matrices.tsv.gz").exists())

    sample_plots_dir.mkdir(parents=True, exist_ok=True)
    outputs: List[Path] = []
    for comp_dir in comp_dirs:
        out_png = sample_plots_dir / f"{comp_dir.name}_samples.png"
        plot_single_comparison_samples(comp_dir, out_png, vmax=vmax, dpi=dpi)
        if out_png.exists():
            outputs.append(out_png)
    return outputs
