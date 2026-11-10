r"""
visualize.py -- Reconstruct 2D WaveQTL posterior effects, rank comparisons, and generate figures.

Outputs generated in <output_dir>/:
  - summary_ranking.tsv: overall logLR (\log \hat{\Lambda}), permutation p-value, and subregion metrics
  - top_domain_subregions.tsv: Top-3 non-overlapping localized square subregions per comparison
  - comparison_<target>.png: 9-column comparison figure(s) showing condition matrices, empirical diff,
    Top-3 wavelet boxes, reconstructed Z-scores, and posterior effect contours (|Z| > 2 and |Z| > 3)
  - summary_ranking.png: 2-panel cross-comparison summary bar chart
"""

from collections import OrderedDict
import json
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import matplotlib
matplotlib.use("Agg")
from matplotlib.collections import LineCollection
import matplotlib.patches as patches
import matplotlib.patheffects as patheffects
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pywt


BOX_COLORS = ["#00e600", "#ffd700", "#00e5ff"]  # #1 Lime Green, #2 Gold, #3 Cyan


def unflatten_to_pywt_list(vec: np.ndarray, decomp_levels: int) -> list:
    """Reconstruct pywt.wavedec2 coefficient list from 1D vector of length 4^J."""
    coeffs = [vec[0:1].reshape((1, 1))]
    pos = 1
    for lvl in range(1, decomp_levels + 1):
        s_dim = 2 ** (lvl - 1)
        n_elem = s_dim * s_dim
        c_h = vec[pos : pos + n_elem].reshape((s_dim, s_dim))
        pos += n_elem
        c_v = vec[pos : pos + n_elem].reshape((s_dim, s_dim))
        pos += n_elem
        c_d = vec[pos : pos + n_elem].reshape((s_dim, s_dim))
        pos += n_elem
        coeffs.append((c_h, c_v, c_d))
    return coeffs


def project_variance_to_dataspace(var_vec: np.ndarray, grid_size: int, decomp_levels: int) -> np.ndarray:
    """Project posterior variance of 2D Haar wavelet coefficients to 2D data space."""
    var_grid = np.zeros((grid_size, grid_size), dtype=np.float64)
    var_grid += var_vec[0] / (grid_size * grid_size)

    pos = 1
    for lvl in range(1, decomp_levels + 1):
        s_dim = 2 ** (lvl - 1)
        block = grid_size // s_dim
        w2_val = 1.0 / (block * block)
        n_elem = s_dim * s_dim

        var_h = var_vec[pos : pos + n_elem].reshape((s_dim, s_dim))
        pos += n_elem
        var_v = var_vec[pos : pos + n_elem].reshape((s_dim, s_dim))
        pos += n_elem
        var_d = var_vec[pos : pos + n_elem].reshape((s_dim, s_dim))
        pos += n_elem

        for u in range(s_dim):
            for v in range(s_dim):
                r0, r1 = u * block, (u + 1) * block
                c0, c1 = v * block, (v + 1) * block
                tot_var = var_h[u, v] + var_v[u, v] + var_d[u, v]
                if tot_var > 0:
                    var_grid[r0:r1, c0:c1] += tot_var * w2_val

    return var_grid


def select_top_domain_boxes(
    meta_df: pd.DataFrame,
    grid_dim: int,
    min_block_size: int = 4,
    max_block_size: int = 16,
    max_boxes: int = 3,
    iou_thresh: float = 0.35,
) -> List[dict]:
    """Select Top-K non-redundant localized square wavelet windows via greedy NMS."""
    cand = meta_df[
        (meta_df["use_flag"] == 1)
        & (~meta_df["is_redundant_symmetry"])
        & ((meta_df["r1"] - meta_df["r0"]) >= min_block_size)
        & ((meta_df["r1"] - meta_df["r0"]) <= max_block_size)
    ].copy()

    if cand.empty:
        return []

    cand["eff_r"] = np.minimum(cand["r1"], grid_dim) - cand["r0"]
    cand["eff_c"] = np.minimum(cand["c1"], grid_dim) - cand["c0"]
    cand = cand[(cand["eff_r"] >= min_block_size) & (cand["eff_c"] >= min_block_size)]
    if cand.empty:
        return []

    cand = cand.sort_values(by="log10_bf", ascending=False)
    selected: List[dict] = []

    for _, row in cand.iterrows():
        r0, r1 = int(row["r0"]), min(grid_dim, int(row["r1"]))
        c0, c1 = int(row["c0"]), min(grid_dim, int(row["c0"] + (row["r1"] - row["r0"])))
        area = (r1 - r0) * (c1 - c0)

        overlap = False
        for s in selected:
            sr0, sr1 = int(s["r0"]), min(grid_dim, int(s["r1"]))
            sc0, sc1 = int(s["c0"]), min(grid_dim, int(s["c0"] + (s["r1"] - s["r0"])))
            sarea = (sr1 - sr0) * (sc1 - sc0)

            inter_r = max(0, min(r1, sr1) - max(r0, sr0))
            inter_c = max(0, min(c1, sc1) - max(c0, sc0))
            inter_area = inter_r * inter_c

            if inter_area > 0:
                iou = inter_area / (area + sarea - inter_area)
                containment = inter_area / min(area, sarea)
                if iou > iou_thresh or containment > 0.6:
                    overlap = True
                    break

        if not overlap:
            selected.append(dict(row))
            if len(selected) >= max_boxes:
                break

    return selected


def draw_wavelet_boxes(
    ax: plt.Axes,
    top_boxes: List[dict],
    start_cell: int,
    end_cell: int,
) -> List[patches.Patch]:
    """Draw neon bounding boxes (plus symmetric mirror boxes) with black halo casing."""
    legend_handles: List[patches.Patch] = []
    for b_idx, box in enumerate(top_boxes):
        color = BOX_COLORS[b_idx % len(BOX_COLORS)]
        r0 = start_cell + int(box["r0"])
        r1 = min(end_cell + 1, start_cell + int(box["r1"]))
        c0 = start_cell + int(box["c0"])
        c1 = min(end_cell + 1, start_cell + int(box["c1"]))
        bf_val = float(box["log10_bf"])
        pip_val = float(box["pip"])

        rect_main = patches.Rectangle(
            (c0 - 0.5, r0 - 0.5),
            c1 - c0,
            r1 - r0,
            linewidth=2.6,
            edgecolor=color,
            facecolor="none",
            linestyle="-",
            zorder=5,
        )
        rect_main.set_path_effects(
            [patheffects.withStroke(linewidth=4.8, foreground="black", alpha=0.9), patheffects.Normal()]
        )
        ax.add_patch(rect_main)

        if r0 != c0 or r1 != c1:
            rect_mirror = patches.Rectangle(
                (r0 - 0.5, c0 - 0.5),
                r1 - r0,
                c1 - c0,
                linewidth=2.6,
                edgecolor=color,
                facecolor="none",
                linestyle="--",
                zorder=5,
            )
            rect_mirror.set_path_effects(
                [patheffects.withStroke(linewidth=4.8, foreground="black", alpha=0.9), patheffects.Normal()]
            )
            ax.add_patch(rect_mirror)

        label_text = f"#{b_idx + 1}: $\\log_{{10}}$BF={bf_val:.1f} (PIP={pip_val * 100:.0f}%)"
        legend_handles.append(patches.Patch(edgecolor=color, facecolor="none", label=label_text, linewidth=2.5))

    return legend_handles


def get_pixel_boundary_segments(mask: np.ndarray) -> List[List[Tuple[float, float]]]:
    """Extract horizontal and vertical segments tracing outer boundaries of True pixels."""
    nr, nc = mask.shape
    segments: List[List[Tuple[float, float]]] = []
    padded = np.pad(mask, 1, mode="constant", constant_values=False)

    for r in range(nr + 1):
        for c in range(nc):
            if padded[r, c + 1] != padded[r + 1, c + 1]:
                segments.append([(c - 0.5, r - 0.5), (c + 0.5, r - 0.5)])

    for c in range(nc + 1):
        for r in range(nr):
            if padded[r + 1, c] != padded[r + 1, c + 1]:
                segments.append([(c - 0.5, r - 0.5), (c - 0.5, r + 0.5)])

    return segments


def draw_stepped_boundaries(
    ax: plt.Axes,
    mask: np.ndarray,
    color: str = "black",
    linewidth: float = 1.8,
    zorder: int = 4,
) -> None:
    """Draw stepped pixel contours around significant cells."""
    if not np.any(mask):
        return
    segments = get_pixel_boundary_segments(mask)
    if segments:
        lc = LineCollection(segments, colors=color, linewidths=linewidth, zorder=zorder)
        ax.add_collection(lc)


def _load_condition_mean_matrices(
    comp_dir: Path,
    comp_meta: dict,
    max_cell: int,
    full_df: Optional[pd.DataFrame] = None,
    full_cond_col: Optional[str] = None,
    row_col: str = "cell_i",
    col_col: str = "cell_j",
    val_col: str = "value",
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Load Condition 1 and Condition 2 reference matrices for Columns 1-5:
      - If full_df is provided (--full-matrix-tsv), extract unpartitioned condition matrices.
      - Otherwise, average across samples in sample_matrices.tsv.gz.
    """
    c1 = str(comp_meta["cond1"])
    c2 = str(comp_meta["cond2"])
    mat1 = np.full((max_cell, max_cell), np.nan, dtype=np.float64)
    mat2 = np.full((max_cell, max_cell), np.nan, dtype=np.float64)

    if full_df is not None:
        sub = full_df
        for g_col, g_val in comp_meta.get("group_values", {}).items():
            if g_col in sub.columns:
                sub = sub[sub[g_col].astype(str) == str(g_val)]

        c_col = full_cond_col if (full_cond_col and full_cond_col in sub.columns) else None
        if c_col is None:
            for candidate in ("condition", "read_set", "run"):
                if candidate in sub.columns:
                    c_col = candidate
                    break
        if c_col is None:
            raise ValueError(f"Could not find condition column in --full-matrix-tsv: {list(full_df.columns)}")

        for cond_label, mat in ((c1, mat1), (c2, mat2)):
            c_sub = sub[sub[c_col].astype(str) == cond_label]
            if not c_sub.empty:
                rows = c_sub[row_col].to_numpy(dtype=np.int64)
                cols = c_sub[col_col].to_numpy(dtype=np.int64)
                vals = c_sub[val_col].to_numpy(dtype=np.float64)
                mask = (rows >= 0) & (rows < max_cell) & (cols >= 0) & (cols < max_cell) & (~np.isnan(vals))
                mat[rows[mask], cols[mask]] = vals[mask]
                mat[cols[mask], rows[mask]] = vals[mask]
    else:
        s_df = pd.read_csv(comp_dir / "sample_matrices.tsv.gz", sep="\t")
        for cond_label, mat in ((c1, mat1), (c2, mat2)):
            c_sub = s_df[s_df["condition"].astype(str) == cond_label]
            if not c_sub.empty:
                agg = c_sub.groupby(["cell_i", "cell_j"], as_index=False)["value"].mean()
                rows = agg["cell_i"].to_numpy(dtype=np.int64)
                cols = agg["cell_j"].to_numpy(dtype=np.int64)
                vals = agg["value"].to_numpy(dtype=np.float64)
                mask = (rows >= 0) & (rows < max_cell) & (cols >= 0) & (cols < max_cell) & (~np.isnan(vals))
                mat[rows[mask], cols[mask]] = vals[mask]
                mat[cols[mask], rows[mask]] = vals[mask]

    return mat1, mat2


def _auto_vmax_pair(mat1: np.ndarray, mat2: np.ndarray, default: float = 0.60) -> float:
    combined = np.concatenate([mat1[~np.isnan(mat1)], mat2[~np.isnan(mat2)]])
    if combined.size == 0:
        return default
    q = float(np.percentile(np.abs(combined), 98.5))
    if q <= 1e-9:
        return default
    return round(max(q, 0.05), 2)


def process_all_outputs(work_dir: Path, output_dir: Path) -> Tuple[pd.DataFrame, pd.DataFrame, List[dict]]:
    """Parse WaveQTL outputs, reconstruct 2D posterior maps, and save summary TSV tables."""
    manifest_path = work_dir / "manifest.json"
    if not manifest_path.exists():
        raise FileNotFoundError(f"Manifest not found at {manifest_path}. Run prepare & execute steps first.")

    with open(manifest_path, "r", encoding="utf-8") as f:
        manifest: List[dict] = json.load(f)

    summary_rows: List[dict] = []
    subregion_rows: List[dict] = []
    valid_manifest: List[dict] = []

    for comp_meta in manifest:
        comp_id = comp_meta["comp_id"]
        comp_dir = work_dir / comp_id
        out_dir = comp_dir / "output"

        loglr_file = out_dir / "waveqtl_out.fph.logLR.txt"
        pval_file = out_dir / "waveqtl_out.fph.pval.txt"
        pi_file = out_dir / "waveqtl_out.fph.pi.txt"
        mean_file = out_dir / "waveqtl_out.fph.mean.txt"
        var_file = out_dir / "waveqtl_out.fph.var.txt"
        meta_file = comp_dir / "coeff_metadata.tsv"

        if not (loglr_file.exists() and mean_file.exists() and var_file.exists() and meta_file.exists()):
            print(f"Skipping {comp_id}: missing WaveQTL output files.")
            continue

        grid_dim = int(comp_meta["grid_dim"])
        grid_size = int(comp_meta["grid_size"])
        decomp_levels = int(comp_meta["decomp_levels"])
        num_coeffs = int(comp_meta["num_coeffs"])
        start_cell = int(comp_meta["start_cell"])
        end_cell = int(comp_meta["end_cell"])

        loglr_df = pd.read_csv(loglr_file, sep=r"\s+", header=None)
        overall_loglr = float(loglr_df.iloc[0, 1])
        wc_log10_bfs = loglr_df.iloc[0, 2 : 2 + num_coeffs].to_numpy(dtype=np.float64)

        p_perm = np.nan
        num_perm = 0
        num_exceed = 0
        if pval_file.exists():
            with open(pval_file, "r", encoding="utf-8") as f:
                p_lines = [line.strip() for line in f.readlines() if line.strip()]
                if len(p_lines) >= 4:
                    num_perm = int(p_lines[1].split()[0])
                    num_exceed = int(p_lines[2].split()[0])
                    p_perm = float(p_lines[3].split()[0])

        pi_df = pd.read_csv(pi_file, sep=r"\s+", header=None)
        pi_vals = pi_df.iloc[0, 1 : 2 + decomp_levels].to_numpy(dtype=np.float64)

        meta_df = pd.read_csv(meta_file, sep="\t")
        meta_df["log10_bf"] = wc_log10_bfs

        pips = []
        for _, row in meta_df.iterrows():
            pi_s = pi_vals[int(row["scale"])]
            lbf = float(row["log10_bf"])
            if lbf > 30:
                pip = 1.0
            elif lbf < -30:
                pip = 0.0
            else:
                bf = 10.0 ** lbf
                pip = (pi_s * bf) / ((1.0 - pi_s) + pi_s * bf)
            pips.append(pip)
        meta_df["pip"] = pips

        mean_vec = pd.read_csv(mean_file, sep=r"\s+", header=None).iloc[0, 1 : 1 + num_coeffs].to_numpy(dtype=np.float64)
        var_vec = pd.read_csv(var_file, sep=r"\s+", header=None).iloc[0, 1 : 1 + num_coeffs].to_numpy(dtype=np.float64)

        sym_mean = mean_vec.copy()
        sym_var = var_vec.copy()
        for _, row in meta_df.iterrows():
            if bool(row["is_redundant_symmetry"]):
                orig_idx = int(row["mirror_coeff_idx"]) - 1
                curr_idx = int(row["coeff_idx"]) - 1
                sym_mean[curr_idx] = sym_mean[orig_idx]
                sym_var[curr_idx] = sym_var[orig_idx]

        pywt_coeffs = unflatten_to_pywt_list(sym_mean, decomp_levels=decomp_levels)
        diff_hat_full_grid = pywt.waverec2(pywt_coeffs, wavelet="haar")
        var_hat_full_grid = project_variance_to_dataspace(sym_var, grid_size=grid_size, decomp_levels=decomp_levels)
        sd_hat_full_grid = np.sqrt(np.maximum(var_hat_full_grid, 1e-15))
        z_hat_full_grid = np.divide(
            diff_hat_full_grid,
            sd_hat_full_grid,
            out=np.zeros_like(diff_hat_full_grid),
            where=(sd_hat_full_grid > 1e-12),
        )

        diff_hat = diff_hat_full_grid[:grid_dim, :grid_dim]
        sd_hat = sd_hat_full_grid[:grid_dim, :grid_dim]
        z_hat = z_hat_full_grid[:grid_dim, :grid_dim]

        frac_z2 = float(np.mean(np.abs(z_hat) > 2.0))
        frac_z3 = float(np.mean(np.abs(z_hat) > 3.0))

        top_boxes = select_top_domain_boxes(meta_df, grid_dim=grid_dim, min_block_size=4, max_block_size=16, max_boxes=3)
        box_bfs = [float(b["log10_bf"]) for b in top_boxes]
        while len(box_bfs) < 3:
            box_bfs.append(np.nan)

        for rank_i, b in enumerate(top_boxes, start=1):
            sub_rec = {"comp_id": comp_id}
            sub_rec.update(comp_meta.get("group_values", {}))
            sub_rec.update(
                {
                    "rank": rank_i,
                    "log10_bf": float(b["log10_bf"]),
                    "pip": float(b["pip"]),
                    "scale": int(b["scale"]),
                    "subband": str(b["subband"]),
                    "cell_r0": start_cell + int(b["r0"]),
                    "cell_r1": min(end_cell + 1, start_cell + int(b["r1"])),
                    "cell_c0": start_cell + int(b["c0"]),
                    "cell_c1": min(end_cell + 1, start_cell + int(b["c1"])),
                }
            )
            subregion_rows.append(sub_rec)

        active_sub = meta_df[meta_df["use_flag"] == 1]
        max_log10_bf = float(active_sub["log10_bf"].max()) if not active_sub.empty else np.nan
        max_pip = float(active_sub["pip"].max()) if not active_sub.empty else np.nan

        sum_rec = {"comp_id": comp_id}
        sum_rec.update(comp_meta.get("group_values", {}))
        sum_rec.update(
            {
                "cond1": comp_meta["cond1"],
                "cond2": comp_meta["cond2"],
                "n_cond1": comp_meta["n_cond1"],
                "n_cond2": comp_meta["n_cond2"],
                "start_cell": start_cell,
                "end_cell": end_cell,
                "log_lr": overall_loglr,
                "p_perm": p_perm,
                "num_perm": num_perm,
                "num_exceed": num_exceed,
                "active_wcs": int(meta_df["use_flag"].sum()),
                "max_log10_bf": max_log10_bf,
                "max_pip": max_pip,
                "top1_domain_bf": box_bfs[0],
                "top2_domain_bf": box_bfs[1],
                "top3_domain_bf": box_bfs[2],
                "frac_cells_z2": frac_z2,
                "frac_cells_z3": frac_z3,
            }
        )
        summary_rows.append(sum_rec)

        np.savez_compressed(
            comp_dir / "reconstructed_diff.npz",
            diff_hat=diff_hat,
            sd_hat=sd_hat,
            z_hat=z_hat,
            top_boxes=np.array(top_boxes, dtype=object),
            log_lr=overall_loglr,
            p_perm=p_perm,
        )
        valid_manifest.append(comp_meta)

    output_dir.mkdir(parents=True, exist_ok=True)
    summary_df = pd.DataFrame(summary_rows)
    if not summary_df.empty:
        summary_df = summary_df.sort_values(by="log_lr", ascending=False).reset_index(drop=True)
    top_sub_df = pd.DataFrame(subregion_rows)

    summary_df.to_csv(output_dir / "summary_ranking.tsv", sep="\t", index=False)
    top_sub_df.to_csv(output_dir / "top_domain_subregions.tsv", sep="\t", index=False)
    return summary_df, top_sub_df, valid_manifest


def plot_9col_group_figure(
    target_name: str,
    comp_metas: List[dict],
    work_dir: Path,
    output_png: Path,
    full_df: Optional[pd.DataFrame] = None,
    full_cond_col: Optional[str] = None,
    row_col: str = "cell_i",
    col_col: str = "cell_j",
    val_col: str = "value",
    vmax: Optional[float] = None,
    dpi: int = 200,
) -> None:
    """Render an (nrows x 9 cols) WaveQTL 2D comparison figure for one primary target."""
    nrows = len(comp_metas)
    fig, axes = plt.subplots(
        nrows=nrows,
        ncols=9,
        figsize=(45, 5.3 * nrows),
        constrained_layout=True,
        squeeze=False,
    )

    # Determine unified display canvas dimension across rows of this target
    max_cell = max(int(m["end_cell"]) + 1 for m in comp_metas)
    if full_df is not None and comp_metas[0].get("group_cols"):
        p_col = comp_metas[0]["group_cols"][0]
        if p_col in full_df.columns:
            p_sub = full_df[full_df[p_col].astype(str) == str(target_name)]
            if not p_sub.empty:
                max_cell = max(max_cell, int(max(p_sub[row_col].max(), p_sub[col_col].max())) + 1)

    cmap = plt.get_cmap("RdBu_r").copy()
    cmap.set_bad(color="#cccccc")

    c1 = str(comp_metas[0]["cond1"])
    c2 = str(comp_metas[0]["cond2"])
    n1 = int(comp_metas[0]["n_cond1"])
    n2 = int(comp_metas[0]["n_cond2"])

    for row_idx, comp_meta in enumerate(comp_metas):
        comp_id = comp_meta["comp_id"]
        comp_dir = work_dir / comp_id
        start_cell = int(comp_meta["start_cell"])
        end_cell = int(comp_meta["end_cell"])

        mat1, mat2 = _load_condition_mean_matrices(
            comp_dir=comp_dir,
            comp_meta=comp_meta,
            max_cell=max_cell,
            full_df=full_df,
            full_cond_col=full_cond_col,
            row_col=row_col,
            col_col=col_col,
            val_col=val_col,
        )
        raw_diff = mat2 - mat1

        npz = np.load(comp_dir / "reconstructed_diff.npz", allow_pickle=True)
        diff_hat = npz["diff_hat"]
        z_hat = npz["z_hat"]
        top_boxes = npz["top_boxes"].tolist()
        llr_val = float(npz["log_lr"])
        p_val = float(npz["p_perm"])

        if np.isnan(p_val):
            p_clause = "p = N/A"
        elif p_val >= 0.001:
            p_clause = f"p = {p_val:.4f}"
        else:
            p_clause = "p < 0.001"

        row_vmax = float(vmax) if vmax is not None else _auto_vmax_pair(mat1, mat2)
        row_vmin = -row_vmax
        vmax_diff = row_vmax * 0.75

        row_label = comp_meta.get("sub_label") or target_name

        # Col 1: Cond 1
        ax0 = axes[row_idx, 0]
        ax0.set_facecolor("#cccccc")
        im0 = ax0.imshow(mat1, cmap=cmap, vmin=row_vmin, vmax=row_vmax, origin="upper")
        ax0.set_title(f"{c1}", fontsize=11, fontweight="bold")
        ax0.set_ylabel(f"{row_label}\nMatrix Bins", fontsize=11, fontweight="bold")
        fig.colorbar(im0, ax=ax0, fraction=0.046, pad=0.04)

        # Col 2: Cond 2
        ax1 = axes[row_idx, 1]
        ax1.set_facecolor("#cccccc")
        im1 = ax1.imshow(mat2, cmap=cmap, vmin=row_vmin, vmax=row_vmax, origin="upper")
        ax1.set_title(f"{c2}", fontsize=11, fontweight="bold")
        fig.colorbar(im1, ax=ax1, fraction=0.046, pad=0.04)

        # Col 3: Cond 1 + Top Boxes
        ax2 = axes[row_idx, 2]
        ax2.set_facecolor("#cccccc")
        im2 = ax2.imshow(mat1, cmap=cmap, vmin=row_vmin, vmax=row_vmax, origin="upper")
        ax2.set_title(f"{c1} + Boxes", fontsize=11, fontweight="bold")
        fig.colorbar(im2, ax=ax2, fraction=0.046, pad=0.04)
        draw_wavelet_boxes(ax2, top_boxes, start_cell, end_cell)

        # Col 4: Cond 2 + Top Boxes + Legend
        ax3 = axes[row_idx, 3]
        ax3.set_facecolor("#cccccc")
        im3 = ax3.imshow(mat2, cmap=cmap, vmin=row_vmin, vmax=row_vmax, origin="upper")
        ax3.set_title(f"{c2} + Boxes", fontsize=11, fontweight="bold")
        fig.colorbar(im3, ax=ax3, fraction=0.046, pad=0.04)
        legend_handles = draw_wavelet_boxes(ax3, top_boxes, start_cell, end_cell)
        if legend_handles:
            ax3.legend(
                handles=legend_handles,
                title="Top Subregions (Diff Certainty):",
                title_fontsize=8.5,
                loc="lower center",
                bbox_to_anchor=(0.5, 1.08),
                ncol=len(legend_handles),
                fontsize=8.0,
                framealpha=0.9,
                edgecolor="black",
                columnspacing=0.8,
                handletextpad=0.4,
            )

        # Col 5: Empirical Diff + Boxes
        ax4 = axes[row_idx, 4]
        ax4.set_facecolor("#cccccc")
        im4 = ax4.imshow(raw_diff, cmap=cmap, vmin=-vmax_diff, vmax=vmax_diff, origin="upper")
        ax4.set_title(
            f"Diff: {c2} - {c1} (log $\\hat{{\\Lambda}}$ = {llr_val:.1f}, {p_clause})",
            fontsize=11,
            fontweight="bold",
        )
        fig.colorbar(im4, ax=ax4, fraction=0.046, pad=0.04)
        draw_wavelet_boxes(ax4, top_boxes, start_cell, end_cell)

        # Col 6: WaveQTL Z-score
        ax5 = axes[row_idx, 5]
        ax5.set_facecolor("#cccccc")
        z_hat_full = np.full((max_cell, max_cell), np.nan, dtype=np.float64)
        z_hat_full[start_cell : end_cell + 1, start_cell : end_cell + 1] = z_hat
        im5 = ax5.imshow(z_hat_full, cmap=cmap, vmin=-3.0, vmax=3.0, origin="upper")
        ax5.set_title("WaveQTL Z-score", fontsize=11, fontweight="bold")
        fig.colorbar(im5, ax=ax5, fraction=0.046, pad=0.04, ticks=[-3, -2, -1, 0, 1, 2, 3])

        # Col 7: WaveQTL Effect (Clean)
        ax6 = axes[row_idx, 6]
        ax6.set_facecolor("#cccccc")
        diff_hat_full = np.full((max_cell, max_cell), np.nan, dtype=np.float64)
        diff_hat_full[start_cell : end_cell + 1, start_cell : end_cell + 1] = diff_hat
        im6 = ax6.imshow(diff_hat_full, cmap=cmap, vmin=-vmax_diff, vmax=vmax_diff, origin="upper")
        ax6.set_title("WaveQTL Effect (Clean)", fontsize=11, fontweight="bold")
        fig.colorbar(im6, ax=ax6, fraction=0.046, pad=0.04)

        # Col 8: WaveQTL Effect + |Z| > 2 Contour
        ax7 = axes[row_idx, 7]
        ax7.set_facecolor("#cccccc")
        im7 = ax7.imshow(diff_hat_full, cmap=cmap, vmin=-vmax_diff, vmax=vmax_diff, origin="upper")
        z_mask_2 = np.abs(np.nan_to_num(z_hat_full, nan=0.0)) > 2.0
        draw_stepped_boundaries(ax7, z_mask_2, color="black", linewidth=1.8, zorder=4)
        frac_sig_2 = float(np.mean(np.abs(z_hat) > 2.0)) * 100.0
        ax7.set_title(f"WaveQTL Effect (|Z| > 2: {frac_sig_2:.1f}%)", fontsize=11, fontweight="bold")
        fig.colorbar(im7, ax=ax7, fraction=0.046, pad=0.04)

        # Col 9: WaveQTL Effect + |Z| > 3 Contour
        ax8 = axes[row_idx, 8]
        ax8.set_facecolor("#cccccc")
        im8 = ax8.imshow(diff_hat_full, cmap=cmap, vmin=-vmax_diff, vmax=vmax_diff, origin="upper")
        z_mask_3 = np.abs(np.nan_to_num(z_hat_full, nan=0.0)) > 3.0
        draw_stepped_boundaries(ax8, z_mask_3, color="black", linewidth=1.8, zorder=4)
        frac_sig_3 = float(np.mean(np.abs(z_hat) > 3.0)) * 100.0
        ax8.set_title(f"WaveQTL Effect (|Z| > 3: {frac_sig_3:.1f}%)", fontsize=11, fontweight="bold")
        fig.colorbar(im8, ax=ax8, fraction=0.046, pad=0.04)

        tick_step = max(5, 10 * max(1, max_cell // 60))
        for ax in (ax0, ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8):
            ax.set_xticks(np.arange(0, max_cell, tick_step))
            ax.set_yticks(np.arange(0, max_cell, tick_step))
            ax.tick_params(axis="both", which="major", labelsize=8.5)

    fig.suptitle(
        f"{target_name} — 2D WaveQTL Association & Subregion Localization (N = {n1} + {n2})",
        fontsize=16,
        fontweight="bold",
        y=1.02,
    )
    output_png.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_png, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved comparison figure -> {output_png.name}")


def plot_summary_ranking_chart(summary_df: pd.DataFrame, output_png: Path, dpi: int = 200) -> None:
    """Render 2-panel horizontal bar chart ranking all comparisons by logLR and Top-1 Domain BF."""
    if summary_df.empty:
        return

    n_bars = len(summary_df)
    fig_h = max(4.5, 0.45 * n_bars + 1.8)
    fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(16, fig_h), constrained_layout=True)

    df_sorted = summary_df.sort_values(by="log_lr", ascending=True).copy()
    y_pos = np.arange(len(df_sorted))
    labels = df_sorted["comp_id"].astype(str).tolist()

    sig_color = "#c93b2b"
    null_color = "#6c757d"
    bar_colors = [
        sig_color if (pd.notna(r["p_perm"]) and float(r["p_perm"]) < 0.05) else null_color
        for _, r in df_sorted.iterrows()
    ]

    bars = ax1.barh(y_pos, df_sorted["log_lr"], color=bar_colors, alpha=0.85, edgecolor="black")
    ax1.set_yticks(y_pos)
    ax1.set_yticklabels(labels, fontsize=10, fontweight="bold")
    ax1.set_xlabel(r"Overall Log Likelihood Ratio ($\log \hat{\Lambda}$)", fontsize=11, fontweight="bold")
    ax1.set_title("Overall Evidence of 2D Matrix Difference\n" + r"(Sorted by $\log \hat{\Lambda}$)", fontsize=12, fontweight="bold")
    ax1.grid(axis="x", linestyle="--", alpha=0.5)

    max_llr = max(1.0, float(df_sorted["log_lr"].max()))
    ax1.set_xlim(0, max_llr * 1.18)

    for bar, (_, r) in zip(bars, df_sorted.iterrows()):
        width = max(0.0, float(bar.get_width()))
        p_val = r["p_perm"]
        if pd.isna(p_val):
            p_text = "p=N/A"
        elif float(p_val) >= 0.001:
            p_text = f"p={float(p_val):.3f}"
        else:
            p_text = "p<0.001"
        ax1.text(width + max_llr * 0.015, bar.get_y() + bar.get_height() / 2, p_text, va="center", fontsize=8.5)

    ax2.barh(y_pos, df_sorted["top1_domain_bf"].fillna(0.0), color=bar_colors, alpha=0.85, edgecolor="black")
    ax2.set_yticks(y_pos)
    ax2.set_yticklabels([])
    ax2.set_xlabel(r"Top-1 Subregion $\log_{10}$ Bayes Factor ($\log_{10}$ BF)", fontsize=11, fontweight="bold")
    ax2.set_title("Strongest Localized Subregion Evidence\n" + r"(Dyadic Square Windows $\geq 4\times 4$ Bins)", fontsize=12, fontweight="bold")
    ax2.axvline(1.0, color="gray", linestyle=":", label=r"Strong Evidence ($\log_{10}$BF=1)")
    ax2.axvline(2.0, color="black", linestyle="--", label=r"Decisive Evidence ($\log_{10}$BF=2)")
    ax2.grid(axis="x", linestyle="--", alpha=0.5)

    legend_elements = [
        patches.Patch(facecolor=sig_color, edgecolor="black", label=r"Significant ($p_{\mathrm{perm}} < 0.05$)"),
        patches.Patch(facecolor=null_color, edgecolor="black", label=r"Not Significant ($p_{\mathrm{perm}} \geq 0.05$ or N/A)"),
    ]
    ax1.legend(handles=legend_elements, loc="lower right", fontsize=9.5, framealpha=0.9, edgecolor="black")
    ax2.legend(loc="lower right", fontsize=9)

    fig.suptitle("WaveQTL-2D Comparison Summary", fontsize=15, fontweight="bold")
    output_png.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_png, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved summary ranking chart -> {output_png.name}")


def summarize_and_visualize_all(
    work_dir: Path,
    output_dir: Path,
    full_matrix_tsv: Optional[Path] = None,
    full_cond_col: Optional[str] = None,
    row_col: str = "cell_i",
    col_col: str = "cell_j",
    val_col: str = "value",
    vmax: Optional[float] = None,
    dpi: int = 200,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Run full post-processing: reconstruct effects, write summary TSVs, and render all comparison plots."""
    summary_df, top_sub_df, valid_manifest = process_all_outputs(work_dir=work_dir, output_dir=output_dir)
    if summary_df.empty:
        print("No completed WaveQTL outputs found to visualize.")
        return summary_df, top_sub_df

    full_df = pd.read_csv(full_matrix_tsv, sep="\t") if full_matrix_tsv is not None else None

    grouped_targets: Dict[str, List[dict]] = OrderedDict()
    for comp_meta in valid_manifest:
        p_target = str(comp_meta.get("primary_target", comp_meta["comp_id"]))
        grouped_targets.setdefault(p_target, []).append(comp_meta)

    for p_target, metas in grouped_targets.items():
        safe_target = "".join(ch if (ch.isalnum() or ch in ("-", "_", ".")) else "_" for ch in p_target)
        out_png = output_dir / f"comparison_{safe_target}.png"
        plot_9col_group_figure(
            target_name=p_target,
            comp_metas=metas,
            work_dir=work_dir,
            output_png=out_png,
            full_df=full_df,
            full_cond_col=full_cond_col,
            row_col=row_col,
            col_col=col_col,
            val_col=val_col,
            vmax=vmax,
            dpi=dpi,
        )

    plot_summary_ranking_chart(summary_df, output_dir / "summary_ranking.png", dpi=dpi)
    return summary_df, top_sub_df
