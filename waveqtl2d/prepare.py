"""
prepare.py -- Transform per-sample 2D matrices into WaveQTL input files.

For each comparison group:
  1. Aligns both conditions to their shared active coordinate range [start_cell, end_cell].
  2. Builds symmetric 2D grids per sample and zero-pads to dyadic size 2^J x 2^J (min 64x64).
  3. Applies a J-level 2D Haar Discrete Wavelet Transform (DWT).
  4. Exports WaveQTL inputs into <work_dir>/<comparison_id>/:
       - WCs.txt: (N1 + N2) x 4^J matrix of raw wavelet coefficients
       - geno.txt: binary condition indicator (N1 zeros, N2 ones)
       - use.txt: binary flags masking symmetric redundancy, zero-padded margins, and zero-variance tiles
       - group.txt: multiscale scale boundaries
       - coeff_metadata.tsv: spatial bounding box and symmetry mapping per coefficient
       - sample_matrices.tsv.gz: standardized sample cells for visualization
       - comparison_meta.json: coordinate bounds, sample counts, and condition labels
"""

import json
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import pywt


def _sort_key_numeric_first(val) -> Tuple[int, object]:
    """Sort integer-like sample IDs numerically, falling back to string order."""
    try:
        return (0, int(val))
    except (ValueError, TypeError):
        return (1, str(val))


def sanitize_name(parts: Sequence[object]) -> str:
    """Create a clean filesystem-safe comparison directory identifier."""
    raw = "_".join(str(p) for p in parts)
    return "".join(ch if (ch.isalnum() or ch in ("-", "_", ".")) else "_" for ch in raw)


def compute_dyadic_grid_size(grid_dim: int, min_grid_size: int = 64) -> Tuple[int, int]:
    """Return (grid_size, decomp_levels) for the smallest 2^J >= max(grid_dim, min_grid_size)."""
    target = max(min_grid_size, 1 << max(1, (grid_dim - 1).bit_length()))
    levels = int(round(np.log2(target)))
    return target, levels


def get_coeff_metadata_and_groups(grid_size: int, decomp_levels: int) -> Tuple[List[dict], List[int]]:
    """
    Build 4^J coefficient metadata (scale, subband, bounding box, symmetry mirror)
    and 1-based WaveQTL scale group start indices.
    """
    metadata: List[dict] = [
        {
            "scale": 0,
            "subband": "A",
            "u": 0,
            "v": 0,
            "r0": 0,
            "r1": grid_size,
            "c0": 0,
            "c1": grid_size,
            "is_redundant_symmetry": False,
            "mirror_key": (0, "A", 0, 0),
        }
    ]

    for lvl in range(1, decomp_levels + 1):
        s_dim = 2 ** (lvl - 1)
        block = grid_size // s_dim
        for sub in ("H", "V", "D"):
            for u in range(s_dim):
                for v in range(s_dim):
                    r0, r1 = u * block, (u + 1) * block
                    c0, c1 = v * block, (v + 1) * block

                    is_redundant = False
                    if sub == "V":
                        # Vertical detail subband is the transpose of horizontal detail subband
                        is_redundant = True
                        mirror_key = (lvl, "H", v, u)
                    elif sub == "D":
                        # Diagonal subband is symmetric; mask strict lower triangle (u > v)
                        if u > v:
                            is_redundant = True
                            mirror_key = (lvl, "D", v, u)
                        else:
                            mirror_key = (lvl, "D", u, v)
                    else:
                        mirror_key = (lvl, "V", v, u)

                    metadata.append(
                        {
                            "scale": lvl,
                            "subband": sub,
                            "u": u,
                            "v": v,
                            "r0": r0,
                            "r1": r1,
                            "c0": c0,
                            "c1": c1,
                            "is_redundant_symmetry": is_redundant,
                            "mirror_key": mirror_key,
                        }
                    )

    key_to_idx = {
        (m["scale"], m["subband"], m["u"], m["v"]): idx + 1
        for idx, m in enumerate(metadata)
    }
    for idx, m in enumerate(metadata):
        m["coeff_idx"] = idx + 1
        m["mirror_coeff_idx"] = key_to_idx[m["mirror_key"]]

    group_starts = [1] + [4 ** (lvl - 1) + 1 for lvl in range(1, decomp_levels + 1)]
    return metadata, group_starts


def flatten_coeffs(coeffs: list) -> np.ndarray:
    """Flatten pywt.wavedec2 output into a 1D vector matching WaveQTL scale ordering."""
    flat = [coeffs[0].ravel()]
    for lvl in range(1, len(coeffs)):
        c_h, c_v, c_d = coeffs[lvl]
        flat.append(c_h.ravel())
        flat.append(c_v.ravel())
        flat.append(c_d.ravel())
    return np.concatenate(flat)


def compute_shared_cell_range(
    df_sub: pd.DataFrame,
    cond_col: str,
    cond1: str,
    cond2: str,
    row_col: str,
    col_col: str,
    val_col: str,
) -> Tuple[int, int]:
    """Compute the shared [start_cell, end_cell] intersection observed in both conditions."""
    valid = df_sub[df_sub[val_col].notna()]
    df_c1 = valid[valid[cond_col].astype(str) == str(cond1)]
    df_c2 = valid[valid[cond_col].astype(str) == str(cond2)]
    if df_c1.empty or df_c2.empty:
        raise ValueError(f"Missing non-null observations for '{cond1}' or '{cond2}'.")

    min_c1 = int(min(df_c1[row_col].min(), df_c1[col_col].min()))
    max_c1 = int(max(df_c1[row_col].max(), df_c1[col_col].max()))
    min_c2 = int(min(df_c2[row_col].min(), df_c2[col_col].min()))
    max_c2 = int(max(df_c2[row_col].max(), df_c2[col_col].max()))

    start_cell = max(min_c1, min_c2)
    end_cell = min(max_c1, max_c2)
    if start_cell > end_cell:
        raise ValueError(
            f"Empty coordinate intersection: {cond1}=[{min_c1}, {max_c1}] vs {cond2}=[{min_c2}, {max_c2}]."
        )
    return start_cell, end_cell


def build_sample_grid(
    sample_df: pd.DataFrame,
    start_cell: int,
    end_cell: int,
    grid_size: int,
    row_col: str,
    col_col: str,
    val_col: str,
) -> np.ndarray:
    """Build a zero-padded (grid_size x grid_size) symmetric matrix for one sample."""
    grid = np.zeros((grid_size, grid_size), dtype=np.float64)
    if sample_df.empty:
        return grid

    rows = sample_df[row_col].to_numpy(dtype=np.int64)
    cols = sample_df[col_col].to_numpy(dtype=np.int64)
    vals = sample_df[val_col].to_numpy(dtype=np.float64)

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

    grid[gi, gj] = v
    grid[gj, gi] = v
    return grid


def resolve_conditions(
    df: pd.DataFrame,
    cond_col: str,
    cond1: Optional[str] = None,
    cond2: Optional[str] = None,
) -> Tuple[str, str]:
    """Determine the two condition labels (cond1 -> 0.0, cond2 -> 1.0)."""
    seen = [str(x) for x in pd.unique(df[cond_col].dropna())]
    if cond1 is not None and cond2 is not None:
        if str(cond1) not in seen or str(cond2) not in seen:
            raise ValueError(f"Requested conditions ({cond1}, {cond2}) not found in '{cond_col}': {seen}")
        return str(cond1), str(cond2)
    if len(seen) != 2:
        raise ValueError(
            f"Expected 2 unique conditions in '{cond_col}', found {len(seen)}: {seen}. "
            f"Please specify --cond1 and --cond2 explicitly."
        )
    return seen[0], seen[1]


def prepare_all_inputs(
    input_tsv: Path,
    work_dir: Path,
    group_by: Optional[Sequence[str]] = None,
    cond_col: str = "condition",
    sample_col: str = "sample",
    row_col: str = "cell_i",
    col_col: str = "cell_j",
    val_col: str = "value",
    cond1: Optional[str] = None,
    cond2: Optional[str] = None,
    min_grid_size: int = 64,
) -> List[dict]:
    """Read input sample 2D matrix TSV(.gz) and prepare WaveQTL directories for all comparisons."""
    df = pd.read_csv(input_tsv, sep="\t")

    group_cols = list(group_by) if group_by else []
    required_cols = group_cols + [cond_col, sample_col, row_col, col_col, val_col]
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns {missing} in {input_tsv}. Available: {list(df.columns)}")

    c1_label, c2_label = resolve_conditions(df, cond_col=cond_col, cond1=cond1, cond2=cond2)
    df = df[df[cond_col].astype(str).isin([c1_label, c2_label])].copy()

    work_dir.mkdir(parents=True, exist_ok=True)

    # Pre-compute shared [start_cell, end_cell] per primary target (group_cols[0]) when multi-key
    # so sub-comparisons of the same target share an aligned coordinate frame.
    primary_ranges: Dict[str, Tuple[int, int]] = {}
    if len(group_cols) >= 1:
        p_col = group_cols[0]
        for p_val, p_df in df.groupby(p_col, sort=False):
            primary_ranges[str(p_val)] = compute_shared_cell_range(
                p_df, cond_col, c1_label, c2_label, row_col, col_col, val_col
            )

    if group_cols:
        grouped_iter = list(df.groupby(group_cols, sort=False))
    else:
        grouped_iter = [(("comparison",), df)]

    metadata_cache: Dict[int, Tuple[List[dict], List[int]]] = {}
    manifest: List[dict] = []

    for key_vals, sub_df in grouped_iter:
        if not isinstance(key_vals, tuple):
            key_tuple = (key_vals,)
        else:
            key_tuple = key_vals

        comp_id = sanitize_name(key_tuple)
        comp_dir = work_dir / comp_id
        comp_dir.mkdir(parents=True, exist_ok=True)

        if group_cols:
            start_cell, end_cell = primary_ranges[str(key_tuple[0])]
        else:
            start_cell, end_cell = compute_shared_cell_range(
                sub_df, cond_col, c1_label, c2_label, row_col, col_col, val_col
            )

        grid_dim = end_cell - start_cell + 1
        grid_size, decomp_levels = compute_dyadic_grid_size(grid_dim, min_grid_size=min_grid_size)

        if grid_size not in metadata_cache:
            metadata_cache[grid_size] = get_coeff_metadata_and_groups(grid_size, decomp_levels)
        metadata_list, group_starts = metadata_cache[grid_size]

        sub_c1 = sub_df[sub_df[cond_col].astype(str) == c1_label]
        sub_c2 = sub_df[sub_df[cond_col].astype(str) == c2_label]

        samples_c1 = sorted(pd.unique(sub_c1[sample_col]), key=_sort_key_numeric_first)
        samples_c2 = sorted(pd.unique(sub_c2[sample_col]), key=_sort_key_numeric_first)
        if not samples_c1 or not samples_c2:
            print(f"Skipping {comp_id}: missing samples in one condition ({len(samples_c1)} vs {len(samples_c2)}).")
            continue

        wcs_rows: List[np.ndarray] = []
        standardized_chunks: List[pd.DataFrame] = []

        for cond_label, cond_df, sample_list in (
            (c1_label, sub_c1, samples_c1),
            (c2_label, sub_c2, samples_c2),
        ):
            for s_idx, s_id in enumerate(sample_list, start=1):
                s_df = cond_df[cond_df[sample_col] == s_id]
                grid = build_sample_grid(
                    s_df, start_cell, end_cell, grid_size, row_col, col_col, val_col
                )
                coeffs = pywt.wavedec2(grid, wavelet="haar", level=decomp_levels)
                wcs_rows.append(flatten_coeffs(coeffs))

                # Save standardized sample records for downstream plotting
                chunk = pd.DataFrame(
                    {
                        "condition": cond_label,
                        "sample_order": s_idx,
                        "sample_id": str(s_id),
                        "cell_i": s_df[row_col].to_numpy(),
                        "cell_j": s_df[col_col].to_numpy(),
                        "value": s_df[val_col].to_numpy(),
                    }
                )
                if "n_obs" in s_df.columns:
                    chunk["n_obs"] = s_df["n_obs"].to_numpy()
                elif "n_reads" in s_df.columns:
                    chunk["n_obs"] = s_df["n_reads"].to_numpy()
                standardized_chunks.append(chunk)

        wcs_mat = np.array(wcs_rows, dtype=np.float64)
        var_wcs = np.var(wcs_mat, axis=0)

        use_flags = []
        for m, var_val in zip(metadata_list, var_wcs):
            if m["is_redundant_symmetry"]:
                use_flags.append(0)
            elif m["r0"] >= grid_dim or m["c0"] >= grid_dim:
                use_flags.append(0)
            elif var_val < 1e-12:
                use_flags.append(0)
            else:
                use_flags.append(1)
        use_arr = np.array(use_flags, dtype=np.int32)

        # 1. WCs.txt
        np.savetxt(comp_dir / "WCs.txt", wcs_mat, fmt="%.8f", delimiter=" ")

        # 2. use.txt
        with open(comp_dir / "use.txt", "w", encoding="utf-8") as f:
            f.write(" ".join(map(str, use_arr)) + "\n")

        # 3. group.txt
        with open(comp_dir / "group.txt", "w", encoding="utf-8") as f:
            f.write(" ".join(map(str, group_starts)) + "\n")

        # 4. geno.txt (BIMBAM mean covariate format: cond1=0.0, cond2=1.0)
        geno_vals = ["0.0"] * len(samples_c1) + ["1.0"] * len(samples_c2)
        with open(comp_dir / "geno.txt", "w", encoding="utf-8") as f:
            f.write(" ".join(["COND", "A", "T"] + geno_vals) + "\n")

        # 5. coeff_metadata.tsv
        meta_records = []
        for m, use_f, var_val in zip(metadata_list, use_arr, var_wcs):
            rec = {k: v for k, v in m.items() if k != "mirror_key"}
            rec["use_flag"] = int(use_f)
            rec["var_wcs"] = float(var_val)
            rec["coord_r0"] = start_cell + m["r0"]
            rec["coord_r1"] = min(end_cell + 1, start_cell + m["r1"])
            rec["coord_c0"] = start_cell + m["c0"]
            rec["coord_c1"] = min(end_cell + 1, start_cell + m["c1"])
            meta_records.append(rec)
        pd.DataFrame(meta_records).to_csv(comp_dir / "coeff_metadata.tsv", sep="\t", index=False)

        # 6. sample_matrices.tsv.gz
        pd.concat(standardized_chunks, ignore_index=True).to_csv(
            comp_dir / "sample_matrices.tsv.gz", sep="\t", index=False, compression="gzip"
        )

        # 7. comparison_meta.json
        group_dict = {col: str(val) for col, val in zip(group_cols, key_tuple)}
        comp_meta = {
            "comp_id": comp_id,
            "group_cols": group_cols,
            "group_values": group_dict,
            "primary_target": str(key_tuple[0]) if group_cols else comp_id,
            "sub_label": " / ".join(str(x) for x in key_tuple[1:]) if len(group_cols) > 1 else "",
            "cond1": c1_label,
            "cond2": c2_label,
            "n_cond1": len(samples_c1),
            "n_cond2": len(samples_c2),
            "start_cell": int(start_cell),
            "end_cell": int(end_cell),
            "grid_dim": int(grid_dim),
            "grid_size": int(grid_size),
            "decomp_levels": int(decomp_levels),
            "num_coeffs": int(grid_size * grid_size),
            "active_coeffs": int(np.sum(use_arr)),
        }
        with open(comp_dir / "comparison_meta.json", "w", encoding="utf-8") as f:
            json.dump(comp_meta, f, indent=2)

        manifest.append(comp_meta)
        print(
            f"[{comp_id}] Prepared: cells=[{start_cell}..{end_cell}] ({grid_dim}x{grid_dim} -> {grid_size}x{grid_size}), "
            f"samples={len(samples_c1)}+{len(samples_c2)}, active WCs={comp_meta['active_coeffs']}/{comp_meta['num_coeffs']}"
        )

    with open(work_dir / "manifest.json", "w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2)

    return manifest
