"""
execute.py -- Locate and execute the compiled WaveQTL C++ binary on prepared comparisons.
"""

import json
from pathlib import Path
import shutil
import subprocess
import time
from typing import List, Optional


def find_waveqtl_binary(custom_path: Optional[Path] = None) -> Path:
    """Locate compiled WaveQTL executable using simple fallback checks."""
    if custom_path is not None:
        resolved = Path(custom_path).expanduser().resolve()
        if resolved.exists():
            return resolved
        raise FileNotFoundError(f"Specified --waveqtl-bin does not exist: {resolved}")

    # 1. Packaged binary inside installed waveqtl2d/bin/WaveQTL
    pkg_bin = Path(__file__).resolve().parent / "bin" / "WaveQTL"
    if pkg_bin.exists():
        return pkg_bin.resolve()

    # 2. Local development checkout (repository root)
    repo_root = Path(__file__).resolve().parent.parent
    if (repo_root / "WaveQTL").exists():
        return (repo_root / "WaveQTL").resolve()
    if (repo_root / "bin" / "WaveQTL").exists():
        return (repo_root / "bin" / "WaveQTL").resolve()

    # 3. Environment variable WAVEQTL_BIN
    import os
    env_bin = os.environ.get("WAVEQTL_BIN")
    if env_bin and Path(env_bin).exists():
        return Path(env_bin).resolve()

    # 4. PATH lookup
    on_path = shutil.which("WaveQTL")
    if on_path:
        return Path(on_path).resolve()

    raise FileNotFoundError(
        "WaveQTL binary not found in package (waveqtl2d/bin/WaveQTL), repo root, or PATH. "
        "Please compile it first by running `make all` in the WaveQTL-2D root directory."
    )


def run_single_comparison(
    comp_dir: Path,
    num_perm: int = 1000,
    waveqtl_bin: Optional[Path] = None,
) -> dict:
    """Execute WaveQTL for a single prepared comparison directory."""
    start_t = time.time()
    out_dir = comp_dir / "output"
    out_dir.mkdir(parents=True, exist_ok=True)

    meta_file = comp_dir / "comparison_meta.json"
    if meta_file.exists():
        with open(meta_file, "r", encoding="utf-8") as f:
            comp_meta = json.load(f)
        num_coeffs = int(comp_meta["num_coeffs"])
    else:
        with open(comp_dir / "use.txt", "r", encoding="utf-8") as f:
            num_coeffs = len(f.readline().split())

    loglr_file = out_dir / "waveqtl_out.fph.logLR.txt"
    mean_file = out_dir / "waveqtl_out.fph.mean.txt"
    pval_file = out_dir / "waveqtl_out.fph.pval.txt"

    # Remove prior output markers so completion check reflects the current run
    for old_f in (loglr_file, mean_file, pval_file):
        if old_f.exists():
            old_f.unlink()

    bin_path = find_waveqtl_binary(waveqtl_bin)
    cmd = [
        str(bin_path),
        "-gmode", "1",
        "-notsnp",
        "-g", "geno.txt",
        "-p", "WCs.txt",
        "-u", "use.txt",
        "-group", "group.txt",
        "-o", "waveqtl_out",
        "-f", str(num_coeffs),
        "-fph", "2" if num_perm > 0 else "1",
    ]
    if num_perm > 0:
        cmd.extend(["-numPerm", str(num_perm)])

    proc = subprocess.run(
        cmd,
        cwd=comp_dir.resolve(),
        capture_output=True,
        text=True,
    )
    elapsed = time.time() - start_t

    # Note: upstream WaveQTL's main() returns 1 on normal completion (src/fp.cpp),
    # so verify fresh output files rather than proc.returncode == 0.
    success = loglr_file.exists() and mean_file.exists() and (loglr_file.stat().st_size > 0)
    loglr_val = None
    pval_val = None

    if success:
        with open(loglr_file, "r", encoding="utf-8") as f:
            parts = f.readline().strip().split()
            if len(parts) >= 2:
                loglr_val = float(parts[1])
        if num_perm > 0 and pval_file.exists():
            with open(pval_file, "r", encoding="utf-8") as f:
                lines = [line.strip() for line in f.readlines() if line.strip()]
                if len(lines) >= 4:
                    pval_val = float(lines[3].split()[0])

    return {
        "comp_id": comp_dir.name,
        "success": success,
        "returncode": proc.returncode,
        "elapsed_sec": elapsed,
        "logLR": loglr_val,
        "pval": pval_val,
        "stderr": proc.stderr[-500:] if proc.stderr else "",
    }


def execute_all_comparisons(
    work_dir: Path,
    num_perm: int = 1000,
    waveqtl_bin: Optional[Path] = None,
) -> List[dict]:
    """Run WaveQTL across all prepared comparison folders in work_dir."""
    manifest_path = work_dir / "manifest.json"
    if manifest_path.exists():
        with open(manifest_path, "r", encoding="utf-8") as f:
            comp_dirs = [work_dir / item["comp_id"] for item in json.load(f)]
    else:
        comp_dirs = sorted(p for p in work_dir.iterdir() if p.is_dir() and (p / "WCs.txt").exists())

    results = []
    for idx, comp_dir in enumerate(comp_dirs, start=1):
        res = run_single_comparison(comp_dir, num_perm=num_perm, waveqtl_bin=waveqtl_bin)
        results.append(res)
        if res["success"]:
            p_str = f", p_perm={res['pval']:.4f}" if res["pval"] is not None else ""
            print(
                f"[{idx}/{len(comp_dirs)}] {comp_dir.name}: OK ({res['elapsed_sec']:.1f}s) "
                f"| logLR={res['logLR']:.3f}{p_str}"
            )
        else:
            print(f"[{idx}/{len(comp_dirs)}] {comp_dir.name}: FAILED (code={res['returncode']}) {res['stderr']}")
    return results
