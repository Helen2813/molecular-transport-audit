from __future__ import annotations

import gc
import gzip
import json
import math
import re
import time
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

try:
    import torch
except Exception as exc:
    raise RuntimeError("PyTorch is required for 04h. Run the CUDA audit first.") from exc

try:
    from scipy.stats import rankdata
except Exception as exc:
    raise RuntimeError("scipy is required for exact Spearman ranking in 04h.") from exc


SCRIPT_VERSION = "04h-run-final-pretarget-audits-gpu-v1-no-cli"

DATA_ROOT = Path(r"D:\paper4_tcbb_data")

CONTRACT = (
    DATA_ROOT
    / "paper4_tcbb_pretarget_confounds_contract_v1"
    / "pretarget_confounds_contract_v1.json"
)

TCGA_EXPR = (
    DATA_ROOT
    / "paper4_tcbb_input_audit_v1"
    / "staged_continuous_inputs"
    / "TCGA_BRCA_PanCanAtlas2018"
    / "data_mrna_seq_v2_rsem.txt"
)
SCANB_EXPR = (
    DATA_ROOT
    / "SCANB_GSE96058"
    / "GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv.gz"
)
METABRIC_EXPR = (
    DATA_ROOT
    / "paper4_tcbb_input_audit_v1"
    / "staged_continuous_inputs"
    / "METABRIC"
    / "data_mrna_illumina_microarray.txt"
)

UNIVERSE = (
    DATA_ROOT
    / "paper4_tcbb_tcga_source_universe_v1"
    / "tcga_source_gene_universe_frozen_v1.tsv"
)
MEMBERSHIP = (
    DATA_ROOT
    / "paper4_tcbb_tcga_source_modules_v1"
    / "tcga_source_module_membership_frozen_v1.tsv"
)
PROGRAM_DIR = DATA_ROOT / "paper4_tcbb_frozen_source_programs_v1"
WEIGHTS = PROGRAM_DIR / "tcga_frozen_source_program_weights_v1.tsv"
PROGRAM_SUMMARY = PROGRAM_DIR / "tcga_frozen_source_program_summary_v1.tsv"

PAIR_DIR = DATA_ROOT / "paper4_tcbb_scanb_pairing_contract_v1"
TECH_PAIRS = PAIR_DIR / "scanb_technical_replicate_pairs_v1.tsv"
PRIMARY = PAIR_DIR / "scanb_primary_profiles_frozen_v1.tsv"

SCANB_PAM50 = (
    DATA_ROOT
    / "paper4_tcbb_pam50_alignment_audit_v1"
    / "scanb_frozen_primary_pam50_alignment_v1.tsv"
)

GTF = (
    DATA_ROOT
    / "SCANB_GSE96058"
    / "GSE96058_UCSC_hg38_knownGenes_22sep2014.gtf.gz"
)

OUT_DIR = DATA_ROOT / "paper4_tcbb_final_pretarget_audits_v1"

BASE_SEED = 20260918
TECH_NULL_NONMATCH_PAIRS = 100_000
IDENTITY_MIN_MARGIN = 0.02
SOURCE_EDGE_BOOTSTRAPS = 500
SOURCE_LOADING_BOOTSTRAPS = 200
SOURCE_MAX_EDGES = 20_000

SOURCE_RECON_TOL = 1e-7
GPU_CPU_TOL = 2e-10
LOADING_FULL_RESIDUAL_TOL = 1e-7
LOADING_BOOT_RESIDUAL_TOL = 1e-7

CANONICAL_CHROMS = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}


def require(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Required file does not exist: {path}")


def canon_symbol(x: object) -> str:
    s = "" if x is None else str(x)
    return " ".join(s.strip().split()).upper()


def ensure_cuda() -> torch.device:
    if not torch.cuda.is_available():
        raise RuntimeError("CUDA is not available in the active Python environment.")
    torch.set_default_dtype(torch.float64)
    return torch.device("cuda:0")


def bh_not_used_notice() -> None:
    # Explicit placeholder to emphasize that this script is not a preservation test.
    return None


def spearman_fast(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    if a.size != b.size or a.size < 3:
        return float("nan")
    ra = rankdata(a, method="average")
    rb = rankdata(b, method="average")
    ra -= ra.mean()
    rb -= rb.mean()
    den = math.sqrt(float(np.dot(ra, ra) * np.dot(rb, rb)))
    if den == 0.0:
        return float("nan")
    return float(np.dot(ra, rb) / den)


def quantile(x: np.ndarray, q: float) -> float:
    return float(np.quantile(np.asarray(x, dtype=np.float64), q, method="linear"))


def zscore_genes_samples(
    gene_by_sample: np.ndarray,
    reference_sample_indices: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    gene_by_sample -> gene-wise z scores.
    Means/SDs can be learned on biological primary samples and then applied to all profiles.
    Sample SD (ddof=1), matching the R scale/correlation convention.
    """
    x = np.asarray(gene_by_sample, dtype=np.float64)
    ref = x if reference_sample_indices is None else x[:, reference_sample_indices]
    mu = np.mean(ref, axis=1)
    sd = np.std(ref, axis=1, ddof=1)
    valid = np.isfinite(mu) & np.isfinite(sd) & (sd > 0)
    z = np.full_like(x, np.nan, dtype=np.float64)
    z[valid, :] = (x[valid, :] - mu[valid, None]) / sd[valid, None]
    return z, valid, sd


def torch_profile_matrix(
    cohort: dict,
    genes: list[str],
    device: torch.device,
) -> torch.Tensor:
    idx = np.fromiter((cohort["gene_index"][g] for g in genes), dtype=np.int64)
    x = cohort["z"][idx, :].T  # samples x genes
    if not np.isfinite(x).all():
        raise RuntimeError(f"Non-finite standardized values in {cohort['name']} fingerprint matrix.")
    t = torch.as_tensor(x, dtype=torch.float64, device=device)
    t = t - t.mean(dim=1, keepdim=True)
    norm = torch.linalg.vector_norm(t, dim=1, keepdim=True)
    if torch.any(norm <= 0):
        raise RuntimeError(f"Zero-norm sample fingerprint in {cohort['name']}.")
    return t / norm


def cpu_profile_matrix(cohort: dict, genes: list[str]) -> np.ndarray:
    idx = np.fromiter((cohort["gene_index"][g] for g in genes), dtype=np.int64)
    x = cohort["z"][idx, :].T.copy()
    x -= x.mean(axis=1, keepdims=True)
    norm = np.linalg.norm(x, axis=1, keepdims=True)
    if np.any(norm <= 0):
        raise RuntimeError(f"Zero-norm CPU fingerprint in {cohort['name']}.")
    return x / norm


def load_source_tcga(
    universe: pd.DataFrame,
) -> dict:
    print("\n[1/7] Reconstructing the exact frozen TCGA 10,000-gene source matrix ...")
    t0 = time.perf_counter()
    raw = pd.read_csv(TCGA_EXPR, sep="\t", low_memory=False)
    if raw.shape[1] != 1084:
        raise RuntimeError(f"Unexpected TCGA expression shape: {raw.shape}")

    row_idx = universe["source_row_index_0based"].astype(int).to_numpy()
    selected = raw.iloc[row_idx, :].copy()

    expected = [canon_symbol(x) for x in universe["Hugo_Symbol"]]
    observed = [canon_symbol(x) for x in selected["Hugo_Symbol"]]
    if expected != observed:
        for i, (e, o) in enumerate(zip(expected, observed)):
            if e != o:
                raise RuntimeError(
                    f"TCGA source-universe row replay mismatch at rank {i+1}: expected {e}, got {o}"
                )
        raise RuntimeError("TCGA source-universe row replay mismatch.")

    sample_ids = list(selected.columns[2:])
    vals = selected.iloc[:, 2:].to_numpy(dtype=np.float64)
    if np.any(vals < 0) or not np.isfinite(vals).all():
        raise RuntimeError("TCGA RSEM replay contains negative or non-finite values.")

    transformed = np.log2(vals + 1.0)
    z, valid, _ = zscore_genes_samples(transformed)
    if not valid.all():
        bad = [expected[i] for i in np.where(~valid)[0][:20]]
        raise RuntimeError(f"Frozen TCGA genes unexpectedly have zero/nonfinite SD: {bad}")

    del raw, selected, vals
    gc.collect()

    print(
        f"  replayed: {transformed.shape[0]:,} genes x {transformed.shape[1]:,} samples "
        f"in {time.perf_counter()-t0:.1f}s"
    )
    return {
        "name": "TCGA_BRCA",
        "genes": expected,
        "gene_index": {g: i for i, g in enumerate(expected)},
        "samples": sample_ids,
        "sample_index": {s: i for i, s in enumerate(sample_ids)},
        "raw_transformed": transformed,
        "z": z,
    }


def load_scanb(universe_set: set[str], primary_titles: list[str]) -> dict:
    print("\n[2/7] Loading SCAN-B frozen-universe genes and replaying the frozen transform ...")
    t0 = time.perf_counter()
    parts = []
    sample_cols = None

    for chunk in pd.read_csv(
        SCANB_EXPR,
        compression="gzip",
        chunksize=1500,
        low_memory=False,
    ):
        symbol_col = chunk.columns[0]
        if sample_cols is None:
            sample_cols = list(chunk.columns[1:])
        symbols = chunk[symbol_col].map(canon_symbol)
        mask = symbols.isin(universe_set)
        if not mask.any():
            continue

        sub = chunk.loc[mask, [symbol_col] + sample_cols].copy()
        sub["__SYMBOL__"] = symbols.loc[mask].to_numpy()
        arr = sub[sample_cols].to_numpy(dtype=np.float64)

        # Published SCAN-B values are log2(FPKM + 0.1).
        fpkm = np.maximum(np.exp2(arr) - 0.1, 0.0)
        arr2 = np.log2(fpkm + 1.0)

        out = pd.DataFrame(arr2, columns=sample_cols)
        out.insert(0, "__SYMBOL__", sub["__SYMBOL__"].to_numpy())
        parts.append(out)

    if not parts:
        raise RuntimeError("No frozen-universe genes found in SCAN-B.")

    df = pd.concat(parts, axis=0, ignore_index=True)
    df = df.groupby("__SYMBOL__", sort=False, as_index=True).mean(numeric_only=True)

    missing_primary = [s for s in primary_titles if s not in df.columns]
    if missing_primary:
        raise RuntimeError(f"SCAN-B frozen primary columns missing: {missing_primary[:20]}")

    primary_idx = np.array([df.columns.get_loc(s) for s in primary_titles], dtype=np.int64)
    genes = list(df.index)
    values = df.to_numpy(dtype=np.float64)
    z, valid, _ = zscore_genes_samples(values, reference_sample_indices=primary_idx)

    genes_valid = [g for g, ok in zip(genes, valid) if ok]
    z = z[valid, :]
    values = values[valid, :]

    print(
        f"  mapped unique frozen-universe genes with nonzero primary SD: {len(genes_valid):,}; "
        f"profiles: {z.shape[1]:,}; elapsed {time.perf_counter()-t0:.1f}s"
    )
    return {
        "name": "SCANB_GSE96058",
        "genes": genes_valid,
        "gene_index": {g: i for i, g in enumerate(genes_valid)},
        "samples": list(df.columns),
        "sample_index": {s: i for i, s in enumerate(df.columns)},
        "raw_transformed": values,
        "z": z,
    }


def load_metabric(universe_set: set[str]) -> dict:
    print("\n[3/7] Loading METABRIC frozen-universe genes on the frozen continuous scale ...")
    t0 = time.perf_counter()
    parts = []
    sample_cols = None

    for chunk in pd.read_csv(
        METABRIC_EXPR,
        sep="\t",
        chunksize=1500,
        low_memory=False,
    ):
        if "Hugo_Symbol" not in chunk.columns:
            raise RuntimeError("METABRIC chunk lacks Hugo_Symbol.")
        if sample_cols is None:
            sample_cols = list(chunk.columns[2:])

        symbols = chunk["Hugo_Symbol"].map(canon_symbol)
        mask = symbols.isin(universe_set)
        if not mask.any():
            continue

        arr = chunk.loc[mask, sample_cols].to_numpy(dtype=np.float64)
        out = pd.DataFrame(arr, columns=sample_cols)
        out.insert(0, "__SYMBOL__", symbols.loc[mask].to_numpy())
        parts.append(out)

    if not parts:
        raise RuntimeError("No frozen-universe genes found in METABRIC.")

    df = pd.concat(parts, axis=0, ignore_index=True)
    # Frozen target rule: average duplicate symbols on the transformed/as-provided scale.
    df = df.groupby("__SYMBOL__", sort=False, as_index=True).mean(numeric_only=True)

    genes = list(df.index)
    values = df.to_numpy(dtype=np.float64)
    z, valid, _ = zscore_genes_samples(values)

    genes_valid = [g for g, ok in zip(genes, valid) if ok]
    z = z[valid, :]
    values = values[valid, :]

    print(
        f"  mapped unique frozen-universe genes with nonzero SD: {len(genes_valid):,}; "
        f"samples: {z.shape[1]:,}; elapsed {time.perf_counter()-t0:.1f}s"
    )
    return {
        "name": "METABRIC",
        "genes": genes_valid,
        "gene_index": {g: i for i, g in enumerate(genes_valid)},
        "samples": list(df.columns),
        "sample_index": {s: i for i, s in enumerate(df.columns)},
        "raw_transformed": values,
        "z": z,
    }


def build_negative_cross_platform_pairs(
    tech_pairs: pd.DataFrame,
    primary_manifest: pd.DataFrame,
    scanb: dict,
    rng: np.random.Generator,
) -> tuple[np.ndarray, np.ndarray]:
    pri_platform = dict(
        zip(primary_manifest["primary_title"], primary_manifest["primary_platform_id"])
    )
    primary_titles = list(primary_manifest["primary_title"])
    by_platform = defaultdict(list)
    for title in primary_titles:
        by_platform[pri_platform[title]].append(title)

    candidates_a = []
    candidates_b = []

    for row in tech_pairs.itertuples(index=False):
        rep_title = row.replicate_title
        rep_platform = row.replicate_platform_id
        true_primary = row.primary_title

        for platform, titles in by_platform.items():
            if platform == rep_platform:
                continue
            for pri_title in titles:
                if pri_title == true_primary:
                    continue
                candidates_a.append(scanb["sample_index"][pri_title])
                candidates_b.append(scanb["sample_index"][rep_title])

    total = len(candidates_a)
    if total < TECH_NULL_NONMATCH_PAIRS:
        raise RuntimeError(
            f"Only {total:,} unique cross-platform nonmatching primary/replicate pairs; "
            f"contract requires {TECH_NULL_NONMATCH_PAIRS:,}."
        )

    choose = rng.choice(total, size=TECH_NULL_NONMATCH_PAIRS, replace=False)
    a = np.asarray(candidates_a, dtype=np.int64)[choose]
    b = np.asarray(candidates_b, dtype=np.int64)[choose]
    return a, b


def paired_profile_dots(
    profile_t: torch.Tensor,
    idx_a: np.ndarray,
    idx_b: np.ndarray,
    chunk: int = 2000,
) -> np.ndarray:
    vals = np.empty(len(idx_a), dtype=np.float64)
    for start in range(0, len(idx_a), chunk):
        end = min(start + chunk, len(idx_a))
        ia = torch.as_tensor(idx_a[start:end], dtype=torch.long, device=profile_t.device)
        ib = torch.as_tensor(idx_b[start:end], dtype=torch.long, device=profile_t.device)
        vals[start:end] = (
            (profile_t.index_select(0, ia) * profile_t.index_select(0, ib))
            .sum(dim=1)
            .detach()
            .cpu()
            .numpy()
        )
    return vals


def second_best_rows(sim: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    # Returns best and second-best values for each row.
    if sim.shape[1] < 2:
        raise RuntimeError("Need at least two comparison samples.")
    part = np.partition(sim, kth=sim.shape[1] - 2, axis=1)
    second = part[:, -2]
    best = part[:, -1]
    return best, second


def identity_pair_audit(
    pair_name: str,
    a: dict,
    b: dict,
    scanb: dict,
    universe_genes: list[str],
    tech_pairs: pd.DataFrame,
    primary_manifest: pd.DataFrame,
    device: torch.device,
    seed: int,
    cpu_validate: bool,
) -> tuple[dict, pd.DataFrame, dict[str, np.ndarray]]:
    common = [
        g for g in universe_genes
        if g in a["gene_index"] and g in b["gene_index"]
    ]
    if len(common) < 1000:
        raise RuntimeError(f"{pair_name}: only {len(common)} common fingerprint genes.")

    # Calibrate on the same gene concept as far as SCAN-B availability permits.
    calib_genes = [g for g in common if g in scanb["gene_index"]]
    if len(calib_genes) < 1000:
        raise RuntimeError(f"{pair_name}: only {len(calib_genes)} SCAN-B calibration genes.")

    print(f"\n  {pair_name}:")
    print(f"    pairwise fingerprint genes: {len(common):,}")
    print(f"    SCAN-B calibration genes:   {len(calib_genes):,}")

    scanb_prof = torch_profile_matrix(scanb, calib_genes, device)

    true = tech_pairs.loc[tech_pairs["cross_platform"].astype(int) == 1].copy()
    if len(true) != 36:
        raise RuntimeError(f"Expected 36 cross-platform technical pairs, found {len(true)}.")

    true_a = np.array(
        [scanb["sample_index"][x] for x in true["primary_title"]],
        dtype=np.int64,
    )
    true_b = np.array(
        [scanb["sample_index"][x] for x in true["replicate_title"]],
        dtype=np.int64,
    )
    true_corr = paired_profile_dots(scanb_prof, true_a, true_b, chunk=100)

    rng = np.random.default_rng(seed)
    neg_a, neg_b = build_negative_cross_platform_pairs(
        tech_pairs, primary_manifest, scanb, rng
    )
    null_corr = paired_profile_dots(scanb_prof, neg_a, neg_b)

    q_true = quantile(true_corr, 0.01)
    q_null = quantile(null_corr, 0.9999)
    separated = bool(q_true > q_null)
    threshold = max(q_true, q_null) if separated else None

    print(
        f"    true cross-platform q01={q_true:.6f}; "
        f"nonmatch q99.99={q_null:.6f}; separated={separated}"
    )
    if separated:
        print(f"    calibrated hard identity threshold={threshold:.6f}")
    else:
        print("    NO hard threshold: rank-based review only")

    del scanb_prof
    torch.cuda.empty_cache()

    pa = torch_profile_matrix(a, common, device)
    pb = torch_profile_matrix(b, common, device)

    if cpu_validate:
        # Deterministic small-submatrix CPU/GPU agreement check.
        na = min(7, pa.shape[0])
        nb = min(11, pb.shape[0])
        gpu_small = (pa[:na] @ pb[:nb].T).detach().cpu().numpy()

        ca = cpu_profile_matrix(a, common)[:na]
        cb = cpu_profile_matrix(b, common)[:nb]
        cpu_small = ca @ cb.T

        max_diff = float(np.max(np.abs(gpu_small - cpu_small)))
        print(f"    CPU↔GPU fingerprint validation max|Δ|={max_diff:.3e}")
        if max_diff > GPU_CPU_TOL:
            raise RuntimeError(
                f"{pair_name}: CPU/GPU fingerprint mismatch {max_diff} > {GPU_CPU_TOL}"
            )
    else:
        max_diff = None

    sim = (pa @ pb.T).detach().cpu().numpy()
    np.clip(sim, -1.0, 1.0, out=sim)

    del pa, pb
    torch.cuda.empty_cache()

    a_best_j = np.argmax(sim, axis=1)
    b_best_i = np.argmax(sim, axis=0)

    a_best, a_second = second_best_rows(sim)
    b_best, b_second = second_best_rows(sim.T)

    rows = []
    flags = []

    for i, j in enumerate(a_best_j):
        if b_best_i[j] != i:
            continue

        s = float(sim[i, j])
        margin_a = float(a_best[i] - a_second[i])
        margin_b = float(b_best[j] - b_second[j])
        above = bool(threshold is not None and s >= threshold)
        flag = bool(
            threshold is not None
            and above
            and margin_a >= IDENTITY_MIN_MARGIN
            and margin_b >= IDENTITY_MIN_MARGIN
        )

        rec = {
            "pair_name": pair_name,
            "sample_a": a["samples"][i],
            "sample_b": b["samples"][j],
            "similarity": s,
            "margin_a": margin_a,
            "margin_b": margin_b,
            "threshold_available": int(threshold is not None),
            "threshold": threshold if threshold is not None else np.nan,
            "above_threshold": int(above),
            "flag_possible_identity_overlap": int(flag),
        }
        rows.append(rec)
        if flag:
            flags.append(rec)

    mutual = pd.DataFrame(rows)
    if len(mutual):
        mutual = mutual.sort_values("similarity", ascending=False).reset_index(drop=True)

    print(
        f"    mutual nearest-neighbor pairs={len(mutual):,}; "
        f"hard flags={len(flags)}"
    )
    if len(mutual):
        top = mutual.iloc[0]
        print(
            f"    top mutual pair: {top['sample_a']} ↔ {top['sample_b']} "
            f"rho={top['similarity']:.6f}"
        )

    summary = {
        "pair_name": pair_name,
        "n_genes_pairwise": len(common),
        "n_genes_scanb_calibration": len(calib_genes),
        "true_cross_platform_pairs": int(len(true_corr)),
        "null_nonmatching_pairs": int(len(null_corr)),
        "true_similarity_min": float(np.min(true_corr)),
        "true_similarity_median": float(np.median(true_corr)),
        "true_similarity_q01": q_true,
        "null_similarity_median": float(np.median(null_corr)),
        "null_similarity_q9999": q_null,
        "reference_separated": separated,
        "hard_threshold": threshold,
        "identity_margin_required": IDENTITY_MIN_MARGIN,
        "mutual_nearest_neighbor_pairs": int(len(mutual)),
        "hard_flag_count": int(len(flags)),
        "cpu_gpu_validation_max_abs_diff": max_diff,
    }
    distributions = {"true": true_corr, "null": null_corr}
    return summary, mutual, distributions


def source_edge_metrics(edges: np.ndarray) -> dict[str, float]:
    e = np.asarray(edges, dtype=np.float64)
    return {
        "source_coherence_mean_abs_r": float(np.mean(np.abs(e))),
        "source_coherence_median_abs_r": float(np.median(np.abs(e))),
        "source_edge_mean_signed_r": float(np.mean(e)),
        "source_edge_sd_signed_r": float(np.std(e, ddof=1)),
        "source_edge_q05": quantile(e, 0.05),
        "source_edge_q25": quantile(e, 0.25),
        "source_edge_median": quantile(e, 0.50),
        "source_edge_q75": quantile(e, 0.75),
        "source_edge_q95": quantile(e, 0.95),
        "source_edge_fraction_positive": float(np.mean(e > 0)),
    }


def standardize_samples_x_genes_torch(x: torch.Tensor) -> torch.Tensor:
    mu = x.mean(dim=0, keepdim=True)
    sd = x.std(dim=0, correction=1, keepdim=True)
    if torch.any(~torch.isfinite(sd)) or torch.any(sd <= 0):
        raise RuntimeError("Zero/nonfinite gene SD during source bootstrap.")
    return (x - mu) / sd


def selected_pair_correlations(
    z: torch.Tensor,
    ii: np.ndarray,
    jj: np.ndarray,
    chunk: int = 5000,
) -> np.ndarray:
    out = np.empty(len(ii), dtype=np.float64)
    n = z.shape[0]
    for start in range(0, len(ii), chunk):
        end = min(start + chunk, len(ii))
        ti = torch.as_tensor(ii[start:end], dtype=torch.long, device=z.device)
        tj = torch.as_tensor(jj[start:end], dtype=torch.long, device=z.device)
        vals = (
            z.index_select(1, ti) * z.index_select(1, tj)
        ).sum(dim=0) / (n - 1)
        out[start:end] = vals.detach().cpu().numpy()
    return out


def leading_loading_power(
    z: torch.Tensor,
    init: torch.Tensor,
    residual_tol: float = LOADING_BOOT_RESIDUAL_TOL,
    max_iter: int = 120,
) -> tuple[torch.Tensor, float, int, bool]:
    v = init / torch.linalg.vector_norm(init)
    converged = False

    for it in range(1, max_iter + 1):
        w = z.T @ (z @ v)
        wn = torch.linalg.vector_norm(w)
        if not torch.isfinite(wn) or wn <= 0:
            raise RuntimeError("Invalid norm in loading power iteration.")
        new_v = w / wn

        # Eigenvector sign is arbitrary.
        if torch.dot(new_v, v) < 0:
            new_v = -new_v

        v = new_v

        if it % 3 == 0 or it == max_iter:
            av = z.T @ (z @ v)
            lam = torch.dot(v, av)
            residual = torch.linalg.vector_norm(av - lam * v) / (
                torch.abs(lam) + torch.finfo(torch.float64).eps
            )
            r = float(residual.detach().cpu())
            if r <= residual_tol:
                converged = True
                return v, r, it, False

    # Deterministic exact fallback for difficult/near-degenerate cases.
    _, _, vh = torch.linalg.svd(z, full_matrices=False)
    v = vh[0, :]
    if torch.dot(v, init) < 0:
        v = -v
    av = z.T @ (z @ v)
    lam = torch.dot(v, av)
    residual = torch.linalg.vector_norm(av - lam * v) / (
        torch.abs(lam) + torch.finfo(torch.float64).eps
    )
    return v, float(residual.detach().cpu()), max_iter, True


def source_stability_audit(
    tcga: dict,
    weights: pd.DataFrame,
    summary: pd.DataFrame,
    device: torch.device,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    print("\n[5/7] Auditing frozen TCGA source-object replay and bootstrap stability ...")
    program_summaries = []
    bootstrap_rows = []
    reconstruction_rows = []

    genes_all = tcga["genes"]
    gene_index = tcga["gene_index"]
    x_all = tcga["raw_transformed"].T  # samples x 10000
    n_source = x_all.shape[0]

    for pnum, program_id in enumerate(sorted(weights["program_id"].unique()), start=1):
        t0 = time.perf_counter()
        wdf = weights.loc[weights["program_id"] == program_id].copy()
        wdf = wdf.sort_values("source_gene_order_within_program")
        genes = [canon_symbol(x) for x in wdf["Hugo_Symbol"]]
        missing = [g for g in genes if g not in gene_index]
        if missing:
            raise RuntimeError(f"{program_id}: missing frozen source genes: {missing[:20]}")

        cols = np.array([gene_index[g] for g in genes], dtype=np.int64)
        x = x_all[:, cols]
        m = x.shape[1]

        if m != len(genes):
            raise RuntimeError(f"{program_id}: module width mismatch.")

        st = summary.loc[summary["program_id"] == program_id]
        if len(st) != 1:
            raise RuntimeError(f"{program_id}: expected one summary row.")
        st = st.iloc[0]

        print(f"\n  [{pnum:02d}/12] {program_id}: n_genes={m:,}")

        xt = torch.as_tensor(x, dtype=torch.float64, device=device)
        z = standardize_samples_x_genes_torch(xt)

        corr = (z.T @ z) / (n_source - 1)
        tri_i, tri_j = np.triu_indices(m, k=1)
        edges = corr[
            torch.as_tensor(tri_i, dtype=torch.long, device=device),
            torch.as_tensor(tri_j, dtype=torch.long, device=device),
        ].detach().cpu().numpy()

        del corr
        torch.cuda.empty_cache()

        if len(edges) != int(st["n_edges"]):
            raise RuntimeError(
                f"{program_id}: n_edges replay {len(edges)} != frozen {int(st['n_edges'])}"
            )

        replay_metrics = source_edge_metrics(edges)
        diffs = {}
        for key, val in replay_metrics.items():
            frozen_val = float(st[key])
            diffs[key] = abs(val - frozen_val)

        max_metric_diff = max(diffs.values())
        print(f"    frozen edge-summary replay max|Δ|={max_metric_diff:.3e}")
        if max_metric_diff > SOURCE_RECON_TOL:
            raise RuntimeError(
                f"{program_id}: source edge replay mismatch {max_metric_diff} > {SOURCE_RECON_TOL}"
            )

        frozen_loading = wdf["source_pc1_loading"].astype(float).to_numpy()
        norm = np.linalg.norm(frozen_loading)
        if not np.isfinite(norm) or norm <= 0:
            raise RuntimeError(f"{program_id}: invalid frozen loading norm.")
        frozen_loading = frozen_loading / norm

        v = torch.as_tensor(frozen_loading, dtype=torch.float64, device=device)
        av = z.T @ (z @ v)
        lam = torch.dot(v, av)
        residual = torch.linalg.vector_norm(av - lam * v) / (
            torch.abs(lam) + torch.finfo(torch.float64).eps
        )
        full_loading_residual = float(residual.detach().cpu())
        explained = float(
            (lam / torch.sum(z * z)).detach().cpu()
        )
        explained_diff = abs(explained - float(st["source_pc1_variance_explained"]))

        print(
            f"    frozen PC1 eigen-residual={full_loading_residual:.3e}; "
            f"variance-explained |Δ|={explained_diff:.3e}"
        )
        if full_loading_residual > LOADING_FULL_RESIDUAL_TOL:
            raise RuntimeError(
                f"{program_id}: frozen loading is not replaying PC1 accurately."
            )
        if explained_diff > SOURCE_RECON_TOL:
            raise RuntimeError(
                f"{program_id}: source PC1 variance replay mismatch {explained_diff}."
            )

        # Fixed uniformly sampled edge subset for all bootstraps.
        rng_edges = np.random.default_rng(BASE_SEED + 10_000 + pnum)
        if len(edges) <= SOURCE_MAX_EDGES:
            chosen = np.arange(len(edges), dtype=np.int64)
        else:
            chosen = np.sort(
                rng_edges.choice(
                    len(edges),
                    size=SOURCE_MAX_EDGES,
                    replace=False,
                )
            )

        ii = tri_i[chosen]
        jj = tri_j[chosen]
        frozen_edge_subset = edges[chosen]
        frozen_edge_ranks = rankdata(frozen_edge_subset, method="average")

        # CPU↔GPU selected-edge validation on full source for deterministic first 100 pairs.
        nval = min(100, len(ii))
        x_cpu = x.copy()
        x_cpu = (x_cpu - x_cpu.mean(axis=0, keepdims=True)) / x_cpu.std(
            axis=0, ddof=1, keepdims=True
        )
        cpu_vals = np.sum(
            x_cpu[:, ii[:nval]] * x_cpu[:, jj[:nval]],
            axis=0,
        ) / (n_source - 1)
        gpu_vals = selected_pair_correlations(z, ii[:nval], jj[:nval], chunk=nval)
        gpu_cpu_diff = float(np.max(np.abs(cpu_vals - gpu_vals)))
        if gpu_cpu_diff > GPU_CPU_TOL:
            raise RuntimeError(
                f"{program_id}: selected-edge CPU/GPU mismatch {gpu_cpu_diff}."
            )

        rng = np.random.default_rng(BASE_SEED + 20_000 + pnum)
        edge_rhos = []
        load_rhos = []
        load_resids = []
        fallback_count = 0

        for b in range(SOURCE_EDGE_BOOTSTRAPS):
            sample_idx = rng.integers(0, n_source, size=n_source, endpoint=False)
            tidx = torch.as_tensor(sample_idx, dtype=torch.long, device=device)
            xb = xt.index_select(0, tidx)
            zb = standardize_samples_x_genes_torch(xb)

            boot_edges = selected_pair_correlations(zb, ii, jj)

            # Spearman against frozen full-source edge subset.
            rb = rankdata(boot_edges, method="average")
            rb -= rb.mean()
            rf = frozen_edge_ranks - frozen_edge_ranks.mean()
            den = math.sqrt(float(np.dot(rb, rb) * np.dot(rf, rf)))
            edge_rho = float(np.dot(rb, rf) / den) if den > 0 else float("nan")
            edge_rhos.append(edge_rho)

            rec = {
                "program_id": program_id,
                "bootstrap_index": b + 1,
                "edge_spearman": edge_rho,
                "loading_spearman": np.nan,
                "loading_eigen_residual": np.nan,
                "loading_exact_svd_fallback": 0,
            }

            if b < SOURCE_LOADING_BOOTSTRAPS:
                init = torch.as_tensor(
                    frozen_loading,
                    dtype=torch.float64,
                    device=device,
                )
                vb, resid, _, fallback = leading_loading_power(zb, init)
                if torch.dot(vb, init) < 0:
                    vb = -vb

                vb_np = vb.detach().cpu().numpy()
                load_rho = spearman_fast(frozen_loading, vb_np)
                load_rhos.append(load_rho)
                load_resids.append(resid)
                fallback_count += int(fallback)

                rec["loading_spearman"] = load_rho
                rec["loading_eigen_residual"] = resid
                rec["loading_exact_svd_fallback"] = int(fallback)

            bootstrap_rows.append(rec)

            del tidx, xb, zb
            if (b + 1) % 50 == 0:
                print(
                    f"    bootstraps {b+1:3d}/{SOURCE_EDGE_BOOTSTRAPS}: "
                    f"edge median-so-far={np.median(edge_rhos):.4f}"
                )

        edge_arr = np.asarray(edge_rhos, dtype=np.float64)
        load_arr = np.asarray(load_rhos, dtype=np.float64)

        ps = {
            "program_id": program_id,
            "n_genes": m,
            "n_edges_total": len(edges),
            "n_edges_bootstrap_subset": len(chosen),
            "source_edge_replay_max_abs_metric_diff": max_metric_diff,
            "source_loading_full_eigen_residual": full_loading_residual,
            "source_pc1_variance_explained_replay_abs_diff": explained_diff,
            "gpu_cpu_selected_edge_max_abs_diff": gpu_cpu_diff,
            "edge_stability_median": float(np.median(edge_arr)),
            "edge_stability_q05": quantile(edge_arr, 0.05),
            "edge_stability_q95": quantile(edge_arr, 0.95),
            "loading_stability_median": float(np.median(load_arr)),
            "loading_stability_q05": quantile(load_arr, 0.05),
            "loading_stability_q95": quantile(load_arr, 0.95),
            "loading_bootstrap_max_eigen_residual": float(np.max(load_resids)),
            "loading_exact_svd_fallback_count": int(fallback_count),
        }
        program_summaries.append(ps)

        reconstruction_rows.append(
            {
                "program_id": program_id,
                **{f"replay_{k}": v for k, v in replay_metrics.items()},
                **{f"absdiff_{k}": v for k, v in diffs.items()},
                "replay_source_pc1_variance_explained": explained,
                "absdiff_source_pc1_variance_explained": explained_diff,
                "frozen_loading_eigen_residual": full_loading_residual,
            }
        )

        print(
            f"    edge stability median [q05,q95] = "
            f"{ps['edge_stability_median']:.4f} "
            f"[{ps['edge_stability_q05']:.4f}, {ps['edge_stability_q95']:.4f}]"
        )
        print(
            f"    loading stability median [q05,q95] = "
            f"{ps['loading_stability_median']:.4f} "
            f"[{ps['loading_stability_q05']:.4f}, {ps['loading_stability_q95']:.4f}]"
        )
        print(
            f"    exact-SVD loading fallbacks={fallback_count}; "
            f"elapsed {time.perf_counter()-t0:.1f}s"
        )

        del xt, z, edges, x_cpu
        torch.cuda.empty_cache()
        gc.collect()

    return (
        pd.DataFrame(program_summaries),
        pd.DataFrame(bootstrap_rows),
        pd.DataFrame(reconstruction_rows),
    )


def genomic_concentration_audit(weights: pd.DataFrame) -> pd.DataFrame:
    print("\n[6/7] Auditing source-module chromosomal concentration (descriptive only) ...")
    counts = defaultdict(Counter)

    gene_pat = re.compile(r'geneSymbol "([^"]+)"')

    with gzip.open(GTF, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue
            chrom = fields[0]
            if chrom not in CANONICAL_CHROMS:
                continue
            m = gene_pat.search(fields[8])
            if not m:
                continue
            gene = canon_symbol(m.group(1))
            if gene:
                counts[gene][chrom] += 1

    gene_to_chrom = {}
    ambiguous = set()
    for gene, cc in counts.items():
        if not cc:
            continue
        max_n = max(cc.values())
        winners = sorted([chrom for chrom, n in cc.items() if n == max_n])
        if len(winners) == 1:
            gene_to_chrom[gene] = winners[0]
        else:
            ambiguous.add(gene)

    rows = []
    for program_id in sorted(weights["program_id"].unique()):
        genes = [
            canon_symbol(x)
            for x in weights.loc[weights["program_id"] == program_id, "Hugo_Symbol"]
        ]
        mapped = [gene_to_chrom[g] for g in genes if g in gene_to_chrom]
        cc = Counter(mapped)
        dominant_chrom = ""
        dominant_n = 0
        if cc:
            dominant_chrom, dominant_n = cc.most_common(1)[0]

        rec = {
            "program_id": program_id,
            "n_genes": len(genes),
            "n_genes_mapped_unique_canonical_chromosome": len(mapped),
            "n_genes_ambiguous_chromosome": sum(g in ambiguous for g in genes),
            "mapped_fraction": len(mapped) / len(genes) if genes else np.nan,
            "dominant_chromosome": dominant_chrom,
            "dominant_chromosome_n": dominant_n,
            "dominant_fraction_all_genes": dominant_n / len(genes) if genes else np.nan,
            "dominant_fraction_mapped_genes": dominant_n / len(mapped) if mapped else np.nan,
        }
        for chrom in sorted(CANONICAL_CHROMS, key=lambda x: (len(x), x)):
            rec[f"n_{chrom}"] = int(cc.get(chrom, 0))
        rows.append(rec)

        print(
            f"  {program_id}: mapped={len(mapped)}/{len(genes)}; "
            f"dominant={dominant_chrom or 'NA'} "
            f"{dominant_n}/{len(genes)} ({rec['dominant_fraction_all_genes']:.3f})"
        )

    return pd.DataFrame(rows)


def platform_pam50_audit() -> pd.DataFrame:
    print("\n[7/7] Reporting SCAN-B PAM50 composition by sequencing platform ...")
    df = pd.read_csv(SCANB_PAM50, sep="\t", dtype=str).fillna("")
    needed = {"platform_id", "pam50_subtype"}
    if not needed.issubset(df.columns):
        raise RuntimeError(f"SCAN-B PAM50 alignment missing columns: {needed-set(df.columns)}")

    ct = pd.crosstab(df["platform_id"], df["pam50_subtype"], dropna=False)
    rows = []
    for platform, row in ct.iterrows():
        total = int(row.sum())
        for subtype, n in row.items():
            rows.append(
                {
                    "platform_id": platform,
                    "pam50_subtype": subtype,
                    "n": int(n),
                    "platform_total": total,
                    "fraction_within_platform": float(n / total) if total else np.nan,
                }
            )
        print(f"  {platform}: n={total}; {dict((k, int(v)) for k, v in row.items())}")
    return pd.DataFrame(rows)


def main() -> None:
    print("=" * 138)
    print("Paper 4 / TCBB - execute FINAL PRE-TARGET identity, source-stability, and source-annotation audits")
    print("=" * 138)
    print(f"Script version: {SCRIPT_VERSION}")
    print()
    print("Scientific guard:")
    print("  Target expression values read:                       YES (identity fingerprint only)")
    print("  Target outcomes/treatment loaded:                    NO")
    print("  Target gene-gene preservation calculated:            NO")
    print("  Target PCA/loading preservation calculated:           NO")
    print("  Mapping-null preservation calculated:                 NO")
    print("  Primary Strong/Partial/No-clear classifier evaluated: NO")
    print("  Source bootstrap stability calculated:                YES")
    print("=" * 138)

    for p in [
        CONTRACT,
        TCGA_EXPR,
        SCANB_EXPR,
        METABRIC_EXPR,
        UNIVERSE,
        MEMBERSHIP,
        WEIGHTS,
        PROGRAM_SUMMARY,
        TECH_PAIRS,
        PRIMARY,
        SCANB_PAM50,
        GTF,
    ]:
        require(p)

    device = ensure_cuda()
    prop = torch.cuda.get_device_properties(0)
    print(
        f"CUDA device: {prop.name}; VRAM={prop.total_memory/1024**3:.2f} GB; "
        f"capability={prop.major}.{prop.minor}"
    )

    contract = json.loads(CONTRACT.read_text(encoding="utf-8"))
    if contract.get("status") != "FROZEN_BEFORE_FIRST_TARGET_PRESERVATION_STATISTIC":
        raise RuntimeError("04f contract is not in the expected frozen state.")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    (OUT_DIR / "identity_calibration_distributions").mkdir(exist_ok=True)
    (OUT_DIR / "identity_mutual_nearest").mkdir(exist_ok=True)

    universe = pd.read_csv(UNIVERSE, sep="\t", dtype=str).fillna("")
    universe["source_gene_rank_by_MAD"] = universe["source_gene_rank_by_MAD"].astype(int)
    universe = universe.sort_values("source_gene_rank_by_MAD").reset_index(drop=True)
    if len(universe) != 10_000:
        raise RuntimeError(f"Expected 10,000 frozen source genes, found {len(universe)}.")

    universe_genes = [canon_symbol(x) for x in universe["Hugo_Symbol"]]
    if len(set(universe_genes)) != 10_000:
        raise RuntimeError("Frozen source universe symbols are not unique after canonicalization.")

    universe_set = set(universe_genes)

    weights = pd.read_csv(WEIGHTS, sep="\t")
    summary = pd.read_csv(PROGRAM_SUMMARY, sep="\t")
    tech_pairs = pd.read_csv(TECH_PAIRS, sep="\t", dtype=str).fillna("")
    tech_pairs["cross_platform"] = tech_pairs["cross_platform"].astype(int)
    tech_pairs["same_platform"] = tech_pairs["same_platform"].astype(int)
    primary_manifest = pd.read_csv(PRIMARY, sep="\t", dtype=str).fillna("")
    primary_titles = list(primary_manifest["primary_title"])

    # Load/replay cohort expression only now.
    tcga = load_source_tcga(universe)
    scanb = load_scanb(universe_set, primary_titles)
    met = load_metabric(universe_set)

    print("\n[4/7] Cross-cohort identity / near-duplicate audit ...")
    pair_defs = [
        ("TCGA_vs_SCANB", tcga, scanb),
        ("TCGA_vs_METABRIC", tcga, met),
        ("SCANB_vs_METABRIC", scanb, met),
    ]

    identity_summaries = []
    all_mutual = []
    any_hard_flag = False

    for idx, (name, a, b) in enumerate(pair_defs, start=1):
        summ, mutual, dist = identity_pair_audit(
            pair_name=name,
            a=a,
            b=b,
            scanb=scanb,
            universe_genes=universe_genes,
            tech_pairs=tech_pairs,
            primary_manifest=primary_manifest,
            device=device,
            seed=BASE_SEED + idx * 100,
            cpu_validate=True,
        )
        identity_summaries.append(summ)
        any_hard_flag = any_hard_flag or bool(summ["hard_flag_count"])

        np.savez_compressed(
            OUT_DIR
            / "identity_calibration_distributions"
            / f"{name}_identity_calibration_v1.npz",
            true_cross_platform=dist["true"],
            null_cross_platform_nonmatch=dist["null"],
        )

        if len(mutual):
            mutual.to_csv(
                OUT_DIR
                / "identity_mutual_nearest"
                / f"{name}_mutual_nearest_v1.tsv",
                sep="\t",
                index=False,
            )
            all_mutual.append(mutual)

    identity_df = pd.DataFrame(identity_summaries)
    identity_df.to_csv(
        OUT_DIR / "cross_cohort_identity_summary_v1.tsv",
        sep="\t",
        index=False,
    )

    if all_mutual:
        pd.concat(all_mutual, ignore_index=True).to_csv(
            OUT_DIR / "cross_cohort_identity_mutual_nearest_all_v1.tsv",
            sep="\t",
            index=False,
        )

    # Source stability is independent of any target preservation result and remains valid
    # even if identity review is required.
    source_summary_df, source_boot_df, source_replay_df = source_stability_audit(
        tcga=tcga,
        weights=weights,
        summary=summary,
        device=device,
    )

    source_summary_df.to_csv(
        OUT_DIR / "source_object_stability_summary_v1.tsv",
        sep="\t",
        index=False,
    )
    source_boot_df.to_csv(
        OUT_DIR / "source_object_stability_bootstraps_v1.tsv",
        sep="\t",
        index=False,
    )
    source_replay_df.to_csv(
        OUT_DIR / "source_object_replay_validation_v1.tsv",
        sep="\t",
        index=False,
    )

    chromosome_df = genomic_concentration_audit(weights)
    chromosome_df.to_csv(
        OUT_DIR / "source_module_chromosome_concentration_v1.tsv",
        sep="\t",
        index=False,
    )

    platform_pam50_df = platform_pam50_audit()
    platform_pam50_df.to_csv(
        OUT_DIR / "scanb_platform_pam50_composition_v1.tsv",
        sep="\t",
        index=False,
    )

    final = {
        "script_version": SCRIPT_VERSION,
        "contract": str(CONTRACT),
        "status": (
            "PAUSE_FOR_IDENTITY_REVIEW"
            if any_hard_flag
            else "PASS_READY_FOR_FIRST_POOLED_TARGET_PRESERVATION"
        ),
        "scientific_guard": {
            "target_expression_values_read_for_identity_only": True,
            "target_outcomes_loaded": False,
            "target_treatment_loaded": False,
            "target_gene_gene_preservation_calculated": False,
            "target_pca_loading_preservation_calculated": False,
            "mapping_null_preservation_calculated": False,
            "primary_classifier_evaluated": False,
        },
        "cuda": {
            "device": prop.name,
            "vram_gb": prop.total_memory / 1024**3,
        },
        "identity_audit": identity_summaries,
        "hard_identity_flag_any": bool(any_hard_flag),
        "source_stability": source_summary_df.to_dict(orient="records"),
        "platform_pam50_composition_file": str(
            OUT_DIR / "scanb_platform_pam50_composition_v1.tsv"
        ),
        "chromosome_concentration_file": str(
            OUT_DIR / "source_module_chromosome_concentration_v1.tsv"
        ),
    }

    json_out = OUT_DIR / "final_pretarget_audit_v1.json"
    json_out.write_text(json.dumps(final, indent=2), encoding="utf-8")

    print("\n" + "=" * 138)
    if any_hard_flag:
        print("04h FINAL PRE-TARGET AUDIT: PAUSE FOR IDENTITY REVIEW")
    else:
        print("04h FINAL PRE-TARGET AUDIT: PASS — READY FOR FIRST POOLED TARGET PRESERVATION")
    print("=" * 138)

    for rec in identity_summaries:
        print(
            f"{rec['pair_name']}: genes={rec['n_genes_pairwise']:,}, "
            f"reference_separated={rec['reference_separated']}, "
            f"hard_flags={rec['hard_flag_count']}"
        )

    print()
    print("Source stability summaries:")
    for row in source_summary_df.itertuples(index=False):
        print(
            f"  {row.program_id}: "
            f"edge median={row.edge_stability_median:.4f} "
            f"[{row.edge_stability_q05:.4f},{row.edge_stability_q95:.4f}]; "
            f"loading median={row.loading_stability_median:.4f} "
            f"[{row.loading_stability_q05:.4f},{row.loading_stability_q95:.4f}]"
        )

    print()
    print("No target preservation statistic was calculated.")
    if any_hard_flag:
        print("DO NOT run pooled preservation until flagged identity candidates are reviewed.")
    else:
        print("No hard cross-cohort identity flag was found; pooled preservation may proceed under the frozen contracts.")
    print()
    print(f"Master output: {json_out}")
    print("=" * 138)


if __name__ == "__main__":
    main()
