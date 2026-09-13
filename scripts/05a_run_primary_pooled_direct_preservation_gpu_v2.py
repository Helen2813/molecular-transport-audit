from __future__ import annotations

import gc
import json
import math
import re
import time
from pathlib import Path

import numpy as np
import pandas as pd

try:
    import torch
except Exception as exc:
    raise RuntimeError("PyTorch is required. Run the CUDA environment audit first.") from exc

try:
    from scipy.stats import rankdata
except Exception as exc:
    raise RuntimeError("scipy is required for exact Spearman statistics.") from exc


SCRIPT_VERSION = "05a-run-primary-pooled-direct-preservation-gpu-v2-no-cli"

DATA_ROOT = Path(r"D:\paper4_tcbb_data")

PRETARGET_AUDIT = (
    DATA_ROOT
    / "paper4_tcbb_final_pretarget_audits_v1"
    / "final_pretarget_audit_v1.json"
)
EVALUABILITY_AMENDMENT = (
    DATA_ROOT
    / "paper4_tcbb_target_statistical_evaluability_amendment_v2"
    / "target_statistical_evaluability_amendment_v2.json"
)
EXACT_VARIATION_CORRECTION = (
    DATA_ROOT
    / "paper4_tcbb_exact_variation_evaluability_correction_v1"
    / "exact_variation_evaluability_correction_v1.json"
)
PRESERVATION_CONTRACT = (
    DATA_ROOT
    / "paper4_tcbb_preservation_operating_contract_v1"
    / "preservation_operating_characteristics_contract_v1.json"
)

UNIVERSE = (
    DATA_ROOT
    / "paper4_tcbb_tcga_source_universe_v1"
    / "tcga_source_gene_universe_frozen_v1.tsv"
)
WEIGHTS = (
    DATA_ROOT
    / "paper4_tcbb_frozen_source_programs_v1"
    / "tcga_frozen_source_program_weights_v1.tsv"
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
SCANB_PRIMARY = (
    DATA_ROOT
    / "paper4_tcbb_scanb_pairing_contract_v1"
    / "scanb_primary_profiles_frozen_v1.tsv"
)
METABRIC_EXPR = (
    DATA_ROOT
    / "paper4_tcbb_input_audit_v1"
    / "staged_continuous_inputs"
    / "METABRIC"
    / "data_mrna_illumina_microarray.txt"
)

OUT_DIR = DATA_ROOT / "paper4_tcbb_primary_pooled_direct_preservation_v2"

DIRECT_PERMUTATIONS = 10_000
DIRECT_ALPHA = 0.05
MIN_EVALUABLE_GENES = 30
MIN_COVERAGE = 0.80
BASE_SEED = 20260916

CPU_GPU_CORR_TOL = 3e-10
PC1_RESIDUAL_TOL = 1e-10
MAX_POWER_ITER = 400

# GPU batch memory guard for the exact full-edge QAP/gene-label permutation.
MAX_EDGE_INDEX_ELEMENTS_PER_BATCH = 20_000_000
LOADING_PERM_BATCH = 1024


def require(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Required file does not exist: {path}")


def canon_symbol(x: object) -> str:
    if x is None:
        return ""
    return " ".join(str(x).strip().split()).upper()


def ensure_cuda() -> torch.device:
    if not torch.cuda.is_available():
        raise RuntimeError("CUDA is not available in the active Python environment.")
    torch.set_default_dtype(torch.float64)
    return torch.device("cuda:0")


def bh_adjust(pvalues: list[float]) -> list[float]:
    p = np.asarray(pvalues, dtype=np.float64)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order]
    qrank = ranked * n / np.arange(1, n + 1, dtype=np.float64)
    qrank = np.minimum.accumulate(qrank[::-1])[::-1]
    qrank = np.minimum(qrank, 1.0)
    q = np.empty_like(qrank)
    q[order] = qrank
    return q.tolist()


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
    if den <= 0:
        return float("nan")
    return float(np.dot(ra, rb) / den)


def standardize_samples_x_genes_np(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=np.float64)
    if not np.isfinite(x).all():
        raise RuntimeError("Non-finite values reached standardization.")
    mu = x.mean(axis=0, keepdims=True)
    sd = x.std(axis=0, ddof=1, keepdims=True)
    if np.any(~np.isfinite(sd)) or np.any(sd <= 0):
        raise RuntimeError("Zero/non-finite gene SD reached standardization.")
    return (x - mu) / sd


def load_tcga_source(universe: pd.DataFrame) -> dict:
    print("\n[1/5] Replaying frozen TCGA source matrix ...")
    t0 = time.perf_counter()
    raw = pd.read_csv(TCGA_EXPR, sep="\t", low_memory=False)

    row_idx = universe["source_row_index_0based"].astype(int).to_numpy()
    selected = raw.iloc[row_idx, :]

    expected = [canon_symbol(x) for x in universe["Hugo_Symbol"]]
    observed = [canon_symbol(x) for x in selected["Hugo_Symbol"]]
    if expected != observed:
        raise RuntimeError("TCGA frozen 10,000-gene universe replay mismatch.")

    values = selected.iloc[:, 2:].to_numpy(dtype=np.float64)
    if np.any(values < 0) or not np.isfinite(values).all():
        raise RuntimeError("TCGA source RSEM contains negative/non-finite values.")

    transformed = np.log2(values + 1.0)
    z = standardize_samples_x_genes_np(transformed.T)  # samples x genes

    out = {
        "name": "TCGA_BRCA",
        "samples": list(selected.columns[2:]),
        "genes": expected,
        "gene_index": {g: i for i, g in enumerate(expected)},
        "z": z,
    }

    del raw, selected, values, transformed
    gc.collect()

    print(f"  {z.shape[0]:,} samples x {z.shape[1]:,} genes; {time.perf_counter()-t0:.1f}s")
    return out


def load_scanb_target(
    universe_set: set[str],
    primary_titles: list[str],
) -> dict:
    print("\n[2/5] Loading frozen SCAN-B primary target matrix ...")
    t0 = time.perf_counter()
    parts = []

    for chunk in pd.read_csv(
        SCANB_EXPR,
        compression="gzip",
        chunksize=1500,
        low_memory=False,
    ):
        symbol_col = chunk.columns[0]
        symbols = chunk[symbol_col].map(canon_symbol)
        mask = symbols.isin(universe_set)
        if not mask.any():
            continue

        sub = chunk.loc[mask, [symbol_col] + primary_titles]
        arr = sub[primary_titles].to_numpy(dtype=np.float64)

        # Frozen preprocessing: published log2(FPKM+0.1) -> FPKM -> log2(FPKM+1).
        fpkm = np.maximum(np.exp2(arr) - 0.1, 0.0)
        transformed = np.log2(fpkm + 1.0)

        tmp = pd.DataFrame(transformed, columns=primary_titles)
        tmp.insert(0, "Hugo_Symbol", symbols.loc[mask].to_numpy())
        parts.append(tmp)

    if not parts:
        raise RuntimeError("No frozen-universe genes found in SCAN-B.")

    df = pd.concat(parts, ignore_index=True)
    df = df.groupby("Hugo_Symbol", sort=False, as_index=True).mean(numeric_only=True)

    x = df.to_numpy(dtype=np.float64).T  # samples x genes
    finite = np.isfinite(x).all(axis=0)
    exact_min = np.full(x.shape[1], np.nan, dtype=np.float64)
    exact_max = np.full(x.shape[1], np.nan, dtype=np.float64)
    if finite.any():
        exact_min[finite] = np.min(x[:, finite], axis=0)
        exact_max[finite] = np.max(x[:, finite], axis=0)
    eligible = finite & (exact_max > exact_min)

    genes_all = list(df.index)
    genes = [g for g, ok in zip(genes_all, eligible) if ok]
    x = x[:, eligible]
    z = standardize_samples_x_genes_np(x)

    if len(genes) != 9220:
        raise RuntimeError(f"Expected 9,220 corrected-evaluable SCAN-B genes; found {len(genes)}.")

    print(f"  {z.shape[0]:,} primary samples x {z.shape[1]:,} evaluable genes; {time.perf_counter()-t0:.1f}s")
    return {
        "name": "SCANB_GSE96058",
        "samples": primary_titles,
        "genes": genes,
        "gene_index": {g: i for i, g in enumerate(genes)},
        "z": z,
    }


def load_metabric_target(universe_set: set[str]) -> dict:
    print("\n[3/5] Loading frozen METABRIC target matrix ...")
    t0 = time.perf_counter()
    parts = []
    sample_cols = None

    for chunk in pd.read_csv(
        METABRIC_EXPR,
        sep="\t",
        chunksize=1500,
        low_memory=False,
    ):
        if sample_cols is None:
            sample_cols = list(chunk.columns[2:])

        symbols = chunk["Hugo_Symbol"].map(canon_symbol)
        mask = symbols.isin(universe_set)
        if not mask.any():
            continue

        arr = chunk.loc[mask, sample_cols].to_numpy(dtype=np.float64)
        tmp = pd.DataFrame(arr, columns=sample_cols)
        tmp.insert(0, "Hugo_Symbol", symbols.loc[mask].to_numpy())
        parts.append(tmp)

    if not parts:
        raise RuntimeError("No frozen-universe genes found in METABRIC.")

    df = pd.concat(parts, ignore_index=True)
    # Frozen target rule: duplicate symbols averaged on the as-provided transformed scale.
    df = df.groupby("Hugo_Symbol", sort=False, as_index=True).mean(numeric_only=True)

    x = df.to_numpy(dtype=np.float64).T
    finite = np.isfinite(x).all(axis=0)
    exact_min = np.full(x.shape[1], np.nan, dtype=np.float64)
    exact_max = np.full(x.shape[1], np.nan, dtype=np.float64)
    if finite.any():
        exact_min[finite] = np.min(x[:, finite], axis=0)
        exact_max[finite] = np.max(x[:, finite], axis=0)
    eligible = finite & (exact_max > exact_min)

    genes_all = list(df.index)
    genes = [g for g, ok in zip(genes_all, eligible) if ok]
    x = x[:, eligible]
    z = standardize_samples_x_genes_np(x)

    if len(genes) != 8485:
        raise RuntimeError(f"Expected 8,485 corrected-evaluable METABRIC genes; found {len(genes)}.")

    print(f"  {z.shape[0]:,} samples x {z.shape[1]:,} evaluable genes; {time.perf_counter()-t0:.1f}s")
    return {
        "name": "METABRIC",
        "samples": sample_cols,
        "genes": genes,
        "gene_index": {g: i for i, g in enumerate(genes)},
        "z": z,
    }


def correlation_matrix_gpu(
    z_np: np.ndarray,
    device: torch.device,
) -> tuple[np.ndarray, float]:
    n = z_np.shape[0]
    zt = torch.as_tensor(z_np, dtype=torch.float64, device=device)
    ct = (zt.T @ zt) / (n - 1)
    c = ct.detach().cpu().numpy()
    del zt, ct
    torch.cuda.empty_cache()
    np.fill_diagonal(c, 1.0)
    return c, float(n)


def cpu_gpu_corr_validation(
    z_np: np.ndarray,
    corr_gpu: np.ndarray,
    max_genes: int = 50,
) -> float:
    k = min(max_genes, z_np.shape[1])
    cpu = (z_np[:, :k].T @ z_np[:, :k]) / (z_np.shape[0] - 1)
    diff = float(np.max(np.abs(cpu - corr_gpu[:k, :k])))
    if diff > CPU_GPU_CORR_TOL:
        raise RuntimeError(
            f"CPU/GPU correlation mismatch {diff:.3e} > {CPU_GPU_CORR_TOL:.3e}"
        )
    return diff


def leading_pc1_loading(
    target_corr: np.ndarray,
    target_z: np.ndarray,
    source_loading: np.ndarray,
    source_sign: np.ndarray,
    device: torch.device,
    seed: int,
) -> tuple[np.ndarray, dict]:
    """
    Leading eigenvector of target correlation matrix with deterministic multi-start
    power iteration, exact-eigh fallback, then orientation by the frozen source-sign score.
    """
    m = target_corr.shape[0]
    ct = torch.as_tensor(target_corr, dtype=torch.float64, device=device)

    rng = np.random.default_rng(seed)
    starts = [
        np.asarray(source_loading, dtype=np.float64),
        np.asarray(source_sign, dtype=np.float64),
        rng.normal(size=m),
    ]

    candidates = []
    eps = torch.finfo(torch.float64).eps

    for start_idx, start in enumerate(starts):
        if not np.isfinite(start).all() or np.linalg.norm(start) == 0:
            continue
        v = torch.as_tensor(start, dtype=torch.float64, device=device)
        v = v / torch.linalg.vector_norm(v)

        residual_value = float("inf")
        for it in range(1, MAX_POWER_ITER + 1):
            w = ct @ v
            wn = torch.linalg.vector_norm(w)
            if not torch.isfinite(wn) or wn <= 0:
                break
            v_new = w / wn
            if torch.dot(v_new, v) < 0:
                v_new = -v_new
            v = v_new

            if it % 5 == 0 or it == MAX_POWER_ITER:
                av = ct @ v
                lam = torch.dot(v, av)
                residual = torch.linalg.vector_norm(av - lam * v) / (
                    torch.abs(lam) + eps
                )
                residual_value = float(residual.detach().cpu())
                if residual_value <= PC1_RESIDUAL_TOL:
                    break

        av = ct @ v
        lam = float(torch.dot(v, av).detach().cpu())
        candidates.append(
            (
                lam,
                residual_value,
                it,
                start_idx,
                v.detach().cpu().numpy(),
            )
        )

    if not candidates:
        raise RuntimeError("No valid power-iteration candidate for target PC1.")

    candidates.sort(key=lambda x: x[0], reverse=True)
    lam, residual, iters, chosen_start, loading = candidates[0]
    exact_fallback = False

    # If the best deterministic multi-start candidate did not converge cleanly,
    # use an exact symmetric eigendecomposition.
    if residual > PC1_RESIDUAL_TOL:
        evals, evecs = torch.linalg.eigh(ct)
        v = evecs[:, -1]
        lam_t = evals[-1]
        av = ct @ v
        residual_t = torch.linalg.vector_norm(av - lam_t * v) / (
            torch.abs(lam_t) + eps
        )
        loading = v.detach().cpu().numpy()
        lam = float(lam_t.detach().cpu())
        residual = float(residual_t.detach().cpu())
        exact_fallback = True

    del ct
    torch.cuda.empty_cache()

    # Frozen orientation rule: PC1 sample scores must correlate positively
    # with mean(sign(source loading) * target z).
    pc_scores = target_z @ loading
    sign_score = np.mean(target_z * source_sign[None, :], axis=1)
    orient_cor = float(np.corrcoef(pc_scores, sign_score)[0, 1])

    if not np.isfinite(orient_cor):
        raise RuntimeError("Target PC1 orientation correlation is non-finite.")
    if orient_cor < 0:
        loading = -loading
        pc_scores = -pc_scores
        orient_cor = -orient_cor

    explained = float(lam / np.trace(target_corr))

    return loading, {
        "pc1_eigen_residual": residual,
        "pc1_exact_eigh_fallback": int(exact_fallback),
        "pc1_power_start_index": chosen_start,
        "pc1_power_iterations": iters,
        "pc1_orientation_correlation_with_frozen_sign_score": orient_cor,
        "pc1_variance_explained": explained,
    }


def exact_edge_label_permutation_test_gpu(
    source_edges: np.ndarray,
    target_corr: np.ndarray,
    observed_rho: float,
    permutations: int,
    seed: int,
    device: torch.device,
    progress_prefix: str,
) -> tuple[float, np.ndarray, dict]:
    """
    Exact implementation of the frozen gene-label permutation scheme on ALL edges.

    Gene labels are permuted as one permutation of target correlation-matrix
    rows/columns. Spearman ranks are computed once because relabeling only
    permutes the target edge values.
    """
    m = target_corr.shape[0]
    tri_i, tri_j = np.triu_indices(m, k=1)
    target_edges = target_corr[tri_i, tri_j]

    if len(source_edges) != len(target_edges):
        raise RuntimeError("Source/target edge-vector length mismatch.")

    source_rank = rankdata(source_edges, method="average").astype(np.float64)
    target_rank = rankdata(target_edges, method="average").astype(np.float64)

    source_rank -= source_rank.mean()
    target_rank -= target_rank.mean()

    source_ss = float(np.dot(source_rank, source_rank))
    target_ss = float(np.dot(target_rank, target_rank))
    denom = math.sqrt(source_ss * target_ss)
    if denom <= 0:
        raise RuntimeError("Degenerate rank variance in edge permutation test.")

    rank_mat = np.zeros((m, m), dtype=np.float64)
    rank_mat[tri_i, tri_j] = target_rank
    rank_mat[tri_j, tri_i] = target_rank

    src_t = torch.as_tensor(source_rank, dtype=torch.float64, device=device)
    rank_t = torch.as_tensor(rank_mat, dtype=torch.float64, device=device)
    ti = torch.as_tensor(tri_i, dtype=torch.long, device=device)
    tj = torch.as_tensor(tri_j, dtype=torch.long, device=device)

    n_edges = len(tri_i)
    batch_size = max(
        1,
        min(
            128,
            int(MAX_EDGE_INDEX_ELEMENTS_PER_BATCH // max(1, n_edges)),
        ),
    )

    null = np.empty(permutations, dtype=np.float64)
    rng = np.random.default_rng(seed)

    done = 0
    while done < permutations:
        k = min(batch_size, permutations - done)

        perms_np = np.empty((k, m), dtype=np.int64)
        for r in range(k):
            perms_np[r, :] = rng.permutation(m)

        p = torch.as_tensor(perms_np, dtype=torch.long, device=device)
        pi = p[:, ti]
        pj = p[:, tj]
        vals = rank_t[pi, pj]
        rhos = torch.matmul(vals, src_t) / denom
        null[done : done + k] = rhos.detach().cpu().numpy()

        del p, pi, pj, vals, rhos
        done += k

        if done == permutations or done % 2000 < batch_size:
            print(
                f"      {progress_prefix} edge permutations: "
                f"{done:,}/{permutations:,}"
            )

    del src_t, rank_t, ti, tj
    torch.cuda.empty_cache()

    extreme = int(np.sum(np.abs(null) >= abs(observed_rho)))
    pvalue = (1 + extreme) / (permutations + 1)

    return pvalue, null, {
        "permutation_batch_size": batch_size,
        "extreme_count": extreme,
        "null_mean": float(np.mean(null)),
        "null_sd": float(np.std(null, ddof=1)),
        "null_q025": float(np.quantile(null, 0.025)),
        "null_q975": float(np.quantile(null, 0.975)),
    }


def exact_loading_label_permutation_test(
    source_loading: np.ndarray,
    target_loading: np.ndarray,
    observed_rho: float,
    permutations: int,
    seed: int,
) -> tuple[float, np.ndarray, dict]:
    sr = rankdata(source_loading, method="average").astype(np.float64)
    tr = rankdata(target_loading, method="average").astype(np.float64)

    sr -= sr.mean()
    tr -= tr.mean()

    denom = math.sqrt(float(np.dot(sr, sr) * np.dot(tr, tr)))
    if denom <= 0:
        raise RuntimeError("Degenerate rank variance in loading permutation test.")

    null = np.empty(permutations, dtype=np.float64)
    rng = np.random.default_rng(seed)

    done = 0
    while done < permutations:
        k = min(LOADING_PERM_BATCH, permutations - done)
        vals = np.empty((k, len(tr)), dtype=np.float64)
        for r in range(k):
            vals[r, :] = tr[rng.permutation(len(tr))]
        null[done : done + k] = vals @ sr / denom
        done += k

    extreme = int(np.sum(np.abs(null) >= abs(observed_rho)))
    pvalue = (1 + extreme) / (permutations + 1)

    return pvalue, null, {
        "extreme_count": extreme,
        "null_mean": float(np.mean(null)),
        "null_sd": float(np.std(null, ddof=1)),
        "null_q025": float(np.quantile(null, 0.025)),
        "null_q975": float(np.quantile(null, 0.975)),
    }


def combo_seed(target_index: int, program_id: str, axis_offset: int) -> int:
    m = re.search(r"M(\d+)$", program_id)
    if not m:
        raise RuntimeError(f"Could not parse module number from {program_id}.")
    pnum = int(m.group(1))
    return BASE_SEED + target_index * 100_000 + pnum * 100 + axis_offset


def program_result_paths(target_name: str, program_id: str) -> tuple[Path, Path, Path]:
    d = OUT_DIR / "per_program" / target_name
    d.mkdir(parents=True, exist_ok=True)
    return (
        d / f"{program_id}_direct_result_v2.json",
        d / f"{program_id}_direct_nulls_v2.npz",
        d / f"{program_id}_evaluable_genes_v2.tsv",
    )


def run_target(
    target: dict,
    target_index: int,
    source: dict,
    weights: pd.DataFrame,
    device: torch.device,
) -> list[dict]:
    print("\n" + "=" * 138)
    print(f"[4/5] PRIMARY DIRECT PRESERVATION — {target['name']}")
    print("=" * 138)

    results = []
    source_gene_index = source["gene_index"]
    target_gene_index = target["gene_index"]

    for p_idx, program_id in enumerate(sorted(weights["program_id"].unique()), start=1):
        wdf = weights.loc[weights["program_id"] == program_id].copy()
        wdf["source_gene_order_within_program"] = (
            wdf["source_gene_order_within_program"].astype(int)
        )
        wdf = wdf.sort_values("source_gene_order_within_program")

        frozen_genes = [canon_symbol(x) for x in wdf["Hugo_Symbol"]]
        eligible_genes = [g for g in frozen_genes if g in target_gene_index]

        n_frozen = len(frozen_genes)
        n_eval = len(eligible_genes)
        coverage = n_eval / n_frozen
        assessable = n_eval >= MIN_EVALUABLE_GENES and coverage >= MIN_COVERAGE

        result_json, null_npz, genes_tsv = program_result_paths(
            target["name"], program_id
        )

        print(
            f"\n  [{p_idx:02d}/12] {program_id}: "
            f"eligible={n_eval}/{n_frozen} coverage={coverage:.3f} "
            f"assessable={'YES' if assessable else 'NO'}"
        )

        if not assessable:
            rec = {
                "script_version": SCRIPT_VERSION,
                "target": target["name"],
                "program_id": program_id,
                "n_frozen_genes": n_frozen,
                "n_evaluable_genes": n_eval,
                "coverage": coverage,
                "primary_assessable": False,
                "reason": "frozen coverage/evaluable-gene guard failed",
            }
            result_json.write_text(json.dumps(rec, indent=2), encoding="utf-8")
            results.append(rec)
            continue

        # Resume only if the exact v1 result + null distributions already exist.
        if result_json.exists() and null_npz.exists() and genes_tsv.exists():
            try:
                old = json.loads(result_json.read_text(encoding="utf-8"))
            except Exception:
                old = {}
            if (
                old.get("script_version") == SCRIPT_VERSION
                and old.get("target") == target["name"]
                and old.get("program_id") == program_id
                and bool(old.get("primary_assessable", False))
                and old.get("direct_permutations") == DIRECT_PERMUTATIONS
            ):
                print("    deterministic checkpoint found — reusing completed result")
                results.append(old)
                continue

        t0 = time.perf_counter()

        source_cols = np.array(
            [source_gene_index[g] for g in eligible_genes],
            dtype=np.int64,
        )
        target_cols = np.array(
            [target_gene_index[g] for g in eligible_genes],
            dtype=np.int64,
        )

        source_z = source["z"][:, source_cols]
        target_z = target["z"][:, target_cols]

        # Frozen source loadings/signs in exactly the same eligible source-gene order.
        weight_by_gene = wdf.set_index(
            wdf["Hugo_Symbol"].map(canon_symbol),
            drop=False,
        )
        source_loading = np.array(
            [
                float(weight_by_gene.loc[g, "source_pc1_loading"])
                for g in eligible_genes
            ],
            dtype=np.float64,
        )
        source_sign = np.array(
            [
                float(weight_by_gene.loc[g, "source_pc1_loading_sign"])
                for g in eligible_genes
            ],
            dtype=np.float64,
        )

        if np.any(~np.isfinite(source_loading)):
            raise RuntimeError(f"{target['name']} {program_id}: non-finite source loading.")

        source_corr, _ = correlation_matrix_gpu(source_z, device)
        target_corr, _ = correlation_matrix_gpu(target_z, device)

        src_gpu_diff = cpu_gpu_corr_validation(source_z, source_corr)
        tgt_gpu_diff = cpu_gpu_corr_validation(target_z, target_corr)

        tri_i, tri_j = np.triu_indices(n_eval, k=1)
        source_edges = source_corr[tri_i, tri_j]
        target_edges = target_corr[tri_i, tri_j]

        rho_edge = spearman_fast(source_edges, target_edges)

        target_loading, pcmeta = leading_pc1_loading(
            target_corr=target_corr,
            target_z=target_z,
            source_loading=source_loading,
            source_sign=source_sign,
            device=device,
            seed=combo_seed(target_index, program_id, 31),
        )
        rho_load = spearman_fast(source_loading, target_loading)

        print(
            f"    OBSERVED rho_edge={rho_edge:.6f}; "
            f"rho_load={rho_load:.6f}; "
            f"target PC1 variance={pcmeta['pc1_variance_explained']:.4f}"
        )

        # First true target inferential calculations: frozen 10,000 direct label permutations.
        p_edge, edge_null, edge_meta = exact_edge_label_permutation_test_gpu(
            source_edges=source_edges,
            target_corr=target_corr,
            observed_rho=rho_edge,
            permutations=DIRECT_PERMUTATIONS,
            seed=combo_seed(target_index, program_id, 41),
            device=device,
            progress_prefix=f"{target['name']} {program_id}",
        )

        p_load, load_null, load_meta = exact_loading_label_permutation_test(
            source_loading=source_loading,
            target_loading=target_loading,
            observed_rho=rho_load,
            permutations=DIRECT_PERMUTATIONS,
            seed=combo_seed(target_index, program_id, 51),
        )

        print(
            f"    direct permutation p_edge={p_edge:.7f}; "
            f"p_load={p_load:.7f}"
        )

        genes_df = pd.DataFrame(
            {
                "source_gene_order_within_evaluable_program": np.arange(
                    1, n_eval + 1
                ),
                "Hugo_Symbol": eligible_genes,
                "source_pc1_loading": source_loading,
                "source_pc1_loading_sign": source_sign,
            }
        )
        genes_df.to_csv(genes_tsv, sep="\t", index=False)

        np.savez_compressed(
            null_npz,
            edge_label_permutation_null=edge_null,
            loading_label_permutation_null=load_null,
        )

        rec = {
            "script_version": SCRIPT_VERSION,
            "target": target["name"],
            "program_id": program_id,
            "n_frozen_genes": n_frozen,
            "n_evaluable_genes": n_eval,
            "coverage": coverage,
            "n_edges": int(len(source_edges)),
            "primary_assessable": True,
            "rho_edge": rho_edge,
            "rho_load": rho_load,
            "direct_edge_permutation_p": p_edge,
            "direct_loading_permutation_p": p_load,
            "direct_permutations": DIRECT_PERMUTATIONS,
            "source_cpu_gpu_corr_validation_max_abs_diff": src_gpu_diff,
            "target_cpu_gpu_corr_validation_max_abs_diff": tgt_gpu_diff,
            **pcmeta,
            "edge_null_summary": edge_meta,
            "loading_null_summary": load_meta,
            "elapsed_seconds": time.perf_counter() - t0,
            "mapping_null_specificity_calculated": False,
            "primary_strong_partial_no_clear_classification_calculated": False,
        }
        result_json.write_text(json.dumps(rec, indent=2), encoding="utf-8")
        results.append(rec)

        print(f"    checkpoint written; elapsed={rec['elapsed_seconds']:.1f}s")

        del (
            source_z,
            target_z,
            source_corr,
            target_corr,
            source_edges,
            target_edges,
            target_loading,
            edge_null,
            load_null,
        )
        torch.cuda.empty_cache()
        gc.collect()

    return results


def aggregate_direct_results(results: list[dict]) -> pd.DataFrame:
    assessable = [r for r in results if bool(r.get("primary_assessable", False))]

    direct_p = []
    keys = []
    for r in assessable:
        direct_p.append(float(r["direct_edge_permutation_p"]))
        keys.append((r["target"], r["program_id"], "edge"))
        direct_p.append(float(r["direct_loading_permutation_p"]))
        keys.append((r["target"], r["program_id"], "loading"))

    # 04a froze BH WITHIN each target, not across targets.
    q_lookup = {}
    for target_name in sorted({r["target"] for r in assessable}):
        target_indices = [
            i for i, k in enumerate(keys) if k[0] == target_name
        ]
        target_p = [direct_p[i] for i in target_indices]
        target_q = bh_adjust(target_p)
        for i, q in zip(target_indices, target_q):
            q_lookup[keys[i]] = q

    rows = []
    for r in results:
        row = dict(r)
        if bool(r.get("primary_assessable", False)):
            q_edge = q_lookup[(r["target"], r["program_id"], "edge")]
            q_load = q_lookup[(r["target"], r["program_id"], "loading")]
            row["direct_edge_bh_q"] = q_edge
            row["direct_loading_bh_q"] = q_load
            row["direct_edge_positive_significant"] = int(
                r["rho_edge"] > 0 and q_edge < DIRECT_ALPHA
            )
            row["direct_loading_positive_significant"] = int(
                r["rho_load"] > 0 and q_load < DIRECT_ALPHA
            )
            row["direct_directionally_discordant"] = int(
                (r["rho_edge"] < 0 and q_edge < DIRECT_ALPHA)
                or (r["rho_load"] < 0 and q_load < DIRECT_ALPHA)
            )
            row["direct_both_positive_significant"] = int(
                row["direct_edge_positive_significant"]
                and row["direct_loading_positive_significant"]
            )
        else:
            row["direct_edge_bh_q"] = np.nan
            row["direct_loading_bh_q"] = np.nan
            row["direct_edge_positive_significant"] = 0
            row["direct_loading_positive_significant"] = 0
            row["direct_directionally_discordant"] = 0
            row["direct_both_positive_significant"] = 0

        # Deliberately no final class before the frozen matched-mapping specificity null.
        row["primary_classification"] = "PENDING_MAPPING_NULL_SPECIFICITY"
        rows.append(row)

    return pd.DataFrame(rows)


def main() -> None:
    print("=" * 140)
    print("Paper 4 / TCBB - corrected pooled direct edge/loading preservation under exact-variation evaluability")
    print("=" * 140)
    print(f"Script version: {SCRIPT_VERSION}")
    print()
    print("Scientific contract:")
    print("  Target outcomes/treatment loaded:                    NO")
    print("  Source program membership frozen:                    YES")
    print("  Target gene evaluability:                           exact finite + max>min correction frozen in 04i3")
    print("  Direct edge effect:                                  ALL evaluable edges")
    print(f"  Edge gene-label permutations:                        {DIRECT_PERMUTATIONS:,}")
    print(f"  Loading gene-label permutations:                     {DIRECT_PERMUTATIONS:,}")
    print("  BH family:                                            edge + loading tests WITHIN target")
    print("  Mapping-null specificity in this script:              NO")
    print("  Final Strong/Partial/No-clear class in this script:   NO")
    print("=" * 140)

    for p in [
        PRETARGET_AUDIT,
        EVALUABILITY_AMENDMENT,
        EXACT_VARIATION_CORRECTION,
        PRESERVATION_CONTRACT,
        UNIVERSE,
        WEIGHTS,
        TCGA_EXPR,
        SCANB_EXPR,
        SCANB_PRIMARY,
        METABRIC_EXPR,
    ]:
        require(p)

    pret = json.loads(PRETARGET_AUDIT.read_text(encoding="utf-8"))
    if pret.get("status") != "PASS_READY_FOR_FIRST_POOLED_TARGET_PRESERVATION":
        raise RuntimeError(f"04h pretarget audit is not PASS: {pret.get('status')}")

    eva = json.loads(EVALUABILITY_AMENDMENT.read_text(encoding="utf-8"))
    if eva.get("status") != "FROZEN_BEFORE_FIRST_TARGET_PRESERVATION_STATISTIC":
        raise RuntimeError("04i v2 evaluability amendment is not frozen as expected.")

    exact_corr = json.loads(EXACT_VARIATION_CORRECTION.read_text(encoding="utf-8"))
    if exact_corr.get("status") != "FROZEN_CORRECTION_BEFORE_RECOMPUTED_PRESERVATION":
        raise RuntimeError("04i3 exact-variation correction has unexpected status.")
    corrected_counts = exact_corr.get("counts", {})
    if int(corrected_counts.get("SCANB_corrected_evaluable", -1)) != 9220:
        raise RuntimeError("04i3 SCAN-B corrected evaluable count is not 9,220.")
    if int(corrected_counts.get("METABRIC_corrected_evaluable", -1)) != 8485:
        raise RuntimeError("04i3 METABRIC corrected evaluable count is not 8,485.")

    pcon = json.loads(PRESERVATION_CONTRACT.read_text(encoding="utf-8"))
    if pcon.get("status") != "FROZEN_BEFORE_FIRST_TARGET_PRESERVATION_STATISTIC":
        raise RuntimeError("04a preservation contract is not frozen as expected.")

    if int(pcon["direct_preservation"]["monte_carlo_permutations"]) != DIRECT_PERMUTATIONS:
        raise RuntimeError("Direct permutation count differs from frozen 04a contract.")

    device = ensure_cuda()
    prop = torch.cuda.get_device_properties(0)
    print(
        f"CUDA device: {prop.name}; VRAM={prop.total_memory/1024**3:.2f} GB; "
        f"capability={prop.major}.{prop.minor}"
    )

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    universe = pd.read_csv(UNIVERSE, sep="\t", dtype=str).fillna("")
    universe["source_gene_rank_by_MAD"] = universe["source_gene_rank_by_MAD"].astype(int)
    universe = universe.sort_values("source_gene_rank_by_MAD").reset_index(drop=True)
    if len(universe) != 10_000:
        raise RuntimeError(f"Expected 10,000 frozen source genes; found {len(universe)}.")

    universe_genes = [canon_symbol(x) for x in universe["Hugo_Symbol"]]
    if len(set(universe_genes)) != 10_000:
        raise RuntimeError("Frozen source universe is not unique after symbol canonicalization.")
    universe_set = set(universe_genes)

    weights = pd.read_csv(WEIGHTS, sep="\t", dtype=str).fillna("")
    if weights["program_id"].nunique() != 12:
        raise RuntimeError(f"Expected 12 frozen programs; found {weights['program_id'].nunique()}.")

    primary = pd.read_csv(SCANB_PRIMARY, sep="\t", dtype=str).fillna("")
    primary_titles = list(primary["primary_title"])
    if len(primary_titles) != 3273 or len(set(primary_titles)) != 3273:
        raise RuntimeError("Frozen SCAN-B primary manifest is not exactly 3,273 unique profiles.")

    source = load_tcga_source(universe)
    scanb = load_scanb_target(universe_set, primary_titles)
    met = load_metabric_target(universe_set)

    all_results = []
    all_results.extend(
        run_target(
            target=scanb,
            target_index=1,
            source=source,
            weights=weights,
            device=device,
        )
    )
    all_results.extend(
        run_target(
            target=met,
            target_index=2,
            source=source,
            weights=weights,
            device=device,
        )
    )

    print("\n[5/5] Applying frozen within-target BH adjustment to direct tests ...")
    final_df = aggregate_direct_results(all_results)

    summary_tsv = OUT_DIR / "primary_pooled_direct_preservation_summary_v2.tsv"
    final_df.to_csv(summary_tsv, sep="\t", index=False)

    # Compact JSON for the next frozen specificity stage.
    compact_cols = [
        "target",
        "program_id",
        "n_frozen_genes",
        "n_evaluable_genes",
        "coverage",
        "n_edges",
        "primary_assessable",
        "rho_edge",
        "rho_load",
        "direct_edge_permutation_p",
        "direct_loading_permutation_p",
        "direct_edge_bh_q",
        "direct_loading_bh_q",
        "direct_edge_positive_significant",
        "direct_loading_positive_significant",
        "direct_directionally_discordant",
        "direct_both_positive_significant",
        "primary_classification",
    ]
    available_cols = [c for c in compact_cols if c in final_df.columns]
    compact = final_df[available_cols].replace({np.nan: None}).to_dict(orient="records")

    master = {
        "script_version": SCRIPT_VERSION,
        "status": "CORRECTED_DIRECT_POOLED_PRESERVATION_COMPLETE_MAPPING_NULL_PENDING",
        "direct_permutations_per_axis": DIRECT_PERMUTATIONS,
        "bh_family": "within each target across assessable direct edge + loading tests",
        "mapping_null_specificity_calculated": False,
        "final_primary_classification_calculated": False,
        "results": compact,
    }
    master_json = OUT_DIR / "primary_pooled_direct_preservation_v2.json"
    master_json.write_text(json.dumps(master, indent=2), encoding="utf-8")

    print("\n" + "=" * 140)
    print("05a v2 CORRECTED PRIMARY POOLED DIRECT PRESERVATION: COMPLETE")
    print("=" * 140)

    for target_name in ["SCANB_GSE96058", "METABRIC"]:
        sub = final_df.loc[final_df["target"] == target_name]
        print(f"\n{target_name}:")
        for r in sub.itertuples(index=False):
            if not bool(r.primary_assessable):
                print(
                    f"  {r.program_id}: NOT ASSESSABLE "
                    f"(eligible={r.n_evaluable_genes}/{r.n_frozen_genes}, coverage={r.coverage:.3f})"
                )
                continue

            print(
                f"  {r.program_id}: "
                f"rho_edge={r.rho_edge:+.4f} q_edge={r.direct_edge_bh_q:.5g}; "
                f"rho_load={r.rho_load:+.4f} q_load={r.direct_loading_bh_q:.5g}; "
                f"both+sig={'YES' if r.direct_both_positive_significant else 'NO'}; "
                f"discordant={'YES' if r.direct_directionally_discordant else 'NO'}"
            )

    print()
    print("IMPORTANT:")
    print("  These are now the first true target preservation results.")
    print("  Do NOT assign Strong/Partial/No-clear yet.")
    print("  The frozen 1,000-panel matched-mapping edge-specificity null is still required.")
    print()
    print(f"Summary TSV:  {summary_tsv}")
    print(f"Master JSON:  {master_json}")
    print("=" * 140)


if __name__ == "__main__":
    main()
