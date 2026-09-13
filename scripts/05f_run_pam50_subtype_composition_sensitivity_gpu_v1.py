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
    raise RuntimeError("PyTorch is required for 05f.") from exc

try:
    from scipy.stats import rankdata
except Exception as exc:
    raise RuntimeError("scipy is required for Spearman statistics.") from exc


SCRIPT_VERSION = "05f-run-pam50-subtype-composition-sensitivity-gpu-v1-no-cli"

DATA_ROOT = Path(r"D:\paper4_tcbb_data")

EXEC_CONTRACT = (
    DATA_ROOT
    / "paper4_tcbb_pam50_sensitivity_execution_contract_v1"
    / "pam50_sensitivity_execution_contract_v1.json"
)
PAM50_CONTRACT = (
    DATA_ROOT
    / "paper4_tcbb_pam50_sensitivity_contract_v1"
    / "pam50_sensitivity_contract_v1.json"
)
TCGA_ALIGN = (
    DATA_ROOT
    / "paper4_tcbb_pam50_alignment_audit_v1"
    / "tcga_frozen_expression_pam50_alignment_v1.tsv"
)
SCANB_ALIGN = (
    DATA_ROOT
    / "paper4_tcbb_pam50_alignment_audit_v1"
    / "scanb_frozen_primary_pam50_alignment_v1.tsv"
)
MET_ALIGN = (
    DATA_ROOT
    / "paper4_tcbb_pam50_alignment_audit_v1"
    / "metabric_frozen_expression_subtype_alignment_v1.tsv"
)

CLASSIFICATION = (
    DATA_ROOT
    / "paper4_tcbb_final_mapping_specificity_null_v2"
    / "primary_pooled_final_classification_v2.tsv"
)
SPECIFICITY_ROOT = DATA_ROOT / "paper4_tcbb_final_mapping_specificity_null_v2"
DIRECT_ROOT = DATA_ROOT / "paper4_tcbb_primary_pooled_direct_preservation_v2"

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

OUT_DIR = DATA_ROOT / "paper4_tcbb_pam50_sensitivity_v1"

BASE_SEED = 20260917
BOOTSTRAP_REPLICATES = 1000
MAX_ATTEMPTS = 10000
CHECKPOINT_EVERY = 50
EDGE_BATCH = 1024
CPU_VALIDATION_EDGES = 25
CPU_GPU_TOL = 3e-10
REPLAY_TOL = 5e-10
PC1_RESIDUAL_TOL = 1e-10
MAX_POWER_ITER = 400

CANONICAL_SUBTYPES = ["LumA", "LumB", "Basal", "Her2", "Normal"]
COMPUTED_SUBTYPES = ["LumA", "LumB", "Basal", "Her2"]
SUBTYPE_CODE = {"LumA": 1, "LumB": 2, "Basal": 3, "Her2": 4, "Normal": 5}
TARGET_INDEX = {"SCANB_GSE96058": 1, "METABRIC": 2}
TIER = {"LumA": "CORE", "LumB": "CORE", "Basal": "CORE", "Her2": "LIMITED_N", "Normal": "NOT_INDIVIDUALLY_ASSESSABLE"}


def require(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Required file does not exist: {path}")


def canon_symbol(x: object) -> str:
    if x is None:
        return ""
    return " ".join(str(x).strip().split()).upper()


def module_number(program_id: str) -> int:
    m = re.search(r"M(\d+)$", str(program_id))
    if not m:
        raise RuntimeError(f"Cannot parse module number from {program_id}.")
    return int(m.group(1))


def boolish(x: object) -> bool:
    if isinstance(x, bool):
        return x
    return str(x).strip().lower() in {"1", "true", "yes", "y"}


def centered_rank_cpu(x: np.ndarray) -> tuple[np.ndarray, float]:
    r = rankdata(np.asarray(x, dtype=np.float64), method="average")
    r -= r.mean()
    ss = float(np.dot(r, r))
    if ss <= 0:
        raise RuntimeError("Degenerate rank vector.")
    return r, math.sqrt(ss)


def spearman_cpu(x: np.ndarray, y: np.ndarray) -> float:
    xr, xn = centered_rank_cpu(x)
    yr, yn = centered_rank_cpu(y)
    return float(np.dot(xr, yr) / (xn * yn))


def gpu_centered_average_ranks(values: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
    sorted_vals, order = torch.sort(values)
    n = sorted_vals.numel()

    new_group = torch.empty(n, dtype=torch.bool, device=values.device)
    new_group[0] = True
    if n > 1:
        new_group[1:] = sorted_vals[1:] != sorted_vals[:-1]

    group_id = torch.cumsum(new_group.to(torch.int64), dim=0) - 1
    counts = torch.bincount(group_id)
    starts0 = torch.cumsum(counts, dim=0) - counts
    avg_rank = starts0.to(torch.float64) + (counts.to(torch.float64) - 1.0) / 2.0 + 1.0
    ranks_sorted = avg_rank[group_id]

    ranks = torch.empty(n, dtype=torch.float64, device=values.device)
    ranks[order] = ranks_sorted
    centered = ranks - (n + 1.0) / 2.0
    norm = torch.linalg.vector_norm(centered)

    del sorted_vals, order, new_group, group_id, counts, starts0, avg_rank, ranks_sorted, ranks
    return centered, norm


def gpu_spearman(x: torch.Tensor, y: torch.Tensor) -> float:
    xr, xn = gpu_centered_average_ranks(x)
    yr, yn = gpu_centered_average_ranks(y)
    rho = torch.dot(xr, yr) / (xn * yn)
    value = float(rho.detach().cpu())
    del xr, xn, yr, yn, rho
    return value


def standardize_np(x: np.ndarray) -> np.ndarray:
    if not np.isfinite(x).all():
        raise RuntimeError("Non-finite values reached standardization.")
    mu = x.mean(axis=0, keepdims=True)
    sd = x.std(axis=0, ddof=1, keepdims=True)
    if np.any(~np.isfinite(sd)) or np.any(sd <= 0):
        raise RuntimeError("Zero/non-finite SD reached standardization.")
    return (x - mu) / sd


def exact_variable_np(x: np.ndarray) -> np.ndarray:
    finite = np.isfinite(x).all(axis=0)
    out = np.zeros(x.shape[1], dtype=bool)
    if finite.any():
        xmin = np.min(x[:, finite], axis=0)
        xmax = np.max(x[:, finite], axis=0)
        out[finite] = xmax > xmin
    return out


def full_corr_from_x(x: np.ndarray, device: torch.device) -> tuple[np.ndarray, float]:
    z = standardize_np(x)
    zt = torch.as_tensor(z, dtype=torch.float64, device=device)
    corr_t = (zt.T @ zt) / (z.shape[0] - 1)
    corr = corr_t.detach().cpu().numpy()

    k = min(30, x.shape[1])
    cpu = (z[:, :k].T @ z[:, :k]) / (z.shape[0] - 1)
    diff = float(np.max(np.abs(cpu - corr[:k, :k])))
    if diff > CPU_GPU_TOL:
        raise RuntimeError(f"CPU/GPU correlation validation failed: max|Δ|={diff:.3e}")

    np.fill_diagonal(corr, 1.0)
    del zt, corr_t
    torch.cuda.empty_cache()
    return corr, diff


def pc1_loading(
    corr: np.ndarray,
    x: np.ndarray,
    source_loading: np.ndarray,
    source_sign: np.ndarray,
    random_start: np.ndarray,
    device: torch.device,
) -> tuple[np.ndarray, dict]:
    z = standardize_np(x)
    ct = torch.as_tensor(corr, dtype=torch.float64, device=device)
    eps = torch.finfo(torch.float64).eps

    starts = [source_loading, source_sign, random_start]
    candidates = []

    for start_idx, start in enumerate(starts):
        v = torch.as_tensor(start, dtype=torch.float64, device=device)
        v = v / torch.linalg.vector_norm(v)
        residual_value = float("inf")

        for it in range(1, MAX_POWER_ITER + 1):
            w = ct @ v
            wn = torch.linalg.vector_norm(w)
            if not torch.isfinite(wn) or wn <= 0:
                break
            nv = w / wn
            if torch.dot(nv, v) < 0:
                nv = -nv
            v = nv

            if it % 5 == 0 or it == MAX_POWER_ITER:
                av = ct @ v
                lam = torch.dot(v, av)
                residual = torch.linalg.vector_norm(av - lam * v) / (torch.abs(lam) + eps)
                residual_value = float(residual.detach().cpu())
                if residual_value <= PC1_RESIDUAL_TOL:
                    break

        av = ct @ v
        lam = float(torch.dot(v, av).detach().cpu())
        candidates.append((lam, residual_value, it, start_idx, v.detach().cpu().numpy()))

    candidates.sort(key=lambda q: q[0], reverse=True)
    lam, residual, iters, start_idx, loading = candidates[0]
    exact_fallback = False

    if residual > PC1_RESIDUAL_TOL:
        evals, evecs = torch.linalg.eigh(ct)
        v = evecs[:, -1]
        lam_t = evals[-1]
        av = ct @ v
        residual_t = torch.linalg.vector_norm(av - lam_t * v) / (torch.abs(lam_t) + eps)
        loading = v.detach().cpu().numpy()
        lam = float(lam_t.detach().cpu())
        residual = float(residual_t.detach().cpu())
        exact_fallback = True

    scores = z @ loading
    sign_score = np.mean(z * source_sign[None, :], axis=1)
    orient = float(np.corrcoef(scores, sign_score)[0, 1])
    if not np.isfinite(orient):
        raise RuntimeError("PC1 orientation correlation is non-finite.")
    if orient < 0:
        loading = -loading
        orient = -orient

    meta = {
        "eigen_residual": residual,
        "exact_eigh_fallback": int(exact_fallback),
        "power_start_index": int(start_idx),
        "power_iterations": int(iters),
        "orientation_correlation": orient,
        "variance_explained": float(lam / np.trace(corr)),
    }

    del ct
    torch.cuda.empty_cache()
    return loading, meta


def load_source_raw(universe: pd.DataFrame) -> dict:
    raw = pd.read_csv(TCGA_EXPR, sep="\t", low_memory=False)
    row_idx = universe["source_row_index_0based"].astype(int).to_numpy()
    selected = raw.iloc[row_idx, :]

    expected = [canon_symbol(x) for x in universe["Hugo_Symbol"]]
    observed = [canon_symbol(x) for x in selected["Hugo_Symbol"]]
    if expected != observed:
        raise RuntimeError("TCGA frozen source-universe replay mismatch.")

    values = selected.iloc[:, 2:].to_numpy(dtype=np.float64)
    if np.any(values < 0) or not np.isfinite(values).all():
        raise RuntimeError("TCGA RSEM contains invalid values.")

    x = np.log2(values + 1.0).T
    return {
        "samples": list(selected.columns[2:]),
        "genes": expected,
        "gene_index": {g: i for i, g in enumerate(expected)},
        "x": x,
    }


def load_scanb_raw(universe_set: set[str], primary_titles: list[str]) -> dict:
    parts = []
    for chunk in pd.read_csv(SCANB_EXPR, compression="gzip", chunksize=1500, low_memory=False):
        symbol_col = chunk.columns[0]
        symbols = chunk[symbol_col].map(canon_symbol)
        mask = symbols.isin(universe_set)
        if not mask.any():
            continue

        sub = chunk.loc[mask, [symbol_col] + primary_titles]
        arr = sub[primary_titles].to_numpy(dtype=np.float64)
        fpkm = np.maximum(np.exp2(arr) - 0.1, 0.0)
        transformed = np.log2(fpkm + 1.0)

        tmp = pd.DataFrame(transformed, columns=primary_titles)
        tmp.insert(0, "Hugo_Symbol", symbols.loc[mask].to_numpy())
        parts.append(tmp)

    df = pd.concat(parts, ignore_index=True)
    df = df.groupby("Hugo_Symbol", sort=False, as_index=True).mean(numeric_only=True)
    x = df.to_numpy(dtype=np.float64).T

    finite = np.isfinite(x).all(axis=0)
    xmin = np.full(x.shape[1], np.nan)
    xmax = np.full(x.shape[1], np.nan)
    xmin[finite] = np.min(x[:, finite], axis=0)
    xmax[finite] = np.max(x[:, finite], axis=0)
    eligible = finite & (xmax > xmin)

    genes_all = list(df.index)
    genes = [g for g, ok in zip(genes_all, eligible) if ok]
    x = x[:, eligible]
    if len(genes) != 9220:
        raise RuntimeError(f"Expected 9,220 corrected SCAN-B genes, found {len(genes)}.")

    return {
        "samples": primary_titles,
        "genes": genes,
        "gene_index": {g: i for i, g in enumerate(genes)},
        "x": x,
    }


def load_metabric_raw(universe_set: set[str]) -> dict:
    parts = []
    sample_cols = None
    for chunk in pd.read_csv(METABRIC_EXPR, sep="\t", chunksize=1500, low_memory=False):
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

    df = pd.concat(parts, ignore_index=True)
    df = df.groupby("Hugo_Symbol", sort=False, as_index=True).mean(numeric_only=True)
    x = df.to_numpy(dtype=np.float64).T

    finite = np.isfinite(x).all(axis=0)
    xmin = np.full(x.shape[1], np.nan)
    xmax = np.full(x.shape[1], np.nan)
    xmin[finite] = np.min(x[:, finite], axis=0)
    xmax[finite] = np.max(x[:, finite], axis=0)
    eligible = finite & (xmax > xmin)

    genes_all = list(df.index)
    genes = [g for g, ok in zip(genes_all, eligible) if ok]
    x = x[:, eligible]
    if len(genes) != 8485:
        raise RuntimeError(f"Expected 8,485 corrected METABRIC genes, found {len(genes)}.")

    return {
        "samples": sample_cols,
        "genes": genes,
        "gene_index": {g: i for i, g in enumerate(genes)},
        "x": x,
    }


def load_labels(samples: list[str], target_name: str) -> np.ndarray:
    if target_name == "TCGA_BRCA":
        df = pd.read_csv(TCGA_ALIGN, sep="\t", dtype=str).fillna("")
        if len(df) != 1082:
            raise RuntimeError("TCGA PAM50 alignment row count is not 1,082.")
        mapping = dict(zip(df["expression_sample_id"], df["pam50_subtype"]))
    elif target_name == "SCANB_GSE96058":
        df = pd.read_csv(SCANB_ALIGN, sep="\t", dtype=str).fillna("")
        if len(df) != 3273:
            raise RuntimeError("SCAN-B PAM50 alignment row count is not 3,273.")
        mapping = dict(zip(df["sample_title"], df["pam50_subtype"]))
    elif target_name == "METABRIC":
        df = pd.read_csv(MET_ALIGN, sep="\t", dtype=str).fillna("")
        if len(df) != 1980:
            raise RuntimeError("METABRIC PAM50 alignment row count is not 1,980.")
        mapping = dict(zip(df["expression_sample_id"], df["subtype_canonical"]))
    else:
        raise RuntimeError(f"Unknown label target: {target_name}")

    labels = np.array([mapping.get(s, "") for s in samples], dtype=object)
    return labels


def validate_label_counts(source_labels: np.ndarray, scanb_labels: np.ndarray, met_labels: np.ndarray) -> None:
    expected = {
        "TCGA_BRCA": {"LumA": 417, "LumB": 188, "Basal": 139, "Her2": 67, "Normal": 23},
        "SCANB_GSE96058": {"LumA": 1657, "LumB": 729, "Basal": 339, "Her2": 327, "Normal": 221},
        "METABRIC": {"LumA": 700, "LumB": 475, "Basal": 209, "Her2": 224, "Normal": 148},
    }
    for name, labels in [("TCGA_BRCA", source_labels), ("SCANB_GSE96058", scanb_labels), ("METABRIC", met_labels)]:
        got = {s: int(np.sum(labels == s)) for s in CANONICAL_SUBTYPES}
        if got != expected[name]:
            raise RuntimeError(f"{name} PAM50 aligned counts mismatch: {got} vs {expected[name]}")


def program_genes(target_name: str, program_id: str) -> list[str]:
    path = DIRECT_ROOT / "per_program" / target_name / f"{program_id}_evaluable_genes_v2.tsv"
    require(path)
    df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    return [canon_symbol(x) for x in df["Hugo_Symbol"]]


def source_loading_signs(weights: pd.DataFrame, program_id: str, genes: list[str]) -> tuple[np.ndarray, np.ndarray]:
    w = weights.loc[weights["program_id"] == program_id].copy()
    w["Hugo_Symbol"] = w["Hugo_Symbol"].map(canon_symbol)
    w = w.set_index("Hugo_Symbol", drop=False)
    loading = np.array([float(w.loc[g, "source_pc1_loading"]) for g in genes], dtype=np.float64)
    sign = np.array([float(w.loc[g, "source_pc1_loading_sign"]) for g in genes], dtype=np.float64)
    return loading, sign


def fixed_edge_subset(target_name: str, program_id: str) -> tuple[np.ndarray, np.ndarray]:
    path = SPECIFICITY_ROOT / "per_program" / target_name / f"{program_id}_mapping_specificity_null_v2.npz"
    require(path)
    z = np.load(path)
    return (
        np.asarray(z["edge_slot_i_0based"], dtype=np.int64),
        np.asarray(z["edge_slot_j_0based"], dtype=np.int64),
    )


def subtype_seed(target_name: str, program_id: str, subtype: str, attempt_id: int) -> np.random.Generator:
    ss = np.random.SeedSequence([
        BASE_SEED,
        TARGET_INDEX[target_name],
        module_number(program_id),
        SUBTYPE_CODE[subtype],
        101,
        int(attempt_id),
    ])
    return np.random.default_rng(ss)


def composition_seed(target_name: str, program_id: str, attempt_id: int) -> np.random.Generator:
    ss = np.random.SeedSequence([
        BASE_SEED,
        TARGET_INDEX[target_name],
        module_number(program_id),
        201,
        int(attempt_id),
    ])
    return np.random.default_rng(ss)


def pc1_random_start(target_name: str, program_id: str, subtype: str, cohort_code: int, m: int) -> np.ndarray:
    ss = np.random.SeedSequence([
        BASE_SEED,
        TARGET_INDEX[target_name],
        module_number(program_id),
        SUBTYPE_CODE[subtype],
        301,
        int(cohort_code),
    ])
    return np.random.default_rng(ss).normal(size=m)


def support_variable_torch(x_t: torch.Tensor, counts_t: torch.Tensor) -> tuple[bool, np.ndarray]:
    support = torch.nonzero(counts_t > 0, as_tuple=False).flatten()
    xs = x_t.index_select(0, support)
    xmin = torch.amin(xs, dim=0)
    xmax = torch.amax(xs, dim=0)
    bad_t = (~torch.isfinite(xmin)) | (~torch.isfinite(xmax)) | (xmax <= xmin)
    bad = bad_t.detach().cpu().numpy().astype(bool)
    ok = not bool(torch.any(bad_t))
    del support, xs, xmin, xmax, bad_t
    return ok, bad


def selected_edge_corr_weighted(
    x_t: torch.Tensor,
    counts_t: torch.Tensor,
    ii_t: torch.Tensor,
    jj_t: torch.Tensor,
) -> tuple[torch.Tensor | None, np.ndarray]:
    ok, bad = support_variable_torch(x_t, counts_t)
    if not ok:
        return None, bad

    total = torch.sum(counts_t)
    w = counts_t / total
    mu = torch.sum(w[:, None] * x_t, dim=0)
    xc = x_t - mu[None, :]
    var = torch.sum(w[:, None] * xc * xc, dim=0)

    if bool(torch.any(~torch.isfinite(var))) or bool(torch.any(var <= 0)):
        bad2 = ((~torch.isfinite(var)) | (var <= 0)).detach().cpu().numpy().astype(bool)
        return None, bad2

    e = ii_t.numel()
    out = torch.empty(e, dtype=torch.float64, device=x_t.device)

    for start in range(0, e, EDGE_BATCH):
        stop = min(e, start + EDGE_BATCH)
        ib = ii_t[start:stop]
        jb = jj_t[start:stop]
        xi = xc.index_select(1, ib)
        xj = xc.index_select(1, jb)
        cov = torch.sum(w[:, None] * xi * xj, dim=0)
        out[start:stop] = cov / torch.sqrt(var[ib] * var[jb])
        del ib, jb, xi, xj, cov

    del w, mu, xc, var
    return out, bad


def selected_edge_corr_residualized(
    x_t: torch.Tensor,
    counts_t: torch.Tensor,
    labels_code_t: torch.Tensor,
    n_groups: int,
    ii_t: torch.Tensor,
    jj_t: torch.Tensor,
) -> tuple[torch.Tensor | None, np.ndarray]:
    means = []
    for g in range(n_groups):
        mask = labels_code_t == g
        cg = counts_t[mask]
        total = torch.sum(cg)
        if total <= 0:
            return None, np.ones(x_t.shape[1], dtype=bool)
        xg = x_t[mask]
        means.append(torch.sum(cg[:, None] * xg, dim=0) / total)

    means_t = torch.stack(means, dim=0)
    resid = x_t - means_t.index_select(0, labels_code_t)
    edges, bad = selected_edge_corr_weighted(resid, counts_t, ii_t, jj_t)
    del means_t, resid
    return edges, bad


def explicit_edges(x: np.ndarray, idx: np.ndarray, ii: np.ndarray, jj: np.ndarray, k: int) -> np.ndarray:
    xb = x[idx]
    out = np.empty(k, dtype=np.float64)
    for q in range(k):
        out[q] = np.corrcoef(xb[:, ii[q]], xb[:, jj[q]])[0, 1]
    return out


def residualize_np(x: np.ndarray, labels: np.ndarray) -> np.ndarray:
    r = x.copy()
    for subtype in CANONICAL_SUBTYPES:
        idx = np.where(labels == subtype)[0]
        if len(idx) == 0:
            raise RuntimeError(f"No {subtype} samples available for residualization.")
        r[idx, :] -= np.mean(x[idx, :], axis=0, keepdims=True)
    return r


def stratified_resample_counts(labels: np.ndarray, rng: np.random.Generator) -> tuple[np.ndarray, np.ndarray]:
    counts = np.zeros(len(labels), dtype=np.float64)
    explicit = []
    for subtype in CANONICAL_SUBTYPES:
        pos = np.where(labels == subtype)[0]
        draw = rng.choice(pos, size=len(pos), replace=True)
        counts += np.bincount(draw, minlength=len(labels)).astype(np.float64)
        explicit.append(draw)
    return counts, np.concatenate(explicit)


def subtype_paths(target_name: str, program_id: str, subtype: str) -> dict[str, Path]:
    d = OUT_DIR / "subtype" / target_name / program_id
    d.mkdir(parents=True, exist_ok=True)
    return {
        "checkpoint": d / f"{subtype}_edge_bootstrap_checkpoint_v1.npz",
        "result": d / f"{subtype}_subtype_sensitivity_result_v1.json",
        "degenerate": d / f"{subtype}_degenerate_genes_v1.tsv",
    }


def composition_paths(target_name: str, program_id: str) -> dict[str, Path]:
    d = OUT_DIR / "composition" / target_name
    d.mkdir(parents=True, exist_ok=True)
    return {
        "checkpoint": d / f"{program_id}_composition_bootstrap_checkpoint_v1.npz",
        "result": d / f"{program_id}_composition_sensitivity_result_v1.json",
        "degenerate": d / f"{program_id}_composition_degenerate_genes_v1.tsv",
    }


def run_subtype_cell(
    target_name: str,
    program_id: str,
    subtype: str,
    source_x: np.ndarray,
    target_x: np.ndarray,
    source_loading: np.ndarray,
    source_sign: np.ndarray,
    ii: np.ndarray,
    jj: np.ndarray,
    device: torch.device,
) -> dict:
    paths = subtype_paths(target_name, program_id, subtype)
    if paths["result"].exists():
        old = json.loads(paths["result"].read_text(encoding="utf-8"))
        if old.get("script_version") == SCRIPT_VERSION and old.get("status") in {
            "COMPLETE", "NOT_ESTIMABLE_EXACT_VARIATION", "BOOTSTRAP_CI_NOT_ESTIMABLE"
        }:
            print("      completed checkpoint found — reusing")
            return old

    bad_source = ~exact_variable_np(source_x)
    bad_target = ~exact_variable_np(target_x)

    if bad_source.any() or bad_target.any():
        pd.DataFrame({
            "gene_index_0based": np.arange(source_x.shape[1]),
            "source_exact_constant": bad_source.astype(int),
            "target_exact_constant": bad_target.astype(int),
        }).loc[lambda d: (d["source_exact_constant"] == 1) | (d["target_exact_constant"] == 1)].to_csv(
            paths["degenerate"], sep="\t", index=False
        )
        result = {
            "script_version": SCRIPT_VERSION,
            "status": "NOT_ESTIMABLE_EXACT_VARIATION",
            "target": target_name,
            "program_id": program_id,
            "subtype": subtype,
            "tier": TIER[subtype],
            "source_n": int(source_x.shape[0]),
            "target_n": int(target_x.shape[0]),
            "n_genes": int(source_x.shape[1]),
            "source_exact_constant_genes": int(bad_source.sum()),
            "target_exact_constant_genes": int(bad_target.sum()),
            "primary_classification_changed": False,
        }
        paths["result"].write_text(json.dumps(result, indent=2), encoding="utf-8")
        return result

    source_corr, source_diff = full_corr_from_x(source_x, device)
    target_corr, target_diff = full_corr_from_x(target_x, device)

    tri_i, tri_j = np.triu_indices(source_x.shape[1], k=1)
    rho_edge = spearman_cpu(
        source_corr[tri_i, tri_j],
        target_corr[tri_i, tri_j],
    )

    m = source_x.shape[1]
    source_pc1, source_meta = pc1_loading(
        source_corr,
        source_x,
        source_loading,
        source_sign,
        pc1_random_start(target_name, program_id, subtype, 1, m),
        device,
    )
    target_pc1, target_meta = pc1_loading(
        target_corr,
        target_x,
        source_loading,
        source_sign,
        pc1_random_start(target_name, program_id, subtype, 2, m),
        device,
    )
    rho_load = spearman_cpu(source_pc1, target_pc1)

    del source_corr, target_corr
    gc.collect()

    source_t = torch.as_tensor(source_x, dtype=torch.float64, device=device)
    target_t = torch.as_tensor(target_x, dtype=torch.float64, device=device)
    ii_t = torch.as_tensor(ii, dtype=torch.long, device=device)
    jj_t = torch.as_tensor(jj, dtype=torch.long, device=device)

    rho_boot = np.full(BOOTSTRAP_REPLICATES, np.nan, dtype=np.float64)
    completed = 0
    attempts = 0
    source_deg = np.zeros(m, dtype=np.int64)
    target_deg = np.zeros(m, dtype=np.int64)

    if paths["checkpoint"].exists():
        ck = np.load(paths["checkpoint"])
        completed = int(ck["completed"])
        attempts = int(ck["attempts"])
        rho_boot[:] = ck["rho_boot"]
        source_deg[:] = ck["source_deg"]
        target_deg[:] = ck["target_deg"]
        print(f"      resuming bootstrap valid={completed}/{BOOTSTRAP_REPLICATES}; attempts={attempts}")

    engine_validated = completed > 0
    t0 = time.perf_counter()

    while completed < BOOTSTRAP_REPLICATES and attempts < MAX_ATTEMPTS:
        attempts += 1
        rng = subtype_seed(target_name, program_id, subtype, attempts)
        src_idx = rng.integers(0, source_x.shape[0], size=source_x.shape[0], endpoint=False)
        tgt_idx = rng.integers(0, target_x.shape[0], size=target_x.shape[0], endpoint=False)
        src_counts = np.bincount(src_idx, minlength=source_x.shape[0]).astype(np.float64)
        tgt_counts = np.bincount(tgt_idx, minlength=target_x.shape[0]).astype(np.float64)

        src_counts_t = torch.as_tensor(src_counts, dtype=torch.float64, device=device)
        tgt_counts_t = torch.as_tensor(tgt_counts, dtype=torch.float64, device=device)

        src_edges_t, src_bad = selected_edge_corr_weighted(source_t, src_counts_t, ii_t, jj_t)
        if src_edges_t is None:
            source_deg += src_bad.astype(np.int64)
            del src_counts_t, tgt_counts_t
            continue

        tgt_edges_t, tgt_bad = selected_edge_corr_weighted(target_t, tgt_counts_t, ii_t, jj_t)
        if tgt_edges_t is None:
            target_deg += tgt_bad.astype(np.int64)
            del src_counts_t, tgt_counts_t, src_edges_t
            continue

        if not engine_validated:
            k = min(CPU_VALIDATION_EDGES, len(ii))
            src_ref = explicit_edges(source_x, src_idx, ii, jj, k)
            tgt_ref = explicit_edges(target_x, tgt_idx, ii, jj, k)
            src_got = src_edges_t[:k].detach().cpu().numpy()
            tgt_got = tgt_edges_t[:k].detach().cpu().numpy()
            diff = max(
                float(np.max(np.abs(src_ref - src_got))),
                float(np.max(np.abs(tgt_ref - tgt_got))),
            )
            if diff > CPU_GPU_TOL:
                raise RuntimeError(
                    f"{target_name} {program_id} {subtype}: bootstrap selected-edge CPU/GPU mismatch {diff:.3e}"
                )
            print(f"      bootstrap selected-edge CPU↔GPU validation max|Δ|={diff:.3e}")
            engine_validated = True

        rho_boot[completed] = gpu_spearman(src_edges_t, tgt_edges_t)
        completed += 1

        del src_counts_t, tgt_counts_t, src_edges_t, tgt_edges_t

        if completed % CHECKPOINT_EVERY == 0 or completed == BOOTSTRAP_REPLICATES:
            np.savez_compressed(
                paths["checkpoint"],
                completed=np.array(completed, dtype=np.int64),
                attempts=np.array(attempts, dtype=np.int64),
                rho_boot=rho_boot,
                source_deg=source_deg,
                target_deg=target_deg,
            )
            print(
                f"      bootstrap {completed:4d}/{BOOTSTRAP_REPLICATES}; attempts={attempts:5d}; "
                f"median rho={np.nanmedian(rho_boot[:completed]):+.4f}; elapsed={time.perf_counter()-t0:.1f}s"
            )

    deg_rows = []
    for g in range(m):
        if source_deg[g] or target_deg[g]:
            deg_rows.append({
                "gene_index_0based": g,
                "source_degenerate_attempts": int(source_deg[g]),
                "target_degenerate_attempts": int(target_deg[g]),
            })
    pd.DataFrame(
        deg_rows,
        columns=["gene_index_0based", "source_degenerate_attempts", "target_degenerate_attempts"],
    ).to_csv(paths["degenerate"], sep="\t", index=False)

    ci_ok = completed == BOOTSTRAP_REPLICATES
    result = {
        "script_version": SCRIPT_VERSION,
        "status": "COMPLETE" if ci_ok else "BOOTSTRAP_CI_NOT_ESTIMABLE",
        "target": target_name,
        "program_id": program_id,
        "subtype": subtype,
        "tier": TIER[subtype],
        "source_n": int(source_x.shape[0]),
        "target_n": int(target_x.shape[0]),
        "n_genes": int(m),
        "n_all_edges": int(len(tri_i)),
        "n_bootstrap_edges": int(len(ii)),
        "rho_edge_all_edges": rho_edge,
        "rho_load_descriptive": rho_load,
        "source_full_corr_cpu_gpu_max_abs_diff": source_diff,
        "target_full_corr_cpu_gpu_max_abs_diff": target_diff,
        "source_pc1": source_meta,
        "target_pc1": target_meta,
        "bootstrap_valid_replicates": int(completed),
        "bootstrap_attempts": int(attempts),
        "bootstrap_invalid_attempts": int(attempts - completed),
        "bootstrap_ci_low": float(np.quantile(rho_boot[:completed], 0.025)) if completed else None,
        "bootstrap_median": float(np.median(rho_boot[:completed])) if completed else None,
        "bootstrap_ci_high": float(np.quantile(rho_boot[:completed], 0.975)) if completed else None,
        "formal_p_value": None,
        "primary_classification_changed": False,
    }
    paths["result"].write_text(json.dumps(result, indent=2), encoding="utf-8")

    del source_t, target_t, ii_t, jj_t
    torch.cuda.empty_cache()
    return result


def composition_point_and_bootstrap(
    target_name: str,
    program_id: str,
    row: dict,
    source_x: np.ndarray,
    source_labels: np.ndarray,
    target_x: np.ndarray,
    target_labels: np.ndarray,
    ii: np.ndarray,
    jj: np.ndarray,
    device: torch.device,
) -> dict:
    paths = composition_paths(target_name, program_id)
    if paths["result"].exists():
        old = json.loads(paths["result"].read_text(encoding="utf-8"))
        if old.get("script_version") == SCRIPT_VERSION and old.get("status") in {
            "COMPLETE", "B_NOT_ESTIMABLE_EXACT_VARIATION", "C_NOT_ESTIMABLE_EXACT_VARIATION", "BOOTSTRAP_CI_NOT_ESTIMABLE"
        }:
            print("      composition checkpoint found — reusing")
            return old

    m = source_x.shape[1]
    bad_B_source = ~exact_variable_np(source_x)
    bad_B_target = ~exact_variable_np(target_x)
    if bad_B_source.any() or bad_B_target.any():
        pd.DataFrame({
            "gene_index_0based": np.arange(m),
            "B_source_exact_constant": bad_B_source.astype(int),
            "B_target_exact_constant": bad_B_target.astype(int),
        }).loc[lambda d: (d["B_source_exact_constant"] == 1) | (d["B_target_exact_constant"] == 1)].to_csv(
            paths["degenerate"], sep="\t", index=False
        )
        result = {
            "script_version": SCRIPT_VERSION,
            "status": "B_NOT_ESTIMABLE_EXACT_VARIATION",
            "target": target_name,
            "program_id": program_id,
            "rho_A_full_primary": float(row["rho_edge"]),
            "B_source_exact_constant_genes": int(bad_B_source.sum()),
            "B_target_exact_constant_genes": int(bad_B_target.sum()),
            "primary_classification_changed": False,
        }
        paths["result"].write_text(json.dumps(result, indent=2), encoding="utf-8")
        return result

    B_source_corr, B_source_diff = full_corr_from_x(source_x, device)
    B_target_corr, B_target_diff = full_corr_from_x(target_x, device)
    tri_i, tri_j = np.triu_indices(m, k=1)
    rho_B = spearman_cpu(B_source_corr[tri_i, tri_j], B_target_corr[tri_i, tri_j])

    source_resid = residualize_np(source_x, source_labels)
    target_resid = residualize_np(target_x, target_labels)
    bad_C_source = ~exact_variable_np(source_resid)
    bad_C_target = ~exact_variable_np(target_resid)

    rho_A = float(row["rho_edge"])

    if bad_C_source.any() or bad_C_target.any():
        pd.DataFrame({
            "gene_index_0based": np.arange(m),
            "C_source_exact_constant": bad_C_source.astype(int),
            "C_target_exact_constant": bad_C_target.astype(int),
        }).loc[lambda d: (d["C_source_exact_constant"] == 1) | (d["C_target_exact_constant"] == 1)].to_csv(
            paths["degenerate"], sep="\t", index=False
        )
        result = {
            "script_version": SCRIPT_VERSION,
            "status": "C_NOT_ESTIMABLE_EXACT_VARIATION",
            "target": target_name,
            "program_id": program_id,
            "rho_A_full_primary": rho_A,
            "rho_B_complete_case": rho_B,
            "B_minus_A": rho_B - rho_A,
            "C_source_exact_constant_genes": int(bad_C_source.sum()),
            "C_target_exact_constant_genes": int(bad_C_target.sum()),
            "primary_classification_changed": False,
        }
        paths["result"].write_text(json.dumps(result, indent=2), encoding="utf-8")
        return result

    C_source_corr, C_source_diff = full_corr_from_x(source_resid, device)
    C_target_corr, C_target_diff = full_corr_from_x(target_resid, device)
    rho_C = spearman_cpu(C_source_corr[tri_i, tri_j], C_target_corr[tri_i, tri_j])

    del B_source_corr, B_target_corr, C_source_corr, C_target_corr, source_resid, target_resid
    gc.collect()

    source_t = torch.as_tensor(source_x, dtype=torch.float64, device=device)
    target_t = torch.as_tensor(target_x, dtype=torch.float64, device=device)
    ii_t = torch.as_tensor(ii, dtype=torch.long, device=device)
    jj_t = torch.as_tensor(jj, dtype=torch.long, device=device)

    source_label_code = np.array([CANONICAL_SUBTYPES.index(x) for x in source_labels], dtype=np.int64)
    target_label_code = np.array([CANONICAL_SUBTYPES.index(x) for x in target_labels], dtype=np.int64)
    source_code_t = torch.as_tensor(source_label_code, dtype=torch.long, device=device)
    target_code_t = torch.as_tensor(target_label_code, dtype=torch.long, device=device)

    rho_B_boot = np.full(BOOTSTRAP_REPLICATES, np.nan, dtype=np.float64)
    rho_C_boot = np.full(BOOTSTRAP_REPLICATES, np.nan, dtype=np.float64)
    completed = 0
    attempts = 0

    if paths["checkpoint"].exists():
        ck = np.load(paths["checkpoint"])
        completed = int(ck["completed"])
        attempts = int(ck["attempts"])
        rho_B_boot[:] = ck["rho_B_boot"]
        rho_C_boot[:] = ck["rho_C_boot"]
        print(f"      resuming composition bootstrap valid={completed}/{BOOTSTRAP_REPLICATES}; attempts={attempts}")

    engine_validated = completed > 0
    t0 = time.perf_counter()

    while completed < BOOTSTRAP_REPLICATES and attempts < MAX_ATTEMPTS:
        attempts += 1
        rng = composition_seed(target_name, program_id, attempts)
        src_counts, src_idx = stratified_resample_counts(source_labels, rng)
        tgt_counts, tgt_idx = stratified_resample_counts(target_labels, rng)

        src_counts_t = torch.as_tensor(src_counts, dtype=torch.float64, device=device)
        tgt_counts_t = torch.as_tensor(tgt_counts, dtype=torch.float64, device=device)

        B_src_t, _ = selected_edge_corr_weighted(source_t, src_counts_t, ii_t, jj_t)
        if B_src_t is None:
            del src_counts_t, tgt_counts_t
            continue
        B_tgt_t, _ = selected_edge_corr_weighted(target_t, tgt_counts_t, ii_t, jj_t)
        if B_tgt_t is None:
            del src_counts_t, tgt_counts_t, B_src_t
            continue

        C_src_t, _ = selected_edge_corr_residualized(source_t, src_counts_t, source_code_t, 5, ii_t, jj_t)
        if C_src_t is None:
            del src_counts_t, tgt_counts_t, B_src_t, B_tgt_t
            continue
        C_tgt_t, _ = selected_edge_corr_residualized(target_t, tgt_counts_t, target_code_t, 5, ii_t, jj_t)
        if C_tgt_t is None:
            del src_counts_t, tgt_counts_t, B_src_t, B_tgt_t, C_src_t
            continue

        if not engine_validated:
            k = min(CPU_VALIDATION_EDGES, len(ii))
            B_src_ref = explicit_edges(source_x, src_idx, ii, jj, k)
            B_tgt_ref = explicit_edges(target_x, tgt_idx, ii, jj, k)

            src_labels_draw = source_labels[src_idx]
            tgt_labels_draw = target_labels[tgt_idx]
            src_r = residualize_np(source_x[src_idx], src_labels_draw)
            tgt_r = residualize_np(target_x[tgt_idx], tgt_labels_draw)
            C_src_ref = explicit_edges(src_r, np.arange(len(src_r)), ii, jj, k)
            C_tgt_ref = explicit_edges(tgt_r, np.arange(len(tgt_r)), ii, jj, k)

            got = [
                B_src_t[:k].detach().cpu().numpy(),
                B_tgt_t[:k].detach().cpu().numpy(),
                C_src_t[:k].detach().cpu().numpy(),
                C_tgt_t[:k].detach().cpu().numpy(),
            ]
            refs = [B_src_ref, B_tgt_ref, C_src_ref, C_tgt_ref]
            diff = max(float(np.max(np.abs(a - b))) for a, b in zip(got, refs))
            if diff > CPU_GPU_TOL:
                raise RuntimeError(
                    f"{target_name} {program_id}: composition bootstrap CPU/GPU mismatch {diff:.3e}"
                )
            print(f"      composition selected-edge CPU↔GPU validation max|Δ|={diff:.3e}")
            engine_validated = True

        rho_B_boot[completed] = gpu_spearman(B_src_t, B_tgt_t)
        rho_C_boot[completed] = gpu_spearman(C_src_t, C_tgt_t)
        completed += 1

        del src_counts_t, tgt_counts_t, B_src_t, B_tgt_t, C_src_t, C_tgt_t

        if completed % CHECKPOINT_EVERY == 0 or completed == BOOTSTRAP_REPLICATES:
            np.savez_compressed(
                paths["checkpoint"],
                completed=np.array(completed, dtype=np.int64),
                attempts=np.array(attempts, dtype=np.int64),
                rho_B_boot=rho_B_boot,
                rho_C_boot=rho_C_boot,
            )
            print(
                f"      composition bootstrap {completed:4d}/{BOOTSTRAP_REPLICATES}; attempts={attempts:5d}; "
                f"B median={np.nanmedian(rho_B_boot[:completed]):+.4f}; "
                f"C median={np.nanmedian(rho_C_boot[:completed]):+.4f}; elapsed={time.perf_counter()-t0:.1f}s"
            )

    ci_ok = completed == BOOTSTRAP_REPLICATES
    delta = rho_C_boot[:completed] - rho_B_boot[:completed]

    result = {
        "script_version": SCRIPT_VERSION,
        "status": "COMPLETE" if ci_ok else "BOOTSTRAP_CI_NOT_ESTIMABLE",
        "target": target_name,
        "program_id": program_id,
        "source_complete_case_n": int(source_x.shape[0]),
        "target_complete_case_n": int(target_x.shape[0]),
        "n_genes": int(m),
        "n_all_edges": int(len(tri_i)),
        "n_bootstrap_edges": int(len(ii)),
        "rho_A_full_primary": rho_A,
        "rho_B_complete_case": rho_B,
        "rho_C_subtype_residualized": rho_C,
        "B_minus_A": rho_B - rho_A,
        "C_minus_B": rho_C - rho_B,
        "B_source_full_corr_cpu_gpu_max_abs_diff": B_source_diff,
        "B_target_full_corr_cpu_gpu_max_abs_diff": B_target_diff,
        "C_source_full_corr_cpu_gpu_max_abs_diff": C_source_diff,
        "C_target_full_corr_cpu_gpu_max_abs_diff": C_target_diff,
        "bootstrap_valid_replicates": int(completed),
        "bootstrap_attempts": int(attempts),
        "bootstrap_invalid_attempts": int(attempts - completed),
        "B_bootstrap_median": float(np.median(rho_B_boot[:completed])) if completed else None,
        "B_bootstrap_ci_low": float(np.quantile(rho_B_boot[:completed], 0.025)) if completed else None,
        "B_bootstrap_ci_high": float(np.quantile(rho_B_boot[:completed], 0.975)) if completed else None,
        "C_bootstrap_median": float(np.median(rho_C_boot[:completed])) if completed else None,
        "C_bootstrap_ci_low": float(np.quantile(rho_C_boot[:completed], 0.025)) if completed else None,
        "C_bootstrap_ci_high": float(np.quantile(rho_C_boot[:completed], 0.975)) if completed else None,
        "C_minus_B_bootstrap_median": float(np.median(delta)) if completed else None,
        "C_minus_B_bootstrap_ci_low": float(np.quantile(delta, 0.025)) if completed else None,
        "C_minus_B_bootstrap_ci_high": float(np.quantile(delta, 0.975)) if completed else None,
        "B_minus_A_bootstrap_ci": None,
        "B_minus_A_ci_reason": "A uses the separate full-cohort primary estimand and does not share the complete-case stratified bootstrap base",
        "primary_classification_changed": False,
    }
    paths["result"].write_text(json.dumps(result, indent=2), encoding="utf-8")

    del source_t, target_t, ii_t, jj_t, source_code_t, target_code_t
    torch.cuda.empty_cache()
    return result


def main() -> None:
    print("=" * 150)
    print("Paper 4 / TCBB - PAM50 subtype heterogeneity + subtype-composition residualization sensitivity")
    print("=" * 150)
    print(f"Script version: {SCRIPT_VERSION}")
    print()
    print("Frozen scientific contract:")
    print("  Target outcomes/treatment loaded:                     NO")
    print("  Subtype tiers:                                        LumA/LumB/Basal CORE; Her2 LIMITED_N; Normal not standalone")
    print("  Subtype point edge effect:                            ALL fixed-gene edges")
    print(f"  Subtype edge bootstrap:                               {BOOTSTRAP_REPLICATES:,} valid source+target resamples")
    print("  Subtype loading:                                      descriptive point estimate only")
    print("  Composition A/B/C:                                    full primary / PAM50-complete / subtype-residualized")
    print(f"  Composition paired B/C bootstrap:                     {BOOTSTRAP_REPLICATES:,} stratified valid replicates")
    print("  Bootstrap edge subset:                                exact corrected 05b v2 fixed subset")
    print("  Exact-constant fixed gene:                            cell NOT ESTIMABLE; no gene dropping")
    print("  Formal new p-values / BH:                             NO")
    print("  Primary classification changed:                       NO")
    print("=" * 150)

    for p in [
        EXEC_CONTRACT,
        PAM50_CONTRACT,
        TCGA_ALIGN,
        SCANB_ALIGN,
        MET_ALIGN,
        CLASSIFICATION,
        UNIVERSE,
        WEIGHTS,
        TCGA_EXPR,
        SCANB_EXPR,
        SCANB_PRIMARY,
        METABRIC_EXPR,
    ]:
        require(p)

    ex = json.loads(EXEC_CONTRACT.read_text(encoding="utf-8"))
    if ex.get("status") != "FROZEN_BEFORE_PAM50_SENSITIVITY_RESULTS":
        raise RuntimeError("05f2 execution contract has unexpected status.")

    pam = json.loads(PAM50_CONTRACT.read_text(encoding="utf-8"))
    if pam.get("status") != "FROZEN_BEFORE_FIRST_TARGET_PRESERVATION_STATISTIC":
        raise RuntimeError("04e PAM50 contract has unexpected status.")

    if not torch.cuda.is_available():
        raise RuntimeError("CUDA is unavailable.")
    device = torch.device("cuda:0")
    prop = torch.cuda.get_device_properties(0)
    print(f"CUDA device: {prop.name}; VRAM={prop.total_memory/1024**3:.2f} GB; capability={prop.major}.{prop.minor}")

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    universe = pd.read_csv(UNIVERSE, sep="\t", dtype=str).fillna("")
    universe["source_gene_rank_by_MAD"] = pd.to_numeric(universe["source_gene_rank_by_MAD"], errors="raise").astype(int)
    universe = universe.sort_values("source_gene_rank_by_MAD").reset_index(drop=True)
    if len(universe) != 10000:
        raise RuntimeError("Frozen source universe is not 10,000 genes.")
    universe_set = {canon_symbol(x) for x in universe["Hugo_Symbol"]}

    weights = pd.read_csv(WEIGHTS, sep="\t", dtype=str).fillna("")
    classification = pd.read_csv(CLASSIFICATION, sep="\t", low_memory=False)

    primary = pd.read_csv(SCANB_PRIMARY, sep="\t", dtype=str).fillna("")
    primary_titles = list(primary["primary_title"])

    print("\n[1/4] Loading transformed source/target matrices ...")
    source = load_source_raw(universe)
    scanb = load_scanb_raw(universe_set, primary_titles)
    met = load_metabric_raw(universe_set)
    print(f"  TCGA:     {source['x'].shape[0]:,} x {source['x'].shape[1]:,}")
    print(f"  SCAN-B:   {scanb['x'].shape[0]:,} x {scanb['x'].shape[1]:,}")
    print(f"  METABRIC: {met['x'].shape[0]:,} x {met['x'].shape[1]:,}")

    source_labels = load_labels(source["samples"], "TCGA_BRCA")
    scanb_labels = load_labels(scanb["samples"], "SCANB_GSE96058")
    met_labels = load_labels(met["samples"], "METABRIC")
    validate_label_counts(source_labels, scanb_labels, met_labels)
    print("  PAM50 alignment/count replay: PASS")

    subtype_results = []
    composition_results = []

    for target_name, target, target_labels in [
        ("SCANB_GSE96058", scanb, scanb_labels),
        ("METABRIC", met, met_labels),
    ]:
        print("\n" + "=" * 150)
        print(f"[2/4] SUBTYPE HETEROGENEITY — {target_name}")
        print("=" * 150)

        subcls = classification.loc[classification["target"] == target_name].copy()

        for i, row in enumerate(subcls.to_dict(orient="records"), start=1):
            program_id = str(row["program_id"])
            if not boolish(row["primary_assessable"]):
                print(f"\n  [{i:02d}/12] {program_id}: primary NOT ASSESSABLE — PAM50 sensitivity skipped")
                continue

            genes = program_genes(target_name, program_id)
            src_cols = np.array([source["gene_index"][g] for g in genes], dtype=np.int64)
            tgt_cols = np.array([target["gene_index"][g] for g in genes], dtype=np.int64)
            src_prog = source["x"][:, src_cols]
            tgt_prog = target["x"][:, tgt_cols]
            source_loading, source_sign = source_loading_signs(weights, program_id, genes)
            ii, jj = fixed_edge_subset(target_name, program_id)

            print(f"\n  [{i:02d}/12] {program_id}: genes={len(genes):,}; CLASS={row['primary_classification']}")

            for subtype in COMPUTED_SUBTYPES:
                sidx = np.where(source_labels == subtype)[0]
                tidx = np.where(target_labels == subtype)[0]
                sx = src_prog[sidx, :]
                tx = tgt_prog[tidx, :]

                print(f"    {subtype:5s} [{TIER[subtype]}]: source n={len(sidx):3d}; target n={len(tidx):4d}")
                result = run_subtype_cell(
                    target_name=target_name,
                    program_id=program_id,
                    subtype=subtype,
                    source_x=sx,
                    target_x=tx,
                    source_loading=source_loading,
                    source_sign=source_sign,
                    ii=ii,
                    jj=jj,
                    device=device,
                )
                subtype_results.append(result)

                if result["status"] == "COMPLETE":
                    print(
                        f"      edge rho={result['rho_edge_all_edges']:+.4f}; "
                        f"95% CI=[{result['bootstrap_ci_low']:+.4f},{result['bootstrap_ci_high']:+.4f}]; "
                        f"loading rho={result['rho_load_descriptive']:+.4f}"
                    )
                else:
                    print(f"      STATUS={result['status']}")

            subtype_results.append({
                "script_version": SCRIPT_VERSION,
                "status": "NOT_INDIVIDUALLY_ASSESSABLE_BY_FROZEN_TIER",
                "target": target_name,
                "program_id": program_id,
                "subtype": "Normal",
                "tier": "NOT_INDIVIDUALLY_ASSESSABLE",
                "source_n": int(np.sum(source_labels == "Normal")),
                "target_n": int(np.sum(target_labels == "Normal")),
                "primary_classification_changed": False,
            })

        print("\n" + "=" * 150)
        print(f"[3/4] PAM50 COMPOSITION A/B/C — {target_name}")
        print("=" * 150)

        src_cc_idx = np.where(np.isin(source_labels, CANONICAL_SUBTYPES))[0]
        if target_name == "SCANB_GSE96058":
            tgt_cc_idx = np.where(np.isin(target_labels, CANONICAL_SUBTYPES))[0]
        else:
            tgt_cc_idx = np.where(np.isin(target_labels, CANONICAL_SUBTYPES))[0]

        src_cc_labels = source_labels[src_cc_idx]
        tgt_cc_labels = target_labels[tgt_cc_idx]

        if len(src_cc_idx) != 834:
            raise RuntimeError(f"Source PAM50 complete-case n is {len(src_cc_idx)}, expected 834.")
        expected_target_n = 3273 if target_name == "SCANB_GSE96058" else 1756
        if len(tgt_cc_idx) != expected_target_n:
            raise RuntimeError(f"{target_name} PAM50 complete-case n is {len(tgt_cc_idx)}, expected {expected_target_n}.")

        for i, row in enumerate(subcls.to_dict(orient="records"), start=1):
            program_id = str(row["program_id"])
            if not boolish(row["primary_assessable"]):
                continue

            genes = program_genes(target_name, program_id)
            src_cols = np.array([source["gene_index"][g] for g in genes], dtype=np.int64)
            tgt_cols = np.array([target["gene_index"][g] for g in genes], dtype=np.int64)
            source_cc_x = source["x"][src_cc_idx][:, src_cols]
            target_cc_x = target["x"][tgt_cc_idx][:, tgt_cols]
            ii, jj = fixed_edge_subset(target_name, program_id)

            print(
                f"\n  [{i:02d}/12] {program_id}: A={float(row['rho_edge']):+.4f}; "
                f"source complete n={len(src_cc_idx)}; target complete n={len(tgt_cc_idx)}"
            )
            result = composition_point_and_bootstrap(
                target_name=target_name,
                program_id=program_id,
                row=row,
                source_x=source_cc_x,
                source_labels=src_cc_labels,
                target_x=target_cc_x,
                target_labels=tgt_cc_labels,
                ii=ii,
                jj=jj,
                device=device,
            )
            composition_results.append(result)

            if result["status"] == "COMPLETE":
                print(
                    f"      B={result['rho_B_complete_case']:+.4f}; "
                    f"C={result['rho_C_subtype_residualized']:+.4f}; "
                    f"B-A={result['B_minus_A']:+.4f}; C-B={result['C_minus_B']:+.4f}; "
                    f"C-B 95% CI=[{result['C_minus_B_bootstrap_ci_low']:+.4f},{result['C_minus_B_bootstrap_ci_high']:+.4f}]"
                )
            else:
                print(f"      STATUS={result['status']}")

        torch.cuda.empty_cache()
        gc.collect()

    print("\n[4/4] Writing PAM50 sensitivity summaries ...")
    subtype_df = pd.DataFrame(subtype_results)
    composition_df = pd.DataFrame(composition_results)

    subtype_path = OUT_DIR / "pam50_subtype_sensitivity_summary_v1.tsv"
    composition_path = OUT_DIR / "pam50_composition_residualization_summary_v1.tsv"
    subtype_df.to_csv(subtype_path, sep="\t", index=False)
    composition_df.to_csv(composition_path, sep="\t", index=False)

    subtype_complete = int((subtype_df["status"] == "COMPLETE").sum())
    subtype_not_estimable = int(subtype_df["status"].isin(["NOT_ESTIMABLE_EXACT_VARIATION", "BOOTSTRAP_CI_NOT_ESTIMABLE"]).sum())
    composition_complete = int((composition_df["status"] == "COMPLETE").sum())
    composition_not_estimable = int((composition_df["status"] != "COMPLETE").sum())

    master = {
        "script_version": SCRIPT_VERSION,
        "status": "PAM50_SUBTYPE_AND_COMPOSITION_SENSITIVITY_COMPLETE",
        "subtype_complete_cells": subtype_complete,
        "subtype_nonestimable_or_ci_nonestimable_cells": subtype_not_estimable,
        "composition_complete_programs": composition_complete,
        "composition_noncomplete_programs": composition_not_estimable,
        "normal_subtype_standalone_computed": False,
        "new_discovery_p_values": False,
        "primary_classification_changed": False,
        "target_outcomes_or_treatment_loaded": False,
        "subtype_summary_file": str(subtype_path),
        "composition_summary_file": str(composition_path),
    }
    master_path = OUT_DIR / "pam50_subtype_composition_sensitivity_v1.json"
    master_path.write_text(json.dumps(master, indent=2), encoding="utf-8")

    print("\n" + "=" * 150)
    print("05f PAM50 SUBTYPE / COMPOSITION SENSITIVITY: COMPLETE")
    print("=" * 150)
    print(f"Subtype COMPLETE cells:                     {subtype_complete}")
    print(f"Subtype non-estimable/CI-non-estimable:     {subtype_not_estimable}")
    print(f"Composition COMPLETE programs:              {composition_complete}")
    print(f"Composition non-complete programs:          {composition_not_estimable}")
    print("Normal standalone subtype conclusion:        NOT COMPUTED (frozen tier)")
    print("Primary classifications changed:             NO")
    print("New subtype discovery p-values:               NO")
    print()
    print(f"Subtype summary:    {subtype_path}")
    print(f"Composition summary:{composition_path}")
    print(f"Master JSON:        {master_path}")
    print("=" * 150)


if __name__ == "__main__":
    main()
