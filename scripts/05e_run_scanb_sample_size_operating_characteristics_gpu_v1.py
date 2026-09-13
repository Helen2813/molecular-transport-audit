from __future__ import annotations

import gc
import json
import math
import re
from pathlib import Path

import numpy as np
import pandas as pd

try:
    import torch
except Exception as exc:
    raise RuntimeError("PyTorch is required for 05e.") from exc

try:
    from scipy.stats import rankdata
except Exception as exc:
    raise RuntimeError("scipy is required for Spearman statistics.") from exc


SCRIPT_VERSION = "05e-run-scanb-sample-size-operating-characteristics-gpu-v1-no-cli"

DATA_ROOT = Path(r"D:\paper4_tcbb_data")

EXEC_CONTRACT = (
    DATA_ROOT
    / "paper4_tcbb_scanb_sample_size_execution_contract_v1"
    / "scanb_sample_size_execution_contract_v1.json"
)
CLASSIFICATION = (
    DATA_ROOT
    / "paper4_tcbb_final_mapping_specificity_null_v2"
    / "primary_pooled_final_classification_v2.tsv"
)
SPECIFICITY_ROOT = DATA_ROOT / "paper4_tcbb_final_mapping_specificity_null_v2"
DIRECT_ROOT = DATA_ROOT / "paper4_tcbb_primary_pooled_direct_preservation_v2"

WEIGHTS = (
    DATA_ROOT
    / "paper4_tcbb_frozen_source_programs_v1"
    / "tcga_frozen_source_program_weights_v1.tsv"
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

OUT_DIR = DATA_ROOT / "paper4_tcbb_scanb_sample_size_operating_characteristics_v1"

SAMPLE_SIZES = [43, 100, 250, 500, 1000, 2000, 3273]
REPEATS_NONFULL = 50
MAX_ATTEMPTS = 5000
BASE_SEED = 20260916

PC1_RESIDUAL_TOL = 1e-10
MAX_POWER_ITER = 400
REPLAY_TOL = 5e-10


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


def centered_rank(x: np.ndarray) -> tuple[np.ndarray, float]:
    r = rankdata(np.asarray(x, dtype=np.float64), method="average")
    r -= r.mean()
    ss = float(np.dot(r, r))
    if ss <= 0:
        raise RuntimeError("Degenerate rank vector.")
    return r, math.sqrt(ss)


def spearman_against_fixed(
    fixed_rank: np.ndarray,
    fixed_norm: float,
    x: np.ndarray,
) -> float:
    r = rankdata(np.asarray(x, dtype=np.float64), method="average")
    r -= r.mean()
    norm = math.sqrt(float(np.dot(r, r)))
    if norm <= 0:
        return float("nan")
    return float(np.dot(fixed_rank, r) / (fixed_norm * norm))


def load_scanb_corrected(
    universe_set: set[str],
    primary_titles: list[str],
) -> dict:
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

        fpkm = np.maximum(np.exp2(arr) - 0.1, 0.0)
        transformed = np.log2(fpkm + 1.0)

        tmp = pd.DataFrame(transformed, columns=primary_titles)
        tmp.insert(0, "Hugo_Symbol", symbols.loc[mask].to_numpy())
        parts.append(tmp)

    if not parts:
        raise RuntimeError("No SCAN-B genes found.")

    df = pd.concat(parts, ignore_index=True)
    df = df.groupby("Hugo_Symbol", sort=False, as_index=True).mean(numeric_only=True)

    x = df.to_numpy(dtype=np.float64).T
    finite = np.isfinite(x).all(axis=0)
    exact_min = np.full(x.shape[1], np.nan)
    exact_max = np.full(x.shape[1], np.nan)
    exact_min[finite] = np.min(x[:, finite], axis=0)
    exact_max[finite] = np.max(x[:, finite], axis=0)
    eligible = finite & (exact_max > exact_min)

    genes_all = list(df.index)
    genes = [g for g, ok in zip(genes_all, eligible) if ok]
    x = x[:, eligible]

    if len(genes) != 9220:
        raise RuntimeError(f"Expected 9,220 corrected SCAN-B genes, got {len(genes)}.")

    return {
        "genes": genes,
        "gene_index": {g: i for i, g in enumerate(genes)},
        "x": x,
    }


def exact_variable_columns(x: np.ndarray) -> np.ndarray:
    finite = np.isfinite(x).all(axis=0)
    out = np.zeros(x.shape[1], dtype=bool)
    if finite.any():
        xmin = np.min(x[:, finite], axis=0)
        xmax = np.max(x[:, finite], axis=0)
        out[finite] = xmax > xmin
    return out


def standardize_np(x: np.ndarray) -> np.ndarray:
    mu = x.mean(axis=0, keepdims=True)
    sd = x.std(axis=0, ddof=1, keepdims=True)
    if np.any(~np.isfinite(sd)) or np.any(sd <= 0):
        raise RuntimeError("Invalid SD reached standardization.")
    return (x - mu) / sd


def leading_pc1(
    corr_t: torch.Tensor,
    z_t: torch.Tensor,
    source_loading_t: torch.Tensor,
    sign_t: torch.Tensor,
    random_start_t: torch.Tensor,
) -> torch.Tensor:
    eps = torch.finfo(torch.float64).eps
    candidates = []

    for start in [source_loading_t, sign_t, random_start_t]:
        v = start / torch.linalg.vector_norm(start)
        residual_value = float("inf")

        for it in range(1, MAX_POWER_ITER + 1):
            w = corr_t @ v
            wn = torch.linalg.vector_norm(w)
            if not torch.isfinite(wn) or wn <= 0:
                break
            nv = w / wn
            if torch.dot(nv, v) < 0:
                nv = -nv
            v = nv

            if it % 5 == 0 or it == MAX_POWER_ITER:
                av = corr_t @ v
                lam = torch.dot(v, av)
                residual = torch.linalg.vector_norm(
                    av - lam * v
                ) / (torch.abs(lam) + eps)
                residual_value = float(residual.detach().cpu())
                if residual_value <= PC1_RESIDUAL_TOL:
                    break

        av = corr_t @ v
        lam = float(torch.dot(v, av).detach().cpu())
        candidates.append((lam, residual_value, v.detach().clone()))

    candidates.sort(key=lambda t: t[0], reverse=True)
    lam, residual, v = candidates[0]

    if residual > PC1_RESIDUAL_TOL:
        evals, evecs = torch.linalg.eigh(corr_t)
        v = evecs[:, -1]

    scores = z_t @ v
    sign_score = torch.mean(z_t * sign_t[None, :], dim=1)

    sc = scores - scores.mean()
    sg = sign_score - sign_score.mean()
    den = torch.linalg.vector_norm(sc) * torch.linalg.vector_norm(sg)
    if den <= 0:
        raise RuntimeError("PC1 orientation score is degenerate.")

    if torch.dot(sc, sg) / den < 0:
        v = -v

    return v


def rng_for_attempt(
    program_id: str,
    sample_size: int,
    attempt_id: int,
) -> np.random.Generator:
    ss = np.random.SeedSequence(
        [
            BASE_SEED,
            module_number(program_id),
            int(sample_size),
            int(attempt_id),
        ]
    )
    return np.random.default_rng(ss)


def random_start_for_program(
    program_id: str,
    m: int,
) -> np.ndarray:
    # Mirrors 05a target-index=1 seed schedule for the target-PC1 random start.
    seed = (
        BASE_SEED
        + 1 * 100_000
        + module_number(program_id) * 100
        + 31
    )
    rng = np.random.default_rng(seed)
    return rng.normal(size=m)


def checkpoint_path(program_id: str) -> Path:
    d = OUT_DIR / "per_program"
    d.mkdir(parents=True, exist_ok=True)
    return d / f"{program_id}_sample_size_replicates_v1.tsv"


def degenerate_path(program_id: str) -> Path:
    d = OUT_DIR / "per_program"
    d.mkdir(parents=True, exist_ok=True)
    return d / f"{program_id}_sample_size_degenerate_gene_counts_v1.tsv"


def run_program(
    target: dict,
    row: dict,
    weights_df: pd.DataFrame,
    device: torch.device,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    program_id = str(row["program_id"])

    gene_file = (
        DIRECT_ROOT
        / "per_program"
        / "SCANB_GSE96058"
        / f"{program_id}_evaluable_genes_v2.tsv"
    )
    require(gene_file)
    eg = pd.read_csv(gene_file, sep="\t", dtype=str).fillna("")
    genes = [canon_symbol(x) for x in eg["Hugo_Symbol"]]
    m = len(genes)

    cols = np.array(
        [target["gene_index"][g] for g in genes],
        dtype=np.int64,
    )
    x = target["x"][:, cols]

    if not exact_variable_columns(x).all():
        raise RuntimeError(
            f"{program_id}: corrected full-size program still contains a constant gene."
        )

    wdf = weights_df.loc[
        weights_df["program_id"] == program_id
    ].copy()
    wdf["Hugo_Symbol"] = wdf["Hugo_Symbol"].map(canon_symbol)
    wdf = wdf.set_index("Hugo_Symbol", drop=False)

    source_loading = np.array(
        [float(wdf.loc[g, "source_pc1_loading"]) for g in genes],
        dtype=np.float64,
    )
    source_sign = np.array(
        [float(wdf.loc[g, "source_pc1_loading_sign"]) for g in genes],
        dtype=np.float64,
    )
    source_loading_rank, source_loading_norm = centered_rank(source_loading)

    spec_npz = (
        SPECIFICITY_ROOT
        / "per_program"
        / "SCANB_GSE96058"
        / f"{program_id}_mapping_specificity_null_v2.npz"
    )
    require(spec_npz)
    nz = np.load(spec_npz)

    source_edge_subset = np.asarray(
        nz["source_edge_subset"],
        dtype=np.float64,
    )
    ii = np.asarray(
        nz["edge_slot_i_0based"],
        dtype=np.int64,
    )
    jj = np.asarray(
        nz["edge_slot_j_0based"],
        dtype=np.int64,
    )
    source_edge_rank, source_edge_norm = centered_rank(source_edge_subset)

    if int(max(ii.max(), jj.max())) >= m:
        raise RuntimeError(f"{program_id}: frozen specificity pair index exceeds gene set.")

    source_loading_t = torch.as_tensor(
        source_loading,
        dtype=torch.float64,
        device=device,
    )
    sign_t = torch.as_tensor(
        source_sign,
        dtype=torch.float64,
        device=device,
    )
    random_start_t = torch.as_tensor(
        random_start_for_program(program_id, m),
        dtype=torch.float64,
        device=device,
    )
    ii_t = torch.as_tensor(ii, dtype=torch.long, device=device)
    jj_t = torch.as_tensor(jj, dtype=torch.long, device=device)

    # Full-size replay once.
    z_full = standardize_np(x)
    zf_t = torch.as_tensor(z_full, dtype=torch.float64, device=device)
    corr_full_t = (zf_t.T @ zf_t) / (z_full.shape[0] - 1)

    edge_full = corr_full_t[ii_t, jj_t].detach().cpu().numpy()
    rho_edge_full = spearman_against_fixed(
        source_edge_rank,
        source_edge_norm,
        edge_full,
    )

    load_full_t = leading_pc1(
        corr_full_t,
        zf_t,
        source_loading_t,
        sign_t,
        random_start_t,
    )
    load_full = load_full_t.detach().cpu().numpy()
    rho_load_full = spearman_against_fixed(
        source_loading_rank,
        source_loading_norm,
        load_full,
    )

    expected_edge = float(row["mapping_specificity_rho"])
    expected_load = float(row["rho_load"])

    if abs(rho_edge_full - expected_edge) > REPLAY_TOL:
        raise RuntimeError(
            f"{program_id}: full-size edge replay mismatch {rho_edge_full} vs {expected_edge}"
        )
    if abs(rho_load_full - expected_load) > REPLAY_TOL:
        raise RuntimeError(
            f"{program_id}: full-size loading replay mismatch {rho_load_full} vs {expected_load}"
        )

    del zf_t, corr_full_t, load_full_t
    torch.cuda.empty_cache()

    rows = [
        {
            "program_id": program_id,
            "sample_size": 3273,
            "repeat_id": 1,
            "attempt_id": 1,
            "rho_edge": rho_edge_full,
            "rho_load": rho_load_full,
            "valid_repeat": 1,
        }
    ]

    degenerate_records = []
    checkpoint = checkpoint_path(program_id)

    for n in SAMPLE_SIZES[:-1]:
        valid = 0
        attempts = 0
        edge_vals = []
        load_vals = []
        deg_counts = np.zeros(m, dtype=np.int64)

        while valid < REPEATS_NONFULL and attempts < MAX_ATTEMPTS:
            attempts += 1
            rng = rng_for_attempt(
                program_id,
                n,
                attempts,
            )
            idx = np.sort(
                rng.choice(
                    x.shape[0],
                    size=n,
                    replace=False,
                )
            )

            xs = x[idx, :]
            variable = exact_variable_columns(xs)
            if not variable.all():
                deg_counts[~variable] += 1
                continue

            z = standardize_np(xs)
            zt = torch.as_tensor(
                z,
                dtype=torch.float64,
                device=device,
            )
            corr_t = (zt.T @ zt) / (n - 1)

            edges = corr_t[ii_t, jj_t].detach().cpu().numpy()
            rho_edge = spearman_against_fixed(
                source_edge_rank,
                source_edge_norm,
                edges,
            )

            load_t = leading_pc1(
                corr_t,
                zt,
                source_loading_t,
                sign_t,
                random_start_t,
            )
            loading = load_t.detach().cpu().numpy()
            rho_load = spearman_against_fixed(
                source_loading_rank,
                source_loading_norm,
                loading,
            )

            valid += 1
            edge_vals.append(rho_edge)
            load_vals.append(rho_load)

            rows.append(
                {
                    "program_id": program_id,
                    "sample_size": n,
                    "repeat_id": valid,
                    "attempt_id": attempts,
                    "rho_edge": rho_edge,
                    "rho_load": rho_load,
                    "valid_repeat": 1,
                }
            )

            del zt, corr_t, load_t

        for gene, count in zip(genes, deg_counts):
            if count:
                degenerate_records.append(
                    {
                        "program_id": program_id,
                        "sample_size": n,
                        "Hugo_Symbol": gene,
                        "degenerate_attempt_count": int(count),
                    }
                )

        status = "ESTIMABLE" if valid == REPEATS_NONFULL else "NOT_ESTIMABLE"

        rows.append(
            {
                "program_id": program_id,
                "sample_size": n,
                "repeat_id": 0,
                "attempt_id": attempts,
                "rho_edge": np.nan,
                "rho_load": np.nan,
                "valid_repeat": 0,
                "cell_status": status,
                "valid_repeats_obtained": valid,
                "invalid_attempts": attempts - valid,
                "valid_attempt_fraction": (
                    valid / attempts if attempts else np.nan
                ),
            }
        )

        print(
            f"      n={n:4d}: valid={valid:2d}/{REPEATS_NONFULL}; "
            f"attempts={attempts:4d}; invalid={attempts-valid:4d}; "
            + (
                f"edge median={np.median(edge_vals):+.4f}; "
                f"load median={np.median(load_vals):+.4f}"
                if valid
                else "NO VALID EFFECT ESTIMATE"
            )
        )

        pd.DataFrame(rows).to_csv(
            checkpoint,
            sep="\t",
            index=False,
        )

    del (
        source_loading_t,
        sign_t,
        random_start_t,
        ii_t,
        jj_t,
    )
    torch.cuda.empty_cache()

    rep_df = pd.DataFrame(rows)
    deg_df = pd.DataFrame(
        degenerate_records,
        columns=[
            "program_id",
            "sample_size",
            "Hugo_Symbol",
            "degenerate_attempt_count",
        ],
    )
    deg_df.to_csv(
        degenerate_path(program_id),
        sep="\t",
        index=False,
    )

    return rep_df, deg_df


def summarize_replicates(reps: pd.DataFrame) -> pd.DataFrame:
    out = []

    for program_id in sorted(reps["program_id"].unique()):
        p = reps.loc[reps["program_id"] == program_id]

        for n in SAMPLE_SIZES:
            cell = p.loc[
                (p["sample_size"] == n)
                & (p["valid_repeat"] == 1)
            ]

            meta = p.loc[
                (p["sample_size"] == n)
                & (p["valid_repeat"] == 0)
            ]

            if n == 3273:
                attempts = 1
                invalid = 0
                status = "FULL_SIZE_PRIMARY_ONCE"
            elif len(meta) == 1:
                attempts = int(meta.iloc[0]["attempt_id"])
                invalid = int(meta.iloc[0]["invalid_attempts"])
                status = str(meta.iloc[0]["cell_status"])
            else:
                attempts = np.nan
                invalid = np.nan
                status = "UNKNOWN"

            if len(cell):
                edge = cell["rho_edge"].to_numpy(dtype=float)
                load = cell["rho_load"].to_numpy(dtype=float)

                out.append(
                    {
                        "program_id": program_id,
                        "sample_size": n,
                        "status": status,
                        "valid_repeats": int(len(cell)),
                        "attempts": attempts,
                        "invalid_attempts": invalid,
                        "valid_attempt_fraction": (
                            len(cell) / attempts
                            if isinstance(attempts, (int, np.integer)) and attempts > 0
                            else 1.0 if n == 3273 else np.nan
                        ),
                        "edge_median": float(np.median(edge)),
                        "edge_q025": float(np.quantile(edge, 0.025)),
                        "edge_q975": float(np.quantile(edge, 0.975)),
                        "loading_median": float(np.median(load)),
                        "loading_q025": float(np.quantile(load, 0.025)),
                        "loading_q975": float(np.quantile(load, 0.975)),
                    }
                )
            else:
                out.append(
                    {
                        "program_id": program_id,
                        "sample_size": n,
                        "status": status,
                        "valid_repeats": 0,
                        "attempts": attempts,
                        "invalid_attempts": invalid,
                        "valid_attempt_fraction": 0.0,
                        "edge_median": np.nan,
                        "edge_q025": np.nan,
                        "edge_q975": np.nan,
                        "loading_median": np.nan,
                        "loading_q025": np.nan,
                        "loading_q975": np.nan,
                    }
                )

    return pd.DataFrame(out)


def main() -> None:
    print("=" * 146)
    print("Paper 4 / TCBB - SCAN-B sample-size operating characteristics")
    print("=" * 146)
    print(f"Script version: {SCRIPT_VERSION}")
    print()
    print("Frozen scientific contract:")
    print(f"  Sample sizes:                                       {SAMPLE_SIZES}")
    print(f"  Valid repeats/nonfull size:                         {REPEATS_NONFULL}")
    print("  Sampling:                                           without replacement")
    print("  Gene set:                                           fixed corrected 05a v2 program genes")
    print("  Degenerate small-n subset:                          discard subset + redraw")
    print(f"  Max attempts/nonfull cell:                          {MAX_ATTEMPTS}")
    print("  Edge statistic:                                     corrected 05b v2 fixed subset")
    print("  Loading statistic:                                  recomputed target PC1 -> rho_load")
    print("  n=3273:                                             primary result replay once")
    print("  Target outcomes/treatment loaded:                   NO")
    print("  Primary classifications changed:                    NO")
    print("=" * 146)

    for p in [
        EXEC_CONTRACT,
        CLASSIFICATION,
        WEIGHTS,
        SCANB_EXPR,
        SCANB_PRIMARY,
    ]:
        require(p)

    ex = json.loads(EXEC_CONTRACT.read_text(encoding="utf-8"))
    if ex.get("status") != "FROZEN_BEFORE_SCANB_SAMPLE_SIZE_RESULTS":
        raise RuntimeError("05e1 execution contract has unexpected status.")

    if not torch.cuda.is_available():
        raise RuntimeError("CUDA is unavailable.")
    device = torch.device("cuda:0")
    prop = torch.cuda.get_device_properties(0)
    print(
        f"CUDA device: {prop.name}; VRAM={prop.total_memory/1024**3:.2f} GB; "
        f"capability={prop.major}.{prop.minor}"
    )

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    classification = pd.read_csv(
        CLASSIFICATION,
        sep="\t",
        low_memory=False,
    )
    scanb_cls = classification.loc[
        classification["target"] == "SCANB_GSE96058"
    ].copy()

    weights_df = pd.read_csv(
        WEIGHTS,
        sep="\t",
        dtype=str,
    ).fillna("")
    weights_df["Hugo_Symbol"] = weights_df["Hugo_Symbol"].map(canon_symbol)

    primary = pd.read_csv(
        SCANB_PRIMARY,
        sep="\t",
        dtype=str,
    ).fillna("")
    primary_titles = list(primary["primary_title"])
    if len(primary_titles) != 3273:
        raise RuntimeError("SCAN-B primary manifest is not 3,273 samples.")

    # The full corrected universe is recoverable from the corrected 05a v2 gene files.
    universe = set()
    for program_id in sorted(scanb_cls["program_id"].unique()):
        if not boolish(
            scanb_cls.loc[
                scanb_cls["program_id"] == program_id,
                "primary_assessable",
            ].iloc[0]
        ):
            continue
        gf = (
            DIRECT_ROOT
            / "per_program"
            / "SCANB_GSE96058"
            / f"{program_id}_evaluable_genes_v2.tsv"
        )
        require(gf)
        gdf = pd.read_csv(gf, sep="\t", dtype=str).fillna("")
        universe.update(canon_symbol(x) for x in gdf["Hugo_Symbol"])

    # Add all 9,220 corrected universe genes through the exact-variation manifest.
    corrected_manifest = (
        DATA_ROOT
        / "paper4_tcbb_exact_variation_evaluability_correction_v1"
        / "scanb_target_exact_variation_evaluability_v1.tsv"
    )
    require(corrected_manifest)
    cm = pd.read_csv(corrected_manifest, sep="\t", dtype=str).fillna("")
    universe = {
        canon_symbol(x)
        for x in cm.loc[
            cm["statistically_evaluable_corrected"].astype(str) == "1",
            "Hugo_Symbol",
        ]
    }
    if len(universe) != 9220:
        raise RuntimeError(f"Expected corrected universe size 9,220, got {len(universe)}.")

    print("\n[1/2] Loading corrected SCAN-B matrix ...")
    target = load_scanb_corrected(
        universe,
        primary_titles,
    )
    print(
        f"  {target['x'].shape[0]:,} samples x "
        f"{target['x'].shape[1]:,} corrected-evaluable genes"
    )

    all_reps = []
    all_deg = []

    for i, row in enumerate(
        scanb_cls.to_dict(orient="records"),
        start=1,
    ):
        program_id = str(row["program_id"])

        if not boolish(row["primary_assessable"]):
            print(f"\n  [{i:02d}/12] {program_id}: NOT ASSESSABLE — skipped")
            continue

        print(
            f"\n  [{i:02d}/12] {program_id}: "
            f"CLASS={row['primary_classification']}; "
            f"full edge={float(row['mapping_specificity_rho']):+.4f}; "
            f"full load={float(row['rho_load']):+.4f}"
        )

        rep_df, deg_df = run_program(
            target=target,
            row=row,
            weights_df=weights_df,
            device=device,
        )
        all_reps.append(rep_df)
        if len(deg_df):
            all_deg.append(deg_df)

    print("\n[2/2] Writing sample-size operating-characteristic outputs ...")

    reps = pd.concat(all_reps, ignore_index=True)
    summary = summarize_replicates(reps)

    deg = (
        pd.concat(all_deg, ignore_index=True)
        if all_deg
        else pd.DataFrame(
            columns=[
                "program_id",
                "sample_size",
                "Hugo_Symbol",
                "degenerate_attempt_count",
            ]
        )
    )

    reps_path = OUT_DIR / "scanb_sample_size_replicates_v1.tsv"
    summary_path = OUT_DIR / "scanb_sample_size_summary_v1.tsv"
    deg_path = OUT_DIR / "scanb_sample_size_degenerate_gene_counts_v1.tsv"

    reps.to_csv(reps_path, sep="\t", index=False)
    summary.to_csv(summary_path, sep="\t", index=False)
    deg.to_csv(deg_path, sep="\t", index=False)

    nonfull = summary.loc[summary["sample_size"] != 3273]
    estimable = int((nonfull["status"] == "ESTIMABLE").sum())
    total_nonfull = len(nonfull)

    master = {
        "script_version": SCRIPT_VERSION,
        "status": "SCANB_SAMPLE_SIZE_OPERATING_CHARACTERISTICS_COMPLETE",
        "sample_sizes": SAMPLE_SIZES,
        "valid_repeats_per_nonfull_size": REPEATS_NONFULL,
        "maximum_attempts_per_nonfull_cell": MAX_ATTEMPTS,
        "nonfull_program_size_cells": int(total_nonfull),
        "estimable_nonfull_cells": estimable,
        "not_estimable_nonfull_cells": int(total_nonfull - estimable),
        "primary_classification_changed": False,
        "target_outcomes_or_treatment_loaded": False,
        "replicate_file": str(reps_path),
        "summary_file": str(summary_path),
        "degenerate_gene_file": str(deg_path),
    }

    master_path = OUT_DIR / "scanb_sample_size_operating_characteristics_v1.json"
    master_path.write_text(
        json.dumps(master, indent=2),
        encoding="utf-8",
    )

    print("\n" + "=" * 146)
    print("05e SCAN-B SAMPLE-SIZE OPERATING CHARACTERISTICS: COMPLETE")
    print("=" * 146)
    print(f"Estimable nonfull program×size cells:     {estimable}/{total_nonfull}")
    print(f"Not-estimable cells under frozen cap:     {total_nonfull-estimable}/{total_nonfull}")
    print("Primary classifications changed:          NO")
    print()
    print(f"Replicates: {reps_path}")
    print(f"Summary:    {summary_path}")
    print(f"Degenerate: {deg_path}")
    print(f"Master:     {master_path}")
    print("=" * 146)


if __name__ == "__main__":
    main()
