from __future__ import annotations

import gc
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

try:
    import torch
except Exception as exc:
    raise RuntimeError("PyTorch is required for 05g.") from exc

try:
    from scipy.stats import rankdata
except Exception as exc:
    raise RuntimeError("scipy is required for Spearman concordance.") from exc


SCRIPT_VERSION = "05g-run-scanb-technical-repeatability-gpu-v1-no-cli"

DATA_ROOT = Path(r"D:\paper4_tcbb_data")

EXEC_CONTRACT = (
    DATA_ROOT
    / "paper4_tcbb_scanb_technical_repeatability_execution_contract_v1"
    / "scanb_technical_repeatability_execution_contract_v1.json"
)
PAIR_MANIFEST = (
    DATA_ROOT
    / "paper4_tcbb_scanb_pairing_contract_v1"
    / "scanb_technical_replicate_pairs_v1.tsv"
)
CLASSIFICATION = (
    DATA_ROOT
    / "paper4_tcbb_final_mapping_specificity_null_v2"
    / "primary_pooled_final_classification_v2.tsv"
)
SPECIFICITY_ROOT = DATA_ROOT / "paper4_tcbb_final_mapping_specificity_null_v2"
DIRECT_ROOT = DATA_ROOT / "paper4_tcbb_primary_pooled_direct_preservation_v2"
SAMPLE_SIZE_REPLICATES = (
    DATA_ROOT
    / "paper4_tcbb_scanb_sample_size_operating_characteristics_v1"
    / "scanb_sample_size_replicates_v1.tsv"
)
SCANB_EXPR = (
    DATA_ROOT
    / "SCANB_GSE96058"
    / "GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv.gz"
)
CORRECTED_MANIFEST = (
    DATA_ROOT
    / "paper4_tcbb_exact_variation_evaluability_correction_v1"
    / "scanb_target_exact_variation_evaluability_v1.tsv"
)

OUT_DIR = DATA_ROOT / "paper4_tcbb_scanb_technical_repeatability_v1"

SAME_PLATFORM_PAIRS = 100
CROSS_PLATFORM_PAIRS = 36
CPU_GPU_TOL = 5e-10


def require(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Required file does not exist: {path}")


def canon_symbol(x: object) -> str:
    if x is None:
        return ""
    return " ".join(str(x).strip().split()).upper()


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


def load_scanb_selected_profiles(
    universe: set[str],
    titles: list[str],
) -> dict:
    parts = []

    for chunk in pd.read_csv(
        SCANB_EXPR,
        compression="gzip",
        chunksize=1500,
        low_memory=False,
    ):
        symbol_col = chunk.columns[0]
        missing_titles = [t for t in titles if t not in chunk.columns]
        if missing_titles:
            raise RuntimeError(
                f"SCAN-B expression file is missing selected titles, e.g. {missing_titles[:5]}"
            )

        symbols = chunk[symbol_col].map(canon_symbol)
        mask = symbols.isin(universe)
        if not mask.any():
            continue

        sub = chunk.loc[mask, [symbol_col] + titles]
        arr = sub[titles].to_numpy(dtype=np.float64)

        # Frozen SCAN-B representation:
        # published y=log2(FPKM+0.1), back-transform, then log2(FPKM+1).
        fpkm = np.maximum(np.exp2(arr) - 0.1, 0.0)
        transformed = np.log2(fpkm + 1.0)

        tmp = pd.DataFrame(transformed, columns=titles)
        tmp.insert(0, "Hugo_Symbol", symbols.loc[mask].to_numpy())
        parts.append(tmp)

    if not parts:
        raise RuntimeError("No corrected-universe SCAN-B genes were loaded.")

    df = pd.concat(parts, ignore_index=True)
    df = df.groupby("Hugo_Symbol", sort=False, as_index=True).mean(numeric_only=True)

    if set(df.index) != universe:
        missing = sorted(universe - set(df.index))
        raise RuntimeError(
            f"Corrected SCAN-B universe mismatch; missing genes: {missing[:20]}"
        )

    return {
        "genes": list(df.index),
        "gene_index": {g: i for i, g in enumerate(df.index)},
        "x": df.to_numpy(dtype=np.float64).T,
        "titles": titles,
        "title_index": {t: i for i, t in enumerate(titles)},
    }


def selected_edge_vector_gpu(
    x: np.ndarray,
    ii: np.ndarray,
    jj: np.ndarray,
    device: torch.device,
) -> tuple[np.ndarray, float]:
    z = standardize_np(x)
    zt = torch.as_tensor(z, dtype=torch.float64, device=device)
    ii_t = torch.as_tensor(ii, dtype=torch.long, device=device)
    jj_t = torch.as_tensor(jj, dtype=torch.long, device=device)

    # Avoid materializing full m x m correlation matrices: compute only selected edges.
    edge = torch.sum(
        zt[:, ii_t] * zt[:, jj_t],
        dim=0,
    ) / (z.shape[0] - 1)
    edge_np = edge.detach().cpu().numpy()

    k = min(200, len(ii))
    cpu = np.sum(
        z[:, ii[:k]] * z[:, jj[:k]],
        axis=0,
    ) / (z.shape[0] - 1)
    diff = float(np.max(np.abs(cpu - edge_np[:k]))) if k else 0.0

    if diff > CPU_GPU_TOL:
        raise RuntimeError(
            f"CPU/GPU selected-edge validation failed: {diff:.3e}"
        )

    del zt, ii_t, jj_t, edge
    torch.cuda.empty_cache()
    return edge_np, diff


def technical_result(
    program_id: str,
    pair_type: str,
    pair_df: pd.DataFrame,
    expr: dict,
    gene_index: dict[str, int],
    ii: np.ndarray,
    jj: np.ndarray,
    device: torch.device,
) -> tuple[dict, pd.DataFrame]:
    primary_titles = pair_df["primary_title"].tolist()
    replicate_titles = pair_df["replicate_title"].tolist()

    p_rows = np.array(
        [expr["title_index"][t] for t in primary_titles],
        dtype=np.int64,
    )
    r_rows = np.array(
        [expr["title_index"][t] for t in replicate_titles],
        dtype=np.int64,
    )

    # Program-local expression matrices are passed through global gene indices.
    cols = np.array(
        [gene_index[g] for g in pair_df.attrs["program_genes"]],
        dtype=np.int64,
    )

    xp = expr["x"][p_rows[:, None], cols[None, :]]
    xr = expr["x"][r_rows[:, None], cols[None, :]]

    p_var = exact_variable_columns(xp)
    r_var = exact_variable_columns(xr)
    jointly_valid = p_var & r_var

    deg_rows = []
    genes = pair_df.attrs["program_genes"]

    for gene, pv, rv in zip(genes, p_var, r_var):
        if not (pv and rv):
            deg_rows.append(
                {
                    "program_id": program_id,
                    "pair_type": pair_type,
                    "Hugo_Symbol": gene,
                    "primary_exact_variable": int(bool(pv)),
                    "replicate_exact_variable": int(bool(rv)),
                }
            )

    if not jointly_valid.all():
        return (
            {
                "program_id": program_id,
                "pair_type": pair_type,
                "n_pairs": int(len(pair_df)),
                "status": "NOT_ESTIMABLE_EXACT_VARIATION",
                "n_fixed_genes": int(len(genes)),
                "n_degenerate_genes": int(np.sum(~jointly_valid)),
                "rho_edge_technical": np.nan,
                "primary_cpu_gpu_selected_edge_max_abs_diff": np.nan,
                "replicate_cpu_gpu_selected_edge_max_abs_diff": np.nan,
            },
            pd.DataFrame(deg_rows),
        )

    p_edges, p_diff = selected_edge_vector_gpu(
        xp,
        ii,
        jj,
        device,
    )
    r_edges, r_diff = selected_edge_vector_gpu(
        xr,
        ii,
        jj,
        device,
    )

    p_rank, p_norm = centered_rank(p_edges)
    rho = spearman_against_fixed(
        p_rank,
        p_norm,
        r_edges,
    )

    return (
        {
            "program_id": program_id,
            "pair_type": pair_type,
            "n_pairs": int(len(pair_df)),
            "status": "COMPLETE",
            "n_fixed_genes": int(len(genes)),
            "n_degenerate_genes": 0,
            "rho_edge_technical": rho,
            "primary_cpu_gpu_selected_edge_max_abs_diff": p_diff,
            "replicate_cpu_gpu_selected_edge_max_abs_diff": r_diff,
        },
        pd.DataFrame(deg_rows),
    )


def main() -> None:
    print("=" * 146)
    print("Paper 4 / TCBB - SCAN-B technical repeatability / matched-n measurement ceiling")
    print("=" * 146)
    print(f"Script version: {SCRIPT_VERSION}")
    print()
    print("Frozen scientific contract:")
    print(f"  Same-platform technical pairs:                    {SAME_PLATFORM_PAIRS}")
    print(f"  Cross-platform technical pairs:                   {CROSS_PLATFORM_PAIRS}")
    print("  Pair selection:                                    exact frozen 01c manifest")
    print("  Gene set:                                          exact corrected 05a v2 genes")
    print("  Technical edge subset:                             exact corrected 05b v2 subset")
    print("  Technical endpoint:                                primary-vs-replicate edge concordance")
    print("  Same-platform matched-n reference:                 exact 05e n=100 50-repeat transport distribution")
    print("  Technical-set-specific gene dropping:              NO")
    print("  Technical pair bootstrap / new p-values:           NO")
    print("  Primary classifications changed:                   NO")
    print("=" * 146)

    for p in [
        EXEC_CONTRACT,
        PAIR_MANIFEST,
        CLASSIFICATION,
        SAMPLE_SIZE_REPLICATES,
        SCANB_EXPR,
        CORRECTED_MANIFEST,
    ]:
        require(p)

    ex = json.loads(EXEC_CONTRACT.read_text(encoding="utf-8"))
    if ex.get("status") != "FROZEN_BEFORE_SCANB_TECHNICAL_REPEATABILITY_RESULTS":
        raise RuntimeError("05g1 execution contract has unexpected status.")

    if not torch.cuda.is_available():
        raise RuntimeError("CUDA is unavailable.")
    device = torch.device("cuda:0")
    prop = torch.cuda.get_device_properties(0)
    print(
        f"CUDA device: {prop.name}; VRAM={prop.total_memory/1024**3:.2f} GB; "
        f"capability={prop.major}.{prop.minor}"
    )

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    pairs = pd.read_csv(
        PAIR_MANIFEST,
        sep="\t",
        dtype=str,
    ).fillna("")

    same = pairs.loc[pairs["same_platform"].astype(str) == "1"].copy()
    cross = pairs.loc[pairs["cross_platform"].astype(str) == "1"].copy()

    if len(same) != SAME_PLATFORM_PAIRS:
        raise RuntimeError(f"Expected 100 same-platform pairs, got {len(same)}.")
    if len(cross) != CROSS_PLATFORM_PAIRS:
        raise RuntimeError(f"Expected 36 cross-platform pairs, got {len(cross)}.")

    if same["pair_id"].duplicated().any() or cross["pair_id"].duplicated().any():
        raise RuntimeError("Duplicate pair_id in frozen technical-pair manifest.")

    all_titles = list(
        dict.fromkeys(
            same["primary_title"].tolist()
            + same["replicate_title"].tolist()
            + cross["primary_title"].tolist()
            + cross["replicate_title"].tolist()
        )
    )

    corrected = pd.read_csv(
        CORRECTED_MANIFEST,
        sep="\t",
        dtype=str,
    ).fillna("")
    universe = {
        canon_symbol(x)
        for x in corrected.loc[
            corrected["statistically_evaluable_corrected"].astype(str) == "1",
            "Hugo_Symbol",
        ]
    }
    if len(universe) != 9220:
        raise RuntimeError(f"Expected 9,220 corrected SCAN-B genes, got {len(universe)}.")

    print("\n[1/3] Loading all frozen technical-pair SCAN-B profiles ...")
    expr = load_scanb_selected_profiles(
        universe,
        all_titles,
    )
    print(
        f"  {len(all_titles):,} unique selected profiles x "
        f"{len(expr['genes']):,} corrected-evaluable genes"
    )

    classification = pd.read_csv(
        CLASSIFICATION,
        sep="\t",
        low_memory=False,
    )
    scanb_cls = classification.loc[
        classification["target"] == "SCANB_GSE96058"
    ].copy()

    ss = pd.read_csv(
        SAMPLE_SIZE_REPLICATES,
        sep="\t",
        low_memory=False,
    )

    technical_rows = []
    degenerate_frames = []

    print("\n[2/3] Technical repeatability by corrected program ...")

    assessable_rows = [
        r for r in scanb_cls.to_dict(orient="records")
        if boolish(r["primary_assessable"])
    ]

    for i, row in enumerate(assessable_rows, start=1):
        program_id = str(row["program_id"])

        gf = (
            DIRECT_ROOT
            / "per_program"
            / "SCANB_GSE96058"
            / f"{program_id}_evaluable_genes_v2.tsv"
        )
        require(gf)
        gdf = pd.read_csv(gf, sep="\t", dtype=str).fillna("")
        program_genes = [canon_symbol(x) for x in gdf["Hugo_Symbol"]]

        npz_path = (
            SPECIFICITY_ROOT
            / "per_program"
            / "SCANB_GSE96058"
            / f"{program_id}_mapping_specificity_null_v2.npz"
        )
        require(npz_path)
        nz = np.load(npz_path)
        ii = np.asarray(nz["edge_slot_i_0based"], dtype=np.int64)
        jj = np.asarray(nz["edge_slot_j_0based"], dtype=np.int64)

        if int(max(ii.max(), jj.max())) >= len(program_genes):
            raise RuntimeError(
                f"{program_id}: frozen edge subset exceeds corrected program gene set."
            )

        print(
            f"\n  [{i:02d}/{len(assessable_rows):02d}] {program_id}: "
            f"genes={len(program_genes):,}; edges={len(ii):,}"
        )

        same_local = same.copy()
        same_local.attrs["program_genes"] = program_genes
        same_result, same_deg = technical_result(
            program_id=program_id,
            pair_type="SAME_PLATFORM",
            pair_df=same_local,
            expr=expr,
            gene_index=expr["gene_index"],
            ii=ii,
            jj=jj,
            device=device,
        )

        cross_local = cross.copy()
        cross_local.attrs["program_genes"] = program_genes
        cross_result, cross_deg = technical_result(
            program_id=program_id,
            pair_type="CROSS_PLATFORM",
            pair_df=cross_local,
            expr=expr,
            gene_index=expr["gene_index"],
            ii=ii,
            jj=jj,
            device=device,
        )

        # Matched-n reference comes from exact 05e repeat-level results.
        transport = ss.loc[
            (ss["program_id"] == program_id)
            & (pd.to_numeric(ss["sample_size"], errors="coerce") == 100)
            & (pd.to_numeric(ss["valid_repeat"], errors="coerce") == 1),
            "rho_edge",
        ].astype(float).to_numpy()

        if len(transport) != 50:
            raise RuntimeError(
                f"{program_id}: expected 50 valid 05e n=100 transport repeats, got {len(transport)}."
            )

        t_med = float(np.median(transport))
        t_q025 = float(np.quantile(transport, 0.025))
        t_q975 = float(np.quantile(transport, 0.975))

        if same_result["status"] == "COMPLETE":
            tech = float(same_result["rho_edge_technical"])
            same_result.update(
                {
                    "matched_n_transport_repeats": 50,
                    "matched_n_transport_median": t_med,
                    "matched_n_transport_q025": t_q025,
                    "matched_n_transport_q975": t_q975,
                    "technical_minus_transport_median": tech - t_med,
                    "technical_empirical_percentile_among_transport_repeats": (
                        100.0 * float(np.mean(transport <= tech))
                    ),
                    "technical_above_transport_median": int(tech >= t_med),
                    "technical_above_transport_q975": int(tech >= t_q975),
                }
            )
        else:
            same_result.update(
                {
                    "matched_n_transport_repeats": 50,
                    "matched_n_transport_median": t_med,
                    "matched_n_transport_q025": t_q025,
                    "matched_n_transport_q975": t_q975,
                    "technical_minus_transport_median": np.nan,
                    "technical_empirical_percentile_among_transport_repeats": np.nan,
                    "technical_above_transport_median": np.nan,
                    "technical_above_transport_q975": np.nan,
                }
            )

        technical_rows.extend([same_result, cross_result])

        if len(same_deg):
            degenerate_frames.append(same_deg)
        if len(cross_deg):
            degenerate_frames.append(cross_deg)

        if same_result["status"] == "COMPLETE":
            print(
                f"      same-platform rho={same_result['rho_edge_technical']:+.4f}; "
                f"transport n=100 median={t_med:+.4f} "
                f"[{t_q025:+.4f},{t_q975:+.4f}]; "
                f"Δtech-median={same_result['technical_minus_transport_median']:+.4f}"
            )
        else:
            print(
                f"      same-platform STATUS={same_result['status']}; "
                f"degenerate genes={same_result['n_degenerate_genes']}"
            )

        if cross_result["status"] == "COMPLETE":
            print(
                f"      cross-platform rho={cross_result['rho_edge_technical']:+.4f}"
            )
        else:
            print(
                f"      cross-platform STATUS={cross_result['status']}; "
                f"degenerate genes={cross_result['n_degenerate_genes']}"
            )

    print("\n[3/3] Writing technical-repeatability outputs ...")

    result_df = pd.DataFrame(technical_rows)
    deg_df = (
        pd.concat(degenerate_frames, ignore_index=True)
        if degenerate_frames
        else pd.DataFrame(
            columns=[
                "program_id",
                "pair_type",
                "Hugo_Symbol",
                "primary_exact_variable",
                "replicate_exact_variable",
            ]
        )
    )

    result_path = OUT_DIR / "scanb_technical_repeatability_summary_v1.tsv"
    deg_path = OUT_DIR / "scanb_technical_repeatability_degenerate_genes_v1.tsv"

    result_df.to_csv(result_path, sep="\t", index=False)
    deg_df.to_csv(deg_path, sep="\t", index=False)

    same_df = result_df.loc[
        result_df["pair_type"] == "SAME_PLATFORM"
    ]
    cross_df = result_df.loc[
        result_df["pair_type"] == "CROSS_PLATFORM"
    ]

    same_complete = int((same_df["status"] == "COMPLETE").sum())
    cross_complete = int((cross_df["status"] == "COMPLETE").sum())

    master = {
        "script_version": SCRIPT_VERSION,
        "status": "SCANB_TECHNICAL_REPEATABILITY_COMPLETE",
        "primary_assessable_programs": int(len(assessable_rows)),
        "same_platform_pairs": SAME_PLATFORM_PAIRS,
        "cross_platform_pairs": CROSS_PLATFORM_PAIRS,
        "same_platform_complete_programs": same_complete,
        "same_platform_not_estimable_programs": int(len(same_df) - same_complete),
        "cross_platform_complete_programs": cross_complete,
        "cross_platform_not_estimable_programs": int(len(cross_df) - cross_complete),
        "matched_n_transport_sample_size": 100,
        "matched_n_transport_repeats_per_program": 50,
        "formal_new_p_values": False,
        "primary_classification_changed": False,
        "target_outcomes_or_treatment_loaded": False,
        "summary_file": str(result_path),
        "degenerate_gene_file": str(deg_path),
    }

    master_path = OUT_DIR / "scanb_technical_repeatability_v1.json"
    master_path.write_text(
        json.dumps(master, indent=2),
        encoding="utf-8",
    )

    print("\n" + "=" * 146)
    print("05g SCAN-B TECHNICAL REPEATABILITY: COMPLETE")
    print("=" * 146)
    print(
        f"Same-platform COMPLETE programs:       "
        f"{same_complete}/{len(same_df)}"
    )
    print(
        f"Cross-platform COMPLETE programs:      "
        f"{cross_complete}/{len(cross_df)}"
    )
    print("Matched-n biological transport:        n=100, exact 05e 50-repeat distribution")
    print("Formal new p-values:                   NO")
    print("Primary classifications changed:       NO")
    print()
    print(f"Summary:    {result_path}")
    print(f"Degenerate: {deg_path}")
    print(f"Master:     {master_path}")
    print("=" * 146)

    del expr
    torch.cuda.empty_cache()
    gc.collect()


if __name__ == "__main__":
    main()
