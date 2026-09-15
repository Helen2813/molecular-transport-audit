from __future__ import annotations

import gc
import importlib.util
import inspect
import json
import math
import re
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from scipy.stats import rankdata, spearmanr


SCRIPT_VERSION = "05j-run-structural-corruption-headtohead-gpu-v1-no-cli"

PROJECT_ROOT = Path(r"C:\Users\olegk\Desktop\molecular-transport-audit")
DATA_ROOT = Path(r"D:\paper4_tcbb_data")

SOURCE_05D = PROJECT_ROOT / "scripts" / "05d_run_controlled_mapping_degradation_gpu_v1.py"

CONTRACT = (
    DATA_ROOT
    / "paper4_tcbb_structural_corruption_headtohead_execution_contract_v1"
    / "structural_corruption_headtohead_execution_contract_v1.json"
)
WGCNA_BASELINE = (
    DATA_ROOT
    / "paper4_tcbb_structural_corruption_headtohead_execution_contract_v1"
    / "wgcna_observed_structural_baseline_v1.tsv"
)

DEGRADATION_REPLICATES = (
    DATA_ROOT
    / "paper4_tcbb_controlled_mapping_degradation_v1"
    / "controlled_degradation_replicates_v1.tsv"
)

NETREP_BASELINE = (
    DATA_ROOT
    / "paper4_tcbb_postprimary_comparator_benchmark_v7"
    / "netrep_all_targets_long_v7.tsv"
)

CLASSIFICATION = (
    DATA_ROOT
    / "paper4_tcbb_final_mapping_specificity_null_v2"
    / "primary_pooled_final_classification_v2.tsv"
)
SPECIFICITY_ROOT = DATA_ROOT / "paper4_tcbb_final_mapping_specificity_null_v2"
DIRECT_ROOT = DATA_ROOT / "paper4_tcbb_primary_pooled_direct_preservation_v2"

SCANB_FEATURES = (
    DATA_ROOT
    / "paper4_tcbb_target_marginal_metrics_v2"
    / "scanb_target_marginal_matching_features_v2.tsv"
)
METABRIC_FEATURES = (
    DATA_ROOT
    / "paper4_tcbb_target_marginal_metrics_v2"
    / "metabric_target_marginal_matching_features_v2.tsv"
)

SCANB_EVAL = (
    DATA_ROOT
    / "paper4_tcbb_exact_variation_evaluability_correction_v1"
    / "scanb_target_exact_variation_evaluability_v1.tsv"
)
METABRIC_EVAL = (
    DATA_ROOT
    / "paper4_tcbb_exact_variation_evaluability_correction_v1"
    / "metabric_target_exact_variation_evaluability_v1.tsv"
)

WEIGHTS = (
    DATA_ROOT
    / "paper4_tcbb_frozen_source_programs_v1"
    / "tcga_frozen_source_program_weights_v1.tsv"
)

SOURCE_UNIVERSE = (
    DATA_ROOT
    / "paper4_tcbb_tcga_source_universe_v1"
    / "tcga_source_gene_universe_frozen_v1.tsv"
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

OUT_DIR = DATA_ROOT / "paper4_tcbb_structural_corruption_headtohead_v1"
COMBO_DIR = OUT_DIR / "per_combo"

FRACTIONS = np.array([0.0, 0.10, 0.25, 0.50, 0.75, 1.0], dtype=float)
BOOTSTRAPS = 2000
BOOTSTRAP_BASE_SEED = 20260920

MTA_REPLAY_TOL = 5e-10
BASELINE_TOL = 5e-8

TARGET_CODE = {
    "SCANB_GSE96058": 1,
    "METABRIC": 2,
}


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


def module_number(program_id: str) -> int:
    m = re.search(r"M(\d+)$", str(program_id))
    if not m:
        raise RuntimeError(f"Cannot parse module number from {program_id}")
    return int(m.group(1))


def import_exact_05d():
    spec = importlib.util.spec_from_file_location("paper4_exact_05d", SOURCE_05D)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot import exact 05d runner: {SOURCE_05D}")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)

    required_functions = [
        "prepare_problem",
        "generate_mapping",
        "quality",
        "rng_for_attempt",
        "round_half_up",
        "centered_rank",
        "spearman_against_fixed",
    ]
    missing = [x for x in required_functions if not hasattr(mod, x)]
    if missing:
        raise RuntimeError(f"Exact 05d runner is missing functions: {missing}")

    if int(getattr(mod, "BASE_SEED")) != 20260916:
        raise RuntimeError("Exact 05d BASE_SEED drift.")
    if list(getattr(mod, "FRACTIONS")) != list(FRACTIONS):
        raise RuntimeError("Exact 05d FRACTIONS drift.")

    return mod


def load_eval_symbols(path: Path) -> list[str]:
    df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    if "statistically_evaluable_corrected" not in df.columns:
        raise RuntimeError(f"Corrected evaluability column missing: {path}")
    good = df["statistically_evaluable_corrected"].map(boolish)
    symbols = [canon_symbol(x) for x in df.loc[good, "Hugo_Symbol"]]
    if len(symbols) != len(set(symbols)):
        raise RuntimeError(f"Duplicate corrected-evaluable symbols in {path}")
    return symbols


def load_features(path: Path, corrected_symbols: list[str]) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", low_memory=False)
    required_cols = {
        "Hugo_Symbol",
        "target_normalized_symbol",
        "source_composite_bin",
        "target_composite_bin",
    }
    if not required_cols.issubset(df.columns):
        raise RuntimeError(
            f"{path}: matching features missing {sorted(required_cols - set(df.columns))}"
        )

    df = df.copy()
    df["Hugo_Symbol"] = df["Hugo_Symbol"].map(canon_symbol)
    df["target_normalized_symbol"] = df["target_normalized_symbol"].map(canon_symbol)

    corrected = set(corrected_symbols)
    df = df.loc[df["Hugo_Symbol"].isin(corrected)].copy()

    if len(df) != len(corrected_symbols):
        missing = sorted(corrected - set(df["Hugo_Symbol"]))
        raise RuntimeError(
            f"{path}: expected {len(corrected_symbols)} corrected feature rows, "
            f"found {len(df)}; missing={missing[:20]}"
        )
    if df["Hugo_Symbol"].duplicated().any():
        raise RuntimeError(f"{path}: duplicate feature symbols after corrected filter.")
    if not (df["Hugo_Symbol"] == df["target_normalized_symbol"]).all():
        raise RuntimeError(f"{path}: exact-symbol mismatch after corrected filter.")

    return df.reset_index(drop=True)


def load_tcga_source() -> dict:
    universe = pd.read_csv(SOURCE_UNIVERSE, sep="\t", dtype=str).fillna("")
    if len(universe) != 10000:
        raise RuntimeError(f"Expected 10,000 frozen source genes, found {len(universe)}")

    raw = pd.read_csv(TCGA_EXPR, sep="\t", low_memory=False)
    rows = universe["source_row_index_0based"].astype(int).to_numpy()
    selected = raw.iloc[rows, :]

    expected = [canon_symbol(x) for x in universe["Hugo_Symbol"]]
    observed = [canon_symbol(x) for x in selected["Hugo_Symbol"]]
    if expected != observed:
        raise RuntimeError("TCGA source-universe row replay mismatch.")

    values = selected.iloc[:, 2:].to_numpy(dtype=np.float64)
    if values.shape[1] != 1082:
        raise RuntimeError(f"Expected 1,082 TCGA samples, found {values.shape[1]}")
    if np.any(values < 0) or not np.isfinite(values).all():
        raise RuntimeError("Invalid TCGA RSEM values.")

    x = np.log2(values + 1.0).T
    return {
        "x": x,
        "genes": expected,
        "gene_index": {g: i for i, g in enumerate(expected)},
    }


def aggregate_duplicate_rows(symbols: list[str], values: np.ndarray) -> tuple[list[str], np.ndarray]:
    df = pd.DataFrame(values)
    df.insert(0, "Hugo_Symbol", symbols)
    agg = df.groupby("Hugo_Symbol", sort=False, as_index=True).mean(numeric_only=True)
    return list(agg.index), agg.to_numpy(dtype=np.float64)


def load_scanb_target(corrected_symbols: list[str]) -> dict:
    primary = pd.read_csv(SCANB_PRIMARY, sep="\t", dtype=str).fillna("")
    titles = primary["primary_title"].tolist()
    if len(titles) != 3273:
        raise RuntimeError("SCAN-B primary manifest is not 3,273 profiles.")

    wanted = set(corrected_symbols)
    parts_symbols: list[str] = []
    parts_values: list[np.ndarray] = []

    for chunk in pd.read_csv(
        SCANB_EXPR,
        compression="gzip",
        chunksize=1200,
        low_memory=False,
    ):
        symbol_col = chunk.columns[0]
        missing_titles = [t for t in titles if t not in chunk.columns]
        if missing_titles:
            raise RuntimeError(
                f"SCAN-B expression missing primary titles, e.g. {missing_titles[:5]}"
            )

        symbols = chunk[symbol_col].map(canon_symbol)
        mask = symbols.isin(wanted)
        if not mask.any():
            continue

        arr = chunk.loc[mask, titles].to_numpy(dtype=np.float64)
        fpkm = np.maximum(np.exp2(arr) - 0.1, 0.0)
        transformed = np.log2(fpkm + 1.0)

        parts_symbols.extend(symbols.loc[mask].tolist())
        parts_values.append(transformed)

    if not parts_values:
        raise RuntimeError("No corrected SCAN-B genes loaded.")

    row_values = np.vstack(parts_values)
    genes, gene_by_sample = aggregate_duplicate_rows(parts_symbols, row_values)

    gene_index_raw = {g: i for i, g in enumerate(genes)}
    missing = [g for g in corrected_symbols if g not in gene_index_raw]
    if missing:
        raise RuntimeError(f"SCAN-B corrected genes missing after aggregation: {missing[:20]}")

    order = np.array([gene_index_raw[g] for g in corrected_symbols], dtype=np.int64)
    x = gene_by_sample[order, :].T

    # Exact corrected variation guard.
    if not np.isfinite(x).all():
        raise RuntimeError("Non-finite values in corrected SCAN-B matrix.")
    if np.any(np.max(x, axis=0) <= np.min(x, axis=0)):
        bad = [
            g for g, ok in zip(corrected_symbols, np.max(x, axis=0) > np.min(x, axis=0))
            if not ok
        ]
        raise RuntimeError(f"Corrected SCAN-B exact-constant genes reappeared: {bad[:20]}")

    return {
        "x": x,
        "genes": corrected_symbols,
        "gene_index": {g: i for i, g in enumerate(corrected_symbols)},
    }


def load_metabric_target(corrected_symbols: list[str]) -> dict:
    raw = pd.read_csv(
        METABRIC_EXPR,
        sep="\t",
        low_memory=False,
        dtype=str,
    ).fillna("")

    if "Hugo_Symbol" not in raw.columns:
        raise RuntimeError("METABRIC Hugo_Symbol column missing.")

    symbols = raw["Hugo_Symbol"].map(canon_symbol)
    wanted = set(corrected_symbols)
    mask = symbols.isin(wanted)

    sample_cols = [
        c for c in raw.columns
        if c not in {"Hugo_Symbol", "Entrez_Gene_Id"}
    ]
    if len(sample_cols) != 1980:
        raise RuntimeError(f"Expected 1,980 METABRIC samples, found {len(sample_cols)}")

    arr = raw.loc[mask, sample_cols].to_numpy(dtype=np.float64)
    genes, gene_by_sample = aggregate_duplicate_rows(symbols.loc[mask].tolist(), arr)

    gene_index_raw = {g: i for i, g in enumerate(genes)}
    missing = [g for g in corrected_symbols if g not in gene_index_raw]
    if missing:
        raise RuntimeError(f"METABRIC corrected genes missing after aggregation: {missing[:20]}")

    order = np.array([gene_index_raw[g] for g in corrected_symbols], dtype=np.int64)
    x = gene_by_sample[order, :].T

    if not np.isfinite(x).all():
        raise RuntimeError("Non-finite values in corrected METABRIC matrix.")
    if np.any(np.max(x, axis=0) <= np.min(x, axis=0)):
        bad = [
            g for g, ok in zip(corrected_symbols, np.max(x, axis=0) > np.min(x, axis=0))
            if not ok
        ]
        raise RuntimeError(f"Corrected METABRIC exact-constant genes reappeared: {bad[:20]}")

    return {
        "x": x,
        "genes": corrected_symbols,
        "gene_index": {g: i for i, g in enumerate(corrected_symbols)},
    }


def corr_gpu(x: np.ndarray, device: torch.device) -> torch.Tensor:
    xt = torch.as_tensor(x, dtype=torch.float64, device=device)
    mu = torch.mean(xt, dim=0, keepdim=True)
    xc = xt - mu
    ss = torch.sum(xc * xc, dim=0)
    if bool(torch.any(~torch.isfinite(ss))) or bool(torch.any(ss <= 0)):
        raise RuntimeError("Invalid variance in GPU correlation input.")
    z = xc / torch.sqrt(ss / (x.shape[0] - 1))
    corr = (z.T @ z) / (x.shape[0] - 1)
    corr = torch.clamp(corr, -1.0, 1.0)
    del xt, mu, xc, ss, z
    torch.cuda.empty_cache()
    return corr


def pearson_against_fixed(
    fixed_centered: torch.Tensor,
    fixed_norm: torch.Tensor,
    values: torch.Tensor,
) -> float:
    vc = values - torch.mean(values)
    vn = torch.linalg.vector_norm(vc)
    if not bool(torch.isfinite(vn)) or float(vn) <= 0:
        return float("nan")
    rho = torch.dot(fixed_centered, vc) / (fixed_norm * vn)
    return float(rho.detach().cpu())


def fixed_pearson_vector(values: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
    vc = values - torch.mean(values)
    vn = torch.linalg.vector_norm(vc)
    if not bool(torch.isfinite(vn)) or float(vn) <= 0:
        raise RuntimeError("Degenerate fixed Pearson vector.")
    return vc, vn


def f50(fracs: np.ndarray, y: np.ndarray) -> float:
    if not np.all(np.isfinite(y)):
        return float("nan")
    for i in range(1, len(fracs)):
        y0 = float(y[i - 1])
        y1 = float(y[i])
        if y0 >= 0.5 and y1 <= 0.5:
            x0 = float(fracs[i - 1])
            x1 = float(fracs[i])
            if y1 == y0:
                return x1
            return x0 + (0.5 - y0) * (x1 - x0) / (y1 - y0)
    return float("nan")


def auc(fracs: np.ndarray, y: np.ndarray) -> float:
    if not np.all(np.isfinite(y)):
        return float("nan")
    return float(np.trapz(y, fracs))


def normalized_curve(raw_medians: np.ndarray) -> np.ndarray:
    floor = float(raw_medians[-1])
    ceiling = float(raw_medians[0])
    denom = ceiling - floor
    if not np.isfinite(denom) or abs(denom) < 1e-12:
        return np.full_like(raw_medians, np.nan, dtype=float)
    return (raw_medians - floor) / denom


def bootstrap_rng(target: str, program_id: str, bootstrap_id: int) -> np.random.Generator:
    ss = np.random.SeedSequence(
        [
            BOOTSTRAP_BASE_SEED,
            TARGET_CODE[target],
            module_number(program_id),
            int(bootstrap_id),
        ]
    )
    return np.random.default_rng(ss)


def compare_float(a: float, b: float, tol: float, label: str) -> None:
    if not (np.isfinite(a) and np.isfinite(b)):
        if np.isnan(a) and np.isnan(b):
            return
        raise RuntimeError(f"{label}: non-finite mismatch {a} vs {b}")
    if abs(a - b) > tol:
        raise RuntimeError(f"{label}: mismatch {a:.16g} vs {b:.16g}; |Δ|={abs(a-b):.3e}")


def target_bundle(target_name: str):
    if target_name == "SCANB_GSE96058":
        corrected = load_eval_symbols(SCANB_EVAL)
        if len(corrected) != 9220:
            raise RuntimeError(f"SCAN-B corrected gene count drift: {len(corrected)}")
        features = load_features(SCANB_FEATURES, corrected)
        target = load_scanb_target(corrected)
    elif target_name == "METABRIC":
        corrected = load_eval_symbols(METABRIC_EVAL)
        if len(corrected) != 8485:
            raise RuntimeError(f"METABRIC corrected gene count drift: {len(corrected)}")
        features = load_features(METABRIC_FEATURES, corrected)
        target = load_metabric_target(corrected)
    else:
        raise RuntimeError(f"Unsupported target: {target_name}")
    return target, features


def main() -> None:
    print("=" * 160)
    print("Paper 4 / TCBB - structural corruption head-to-head: MTA vs NetRep vs WGCNA")
    print("=" * 160)
    print(f"Script version: {SCRIPT_VERSION}")
    print()
    print("Frozen execution:")
    print("  Exact 05d generator imported/reused:                YES")
    print("  All 11,022 saved 05d rows replayed:                 YES")
    print("  MTA replay required before comparator acceptance:   YES")
    print("  NetRep endpoint:                                    cor.cor")
    print("  WGCNA endpoints:                                    cor.cor + cor.kIM")
    print("  Comparator p-values/permutations in ladder:         NO")
    print("  Paired bootstrap replicates:                        2,000")
    print("  Primary MTA classifications changed:                NO")
    print("=" * 160)

    for p in [
        SOURCE_05D,
        CONTRACT,
        WGCNA_BASELINE,
        DEGRADATION_REPLICATES,
        NETREP_BASELINE,
        CLASSIFICATION,
        WEIGHTS,
        SOURCE_UNIVERSE,
        TCGA_EXPR,
    ]:
        require(p)

    c = json.loads(CONTRACT.read_text(encoding="utf-8"))
    if c.get("status") != "FROZEN_BEFORE_FIRST_CORRUPTION_HEADTOHEAD_COMPARATOR_STATISTIC":
        raise RuntimeError("05j2 execution contract has unexpected status.")

    exact05d = import_exact_05d()

    if not torch.cuda.is_available():
        raise RuntimeError("CUDA is unavailable.")
    device = torch.device("cuda:0")
    prop = torch.cuda.get_device_properties(0)
    print(
        f"CUDA device: {prop.name}; VRAM={prop.total_memory/1024**3:.2f} GB; "
        f"capability={prop.major}.{prop.minor}"
    )

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    COMBO_DIR.mkdir(parents=True, exist_ok=True)

    reps = pd.read_csv(DEGRADATION_REPLICATES, sep="\t", low_memory=False)
    if len(reps) != 11022:
        raise RuntimeError(f"Expected 11,022 saved 05d rows, found {len(reps)}")

    netrep = pd.read_csv(NETREP_BASELINE, sep="\t", low_memory=False)
    netrep = netrep.loc[netrep["statistic"] == "cor.cor"].copy()

    wgcna = pd.read_csv(WGCNA_BASELINE, sep="\t", low_memory=False)

    classification = pd.read_csv(CLASSIFICATION, sep="\t", low_memory=False)
    classification = classification.loc[
        classification["primary_assessable"].map(boolish)
    ].copy()
    if len(classification) != 22:
        raise RuntimeError(f"Expected 22 primary-assessable pairs, found {len(classification)}")

    weights = pd.read_csv(WEIGHTS, sep="\t", dtype=str).fillna("")
    weights["Hugo_Symbol"] = weights["Hugo_Symbol"].map(canon_symbol)

    source = load_tcga_source()

    all_raw = []
    all_fraction = []
    all_metrics = []
    all_bootstrap = []

    for target_name in ["SCANB_GSE96058", "METABRIC"]:
        print("\n" + "=" * 160)
        print(f"TARGET: {target_name}")
        print("=" * 160)

        target, features = target_bundle(target_name)
        print(
            f"  target matrix: {target['x'].shape[0]:,} samples x "
            f"{target['x'].shape[1]:,} corrected genes"
        )

        print("  computing full target Pearson correlation on GPU ...")
        target_corr = corr_gpu(target["x"], device)
        del target["x"]
        gc.collect()
        torch.cuda.empty_cache()

        rows_target = classification.loc[classification["target"] == target_name].copy()

        for ci, class_row in enumerate(rows_target.to_dict(orient="records"), start=1):
            program_id = str(class_row["program_id"])
            final_path = COMBO_DIR / f"{target_name}__{program_id}__structural_headtohead_v1.tsv"
            meta_path = COMBO_DIR / f"{target_name}__{program_id}__structural_headtohead_v1.json"

            print(
                f"\n  [{ci:02d}/{len(rows_target):02d}] {program_id}"
            )

            # Completed combo replay can be reused exactly.
            if final_path.exists() and meta_path.exists():
                meta = json.loads(meta_path.read_text(encoding="utf-8"))
                if meta.get("status") == "COMPLETE" and meta.get("script_version") == SCRIPT_VERSION:
                    print("      completed combo found — reusing")
                    combo_df = pd.read_csv(final_path, sep="\t", low_memory=False)
                    all_raw.append(combo_df)
                    continue

            gene_file = (
                DIRECT_ROOT
                / "per_program"
                / target_name
                / f"{program_id}_evaluable_genes_v2.tsv"
            )
            require(gene_file)
            eg = pd.read_csv(gene_file, sep="\t", dtype=str).fillna("")
            slot_genes = [canon_symbol(x) for x in eg["Hugo_Symbol"]]
            m = len(slot_genes)

            full_program_genes = {
                canon_symbol(x)
                for x in weights.loc[
                    weights["program_id"] == program_id,
                    "Hugo_Symbol",
                ]
            }

            # Source expression and full module correlation.
            source_cols = np.array(
                [source["gene_index"][g] for g in slot_genes],
                dtype=np.int64,
            )
            source_x = source["x"][:, source_cols]
            source_corr = corr_gpu(source_x, device)
            del source_x

            tri_i_np, tri_j_np = np.triu_indices(m, k=1)
            tri_i_t = torch.as_tensor(tri_i_np, dtype=torch.long, device=device)
            tri_j_t = torch.as_tensor(tri_j_np, dtype=torch.long, device=device)

            source_all_edges = source_corr[tri_i_t, tri_j_t]
            source_edge_centered, source_edge_norm = fixed_pearson_vector(source_all_edges)

            # WGCNA signed modulePreservation adjacency uses power 12 internally.
            source_adj12 = torch.pow((1.0 + source_corr) / 2.0, 12.0)
            source_kim = torch.sum(source_adj12, dim=0)
            source_kim_centered, source_kim_norm = fixed_pearson_vector(source_kim)

            null_path = (
                SPECIFICITY_ROOT
                / "per_program"
                / target_name
                / f"{program_id}_mapping_specificity_null_v2.npz"
            )
            require(null_path)
            nz = np.load(null_path)

            ii = np.asarray(nz["edge_slot_i_0based"], dtype=np.int64)
            jj = np.asarray(nz["edge_slot_j_0based"], dtype=np.int64)
            source_edges_subset = np.asarray(nz["source_edge_subset"], dtype=np.float64)
            source_rank, source_rank_norm = exact05d.centered_rank(source_edges_subset)

            ii_t = torch.as_tensor(ii, dtype=torch.long, device=device)
            jj_t = torch.as_tensor(jj, dtype=torch.long, device=device)

            correct_target_idx = np.array(
                [target["gene_index"][g] for g in slot_genes],
                dtype=np.int64,
            )

            # Baseline references.
            nr = netrep.loc[
                (netrep["target"] == target_name)
                & (netrep["program_id"] == program_id),
                "observed",
            ]
            if len(nr) != 1:
                raise RuntimeError(f"{target_name} {program_id}: NetRep baseline row count={len(nr)}")
            expected_netrep = float(nr.iloc[0])

            wr = wgcna.loc[
                (wgcna["target"] == target_name)
                & (wgcna["program_id"] == program_id)
            ]
            if len(wr) != 1:
                raise RuntimeError(f"{target_name} {program_id}: WGCNA baseline row count={len(wr)}")
            expected_wgcna_corcor = float(wr["cor.cor"].iloc[0])
            expected_wgcna_corkim = float(wr["cor.kIM"].iloc[0])

            saved = reps.loc[
                (reps["target"] == target_name)
                & (reps["program_id"] == program_id)
            ].copy()
            saved = saved.sort_values(["fraction", "replicate_id"]).reset_index(drop=True)
            if len(saved) != 501:
                raise RuntimeError(f"{target_name} {program_id}: expected 501 05d rows, got {len(saved)}")

            out_rows = []
            max_mta_diff = 0.0

            for ri, row in enumerate(saved.to_dict(orient="records"), start=1):
                fraction = float(row["fraction"])
                replicate_id = int(row["replicate_id"])
                attempt_id = int(row["attempt_id"])

                if fraction == 0.0:
                    mapped = correct_target_idx.copy()
                else:
                    k = max(1, min(m, exact05d.round_half_up(fraction * m)))
                    if k != int(row["n_corrupted_slots"]):
                        raise RuntimeError(
                            f"{target_name} {program_id} f={fraction} rep={replicate_id}: "
                            f"corrupted count replay mismatch {k} vs {row['n_corrupted_slots']}"
                        )

                    rng = exact05d.rng_for_attempt(
                        target_name,
                        program_id,
                        fraction,
                        attempt_id,
                    )
                    corrupted = np.sort(
                        rng.choice(
                            m,
                            size=k,
                            replace=False,
                        )
                    )
                    problem = exact05d.prepare_problem(
                        slot_genes=slot_genes,
                        corrupted_slot_idx=corrupted,
                        full_program_genes=full_program_genes,
                        features=features,
                        target_gene_index=target["gene_index"],
                    )
                    mapping_rows, distances, success, reason = exact05d.generate_mapping(
                        problem,
                        rng,
                    )
                    if not success:
                        raise RuntimeError(
                            f"{target_name} {program_id} f={fraction} rep={replicate_id}: "
                            f"saved valid attempt replay failed: {reason}"
                        )
                    q = exact05d.quality(distances)
                    if not q["panel_valid"]:
                        raise RuntimeError(
                            f"{target_name} {program_id} f={fraction} rep={replicate_id}: "
                            "saved valid attempt replayed as invalid"
                        )

                    compare_float(
                        float(q["fraction_distance_0"]),
                        float(row["fraction_distance_0"]),
                        1e-12,
                        "fraction_distance_0",
                    )
                    compare_float(
                        float(q["fraction_distance_le1"]),
                        float(row["fraction_distance_le1"]),
                        1e-12,
                        "fraction_distance_le1",
                    )
                    if int(q["maximum_distance"]) != int(row["maximum_distance"]):
                        raise RuntimeError("maximum_distance replay mismatch")
                    compare_float(
                        float(q["mean_distance"]),
                        float(row["mean_distance"]),
                        1e-12,
                        "mean_distance",
                    )

                    candidate_rows = problem["candidate"].iloc[mapping_rows]
                    replacement_target_idx = candidate_rows[
                        "target_corr_index"
                    ].to_numpy(dtype=np.int64)

                    mapped = correct_target_idx.copy()
                    mapped[corrupted] = replacement_target_idx

                if len(np.unique(mapped)) != len(mapped):
                    raise RuntimeError(
                        f"{target_name} {program_id} f={fraction} rep={replicate_id}: "
                        "duplicate mapped target genes"
                    )

                mapped_t = torch.as_tensor(mapped, dtype=torch.long, device=device)
                sub = target_corr.index_select(0, mapped_t).index_select(1, mapped_t)

                # Exact 05d MTA replay on frozen subset.
                mta_edges = sub[ii_t, jj_t].detach().cpu().numpy()
                mta_rho = exact05d.spearman_against_fixed(
                    source_rank,
                    source_rank_norm,
                    mta_edges,
                )
                saved_rho = float(row["rho_edge"])
                diff = abs(mta_rho - saved_rho)
                max_mta_diff = max(max_mta_diff, diff)
                if diff > MTA_REPLAY_TOL:
                    raise RuntimeError(
                        f"{target_name} {program_id} f={fraction} rep={replicate_id}: "
                        f"MTA replay mismatch {mta_rho:.16g} vs {saved_rho:.16g}; |Δ|={diff:.3e}"
                    )

                # NetRep cor.cor and WGCNA cor.cor are both Pearson concordance
                # of the full within-module correlation edge vectors.
                target_all_edges = sub[tri_i_t, tri_j_t]
                corcor = pearson_against_fixed(
                    source_edge_centered,
                    source_edge_norm,
                    target_all_edges,
                )

                # WGCNA cor.kIM: signed adjacency ((1+r)/2)^12.
                target_adj12 = torch.pow((1.0 + sub) / 2.0, 12.0)
                target_kim = torch.sum(target_adj12, dim=0)
                corkim = pearson_against_fixed(
                    source_kim_centered,
                    source_kim_norm,
                    target_kim,
                )

                if fraction == 0.0:
                    compare_float(
                        corcor,
                        expected_netrep,
                        BASELINE_TOL,
                        f"{target_name} {program_id} NetRep fraction0 cor.cor",
                    )
                    compare_float(
                        corcor,
                        expected_wgcna_corcor,
                        BASELINE_TOL,
                        f"{target_name} {program_id} WGCNA fraction0 cor.cor",
                    )
                    compare_float(
                        corkim,
                        expected_wgcna_corkim,
                        BASELINE_TOL,
                        f"{target_name} {program_id} WGCNA fraction0 cor.kIM",
                    )

                out_rows.append(
                    {
                        "target": target_name,
                        "program_id": program_id,
                        "fraction": fraction,
                        "replicate_id": replicate_id,
                        "attempt_id": attempt_id,
                        "n_slots": m,
                        "n_corrupted_slots": int(row["n_corrupted_slots"]),
                        "mta_rho_edge_saved": saved_rho,
                        "mta_rho_edge_replayed": mta_rho,
                        "netrep_cor_cor": corcor,
                        "wgcna_cor_cor": corcor,
                        "wgcna_cor_kIM": corkim,
                    }
                )

                del mapped_t, sub, target_all_edges, target_adj12, target_kim

                if ri % 50 == 0 or ri == len(saved):
                    print(
                        f"      replay/comparator {ri:3d}/{len(saved)}; "
                        f"max MTA |Δ|={max_mta_diff:.2e}"
                    )

            combo_df = pd.DataFrame(out_rows)
            combo_df.to_csv(final_path, sep="\t", index=False)

            meta = {
                "script_version": SCRIPT_VERSION,
                "status": "COMPLETE",
                "target": target_name,
                "program_id": program_id,
                "rows": int(len(combo_df)),
                "max_mta_replay_abs_diff": float(max_mta_diff),
                "netrep_fraction0_reference": expected_netrep,
                "wgcna_corcor_fraction0_reference": expected_wgcna_corcor,
                "wgcna_corkim_fraction0_reference": expected_wgcna_corkim,
                "result_file": str(final_path),
            }
            meta_path.write_text(json.dumps(meta, indent=2), encoding="utf-8")
            all_raw.append(combo_df)

            del (
                source_corr,
                tri_i_t,
                tri_j_t,
                source_all_edges,
                source_edge_centered,
                source_edge_norm,
                source_adj12,
                source_kim,
                source_kim_centered,
                source_kim_norm,
                ii_t,
                jj_t,
            )
            torch.cuda.empty_cache()
            gc.collect()

        del target_corr, target, features
        torch.cuda.empty_cache()
        gc.collect()

    raw = pd.concat(all_raw, ignore_index=True)
    raw_path = OUT_DIR / "structural_headtohead_replicates_v1.tsv"
    raw.to_csv(raw_path, sep="\t", index=False)

    print("\n" + "=" * 160)
    print("Replay and raw structural comparator calculations COMPLETE.")
    print("Computing frozen normalized operating summaries + paired bootstrap ...")
    print("=" * 160)

    methods = [
        ("MTA", "mta_rho_edge_saved"),
        ("NetRep_cor.cor", "netrep_cor_cor"),
        ("WGCNA_cor.cor", "wgcna_cor_cor"),
        ("WGCNA_cor.kIM", "wgcna_cor_kIM"),
    ]

    fraction_rows = []
    metric_rows = []
    bootstrap_rows = []

    for (target_name, program_id), g in raw.groupby(["target", "program_id"], sort=True):
        # Raw medians and method-specific normalization.
        method_medians: dict[str, np.ndarray] = {}
        method_norm: dict[str, np.ndarray] = {}

        for method, col in methods:
            meds = np.array(
                [
                    float(np.median(g.loc[g["fraction"] == f, col].to_numpy(dtype=float)))
                    for f in FRACTIONS
                ],
                dtype=float,
            )
            norm = normalized_curve(meds)
            method_medians[method] = meds
            method_norm[method] = norm

            trend = float(spearmanr(FRACTIONS, meds).statistic)
            metric_rows.append(
                {
                    "target": target_name,
                    "program_id": program_id,
                    "method": method,
                    "auc_normalized": auc(FRACTIONS, norm),
                    "f50": f50(FRACTIONS, norm),
                    "spearman_fraction_vs_median_raw": trend,
                    "raw_zero": meds[0],
                    "raw_full_corruption_median": meds[-1],
                }
            )

            for f, raw_med, norm_med in zip(FRACTIONS, meds, norm):
                fraction_rows.append(
                    {
                        "target": target_name,
                        "program_id": program_id,
                        "method": method,
                        "fraction": float(f),
                        "raw_median": float(raw_med),
                        "normalized_median": float(norm_med),
                    }
                )

        # Paired bootstrap against MTA.
        nonzero_groups = {
            float(f): g.loc[g["fraction"] == f].sort_values("replicate_id").reset_index(drop=True)
            for f in FRACTIONS[1:]
        }
        for f, gg in nonzero_groups.items():
            if len(gg) != 100:
                raise RuntimeError(
                    f"{target_name} {program_id}: fraction {f} has {len(gg)} rows, expected 100"
                )

        boot_metrics: dict[str, list[tuple[float, float]]] = {
            method: [] for method, _ in methods
        }

        zero_row = g.loc[g["fraction"] == 0.0]
        if len(zero_row) != 1:
            raise RuntimeError(f"{target_name} {program_id}: fraction 0 row count != 1")

        for b in range(1, BOOTSTRAPS + 1):
            rng = bootstrap_rng(target_name, program_id, b)
            sampled = {
                float(f): rng.integers(0, 100, size=100, endpoint=False)
                for f in FRACTIONS[1:]
            }

            for method, col in methods:
                meds = [float(zero_row[col].iloc[0])]
                for f in FRACTIONS[1:]:
                    vals = nonzero_groups[float(f)][col].to_numpy(dtype=float)
                    meds.append(float(np.median(vals[sampled[float(f)]])))
                meds_arr = np.asarray(meds, dtype=float)
                norm = normalized_curve(meds_arr)
                boot_metrics[method].append(
                    (
                        auc(FRACTIONS, norm),
                        f50(FRACTIONS, norm),
                    )
                )

        mta_boot = np.asarray(boot_metrics["MTA"], dtype=float)

        for comp_method in ["NetRep_cor.cor", "WGCNA_cor.cor", "WGCNA_cor.kIM"]:
            comp_boot = np.asarray(boot_metrics[comp_method], dtype=float)

            delta_auc = mta_boot[:, 0] - comp_boot[:, 0]
            finite_auc = np.isfinite(delta_auc)

            delta_f50 = mta_boot[:, 1] - comp_boot[:, 1]
            finite_f50 = np.isfinite(delta_f50)

            observed_mta = next(
                r for r in metric_rows
                if r["target"] == target_name
                and r["program_id"] == program_id
                and r["method"] == "MTA"
            )
            observed_comp = next(
                r for r in metric_rows
                if r["target"] == target_name
                and r["program_id"] == program_id
                and r["method"] == comp_method
            )

            row = {
                "target": target_name,
                "program_id": program_id,
                "comparator": comp_method,
                "delta_auc_mta_minus_comparator": (
                    float(observed_mta["auc_normalized"] - observed_comp["auc_normalized"])
                ),
                "delta_auc_boot_q025": (
                    float(np.quantile(delta_auc[finite_auc], 0.025))
                    if finite_auc.any() else np.nan
                ),
                "delta_auc_boot_q975": (
                    float(np.quantile(delta_auc[finite_auc], 0.975))
                    if finite_auc.any() else np.nan
                ),
                "delta_auc_boot_valid": int(finite_auc.sum()),
                "delta_f50_mta_minus_comparator": (
                    float(observed_mta["f50"] - observed_comp["f50"])
                    if np.isfinite(observed_mta["f50"]) and np.isfinite(observed_comp["f50"])
                    else np.nan
                ),
                "delta_f50_boot_q025": (
                    float(np.quantile(delta_f50[finite_f50], 0.025))
                    if finite_f50.any() else np.nan
                ),
                "delta_f50_boot_q975": (
                    float(np.quantile(delta_f50[finite_f50], 0.975))
                    if finite_f50.any() else np.nan
                ),
                "delta_f50_boot_valid": int(finite_f50.sum()),
            }

            if (
                np.isfinite(row["delta_auc_boot_q025"])
                and np.isfinite(row["delta_auc_boot_q975"])
            ):
                if row["delta_auc_boot_q975"] < 0:
                    row["auc_interpretation"] = "MTA_MORE_RESPONSIVE"
                elif row["delta_auc_boot_q025"] > 0:
                    row["auc_interpretation"] = "COMPARATOR_MORE_RESPONSIVE"
                else:
                    row["auc_interpretation"] = "NO_CLEAR_DIFFERENCE"
            else:
                row["auc_interpretation"] = "NOT_ESTIMABLE"

            bootstrap_rows.append(row)

    frac_df = pd.DataFrame(fraction_rows)
    metrics_df = pd.DataFrame(metric_rows)
    boot_df = pd.DataFrame(bootstrap_rows)

    frac_path = OUT_DIR / "structural_headtohead_fraction_summary_v1.tsv"
    metrics_path = OUT_DIR / "structural_headtohead_operating_metrics_v1.tsv"
    boot_path = OUT_DIR / "structural_headtohead_paired_contrasts_v1.tsv"

    frac_df.to_csv(frac_path, sep="\t", index=False)
    metrics_df.to_csv(metrics_path, sep="\t", index=False)
    boot_df.to_csv(boot_path, sep="\t", index=False)

    count_summary = (
        boot_df.groupby(["comparator", "auc_interpretation"])
        .size()
        .reset_index(name="n_target_program_pairs")
    )
    count_path = OUT_DIR / "structural_headtohead_auc_interpretation_counts_v1.tsv"
    count_summary.to_csv(count_path, sep="\t", index=False)

    master = {
        "script_version": SCRIPT_VERSION,
        "status": "STRUCTURAL_CORRUPTION_HEADTOHEAD_COMPLETE",
        "assessable_target_program_pairs": 22,
        "saved_05d_rows_replayed": int(len(raw)),
        "fractions": FRACTIONS.tolist(),
        "replicates_per_nonzero_fraction": 100,
        "paired_bootstrap_replicates": BOOTSTRAPS,
        "methods": [x[0] for x in methods],
        "primary_operating_contrast": "AUC_MTA - AUC_comparator; negative favors MTA",
        "primary_mta_classification_changed": False,
        "replicate_file": str(raw_path),
        "fraction_summary_file": str(frac_path),
        "operating_metrics_file": str(metrics_path),
        "paired_contrasts_file": str(boot_path),
        "interpretation_counts_file": str(count_path),
    }
    master_path = OUT_DIR / "structural_corruption_headtohead_v1.json"
    master_path.write_text(json.dumps(master, indent=2), encoding="utf-8")

    print("\n" + "=" * 160)
    print("05j STRUCTURAL CORRUPTION HEAD-TO-HEAD: COMPLETE")
    print("=" * 160)
    print(f"05d rows exactly replayed:                 {len(raw):,}/11,022")
    print("Fraction-0 NetRep/WGCNA validation:        PASS for all completed pairs")
    print("Paired bootstrap:                          2,000/pair")
    print("Primary classifications changed:           NO")
    print()
    print("AUC interpretation counts:")
    print(count_summary.to_string(index=False))
    print()
    print(f"Replicates:  {raw_path}")
    print(f"Fractions:   {frac_path}")
    print(f"Metrics:     {metrics_path}")
    print(f"Contrasts:   {boot_path}")
    print(f"Master:      {master_path}")
    print("=" * 160)


if __name__ == "__main__":
    main()
