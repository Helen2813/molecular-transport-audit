from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pandas as pd


SCRIPT_VERSION = "01b-freeze-expression-input-contract-v1-no-cli"

DATA_ROOT = Path(r"D:\paper4_tcbb_data")
AUDIT_ROOT = DATA_ROOT / "paper4_tcbb_input_audit_v1"
STAGED_ROOT = AUDIT_ROOT / "staged_continuous_inputs"

TCGA_EXPR = (
    STAGED_ROOT
    / "TCGA_BRCA_PanCanAtlas2018"
    / "data_mrna_seq_v2_rsem.txt"
)
METABRIC_EXPR = (
    STAGED_ROOT
    / "METABRIC"
    / "data_mrna_illumina_microarray.txt"
)
SCANB_EXPR = (
    DATA_ROOT
    / "SCANB_GSE96058"
    / "GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv.gz"
)
SCANB_TECH = AUDIT_ROOT / "scanb_technical_inventory_v1.tsv"

OUT_DIR = DATA_ROOT / "paper4_tcbb_expression_contract_v1"

CHUNK_ROWS = 256
EPS = 1e-12


def require(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Required file does not exist: {path}")


def banner() -> None:
    print("=" * 116)
    print("Paper 4 / TCBB - freeze expression-input representation BEFORE target preservation analysis")
    print("=" * 116)
    print(f"Script version: {SCRIPT_VERSION}")
    print()
    print("Scientific guard:")
    print("  Target outcomes loaded:                              NO")
    print("  Target clinical characteristics loaded:              NO")
    print("  Target preservation statistics calculated:           NO")
    print("  Target correlation matrices calculated:              NO")
    print("  Module selection/tuning from target data:             NO")
    print("  Operations here: deterministic expression-scale QC and contract freeze only")
    print("=" * 116)


def header_info(path: Path, sep: str, compression: str | None, id_columns: int) -> dict:
    df = pd.read_csv(path, sep=sep, compression=compression, nrows=0, low_memory=False)
    cols = list(df.columns)
    return {
        "id_columns": cols[:id_columns],
        "sample_columns": len(cols) - id_columns,
        "first_sample_columns": cols[id_columns:id_columns + 5],
    }


def verify_tcga_log2p1(path: Path) -> dict:
    min_raw = math.inf
    max_raw = -math.inf
    min_tx = math.inf
    max_tx = -math.inf
    finite_n = 0
    negative_n = 0

    hdr = pd.read_csv(path, sep="\t", nrows=0, low_memory=False)
    sample_cols = list(hdr.columns)[2:]

    for df in pd.read_csv(path, sep="\t", chunksize=CHUNK_ROWS, low_memory=False):
        arr = df[sample_cols].apply(pd.to_numeric, errors="coerce").to_numpy(float)
        finite = arr[np.isfinite(arr)]
        if finite.size == 0:
            continue
        finite_n += int(finite.size)
        negative_n += int((finite < 0).sum())
        min_raw = min(min_raw, float(finite.min()))
        max_raw = max(max_raw, float(finite.max()))

        if np.any(finite < 0):
            raise RuntimeError("TCGA contains negative values; log2(x+1) contract cannot be frozen.")

        tx = np.log2(finite + 1.0)
        min_tx = min(min_tx, float(tx.min()))
        max_tx = max(max_tx, float(tx.max()))

    return {
        "raw_min": float(min_raw),
        "raw_max": float(max_raw),
        "raw_negative_fraction": float(negative_n / finite_n) if finite_n else None,
        "frozen_transform": "log2(x + 1)",
        "transformed_min": float(min_tx),
        "transformed_max": float(max_tx),
    }


def verify_scanb_relog1(path: Path) -> dict:
    """
    Published SCAN-B file is treated as y = log2(FPKM + 0.1).
    Reconstruct FPKM deterministically and harmonize the RNA-seq pseudocount:
        FPKM = max(2**y - 0.1, 0)
        primary = log2(FPKM + 1)
    No target preservation statistic is computed.
    """
    hdr = pd.read_csv(path, compression="gzip", nrows=0, low_memory=False)
    sample_cols = list(hdr.columns)[1:]

    min_y = math.inf
    max_y = -math.inf
    min_fpkm = math.inf
    max_fpkm = -math.inf
    min_tx = math.inf
    max_tx = -math.inf
    finite_n = 0
    below_floor_n = 0
    floor_hits = 0

    theoretical_floor = math.log2(0.1)

    for df in pd.read_csv(
        path, compression="gzip", chunksize=CHUNK_ROWS, low_memory=False
    ):
        arr = df[sample_cols].apply(pd.to_numeric, errors="coerce").to_numpy(float)
        y = arr[np.isfinite(arr)]
        if y.size == 0:
            continue

        finite_n += int(y.size)
        min_y = min(min_y, float(y.min()))
        max_y = max(max_y, float(y.max()))
        below_floor_n += int((y < theoretical_floor - 1e-9).sum())
        floor_hits += int((np.abs(y - theoretical_floor) <= 1e-12).sum())

        fpkm = np.exp2(y) - 0.1
        # Tiny negative values can occur from roundoff at the encoded zero floor.
        if np.nanmin(fpkm) < -1e-8:
            raise RuntimeError(
                "SCAN-B back-transformation produced values materially below zero; "
                "the assumed log2(FPKM + 0.1) representation needs review."
            )
        fpkm = np.maximum(fpkm, 0.0)

        tx = np.log2(fpkm + 1.0)

        min_fpkm = min(min_fpkm, float(fpkm.min()))
        max_fpkm = max(max_fpkm, float(fpkm.max()))
        min_tx = min(min_tx, float(tx.min()))
        max_tx = max(max_tx, float(tx.max()))

    return {
        "published_scale_min": float(min_y),
        "published_scale_max": float(max_y),
        "theoretical_log2_0p1_floor": theoretical_floor,
        "values_below_floor_tolerance": int(below_floor_n),
        "exact_floor_hits_tolerance": int(floor_hits),
        "finite_values": int(finite_n),
        "assumed_published_transform": "y = log2(FPKM + 0.1)",
        "back_transform": "FPKM = max(2**y - 0.1, 0)",
        "back_transformed_fpkm_min": float(min_fpkm),
        "back_transformed_fpkm_max": float(max_fpkm),
        "frozen_primary_transform": "log2(FPKM + 1)",
        "primary_transformed_min": float(min_tx),
        "primary_transformed_max": float(max_tx),
    }


def verify_metabric_as_provided(path: Path) -> dict:
    hdr = pd.read_csv(path, sep="\t", nrows=0, low_memory=False)
    sample_cols = list(hdr.columns)[2:]

    min_v = math.inf
    max_v = -math.inf
    finite_n = 0
    negative_n = 0

    for df in pd.read_csv(path, sep="\t", chunksize=CHUNK_ROWS, low_memory=False):
        arr = df[sample_cols].apply(pd.to_numeric, errors="coerce").to_numpy(float)
        x = arr[np.isfinite(arr)]
        if x.size == 0:
            continue
        finite_n += int(x.size)
        negative_n += int((x < 0).sum())
        min_v = min(min_v, float(x.min()))
        max_v = max(max_v, float(x.max()))

    return {
        "raw_min": float(min_v),
        "raw_max": float(max_v),
        "negative_fraction": float(negative_n / finite_n) if finite_n else None,
        "frozen_primary_transform": "as provided: log2 Illumina HT-12 v3 intensity",
    }


def scanb_role_summary(path: Path) -> tuple[pd.DataFrame, dict]:
    df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    title = df["sample_title"].str.lower()

    # This is only a technical-role inventory; no clinical fields are read.
    is_rep = title.str.contains(r"repl|replicate|repeat|resequenc", regex=True)
    out = df.copy()
    out["is_apparent_technical_replicate"] = is_rep.astype(int)

    counts = (
        out.groupby(["platform_id", "is_apparent_technical_replicate"], dropna=False)
        .size()
        .reset_index(name="n_profiles")
    )

    info = {
        "total_profiles": int(len(out)),
        "apparent_technical_replicates": int(is_rep.sum()),
        "apparent_primary_profiles": int((~is_rep).sum()),
        "platform_by_role": counts.to_dict(orient="records"),
        "selection_status": (
            "replicate keyword identifies the expected 136 profiles, but exact "
            "pairing/which member is primary remains TO_FREEZE after title-pattern review"
        ),
    }
    return out, info


def main() -> None:
    banner()
    for p in [TCGA_EXPR, METABRIC_EXPR, SCANB_EXPR, SCANB_TECH]:
        require(p)

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print("\n[1/4] Inspecting input schemas ...")
    schemas = {
        "TCGA_BRCA": header_info(TCGA_EXPR, "\t", None, 2),
        "SCANB_GSE96058": header_info(SCANB_EXPR, ",", "gzip", 1),
        "METABRIC": header_info(METABRIC_EXPR, "\t", None, 2),
    }
    for k, v in schemas.items():
        print(f"  {k}: id_columns={v['id_columns']} sample_columns={v['sample_columns']}")

    print("\n[2/4] Verifying deterministic expression transforms ...")
    tcga = verify_tcga_log2p1(TCGA_EXPR)
    print(
        f"  TCGA-BRCA: raw [{tcga['raw_min']}, {tcga['raw_max']}] -> "
        f"log2(x+1) [{tcga['transformed_min']:.6g}, {tcga['transformed_max']:.6g}]"
    )

    scanb = verify_scanb_relog1(SCANB_EXPR)
    print(
        f"  SCAN-B: published [{scanb['published_scale_min']}, {scanb['published_scale_max']}] ; "
        f"values below log2(0.1) tolerance={scanb['values_below_floor_tolerance']}"
    )
    print(
        f"          reconstructed FPKM [{scanb['back_transformed_fpkm_min']:.6g}, "
        f"{scanb['back_transformed_fpkm_max']:.6g}] -> "
        f"log2(FPKM+1) [{scanb['primary_transformed_min']:.6g}, "
        f"{scanb['primary_transformed_max']:.6g}]"
    )

    metabric = verify_metabric_as_provided(METABRIC_EXPR)
    print(
        f"  METABRIC: keep published log2 intensity scale "
        f"[{metabric['raw_min']}, {metabric['raw_max']}]"
    )

    print("\n[3/4] Auditing SCAN-B technical-profile roles ...")
    role_df, role_info = scanb_role_summary(SCANB_TECH)
    role_out = OUT_DIR / "scanb_profile_roles_pre_freeze_v1.tsv"
    role_df.to_csv(role_out, sep="\t", index=False)

    repl_examples = role_df[
        role_df["is_apparent_technical_replicate"].astype(int) == 1
    ][["geo_accession", "sample_title", "platform_id"]].head(40)
    repl_examples_out = OUT_DIR / "scanb_replicate_title_examples_v1.tsv"
    repl_examples.to_csv(repl_examples_out, sep="\t", index=False)

    print(
        f"  profiles={role_info['total_profiles']} "
        f"primary-like={role_info['apparent_primary_profiles']} "
        f"replicate-like={role_info['apparent_technical_replicates']}"
    )
    print(f"  platform x role: {role_info['platform_by_role']}")

    print("\n[4/4] Freezing expression-input contract ...")
    contract = {
        "contract_id": "paper4-tcbb-expression-input-contract-v1",
        "script_version": SCRIPT_VERSION,
        "status": "FROZEN_FOR_EXPRESSION_REPRESENTATION_ONLY",
        "scientific_guard": {
            "target_outcomes_loaded": False,
            "target_clinical_characteristics_loaded": False,
            "target_preservation_statistics_calculated": False,
            "target_correlation_matrices_calculated": False,
            "target_driven_tuning": False,
        },
        "schemas": schemas,
        "expression_representation": {
            "TCGA_BRCA": tcga,
            "SCANB_GSE96058": scanb,
            "METABRIC": metabric,
            "post_transform_gene_standardization": (
                "z-score each evaluable gene across samples within each cohort"
            ),
        },
        "correlation_contract": {
            "primary_within_cohort_gene_gene_correlation": "Pearson",
            "primary_rationale": (
                "Pearson correlation is computed after prespecified log-scale handling "
                "and within-cohort gene standardization, preserving magnitude-sensitive "
                "coexpression structure and compatibility with PCA/loadings."
            ),
            "prespecified_sensitivity": "Spearman",
            "sensitivity_rationale": (
                "Spearman checks robustness to remaining monotone scale/platform effects."
            ),
            "edge_vector_cross_cohort_concordance": "Spearman",
        },
        "scanb_profiles": role_info,
        "still_to_freeze": [
            "exact SCAN-B technical-replicate pairing",
            "deterministic choice of the primary member of each replicate group",
            "TCGA source-only gene filter",
            "TCGA WGCNA/module construction parameters",
            "source PC1 orientation rule",
            "assessability threshold",
            "coherence-matched null details",
            "permutation/bootstrapping counts",
        ],
    }

    json_out = OUT_DIR / "expression_input_contract_v1.json"
    json_out.write_text(json.dumps(contract, indent=2), encoding="utf-8")

    md = f"""# Paper 4 / TCBB — Frozen Expression-Input Contract v1

Status: **frozen for expression representation only**

No target outcomes, clinical characteristics, preservation statistics, or target correlation matrices were used.

## Primary representations

- **TCGA-BRCA:** `log2(RSEM + 1)` from the continuous cBioPortal RSEM matrix.
- **SCAN-B GSE96058:** published values are treated as `log2(FPKM + 0.1)`, deterministically back-transformed to FPKM, then represented as `log2(FPKM + 1)` for the primary matched-modality analysis.
- **METABRIC:** published continuous log2 Illumina HT-12 v3 intensity values are used as provided.
- After cohort-specific representation, each evaluable gene is standardized across samples within its cohort.

## Correlation contract

- Primary within-cohort gene–gene correlation: **Pearson**.
- Prespecified sensitivity: **Spearman**.
- Cross-cohort concordance of edge vectors: **Spearman**.

## SCAN-B sample-role status

- Total profiles: {role_info['total_profiles']}
- Apparent primary profiles: {role_info['apparent_primary_profiles']}
- Apparent technical replicates: {role_info['apparent_technical_replicates']}

Exact replicate pairing and the deterministic primary-member rule remain to be frozen before any preservation analysis.
"""
    md_out = OUT_DIR / "expression_input_contract_v1.md"
    md_out.write_text(md, encoding="utf-8")

    print("\n" + "=" * 116)
    print("01b EXPRESSION-INPUT CONTRACT: PASS")
    print("=" * 116)
    print("Frozen primary expression handling:")
    print("  TCGA-BRCA : log2(RSEM + 1)")
    print("  SCAN-B    : log2(max(2**published - 0.1, 0) + 1)")
    print("  METABRIC  : published continuous log2 microarray intensity")
    print("  Per gene  : within-cohort z-standardization after the above representation")
    print()
    print("Frozen correlation handling:")
    print("  Primary within-cohort gene-gene correlation: Pearson")
    print("  Prespecified sensitivity:                    Spearman")
    print("  Cross-cohort edge-vector concordance:         Spearman")
    print()
    print("SCAN-B exact replicate pairing / primary-member rule: NOT YET FROZEN")
    print("No preservation statistic was calculated.")
    print()
    print("Outputs:")
    print(f"  {json_out}")
    print(f"  {md_out}")
    print(f"  {role_out}")
    print(f"  {repl_examples_out}")
    print("=" * 116)


if __name__ == "__main__":
    main()
