from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np
import pandas as pd


SCRIPT_VERSION = "02b-freeze-tcga-source-universe-v1-no-cli"

DATA_ROOT = Path(r"D:\paper4_tcbb_data")
AUDIT_ROOT = DATA_ROOT / "paper4_tcbb_tcga_source_filter_audit_v1"
INPUT_ROOT = (
    DATA_ROOT
    / "paper4_tcbb_input_audit_v1"
    / "staged_continuous_inputs"
    / "TCGA_BRCA_PanCanAtlas2018"
)

TCGA_EXPR = INPUT_ROOT / "data_mrna_seq_v2_rsem.txt"
GENE_METRICS = AUDIT_ROOT / "tcga_source_gene_metrics_v1.tsv"

OUT_DIR = DATA_ROOT / "paper4_tcbb_tcga_source_universe_v1"

# ---------------------------------------------------------------------
# FROZEN SOURCE-ONLY GENE-UNIVERSE RULE
# ---------------------------------------------------------------------
MIN_FRAC_RSEM_GE1 = 0.20
MAX_MISSING_FRACTION = 0.05
N_TOP_BY_MAD = 10_000

# FROZEN SOURCE SAMPLE RULE
PRIMARY_TUMOR_SAMPLE_TYPE = "01"


def require(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Required file does not exist: {path}")


def banner() -> None:
    print("=" * 120)
    print("Paper 4 / TCBB - freeze TCGA-BRCA source sample set and source gene universe")
    print("=" * 120)
    print(f"Script version: {SCRIPT_VERSION}")
    print()
    print("Scientific guard:")
    print("  SCAN-B opened:                                      NO")
    print("  METABRIC opened:                                    NO")
    print("  GSE239948 opened:                                   NO")
    print("  Clinical outcomes loaded:                           NO")
    print("  Target coverage used for gene selection:            NO")
    print("  Target preservation statistics calculated:          NO")
    print("  WGCNA/modules calculated:                           NO")
    print("=" * 120)


def parse_tcga_sample_id(sample_id: str) -> tuple[str, str]:
    """
    Parse canonical TCGA sample IDs of the form TCGA-XX-YYYY-ZZ...
    Returns (patient_id, two-digit sample-type code).
    """
    s = str(sample_id).strip()
    parts = s.split("-")
    if len(parts) < 4 or parts[0].upper() != "TCGA":
        raise ValueError(f"Unrecognized TCGA sample ID: {s}")
    patient_id = "-".join(parts[:3])
    m = re.match(r"^(\d{2})", parts[3])
    if not m:
        raise ValueError(f"Cannot parse TCGA sample-type code from: {s}")
    return patient_id, m.group(1)


def build_sample_manifest(sample_cols: list[str]) -> tuple[pd.DataFrame, dict]:
    rows = []
    for sid in sample_cols:
        patient_id, sample_type = parse_tcga_sample_id(sid)
        rows.append(
            {
                "sample_id": sid,
                "patient_id": patient_id,
                "sample_type_code": sample_type,
            }
        )
    df = pd.DataFrame(rows)

    type_counts = {
        str(k): int(v)
        for k, v in df["sample_type_code"].value_counts().sort_index().to_dict().items()
    }

    # Frozen design: primary solid tumour (TCGA code 01) only.
    primary = df.loc[df["sample_type_code"] == PRIMARY_TUMOR_SAMPLE_TYPE].copy()
    nonprimary = df.loc[df["sample_type_code"] != PRIMARY_TUMOR_SAMPLE_TYPE].copy()

    # If multiple source samples exist for one patient, keep a deterministic
    # single profile: lexicographically first sample ID. No target information
    # and no outcome/quality criterion is used.
    primary = primary.sort_values(["patient_id", "sample_id"], kind="mergesort")
    primary["patient_sample_rank"] = (
        primary.groupby("patient_id").cumcount() + 1
    )
    primary["include_in_source_network"] = (
        primary["patient_sample_rank"] == 1
    ).astype(int)

    selected = primary.loc[primary["include_in_source_network"] == 1].copy()

    info = {
        "expression_columns_total": int(len(df)),
        "sample_type_counts": type_counts,
        "primary_type01_columns": int(len(primary)),
        "nonprimary_columns_excluded": int(len(nonprimary)),
        "unique_primary_patients": int(primary["patient_id"].nunique()),
        "duplicate_primary_profiles_beyond_first": int(
            (primary["patient_sample_rank"] > 1).sum()
        ),
        "frozen_selected_source_profiles": int(len(selected)),
        "frozen_rule": (
            "retain TCGA sample-type code 01 only; if more than one code-01 "
            "profile exists for a patient, retain the lexicographically first sample ID"
        ),
    }
    return df.merge(
        primary[
            ["sample_id", "patient_sample_rank", "include_in_source_network"]
        ],
        on="sample_id",
        how="left",
    ).fillna(
        {
            "patient_sample_rank": 0,
            "include_in_source_network": 0,
        }
    ), info


def numeric_entrez(x: str) -> float:
    try:
        return float(x)
    except Exception:
        return np.inf


def main() -> None:
    banner()
    require(TCGA_EXPR)
    require(GENE_METRICS)
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    header = pd.read_csv(TCGA_EXPR, sep="\t", nrows=0, low_memory=False)
    cols = list(header.columns)
    if cols[:2] != ["Hugo_Symbol", "Entrez_Gene_Id"]:
        raise RuntimeError(f"Unexpected TCGA identifier columns: {cols[:2]}")
    sample_cols = cols[2:]

    print("\n[1/3] Freezing source sample set from TCGA sample IDs only ...")
    sample_manifest, sample_info = build_sample_manifest(sample_cols)
    selected_samples = sample_manifest.loc[
        sample_manifest["include_in_source_network"].astype(int) == 1,
        "sample_id",
    ].tolist()

    print(f"  expression columns:        {sample_info['expression_columns_total']:,}")
    print(f"  sample type counts:        {sample_info['sample_type_counts']}")
    print(f"  primary type-01 columns:   {sample_info['primary_type01_columns']:,}")
    print(f"  duplicate primary profiles beyond first: "
          f"{sample_info['duplicate_primary_profiles_beyond_first']:,}")
    print(f"  frozen source profiles:    {sample_info['frozen_selected_source_profiles']:,}")

    print("\n[2/3] Freezing source gene universe ...")
    metrics = pd.read_csv(GENE_METRICS, sep="\t", dtype={"Hugo_Symbol": str, "Entrez_Gene_Id": str})
    metrics = metrics.reset_index(drop=True)
    metrics["source_row_index_0based"] = np.arange(len(metrics), dtype=int)
    metrics["Hugo_Symbol"] = metrics["Hugo_Symbol"].fillna("").str.strip()
    metrics["Entrez_Gene_Id"] = metrics["Entrez_Gene_Id"].fillna("").str.strip()

    if len(metrics) != 20531:
        raise RuntimeError(f"Expected 20,531 source gene rows; found {len(metrics)}")

    measurable = (
        (metrics["missing_fraction"] <= MAX_MISSING_FRACTION)
        & (metrics["frac_rsem_ge1"] >= MIN_FRAC_RSEM_GE1)
        & (metrics["mad_log2_rsem1"] > 0)
        & metrics["Hugo_Symbol"].ne("")
    )
    candidates = metrics.loc[measurable].copy()

    # Deterministic duplicate-symbol resolution:
    # highest MAD, then variance, then measured fraction, then smaller numeric
    # Entrez ID, then earlier source-row index.
    candidates["_entrez_numeric"] = candidates["Entrez_Gene_Id"].map(numeric_entrez)
    candidates = candidates.sort_values(
        [
            "Hugo_Symbol",
            "mad_log2_rsem1",
            "var_log2_rsem1",
            "frac_rsem_ge1",
            "_entrez_numeric",
            "source_row_index_0based",
        ],
        ascending=[True, False, False, False, True, True],
        kind="mergesort",
    )
    before_dedup = len(candidates)
    candidates = candidates.drop_duplicates(subset=["Hugo_Symbol"], keep="first").copy()
    duplicate_rows_removed = before_dedup - len(candidates)

    # Frozen robust-variability ranking and fixed cap.
    candidates = candidates.sort_values(
        [
            "mad_log2_rsem1",
            "var_log2_rsem1",
            "Hugo_Symbol",
            "source_row_index_0based",
        ],
        ascending=[False, False, True, True],
        kind="mergesort",
    )

    if len(candidates) < N_TOP_BY_MAD:
        raise RuntimeError(
            f"Only {len(candidates)} unique measurable genes remain; "
            f"cannot freeze top {N_TOP_BY_MAD}."
        )

    frozen = candidates.iloc[:N_TOP_BY_MAD].copy()
    frozen.insert(0, "source_gene_rank_by_MAD", np.arange(1, len(frozen) + 1))
    frozen["include_in_source_network"] = 1
    frozen["frozen_filter_rule"] = (
        "missing<=0.05; RSEM>=1 in >=20% source samples; MAD(log2(RSEM+1))>0; "
        "nonblank Hugo symbol; duplicate symbols resolved source-only; "
        "top 10000 unique genes by MAD"
    )
    frozen = frozen.drop(columns=["_entrez_numeric"])

    threshold_mad = float(frozen["mad_log2_rsem1"].iloc[-1])

    print(f"  measurable nonblank rows before symbol dedup: {before_dedup:,}")
    print(f"  duplicate-symbol rows removed:                {duplicate_rows_removed:,}")
    print(f"  unique measurable symbols:                    {len(candidates):,}")
    print(f"  frozen top-by-MAD universe:                    {len(frozen):,}")
    print(f"  MAD at rank 10,000:                            {threshold_mad:.6g}")

    print("\n[3/3] Writing frozen manifests ...")
    sample_out = OUT_DIR / "tcga_source_samples_frozen_v1.tsv"
    gene_out = OUT_DIR / "tcga_source_gene_universe_frozen_v1.tsv"
    sample_manifest.to_csv(sample_out, sep="\t", index=False)
    frozen.to_csv(gene_out, sep="\t", index=False)

    contract = {
        "contract_id": "paper4-tcbb-tcga-source-universe-v1",
        "script_version": SCRIPT_VERSION,
        "status": "FROZEN",
        "scientific_guard": {
            "scanb_opened": False,
            "metabric_opened": False,
            "gse239948_opened": False,
            "outcomes_loaded": False,
            "target_coverage_used_for_selection": False,
            "target_preservation_statistics_calculated": False,
            "wgcna_modules_calculated": False,
        },
        "source_samples": sample_info,
        "source_gene_universe": {
            "starting_gene_rows": int(len(metrics)),
            "max_missing_fraction": MAX_MISSING_FRACTION,
            "minimum_fraction_source_samples_with_RSEM_ge_1": MIN_FRAC_RSEM_GE1,
            "require_positive_MAD_log2_RSEM_plus_1": True,
            "require_nonblank_Hugo_symbol": True,
            "duplicate_symbol_resolution": (
                "highest source MAD; then highest source variance; then highest "
                "fraction RSEM>=1; then smallest numeric Entrez ID; then earliest source row"
            ),
            "unique_measurable_symbols_after_dedup": int(len(candidates)),
            "final_selection": f"top {N_TOP_BY_MAD} by source MAD",
            "frozen_gene_count": int(len(frozen)),
            "MAD_threshold_at_final_rank": threshold_mad,
            "random_panel_background_rule": (
                "primary MTA random/control panels for TCGA-derived modules must "
                "draw from this same frozen 10,000-gene source universe, after "
                "the separately frozen target-availability intersection and exclusions"
            ),
        },
        "next_step": (
            "Perform TCGA-source-only soft-threshold/network diagnostics using only "
            "the frozen source sample manifest and frozen 10,000-gene universe."
        ),
    }

    json_out = OUT_DIR / "tcga_source_universe_contract_v1.json"
    json_out.write_text(json.dumps(contract, indent=2), encoding="utf-8")

    md = f"""# Paper 4 / TCBB — Frozen TCGA-BRCA Source Universe v1

Status: **FROZEN**

No SCAN-B, METABRIC, GSE239948, target coverage, target preservation result, or
clinical outcome was used.

## Source samples

Frozen rule:
- retain TCGA sample-type code **01** (primary solid tumour);
- if more than one code-01 profile exists for a patient, retain the
  lexicographically first sample ID;
- no outcome or expression-quality ranking is used to choose among profiles.

Selected source profiles: **{sample_info['frozen_selected_source_profiles']}**

## Source gene universe

A row is eligible if:
- missing fraction <= **{MAX_MISSING_FRACTION}**;
- RSEM >= 1 in at least **{MIN_FRAC_RSEM_GE1:.0%}** of source samples;
- MAD of `log2(RSEM+1)` > 0;
- Hugo symbol is nonblank.

Duplicate Hugo symbols are resolved using source-only statistics:
highest MAD, then highest variance, then highest measured fraction, then
smallest numeric Entrez ID, then earliest source row.

From the unique eligible symbols, retain the **top {N_TOP_BY_MAD} genes by source MAD**.

The final gene at rank {N_TOP_BY_MAD} has MAD **{threshold_mad:.6g}**.

This same frozen source universe is the starting candidate universe for
TCGA-derived matched random/control panels; the null must not draw unrelated
genes from outside the network-discovery universe.
"""
    md_out = OUT_DIR / "tcga_source_universe_contract_v1.md"
    md_out.write_text(md, encoding="utf-8")

    print("\n" + "=" * 120)
    print("02b TCGA SOURCE UNIVERSE CONTRACT: PASS")
    print("=" * 120)
    print(f"Frozen source profiles: {sample_info['frozen_selected_source_profiles']:,}")
    print(f"Frozen source genes:    {len(frozen):,}")
    print(f"MAD threshold rank 10,000: {threshold_mad:.6g}")
    print()
    print("Random-panel background is now constrained to this frozen source universe.")
    print("No target dataset was opened.")
    print("No WGCNA modules or preservation statistics were calculated.")
    print()
    print("Outputs:")
    print(f"  {sample_out}")
    print(f"  {gene_out}")
    print(f"  {json_out}")
    print(f"  {md_out}")
    print("=" * 120)


if __name__ == "__main__":
    main()
