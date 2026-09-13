from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd


SCRIPT_VERSION = "04i-freeze-target-statistical-evaluability-amendment-v1-no-cli"

DATA_ROOT = Path(r"D:\paper4_tcbb_data")

PRETARGET_AUDIT = (
    DATA_ROOT
    / "paper4_tcbb_final_pretarget_audits_v1"
    / "final_pretarget_audit_v1.json"
)
SOURCE_UNIVERSE = (
    DATA_ROOT
    / "paper4_tcbb_tcga_source_universe_v1"
    / "tcga_source_gene_universe_frozen_v1.tsv"
)
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
METABRIC_EXPR = (
    DATA_ROOT
    / "paper4_tcbb_input_audit_v1"
    / "staged_continuous_inputs"
    / "METABRIC"
    / "data_mrna_illumina_microarray.txt"
)

OUT_DIR = DATA_ROOT / "paper4_tcbb_target_statistical_evaluability_amendment_v1"

EXPECTED_IDENTIFIER_MAPPED = {
    "SCANB_GSE96058": 9225,
    "METABRIC": 8490,
}
EXPECTED_NONZERO_SD_FROM_04H = {
    "SCANB_GSE96058": 9225,
    "METABRIC": 8485,
}

MIN_EVALUABLE_GENES = 30
MIN_COVERAGE = 0.80


def require(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Required file does not exist: {path}")


def canon_symbol(x: object) -> str:
    if x is None:
        return ""
    return " ".join(str(x).strip().split()).upper()


def load_scanb_target_stats(
    universe: set[str],
    primary_titles: list[str],
) -> pd.DataFrame:
    parts = []

    for chunk in pd.read_csv(
        SCANB_EXPR,
        compression="gzip",
        chunksize=1500,
        low_memory=False,
    ):
        symbol_col = chunk.columns[0]
        symbols = chunk[symbol_col].map(canon_symbol)
        mask = symbols.isin(universe)
        if not mask.any():
            continue

        missing_cols = [c for c in primary_titles if c not in chunk.columns]
        if missing_cols:
            raise RuntimeError(
                f"SCAN-B primary columns missing from expression matrix: {missing_cols[:20]}"
            )

        sub = chunk.loc[mask, [symbol_col] + primary_titles].copy()
        arr = sub[primary_titles].to_numpy(dtype=np.float64)

        # Frozen transform: published log2(FPKM+0.1) -> FPKM -> log2(FPKM+1).
        fpkm = np.maximum(np.exp2(arr) - 0.1, 0.0)
        transformed = np.log2(fpkm + 1.0)

        out = pd.DataFrame(transformed, columns=primary_titles)
        out.insert(0, "Hugo_Symbol", symbols.loc[mask].to_numpy())
        parts.append(out)

    if not parts:
        raise RuntimeError("No frozen-universe genes found in SCAN-B.")

    df = pd.concat(parts, ignore_index=True)
    # Frozen duplicate rule: average duplicate symbols on transformed scale.
    df = df.groupby("Hugo_Symbol", sort=False, as_index=True).mean(numeric_only=True)

    x = df.to_numpy(dtype=np.float64)
    finite = np.isfinite(x).all(axis=1)
    sd = np.full(len(df), np.nan, dtype=np.float64)
    if finite.any():
        sd[finite] = np.std(x[finite], axis=1, ddof=1)

    return pd.DataFrame(
        {
            "Hugo_Symbol": df.index.to_numpy(),
            "all_finite": finite.astype(int),
            "target_sample_sd": sd,
            "statistically_evaluable": (finite & np.isfinite(sd) & (sd > 0)).astype(int),
        }
    )


def load_metabric_target_stats(universe: set[str]) -> pd.DataFrame:
    parts = []
    sample_cols = None

    for chunk in pd.read_csv(
        METABRIC_EXPR,
        sep="\t",
        chunksize=1500,
        low_memory=False,
    ):
        if "Hugo_Symbol" not in chunk.columns:
            raise RuntimeError("METABRIC expression chunk lacks Hugo_Symbol.")
        if sample_cols is None:
            sample_cols = list(chunk.columns[2:])

        symbols = chunk["Hugo_Symbol"].map(canon_symbol)
        mask = symbols.isin(universe)
        if not mask.any():
            continue

        arr = chunk.loc[mask, sample_cols].to_numpy(dtype=np.float64)
        out = pd.DataFrame(arr, columns=sample_cols)
        out.insert(0, "Hugo_Symbol", symbols.loc[mask].to_numpy())
        parts.append(out)

    if not parts:
        raise RuntimeError("No frozen-universe genes found in METABRIC.")

    df = pd.concat(parts, ignore_index=True)
    # Frozen duplicate rule: average duplicate symbols on as-provided log2 scale.
    df = df.groupby("Hugo_Symbol", sort=False, as_index=True).mean(numeric_only=True)

    x = df.to_numpy(dtype=np.float64)
    finite = np.isfinite(x).all(axis=1)
    sd = np.full(len(df), np.nan, dtype=np.float64)
    if finite.any():
        sd[finite] = np.std(x[finite], axis=1, ddof=1)

    return pd.DataFrame(
        {
            "Hugo_Symbol": df.index.to_numpy(),
            "all_finite": finite.astype(int),
            "target_sample_sd": sd,
            "statistically_evaluable": (finite & np.isfinite(sd) & (sd > 0)).astype(int),
        }
    )


def program_impact(
    target_name: str,
    target_stats: pd.DataFrame,
    weights: pd.DataFrame,
) -> pd.DataFrame:
    represented = set(target_stats["Hugo_Symbol"])
    eligible = set(
        target_stats.loc[
            target_stats["statistically_evaluable"].astype(int) == 1,
            "Hugo_Symbol",
        ]
    )

    rows = []
    for program_id in sorted(weights["program_id"].unique()):
        genes = [
            canon_symbol(x)
            for x in weights.loc[weights["program_id"] == program_id, "Hugo_Symbol"]
        ]
        n_frozen = len(genes)
        n_mapped = sum(g in represented for g in genes)
        n_eligible = sum(g in eligible for g in genes)

        mapped_coverage = n_mapped / n_frozen
        eligible_coverage = n_eligible / n_frozen
        primary_assessable = (
            n_eligible >= MIN_EVALUABLE_GENES
            and eligible_coverage >= MIN_COVERAGE
        )

        rows.append(
            {
                "target": target_name,
                "program_id": program_id,
                "n_frozen_genes": n_frozen,
                "n_identifier_mapped": n_mapped,
                "identifier_coverage": mapped_coverage,
                "n_statistically_evaluable": n_eligible,
                "statistical_coverage": eligible_coverage,
                "n_mapped_but_ineligible": n_mapped - n_eligible,
                "primary_assessable_after_amendment": int(primary_assessable),
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    print("=" * 136)
    print("Paper 4 / TCBB - freeze target statistical-evaluability amendment before first preservation")
    print("=" * 136)
    print(f"Script version: {SCRIPT_VERSION}")
    print()
    print("Reason for amendment:")
    print("  04h confirmed that identifier mapping and statistical evaluability are not identical in METABRIC.")
    print("  Correlation/PCA require finite, nonzero-variance target genes.")
    print()
    print("Scientific guard:")
    print("  Target expression values read:                       YES (variance/evaluability only)")
    print("  Target outcomes/treatment loaded:                    NO")
    print("  Target gene-gene preservation calculated:            NO")
    print("  Target PCA/loading preservation calculated:           NO")
    print("  Mapping-null preservation calculated:                 NO")
    print("  Primary classifier evaluated:                         NO")
    print("=" * 136)

    for p in [
        PRETARGET_AUDIT,
        SOURCE_UNIVERSE,
        WEIGHTS,
        SCANB_EXPR,
        SCANB_PRIMARY,
        METABRIC_EXPR,
    ]:
        require(p)

    pret = json.loads(PRETARGET_AUDIT.read_text(encoding="utf-8"))
    if pret.get("status") != "PASS_READY_FOR_FIRST_POOLED_TARGET_PRESERVATION":
        raise RuntimeError(
            f"04h did not finish in the expected PASS state: {pret.get('status')}"
        )

    universe_df = pd.read_csv(SOURCE_UNIVERSE, sep="\t", dtype=str).fillna("")
    universe = {canon_symbol(x) for x in universe_df["Hugo_Symbol"]}
    if len(universe) != 10_000:
        raise RuntimeError(f"Expected 10,000 unique frozen source genes; found {len(universe)}.")

    weights = pd.read_csv(WEIGHTS, sep="\t", dtype=str).fillna("")
    weights["Hugo_Symbol"] = weights["Hugo_Symbol"].map(canon_symbol)

    primary = pd.read_csv(SCANB_PRIMARY, sep="\t", dtype=str).fillna("")
    primary_titles = list(primary["primary_title"])
    if len(primary_titles) != 3273:
        raise RuntimeError(f"Expected 3,273 frozen SCAN-B primary profiles; found {len(primary_titles)}.")

    print("\n[1/3] SCAN-B statistical evaluability ...")
    scanb = load_scanb_target_stats(universe)
    scanb_mapped = len(scanb)
    scanb_eval = int(scanb["statistically_evaluable"].sum())
    print(f"  identifier-mapped frozen-universe genes: {scanb_mapped:,}")
    print(f"  finite + nonzero-SD evaluable genes:      {scanb_eval:,}")

    print("\n[2/3] METABRIC statistical evaluability ...")
    met = load_metabric_target_stats(universe)
    met_mapped = len(met)
    met_eval = int(met["statistically_evaluable"].sum())
    print(f"  identifier-mapped frozen-universe genes: {met_mapped:,}")
    print(f"  finite + nonzero-SD evaluable genes:      {met_eval:,}")

    observed = {
        "SCANB_GSE96058": (scanb_mapped, scanb_eval),
        "METABRIC": (met_mapped, met_eval),
    }
    for target, (mapped, evaluable) in observed.items():
        if mapped != EXPECTED_IDENTIFIER_MAPPED[target]:
            raise RuntimeError(
                f"{target}: mapped count {mapped} != expected {EXPECTED_IDENTIFIER_MAPPED[target]}"
            )
        if evaluable != EXPECTED_NONZERO_SD_FROM_04H[target]:
            raise RuntimeError(
                f"{target}: evaluable count {evaluable} != 04h count "
                f"{EXPECTED_NONZERO_SD_FROM_04H[target]}"
            )

    scanb_ineligible = scanb.loc[
        scanb["statistically_evaluable"].astype(int) == 0
    ].copy()
    scanb_ineligible.insert(0, "target", "SCANB_GSE96058")

    met_ineligible = met.loc[
        met["statistically_evaluable"].astype(int) == 0
    ].copy()
    met_ineligible.insert(0, "target", "METABRIC")

    affected = pd.concat(
        [scanb_ineligible, met_ineligible],
        ignore_index=True,
    )

    print("\nMapped but statistically ineligible genes:")
    if affected.empty:
        print("  NONE")
    else:
        program_by_gene = defaultdict(list)
        for row in weights.itertuples(index=False):
            program_by_gene[canon_symbol(row.Hugo_Symbol)].append(row.program_id)

        for row in affected.itertuples(index=False):
            progs = program_by_gene.get(row.Hugo_Symbol, [])
            print(
                f"  {row.target}: {row.Hugo_Symbol} "
                f"SD={row.target_sample_sd:.12g}; programs={progs if progs else ['grey/not-tested']}"
            )

    print("\n[3/3] Recomputing the frozen assessability guard using statistically evaluable genes ...")
    impact = pd.concat(
        [
            program_impact("SCANB_GSE96058", scanb, weights),
            program_impact("METABRIC", met, weights),
        ],
        ignore_index=True,
    )

    for target in ["SCANB_GSE96058", "METABRIC"]:
        sub = impact.loc[impact["target"] == target]
        print(f"  {target}:")
        for r in sub.itertuples(index=False):
            print(
                f"    {r.program_id}: mapped={r.n_identifier_mapped:4d}/{r.n_frozen_genes:4d} "
                f"eligible={r.n_statistically_evaluable:4d}/{r.n_frozen_genes:4d} "
                f"coverage={r.statistical_coverage:.3f} "
                f"assessable={'YES' if r.primary_assessable_after_amendment else 'NO'}"
            )

    scanb_assess = int(
        impact.loc[
            impact["target"] == "SCANB_GSE96058",
            "primary_assessable_after_amendment",
        ].sum()
    )
    met_assess = int(
        impact.loc[
            impact["target"] == "METABRIC",
            "primary_assessable_after_amendment",
        ].sum()
    )

    if scanb_assess != 12:
        raise RuntimeError(f"Unexpected SCAN-B assessable-program count after amendment: {scanb_assess}")
    if met_assess != 10:
        raise RuntimeError(f"Unexpected METABRIC assessable-program count after amendment: {met_assess}")

    contract = {
        "contract_id": "paper4-tcbb-target-statistical-evaluability-amendment-v1",
        "script_version": SCRIPT_VERSION,
        "status": "FROZEN_BEFORE_FIRST_TARGET_PRESERVATION_STATISTIC",
        "trigger": (
            "04h observed 8,485 finite/nonzero-SD frozen-universe genes in METABRIC "
            "versus 8,490 identifier-mapped genes from the pre-expression mapping audit."
        ),
        "scientific_guard": {
            "target_outcomes_loaded": False,
            "target_treatment_loaded": False,
            "target_gene_gene_preservation_calculated": False,
            "target_pca_loading_preservation_calculated": False,
            "mapping_null_preservation_calculated": False,
            "primary_classifier_evaluated": False,
        },
        "statistical_evaluability_rule": {
            "identifier_mapping": "retain the frozen exact canonical-symbol mapping and duplicate averaging rule",
            "target_gene_required": (
                "mapped target gene must have finite transformed values across the analysis sample set "
                "and strictly positive sample standard deviation on that transformed scale"
            ),
            "reason": (
                "Pearson correlation, z-standardization, and PCA are undefined for zero-variance genes"
            ),
            "program_coverage_numerator": "number of statistically evaluable frozen program genes",
            "program_coverage_denominator": "all genes in the frozen source program",
            "assessability_guard": (
                f"n_evaluable >= {MIN_EVALUABLE_GENES} AND "
                f"n_evaluable/n_frozen >= {MIN_COVERAGE:.2f}"
            ),
            "source_edge_and_loading_restriction": (
                "for a target/program comparison, restrict the frozen source edge/loading object "
                "to exactly the same statistically evaluable frozen genes used in the target"
            ),
        },
        "mapping_null_amendment": {
            "candidate_pool": (
                "the frozen 10,000-gene source universe intersected with target-available genes "
                "that are statistically evaluable under the rule above, excluding all tested-program genes"
            ),
            "observed_slots": (
                "only statistically evaluable frozen program genes enter the observed specificity statistic"
            ),
            "matching_features": (
                "retain the already frozen source/target mean-MAD percentile features and 5x5 matching strata "
                "for genes that remain eligible"
            ),
            "generator": "otherwise unchanged from frozen 03d/03h",
            "panel_validity": "unchanged from frozen 03d",
            "attempt_cap": "unchanged: 10,000 attempts for 1,000 valid mappings",
            "reason_no_posthoc_rescue": (
                "zero-variance genes are removed because the planned correlation statistic is mathematically "
                "undefined, not because of any observed preservation result"
            ),
        },
        "observed_counts": {
            "SCANB_GSE96058": {
                "identifier_mapped": scanb_mapped,
                "statistically_evaluable": scanb_eval,
                "mapped_but_ineligible": scanb_mapped - scanb_eval,
                "primary_assessable_programs": scanb_assess,
            },
            "METABRIC": {
                "identifier_mapped": met_mapped,
                "statistically_evaluable": met_eval,
                "mapped_but_ineligible": met_mapped - met_eval,
                "primary_assessable_programs": met_assess,
            },
        },
        "primary_assessability_changed_from_03a": False,
        "note": (
            "The amendment changes only mathematical gene eligibility. It does not alter any frozen "
            "module membership, coverage threshold, classifier threshold, permutation rule, or target outcome firewall."
        ),
    }

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    scanb.to_csv(
        OUT_DIR / "scanb_target_statistical_evaluability_v1.tsv",
        sep="\t",
        index=False,
    )
    met.to_csv(
        OUT_DIR / "metabric_target_statistical_evaluability_v1.tsv",
        sep="\t",
        index=False,
    )
    affected.to_csv(
        OUT_DIR / "mapped_but_statistically_ineligible_genes_v1.tsv",
        sep="\t",
        index=False,
    )
    impact.to_csv(
        OUT_DIR / "program_assessability_after_statistical_evaluability_v1.tsv",
        sep="\t",
        index=False,
    )

    json_out = OUT_DIR / "target_statistical_evaluability_amendment_v1.json"
    json_out.write_text(json.dumps(contract, indent=2), encoding="utf-8")

    print("\n" + "=" * 136)
    print("04i TARGET STATISTICAL-EVALUABILITY AMENDMENT: PASS")
    print("=" * 136)
    print(f"SCAN-B:   {scanb_eval:,}/{scanb_mapped:,} mapped genes statistically evaluable; 12/12 programs assessable")
    print(f"METABRIC: {met_eval:,}/{met_mapped:,} mapped genes statistically evaluable; 10/12 programs assessable")
    print("Primary assessability changed from 03a: NO")
    print("Preservation statistics seen:             NO")
    print("Mapping-null generator:                    unchanged except ineligible genes cannot enter")
    print()
    print("Outputs:")
    print(f"  {json_out}")
    print(f"  {OUT_DIR / 'mapped_but_statistically_ineligible_genes_v1.tsv'}")
    print(f"  {OUT_DIR / 'program_assessability_after_statistical_evaluability_v1.tsv'}")
    print("=" * 136)


if __name__ == "__main__":
    main()
