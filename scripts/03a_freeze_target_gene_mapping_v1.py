from __future__ import annotations

import json
import re
from pathlib import Path

import pandas as pd


SCRIPT_VERSION = "03a-freeze-target-gene-mapping-v1-no-cli"

DATA_ROOT = Path(r"D:\paper4_tcbb_data")

SOURCE_WEIGHTS = (
    DATA_ROOT
    / "paper4_tcbb_frozen_source_programs_v1"
    / "tcga_frozen_source_program_weights_v1.tsv"
)
SOURCE_UNIVERSE = (
    DATA_ROOT
    / "paper4_tcbb_tcga_source_universe_v1"
    / "tcga_source_gene_universe_frozen_v1.tsv"
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

OUT_DIR = DATA_ROOT / "paper4_tcbb_target_mapping_v1"

MIN_EVALUABLE_GENES = 30
MIN_TARGET_COVERAGE = 0.80


def require(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Required file does not exist: {path}")


def normalize_symbol(x: object) -> str:
    if pd.isna(x):
        return ""
    s = str(x).strip()
    s = re.sub(r"\s+", " ", s)
    return s.upper()


def banner() -> None:
    print("=" * 122)
    print("Paper 4 / TCBB - freeze target gene-symbol mapping and coverage BEFORE preservation analysis")
    print("=" * 122)
    print(f"Script version: {SCRIPT_VERSION}")
    print()
    print("Scientific guard:")
    print("  Target clinical outcomes loaded:                     NO")
    print("  Target clinical characteristics loaded:              NO")
    print("  Target expression VALUES loaded:                     NO")
    print("  Target correlations calculated:                      NO")
    print("  Target preservation statistics calculated:           NO")
    print("  Target identifiers/gene coverage inspected:          YES (allowed input QC)")
    print("=" * 122)


def target_identifier_table_scanb(path: Path) -> pd.DataFrame:
    df = pd.read_csv(
        path,
        compression="gzip",
        usecols=[0],
        dtype=str,
        low_memory=False,
    )
    col = df.columns[0]
    out = pd.DataFrame(
        {
            "target_row_index_0based": range(len(df)),
            "target_original_symbol": df[col].fillna("").astype(str),
        }
    )
    out["target_normalized_symbol"] = out["target_original_symbol"].map(normalize_symbol)
    return out


def target_identifier_table_metabric(path: Path) -> pd.DataFrame:
    df = pd.read_csv(
        path,
        sep="\t",
        usecols=["Hugo_Symbol", "Entrez_Gene_Id"],
        dtype=str,
        low_memory=False,
    )
    out = pd.DataFrame(
        {
            "target_row_index_0based": range(len(df)),
            "target_original_symbol": df["Hugo_Symbol"].fillna("").astype(str),
            "target_entrez_gene_id": df["Entrez_Gene_Id"].fillna("").astype(str),
        }
    )
    out["target_normalized_symbol"] = out["target_original_symbol"].map(normalize_symbol)
    return out


def summarize_identifier_table(df: pd.DataFrame) -> dict:
    nonblank = df["target_normalized_symbol"].ne("")
    dup_nonblank = nonblank & df["target_normalized_symbol"].duplicated(keep=False)
    counts = df.loc[nonblank, "target_normalized_symbol"].value_counts()
    return {
        "rows": int(len(df)),
        "blank_normalized_symbols": int((~nonblank).sum()),
        "unique_nonblank_normalized_symbols": int(counts.size),
        "rows_in_duplicate_normalized_symbol_groups": int(dup_nonblank.sum()),
        "duplicate_normalized_symbols": int((counts > 1).sum()),
        "maximum_rows_per_normalized_symbol": int(counts.max()) if len(counts) else 0,
    }


def build_unique_mapping(df: pd.DataFrame) -> pd.DataFrame:
    nonblank = df.loc[df["target_normalized_symbol"].ne("")].copy()

    grouped = (
        nonblank.groupby("target_normalized_symbol", sort=True)
        .agg(
            target_row_indices_0based=(
                "target_row_index_0based",
                lambda x: ",".join(str(int(v)) for v in x),
            ),
            target_original_symbols=(
                "target_original_symbol",
                lambda x: "|".join(sorted(set(str(v) for v in x))),
            ),
            target_row_count=("target_row_index_0based", "size"),
        )
        .reset_index()
    )
    grouped["duplicate_expression_rule"] = grouped["target_row_count"].map(
        lambda n: (
            "single_row"
            if int(n) == 1
            else "arithmetic_mean_of_duplicate_rows_on_frozen_transformed_scale"
        )
    )
    return grouped


def coverage_table(
    weights: pd.DataFrame,
    mapping_symbols: set[str],
    target_name: str,
) -> pd.DataFrame:
    rows = []
    for pid, grp in weights.groupby("program_id", sort=True):
        genes = [normalize_symbol(x) for x in grp["Hugo_Symbol"]]
        if len(genes) != len(set(genes)):
            raise RuntimeError(f"Frozen source program {pid} has duplicate normalized symbols.")

        evaluable = [g for g in genes if g in mapping_symbols]
        n_frozen = len(genes)
        n_eval = len(evaluable)
        cov = n_eval / n_frozen if n_frozen else 0.0
        assessable = n_eval >= MIN_EVALUABLE_GENES and cov >= MIN_TARGET_COVERAGE

        rows.append(
            {
                "target": target_name,
                "program_id": pid,
                "frozen_genes": n_frozen,
                "evaluable_genes": n_eval,
                "coverage": cov,
                "minimum_evaluable_genes_rule": MIN_EVALUABLE_GENES,
                "minimum_coverage_rule": MIN_TARGET_COVERAGE,
                "primary_assessable": int(assessable),
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    banner()

    for p in [SOURCE_WEIGHTS, SOURCE_UNIVERSE, SCANB_EXPR, METABRIC_EXPR]:
        require(p)
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print("\n[1/4] Loading frozen source identifiers ...")
    weights = pd.read_csv(SOURCE_WEIGHTS, sep="\t", dtype=str)
    universe = pd.read_csv(SOURCE_UNIVERSE, sep="\t", dtype=str)

    weights["source_normalized_symbol"] = weights["Hugo_Symbol"].map(normalize_symbol)
    universe["source_normalized_symbol"] = universe["Hugo_Symbol"].map(normalize_symbol)

    if weights["source_normalized_symbol"].eq("").any():
        raise RuntimeError("Blank normalized source-program symbol.")
    if universe["source_normalized_symbol"].eq("").any():
        raise RuntimeError("Blank normalized source-universe symbol.")
    if universe["source_normalized_symbol"].duplicated().any():
        raise RuntimeError("Frozen 10,000-gene source universe is not unique after normalization.")

    source_universe_symbols = set(universe["source_normalized_symbol"])
    program_symbols = set(weights["source_normalized_symbol"])

    print(f"  frozen source universe: {len(source_universe_symbols):,} symbols")
    print(f"  frozen non-grey program genes: {len(program_symbols):,} symbols")

    print("\n[2/4] Reading TARGET IDENTIFIERS ONLY ...")
    scanb_ids = target_identifier_table_scanb(SCANB_EXPR)
    metabric_ids = target_identifier_table_metabric(METABRIC_EXPR)

    scanb_info = summarize_identifier_table(scanb_ids)
    metabric_info = summarize_identifier_table(metabric_ids)

    print(
        "  SCAN-B:   "
        f"rows={scanb_info['rows']:,}, "
        f"unique symbols={scanb_info['unique_nonblank_normalized_symbols']:,}, "
        f"duplicate symbols={scanb_info['duplicate_normalized_symbols']:,}"
    )
    print(
        "  METABRIC: "
        f"rows={metabric_info['rows']:,}, "
        f"unique symbols={metabric_info['unique_nonblank_normalized_symbols']:,}, "
        f"duplicate symbols={metabric_info['duplicate_normalized_symbols']:,}"
    )

    scanb_map = build_unique_mapping(scanb_ids)
    metabric_map = build_unique_mapping(metabric_ids)

    scanb_symbols = set(scanb_map["target_normalized_symbol"])
    metabric_symbols = set(metabric_map["target_normalized_symbol"])

    scanb_shared_universe = source_universe_symbols & scanb_symbols
    metabric_shared_universe = source_universe_symbols & metabric_symbols

    print("\n[3/4] Computing identifier coverage only ...")
    print(
        f"  frozen 10k universe represented in SCAN-B:   "
        f"{len(scanb_shared_universe):,} / 10,000 "
        f"({len(scanb_shared_universe)/10000:.1%})"
    )
    print(
        f"  frozen 10k universe represented in METABRIC: "
        f"{len(metabric_shared_universe):,} / 10,000 "
        f"({len(metabric_shared_universe)/10000:.1%})"
    )

    scanb_cov = coverage_table(weights, scanb_symbols, "SCANB_GSE96058")
    metabric_cov = coverage_table(weights, metabric_symbols, "METABRIC")
    coverage = pd.concat([scanb_cov, metabric_cov], ignore_index=True)

    for target, grp in coverage.groupby("target", sort=False):
        print(f"\n  {target}:")
        for _, row in grp.iterrows():
            print(
                f"    {row['program_id']:<10} "
                f"{int(row['evaluable_genes']):>4}/{int(row['frozen_genes']):<4} "
                f"coverage={row['coverage']:.3f} "
                f"assessable={'YES' if int(row['primary_assessable']) else 'NO'}"
            )

    print("\n[4/4] Freezing mapping contract ...")

    scanb_map_out = OUT_DIR / "scanb_gene_mapping_frozen_v1.tsv"
    metabric_map_out = OUT_DIR / "metabric_gene_mapping_frozen_v1.tsv"
    coverage_out = OUT_DIR / "target_program_gene_coverage_v1.tsv"

    scanb_map.to_csv(scanb_map_out, sep="\t", index=False)
    metabric_map.to_csv(metabric_map_out, sep="\t", index=False)
    coverage.to_csv(coverage_out, sep="\t", index=False)

    contract = {
        "contract_id": "paper4-tcbb-target-gene-mapping-v1",
        "script_version": SCRIPT_VERSION,
        "status": "FROZEN_BEFORE_TARGET_EXPRESSION_OR_PRESERVATION",
        "scientific_guard": {
            "target_outcomes_loaded": False,
            "target_clinical_characteristics_loaded": False,
            "target_expression_values_loaded": False,
            "target_correlations_calculated": False,
            "target_preservation_statistics_calculated": False,
        },
        "normalization": {
            "rule": "trim whitespace; collapse internal whitespace; uppercase",
            "alias_substitution": False,
            "fuzzy_matching": False,
            "synonym_rescue": False,
            "ortholog_mapping": False,
        },
        "duplicate_target_symbol_rule": (
            "if multiple target rows map to the same normalized symbol, "
            "later expression is the arithmetic mean of those rows on the already-frozen "
            "target transformed scale"
        ),
        "assessability": {
            "minimum_evaluable_genes": MIN_EVALUABLE_GENES,
            "minimum_target_coverage": MIN_TARGET_COVERAGE,
            "rule": "both thresholds must be met; otherwise label not assessable",
        },
        "SCANB_GSE96058": {
            "identifier_summary": scanb_info,
            "source_universe_shared_symbols": len(scanb_shared_universe),
        },
        "METABRIC": {
            "identifier_summary": metabric_info,
            "source_universe_shared_symbols": len(metabric_shared_universe),
        },
        "random_panel_candidate_universe_rule": (
            "For each target, matched random/control panels start from the frozen "
            "10,000-gene TCGA source universe intersected with that target's frozen "
            "mapped-symbol universe. The tested program's own genes are excluded "
            "from its panel candidate pool. Genes from other frozen modules are allowed, "
            "because nearly the entire source universe belongs to frozen modules."
        ),
    }

    json_out = OUT_DIR / "target_gene_mapping_contract_v1.json"
    json_out.write_text(json.dumps(contract, indent=2), encoding="utf-8")

    print("\n" + "=" * 122)
    print("03a TARGET GENE-MAPPING CONTRACT: PASS")
    print("=" * 122)
    print(
        f"SCAN-B shared frozen source universe:   {len(scanb_shared_universe):,}/10,000"
    )
    print(
        f"METABRIC shared frozen source universe: {len(metabric_shared_universe):,}/10,000"
    )
    print(
        f"SCAN-B assessable programs:   "
        f"{int(scanb_cov['primary_assessable'].sum())}/12"
    )
    print(
        f"METABRIC assessable programs: "
        f"{int(metabric_cov['primary_assessable'].sum())}/12"
    )
    print()
    print("No target expression values were loaded.")
    print("No target correlation/preservation statistic was calculated.")
    print()
    print("Outputs:")
    print(f"  {scanb_map_out}")
    print(f"  {metabric_map_out}")
    print(f"  {coverage_out}")
    print(f"  {json_out}")
    print("=" * 122)


if __name__ == "__main__":
    main()
