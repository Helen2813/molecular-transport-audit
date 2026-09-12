from __future__ import annotations

import json
import re
from pathlib import Path

import pandas as pd


SCRIPT_VERSION = "01c-freeze-scanb-replicate-pairing-v1-no-cli"

DATA_ROOT = Path(r"D:\paper4_tcbb_data")
AUDIT_ROOT = DATA_ROOT / "paper4_tcbb_input_audit_v1"

TECH_INVENTORY = AUDIT_ROOT / "scanb_technical_inventory_v1.tsv"
SCANB_EXPR = (
    DATA_ROOT
    / "SCANB_GSE96058"
    / "GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv.gz"
)

OUT_DIR = DATA_ROOT / "paper4_tcbb_scanb_pairing_contract_v1"


def require(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Required file does not exist: {path}")


def banner() -> None:
    print("=" * 118)
    print("Paper 4 / TCBB - freeze SCAN-B technical-replicate pairing and primary-profile rule")
    print("=" * 118)
    print(f"Script version: {SCRIPT_VERSION}")
    print()
    print("Scientific guard:")
    print("  Target outcomes loaded:                              NO")
    print("  Target clinical characteristics loaded:              NO")
    print("  Target preservation statistics calculated:           NO")
    print("  Target correlation matrices calculated:              NO")
    print("  Expression values read:                              NO")
    print("  Only sample titles/accessions/platforms + expression header are used")
    print("=" * 118)


def normalize_title(x: str) -> str:
    return str(x).strip()


def is_replicate_title(title: str) -> bool:
    return bool(re.search(r"repl$", title, flags=re.IGNORECASE))


def primary_title_from_replicate(title: str) -> str:
    if not is_replicate_title(title):
        raise ValueError(f"Not a terminal-'repl' title: {title}")
    return re.sub(r"repl$", "", title, flags=re.IGNORECASE)


def main() -> None:
    banner()
    require(TECH_INVENTORY)
    require(SCANB_EXPR)

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    tech = pd.read_csv(TECH_INVENTORY, sep="\t", dtype=str).fillna("")
    required_cols = {
        "geo_accession",
        "sample_title",
        "platform_id",
        "source_series_matrix",
    }
    missing = required_cols.difference(tech.columns)
    if missing:
        raise RuntimeError(f"Technical inventory missing columns: {sorted(missing)}")

    tech["sample_title"] = tech["sample_title"].map(normalize_title)
    tech["geo_accession"] = tech["geo_accession"].str.strip()
    tech["platform_id"] = tech["platform_id"].str.strip()

    if tech["geo_accession"].duplicated().any():
        dup = tech.loc[tech["geo_accession"].duplicated(keep=False), "geo_accession"].tolist()
        raise RuntimeError(f"Duplicate GEO accessions in technical inventory: {dup[:20]}")

    # Read only the CSV header; no expression values are opened.
    expr_header = pd.read_csv(SCANB_EXPR, compression="gzip", nrows=0).columns.tolist()
    expr_sample_titles = [str(x).strip() for x in expr_header[1:]]
    expr_title_set = set(expr_sample_titles)

    if len(expr_sample_titles) != 3409:
        raise RuntimeError(
            f"Expected 3409 SCAN-B expression profiles; found {len(expr_sample_titles)}"
        )
    if len(expr_title_set) != len(expr_sample_titles):
        raise RuntimeError("Duplicate sample-title columns in SCAN-B expression matrix.")

    inventory_title_set = set(tech["sample_title"])
    only_expr = sorted(expr_title_set - inventory_title_set)
    only_inventory = sorted(inventory_title_set - expr_title_set)
    if only_expr or only_inventory:
        raise RuntimeError(
            "Expression-header titles do not match technical inventory. "
            f"Only in expression={only_expr[:10]}, only in inventory={only_inventory[:10]}"
        )

    tech["is_technical_replicate"] = tech["sample_title"].map(is_replicate_title)

    replicates = tech.loc[tech["is_technical_replicate"]].copy()
    primaries = tech.loc[~tech["is_technical_replicate"]].copy()

    if len(replicates) != 136:
        raise RuntimeError(f"Expected 136 terminal-'repl' profiles; found {len(replicates)}")
    if len(primaries) != 3273:
        raise RuntimeError(f"Expected 3273 primary profiles; found {len(primaries)}")

    if primaries["sample_title"].duplicated().any():
        d = primaries.loc[
            primaries["sample_title"].duplicated(keep=False), "sample_title"
        ].tolist()
        raise RuntimeError(f"Duplicate primary titles: {d[:20]}")

    pair_rows = []
    unmatched = []
    ambiguous = []

    for _, rep in replicates.iterrows():
        rep_title = rep["sample_title"]
        primary_title = primary_title_from_replicate(rep_title)

        matches = primaries.loc[primaries["sample_title"] == primary_title]
        if len(matches) == 0:
            unmatched.append((rep_title, primary_title))
            continue
        if len(matches) != 1:
            ambiguous.append((rep_title, primary_title, len(matches)))
            continue

        pri = matches.iloc[0]

        pair_rows.append(
            {
                "pair_id": primary_title,
                "primary_title": primary_title,
                "primary_geo_accession": pri["geo_accession"],
                "primary_platform_id": pri["platform_id"],
                "replicate_title": rep_title,
                "replicate_geo_accession": rep["geo_accession"],
                "replicate_platform_id": rep["platform_id"],
                "same_platform": int(pri["platform_id"] == rep["platform_id"]),
                "cross_platform": int(pri["platform_id"] != rep["platform_id"]),
                "primary_selection_rule": "title_without_terminal_repl",
                "replicate_pairing_rule": "strip_terminal_repl_from_replicate_title",
            }
        )

    if unmatched:
        raise RuntimeError(f"Unmatched replicate titles: {unmatched[:20]}")
    if ambiguous:
        raise RuntimeError(f"Ambiguous replicate titles: {ambiguous[:20]}")
    if len(pair_rows) != 136:
        raise RuntimeError(f"Expected 136 complete pairs; built {len(pair_rows)}")

    pairs = pd.DataFrame(pair_rows)

    if pairs["pair_id"].duplicated().any():
        d = pairs.loc[pairs["pair_id"].duplicated(keep=False), "pair_id"].tolist()
        raise RuntimeError(f"A primary title is paired more than once: {d[:20]}")

    primary_manifest = primaries[
        ["sample_title", "geo_accession", "platform_id", "source_series_matrix"]
    ].copy()
    primary_manifest = primary_manifest.rename(
        columns={
            "sample_title": "primary_title",
            "geo_accession": "primary_geo_accession",
            "platform_id": "primary_platform_id",
        }
    )
    primary_manifest["include_in_primary_scanb_analysis"] = 1
    primary_manifest["selection_rule"] = "exclude_titles_ending_case_insensitive_terminal_repl"

    primary_title_set = set(primary_manifest["primary_title"])
    excluded_title_set = set(replicates["sample_title"])
    if len(primary_title_set) != 3273 or len(excluded_title_set) != 136:
        raise RuntimeError("Primary/replicate set cardinality validation failed.")
    if primary_title_set.intersection(excluded_title_set):
        raise RuntimeError("Primary and replicate title sets overlap.")
    if primary_title_set.union(excluded_title_set) != expr_title_set:
        raise RuntimeError("Primary + replicate titles do not exactly partition expression columns.")

    pair_out = OUT_DIR / "scanb_technical_replicate_pairs_v1.tsv"
    primary_out = OUT_DIR / "scanb_primary_profiles_frozen_v1.tsv"
    pairs.to_csv(pair_out, sep="\t", index=False)
    primary_manifest.to_csv(primary_out, sep="\t", index=False)

    combo = (
        pairs.groupby(
            ["primary_platform_id", "replicate_platform_id", "same_platform", "cross_platform"],
            dropna=False,
        )
        .size()
        .reset_index(name="n_pairs")
        .sort_values(["primary_platform_id", "replicate_platform_id"])
    )
    combo_out = OUT_DIR / "scanb_replicate_platform_pair_counts_v1.tsv"
    combo.to_csv(combo_out, sep="\t", index=False)

    contract = {
        "contract_id": "paper4-tcbb-scanb-replicate-pairing-contract-v1",
        "script_version": SCRIPT_VERSION,
        "status": "FROZEN",
        "scientific_guard": {
            "target_outcomes_loaded": False,
            "target_clinical_characteristics_loaded": False,
            "target_preservation_statistics_calculated": False,
            "target_correlation_matrices_calculated": False,
            "expression_values_read": False,
        },
        "source_basis": {
            "official_GEO_naming_rule": (
                "primary sample titles have base names such as F30; "
                "technical replicate titles append terminal 'repl', e.g. F30repl"
            ),
            "local_validation": {
                "expression_profiles": len(expr_sample_titles),
                "primary_profiles": len(primary_manifest),
                "technical_replicate_profiles": len(replicates),
                "complete_pairs": len(pairs),
            },
        },
        "frozen_rules": {
            "technical_replicate_definition": (
                "sample_title ends with terminal 'repl' (case-insensitive)"
            ),
            "pairing_rule": (
                "remove terminal 'repl' from replicate title and require exactly one "
                "non-replicate title with the resulting base title"
            ),
            "primary_profile_rule": (
                "for the main SCAN-B biological validation include the non-'repl' profile; "
                "exclude the corresponding 'repl' profile"
            ),
            "quality_based_profile_selection": False,
            "technical_replicate_use": (
                "excluded from primary biological target; retained only for "
                "technical-repeatability / measurement-ceiling analyses"
            ),
        },
        "platform_pair_counts": combo.to_dict(orient="records"),
        "remaining_to_freeze": [
            "technical-repeatability subsampling schedule at matched n",
            "same-platform versus cross-platform ceiling reporting details",
        ],
    }

    json_out = OUT_DIR / "scanb_pairing_contract_v1.json"
    json_out.write_text(json.dumps(contract, indent=2), encoding="utf-8")

    md = f"""# Paper 4 / TCBB — SCAN-B Replicate Pairing Contract v1

Status: **FROZEN**

No target outcomes, clinical characteristics, expression values, correlation matrices,
or preservation statistics were used.

## Frozen naming and pairing rule

The official GSE96058 convention identifies a primary sample by a base title such as
`F30` and its technical replicate by appending `repl`, e.g. `F30repl`.

Primary analysis rule:

- include every non-`repl` SCAN-B profile;
- exclude every title ending in terminal `repl` (case-insensitive);
- pair each replicate by stripping the terminal `repl`;
- require exactly one corresponding non-replicate primary title;
- never choose between pair members using quality metrics or preservation results.

## Local validation

- Expression profiles: **{len(expr_sample_titles)}**
- Primary profiles: **{len(primary_manifest)}**
- Technical replicate profiles: **{len(replicates)}**
- Complete technical pairs: **{len(pairs)}**

The replicate set is reserved for technical-repeatability analyses and is not part of the
primary TCGA-BRCA → SCAN-B biological transport cohort.
"""
    md_out = OUT_DIR / "scanb_pairing_contract_v1.md"
    md_out.write_text(md, encoding="utf-8")

    print("\n" + "=" * 118)
    print("01c SCAN-B PAIRING CONTRACT: PASS")
    print("=" * 118)
    print(f"Expression profiles:             {len(expr_sample_titles)}")
    print(f"Frozen primary SCAN-B profiles:  {len(primary_manifest)}")
    print(f"Technical replicate profiles:    {len(replicates)}")
    print(f"Complete replicate pairs:        {len(pairs)}")
    print()
    print("Platform combinations:")
    for _, row in combo.iterrows():
        print(
            f"  primary={row['primary_platform_id']} -> "
            f"replicate={row['replicate_platform_id']}: "
            f"{int(row['n_pairs'])} pairs"
        )
    print()
    print("Frozen primary-member rule: NON-'repl' title")
    print("No quality-based selection was used.")
    print("No expression values or preservation statistics were read/calculated.")
    print()
    print("Outputs:")
    print(f"  {pair_out}")
    print(f"  {primary_out}")
    print(f"  {combo_out}")
    print(f"  {json_out}")
    print(f"  {md_out}")
    print("=" * 118)


if __name__ == "__main__":
    main()
