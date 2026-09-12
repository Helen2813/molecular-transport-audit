from __future__ import annotations

import json
from pathlib import Path


SCRIPT_VERSION = "03d-freeze-structure-preserving-mapping-null-contract-v1-no-cli"

DATA_ROOT = Path(r"D:\paper4_tcbb_data")
OUT_DIR = DATA_ROOT / "paper4_tcbb_mapping_null_contract_v1"

FINAL_NULL_PANELS = 1000
PILOT_PANELS = 100
MAX_PANEL_ATTEMPTS_FACTOR = 10

N_BINS_SOURCE_COMPOSITE = 5
N_BINS_TARGET_COMPOSITE = 5

VALID_PANEL_MIN_FRACTION_WITHIN_DISTANCE_1 = 0.90
VALID_PANEL_MAX_STRATUM_DISTANCE = 2

BASE_SEED = 20260913


def main() -> None:
    print("=" * 124)
    print("Paper 4 / TCBB - freeze source-structure-preserving matched-mapping specificity-null contract")
    print("=" * 124)
    print(f"Script version: {SCRIPT_VERSION}")
    print()
    print("Scientific guard:")
    print("  Target outcomes loaded:                              NO")
    print("  Target clinical characteristics loaded:              NO")
    print("  Target preservation statistics calculated:           NO")
    print("  Final null statistics calculated by this script:     NO")
    print("=" * 124)

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    contract = {
        "contract_id": "paper4-tcbb-structure-preserving-mapping-null-v1",
        "script_version": SCRIPT_VERSION,
        "status": "FROZEN_BEFORE_TARGET_MARGINAL_METRICS_AND_PRESERVATION",
        "motivation": (
            "The predeclared 03c source-only feasibility audit showed that ordinary "
            "same-size random panels cannot match the source coherence of 11/12 frozen "
            "TCGA modules in either target. The final specificity null therefore keeps "
            "the tested source program structure exactly fixed and randomizes only the "
            "source-to-target gene mapping."
        ),
        "scientific_guard": {
            "target_outcomes_loaded": False,
            "target_clinical_characteristics_loaded": False,
            "target_preservation_statistics_calculated": False,
        },
        "null_object": {
            "source_program_membership": "fixed",
            "source_edge_vector": "fixed exactly on the target-evaluable frozen source genes",
            "source_loading_vector": "fixed exactly on the target-evaluable frozen source genes",
            "source_loading_signs": "fixed",
            "randomized_component": (
                "one-to-one assignment of unrelated target-available background genes "
                "to the frozen source-gene slots"
            ),
            "interpretation": (
                "tests whether preservation depends on the correct biological gene mapping "
                "rather than merely on the source module's size/coherence"
            ),
        },
        "candidate_pool": {
            "base": (
                "frozen TCGA 10,000-gene source universe intersected with the target's "
                "frozen mapped-symbol universe"
            ),
            "exclusion": "all genes belonging to the tested frozen source program",
            "replacement": False,
        },
        "marginal_matching": {
            "source_features": [
                "source mean log2(RSEM+1)",
                "source MAD log2(RSEM+1)",
            ],
            "target_features": [
                "target mean on frozen transformed scale",
                "target MAD on frozen transformed scale",
            ],
            "percentile_ranks": (
                "computed separately within each target-available frozen 10k universe"
            ),
            "source_composite_score": (
                "mean(percentile_rank(source mean), percentile_rank(source MAD))"
            ),
            "target_composite_score": (
                "mean(percentile_rank(target mean), percentile_rank(target MAD))"
            ),
            "strata": (
                f"{N_BINS_SOURCE_COMPOSITE} source-composite quantile bins x "
                f"{N_BINS_TARGET_COMPOSITE} target-composite quantile bins"
            ),
            "assignment_rule": (
                "for each null panel, randomize source-slot order; for each slot sample "
                "uniformly without replacement from unused candidates in the same 2D stratum; "
                "if empty, expand to strata in increasing Manhattan distance; ties between "
                "equidistant strata are randomized"
            ),
        },
        "panel_quality_guard": {
            "minimum_fraction_slots_with_stratum_distance_le_1":
                VALID_PANEL_MIN_FRACTION_WITHIN_DISTANCE_1,
            "maximum_allowed_stratum_distance": VALID_PANEL_MAX_STRATUM_DISTANCE,
            "invalid_panel_action": "discard panel and generate a new attempt",
            "not_estimable_rule": (
                f"if {FINAL_NULL_PANELS} valid panels cannot be obtained within "
                f"{MAX_PANEL_ATTEMPTS_FACTOR}x that many attempts, specificity for that "
                "target/program is labelled not estimable"
            ),
        },
        "panel_counts": {
            "pilot_match_quality_panels": PILOT_PANELS,
            "final_specificity_panels": FINAL_NULL_PANELS,
        },
        "statistics": {
            "edge_null": (
                "Spearman concordance between the fixed source edge vector and the target "
                "edge vector of randomly mapped background genes in source-slot order"
            ),
            "loading_null": (
                "Spearman concordance between fixed source loadings and target PC1 loadings "
                "of randomly mapped genes, with target PC1 oriented by the frozen source-sign score"
            ),
            "empirical_p": "(1 + count(|null| >= |observed|)) / (K + 1)",
            "multiplicity": (
                "within each target, all assessable mapping-null edge and loading p-values "
                "form one Benjamini-Hochberg family"
            ),
        },
        "seed": {
            "base_seed": BASE_SEED,
            "schedule": "deterministic target/program/panel offsets from base seed",
        },
        "important_distinction": (
            "This null is distinct from the direct within-program gene-label permutation test. "
            "The direct test asks whether the correct correspondence within the measured program "
            "matters; this specificity null asks whether the exact source program maps unusually "
            "well to its true target genes compared with unrelated but marginally comparable genes."
        ),
    }

    json_out = OUT_DIR / "structure_preserving_mapping_null_contract_v1.json"
    md_out = OUT_DIR / "structure_preserving_mapping_null_contract_v1.md"

    json_out.write_text(json.dumps(contract, indent=2), encoding="utf-8")

    md_out.write_text(
        f"""# Paper 4 / TCBB — Structure-Preserving Matched-Mapping Null v1

Status: **FROZEN BEFORE TARGET MARGINAL METRICS OR PRESERVATION**

The 03c source-only audit showed that ordinary same-size random panels do not
match the source coherence of 11/12 frozen TCGA modules. The specificity null
therefore keeps each tested source program **exactly fixed** and randomizes only
its mapping to unrelated target-available genes.

## Fixed source object

For every assessable target/program:

- frozen source membership is fixed;
- the source edge vector is fixed;
- source PC1 loadings and signs are fixed;
- only the target-gene identity assigned to each source slot is randomized.

## Marginal matching

Background genes are matched using source and target measurability, not target
preservation:

- source composite = average percentile rank of source mean and source MAD;
- target composite = average percentile rank of target mean and target MAD;
- a {N_BINS_SOURCE_COMPOSITE} x {N_BINS_TARGET_COMPOSITE} quantile grid is used.

Assignments use the same stratum first and expand by Manhattan distance only
when needed, always without replacement.

## Quality guard

A null panel is valid only if:

- at least {VALID_PANEL_MIN_FRACTION_WITHIN_DISTANCE_1:.0%} of assigned slots have stratum distance <= 1; and
- no assigned slot has stratum distance > {VALID_PANEL_MAX_STRATUM_DISTANCE}.

Invalid panels are discarded. If {FINAL_NULL_PANELS} valid panels cannot be
generated within {MAX_PANEL_ATTEMPTS_FACTOR * FINAL_NULL_PANELS} attempts,
specificity is **not estimable** for that target/program.

## Final null size

- pilot matching-quality audit: {PILOT_PANELS} panels;
- final specificity analysis: {FINAL_NULL_PANELS} valid panels.

No target outcome or preservation result may modify this contract.
""",
        encoding="utf-8",
    )

    print("=" * 124)
    print("03d STRUCTURE-PRESERVING MAPPING NULL CONTRACT: PASS")
    print("=" * 124)
    print("Source edges/loadings: fixed exactly")
    print("Randomized component: target-gene mapping only")
    print(
        f"Marginal matching: {N_BINS_SOURCE_COMPOSITE}x"
        f"{N_BINS_TARGET_COMPOSITE} source/target composite strata"
    )
    print(
        f"Valid panel guard: >= {VALID_PANEL_MIN_FRACTION_WITHIN_DISTANCE_1:.0%} "
        f"slots within distance <=1; max distance <= {VALID_PANEL_MAX_STRATUM_DISTANCE}"
    )
    print(f"Pilot panels: {PILOT_PANELS}")
    print(f"Final panels: {FINAL_NULL_PANELS}")
    print()
    print("No target marginal metric or preservation statistic was calculated.")
    print("Outputs:")
    print(f"  {json_out}")
    print(f"  {md_out}")
    print("=" * 124)


if __name__ == "__main__":
    main()
