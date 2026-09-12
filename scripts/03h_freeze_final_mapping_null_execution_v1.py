from __future__ import annotations

import json
from pathlib import Path

SCRIPT_VERSION = "03h-freeze-final-mapping-null-execution-v1-no-cli"

DATA_ROOT = Path(r"D:\paper4_tcbb_data")
BASE_NULL_CONTRACT = (
    DATA_ROOT
    / "paper4_tcbb_mapping_null_contract_v1"
    / "structure_preserving_mapping_null_contract_v1.json"
)
PILOT_AUDIT = (
    DATA_ROOT
    / "paper4_tcbb_mapping_null_pilot_audit_v1"
    / "mapping_null_pilot_audit_v1.json"
)
OUT_DIR = DATA_ROOT / "paper4_tcbb_final_mapping_null_execution_v1"

FINAL_VALID_PANELS = 1000
MAX_ATTEMPTS = 10000
MAX_SPECIFICITY_EDGES = 20000
EDGE_SUBSET_BASE_SEED = 20260914
NULL_MAPPING_BASE_SEED = 20260915


def require(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Required file does not exist: {path}")


def main() -> None:
    print("=" * 126)
    print("Paper 4 / TCBB - freeze FINAL execution details for the matched-mapping specificity null")
    print("=" * 126)
    print(f"Script version: {SCRIPT_VERSION}")
    print()
    print("Scientific guard:")
    print("  Target outcomes loaded:                              NO")
    print("  Target clinical characteristics loaded:              NO")
    print("  Target gene-gene correlations calculated:            NO")
    print("  Target PCA/loadings calculated:                       NO")
    print("  Target preservation statistics calculated:           NO")
    print("=" * 126)

    for p in [BASE_NULL_CONTRACT, PILOT_AUDIT]:
        require(p)

    base = json.loads(BASE_NULL_CONTRACT.read_text(encoding="utf-8"))
    pilot = json.loads(PILOT_AUDIT.read_text(encoding="utf-8"))

    if not bool(pilot.get("overall_pass", False)):
        raise RuntimeError(
            "03g pilot did not pass; final mapping-null execution cannot be frozen."
        )

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    contract = {
        "contract_id": "paper4-tcbb-final-mapping-null-execution-v1",
        "script_version": SCRIPT_VERSION,
        "status": "FROZEN_BEFORE_FIRST_TARGET_CORRELATION_OR_PRESERVATION",
        "base_generator_contract": str(BASE_NULL_CONTRACT),
        "pilot_audit": str(PILOT_AUDIT),
        "pilot_result": "PASS_ALL_PRIMARY_ASSESSABLE_COMBINATIONS",
        "scientific_guard": {
            "target_outcomes_loaded": False,
            "target_clinical_characteristics_loaded": False,
            "target_gene_gene_correlations_calculated": False,
            "target_pca_calculated": False,
            "target_preservation_statistics_calculated": False,
        },
        "final_edge_specificity_null": {
            "valid_panels_required": FINAL_VALID_PANELS,
            "maximum_panel_attempts": MAX_ATTEMPTS,
            "mapping_generator": "exactly the frozen 03d generator validated by 03g",
            "edge_statistic": "Spearman concordance of source and target edge vectors",
            "observed_specificity_statistic": (
                "computed on exactly the same fixed edge subset used by every null mapping"
            ),
            "maximum_edges_per_target_program": MAX_SPECIFICITY_EDGES,
            "edge_subset_rule": (
                "if the assessable program has <=20,000 nonredundant edges, use all; "
                "otherwise sample exactly 20,000 unordered source-gene pairs uniformly "
                "without replacement once, using a deterministic target/program seed, "
                "BEFORE any target edge value is inspected; reuse that identical pair set "
                "for the observed specificity statistic and all 1,000 null mappings"
            ),
            "edge_subset_seed_base": EDGE_SUBSET_BASE_SEED,
            "null_mapping_seed_base": NULL_MAPPING_BASE_SEED,
            "empirical_p": "(1 + count(|rho_null| >= |rho_observed_specificity|)) / 1001",
            "reason_for_edge_subsample": (
                "keeps the repeated 1,000-panel specificity analysis computationally "
                "tractable for modules with >1 million edges while retaining a large, "
                "predeclared and unbiased edge sample"
            ),
        },
        "primary_direct_edge_effect": {
            "use_all_edges": True,
            "note": (
                "The reported direct edge-preservation effect size uses ALL evaluable edges. "
                "Only the repeated specificity-null statistic is computed on the fixed "
                "up-to-20,000-edge audit subset."
            ),
        },
        "loading_specificity_amendment": {
            "random_mapping_loading_null_primary": False,
            "direct_loading_label_permutation_retained": True,
            "reason": (
                "A 1,000-panel random-mapping loading null would require refitting target "
                "PC1 thousands of times for modules up to ~1,700 genes. It is computationally "
                "disproportionate and partly redundant with the already frozen direct loading "
                "gene-label permutation test. Therefore mapping-null specificity is defined "
                "for the structural edge axis; loading remains an independent direct "
                "preservation axis with its own permutation inference."
            ),
            "timing": (
                "This amendment is frozen before any target correlation matrix, target PCA, "
                "loading preservation, edge preservation, or target outcome was calculated."
            ),
        },
        "matching_balance_diagnostics": {
            "report_for_final_panels": [
                "fraction exact 2D-stratum assignments",
                "fraction assignments at Manhattan distance <=1",
                "maximum stratum distance",
                "mean stratum distance",
                "mean absolute difference in source-mean percentile",
                "mean absolute difference in source-MAD percentile",
                "mean absolute difference in target-mean percentile",
                "mean absolute difference in target-MAD percentile",
            ],
            "use_for_posthoc_panel_rejection": False,
            "note": (
                "The frozen 03d validity rule controls panel acceptance. Continuous balance "
                "diagnostics are descriptive and may not be used to tune or discard panels "
                "after preservation results are seen."
            ),
        },
        "multiplicity": {
            "edge_mapping_null_p_values": (
                "BH-adjust all primary-assessable edge mapping-null p-values within each target"
            ),
            "loading_mapping_null_p_values": "not applicable",
        },
        "amendment_provenance": {
            "supersedes": (
                "only the loading-null execution and repeated-edge computation details "
                "in 03d; the 03d mapping generator/marginal matching/quality guard remain unchanged"
            ),
            "reason": (
                "computational scalability and removal of redundant loading-null refits; "
                "decision made before any target preservation statistic was calculated"
            ),
        },
    }

    json_out = OUT_DIR / "final_mapping_null_execution_contract_v1.json"
    md_out = OUT_DIR / "final_mapping_null_execution_contract_v1.md"

    json_out.write_text(json.dumps(contract, indent=2), encoding="utf-8")
    md_out.write_text(
        f"""# Paper 4 / TCBB — Final Mapping-Null Execution Contract v1

Status: **FROZEN BEFORE FIRST TARGET CORRELATION/PRESERVATION**

The 03g pilot passed for every primary-assessable target/program combination
with 100/100 valid pilot mappings, so the 03d mapping generator is retained
unchanged.

## Edge specificity

- Final valid mappings per target/program: **{FINAL_VALID_PANELS}**
- Maximum attempts: **{MAX_ATTEMPTS}**
- Direct edge-preservation effect size: **all evaluable edges**
- Repeated mapping-null specificity statistic:
  - all edges when there are <= {MAX_SPECIFICITY_EDGES:,};
  - otherwise one deterministic uniform sample of exactly
    **{MAX_SPECIFICITY_EDGES:,}** unordered source-gene pairs;
  - the same pair set is used for the correct mapping and all null mappings.
- Empirical p-value uses the +1 correction.

## Loading axis

The 1,000-panel random-mapping loading null is **not primary**. Loading
preservation retains the already frozen direct gene-label permutation test.
This amendment is made before any target PCA or preservation statistic and
avoids thousands of redundant high-dimensional PC1 refits.

## Matching diagnostics

Continuous source/target mean/MAD percentile balance is reported
descriptively for the final mappings. It cannot be used for post-hoc panel
rejection or tuning.
""",
        encoding="utf-8",
    )

    print("=" * 126)
    print("03h FINAL MAPPING-NULL EXECUTION CONTRACT: PASS")
    print("=" * 126)
    print(f"Final valid edge-null mappings: {FINAL_VALID_PANELS}")
    print(f"Maximum specificity edges:      {MAX_SPECIFICITY_EDGES:,}")
    print("Direct edge effect:              ALL evaluable edges")
    print("Mapping-null edge statistic:     fixed <=20,000-edge subset")
    print("Random-mapping loading null:     NOT primary")
    print("Direct loading permutation:      retained")
    print()
    print("No target correlation/PCA/preservation statistic was calculated.")
    print("Outputs:")
    print(f"  {json_out}")
    print(f"  {md_out}")
    print("=" * 126)


if __name__ == "__main__":
    main()
