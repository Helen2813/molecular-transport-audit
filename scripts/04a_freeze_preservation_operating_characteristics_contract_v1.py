from __future__ import annotations

import json
from pathlib import Path


SCRIPT_VERSION = "04a-freeze-preservation-operating-characteristics-contract-v1-no-cli"

DATA_ROOT = Path(r"D:\paper4_tcbb_data")
OUT_DIR = DATA_ROOT / "paper4_tcbb_preservation_operating_contract_v1"

# -----------------------------------------------------------------------------
# DIRECT PRESERVATION INFERENCE
# -----------------------------------------------------------------------------
DIRECT_PERMUTATIONS = 10_000
DIRECT_ALPHA = 0.05

# -----------------------------------------------------------------------------
# TARGET-SAMPLING UNCERTAINTY
# -----------------------------------------------------------------------------
BOOTSTRAP_REPLICATES = 1_000
BOOTSTRAP_CI = 0.95

# -----------------------------------------------------------------------------
# RELIABILITY
# -----------------------------------------------------------------------------
SPLIT_HALF_REPEATS = 2_000

# -----------------------------------------------------------------------------
# CONTROLLED MAPPING DEGRADATION
# -----------------------------------------------------------------------------
DEGRADATION_FRACTIONS = [0.0, 0.10, 0.25, 0.50, 0.75, 1.0]
DEGRADATION_REPLICATES = 100
DEGRADATION_TREND_SPEARMAN_MAX = -0.80
FULL_CORRUPTION_MAX_FALSE_SPECIFIC_FRACTION = 0.10

# -----------------------------------------------------------------------------
# SCAN-B SAMPLE-SIZE OPERATING CHARACTERISTICS
# -----------------------------------------------------------------------------
SCANB_SUBSAMPLE_SIZES = [43, 100, 250, 500, 1000, 2000, 3273]
SCANB_SUBSAMPLE_REPEATS = 50

# -----------------------------------------------------------------------------
# TECHNICAL-REPLICATE CEILING
# -----------------------------------------------------------------------------
SAME_PLATFORM_TECH_PAIRS = 100
CROSS_PLATFORM_TECH_PAIRS = 36

BASE_SEED = 20260916


def main() -> None:
    print("=" * 128)
    print("Paper 4 / TCBB - freeze target preservation, uncertainty, degradation, and operating-characteristic contract")
    print("=" * 128)
    print(f"Script version: {SCRIPT_VERSION}")
    print()
    print("Scientific guard:")
    print("  Target outcomes loaded:                              NO")
    print("  Target clinical characteristics loaded:              NO")
    print("  Target edge/loading preservation calculated:          NO")
    print("  Target bootstrap/permutation results calculated:      NO")
    print("  Controlled degradation results calculated:            NO")
    print("=" * 128)

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    contract = {
        "contract_id": "paper4-tcbb-preservation-operating-characteristics-v1",
        "script_version": SCRIPT_VERSION,
        "status": "FROZEN_BEFORE_FIRST_TARGET_PRESERVATION_STATISTIC",
        "scientific_guard": {
            "target_outcomes_loaded": False,
            "target_clinical_characteristics_loaded": False,
            "target_preservation_statistics_calculated": False,
            "target_bootstrap_results_calculated": False,
            "controlled_degradation_results_calculated": False,
        },

        "direct_preservation": {
            "edge": {
                "within_cohort_gene_gene_correlation": "Pearson",
                "effect": (
                    "Spearman correlation between frozen source and target "
                    "nonredundant edge vectors, using all evaluable edges"
                ),
                "positive_direction": "rho_edge > 0",
                "permutation": (
                    "permute target gene labels relative to frozen source labels; "
                    "apply the same permutation to target correlation-matrix rows/columns"
                ),
            },
            "loading": {
                "target_representation": (
                    "target PC1 loadings on gene-standardized evaluable frozen genes"
                ),
                "orientation": (
                    "orient target PC1 using the previously frozen source-sign target score"
                ),
                "effect": (
                    "Spearman correlation between frozen source loadings and "
                    "oriented target PC1 loadings"
                ),
                "positive_direction": "rho_load > 0",
                "permutation": (
                    "permute oriented target PC1 loading values across frozen gene labels"
                ),
            },
            "monte_carlo_permutations": DIRECT_PERMUTATIONS,
            "permutation_p_value": "(1 + count(|null| >= |observed|)) / (B + 1)",
            "multiplicity": (
                "within each target, ALL primary-assessable direct edge and loading "
                "permutation p-values form one Benjamini-Hochberg family"
            ),
            "alpha": DIRECT_ALPHA,
        },

        "specificity": {
            "primary_axis": "edge",
            "null": (
                "the frozen 03d/03h structure-preserving matched target-mapping null"
            ),
            "valid_mappings": 1000,
            "multiplicity": (
                "within each target, all primary-assessable edge mapping-null "
                "p-values form one Benjamini-Hochberg family"
            ),
            "alpha": 0.05,
        },

        "primary_classification": {
            "not_assessable": (
                "apply the already frozen coverage/evaluable-gene guard before classification"
            ),
            "directionally_discordant": (
                "either direct edge or direct loading test has BH q<0.05 "
                "and the corresponding observed rho is negative; this takes precedence"
            ),
            "strong": (
                "positive direct edge q<0.05 AND positive direct loading q<0.05 "
                "AND edge mapping-specificity q<0.05"
            ),
            "partial": (
                "edge mapping-specificity q<0.05 AND exactly one positive direct axis "
                "has q<0.05, with no significantly negative direct axis"
            ),
            "no_clear": "all remaining primary-assessable cases",
            "important_change_from_v1": (
                "specificity is part of the primary class definition rather than a side annotation; "
                "the old universal split-half >=0.60 threshold is retired from classification"
            ),
        },

        "reliability": {
            "split_half_repeats": SPLIT_HALF_REPEATS,
            "split_definition": (
                "randomly partition evaluable frozen genes into two nonoverlapping halves; "
                "compute frozen-sign target scores for both halves"
            ),
            "raw_statistic": "Pearson correlation of the two half scores across target samples",
            "spearman_brown_corrected": "2*r/(1+r), reported per repeat when finite",
            "summary": "median, 5th, 95th percentiles for raw and corrected reliability",
            "classifier_threshold": None,
            "leave_one_out": (
                "retain minimum and median correlation of leave-one-gene-out frozen-sign "
                "score with the full frozen-sign score as descriptive diagnostics"
            ),
        },

        "target_sampling_uncertainty": {
            "bootstrap_replicates": BOOTSTRAP_REPLICATES,
            "resampled_unit": "target biological samples only",
            "source_object": "fixed; source membership, edges, and loadings are not bootstrapped",
            "edge_bootstrap": (
                "recompute target Pearson edge matrix and full-edge rho_edge on each bootstrap"
            ),
            "loading_bootstrap": (
                "recompute target PC1/orientation and rho_load on each bootstrap"
            ),
            "ci": f"{int(BOOTSTRAP_CI * 100)}% percentile interval",
            "reason_source_not_bootstrapped": (
                "the estimand is transport relative to a frozen source program; "
                "bootstrap uncertainty quantifies finite target-sample variability"
            ),
        },

        "controlled_mapping_degradation": {
            "role": "primary known-answer operating-characteristic analysis",
            "fractions": DEGRADATION_FRACTIONS,
            "replicates_per_nonzero_fraction": DEGRADATION_REPLICATES,
            "zero_fraction": "the unique correct frozen target mapping",
            "corrupted_slot_selection": (
                "choose exactly round(fraction * m) evaluable source-gene slots uniformly "
                "without replacement using a deterministic seed schedule"
            ),
            "replacement_rule": (
                "for corrupted slots replace the true target gene with an unrelated "
                "marginally matched background gene drawn by the frozen 03d mapping generator; "
                "keep source slot identity, source edge structure, source loading, and source sign fixed"
            ),
            "uncorrupted_slots": "retain the correct target gene mapping",
            "no_duplicate_target_gene_assignments": True,
            "primary_degradation_statistic": (
                "edge rho on the same fixed up-to-20,000-edge subset used by the 03h specificity audit"
            ),
            "formal_trend_gate_scope": (
                "apply only to target/program combinations whose uncorrupted edge result "
                "is positive and mapping-specificity-supported"
            ),
            "formal_trend_gate": (
                f"Spearman correlation between corruption fraction and median rho_edge "
                f"across fractions must be <= {DEGRADATION_TREND_SPEARMAN_MAX}"
            ),
            "full_corruption_calibration": {
                "expected_distribution": (
                    "100% corruption is generated by the same matched background-mapping mechanism "
                    "as the frozen specificity null"
                ),
                "median_check": (
                    "100%-corruption median rho_edge must fall inside the frozen final mapping-null "
                    "2.5th-97.5th percentile interval"
                ),
                "false_specificity_check": (
                    f"among 100 independent 100%-corruption replicates, no more than "
                    f"{int(FULL_CORRUPTION_MAX_FALSE_SPECIFIC_FRACTION * 100)}% may have "
                    "nominal mapping-null p<0.05"
                ),
            },
            "loading_corruption": (
                "not a primary degradation endpoint; weight/sign corruption is Supplement-only "
                "because edge preservation is weight-invariant and loading corruption partly "
                "reproduces the loading permutation null by construction"
            ),
        },

        "scanb_sample_size_operating_characteristics": {
            "role": (
                "quantify dependence of preservation precision/effect on target sample size"
            ),
            "sample_sizes": SCANB_SUBSAMPLE_SIZES,
            "repeats_per_nonfull_size": SCANB_SUBSAMPLE_REPEATS,
            "sampling": "without replacement from the frozen 3,273 primary SCAN-B profiles",
            "statistics": [
                "edge rho on the fixed specificity edge subset",
                "loading rho",
            ],
            "full_size": (
                "n=3273 is evaluated once as the primary SCAN-B result rather than repeatedly"
            ),
            "no_threshold_tuning": True,
        },

        "scanb_technical_repeatability": {
            "same_platform_primary_ceiling": {
                "pairs": SAME_PLATFORM_TECH_PAIRS,
                "definition": (
                    "all 100 frozen technical pairs whose primary and replicate profiles "
                    "use the same sequencing platform"
                ),
                "statistic": (
                    "edge concordance between primary-profile and replicate-profile "
                    "correlation structures, computed on the same frozen module genes"
                ),
            },
            "cross_platform_secondary": {
                "pairs": CROSS_PLATFORM_TECH_PAIRS,
                "definition": (
                    "all 36 frozen technical pairs whose primary and replicate profiles "
                    "use different sequencing platforms"
                ),
                "role": "secondary platform/resequencing sensitivity",
            },
            "matched_n_transport_comparison": (
                "compare same-platform technical repeatability at n=100 with the frozen "
                "SCAN-B n=100 source-to-target subsampling distribution"
            ),
            "interpretation": (
                "technical repeatability is a measurement ceiling at matched n, not a biological transport null"
            ),
        },

        "computation": {
            "gpu_permitted_and_preferred_for_heavy_matrix_operations": True,
            "numerical_contract": (
                "scientific statistics must be reproducible to a predeclared validation tolerance "
                "against CPU reference calculations on deterministic subsets before GPU results are accepted"
            ),
            "gpu_does_not_change_statistical_rules": True,
            "base_seed": BASE_SEED,
        },

        "historical_dog2_gse239948": {
            "status": "separate later historical-known-target reanalysis",
            "not_used_to_tune_this_contract": True,
        },
    }

    json_out = OUT_DIR / "preservation_operating_characteristics_contract_v1.json"
    md_out = OUT_DIR / "preservation_operating_characteristics_contract_v1.md"

    json_out.write_text(json.dumps(contract, indent=2), encoding="utf-8")
    md_out.write_text(
        f"""# Paper 4 / TCBB — Preservation & Operating-Characteristics Contract v1

Status: **FROZEN BEFORE FIRST TARGET PRESERVATION STATISTIC**

## Primary preservation classes

- **Strong:** positive edge q<0.05 + positive loading q<0.05 +
  edge mapping-specificity q<0.05.
- **Partial:** mapping-specificity q<0.05 + exactly one positive direct
  axis q<0.05, with no significantly negative axis.
- **Directionally discordant:** any direct axis is significantly negative;
  this takes precedence.
- **No clear:** remaining assessable cases.
- **Not assessable:** frozen coverage/gene-count guard fails.

The old fixed split-half >=0.60 cutoff is not part of classification.

## Direct inference

- {DIRECT_PERMUTATIONS:,} gene-label permutations per direct edge/loading test.
- BH correction across all assessable direct edge+loading tests within each target.

## Uncertainty

- {BOOTSTRAP_REPLICATES:,} target-sample bootstraps.
- Frozen source program remains fixed.
- {int(BOOTSTRAP_CI * 100)}% percentile CIs for edge and loading effects.

## Controlled mapping degradation

Fractions: {", ".join(str(x) for x in DEGRADATION_FRACTIONS)}.

At each nonzero fraction, {DEGRADATION_REPLICATES} deterministic corruption
replicates replace the selected correct target mappings with unrelated,
marginally matched genes while preserving source slots/edges/weights.

Formal trend testing is applied only when the uncorrupted edge result is
positive and specificity-supported.

## SCAN-B operating regimes

Target sample sizes:
{", ".join(str(x) for x in SCANB_SUBSAMPLE_SIZES)}.

For n < 3273, {SCANB_SUBSAMPLE_REPEATS} subsamples are used.

Technical-repeatability ceiling:
- 100 same-platform replicate pairs (primary ceiling)
- 36 cross-platform replicate pairs (secondary sensitivity)

## Computation

GPU acceleration is permitted and preferred for heavy matrix operations, but
must pass deterministic CPU-reference numerical validation before scientific
results are accepted.
""",
        encoding="utf-8",
    )

    print("=" * 128)
    print("04a PRESERVATION / OPERATING-CHARACTERISTICS CONTRACT: PASS")
    print("=" * 128)
    print(f"Direct permutations:          {DIRECT_PERMUTATIONS:,}")
    print(f"Target bootstraps:            {BOOTSTRAP_REPLICATES:,}")
    print(f"Split-half repeats:           {SPLIT_HALF_REPEATS:,}")
    print(f"Degradation reps/fraction:    {DEGRADATION_REPLICATES}")
    print(
        "SCAN-B subsample sizes:       "
        + ", ".join(str(x) for x in SCANB_SUBSAMPLE_SIZES)
    )
    print("Technical same-platform pairs:100")
    print("Technical cross-platform pairs:36")
    print("GPU heavy computation:        permitted + preferred with CPU validation")
    print()
    print("No target preservation statistic was calculated.")
    print("Outputs:")
    print(f"  {json_out}")
    print(f"  {md_out}")
    print("=" * 128)


if __name__ == "__main__":
    main()
