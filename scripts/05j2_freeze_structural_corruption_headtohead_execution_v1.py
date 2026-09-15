from __future__ import annotations

import json
from pathlib import Path


SCRIPT_VERSION = "05j2-freeze-structural-corruption-headtohead-execution-v1-no-cli"

DATA_ROOT = Path(r"D:\paper4_tcbb_data")

INTERPRETATION_CONTRACT = (
    DATA_ROOT
    / "paper4_tcbb_comparator_interpretation_headtohead_contract_v1"
    / "comparator_interpretation_headtohead_contract_v1.json"
)
DEGRADATION_MASTER = (
    DATA_ROOT
    / "paper4_tcbb_controlled_mapping_degradation_v1"
    / "controlled_mapping_degradation_v1.json"
)
COMPARATOR_MASTER = (
    DATA_ROOT
    / "paper4_tcbb_postprimary_comparator_benchmark_v7"
    / "postprimary_comparator_benchmark_v7.json"
)

OUT_DIR = (
    DATA_ROOT
    / "paper4_tcbb_structural_corruption_headtohead_execution_contract_v1"
)
OUT_JSON = OUT_DIR / "structural_corruption_headtohead_execution_contract_v1.json"

FRACTIONS = [0.0, 0.10, 0.25, 0.50, 0.75, 1.0]
NONZERO_REPLICATES = 100
BOOTSTRAP_REPLICATES = 2000
BOOTSTRAP_BASE_SEED = 20260920

MTA_REPLAY_TOL = 5e-10
NETREP_BASELINE_TOL = 5e-8
WGCNA_BASELINE_TOL = 5e-8


def require(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Required file does not exist: {path}")


def main() -> None:
    print("=" * 158)
    print("Paper 4 / TCBB - freeze structural known-answer head-to-head execution")
    print("=" * 158)
    print(f"Script version: {SCRIPT_VERSION}")
    print()
    print("Scientific provenance:")
    print("  Baseline MTA results observed:                         YES")
    print("  Baseline NetRep/WGCNA comparator results observed:    YES")
    print("  Corruption head-to-head results observed:              NO")
    print("  Exact 05d generator source inspected/recovered:        YES")
    print("  Head-to-head outcome thresholds tuned here:            NO")
    print("=" * 158)

    for p in [
        INTERPRETATION_CONTRACT,
        DEGRADATION_MASTER,
        COMPARATOR_MASTER,
    ]:
        require(p)

    interp = json.loads(INTERPRETATION_CONTRACT.read_text(encoding="utf-8"))
    deg = json.loads(DEGRADATION_MASTER.read_text(encoding="utf-8"))
    comp = json.loads(COMPARATOR_MASTER.read_text(encoding="utf-8"))

    if interp.get("status") != "FROZEN_BEFORE_FIRST_NETREP_OR_WGCNA_COMPARATOR_RESULT":
        raise RuntimeError("05i3 interpretation/head-to-head contract has unexpected status.")
    if deg.get("status") != "CONTROLLED_DEGRADATION_ALL_FROZEN_GATES_PASS":
        raise RuntimeError("05d degradation master has unexpected status.")
    if comp.get("status") != "POSTPRIMARY_COMPARATOR_BENCHMARK_COMPLETE":
        raise RuntimeError("05i v7 comparator master has unexpected status.")

    if [float(x) for x in deg["fractions"]] != FRACTIONS:
        raise RuntimeError("05d fractions differ from frozen head-to-head fractions.")
    if int(deg["replicates_per_nonzero_fraction"]) != NONZERO_REPLICATES:
        raise RuntimeError("05d replicate count differs from frozen head-to-head count.")

    contract = {
        "contract_id": "paper4-tcbb-structural-corruption-headtohead-execution-v1",
        "script_version": SCRIPT_VERSION,
        "status": "FROZEN_BEFORE_FIRST_CORRUPTION_HEADTOHEAD_COMPARATOR_STATISTIC",

        "provenance": {
            "primary_mta_results_already_observed": True,
            "baseline_comparator_results_already_observed": True,
            "corruption_headtohead_results_observed_before_this_freeze": False,
            "inherits_05i3_falsifiable_interpretation": True,
            "exact_05d_generator_recovered_from_source": True,
        },

        "scope_amendment_relative_to_05i3": {
            "status": "EXECUTION_CLARIFICATION_BEFORE_HEADTOHEAD_RESULTS",
            "reason": (
                "05d deliberately degraded only the source-to-target gene correspondence "
                "for the structural edge endpoint; loading degradation was explicitly not "
                "a 05d endpoint. The head-to-head operating comparison is therefore confined "
                "to structural preservation statistics. This prevents mixing the known-answer "
                "structural corruption experiment with unrelated density/coherence/loading axes."
            ),
            "supersedes_only": (
                "the overly broad 05i3 instruction to retain every observed WGCNA/NetRep "
                "statistic in the corruption ladder"
            ),
            "does_not_change": [
                "05d corruption mappings",
                "MTA rho_edge endpoint",
                "NetRep primary cor.cor endpoint",
                "05i3 falsifiable interpretation",
                "baseline NetRep/WGCNA comparator results",
            ],
        },

        "exact_mapping_replay": {
            "source": (
                "dynamically import the exact local "
                "05d_run_controlled_mapping_degradation_gpu_v1.py generator"
            ),
            "rows": (
                "replay every saved 05d (target, program_id, fraction, replicate_id, "
                "attempt_id) row; fraction 0 uses the correct mapping"
            ),
            "seed_schedule": (
                "unchanged 05d SeedSequence([20260916,target_index,module_number,"
                "fraction_code,attempt_id])"
            ),
            "generator_functions": [
                "prepare_problem",
                "generate_mapping",
                "quality",
                "rng_for_attempt",
                "round_half_up",
            ],
            "mta_replay": (
                "recompute frozen-subset MTA rho_edge from each replayed mapping and require "
                "agreement with controlled_degradation_replicates_v1.tsv"
            ),
            "mta_absolute_tolerance": MTA_REPLAY_TOL,
            "quality_fields_replayed": [
                "n_corrupted_slots",
                "realized_corruption_fraction",
                "fraction_distance_0",
                "fraction_distance_le1",
                "maximum_distance",
                "mean_distance",
            ],
            "failure_rule": (
                "any replay mismatch stops the entire head-to-head before accepting "
                "comparator corruption statistics"
            ),
        },

        "structural_methods": {
            "MTA": {
                "statistic": "rho_edge",
                "definition": (
                    "existing 05d Spearman concordance on exact corrected 05b v2 "
                    "fixed edge subset"
                ),
                "larger_means_stronger_preservation": True,
            },
            "NetRep": {
                "statistic": "cor.cor",
                "definition": (
                    "Pearson concordance between all unique within-module source and "
                    "mapped-target Pearson correlation coefficients"
                ),
                "implementation": (
                    "package-equivalent direct calculation, accepted only if fraction-0 "
                    "value reproduces completed NetRep v7 observed cor.cor"
                ),
                "fraction0_absolute_tolerance": NETREP_BASELINE_TOL,
                "larger_means_stronger_preservation": True,
            },
            "WGCNA": {
                "statistics": ["cor.cor", "cor.kIM"],
                "cor_cor_definition": (
                    "Pearson concordance between all unique within-module source and "
                    "mapped-target Pearson correlation coefficients"
                ),
                "cor_kIM_definition": (
                    "Pearson concordance of intramodular connectivity vectors using "
                    "WGCNA modulePreservation signed adjacency ((1+r)/2)^12"
                ),
                "implementation": (
                    "package-equivalent direct structural calculation, accepted only if "
                    "fraction-0 values reproduce raw observed WGCNA v7 preservation values"
                ),
                "fraction0_absolute_tolerance": WGCNA_BASELINE_TOL,
                "larger_means_stronger_preservation": True,
            },
        },

        "baseline_validation": {
            "NetRep_reference": (
                "paper4_tcbb_postprimary_comparator_benchmark_v7/"
                "netrep_all_targets_long_v7.tsv"
            ),
            "WGCNA_reference": (
                "raw preservation$observed values extracted from the already completed "
                "v7 WGCNA RDS objects; extraction performs no scientific recalculation"
            ),
            "rule": (
                "no nonzero corruption comparator statistic is accepted for a "
                "target/program until its fraction-0 direct calculation validates"
            ),
        },

        "normalization": {
            "formula": "(S(f)-median[S(1.0)])/(S(0)-median[S(1.0)])",
            "performed_separately_for_each": "target/program/method",
            "no_clipping": True,
            "interpretation": (
                "1 is the uncorrupted observed level and 0 is the method's own median "
                "full-corruption level; lower response at a given corruption fraction "
                "means greater sensitivity to deliberately wrong correspondence"
            ),
        },

        "operating_summaries": {
            "fractions": FRACTIONS,
            "median_normalized_response_by_fraction": True,
            "auc": {
                "definition": "trapezoidal area under median normalized response versus fraction",
                "direction": "lower AUC = faster loss under known-wrong mapping",
            },
            "f50": {
                "definition": (
                    "first corruption fraction where the median normalized response crosses "
                    "0.5, using linear interpolation between adjacent frozen fractions"
                ),
                "direction": "lower f50 = faster loss under known-wrong mapping",
                "if_no_crossing": "NA; do not invent extrapolation beyond fraction 1",
            },
            "trend": (
                "Spearman correlation between frozen corruption fractions and median raw statistic"
            ),
        },

        "paired_bootstrap": {
            "replicates": BOOTSTRAP_REPLICATES,
            "base_seed": BOOTSTRAP_BASE_SEED,
            "unit": "saved valid corruption replicate within each nonzero fraction",
            "pairing": (
                "the same resampled replicate indices are used for MTA and every comparator"
            ),
            "fraction_zero": "fixed unique uncorrupted observation",
            "normalization_within_bootstrap": (
                "each bootstrap draw uses that draw's own fraction-1 median as the method-specific floor"
            ),
            "interval": "95% percentile",
            "primary_contrast": {
                "definition": "AUC_MTA - AUC_comparator",
                "negative": "MTA more responsive to known mapping corruption",
                "positive": "comparator more responsive",
            },
            "secondary_contrast": {
                "definition": "f50_MTA - f50_comparator",
                "negative": "MTA crosses 0.5 earlier",
                "positive": "comparator crosses 0.5 earlier",
            },
        },

        "interpretation": {
            "MTA_more_responsive": (
                "paired AUC contrast 95% interval entirely below zero supports greater "
                "MTA sensitivity to mapping corruption for this target/program only"
            ),
            "comparator_more_responsive": (
                "paired AUC contrast 95% interval entirely above zero challenges any "
                "claim of greater MTA corruption discrimination"
            ),
            "indistinguishable": (
                "interval overlapping zero means no incremental operating advantage demonstrated"
            ),
            "overall_reporting": (
                "report counts and distributions across the 22 assessable pairs; no new "
                "post-hoc majority threshold defines overall superiority"
            ),
        },

        "fractions": FRACTIONS,
        "replicates_per_nonzero_fraction": NONZERO_REPLICATES,
        "primary_mta_classification_changed": False,
        "target_outcomes_or_treatment_loaded": False,
    }

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    OUT_JSON.write_text(json.dumps(contract, indent=2), encoding="utf-8")

    print()
    print("=" * 158)
    print("05j2 STRUCTURAL CORRUPTION HEAD-TO-HEAD EXECUTION CONTRACT: PASS")
    print("=" * 158)
    print("Exact 05d mapping replay required:             YES")
    print(f"MTA replay tolerance:                          {MTA_REPLAY_TOL:g}")
    print("Primary comparator:                            NetRep cor.cor")
    print("WGCNA structural descriptors:                  cor.cor + cor.kIM")
    print("Non-structural comparator ladder statistics:   NOT USED")
    print(f"Paired bootstrap replicates:                    {BOOTSTRAP_REPLICATES}")
    print("Negative AUC_MTA-AUC_comparator favors MTA:     YES")
    print("New overall superiority threshold:              NO")
    print("Primary classifications changed:                NO")
    print()
    print(f"Output: {OUT_JSON}")
    print("=" * 158)


if __name__ == "__main__":
    main()
