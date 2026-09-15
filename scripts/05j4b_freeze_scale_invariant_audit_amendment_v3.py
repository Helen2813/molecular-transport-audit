from __future__ import annotations

import json
from pathlib import Path


SCRIPT_VERSION = "05j4b-freeze-scale-invariant-audit-amendment-v3-no-cli"

DATA_ROOT = Path(r"D:\paper4_tcbb_data")

V1_CONTRACT = (
    DATA_ROOT
    / "paper4_tcbb_scale_invariant_corruption_audit_contract_v1"
    / "scale_invariant_corruption_audit_contract_v1.json"
)
HEADTOHEAD_MASTER = (
    DATA_ROOT
    / "paper4_tcbb_structural_corruption_headtohead_v1"
    / "structural_corruption_headtohead_v1.json"
)
INTERPRETATION_CONTRACT = (
    DATA_ROOT
    / "paper4_tcbb_comparator_interpretation_headtohead_contract_v1"
    / "comparator_interpretation_headtohead_contract_v1.json"
)

V1_RESULT_DIR = DATA_ROOT / "paper4_tcbb_scale_invariant_corruption_audit_v1"
V2_RESULT_DIR = DATA_ROOT / "paper4_tcbb_scale_invariant_corruption_audit_v2"

OUT_DIR = DATA_ROOT / "paper4_tcbb_scale_invariant_corruption_audit_contract_v3"
OUT_JSON = OUT_DIR / "scale_invariant_corruption_audit_contract_v3.json"

ALPHA = 0.05
FULL_NULL_VALID = 1000
CAL_N = 800
HOLDOUT_N = 200
PER_PARTIAL_FRACTION = 100
MAX_ATTEMPTS_PER_NEW_FRACTION = 1000

FRACTIONS = [0.0, 0.05, 0.10, 0.15, 0.25, 0.50, 0.75, 1.0]
PERFORMANCE_FRACTIONS = [0.0, 0.05, 0.10, 0.15, 0.25, 0.50, 0.75]

PAIR_BOOTSTRAPS = 2000
CLUSTER_BOOTSTRAPS = 2000
BASE_SEED = 20260922


def require(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Required file does not exist: {path}")


def ensure_no_results(path: Path, label: str) -> None:
    if not path.exists():
        return
    files = [p for p in path.rglob("*") if p.is_file()]
    if files:
        raise RuntimeError(
            f"{label} contains result files, but v3 assumes no scale-invariant "
            "scientific result has been produced yet:\n"
            + "\n".join(str(p) for p in files)
        )


def main() -> None:
    print("=" * 164)
    print("Paper 4 / TCBB - freeze scale-invariant mapping-error audit v3")
    print("=" * 164)
    print(f"Script version: {SCRIPT_VERSION}")
    print()
    print("Provenance:")
    print("  05j normalized AUC/f50 already observed:              YES")
    print("  05j4 v1 scale-invariant scientific results observed:  NO")
    print("  05j4 v2 scale-invariant scientific results observed:  NO")
    print("  v3 frozen before any scale-invariant result:           YES")
    print("  Primary MTA statistic/classifier changed:              NO")
    print("=" * 164)

    for p in [V1_CONTRACT, HEADTOHEAD_MASTER, INTERPRETATION_CONTRACT]:
        require(p)

    v1 = json.loads(V1_CONTRACT.read_text(encoding="utf-8"))
    h = json.loads(HEADTOHEAD_MASTER.read_text(encoding="utf-8"))
    i = json.loads(INTERPRETATION_CONTRACT.read_text(encoding="utf-8"))

    if v1.get("status") != "FROZEN_POST_05J3_BEFORE_SCALE_INVARIANT_AUDIT_RESULTS":
        raise RuntimeError("05j4 v1 contract has unexpected status.")
    if h.get("status") != "STRUCTURAL_CORRUPTION_HEADTOHEAD_COMPLETE":
        raise RuntimeError("05j head-to-head master has unexpected status.")
    if i.get("status") != "FROZEN_BEFORE_FIRST_NETREP_OR_WGCNA_COMPARATOR_RESULT":
        raise RuntimeError("05i3 interpretation contract has unexpected status.")

    ensure_no_results(V1_RESULT_DIR, "05j4 v1 result directory")
    ensure_no_results(V2_RESULT_DIR, "05j4 v2 result directory")

    contract = {
        "contract_id": "paper4-tcbb-scale-invariant-corruption-audit-v3",
        "script_version": SCRIPT_VERSION,
        "status": "FROZEN_POST_05J3_BEFORE_SCALE_INVARIANT_AUDIT_V3_RESULTS",

        "provenance": {
            "normalized_auc_f50_05j_results_already_observed": True,
            "scale_invariant_followup_is_post_hoc_relative_to_05j3": True,
            "scale_invariant_v1_results_observed": False,
            "scale_invariant_v2_results_observed": False,
            "scale_invariant_v3_results_observed_before_freeze": False,
            "primary_mta_statistic_changed": False,
            "primary_mta_classifier_changed": False,
        },

        "why_v3": [
            "replace pooled-across-fractions redundancy correlation with within-fraction correlation",
            "extend 100%-corruption null from 100 to 1000 valid mappings",
            "split full-null mappings into 800 calibration and 200 held-out mappings",
            "add corruption fractions 0.05 and 0.15 to improve f50 resolution",
            "predeclare censoring/interpolation rules for f50",
            "add same-edge Pearson diagnostic to separate correlation functional form from edge-set scope",
            "avoid treating 22 target/program pairs as independent trials",
            "freeze Plan B before scale-invariant results",
        ],

        "known_null_extension": {
            "generator": "exact frozen 05d mapping-corruption generator",
            "valid_full_corruption_mappings_per_pair": FULL_NULL_VALID,
            "maximum_attempts": 10000,
            "first_100_replay_guard": (
                "first 100 valid full-corruption mappings must reproduce existing "
                "05d/05j attempt IDs and MTA/NetRep/WGCNA statistics"
            ),
            "calibration": {
                "valid_ids": "1-800",
                "n": CAL_N,
            },
            "heldout_validation": {
                "valid_ids": "801-1000",
                "n": HOLDOUT_N,
                "scope_statement": (
                    "held-out false-preservation is validated only on independent mappings "
                    "from the same frozen full-corruption generator, not on other mapping-error mechanisms"
                ),
            },
        },

        "fraction_grid": {
            "fractions": FRACTIONS,
            "existing_05j_fractions_reused": [0.0, 0.10, 0.25, 0.50, 0.75],
            "new_post_05j_fractions": [0.05, 0.15],
            "new_fraction_generator": "same exact 05d generator and quality rule",
            "valid_replicates_per_new_fraction": PER_PARTIAL_FRACTION,
            "max_attempts_per_new_fraction": MAX_ATTEMPTS_PER_NEW_FRACTION,
            "reason": (
                "increase horizontal resolution around the expected f50 crossing; "
                "these fractions are explicitly post-05j robustness extensions"
            ),
        },

        "method_decomposition": {
            "MTA": {
                "correlation": "Spearman",
                "edges": "exact corrected 05b deterministic <=20k subset",
            },
            "Pearson_subset_matched": {
                "correlation": "Pearson",
                "edges": "the exact same <=20k subset as MTA",
                "purpose": "isolate rank-correlation versus Pearson functional-form effect",
                "published_comparator": False,
            },
            "NetRep_cor_cor": {
                "correlation": "Pearson",
                "edges": "all unique within-module edges",
                "purpose": "published structural comparator",
            },
            "sequential_factor_decomposition": [
                "MTA versus Pearson_subset_matched: correlation functional form at identical edges",
                "Pearson_subset_matched versus NetRep_cor_cor: edge-set scope at identical Pearson functional form",
            ],
            "Spearman_all_edges_not_computed": {
                "reason": (
                    "the two main design differences are already isolated sequentially; "
                    "an all-edge Spearman interaction diagnostic would require ranking up to "
                    "millions of edges for thousands of mappings and is not required for the "
                    "predeclared primary decomposition"
                ),
                "limitation": (
                    "a correlation-form-by-edge-scope interaction is therefore not separately estimated"
                ),
            },
            "WGCNA_cor_kIM": {
                "role": "secondary connectivity/topology comparator",
            },
        },

        "equal_false_preservation_calibration": {
            "alpha": ALPHA,
            "positive_call": "PRESERVED",
            "direction": "larger statistic = stronger preservation",
            "threshold_data": "800 calibration full-corruption mappings",
            "threshold_rule": (
                "conservative empirical 95th-percentile order-statistic rule; "
                "PRESERVED only for statistic > threshold"
            ),
            "threshold_bootstrap_interval": "95% percentile",
            "heldout_false_preservation": (
                "measured on 200 independent full-corruption mappings from the same generator"
            ),
            "strictly_monotone_rescaling_invariant": True,
        },

        "within_fraction_redundancy": {
            "computed_separately_at": [0.05, 0.10, 0.15, 0.25, 0.50, 0.75, 1.0],
            "statistic": "Spearman across mapping replicates at fixed corruption fraction",
            "n_partial": 100,
            "n_fraction_1": 200,
            "pooled_across_fraction_correlation_not_used": True,
        },

        "f50": {
            "headline_resolution_metric": True,
            "definition": (
                "first corruption fraction where preservation-call rate crosses from >0.5 "
                "to <=0.5"
            ),
            "interpolation": "piecewise linear between adjacent frozen grid points",
            "if_fraction_zero_rate_le_0p5": {
                "value": 0.0,
                "status": "FAILED_AT_UNCORRUPTED_MAPPING",
            },
            "if_no_crossing_by_fraction_1": {
                "value": None,
                "status": "RIGHT_CENSORED_GT_1",
                "numeric_extrapolation": False,
            },
            "pairwise_delta_rule": (
                "numeric delta-f50 is reported only when both methods have finite f50; "
                "if either method is right-censored, report censoring statuses and direction "
                "when logically determined, but do not invent a numeric difference"
            ),
            "primary_text_use": (
                "f50 is the headline interpretable resolution boundary, with explicit censoring"
            ),
        },

        "retention_auc": {
            "primary_inferential_scale_invariant_contrast": True,
            "fractions": PERFORMANCE_FRACTIONS,
            "fraction_1_excluded": True,
            "definition": "trapezoidal area under preservation-call-rate curve",
            "reason": (
                "always estimable even when f50 is right-censored; used for paired inferential comparison"
            ),
        },

        "uncertainty": {
            "within_pair_bootstrap": {
                "replicates": PAIR_BOOTSTRAPS,
                "base_seed": BASE_SEED,
                "reestimate_threshold_each_draw": True,
                "paired_mapping_indices_across_methods": True,
                "resample_calibration_null": True,
                "resample_heldout_null": True,
                "resample_partial_fraction_replicates": True,
            },
            "aggregate_module_cluster_bootstrap": {
                "replicates": CLUSTER_BOOTSTRAPS,
                "cluster": "source program_id/module",
                "sampling": (
                    "sample 12 source modules with replacement and carry all available target-cohort "
                    "pairs for each selected module; within each selected pair draw one stored paired "
                    "mapping-bootstrap AUC contrast"
                ),
                "primary_aggregate": "median MTA-minus-comparator retention AUC",
                "reason": (
                    "the 22 target/program pairs are not treated as independent trials because the "
                    "same source modules recur across target cohorts"
                ),
            },
            "no_binomial_test_over_22_pairs": True,
        },

        "interpretation": {
            "mapping_failure_domain_only": (
                "mapping corruption is explicit in MTA but not the primary design target of "
                "NetRep/WGCNA; any surviving advantage is domain-specific mapping-error sensitivity, "
                "not general preservation-method superiority"
            ),
            "if_same_edge_diagnostic_erases_difference": (
                "do not claim intrinsic MTA advantage; attribute the apparent difference primarily "
                "to edge-set scope and/or design choice"
            ),
            "if_equal_fpr_auc_advantage_survives": (
                "supports earlier loss of preservation calls under known mapping error at matched "
                "held-out false-preservation calibration"
            ),
            "if_no_advantage": (
                "normalized-AUC differences from 05j3 are treated as scale-dependent; "
                "no operating superiority claim is made"
            ),
        },

        "plan_B_frozen_before_v3_results": {
            "trigger": (
                "no coherent scale-invariant MTA advantage, comparator advantage, or apparent "
                "advantage explained by same-edge diagnostic"
            ),
            "paper_reframe": (
                "evaluation/benchmark paper: empirical redundancy and saturation of preservation "
                "statistics on real cohorts, known-answer resolution limits, sample-size dependence, "
                "technical ceilings, and executable provenance"
            ),
            "claims_removed": [
                "MTA is a superior preservation statistic",
                "MTA has general operating superiority over NetRep/WGCNA",
            ],
            "sign_or_weight_corruption_cannot_rescue_superiority_claim": True,
        },

        "alpha": ALPHA,
        "fractions": FRACTIONS,
        "performance_fractions": PERFORMANCE_FRACTIONS,
        "pair_bootstraps": PAIR_BOOTSTRAPS,
        "cluster_bootstraps": CLUSTER_BOOTSTRAPS,
        "primary_mta_classification_changed": False,
    }

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    OUT_JSON.write_text(json.dumps(contract, indent=2), encoding="utf-8")

    print()
    print("=" * 164)
    print("05j4b SCALE-INVARIANT AUDIT v3 CONTRACT: PASS")
    print("=" * 164)
    print("Full-corruption valid mappings/pair:     1,000")
    print("Calibration / held-out null:             800 / 200")
    print("Corruption grid:                         0,.05,.10,.15,.25,.50,.75,1")
    print("Headline resolution metric:              f50 with censoring")
    print("Primary inferential contrast:            retention-call AUC")
    print("Same-edge Pearson diagnostic:            YES")
    print("All-edge Spearman diagnostic:            NO (interaction not estimated)")
    print("Within-fraction correlations:            YES")
    print("Module-cluster aggregate bootstrap:      YES")
    print("Plan B frozen before results:            YES")
    print("Primary MTA method changed:              NO")
    print()
    print(f"Output: {OUT_JSON}")
    print("=" * 164)


if __name__ == "__main__":
    main()
