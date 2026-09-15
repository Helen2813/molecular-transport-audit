from __future__ import annotations

import json
from pathlib import Path


SCRIPT_VERSION = "05j4a-freeze-scale-invariant-audit-amendment-v2-no-cli"

DATA_ROOT = Path(r"D:\paper4_tcbb_data")

V1_CONTRACT = (
    DATA_ROOT
    / "paper4_tcbb_scale_invariant_corruption_audit_contract_v1"
    / "scale_invariant_corruption_audit_contract_v1.json"
)
V1_RESULT_DIR = DATA_ROOT / "paper4_tcbb_scale_invariant_corruption_audit_v1"

HEADTOHEAD_MASTER = (
    DATA_ROOT
    / "paper4_tcbb_structural_corruption_headtohead_v1"
    / "structural_corruption_headtohead_v1.json"
)
HEADTOHEAD_REPLICATES = (
    DATA_ROOT
    / "paper4_tcbb_structural_corruption_headtohead_v1"
    / "structural_headtohead_replicates_v1.tsv"
)
INTERPRETATION_CONTRACT = (
    DATA_ROOT
    / "paper4_tcbb_comparator_interpretation_headtohead_contract_v1"
    / "comparator_interpretation_headtohead_contract_v1.json"
)

OUT_DIR = DATA_ROOT / "paper4_tcbb_scale_invariant_corruption_audit_contract_v2"
OUT_JSON = OUT_DIR / "scale_invariant_corruption_audit_contract_v2.json"

ALPHA = 0.05
FULL_NULL_VALID_REPLICATES = 1000
CALIBRATION_NULL_REPLICATES = 800
HOLDOUT_NULL_REPLICATES = 200
MAX_ATTEMPTS_FULL_NULL = 10000

FRACTIONS = [0.0, 0.10, 0.25, 0.50, 0.75, 1.0]
PERFORMANCE_FRACTIONS = [0.0, 0.10, 0.25, 0.50, 0.75]

BOOTSTRAPS = 2000
BASE_SEED = 20260922


def require(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Required file does not exist: {path}")


def main() -> None:
    print("=" * 162)
    print("Paper 4 / TCBB - freeze amended scale-invariant mapping-error audit v2")
    print("=" * 162)
    print(f"Script version: {SCRIPT_VERSION}")
    print()
    print("Provenance:")
    print("  05j normalized-AUC/f50 results already observed:       YES")
    print("  05j4 v1 scale-invariant results observed:               NO")
    print("  05j4 v1 runner failed before first result (key-name bug): REQUIRED")
    print("  v2 changes motivated by methodological critique:        YES")
    print("  Primary MTA statistic/classifier changed:               NO")
    print("=" * 162)

    for p in [
        V1_CONTRACT,
        HEADTOHEAD_MASTER,
        HEADTOHEAD_REPLICATES,
        INTERPRETATION_CONTRACT,
    ]:
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

    if V1_RESULT_DIR.exists():
        scientific = [
            p for p in V1_RESULT_DIR.rglob("*")
            if p.is_file()
        ]
        if scientific:
            raise RuntimeError(
                "05j4 v1 result files already exist. This v2 amendment assumes the "
                "v1 runner failed before writing any scientific result:\n"
                + "\n".join(str(p) for p in scientific)
            )

    contract = {
        "contract_id": "paper4-tcbb-scale-invariant-corruption-audit-v2",
        "script_version": SCRIPT_VERSION,
        "status": "FROZEN_POST_05J3_BEFORE_SCALE_INVARIANT_AUDIT_V2_RESULTS",

        "supersedes": {
            "contract": str(V1_CONTRACT),
            "reason": (
                "05j4 v1 runner stopped on the first target/program before producing "
                "any scale-invariant result because of an internal method-key naming bug. "
                "Before rerunning, v2 incorporates methodological robustness points raised "
                "after 05j3 but before any 05j4 result."
            ),
            "v1_scientific_results_calculated": False,
        },

        "provenance": {
            "normalized_auc_f50_05j_results_already_observed": True,
            "scale_invariant_audit_is_post_hoc_relative_to_05j3": True,
            "scale_invariant_v2_results_observed_before_this_freeze": False,
            "primary_mta_statistics_changed": False,
            "primary_mta_classifier_changed": False,
        },

        "known_null_extension": {
            "reason": (
                "A 5% tail threshold estimated from only 100 full-corruption replicates "
                "is coarse. v2 extends the same exact frozen 05d full-corruption generator "
                "to 1,000 valid full-corruption mappings per target/program."
            ),
            "generator": (
                "exact 05d SeedSequence/matched-background generator; valid mappings are "
                "accepted in attempt order"
            ),
            "valid_full_corruption_mappings": FULL_NULL_VALID_REPLICATES,
            "maximum_attempts": MAX_ATTEMPTS_FULL_NULL,
            "first_100_replay_guard": (
                "the first 100 valid mappings MUST reproduce the existing 05d/05j "
                "fraction=1 replicate_id, attempt_id, MTA rho, NetRep cor.cor, and "
                "WGCNA cor.kIM before mappings 101-1000 are accepted"
            ),
            "split": {
                "calibration": (
                    "valid full-corruption mappings 1-800; used only to estimate thresholds"
                ),
                "holdout": (
                    "valid full-corruption mappings 801-1000; used to estimate held-out "
                    "false-preservation rate and the fraction=1 endpoint for f50"
                ),
                "calibration_n": CALIBRATION_NULL_REPLICATES,
                "holdout_n": HOLDOUT_NULL_REPLICATES,
            },
        },

        "methods": {
            "MTA": {
                "statistic": "Spearman edge concordance on exact corrected 05b <=20k edge subset",
                "role": "primary MTA structural endpoint",
            },
            "Pearson_subset_matched": {
                "statistic": (
                    "Pearson concordance on the exact same corrected 05b <=20k edge subset"
                ),
                "role": (
                    "matched-input diagnostic isolating rank-correlation functional form "
                    "from the fixed-edge-subset design"
                ),
                "not_an_independent_published_comparator": True,
            },
            "NetRep_cor_cor": {
                "statistic": "NetRep-style Pearson cor.cor on all unique within-module edges",
                "role": "primary published structural comparator",
            },
            "WGCNA_cor_kIM": {
                "statistic": "WGCNA signed-adjacency intramodular-connectivity concordance",
                "role": "secondary topology comparator",
            },
            "WGCNA_cor_cor": (
                "not counted separately because it is numerically the same cor.cor axis "
                "as NetRep in this benchmark"
            ),
        },

        "edge_set_asymmetry": {
            "acknowledged": True,
            "MTA_primary_edges": "fixed deterministic <=20,000 subset",
            "NetRep_cor_cor_edges": "all unique module edges",
            "diagnostic": (
                "Pearson_subset_matched uses the same edge subset as MTA. "
                "This does not redefine MTA; it diagnoses whether any observed difference "
                "is attributable primarily to rank-vs-Pearson statistic or to edge-set scope."
            ),
            "interpretation_rule": (
                "No claim that MTA intrinsically outperforms cor.cor is permitted if the "
                "scale-invariant difference disappears under the same-edge diagnostic."
            ),
        },

        "threshold_calibration": {
            "known_null": "800 calibration full-corruption mappings",
            "alpha": ALPHA,
            "direction": "larger statistic = stronger preservation",
            "threshold_rule": (
                "sort the 800 calibration scores; k=floor(alpha*n)=40; threshold is "
                "the (n-k)th order statistic and PRESERVED requires score > threshold"
            ),
            "heldout_validation": (
                "report false-preservation rate on the 200 independent holdout "
                "full-corruption mappings"
            ),
            "strict_monotone_invariance": True,
            "threshold_uncertainty": (
                "bootstrap calibration mappings and report 95% percentile interval "
                "for each method-specific threshold"
            ),
        },

        "within_fraction_redundancy": {
            "statistic": "Spearman correlation between methods across mapping replicates",
            "computed_separately_at": [0.10, 0.25, 0.50, 0.75, 1.0],
            "fraction_1_data": "200 held-out full-corruption mappings",
            "reason": (
                "pooled correlation across all corruption fractions is dominated by the "
                "common corruption gradient and is not a valid redundancy diagnostic"
            ),
        },

        "primary_operating_endpoint": {
            "name": "f50",
            "definition": (
                "corruption fraction at which preservation-call rate first crosses 0.5, "
                "using linear interpolation across [0,.10,.25,.50,.75,1.0]"
            ),
            "fraction_1_call_rate": "estimated on 200 held-out full-corruption mappings",
            "if_uncorrupted_call_rate_le_0p5": "f50=0",
            "if_no_crossing_by_fraction_1": "f50=NA; no extrapolation",
            "direction": "lower f50 = earlier detection of mapping failure",
            "reason": "interpretable resolution boundary in corruption-fraction units",
        },

        "secondary_operating_endpoint": {
            "name": "retention_call_auc",
            "definition": (
                "trapezoidal area under preservation-call-rate curve on "
                "[0,.10,.25,.50,.75]; fraction 1 is excluded because it supplies the null endpoint"
            ),
            "direction": "lower AUC = preservation calls disappear earlier",
        },

        "paired_bootstrap": {
            "replicates": BOOTSTRAPS,
            "base_seed": BASE_SEED,
            "paired_indices_across_methods": True,
            "resample": {
                "calibration_full_corruption": 800,
                "holdout_full_corruption": 200,
                "each_partial_nonzero_fraction": 100,
                "fraction_zero": "fixed unique uncorrupted mapping",
            },
            "reestimate_threshold_each_draw": True,
            "primary_contrast": "f50_MTA - f50_comparator",
            "negative_primary_contrast": "MTA detects mapping failure earlier",
            "secondary_contrast": "retention_AUC_MTA - retention_AUC_comparator",
            "interval": "95% percentile",
        },

        "reporting_independence": {
            "22_pairs_are_not_treated_as_independent_replicates": True,
            "no_binomial_or_majority_test_over_22_pairs": True,
            "report": [
                "per-pair f50 and paired interval",
                "per-pair AUC and paired interval",
                "distribution summaries by target",
                "module-linked paired table across targets",
            ],
        },

        "interpretation": {
            "domain_specific_only": (
                "Mapping corruption is an explicit MTA failure mode but is not the primary "
                "design target of NetRep/WGCNA. Any advantage is therefore described as "
                "mapping-error sensitivity, not general preservation-method superiority."
            ),
            "if_scale_invariant_advantage_survives": (
                "A negative paired f50 contrast interval supports earlier mapping-failure "
                "detection for that target/program at matched held-out false-preservation calibration."
            ),
            "if_difference_disappears": (
                "The normalized-AUC difference is treated as scale-dependent; no operating "
                "advantage claim is made."
            ),
            "if_comparator_is_better": (
                "This directly challenges any MTA mapping-error sensitivity advantage and "
                "is reported without reinterpretation."
            ),
        },

        "plan_B_frozen_before_v2_results": {
            "trigger": (
                "scale-invariant f50/call-curve differences show no coherent MTA advantage "
                "or are explained by edge-set asymmetry"
            ),
            "paper_reframe": (
                "evaluation/benchmark study of molecular-program preservation statistics: "
                "strong real-cohort redundancy/saturation, known-answer resolution limits, "
                "sample-size dependence, technical ceilings, and reproducible audit provenance"
            ),
            "claims_removed": [
                "MTA is a superior preservation statistic",
                "MTA has general operating superiority over NetRep/WGCNA",
            ],
            "sign_corruption_not_used_as_rescue": (
                "a later sign/weight perturbation may only illustrate a distinct scope question "
                "(frozen weighted-definition transport) and cannot rescue a failed general "
                "superiority claim"
            ),
        },

        "fractions": FRACTIONS,
        "performance_fractions": PERFORMANCE_FRACTIONS,
        "bootstrap_replicates": BOOTSTRAPS,
        "primary_mta_classification_changed": False,
    }

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    OUT_JSON.write_text(json.dumps(contract, indent=2), encoding="utf-8")

    print()
    print("=" * 162)
    print("05j4a SCALE-INVARIANT AUDIT AMENDMENT v2: PASS")
    print("=" * 162)
    print("Full-corruption valid mappings/pair:       1,000")
    print("Calibration / holdout split:               800 / 200")
    print("Primary endpoint:                          f50")
    print("Secondary endpoint:                        retention-call AUC")
    print("Within-fraction correlation diagnostic:    YES")
    print("Matched-edge Pearson diagnostic:           YES")
    print("22 pairs treated as independent trials:    NO")
    print("Plan B frozen before results:              YES")
    print("Primary MTA method changed:                NO")
    print()
    print(f"Output: {OUT_JSON}")
    print("=" * 162)


if __name__ == "__main__":
    main()
