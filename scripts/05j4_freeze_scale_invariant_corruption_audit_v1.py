from __future__ import annotations

import json
from pathlib import Path


SCRIPT_VERSION = "05j4-freeze-scale-invariant-corruption-audit-v1-no-cli"

DATA_ROOT = Path(r"D:\paper4_tcbb_data")

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

OUT_DIR = DATA_ROOT / "paper4_tcbb_scale_invariant_corruption_audit_contract_v1"
OUT_JSON = OUT_DIR / "scale_invariant_corruption_audit_contract_v1.json"

ALPHA = 0.05
FRACTIONS = [0.0, 0.10, 0.25, 0.50, 0.75, 1.0]
PERFORMANCE_FRACTIONS = [0.0, 0.10, 0.25, 0.50, 0.75]
BOOTSTRAPS = 2000
BASE_SEED = 20260921


def require(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Required file does not exist: {path}")


def main() -> None:
    print("=" * 158)
    print("Paper 4 / TCBB - freeze scale-invariant corruption-detection audit")
    print("=" * 158)
    print(f"Script version: {SCRIPT_VERSION}")
    print()
    print("Provenance:")
    print("  05j normalized-AUC results already observed:          YES")
    print("  This audit was motivated by scale/reparameterization critique: YES")
    print("  Results of this scale-invariant audit observed:        NO")
    print("  Primary MTA statistics/classification changed:         NO")
    print("=" * 158)

    for p in [
        HEADTOHEAD_MASTER,
        HEADTOHEAD_REPLICATES,
        INTERPRETATION_CONTRACT,
    ]:
        require(p)

    h = json.loads(HEADTOHEAD_MASTER.read_text(encoding="utf-8"))
    i = json.loads(INTERPRETATION_CONTRACT.read_text(encoding="utf-8"))

    if h.get("status") != "STRUCTURAL_CORRUPTION_HEADTOHEAD_COMPLETE":
        raise RuntimeError("05j head-to-head master has unexpected status.")
    if i.get("status") != "FROZEN_BEFORE_FIRST_NETREP_OR_WGCNA_COMPARATOR_RESULT":
        raise RuntimeError("05i3 interpretation contract has unexpected status.")

    contract = {
        "contract_id": "paper4-tcbb-scale-invariant-corruption-audit-v1",
        "script_version": SCRIPT_VERSION,
        "status": "FROZEN_POST_05J3_BEFORE_SCALE_INVARIANT_AUDIT_RESULTS",

        "provenance": {
            "normalized_auc_results_already_observed": True,
            "motivation": (
                "The 05j3 normalized response AUC can differ under monotone "
                "reparameterization. This follow-up asks whether the apparent operating "
                "difference remains after calibrating each method to the same known-null "
                "false-preservation rate."
            ),
            "post_hoc_relative_to_05j3": True,
            "thresholds_or_interpretation_tuned_after_this_audit": False,
        },

        "important_clarification": {
            "auc_f50_predeclared_before_05j_results": True,
            "where": (
                "05i3 comparator interpretation/head-to-head contract and "
                "05j2 structural head-to-head execution contract"
            ),
            "therefore": (
                "05j AUC/f50 were not invented after seeing corruption curves, but this "
                "scale-invariant detector audit is a new robustness analysis motivated "
                "after observing 05j3 effect sizes."
            ),
        },

        "task_definition": {
            "positive_class": "preserved/transport-supported mapping",
            "known_null": (
                "100% corruption replicates from the exact frozen 05d/05j mapping generator"
            ),
            "type_I_error": (
                "false preservation call under fully corrupted mapping"
            ),
            "alpha": ALPHA,
            "reason_for_full_corruption_null": (
                "0% corruption is the true uncorrupted mapping and is therefore not a null "
                "for preservation. Full corruption is the known-answer non-transport mapping "
                "and was already calibrated against the frozen matched mapping null in 05d."
            ),
        },

        "method_specific_threshold": {
            "direction": "larger statistic = stronger preservation for every method used",
            "rule": (
                "Within each target/program/method, sort the 100 full-corruption scores. "
                "Let k=floor(alpha*n). Use the (n-k)th order statistic as threshold and "
                "call PRESERVED only when score > threshold. With n=100 and alpha=.05, "
                "this guarantees at most 5/100 false-preservation calls, modulo ties "
                "(ties make the rule more conservative)."
            ),
            "monotone_invariance": (
                "Any strictly increasing transformation of a method's statistic preserves "
                "both the null quantile threshold ordering and all resulting calls."
            ),
        },

        "methods": {
            "MTA": "mta_rho_edge_saved",
            "NetRep_cor.cor": "netrep_cor_cor",
            "WGCNA_cor.kIM": "wgcna_cor_kIM",
            "WGCNA_cor.cor": (
                "retained as an identity/redundancy audit only; not counted as an "
                "independent comparator because it is numerically the same cor.cor axis"
            ),
        },

        "outputs": {
            "raw_monotone_redundancy": (
                "Spearman correlation between MTA and each comparator across all 501 "
                "corruption rows within each target/program"
            ),
            "preservation_call_rate": (
                "fraction of replicates exceeding the method-specific full-corruption "
                "threshold at each frozen corruption fraction"
            ),
            "failure_detection_rate": "1 - preservation_call_rate",
            "primary_scale_invariant_summary": (
                "area under preservation-call-rate curve from corruption fraction 0 through "
                "0.75 only; fraction 1.0 is the calibration set and is not used as a "
                "performance endpoint"
            ),
            "direction": (
                "lower retention AUC means preservation calls disappear earlier as mapping "
                "corruption increases"
            ),
        },

        "paired_bootstrap": {
            "replicates": BOOTSTRAPS,
            "base_seed": BASE_SEED,
            "pairing": (
                "resample the same replicate indices across methods within each fraction; "
                "also bootstrap the 100% corruption calibration rows and re-estimate each "
                "method's threshold inside every bootstrap"
            ),
            "fraction_zero": "fixed unique uncorrupted mapping",
            "primary_contrast": "retention_AUC_MTA - retention_AUC_comparator",
            "negative": "MTA loses preservation calls earlier at equal null false-preservation rate",
            "positive": "comparator loses preservation calls earlier",
            "interval": "95% percentile",
        },

        "interpretation": {
            "mta_advantage_survives": (
                "95% paired-bootstrap interval for MTA-minus-comparator retention AUC is "
                "entirely below zero"
            ),
            "comparator_advantage": (
                "interval entirely above zero"
            ),
            "no_clear_difference": (
                "interval includes zero"
            ),
            "critical_claim_rule": (
                "Only this scale-invariant calibrated-call comparison may support an "
                "operating sensitivity claim that is robust to monotone rescaling. "
                "The earlier normalized-AUC result remains descriptive."
            ),
        },

        "scope_caveat": {
            "mapping_corruption": (
                "The experiment targets correspondence/mapping failure, a failure mode "
                "explicitly modeled by MTA but not the primary design target of NetRep/WGCNA. "
                "Therefore any surviving advantage is reported as domain-specific mapping-error "
                "sensitivity, not general method superiority."
            )
        },

        "alpha": ALPHA,
        "fractions": FRACTIONS,
        "performance_fractions": PERFORMANCE_FRACTIONS,
        "bootstrap_replicates": BOOTSTRAPS,
        "primary_mta_classification_changed": False,
    }

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    OUT_JSON.write_text(json.dumps(contract, indent=2), encoding="utf-8")

    print()
    print("=" * 158)
    print("05j4 SCALE-INVARIANT CORRUPTION AUDIT CONTRACT: PASS")
    print("=" * 158)
    print("Known-null calibration:                  100% corruption")
    print("False-preservation target:               <=5%")
    print("Threshold rule:                          exact conservative order statistic")
    print("Performance fractions:                   0, .10, .25, .50, .75")
    print("Primary summary:                         retention-call AUC")
    print("Strict monotone-rescaling invariant:     YES")
    print(f"Paired bootstrap replicates:              {BOOTSTRAPS}")
    print("General superiority claim permitted:      NO")
    print("Domain-specific mapping-error claim only: YES")
    print()
    print(f"Output: {OUT_JSON}")
    print("=" * 158)


if __name__ == "__main__":
    main()
