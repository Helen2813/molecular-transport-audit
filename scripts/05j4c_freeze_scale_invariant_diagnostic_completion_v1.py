#!/usr/bin/env python
from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

SCRIPT_VERSION = "05j4c-freeze-scale-invariant-diagnostic-completion-v1-no-cli"

DATA_ROOT = Path(r"D:\paper4_tcbb_data")

V3_CONTRACT = (
    DATA_ROOT
    / "paper4_tcbb_scale_invariant_corruption_audit_contract_v3"
    / "scale_invariant_corruption_audit_contract_v3.json"
)
V3_INPUT_MASTER = (
    DATA_ROOT
    / "paper4_tcbb_scale_invariant_corruption_inputs_v3"
    / "scale_invariant_inputs_v3.json"
)
V3_PERFORMANCE = (
    DATA_ROOT
    / "paper4_tcbb_scale_invariant_corruption_inputs_v3"
    / "scale_invariant_performance_grid_v3.tsv"
)
V3_FULLNULL = (
    DATA_ROOT
    / "paper4_tcbb_scale_invariant_corruption_inputs_v3"
    / "scale_invariant_fullnull_1000_v3.tsv"
)
V3_AUDIT_MASTER = (
    DATA_ROOT
    / "paper4_tcbb_scale_invariant_corruption_audit_v3"
    / "scale_invariant_corruption_audit_v3.json"
)
V3_THRESHOLDS = (
    DATA_ROOT
    / "paper4_tcbb_scale_invariant_corruption_audit_v3"
    / "scale_invariant_thresholds_v3.tsv"
)

OUT_DIR = DATA_ROOT / "paper4_tcbb_scale_invariant_diagnostic_completion_contract_v1"
OUT_JSON = OUT_DIR / "scale_invariant_diagnostic_completion_contract_v1.json"

SEP = "=" * 168

def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)

def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()

def main() -> None:
    print(SEP)
    print("Paper 4 / TCBB - freeze scale-invariant diagnostic completion")
    print(SEP)
    print(f"Script version: {SCRIPT_VERSION}")
    print()
    print("Scientific provenance:")
    print("  05j4 v3 scale-invariant scientific results observed:  YES")
    print("  05j4c held-out FPR results observed:                  NO")
    print("  05j4c all-edge Spearman results observed:             NO")
    print("  Existing 05j4 v3 thresholds/results changed here:     NO")
    print("  Primary MTA statistic/classifier changed:              NO")
    print()

    required = [
        V3_CONTRACT,
        V3_INPUT_MASTER,
        V3_PERFORMANCE,
        V3_FULLNULL,
        V3_AUDIT_MASTER,
        V3_THRESHOLDS,
    ]
    for p in required:
        require(p.exists(), f"Required frozen v3 artifact does not exist: {p}")

    # Read only frozen master/contract metadata needed to verify identity.
    # Scientific values from the new 05j4c analyses are not accessed here.
    v3_contract = json.loads(V3_CONTRACT.read_text(encoding="utf-8"))
    v3_input_master = json.loads(V3_INPUT_MASTER.read_text(encoding="utf-8"))
    v3_audit_master = json.loads(V3_AUDIT_MASTER.read_text(encoding="utf-8"))

    contract = {
        "script_version": SCRIPT_VERSION,
        "frozen_at_utc": datetime.now(timezone.utc).isoformat(),
        "project": "Paper 4 / TCBB / molecular-transport-audit",
        "purpose": (
            "Complete the already-observed 05j4 v3 scale-invariant audit with "
            "(A) a descriptive held-out false-preservation calibration audit and "
            "(B) the pre-specified missing all-edge Spearman diagnostic."
        ),
        "scientific_provenance": {
            "05j4_v3_results_already_observed": True,
            "05j4c_heldout_fpr_results_observed_before_this_freeze": False,
            "05j4c_all_edge_spearman_results_observed_before_this_freeze": False,
            "primary_mta_statistic_changed": False,
            "primary_mta_classifier_changed": False,
            "v3_thresholds_reestimated": False,
            "v3_results_recomputed_or_replaced": False,
        },
        "locked_existing_hierarchy": {
            "headline_resolution_metric": "f50 with censoring",
            "primary_inferential_contrast": "retention-call AUC",
            "f50_role_in_05j4c": (
                "descriptive effect-size reporting only; no new inferential bootstrap "
                "or threshold is introduced after observing v3"
            ),
            "auc_role_in_05j4c": (
                "existing v3 inferential result remains authoritative; 05j4c does not "
                "replace or re-tune it"
            ),
        },
        "heldout_false_preservation_audit": {
            "status": "FROZEN_BEFORE_05J4C_HELDOUT_RESULTS",
            "full_corruption_valid_mappings_per_pair": 1000,
            "calibration_mappings": 800,
            "heldout_mappings": 200,
            "split_rule": (
                "reuse the exact frozen v3 ordered split: first 800 valid full-corruption "
                "mappings for calibration, final 200 for held-out evaluation"
            ),
            "threshold_rule": (
                "reuse the exact frozen threshold written by 05j4 v3 for every "
                "target-program-method; do not re-estimate, shift, clip, or otherwise "
                "modify thresholds after held-out inspection"
            ),
            "preservation_call_rule": (
                "reuse the v3 score orientation. The audit implementation must fail "
                "rather than guess if the frozen artifacts do not identify the score "
                "and threshold unambiguously."
            ),
            "primary_descriptive_comparison": (
                "paired within target-program pair difference in held-out "
                "false-preservation rate: MTA minus comparator"
            ),
            "comparators": [
                "NetRep_cor.cor",
                "WGCNA_cor.kIM",
                "Pearson_subset_matched",
            ],
            "report": [
                "calibration false-preservation count/rate under the frozen threshold",
                "held-out false-preservation count/rate",
                "held-out minus calibration rate",
                "paired MTA-minus-comparator held-out rate per target-program pair",
                "median, IQR, range, and direction counts of paired rate differences",
            ],
            "interpretation_guard": (
                "A held-out rate above 0.05 for an individual n=200 pair is not itself "
                "declared calibration failure. The diagnostic target is systematic "
                "paired between-method imbalance."
            ),
            "if_systematic_imbalance_is_seen": (
                "report as a limitation/sensitivity finding; do not recalculate thresholds "
                "or introduce a new calibration target after seeing the result"
            ),
            "new_p_values_or_new_inferential_claims": False,
        },
        "all_edge_spearman_diagnostic": {
            "status": "FROZEN_BEFORE_05J4C_ALL_EDGE_SPEARMAN_RESULTS",
            "scientific_question": (
                "complete the 2x2 decomposition of correlation type "
                "(Spearman/Pearson) by edge universe (matched subset/all edges)"
            ),
            "cells": {
                "Spearman_x_matched_subset": "MTA",
                "Pearson_x_matched_subset": "Pearson_subset_matched",
                "Pearson_x_all_edges": "NetRep.cor.cor code path used in 05j4",
                "Spearman_x_all_edges": "05j4c missing diagnostic cell",
            },
            "implementation_definition": (
                "Use the exact same source/target objects, exact same mapped gene ordering, "
                "exact same all-edge eligibility/universe, exact same corruption mappings, "
                "and exact same missing/finite-value handling as the existing all-edge "
                "Pearson cor.cor code path. Change only the correlation across the paired "
                "edge vectors from Pearson to Spearman rank correlation."
            ),
            "source_binding_rule": (
                "Before scientific execution, inspect the existing 05j/05j4b source code "
                "and bind the implementation to that exact all-edge cor.cor path. "
                "No scientifically meaningful fallback implementation is allowed."
            ),
            "reuse": {
                "corruption_grid": [0.0, 0.05, 0.10, 0.15, 0.25, 0.50, 0.75, 1.0],
                "replicates_per_nonzero_performance_fraction": 100,
                "full_corruption_valid_mappings_per_pair": 1000,
                "calibration_mappings": 800,
                "heldout_mappings": 200,
                "existing_mapping_identity_required": True,
                "new_corruption_generation": False,
                "new_target_or_program_exclusions": False,
            },
            "report": [
                "descriptive f50 with the same censoring convention as v3",
                "retention-call AUC under the same calibration/call convention as v3",
                "within-fraction rank correlation with MTA",
                "2x2 descriptive component contrasts",
            ],
            "component_interpretation_guard": (
                "non-detectability or small/unstable decomposition components are valid "
                "outcomes and will not trigger alternative edge subsets, thresholds, "
                "fraction grids, transformations, or post-hoc slicing"
            ),
            "general_superiority_claim": False,
            "existing_v3_primary_inference_changed": False,
        },
        "f50_reporting_guard": {
            "all_22_pairs_were_finite_in_observed_v3": True,
            "planned_descriptive_reporting": [
                "overall median delta-f50",
                "IQR delta-f50",
                "target-stratified median/IQR",
                "direction counts across 22 target-program pairs",
                "translate delta-f50 into percentage points of corrupted correspondence",
            ],
            "new_f50_bootstrap": False,
        },
        "locked_input_artifacts": {
            str(p): {
                "sha256": sha256_file(p),
                "bytes": p.stat().st_size,
            }
            for p in required
        },
        "observed_v3_metadata_snapshot": {
            "v3_contract_script_version": v3_contract.get("script_version"),
            "v3_input_script_version": v3_input_master.get("script_version"),
            "v3_audit_script_version": v3_audit_master.get("script_version"),
        },
        "status": "FROZEN_05J4C_DIAGNOSTIC_COMPLETION_CONTRACT",
    }

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    require(not OUT_JSON.exists(), f"Refusing to overwrite existing frozen contract: {OUT_JSON}")
    OUT_JSON.write_text(json.dumps(contract, indent=2) + "\n", encoding="utf-8")

    print(SEP)
    print("05j4c DIAGNOSTIC COMPLETION CONTRACT: PASS")
    print(SEP)
    print("Existing v3 thresholds reused unchanged:          YES")
    print("Held-out comparison paired within module/cohort: YES")
    print("Held-out threshold re-fitting allowed:           NO")
    print("All-edge Spearman exact mapping reuse required:  YES")
    print("All-edge Spearman post-hoc re-slicing allowed:   NO")
    print("New f50 inferential bootstrap:                   NO")
    print("Primary v3 AUC inference changed:                NO")
    print()
    print(f"Output: {OUT_JSON}")
    print(SEP)

if __name__ == "__main__":
    main()
