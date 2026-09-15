#!/usr/bin/env python
from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

SCRIPT_VERSION = "05j4c-freeze-all-edge-spearman-interpretation-v1-no-cli"

PROJECT_ROOT = Path(r"C:\Users\olegk\Desktop\molecular-transport-audit")
DATA_ROOT = Path(r"D:\paper4_tcbb_data")

SOURCE_PREP = PROJECT_ROOT / "scripts" / "05j4b_prepare_scale_invariant_inputs_gpu_v2.py"
SOURCE_AUDIT = PROJECT_ROOT / "scripts" / "05j4_run_scale_invariant_corruption_audit_v3.py"
SOURCE_05J = PROJECT_ROOT / "scripts" / "05j_run_structural_corruption_headtohead_gpu_v1.py"

COMPLETION_CONTRACT = (
    DATA_ROOT
    / "paper4_tcbb_scale_invariant_diagnostic_completion_contract_v1"
    / "scale_invariant_diagnostic_completion_contract_v1.json"
)
PROBE_JSON = (
    DATA_ROOT
    / "paper4_tcbb_scale_invariant_diagnostic_completion_contract_v1"
    / "all_edge_spearman_implementation_probe_v1.json"
)

OUT_DIR = DATA_ROOT / "paper4_tcbb_scale_invariant_diagnostic_completion_contract_v1"
OUT_JSON = OUT_DIR / "all_edge_spearman_interpretation_contract_v1.json"

EXPECTED_SHA256 = {
    str(SOURCE_PREP): "64580c8b1aa3f2e7382c9c4fc7433ca0e0af87b26c6c56f7f2f9789ea1d7d5e4",
    str(SOURCE_AUDIT): "37e2286612cf7a708ef99e6de5667073cc8a78314ba5bbc6402bca86cc24792f",
    str(SOURCE_05J): "e25583122628596eb52d99a7b837031564a5c1764c0e6e371dc2743c7c4ae2e1",
}

SEP = "=" * 168

def require(cond: bool, msg: str) -> None:
    if not cond:
        raise RuntimeError(msg)

def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()

def main() -> None:
    print(SEP)
    print("Paper 4 / TCBB - freeze all-edge Spearman implementation + interpretation")
    print(SEP)
    print(f"Script version: {SCRIPT_VERSION}")
    print()
    print("Provenance:")
    print("  05j4 v3 scientific results already observed:       YES")
    print("  held-out calibration audit already observed:       YES")
    print("  source-only implementation probe already observed: YES")
    print("  Spearman_all_edges scientific results observed:    NO")
    print("  2x2 interpretation frozen before its result:       YES")
    print()

    for p in [SOURCE_PREP, SOURCE_AUDIT, SOURCE_05J, COMPLETION_CONTRACT, PROBE_JSON]:
        require(p.exists(), f"Required file missing: {p}")

    actual_hashes = {}
    for p in [SOURCE_PREP, SOURCE_AUDIT, SOURCE_05J]:
        actual = sha256_file(p)
        expected = EXPECTED_SHA256[str(p)]
        require(actual == expected, f"Historical source hash drift: {p}\nexpected={expected}\nactual={actual}")
        actual_hashes[str(p)] = actual

    completion = json.loads(COMPLETION_CONTRACT.read_text(encoding="utf-8"))
    require(
        completion.get("status") == "FROZEN_05J4C_DIAGNOSTIC_COMPLETION_CONTRACT",
        "Unexpected completion-contract status.",
    )

    contract = {
        "script_version": SCRIPT_VERSION,
        "frozen_at_utc": datetime.now(timezone.utc).isoformat(),
        "project": "Paper 4 / TCBB / molecular-transport-audit",
        "provenance": {
            "05j4_v3_results_observed": True,
            "heldout_calibration_results_observed": True,
            "source_only_probe_observed": True,
            "spearman_all_edges_results_observed_before_freeze": False,
        },
        "exact_all_edge_cell_definition": {
            "source_gene_set": "exact same slot_genes as 05j4b/NetRep cor.cor",
            "target_mapped_gene_set": "exact same mapped indices reconstructed from saved attempt_id",
            "edge_universe": "all strict upper-triangle within-module edges via np.triu_indices(m, k=1)",
            "source_vector": "source_corr[tri_i, tri_j]",
            "target_vector": "mapped target sub-correlation matrix [tri_i, tri_j]",
            "only_scientific_change_from_netrep_cor_cor": "Pearson correlation -> Spearman rank correlation",
            "rank_ties": "average ranks, validated against scipy.stats.rankdata(method='average')",
            "mapping_generation_changed": False,
            "edge_filtering_changed": False,
            "missing_value_rule_changed": False,
        },
        "threshold_and_call_rule": {
            "calibration_n": 800,
            "holdout_n": 200,
            "alpha": 0.05,
            "threshold": "sort finite calibration scores; k=floor(alpha*n); threshold=xs[n-k-1]",
            "preserved_call": "score > threshold",
            "guarantee_on_calibration_set": "false-preservation rate <= 0.05",
            "threshold_retuning_after_result": False,
        },
        "corruption_reuse": {
            "fractions": [0.0, 0.05, 0.10, 0.15, 0.25, 0.50, 0.75, 1.0],
            "performance_rows_per_pair": 601,
            "full_null_rows_per_pair": 1000,
            "target_program_pairs": 22,
            "reuse_saved_attempt_ids": True,
            "new_mapping_generation": False,
        },
        "pre_result_2x2_interpretation": {
            "cell_labels": {
                "A": "MTA = Spearman x matched/frozen subset",
                "B": "Pearson_subset_matched = Pearson x matched/frozen subset",
                "C": "Spearman_all_edges = Spearman x all edges (new cell)",
                "D": "NetRep_cor.cor = Pearson x all edges",
            },
            "observed_reference_effects_before_new_cell": {
                "A_minus_D_median_delta_f50_approx": -0.0388,
                "A_minus_B_median_delta_f50_approx": -0.0149,
            },
            "one_factor_effects_to_report_pairwise": {
                "correlation_type_at_subset": "A - B",
                "correlation_type_at_all_edges": "C - D",
                "edge_universe_under_spearman": "A - C",
                "edge_universe_under_pearson": "B - D",
                "interaction_difference_of_differences": "(A - B) - (C - D)",
            },
            "outcome_A_edge_universe_small": (
                "If A-C is near zero, changing from the frozen matched subset to all "
                "edges has little effect under Spearman. The remaining A-D gap must "
                "then arise from correlation-type dependence on the all-edge universe "
                "and/or interaction; do not attribute it to edge selection."
            ),
            "outcome_B_similar_one_factor_effects": (
                "If A-C is close to the already-observed A-B (~-0.015), the direct "
                "one-factor edge-universe effect under Spearman and correlation-type "
                "effect on the matched subset are of similar magnitude. The larger "
                "combined A-D gap (~-0.039) can be described only after inspecting the "
                "pairwise additive/interaction decomposition."
            ),
            "outcome_C_edge_universe_dominant": (
                "If A-C is close to the already-observed A-D (~-0.039), changing the "
                "edge universe alone nearly reproduces the full MTA-vs-NetRep gap; "
                "this supports edge-universe selection as the dominant contributor, "
                "with C-D expected to be comparatively small if the decomposition is coherent."
            ),
            "outcome_D_mixed_or_unresolved": (
                "If A-C is intermediate, heterogeneous, or the difference-of-differences "
                "is unstable, report a mixed or unresolved decomposition. Do not search "
                "alternative subsets, thresholds, corruption grids, or post-hoc strata."
            ),
            "aggregate_medians_not_assumed_additive": True,
            "pairwise_decomposition_preferred": True,
            "no_general_superiority_claim": True,
            "non_detectability_is_valid_outcome": True,
        },
        "discussion_mechanism_status": {
            "signal_dilution_hypothesis": (
                "A fixed informative edge subset may react earlier because the all-edge "
                "universe can contain many weak/noisy edges that dilute mapping-error signal."
            ),
            "status": "plausible mechanism to discuss, not established causal mechanism",
            "wording_guard": "describe as a hypothesis/possible explanation unless directly tested",
        },
        "historical_source_sha256": actual_hashes,
        "status": "FROZEN_BEFORE_ALL_EDGE_SPEARMAN_RESULTS",
    }

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    require(not OUT_JSON.exists(), f"Refusing to overwrite existing contract: {OUT_JSON}")
    OUT_JSON.write_text(json.dumps(contract, indent=2) + "\n", encoding="utf-8")

    print(SEP)
    print("05j4c ALL-EDGE SPEARMAN INTERPRETATION CONTRACT: PASS")
    print(SEP)
    print("Exact all-edge universe locked:             YES")
    print("Only Pearson -> Spearman change permitted:  YES")
    print("Saved attempt IDs/mappings reused:          YES")
    print("Strict score > threshold convention locked: YES")
    print("Three interpretation branches pre-specified:YES")
    print("Post-hoc re-slicing permitted:              NO")
    print()
    print(f"Output: {OUT_JSON}")
    print(SEP)

if __name__ == "__main__":
    main()
