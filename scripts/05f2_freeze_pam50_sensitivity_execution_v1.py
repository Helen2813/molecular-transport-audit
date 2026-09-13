from __future__ import annotations

import json
from pathlib import Path


SCRIPT_VERSION = "05f2-freeze-pam50-sensitivity-execution-v1-no-cli"

DATA_ROOT = Path(r"D:\paper4_tcbb_data")

PAM50_CONTRACT = (
    DATA_ROOT
    / "paper4_tcbb_pam50_sensitivity_contract_v1"
    / "pam50_sensitivity_contract_v1.json"
)
ALIGNMENT_AUDIT = (
    DATA_ROOT
    / "paper4_tcbb_pam50_alignment_audit_v1"
    / "pam50_alignment_audit_v1.json"
)
EXACT_CORRECTION = (
    DATA_ROOT
    / "paper4_tcbb_exact_variation_evaluability_correction_v1"
    / "exact_variation_evaluability_correction_v1.json"
)
FINAL_CLASSIFICATION = (
    DATA_ROOT
    / "paper4_tcbb_final_mapping_specificity_null_v2"
    / "final_mapping_specificity_and_classification_v2.json"
)
BOOTSTRAP_RESULT = (
    DATA_ROOT
    / "paper4_tcbb_corrected_classical_bootstrap_reliability_v4"
    / "primary_pooled_corrected_classical_bootstrap_reliability_v4.json"
)
DEGRADATION_RESULT = (
    DATA_ROOT
    / "paper4_tcbb_controlled_mapping_degradation_v1"
    / "controlled_mapping_degradation_v1.json"
)
SAMPLE_SIZE_RESULT = (
    DATA_ROOT
    / "paper4_tcbb_scanb_sample_size_operating_characteristics_v1"
    / "scanb_sample_size_operating_characteristics_v1.json"
)

OUT_DIR = DATA_ROOT / "paper4_tcbb_pam50_sensitivity_execution_contract_v1"

BASE_SEED = 20260917
BOOTSTRAP_REPLICATES = 1000
MAX_BOOTSTRAP_ATTEMPTS = 10000
MAX_BOOTSTRAP_EDGES = 20000
CORE = ["LumA", "LumB", "Basal"]
LIMITED_N = ["Her2"]
NOT_INDIVIDUALLY_ASSESSABLE = ["Normal"]


def require(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Required file does not exist: {path}")


def main() -> None:
    print("=" * 146)
    print("Paper 4 / TCBB - freeze PAM50 subtype/composition execution details")
    print("=" * 146)
    print(f"Script version: {SCRIPT_VERSION}")
    print()
    print("This step freezes implementation details not already explicit in 04e.")
    print("No expression matrix, subtype preservation effect, or composition effect is calculated here.")
    print("=" * 146)

    for p in [
        PAM50_CONTRACT,
        ALIGNMENT_AUDIT,
        EXACT_CORRECTION,
        FINAL_CLASSIFICATION,
        BOOTSTRAP_RESULT,
        DEGRADATION_RESULT,
        SAMPLE_SIZE_RESULT,
    ]:
        require(p)

    pam = json.loads(PAM50_CONTRACT.read_text(encoding="utf-8"))
    ali = json.loads(ALIGNMENT_AUDIT.read_text(encoding="utf-8"))
    exact = json.loads(EXACT_CORRECTION.read_text(encoding="utf-8"))
    cls = json.loads(FINAL_CLASSIFICATION.read_text(encoding="utf-8"))
    boot = json.loads(BOOTSTRAP_RESULT.read_text(encoding="utf-8"))
    degr = json.loads(DEGRADATION_RESULT.read_text(encoding="utf-8"))
    ssize = json.loads(SAMPLE_SIZE_RESULT.read_text(encoding="utf-8"))

    if pam.get("status") != "FROZEN_BEFORE_FIRST_TARGET_PRESERVATION_STATISTIC":
        raise RuntimeError("04e PAM50 contract has unexpected status.")
    if ali.get("status") != "ALIGNMENT_AUDIT_ONLY_NOT_YET_SUBTYPE_ANALYSIS_CONTRACT":
        raise RuntimeError("04d PAM50 alignment audit has unexpected status.")
    if exact.get("status") != "FROZEN_CORRECTION_BEFORE_RECOMPUTED_PRESERVATION":
        raise RuntimeError("04i3 exact-variation correction has unexpected status.")
    if cls.get("status") != "CORRECTED_FINAL_PRIMARY_POOLED_CLASSIFICATION_COMPLETE":
        raise RuntimeError("05b v2 classification has unexpected status.")
    if boot.get("status") != "CORRECTED_PRIMARY_POOLED_CLASSICAL_BOOTSTRAP_RELIABILITY_COMPLETE":
        raise RuntimeError("05c v4 bootstrap/reliability has unexpected status.")
    if degr.get("status") not in {
        "CONTROLLED_DEGRADATION_ALL_FROZEN_GATES_PASS",
        "CONTROLLED_DEGRADATION_ONE_OR_MORE_FROZEN_GATES_FAIL",
    }:
        raise RuntimeError("05d degradation result has unexpected status.")
    if ssize.get("status") != "SCANB_SAMPLE_SIZE_OPERATING_CHARACTERISTICS_COMPLETE":
        raise RuntimeError("05e sample-size result has unexpected status.")

    if int(pam["computation"]["base_seed"]) != BASE_SEED:
        raise RuntimeError("04e base seed is not 20260917.")
    if int(pam["subtype_stratified_sensitivity"]["bootstrap"]["replicates"]) != BOOTSTRAP_REPLICATES:
        raise RuntimeError("Subtype bootstrap count differs from 04e.")
    if int(pam["composition_residualization_sensitivity"]["uncertainty"]["bootstrap_replicates"]) != BOOTSTRAP_REPLICATES:
        raise RuntimeError("Composition bootstrap count differs from 04e.")

    tiers = pam["subtype_reporting_tiers"]["expected_result_from_counts_only"]
    if tiers["core"] != CORE or tiers["limited_n"] != LIMITED_N or tiers["not_individually_assessable"] != NOT_INDIVIDUALLY_ASSESSABLE:
        raise RuntimeError("Frozen PAM50 reporting tiers differ from expected 04e values.")

    contract = {
        "contract_id": "paper4-tcbb-pam50-sensitivity-execution-v1",
        "script_version": SCRIPT_VERSION,
        "status": "FROZEN_BEFORE_PAM50_SENSITIVITY_RESULTS",
        "inherits": {
            "04e": str(PAM50_CONTRACT),
            "04d": str(ALIGNMENT_AUDIT),
            "04i3": str(EXACT_CORRECTION),
            "05b_v2": str(FINAL_CLASSIFICATION),
            "05c_v4": str(BOOTSTRAP_RESULT),
            "05d": str(DEGRADATION_RESULT),
            "05e": str(SAMPLE_SIZE_RESULT),
        },
        "reporting_tiers": {
            "CORE": CORE,
            "LIMITED_N": LIMITED_N,
            "NOT_INDIVIDUALLY_ASSESSABLE": NOT_INDIVIDUALLY_ASSESSABLE,
        },
        "fixed_gene_set_rule": (
            "for each target/program use exactly the corrected full-target evaluable gene set from 05a v2; "
            "do not drop or replace genes within subtype or composition sensitivity analyses"
        ),
        "exact_variation_execution_guard": {
            "point_estimate": (
                "if any fixed program gene is nonfinite or exactly constant in a required source/target subtype, "
                "complete-case cohort, or residualized cohort, that specific sensitivity cell is NOT_ESTIMABLE; "
                "do not drop genes, add jitter, add epsilon, or change module membership"
            ),
            "bootstrap_attempt": (
                "if a source or target bootstrap resample makes any fixed gene nonfinite/exactly constant, "
                "or makes a required residualized gene exactly constant, discard that bootstrap attempt and redraw"
            ),
            "valid_replicates_required": BOOTSTRAP_REPLICATES,
            "maximum_attempts_per_cell": MAX_BOOTSTRAP_ATTEMPTS,
            "if_cap_not_reached": "retain the point estimate but label the bootstrap CI NOT_ESTIMABLE",
        },
        "bootstrap_edge_subset": {
            "rule": (
                "reuse the exact corrected 05b v2 per-target/program specificity edge subset. "
                "It is already a deterministic uniform up-to-20,000-edge subset frozen before target results, "
                "and therefore satisfies the 04e deterministic bootstrap-edge requirement without introducing a new post-result pair sample"
            ),
            "maximum_edges": MAX_BOOTSTRAP_EDGES,
            "point_estimates": "still use all nonredundant fixed-gene edges as required by 04e",
        },
        "subtype_bootstrap": {
            "replicates": BOOTSTRAP_REPLICATES,
            "source_and_target_resampling": "independent ordinary sample-with-replacement bootstrap within the same PAM50 subtype",
            "ci": "95% percentile interval of edge rho",
            "loading_bootstrap": False,
            "formal_p_values": False,
            "normal_subtype_computed": False,
            "normal_reason": "04e tier is NOT_INDIVIDUALLY_ASSESSABLE",
        },
        "composition_execution": {
            "A": "corrected 05a v2 full-cohort all-edge rho_edge",
            "B": "PAM50-complete-case pooled all-edge rho_edge on identical fixed genes",
            "C": "same complete cases after within-cohort PAM50 mean residualization and gene-wise residual z-standardization",
            "METABRIC_complete_case_labels": ["LumA", "LumB", "Basal", "Her2", "Normal"],
            "METABRIC_excluded": ["Claudin-low", "NC"],
            "bootstrap": (
                "paired B/C stratified bootstrap: independently resample source and target WITHIN each canonical PAM50 subtype, "
                "preserving the frozen subtype sample counts; refit B correlations and C subtype means/correlations in every replicate"
            ),
            "bootstrap_outputs": ["B rho", "C rho", "C minus B"],
            "B_minus_A_uncertainty": (
                "point estimate only in this block because A is the separate full-cohort primary estimand with its own 05c uncertainty "
                "and does not share the PAM50-complete-case bootstrap sample base"
            ),
            "no_robustness_retention_threshold": True,
        },
        "randomization": {
            "base_seed": BASE_SEED,
            "generator": "numpy.random.default_rng with numpy SeedSequence entropy vectors",
            "target_index": {"SCANB_GSE96058": 1, "METABRIC": 2},
            "subtype_code": {"LumA": 1, "LumB": 2, "Basal": 3, "Her2": 4, "Normal": 5},
            "subtype_bootstrap_attempt": (
                "SeedSequence([20260917, target_index, module_number, subtype_code, 101, attempt_id])"
            ),
            "composition_bootstrap_attempt": (
                "SeedSequence([20260917, target_index, module_number, 201, attempt_id])"
            ),
            "pc1_random_start": (
                "SeedSequence([20260917, target_index, module_number, subtype_code, 301, cohort_code]), cohort_code 1=source, 2=target"
            ),
        },
        "cpu_reference_validation": {
            "point_full_correlation": "validate a fixed leading gene block CPU versus GPU for every computed full correlation matrix",
            "bootstrap_selected_edges": "on the first valid bootstrap attempt in each cell, compare a fixed leading selected-edge block against explicit CPU resampling",
        },
        "claim_guard": {
            "new_discovery_q_values": False,
            "changes_primary_classification": False,
            "posthoc_threshold_tuning": False,
        },
        "target_outcomes_or_treatment_loaded": False,
    }

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out = OUT_DIR / "pam50_sensitivity_execution_contract_v1.json"
    out.write_text(json.dumps(contract, indent=2), encoding="utf-8")

    print()
    print("=" * 146)
    print("05f2 PAM50 SENSITIVITY EXECUTION CONTRACT: PASS")
    print("=" * 146)
    print(f"Subtype bootstrap valid replicates:     {BOOTSTRAP_REPLICATES:,}")
    print(f"Composition bootstrap valid replicates: {BOOTSTRAP_REPLICATES:,}")
    print(f"Max attempts/cell:                      {MAX_BOOTSTRAP_ATTEMPTS:,}")
    print("Bootstrap edge subset:                  exact corrected 05b v2 subset")
    print("Point estimates:                        ALL fixed-gene edges")
    print("Normal standalone subtype conclusion:   NO")
    print("Primary classification changed:         NO")
    print()
    print(f"Output: {out}")
    print("=" * 146)


if __name__ == "__main__":
    main()
