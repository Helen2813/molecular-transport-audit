from __future__ import annotations

import json
from pathlib import Path


SCRIPT_VERSION = "05d1-freeze-controlled-degradation-execution-v1-no-cli"

DATA_ROOT = Path(r"D:\paper4_tcbb_data")

OPERATING_CONTRACT = (
    DATA_ROOT
    / "paper4_tcbb_preservation_operating_contract_v1"
    / "preservation_operating_characteristics_contract_v1.json"
)
MAPPING_CONTRACT = (
    DATA_ROOT
    / "paper4_tcbb_mapping_null_contract_v1"
    / "structure_preserving_mapping_null_contract_v1.json"
)
FINAL_MAPPING_EXECUTION = (
    DATA_ROOT
    / "paper4_tcbb_final_mapping_null_execution_v1"
    / "final_mapping_null_execution_contract_v1.json"
)
EXACT_VARIATION_CORRECTION = (
    DATA_ROOT
    / "paper4_tcbb_exact_variation_evaluability_correction_v1"
    / "exact_variation_evaluability_correction_v1.json"
)
CORRECTED_CLASSIFICATION = (
    DATA_ROOT
    / "paper4_tcbb_final_mapping_specificity_null_v2"
    / "final_mapping_specificity_and_classification_v2.json"
)
CORRECTED_BOOTSTRAP = (
    DATA_ROOT
    / "paper4_tcbb_corrected_classical_bootstrap_reliability_v4"
    / "primary_pooled_corrected_classical_bootstrap_reliability_v4.json"
)

OUT_DIR = DATA_ROOT / "paper4_tcbb_controlled_degradation_execution_contract_v1"

BASE_SEED = 20260916
FRACTIONS = [0.0, 0.10, 0.25, 0.50, 0.75, 1.0]
REPLICATES = 100
MAX_ATTEMPTS_PER_NONZERO_FRACTION = 1000
MAX_EDGES = 20000


def require(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Required file does not exist: {path}")


def main() -> None:
    print("=" * 142)
    print("Paper 4 / TCBB - freeze exact execution details for controlled mapping degradation")
    print("=" * 142)
    print(f"Script version: {SCRIPT_VERSION}")
    print()
    print("This file freezes only execution details left open by 04a.")
    print("No target expression, degradation statistic, or gate result is calculated here.")
    print("=" * 142)

    for p in [
        OPERATING_CONTRACT,
        MAPPING_CONTRACT,
        FINAL_MAPPING_EXECUTION,
        EXACT_VARIATION_CORRECTION,
        CORRECTED_CLASSIFICATION,
        CORRECTED_BOOTSTRAP,
    ]:
        require(p)

    op = json.loads(OPERATING_CONTRACT.read_text(encoding="utf-8"))
    mp = json.loads(MAPPING_CONTRACT.read_text(encoding="utf-8"))
    ex = json.loads(FINAL_MAPPING_EXECUTION.read_text(encoding="utf-8"))
    corr = json.loads(EXACT_VARIATION_CORRECTION.read_text(encoding="utf-8"))
    cls = json.loads(CORRECTED_CLASSIFICATION.read_text(encoding="utf-8"))
    boot = json.loads(CORRECTED_BOOTSTRAP.read_text(encoding="utf-8"))

    if op.get("status") != "FROZEN_BEFORE_FIRST_TARGET_PRESERVATION_STATISTIC":
        raise RuntimeError("04a operating-characteristics contract has unexpected status.")
    if mp.get("status") != "FROZEN_BEFORE_TARGET_MARGINAL_METRICS_AND_PRESERVATION":
        raise RuntimeError("03d mapping contract has unexpected status.")
    if ex.get("status") != "FROZEN_BEFORE_FIRST_TARGET_CORRELATION_OR_PRESERVATION":
        raise RuntimeError("03h mapping execution contract has unexpected status.")
    if corr.get("status") != "FROZEN_CORRECTION_BEFORE_RECOMPUTED_PRESERVATION":
        raise RuntimeError("04i3 exact-variation correction has unexpected status.")
    if cls.get("status") != "CORRECTED_FINAL_PRIMARY_POOLED_CLASSIFICATION_COMPLETE":
        raise RuntimeError("05b v2 corrected classification has unexpected status.")
    if boot.get("status") != "CORRECTED_PRIMARY_POOLED_CLASSICAL_BOOTSTRAP_RELIABILITY_COMPLETE":
        raise RuntimeError("05c v4 corrected bootstrap/reliability has unexpected status.")

    frozen = op["controlled_mapping_degradation"]
    if [float(x) for x in frozen["fractions"]] != FRACTIONS:
        raise RuntimeError("Fractions differ from 04a.")
    if int(frozen["replicates_per_nonzero_fraction"]) != REPLICATES:
        raise RuntimeError("Replicate count differs from 04a.")
    if int(op["computation"]["base_seed"]) != BASE_SEED:
        raise RuntimeError("Base seed differs from 04a.")
    if int(ex["final_edge_specificity_null"]["maximum_edges_per_target_program"]) != MAX_EDGES:
        raise RuntimeError("Maximum degradation edge subset differs from frozen 03h specificity subset.")

    contract = {
        "contract_id": "paper4-tcbb-controlled-degradation-execution-v1",
        "script_version": SCRIPT_VERSION,
        "status": "FROZEN_BEFORE_CONTROLLED_DEGRADATION_RESULTS",
        "inherits": {
            "04a": str(OPERATING_CONTRACT),
            "03d": str(MAPPING_CONTRACT),
            "03h": str(FINAL_MAPPING_EXECUTION),
            "04i3": str(EXACT_VARIATION_CORRECTION),
            "05b_v2": str(CORRECTED_CLASSIFICATION),
            "05c_v4": str(CORRECTED_BOOTSTRAP),
        },
        "fractions": FRACTIONS,
        "replicates_per_nonzero_fraction": REPLICATES,
        "zero_fraction": "one unique correct frozen mapping; no repeated zero-fraction pseudo-replicates",
        "corrupted_slot_count": {
            "formula": "floor(fraction * m + 0.5)",
            "interpretation": "nearest integer with .5 rounded upward",
            "minimum_for_nonzero_fraction": 1,
            "maximum_at_fraction_1": "m",
        },
        "randomization": {
            "base_seed": BASE_SEED,
            "generator": "numpy.random.default_rng",
            "seed_schedule": (
                "numpy SeedSequence entropy vector "
                "[20260916, target_index, module_number, fraction_code, attempt_id]"
            ),
            "target_index": {
                "SCANB_GSE96058": 1,
                "METABRIC": 2,
            },
            "fraction_code": "round(1000*fraction): 100,250,500,750,1000",
            "attempt_id": "1-based attempt index within target/program/fraction",
            "same_rng_use": (
                "within an attempt, the seeded RNG first selects corrupted slots uniformly "
                "without replacement and then executes the frozen 03d matched-background assignment"
            ),
        },
        "matched_background_assignment": {
            "candidate_pool": (
                "corrected target-evaluable frozen-universe genes excluding ALL genes "
                "in the tested frozen source program"
            ),
            "slot_features": (
                "use the frozen 03e source/target marginal percentile features for the "
                "selected corrupted source slots"
            ),
            "generator": (
                "exact 03d rule: randomize corrupted-slot order; sample uniformly without "
                "replacement in the same 5x5 stratum; if empty expand by increasing "
                "Manhattan stratum distance with randomized ties"
            ),
            "no_duplicate_target_gene_assignments": True,
            "uncorrupted_slots": "retain correct target genes",
            "candidate_program_gene_exclusion_guarantees_no_collision_with_uncorrupted_true_program_genes": True,
        },
        "panel_quality": {
            "minimum_fraction_distance_le1": 0.90,
            "maximum_distance": 2,
            "invalid_attempt": "discard and generate another attempt",
            "valid_replicates_required_per_nonzero_fraction": REPLICATES,
            "maximum_attempts_per_nonzero_fraction": MAX_ATTEMPTS_PER_NONZERO_FRACTION,
            "attempt_ratio": "10x, matching the frozen 03d validity-cap ratio",
        },
        "edge_statistic": {
            "source": "reuse the exact corrected 05b v2 per-program source_edge_subset",
            "pair_subset": "reuse exact corrected 05b v2 edge_slot_i/edge_slot_j arrays",
            "maximum_edges": MAX_EDGES,
            "observed_zero_fraction_guard": (
                "recomputed zero-fraction rho must reproduce corrected 05b v2 "
                "mapping_specificity_rho to numerical tolerance"
            ),
        },
        "formal_gate": {
            "scope": (
                "only corrected primary-assessable target/program combinations with "
                "positive uncorrupted edge result and corrected mapping-specificity support"
            ),
            "trend": (
                "Spearman correlation between nominal corruption fractions "
                "[0,.1,.25,.5,.75,1] and median rho_edge at each fraction <= -0.8"
            ),
        },
        "full_corruption_calibration": {
            "frozen_null_reference": "the corrected 05b v2 1000-panel mapping-null rho distribution",
            "median_check": "100%-corruption median must lie within frozen null q2.5-q97.5",
            "nominal_p_per_full_corruption_replicate": (
                "(1 + count(abs(frozen_null_rho) >= abs(rho_rep))) / 1001"
            ),
            "false_specificity_gate": "fraction of 100 full-corruption replicates with nominal p<0.05 <= 0.10",
        },
        "loading": "not a primary controlled-degradation endpoint",
        "target_outcomes_or_treatment_loaded": False,
    }

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out = OUT_DIR / "controlled_degradation_execution_contract_v1.json"
    out.write_text(json.dumps(contract, indent=2), encoding="utf-8")

    print()
    print("=" * 142)
    print("05d1 CONTROLLED-DEGRADATION EXECUTION CONTRACT: PASS")
    print("=" * 142)
    print(f"Fractions:                         {FRACTIONS}")
    print(f"Replicates/nonzero fraction:       {REPLICATES}")
    print(f"Max attempts/nonzero fraction:     {MAX_ATTEMPTS_PER_NONZERO_FRACTION}")
    print("Corrupted-slot rounding:            floor(f*m + 0.5)")
    print("Edge subset:                        reuse exact corrected 05b v2 subset")
    print("100% calibration null:              corrected 05b v2 mapping-null distribution")
    print("Loading degradation endpoint:       NO")
    print("Primary classifications changed:    NO")
    print()
    print(f"Output: {out}")
    print("=" * 142)


if __name__ == "__main__":
    main()
