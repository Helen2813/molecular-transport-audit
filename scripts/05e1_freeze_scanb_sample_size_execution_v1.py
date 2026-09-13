from __future__ import annotations

import json
from pathlib import Path


SCRIPT_VERSION = "05e1-freeze-scanb-sample-size-execution-v1-no-cli"

DATA_ROOT = Path(r"D:\paper4_tcbb_data")

OPERATING_CONTRACT = (
    DATA_ROOT
    / "paper4_tcbb_preservation_operating_contract_v1"
    / "preservation_operating_characteristics_contract_v1.json"
)
EXACT_VARIATION_CORRECTION = (
    DATA_ROOT
    / "paper4_tcbb_exact_variation_evaluability_correction_v1"
    / "exact_variation_evaluability_correction_v1.json"
)
CORRECTED_DIRECT = (
    DATA_ROOT
    / "paper4_tcbb_primary_pooled_direct_preservation_v2"
    / "primary_pooled_direct_preservation_v2.json"
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
CONTROLLED_DEGRADATION = (
    DATA_ROOT
    / "paper4_tcbb_controlled_mapping_degradation_v1"
    / "controlled_mapping_degradation_v1.json"
)

OUT_DIR = DATA_ROOT / "paper4_tcbb_scanb_sample_size_execution_contract_v1"

SAMPLE_SIZES = [43, 100, 250, 500, 1000, 2000, 3273]
REPEATS_NONFULL = 50
MAX_ATTEMPTS_PER_NONFULL_SIZE = 5000
BASE_SEED = 20260916


def require(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Required file does not exist: {path}")


def main() -> None:
    print("=" * 144)
    print("Paper 4 / TCBB - freeze exact SCAN-B sample-size operating-characteristic execution")
    print("=" * 144)
    print(f"Script version: {SCRIPT_VERSION}")
    print()
    print("No sample-size preservation result is calculated by this contract.")
    print("=" * 144)

    for p in [
        OPERATING_CONTRACT,
        EXACT_VARIATION_CORRECTION,
        CORRECTED_DIRECT,
        CORRECTED_CLASSIFICATION,
        CORRECTED_BOOTSTRAP,
        CONTROLLED_DEGRADATION,
    ]:
        require(p)

    op = json.loads(OPERATING_CONTRACT.read_text(encoding="utf-8"))
    exact = json.loads(EXACT_VARIATION_CORRECTION.read_text(encoding="utf-8"))
    direct = json.loads(CORRECTED_DIRECT.read_text(encoding="utf-8"))
    cls = json.loads(CORRECTED_CLASSIFICATION.read_text(encoding="utf-8"))
    boot = json.loads(CORRECTED_BOOTSTRAP.read_text(encoding="utf-8"))
    degr = json.loads(CONTROLLED_DEGRADATION.read_text(encoding="utf-8"))

    if op.get("status") != "FROZEN_BEFORE_FIRST_TARGET_PRESERVATION_STATISTIC":
        raise RuntimeError("04a operating contract has unexpected status.")
    if exact.get("status") != "FROZEN_CORRECTION_BEFORE_RECOMPUTED_PRESERVATION":
        raise RuntimeError("04i3 exact-variation correction has unexpected status.")
    if direct.get("status") != "CORRECTED_DIRECT_POOLED_PRESERVATION_COMPLETE_MAPPING_NULL_PENDING":
        raise RuntimeError("05a v2 corrected direct result has unexpected status.")
    if cls.get("status") != "CORRECTED_FINAL_PRIMARY_POOLED_CLASSIFICATION_COMPLETE":
        raise RuntimeError("05b v2 corrected classification has unexpected status.")
    if boot.get("status") != "CORRECTED_PRIMARY_POOLED_CLASSICAL_BOOTSTRAP_RELIABILITY_COMPLETE":
        raise RuntimeError("05c v4 corrected bootstrap has unexpected status.")
    if degr.get("status") not in {
        "CONTROLLED_DEGRADATION_ALL_FROZEN_GATES_PASS",
        "CONTROLLED_DEGRADATION_ONE_OR_MORE_FROZEN_GATES_FAIL",
    }:
        raise RuntimeError("05d controlled-degradation result has unexpected status.")

    frozen = op["scanb_sample_size_operating_characteristics"]
    if [int(x) for x in frozen["sample_sizes"]] != SAMPLE_SIZES:
        raise RuntimeError("Sample sizes differ from frozen 04a.")
    if int(frozen["repeats_per_nonfull_size"]) != REPEATS_NONFULL:
        raise RuntimeError("Repeat count differs from frozen 04a.")
    if int(op["computation"]["base_seed"]) != BASE_SEED:
        raise RuntimeError("Base seed differs from frozen 04a.")

    contract = {
        "contract_id": "paper4-tcbb-scanb-sample-size-execution-v1",
        "script_version": SCRIPT_VERSION,
        "status": "FROZEN_BEFORE_SCANB_SAMPLE_SIZE_RESULTS",
        "inherits": {
            "04a": str(OPERATING_CONTRACT),
            "04i3": str(EXACT_VARIATION_CORRECTION),
            "05a_v2": str(CORRECTED_DIRECT),
            "05b_v2": str(CORRECTED_CLASSIFICATION),
            "05c_v4": str(CORRECTED_BOOTSTRAP),
            "05d": str(CONTROLLED_DEGRADATION),
        },
        "sample_sizes": SAMPLE_SIZES,
        "nonfull_repeats_required": REPEATS_NONFULL,
        "sampling": "without replacement from the frozen 3,273 SCAN-B primary biological profiles",
        "full_size": "n=3273 evaluated exactly once",
        "fixed_gene_set": (
            "for each program use exactly the corrected full-SCAN-B evaluable gene set "
            "from 05a v2; do not change genes across sample-size repeats"
        ),
        "small_n_degeneracy_rule": {
            "valid_repeat": (
                "every fixed program gene remains exactly variable in the sampled profiles "
                "(finite values and exact max>min)"
            ),
            "invalid_repeat": (
                "if any fixed program gene becomes constant/nonfinite, discard that sample "
                "subset and draw another subset; do not drop genes, add jitter, add epsilon, "
                "or change the source program"
            ),
            "maximum_attempts_per_nonfull_size_program": MAX_ATTEMPTS_PER_NONFULL_SIZE,
            "if_fewer_than_50_valid_within_cap": (
                "label that sample-size/program operating-characteristic cell NOT ESTIMABLE; "
                "report validity fraction and degenerate-gene counts; do not tune the rule"
            ),
            "interpretation": (
                "effect distributions are conditional on mathematical estimability of the "
                "fixed full-size program; the attempt/validity rate is separately reported "
                "as part of the small-n operating characteristic"
            ),
        },
        "statistics": {
            "edge": (
                "rho_edge on exactly the corrected 05b v2 fixed specificity edge subset "
                "for that program"
            ),
            "loading": (
                "recompute target PC1 on all fixed corrected program genes in the sampled "
                "profiles, orient by the frozen source-sign score, then Spearman rho with "
                "the frozen source loading vector"
            ),
        },
        "full_size_replay_guards": {
            "edge": "n=3273 must reproduce corrected 05b v2 mapping_specificity_rho",
            "loading": "n=3273 must reproduce corrected 05a v2 rho_load",
        },
        "randomization": {
            "base_seed": BASE_SEED,
            "generator": "numpy.random.default_rng",
            "seed_schedule": (
                "SeedSequence([20260916, module_number, sample_size, attempt_id])"
            ),
            "attempt_id": "1-based within program/sample-size",
        },
        "reporting": {
            "per_repeat_edge_rho": True,
            "per_repeat_loading_rho": True,
            "validity_rate": True,
            "invalid_attempt_count": True,
            "degenerate_gene_counts": True,
            "median_q025_q975": True,
            "no_classification_change": True,
            "no_threshold_tuning": True,
        },
        "target_outcomes_or_treatment_loaded": False,
    }

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out = OUT_DIR / "scanb_sample_size_execution_contract_v1.json"
    out.write_text(json.dumps(contract, indent=2), encoding="utf-8")

    print()
    print("=" * 144)
    print("05e1 SCAN-B SAMPLE-SIZE EXECUTION CONTRACT: PASS")
    print("=" * 144)
    print(f"Sample sizes:                    {SAMPLE_SIZES}")
    print(f"Valid repeats/nonfull size:      {REPEATS_NONFULL}")
    print(f"Max attempts/nonfull cell:       {MAX_ATTEMPTS_PER_NONFULL_SIZE}")
    print("Gene set across repeats:          FIXED corrected 05a v2 set")
    print("Undefined subset handling:        discard subset + redraw")
    print("Effect if cap not reached:        NOT ESTIMABLE")
    print("Primary classification changed:  NO")
    print()
    print(f"Output: {out}")
    print("=" * 144)


if __name__ == "__main__":
    main()
