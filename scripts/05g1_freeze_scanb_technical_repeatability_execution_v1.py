from __future__ import annotations

import json
from pathlib import Path


SCRIPT_VERSION = "05g1-freeze-scanb-technical-repeatability-execution-v1-no-cli"

DATA_ROOT = Path(r"D:\paper4_tcbb_data")

OPERATING_CONTRACT = (
    DATA_ROOT
    / "paper4_tcbb_preservation_operating_contract_v1"
    / "preservation_operating_characteristics_contract_v1.json"
)
PAIRING_CONTRACT = (
    DATA_ROOT
    / "paper4_tcbb_scanb_pairing_contract_v1"
    / "scanb_pairing_contract_v1.json"
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
SAMPLE_SIZE_MASTER = (
    DATA_ROOT
    / "paper4_tcbb_scanb_sample_size_operating_characteristics_v1"
    / "scanb_sample_size_operating_characteristics_v1.json"
)
PAM50_MASTER = (
    DATA_ROOT
    / "paper4_tcbb_pam50_sensitivity_v1"
    / "pam50_subtype_composition_sensitivity_v1.json"
)

PAIR_MANIFEST = (
    DATA_ROOT
    / "paper4_tcbb_scanb_pairing_contract_v1"
    / "scanb_technical_replicate_pairs_v1.tsv"
)
SAMPLE_SIZE_REPLICATES = (
    DATA_ROOT
    / "paper4_tcbb_scanb_sample_size_operating_characteristics_v1"
    / "scanb_sample_size_replicates_v1.tsv"
)

OUT_DIR = (
    DATA_ROOT
    / "paper4_tcbb_scanb_technical_repeatability_execution_contract_v1"
)

SAME_PLATFORM_PAIRS = 100
CROSS_PLATFORM_PAIRS = 36


def require(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Required file does not exist: {path}")


def main() -> None:
    print("=" * 146)
    print("Paper 4 / TCBB - freeze SCAN-B technical-repeatability execution details")
    print("=" * 146)
    print(f"Script version: {SCRIPT_VERSION}")
    print()
    print("No expression values or technical-repeatability statistic are calculated here.")
    print("=" * 146)

    for p in [
        OPERATING_CONTRACT,
        PAIRING_CONTRACT,
        EXACT_VARIATION_CORRECTION,
        CORRECTED_CLASSIFICATION,
        SAMPLE_SIZE_MASTER,
        PAM50_MASTER,
        PAIR_MANIFEST,
        SAMPLE_SIZE_REPLICATES,
    ]:
        require(p)

    op = json.loads(OPERATING_CONTRACT.read_text(encoding="utf-8"))
    pairing = json.loads(PAIRING_CONTRACT.read_text(encoding="utf-8"))
    exact = json.loads(EXACT_VARIATION_CORRECTION.read_text(encoding="utf-8"))
    cls = json.loads(CORRECTED_CLASSIFICATION.read_text(encoding="utf-8"))
    ss = json.loads(SAMPLE_SIZE_MASTER.read_text(encoding="utf-8"))
    pam = json.loads(PAM50_MASTER.read_text(encoding="utf-8"))

    if op.get("status") != "FROZEN_BEFORE_FIRST_TARGET_PRESERVATION_STATISTIC":
        raise RuntimeError("04a operating-characteristics contract has unexpected status.")
    if pairing.get("status") != "FROZEN":
        raise RuntimeError("01c SCAN-B pairing contract has unexpected status.")
    if exact.get("status") != "FROZEN_CORRECTION_BEFORE_RECOMPUTED_PRESERVATION":
        raise RuntimeError("04i3 exact-variation correction has unexpected status.")
    if cls.get("status") != "CORRECTED_FINAL_PRIMARY_POOLED_CLASSIFICATION_COMPLETE":
        raise RuntimeError("05b v2 corrected classification has unexpected status.")
    if ss.get("status") != "SCANB_SAMPLE_SIZE_OPERATING_CHARACTERISTICS_COMPLETE":
        raise RuntimeError("05e sample-size master has unexpected status.")
    if pam.get("status") != "PAM50_SUBTYPE_COMPOSITION_SENSITIVITY_COMPLETE":
        raise RuntimeError("05f PAM50 master has unexpected status.")

    frozen = op["scanb_technical_repeatability"]
    if int(frozen["same_platform_primary_ceiling"]["pairs"]) != SAME_PLATFORM_PAIRS:
        raise RuntimeError("Same-platform pair count differs from frozen 04a.")
    if int(frozen["cross_platform_secondary"]["pairs"]) != CROSS_PLATFORM_PAIRS:
        raise RuntimeError("Cross-platform pair count differs from frozen 04a.")

    contract = {
        "contract_id": "paper4-tcbb-scanb-technical-repeatability-execution-v1",
        "script_version": SCRIPT_VERSION,
        "status": "FROZEN_BEFORE_SCANB_TECHNICAL_REPEATABILITY_RESULTS",
        "inherits": {
            "04a": str(OPERATING_CONTRACT),
            "01c": str(PAIRING_CONTRACT),
            "04i3": str(EXACT_VARIATION_CORRECTION),
            "05b_v2": str(CORRECTED_CLASSIFICATION),
            "05e": str(SAMPLE_SIZE_MASTER),
            "05f": str(PAM50_MASTER),
        },
        "pair_manifest": {
            "path": str(PAIR_MANIFEST),
            "pair_id": "pair_id",
            "primary_title": "primary_title",
            "replicate_title": "replicate_title",
            "primary_platform": "primary_platform_id",
            "replicate_platform": "replicate_platform_id",
            "same_platform_flag": "same_platform",
            "cross_platform_flag": "cross_platform",
            "same_platform_pairs_required": SAME_PLATFORM_PAIRS,
            "cross_platform_pairs_required": CROSS_PLATFORM_PAIRS,
            "pairing_is_frozen_and_not_reselected": True,
        },
        "same_platform_primary_ceiling": {
            "pairs": SAME_PLATFORM_PAIRS,
            "use": "all frozen same-platform pairs exactly once",
            "primary_matrix": "100 paired primary profiles",
            "replicate_matrix": "the corresponding 100 technical replicate profiles",
            "statistic": (
                "Spearman concordance between primary-profile and replicate-profile "
                "gene-gene Pearson edge vectors"
            ),
            "edge_subset": (
                "exact corrected 05b v2 fixed specificity edge subset for the same "
                "target/program; this is required for a like-for-like comparison with "
                "the frozen 05e SCAN-B n=100 transport distribution"
            ),
            "role": "measurement ceiling at matched n, not a biological null",
        },
        "cross_platform_secondary": {
            "pairs": CROSS_PLATFORM_PAIRS,
            "use": "all frozen cross-platform pairs exactly once",
            "statistic": (
                "same edge-concordance statistic on the same corrected 05b v2 edge subset"
            ),
            "role": "secondary platform/resequencing sensitivity; no formal gate",
        },
        "fixed_gene_evaluability": {
            "gene_set": (
                "use exactly the corrected SCAN-B evaluable program genes from 05a v2; "
                "no technical-set-specific gene dropping"
            ),
            "estimability_rule": (
                "for a technical comparison, every fixed program gene must be finite and "
                "exactly variable (max>min) in BOTH the primary and replicate matrices"
            ),
            "if_not_estimable": (
                "label that target/program technical comparison NOT_ESTIMABLE_EXACT_VARIATION; "
                "report the exact degenerate genes separately; do not drop genes, jitter, "
                "epsilon-shift, or alter the frozen edge subset"
            ),
        },
        "matched_n_transport_comparison": {
            "reference": str(SAMPLE_SIZE_REPLICATES),
            "transport_sample_size": 100,
            "reference_repeats": 50,
            "same_statistic": (
                "05e edge rho on the same corrected 05b v2 fixed specificity edge subset"
            ),
            "report": [
                "technical same-platform rho",
                "05e n=100 transport median",
                "05e n=100 transport 2.5th-97.5th percentile interval",
                "technical minus transport median",
                "empirical percentile rank of technical rho among the 50 transport repeats",
            ],
            "formal_p_value": False,
            "formal_gate": False,
            "interpretation": (
                "descriptive matched-n measurement-ceiling comparison only; "
                "technical repeatability is not used to reclassify biological transport"
            ),
        },
        "uncertainty": {
            "technical_pair_bootstrap": False,
            "reason": (
                "04a froze the all-pairs ceiling statistic but did not prespecify a "
                "technical-pair bootstrap; no post-result uncertainty branch is introduced"
            ),
        },
        "loading_axis": {
            "computed": False,
            "reason": "04a technical-repeatability endpoint is edge concordance only",
        },
        "formal_new_p_values": False,
        "changes_primary_classification": False,
        "target_outcomes_or_treatment_loaded": False,
    }

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out = OUT_DIR / "scanb_technical_repeatability_execution_contract_v1.json"
    out.write_text(json.dumps(contract, indent=2), encoding="utf-8")

    print()
    print("=" * 146)
    print("05g1 SCAN-B TECHNICAL-REPEATABILITY EXECUTION CONTRACT: PASS")
    print("=" * 146)
    print(f"Same-platform pairs:                 {SAME_PLATFORM_PAIRS}")
    print(f"Cross-platform pairs:                {CROSS_PLATFORM_PAIRS}")
    print("Pair selection:                       frozen 01c manifest; no reselection")
    print("Technical edge subset:                exact corrected 05b v2 subset")
    print("Matched-n transport reference:        exact 05e n=100 50-repeat distribution")
    print("Technical-set-specific gene dropping: NO")
    print("Technical pair bootstrap:             NO (not prespecified)")
    print("Primary classification changed:       NO")
    print()
    print(f"Output: {out}")
    print("=" * 146)


if __name__ == "__main__":
    main()
