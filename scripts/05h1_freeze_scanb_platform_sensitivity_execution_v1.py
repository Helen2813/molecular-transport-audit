from __future__ import annotations

import json
from pathlib import Path


SCRIPT_VERSION = "05h1-freeze-scanb-platform-sensitivity-execution-v1-no-cli"

DATA_ROOT = Path(r"D:\paper4_tcbb_data")

PRETARGET_CONTRACT = (
    DATA_ROOT
    / "paper4_tcbb_pretarget_confounds_contract_v1"
    / "pretarget_confounds_contract_v1.json"
)
TECHNICAL_MASTER = (
    DATA_ROOT
    / "paper4_tcbb_scanb_technical_repeatability_v2"
    / "scanb_technical_repeatability_v2.json"
)
CLASSIFICATION = (
    DATA_ROOT
    / "paper4_tcbb_final_mapping_specificity_null_v2"
    / "primary_pooled_final_classification_v2.tsv"
)
PLATFORM_MANIFEST = (
    DATA_ROOT
    / "paper4_tcbb_pam50_alignment_audit_v1"
    / "scanb_frozen_primary_pam50_alignment_v1.tsv"
)
PLATFORM_COMPOSITION = (
    DATA_ROOT
    / "paper4_tcbb_final_pretarget_audits_v1"
    / "scanb_platform_pam50_composition_v1.tsv"
)

OUT_DIR = DATA_ROOT / "paper4_tcbb_scanb_platform_sensitivity_execution_contract_v1"

PLATFORM_COUNTS = {
    "GPL11154": 2969,
    "GPL18573": 304,
}
BOOTSTRAP_REPLICATES = 1000
MAX_ATTEMPTS_PER_CELL = 10000
BASE_SEED = 20260918


def require(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Required file does not exist: {path}")


def main() -> None:
    print("=" * 148)
    print("Paper 4 / TCBB - freeze SCAN-B platform-sensitivity execution details")
    print("=" * 148)
    print(f"Script version: {SCRIPT_VERSION}")
    print()
    print("This freezes only execution details left implicit by 04f.")
    print("No expression matrix or platform-specific preservation result is calculated here.")
    print("=" * 148)

    for p in [
        PRETARGET_CONTRACT,
        TECHNICAL_MASTER,
        CLASSIFICATION,
        PLATFORM_MANIFEST,
        PLATFORM_COMPOSITION,
    ]:
        require(p)

    pre = json.loads(PRETARGET_CONTRACT.read_text(encoding="utf-8"))
    tech = json.loads(TECHNICAL_MASTER.read_text(encoding="utf-8"))

    if pre.get("status") != "FROZEN_BEFORE_FIRST_TARGET_PRESERVATION_STATISTIC":
        raise RuntimeError("04f pretarget-confounds contract has unexpected status.")
    if tech.get("status") != "SCANB_TECHNICAL_REPEATABILITY_COMPLETE":
        raise RuntimeError("05g technical-repeatability master has unexpected status.")

    frozen = pre["scanb_platform_sensitivity"]
    frozen_counts = {
        str(k): int(v)
        for k, v in frozen["platform_groups"].items()
    }
    if frozen_counts != PLATFORM_COUNTS:
        raise RuntimeError(
            f"04f platform counts mismatch: {frozen_counts} vs {PLATFORM_COUNTS}"
        )
    if int(frozen["uncertainty"]["target_sample_bootstraps"]) != BOOTSTRAP_REPLICATES:
        raise RuntimeError("04f platform-bootstrap count mismatch.")
    if int(pre["computation"]["base_seed"]) != BASE_SEED:
        raise RuntimeError("04f base seed mismatch.")

    contract = {
        "contract_id": "paper4-tcbb-scanb-platform-sensitivity-execution-v1",
        "script_version": SCRIPT_VERSION,
        "status": "FROZEN_BEFORE_SCANB_PLATFORM_SENSITIVITY_RESULTS",
        "inherits": {
            "04f": str(PRETARGET_CONTRACT),
            "05g_v2": str(TECHNICAL_MASTER),
        },
        "platform_groups": PLATFORM_COUNTS,
        "platform_manifest": {
            "path": str(PLATFORM_MANIFEST),
            "sample_id_column": "sample_title",
            "platform_column": "platform_id",
            "pam50_column": "pam50_subtype",
            "all_3273_primary_profiles_required": True,
        },
        "source": {
            "cohort": "TCGA_BRCA",
            "samples": 1082,
            "rule": "full frozen TCGA source object remains unchanged",
        },
        "program_gene_set": {
            "rule": (
                "for each program use exactly the corrected full-SCAN-B evaluable gene set "
                "from corrected 05a v2; do not redefine genes within either platform"
            ),
            "platform_specific_gene_dropping": False,
            "exact_variation_rule": (
                "every fixed program gene must be finite and exact-variable (max>min) "
                "within the platform-specific target group"
            ),
            "if_not_estimable": (
                "label target/program/platform NOT_ESTIMABLE_EXACT_VARIATION and report "
                "the exact constant genes; do not drop genes, jitter, epsilon-shift, or "
                "alter module membership"
            ),
        },
        "point_estimates": {
            "rho_edge": (
                "Spearman concordance between full frozen-source Pearson edge vector and "
                "platform-specific target Pearson edge vector on ALL nonredundant fixed-gene edges"
            ),
            "rho_load": (
                "recompute platform-specific target PC1 on all fixed genes, orient by the "
                "frozen source-sign score, and Spearman-correlate with the frozen full-source "
                "loading vector"
            ),
        },
        "uncertainty": {
            "edge_bootstrap_valid_replicates": BOOTSTRAP_REPLICATES,
            "maximum_attempts_per_cell": MAX_ATTEMPTS_PER_CELL,
            "resampling": (
                "ordinary nonparametric target-sample bootstrap with replacement within "
                "the fixed platform group; source object remains fixed"
            ),
            "undefined_bootstrap_draw": (
                "discard and redraw when any fixed target gene becomes nonfinite or "
                "exact-constant in the resampled support; no gene dropping/jitter/epsilon"
            ),
            "edge_subset": (
                "reuse the exact corrected 05b v2 deterministic up-to-20,000-edge subset "
                "for that SCAN-B program"
            ),
            "edge_interval": "95% percentile interval from 1000 valid target bootstraps",
            "loading_bootstrap": False,
            "loading_bootstrap_reason": (
                "04f explicitly prespecified the edge-bootstrap subset but did not specify "
                "a loading-bootstrap algorithm for this secondary platform sensitivity; "
                "rho_load is therefore retained as a descriptive point estimate rather than "
                "adding a post-result inference branch"
            ),
        },
        "randomization": {
            "base_seed": BASE_SEED,
            "generator": "numpy.random.default_rng via SeedSequence",
            "platform_code": {
                "GPL11154": 1,
                "GPL18573": 2,
            },
            "bootstrap_seed_schedule": (
                "SeedSequence([20260918, platform_code, module_number, 101, attempt_id])"
            ),
            "pc1_random_start_schedule": (
                "SeedSequence([20260918, platform_code, module_number, 301])"
            ),
        },
        "cpu_gpu_validation_required": True,
        "additional_reporting": {
            "pam50_composition_file": str(PLATFORM_COMPOSITION),
            "purpose": (
                "report frozen platform-specific PAM50 composition alongside platform effects "
                "to distinguish platform imbalance from subtype-composition imbalance"
            ),
        },
        "formal_new_discovery_p_values": False,
        "changes_primary_classifier": False,
        "target_outcomes_or_treatment_loaded": False,
    }

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out = OUT_DIR / "scanb_platform_sensitivity_execution_contract_v1.json"
    out.write_text(json.dumps(contract, indent=2), encoding="utf-8")

    print()
    print("=" * 148)
    print("05h1 SCAN-B PLATFORM-SENSITIVITY EXECUTION CONTRACT: PASS")
    print("=" * 148)
    print(f"GPL11154 primary profiles:             {PLATFORM_COUNTS['GPL11154']:,}")
    print(f"GPL18573 primary profiles:             {PLATFORM_COUNTS['GPL18573']:,}")
    print(f"Target bootstraps/platform/program:    {BOOTSTRAP_REPLICATES:,}")
    print(f"Maximum attempts/cell:                 {MAX_ATTEMPTS_PER_CELL:,}")
    print("Point rho_edge:                         ALL fixed-gene edges")
    print("Bootstrap edge subset:                  exact corrected 05b v2 subset")
    print("rho_load:                               descriptive point estimate")
    print("Platform-specific gene dropping:        NO")
    print("New discovery p-values:                 NO")
    print("Primary classification changed:         NO")
    print()
    print(f"Output: {out}")
    print("=" * 148)


if __name__ == "__main__":
    main()
