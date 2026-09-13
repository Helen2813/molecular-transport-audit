from __future__ import annotations

import json
from pathlib import Path

SCRIPT_VERSION = "04f-freeze-pretarget-independence-platform-source-stability-contract-v1-no-cli"

DATA_ROOT = Path(r"D:\paper4_tcbb_data")
PAM50_CONTRACT = (
    DATA_ROOT
    / "paper4_tcbb_pam50_sensitivity_contract_v1"
    / "pam50_sensitivity_contract_v1.json"
)
GPU_AUDIT = (
    DATA_ROOT
    / "paper4_tcbb_gpu_audit_v1"
    / "gpu_compute_stack_audit_v1.json"
)
OUT_DIR = DATA_ROOT / "paper4_tcbb_pretarget_confounds_contract_v1"

# Cross-cohort identity audit.
TECH_NULL_NONMATCH_PAIRS = 100_000
TECH_TRUE_PERCENTILE = 0.01
TECH_NULL_PERCENTILE = 0.9999
IDENTITY_MIN_MARGIN_TO_SECOND_BEST = 0.02

# Source-object stability.
SOURCE_EDGE_BOOTSTRAPS = 500
SOURCE_LOADING_BOOTSTRAPS = 200
SOURCE_STABILITY_MAX_EDGES = 20_000

# SCAN-B platform sensitivity.
PLATFORM_BOOTSTRAPS = 1_000

BASE_SEED = 20260918


def require(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Required file does not exist: {path}")


def main() -> None:
    print("=" * 134)
    print("Paper 4 / TCBB - freeze final pre-target independence/platform/source-stability contract")
    print("=" * 134)
    print(f"Script version: {SCRIPT_VERSION}")
    print()
    print("Scientific guard:")
    print("  Target outcomes loaded:                              NO")
    print("  Target treatment variables loaded:                   NO")
    print("  Target gene-gene preservation calculated:            NO")
    print("  Target PCA/loading preservation calculated:           NO")
    print("  Primary classifier evaluated:                         NO")
    print("=" * 134)

    require(PAM50_CONTRACT)
    require(GPU_AUDIT)

    pam50 = json.loads(PAM50_CONTRACT.read_text(encoding="utf-8"))
    gpu = json.loads(GPU_AUDIT.read_text(encoding="utf-8"))

    if pam50.get("status") != "FROZEN_BEFORE_FIRST_TARGET_PRESERVATION_STATISTIC":
        raise RuntimeError("04e PAM50 contract is not in the expected frozen state.")
    if not bool(gpu.get("cuda_available", False)):
        raise RuntimeError("GPU audit does not report CUDA available.")

    contract = {
        "contract_id": "paper4-tcbb-pretarget-confounds-v1",
        "script_version": SCRIPT_VERSION,
        "status": "FROZEN_BEFORE_FIRST_TARGET_PRESERVATION_STATISTIC",
        "scientific_guard": {
            "target_outcomes_loaded": False,
            "target_treatment_loaded": False,
            "target_gene_gene_preservation_calculated": False,
            "target_pca_loading_preservation_calculated": False,
            "primary_classifier_evaluated": False,
        },

        "cross_cohort_identity_audit": {
            "role": (
                "rule out exact/near-duplicate biological specimens across nominally independent "
                "TCGA-BRCA, SCAN-B, and METABRIC cohorts before transport results are interpreted"
            ),
            "gene_universe": (
                "pairwise intersection of the frozen TCGA 10,000-gene source universe and genes "
                "available in both cohorts being compared"
            ),
            "expression_scale": (
                "use the already frozen transformed expression scale in each cohort; then, within "
                "each cohort separately, center and scale each gene across biological samples"
            ),
            "sample_fingerprint": (
                "Pearson correlation across the common standardized genes between cross-cohort sample profiles"
            ),
            "threshold_calibration": {
                "positive_reference": (
                    "the 36 frozen SCAN-B cross-platform technical replicate pairs"
                ),
                "negative_reference": (
                    f"{TECH_NULL_NONMATCH_PAIRS:,} deterministic random SCAN-B cross-platform "
                    "NON-matching primary/replicate pairs"
                ),
                "candidate_threshold": (
                    f"max({100*TECH_TRUE_PERCENTILE:.0f}th percentile of true cross-platform replicate "
                    f"similarities, {100*TECH_NULL_PERCENTILE:.2f}th percentile of non-match similarities)"
                ),
                "separation_guard": (
                    "a hard automated threshold is accepted only if the true-replicate lower-tail "
                    "reference exceeds the non-match upper-tail reference; otherwise the audit is "
                    "rank-based and no sample is automatically called duplicate"
                ),
            },
            "flag_rule": (
                "flag a possible cross-cohort identity overlap only when the pair is mutual nearest "
                f"neighbors, exceeds the calibrated threshold when available, and each member's "
                f"similarity margin over its second-best cross-cohort match is >= "
                f"{IDENTITY_MIN_MARGIN_TO_SECOND_BEST:.2f}"
            ),
            "automatic_exclusion": False,
            "if_flagged": (
                "PAUSE before preservation; manually review identifiers/provenance. No flagged sample "
                "may be removed solely because removal improves preservation."
            ),
            "pairs_audited": [
                "TCGA-BRCA vs SCAN-B",
                "TCGA-BRCA vs METABRIC",
                "SCAN-B vs METABRIC",
            ],
        },

        "scanb_platform_sensitivity": {
            "role": "secondary technical heterogeneity sensitivity",
            "primary_pooled_target_remains": "all 3,273 frozen SCAN-B primary profiles",
            "platform_groups": {
                "GPL11154": 2969,
                "GPL18573": 304,
            },
            "frozen_source": "full 1,082-sample TCGA source object/program membership remains unchanged",
            "statistics": [
                "all-edge rho_edge point estimate",
                "rho_load point estimate",
            ],
            "uncertainty": {
                "target_sample_bootstraps": PLATFORM_BOOTSTRAPS,
                "edge_bootstrap_subset": "same deterministic up-to-20,000-edge rule used elsewhere",
                "ci": "95% percentile interval",
            },
            "formal_new_discovery_p_values": False,
            "changes_primary_classifier": False,
            "additional_reporting": (
                "report PAM50 composition within each SCAN-B platform to separate platform imbalance "
                "from subtype-composition imbalance"
            ),
        },

        "source_object_stability": {
            "role": (
                "characterize finite-sample stability of the already frozen TCGA source representation "
                "without rediscovering or reselecting modules"
            ),
            "module_membership": "fixed exactly to the 12 frozen source modules",
            "source_samples": "bootstrap the 1,082 frozen TCGA primary tumors with replacement",
            "edge_stability": {
                "replicates": SOURCE_EDGE_BOOTSTRAPS,
                "statistic": (
                    "Spearman concordance between the frozen full-source edge vector and each bootstrap "
                    "source edge vector on one deterministic fixed edge subset"
                ),
                "max_edges_per_module": SOURCE_STABILITY_MAX_EDGES,
                "summary": "median, 5th percentile, 95th percentile",
            },
            "loading_stability": {
                "replicates": SOURCE_LOADING_BOOTSTRAPS,
                "statistic": (
                    "Spearman concordance between frozen source PC1 loadings and bootstrap source PC1 "
                    "loadings after deterministic sign alignment to the frozen loading vector"
                ),
                "summary": "median, 5th percentile, 95th percentile",
            },
            "module_exclusion_based_on_stability": False,
            "posthoc_threshold": None,
            "interpretation": (
                "source stability is a descriptive property of the frozen estimand and cannot be used "
                "to remove inconvenient modules after target results are observed"
            ),
        },

        "bulk_composition_followup": {
            "status": "PRESPECIFIED_SECONDARY_SENSITIVITY_NOT_REQUIRED_BEFORE_FIRST_POOLED_RESULT",
            "reason": (
                "bulk breast-tumor coexpression may reflect immune/stromal admixture in addition to "
                "cell-intrinsic regulation"
            ),
            "planned_method": (
                "derive the same published ESTIMATE immune/stromal scores in each cohort; on identical "
                "samples, residualize each transformed gene on standardized immune and stromal scores, "
                "then recompute edge preservation"
            ),
            "comparison_rule": (
                "report unadjusted complete-case edge preservation beside composition-residualized "
                "preservation; no arbitrary retention threshold"
            ),
            "changes_primary_classifier": False,
            "can_be_tuned_after_primary_result": False,
        },

        "source_genomic_concentration_annotation": {
            "status": "PRESPECIFIED_DESCRIPTIVE_SOURCE_ONLY_AUDIT",
            "purpose": (
                "identify modules whose coherence may be strongly influenced by chromosomal/genomic "
                "concentration such as recurrent copy-number regions"
            ),
            "mapping": "hg38 gene annotation already downloaded with SCAN-B inputs",
            "report": [
                "fraction of module genes on the most represented chromosome",
                "chromosome-level gene counts",
            ],
            "module_exclusion": False,
            "changes_primary_classifier": False,
        },

        "computation": {
            "gpu_available": True,
            "gpu_model_from_audit": (
                gpu.get("devices", [{}])[0].get("name", "CUDA device")
                if gpu.get("devices")
                else "CUDA device"
            ),
            "gpu_preferred_for_cross_sample_similarity_and_bootstrap_correlations": True,
            "cpu_reference_validation_required": True,
            "base_seed": BASE_SEED,
        },
    }

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    json_out = OUT_DIR / "pretarget_confounds_contract_v1.json"
    md_out = OUT_DIR / "pretarget_confounds_contract_v1.md"

    json_out.write_text(json.dumps(contract, indent=2), encoding="utf-8")
    md_out.write_text(
        f"""# Paper 4 / TCBB — Final Pre-target Confounds Contract v1

Status: **FROZEN BEFORE FIRST TARGET PRESERVATION STATISTIC**

## Cross-cohort identity audit

Cross-cohort sample fingerprints will be computed on the pairwise intersection
of the frozen 10,000-gene source universe after cohort-wise gene
standardization.

The duplicate-screen threshold is calibrated from:
- 36 true SCAN-B cross-platform technical replicate pairs; and
- {TECH_NULL_NONMATCH_PAIRS:,} deterministic non-matching cross-platform pairs.

No sample is automatically excluded. Any credible overlap pauses the analysis
for provenance review.

## SCAN-B platform sensitivity

The pooled 3,273-sample SCAN-B result remains primary.

Secondary target strata:
- GPL11154: 2,969
- GPL18573: 304

Platform-specific effects and {PLATFORM_BOOTSTRAPS:,} target bootstraps are
reported without creating a new discovery family or changing classification.

## Frozen-source stability

Memberships are never rediscovered.

- edge stability: {SOURCE_EDGE_BOOTSTRAPS} TCGA sample bootstraps, fixed
  up-to-{SOURCE_STABILITY_MAX_EDGES:,}-edge subset;
- loading stability: {SOURCE_LOADING_BOOTSTRAPS} TCGA sample bootstraps.

No source-stability cutoff can exclude a module.

## Prespecified later sensitivities

Immune/stromal composition residualization and source chromosome-concentration
annotation are prespecified now but do not block the first pooled preservation
result. Neither can alter the primary classifier.
""",
        encoding="utf-8",
    )

    print()
    print("=" * 134)
    print("04f FINAL PRE-TARGET CONFOUNDS CONTRACT: PASS")
    print("=" * 134)
    print(f"Identity null nonmatches:        {TECH_NULL_NONMATCH_PAIRS:,}")
    print(f"Source edge bootstraps:          {SOURCE_EDGE_BOOTSTRAPS}")
    print(f"Source loading bootstraps:       {SOURCE_LOADING_BOOTSTRAPS}")
    print(f"SCAN-B platform bootstraps:      {PLATFORM_BOOTSTRAPS:,}")
    print("Bulk composition sensitivity:    prespecified; secondary")
    print("Genomic concentration audit:     prespecified; descriptive")
    print("Primary classifier changed:      NO")
    print("Preservation statistics seen:    NO")
    print()
    print("Outputs:")
    print(f"  {json_out}")
    print(f"  {md_out}")
    print("=" * 134)


if __name__ == "__main__":
    main()
