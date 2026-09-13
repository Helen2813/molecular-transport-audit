from __future__ import annotations

import json
from pathlib import Path

SCRIPT_VERSION = "04e-freeze-pam50-sensitivity-contract-v1-no-cli"

DATA_ROOT = Path(r"D:\paper4_tcbb_data")
ALIGNMENT_AUDIT = (
    DATA_ROOT
    / "paper4_tcbb_pam50_alignment_audit_v1"
    / "pam50_alignment_audit_v1.json"
)
OUT_DIR = DATA_ROOT / "paper4_tcbb_pam50_sensitivity_contract_v1"

CANONICAL_SUBTYPES = ["LumA", "LumB", "Basal", "Her2", "Normal"]

CORE_MIN_SOURCE_N = 100
CORE_MIN_TARGET_N = 100
LIMITED_MIN_SOURCE_N = 50
LIMITED_MIN_TARGET_N = 100

BOOTSTRAP_REPLICATES = 1000
BOOTSTRAP_CI = 0.95
MAX_BOOTSTRAP_EDGES = 20000
BASE_SEED = 20260917


def require(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Required file does not exist: {path}")


def classify_tier(source_n: int, target_n: int) -> str:
    if source_n >= CORE_MIN_SOURCE_N and target_n >= CORE_MIN_TARGET_N:
        return "CORE"
    if source_n >= LIMITED_MIN_SOURCE_N and target_n >= LIMITED_MIN_TARGET_N:
        return "LIMITED_N"
    return "NOT_INDIVIDUALLY_ASSESSABLE"


def main() -> None:
    print("=" * 132)
    print("Paper 4 / TCBB - freeze PAM50 subtype-stratified and subtype-residualized sensitivity contract")
    print("=" * 132)
    print(f"Script version: {SCRIPT_VERSION}")
    print()
    print("Scientific guard:")
    print("  Target outcomes loaded:                              NO")
    print("  Target treatment variables loaded:                   NO")
    print("  Target gene-gene correlations calculated:            NO")
    print("  Target PCA/loadings calculated:                      NO")
    print("  Preservation statistics calculated:                  NO")
    print("  Inputs used here: aligned PAM50/subtype COUNTS only")
    print("=" * 132)

    require(ALIGNMENT_AUDIT)
    audit = json.loads(ALIGNMENT_AUDIT.read_text(encoding="utf-8"))

    tcga = audit["TCGA_BRCA"]
    scanb = audit["SCANB_GSE96058"]
    met = audit["METABRIC"]

    source_counts = tcga["matched_counts"]
    scanb_counts = scanb["pam50_counts"]
    met_counts = met["pam50_compatible_counts"]

    expected_tcga = {"LumA": 417, "LumB": 188, "Basal": 139, "Her2": 67, "Normal": 23}
    expected_scanb = {"LumA": 1657, "LumB": 729, "Basal": 339, "Her2": 327, "Normal": 221}
    expected_met = {"LumA": 700, "LumB": 475, "Basal": 209, "Her2": 224, "Normal": 148}

    if source_counts != expected_tcga:
        raise RuntimeError(f"Unexpected TCGA aligned PAM50 counts: {source_counts}")
    if scanb_counts != expected_scanb:
        raise RuntimeError(f"Unexpected SCAN-B aligned PAM50 counts: {scanb_counts}")
    if met_counts != expected_met:
        raise RuntimeError(f"Unexpected METABRIC PAM50-compatible counts: {met_counts}")

    tiers = {}
    for target_name, target_counts in [
        ("SCANB_GSE96058", scanb_counts),
        ("METABRIC", met_counts),
    ]:
        tiers[target_name] = {}
        for subtype in CANONICAL_SUBTYPES:
            s_n = int(source_counts.get(subtype, 0))
            t_n = int(target_counts.get(subtype, 0))
            tiers[target_name][subtype] = {
                "source_n": s_n,
                "target_n": t_n,
                "tier": classify_tier(s_n, t_n),
            }

    core = {
        target: [s for s, x in by_subtype.items() if x["tier"] == "CORE"]
        for target, by_subtype in tiers.items()
    }
    limited = {
        target: [s for s, x in by_subtype.items() if x["tier"] == "LIMITED_N"]
        for target, by_subtype in tiers.items()
    }
    not_assessable = {
        target: [
            s for s, x in by_subtype.items()
            if x["tier"] == "NOT_INDIVIDUALLY_ASSESSABLE"
        ]
        for target, by_subtype in tiers.items()
    }

    contract = {
        "contract_id": "paper4-tcbb-pam50-sensitivity-v1",
        "script_version": SCRIPT_VERSION,
        "status": "FROZEN_BEFORE_FIRST_TARGET_PRESERVATION_STATISTIC",
        "alignment_audit": str(ALIGNMENT_AUDIT),
        "scientific_guard": {
            "target_outcomes_loaded": False,
            "target_treatment_loaded": False,
            "target_gene_gene_correlations_calculated": False,
            "target_pca_calculated": False,
            "preservation_statistics_calculated": False,
        },
        "aligned_counts": {
            "TCGA_BRCA_source": source_counts,
            "SCANB_GSE96058_target": scanb_counts,
            "METABRIC_target_pam50_compatible": met_counts,
            "TCGA_source_pam50_labeled_total": int(tcga["pam50_matched_profiles"]),
            "TCGA_source_full_total": int(tcga["frozen_expression_profiles"]),
            "SCANB_pam50_labeled_total": int(scanb["pam50_matched_profiles"]),
            "METABRIC_pam50_compatible_total": int(met["pam50_compatible_profiles"]),
        },
        "subtype_reporting_tiers": {
            "purpose": (
                "prevent over-interpretation of small within-subtype source cohorts; "
                "these are reporting tiers, not post-hoc significance thresholds"
            ),
            "core_rule": (
                f"source n >= {CORE_MIN_SOURCE_N} AND target n >= {CORE_MIN_TARGET_N}"
            ),
            "limited_n_rule": (
                f"{LIMITED_MIN_SOURCE_N} <= source n < {CORE_MIN_SOURCE_N} "
                f"AND target n >= {LIMITED_MIN_TARGET_N}"
            ),
            "not_individually_assessable_rule": "all remaining subtype/target pairs",
            "tiers": tiers,
            "core_subtypes": core,
            "limited_n_subtypes": limited,
            "not_individually_assessable_subtypes": not_assessable,
            "expected_result_from_counts_only": {
                "core": ["LumA", "LumB", "Basal"],
                "limited_n": ["Her2"],
                "not_individually_assessable": ["Normal"],
            },
        },
        "subtype_stratified_sensitivity": {
            "role": (
                "secondary heterogeneity analysis; it does not replace the pooled primary "
                "transport estimand and does not alter the primary classifier"
            ),
            "membership": "retain the 12 frozen TCGA-derived module memberships unchanged",
            "gene_mapping": "retain the frozen target gene mapping/coverage rules",
            "comparison": (
                "for each assessable PAM50 subtype, compare TCGA source samples of that subtype "
                "with target samples of the same subtype"
            ),
            "source_edges": (
                "recompute Pearson gene-gene correlations within the source subtype on fixed "
                "evaluable module genes"
            ),
            "target_edges": (
                "recompute Pearson gene-gene correlations within the matching target subtype "
                "on the same fixed evaluable genes"
            ),
            "gene_standardization": (
                "within each cohort/subtype, z-standardize each transformed gene across included "
                "subtype samples before Pearson correlations"
            ),
            "primary_effect": (
                "Spearman concordance between source-subtype and target-subtype edge vectors"
            ),
            "point_estimate_edges": "all nonredundant evaluable module edges",
            "uncertainty_edges": (
                f"one deterministic uniform edge subset of up to {MAX_BOOTSTRAP_EDGES:,} "
                "edges per target/program, reused across bootstrap replicates"
            ),
            "bootstrap": {
                "replicates": BOOTSTRAP_REPLICATES,
                "resampling": (
                    "independently resample source and target biological samples with replacement "
                    "within subtype; recompute correlations each replicate"
                ),
                "ci": f"{int(BOOTSTRAP_CI * 100)}% percentile interval",
                "formal_p_values": False,
                "bh_adjustment": False,
            },
            "loading_axis": {
                "status": "secondary descriptive point estimate only",
                "definition": (
                    "recompute source-subtype and target-subtype PC1 loadings on fixed evaluable "
                    "genes; orient each PC1 to its within-subtype frozen-sign score; report Spearman rho"
                ),
                "bootstrap": False,
                "classification_use": False,
            },
            "interpretation": {
                "CORE": "full subtype-sensitivity interpretation with bootstrap CI",
                "LIMITED_N": (
                    "report effect and bootstrap CI but explicitly label limited source sample size; "
                    "do not use for a strong subtype-wide robustness/failure claim"
                ),
                "NOT_INDIVIDUALLY_ASSESSABLE": (
                    "do not compute a standalone subtype transport conclusion"
                ),
            },
        },
        "composition_residualization_sensitivity": {
            "role": (
                "primary sensitivity for whether pooled transport is driven mainly by between-PAM50 "
                "subtype composition rather than within-subtype covariance"
            ),
            "complete_case_sets": {
                "TCGA_BRCA_source": int(tcga["pam50_matched_profiles"]),
                "SCANB_GSE96058": int(scanb["pam50_matched_profiles"]),
                "METABRIC": int(met["pam50_compatible_profiles"]),
            },
            "metabric_rule": (
                "use only canonical PAM50-compatible LumA/LumB/Basal/Her2/Normal samples; "
                "exclude Claudin-low and NC from this sensitivity"
            ),
            "required_two_step_comparison": {
                "step_B_complete_case_pooled": (
                    "recompute pooled source-target edge preservation on exactly the PAM50-complete "
                    "samples, with no subtype adjustment"
                ),
                "step_C_subtype_residualized": (
                    "on those identical samples, subtract each gene's within-cohort PAM50 subtype mean; "
                    "then z-standardize each gene residual vector and recompute edge preservation"
                ),
                "reason": (
                    "B versus C isolates subtype-composition removal; A(full primary) versus B "
                    "separately shows the effect of PAM50 label availability/sample restriction"
                ),
            },
            "residualization_model": (
                "for every gene and cohort: transformed_expression = intercept + categorical PAM50 "
                "subtype effect + residual; operationally subtract the within-subtype gene mean, "
                "then gene-wise z-standardize residuals"
            ),
            "module_membership": "frozen pooled TCGA memberships remain unchanged",
            "gene_mapping": "frozen evaluable-gene mapping remains unchanged",
            "primary_effect": "full-edge Spearman source-target edge concordance",
            "uncertainty": {
                "bootstrap_replicates": BOOTSTRAP_REPLICATES,
                "bootstrap_edges": f"fixed up-to-{MAX_BOOTSTRAP_EDGES:,}-edge subset",
                "resampling": (
                    "independently resample samples within PAM50 subtype in source and target; "
                    "refit subtype means and edge correlations each replicate"
                ),
                "ci": f"{int(BOOTSTRAP_CI * 100)}% percentile interval",
            },
            "reported_contrasts": [
                "A: full-cohort primary pooled rho_edge",
                "B: PAM50-complete-case pooled rho_edge",
                "C: PAM50-residualized rho_edge",
                "B minus A (sample-restriction effect)",
                "C minus B (subtype-composition effect)",
            ],
            "binary_classifier_change": False,
            "posthoc_threshold_for_robustness": None,
            "note": (
                "No arbitrary percent-retention cutoff is allowed. Report effects and CIs and "
                "interpret attenuation quantitatively."
            ),
        },
        "multiplicity_and_claim_guard": {
            "subtype_sensitivity_generates_new_discovery_q_values": False,
            "subtype_sensitivity_changes_primary_classification": False,
            "composition_sensitivity_changes_primary_classification": False,
            "claim_rule": (
                "PAM50 analyses may refine interpretation as broadly within-subtype, subtype-dependent, "
                "or composition-associated, but may not retroactively change the prespecified pooled "
                "primary statistical result"
            ),
        },
        "computation": {
            "gpu_permitted": True,
            "gpu_preferred_for_bootstrap_edge_recomputation": True,
            "cpu_reference_validation_required": True,
            "base_seed": BASE_SEED,
        },
    }

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    json_out = OUT_DIR / "pam50_sensitivity_contract_v1.json"
    md_out = OUT_DIR / "pam50_sensitivity_contract_v1.md"

    json_out.write_text(json.dumps(contract, indent=2), encoding="utf-8")

    md_out.write_text(
        f"""# Paper 4 / TCBB — PAM50 Sensitivity Contract v1

Status: **FROZEN BEFORE FIRST TARGET PRESERVATION STATISTIC**

## Aligned PAM50 counts

- TCGA frozen source: {source_counts}
- SCAN-B frozen primary target: {scanb_counts}
- METABRIC PAM50-compatible target: {met_counts}

## Prespecified reporting tiers

- **CORE:** source n >= {CORE_MIN_SOURCE_N} and target n >= {CORE_MIN_TARGET_N}
  - LumA, LumB, Basal
- **LIMITED_N:** {LIMITED_MIN_SOURCE_N} <= source n < {CORE_MIN_SOURCE_N}
  and target n >= {LIMITED_MIN_TARGET_N}
  - Her2
- **NOT_INDIVIDUALLY_ASSESSABLE:** remaining subtype pairs
  - Normal

These are sensitivity-reporting guards, not significance thresholds.

## Subtype-stratified transport

Frozen module memberships and target mappings remain unchanged. Source and target
correlation structures are recomputed within the same PAM50 subtype. Point
estimates use all evaluable edges. Uncertainty uses {BOOTSTRAP_REPLICATES}
independent source/target sample bootstraps on one fixed up-to-{MAX_BOOTSTRAP_EDGES:,}-edge
subset.

No subtype-specific discovery p-values are generated.

## PAM50 composition residualization

Three quantities are required:

1. **A:** original full-cohort pooled preservation.
2. **B:** pooled preservation restricted to PAM50-complete samples.
3. **C:** preservation after subtracting cohort-specific PAM50 subtype means
   from every gene on exactly the same samples as B.

Therefore B-A measures sample-restriction effects, whereas C-B isolates
attenuation attributable to between-subtype composition.

## METABRIC

Only LumA/LumB/Basal/Her2/Normal enter the PAM50-composition analysis.
Claudin-low and NC remain outside this sensitivity.

## Primary-analysis guard

PAM50 sensitivity results do not alter the prespecified pooled
Strong/Partial/No-clear classification. They refine interpretation.
""",
        encoding="utf-8",
    )

    print()
    print("Prespecified subtype reporting tiers:")
    for target, by_subtype in tiers.items():
        print(f"  {target}:")
        for subtype in CANONICAL_SUBTYPES:
            x = by_subtype[subtype]
            print(
                f"    {subtype:6s} source_n={x['source_n']:4d} "
                f"target_n={x['target_n']:4d} -> {x['tier']}"
            )

    print()
    print("=" * 132)
    print("04e PAM50 SENSITIVITY CONTRACT: PASS")
    print("=" * 132)
    print("Core subtype sensitivity:        LumA, LumB, Basal")
    print("Limited-n subtype sensitivity:   Her2")
    print("Not individually assessable:     Normal")
    print(f"Subtype/composition bootstraps:  {BOOTSTRAP_REPLICATES:,}")
    print("Residualization comparator:      A(full) -> B(complete-case) -> C(PAM50-residualized)")
    print("Primary classifier changed:      NO")
    print("Preservation statistics seen:    NO")
    print()
    print("Outputs:")
    print(f"  {json_out}")
    print(f"  {md_out}")
    print("=" * 132)


if __name__ == "__main__":
    main()
