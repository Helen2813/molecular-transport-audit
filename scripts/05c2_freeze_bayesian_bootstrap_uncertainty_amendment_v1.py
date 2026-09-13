from __future__ import annotations

import json
from pathlib import Path


SCRIPT_VERSION = "05c2-freeze-bayesian-bootstrap-uncertainty-amendment-v1-no-cli"

DATA_ROOT = Path(r"D:\paper4_tcbb_data")

PRESERVATION_CONTRACT = (
    DATA_ROOT
    / "paper4_tcbb_preservation_operating_contract_v1"
    / "preservation_operating_characteristics_contract_v1.json"
)
DEGENERACY_AMENDMENT = (
    DATA_ROOT
    / "paper4_tcbb_bootstrap_degeneracy_amendment_v1"
    / "bootstrap_degeneracy_amendment_v1.json"
)
V2_CHECKPOINT = (
    DATA_ROOT
    / "paper4_tcbb_primary_bootstrap_reliability_v2"
    / "per_program"
    / "SCANB_GSE96058"
    / "TCGA_M001_bootstrap_checkpoint_v2.npz"
)
FINAL_CLASSIFICATION = (
    DATA_ROOT
    / "paper4_tcbb_final_mapping_specificity_null_v1"
    / "final_mapping_specificity_and_classification_v1.json"
)

OUT_DIR = DATA_ROOT / "paper4_tcbb_bayesian_bootstrap_amendment_v1"

BAYESIAN_BOOTSTRAP_DRAWS = 1000
DIRICHLET_ALPHA = 1.0
INTERVAL = 0.95


def require(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Required file does not exist: {path}")


def main() -> None:
    import numpy as np

    print("=" * 138)
    print("Paper 4 / TCBB - freeze Bayesian-bootstrap replacement for pooled target uncertainty")
    print("=" * 138)
    print(f"Script version: {SCRIPT_VERSION}")
    print()
    print("Trigger:")
    print("  The frozen ordinary fixed-gene target bootstrap was empirically non-estimable:")
    print("  SCAN-B M001 produced 0 valid draws in the full 10,000-attempt cap.")
    print()
    print("Scientific guard:")
    print("  Ordinary-bootstrap preservation effect/CI obtained:  NO")
    print("  Primary pooled classifications already frozen:       YES")
    print("  Primary classifications may change here:             NO")
    print("  New uncertainty estimator chosen using effect size:   NO")
    print("=" * 138)

    for p in [
        PRESERVATION_CONTRACT,
        DEGENERACY_AMENDMENT,
        V2_CHECKPOINT,
        FINAL_CLASSIFICATION,
    ]:
        require(p)

    pcon = json.loads(PRESERVATION_CONTRACT.read_text(encoding="utf-8"))
    if pcon.get("status") != "FROZEN_BEFORE_FIRST_TARGET_PRESERVATION_STATISTIC":
        raise RuntimeError("04a preservation contract has unexpected status.")

    dcon = json.loads(DEGENERACY_AMENDMENT.read_text(encoding="utf-8"))
    if dcon.get("status") != "FROZEN_BEFORE_ANY_BOOTSTRAP_PRESERVATION_RESULT":
        raise RuntimeError("05c0 degeneracy amendment has unexpected status.")

    final = json.loads(FINAL_CLASSIFICATION.read_text(encoding="utf-8"))
    if final.get("status") != "FINAL_PRIMARY_POOLED_CLASSIFICATION_COMPLETE":
        raise RuntimeError("05b final classification is not complete.")

    ck = np.load(V2_CHECKPOINT)
    completed = int(ck["completed"])
    attempts = int(ck["attempts"])
    invalid = int(ck["invalid_attempts"])

    if (completed, attempts, invalid) != (0, 10_000, 10_000):
        raise RuntimeError(
            "The Bayesian-bootstrap amendment trigger does not match the "
            f"recorded 05c v2 state: {(completed, attempts, invalid)}"
        )

    contract = {
        "contract_id": "paper4-tcbb-bayesian-bootstrap-uncertainty-amendment-v1",
        "script_version": SCRIPT_VERSION,
        "status": "FROZEN_BEFORE_ANY_VALID_TARGET_BOOTSTRAP_PRESERVATION_RESULT",
        "supersedes_for_uncertainty_only": [
            "04a ordinary target-sample resampling-with-replacement bootstrap",
            "05c0 redraw rule for undefined ordinary bootstrap attempts",
        ],
        "trigger": {
            "target": "SCANB_GSE96058",
            "program_id": "TCGA_M001",
            "ordinary_bootstrap_valid_draws": completed,
            "ordinary_bootstrap_attempts": attempts,
            "ordinary_bootstrap_invalid_draws": invalid,
            "effect_or_ci_seen_before_amendment": False,
        },
        "unchanged_primary_analysis": {
            "05a_direct_effects": "unchanged",
            "05b_mapping_specificity": "unchanged",
            "primary_classification": "unchanged",
            "source_object": "fixed",
            "module_gene_sets": "fixed",
        },
        "bayesian_bootstrap": {
            "draws": BAYESIAN_BOOTSTRAP_DRAWS,
            "weight_distribution": (
                "Dirichlet(1,...,1) over all original target biological samples, "
                "implemented as iid Exp(1) weights normalized to sum 1"
            ),
            "why_no_zero_variance_artifact": (
                "all original samples receive strictly positive weight in every draw; "
                "an originally nonconstant finite gene therefore retains positive "
                "weighted variance except for numerical failure"
            ),
            "weighted_gene_mean": "sum_i w_i x_ig",
            "weighted_gene_variance": "sum_i w_i (x_ig - weighted_mean_g)^2",
            "weighted_gene_standardization": (
                "(x_ig-weighted_mean_g)/sqrt(weighted_variance_g)"
            ),
            "weighted_gene_correlation": (
                "Z' diag(w) Z after weighted standardization"
            ),
            "edge_statistic": (
                "Spearman correlation between frozen source full-edge vector and "
                "the weighted target full-edge correlation vector"
            ),
            "loading_statistic": (
                "leading eigenvector of the weighted target gene-correlation matrix, "
                "oriented by weighted correlation with the frozen-sign target score; "
                "then Spearman rho with frozen source loadings"
            ),
            "interval": f"{int(INTERVAL*100)}% percentile interval over Bayesian-bootstrap draws",
            "interpretation": (
                "nonparametric weighted-sample uncertainty interval; report explicitly "
                "as a Bayesian-bootstrap interval, not as the failed classical bootstrap CI"
            ),
        },
        "uniform_weight_replay_guard": {
            "required": True,
            "rule": (
                "weights w_i=1/n must reproduce the observed 05a rho_edge and rho_load "
                "to numerical tolerance before Bayesian-bootstrap draws are accepted"
            ),
        },
        "reliability": {
            "split_half_repeats": 2000,
            "spearman_brown": True,
            "leave_one_out": True,
            "unchanged_from_04a": True,
        },
        "reporting": {
            "classical_bootstrap_failure_reported": True,
            "bayesian_bootstrap_used_for_uncertainty": True,
            "classification_changed": False,
        },
    }

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    json_out = OUT_DIR / "bayesian_bootstrap_uncertainty_amendment_v1.json"
    md_out = OUT_DIR / "bayesian_bootstrap_uncertainty_amendment_v1.md"

    json_out.write_text(json.dumps(contract, indent=2), encoding="utf-8")
    md_out.write_text(
        f"""# Paper 4 / TCBB — Bayesian-Bootstrap Uncertainty Amendment v1

Status: **FROZEN BEFORE ANY VALID TARGET BOOTSTRAP PRESERVATION RESULT**

The prespecified ordinary sample-with-replacement bootstrap was not estimable
on the fixed high-dimensional gene set: SCAN-B M001 yielded **0 valid draws
in 10,000 attempts** because at least one fixed gene became constant in every
resample.

No bootstrap preservation effect or interval had been obtained.

For pooled target uncertainty only, the ordinary bootstrap is replaced by a
Bayesian bootstrap with **{BAYESIAN_BOOTSTRAP_DRAWS} Dirichlet(1,...,1) sample-weight
draws**. Every original target sample receives positive weight in every draw,
so the frozen gene set remains fixed and originally nonconstant genes retain
positive weighted variance.

Weighted means, variances, correlations and target PC1 are recomputed in each
draw. The reported uncertainty is the **{int(INTERVAL*100)}% Bayesian-bootstrap
percentile interval**.

A uniform-weight draw must exactly replay the observed 05a edge/loading
statistics before weighted draws are accepted.

Primary effects, specificity results, and Strong/Partial/No-clear
classifications are unchanged.
""",
        encoding="utf-8",
    )

    print()
    print("=" * 138)
    print("05c2 BAYESIAN-BOOTSTRAP UNCERTAINTY AMENDMENT: PASS")
    print("=" * 138)
    print(f"Bayesian-bootstrap draws:       {BAYESIAN_BOOTSTRAP_DRAWS:,}")
    print(f"Dirichlet concentration:       alpha={DIRICHLET_ALPHA:g} per sample")
    print(f"Uncertainty interval:           {int(INTERVAL*100)}% percentile")
    print("Uniform-weight replay required: YES")
    print("Primary classification changed: NO")
    print("Valid ordinary-bootstrap result seen before amendment: NO")
    print()
    print(f"Output: {json_out}")
    print("=" * 138)


if __name__ == "__main__":
    main()
