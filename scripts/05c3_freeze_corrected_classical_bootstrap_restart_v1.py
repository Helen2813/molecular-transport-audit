from __future__ import annotations

import json
from pathlib import Path


SCRIPT_VERSION = "05c3-freeze-corrected-classical-bootstrap-restart-v1-no-cli"

DATA_ROOT = Path(r"D:\paper4_tcbb_data")

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
DEGENERACY_RULE = (
    DATA_ROOT
    / "paper4_tcbb_bootstrap_degeneracy_amendment_v1"
    / "bootstrap_degeneracy_amendment_v1.json"
)
BAYESIAN_AMENDMENT = (
    DATA_ROOT
    / "paper4_tcbb_bayesian_bootstrap_amendment_v1"
    / "bayesian_bootstrap_uncertainty_amendment_v1.json"
)

OUT_DIR = DATA_ROOT / "paper4_tcbb_corrected_classical_bootstrap_restart_v1"

VALID_DRAWS = 1000
MAX_ATTEMPTS = 10000
CI = 0.95
SPLIT_HALF_REPEATS = 2000


def require(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Required file does not exist: {path}")


def main() -> None:
    print("=" * 140)
    print("Paper 4 / TCBB - freeze restart of ORIGINAL classical target bootstrap after exact-variation correction")
    print("=" * 140)
    print(f"Script version: {SCRIPT_VERSION}")
    print()
    print("Rationale:")
    print("  05c1 showed the universal failure was caused by five exact-constant SCAN-B genes")
    print("  that had been admitted by the old floating-point SD>0 screen.")
    print("  04i3 removed exact-constant genes universally using finite AND exact max>min.")
    print("  Corrected 05a v2 / 05b v2 have now been rerun.")
    print()
    print("Frozen restart:")
    print("  Return to the originally planned nonparametric target-sample bootstrap.")
    print("  The Bayesian-bootstrap amendment remains archival and is NOT used.")
    print("  Rare undefined resamples, if any, follow the already frozen 05c0 redraw rule.")
    print("=" * 140)

    for p in [
        EXACT_CORRECTION,
        FINAL_CLASSIFICATION,
        DEGENERACY_RULE,
        BAYESIAN_AMENDMENT,
    ]:
        require(p)

    exact = json.loads(EXACT_CORRECTION.read_text(encoding="utf-8"))
    if exact.get("status") != "FROZEN_CORRECTION_BEFORE_RECOMPUTED_PRESERVATION":
        raise RuntimeError("04i3 exact-variation correction has unexpected status.")

    counts = exact.get("counts", {})
    if int(counts.get("SCANB_corrected_evaluable", -1)) != 9220:
        raise RuntimeError("04i3 SCAN-B corrected count is not 9,220.")
    if int(counts.get("METABRIC_corrected_evaluable", -1)) != 8485:
        raise RuntimeError("04i3 METABRIC corrected count is not 8,485.")

    final = json.loads(FINAL_CLASSIFICATION.read_text(encoding="utf-8"))
    if final.get("status") != "CORRECTED_FINAL_PRIMARY_POOLED_CLASSIFICATION_COMPLETE":
        raise RuntimeError("Corrected 05b v2 final classification is not complete.")

    dcon = json.loads(DEGENERACY_RULE.read_text(encoding="utf-8"))
    if dcon.get("status") != "FROZEN_BEFORE_ANY_BOOTSTRAP_PRESERVATION_RESULT":
        raise RuntimeError("05c0 degeneracy/redraw rule has unexpected status.")

    bayes = json.loads(BAYESIAN_AMENDMENT.read_text(encoding="utf-8"))
    if bayes.get("status") != "FROZEN_BEFORE_ANY_VALID_TARGET_BOOTSTRAP_PRESERVATION_RESULT":
        raise RuntimeError("05c2 Bayesian amendment has unexpected status.")

    contract = {
        "contract_id": "paper4-tcbb-corrected-classical-bootstrap-restart-v1",
        "script_version": SCRIPT_VERSION,
        "status": "FROZEN_BEFORE_CORRECTED_CLASSICAL_BOOTSTRAP_RESULTS",
        "reason": (
            "The original classical bootstrap failed because exact-constant target genes "
            "were admitted by the prior floating SD>0 screen. The universal 04i3 "
            "finite-and-exact-max>min correction removed those genes before corrected "
            "05a/05b were rerun."
        ),
        "uncertainty_estimator": {
            "type": "ordinary nonparametric target-sample bootstrap",
            "sampling": "sample target biological samples with replacement at original n",
            "source_object": "fixed",
            "corrected_target_gene_set": "fixed within target/program",
            "valid_draws_required": VALID_DRAWS,
            "maximum_total_attempts": MAX_ATTEMPTS,
            "undefined_draw_rule": (
                "if any fixed corrected module gene has zero/nonfinite variance in a "
                "particular resample, discard that attempt and redraw; do not alter genes"
            ),
            "edge_statistic": "full-evaluable-edge rho_edge",
            "loading_statistic": "recompute target PC1, then rho_load",
            "interval": f"{int(CI*100)}% percentile interval",
        },
        "implementation_equivalence": {
            "allowed": (
                "represent a sample-with-replacement resample by integer multiplicity "
                "counts and compute the exactly equivalent count-weighted means, "
                "variances, correlations, and orientation"
            ),
            "required_validation": (
                "uniform count=1 replay must reproduce corrected 05a v2 observed "
                "rho_edge/rho_load; first bootstrap count-weight engine must also pass "
                "an explicit-resample correlation reference check on a fixed gene subset"
            ),
        },
        "reliability": {
            "split_half_repeats": SPLIT_HALF_REPEATS,
            "spearman_brown_correction": True,
            "leave_one_gene_out": True,
            "classification_threshold": None,
        },
        "bayesian_bootstrap": {
            "status": "ARCHIVAL_NOT_USED",
            "reason": (
                "it was proposed before the root cause was identified as exact-constant "
                "genes admitted by the old evaluability screen"
            ),
        },
        "primary_classification": {
            "source": "corrected 05b v2",
            "may_change_here": False,
        },
        "target_outcomes_or_treatment_loaded": False,
    }

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out = OUT_DIR / "corrected_classical_bootstrap_restart_v1.json"
    out.write_text(json.dumps(contract, indent=2), encoding="utf-8")

    print()
    print("=" * 140)
    print("05c3 CORRECTED CLASSICAL-BOOTSTRAP RESTART: PASS")
    print("=" * 140)
    print(f"Valid draws required:            {VALID_DRAWS:,}")
    print(f"Maximum attempts:                {MAX_ATTEMPTS:,}")
    print(f"Percentile interval:             {int(CI*100)}%")
    print("Exact-constant genes excluded:   YES (04i3)")
    print("Bayesian bootstrap used:         NO")
    print("Rare undefined draw -> redraw:   YES (05c0)")
    print("Primary classification changed: NO")
    print()
    print(f"Output: {out}")
    print("=" * 140)


if __name__ == "__main__":
    main()
