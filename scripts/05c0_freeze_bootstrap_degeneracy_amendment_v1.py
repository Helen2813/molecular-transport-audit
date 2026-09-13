from __future__ import annotations

import json
from pathlib import Path


SCRIPT_VERSION = "05c0-freeze-bootstrap-degeneracy-amendment-v1-no-cli"

DATA_ROOT = Path(r"D:\paper4_tcbb_data")
PRESERVATION_CONTRACT = (
    DATA_ROOT
    / "paper4_tcbb_preservation_operating_contract_v1"
    / "preservation_operating_characteristics_contract_v1.json"
)
FINAL_CLASSIFICATION = (
    DATA_ROOT
    / "paper4_tcbb_final_mapping_specificity_null_v1"
    / "final_mapping_specificity_and_classification_v1.json"
)
OUT_DIR = DATA_ROOT / "paper4_tcbb_bootstrap_degeneracy_amendment_v1"

VALID_BOOTSTRAPS_REQUIRED = 1000
MAX_BOOTSTRAP_ATTEMPTS = 10000


def require(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Required file does not exist: {path}")


def main() -> None:
    print("=" * 136)
    print("Paper 4 / TCBB - freeze handling of degenerate target-sample bootstrap draws")
    print("=" * 136)
    print(f"Script version: {SCRIPT_VERSION}")
    print()
    print("Trigger:")
    print("  05c v1 stopped on the FIRST SCAN-B M001 bootstrap draw because at least one")
    print("  originally evaluable target gene had zero/nonfinite SD in that resample.")
    print()
    print("Scientific guard:")
    print("  Bootstrap preservation effect recorded before this amendment: NO")
    print("  Bootstrap CI recorded before this amendment:                  NO")
    print("  Primary pooled classifications already frozen:                YES")
    print("  Primary classifications may change here:                      NO")
    print("=" * 136)

    for p in [PRESERVATION_CONTRACT, FINAL_CLASSIFICATION]:
        require(p)

    pcon = json.loads(PRESERVATION_CONTRACT.read_text(encoding="utf-8"))
    if pcon.get("status") != "FROZEN_BEFORE_FIRST_TARGET_PRESERVATION_STATISTIC":
        raise RuntimeError("04a preservation contract has unexpected status.")

    final = json.loads(FINAL_CLASSIFICATION.read_text(encoding="utf-8"))
    if final.get("status") != "FINAL_PRIMARY_POOLED_CLASSIFICATION_COMPLETE":
        raise RuntimeError("05b final pooled classification is not complete.")

    if int(pcon["target_sampling_uncertainty"]["bootstrap_replicates"]) != VALID_BOOTSTRAPS_REQUIRED:
        raise RuntimeError("Required bootstrap count differs from the frozen 04a contract.")

    contract = {
        "contract_id": "paper4-tcbb-bootstrap-degeneracy-amendment-v1",
        "script_version": SCRIPT_VERSION,
        "status": "FROZEN_BEFORE_ANY_BOOTSTRAP_PRESERVATION_RESULT",
        "trigger": (
            "05c v1 failed on its first SCAN-B TCGA_M001 bootstrap attempt because "
            "sample resampling with replacement can make an originally nonconstant gene "
            "constant or nonfinite within a particular bootstrap draw."
        ),
        "scientific_guard": {
            "bootstrap_preservation_effect_seen_before_amendment": False,
            "bootstrap_ci_seen_before_amendment": False,
            "primary_classification_already_frozen": True,
            "primary_classification_can_change": False,
        },
        "unchanged_bootstrap_estimand": {
            "resampled_unit": "target biological samples",
            "sampling": "nonparametric sampling with replacement at the original target sample size",
            "source_object": "fixed",
            "module_gene_set": "fixed to the statistically evaluable genes used in the observed target/program result",
            "edge_statistic": "full-evaluable-edge rho_edge",
            "loading_statistic": "target PC1 recomputed on the same valid resample, then rho_load",
        },
        "degenerate_draw_rule": {
            "valid_draw": (
                "after resampling target samples, every fixed module gene must have finite values "
                "and strictly positive sample SD in that bootstrap draw"
            ),
            "invalid_draw": (
                "if any fixed module gene has nonfinite or zero SD, the statistic is mathematically "
                "undefined on the frozen gene set; discard that bootstrap ATTEMPT and draw a new "
                "sample resample"
            ),
            "required_valid_draws": VALID_BOOTSTRAPS_REQUIRED,
            "maximum_total_attempts": MAX_BOOTSTRAP_ATTEMPTS,
            "if_cap_not_reached": (
                "bootstrap uncertainty is declared not estimable under the frozen nonparametric "
                "scheme for that target/program; do not drop genes, add jitter, add epsilon, "
                "or change the module after seeing the failure rate"
            ),
        },
        "explicitly_prohibited": [
            "dropping a degenerate gene only within an individual bootstrap replicate",
            "setting undefined correlations to zero",
            "adding random or deterministic jitter",
            "adding an epsilon to the resampled gene SD",
            "changing module membership or the observed target result",
            "changing the already frozen Strong/Partial/No-clear classification",
        ],
        "reporting": {
            "report_total_attempts": True,
            "report_invalid_attempts": True,
            "report_valid_attempt_fraction": True,
            "report_per_gene_degenerate_attempt_counts": True,
            "same_valid_resample_for_edge_and_loading": True,
        },
        "rationale": (
            "Redrawing an undefined bootstrap statistic preserves the fixed observed gene set "
            "and the prespecified nonparametric target-sample bootstrap. It avoids silently "
            "changing dimensionality across replicates or manufacturing correlations for "
            "zero-variance genes."
        ),
    }

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    json_out = OUT_DIR / "bootstrap_degeneracy_amendment_v1.json"
    md_out = OUT_DIR / "bootstrap_degeneracy_amendment_v1.md"

    json_out.write_text(json.dumps(contract, indent=2), encoding="utf-8")
    md_out.write_text(
        f"""# Paper 4 / TCBB — Bootstrap Degeneracy Amendment v1

Status: **FROZEN BEFORE ANY BOOTSTRAP PRESERVATION RESULT**

05c v1 stopped on its first SCAN-B M001 resample before an edge/loading
bootstrap statistic or CI was produced.

The target-sample bootstrap remains nonparametric sampling with replacement.
The frozen observed module gene set is not allowed to change.

A bootstrap attempt is valid only when every fixed module gene has finite
values and strictly positive sample SD in that resample. Undefined attempts
are discarded and redrawn until **{VALID_BOOTSTRAPS_REQUIRED} valid draws**
are obtained, with a maximum of **{MAX_BOOTSTRAP_ATTEMPTS} total attempts**.

If the cap is not sufficient, uncertainty is reported as not estimable under
this frozen scheme. No gene dropping, jitter, epsilon regularization, or
post-hoc module modification is permitted.

This amendment cannot change the already frozen primary classification.
""",
        encoding="utf-8",
    )

    print()
    print("=" * 136)
    print("05c0 BOOTSTRAP DEGENERACY AMENDMENT: PASS")
    print("=" * 136)
    print(f"Valid bootstrap draws required: {VALID_BOOTSTRAPS_REQUIRED:,}")
    print(f"Maximum total attempts:         {MAX_BOOTSTRAP_ATTEMPTS:,}")
    print("Per-replicate gene dropping:    PROHIBITED")
    print("Jitter / epsilon rescue:        PROHIBITED")
    print("Primary classification changed: NO")
    print("Bootstrap preservation seen:    NO")
    print()
    print(f"Output: {json_out}")
    print("=" * 136)


if __name__ == "__main__":
    main()
