from __future__ import annotations

import json
from pathlib import Path


SCRIPT_VERSION = "03f-freeze-mapping-null-pilot-contract-v1-no-cli"

DATA_ROOT = Path(r"D:\paper4_tcbb_data")
OUT_DIR = DATA_ROOT / "paper4_tcbb_mapping_null_pilot_contract_v1"

PILOT_PANELS = 100
MIN_VALID_PANELS_TO_PROCEED = 25
MIN_VALID_FRACTION_TO_PROCEED = MIN_VALID_PANELS_TO_PROCEED / PILOT_PANELS

BASE_SEED = 20260913

VALID_PANEL_MIN_FRACTION_WITHIN_DISTANCE_1 = 0.90
VALID_PANEL_MAX_STRATUM_DISTANCE = 2


def main() -> None:
    print("=" * 124)
    print("Paper 4 / TCBB - freeze mapping-null PILOT feasibility contract")
    print("=" * 124)
    print(f"Script version: {SCRIPT_VERSION}")
    print()
    print("Scientific guard:")
    print("  Target outcomes loaded:                              NO")
    print("  Target clinical characteristics loaded:              NO")
    print("  Target gene-gene correlations calculated:            NO")
    print("  Target PCA/loadings calculated:                       NO")
    print("  Target preservation statistics calculated:           NO")
    print("  Pilot operation: mapping quality only")
    print("=" * 124)

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    contract = {
        "contract_id": "paper4-tcbb-mapping-null-pilot-feasibility-v1",
        "script_version": SCRIPT_VERSION,
        "status": "FROZEN_BEFORE_PILOT_MAPPINGS",
        "pilot_panels_per_target_program": PILOT_PANELS,
        "scientific_guard": {
            "target_outcomes_loaded": False,
            "target_clinical_characteristics_loaded": False,
            "target_gene_gene_correlations_calculated": False,
            "target_pca_calculated": False,
            "target_preservation_statistics_calculated": False,
        },
        "panel_validity": {
            "minimum_fraction_assignments_with_stratum_distance_le_1":
                VALID_PANEL_MIN_FRACTION_WITHIN_DISTANCE_1,
            "maximum_assignment_stratum_distance":
                VALID_PANEL_MAX_STRATUM_DISTANCE,
        },
        "pilot_feasibility_gate": {
            "minimum_valid_panels": MIN_VALID_PANELS_TO_PROCEED,
            "minimum_valid_fraction": MIN_VALID_FRACTION_TO_PROCEED,
            "scope": (
                "gate applies to every primary-assessable target/program combination; "
                "non-assessable combinations are reported but do not block the design"
            ),
            "action_if_pass": (
                "retain the 03d mapping-null generator unchanged for final 1,000-panel specificity analysis"
            ),
            "action_if_fail": (
                "do not compute preservation; revise the mapping generator in a newly versioned "
                "source/marginal-only contract before any target preservation statistic is opened"
            ),
            "rationale": (
                "25% valid-pilot yield gives substantial margin relative to the already frozen "
                "10x-attempt cap for obtaining 1,000 final valid panels"
            ),
        },
        "seed": {
            "base_seed": BASE_SEED,
            "schedule": "deterministic target/program/panel offsets",
        },
    }

    json_out = OUT_DIR / "mapping_null_pilot_contract_v1.json"
    md_out = OUT_DIR / "mapping_null_pilot_contract_v1.md"

    json_out.write_text(json.dumps(contract, indent=2), encoding="utf-8")
    md_out.write_text(
        f"""# Paper 4 / TCBB — Mapping-Null Pilot Contract v1

Status: **FROZEN BEFORE PILOT MAPPINGS**

The pilot evaluates only whether the already frozen 03d target-gene mapping
generator can produce sufficiently close marginal matches.

- Pilot panels per target/program: **{PILOT_PANELS}**
- Valid panel: at least **{VALID_PANEL_MIN_FRACTION_WITHIN_DISTANCE_1:.0%}**
  of assignments have stratum distance <= 1, and no assignment has distance > {VALID_PANEL_MAX_STRATUM_DISTANCE}.
- Feasibility gate: at least **{MIN_VALID_PANELS_TO_PROCEED}/{PILOT_PANELS}**
  valid panels for every **primary-assessable** target/program combination.

If the gate fails, no preservation statistic is calculated. A newly versioned
mapping-generator contract must be frozen first.
""",
        encoding="utf-8",
    )

    print("=" * 124)
    print("03f MAPPING-NULL PILOT CONTRACT: PASS")
    print("=" * 124)
    print(f"Pilot panels per target/program: {PILOT_PANELS}")
    print(
        f"Proceed gate: >= {MIN_VALID_PANELS_TO_PROCEED}/{PILOT_PANELS} "
        f"valid panels for every primary-assessable combination"
    )
    print(
        f"Valid panel: >= {VALID_PANEL_MIN_FRACTION_WITHIN_DISTANCE_1:.0%} "
        f"assignments at distance <=1; max distance <= {VALID_PANEL_MAX_STRATUM_DISTANCE}"
    )
    print()
    print("No target correlation/PCA/preservation statistic was calculated.")
    print("Outputs:")
    print(f"  {json_out}")
    print(f"  {md_out}")
    print("=" * 124)


if __name__ == "__main__":
    main()
