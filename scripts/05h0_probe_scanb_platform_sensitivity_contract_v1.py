from __future__ import annotations

import json
from pathlib import Path

import pandas as pd


SCRIPT_VERSION = "05h0-probe-scanb-platform-sensitivity-contract-v1-no-cli"

DATA_ROOT = Path(r"D:\paper4_tcbb_data")

PRETARGET_CONTRACT = (
    DATA_ROOT
    / "paper4_tcbb_pretarget_confounds_contract_v1"
    / "pretarget_confounds_contract_v1.json"
)
PAM50_ALIGNMENT = (
    DATA_ROOT
    / "paper4_tcbb_pam50_alignment_audit_v1"
    / "scanb_frozen_primary_pam50_alignment_v1.tsv"
)
PLATFORM_COMPOSITION = (
    DATA_ROOT
    / "paper4_tcbb_final_pretarget_audits_v1"
    / "scanb_platform_pam50_composition_v1.tsv"
)
TECHNICAL_MASTER = (
    DATA_ROOT
    / "paper4_tcbb_scanb_technical_repeatability_v2"
    / "scanb_technical_repeatability_v2.json"
)


def require(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Required path does not exist: {path}")


def main() -> None:
    print("=" * 146)
    print("Paper 4 / TCBB - probe exact frozen SCAN-B platform-sensitivity contract")
    print("=" * 146)
    print(f"Script version: {SCRIPT_VERSION}")
    print()
    print("Scientific guard:")
    print("  Expression matrices read:                       NO")
    print("  Platform-specific preservation calculated:      NO")
    print("  Platform bootstrap calculated:                  NO")
    print("  Operation: inspect frozen 04f contract + existing manifests only")
    print("=" * 146)

    for p in [
        PRETARGET_CONTRACT,
        PAM50_ALIGNMENT,
        PLATFORM_COMPOSITION,
        TECHNICAL_MASTER,
    ]:
        require(p)

    contract = json.loads(PRETARGET_CONTRACT.read_text(encoding="utf-8"))
    tech = json.loads(TECHNICAL_MASTER.read_text(encoding="utf-8"))

    print("\n" + "#" * 146)
    print("FULL 04f PRETARGET-CONFOUNDS CONTRACT")
    print(f"FILE: {PRETARGET_CONTRACT}")
    print("#" * 146)
    print(json.dumps(contract, indent=2, ensure_ascii=False))

    print("\n" + "#" * 146)
    print("SCAN-B PRIMARY PAM50 / PLATFORM MANIFEST")
    print(f"FILE: {PAM50_ALIGNMENT}")
    print("#" * 146)
    df = pd.read_csv(PAM50_ALIGNMENT, sep="\t", dtype=str).fillna("")
    print(f"columns ({len(df.columns)}): {list(df.columns)}")
    print(f"rows: {len(df)}")
    print("platform counts:")
    print(df["platform_id"].value_counts(dropna=False).to_string())
    print("\nfirst rows:")
    print(df.head(8).to_string(index=False))

    print("\n" + "#" * 146)
    print("PRETARGET PLATFORM × PAM50 COMPOSITION")
    print(f"FILE: {PLATFORM_COMPOSITION}")
    print("#" * 146)
    comp = pd.read_csv(PLATFORM_COMPOSITION, sep="\t", dtype=str).fillna("")
    print(comp.to_string(index=False))

    print("\n" + "#" * 146)
    print("05g TECHNICAL-REPEATABILITY MASTER")
    print(f"FILE: {TECHNICAL_MASTER}")
    print("#" * 146)
    print(json.dumps(tech, indent=2, ensure_ascii=False))

    print("\n" + "=" * 146)
    print("05h0 SCAN-B PLATFORM-SENSITIVITY CONTRACT PROBE: COMPLETE")
    print("=" * 146)
    print("No expression matrix or platform-specific scientific statistic was calculated.")
    print("Paste this output back so the exact bootstrap scope, seed, and execution details")
    print("can be frozen without guessing.")
    print("=" * 146)


if __name__ == "__main__":
    main()
