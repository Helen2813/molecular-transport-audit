from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import numpy as np


SCRIPT_VERSION = "05j0-probe-headtohead-corruption-artifacts-v1-no-cli"

DATA_ROOT = Path(r"D:\paper4_tcbb_data")
DEGRADATION_ROOT = (
    DATA_ROOT
    / "paper4_tcbb_controlled_mapping_degradation_v1"
)
INTERPRETATION_CONTRACT = (
    DATA_ROOT
    / "paper4_tcbb_comparator_interpretation_headtohead_contract_v1"
    / "comparator_interpretation_headtohead_contract_v1.json"
)


def print_table(path: Path) -> None:
    try:
        df = pd.read_csv(path, sep="\t", low_memory=False)
    except Exception as exc:
        print(f"  TSV READ FAILED: {exc}")
        return

    print(f"  rows={len(df):,}; columns ({len(df.columns)}): {list(df.columns)}")
    if len(df):
        print("  first rows:")
        print(df.head(8).to_string(index=False))


def print_npz(path: Path) -> None:
    try:
        z = np.load(path, allow_pickle=False)
    except Exception as exc:
        print(f"  NPZ READ FAILED: {exc}")
        return

    print(f"  arrays: {list(z.files)}")
    for key in z.files:
        arr = z[key]
        print(
            f"    {key}: shape={arr.shape}, dtype={arr.dtype}"
        )


def main() -> None:
    print("=" * 154)
    print("Paper 4 / TCBB - probe exact 05d corruption artifacts for mandatory comparator head-to-head")
    print("=" * 154)
    print(f"Script version: {SCRIPT_VERSION}")
    print()
    print("Scientific guard:")
    print("  New target expression statistic calculated:         NO")
    print("  NetRep/WGCNA comparator statistic calculated:       NO")
    print("  05d degradation statistic recalculated:             NO")
    print("  Operation: inspect serialized 05d artifacts/schema only")
    print("=" * 154)

    if not INTERPRETATION_CONTRACT.exists():
        raise FileNotFoundError(
            f"05i3 interpretation contract missing: {INTERPRETATION_CONTRACT}"
        )
    c = json.loads(INTERPRETATION_CONTRACT.read_text(encoding="utf-8"))
    if c.get("status") != "FROZEN_BEFORE_FIRST_NETREP_OR_WGCNA_COMPARATOR_RESULT":
        raise RuntimeError("05i3 interpretation contract has unexpected status.")

    if not DEGRADATION_ROOT.exists():
        raise FileNotFoundError(
            f"05d degradation root does not exist: {DEGRADATION_ROOT}"
        )

    files = sorted(
        p for p in DEGRADATION_ROOT.rglob("*")
        if p.is_file()
    )

    print(f"\n05d root: {DEGRADATION_ROOT}")
    print(f"Files found: {len(files)}")

    for path in files:
        rel = path.relative_to(DEGRADATION_ROOT)
        print("\n" + "#" * 154)
        print(f"FILE: {rel}")
        print("#" * 154)

        suffix = path.suffix.lower()
        if suffix == ".json":
            try:
                obj = json.loads(path.read_text(encoding="utf-8"))
                if isinstance(obj, dict):
                    print(f"  top-level keys: {list(obj.keys())}")
                    if "status" in obj:
                        print(f"  status: {obj['status']}")
                    if "script_version" in obj:
                        print(f"  script_version: {obj['script_version']}")
                text = json.dumps(obj, indent=2)
                print(text[:12000])
                if len(text) > 12000:
                    print("... <JSON truncated>")
            except Exception as exc:
                print(f"  JSON READ FAILED: {exc}")

        elif suffix in {".tsv", ".txt"}:
            if suffix == ".tsv":
                print_table(path)
            else:
                try:
                    text = path.read_text(encoding="utf-8", errors="ignore")
                    print(text[:12000])
                    if len(text) > 12000:
                        print("... <TEXT truncated>")
                except Exception as exc:
                    print(f"  TEXT READ FAILED: {exc}")

        elif suffix == ".npz":
            print_npz(path)

        else:
            print(f"  size={path.stat().st_size:,} bytes")

    print("\n" + "=" * 154)
    print("05j0 HEAD-TO-HEAD CORRUPTION ARTIFACT PROBE: COMPLETE")
    print("=" * 154)
    print("Paste this output back. The 05j execution runner will then either:")
    print("  (A) reuse exact serialized 05d corrupted mappings, or")
    print("  (B) replay the exact frozen 05d generator/seeds and verify rho_edge against saved 05d results")
    print("before calculating any comparator corruption statistic.")
    print("=" * 154)


if __name__ == "__main__":
    main()
