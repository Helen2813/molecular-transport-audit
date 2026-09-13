from __future__ import annotations

import json
from pathlib import Path

import pandas as pd


SCRIPT_VERSION = "05f1-probe-pam50-input-manifests-v1-no-cli"

DATA_ROOT = Path(r"D:\paper4_tcbb_data")

PAM50_CONTRACT = (
    DATA_ROOT
    / "paper4_tcbb_pam50_sensitivity_contract_v1"
    / "pam50_sensitivity_contract_v1.json"
)
ALIGNMENT_AUDIT = (
    DATA_ROOT
    / "paper4_tcbb_pam50_alignment_audit_v1"
    / "pam50_alignment_audit_v1.json"
)

INSPECT_DIRS = [
    DATA_ROOT / "paper4_tcbb_pam50_alignment_audit_v1",
    DATA_ROOT / "paper4_tcbb_breast_subtype_audit_v1",
    DATA_ROOT / "paper4_tcbb_final_pretarget_audits_v1",
]


def require(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Required path does not exist: {path}")


def inspect_table(path: Path) -> None:
    suffix = path.suffix.lower()

    try:
        if suffix in {".tsv", ".txt"}:
            df = pd.read_csv(path, sep="\t", dtype=str, nrows=5).fillna("")
        elif suffix == ".csv":
            df = pd.read_csv(path, dtype=str, nrows=5).fillna("")
        else:
            return
    except Exception as exc:
        print(f"  READ FAILED: {exc}")
        return

    print(f"  columns ({len(df.columns)}): {list(df.columns)}")
    print("  first rows:")
    if len(df):
        print(df.to_string(index=False))
    else:
        print("  <empty>")

    try:
        if suffix in {".tsv", ".txt"}:
            count = sum(1 for _ in path.open("r", encoding="utf-8", errors="ignore")) - 1
        elif suffix == ".csv":
            count = sum(1 for _ in path.open("r", encoding="utf-8", errors="ignore")) - 1
        else:
            count = None
        if count is not None:
            print(f"  approximate data rows: {max(count, 0)}")
    except Exception:
        pass


def main() -> None:
    print("=" * 144)
    print("Paper 4 / TCBB - probe exact PAM50 input manifests, schemas, and frozen computation fields")
    print("=" * 144)
    print(f"Script version: {SCRIPT_VERSION}")
    print()
    print("Scientific guard:")
    print("  Expression matrices read:                         NO")
    print("  Subtype preservation calculated:                  NO")
    print("  Composition residualization calculated:           NO")
    print("  Operation: inspect frozen JSON contracts + small metadata/manifests only")
    print("=" * 144)

    require(PAM50_CONTRACT)
    require(ALIGNMENT_AUDIT)

    contract = json.loads(PAM50_CONTRACT.read_text(encoding="utf-8"))
    alignment = json.loads(ALIGNMENT_AUDIT.read_text(encoding="utf-8"))

    print("\n" + "#" * 144)
    print("FULL 04e PAM50 SENSITIVITY CONTRACT")
    print(f"FILE: {PAM50_CONTRACT}")
    print("#" * 144)
    print(json.dumps(contract, indent=2, ensure_ascii=False))

    print("\n" + "#" * 144)
    print("FULL 04d PAM50 ALIGNMENT AUDIT")
    print(f"FILE: {ALIGNMENT_AUDIT}")
    print("#" * 144)
    print(json.dumps(alignment, indent=2, ensure_ascii=False))

    for directory in INSPECT_DIRS:
        require(directory)

        print("\n" + "#" * 144)
        print(f"DIRECTORY: {directory}")
        print("#" * 144)

        files = sorted(
            p for p in directory.iterdir()
            if p.is_file()
        )

        for path in files:
            print(f"\nFILE: {path.name}")
            print(f"  full path: {path}")

            if path.suffix.lower() == ".json":
                try:
                    obj = json.loads(path.read_text(encoding="utf-8"))
                    if isinstance(obj, dict):
                        print(f"  JSON top-level keys: {list(obj.keys())}")
                        if "status" in obj:
                            print(f"  status: {obj['status']}")
                except Exception as exc:
                    print(f"  JSON READ FAILED: {exc}")
            else:
                inspect_table(path)

    print("\n" + "=" * 144)
    print("05f1 PAM50 INPUT-MANIFEST PROBE: COMPLETE")
    print("=" * 144)
    print("No expression matrix or new scientific statistic was calculated.")
    print("Paste this output back; it is sufficient to build the exact subtype/composition runner")
    print("without guessing metadata filenames, sample-ID columns, or seed/computation rules.")
    print("=" * 144)


if __name__ == "__main__":
    main()
