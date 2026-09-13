from __future__ import annotations

import json
from pathlib import Path

import pandas as pd


SCRIPT_VERSION = "05g0-probe-scanb-technical-repeatability-inputs-v1-no-cli"

DATA_ROOT = Path(r"D:\paper4_tcbb_data")

OPERATING_CONTRACT = (
    DATA_ROOT
    / "paper4_tcbb_preservation_operating_contract_v1"
    / "preservation_operating_characteristics_contract_v1.json"
)
PAIRING_CONTRACT = (
    DATA_ROOT
    / "paper4_tcbb_scanb_pairing_contract_v1"
    / "scanb_pairing_contract_v1.json"
)
EXPRESSION_CONTRACT = (
    DATA_ROOT
    / "paper4_tcbb_expression_contract_v1"
    / "expression_input_contract_v1.json"
)
SAMPLE_SIZE_MASTER = (
    DATA_ROOT
    / "paper4_tcbb_scanb_sample_size_operating_characteristics_v1"
    / "scanb_sample_size_operating_characteristics_v1.json"
)
SAMPLE_SIZE_SUMMARY = (
    DATA_ROOT
    / "paper4_tcbb_scanb_sample_size_operating_characteristics_v1"
    / "scanb_sample_size_summary_v1.tsv"
)

PAIRING_DIR = DATA_ROOT / "paper4_tcbb_scanb_pairing_contract_v1"
INPUT_AUDIT_DIR = DATA_ROOT / "paper4_tcbb_input_audit_v1"

TABLE_SUFFIXES = {".tsv", ".csv", ".txt"}


def require(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Required path does not exist: {path}")


def inspect_table(path: Path, nrows: int = 8) -> None:
    try:
        if path.suffix.lower() in {".tsv", ".txt"}:
            df = pd.read_csv(path, sep="\t", dtype=str, nrows=nrows).fillna("")
        elif path.suffix.lower() == ".csv":
            df = pd.read_csv(path, dtype=str, nrows=nrows).fillna("")
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


def inspect_dir(directory: Path) -> None:
    require(directory)
    print("\n" + "#" * 144)
    print(f"DIRECTORY: {directory}")
    print("#" * 144)

    for path in sorted(p for p in directory.iterdir() if p.is_file()):
        print(f"\nFILE: {path.name}")
        print(f"  full path: {path}")

        if path.suffix.lower() == ".json":
            try:
                obj = json.loads(path.read_text(encoding="utf-8"))
                if isinstance(obj, dict):
                    print(f"  JSON top-level keys: {list(obj.keys())}")
                    if "status" in obj:
                        print(f"  status: {obj['status']}")
                    if "script_version" in obj:
                        print(f"  script_version: {obj['script_version']}")
            except Exception as exc:
                print(f"  JSON READ FAILED: {exc}")

        elif path.suffix.lower() in TABLE_SUFFIXES:
            inspect_table(path)


def main() -> None:
    print("=" * 144)
    print("Paper 4 / TCBB - probe exact SCAN-B technical-repeatability manifests and matched-n comparison inputs")
    print("=" * 144)
    print(f"Script version: {SCRIPT_VERSION}")
    print()
    print("Scientific guard:")
    print("  Expression matrices read:                         NO")
    print("  Technical-repeatability statistic calculated:     NO")
    print("  Platform-sensitivity statistic calculated:        NO")
    print("  Operation: inspect frozen contracts/manifests and 05e n=100 summary only")
    print("=" * 144)

    for p in [
        OPERATING_CONTRACT,
        PAIRING_CONTRACT,
        EXPRESSION_CONTRACT,
        SAMPLE_SIZE_MASTER,
        SAMPLE_SIZE_SUMMARY,
    ]:
        require(p)

    op = json.loads(OPERATING_CONTRACT.read_text(encoding="utf-8"))
    pairing = json.loads(PAIRING_CONTRACT.read_text(encoding="utf-8"))
    expr = json.loads(EXPRESSION_CONTRACT.read_text(encoding="utf-8"))
    ss = json.loads(SAMPLE_SIZE_MASTER.read_text(encoding="utf-8"))

    print("\n" + "#" * 144)
    print("04a TECHNICAL-REPEATABILITY CONTRACT BLOCK")
    print("#" * 144)
    print(json.dumps(op["scanb_technical_repeatability"], indent=2, ensure_ascii=False))
    print("\n04a computation:")
    print(json.dumps(op.get("computation", {}), indent=2, ensure_ascii=False))

    print("\n" + "#" * 144)
    print("01c FULL SCAN-B PAIRING CONTRACT")
    print("#" * 144)
    print(json.dumps(pairing, indent=2, ensure_ascii=False))

    print("\n" + "#" * 144)
    print("01b SCAN-B EXPRESSION/PROFILE CONTRACT FIELDS")
    print("#" * 144)
    print(json.dumps(expr.get("scanb_profiles", {}), indent=2, ensure_ascii=False))
    print(json.dumps(expr.get("expression_representation", {}).get("SCANB_GSE96058", {}), indent=2, ensure_ascii=False))

    print("\n" + "#" * 144)
    print("05e SAMPLE-SIZE MASTER")
    print("#" * 144)
    print(json.dumps(ss, indent=2, ensure_ascii=False))

    print("\n" + "#" * 144)
    print("05e n=100 MATCHED TRANSPORT SUMMARY")
    print("#" * 144)
    sdf = pd.read_csv(SAMPLE_SIZE_SUMMARY, sep="\t", dtype=str).fillna("")
    print(f"columns ({len(sdf.columns)}): {list(sdf.columns)}")
    n100 = sdf.loc[sdf["sample_size"].astype(str) == "100"].copy()
    print(n100.to_string(index=False))

    inspect_dir(PAIRING_DIR)

    print("\n" + "#" * 144)
    print(f"SELECTED INPUT-AUDIT FILES: {INPUT_AUDIT_DIR}")
    print("#" * 144)
    require(INPUT_AUDIT_DIR)

    keywords = ("scanb", "technical", "pair", "replicate", "inventory")
    for path in sorted(p for p in INPUT_AUDIT_DIR.iterdir() if p.is_file()):
        if not any(k in path.name.lower() for k in keywords):
            continue

        print(f"\nFILE: {path.name}")
        print(f"  full path: {path}")
        if path.suffix.lower() in TABLE_SUFFIXES:
            inspect_table(path)
        elif path.suffix.lower() == ".json":
            try:
                obj = json.loads(path.read_text(encoding="utf-8"))
                print(json.dumps(obj, indent=2, ensure_ascii=False))
            except Exception as exc:
                print(f"  JSON READ FAILED: {exc}")

    print("\n" + "=" * 144)
    print("05g0 TECHNICAL-REPEATABILITY INPUT PROBE: COMPLETE")
    print("=" * 144)
    print("No expression matrix or technical-repeatability statistic was calculated.")
    print("Paste this output back so the execution contract can freeze exact pair-manifest columns,")
    print("fixed-gene estimability handling, and matched-n reporting without guessing.")
    print("=" * 144)


if __name__ == "__main__":
    main()
