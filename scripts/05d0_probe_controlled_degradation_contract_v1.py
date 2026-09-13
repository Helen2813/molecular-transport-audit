from __future__ import annotations

import json
from pathlib import Path


SCRIPT_VERSION = "05d0-probe-controlled-degradation-contract-v1-no-cli"

DATA_ROOT = Path(r"D:\paper4_tcbb_data")

FILES = {
    "04a preservation / operating-characteristics contract": (
        DATA_ROOT
        / "paper4_tcbb_preservation_operating_contract_v1"
        / "preservation_operating_characteristics_contract_v1.json"
    ),
    "03d matched-mapping generator contract": (
        DATA_ROOT
        / "paper4_tcbb_mapping_null_contract_v1"
        / "structure_preserving_mapping_null_contract_v1.json"
    ),
    "03h final mapping-null execution contract": (
        DATA_ROOT
        / "paper4_tcbb_final_mapping_null_execution_v1"
        / "final_mapping_null_execution_contract_v1.json"
    ),
    "04i3 exact-variation evaluability correction": (
        DATA_ROOT
        / "paper4_tcbb_exact_variation_evaluability_correction_v1"
        / "exact_variation_evaluability_correction_v1.json"
    ),
    "05b v2 corrected final classification": (
        DATA_ROOT
        / "paper4_tcbb_final_mapping_specificity_null_v2"
        / "final_mapping_specificity_and_classification_v2.json"
    ),
    "05c v4 corrected bootstrap / reliability": (
        DATA_ROOT
        / "paper4_tcbb_corrected_classical_bootstrap_reliability_v4"
        / "primary_pooled_corrected_classical_bootstrap_reliability_v4.json"
    ),
}

KEYWORDS = (
    "degrad",
    "mapping",
    "corrupt",
    "fraction",
    "operating",
    "sample_size",
    "sample-size",
    "monot",
    "calibr",
    "seed",
    "edge_subset",
    "edge subset",
    "panel",
    "specific",
)


def require(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Required file does not exist: {path}")


def key_relevant(key: str) -> bool:
    k = key.lower()
    return any(term in k for term in KEYWORDS)


def value_relevant(value) -> bool:
    if isinstance(value, str):
        v = value.lower()
        return any(term in v for term in KEYWORDS)
    return False


def print_relevant(obj, path: str = "root") -> None:
    if isinstance(obj, dict):
        for k, v in obj.items():
            here = f"{path}.{k}"
            if key_relevant(str(k)) or value_relevant(v):
                print(f"\n[{here}]")
                print(json.dumps(v, indent=2, ensure_ascii=False))
            elif isinstance(v, (dict, list)):
                print_relevant(v, here)

    elif isinstance(obj, list):
        for i, v in enumerate(obj):
            if isinstance(v, (dict, list)):
                print_relevant(v, f"{path}[{i}]")
            elif value_relevant(v):
                print(f"\n[{path}[{i}]]")
                print(json.dumps(v, indent=2, ensure_ascii=False))


def main() -> None:
    print("=" * 138)
    print("Paper 4 / TCBB - probe exact frozen contract fields for controlled mapping degradation")
    print("=" * 138)
    print(f"Script version: {SCRIPT_VERSION}")
    print()
    print("Scientific guard:")
    print("  Target expression values read:                    NO")
    print("  New preservation/degradation statistic calculated:NO")
    print("  Frozen results modified:                          NO")
    print("  Operation: JSON contract/schema inspection only")
    print("=" * 138)

    for label, path in FILES.items():
        require(path)

        print("\n" + "#" * 138)
        print(label)
        print(f"FILE: {path}")
        print("#" * 138)

        obj = json.loads(path.read_text(encoding="utf-8"))

        if isinstance(obj, dict):
            print(f"Top-level keys: {list(obj.keys())}")
            if "status" in obj:
                print(f"Status: {obj['status']}")
            if "script_version" in obj:
                print(f"Script version: {obj['script_version']}")

        print_relevant(obj)

    print("\n" + "=" * 138)
    print("05d0 CONTROLLED-DEGRADATION CONTRACT PROBE: COMPLETE")
    print("=" * 138)
    print("No target expression or new scientific statistic was calculated.")
    print("Paste this output back; it will be used to build the exact 05d runner without inventing any seed or gate.")
    print("=" * 138)


if __name__ == "__main__":
    main()
