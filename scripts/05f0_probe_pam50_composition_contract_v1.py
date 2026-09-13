from __future__ import annotations

import json
from pathlib import Path


SCRIPT_VERSION = "05f0-probe-pam50-composition-contract-v1-no-cli"

DATA_ROOT = Path(r"D:\paper4_tcbb_data")

KEYWORDS = (
    "pam50",
    "subtype",
    "composition",
    "residual",
    "luminal",
    "basal",
    "her2",
    "normal",
    "claudin",
    "platform",
    "technical",
    "repeatability",
    "estimate",
    "immune",
    "stromal",
)

# Restrict to JSON contracts/results already produced for Paper 4.
SEARCH_ROOTS = [
    DATA_ROOT,
]


def relevant_text(text: str) -> bool:
    t = text.lower()
    return any(k in t for k in KEYWORDS)


def print_relevant(obj, path: str = "root") -> None:
    if isinstance(obj, dict):
        for k, v in obj.items():
            here = f"{path}.{k}"
            hit = relevant_text(str(k))
            if isinstance(v, str):
                hit = hit or relevant_text(v)

            if hit:
                print(f"\n[{here}]")
                try:
                    print(json.dumps(v, indent=2, ensure_ascii=False))
                except Exception:
                    print(repr(v))
            elif isinstance(v, (dict, list)):
                print_relevant(v, here)

    elif isinstance(obj, list):
        for i, v in enumerate(obj):
            here = f"{path}[{i}]"
            if isinstance(v, str) and relevant_text(v):
                print(f"\n[{here}]")
                print(json.dumps(v, ensure_ascii=False))
            elif isinstance(v, (dict, list)):
                print_relevant(v, here)


def main() -> None:
    print("=" * 142)
    print("Paper 4 / TCBB - probe frozen PAM50/subtype/composition/technical-repeatability contract fields")
    print("=" * 142)
    print(f"Script version: {SCRIPT_VERSION}")
    print()
    print("Scientific guard:")
    print("  Target expression values read:                    NO")
    print("  PAM50/subtype result calculated:                  NO")
    print("  Composition-residualized result calculated:       NO")
    print("  Technical-repeatability result calculated:        NO")
    print("  Operation: inspect existing JSON contracts/results only")
    print("=" * 142)

    found_files = []
    for root in SEARCH_ROOTS:
        if not root.exists():
            raise FileNotFoundError(f"Data root does not exist: {root}")

        for path in root.rglob("*.json"):
            try:
                text = path.read_text(encoding="utf-8")
            except Exception:
                continue

            if relevant_text(path.name) or relevant_text(text):
                found_files.append(path)

    if not found_files:
        raise RuntimeError(
            "No JSON files containing PAM50/subtype/composition/platform/"
            "technical-repeatability/ESTIMATE keywords were found."
        )

    print(f"\nMatched JSON files: {len(found_files)}")

    for path in sorted(found_files):
        print("\n" + "#" * 142)
        print(f"FILE: {path}")
        print("#" * 142)

        try:
            obj = json.loads(path.read_text(encoding="utf-8"))
        except Exception as exc:
            print(f"Could not parse JSON: {exc}")
            continue

        if isinstance(obj, dict):
            if "status" in obj:
                print(f"Status: {obj['status']}")
            if "script_version" in obj:
                print(f"Script version: {obj['script_version']}")
            print(f"Top-level keys: {list(obj.keys())}")

        print_relevant(obj)

    print("\n" + "=" * 142)
    print("05f0 PAM50 / COMPOSITION CONTRACT PROBE: COMPLETE")
    print("=" * 142)
    print("No scientific statistic was calculated.")
    print("Paste this output back so the next runner can reuse the exact frozen tiers,")
    print("sample restrictions, residualization rules, seed schedule, and any predeclared platform controls.")
    print("=" * 142)


if __name__ == "__main__":
    main()
