from __future__ import annotations

import json
from pathlib import Path


SCRIPT_VERSION = "05i0-probe-comparator-benchmark-contract-v1-no-cli"

DATA_ROOT = Path(r"D:\paper4_tcbb_data")
PROJECT_ROOT = Path(r"C:\Users\olegk\Desktop\molecular-transport-audit")

KEYWORDS = (
    "netrep",
    "modulepreservation",
    "module preservation",
    "comparator",
    "benchmark",
    "wgcna",
    "preservation benchmark",
    "zsummary",
    "medianrank",
)

SEARCH_ROOTS = [
    DATA_ROOT,
    PROJECT_ROOT / "scripts",
]


def relevant_text(text: str) -> bool:
    t = text.lower()
    return any(k in t for k in KEYWORDS)


def print_json_hits(obj, path: str = "root") -> None:
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
                print_json_hits(v, here)

    elif isinstance(obj, list):
        for i, v in enumerate(obj):
            here = f"{path}[{i}]"
            if isinstance(v, str) and relevant_text(v):
                print(f"\n[{here}]")
                print(json.dumps(v, ensure_ascii=False))
            elif isinstance(v, (dict, list)):
                print_json_hits(v, here)


def main() -> None:
    print("=" * 148)
    print("Paper 4 / TCBB - probe frozen comparator-benchmarking contract and existing implementation files")
    print("=" * 148)
    print(f"Script version: {SCRIPT_VERSION}")
    print()
    print("Scientific guard:")
    print("  Target expression values read:                    NO")
    print("  NetRep statistics calculated:                     NO")
    print("  WGCNA modulePreservation calculated:              NO")
    print("  Operation: inspect existing JSON/text/R/Python files only")
    print("=" * 148)

    matched_json = []
    matched_text = []

    for root in SEARCH_ROOTS:
        if not root.exists():
            print(f"\nWARNING: search root does not exist: {root}")
            continue

        for path in root.rglob("*"):
            if not path.is_file():
                continue

            suffix = path.suffix.lower()
            if suffix not in {".json", ".md", ".txt", ".py", ".r"}:
                continue

            try:
                text = path.read_text(encoding="utf-8", errors="ignore")
            except Exception:
                continue

            if not (relevant_text(path.name) or relevant_text(text)):
                continue

            if suffix == ".json":
                matched_json.append(path)
            else:
                matched_text.append(path)

    print(f"\nMatched JSON files: {len(matched_json)}")
    print(f"Matched text/code files: {len(matched_text)}")

    for path in sorted(matched_json):
        print("\n" + "#" * 148)
        print(f"JSON FILE: {path}")
        print("#" * 148)

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

        print_json_hits(obj)

    for path in sorted(matched_text):
        print("\n" + "#" * 148)
        print(f"TEXT/CODE FILE: {path}")
        print("#" * 148)

        try:
            lines = path.read_text(
                encoding="utf-8",
                errors="ignore",
            ).splitlines()
        except Exception as exc:
            print(f"READ FAILED: {exc}")
            continue

        hit_lines = []
        for i, line in enumerate(lines, start=1):
            if relevant_text(line):
                lo = max(1, i - 3)
                hi = min(len(lines), i + 8)
                hit_lines.append((i, lo, hi))

        merged = []
        for _, lo, hi in hit_lines:
            if not merged or lo > merged[-1][1] + 1:
                merged.append([lo, hi])
            else:
                merged[-1][1] = max(merged[-1][1], hi)

        for lo, hi in merged[:30]:
            print(f"\n--- lines {lo}-{hi} ---")
            for j in range(lo, hi + 1):
                print(f"{j:05d}: {lines[j-1]}")

    print("\n" + "=" * 148)
    print("05i0 COMPARATOR-BENCHMARK CONTRACT PROBE: COMPLETE")
    print("=" * 148)
    print("No target matrix or comparator statistic was calculated.")
    print("Paste this output back so the NetRep/WGCNA benchmark can reuse exact frozen")
    print("module sets, target mappings, sample restrictions, permutation counts, seeds,")
    print("and reporting rules without introducing post-result choices.")
    print("=" * 148)


if __name__ == "__main__":
    main()
