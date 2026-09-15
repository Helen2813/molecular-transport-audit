from __future__ import annotations

import ast
import json
import re
from pathlib import Path


SCRIPT_VERSION = "05j1-probe-controlled-degradation-generator-v1-no-cli"

PROJECT_ROOT = Path(r"C:\Users\olegk\Desktop\molecular-transport-audit")
DATA_ROOT = Path(r"D:\paper4_tcbb_data")

RUNNER = PROJECT_ROOT / "scripts" / "05d_run_controlled_mapping_degradation_gpu_v1.py"
FREEZE = PROJECT_ROOT / "scripts" / "05d1_freeze_controlled_degradation_execution_v1.py"

DEGRADATION_ROOT = DATA_ROOT / "paper4_tcbb_controlled_mapping_degradation_v1"
REPLICATES = DEGRADATION_ROOT / "controlled_degradation_replicates_v1.tsv"
MASTER = DEGRADATION_ROOT / "controlled_mapping_degradation_v1.json"

INTERPRETATION_CONTRACT = (
    DATA_ROOT
    / "paper4_tcbb_comparator_interpretation_headtohead_contract_v1"
    / "comparator_interpretation_headtohead_contract_v1.json"
)

KEYWORDS = (
    "seed",
    "rng",
    "random",
    "corrupt",
    "replacement",
    "candidate",
    "mapping",
    "distance",
    "strata",
    "panel",
    "attempt",
    "edge",
    "specificity",
    "source_to_target",
    "slot",
)


def require(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Required file does not exist: {path}")


def line_interesting(line: str) -> bool:
    x = line.lower()
    return any(k in x for k in KEYWORDS)


def print_relevant_source(path: Path, max_context: int = 8) -> None:
    text = path.read_text(encoding="utf-8", errors="ignore")
    lines = text.splitlines()

    print("\n" + "#" * 156)
    print(f"FILE: {path}")
    print("#" * 156)

    print(f"Total lines: {len(lines)}")

    # Parse function names so we can identify generator helpers explicitly.
    try:
        tree = ast.parse(text)
        funcs = [
            node.name
            for node in tree.body
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
        ]
        print(f"Top-level functions ({len(funcs)}): {funcs}")
    except Exception as exc:
        print(f"AST parse warning: {exc}")

    hits = []
    for i, line in enumerate(lines, start=1):
        if line_interesting(line):
            lo = max(1, i - max_context)
            hi = min(len(lines), i + max_context)
            hits.append((lo, hi))

    # Merge overlapping ranges.
    merged: list[list[int]] = []
    for lo, hi in hits:
        if not merged or lo > merged[-1][1] + 1:
            merged.append([lo, hi])
        else:
            merged[-1][1] = max(merged[-1][1], hi)

    print(f"Relevant merged source ranges: {len(merged)}")

    for lo, hi in merged:
        print(f"\n--- lines {lo}-{hi} ---")
        for j in range(lo, hi + 1):
            print(f"{j:05d}: {lines[j - 1]}")

    # Also print top-level constants/assignments likely needed for exact replay.
    print("\n--- top-level assignments containing generator-related names ---")
    try:
        tree = ast.parse(text)
        for node in tree.body:
            if isinstance(node, (ast.Assign, ast.AnnAssign)):
                seg = ast.get_source_segment(text, node)
                if seg and line_interesting(seg):
                    print(seg)
    except Exception:
        pass


def main() -> None:
    print("=" * 156)
    print("Paper 4 / TCBB - probe exact 05d controlled-degradation generator for comparator head-to-head replay")
    print("=" * 156)
    print(f"Script version: {SCRIPT_VERSION}")
    print()
    print("Scientific guard:")
    print("  New target expression statistic calculated:         NO")
    print("  New corruption mapping generated:                    NO")
    print("  NetRep/WGCNA comparator statistic calculated:       NO")
    print("  Operation: inspect exact 05d source code/contracts only")
    print("=" * 156)

    for p in [
        RUNNER,
        FREEZE,
        REPLICATES,
        MASTER,
        INTERPRETATION_CONTRACT,
    ]:
        require(p)

    master = json.loads(MASTER.read_text(encoding="utf-8"))
    interp = json.loads(INTERPRETATION_CONTRACT.read_text(encoding="utf-8"))

    if master.get("status") != "CONTROLLED_DEGRADATION_ALL_FROZEN_GATES_PASS":
        raise RuntimeError("05d degradation master has unexpected status.")
    if interp.get("status") != "FROZEN_BEFORE_FIRST_NETREP_OR_WGCNA_COMPARATOR_RESULT":
        raise RuntimeError("05i3 interpretation contract has unexpected status.")

    print("\n05d master:")
    print(json.dumps(master, indent=2))

    print_relevant_source(FREEZE)
    print_relevant_source(RUNNER)

    print("\n" + "=" * 156)
    print("05j1 CONTROLLED-DEGRADATION GENERATOR PROBE: COMPLETE")
    print("=" * 156)
    print("No scientific result was calculated.")
    print("Paste this output back. The next runner will copy/reuse the exact 05d")
    print("generator logic, replay the saved (target, program, fraction, replicate_id, attempt_id)")
    print("rows, and REQUIRE rho_edge reproduction before any comparator corruption statistic")
    print("is accepted.")
    print("=" * 156)


if __name__ == "__main__":
    main()
