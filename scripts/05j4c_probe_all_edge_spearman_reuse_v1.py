#!/usr/bin/env python
from __future__ import annotations

import ast
import hashlib
import json
import re
from pathlib import Path

SCRIPT_VERSION = "05j4c-probe-all-edge-spearman-reuse-v1-no-cli"

ROOT = Path(__file__).resolve().parents[1]
DATA_ROOT = Path(r"D:\paper4_tcbb_data")

CONTRACT = (
    DATA_ROOT
    / "paper4_tcbb_scale_invariant_diagnostic_completion_contract_v1"
    / "scale_invariant_diagnostic_completion_contract_v1.json"
)

CANDIDATE_SOURCES = [
    ROOT / "scripts" / "05j4b_prepare_scale_invariant_inputs_gpu_v2.py",
    ROOT / "scripts" / "05j4_run_scale_invariant_corruption_audit_v3.py",
    ROOT / "scripts" / "05j_run_structural_corruption_headtohead_gpu_v1.py",
    ROOT / "scripts" / "05j2_extract_wgcna_observed_structural_v1.R",
]

OUT_DIR = DATA_ROOT / "paper4_tcbb_scale_invariant_diagnostic_completion_contract_v1"
OUT_JSON = OUT_DIR / "all_edge_spearman_implementation_probe_v1.json"
OUT_TXT = OUT_DIR / "all_edge_spearman_implementation_probe_v1.txt"

SEP = "=" * 168

TOKENS = [
    "netrep", "cor.cor", "pearson", "corrcoef", "correlation",
    "edge", "upper", "triu", "mapping", "mapped", "permutation",
    "source_corr", "target_corr", "source_cor", "target_cor",
    "pearson_subset_matched", "fullnull", "full_null",
    "threshold", "preserv", "retention", "rankdata", "spearman",
]

def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)

def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()

def python_function_ranges(text: str):
    try:
        tree = ast.parse(text)
    except SyntaxError:
        return []
    out = []
    for node in ast.walk(tree):
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            out.append({
                "name": node.name,
                "lineno": node.lineno,
                "end_lineno": getattr(node, "end_lineno", node.lineno),
            })
    return sorted(out, key=lambda x: x["lineno"])

def merged_relevant_ranges(lines: list[str], radius: int = 12):
    hits = []
    low_tokens = [t.lower() for t in TOKENS]
    for i, line in enumerate(lines, start=1):
        ll = line.lower()
        if any(t in ll for t in low_tokens):
            hits.append((max(1, i - radius), min(len(lines), i + radius)))
    if not hits:
        return []
    hits.sort()
    merged = [list(hits[0])]
    for a, b in hits[1:]:
        if a <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], b)
        else:
            merged.append([a, b])
    return [tuple(x) for x in merged]

def main() -> None:
    print(SEP)
    print("Paper 4 / TCBB - source-only probe for the frozen all-edge Pearson code path")
    print(SEP)
    print(f"Script version: {SCRIPT_VERSION}")
    print()
    print("Scientific guard:")
    print("  New target-expression statistic calculated: NO")
    print("  New corruption mapping generated:           NO")
    print("  Spearman_all_edges calculated:              NO")
    print("  Existing scientific TSV/JSON results read:  NO")
    print("  Operation: inspect source code only")
    print()

    require(CONTRACT.exists(), f"Run the 05j4c freeze first; missing: {CONTRACT}")
    contract = json.loads(CONTRACT.read_text(encoding="utf-8"))
    require(
        contract.get("status") == "FROZEN_05J4C_DIAGNOSTIC_COMPLETION_CONTRACT",
        "05j4c contract is not frozen/authoritative.",
    )

    present = [p for p in CANDIDATE_SOURCES if p.exists()]
    require(
        len(present) >= 2,
        "Too few expected historical source files are present. "
        "This probe refuses to reconstruct the implementation from results.",
    )

    report = {
        "script_version": SCRIPT_VERSION,
        "status": "SOURCE_ONLY_PROBE_COMPLETE",
        "scientific_results_accessed": False,
        "spearman_all_edges_calculated": False,
        "files": [],
    }
    text_blocks = []

    for path in present:
        text = path.read_text(encoding="utf-8", errors="replace")
        lines = text.splitlines()
        suffix = path.suffix.lower()
        functions = python_function_ranges(text) if suffix == ".py" else []
        ranges = merged_relevant_ranges(lines)

        file_rec = {
            "path": str(path),
            "sha256": sha256_file(path),
            "line_count": len(lines),
            "python_functions": functions,
            "relevant_ranges": [{"start": a, "end": b} for a, b in ranges],
        }
        report["files"].append(file_rec)

        text_blocks.append("#" * 168)
        text_blocks.append(f"FILE: {path}")
        text_blocks.append(f"SHA256: {file_rec['sha256']}")
        text_blocks.append(f"Total lines: {len(lines)}")
        if functions:
            text_blocks.append(
                "Top-level/nested Python functions: "
                + repr([f["name"] for f in functions])
            )
        text_blocks.append(
            "Relevant source ranges (source-only; no scientific result files opened):"
        )

        for a, b in ranges:
            text_blocks.append("")
            text_blocks.append(f"--- lines {a}-{b} ---")
            for lineno in range(a, b + 1):
                text_blocks.append(f"{lineno:05d}: {lines[lineno - 1]}")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    require(not OUT_JSON.exists(), f"Refusing to overwrite probe JSON: {OUT_JSON}")
    require(not OUT_TXT.exists(), f"Refusing to overwrite probe text: {OUT_TXT}")
    OUT_JSON.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    OUT_TXT.write_text("\n".join(text_blocks) + "\n", encoding="utf-8")

    print("\n".join(text_blocks))
    print()
    print(SEP)
    print("05j4c ALL-EDGE SPEARMAN SOURCE-ONLY PROBE: COMPLETE")
    print(SEP)
    print("Scientific results accessed:       NO")
    print("Spearman_all_edges calculated:     NO")
    print("Existing mappings changed:         NO")
    print("Ready for exact implementation binding after source inspection: YES")
    print()
    print(f"Probe text: {OUT_TXT}")
    print(f"Probe JSON: {OUT_JSON}")
    print(SEP)

if __name__ == "__main__":
    main()
