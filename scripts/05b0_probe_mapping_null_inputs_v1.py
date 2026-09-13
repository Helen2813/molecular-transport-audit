from __future__ import annotations

import json
from pathlib import Path

import pandas as pd


SCRIPT_VERSION = "05b0-probe-mapping-null-inputs-v1-no-cli"

DATA_ROOT = Path(r"D:\paper4_tcbb_data")

DIRS = {
    "03d mapping-null contract": DATA_ROOT / "paper4_tcbb_mapping_null_contract_v1",
    "03e target marginal metrics": DATA_ROOT / "paper4_tcbb_target_marginal_metrics_v2",
    "03g pilot audit": DATA_ROOT / "paper4_tcbb_mapping_null_pilot_audit_v1",
    "03h final execution contract": DATA_ROOT / "paper4_tcbb_final_mapping_null_execution_v1",
    "04i evaluability amendment": DATA_ROOT / "paper4_tcbb_target_statistical_evaluability_amendment_v2",
    "05a direct preservation": DATA_ROOT / "paper4_tcbb_primary_pooled_direct_preservation_v1",
}


def human_size(n: int) -> str:
    if n >= 1024**3:
        return f"{n / 1024**3:.2f} GiB"
    if n >= 1024**2:
        return f"{n / 1024**2:.2f} MiB"
    if n >= 1024:
        return f"{n / 1024:.2f} KiB"
    return f"{n} B"


def preview_tsv(path: Path) -> None:
    try:
        df = pd.read_csv(path, sep="\t", dtype=str, nrows=3).fillna("")
    except Exception as exc:
        print(f"      TSV preview FAILED: {exc}")
        return
    print(f"      columns ({len(df.columns)}): {list(df.columns)}")
    for rec in df.to_dict(orient="records"):
        print(f"        {rec}")


def preview_json(path: Path) -> None:
    try:
        obj = json.loads(path.read_text(encoding="utf-8"))
    except Exception as exc:
        print(f"      JSON preview FAILED: {exc}")
        return

    if isinstance(obj, dict):
        print(f"      top-level keys: {list(obj.keys())}")
        # Print small scalar/short-list values only.
        for k, v in obj.items():
            if isinstance(v, (str, int, float, bool)) or v is None:
                print(f"        {k}: {v}")
            elif isinstance(v, list) and len(v) <= 10 and all(
                isinstance(x, (str, int, float, bool)) or x is None for x in v
            ):
                print(f"        {k}: {v}")
    else:
        print(f"      JSON root type: {type(obj).__name__}")


def main() -> None:
    print("=" * 132)
    print("Paper 4 / TCBB - probe exact inputs for frozen 1,000-panel matched-mapping specificity null")
    print("=" * 132)
    print(f"Script version: {SCRIPT_VERSION}")
    print()
    print("Scientific guard:")
    print("  Target expression values read:                       NO")
    print("  Target preservation statistics newly calculated:     NO")
    print("  Mapping-null panels generated:                        NO")
    print("  Operation: filesystem/schema inspection only")
    print("=" * 132)

    for label, d in DIRS.items():
        print()
        print(label)
        print(f"  directory: {d}")

        if not d.exists():
            print("  MISSING")
            continue

        files = sorted([p for p in d.rglob("*") if p.is_file()])
        print(f"  files: {len(files)}")

        for p in files:
            rel = p.relative_to(d)
            print(f"    {rel}  [{human_size(p.stat().st_size)}]")

            suffix = p.suffix.lower()
            if suffix in {".tsv", ".csv"}:
                if suffix == ".tsv":
                    preview_tsv(p)
                else:
                    try:
                        df = pd.read_csv(p, dtype=str, nrows=3).fillna("")
                        print(f"      columns ({len(df.columns)}): {list(df.columns)}")
                        for rec in df.to_dict(orient="records"):
                            print(f"        {rec}")
                    except Exception as exc:
                        print(f"      CSV preview FAILED: {exc}")
            elif suffix == ".json":
                preview_json(p)

    print()
    print("=" * 132)
    print("05b0 MAPPING-NULL INPUT PROBE: COMPLETE")
    print("=" * 132)
    print("No new scientific statistic was calculated.")
    print("Paste this output back; it is sufficient to build the exact 05b specificity-null runner.")
    print("=" * 132)


if __name__ == "__main__":
    main()
