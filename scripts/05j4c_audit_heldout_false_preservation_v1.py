#!/usr/bin/env python
from __future__ import annotations

import hashlib
import json
import math
import re
from pathlib import Path

import numpy as np
import pandas as pd

SCRIPT_VERSION = "05j4c-audit-heldout-false-preservation-v1-no-cli"

DATA_ROOT = Path(r"D:\paper4_tcbb_data")

CONTRACT = (
    DATA_ROOT
    / "paper4_tcbb_scale_invariant_diagnostic_completion_contract_v1"
    / "scale_invariant_diagnostic_completion_contract_v1.json"
)
FULLNULL = (
    DATA_ROOT
    / "paper4_tcbb_scale_invariant_corruption_inputs_v3"
    / "scale_invariant_fullnull_1000_v3.tsv"
)
THRESHOLDS = (
    DATA_ROOT
    / "paper4_tcbb_scale_invariant_corruption_audit_v3"
    / "scale_invariant_thresholds_v3.tsv"
)

OUT_DIR = DATA_ROOT / "paper4_tcbb_scale_invariant_heldout_calibration_audit_v1"
BY_PAIR = OUT_DIR / "heldout_false_preservation_by_pair_v1.tsv"
PAIRED = OUT_DIR / "heldout_false_preservation_paired_differences_v1.tsv"
METHOD_SUMMARY = OUT_DIR / "heldout_false_preservation_method_summary_v1.tsv"
PAIRED_SUMMARY = OUT_DIR / "heldout_false_preservation_paired_summary_v1.tsv"
SCHEMA_REPORT = OUT_DIR / "heldout_false_preservation_schema_binding_v1.json"
MASTER = OUT_DIR / "heldout_false_preservation_audit_v1.json"

SEP = "=" * 168

CAL_N = 800
HOLD_N = 200
EXPECTED_FULL_N = 1000
TARGET_FPR = 0.05

METHOD_CANONICAL = {
    "MTA": {
        "fullnull_aliases": [
            "MTA", "MTA_rho_edge", "rho_edge", "mta_score", "mta_rho",
            "mta", "mta.rho_edge",
        ],
        "threshold_method_aliases": [
            "MTA", "MTA_rho_edge", "rho_edge", "mta_score", "mta_rho",
        ],
    },
    "NetRep_cor.cor": {
        "fullnull_aliases": [
            "NetRep_cor.cor", "NetRep_cor_cor", "netrep_cor_cor",
            "netrep.cor.cor", "cor.cor", "netrep_corcor",
        ],
        "threshold_method_aliases": [
            "NetRep_cor.cor", "NetRep_cor_cor", "netrep_cor_cor",
            "netrep.cor.cor", "cor.cor", "netrep_corcor",
        ],
    },
    "WGCNA_cor.kIM": {
        "fullnull_aliases": [
            "WGCNA_cor.kIM", "WGCNA_cor_kIM", "wgcna_cor_kim",
            "wgcna.cor.kIM", "cor.kIM", "wgcna_corkim",
        ],
        "threshold_method_aliases": [
            "WGCNA_cor.kIM", "WGCNA_cor_kIM", "wgcna_cor_kim",
            "wgcna.cor.kIM", "cor.kIM", "wgcna_corkim",
        ],
    },
    "Pearson_subset_matched": {
        "fullnull_aliases": [
            "Pearson_subset_matched", "pearson_subset_matched",
            "pearson_matched", "same_edge_pearson", "matched_pearson",
        ],
        "threshold_method_aliases": [
            "Pearson_subset_matched", "pearson_subset_matched",
            "pearson_matched", "same_edge_pearson", "matched_pearson",
        ],
    },
}

PAIR_KEY_ALIASES = {
    "target": ["target", "target_name", "target_dataset", "cohort", "dataset"],
    "program_id": ["program_id", "module_id", "program", "module"],
}
ORDER_ALIASES = [
    "replicate", "replicate_id", "mapping_id", "fullnull_id", "null_id",
    "valid_index", "valid_mapping_index", "draw", "iteration",
]
LONG_METHOD_ALIASES = ["method", "statistic", "comparator", "score_name"]
LONG_SCORE_ALIASES = ["score", "value", "statistic_value", "observed", "raw_score"]
THRESHOLD_VALUE_ALIASES = [
    "threshold", "threshold_value", "preservation_threshold", "cutoff", "cut_point",
]
THRESHOLD_METHOD_ALIASES = ["method", "statistic", "comparator", "score_name"]

def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)

def norm(s: object) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(s).lower())

def find_unique_col(columns, aliases, label: str, required: bool = True):
    by_norm = {}
    for c in columns:
        by_norm.setdefault(norm(c), []).append(c)
    hits = []
    for a in aliases:
        hits.extend(by_norm.get(norm(a), []))
    hits = list(dict.fromkeys(hits))
    if len(hits) == 1:
        return hits[0]
    if not hits and not required:
        return None
    if not hits:
        raise RuntimeError(
            f"Could not bind required {label} column. "
            f"Expected one of {aliases}; observed columns={list(columns)}"
        )
    raise RuntimeError(f"Ambiguous {label} binding: {hits}")

def bind_pair_keys(df: pd.DataFrame):
    return {
        k: find_unique_col(df.columns, aliases, f"{k}")
        for k, aliases in PAIR_KEY_ALIASES.items()
    }

def canonical_method(value: object) -> str | None:
    nv = norm(value)
    for canonical, cfg in METHOD_CANONICAL.items():
        aliases = cfg["threshold_method_aliases"] + cfg["fullnull_aliases"] + [canonical]
        if nv in {norm(x) for x in aliases}:
            return canonical
    return None

def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()

def quantile_lower(values: np.ndarray, q: float) -> float:
    # This helper is NOT used to re-estimate thresholds.
    # It is retained only for optional diagnostic display if needed later.
    try:
        return float(np.quantile(values, q, method="lower"))
    except TypeError:
        return float(np.quantile(values, q, interpolation="lower"))

def load_thresholds() -> tuple[pd.DataFrame, dict]:
    t = pd.read_csv(THRESHOLDS, sep="\t")
    pair = bind_pair_keys(t)
    method_col = find_unique_col(t.columns, THRESHOLD_METHOD_ALIASES, "threshold method")
    threshold_col = find_unique_col(t.columns, THRESHOLD_VALUE_ALIASES, "threshold value")

    t = t.copy()
    t["_method"] = t[method_col].map(canonical_method)
    unknown = sorted(set(t.loc[t["_method"].isna(), method_col].astype(str)))
    require(
        not unknown,
        "Threshold table contains method labels that cannot be bound safely: "
        + repr(unknown),
    )
    t["_target"] = t[pair["target"]].astype(str)
    t["_program_id"] = t[pair["program_id"]].astype(str)
    t["_threshold"] = pd.to_numeric(t[threshold_col], errors="coerce")
    require(t["_threshold"].notna().all(), "Non-numeric frozen threshold detected.")

    dups = t.duplicated(["_target", "_program_id", "_method"], keep=False)
    require(
        not dups.any(),
        "Frozen threshold table is not unique by target/program/method:\n"
        + t.loc[dups, ["_target", "_program_id", "_method", "_threshold"]].to_string(index=False),
    )
    return t, {
        "target_col": pair["target"],
        "program_id_col": pair["program_id"],
        "method_col": method_col,
        "threshold_col": threshold_col,
    }

def load_fullnull_long() -> tuple[pd.DataFrame, dict]:
    f = pd.read_csv(FULLNULL, sep="\t")
    pair = bind_pair_keys(f)
    order_col = find_unique_col(f.columns, ORDER_ALIASES, "full-null replicate/order")

    method_col = find_unique_col(f.columns, LONG_METHOD_ALIASES, "long-form method", required=False)
    score_col = find_unique_col(f.columns, LONG_SCORE_ALIASES, "long-form score", required=False)

    binding = {
        "target_col": pair["target"],
        "program_id_col": pair["program_id"],
        "order_col": order_col,
        "layout": None,
        "method_col": method_col,
        "score_col": score_col,
        "wide_score_columns": {},
    }

    if method_col is not None and score_col is not None:
        x = f[[pair["target"], pair["program_id"], order_col, method_col, score_col]].copy()
        x.columns = ["_target", "_program_id", "_order", "_method_raw", "_score"]
        x["_method"] = x["_method_raw"].map(canonical_method)
        unknown = sorted(set(x.loc[x["_method"].isna(), "_method_raw"].astype(str)))
        require(
            not unknown,
            "Full-null long table contains method labels that cannot be bound safely: "
            + repr(unknown),
        )
        x["_score"] = pd.to_numeric(x["_score"], errors="coerce")
        binding["layout"] = "long"
    else:
        wide_cols = {}
        for canonical, cfg in METHOD_CANONICAL.items():
            col = find_unique_col(
                f.columns,
                cfg["fullnull_aliases"] + [canonical],
                f"wide full-null score for {canonical}",
            )
            wide_cols[canonical] = col

        pieces = []
        for canonical, score_name in wide_cols.items():
            p = f[[pair["target"], pair["program_id"], order_col, score_name]].copy()
            p.columns = ["_target", "_program_id", "_order", "_score"]
            p["_method"] = canonical
            pieces.append(p)
        x = pd.concat(pieces, ignore_index=True)
        x["_score"] = pd.to_numeric(x["_score"], errors="coerce")
        binding["layout"] = "wide"
        binding["wide_score_columns"] = wide_cols

    x["_target"] = x["_target"].astype(str)
    x["_program_id"] = x["_program_id"].astype(str)
    require(x["_score"].notna().all(), "Full-null table contains non-numeric/NA scores.")

    # Preserve numeric or lexical mapping identity, but never silently use file row order.
    numeric_order = pd.to_numeric(x["_order"], errors="coerce")
    if numeric_order.notna().all():
        x["_order_sort"] = numeric_order.astype(float)
        binding["order_type"] = "numeric"
    else:
        x["_order_sort"] = x["_order"].astype(str)
        binding["order_type"] = "lexical"

    dups = x.duplicated(["_target", "_program_id", "_method", "_order"], keep=False)
    require(
        not dups.any(),
        "Duplicate target/program/method/full-null mapping identities detected; refusing to infer order.",
    )

    return x, binding

def summarize_numeric(a: pd.Series) -> dict:
    z = pd.to_numeric(a, errors="coerce").dropna().to_numpy(float)
    require(len(z) > 0, "Cannot summarize an empty numeric vector.")
    return {
        "n": int(len(z)),
        "median": float(np.median(z)),
        "q25": float(np.quantile(z, 0.25)),
        "q75": float(np.quantile(z, 0.75)),
        "min": float(np.min(z)),
        "max": float(np.max(z)),
    }

def main() -> None:
    print(SEP)
    print("Paper 4 / TCBB - held-out false-preservation audit for frozen 05j4 v3 thresholds")
    print(SEP)
    print(f"Script version: {SCRIPT_VERSION}")
    print()
    print("Execution contract:")
    print("  Existing 05j4 v3 thresholds reused unchanged: YES")
    print("  First 800 / final 200 split reused:            YES")
    print("  Threshold recalibration here:                  NO")
    print("  Pairwise within module/cohort comparisons:     YES")
    print("  New inferential p-value/bootstrap:             NO")
    print()

    for p in [CONTRACT, FULLNULL, THRESHOLDS]:
        require(p.exists(), f"Missing required frozen artifact: {p}")

    c = json.loads(CONTRACT.read_text(encoding="utf-8"))
    require(
        c.get("status") == "FROZEN_05J4C_DIAGNOSTIC_COMPLETION_CONTRACT",
        "05j4c freeze contract status is not authoritative.",
    )
    require(
        c["heldout_false_preservation_audit"]["calibration_mappings"] == CAL_N
        and c["heldout_false_preservation_audit"]["heldout_mappings"] == HOLD_N,
        "05j4c contract split does not match this runner.",
    )

    thresholds, threshold_binding = load_thresholds()
    fullnull, fullnull_binding = load_fullnull_long()

    expected_methods = set(METHOD_CANONICAL)
    observed_methods = set(fullnull["_method"].unique())
    require(
        expected_methods <= observed_methods,
        f"Full-null methods missing: {sorted(expected_methods - observed_methods)}",
    )

    rows = []
    grouped = fullnull.groupby(["_target", "_program_id", "_method"], sort=True)
    for (target, program_id, method), g in grouped:
        if method not in expected_methods:
            continue
        g = g.sort_values("_order_sort", kind="mergesort").reset_index(drop=True)
        require(
            len(g) == EXPECTED_FULL_N,
            f"{target}/{program_id}/{method}: expected {EXPECTED_FULL_N} full-null mappings, got {len(g)}",
        )

        tr = thresholds[
            (thresholds["_target"] == target)
            & (thresholds["_program_id"] == program_id)
            & (thresholds["_method"] == method)
        ]
        require(
            len(tr) == 1,
            f"{target}/{program_id}/{method}: expected exactly one frozen threshold, got {len(tr)}",
        )
        threshold = float(tr.iloc[0]["_threshold"])

        cal = g.iloc[:CAL_N]
        hold = g.iloc[CAL_N:CAL_N + HOLD_N]
        require(len(cal) == CAL_N and len(hold) == HOLD_N, "Split length failure.")

        # All four preservation statistics in this audit are positive-direction
        # structural-preservation scores: larger = more preserved.
        # We still guard this explicitly rather than auto-flipping from observed data.
        cal_calls = cal["_score"].to_numpy(float) >= threshold
        hold_calls = hold["_score"].to_numpy(float) >= threshold

        cal_n_pres = int(np.sum(cal_calls))
        hold_n_pres = int(np.sum(hold_calls))
        cal_rate = cal_n_pres / CAL_N
        hold_rate = hold_n_pres / HOLD_N

        rows.append({
            "target": target,
            "program_id": program_id,
            "method": method,
            "threshold_frozen_v3": threshold,
            "score_orientation": "higher_is_more_preserved",
            "calibration_n": CAL_N,
            "calibration_false_preserved_n": cal_n_pres,
            "calibration_false_preservation_rate": cal_rate,
            "heldout_n": HOLD_N,
            "heldout_false_preserved_n": hold_n_pres,
            "heldout_false_preservation_rate": hold_rate,
            "heldout_minus_calibration_rate": hold_rate - cal_rate,
            "heldout_minus_nominal_0p05": hold_rate - TARGET_FPR,
        })

    by_pair = pd.DataFrame(rows).sort_values(["target", "program_id", "method"]).reset_index(drop=True)

    # Ensure each target-program pair has exactly the four predeclared methods.
    counts = by_pair.groupby(["target", "program_id"])["method"].nunique()
    require(
        (counts == len(expected_methods)).all(),
        "One or more target-program pairs do not contain all four predeclared v3 methods.",
    )
    require(len(counts) == 22, f"Expected 22 target-program pairs, found {len(counts)}.")

    # Descriptive method summaries.
    method_summary_rows = []
    for method, g in by_pair.groupby("method", sort=True):
        s = summarize_numeric(g["heldout_false_preservation_rate"])
        method_summary_rows.append({
            "method": method,
            "n_pairs": len(g),
            "median_heldout_fpr": s["median"],
            "q25_heldout_fpr": s["q25"],
            "q75_heldout_fpr": s["q75"],
            "min_heldout_fpr": s["min"],
            "max_heldout_fpr": s["max"],
            "median_calibration_fpr": float(np.median(g["calibration_false_preservation_rate"])),
            "pairs_heldout_gt_0p05": int((g["heldout_false_preservation_rate"] > TARGET_FPR).sum()),
            "pairs_heldout_eq_0p05": int(np.isclose(g["heldout_false_preservation_rate"], TARGET_FPR).sum()),
            "pairs_heldout_lt_0p05": int((g["heldout_false_preservation_rate"] < TARGET_FPR).sum()),
        })
    method_summary = pd.DataFrame(method_summary_rows)

    # Paired MTA-minus-comparator held-out differences.
    wide = by_pair.pivot(
        index=["target", "program_id"],
        columns="method",
        values="heldout_false_preservation_rate",
    ).reset_index()

    paired_rows = []
    for comparator in ["NetRep_cor.cor", "WGCNA_cor.kIM", "Pearson_subset_matched"]:
        for _, r in wide.iterrows():
            diff = float(r["MTA"] - r[comparator])
            paired_rows.append({
                "target": r["target"],
                "program_id": r["program_id"],
                "comparator": comparator,
                "mta_heldout_fpr": float(r["MTA"]),
                "comparator_heldout_fpr": float(r[comparator]),
                "paired_difference_mta_minus_comparator": diff,
            })
    paired = pd.DataFrame(paired_rows).sort_values(["comparator", "target", "program_id"]).reset_index(drop=True)

    paired_summary_rows = []
    for comparator, g in paired.groupby("comparator", sort=True):
        d = g["paired_difference_mta_minus_comparator"].to_numpy(float)
        paired_summary_rows.append({
            "comparator": comparator,
            "n_pairs": len(g),
            "median_mta_minus_comparator_fpr": float(np.median(d)),
            "q25_mta_minus_comparator_fpr": float(np.quantile(d, 0.25)),
            "q75_mta_minus_comparator_fpr": float(np.quantile(d, 0.75)),
            "min_mta_minus_comparator_fpr": float(np.min(d)),
            "max_mta_minus_comparator_fpr": float(np.max(d)),
            "mta_higher_fpr_pairs": int(np.sum(d > 0)),
            "equal_fpr_pairs": int(np.sum(np.isclose(d, 0.0))),
            "mta_lower_fpr_pairs": int(np.sum(d < 0)),
            "interpretation": "DESCRIPTIVE_ONLY_NO_RECALIBRATION",
        })
    paired_summary = pd.DataFrame(paired_summary_rows)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    for p in [BY_PAIR, PAIRED, METHOD_SUMMARY, PAIRED_SUMMARY, SCHEMA_REPORT, MASTER]:
        require(not p.exists(), f"Refusing to overwrite existing 05j4c result: {p}")

    by_pair.to_csv(BY_PAIR, sep="\t", index=False)
    paired.to_csv(PAIRED, sep="\t", index=False)
    method_summary.to_csv(METHOD_SUMMARY, sep="\t", index=False)
    paired_summary.to_csv(PAIRED_SUMMARY, sep="\t", index=False)

    schema = {
        "script_version": SCRIPT_VERSION,
        "fullnull_binding": fullnull_binding,
        "threshold_binding": threshold_binding,
        "score_orientation": "higher_is_more_preserved",
        "calibration_n": CAL_N,
        "heldout_n": HOLD_N,
        "input_sha256": {
            str(CONTRACT): sha256_file(CONTRACT),
            str(FULLNULL): sha256_file(FULLNULL),
            str(THRESHOLDS): sha256_file(THRESHOLDS),
        },
    }
    SCHEMA_REPORT.write_text(json.dumps(schema, indent=2) + "\n", encoding="utf-8")

    master = {
        "script_version": SCRIPT_VERSION,
        "status": "HELDOUT_FALSE_PRESERVATION_AUDIT_COMPLETE",
        "scientific_role": "descriptive calibration validation; no threshold recalibration",
        "n_target_program_pairs": int(len(counts)),
        "methods": sorted(expected_methods),
        "calibration_n_per_pair_method": CAL_N,
        "heldout_n_per_pair_method": HOLD_N,
        "nominal_false_preservation_target": TARGET_FPR,
        "thresholds_reestimated": False,
        "new_inferential_test_performed": False,
        "outputs": {
            "by_pair": str(BY_PAIR),
            "paired_differences": str(PAIRED),
            "method_summary": str(METHOD_SUMMARY),
            "paired_summary": str(PAIRED_SUMMARY),
            "schema_binding": str(SCHEMA_REPORT),
        },
        "input_sha256": schema["input_sha256"],
    }
    MASTER.write_text(json.dumps(master, indent=2) + "\n", encoding="utf-8")

    print(SEP)
    print("05j4c HELD-OUT FALSE-PRESERVATION AUDIT: COMPLETE")
    print(SEP)
    print()
    print("Method-level held-out FPR summaries:")
    print(method_summary.to_string(index=False))
    print()
    print("Paired held-out FPR differences (MTA - comparator):")
    print(paired_summary.to_string(index=False))
    print()
    print("Scientific guard:")
    print("  Frozen thresholds re-estimated: NO")
    print("  New p-value/bootstrap computed: NO")
    print("  Any imbalance must be reported, not tuned away: YES")
    print()
    print(f"By-pair:       {BY_PAIR}")
    print(f"Paired:        {PAIRED}")
    print(f"Method summary:{METHOD_SUMMARY}")
    print(f"Pair summary:  {PAIRED_SUMMARY}")
    print(f"Master:        {MASTER}")
    print(SEP)

if __name__ == "__main__":
    main()
