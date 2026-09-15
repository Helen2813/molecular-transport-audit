#!/usr/bin/env python
from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

SCRIPT_VERSION = "05j4c-run-all-edge-spearman-diagnostic-v1-no-cli"

PROJECT_ROOT = Path(r"C:\Users\olegk\Desktop\molecular-transport-audit")
DATA_ROOT = Path(r"D:\paper4_tcbb_data")

SOURCE_05J4 = PROJECT_ROOT / "scripts" / "05j4_run_scale_invariant_corruption_audit_v3.py"
EXPECTED_05J4_SHA256 = "37e2286612cf7a708ef99e6de5667073cc8a78314ba5bbc6402bca86cc24792f"

CONTRACT = (
    DATA_ROOT
    / "paper4_tcbb_scale_invariant_diagnostic_completion_contract_v1"
    / "all_edge_spearman_interpretation_contract_v1.json"
)
OLD_PERF = (
    DATA_ROOT
    / "paper4_tcbb_scale_invariant_corruption_inputs_v3"
    / "scale_invariant_performance_grid_v3.tsv"
)
OLD_NULL = (
    DATA_ROOT
    / "paper4_tcbb_scale_invariant_corruption_inputs_v3"
    / "scale_invariant_fullnull_1000_v3.tsv"
)
OLD_METRICS = (
    DATA_ROOT
    / "paper4_tcbb_scale_invariant_corruption_audit_v3"
    / "scale_invariant_operating_metrics_v3.tsv"
)
NEW_PERF = (
    DATA_ROOT
    / "paper4_tcbb_all_edge_spearman_inputs_v1"
    / "all_edge_spearman_performance_grid_v1.tsv"
)
NEW_NULL = (
    DATA_ROOT
    / "paper4_tcbb_all_edge_spearman_inputs_v1"
    / "all_edge_spearman_fullnull_1000_v1.tsv"
)
NEW_MASTER = (
    DATA_ROOT
    / "paper4_tcbb_all_edge_spearman_inputs_v1"
    / "all_edge_spearman_inputs_v1.json"
)

OUT_DIR = DATA_ROOT / "paper4_tcbb_all_edge_spearman_diagnostic_v1"
METRICS_OUT = OUT_DIR / "all_edge_spearman_operating_metrics_v1.tsv"
CONTRASTS_OUT = OUT_DIR / "all_edge_spearman_paired_contrasts_v1.tsv"
WITHIN_OUT = OUT_DIR / "all_edge_spearman_within_fraction_correlations_v1.tsv"
DECOMP_OUT = OUT_DIR / "all_edge_spearman_2x2_decomposition_v1.tsv"
SUMMARY_OUT = OUT_DIR / "all_edge_spearman_2x2_summary_v1.tsv"
MASTER_OUT = OUT_DIR / "all_edge_spearman_diagnostic_v1.json"

SEP = "=" * 168
REPLAY_TOL = 1e-12

def require(cond: bool, msg: str) -> None:
    if not cond:
        raise RuntimeError(msg)

def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()

def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, f"Cannot import {path}")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod

def q(x, p):
    z = np.asarray(x, dtype=float)
    z = z[np.isfinite(z)]
    return float(np.quantile(z, p)) if len(z) else np.nan

def cell_metrics(exact05j4, g, n, col: str) -> dict:
    cal = n.iloc[: exact05j4.CAL_N]
    hold = n.iloc[exact05j4.CAL_N :]
    threshold = exact05j4.threshold_from_calibration(cal[col].to_numpy(dtype=float))
    cal_fpr = float(np.mean(cal[col].to_numpy(dtype=float) > threshold))
    hold_fpr = float(np.mean(hold[col].to_numpy(dtype=float) > threshold))
    require(cal_fpr <= exact05j4.ALPHA + 1e-15, f"{col}: calibration FPR > nominal alpha.")

    rates = []
    for f in exact05j4.FRACTIONS[:-1]:
        gg = g.loc[g["fraction"] == f]
        expected = 1 if f == 0.0 else 100
        require(len(gg) == expected, f"{col}/f={f}: row-count drift.")
        rates.append(float(np.mean(gg[col].to_numpy(dtype=float) > threshold)))
    rates.append(hold_fpr)
    rates = np.asarray(rates, dtype=float)

    f50, fstatus = exact05j4.f50_with_status(exact05j4.FRACTIONS, rates)
    auc = exact05j4.auc_perf(rates)
    return {
        "threshold": float(threshold),
        "calibration_fpr": cal_fpr,
        "holdout_fpr": hold_fpr,
        "rates": rates,
        "f50": float(f50) if np.isfinite(f50) else np.nan,
        "f50_status": fstatus,
        "auc": float(auc),
    }

def main() -> None:
    print(SEP)
    print("Paper 4 / TCBB - all-edge Spearman diagnostic")
    print(SEP)
    print(f"Script version: {SCRIPT_VERSION}")
    print()
    print("Frozen interpretation:")
    print("  Exact original 05j4 threshold/f50/AUC functions reused: YES")
    print("  Pairwise 2x2 decomposition:                         YES")
    print("  Aggregate medians assumed additive:                 NO")
    print("  New inferential bootstrap:                          NO")
    print("  Post-hoc re-slicing permitted:                      NO")
    print("  Existing 05j4 primary inference changed:            NO")
    print()

    for p in [SOURCE_05J4, CONTRACT, OLD_PERF, OLD_NULL, OLD_METRICS, NEW_PERF, NEW_NULL, NEW_MASTER]:
        require(p.exists(), f"Missing required file: {p}")
    require(sha256_file(SOURCE_05J4) == EXPECTED_05J4_SHA256, "05j4 source hash drift.")

    c = json.loads(CONTRACT.read_text(encoding="utf-8"))
    require(c.get("status") == "FROZEN_BEFORE_ALL_EDGE_SPEARMAN_RESULTS", "Unexpected contract status.")
    nm = json.loads(NEW_MASTER.read_text(encoding="utf-8"))
    require(nm.get("status") == "ALL_EDGE_SPEARMAN_INPUTS_COMPLETE", "New input master incomplete.")

    exact = load_module(SOURCE_05J4, "paper4_exact_05j4_for_spearman_diag")

    oldp = pd.read_csv(OLD_PERF, sep="\t")
    oldn = pd.read_csv(OLD_NULL, sep="\t")
    old_metrics = pd.read_csv(OLD_METRICS, sep="\t")
    newp = pd.read_csv(NEW_PERF, sep="\t")
    newn = pd.read_csv(NEW_NULL, sep="\t")

    kp = ["target", "program_id", "fraction", "replicate_id", "attempt_id"]
    kn = ["target", "program_id", "valid_null_id", "attempt_id"]
    require(
        oldp[kp].sort_values(kp).reset_index(drop=True).equals(
            newp[kp].sort_values(kp).reset_index(drop=True)
        ),
        "Performance mapping identity mismatch.",
    )
    require(
        oldn[kn].sort_values(kn).reset_index(drop=True).equals(
            newn[kn].sort_values(kn).reset_index(drop=True)
        ),
        "Full-null mapping identity mismatch.",
    )

    perf = oldp.merge(newp[kp + ["spearman_all_edges"]], on=kp, how="inner", validate="one_to_one")
    full = oldn.merge(newn[kn + ["spearman_all_edges"]], on=kn, how="inner", validate="one_to_one")
    require(len(perf) == 22 * 601, "Merged performance row-count drift.")
    require(len(full) == 22 * 1000, "Merged full-null row-count drift.")

    # A = Spearman subset (MTA), B = Pearson subset, C = Spearman all, D = Pearson all.
    cells = {
        "A_MTA_Spearman_subset": "mta_rho_edge",
        "B_Pearson_subset": "pearson_subset_matched",
        "C_Spearman_all": "spearman_all_edges",
        "D_NetRep_Pearson_all": "netrep_cor_cor",
    }

    metric_rows = []
    contrast_rows = []
    within_rows = []
    decomp_rows = []

    for (target, program_id), g in perf.groupby(["target", "program_id"], sort=True):
        g = g.sort_values(["fraction", "replicate_id"]).reset_index(drop=True)
        n = full.loc[
            (full["target"] == target) & (full["program_id"] == program_id)
        ].sort_values("valid_null_id").reset_index(drop=True)
        require(len(n) == 1000, f"{target}/{program_id}: full-null row-count drift.")

        m = {label: cell_metrics(exact, g, n, col) for label, col in cells.items()}

        # Replay the three already-observed v3 cells exactly against stored metrics.
        for label, method_name in [
            ("A_MTA_Spearman_subset", "MTA"),
            ("B_Pearson_subset", "Pearson_subset_matched"),
            ("D_NetRep_Pearson_all", "NetRep_cor.cor"),
        ]:
            r = old_metrics.loc[
                (old_metrics["target"] == target)
                & (old_metrics["program_id"] == program_id)
                & (old_metrics["method"] == method_name)
            ]
            require(len(r) == 1, f"Stored old metric row missing/duplicated: {target}/{program_id}/{method_name}")
            rr = r.iloc[0]
            require(abs(m[label]["threshold"] - float(rr["threshold"])) <= REPLAY_TOL,
                    f"{target}/{program_id}/{method_name}: threshold replay mismatch.")
            require(abs(m[label]["holdout_fpr"] - float(rr["holdout_false_preservation_rate"])) <= REPLAY_TOL,
                    f"{target}/{program_id}/{method_name}: holdout-FPR replay mismatch.")
            stored_f50 = float(rr["f50"])
            if np.isfinite(stored_f50) or np.isfinite(m[label]["f50"]):
                require(abs(m[label]["f50"] - stored_f50) <= REPLAY_TOL,
                        f"{target}/{program_id}/{method_name}: f50 replay mismatch.")
            require(abs(m[label]["auc"] - float(rr["retention_auc_0_to_0p75"])) <= REPLAY_TOL,
                    f"{target}/{program_id}/{method_name}: AUC replay mismatch.")

        for label, mm in m.items():
            metric_rows.append({
                "target": target,
                "program_id": program_id,
                "cell": label,
                "threshold": mm["threshold"],
                "calibration_false_preservation_rate": mm["calibration_fpr"],
                "holdout_false_preservation_rate": mm["holdout_fpr"],
                "f50": mm["f50"],
                "f50_status": mm["f50_status"],
                "retention_auc_0_to_0p75": mm["auc"],
            })

        A, B, C, D = (m["A_MTA_Spearman_subset"], m["B_Pearson_subset"],
                      m["C_Spearman_all"], m["D_NetRep_Pearson_all"])

        def diff(x, y, key):
            xv, yv = x[key], y[key]
            if key == "f50" and (not np.isfinite(xv) or not np.isfinite(yv)):
                return np.nan
            return float(xv - yv)

        row = {
            "target": target,
            "program_id": program_id,
            # One-factor effects:
            "corr_type_effect_subset_f50_A_minus_B": diff(A, B, "f50"),
            "corr_type_effect_all_f50_C_minus_D": diff(C, D, "f50"),
            "edge_universe_effect_spearman_f50_A_minus_C": diff(A, C, "f50"),
            "edge_universe_effect_pearson_f50_B_minus_D": diff(B, D, "f50"),
            "interaction_f50": diff(A, B, "f50") - diff(C, D, "f50")
                if np.isfinite(diff(A, B, "f50")) and np.isfinite(diff(C, D, "f50")) else np.nan,
            "combined_A_minus_D_f50": diff(A, D, "f50"),
            "corr_type_effect_subset_auc_A_minus_B": diff(A, B, "auc"),
            "corr_type_effect_all_auc_C_minus_D": diff(C, D, "auc"),
            "edge_universe_effect_spearman_auc_A_minus_C": diff(A, C, "auc"),
            "edge_universe_effect_pearson_auc_B_minus_D": diff(B, D, "auc"),
            "interaction_auc": diff(A, B, "auc") - diff(C, D, "auc"),
            "combined_A_minus_D_auc": diff(A, D, "auc"),
        }
        decomp_rows.append(row)

        contrast_rows.extend([
            {"target": target, "program_id": program_id, "contrast": "A_minus_C_MTA_vs_Spearman_all",
             "delta_f50": row["edge_universe_effect_spearman_f50_A_minus_C"],
             "delta_auc": row["edge_universe_effect_spearman_auc_A_minus_C"]},
            {"target": target, "program_id": program_id, "contrast": "A_minus_B_MTA_vs_Pearson_subset",
             "delta_f50": row["corr_type_effect_subset_f50_A_minus_B"],
             "delta_auc": row["corr_type_effect_subset_auc_A_minus_B"]},
            {"target": target, "program_id": program_id, "contrast": "A_minus_D_MTA_vs_NetRep",
             "delta_f50": row["combined_A_minus_D_f50"],
             "delta_auc": row["combined_A_minus_D_auc"]},
        ])

        # Within-fraction rank concordance with MTA for the new cell.
        for f in [0.05, 0.10, 0.15, 0.25, 0.50, 0.75]:
            gg = g.loc[g["fraction"] == f]
            rho = float(spearmanr(
                gg["mta_rho_edge"].to_numpy(float),
                gg["spearman_all_edges"].to_numpy(float),
            ).statistic)
            within_rows.append({
                "target": target, "program_id": program_id, "fraction": f,
                "spearman_mta_vs_spearman_all_edges": rho, "n": 100,
            })
        hold = n.iloc[exact.CAL_N:]
        within_rows.append({
            "target": target, "program_id": program_id, "fraction": 1.0,
            "spearman_mta_vs_spearman_all_edges": float(spearmanr(
                hold["mta_rho_edge"].to_numpy(float),
                hold["spearman_all_edges"].to_numpy(float),
            ).statistic),
            "n": exact.HOLDOUT_N,
        })

    metrics = pd.DataFrame(metric_rows)
    contrasts = pd.DataFrame(contrast_rows)
    within = pd.DataFrame(within_rows)
    decomp = pd.DataFrame(decomp_rows)

    summary_rows = []
    for col in [
        "corr_type_effect_subset_f50_A_minus_B",
        "corr_type_effect_all_f50_C_minus_D",
        "edge_universe_effect_spearman_f50_A_minus_C",
        "edge_universe_effect_pearson_f50_B_minus_D",
        "interaction_f50",
        "combined_A_minus_D_f50",
        "corr_type_effect_subset_auc_A_minus_B",
        "corr_type_effect_all_auc_C_minus_D",
        "edge_universe_effect_spearman_auc_A_minus_C",
        "edge_universe_effect_pearson_auc_B_minus_D",
        "interaction_auc",
        "combined_A_minus_D_auc",
    ]:
        z = decomp[col].to_numpy(float)
        finite = z[np.isfinite(z)]
        summary_rows.append({
            "quantity": col,
            "n_pairs": len(z),
            "finite_pairs": len(finite),
            "median": float(np.median(finite)) if len(finite) else np.nan,
            "q25": q(finite, 0.25),
            "q75": q(finite, 0.75),
            "negative_pairs": int(np.sum(finite < 0)),
            "zero_pairs": int(np.sum(np.isclose(finite, 0.0))),
            "positive_pairs": int(np.sum(finite > 0)),
        })
    summary = pd.DataFrame(summary_rows)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    for path in [METRICS_OUT, CONTRASTS_OUT, WITHIN_OUT, DECOMP_OUT, SUMMARY_OUT, MASTER_OUT]:
        require(not path.exists(), f"Refusing to overwrite existing result: {path}")

    metrics.to_csv(METRICS_OUT, sep="\t", index=False)
    contrasts.to_csv(CONTRASTS_OUT, sep="\t", index=False)
    within.to_csv(WITHIN_OUT, sep="\t", index=False)
    decomp.to_csv(DECOMP_OUT, sep="\t", index=False)
    summary.to_csv(SUMMARY_OUT, sep="\t", index=False)

    master = {
        "script_version": SCRIPT_VERSION,
        "status": "ALL_EDGE_SPEARMAN_DIAGNOSTIC_COMPLETE",
        "scientific_role": "pre-specified descriptive completion of the 2x2 decomposition",
        "cell_labels": cells,
        "old_three_cells_exact_metric_replay": "PASS",
        "new_inferential_bootstrap": False,
        "existing_05j4_primary_inference_changed": False,
        "outputs": {
            "metrics": str(METRICS_OUT),
            "contrasts": str(CONTRASTS_OUT),
            "within_fraction": str(WITHIN_OUT),
            "decomposition": str(DECOMP_OUT),
            "summary": str(SUMMARY_OUT),
        },
    }
    MASTER_OUT.write_text(json.dumps(master, indent=2) + "\n", encoding="utf-8")

    print(SEP)
    print("05j4c ALL-EDGE SPEARMAN DIAGNOSTIC: COMPLETE")
    print(SEP)
    print("\n2x2 descriptive summaries:")
    print(summary.to_string(index=False))
    print("\nNew C-cell operating metrics:")
    print(
        metrics.loc[metrics["cell"] == "C_Spearman_all",
                    ["target", "program_id", "calibration_false_preservation_rate",
                     "holdout_false_preservation_rate", "f50",
                     "retention_auc_0_to_0p75"]].to_string(index=False)
    )
    print()
    print("Old A/B/D cells exact replay: PASS")
    print("New inferential bootstrap: NO")
    print("Existing 05j4 primary inference changed: NO")
    print(f"Summary:       {SUMMARY_OUT}")
    print(f"Decomposition: {DECOMP_OUT}")
    print(f"Master:        {MASTER_OUT}")
    print(SEP)

if __name__ == "__main__":
    main()
