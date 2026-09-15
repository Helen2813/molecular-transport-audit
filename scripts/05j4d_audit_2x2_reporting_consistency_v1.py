#!/usr/bin/env python
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr, skew

SCRIPT_VERSION = "05j4d-audit-2x2-reporting-consistency-v1-no-cli"

DATA_ROOT = Path(r"D:\paper4_tcbb_data")

INPUT_DIR = DATA_ROOT / "paper4_tcbb_all_edge_spearman_diagnostic_v1"
METRICS = INPUT_DIR / "all_edge_spearman_operating_metrics_v1.tsv"
DECOMP = INPUT_DIR / "all_edge_spearman_2x2_decomposition_v1.tsv"
DIAG_MASTER = INPUT_DIR / "all_edge_spearman_diagnostic_v1.json"

OUT_DIR = DATA_ROOT / "paper4_tcbb_2x2_reporting_audit_v1"
PAIRWISE_OUT = OUT_DIR / "pairwise_2x2_reconstruction_v1.tsv"
SUMMARY_OUT = OUT_DIR / "2x2_reporting_summary_v1.tsv"
TARGET_OUT = OUT_DIR / "2x2_reporting_summary_by_target_v1.tsv"
IDENTITY_OUT = OUT_DIR / "2x2_mean_identity_audit_v1.tsv"
DEPENDENCE_OUT = OUT_DIR / "2x2_component_dependence_v1.tsv"
MASTER_OUT = OUT_DIR / "2x2_reporting_audit_v1.json"

SEP = "=" * 168
TOL = 1e-12

# A = Spearman x frozen matched subset (MTA)
# B = Pearson x frozen matched subset
# C = Spearman x all edges
# D = Pearson x all edges, cor.cor-equivalent common-benchmark implementation
CELL_A = "A_MTA_Spearman_subset"
CELL_B = "B_Pearson_subset"
CELL_C = "C_Spearman_all"
CELL_D = "D_NetRep_Pearson_all"

METRIC_COLUMNS = {
    "f50": "f50",
    "auc": "retention_auc_0_to_0p75",
}

DECOMP_EXPECTED = {
    "f50": {
        "AB": "corr_type_effect_subset_f50_A_minus_B",
        "CD": "corr_type_effect_all_f50_C_minus_D",
        "AC": "edge_universe_effect_spearman_f50_A_minus_C",
        "BD": "edge_universe_effect_pearson_f50_B_minus_D",
        "interaction": "interaction_f50",
        "AD": "combined_A_minus_D_f50",
    },
    "auc": {
        "AB": "corr_type_effect_subset_auc_A_minus_B",
        "CD": "corr_type_effect_all_auc_C_minus_D",
        "AC": "edge_universe_effect_spearman_auc_A_minus_C",
        "BD": "edge_universe_effect_pearson_auc_B_minus_D",
        "interaction": "interaction_auc",
        "AD": "combined_A_minus_D_auc",
    },
}

def require(cond: bool, msg: str) -> None:
    if not cond:
        raise RuntimeError(msg)

def q(x: np.ndarray, p: float) -> float:
    z = np.asarray(x, dtype=float)
    z = z[np.isfinite(z)]
    return float(np.quantile(z, p)) if len(z) else np.nan

def summarize(values: np.ndarray) -> dict:
    z = np.asarray(values, dtype=float)
    z = z[np.isfinite(z)]
    require(len(z) > 0, "Cannot summarize empty vector.")
    return {
        "n": int(len(z)),
        "mean": float(np.mean(z)),
        "median": float(np.median(z)),
        "sd": float(np.std(z, ddof=1)) if len(z) > 1 else 0.0,
        "q25": q(z, 0.25),
        "q75": q(z, 0.75),
        "min": float(np.min(z)),
        "max": float(np.max(z)),
        "skewness_descriptive": float(skew(z, bias=False)) if len(z) >= 3 else np.nan,
        "negative_pairs": int(np.sum(z < 0)),
        "zero_pairs": int(np.sum(np.isclose(z, 0.0, atol=TOL, rtol=0.0))),
        "positive_pairs": int(np.sum(z > 0)),
    }

def finite_corr(x: np.ndarray, y: np.ndarray) -> tuple[float, float, int]:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    xx, yy = x[mask], y[mask]
    require(len(xx) >= 3, "Too few finite pairs for descriptive component correlation.")
    return (
        float(pearsonr(xx, yy).statistic),
        float(spearmanr(xx, yy).statistic),
        int(len(xx)),
    )

def main() -> None:
    print(SEP)
    print("Paper 4 / TCBB - deterministic 2x2 reporting consistency audit")
    print(SEP)
    print(f"Script version: {SCRIPT_VERSION}")
    print()
    print("Scientific guard:")
    print("  New biological/statistical experiment:            NO")
    print("  New threshold or corruption grid:                 NO")
    print("  New p-value / CI / bootstrap:                     NO")
    print("  Operation: algebraic/reporting audit of frozen 2x2 results")
    print("  Aggregate medians assumed additive:               NO")
    print("  Pairwise A/B/C/D identities required:             YES")
    print()

    for p in [METRICS, DECOMP, DIAG_MASTER]:
        require(p.exists(), f"Missing required frozen diagnostic artifact: {p}")

    master = json.loads(DIAG_MASTER.read_text(encoding="utf-8"))
    require(
        master.get("status") == "ALL_EDGE_SPEARMAN_DIAGNOSTIC_COMPLETE",
        "Unexpected all-edge Spearman diagnostic master status.",
    )

    metrics = pd.read_csv(METRICS, sep="\t")
    decomp_saved = pd.read_csv(DECOMP, sep="\t")

    required_metric_cols = {
        "target", "program_id", "cell", "f50", "retention_auc_0_to_0p75"
    }
    require(
        required_metric_cols.issubset(metrics.columns),
        f"Metrics missing columns: {sorted(required_metric_cols - set(metrics.columns))}",
    )
    require(
        len(metrics) == 22 * 4,
        f"Expected 88 cell-metric rows (22x4), found {len(metrics)}.",
    )

    expected_cells = {CELL_A, CELL_B, CELL_C, CELL_D}
    observed_cells = set(metrics["cell"].astype(str))
    require(
        observed_cells == expected_cells,
        f"Unexpected 2x2 cells: observed={sorted(observed_cells)}",
    )

    counts = metrics.groupby(["target", "program_id"])["cell"].nunique()
    require((counts == 4).all(), "Not every target/program pair has exactly four 2x2 cells.")
    require(len(counts) == 22, f"Expected 22 target/program pairs, found {len(counts)}.")

    # Build one row per pair with raw cell values.
    pair_keys = ["target", "program_id"]
    pairwise = None

    for metric_name, value_col in METRIC_COLUMNS.items():
        p = (
            metrics.pivot(index=pair_keys, columns="cell", values=value_col)
            .reset_index()
        )
        require(not p[list(expected_cells)].isna().any().any(),
                f"Missing/non-finite cell value in {metric_name} pivot.")

        A = p[CELL_A].to_numpy(float)
        B = p[CELL_B].to_numpy(float)
        C = p[CELL_C].to_numpy(float)
        D = p[CELL_D].to_numpy(float)

        out = p[pair_keys].copy()
        out[f"{metric_name}_A"] = A
        out[f"{metric_name}_B"] = B
        out[f"{metric_name}_C"] = C
        out[f"{metric_name}_D"] = D

        out[f"{metric_name}_AB"] = A - B
        out[f"{metric_name}_CD"] = C - D
        out[f"{metric_name}_AC"] = A - C
        out[f"{metric_name}_BD"] = B - D
        out[f"{metric_name}_AD"] = A - D

        # Difference-of-differences interaction, equivalent along either path.
        out[f"{metric_name}_interaction_path1"] = (A - B) - (C - D)
        out[f"{metric_name}_interaction_path2"] = (A - C) - (B - D)

        # Two exact pairwise path identities.
        out[f"{metric_name}_identity_residual_via_B"] = (
            (A - D) - ((A - B) + (B - D))
        )
        out[f"{metric_name}_identity_residual_via_C"] = (
            (A - D) - ((A - C) + (C - D))
        )
        out[f"{metric_name}_interaction_identity_residual"] = (
            out[f"{metric_name}_interaction_path1"]
            - out[f"{metric_name}_interaction_path2"]
        )

        # Symmetric/Shapley accounting of the combined A-D gap.
        # Correlation-type contribution averages the correlation switch at
        # subset and all-edge representations.
        out[f"{metric_name}_shapley_corr_type"] = 0.5 * (
            (A - B) + (C - D)
        )
        # Edge-universe contribution averages the subset->all switch under
        # Spearman and Pearson.
        out[f"{metric_name}_shapley_edge_universe"] = 0.5 * (
            (A - C) + (B - D)
        )
        out[f"{metric_name}_shapley_identity_residual"] = (
            (A - D)
            - (
                out[f"{metric_name}_shapley_corr_type"]
                + out[f"{metric_name}_shapley_edge_universe"]
            )
        )

        if pairwise is None:
            pairwise = out
        else:
            pairwise = pairwise.merge(out, on=pair_keys, how="inner", validate="one_to_one")

    require(pairwise is not None and len(pairwise) == 22, "Pairwise reconstruction failed.")

    # Pairwise algebra must be exact up to machine precision.
    residual_cols = [
        c for c in pairwise.columns
        if "identity_residual" in c
    ]
    max_abs_residual = 0.0
    for col in residual_cols:
        val = float(np.max(np.abs(pairwise[col].to_numpy(float))))
        max_abs_residual = max(max_abs_residual, val)
        require(
            val <= TOL,
            f"Algebraic identity failure in {col}: max |residual|={val:.3e}",
        )

    # Verify against the previously saved decomposition TSV.
    saved = decomp_saved.copy()
    require(len(saved) == 22, f"Expected 22 saved decomposition rows, found {len(saved)}.")
    check = pairwise.merge(saved, on=pair_keys, how="inner", validate="one_to_one")
    require(len(check) == 22, "Saved-decomposition merge did not preserve all 22 pairs.")

    max_saved_replay_diff = 0.0
    for metric_name in ["f50", "auc"]:
        mapping = DECOMP_EXPECTED[metric_name]
        reconstructed = {
            "AB": f"{metric_name}_AB",
            "CD": f"{metric_name}_CD",
            "AC": f"{metric_name}_AC",
            "BD": f"{metric_name}_BD",
            "interaction": f"{metric_name}_interaction_path1",
            "AD": f"{metric_name}_AD",
        }
        for label, saved_col in mapping.items():
            require(saved_col in check.columns, f"Saved decomposition missing: {saved_col}")
            diff = np.abs(
                check[reconstructed[label]].to_numpy(float)
                - check[saved_col].to_numpy(float)
            )
            md = float(np.max(diff))
            max_saved_replay_diff = max(max_saved_replay_diff, md)
            require(
                md <= TOL,
                f"Saved decomposition replay mismatch {metric_name}/{label}: max |Δ|={md:.3e}",
            )

    # Reporting summaries: mean + median/IQR/range together.
    summary_rows = []
    quantities = {
        "AB_corr_type_subset": "AB",
        "CD_corr_type_all": "CD",
        "AC_edge_universe_spearman": "AC",
        "BD_edge_universe_pearson": "BD",
        "interaction": "interaction_path1",
        "AD_combined": "AD",
        "shapley_corr_type": "shapley_corr_type",
        "shapley_edge_universe": "shapley_edge_universe",
    }

    for metric_name in ["f50", "auc"]:
        for label, suffix in quantities.items():
            s = summarize(pairwise[f"{metric_name}_{suffix}"].to_numpy(float))
            summary_rows.append({
                "metric": metric_name,
                "quantity": label,
                **s,
            })

    summary = pd.DataFrame(summary_rows)

    # By-target descriptive replication; still no inferential test.
    target_rows = []
    for target, gg in pairwise.groupby("target", sort=True):
        for metric_name in ["f50", "auc"]:
            for label, suffix in quantities.items():
                s = summarize(gg[f"{metric_name}_{suffix}"].to_numpy(float))
                target_rows.append({
                    "target": target,
                    "metric": metric_name,
                    "quantity": label,
                    **s,
                })
    target_summary = pd.DataFrame(target_rows)

    # Mean identities: these MUST add exactly because expectation/mean is linear.
    identity_rows = []
    for metric_name in ["f50", "auc"]:
        AB = pairwise[f"{metric_name}_AB"].to_numpy(float)
        CD = pairwise[f"{metric_name}_CD"].to_numpy(float)
        AC = pairwise[f"{metric_name}_AC"].to_numpy(float)
        BD = pairwise[f"{metric_name}_BD"].to_numpy(float)
        AD = pairwise[f"{metric_name}_AD"].to_numpy(float)
        SC = pairwise[f"{metric_name}_shapley_corr_type"].to_numpy(float)
        SE = pairwise[f"{metric_name}_shapley_edge_universe"].to_numpy(float)

        mean_AD = float(np.mean(AD))
        via_B = float(np.mean(AB) + np.mean(BD))
        via_C = float(np.mean(AC) + np.mean(CD))
        via_shapley = float(np.mean(SC) + np.mean(SE))

        identity_rows.extend([
            {
                "metric": metric_name,
                "identity": "mean_AD_equals_mean_AB_plus_mean_BD",
                "lhs": mean_AD,
                "rhs": via_B,
                "residual": mean_AD - via_B,
            },
            {
                "metric": metric_name,
                "identity": "mean_AD_equals_mean_AC_plus_mean_CD",
                "lhs": mean_AD,
                "rhs": via_C,
                "residual": mean_AD - via_C,
            },
            {
                "metric": metric_name,
                "identity": "mean_AD_equals_mean_shapley_corr_plus_edge",
                "lhs": mean_AD,
                "rhs": via_shapley,
                "residual": mean_AD - via_shapley,
            },
        ])

    identity = pd.DataFrame(identity_rows)
    require(
        float(np.max(np.abs(identity["residual"].to_numpy(float)))) <= TOL,
        "Mean linearity identity failed beyond numerical tolerance.",
    )

    # Descriptive joint-dependence diagnostics.
    # A large median non-additivity does NOT by itself prove marginal skewness;
    # dependence between pair-level components can also drive it.
    dep_rows = []
    for metric_name in ["f50", "auc"]:
        for path_label, x_suffix, y_suffix in [
            ("via_B_AB_plus_BD", "AB", "BD"),
            ("via_C_AC_plus_CD", "AC", "CD"),
        ]:
            x = pairwise[f"{metric_name}_{x_suffix}"].to_numpy(float)
            y = pairwise[f"{metric_name}_{y_suffix}"].to_numpy(float)
            pr, sr, n = finite_corr(x, y)
            dep_rows.append({
                "metric": metric_name,
                "path": path_label,
                "n_pairs": n,
                "pearson_component_correlation_descriptive": pr,
                "spearman_component_correlation_descriptive": sr,
                "p_values_computed": False,
            })
    dependence = pd.DataFrame(dep_rows)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    for p in [PAIRWISE_OUT, SUMMARY_OUT, TARGET_OUT, IDENTITY_OUT, DEPENDENCE_OUT, MASTER_OUT]:
        require(not p.exists(), f"Refusing to overwrite existing audit output: {p}")

    pairwise.to_csv(PAIRWISE_OUT, sep="\t", index=False)
    summary.to_csv(SUMMARY_OUT, sep="\t", index=False)
    target_summary.to_csv(TARGET_OUT, sep="\t", index=False)
    identity.to_csv(IDENTITY_OUT, sep="\t", index=False)
    dependence.to_csv(DEPENDENCE_OUT, sep="\t", index=False)

    master_out = {
        "script_version": SCRIPT_VERSION,
        "status": "TWO_BY_TWO_REPORTING_AUDIT_COMPLETE",
        "scientific_role": "deterministic reporting/algebra audit of already-observed frozen 2x2 results",
        "new_experiment": False,
        "new_inferential_test": False,
        "p_values_computed": False,
        "bootstrap_computed": False,
        "pairwise_algebraic_identities": "PASS",
        "saved_decomposition_exact_replay": "PASS",
        "mean_linearity_identities": "PASS",
        "max_abs_pairwise_identity_residual": max_abs_residual,
        "max_abs_saved_decomposition_replay_diff": max_saved_replay_diff,
        "cell_D_reporting_name_recommendation": (
            "all-edge Pearson cor.cor-equivalent common-benchmark statistic; "
            "validated at zero corruption against package output in historical 05j"
        ),
        "outputs": {
            "pairwise": str(PAIRWISE_OUT),
            "summary": str(SUMMARY_OUT),
            "summary_by_target": str(TARGET_OUT),
            "mean_identity_audit": str(IDENTITY_OUT),
            "component_dependence": str(DEPENDENCE_OUT),
        },
    }
    MASTER_OUT.write_text(json.dumps(master_out, indent=2) + "\n", encoding="utf-8")

    print(SEP)
    print("05j4d 2x2 REPORTING CONSISTENCY AUDIT: COMPLETE")
    print(SEP)
    print(f"Pairwise algebra identities:          PASS (max |residual|={max_abs_residual:.3e})")
    print(f"Saved decomposition exact replay:     PASS (max |Δ|={max_saved_replay_diff:.3e})")
    print("Mean linearity identities:            PASS")
    print("New p-values / bootstrap:             NO")
    print()

    print("Mean identity audit:")
    print(identity.to_string(index=False))
    print()

    print("Overall reporting summary:")
    print(
        summary[
            ["metric", "quantity", "n", "mean", "median", "q25", "q75",
             "min", "max", "skewness_descriptive",
             "negative_pairs", "zero_pairs", "positive_pairs"]
        ].to_string(index=False)
    )
    print()

    print("Descriptive component dependence:")
    print(dependence.to_string(index=False))
    print()

    print("Recommended D-cell manuscript label:")
    print("  all-edge Pearson cor.cor-equivalent common-benchmark statistic")
    print("  (historically validated at zero corruption against package output)")
    print()
    print(f"Pairwise:   {PAIRWISE_OUT}")
    print(f"Summary:    {SUMMARY_OUT}")
    print(f"By target:  {TARGET_OUT}")
    print(f"Identity:   {IDENTITY_OUT}")
    print(f"Dependence: {DEPENDENCE_OUT}")
    print(f"Master:     {MASTER_OUT}")
    print(SEP)

if __name__ == "__main__":
    main()
