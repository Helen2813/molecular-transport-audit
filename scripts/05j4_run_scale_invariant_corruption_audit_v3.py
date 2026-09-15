from __future__ import annotations

import json
import math
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr


SCRIPT_VERSION = "05j4-run-scale-invariant-corruption-audit-v3-no-cli"

DATA_ROOT = Path(r"D:\paper4_tcbb_data")

CONTRACT = (
    DATA_ROOT
    / "paper4_tcbb_scale_invariant_corruption_audit_contract_v3"
    / "scale_invariant_corruption_audit_contract_v3.json"
)
INPUT_MASTER = (
    DATA_ROOT
    / "paper4_tcbb_scale_invariant_corruption_inputs_v3"
    / "scale_invariant_inputs_v3.json"
)
PERFORMANCE = (
    DATA_ROOT
    / "paper4_tcbb_scale_invariant_corruption_inputs_v3"
    / "scale_invariant_performance_grid_v3.tsv"
)
FULL_NULL = (
    DATA_ROOT
    / "paper4_tcbb_scale_invariant_corruption_inputs_v3"
    / "scale_invariant_fullnull_1000_v3.tsv"
)

OUT_DIR = DATA_ROOT / "paper4_tcbb_scale_invariant_corruption_audit_v3"

ALPHA = 0.05
CAL_N = 800
HOLDOUT_N = 200
PAIR_BOOTSTRAPS = 2000
CLUSTER_BOOTSTRAPS = 2000
BASE_SEED = 20260922

FRACTIONS = np.array([0.0, 0.05, 0.10, 0.15, 0.25, 0.50, 0.75, 1.0], dtype=float)
PERF_FRACTIONS = np.array([0.0, 0.05, 0.10, 0.15, 0.25, 0.50, 0.75], dtype=float)

METHODS = {
    "MTA": "mta_rho_edge",
    "Pearson_subset_matched": "pearson_subset_matched",
    "NetRep_cor.cor": "netrep_cor_cor",
    "WGCNA_cor.kIM": "wgcna_cor_kIM",
}
PRIMARY_COMPARATORS = ["NetRep_cor.cor", "WGCNA_cor.kIM"]
DIAGNOSTIC_COMPARATORS = ["Pearson_subset_matched"]

TARGET_CODE = {
    "SCANB_GSE96058": 1,
    "METABRIC": 2,
}


def require(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Required file does not exist: {path}")


def module_number(program_id: str) -> int:
    m = re.search(r"M(\d+)$", str(program_id))
    if not m:
        raise RuntimeError(f"Cannot parse module number from {program_id}")
    return int(m.group(1))


def threshold_from_calibration(x: np.ndarray) -> float:
    vals = np.asarray(x, dtype=float)
    vals = vals[np.isfinite(vals)]
    if len(vals) != CAL_N:
        raise RuntimeError(f"Expected {CAL_N} calibration scores, got {len(vals)}")
    k = int(math.floor(ALPHA * len(vals)))
    xs = np.sort(vals)
    return float(xs[len(xs) - k - 1])


def f50_with_status(fracs: np.ndarray, rates: np.ndarray) -> tuple[float, str]:
    y = np.asarray(rates, dtype=float)
    if not np.all(np.isfinite(y)):
        return np.nan, "NOT_ESTIMABLE"
    if y[0] <= 0.5:
        return 0.0, "FAILED_AT_UNCORRUPTED_MAPPING"

    for i in range(1, len(fracs)):
        if y[i - 1] > 0.5 and y[i] <= 0.5:
            x0, x1 = float(fracs[i - 1]), float(fracs[i])
            y0, y1 = float(y[i - 1]), float(y[i])
            if y1 == y0:
                return x1, "FINITE"
            val = x0 + (0.5 - y0) * (x1 - x0) / (y1 - y0)
            return float(val), "FINITE"

    return np.nan, "RIGHT_CENSORED_GT_1"


def auc_perf(rates: np.ndarray) -> float:
    y = np.asarray(rates[: len(PERF_FRACTIONS)], dtype=float)
    if not np.all(np.isfinite(y)):
        return np.nan
    return float(np.trapezoid(y, PERF_FRACTIONS))


def pair_boot_rng(target: str, program_id: str, b: int) -> np.random.Generator:
    ss = np.random.SeedSequence(
        [BASE_SEED, TARGET_CODE[target], module_number(program_id), int(b)]
    )
    return np.random.default_rng(ss)


def cluster_rng(b: int) -> np.random.Generator:
    ss = np.random.SeedSequence([BASE_SEED, 999, int(b)])
    return np.random.default_rng(ss)


def ci_interpret(lo: float, hi: float) -> str:
    if np.isfinite(hi) and hi < 0:
        return "MTA_EARLIER"
    if np.isfinite(lo) and lo > 0:
        return "COMPARATOR_EARLIER"
    return "NO_CLEAR_DIFFERENCE"


def quantile_or_nan(x: np.ndarray, p: float) -> float:
    z = np.asarray(x, dtype=float)
    z = z[np.isfinite(z)]
    return float(np.quantile(z, p)) if len(z) else np.nan


def main() -> None:
    print("=" * 164)
    print("Paper 4 / TCBB - scale-invariant mapping-error audit v3")
    print("=" * 164)
    print(f"Script version: {SCRIPT_VERSION}")
    print()
    print("Frozen analysis:")
    print("  Corruption grid:                         0,.05,.10,.15,.25,.50,.75,1")
    print("  Full-null calibration / holdout:         800 / 200")
    print("  Equal false-preservation target:         <=5%")
    print("  Headline resolution metric:              f50 with censoring")
    print("  Primary inferential contrast:            retention-call AUC")
    print("  Same-edge Pearson diagnostic:            YES")
    print("  Within-fraction correlations:            YES")
    print("  Module-cluster aggregate bootstrap:      2,000")
    print("  General superiority claim:               NO")
    print("=" * 164)

    for p in [CONTRACT, INPUT_MASTER, PERFORMANCE, FULL_NULL]:
        require(p)

    contract = json.loads(CONTRACT.read_text(encoding="utf-8"))
    inp = json.loads(INPUT_MASTER.read_text(encoding="utf-8"))

    if contract.get("status") != "FROZEN_POST_05J3_BEFORE_SCALE_INVARIANT_AUDIT_V3_RESULTS":
        raise RuntimeError("05j4b v3 contract has unexpected status.")
    if inp.get("status") != "SCALE_INVARIANT_INPUTS_V3_COMPLETE":
        raise RuntimeError("05j4b v3 input master has unexpected status.")

    perf = pd.read_csv(PERFORMANCE, sep="\t", low_memory=False)
    fullnull = pd.read_csv(FULL_NULL, sep="\t", low_memory=False)

    if len(perf) != 22 * 601:
        raise RuntimeError(f"Performance-grid row count={len(perf)}, expected {22*601}")
    if len(fullnull) != 22000:
        raise RuntimeError(f"Full-null row count={len(fullnull)}, expected 22000")

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    threshold_rows = []
    curve_rows = []
    within_rows = []
    metric_rows = []
    contrast_rows = []

    # Store pair-level bootstrap AUC contrasts for cluster bootstrap.
    pair_boot_auc: dict[tuple[str, str, str], np.ndarray] = {}

    pairs = list(perf.groupby(["target", "program_id"], sort=True))
    if len(pairs) != 22:
        raise RuntimeError(f"Expected 22 pairs, found {len(pairs)}")

    for pi, ((target, program_id), g) in enumerate(pairs, start=1):
        print(f"\n[{pi:02d}/22] {target} {program_id}")

        g = g.sort_values(["fraction", "replicate_id"]).reset_index(drop=True)
        n = fullnull.loc[
            (fullnull["target"] == target)
            & (fullnull["program_id"] == program_id)
        ].sort_values("valid_null_id").reset_index(drop=True)

        if len(n) != 1000:
            raise RuntimeError(f"{target} {program_id}: extended null n={len(n)}")

        cal = n.iloc[:CAL_N].copy()
        hold = n.iloc[CAL_N:].copy()

        # Within-fraction correlations only.
        for comp in PRIMARY_COMPARATORS + DIAGNOSTIC_COMPARATORS:
            for fraction in [0.05, 0.10, 0.15, 0.25, 0.50, 0.75]:
                gg = g.loc[g["fraction"] == fraction]
                if len(gg) != 100:
                    raise RuntimeError("Partial fraction row-count drift.")
                rho = float(
                    spearmanr(
                        gg[METHODS["MTA"]].to_numpy(dtype=float),
                        gg[METHODS[comp]].to_numpy(dtype=float),
                    ).statistic
                )
                within_rows.append(
                    {
                        "target": target,
                        "program_id": program_id,
                        "comparator": comp,
                        "fraction": fraction,
                        "spearman_within_fraction": rho,
                        "n": 100,
                    }
                )

            rho1 = float(
                spearmanr(
                    hold[METHODS["MTA"]].to_numpy(dtype=float),
                    hold[METHODS[comp]].to_numpy(dtype=float),
                ).statistic
            )
            within_rows.append(
                {
                    "target": target,
                    "program_id": program_id,
                    "comparator": comp,
                    "fraction": 1.0,
                    "spearman_within_fraction": rho1,
                    "n": HOLDOUT_N,
                }
            )

        observed = {}
        threshold_boot_store = {
            method: np.empty(PAIR_BOOTSTRAPS, dtype=float)
            for method in METHODS
        }

        for method, col in METHODS.items():
            threshold = threshold_from_calibration(cal[col].to_numpy(dtype=float))
            hold_fp = float(np.mean(hold[col].to_numpy(dtype=float) > threshold))

            rates = []
            for fraction in FRACTIONS[:-1]:
                gg = g.loc[g["fraction"] == fraction]
                expected = 1 if fraction == 0.0 else 100
                if len(gg) != expected:
                    raise RuntimeError(
                        f"{target} {program_id} fraction {fraction}: "
                        f"{len(gg)} rows, expected {expected}"
                    )
                rates.append(float(np.mean(gg[col].to_numpy(dtype=float) > threshold)))
            rates.append(hold_fp)
            rates = np.asarray(rates, dtype=float)

            f50_value, f50_status = f50_with_status(FRACTIONS, rates)
            auc_value = auc_perf(rates)

            observed[method] = {
                "threshold": threshold,
                "holdout_fp": hold_fp,
                "rates": rates,
                "f50": f50_value,
                "f50_status": f50_status,
                "auc": auc_value,
            }

            for fraction, rate in zip(FRACTIONS, rates):
                curve_rows.append(
                    {
                        "target": target,
                        "program_id": program_id,
                        "method": method,
                        "fraction": float(fraction),
                        "preservation_call_rate": float(rate),
                        "failure_detection_rate": float(1.0 - rate),
                        "threshold": threshold,
                    }
                )

            metric_rows.append(
                {
                    "target": target,
                    "program_id": program_id,
                    "method": method,
                    "threshold": threshold,
                    "holdout_false_preservation_rate": hold_fp,
                    "f50": f50_value,
                    "f50_status": f50_status,
                    "retention_auc_0_to_0p75": auc_value,
                }
            )

        # Pair-level paired bootstrap.
        comp_list = PRIMARY_COMPARATORS + DIAGNOSTIC_COMPARATORS
        boot_delta_auc = {
            comp: np.full(PAIR_BOOTSTRAPS, np.nan, dtype=float)
            for comp in comp_list
        }
        boot_delta_f50 = {
            comp: np.full(PAIR_BOOTSTRAPS, np.nan, dtype=float)
            for comp in comp_list
        }

        zero = g.loc[g["fraction"] == 0.0]
        partial_by_fraction = {
            float(f): g.loc[g["fraction"] == f]
            .sort_values("replicate_id").reset_index(drop=True)
            for f in [0.05, 0.10, 0.15, 0.25, 0.50, 0.75]
        }

        for b in range(1, PAIR_BOOTSTRAPS + 1):
            rng = pair_boot_rng(target, program_id, b)

            cal_idx = rng.integers(0, CAL_N, size=CAL_N, endpoint=False)
            hold_idx = rng.integers(0, HOLDOUT_N, size=HOLDOUT_N, endpoint=False)
            part_idx = {
                f: rng.integers(0, 100, size=100, endpoint=False)
                for f in [0.05, 0.10, 0.15, 0.25, 0.50, 0.75]
            }

            boot = {}
            for method, col in METHODS.items():
                threshold_b = threshold_from_calibration(
                    cal[col].to_numpy(dtype=float)[cal_idx]
                )
                threshold_boot_store[method][b - 1] = threshold_b

                rates_b = [
                    float(float(zero[col].iloc[0]) > threshold_b)
                ]
                for f in [0.05, 0.10, 0.15, 0.25, 0.50, 0.75]:
                    vals = partial_by_fraction[f][col].to_numpy(dtype=float)
                    rates_b.append(float(np.mean(vals[part_idx[f]] > threshold_b)))
                hold_vals = hold[col].to_numpy(dtype=float)
                rates_b.append(float(np.mean(hold_vals[hold_idx] > threshold_b)))
                rates_b = np.asarray(rates_b, dtype=float)

                fv, fs = f50_with_status(FRACTIONS, rates_b)
                boot[method] = {
                    "f50": fv,
                    "f50_status": fs,
                    "auc": auc_perf(rates_b),
                }

            for comp in comp_list:
                boot_delta_auc[comp][b - 1] = boot["MTA"]["auc"] - boot[comp]["auc"]
                if np.isfinite(boot["MTA"]["f50"]) and np.isfinite(boot[comp]["f50"]):
                    boot_delta_f50[comp][b - 1] = boot["MTA"]["f50"] - boot[comp]["f50"]

        for method in METHODS:
            tb = threshold_boot_store[method]
            threshold_rows.append(
                {
                    "target": target,
                    "program_id": program_id,
                    "method": method,
                    "threshold_observed": observed[method]["threshold"],
                    "threshold_boot_q025": float(np.quantile(tb, 0.025)),
                    "threshold_boot_q975": float(np.quantile(tb, 0.975)),
                    "holdout_false_preservation_rate": observed[method]["holdout_fp"],
                    "calibration_n": CAL_N,
                    "holdout_n": HOLDOUT_N,
                }
            )

        for comp in comp_list:
            da = boot_delta_auc[comp]
            df50 = boot_delta_f50[comp]

            a_lo = quantile_or_nan(da, 0.025)
            a_hi = quantile_or_nan(da, 0.975)
            f_lo = quantile_or_nan(df50, 0.025)
            f_hi = quantile_or_nan(df50, 0.975)

            obs_da = observed["MTA"]["auc"] - observed[comp]["auc"]
            obs_df50 = (
                observed["MTA"]["f50"] - observed[comp]["f50"]
                if np.isfinite(observed["MTA"]["f50"])
                and np.isfinite(observed[comp]["f50"])
                else np.nan
            )

            pair_boot_auc[(target, program_id, comp)] = da.copy()

            contrast_rows.append(
                {
                    "target": target,
                    "program_id": program_id,
                    "comparator": comp,
                    "mta_f50": observed["MTA"]["f50"],
                    "mta_f50_status": observed["MTA"]["f50_status"],
                    "comparator_f50": observed[comp]["f50"],
                    "comparator_f50_status": observed[comp]["f50_status"],
                    "delta_f50_mta_minus_comparator": obs_df50,
                    "delta_f50_boot_q025": f_lo,
                    "delta_f50_boot_q975": f_hi,
                    "delta_f50_boot_valid": int(np.isfinite(df50).sum()),
                    "f50_interpretation": (
                        ci_interpret(f_lo, f_hi)
                        if np.isfinite(f_lo) and np.isfinite(f_hi)
                        else "CENSORED_OR_NOT_ESTIMABLE"
                    ),
                    "delta_auc_mta_minus_comparator": obs_da,
                    "delta_auc_boot_q025": a_lo,
                    "delta_auc_boot_q975": a_hi,
                    "delta_auc_boot_valid": int(np.isfinite(da).sum()),
                    "auc_interpretation": ci_interpret(a_lo, a_hi),
                }
            )

            print(
                f"  {comp}: Δf50={obs_df50 if np.isfinite(obs_df50) else float('nan'):+.4f} "
                f"[{f_lo:+.4f},{f_hi:+.4f}] | "
                f"ΔAUC={obs_da:+.4f} [{a_lo:+.4f},{a_hi:+.4f}]"
            )

    thresholds = pd.DataFrame(threshold_rows)
    curves = pd.DataFrame(curve_rows)
    within = pd.DataFrame(within_rows)
    metrics = pd.DataFrame(metric_rows)
    contrasts = pd.DataFrame(contrast_rows)

    # By-target descriptive distributions; no independence claim.
    summary_rows = []
    for (target, comp), gg in contrasts.groupby(["target", "comparator"], sort=True):
        fvals = pd.to_numeric(gg["delta_f50_mta_minus_comparator"], errors="coerce")
        avals = pd.to_numeric(gg["delta_auc_mta_minus_comparator"], errors="coerce")
        summary_rows.append(
            {
                "target": target,
                "comparator": comp,
                "n_pairs": int(len(gg)),
                "finite_delta_f50_pairs": int(np.isfinite(fvals).sum()),
                "median_delta_f50": float(np.nanmedian(fvals)) if np.isfinite(fvals).any() else np.nan,
                "q25_delta_f50": quantile_or_nan(fvals.to_numpy(dtype=float), 0.25),
                "q75_delta_f50": quantile_or_nan(fvals.to_numpy(dtype=float), 0.75),
                "median_delta_auc": float(np.nanmedian(avals)),
                "q25_delta_auc": quantile_or_nan(avals.to_numpy(dtype=float), 0.25),
                "q75_delta_auc": quantile_or_nan(avals.to_numpy(dtype=float), 0.75),
            }
        )
    summary = pd.DataFrame(summary_rows)

    within_summary = (
        within.groupby(["comparator", "fraction"])["spearman_within_fraction"]
        .agg(["median", "min", "max"])
        .reset_index()
    )

    # Aggregate cluster bootstrap over source modules, carrying both target-cohort pairs.
    cluster_rows = []
    module_ids = sorted(contrasts["program_id"].unique())
    if len(module_ids) != 12:
        raise RuntimeError(f"Expected 12 source-module clusters, found {len(module_ids)}")

    for comp in PRIMARY_COMPARATORS + DIAGNOSTIC_COMPARATORS:
        obs_vals = contrasts.loc[
            contrasts["comparator"] == comp,
            "delta_auc_mta_minus_comparator",
        ].to_numpy(dtype=float)
        observed_median = float(np.median(obs_vals))

        boot_med = np.full(CLUSTER_BOOTSTRAPS, np.nan, dtype=float)

        for b in range(1, CLUSTER_BOOTSTRAPS + 1):
            rng = cluster_rng(b)
            sampled_modules = rng.choice(module_ids, size=len(module_ids), replace=True)

            vals = []
            for module_id in sampled_modules:
                pair_rows = contrasts.loc[
                    (contrasts["program_id"] == module_id)
                    & (contrasts["comparator"] == comp)
                ]
                for pr in pair_rows.to_dict(orient="records"):
                    key = (str(pr["target"]), str(module_id), comp)
                    arr = pair_boot_auc[key]
                    finite = arr[np.isfinite(arr)]
                    if len(finite) == 0:
                        continue
                    vals.append(float(finite[int(rng.integers(0, len(finite)))]))

            if vals:
                boot_med[b - 1] = float(np.median(vals))

        finite = boot_med[np.isfinite(boot_med)]
        cluster_rows.append(
            {
                "comparator": comp,
                "observed_median_delta_auc": observed_median,
                "cluster_boot_median_q025": float(np.quantile(finite, 0.025)),
                "cluster_boot_median_q975": float(np.quantile(finite, 0.975)),
                "cluster_boot_valid": int(len(finite)),
                "interpretation": ci_interpret(
                    float(np.quantile(finite, 0.025)),
                    float(np.quantile(finite, 0.975)),
                ),
            }
        )

    cluster = pd.DataFrame(cluster_rows)

    thresholds_path = OUT_DIR / "scale_invariant_thresholds_v3.tsv"
    curves_path = OUT_DIR / "scale_invariant_call_curves_v3.tsv"
    within_path = OUT_DIR / "within_fraction_method_correlations_v3.tsv"
    within_summary_path = OUT_DIR / "within_fraction_method_correlations_summary_v3.tsv"
    metrics_path = OUT_DIR / "scale_invariant_operating_metrics_v3.tsv"
    contrasts_path = OUT_DIR / "scale_invariant_paired_contrasts_v3.tsv"
    summary_path = OUT_DIR / "scale_invariant_effect_distribution_by_target_v3.tsv"
    cluster_path = OUT_DIR / "module_cluster_bootstrap_auc_v3.tsv"

    thresholds.to_csv(thresholds_path, sep="\t", index=False)
    curves.to_csv(curves_path, sep="\t", index=False)
    within.to_csv(within_path, sep="\t", index=False)
    within_summary.to_csv(within_summary_path, sep="\t", index=False)
    metrics.to_csv(metrics_path, sep="\t", index=False)
    contrasts.to_csv(contrasts_path, sep="\t", index=False)
    summary.to_csv(summary_path, sep="\t", index=False)
    cluster.to_csv(cluster_path, sep="\t", index=False)

    master = {
        "script_version": SCRIPT_VERSION,
        "status": "SCALE_INVARIANT_CORRUPTION_AUDIT_V3_COMPLETE",
        "post_hoc_relative_to_05j3": True,
        "full_null_valid_per_pair": 1000,
        "calibration_n": CAL_N,
        "holdout_n": HOLDOUT_N,
        "fractions": FRACTIONS.tolist(),
        "headline_metric": "f50_with_censoring",
        "primary_inferential_contrast": "retention_call_auc_0_to_0.75",
        "same_edge_pearson_diagnostic": True,
        "all_edge_spearman_interaction_diagnostic": False,
        "within_fraction_correlations": True,
        "module_cluster_bootstrap": True,
        "general_superiority_claim_permitted": False,
        "domain_specific_mapping_error_claim_only": True,
        "thresholds_file": str(thresholds_path),
        "curves_file": str(curves_path),
        "within_fraction_correlations_file": str(within_path),
        "within_fraction_summary_file": str(within_summary_path),
        "metrics_file": str(metrics_path),
        "paired_contrasts_file": str(contrasts_path),
        "effect_distribution_file": str(summary_path),
        "module_cluster_bootstrap_file": str(cluster_path),
    }
    master_path = OUT_DIR / "scale_invariant_corruption_audit_v3.json"
    master_path.write_text(json.dumps(master, indent=2), encoding="utf-8")

    print("\n" + "=" * 164)
    print("05j4 SCALE-INVARIANT CORRUPTION AUDIT v3: COMPLETE")
    print("=" * 164)
    print("\nBy-target effect distributions:")
    print(summary.to_string(index=False))
    print("\nWithin-fraction rank-correlation summaries:")
    print(within_summary.to_string(index=False))
    print("\nModule-cluster aggregate AUC bootstrap:")
    print(cluster.to_string(index=False))
    print()
    print(f"Contrasts:      {contrasts_path}")
    print(f"Thresholds:     {thresholds_path}")
    print(f"Curves:         {curves_path}")
    print(f"Within-fraction:{within_path}")
    print(f"Cluster:        {cluster_path}")
    print(f"Master:         {master_path}")
    print("=" * 164)


if __name__ == "__main__":
    main()
