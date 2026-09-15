from __future__ import annotations

import json
import math
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr


SCRIPT_VERSION = "05j4-run-scale-invariant-corruption-audit-v2-no-cli"

DATA_ROOT = Path(r"D:\paper4_tcbb_data")

CONTRACT = (
    DATA_ROOT
    / "paper4_tcbb_scale_invariant_corruption_audit_contract_v2"
    / "scale_invariant_corruption_audit_contract_v2.json"
)
INPUT_MASTER = (
    DATA_ROOT
    / "paper4_tcbb_scale_invariant_corruption_inputs_v2"
    / "scale_invariant_inputs_v2.json"
)
PARTIAL = (
    DATA_ROOT
    / "paper4_tcbb_scale_invariant_corruption_inputs_v2"
    / "scale_invariant_partial_augmented_v2.tsv"
)
FULL_NULL = (
    DATA_ROOT
    / "paper4_tcbb_scale_invariant_corruption_inputs_v2"
    / "scale_invariant_fullnull_1000_v2.tsv"
)

OUT_DIR = DATA_ROOT / "paper4_tcbb_scale_invariant_corruption_audit_v2"

ALPHA = 0.05
CAL_N = 800
HOLDOUT_N = 200
BOOTSTRAPS = 2000
BASE_SEED = 20260922

FRACTIONS = np.array([0.0, 0.10, 0.25, 0.50, 0.75, 1.0], dtype=float)
PERF_FRACTIONS = np.array([0.0, 0.10, 0.25, 0.50, 0.75], dtype=float)

METHODS = {
    "MTA": "mta_rho_edge",
    "Pearson_subset_matched": "pearson_subset_matched",
    "NetRep_cor.cor": "netrep_cor_cor",
    "WGCNA_cor.kIM": "wgcna_cor_kIM",
}
PRIMARY_COMPARATORS = [
    "NetRep_cor.cor",
    "WGCNA_cor.kIM",
]
DIAGNOSTIC_COMPARATORS = [
    "Pearson_subset_matched",
]

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


def f50(fracs: np.ndarray, rates: np.ndarray) -> float:
    y = np.asarray(rates, dtype=float)
    if not np.all(np.isfinite(y)):
        return np.nan
    if y[0] <= 0.5:
        return 0.0
    for i in range(1, len(fracs)):
        if y[i - 1] >= 0.5 and y[i] <= 0.5:
            x0, x1 = float(fracs[i - 1]), float(fracs[i])
            y0, y1 = float(y[i - 1]), float(y[i])
            if y1 == y0:
                return x1
            return x0 + (0.5 - y0) * (x1 - x0) / (y1 - y0)
    return np.nan


def auc_perf(rates: np.ndarray) -> float:
    y = np.asarray(rates[:5], dtype=float)
    if not np.all(np.isfinite(y)):
        return np.nan
    return float(np.trapezoid(y, PERF_FRACTIONS))


def boot_rng(target: str, program_id: str, b: int) -> np.random.Generator:
    ss = np.random.SeedSequence(
        [BASE_SEED, TARGET_CODE[target], module_number(program_id), int(b)]
    )
    return np.random.default_rng(ss)


def interp(lo: float, hi: float, metric: str) -> str:
    # For both f50 and AUC, negative MTA-comparator means MTA detects failure earlier.
    if np.isfinite(hi) and hi < 0:
        return "MTA_EARLIER"
    if np.isfinite(lo) and lo > 0:
        return "COMPARATOR_EARLIER"
    return "NO_CLEAR_DIFFERENCE"


def main() -> None:
    print("=" * 162)
    print("Paper 4 / TCBB - scale-invariant mapping-error audit v2")
    print("=" * 162)
    print(f"Script version: {SCRIPT_VERSION}")
    print()
    print("Frozen analysis:")
    print("  Extended full-corruption null/pair:       1,000")
    print("  Calibration / held-out null split:        800 / 200")
    print("  Equal false-preservation calibration:     <=5%")
    print("  Primary endpoint:                         f50")
    print("  Secondary endpoint:                       retention-call AUC")
    print("  Same-edge Pearson diagnostic:             YES")
    print("  Within-fraction method correlation:       YES")
    print("  Paired bootstrap:                         2,000")
    print("  General superiority claim:                NO")
    print("=" * 162)

    for p in [CONTRACT, INPUT_MASTER, PARTIAL, FULL_NULL]:
        require(p)

    contract = json.loads(CONTRACT.read_text(encoding="utf-8"))
    inp = json.loads(INPUT_MASTER.read_text(encoding="utf-8"))

    if contract.get("status") != "FROZEN_POST_05J3_BEFORE_SCALE_INVARIANT_AUDIT_V2_RESULTS":
        raise RuntimeError("05j4a v2 contract has unexpected status.")
    if inp.get("status") != "SCALE_INVARIANT_INPUTS_V2_COMPLETE":
        raise RuntimeError("05j4a input master has unexpected status.")

    partial = pd.read_csv(PARTIAL, sep="\t", low_memory=False)
    fullnull = pd.read_csv(FULL_NULL, sep="\t", low_memory=False)

    if len(partial) != 11022:
        raise RuntimeError(f"Expected 11,022 partial/head-to-head rows, found {len(partial)}")
    if len(fullnull) != 22000:
        raise RuntimeError(f"Expected 22,000 extended full-null rows, found {len(fullnull)}")

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    threshold_rows = []
    curve_rows = []
    within_rows = []
    metric_rows = []
    contrast_rows = []
    target_summary_rows = []

    pairs = list(partial.groupby(["target", "program_id"], sort=True))
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
        if len(hold) != HOLDOUT_N:
            raise RuntimeError("Holdout null size drift.")

        # Within-fraction correlations: do not pool across corruption levels.
        for comp in PRIMARY_COMPARATORS + DIAGNOSTIC_COMPARATORS:
            for fraction in [0.10, 0.25, 0.50, 0.75]:
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

        observed_metrics = {}
        threshold_boot_store = {method: np.empty(BOOTSTRAPS) for method in METHODS}

        for method, col in METHODS.items():
            threshold = threshold_from_calibration(cal[col].to_numpy(dtype=float))
            hold_fp = float(np.mean(hold[col].to_numpy(dtype=float) > threshold))

            rates = []
            for fraction in FRACTIONS[:-1]:
                gg = g.loc[g["fraction"] == fraction]
                expected = 1 if fraction == 0.0 else 100
                if len(gg) != expected:
                    raise RuntimeError(
                        f"{target} {program_id} {fraction}: rows={len(gg)}, expected={expected}"
                    )
                rate = float(np.mean(gg[col].to_numpy(dtype=float) > threshold))
                rates.append(rate)
            rates.append(hold_fp)
            rates_arr = np.asarray(rates, dtype=float)

            obs_f50 = f50(FRACTIONS, rates_arr)
            obs_auc = auc_perf(rates_arr)

            observed_metrics[method] = {
                "threshold": threshold,
                "holdout_fp": hold_fp,
                "rates": rates_arr,
                "f50": obs_f50,
                "auc": obs_auc,
            }

            for fraction, rate in zip(FRACTIONS, rates_arr):
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
                    "f50": obs_f50,
                    "retention_auc_0_to_0p75": obs_auc,
                }
            )

        # Paired bootstrap: re-estimate threshold and holdout performance.
        comp_list = PRIMARY_COMPARATORS + DIAGNOSTIC_COMPARATORS
        delta_f50 = {comp: np.full(BOOTSTRAPS, np.nan) for comp in comp_list}
        delta_auc = {comp: np.full(BOOTSTRAPS, np.nan) for comp in comp_list}

        partial_by_fraction = {
            float(f): g.loc[g["fraction"] == f]
            .sort_values("replicate_id").reset_index(drop=True)
            for f in FRACTIONS[1:-1]
        }
        zero = g.loc[g["fraction"] == 0.0]
        if len(zero) != 1:
            raise RuntimeError("Fraction-zero row-count drift.")

        for b in range(1, BOOTSTRAPS + 1):
            rng = boot_rng(target, program_id, b)

            cal_idx = rng.integers(0, CAL_N, size=CAL_N, endpoint=False)
            hold_idx = rng.integers(0, HOLDOUT_N, size=HOLDOUT_N, endpoint=False)
            part_idx = {
                f: rng.integers(0, 100, size=100, endpoint=False)
                for f in [0.10, 0.25, 0.50, 0.75]
            }

            boot = {}
            for method, col in METHODS.items():
                cal_vals = cal[col].to_numpy(dtype=float)[cal_idx]
                threshold_b = threshold_from_calibration(cal_vals)
                threshold_boot_store[method][b - 1] = threshold_b

                rates_b = [
                    float(float(zero[col].iloc[0]) > threshold_b)
                ]
                for f in [0.10, 0.25, 0.50, 0.75]:
                    vals = partial_by_fraction[f][col].to_numpy(dtype=float)
                    rates_b.append(float(np.mean(vals[part_idx[f]] > threshold_b)))

                hold_vals = hold[col].to_numpy(dtype=float)
                rates_b.append(float(np.mean(hold_vals[hold_idx] > threshold_b)))
                rates_b = np.asarray(rates_b, dtype=float)

                boot[method] = {
                    "f50": f50(FRACTIONS, rates_b),
                    "auc": auc_perf(rates_b),
                }

            for comp in comp_list:
                if np.isfinite(boot["MTA"]["f50"]) and np.isfinite(boot[comp]["f50"]):
                    delta_f50[comp][b - 1] = boot["MTA"]["f50"] - boot[comp]["f50"]
                if np.isfinite(boot["MTA"]["auc"]) and np.isfinite(boot[comp]["auc"]):
                    delta_auc[comp][b - 1] = boot["MTA"]["auc"] - boot[comp]["auc"]

        for method in METHODS:
            tb = threshold_boot_store[method]
            threshold_rows.append(
                {
                    "target": target,
                    "program_id": program_id,
                    "method": method,
                    "threshold_observed": observed_metrics[method]["threshold"],
                    "threshold_boot_q025": float(np.quantile(tb, 0.025)),
                    "threshold_boot_q975": float(np.quantile(tb, 0.975)),
                    "holdout_false_preservation_rate": observed_metrics[method]["holdout_fp"],
                    "calibration_n": CAL_N,
                    "holdout_n": HOLDOUT_N,
                }
            )

        for comp in comp_list:
            df50 = delta_f50[comp]
            da = delta_auc[comp]
            fmask = np.isfinite(df50)
            amask = np.isfinite(da)

            f_lo = float(np.quantile(df50[fmask], 0.025)) if fmask.any() else np.nan
            f_hi = float(np.quantile(df50[fmask], 0.975)) if fmask.any() else np.nan
            a_lo = float(np.quantile(da[amask], 0.025)) if amask.any() else np.nan
            a_hi = float(np.quantile(da[amask], 0.975)) if amask.any() else np.nan

            obs_df50 = (
                observed_metrics["MTA"]["f50"] - observed_metrics[comp]["f50"]
                if np.isfinite(observed_metrics["MTA"]["f50"])
                and np.isfinite(observed_metrics[comp]["f50"])
                else np.nan
            )
            obs_da = observed_metrics["MTA"]["auc"] - observed_metrics[comp]["auc"]

            contrast_rows.append(
                {
                    "target": target,
                    "program_id": program_id,
                    "comparator": comp,
                    "delta_f50_mta_minus_comparator": obs_df50,
                    "delta_f50_boot_q025": f_lo,
                    "delta_f50_boot_q975": f_hi,
                    "delta_f50_boot_valid": int(fmask.sum()),
                    "f50_interpretation": (
                        "MTA_EARLIER" if np.isfinite(f_hi) and f_hi < 0
                        else "COMPARATOR_EARLIER" if np.isfinite(f_lo) and f_lo > 0
                        else "NO_CLEAR_DIFFERENCE"
                    ),
                    "delta_auc_mta_minus_comparator": obs_da,
                    "delta_auc_boot_q025": a_lo,
                    "delta_auc_boot_q975": a_hi,
                    "delta_auc_boot_valid": int(amask.sum()),
                    "auc_interpretation": (
                        "MTA_EARLIER" if np.isfinite(a_hi) and a_hi < 0
                        else "COMPARATOR_EARLIER" if np.isfinite(a_lo) and a_lo > 0
                        else "NO_CLEAR_DIFFERENCE"
                    ),
                }
            )

            print(
                f"  {comp}: Δf50={obs_df50:+.4f} "
                f"[{f_lo:+.4f},{f_hi:+.4f}] | "
                f"ΔAUC={obs_da:+.4f} [{a_lo:+.4f},{a_hi:+.4f}]"
            )

    thresholds = pd.DataFrame(threshold_rows)
    curves = pd.DataFrame(curve_rows)
    within = pd.DataFrame(within_rows)
    metrics = pd.DataFrame(metric_rows)
    contrasts = pd.DataFrame(contrast_rows)

    # Distribution summaries, not inferential counts.
    summary_rows = []
    for (target, comp), gg in contrasts.groupby(["target", "comparator"], sort=True):
        fvals = pd.to_numeric(gg["delta_f50_mta_minus_comparator"], errors="coerce")
        avals = pd.to_numeric(gg["delta_auc_mta_minus_comparator"], errors="coerce")
        summary_rows.append(
            {
                "target": target,
                "comparator": comp,
                "n_pairs": len(gg),
                "median_delta_f50": float(np.nanmedian(fvals)),
                "q25_delta_f50": float(np.nanquantile(fvals, 0.25)),
                "q75_delta_f50": float(np.nanquantile(fvals, 0.75)),
                "median_delta_auc": float(np.nanmedian(avals)),
                "q25_delta_auc": float(np.nanquantile(avals, 0.25)),
                "q75_delta_auc": float(np.nanquantile(avals, 0.75)),
            }
        )
    summary = pd.DataFrame(summary_rows)

    within_summary = (
        within.groupby(["comparator", "fraction"])["spearman_within_fraction"]
        .agg(["median", "min", "max"])
        .reset_index()
    )

    thresholds_path = OUT_DIR / "scale_invariant_thresholds_v2.tsv"
    curves_path = OUT_DIR / "scale_invariant_call_curves_v2.tsv"
    within_path = OUT_DIR / "within_fraction_method_correlations_v2.tsv"
    within_summary_path = OUT_DIR / "within_fraction_method_correlations_summary_v2.tsv"
    metrics_path = OUT_DIR / "scale_invariant_operating_metrics_v2.tsv"
    contrasts_path = OUT_DIR / "scale_invariant_paired_contrasts_v2.tsv"
    summary_path = OUT_DIR / "scale_invariant_effect_distribution_by_target_v2.tsv"

    thresholds.to_csv(thresholds_path, sep="\t", index=False)
    curves.to_csv(curves_path, sep="\t", index=False)
    within.to_csv(within_path, sep="\t", index=False)
    within_summary.to_csv(within_summary_path, sep="\t", index=False)
    metrics.to_csv(metrics_path, sep="\t", index=False)
    contrasts.to_csv(contrasts_path, sep="\t", index=False)
    summary.to_csv(summary_path, sep="\t", index=False)

    master = {
        "script_version": SCRIPT_VERSION,
        "status": "SCALE_INVARIANT_CORRUPTION_AUDIT_V2_COMPLETE",
        "post_hoc_relative_to_05j3": True,
        "full_null_valid_per_pair": 1000,
        "calibration_n": CAL_N,
        "holdout_n": HOLDOUT_N,
        "primary_endpoint": "f50",
        "secondary_endpoint": "retention_call_auc_0_to_0.75",
        "general_superiority_claim_permitted": False,
        "domain_specific_mapping_error_claim_only": True,
        "thresholds_file": str(thresholds_path),
        "curves_file": str(curves_path),
        "within_fraction_correlations_file": str(within_path),
        "within_fraction_summary_file": str(within_summary_path),
        "metrics_file": str(metrics_path),
        "paired_contrasts_file": str(contrasts_path),
        "effect_distribution_file": str(summary_path),
    }
    master_path = OUT_DIR / "scale_invariant_corruption_audit_v2.json"
    master_path.write_text(json.dumps(master, indent=2), encoding="utf-8")

    print("\n" + "=" * 162)
    print("05j4 SCALE-INVARIANT CORRUPTION AUDIT v2: COMPLETE")
    print("=" * 162)
    print("\nEffect distributions by target:")
    print(summary.to_string(index=False))
    print("\nWithin-fraction rank-correlation summaries:")
    print(within_summary.to_string(index=False))
    print()
    print(f"Contrasts: {contrasts_path}")
    print(f"Thresholds: {thresholds_path}")
    print(f"Curves: {curves_path}")
    print(f"Within-fraction: {within_path}")
    print(f"Master: {master_path}")
    print("=" * 162)


if __name__ == "__main__":
    main()
