from __future__ import annotations

import json
import math
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr


SCRIPT_VERSION = "05j4-run-scale-invariant-corruption-audit-v1-no-cli"

DATA_ROOT = Path(r"D:\paper4_tcbb_data")

CONTRACT = (
    DATA_ROOT
    / "paper4_tcbb_scale_invariant_corruption_audit_contract_v1"
    / "scale_invariant_corruption_audit_contract_v1.json"
)
REPLICATES = (
    DATA_ROOT
    / "paper4_tcbb_structural_corruption_headtohead_v1"
    / "structural_headtohead_replicates_v1.tsv"
)
MASTER = (
    DATA_ROOT
    / "paper4_tcbb_structural_corruption_headtohead_v1"
    / "structural_corruption_headtohead_v1.json"
)

OUT_DIR = DATA_ROOT / "paper4_tcbb_scale_invariant_corruption_audit_v1"

ALPHA = 0.05
FRACTIONS = np.array([0.0, 0.10, 0.25, 0.50, 0.75, 1.0], dtype=float)
PERFORMANCE_FRACTIONS = np.array([0.0, 0.10, 0.25, 0.50, 0.75], dtype=float)
BOOTSTRAPS = 2000
BASE_SEED = 20260921

METHODS = {
    "MTA": "mta_rho_edge_saved",
    "NetRep_cor.cor": "netrep_cor_cor",
    "WGCNA.cor.cor": "wgcna_cor_cor",
    "WGCNA.cor.kIM": "wgcna_cor_kIM",
}

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


def threshold_from_null(values: np.ndarray, alpha: float = ALPHA) -> float:
    x = np.asarray(values, dtype=float)
    x = x[np.isfinite(x)]
    if len(x) != 100:
        raise RuntimeError(f"Expected exactly 100 full-corruption null values, got {len(x)}")
    k = int(math.floor(alpha * len(x)))
    if k < 1:
        return float(np.max(x))
    xs = np.sort(x)
    idx = len(xs) - k - 1
    return float(xs[idx])


def retained(scores: np.ndarray, threshold: float) -> np.ndarray:
    return np.asarray(scores, dtype=float) > float(threshold)


def retention_auc(call_rates: np.ndarray) -> float:
    y = np.asarray(call_rates, dtype=float)
    if not np.all(np.isfinite(y)):
        return np.nan
    return float(np.trapezoid(y, PERFORMANCE_FRACTIONS))


def rng_for_boot(target: str, program_id: str, boot_id: int) -> np.random.Generator:
    ss = np.random.SeedSequence(
        [
            BASE_SEED,
            TARGET_CODE[target],
            module_number(program_id),
            int(boot_id),
        ]
    )
    return np.random.default_rng(ss)


def pair_interpretation(lo: float, hi: float) -> str:
    if np.isfinite(hi) and hi < 0:
        return "MTA_MORE_RESPONSIVE"
    if np.isfinite(lo) and lo > 0:
        return "COMPARATOR_MORE_RESPONSIVE"
    return "NO_CLEAR_DIFFERENCE"


def main() -> None:
    print("=" * 158)
    print("Paper 4 / TCBB - scale-invariant corruption-detection audit")
    print("=" * 158)
    print(f"Script version: {SCRIPT_VERSION}")
    print()
    print("Frozen analysis:")
    print("  Full-corruption known-null calibration:          YES")
    print("  Target false-preservation rate:                  <=5%")
    print("  Strict-monotone-rescaling invariant calls:       YES")
    print("  Fraction 1.0 used for calibration, not AUC:      YES")
    print("  Performance AUC fractions:                       0,.10,.25,.50,.75")
    print("  Paired threshold+curve bootstrap:                2,000")
    print("  General superiority claim:                       NO")
    print("=" * 158)

    for p in [CONTRACT, REPLICATES, MASTER]:
        require(p)

    contract = json.loads(CONTRACT.read_text(encoding="utf-8"))
    master = json.loads(MASTER.read_text(encoding="utf-8"))

    if contract.get("status") != "FROZEN_POST_05J3_BEFORE_SCALE_INVARIANT_AUDIT_RESULTS":
        raise RuntimeError("05j4 contract has unexpected status.")
    if master.get("status") != "STRUCTURAL_CORRUPTION_HEADTOHEAD_COMPLETE":
        raise RuntimeError("05j head-to-head master has unexpected status.")

    df = pd.read_csv(REPLICATES, sep="\t", low_memory=False)
    if len(df) != 11022:
        raise RuntimeError(f"Expected 11,022 head-to-head rows, found {len(df)}")

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    threshold_rows = []
    curve_rows = []
    monotone_rows = []
    contrast_rows = []

    pair_groups = list(df.groupby(["target", "program_id"], sort=True))
    if len(pair_groups) != 22:
        raise RuntimeError(f"Expected 22 target/program pairs, found {len(pair_groups)}")

    for pair_idx, ((target, program_id), g) in enumerate(pair_groups, start=1):
        g = g.sort_values(["fraction", "replicate_id"]).reset_index(drop=True)
        print(f"\n[{pair_idx:02d}/22] {target} {program_id}")

        zero = g.loc[g["fraction"] == 0.0]
        full = g.loc[g["fraction"] == 1.0].sort_values("replicate_id").reset_index(drop=True)
        if len(zero) != 1:
            raise RuntimeError(f"{target} {program_id}: fraction-0 row count={len(zero)}")
        if len(full) != 100:
            raise RuntimeError(f"{target} {program_id}: full-corruption row count={len(full)}")

        per_method_auc: dict[str, float] = {}

        for comp_method in ["NetRep.cor.cor", "WGCNA.cor.kIM"]:
            rho = float(
                spearmanr(
                    pd.to_numeric(g[METHODS["MTA"]], errors="coerce"),
                    pd.to_numeric(g[METHODS[comp_method]], errors="coerce"),
                ).statistic
            )
            monotone_rows.append(
                {
                    "target": target,
                    "program_id": program_id,
                    "comparator": comp_method,
                    "spearman_all_501_rows": rho,
                }
            )

        for method, col in METHODS.items():
            full_scores = pd.to_numeric(full[col], errors="coerce").to_numpy(dtype=float)
            threshold = threshold_from_null(full_scores, ALPHA)
            full_fp = float(np.mean(retained(full_scores, threshold)))

            threshold_rows.append(
                {
                    "target": target,
                    "program_id": program_id,
                    "method": method,
                    "threshold": threshold,
                    "full_corruption_false_preservation_rate": full_fp,
                    "null_n": 100,
                    "alpha_target": ALPHA,
                }
            )

            rates = []
            for fraction in FRACTIONS:
                s = g.loc[g["fraction"] == fraction]
                vals = pd.to_numeric(s[col], errors="coerce").to_numpy(dtype=float)
                if fraction == 0.0 and len(vals) != 1:
                    raise RuntimeError("Fraction-zero cardinality drift.")
                if fraction != 0.0 and len(vals) != 100:
                    raise RuntimeError(
                        f"{target} {program_id} {method}: fraction {fraction} has {len(vals)} rows"
                    )

                rate = float(np.mean(retained(vals, threshold)))
                rates.append(rate)

                curve_rows.append(
                    {
                        "target": target,
                        "program_id": program_id,
                        "method": method,
                        "fraction": float(fraction),
                        "preservation_call_rate": rate,
                        "failure_detection_rate": 1.0 - rate,
                        "threshold": threshold,
                    }
                )

            per_method_auc[method] = retention_auc(np.array(rates[:5], dtype=float))

        comp_methods = ["NetRep.cor.cor", "WGCNA.cor.kIM"]
        boot_delta = {x: np.empty(BOOTSTRAPS, dtype=float) for x in comp_methods}

        by_fraction = {
            float(f): g.loc[g["fraction"] == f]
            .sort_values("replicate_id")
            .reset_index(drop=True)
            for f in FRACTIONS[1:]
        }
        for f, gg in by_fraction.items():
            if len(gg) != 100:
                raise RuntimeError(f"{target} {program_id}: fraction {f} != 100 rows")

        for b in range(1, BOOTSTRAPS + 1):
            rng = rng_for_boot(target, program_id, b)
            idx = {
                float(f): rng.integers(0, 100, size=100, endpoint=False)
                for f in FRACTIONS[1:]
            }

            boot_auc = {}
            for method, col in METHODS.items():
                full_vals = pd.to_numeric(
                    by_fraction[1.0][col], errors="coerce"
                ).to_numpy(dtype=float)
                full_boot = full_vals[idx[1.0]]
                threshold_b = threshold_from_null(full_boot, ALPHA)

                rates_b = [
                    float(float(zero[col].iloc[0]) > threshold_b)
                ]
                for f in PERFORMANCE_FRACTIONS[1:]:
                    vals = pd.to_numeric(
                        by_fraction[float(f)][col], errors="coerce"
                    ).to_numpy(dtype=float)
                    vb = vals[idx[float(f)]]
                    rates_b.append(float(np.mean(vb > threshold_b)))

                boot_auc[method] = retention_auc(np.array(rates_b, dtype=float))

            for comp in comp_methods:
                boot_delta[comp][b - 1] = boot_auc["MTA"] - boot_auc[comp]

        for comp in comp_methods:
            d = boot_delta[comp]
            lo = float(np.quantile(d, 0.025))
            hi = float(np.quantile(d, 0.975))
            observed = float(per_method_auc["MTA"] - per_method_auc[comp])

            contrast_rows.append(
                {
                    "target": target,
                    "program_id": program_id,
                    "comparator": comp,
                    "retention_auc_mta": per_method_auc["MTA"],
                    "retention_auc_comparator": per_method_auc[comp],
                    "delta_auc_mta_minus_comparator": observed,
                    "delta_auc_boot_q025": lo,
                    "delta_auc_boot_q975": hi,
                    "interpretation": pair_interpretation(lo, hi),
                }
            )

            print(
                f"  {comp}: Δretention-AUC={observed:+.4f} "
                f"[{lo:+.4f},{hi:+.4f}] -> {pair_interpretation(lo,hi)}"
            )

    thresholds = pd.DataFrame(threshold_rows)
    curves = pd.DataFrame(curve_rows)
    monotone = pd.DataFrame(monotone_rows)
    contrasts = pd.DataFrame(contrast_rows)

    cc = curves.loc[curves["method"].isin(["NetRep.cor.cor", "WGCNA.cor.cor"])].copy()
    wide = cc.pivot_table(
        index=["target", "program_id", "fraction"],
        columns="method",
        values="preservation_call_rate",
        aggfunc="first",
    ).reset_index()
    if not np.allclose(
        wide["NetRep.cor.cor"].to_numpy(dtype=float),
        wide["WGCNA.cor.cor"].to_numpy(dtype=float),
        atol=0,
        rtol=0,
    ):
        raise RuntimeError("NetRep cor.cor and WGCNA cor.cor calibrated call curves differ unexpectedly.")

    summary = (
        contrasts.groupby(["comparator", "interpretation"])
        .size()
        .reset_index(name="n_target_program_pairs")
    )
    by_target = (
        contrasts.groupby(["target", "comparator", "interpretation"])
        .size()
        .reset_index(name="n_target_program_pairs")
    )

    thresholds_path = OUT_DIR / "scale_invariant_thresholds_v1.tsv"
    curves_path = OUT_DIR / "scale_invariant_call_curves_v1.tsv"
    monotone_path = OUT_DIR / "raw_monotone_redundancy_v1.tsv"
    contrasts_path = OUT_DIR / "scale_invariant_paired_contrasts_v1.tsv"
    summary_path = OUT_DIR / "scale_invariant_interpretation_counts_v1.tsv"
    by_target_path = OUT_DIR / "scale_invariant_interpretation_counts_by_target_v1.tsv"

    thresholds.to_csv(thresholds_path, sep="\t", index=False)
    curves.to_csv(curves_path, sep="\t", index=False)
    monotone.to_csv(monotone_path, sep="\t", index=False)
    contrasts.to_csv(contrasts_path, sep="\t", index=False)
    summary.to_csv(summary_path, sep="\t", index=False)
    by_target.to_csv(by_target_path, sep="\t", index=False)

    audit = {
        "script_version": SCRIPT_VERSION,
        "status": "SCALE_INVARIANT_CORRUPTION_AUDIT_COMPLETE",
        "post_hoc_relative_to_05j3": True,
        "false_preservation_alpha": ALPHA,
        "calibration_fraction": 1.0,
        "performance_fractions": PERFORMANCE_FRACTIONS.tolist(),
        "bootstrap_replicates": BOOTSTRAPS,
        "general_superiority_claim_permitted": False,
        "domain_specific_mapping_error_claim_only": True,
        "thresholds_file": str(thresholds_path),
        "curves_file": str(curves_path),
        "monotone_redundancy_file": str(monotone_path),
        "paired_contrasts_file": str(contrasts_path),
        "interpretation_counts_file": str(summary_path),
        "interpretation_counts_by_target_file": str(by_target_path),
    }
    master_path = OUT_DIR / "scale_invariant_corruption_audit_v1.json"
    master_path.write_text(json.dumps(audit, indent=2), encoding="utf-8")

    print("\n" + "=" * 158)
    print("05j4 SCALE-INVARIANT CORRUPTION AUDIT: COMPLETE")
    print("=" * 158)
    print("Overall interpretation counts:")
    print(summary.to_string(index=False))
    print("\nBy target:")
    print(by_target.to_string(index=False))
    print("\nRaw-score monotone redundancy summary:")
    print(
        monotone.groupby("comparator")["spearman_all_501_rows"]
        .agg(["median", "min", "max"])
        .to_string()
    )
    print()
    print(f"Contrasts: {contrasts_path}")
    print(f"Call curves: {curves_path}")
    print(f"Monotone audit: {monotone_path}")
    print(f"Master: {master_path}")
    print("=" * 158)


if __name__ == "__main__":
    main()
