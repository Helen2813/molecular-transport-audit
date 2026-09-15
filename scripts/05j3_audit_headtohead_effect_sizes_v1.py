from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd


SCRIPT_VERSION = "05j3-audit-headtohead-effect-sizes-v1-no-cli"

DATA_ROOT = Path(r"D:\paper4_tcbb_data")
ROOT = DATA_ROOT / "paper4_tcbb_structural_corruption_headtohead_v1"

CONTRASTS = ROOT / "structural_headtohead_paired_contrasts_v1.tsv"
METRICS = ROOT / "structural_headtohead_operating_metrics_v1.tsv"
FRACTIONS = ROOT / "structural_headtohead_fraction_summary_v1.tsv"
MASTER = ROOT / "structural_corruption_headtohead_v1.json"

OUT_DIR = DATA_ROOT / "paper4_tcbb_structural_corruption_headtohead_effect_audit_v1"


def require(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Required file does not exist: {path}")


def q(x: pd.Series, p: float) -> float:
    arr = pd.to_numeric(x, errors="coerce").to_numpy(dtype=float)
    arr = arr[np.isfinite(arr)]
    return float(np.quantile(arr, p)) if arr.size else np.nan


def main() -> None:
    print("=" * 154)
    print("Paper 4 / TCBB - audit structural head-to-head effect sizes")
    print("=" * 154)
    print(f"Script version: {SCRIPT_VERSION}")
    print()
    print("Scientific guard:")
    print("  New corruption mapping generated:                  NO")
    print("  New preservation statistic calculated:             NO")
    print("  New threshold introduced:                          NO")
    print("  Operation: summarize already-frozen 05j AUC/f50 contrasts only")
    print("=" * 154)

    for p in [CONTRASTS, METRICS, FRACTIONS, MASTER]:
        require(p)

    master = json.loads(MASTER.read_text(encoding="utf-8"))
    if master.get("status") != "STRUCTURAL_CORRUPTION_HEADTOHEAD_COMPLETE":
        raise RuntimeError("05j head-to-head master is not COMPLETE.")

    c = pd.read_csv(CONTRASTS, sep="\t", low_memory=False)
    m = pd.read_csv(METRICS, sep="\t", low_memory=False)
    f = pd.read_csv(FRACTIONS, sep="\t", low_memory=False)

    print("\n[1/4] Overall AUC effect-size summaries")

    overall_rows = []
    for comp, g in c.groupby("comparator", sort=True):
        d = pd.to_numeric(
            g["delta_auc_mta_minus_comparator"],
            errors="coerce",
        )
        finite = np.isfinite(d)
        gg = g.loc[finite].copy()
        d = d.loc[finite]

        overall_rows.append(
            {
                "comparator": comp,
                "n_pairs": int(len(d)),
                "median_delta_auc": float(np.median(d)),
                "mean_delta_auc": float(np.mean(d)),
                "q25_delta_auc": q(d, 0.25),
                "q75_delta_auc": q(d, 0.75),
                "min_delta_auc": float(np.min(d)),
                "max_delta_auc": float(np.max(d)),
                "mta_more_responsive_n": int(
                    (gg["auc_interpretation"] == "MTA_MORE_RESPONSIVE").sum()
                ),
                "no_clear_difference_n": int(
                    (gg["auc_interpretation"] == "NO_CLEAR_DIFFERENCE").sum()
                ),
                "comparator_more_responsive_n": int(
                    (gg["auc_interpretation"] == "COMPARATOR_MORE_RESPONSIVE").sum()
                ),
            }
        )

    overall = pd.DataFrame(overall_rows)
    print(overall.to_string(index=False))

    print("\n[2/4] By-target AUC summaries")

    target_rows = []
    for (target, comp), g in c.groupby(["target", "comparator"], sort=True):
        d = pd.to_numeric(g["delta_auc_mta_minus_comparator"], errors="coerce")
        finite = np.isfinite(d)
        gg = g.loc[finite].copy()
        d = d.loc[finite]

        target_rows.append(
            {
                "target": target,
                "comparator": comp,
                "n_pairs": int(len(d)),
                "median_delta_auc": float(np.median(d)),
                "q25_delta_auc": q(d, 0.25),
                "q75_delta_auc": q(d, 0.75),
                "mta_more_responsive_n": int(
                    (gg["auc_interpretation"] == "MTA_MORE_RESPONSIVE").sum()
                ),
                "no_clear_difference_n": int(
                    (gg["auc_interpretation"] == "NO_CLEAR_DIFFERENCE").sum()
                ),
                "comparator_more_responsive_n": int(
                    (gg["auc_interpretation"] == "COMPARATOR_MORE_RESPONSIVE").sum()
                ),
            }
        )

    by_target = pd.DataFrame(target_rows)
    print(by_target.to_string(index=False))

    print("\n[3/4] Exact pairs where MTA does NOT clearly outperform cor.cor")

    corcor = c.loc[c["comparator"] == "NetRep_cor.cor"].copy()
    nonwin = corcor.loc[
        corcor["auc_interpretation"] != "MTA_MORE_RESPONSIVE"
    ].copy()

    show_cols = [
        "target",
        "program_id",
        "comparator",
        "delta_auc_mta_minus_comparator",
        "delta_auc_boot_q025",
        "delta_auc_boot_q975",
        "delta_f50_mta_minus_comparator",
        "delta_f50_boot_q025",
        "delta_f50_boot_q975",
        "auc_interpretation",
    ]
    show_cols = [x for x in show_cols if x in nonwin.columns]
    print(nonwin[show_cols].sort_values(
        ["auc_interpretation", "target", "program_id"]
    ).to_string(index=False))

    print("\n[4/4] Per-pair compact effect table")

    # Join observed method AUC/f50 so the size of the contrast is interpretable.
    mwide = m.pivot_table(
        index=["target", "program_id"],
        columns="method",
        values=["auc_normalized", "f50"],
        aggfunc="first",
    )
    mwide.columns = [f"{a}.{b}" for a, b in mwide.columns]
    mwide = mwide.reset_index()

    compact = c.merge(
        mwide,
        on=["target", "program_id"],
        how="left",
    )

    compact = compact.sort_values(
        ["comparator", "target", "program_id"]
    )

    compact_path = OUT_DIR / "headtohead_effect_size_compact_v1.tsv"
    overall_path = OUT_DIR / "headtohead_effect_size_overall_v1.tsv"
    target_path = OUT_DIR / "headtohead_effect_size_by_target_v1.tsv"
    nonwin_path = OUT_DIR / "headtohead_corcor_nonwins_v1.tsv"

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    compact.to_csv(compact_path, sep="\t", index=False)
    overall.to_csv(overall_path, sep="\t", index=False)
    by_target.to_csv(target_path, sep="\t", index=False)
    nonwin.to_csv(nonwin_path, sep="\t", index=False)

    # A deliberately conservative descriptive flag; not a scientific threshold.
    # We do NOT use it for claims, only print absolute contrast scale.
    print("\nAbsolute delta-AUC magnitude distribution for NetRep cor.cor:")
    dabs = np.abs(
        pd.to_numeric(
            corcor["delta_auc_mta_minus_comparator"],
            errors="coerce",
        ).to_numpy(dtype=float)
    )
    dabs = dabs[np.isfinite(dabs)]
    print(
        f"  median |ΔAUC|={np.median(dabs):.4f}; "
        f"IQR=[{np.quantile(dabs,0.25):.4f},{np.quantile(dabs,0.75):.4f}]; "
        f"range=[{np.min(dabs):.4f},{np.max(dabs):.4f}]"
    )

    audit = {
        "script_version": SCRIPT_VERSION,
        "status": "HEADTOHEAD_EFFECT_SIZE_AUDIT_COMPLETE",
        "new_scientific_statistic_calculated": False,
        "new_threshold_introduced": False,
        "overall_table": str(overall_path),
        "by_target_table": str(target_path),
        "corcor_nonwins_table": str(nonwin_path),
        "compact_table": str(compact_path),
    }
    master_path = OUT_DIR / "headtohead_effect_size_audit_v1.json"
    master_path.write_text(json.dumps(audit, indent=2), encoding="utf-8")

    print("\n" + "=" * 154)
    print("05j3 HEAD-TO-HEAD EFFECT-SIZE AUDIT: COMPLETE")
    print("=" * 154)
    print("No new scientific statistic or decision threshold was introduced.")
    print(f"Overall: {overall_path}")
    print(f"By target: {target_path}")
    print(f"Non-wins: {nonwin_path}")
    print(f"Compact: {compact_path}")
    print("=" * 154)


if __name__ == "__main__":
    main()
