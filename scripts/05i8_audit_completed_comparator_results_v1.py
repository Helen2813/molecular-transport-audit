from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

SCRIPT_VERSION = "05i8-audit-completed-comparator-results-v1-no-cli"
DATA_ROOT = Path(r"D:\paper4_tcbb_data")

V7_ROOT = DATA_ROOT / "paper4_tcbb_postprimary_comparator_benchmark_v7"
V7_MASTER = V7_ROOT / "postprimary_comparator_benchmark_v7.json"
V7_NETREP = V7_ROOT / "netrep_all_targets_long_v7.tsv"
V7_WGCNA = V7_ROOT / "wgcna_modulePreservation_all_targets_v7.tsv"
CLASSIFICATION = DATA_ROOT / "paper4_tcbb_final_mapping_specificity_null_v2" / "primary_pooled_final_classification_v2.tsv"
INTERPRETATION_CONTRACT = DATA_ROOT / "paper4_tcbb_comparator_interpretation_headtohead_contract_v1" / "comparator_interpretation_headtohead_contract_v1.json"
CONTINUATION_CONTRACT = DATA_ROOT / "paper4_tcbb_partial_comparator_continuation_contract_v1" / "partial_comparator_continuation_contract_v1.json"

V4_SCANB_RDS = DATA_ROOT / "paper4_tcbb_postprimary_comparator_benchmark_v4" / "SCANB_GSE96058" / "netrep_modulePreservation_v4.rds"
V4_SCANB_TSV = DATA_ROOT / "paper4_tcbb_postprimary_comparator_benchmark_v4" / "SCANB_GSE96058" / "netrep_modulePreservation_long_v4.tsv"
V7_SCANB_RDS = V7_ROOT / "SCANB_GSE96058" / "netrep_modulePreservation_v7.rds"
V7_SCANB_TSV = V7_ROOT / "SCANB_GSE96058" / "netrep_modulePreservation_long_v7.tsv"
V4_METABRIC_RDS = DATA_ROOT / "paper4_tcbb_postprimary_comparator_benchmark_v4" / "METABRIC" / "netrep_modulePreservation_v4.rds"
V4_METABRIC_TSV = DATA_ROOT / "paper4_tcbb_postprimary_comparator_benchmark_v4" / "METABRIC" / "netrep_modulePreservation_long_v4.tsv"
OUT_DIR = DATA_ROOT / "paper4_tcbb_comparator_result_audit_v1"


def require(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Required file does not exist: {path}")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def boolish(x: object) -> bool:
    if isinstance(x, bool):
        return x
    return str(x).strip().lower() in {"1", "true", "yes", "y"}


def finite_spearman(x: pd.Series, y: pd.Series) -> tuple[float, int]:
    a = pd.to_numeric(x, errors="coerce").to_numpy(dtype=float)
    b = pd.to_numeric(y, errors="coerce").to_numpy(dtype=float)
    keep = np.isfinite(a) & np.isfinite(b)
    if keep.sum() < 3:
        return np.nan, int(keep.sum())
    return float(spearmanr(a[keep], b[keep]).statistic), int(keep.sum())


def main() -> None:
    print("=" * 156)
    print("Paper 4 / TCBB - audit completed NetRep/WGCNA comparator results under frozen 05i3 interpretation")
    print("=" * 156)
    print(f"Script version: {SCRIPT_VERSION}")
    print("\nScientific guard:")
    print("  New NetRep statistic calculated:                    NO")
    print("  New WGCNA preservation statistic calculated:        NO")
    print("  Comparator thresholds changed:                      NO")
    print("  Primary MTA classification changed:                 NO")
    print("=" * 156)

    for p in [V7_MASTER, V7_NETREP, V7_WGCNA, CLASSIFICATION, INTERPRETATION_CONTRACT,
              CONTINUATION_CONTRACT, V4_SCANB_RDS, V4_SCANB_TSV, V7_SCANB_RDS, V7_SCANB_TSV]:
        require(p)

    master = json.loads(V7_MASTER.read_text(encoding="utf-8"))
    interp = json.loads(INTERPRETATION_CONTRACT.read_text(encoding="utf-8"))
    cont = json.loads(CONTINUATION_CONTRACT.read_text(encoding="utf-8"))

    if master.get("status") != "POSTPRIMARY_COMPARATOR_BENCHMARK_COMPLETE":
        raise RuntimeError("v7 comparator master is not COMPLETE.")
    if interp.get("status") != "FROZEN_BEFORE_FIRST_NETREP_OR_WGCNA_COMPARATOR_RESULT":
        raise RuntimeError("05i3 interpretation contract has unexpected status.")
    if cont.get("status") != "FROZEN_AFTER_SCANB_NETREP_BYTES_EXIST_BEFORE_NUMERICAL_INSPECTION":
        raise RuntimeError("05i5a continuation contract has unexpected status.")

    print("\n[1/5] Provenance / byte replay")
    frozen_rds = cont["authoritative_completed_artifact"]["rds"]["sha256"]
    frozen_tsv = cont["authoritative_completed_artifact"]["tsv"]["sha256"]
    hashes = {
        "v4_scanb_rds": sha256(V4_SCANB_RDS),
        "v4_scanb_tsv": sha256(V4_SCANB_TSV),
        "v7_scanb_rds": sha256(V7_SCANB_RDS),
        "v7_scanb_tsv": sha256(V7_SCANB_TSV),
    }
    if hashes["v4_scanb_rds"] != frozen_rds or hashes["v7_scanb_rds"] != frozen_rds:
        raise RuntimeError("SCAN-B NetRep RDS byte identity failed.")
    if hashes["v4_scanb_tsv"] != frozen_tsv or hashes["v7_scanb_tsv"] != frozen_tsv:
        raise RuntimeError("SCAN-B NetRep TSV byte identity failed.")
    print("  SCAN-B NetRep v4 -> v7 byte identity: PASS")

    accidental = {"rds_exists": V4_METABRIC_RDS.exists(), "tsv_exists": V4_METABRIC_TSV.exists()}
    print(f"  Accidental old v4 METABRIC serialized result: {accidental}")
    if any(accidental.values()):
        print("  NOTE: preserve those files; do not delete or substitute them for v7.")

    print("\n[2/5] Loading completed outputs")
    net = pd.read_csv(V7_NETREP, sep="\t", low_memory=False)
    wgc = pd.read_csv(V7_WGCNA, sep="\t", low_memory=False)
    cls = pd.read_csv(CLASSIFICATION, sep="\t", low_memory=False)

    required_net = {"target", "program_id", "statistic", "observed", "raw_p", "bh_q"}
    if not required_net.issubset(net.columns):
        raise RuntimeError(f"NetRep missing columns: {sorted(required_net - set(net.columns))}")

    cls["primary_assessable_bool"] = cls["primary_assessable"].map(boolish)
    cls_assess = cls.loc[cls["primary_assessable_bool"]].copy()
    stats = sorted(net["statistic"].astype(str).unique())
    print(f"  NetRep statistics present ({len(stats)}): {stats}")
    print(f"  Primary-assessable MTA pairs: {len(cls_assess)}")

    net_assess = net.merge(cls_assess[["target", "program_id"]], on=["target", "program_id"], how="inner")
    rows = []
    for target in sorted(net_assess["target"].unique()):
        for stat in stats:
            s = net_assess[(net_assess["target"] == target) & (net_assess["statistic"] == stat)]
            rows.append({
                "target": target,
                "statistic": stat,
                "n": int(len(s)),
                "raw_p_lt_0p05": int((pd.to_numeric(s["raw_p"], errors="coerce") < 0.05).sum()),
                "bh_q_lt_0p05": int((pd.to_numeric(s["bh_q"], errors="coerce") < 0.05).sum()),
            })
    support = pd.DataFrame(rows)
    print("\n[3/5] NetRep support counts")
    print(support.to_string(index=False))

    piv = net_assess.pivot_table(index=["target", "program_id"], columns="statistic", values=["observed", "bh_q"], aggfunc="first")
    piv.columns = [f"{a}.{b}" for a, b in piv.columns]
    piv = piv.reset_index()

    for c in ["bh_q.cor.cor", "bh_q.avg.contrib", "bh_q.cor.contrib"]:
        if c not in piv.columns:
            piv[c] = np.nan
    piv["netrep_structural_support"] = pd.to_numeric(piv["bh_q.cor.cor"], errors="coerce") < 0.05
    piv["netrep_loading_support"] = (
        (pd.to_numeric(piv["bh_q.avg.contrib"], errors="coerce") < 0.05)
        & (pd.to_numeric(piv["bh_q.cor.contrib"], errors="coerce") < 0.05)
    )

    for target in sorted(piv["target"].unique()):
        s = piv[piv["target"] == target]
        print(f"\n  {target}: structural support {int(s['netrep_structural_support'].sum())}/{len(s)}; loading support {int(s['netrep_loading_support'].sum())}/{len(s)}")

    if "Zsummary.pres" not in wgc.columns or "medianRank.pres" not in wgc.columns:
        raise RuntimeError("WGCNA table lacks Zsummary.pres / medianRank.pres.")
    wgc["Zsummary.pres"] = pd.to_numeric(wgc["Zsummary.pres"], errors="coerce")
    wgc["medianRank.pres"] = pd.to_numeric(wgc["medianRank.pres"], errors="coerce")
    wgc["Zsummary_category"] = np.select(
        [wgc["Zsummary.pres"] < 2, (wgc["Zsummary.pres"] >= 2) & (wgc["Zsummary.pres"] < 10), wgc["Zsummary.pres"] >= 10],
        ["NO_EVIDENCE_LT2", "MODERATE_2_TO_LT10", "STRONG_GE10"],
        default="NOT_FINITE",
    )
    wgc_assess = wgc.merge(cls_assess[["target", "program_id"]], on=["target", "program_id"], how="inner")
    print("\nWGCNA Zsummary categories:")
    for target in sorted(wgc_assess["target"].unique()):
        s = wgc_assess[wgc_assess["target"] == target]
        print(f"  {target}: {s['Zsummary_category'].value_counts().to_dict()}")

    print("\n[4/5] Predeclared continuous rank concordance")
    compact = cls_assess[["target", "program_id", "primary_classification", "rho_edge", "rho_load", "coverage", "n_evaluable_genes"]].copy()
    compact = compact.merge(piv, on=["target", "program_id"], how="left")
    keep_w = [c for c in ["target", "program_id", "Zsummary.pres", "medianRank.pres", "Zdensity.pres", "Zconnectivity.pres", "Zsummary_category"] if c in wgc.columns]
    compact = compact.merge(wgc[keep_w], on=["target", "program_id"], how="left")

    con_rows = []
    comps = [
        ("MTA rho_edge vs NetRep cor.cor", "rho_edge", "observed.cor.cor", False),
        ("MTA rho_load vs NetRep cor.contrib", "rho_load", "observed.cor.contrib", False),
        ("MTA rho_edge vs WGCNA Zsummary", "rho_edge", "Zsummary.pres", False),
        ("MTA rho_edge vs -WGCNA medianRank", "rho_edge", "medianRank.pres", True),
    ]
    for target in sorted(compact["target"].unique()):
        s = compact[compact["target"] == target]
        for label, xc, yc, negate in comps:
            if yc not in s.columns:
                continue
            y = -pd.to_numeric(s[yc], errors="coerce") if negate else s[yc]
            rho, n = finite_spearman(s[xc], y)
            con_rows.append({"target": target, "comparison": label, "spearman_rho": rho, "n_programs": n})
    concordance = pd.DataFrame(con_rows)
    print(concordance.to_string(index=False))

    print("\n[5/5] Compact novelty-risk table")
    show = ["target", "program_id", "primary_classification", "rho_edge", "rho_load"]
    for c in ["observed.cor.cor", "bh_q.cor.cor", "observed.avg.contrib", "bh_q.avg.contrib", "observed.cor.contrib", "bh_q.cor.contrib", "netrep_structural_support", "netrep_loading_support", "Zsummary.pres", "Zsummary_category", "medianRank.pres"]:
        if c in compact.columns:
            show.append(c)
    print(compact[show].to_string(index=False))

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    support_path = OUT_DIR / "netrep_statistic_support_counts_v1.tsv"
    compact_path = OUT_DIR / "comparator_compact_primary_assessable_v1.tsv"
    concordance_path = OUT_DIR / "cross_method_rank_concordance_v1.tsv"
    support.to_csv(support_path, sep="\t", index=False)
    compact.to_csv(compact_path, sep="\t", index=False)
    concordance.to_csv(concordance_path, sep="\t", index=False)
    (OUT_DIR / "comparator_result_audit_v1.json").write_text(json.dumps({
        "script_version": SCRIPT_VERSION,
        "status": "COMPLETED_COMPARATOR_RESULT_AUDIT_COMPLETE",
        "scientific_statistics_recomputed": False,
        "scanb_netrep_v4_to_v7_byte_identity": True,
        "accidental_v4_metabric_serialized_artifacts": accidental,
        "primary_assessable_pairs": int(len(cls_assess)),
        "primary_mta_classification_changed": False,
    }, indent=2), encoding="utf-8")

    print("\n" + "=" * 156)
    print("05i8 COMPLETED COMPARATOR RESULT AUDIT: COMPLETE")
    print("=" * 156)
    print("No comparator statistic was recomputed.")
    print("This output is sufficient for the frozen GREEN/YELLOW/RED novelty assessment.")
    print("=" * 156)


if __name__ == "__main__":
    main()
