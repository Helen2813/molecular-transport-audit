from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd


SCRIPT_VERSION = "04i3-freeze-exact-variation-evaluability-correction-v1-no-cli"

DATA_ROOT = Path(r"D:\paper4_tcbb_data")

PRIOR_EVAL_DIR = DATA_ROOT / "paper4_tcbb_target_statistical_evaluability_amendment_v2"
PRIOR_SCANB_EVAL = PRIOR_EVAL_DIR / "scanb_target_statistical_evaluability_v2.tsv"
PRIOR_MET_EVAL = PRIOR_EVAL_DIR / "metabric_target_statistical_evaluability_v2.tsv"

CLASSICAL_DEGENERACY_AUDIT = (
    DATA_ROOT
    / "paper4_tcbb_classical_bootstrap_degeneracy_audit_v1"
    / "classical_bootstrap_degeneracy_audit_v1.json"
)
BAYESIAN_AMENDMENT = (
    DATA_ROOT
    / "paper4_tcbb_bayesian_bootstrap_amendment_v1"
    / "bayesian_bootstrap_uncertainty_amendment_v1.json"
)

SOURCE_UNIVERSE = (
    DATA_ROOT
    / "paper4_tcbb_tcga_source_universe_v1"
    / "tcga_source_gene_universe_frozen_v1.tsv"
)
WEIGHTS = (
    DATA_ROOT
    / "paper4_tcbb_frozen_source_programs_v1"
    / "tcga_frozen_source_program_weights_v1.tsv"
)

SCANB_EXPR = (
    DATA_ROOT
    / "SCANB_GSE96058"
    / "GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv.gz"
)
SCANB_PRIMARY = (
    DATA_ROOT
    / "paper4_tcbb_scanb_pairing_contract_v1"
    / "scanb_primary_profiles_frozen_v1.tsv"
)
METABRIC_EXPR = (
    DATA_ROOT
    / "paper4_tcbb_input_audit_v1"
    / "staged_continuous_inputs"
    / "METABRIC"
    / "data_mrna_illumina_microarray.txt"
)

OUT_DIR = DATA_ROOT / "paper4_tcbb_exact_variation_evaluability_correction_v1"

MIN_EVALUABLE_GENES = 30
MIN_COVERAGE = 0.80


def require(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Required file does not exist: {path}")


def canon_symbol(x: object) -> str:
    if x is None:
        return ""
    return " ".join(str(x).strip().split()).upper()


def load_scanb_exact_variation(
    universe: set[str],
    primary_titles: list[str],
) -> pd.DataFrame:
    parts = []

    for chunk in pd.read_csv(
        SCANB_EXPR,
        compression="gzip",
        chunksize=1500,
        low_memory=False,
    ):
        symbol_col = chunk.columns[0]
        symbols = chunk[symbol_col].map(canon_symbol)
        mask = symbols.isin(universe)
        if not mask.any():
            continue

        missing = [c for c in primary_titles if c not in chunk.columns]
        if missing:
            raise RuntimeError(f"Missing SCAN-B primary columns: {missing[:20]}")

        sub = chunk.loc[mask, [symbol_col] + primary_titles]
        arr = sub[primary_titles].to_numpy(dtype=np.float64)

        # Frozen preprocessing.
        fpkm = np.maximum(np.exp2(arr) - 0.1, 0.0)
        transformed = np.log2(fpkm + 1.0)

        tmp = pd.DataFrame(transformed, columns=primary_titles)
        tmp.insert(0, "Hugo_Symbol", symbols.loc[mask].to_numpy())
        parts.append(tmp)

    if not parts:
        raise RuntimeError("No SCAN-B frozen-universe genes found.")

    df = pd.concat(parts, ignore_index=True)
    df = df.groupby("Hugo_Symbol", sort=False, as_index=True).mean(numeric_only=True)

    x = df.to_numpy(dtype=np.float64)
    finite = np.isfinite(x).all(axis=1)

    exact_min = np.full(len(df), np.nan, dtype=np.float64)
    exact_max = np.full(len(df), np.nan, dtype=np.float64)
    exact_range = np.full(len(df), np.nan, dtype=np.float64)
    sample_sd = np.full(len(df), np.nan, dtype=np.float64)
    n_unique = np.full(len(df), -1, dtype=np.int64)

    for i in np.where(finite)[0]:
        v = x[i]
        vmin = float(np.min(v))
        vmax = float(np.max(v))
        exact_min[i] = vmin
        exact_max[i] = vmax
        exact_range[i] = vmax - vmin
        sample_sd[i] = float(np.std(v, ddof=1))
        n_unique[i] = int(np.unique(v).size)

    exact_variable = finite & (exact_max > exact_min)

    return pd.DataFrame(
        {
            "Hugo_Symbol": df.index.to_numpy(),
            "all_finite": finite.astype(int),
            "exact_min": exact_min,
            "exact_max": exact_max,
            "exact_range": exact_range,
            "n_unique_exact_values": n_unique,
            "numpy_sample_sd_diagnostic": sample_sd,
            "exactly_variable": exact_variable.astype(int),
            "statistically_evaluable_corrected": exact_variable.astype(int),
        }
    )


def load_metabric_exact_variation(universe: set[str]) -> pd.DataFrame:
    parts = []
    sample_cols = None

    for chunk in pd.read_csv(
        METABRIC_EXPR,
        sep="\t",
        chunksize=1500,
        low_memory=False,
    ):
        if sample_cols is None:
            sample_cols = list(chunk.columns[2:])

        symbols = chunk["Hugo_Symbol"].map(canon_symbol)
        mask = symbols.isin(universe)
        if not mask.any():
            continue

        arr = chunk.loc[mask, sample_cols].to_numpy(dtype=np.float64)
        tmp = pd.DataFrame(arr, columns=sample_cols)
        tmp.insert(0, "Hugo_Symbol", symbols.loc[mask].to_numpy())
        parts.append(tmp)

    if not parts:
        raise RuntimeError("No METABRIC frozen-universe genes found.")

    df = pd.concat(parts, ignore_index=True)
    df = df.groupby("Hugo_Symbol", sort=False, as_index=True).mean(numeric_only=True)

    x = df.to_numpy(dtype=np.float64)
    finite = np.isfinite(x).all(axis=1)

    exact_min = np.full(len(df), np.nan, dtype=np.float64)
    exact_max = np.full(len(df), np.nan, dtype=np.float64)
    exact_range = np.full(len(df), np.nan, dtype=np.float64)
    sample_sd = np.full(len(df), np.nan, dtype=np.float64)
    n_unique = np.full(len(df), -1, dtype=np.int64)

    for i in np.where(finite)[0]:
        v = x[i]
        vmin = float(np.min(v))
        vmax = float(np.max(v))
        exact_min[i] = vmin
        exact_max[i] = vmax
        exact_range[i] = vmax - vmin
        sample_sd[i] = float(np.std(v, ddof=1))
        n_unique[i] = int(np.unique(v).size)

    exact_variable = finite & (exact_max > exact_min)

    return pd.DataFrame(
        {
            "Hugo_Symbol": df.index.to_numpy(),
            "all_finite": finite.astype(int),
            "exact_min": exact_min,
            "exact_max": exact_max,
            "exact_range": exact_range,
            "n_unique_exact_values": n_unique,
            "numpy_sample_sd_diagnostic": sample_sd,
            "exactly_variable": exact_variable.astype(int),
            "statistically_evaluable_corrected": exact_variable.astype(int),
        }
    )


def program_impact(
    target_name: str,
    target_stats: pd.DataFrame,
    weights: pd.DataFrame,
) -> pd.DataFrame:
    represented = set(target_stats["Hugo_Symbol"])
    eligible = set(
        target_stats.loc[
            target_stats["statistically_evaluable_corrected"].astype(int) == 1,
            "Hugo_Symbol",
        ]
    )

    rows = []
    for program_id in sorted(weights["program_id"].unique()):
        genes = [
            canon_symbol(x)
            for x in weights.loc[weights["program_id"] == program_id, "Hugo_Symbol"]
        ]
        n_frozen = len(genes)
        n_mapped = sum(g in represented for g in genes)
        n_eval = sum(g in eligible for g in genes)
        coverage = n_eval / n_frozen

        rows.append(
            {
                "target": target_name,
                "program_id": program_id,
                "n_frozen_genes": n_frozen,
                "n_identifier_mapped": n_mapped,
                "n_corrected_evaluable": n_eval,
                "corrected_coverage": coverage,
                "n_mapped_but_ineligible": n_mapped - n_eval,
                "primary_assessable_corrected": int(
                    n_eval >= MIN_EVALUABLE_GENES
                    and coverage >= MIN_COVERAGE
                ),
            }
        )

    return pd.DataFrame(rows)


def main() -> None:
    print("=" * 142)
    print("Paper 4 / TCBB - freeze exact-variation correction to target statistical evaluability")
    print("=" * 142)
    print(f"Script version: {SCRIPT_VERSION}")
    print()
    print("Trigger:")
    print("  05c1 showed that SCAN-B M001 genes LOC391322 and GSTT1 are EXACTLY constant")
    print("  in the original frozen target samples, despite having passed the prior SD>0 test.")
    print("  The prior SD test was therefore vulnerable to floating-point reduction noise.")
    print()
    print("Scientific correction frozen BEFORE any recomputed preservation result:")
    print("  A target gene is evaluable iff:")
    print("    (1) all transformed values on the analysis sample set are finite; AND")
    print("    (2) exact max(value) > exact min(value).")
    print("  No arbitrary variance tolerance is introduced.")
    print()
    print("Consequences:")
    print("  Prior 05a/05b v1 effects/classes remain in the audit trail but are PROVISIONAL")
    print("  until rerun with the corrected universal gene-eligibility rule.")
    print("  The Bayesian-bootstrap replacement is NOT used at this stage.")
    print("=" * 142)

    for p in [
        PRIOR_SCANB_EVAL,
        PRIOR_MET_EVAL,
        CLASSICAL_DEGENERACY_AUDIT,
        BAYESIAN_AMENDMENT,
        SOURCE_UNIVERSE,
        WEIGHTS,
        SCANB_EXPR,
        SCANB_PRIMARY,
        METABRIC_EXPR,
    ]:
        require(p)

    deg = json.loads(CLASSICAL_DEGENERACY_AUDIT.read_text(encoding="utf-8"))
    if deg.get("status") != "CLASSICAL_FIXED_GENE_BOOTSTRAP_DEGENERACY_CONFIRMED":
        raise RuntimeError("05c1 degeneracy audit has unexpected status.")

    if int(deg.get("valid_attempts", -1)) != 0:
        raise RuntimeError("05c1 unexpectedly reports a valid classical bootstrap draw.")

    bayes = json.loads(BAYESIAN_AMENDMENT.read_text(encoding="utf-8"))
    if bayes.get("status") != "FROZEN_BEFORE_ANY_VALID_TARGET_BOOTSTRAP_PRESERVATION_RESULT":
        raise RuntimeError("05c2 Bayesian amendment has unexpected status.")

    universe_df = pd.read_csv(SOURCE_UNIVERSE, sep="\t", dtype=str).fillna("")
    universe = {canon_symbol(x) for x in universe_df["Hugo_Symbol"]}
    if len(universe) != 10_000:
        raise RuntimeError(f"Expected 10,000 frozen source genes, found {len(universe)}.")

    weights = pd.read_csv(WEIGHTS, sep="\t", dtype=str).fillna("")
    weights["Hugo_Symbol"] = weights["Hugo_Symbol"].map(canon_symbol)

    primary = pd.read_csv(SCANB_PRIMARY, sep="\t", dtype=str).fillna("")
    primary_titles = list(primary["primary_title"])
    if len(primary_titles) != 3273 or len(set(primary_titles)) != 3273:
        raise RuntimeError("SCAN-B frozen primary manifest is not exactly 3,273 unique samples.")

    print("\n[1/3] Re-auditing SCAN-B with exact-variation rule ...")
    scanb = load_scanb_exact_variation(universe, primary_titles)

    print("\n[2/3] Re-auditing METABRIC with exact-variation rule ...")
    met = load_metabric_exact_variation(universe)

    prior_scanb = pd.read_csv(PRIOR_SCANB_EVAL, sep="\t", dtype=str).fillna("")
    prior_met = pd.read_csv(PRIOR_MET_EVAL, sep="\t", dtype=str).fillna("")

    prior_scanb["Hugo_Symbol"] = prior_scanb["Hugo_Symbol"].map(canon_symbol)
    prior_met["Hugo_Symbol"] = prior_met["Hugo_Symbol"].map(canon_symbol)

    prior_scanb_eval = set(
        prior_scanb.loc[
            prior_scanb["statistically_evaluable"].astype(str) == "1",
            "Hugo_Symbol",
        ]
    )
    prior_met_eval = set(
        prior_met.loc[
            prior_met["statistically_evaluable"].astype(str) == "1",
            "Hugo_Symbol",
        ]
    )

    new_scanb_eval = set(
        scanb.loc[
            scanb["statistically_evaluable_corrected"].astype(int) == 1,
            "Hugo_Symbol",
        ]
    )
    new_met_eval = set(
        met.loc[
            met["statistically_evaluable_corrected"].astype(int) == 1,
            "Hugo_Symbol",
        ]
    )

    scanb_newly_ineligible = sorted(prior_scanb_eval - new_scanb_eval)
    met_newly_ineligible = sorted(prior_met_eval - new_met_eval)

    print(
        f"  SCAN-B prior evaluable={len(prior_scanb_eval):,}; "
        f"corrected evaluable={len(new_scanb_eval):,}; "
        f"newly ineligible={len(scanb_newly_ineligible)}"
    )
    print(
        f"  METABRIC prior evaluable={len(prior_met_eval):,}; "
        f"corrected evaluable={len(new_met_eval):,}; "
        f"newly ineligible={len(met_newly_ineligible)}"
    )

    print("\nNewly ineligible after exact-variation correction:")
    if not scanb_newly_ineligible and not met_newly_ineligible:
        print("  NONE")
    else:
        for g in scanb_newly_ineligible:
            print(f"  SCANB_GSE96058: {g}")
        for g in met_newly_ineligible:
            print(f"  METABRIC: {g}")

    impact = pd.concat(
        [
            program_impact("SCANB_GSE96058", scanb, weights),
            program_impact("METABRIC", met, weights),
        ],
        ignore_index=True,
    )

    print("\n[3/3] Corrected program assessability:")
    for target in ["SCANB_GSE96058", "METABRIC"]:
        print(f"  {target}:")
        sub = impact.loc[impact["target"] == target]
        for r in sub.itertuples(index=False):
            print(
                f"    {r.program_id}: corrected={r.n_corrected_evaluable:4d}/"
                f"{r.n_frozen_genes:4d} coverage={r.corrected_coverage:.3f} "
                f"assessable={'YES' if r.primary_assessable_corrected else 'NO'}"
            )

    program_by_gene = defaultdict(list)
    for r in weights.itertuples(index=False):
        program_by_gene[canon_symbol(r.Hugo_Symbol)].append(r.program_id)

    newly_rows = []
    for target, genes, stats in [
        ("SCANB_GSE96058", scanb_newly_ineligible, scanb),
        ("METABRIC", met_newly_ineligible, met),
    ]:
        stat_by_gene = stats.set_index("Hugo_Symbol")
        for g in genes:
            s = stat_by_gene.loc[g]
            newly_rows.append(
                {
                    "target": target,
                    "Hugo_Symbol": g,
                    "exact_min": s["exact_min"],
                    "exact_max": s["exact_max"],
                    "exact_range": s["exact_range"],
                    "n_unique_exact_values": int(s["n_unique_exact_values"]),
                    "numpy_sample_sd_diagnostic": s["numpy_sample_sd_diagnostic"],
                    "programs": ";".join(program_by_gene.get(g, [])),
                }
            )

    newly_df = pd.DataFrame(newly_rows)

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    scanb_out = OUT_DIR / "scanb_target_exact_variation_evaluability_v1.tsv"
    met_out = OUT_DIR / "metabric_target_exact_variation_evaluability_v1.tsv"
    new_out = OUT_DIR / "newly_ineligible_exact_constant_genes_v1.tsv"
    impact_out = OUT_DIR / "program_assessability_after_exact_variation_correction_v1.tsv"

    scanb.to_csv(scanb_out, sep="\t", index=False)
    met.to_csv(met_out, sep="\t", index=False)
    newly_df.to_csv(new_out, sep="\t", index=False)
    impact.to_csv(impact_out, sep="\t", index=False)

    contract = {
        "contract_id": "paper4-tcbb-exact-variation-evaluability-correction-v1",
        "script_version": SCRIPT_VERSION,
        "status": "FROZEN_CORRECTION_BEFORE_RECOMPUTED_PRESERVATION",
        "trigger": (
            "05c1 identified exact-constant SCAN-B genes admitted by the prior "
            "floating-point SD>0 evaluability test."
        ),
        "corrected_rule": {
            "finite": "all transformed values on the analysis sample set must be finite",
            "exact_variation": "exact max(value) > exact min(value)",
            "arbitrary_tolerance": None,
            "reason": (
                "Pearson correlation/z-standardization/PCA require mathematical "
                "variation; exact constancy must not be inferred from a floating "
                "variance reduction that can become spuriously positive."
            ),
        },
        "universal_application": (
            "apply to every frozen-universe target gene in SCAN-B and METABRIC, "
            "not only genes implicated by the bootstrap failure"
        ),
        "prior_results": {
            "05a_v1": "provisional; must be rerun under corrected eligibility",
            "05b_v1": "provisional; must be rerun under corrected eligibility",
            "05c_v1_v2": "failed before valid bootstrap uncertainty",
            "05c2_bayesian_amendment": (
                "archival only; not used unless classical bootstrap remains "
                "non-estimable after corrected exact-variation eligibility"
            ),
        },
        "bootstrap_plan_after_correction": (
            "return first to the original 04a ordinary target-sample bootstrap; "
            "Bayesian bootstrap is not the default after removal of exact-constant genes"
        ),
        "primary_thresholds_changed": False,
        "module_membership_changed": False,
        "outcome_firewall_changed": False,
        "counts": {
            "SCANB_prior_evaluable": len(prior_scanb_eval),
            "SCANB_corrected_evaluable": len(new_scanb_eval),
            "SCANB_newly_ineligible": len(scanb_newly_ineligible),
            "METABRIC_prior_evaluable": len(prior_met_eval),
            "METABRIC_corrected_evaluable": len(new_met_eval),
            "METABRIC_newly_ineligible": len(met_newly_ineligible),
        },
    }

    json_out = OUT_DIR / "exact_variation_evaluability_correction_v1.json"
    json_out.write_text(json.dumps(contract, indent=2), encoding="utf-8")

    print("\n" + "=" * 142)
    print("04i3 EXACT-VARIATION EVALUABILITY CORRECTION: PASS")
    print("=" * 142)
    print(f"SCAN-B corrected evaluable genes:   {len(new_scanb_eval):,}")
    print(f"METABRIC corrected evaluable genes: {len(new_met_eval):,}")
    print(f"SCAN-B newly ineligible:             {len(scanb_newly_ineligible)}")
    print(f"METABRIC newly ineligible:           {len(met_newly_ineligible)}")
    print("Prior 05a/05b v1 final-use status:   PROVISIONAL / MUST RERUN")
    print("Bayesian bootstrap current status:   ARCHIVAL / NOT USED")
    print("Classical bootstrap after correction:RETRY ORIGINAL 04a SCHEME")
    print()
    print(f"Contract: {json_out}")
    print(f"Newly ineligible genes: {new_out}")
    print(f"Program impact: {impact_out}")
    print("=" * 142)


if __name__ == "__main__":
    main()
