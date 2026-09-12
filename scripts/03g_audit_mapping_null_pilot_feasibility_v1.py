from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd


SCRIPT_VERSION = "03g-audit-mapping-null-pilot-feasibility-v1-no-cli"

DATA_ROOT = Path(r"D:\paper4_tcbb_data")

WEIGHTS = (
    DATA_ROOT
    / "paper4_tcbb_frozen_source_programs_v1"
    / "tcga_frozen_source_program_weights_v1.tsv"
)
COVERAGE = (
    DATA_ROOT
    / "paper4_tcbb_target_mapping_v1"
    / "target_program_gene_coverage_v1.tsv"
)
MARGINAL_ROOT = DATA_ROOT / "paper4_tcbb_target_marginal_metrics_v2"
SCANB_FEATURES = MARGINAL_ROOT / "scanb_target_marginal_matching_features_v2.tsv"
METABRIC_FEATURES = MARGINAL_ROOT / "metabric_target_marginal_matching_features_v2.tsv"

NULL_CONTRACT = (
    DATA_ROOT
    / "paper4_tcbb_mapping_null_contract_v1"
    / "structure_preserving_mapping_null_contract_v1.json"
)
PILOT_CONTRACT = (
    DATA_ROOT
    / "paper4_tcbb_mapping_null_pilot_contract_v1"
    / "mapping_null_pilot_contract_v1.json"
)

OUT_DIR = DATA_ROOT / "paper4_tcbb_mapping_null_pilot_audit_v1"


def require(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Required file does not exist: {path}")


def norm(x: object) -> str:
    if pd.isna(x):
        return ""
    return str(x).strip().upper()


def manhattan(a_s: int, a_t: int, b_s: int, b_t: int) -> int:
    return abs(a_s - b_s) + abs(a_t - b_t)


def sample_mapping(
    slots: pd.DataFrame,
    candidates: pd.DataFrame,
    rng: np.random.Generator,
    max_search_distance: int = 2,
) -> dict:
    """
    Frozen 03d assignment rule:
      - randomize slot order;
      - without replacement;
      - exact stratum first;
      - if empty, expand by increasing Manhattan distance;
      - ties among candidates at the same minimum distance are random.
    """
    if len(candidates) < len(slots):
        return {
            "success": False,
            "failure_reason": "candidate_pool_smaller_than_panel",
        }

    cand = candidates[
        [
            "Hugo_Symbol",
            "source_composite_bin",
            "target_composite_bin",
        ]
    ].copy().reset_index(drop=True)

    unused = np.ones(len(cand), dtype=bool)
    order = rng.permutation(len(slots))

    dists = np.empty(len(slots), dtype=int)
    assigned_symbols = [None] * len(slots)

    for slot_idx in order:
        sb = int(slots.iloc[slot_idx]["source_composite_bin"])
        tb = int(slots.iloc[slot_idx]["target_composite_bin"])

        available_idx = np.flatnonzero(unused)
        if available_idx.size == 0:
            return {
                "success": False,
                "failure_reason": "candidate_pool_exhausted",
            }

        cs = cand.loc[available_idx, "source_composite_bin"].to_numpy(dtype=int)
        ct = cand.loc[available_idx, "target_composite_bin"].to_numpy(dtype=int)
        dist = np.abs(cs - sb) + np.abs(ct - tb)

        min_dist = int(dist.min())
        if min_dist > max_search_distance:
            return {
                "success": False,
                "failure_reason": f"minimum_available_distance_{min_dist}_gt_{max_search_distance}",
            }

        nearest = available_idx[dist == min_dist]
        chosen = int(rng.choice(nearest))

        unused[chosen] = False
        dists[slot_idx] = min_dist
        assigned_symbols[slot_idx] = cand.iloc[chosen]["Hugo_Symbol"]

    return {
        "success": True,
        "distances": dists,
        "assigned_symbols": assigned_symbols,
    }


def summarize_panel(
    target: str,
    program_id: str,
    panel_id: int,
    result: dict,
    min_frac_le1: float,
    max_dist_allowed: int,
    primary_assessable: int,
) -> dict:
    if not result["success"]:
        return {
            "target": target,
            "program_id": program_id,
            "panel_id": panel_id,
            "assignment_success": 0,
            "failure_reason": result["failure_reason"],
            "fraction_distance_0": np.nan,
            "fraction_distance_le1": np.nan,
            "maximum_distance": np.nan,
            "mean_distance": np.nan,
            "panel_valid": 0,
            "primary_assessable": primary_assessable,
        }

    d = result["distances"]
    frac0 = float(np.mean(d == 0))
    frac_le1 = float(np.mean(d <= 1))
    maxd = int(d.max())
    meand = float(d.mean())

    valid = int(frac_le1 >= min_frac_le1 and maxd <= max_dist_allowed)

    return {
        "target": target,
        "program_id": program_id,
        "panel_id": panel_id,
        "assignment_success": 1,
        "failure_reason": "",
        "fraction_distance_0": frac0,
        "fraction_distance_le1": frac_le1,
        "maximum_distance": maxd,
        "mean_distance": meand,
        "panel_valid": valid,
        "primary_assessable": primary_assessable,
    }


def main() -> None:
    print("=" * 126)
    print("Paper 4 / TCBB - audit pilot feasibility of the frozen 03d matched-mapping null")
    print("=" * 126)
    print(f"Script version: {SCRIPT_VERSION}")
    print()
    print("Scientific guard:")
    print("  Target outcomes loaded:                              NO")
    print("  Target clinical characteristics loaded:              NO")
    print("  Target gene-gene correlations calculated:            NO")
    print("  Target PCA/loadings calculated:                       NO")
    print("  Target preservation statistics calculated:           NO")
    print("  Target expression values read here:                  NO (uses frozen marginal tables)")
    print("  Operation: mapping-distance feasibility only")
    print("=" * 126)

    for p in [
        WEIGHTS,
        COVERAGE,
        SCANB_FEATURES,
        METABRIC_FEATURES,
        NULL_CONTRACT,
        PILOT_CONTRACT,
    ]:
        require(p)

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    null_contract = json.loads(NULL_CONTRACT.read_text(encoding="utf-8"))
    pilot_contract = json.loads(PILOT_CONTRACT.read_text(encoding="utf-8"))

    n_panels = int(pilot_contract["pilot_panels_per_target_program"])
    min_valid = int(
        pilot_contract["pilot_feasibility_gate"]["minimum_valid_panels"]
    )
    min_frac_le1 = float(
        pilot_contract["panel_validity"][
            "minimum_fraction_assignments_with_stratum_distance_le_1"
        ]
    )
    max_dist_allowed = int(
        pilot_contract["panel_validity"]["maximum_assignment_stratum_distance"]
    )
    base_seed = int(pilot_contract["seed"]["base_seed"])

    weights = pd.read_csv(WEIGHTS, sep="\t", dtype=str)
    coverage = pd.read_csv(COVERAGE, sep="\t")
    features_by_target = {
        "SCANB_GSE96058": pd.read_csv(SCANB_FEATURES, sep="\t"),
        "METABRIC": pd.read_csv(METABRIC_FEATURES, sep="\t"),
    }

    weights["GENE"] = weights["Hugo_Symbol"].map(norm)

    program_genes = {
        pid: grp["GENE"].tolist()
        for pid, grp in weights.groupby("program_id", sort=True)
    }

    panel_rows = []
    target_program_rows = []

    print("\n[1/2] Generating predeclared pilot mappings ...")

    for t_idx, (target, feat) in enumerate(features_by_target.items(), start=1):
        feat = feat.copy()
        feat["Hugo_Symbol"] = feat["Hugo_Symbol"].map(norm)

        if feat["Hugo_Symbol"].duplicated().any():
            raise RuntimeError(f"Duplicate feature-table symbols for {target}")

        feat_indexed = feat.set_index("Hugo_Symbol", drop=False)

        print(f"\n  {target}:")

        for p_idx, (pid, frozen_program_genes) in enumerate(
            sorted(program_genes.items()), start=1
        ):
            cov_row = coverage[
                (coverage["target"] == target)
                & (coverage["program_id"] == pid)
            ]
            if len(cov_row) != 1:
                raise RuntimeError(f"Coverage row missing/duplicated for {target} {pid}")
            primary_assessable = int(cov_row.iloc[0]["primary_assessable"])
            expected_evaluable = int(cov_row.iloc[0]["evaluable_genes"])

            eval_genes = [
                g for g in frozen_program_genes if g in feat_indexed.index
            ]
            if len(eval_genes) != expected_evaluable:
                raise RuntimeError(
                    f"Evaluable-gene replay mismatch for {target} {pid}: "
                    f"{len(eval_genes)} vs {expected_evaluable}"
                )

            slots = feat_indexed.loc[
                eval_genes,
                [
                    "Hugo_Symbol",
                    "source_composite_bin",
                    "target_composite_bin",
                ],
            ].reset_index(drop=True)

            frozen_program_set = set(frozen_program_genes)
            candidates = feat[
                ~feat["Hugo_Symbol"].isin(frozen_program_set)
            ][
                [
                    "Hugo_Symbol",
                    "source_composite_bin",
                    "target_composite_bin",
                ]
            ].copy()

            seed = base_seed + (t_idx * 100000) + (p_idx * 1000)
            rng = np.random.default_rng(seed)

            this_panels = []
            for panel_id in range(1, n_panels + 1):
                result = sample_mapping(
                    slots=slots,
                    candidates=candidates,
                    rng=rng,
                    max_search_distance=max_dist_allowed,
                )
                row = summarize_panel(
                    target=target,
                    program_id=pid,
                    panel_id=panel_id,
                    result=result,
                    min_frac_le1=min_frac_le1,
                    max_dist_allowed=max_dist_allowed,
                    primary_assessable=primary_assessable,
                )
                panel_rows.append(row)
                this_panels.append(row)

            this_df = pd.DataFrame(this_panels)
            n_success = int(this_df["assignment_success"].sum())
            n_valid = int(this_df["panel_valid"].sum())
            feasible = int(
                (primary_assessable == 0) or (n_valid >= min_valid)
            )

            success_df = this_df[this_df["assignment_success"] == 1]

            target_program_rows.append(
                {
                    "target": target,
                    "program_id": pid,
                    "evaluable_genes": len(eval_genes),
                    "candidate_pool_size": len(candidates),
                    "primary_assessable": primary_assessable,
                    "pilot_panels": n_panels,
                    "assignment_success_panels": n_success,
                    "valid_panels": n_valid,
                    "valid_fraction": n_valid / n_panels,
                    "median_fraction_distance_0": (
                        float(success_df["fraction_distance_0"].median())
                        if len(success_df)
                        else np.nan
                    ),
                    "median_fraction_distance_le1": (
                        float(success_df["fraction_distance_le1"].median())
                        if len(success_df)
                        else np.nan
                    ),
                    "worst_maximum_distance": (
                        int(success_df["maximum_distance"].max())
                        if len(success_df)
                        else np.nan
                    ),
                    "median_mean_distance": (
                        float(success_df["mean_distance"].median())
                        if len(success_df)
                        else np.nan
                    ),
                    "pilot_gate_pass": feasible,
                    "seed": seed,
                }
            )

            print(
                f"    {pid:<10} "
                f"m={len(eval_genes):>4} "
                f"valid={n_valid:>3}/{n_panels} "
                f"medianExact="
                f"{(success_df['fraction_distance_0'].median() if len(success_df) else float('nan')):.3f} "
                f"median<=1="
                f"{(success_df['fraction_distance_le1'].median() if len(success_df) else float('nan')):.3f} "
                f"worstMaxDist="
                f"{(int(success_df['maximum_distance'].max()) if len(success_df) else -1)} "
                f"gate={'PASS' if feasible else 'FAIL'} "
                f"{'(not primary-assessable)' if not primary_assessable else ''}"
            )

    panel_df = pd.DataFrame(panel_rows)
    summary_df = pd.DataFrame(target_program_rows)

    print("\n[2/2] Applying frozen pilot feasibility gate ...")
    blocking = summary_df[
        (summary_df["primary_assessable"] == 1)
        & (summary_df["pilot_gate_pass"] == 0)
    ].copy()

    overall_pass = len(blocking) == 0

    panel_out = OUT_DIR / "mapping_null_pilot_panel_quality_v1.tsv"
    summary_out = OUT_DIR / "mapping_null_pilot_summary_v1.tsv"
    panel_df.to_csv(panel_out, sep="\t", index=False)
    summary_df.to_csv(summary_out, sep="\t", index=False)

    result = {
        "result_id": "paper4-tcbb-mapping-null-pilot-feasibility-v1",
        "script_version": SCRIPT_VERSION,
        "status": (
            "PASS_FINAL_NULL_GENERATOR_UNCHANGED"
            if overall_pass
            else "FAIL_REVISE_BEFORE_PRESERVATION"
        ),
        "scientific_guard": {
            "target_outcomes_loaded": False,
            "target_clinical_characteristics_loaded": False,
            "target_gene_gene_correlations_calculated": False,
            "target_pca_calculated": False,
            "target_preservation_statistics_calculated": False,
            "target_expression_values_read_by_this_script": False,
        },
        "pilot_contract": str(PILOT_CONTRACT),
        "null_contract": str(NULL_CONTRACT),
        "overall_pass": overall_pass,
        "blocking_combinations": blocking[
            ["target", "program_id", "valid_panels", "pilot_panels"]
        ].to_dict(orient="records"),
    }

    json_out = OUT_DIR / "mapping_null_pilot_audit_v1.json"
    json_out.write_text(json.dumps(result, indent=2), encoding="utf-8")

    print("\n" + "=" * 126)
    if overall_pass:
        print("03g MAPPING-NULL PILOT FEASIBILITY: PASS")
    else:
        print("03g MAPPING-NULL PILOT FEASIBILITY: FAIL")
    print("=" * 126)

    for target, grp in summary_df.groupby("target", sort=False):
        primary = grp[grp["primary_assessable"] == 1]
        print(f"{target}:")
        print(
            f"  primary-assessable combinations passing gate: "
            f"{int(primary['pilot_gate_pass'].sum())}/{len(primary)}"
        )
        print(
            f"  valid-panel fraction range: "
            f"{primary['valid_fraction'].min():.3f} - "
            f"{primary['valid_fraction'].max():.3f}"
        )

    if overall_pass:
        print()
        print(
            "Frozen 03d mapping generator is practically feasible and may be "
            "retained unchanged for the later 1,000-panel specificity null."
        )
    else:
        print()
        print("BLOCKING combinations:")
        for _, row in blocking.iterrows():
            print(
                f"  {row['target']} {row['program_id']}: "
                f"{int(row['valid_panels'])}/{int(row['pilot_panels'])} valid"
            )
        print(
            "Do NOT calculate preservation. Freeze a newly versioned mapping "
            "generator first."
        )

    print()
    print("No target correlation/PCA/preservation statistic was calculated.")
    print("Outputs:")
    print(f"  {panel_out}")
    print(f"  {summary_out}")
    print(f"  {json_out}")
    print("=" * 126)


if __name__ == "__main__":
    main()
