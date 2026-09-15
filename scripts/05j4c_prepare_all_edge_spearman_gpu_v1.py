#!/usr/bin/env python
from __future__ import annotations

import gc
import hashlib
import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from scipy.stats import rankdata

SCRIPT_VERSION = "05j4c-prepare-all-edge-spearman-gpu-v1-no-cli"

PROJECT_ROOT = Path(r"C:\Users\olegk\Desktop\molecular-transport-audit")
DATA_ROOT = Path(r"D:\paper4_tcbb_data")

SOURCE_PREP = PROJECT_ROOT / "scripts" / "05j4b_prepare_scale_invariant_inputs_gpu_v2.py"
SOURCE_05D = PROJECT_ROOT / "scripts" / "05d_run_controlled_mapping_degradation_gpu_v1.py"
SOURCE_05J = PROJECT_ROOT / "scripts" / "05j_run_structural_corruption_headtohead_gpu_v1.py"

INTERP_CONTRACT = (
    DATA_ROOT
    / "paper4_tcbb_scale_invariant_diagnostic_completion_contract_v1"
    / "all_edge_spearman_interpretation_contract_v1.json"
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
CLASSIFICATION = (
    DATA_ROOT
    / "paper4_tcbb_final_mapping_specificity_null_v2"
    / "primary_pooled_final_classification_v2.tsv"
)
SPECIFICITY_ROOT = DATA_ROOT / "paper4_tcbb_final_mapping_specificity_null_v2"
DIRECT_ROOT = DATA_ROOT / "paper4_tcbb_primary_pooled_direct_preservation_v2"

SCANB_FEATURES = (
    DATA_ROOT / "paper4_tcbb_target_marginal_metrics_v2"
    / "scanb_target_marginal_matching_features_v2.tsv"
)
METABRIC_FEATURES = (
    DATA_ROOT / "paper4_tcbb_target_marginal_metrics_v2"
    / "metabric_target_marginal_matching_features_v2.tsv"
)
SCANB_EVAL = (
    DATA_ROOT / "paper4_tcbb_exact_variation_evaluability_correction_v1"
    / "scanb_target_exact_variation_evaluability_v1.tsv"
)
METABRIC_EVAL = (
    DATA_ROOT / "paper4_tcbb_exact_variation_evaluability_correction_v1"
    / "metabric_target_exact_variation_evaluability_v1.tsv"
)

# These names are read from the historical prep module if present; keep local
# fallback only as a hard failure aid, never to silently choose another file.
WEIGHTS = DATA_ROOT / "paper4_tcbb_tcga_source_module_membership_frozen_v1" / "tcga_source_module_membership_frozen_v1.tsv"

OUT_DIR = DATA_ROOT / "paper4_tcbb_all_edge_spearman_inputs_v1"
COMBO_DIR = OUT_DIR / "per_pair"
PERF_OUT = OUT_DIR / "all_edge_spearman_performance_grid_v1.tsv"
FULL_OUT = OUT_DIR / "all_edge_spearman_fullnull_1000_v1.tsv"
MASTER_OUT = OUT_DIR / "all_edge_spearman_inputs_v1.json"

EXPECTED_PREP_SHA256 = "64580c8b1aa3f2e7382c9c4fc7433ca0e0af87b26c6c56f7f2f9789ea1d7d5e4"
EXPECTED_05J_SHA256 = "e25583122628596eb52d99a7b837031564a5c1764c0e6e371dc2743c7c4ae2e1"

SEP = "=" * 168
TOL = 1e-10

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

def boolish(x: object) -> bool:
    if isinstance(x, bool):
        return x
    return str(x).strip().lower() in {"1", "true", "yes", "y"}

def canon_symbol(x: object) -> str:
    if x is None:
        return ""
    return " ".join(str(x).strip().split()).upper()

def gpu_average_ranks(values: torch.Tensor) -> torch.Tensor:
    """Exact average ranks for equal float64 values; 1-based ranks."""
    require(values.ndim == 1, "gpu_average_ranks expects a 1-D tensor.")
    require(values.dtype == torch.float64, "gpu_average_ranks requires float64.")
    require(bool(torch.all(torch.isfinite(values))), "Non-finite value in rank input.")

    sorted_vals, order = torch.sort(values, stable=True)
    n = int(values.numel())
    new_group = torch.ones(n, dtype=torch.bool, device=values.device)
    if n > 1:
        new_group[1:] = sorted_vals[1:] != sorted_vals[:-1]
    gid = torch.cumsum(new_group.to(torch.int64), dim=0) - 1
    ng = int(gid[-1].item()) + 1

    counts = torch.bincount(gid, minlength=ng)
    ends = torch.cumsum(counts, dim=0)
    starts = ends - counts
    avg = (starts.to(torch.float64) + 1.0 + ends.to(torch.float64)) / 2.0
    ranks_sorted = avg[gid]

    ranks = torch.empty(n, dtype=torch.float64, device=values.device)
    ranks[order] = ranks_sorted
    return ranks

def fixed_rank_vector(values: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
    r = gpu_average_ranks(values)
    rc = r - torch.mean(r)
    rn = torch.linalg.vector_norm(rc)
    require(bool(torch.isfinite(rn)) and float(rn) > 0.0, "Degenerate fixed rank vector.")
    return rc, rn

def spearman_against_fixed_gpu(
    fixed_rank_centered: torch.Tensor,
    fixed_rank_norm: torch.Tensor,
    values: torch.Tensor,
) -> float:
    r = gpu_average_ranks(values)
    rc = r - torch.mean(r)
    rn = torch.linalg.vector_norm(rc)
    if (not bool(torch.isfinite(rn))) or float(rn) <= 0.0:
        return float("nan")
    rho = torch.dot(fixed_rank_centered, rc) / (fixed_rank_norm * rn)
    return float(rho.detach().cpu())

def rank_self_test(device: torch.device) -> None:
    tests = [
        np.array([3.0, 1.0, 2.0, 4.0], dtype=np.float64),
        np.array([1.0, 1.0, 2.0, 4.0, 4.0, 4.0], dtype=np.float64),
        np.array([-2.0, 0.0, -2.0, 7.0, 1.0, 1.0], dtype=np.float64),
    ]
    for a in tests:
        got = gpu_average_ranks(torch.as_tensor(a, dtype=torch.float64, device=device)).cpu().numpy()
        exp = rankdata(a, method="average").astype(np.float64)
        require(np.array_equal(got, exp), f"GPU average-rank self-test failed: got={got}, exp={exp}")

def main() -> None:
    print(SEP)
    print("Paper 4 / TCBB - prepare exact all-edge Spearman diagnostic on frozen mappings")
    print(SEP)
    print(f"Script version: {SCRIPT_VERSION}")
    print()
    print("Execution contract:")
    print("  Spearman_all_edges results previously observed: NO")
    print("  Exact 05j4b mapping generator reused:          YES")
    print("  Saved attempt_id for every mapping reused:      YES")
    print("  All strict upper-triangle edges reused:         YES")
    print("  Only Pearson -> Spearman changed:               YES")
    print("  New mappings generated:                         NO")
    print("  Existing MTA/NetRep/WGCNA results changed:      NO")
    print()

    for p in [SOURCE_PREP, SOURCE_05D, SOURCE_05J, INTERP_CONTRACT, PERFORMANCE, FULL_NULL, CLASSIFICATION]:
        require(p.exists(), f"Missing required file: {p}")

    require(sha256_file(SOURCE_PREP) == EXPECTED_PREP_SHA256, "05j4b source hash drift.")
    require(sha256_file(SOURCE_05J) == EXPECTED_05J_SHA256, "05j source hash drift.")

    c = json.loads(INTERP_CONTRACT.read_text(encoding="utf-8"))
    require(c.get("status") == "FROZEN_BEFORE_ALL_EDGE_SPEARMAN_RESULTS", "Interpretation contract not frozen.")

    prep = load_module(SOURCE_PREP, "paper4_05j4b_frozen")
    exact05d = prep.load_module(SOURCE_05D, "paper4_05d_exact_for_05j4c")
    exact05j = prep.load_module(SOURCE_05J, "paper4_05j_exact_for_05j4c")

    require(torch.cuda.is_available(), "CUDA GPU required for exact all-edge Spearman run.")
    device = torch.device("cuda")
    prop = torch.cuda.get_device_properties(device)
    print(f"CUDA device: {prop.name}; VRAM={prop.total_memory/1024**3:.2f} GB")

    rank_self_test(device)
    print("GPU average-rank tie-handling self-test vs scipy: PASS")

    perf = pd.read_csv(PERFORMANCE, sep="\t", low_memory=False)
    full = pd.read_csv(FULL_NULL, sep="\t", low_memory=False)
    require(len(perf) == 22 * 601, f"Unexpected performance rows: {len(perf)}")
    require(len(full) == 22 * 1000, f"Unexpected full-null rows: {len(full)}")

    classification = pd.read_csv(CLASSIFICATION, sep="\t", low_memory=False)
    classification = classification.loc[classification["primary_assessable"].map(boolish)].copy()
    require(len(classification) == 22, f"Expected 22 assessable pairs, found {len(classification)}")

    # Use the exact historical source loaders and the exact weights path used by 05j4b.
    weights_path = getattr(prep, "WEIGHTS", None)
    require(weights_path is not None and Path(weights_path).exists(), "Historical 05j4b WEIGHTS path unavailable.")
    weights = pd.read_csv(weights_path, sep="\t", dtype=str).fillna("")
    weights["Hugo_Symbol"] = weights["Hugo_Symbol"].map(canon_symbol)

    source = exact05j.load_tcga_source()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    COMBO_DIR.mkdir(parents=True, exist_ok=True)

    all_perf = []
    all_full = []

    for target_name in ["SCANB_GSE96058", "METABRIC"]:
        print("\n" + SEP)
        print(f"TARGET: {target_name}")
        print(SEP)

        if target_name == "SCANB_GSE96058":
            corrected = exact05j.load_eval_symbols(SCANB_EVAL)
            features = exact05j.load_features(SCANB_FEATURES, corrected)
            target = exact05j.load_scanb_target(corrected)
        else:
            corrected = exact05j.load_eval_symbols(METABRIC_EVAL)
            features = exact05j.load_features(METABRIC_FEATURES, corrected)
            target = exact05j.load_metabric_target(corrected)

        print("  computing full target Pearson correlation on GPU ...")
        target_corr = exact05j.corr_gpu(target["x"], device)
        del target["x"]
        torch.cuda.empty_cache()
        gc.collect()

        rows_target = classification.loc[classification["target"] == target_name].copy()

        for combo_i, class_row in enumerate(rows_target.to_dict(orient="records"), start=1):
            program_id = str(class_row["program_id"])
            print(f"\n  [{combo_i:02d}/{len(rows_target):02d}] {program_id}")

            perf_path = COMBO_DIR / f"{target_name}__{program_id}__spearman_all_edges_performance_v1.tsv"
            full_path = COMBO_DIR / f"{target_name}__{program_id}__spearman_all_edges_fullnull_v1.tsv"
            meta_path = COMBO_DIR / f"{target_name}__{program_id}__spearman_all_edges_meta_v1.json"

            if perf_path.exists() and full_path.exists() and meta_path.exists():
                meta = json.loads(meta_path.read_text(encoding="utf-8"))
                if meta.get("status") == "COMPLETE":
                    print("      completed pair found — reusing")
                    all_perf.append(pd.read_csv(perf_path, sep="\t"))
                    all_full.append(pd.read_csv(full_path, sep="\t"))
                    continue

            gene_file = (
                Path(getattr(prep, "DIRECT_ROOT"))
                / "per_program" / target_name
                / f"{program_id}_evaluable_genes_v2.tsv"
            )
            require(gene_file.exists(), f"Missing evaluable-gene file: {gene_file}")
            eg = pd.read_csv(gene_file, sep="\t", dtype=str).fillna("")
            slot_genes = [canon_symbol(x) for x in eg["Hugo_Symbol"]]
            m = len(slot_genes)

            full_program_genes = {
                canon_symbol(x)
                for x in weights.loc[weights["program_id"] == program_id, "Hugo_Symbol"]
            }

            correct_target_idx = np.array(
                [target["gene_index"][g] for g in slot_genes],
                dtype=np.int64,
            )
            source_cols = np.array(
                [source["gene_index"][g] for g in slot_genes],
                dtype=np.int64,
            )
            source_corr = exact05j.corr_gpu(source["x"][:, source_cols], device)

            tri_i_np, tri_j_np = np.triu_indices(m, k=1)
            tri_i_t = torch.as_tensor(tri_i_np, dtype=torch.long, device=device)
            tri_j_t = torch.as_tensor(tri_j_np, dtype=torch.long, device=device)
            source_all_edges = source_corr[tri_i_t, tri_j_t]
            source_rank_centered, source_rank_norm = fixed_rank_vector(source_all_edges)
            source_pearson_centered, source_pearson_norm = exact05j.fixed_pearson_vector(source_all_edges)

            n_edges = int(source_all_edges.numel())
            print(f"      module genes={m:,}; all edges={n_edges:,}")

            def score_attempt(fraction: float, attempt_id: int) -> tuple[float, float]:
                if fraction == 0.0:
                    mapped = correct_target_idx.copy()
                else:
                    mapped, _ = prep.build_mapping(
                        exact05d,
                        target_name,
                        program_id,
                        float(fraction),
                        int(attempt_id),
                        slot_genes,
                        full_program_genes,
                        features,
                        target["gene_index"],
                        correct_target_idx,
                    )
                mapped_t = torch.as_tensor(mapped, dtype=torch.long, device=device)
                sub = target_corr.index_select(0, mapped_t).index_select(1, mapped_t)
                target_all_edges = sub[tri_i_t, tri_j_t]
                rho_s = spearman_against_fixed_gpu(
                    source_rank_centered,
                    source_rank_norm,
                    target_all_edges,
                )
                rho_p = exact05j.pearson_against_fixed(
                    source_pearson_centered,
                    source_pearson_norm,
                    target_all_edges,
                )
                del mapped_t, sub, target_all_edges
                return rho_s, rho_p

            pg = perf.loc[
                (perf["target"] == target_name)
                & (perf["program_id"] == program_id)
            ].copy().sort_values(["fraction", "replicate_id"]).reset_index(drop=True)
            require(len(pg) == 601, f"{target_name}/{program_id}: performance rows={len(pg)}")

            perf_rows = []
            for ri, row in enumerate(pg.to_dict(orient="records"), start=1):
                fraction = float(row["fraction"])
                attempt_id = int(row["attempt_id"])
                rho, pearson_replay = score_attempt(fraction, attempt_id)
                pearson_saved = float(row["netrep_cor_cor"])
                diff = abs(pearson_replay - pearson_saved)
                require(
                    diff <= TOL,
                    (
                        f"{target_name}/{program_id}/f={fraction}/rep={int(row['replicate_id'])}: "
                        f"all-edge Pearson replay mismatch {pearson_replay:.16g} vs "
                        f"{pearson_saved:.16g}; |Δ|={diff:.3e}"
                    ),
                )
                perf_rows.append({
                    "target": target_name,
                    "program_id": program_id,
                    "fraction": fraction,
                    "replicate_id": int(row["replicate_id"]),
                    "attempt_id": attempt_id,
                    "spearman_all_edges": rho,
                    "netrep_cor_cor_replayed_guard": pearson_replay,
                    "netrep_cor_cor_saved_guard": pearson_saved,
                    "netrep_replay_abs_diff": diff,
                })
                if ri % 100 == 0 or ri == len(pg):
                    print(f"      performance {ri:3d}/{len(pg)}")

            ng = full.loc[
                (full["target"] == target_name)
                & (full["program_id"] == program_id)
            ].copy().sort_values("valid_null_id").reset_index(drop=True)
            require(len(ng) == 1000, f"{target_name}/{program_id}: full-null rows={len(ng)}")

            full_rows = []
            for ri, row in enumerate(ng.to_dict(orient="records"), start=1):
                rho, pearson_replay = score_attempt(1.0, int(row["attempt_id"]))
                pearson_saved = float(row["netrep_cor_cor"])
                diff = abs(pearson_replay - pearson_saved)
                require(
                    diff <= TOL,
                    (
                        f"{target_name}/{program_id}/null={int(row['valid_null_id'])}: "
                        f"all-edge Pearson replay mismatch {pearson_replay:.16g} vs "
                        f"{pearson_saved:.16g}; |Δ|={diff:.3e}"
                    ),
                )
                full_rows.append({
                    "target": target_name,
                    "program_id": program_id,
                    "valid_null_id": int(row["valid_null_id"]),
                    "attempt_id": int(row["attempt_id"]),
                    "spearman_all_edges": rho,
                    "netrep_cor_cor_replayed_guard": pearson_replay,
                    "netrep_cor_cor_saved_guard": pearson_saved,
                    "netrep_replay_abs_diff": diff,
                })
                if ri % 100 == 0 or ri == len(ng):
                    print(f"      full-null {ri:4d}/{len(ng)}")

            pdf = pd.DataFrame(perf_rows)
            ndf = pd.DataFrame(full_rows)
            require(np.isfinite(pdf["spearman_all_edges"]).all(), "Non-finite performance Spearman.")
            require(np.isfinite(ndf["spearman_all_edges"]).all(), "Non-finite full-null Spearman.")

            pdf.to_csv(perf_path, sep="\t", index=False)
            ndf.to_csv(full_path, sep="\t", index=False)

            meta = {
                "script_version": SCRIPT_VERSION,
                "status": "COMPLETE",
                "target": target_name,
                "program_id": program_id,
                "module_genes": m,
                "all_edges": n_edges,
                "performance_rows": len(pdf),
                "full_null_rows": len(ndf),
                "max_netrep_pearson_replay_abs_diff": float(max(
                    pdf["netrep_replay_abs_diff"].max(),
                    ndf["netrep_replay_abs_diff"].max(),
                )),
                "edge_universe": "strict_upper_triangle_all_module_edges",
                "only_change_from_netrep_cor_cor": "Pearson_to_Spearman",
                "result_performance_file": str(perf_path),
                "result_fullnull_file": str(full_path),
            }
            meta_path.write_text(json.dumps(meta, indent=2) + "\n", encoding="utf-8")

            all_perf.append(pdf)
            all_full.append(ndf)

            del (
                source_corr,
                tri_i_t,
                tri_j_t,
                source_all_edges,
                source_rank_centered,
                source_rank_norm,
                source_pearson_centered,
                source_pearson_norm,
            )
            torch.cuda.empty_cache()
            gc.collect()

        del target_corr, target, features
        torch.cuda.empty_cache()
        gc.collect()

    perf_all = pd.concat(all_perf, ignore_index=True)
    full_all = pd.concat(all_full, ignore_index=True)
    require(len(perf_all) == 22 * 601, "Combined performance row-count drift.")
    require(len(full_all) == 22 * 1000, "Combined full-null row-count drift.")

    require(not PERF_OUT.exists(), f"Refusing to overwrite: {PERF_OUT}")
    require(not FULL_OUT.exists(), f"Refusing to overwrite: {FULL_OUT}")
    require(not MASTER_OUT.exists(), f"Refusing to overwrite: {MASTER_OUT}")

    perf_all.to_csv(PERF_OUT, sep="\t", index=False)
    full_all.to_csv(FULL_OUT, sep="\t", index=False)

    master = {
        "script_version": SCRIPT_VERSION,
        "status": "ALL_EDGE_SPEARMAN_INPUTS_COMPLETE",
        "target_program_pairs": 22,
        "performance_rows": len(perf_all),
        "full_null_rows": len(full_all),
        "mapping_attempt_ids_reused": True,
        "new_mappings_generated": False,
        "edge_universe": "exact all-edge upper triangle used by NetRep cor.cor",
        "only_scientific_change": "Pearson_to_Spearman",
        "performance_file": str(PERF_OUT),
        "full_null_file": str(FULL_OUT),
    }
    MASTER_OUT.write_text(json.dumps(master, indent=2) + "\n", encoding="utf-8")

    print("\n" + SEP)
    print("05j4c ALL-EDGE SPEARMAN INPUT PREPARATION: COMPLETE")
    print(SEP)
    print(f"Performance rows: {len(perf_all):,}")
    print(f"Full-null rows:   {len(full_all):,}")
    print("New mappings generated: NO")
    print("Only Pearson -> Spearman changed: YES")
    print(f"Performance: {PERF_OUT}")
    print(f"Full null:   {FULL_OUT}")
    print(f"Master:      {MASTER_OUT}")
    print(SEP)

if __name__ == "__main__":
    main()
