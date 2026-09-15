#!/usr/bin/env python
from __future__ import annotations

import gc
import hashlib
import importlib.util
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch

SCRIPT_VERSION = "05j4f-run-backend-matched-kernel-benchmark-v1-no-cli"

PROJECT_ROOT = Path(r"C:\Users\olegk\Desktop\molecular-transport-audit")
DATA_ROOT = Path(r"D:\paper4_tcbb_data")

CONTRACT = (
    DATA_ROOT
    / "paper4_tcbb_backend_matched_kernel_benchmark_contract_v1"
    / "backend_matched_kernel_benchmark_contract_v1.json"
)

SOURCE_05J4B = PROJECT_ROOT / "scripts" / "05j4b_prepare_scale_invariant_inputs_gpu_v2.py"
SOURCE_05J = PROJECT_ROOT / "scripts" / "05j_run_structural_corruption_headtohead_gpu_v1.py"
SOURCE_05D = PROJECT_ROOT / "scripts" / "05d_run_controlled_mapping_degradation_gpu_v1.py"
SOURCE_SPEARMAN = PROJECT_ROOT / "scripts" / "05j4c_prepare_all_edge_spearman_gpu_v1.py"

FULL_NULL = (
    DATA_ROOT
    / "paper4_tcbb_scale_invariant_corruption_inputs_v3"
    / "scale_invariant_fullnull_1000_v3.tsv"
)
ALL_EDGE_FULL_NULL = (
    DATA_ROOT
    / "paper4_tcbb_all_edge_spearman_inputs_v1"
    / "all_edge_spearman_fullnull_1000_v1.tsv"
)
CLASSIFICATION = (
    DATA_ROOT
    / "paper4_tcbb_final_mapping_specificity_null_v2"
    / "primary_pooled_final_classification_v2.tsv"
)
SPECIFICITY_ROOT = DATA_ROOT / "paper4_tcbb_final_mapping_specificity_null_v2"

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

OUT_DIR = DATA_ROOT / "paper4_tcbb_backend_matched_kernel_benchmark_v1"
REPEATS_OUT = OUT_DIR / "backend_matched_kernel_repeat_timings_v1.tsv"
MAPPING_OUT = OUT_DIR / "backend_matched_kernel_by_mapping_v1.tsv"
PAIR_OUT = OUT_DIR / "backend_matched_kernel_by_pair_v1.tsv"
SUMMARY_OUT = OUT_DIR / "backend_matched_kernel_summary_v1.tsv"
MASTER_OUT = OUT_DIR / "backend_matched_kernel_benchmark_v1.json"

SEP = "=" * 168

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

def desc(values) -> dict:
    z = np.asarray(values, dtype=float)
    require(len(z) > 0 and np.isfinite(z).all(), "Cannot summarize empty/non-finite timings.")
    return {
        "mean": float(np.mean(z)),
        "median": float(np.median(z)),
        "q25": float(np.quantile(z, 0.25)),
        "q75": float(np.quantile(z, 0.75)),
        "min": float(np.min(z)),
        "max": float(np.max(z)),
    }

def timed_kernel(fn, device: torch.device) -> dict:
    torch.cuda.synchronize(device)
    baseline_alloc = float(torch.cuda.memory_allocated(device))
    baseline_reserved = float(torch.cuda.memory_reserved(device))
    torch.cuda.reset_peak_memory_stats(device)
    t0 = time.perf_counter()
    value = fn()
    torch.cuda.synchronize(device)
    dt = time.perf_counter() - t0
    peak_alloc = float(torch.cuda.max_memory_allocated(device))
    peak_reserved = float(torch.cuda.max_memory_reserved(device))
    return {
        "value": float(value),
        "seconds": float(dt),
        "gpu_baseline_allocated_bytes": baseline_alloc,
        "gpu_peak_allocated_bytes": peak_alloc,
        "gpu_peak_incremental_allocated_bytes": max(0.0, peak_alloc - baseline_alloc),
        "gpu_baseline_reserved_bytes": baseline_reserved,
        "gpu_peak_reserved_bytes": peak_reserved,
        "gpu_peak_incremental_reserved_bytes": max(0.0, peak_reserved - baseline_reserved),
    }

def main() -> None:
    print(SEP)
    print("Paper 4 / TCBB - backend-matched Spearman kernel benchmark")
    print(SEP)
    print(f"Script version: {SCRIPT_VERSION}")
    print()

    required = [
        CONTRACT,
        SOURCE_05J4B,
        SOURCE_05J,
        SOURCE_05D,
        SOURCE_SPEARMAN,
        FULL_NULL,
        ALL_EDGE_FULL_NULL,
        CLASSIFICATION,
    ]
    for p in required:
        require(p.exists(), f"Missing required file: {p}")

    contract = json.loads(CONTRACT.read_text(encoding="utf-8"))
    require(
        contract.get("status") == "FROZEN_BEFORE_BACKEND_MATCHED_KERNEL_RESULTS",
        "Benchmark contract is not frozen/authoritative.",
    )

    for path_s, expected in contract["input_sha256"].items():
        p = Path(path_s)
        require(p.exists(), f"Frozen input missing: {p}")
        actual = sha256_file(p)
        require(actual == expected, f"Frozen input hash drift: {p}")

    require(torch.cuda.is_available(), "Frozen benchmark requires CUDA.")
    device = torch.device("cuda")
    prop = torch.cuda.get_device_properties(device)
    hw = contract["hardware_lock"]
    require(prop.name == hw["gpu_name"], f"GPU name drift: {prop.name} vs {hw['gpu_name']}")
    require(int(prop.total_memory) == int(hw["gpu_total_memory_bytes"]), "GPU memory drift.")

    valid_ids = [int(x) for x in contract["fixed_workload"]["valid_null_ids"]]
    warmups = int(contract["fixed_workload"]["warmups_per_mapping"])
    repeats = int(contract["fixed_workload"]["timed_repeats_per_mapping_per_method"])
    replay_tol = float(contract["score_replay_guards"]["tolerance"])

    prep = load_module(SOURCE_05J4B, "paper4_05j4b_for_05j4f")
    exact05d = prep.load_module(SOURCE_05D, "paper4_05d_for_05j4f")
    exact05j = prep.load_module(SOURCE_05J, "paper4_05j_for_05j4f")
    spearmod = load_module(SOURCE_SPEARMAN, "paper4_05j4c_spearman_for_05j4f")

    spearmod.rank_self_test(device)

    print(f"CUDA device: {prop.name}; VRAM={prop.total_memory/1024**3:.2f} GB")
    print("GPU average-rank tie-handling self-test: PASS")
    print(f"Fixed valid_null_id: {valid_ids}")
    print(f"Warmups/mapping: {warmups}; timed repeats/mapping/method: {repeats}")
    print()

    full = pd.read_csv(FULL_NULL, sep="\t", low_memory=False)
    all_saved = pd.read_csv(ALL_EDGE_FULL_NULL, sep="\t", low_memory=False)

    classification = pd.read_csv(CLASSIFICATION, sep="\t", low_memory=False)
    classification = classification.loc[classification["primary_assessable"].map(boolish)].copy()
    require(len(classification) == 22, f"Expected 22 assessable pairs, found {len(classification)}")

    weights_path = Path(getattr(prep, "WEIGHTS"))
    require(weights_path.exists(), f"Historical weights path missing: {weights_path}")
    weights = pd.read_csv(weights_path, sep="\t", dtype=str).fillna("")
    weights["Hugo_Symbol"] = weights["Hugo_Symbol"].map(canon_symbol)

    source = exact05j.load_tcga_source()

    repeat_rows = []
    mapping_rows = []
    pair_rows = []

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

        print("  computing shared full target correlation matrix (NOT TIMED) ...")
        target_corr = exact05j.corr_gpu(target["x"], device)
        del target["x"]
        torch.cuda.empty_cache()
        gc.collect()

        rows_target = classification.loc[classification["target"] == target_name].copy()

        for combo_i, class_row in enumerate(rows_target.to_dict(orient="records"), start=1):
            program_id = str(class_row["program_id"])
            print(f"\n  [{combo_i:02d}/{len(rows_target):02d}] {program_id}")

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
                for x in weights.loc[
                    weights["program_id"] == program_id,
                    "Hugo_Symbol",
                ]
            }

            correct_target_idx = np.array(
                [target["gene_index"][g] for g in slot_genes],
                dtype=np.int64,
            )

            null_path = (
                SPECIFICITY_ROOT
                / "per_program"
                / target_name
                / f"{program_id}_mapping_specificity_null_v2.npz"
            )
            require(null_path.exists(), f"Missing frozen subset NPZ: {null_path}")
            nz = np.load(null_path)

            ii = np.asarray(nz["edge_slot_i_0based"], dtype=np.int64)
            jj = np.asarray(nz["edge_slot_j_0based"], dtype=np.int64)
            source_subset_np = np.asarray(nz["source_edge_subset"], dtype=np.float64)
            subset_edges = int(len(source_subset_np))
            require(len(ii) == len(jj) == subset_edges, "Frozen subset length mismatch.")

            ii_t = torch.as_tensor(ii, dtype=torch.long, device=device)
            jj_t = torch.as_tensor(jj, dtype=torch.long, device=device)
            source_subset_t = torch.as_tensor(source_subset_np, dtype=torch.float64, device=device)
            subset_rank_centered, subset_rank_norm = spearmod.fixed_rank_vector(source_subset_t)

            source_cols = np.array(
                [source["gene_index"][g] for g in slot_genes],
                dtype=np.int64,
            )
            source_corr = exact05j.corr_gpu(source["x"][:, source_cols], device)

            tri_i_np, tri_j_np = np.triu_indices(m, k=1)
            tri_i_t = torch.as_tensor(tri_i_np, dtype=torch.long, device=device)
            tri_j_t = torch.as_tensor(tri_j_np, dtype=torch.long, device=device)
            source_all_edges_t = source_corr[tri_i_t, tri_j_t]
            all_edges = int(source_all_edges_t.numel())
            all_rank_centered, all_rank_norm = spearmod.fixed_rank_vector(source_all_edges_t)

            require(
                all_edges == m * (m - 1) // 2,
                f"{target_name}/{program_id}: all-edge count identity failed.",
            )
            require(
                subset_edges <= all_edges,
                f"{target_name}/{program_id}: subset edge count exceeds all-edge count.",
            )

            frozen = full.loc[
                (full["target"] == target_name)
                & (full["program_id"] == program_id)
                & (full["valid_null_id"].astype(int).isin(valid_ids))
            ].copy().sort_values("valid_null_id").reset_index(drop=True)
            require(
                frozen["valid_null_id"].astype(int).tolist() == valid_ids,
                f"{target_name}/{program_id}: frozen valid_null_id workload mismatch.",
            )

            saved_all = all_saved.loc[
                (all_saved["target"] == target_name)
                & (all_saved["program_id"] == program_id)
                & (all_saved["valid_null_id"].astype(int).isin(valid_ids))
            ].copy().sort_values("valid_null_id").reset_index(drop=True)
            require(
                saved_all["valid_null_id"].astype(int).tolist() == valid_ids,
                f"{target_name}/{program_id}: all-edge saved workload mismatch.",
            )

            # attempt_id identity must agree across the two saved inputs.
            require(
                frozen["attempt_id"].astype(int).tolist()
                == saved_all["attempt_id"].astype(int).tolist(),
                f"{target_name}/{program_id}: saved attempt_id disagreement.",
            )

            pair_subset_times = []
            pair_all_times = []
            pair_subset_mem = []
            pair_all_mem = []

            for map_i, row in enumerate(frozen.to_dict(orient="records"), start=1):
                valid_null_id = int(row["valid_null_id"])
                attempt_id = int(row["attempt_id"])

                mapped, _ = prep.build_mapping(
                    exact05d,
                    target_name,
                    program_id,
                    1.0,
                    attempt_id,
                    slot_genes,
                    full_program_genes,
                    features,
                    target["gene_index"],
                    correct_target_idx,
                )
                mapped_t = torch.as_tensor(mapped, dtype=torch.long, device=device)

                # CRITICAL boundary: mapped target submatrix is constructed before timing
                # and retained as the common baseline allocation for both methods.
                sub = target_corr.index_select(0, mapped_t).index_select(1, mapped_t)
                del mapped_t
                torch.cuda.synchronize(device)

                def subset_kernel():
                    vals = sub[ii_t, jj_t]
                    return spearmod.spearman_against_fixed_gpu(
                        subset_rank_centered,
                        subset_rank_norm,
                        vals,
                    )

                def all_kernel():
                    vals = sub[tri_i_t, tri_j_t]
                    return spearmod.spearman_against_fixed_gpu(
                        all_rank_centered,
                        all_rank_norm,
                        vals,
                    )

                # Scientific score replay before timing.
                subset_score = float(subset_kernel())
                all_score = float(all_kernel())

                saved_subset_score = float(row["mta_rho_edge"])
                all_row = saved_all.loc[
                    saved_all["valid_null_id"].astype(int) == valid_null_id
                ]
                require(len(all_row) == 1, "Saved all-edge score row missing/duplicated.")
                saved_all_score = float(all_row.iloc[0]["spearman_all_edges"])

                subset_diff = abs(subset_score - saved_subset_score)
                all_diff = abs(all_score - saved_all_score)
                require(
                    subset_diff <= replay_tol,
                    (
                        f"{target_name}/{program_id}/null={valid_null_id}: "
                        f"subset GPU replay mismatch |Δ|={subset_diff:.3e}"
                    ),
                )
                require(
                    all_diff <= replay_tol,
                    (
                        f"{target_name}/{program_id}/null={valid_null_id}: "
                        f"all-edge GPU replay mismatch |Δ|={all_diff:.3e}"
                    ),
                )

                # Warmups are outside timed records.
                for _ in range(warmups):
                    subset_kernel()
                    all_kernel()
                    torch.cuda.synchronize(device)

                map_times = {"subset": [], "all": []}
                map_mem = {"subset": [], "all": []}

                for rep in range(repeats):
                    order = ["subset", "all"] if rep % 2 == 0 else ["all", "subset"]
                    for method in order:
                        mm = timed_kernel(
                            subset_kernel if method == "subset" else all_kernel,
                            device,
                        )
                        map_times[method].append(mm["seconds"])
                        map_mem[method].append(mm["gpu_peak_incremental_allocated_bytes"])
                        repeat_rows.append({
                            "target": target_name,
                            "program_id": program_id,
                            "module_genes": m,
                            "valid_null_id": valid_null_id,
                            "attempt_id": attempt_id,
                            "method": method,
                            "repeat": rep + 1,
                            "edge_count": subset_edges if method == "subset" else all_edges,
                            "seconds": mm["seconds"],
                            "gpu_peak_incremental_allocated_bytes": mm[
                                "gpu_peak_incremental_allocated_bytes"
                            ],
                            "gpu_peak_incremental_reserved_bytes": mm[
                                "gpu_peak_incremental_reserved_bytes"
                            ],
                        })

                ss = desc(map_times["subset"])
                aa = desc(map_times["all"])
                ratio = aa["median"] / ss["median"]

                pair_subset_times.extend(map_times["subset"])
                pair_all_times.extend(map_times["all"])
                pair_subset_mem.extend(map_mem["subset"])
                pair_all_mem.extend(map_mem["all"])

                mapping_rows.append({
                    "target": target_name,
                    "program_id": program_id,
                    "module_genes": m,
                    "valid_null_id": valid_null_id,
                    "attempt_id": attempt_id,
                    "subset_edges": subset_edges,
                    "all_edges": all_edges,
                    "all_over_subset_edge_count_ratio": all_edges / subset_edges,
                    "subset_seconds_median": ss["median"],
                    "subset_seconds_q25": ss["q25"],
                    "subset_seconds_q75": ss["q75"],
                    "all_seconds_median": aa["median"],
                    "all_seconds_q25": aa["q25"],
                    "all_seconds_q75": aa["q75"],
                    "all_over_subset_time_ratio": ratio,
                    "subset_score_replay_abs_diff": subset_diff,
                    "all_score_replay_abs_diff": all_diff,
                    "subset_gpu_peak_incremental_allocated_bytes_max": max(map_mem["subset"]),
                    "all_gpu_peak_incremental_allocated_bytes_max": max(map_mem["all"]),
                })

                del sub
                torch.cuda.empty_cache()

            ps = desc(pair_subset_times)
            pa = desc(pair_all_times)
            pair_ratio = pa["median"] / ps["median"]

            pair_rows.append({
                "target": target_name,
                "program_id": program_id,
                "module_genes": m,
                "subset_edges": subset_edges,
                "all_edges": all_edges,
                "subset_edge_fraction": subset_edges / all_edges,
                "edge_count_reduction_fraction": 1.0 - subset_edges / all_edges,
                "all_over_subset_edge_count_ratio": all_edges / subset_edges,
                "subset_raw_float64_edge_bytes": subset_edges * 8,
                "all_raw_float64_edge_bytes": all_edges * 8,
                "n_fixed_mappings": len(valid_ids),
                "timed_repeats_per_mapping": repeats,
                "n_timed_observations_per_method": len(pair_subset_times),
                "subset_kernel_seconds_median": ps["median"],
                "subset_kernel_seconds_q25": ps["q25"],
                "subset_kernel_seconds_q75": ps["q75"],
                "subset_kernel_seconds_min": ps["min"],
                "subset_kernel_seconds_max": ps["max"],
                "all_kernel_seconds_median": pa["median"],
                "all_kernel_seconds_q25": pa["q25"],
                "all_kernel_seconds_q75": pa["q75"],
                "all_kernel_seconds_min": pa["min"],
                "all_kernel_seconds_max": pa["max"],
                "all_over_subset_kernel_time_ratio": pair_ratio,
                "subset_gpu_peak_incremental_allocated_bytes_max": max(pair_subset_mem),
                "all_gpu_peak_incremental_allocated_bytes_max": max(pair_all_mem),
                "backend": "identical_05j4c_gpu_rank_spearman",
                "mapped_submatrix_timing_excluded": True,
            })

            print(
                f"      genes={m:,}; edges={subset_edges:,}/{all_edges:,} "
                f"({all_edges/subset_edges:.1f}x); "
                f"backend-matched kernel all/subset={pair_ratio:.2f}x"
            )

            del (
                source_corr,
                source_subset_t,
                subset_rank_centered,
                subset_rank_norm,
                source_all_edges_t,
                all_rank_centered,
                all_rank_norm,
                ii_t,
                jj_t,
                tri_i_t,
                tri_j_t,
            )
            torch.cuda.empty_cache()
            gc.collect()

        del target_corr, target, features
        torch.cuda.empty_cache()
        gc.collect()

    repeats_df = pd.DataFrame(repeat_rows)
    mapping_df = pd.DataFrame(mapping_rows)
    pair_df = pd.DataFrame(pair_rows)

    require(len(pair_df) == 22, f"Expected 22 pair rows, found {len(pair_df)}.")
    require(
        len(mapping_df) == 22 * len(valid_ids),
        f"Expected {22 * len(valid_ids)} mapping rows, found {len(mapping_df)}.",
    )
    require(
        len(repeats_df) == 22 * len(valid_ids) * repeats * 2,
        "Repeat timing row-count drift.",
    )

    quantities = [
        "subset_edge_fraction",
        "edge_count_reduction_fraction",
        "all_over_subset_edge_count_ratio",
        "all_over_subset_kernel_time_ratio",
        "subset_kernel_seconds_median",
        "all_kernel_seconds_median",
        "subset_gpu_peak_incremental_allocated_bytes_max",
        "all_gpu_peak_incremental_allocated_bytes_max",
    ]
    summary_rows = []
    for col in quantities:
        s = desc(pair_df[col].to_numpy(float))
        summary_rows.append({
            "quantity": col,
            "n_pairs": len(pair_df),
            **s,
        })

    # Descriptive scaling association only; no p-value.
    x = np.log10(pair_df["all_over_subset_edge_count_ratio"].to_numpy(float))
    y = np.log10(pair_df["all_over_subset_kernel_time_ratio"].to_numpy(float))
    if np.std(x) > 0 and np.std(y) > 0:
        scaling_corr = float(np.corrcoef(x, y)[0, 1])
    else:
        scaling_corr = np.nan

    summary_df = pd.DataFrame(summary_rows)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    for p in [REPEATS_OUT, MAPPING_OUT, PAIR_OUT, SUMMARY_OUT, MASTER_OUT]:
        require(not p.exists(), f"Refusing to overwrite benchmark output: {p}")

    repeats_df.to_csv(REPEATS_OUT, sep="\t", index=False)
    mapping_df.to_csv(MAPPING_OUT, sep="\t", index=False)
    pair_df.to_csv(PAIR_OUT, sep="\t", index=False)
    summary_df.to_csv(SUMMARY_OUT, sep="\t", index=False)

    master = {
        "script_version": SCRIPT_VERSION,
        "status": "BACKEND_MATCHED_KERNEL_BENCHMARK_COMPLETE",
        "scientific_results_changed": False,
        "new_inferential_test": False,
        "p_values_computed": False,
        "same_gpu_backend_for_both_methods": True,
        "mapped_submatrix_timing_excluded": True,
        "score_replay_guards": "PASS",
        "fixed_valid_null_ids": valid_ids,
        "descriptive_log10_edge_ratio_vs_log10_time_ratio_pearson": scaling_corr,
        "interpretation_guard": (
            "kernel-level timing only; does not replace 05j4e current-implementation "
            "or whole-workflow result"
        ),
        "outputs": {
            "repeat_timings": str(REPEATS_OUT),
            "by_mapping": str(MAPPING_OUT),
            "by_pair": str(PAIR_OUT),
            "summary": str(SUMMARY_OUT),
        },
    }
    MASTER_OUT.write_text(json.dumps(master, indent=2) + "\n", encoding="utf-8")

    print("\n" + SEP)
    print("05j4f BACKEND-MATCHED KERNEL BENCHMARK: COMPLETE")
    print(SEP)
    print("\nOverall descriptive summary:")
    print(summary_df.to_string(index=False))
    print()
    print(
        "Descriptive corr(log10 edge-count ratio, log10 kernel-time ratio): "
        f"{scaling_corr:.4f}"
    )
    print()
    print("Guards:")
    print("  Identical GPU Spearman backend:       YES")
    print("  Mapped submatrix constructed pre-time:YES")
    print("  Frozen subset score replay:           PASS")
    print("  Frozen all-edge score replay:         PASS")
    print("  New scientific inference:             NO")
    print("  Whole-workflow speedup claim:         NO")
    print()
    print(f"Repeats: {REPEATS_OUT}")
    print(f"Mapping: {MAPPING_OUT}")
    print(f"Pairs:   {PAIR_OUT}")
    print(f"Summary: {SUMMARY_OUT}")
    print(f"Master:  {MASTER_OUT}")
    print(SEP)

if __name__ == "__main__":
    main()
