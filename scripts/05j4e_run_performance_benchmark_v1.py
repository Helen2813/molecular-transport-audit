#!/usr/bin/env python
from __future__ import annotations

import gc
import hashlib
import importlib.util
import json
import math
import statistics
import threading
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch

try:
    import psutil  # optional; CPU RSS benchmark degrades gracefully if absent
except Exception:
    psutil = None

SCRIPT_VERSION = "05j4e-run-performance-benchmark-v1-no-cli"

PROJECT_ROOT = Path(r"C:\Users\olegk\Desktop\molecular-transport-audit")
DATA_ROOT = Path(r"D:\paper4_tcbb_data")

CONTRACT = (
    DATA_ROOT
    / "paper4_tcbb_performance_benchmark_contract_v1"
    / "performance_benchmark_contract_v1.json"
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

OUT_DIR = DATA_ROOT / "paper4_tcbb_performance_benchmark_v1"
PAIR_OUT = OUT_DIR / "performance_benchmark_by_pair_v1.tsv"
REPEATS_OUT = OUT_DIR / "performance_benchmark_repeat_timings_v1.tsv"
TARGET_OUT = OUT_DIR / "performance_benchmark_shared_target_preprocessing_v1.tsv"
SUMMARY_OUT = OUT_DIR / "performance_benchmark_summary_v1.tsv"
MASTER_OUT = OUT_DIR / "performance_benchmark_v1.json"

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

def ssum(x):
    z = np.asarray(x, dtype=float)
    return {
        "median": float(np.median(z)),
        "q25": float(np.quantile(z, 0.25)),
        "q75": float(np.quantile(z, 0.75)),
        "min": float(np.min(z)),
        "max": float(np.max(z)),
        "mean": float(np.mean(z)),
    }

class RSSSampler:
    def __init__(self, interval_s: float = 0.005):
        self.interval_s = interval_s
        self.available = psutil is not None
        self._stop = threading.Event()
        self._thread = None
        self.baseline = np.nan
        self.peak = np.nan

    def __enter__(self):
        if not self.available:
            return self
        proc = psutil.Process()
        self.baseline = float(proc.memory_info().rss)
        self.peak = self.baseline

        def run():
            while not self._stop.is_set():
                try:
                    self.peak = max(self.peak, float(proc.memory_info().rss))
                except Exception:
                    pass
                time.sleep(self.interval_s)

        self._thread = threading.Thread(target=run, daemon=True)
        self._thread.start()
        return self

    def __exit__(self, exc_type, exc, tb):
        if self.available:
            self._stop.set()
            if self._thread is not None:
                self._thread.join(timeout=1.0)

def timed_call(fn, device: torch.device):
    torch.cuda.synchronize(device)
    gc.collect()
    torch.cuda.empty_cache()
    baseline_alloc = float(torch.cuda.memory_allocated(device))
    baseline_reserved = float(torch.cuda.memory_reserved(device))
    torch.cuda.reset_peak_memory_stats(device)

    with RSSSampler() as rss:
        t0 = time.perf_counter()
        out = fn()
        torch.cuda.synchronize(device)
        dt = time.perf_counter() - t0

    peak_alloc = float(torch.cuda.max_memory_allocated(device))
    peak_reserved = float(torch.cuda.max_memory_reserved(device))
    return {
        "result": out,
        "seconds": dt,
        "gpu_baseline_allocated_bytes": baseline_alloc,
        "gpu_peak_allocated_bytes": peak_alloc,
        "gpu_peak_incremental_allocated_bytes": max(0.0, peak_alloc - baseline_alloc),
        "gpu_baseline_reserved_bytes": baseline_reserved,
        "gpu_peak_reserved_bytes": peak_reserved,
        "gpu_peak_incremental_reserved_bytes": max(0.0, peak_reserved - baseline_reserved),
        "cpu_rss_baseline_bytes": float(rss.baseline) if rss.available else np.nan,
        "cpu_rss_peak_bytes": float(rss.peak) if rss.available else np.nan,
        "cpu_rss_peak_incremental_bytes": (
            max(0.0, float(rss.peak) - float(rss.baseline))
            if rss.available else np.nan
        ),
    }

def main() -> None:
    print(SEP)
    print("Paper 4 / TCBB - computational performance benchmark")
    print(SEP)
    print(f"Script version: {SCRIPT_VERSION}")
    print()

    required = [
        CONTRACT, SOURCE_05J4B, SOURCE_05J, SOURCE_05D, SOURCE_SPEARMAN,
        FULL_NULL, CLASSIFICATION,
    ]
    for p in required:
        require(p.exists(), f"Missing required file: {p}")

    contract = json.loads(CONTRACT.read_text(encoding="utf-8"))
    require(
        contract.get("status") == "FROZEN_BEFORE_PERFORMANCE_BENCHMARK_RESULTS",
        "Benchmark contract is not frozen/authoritative.",
    )

    # Source/input immutability.
    for path_s, expected in contract["input_sha256"].items():
        p = Path(path_s)
        require(p.exists(), f"Frozen contract input disappeared: {p}")
        actual = sha256_file(p)
        require(actual == expected, f"Frozen input hash drift: {p}")

    require(torch.cuda.is_available(), "Frozen benchmark requires CUDA.")
    device = torch.device("cuda")
    prop = torch.cuda.get_device_properties(device)
    hw = contract["hardware_lock"]
    require(prop.name == hw["gpu_name"], f"GPU name drift: {prop.name} vs {hw['gpu_name']}")
    require(int(prop.total_memory) == int(hw["gpu_total_memory_bytes"]), "GPU VRAM drift.")

    bench_ids = list(map(int, contract["fixed_benchmark_workload"]["valid_null_ids_per_pair"]))
    nmap = int(contract["fixed_benchmark_workload"]["mappings_per_pair"])
    sw = int(contract["measurement_tiers"]["tier_3_per_mapping_scoring"]["warmups"])
    sr = int(contract["measurement_tiers"]["tier_3_per_mapping_scoring"]["timed_repeats"])
    kw = int(contract["measurement_tiers"]["tier_2_concordance_kernel"]["warmups"])
    kr = int(contract["measurement_tiers"]["tier_2_concordance_kernel"]["timed_repeats"])
    tw = int(contract["measurement_tiers"]["tier_1_shared_target_preprocessing"]["warmups"])
    tr = int(contract["measurement_tiers"]["tier_1_shared_target_preprocessing"]["timed_repeats"])
    require(len(bench_ids) == nmap, "Frozen mapping-ID workload length drift.")

    prep = load_module(SOURCE_05J4B, "paper4_05j4b_for_perf")
    exact05d = prep.load_module(SOURCE_05D, "paper4_05d_for_perf")
    exact05j = prep.load_module(SOURCE_05J, "paper4_05j_for_perf")
    spearmod = load_module(SOURCE_SPEARMAN, "paper4_05j4c_spearman_for_perf")

    print(f"CUDA device: {prop.name}; VRAM={prop.total_memory/1024**3:.2f} GB")
    spearmod.rank_self_test(device)
    print("GPU rank tie-handling self-test: PASS")
    print(f"Fixed mappings/pair: {nmap}; scoring repeats: {sr}; kernel repeats: {kr}")

    full = pd.read_csv(FULL_NULL, sep="\t", low_memory=False)
    classification = pd.read_csv(CLASSIFICATION, sep="\t", low_memory=False)
    classification = classification.loc[classification["primary_assessable"].map(boolish)].copy()
    require(len(classification) == 22, f"Expected 22 assessable pairs, found {len(classification)}")

    weights_path = Path(getattr(prep, "WEIGHTS"))
    require(weights_path.exists(), f"Historical weights path missing: {weights_path}")
    weights = pd.read_csv(weights_path, sep="\t", dtype=str).fillna("")
    weights["Hugo_Symbol"] = weights["Hugo_Symbol"].map(canon_symbol)
    source = exact05j.load_tcga_source()

    target_rows = []
    pair_rows = []
    repeat_rows = []

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

        target_x = target["x"]
        n_samples, n_genes = target_x.shape

        # Tier 1: shared target-correlation preprocessing, measured independently.
        for _ in range(tw):
            tmp = exact05j.corr_gpu(target_x, device)
            torch.cuda.synchronize(device)
            del tmp
            torch.cuda.empty_cache()

        corr_measurements = []
        target_corr = None
        for r in range(tr):
            m = timed_call(lambda: exact05j.corr_gpu(target_x, device), device)
            corr_measurements.append(m)
            if target_corr is not None:
                del target_corr
                torch.cuda.empty_cache()
            target_corr = m["result"]
            m["result"] = None
            print(f"  target-corr timed repeat {r+1}/{tr}: {m['seconds']:.4f}s")

        corr_times = [x["seconds"] for x in corr_measurements]
        corr_sum = ssum(corr_times)
        target_rows.append({
            "target": target_name,
            "n_samples": n_samples,
            "n_corrected_genes": n_genes,
            "target_corr_warmups": tw,
            "target_corr_repeats": tr,
            "target_corr_seconds_mean": corr_sum["mean"],
            "target_corr_seconds_median": corr_sum["median"],
            "target_corr_seconds_q25": corr_sum["q25"],
            "target_corr_seconds_q75": corr_sum["q75"],
            "target_corr_seconds_min": corr_sum["min"],
            "target_corr_seconds_max": corr_sum["max"],
            "target_corr_gpu_peak_incremental_allocated_bytes_max": max(
                x["gpu_peak_incremental_allocated_bytes"] for x in corr_measurements
            ),
            "target_corr_cpu_rss_peak_incremental_bytes_max": np.nanmax(
                [x["cpu_rss_peak_incremental_bytes"] for x in corr_measurements]
            ) if psutil is not None else np.nan,
        })

        del target["x"]

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
            module_genes = len(slot_genes)

            full_program_genes = {
                canon_symbol(x)
                for x in weights.loc[weights["program_id"] == program_id, "Hugo_Symbol"]
            }
            correct_target_idx = np.array(
                [target["gene_index"][g] for g in slot_genes],
                dtype=np.int64,
            )

            null_path = (
                SPECIFICITY_ROOT / "per_program" / target_name
                / f"{program_id}_mapping_specificity_null_v2.npz"
            )
            require(null_path.exists(), f"Missing frozen subset NPZ: {null_path}")
            nz = np.load(null_path)
            ii = np.asarray(nz["edge_slot_i_0based"], dtype=np.int64)
            jj = np.asarray(nz["edge_slot_j_0based"], dtype=np.int64)
            source_subset = np.asarray(nz["source_edge_subset"], dtype=np.float64)
            require(len(ii) == len(jj) == len(source_subset), "Frozen subset length mismatch.")
            subset_edges = int(len(ii))
            source_rank, source_rank_norm = exact05d.centered_rank(source_subset)
            ii_t = torch.as_tensor(ii, dtype=torch.long, device=device)
            jj_t = torch.as_tensor(jj, dtype=torch.long, device=device)

            source_cols = np.array(
                [source["gene_index"][g] for g in slot_genes],
                dtype=np.int64,
            )
            source_corr = exact05j.corr_gpu(source["x"][:, source_cols], device)
            tri_i_np, tri_j_np = np.triu_indices(module_genes, k=1)
            tri_i_t = torch.as_tensor(tri_i_np, dtype=torch.long, device=device)
            tri_j_t = torch.as_tensor(tri_j_np, dtype=torch.long, device=device)
            source_all_edges = source_corr[tri_i_t, tri_j_t]
            source_all_rank_centered, source_all_rank_norm = spearmod.fixed_rank_vector(source_all_edges)
            all_edges = int(source_all_edges.numel())
            require(all_edges == module_genes * (module_genes - 1) // 2, "All-edge count drift.")

            # Fixed workload from already-generated full-null mappings.
            ng = full.loc[
                (full["target"] == target_name)
                & (full["program_id"] == program_id)
                & (full["valid_null_id"].isin(bench_ids))
            ].sort_values("valid_null_id").reset_index(drop=True)
            require(
                ng["valid_null_id"].astype(int).tolist() == bench_ids,
                f"{target_name}/{program_id}: fixed valid_null_id workload mismatch.",
            )

            # Common mapping reconstruction measured once; mappings are then reused
            # outside all method-specific timing regions.
            def reconstruct_all():
                mapped_list = []
                for row in ng.to_dict(orient="records"):
                    mapped, _ = prep.build_mapping(
                        exact05d,
                        target_name,
                        program_id,
                        1.0,
                        int(row["attempt_id"]),
                        slot_genes,
                        full_program_genes,
                        features,
                        target["gene_index"],
                        correct_target_idx,
                    )
                    mapped_list.append(mapped)
                return mapped_list

            map_measure = timed_call(reconstruct_all, device)
            mapped_list = map_measure["result"]
            require(len(mapped_list) == nmap, "Mapped workload reconstruction count drift.")
            mapping_seconds = float(map_measure["seconds"])

            def score_subset_batch():
                acc = 0.0
                for mapped in mapped_list:
                    mapped_t = torch.as_tensor(mapped, dtype=torch.long, device=device)
                    sub = target_corr.index_select(0, mapped_t).index_select(1, mapped_t)
                    vals = sub[ii_t, jj_t].detach().cpu().numpy()
                    rho = exact05d.spearman_against_fixed(source_rank, source_rank_norm, vals)
                    acc += float(rho)
                    del mapped_t, sub, vals
                return acc

            def score_all_batch():
                acc = 0.0
                for mapped in mapped_list:
                    mapped_t = torch.as_tensor(mapped, dtype=torch.long, device=device)
                    sub = target_corr.index_select(0, mapped_t).index_select(1, mapped_t)
                    vals = sub[tri_i_t, tri_j_t]
                    rho = spearmod.spearman_against_fixed_gpu(
                        source_all_rank_centered,
                        source_all_rank_norm,
                        vals,
                    )
                    acc += float(rho)
                    del mapped_t, sub, vals
                return acc

            # Kernel-only representative mapping: submatrix preconstructed.
            rep_mapped = mapped_list[0]
            rep_t = torch.as_tensor(rep_mapped, dtype=torch.long, device=device)
            rep_sub = target_corr.index_select(0, rep_t).index_select(1, rep_t)
            del rep_t

            def kernel_subset():
                vals = rep_sub[ii_t, jj_t].detach().cpu().numpy()
                rho = exact05d.spearman_against_fixed(source_rank, source_rank_norm, vals)
                del vals
                return float(rho)

            def kernel_all():
                vals = rep_sub[tri_i_t, tri_j_t]
                rho = spearmod.spearman_against_fixed_gpu(
                    source_all_rank_centered,
                    source_all_rank_norm,
                    vals,
                )
                del vals
                return float(rho)

            # Warmups, balanced across methods.
            for _ in range(sw):
                score_subset_batch()
                torch.cuda.synchronize(device)
                score_all_batch()
                torch.cuda.synchronize(device)

            score_meas = {"subset": [], "all": []}
            for r in range(sr):
                order = ["subset", "all"] if r % 2 == 0 else ["all", "subset"]
                for method in order:
                    fn = score_subset_batch if method == "subset" else score_all_batch
                    mm = timed_call(fn, device)
                    score_meas[method].append(mm)
                    repeat_rows.append({
                        "target": target_name,
                        "program_id": program_id,
                        "tier": "per_mapping_scoring_batch",
                        "method": method,
                        "repeat": r + 1,
                        "n_mappings": nmap,
                        "seconds": mm["seconds"],
                        "seconds_per_mapping": mm["seconds"] / nmap,
                        "gpu_peak_incremental_allocated_bytes": mm["gpu_peak_incremental_allocated_bytes"],
                        "gpu_peak_incremental_reserved_bytes": mm["gpu_peak_incremental_reserved_bytes"],
                        "cpu_rss_peak_incremental_bytes": mm["cpu_rss_peak_incremental_bytes"],
                    })

            for _ in range(kw):
                kernel_subset()
                torch.cuda.synchronize(device)
                kernel_all()
                torch.cuda.synchronize(device)

            kernel_meas = {"subset": [], "all": []}
            for r in range(kr):
                order = ["subset", "all"] if r % 2 == 0 else ["all", "subset"]
                for method in order:
                    fn = kernel_subset if method == "subset" else kernel_all
                    mm = timed_call(fn, device)
                    kernel_meas[method].append(mm)
                    repeat_rows.append({
                        "target": target_name,
                        "program_id": program_id,
                        "tier": "concordance_kernel",
                        "method": method,
                        "repeat": r + 1,
                        "n_mappings": 1,
                        "seconds": mm["seconds"],
                        "seconds_per_mapping": mm["seconds"],
                        "gpu_peak_incremental_allocated_bytes": mm["gpu_peak_incremental_allocated_bytes"],
                        "gpu_peak_incremental_reserved_bytes": mm["gpu_peak_incremental_reserved_bytes"],
                        "cpu_rss_peak_incremental_bytes": mm["cpu_rss_peak_incremental_bytes"],
                    })

            ss = ssum([x["seconds"] for x in score_meas["subset"]])
            sa = ssum([x["seconds"] for x in score_meas["all"]])
            ks = ssum([x["seconds"] for x in kernel_meas["subset"]])
            ka = ssum([x["seconds"] for x in kernel_meas["all"]])

            shared_corr_median = next(
                row["target_corr_seconds_median"]
                for row in target_rows
                if row["target"] == target_name
            )
            derived_subset = shared_corr_median + mapping_seconds + ss["median"]
            derived_all = shared_corr_median + mapping_seconds + sa["median"]

            pair_rows.append({
                "target": target_name,
                "program_id": program_id,
                "module_genes": module_genes,
                "subset_edges": subset_edges,
                "all_edges": all_edges,
                "subset_edge_fraction": subset_edges / all_edges,
                "edge_count_reduction_fraction": 1.0 - (subset_edges / all_edges),
                "subset_raw_float64_edge_vector_bytes": subset_edges * 8,
                "all_raw_float64_edge_vector_bytes": all_edges * 8,
                "common_mapping_reconstruction_seconds": mapping_seconds,
                "scoring_mappings_per_repeat": nmap,
                "subset_scoring_seconds_median": ss["median"],
                "subset_scoring_seconds_q25": ss["q25"],
                "subset_scoring_seconds_q75": ss["q75"],
                "all_scoring_seconds_median": sa["median"],
                "all_scoring_seconds_q25": sa["q25"],
                "all_scoring_seconds_q75": sa["q75"],
                "scoring_speed_ratio_all_over_subset": sa["median"] / ss["median"],
                "subset_seconds_per_mapping_median": ss["median"] / nmap,
                "all_seconds_per_mapping_median": sa["median"] / nmap,
                "subset_kernel_seconds_median": ks["median"],
                "all_kernel_seconds_median": ka["median"],
                "kernel_speed_ratio_all_over_subset": ka["median"] / ks["median"],
                "subset_scoring_gpu_peak_incremental_allocated_bytes_max": max(
                    x["gpu_peak_incremental_allocated_bytes"] for x in score_meas["subset"]
                ),
                "all_scoring_gpu_peak_incremental_allocated_bytes_max": max(
                    x["gpu_peak_incremental_allocated_bytes"] for x in score_meas["all"]
                ),
                "subset_kernel_gpu_peak_incremental_allocated_bytes_max": max(
                    x["gpu_peak_incremental_allocated_bytes"] for x in kernel_meas["subset"]
                ),
                "all_kernel_gpu_peak_incremental_allocated_bytes_max": max(
                    x["gpu_peak_incremental_allocated_bytes"] for x in kernel_meas["all"]
                ),
                "shared_target_corr_seconds_median": shared_corr_median,
                "derived_eval_workflow_subset_seconds": derived_subset,
                "derived_eval_workflow_all_seconds": derived_all,
                "derived_eval_workflow_speed_ratio_all_over_subset": derived_all / derived_subset,
                "derived_eval_workflow_is_independently_timed": False,
                "upstream_subset_selection_cost_included": False,
            })

            print(
                f"      genes={module_genes:,}; subset/all edges={subset_edges:,}/{all_edges:,}; "
                f"score ratio all/subset={sa['median']/ss['median']:.2f}x; "
                f"kernel ratio={ka['median']/ks['median']:.2f}x"
            )

            del (
                source_corr, tri_i_t, tri_j_t, source_all_edges,
                source_all_rank_centered, source_all_rank_norm,
                ii_t, jj_t, rep_sub, mapped_list,
            )
            torch.cuda.empty_cache()
            gc.collect()

        del target_corr, target, features, target_x
        torch.cuda.empty_cache()
        gc.collect()

    pair = pd.DataFrame(pair_rows)
    repeats = pd.DataFrame(repeat_rows)
    targets = pd.DataFrame(target_rows)

    require(len(pair) == 22, f"Expected 22 pair benchmark rows, found {len(pair)}.")

    summary_rows = []
    summary_cols = [
        "subset_edge_fraction",
        "edge_count_reduction_fraction",
        "scoring_speed_ratio_all_over_subset",
        "kernel_speed_ratio_all_over_subset",
        "derived_eval_workflow_speed_ratio_all_over_subset",
        "subset_seconds_per_mapping_median",
        "all_seconds_per_mapping_median",
        "subset_scoring_gpu_peak_incremental_allocated_bytes_max",
        "all_scoring_gpu_peak_incremental_allocated_bytes_max",
    ]
    for col in summary_cols:
        s = ssum(pair[col].to_numpy(float))
        summary_rows.append({"quantity": col, "n_pairs": len(pair), **s})
    summary = pd.DataFrame(summary_rows)

    # Target-level derived batch totals: shared correlation once + common mapping
    # reconstruction + sum of pair-specific method scoring medians.
    target_batch_rows = []
    for target_name, gg in pair.groupby("target", sort=True):
        corr = float(gg["shared_target_corr_seconds_median"].iloc[0])
        common_map = float(gg["common_mapping_reconstruction_seconds"].sum())
        subset_score = float(gg["subset_scoring_seconds_median"].sum())
        all_score = float(gg["all_scoring_seconds_median"].sum())
        subset_total = corr + common_map + subset_score
        all_total = corr + common_map + all_score
        target_batch_rows.append({
            "target": target_name,
            "shared_target_corr_seconds_median": corr,
            "common_mapping_reconstruction_seconds_sum": common_map,
            "subset_scoring_seconds_sum_of_pair_medians": subset_score,
            "all_scoring_seconds_sum_of_pair_medians": all_score,
            "derived_subset_batch_seconds": subset_total,
            "derived_all_batch_seconds": all_total,
            "derived_batch_speed_ratio_all_over_subset": all_total / subset_total,
            "derived_not_independently_timed": True,
            "upstream_subset_selection_excluded": True,
        })
    target_batch = pd.DataFrame(target_batch_rows)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    for p in [PAIR_OUT, REPEATS_OUT, TARGET_OUT, SUMMARY_OUT, MASTER_OUT]:
        require(not p.exists(), f"Refusing to overwrite benchmark output: {p}")

    pair.to_csv(PAIR_OUT, sep="\t", index=False)
    repeats.to_csv(REPEATS_OUT, sep="\t", index=False)

    # Append target-batch rows below shared preprocessing rows in JSON-friendly way:
    # save a single target TSV with row_type discriminator.
    t1 = targets.copy()
    t1.insert(0, "row_type", "shared_target_preprocessing")
    t2 = target_batch.copy()
    t2.insert(0, "row_type", "derived_target_batch")
    target_combined = pd.concat([t1, t2], ignore_index=True, sort=False)
    target_combined.to_csv(TARGET_OUT, sep="\t", index=False)

    summary.to_csv(SUMMARY_OUT, sep="\t", index=False)

    master = {
        "script_version": SCRIPT_VERSION,
        "status": "PERFORMANCE_BENCHMARK_COMPLETE",
        "scientific_results_changed": False,
        "new_inferential_test": False,
        "benchmark_scope": (
            "controlled-corruption repeated scoring after subset/reference establishment"
        ),
        "whole_pipeline_speedup_claim_permitted": False,
        "upstream_subset_selection_cost_measured": False,
        "cpu_rss_available": psutil is not None,
        "outputs": {
            "by_pair": str(PAIR_OUT),
            "repeat_timings": str(REPEATS_OUT),
            "target_preprocessing_and_derived_batches": str(TARGET_OUT),
            "summary": str(SUMMARY_OUT),
        },
    }
    MASTER_OUT.write_text(json.dumps(master, indent=2) + "\n", encoding="utf-8")

    print("\n" + SEP)
    print("05j4e PERFORMANCE BENCHMARK: COMPLETE")
    print(SEP)
    print("\nOverall descriptive summary:")
    print(summary.to_string(index=False))
    print("\nDerived target-batch accounting:")
    print(target_batch.to_string(index=False))
    print()
    print("Interpretation guards:")
    print("  Shared target-correlation cost attributed to subset speedup: NO")
    print("  Upstream subset-selection cost measured:                  NO")
    print("  Whole MTA pipeline speedup claim permitted:               NO")
    print("  Repeated scoring-stage efficiency claim permitted:        YES")
    print("  New scientific inference:                                 NO")
    print()
    print(f"By pair:  {PAIR_OUT}")
    print(f"Repeats:  {REPEATS_OUT}")
    print(f"Targets:  {TARGET_OUT}")
    print(f"Summary:  {SUMMARY_OUT}")
    print(f"Master:   {MASTER_OUT}")
    print(SEP)

if __name__ == "__main__":
    main()
