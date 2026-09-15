#!/usr/bin/env python
from __future__ import annotations

import hashlib
import json
import platform
from datetime import datetime, timezone
from pathlib import Path

import torch

SCRIPT_VERSION = "05j4e-freeze-performance-benchmark-v1-no-cli"

PROJECT_ROOT = Path(r"C:\Users\olegk\Desktop\molecular-transport-audit")
DATA_ROOT = Path(r"D:\paper4_tcbb_data")

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

OUT_DIR = DATA_ROOT / "paper4_tcbb_performance_benchmark_contract_v1"
OUT_JSON = OUT_DIR / "performance_benchmark_contract_v1.json"

BENCH_MAPPINGS_PER_PAIR = 20
BENCH_VALID_NULL_IDS = list(range(1, BENCH_MAPPINGS_PER_PAIR + 1))
SCORING_WARMUPS = 2
SCORING_REPEATS = 5
KERNEL_WARMUPS = 3
KERNEL_REPEATS = 10
TARGET_CORR_WARMUPS = 1
TARGET_CORR_REPEATS = 3

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

def main() -> None:
    print(SEP)
    print("Paper 4 / TCBB - freeze computational performance benchmark")
    print(SEP)
    print(f"Script version: {SCRIPT_VERSION}")
    print()
    print("Scientific provenance:")
    print("  05j4/05j4c scientific results already observed: YES")
    print("  Performance benchmark results observed:         NO")
    print("  Scientific thresholds/statistics changed here:  NO")
    print()

    required = [SOURCE_05J4B, SOURCE_05J, SOURCE_05D, SOURCE_SPEARMAN, FULL_NULL, CLASSIFICATION]
    for p in required:
        require(p.exists(), f"Required file missing: {p}")

    require(torch.cuda.is_available(), "CUDA device required for the frozen benchmark.")
    dev = torch.device("cuda")
    prop = torch.cuda.get_device_properties(dev)

    contract = {
        "script_version": SCRIPT_VERSION,
        "frozen_at_utc": datetime.now(timezone.utc).isoformat(),
        "project": "Paper 4 / TCBB / molecular-transport-audit",
        "status": "FROZEN_BEFORE_PERFORMANCE_BENCHMARK_RESULTS",
        "purpose": (
            "Quantify the computational effect of using the frozen matched edge subset "
            "versus the complete within-module edge universe without making a new "
            "scientific superiority claim."
        ),
        "scope_guard": {
            "benchmark_claim_permitted": (
                "computational efficiency of the repeated controlled-corruption scoring "
                "stage after source references/subset have been established"
            ),
            "total_MTA_design_pipeline_speedup_claim_permitted": False,
            "upstream_subset_selection_cost_measured": False,
            "reason": (
                "the frozen subset was selected upstream; this benchmark must not treat "
                "that selection/design cost as zero"
            ),
        },
        "methods": {
            "subset_method": "A: MTA Spearman on frozen matched edge subset",
            "all_edge_method": "C: Spearman on all strict upper-triangle module edges",
            "scientific_difference": "edge universe only",
            "same_mappings": True,
            "same_target_correlation_matrix": True,
        },
        "fixed_benchmark_workload": {
            "target_program_pairs": 22,
            "valid_null_ids_per_pair": BENCH_VALID_NULL_IDS,
            "mappings_per_pair": BENCH_MAPPINGS_PER_PAIR,
            "mapping_source": "existing frozen full-null attempt_id; no new mapping identity",
            "corruption_fraction": 1.0,
            "corruption_fraction_rationale": (
                "scoring computational complexity is independent of corruption dose; "
                "fraction=1.0 provides an already-frozen deterministic mapping set"
            ),
        },
        "measurement_tiers": {
            "tier_0_edge_census": {
                "metrics": [
                    "subset edge count",
                    "all-edge count",
                    "subset/all edge fraction",
                    "edge-count reduction",
                    "float64 raw vector bytes",
                ],
                "timing": False,
            },
            "tier_1_shared_target_preprocessing": {
                "operation": "full target Pearson correlation matrix on GPU",
                "shared_by_methods": True,
                "warmups": TARGET_CORR_WARMUPS,
                "timed_repeats": TARGET_CORR_REPEATS,
                "claim_guard": "reported separately; not attributed as a subset-method speedup",
            },
            "tier_2_concordance_kernel": {
                "operation": (
                    "with mapped target sub-correlation matrix already constructed, "
                    "extract required edges and compute Spearman concordance"
                ),
                "warmups": KERNEL_WARMUPS,
                "timed_repeats": KERNEL_REPEATS,
                "representative_mapping": "lowest frozen valid_null_id in the fixed workload",
            },
            "tier_3_per_mapping_scoring": {
                "operation": (
                    "for the fixed batch of mapped indices: construct mapped target "
                    "sub-correlation matrix + extract edges + compute Spearman"
                ),
                "warmups": SCORING_WARMUPS,
                "timed_repeats": SCORING_REPEATS,
                "method_order": (
                    "balanced by repeat: subset->all on even repeats; "
                    "all->subset on odd repeats"
                ),
            },
            "tier_4_derived_evaluation_workflow": {
                "definition": (
                    "shared measured target-correlation time + common measured mapping "
                    "reconstruction time + method-specific median repeated-scoring time"
                ),
                "label": "derived, not independently timed end-to-end",
                "upstream_subset_selection_excluded": True,
            },
        },
        "timing_rules": {
            "clock": "time.perf_counter",
            "cuda_synchronize_before_and_after_timed_region": True,
            "cuda_warmup_required": True,
            "same_process_same_device": True,
            "no_concurrent_benchmark_jobs_expected": True,
            "summary": "median, IQR, min, max of frozen repeats; no p-values",
        },
        "memory_rules": {
            "gpu": (
                "report baseline allocated/reserved and peak incremental allocated/reserved "
                "for each timed method block"
            ),
            "cpu": (
                "report sampled process RSS peak increment if psutil is available; "
                "otherwise NA without failing the benchmark"
            ),
            "analytical_edge_vector_memory": True,
        },
        "interpretation_rules": {
            "scoring_speedup_if_large_but_derived_workflow_small": (
                "claim efficiency specifically for repeated concordance scoring, not the "
                "whole pipeline"
            ),
            "derived_workflow_speedup_if_material": (
                "may report practical savings for the controlled-corruption evaluation stage "
                "after subset/reference construction"
            ),
            "no_scientific_sensitivity_reanalysis": True,
            "no_new_thresholds_or_corruption_grids": True,
        },
        "hardware_lock": {
            "platform": platform.platform(),
            "python": platform.python_version(),
            "torch": torch.__version__,
            "cuda_runtime": torch.version.cuda,
            "gpu_name": prop.name,
            "gpu_total_memory_bytes": int(prop.total_memory),
            "gpu_capability": [int(prop.major), int(prop.minor)],
        },
        "input_sha256": {str(p): sha256_file(p) for p in required},
    }

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    require(not OUT_JSON.exists(), f"Refusing to overwrite frozen benchmark contract: {OUT_JSON}")
    OUT_JSON.write_text(json.dumps(contract, indent=2) + "\n", encoding="utf-8")

    print(SEP)
    print("05j4e PERFORMANCE BENCHMARK CONTRACT: PASS")
    print(SEP)
    print(f"GPU: {prop.name}; VRAM={prop.total_memory/1024**3:.2f} GB")
    print(f"Fixed mappings/pair:              {BENCH_MAPPINGS_PER_PAIR}")
    print("Shared target correlation separated: YES")
    print("Scoring-only benchmark:              YES")
    print("Concordance-kernel benchmark:        YES")
    print("Derived workflow benchmark:          YES")
    print("Upstream subset-selection cost measured: NO")
    print("Whole MTA pipeline speedup claim allowed: NO")
    print("New scientific inference:            NO")
    print()
    print(f"Output: {OUT_JSON}")
    print(SEP)

if __name__ == "__main__":
    main()
