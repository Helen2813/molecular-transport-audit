#!/usr/bin/env python
from __future__ import annotations

import hashlib
import json
import platform
from datetime import datetime, timezone
from pathlib import Path

import torch

SCRIPT_VERSION = "05j4f-freeze-backend-matched-kernel-benchmark-v1-no-cli"

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
PERF_MASTER = (
    DATA_ROOT
    / "paper4_tcbb_performance_benchmark_v1"
    / "performance_benchmark_v1.json"
)

OUT_DIR = DATA_ROOT / "paper4_tcbb_backend_matched_kernel_benchmark_contract_v1"
OUT_JSON = OUT_DIR / "backend_matched_kernel_benchmark_contract_v1.json"

# Fixed before any 05j4f timing result.
VALID_NULL_IDS = [1, 5, 10, 15, 20]
WARMUPS_PER_MAPPING = 3
TIMED_REPEATS_PER_MAPPING = 10

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
    print("Paper 4 / TCBB - freeze backend-matched Spearman kernel benchmark")
    print(SEP)
    print(f"Script version: {SCRIPT_VERSION}")
    print()
    print("Provenance:")
    print("  05j4e implementation-level performance results observed: YES")
    print("  05j4f backend-matched timing results observed:           NO")
    print("  Scientific operating results changed:                    NO")
    print()

    required = [
        SOURCE_05J4B,
        SOURCE_05J,
        SOURCE_05D,
        SOURCE_SPEARMAN,
        FULL_NULL,
        ALL_EDGE_FULL_NULL,
        CLASSIFICATION,
        PERF_MASTER,
    ]
    for p in required:
        require(p.exists(), f"Required file missing: {p}")

    pm = json.loads(PERF_MASTER.read_text(encoding="utf-8"))
    require(
        pm.get("status") == "PERFORMANCE_BENCHMARK_COMPLETE",
        "05j4e performance master is not complete.",
    )

    require(torch.cuda.is_available(), "CUDA GPU required for frozen benchmark.")
    device = torch.device("cuda")
    prop = torch.cuda.get_device_properties(device)

    contract = {
        "script_version": SCRIPT_VERSION,
        "frozen_at_utc": datetime.now(timezone.utc).isoformat(),
        "project": "Paper 4 / TCBB / molecular-transport-audit",
        "status": "FROZEN_BEFORE_BACKEND_MATCHED_KERNEL_RESULTS",
        "scientific_role": (
            "computational microbenchmark motivated by 05j4e implementation asymmetry; "
            "not a new preservation experiment"
        ),
        "motivation_locked_before_result": (
            "05j4e compared current implementations in which subset Spearman crossed "
            "GPU->CPU and all-edge Spearman remained on GPU. 05j4f isolates intrinsic "
            "edge-universe kernel cost by using the identical GPU Spearman backend for both."
        ),
        "methods": {
            "subset": "Spearman on frozen matched subset",
            "all_edges": "Spearman on all strict upper-triangle module edges",
            "identical_backend": (
                "05j4c gpu_average_ranks + fixed_rank_vector + "
                "spearman_against_fixed_gpu for both methods"
            ),
            "mapped_target_submatrix_built_before_timing": True,
            "mapping_generation_inside_timing": False,
            "target_correlation_inside_timing": False,
            "source_reference_construction_inside_timing": False,
            "edge_extraction_inside_timing": True,
            "rank_and_spearman_inside_timing": True,
        },
        "fixed_workload": {
            "target_program_pairs": 22,
            "corruption_fraction": 1.0,
            "valid_null_ids": VALID_NULL_IDS,
            "mappings_per_pair": len(VALID_NULL_IDS),
            "mapping_identity": "reconstruct exact saved attempt_id; no new mapping identity",
            "warmups_per_mapping": WARMUPS_PER_MAPPING,
            "timed_repeats_per_mapping_per_method": TIMED_REPEATS_PER_MAPPING,
            "method_order": (
                "balanced by timed repeat; subset->all on even-index repeats, "
                "all->subset on odd-index repeats"
            ),
        },
        "score_replay_guards": {
            "subset_score": (
                "backend-matched GPU subset score must reproduce saved frozen "
                "mta_rho_edge for each fixed mapping"
            ),
            "all_edge_score": (
                "backend-matched GPU all-edge score must reproduce saved "
                "spearman_all_edges for each fixed mapping"
            ),
            "tolerance": 1e-10,
        },
        "timing": {
            "clock": "time.perf_counter",
            "cuda_synchronize_before_after": True,
            "empty_cache_inside_timed_region": False,
            "submatrix_allocation_excluded": True,
            "summary": (
                "per mapping and per pair: median/IQR/min/max timing; "
                "all/subset time ratio; no p-values"
            ),
        },
        "memory": {
            "submatrix_memory_excluded_from_incremental_kernel_memory": True,
            "gpu_peak_incremental_allocated_bytes": True,
            "analytical_raw_edge_vector_bytes": True,
        },
        "pre_result_interpretation": {
            "if_all_edge_slower": (
                "the restricted edge universe reduces intrinsic concordance-kernel cost; "
                "this does not override the 05j4e finding that the current full evaluation "
                "workflow showed negligible end-to-end speedup"
            ),
            "if_ratio_near_one": (
                "edge-count reduction does not materially change intrinsic kernel runtime "
                "on this GPU/backend over the tested module sizes"
            ),
            "if_all_edge_faster": (
                "GPU execution efficiency dominates edge-count savings in this implementation; "
                "make no subset-efficiency claim"
            ),
            "whole_pipeline_speedup_claim_permitted": False,
            "current_implementation_speedup_claim_replaced": False,
            "no_post_hoc_workload_or_backend_change": True,
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
    require(not OUT_JSON.exists(), f"Refusing to overwrite frozen contract: {OUT_JSON}")
    OUT_JSON.write_text(json.dumps(contract, indent=2) + "\n", encoding="utf-8")

    print(SEP)
    print("05j4f BACKEND-MATCHED KERNEL BENCHMARK CONTRACT: PASS")
    print(SEP)
    print(f"GPU: {prop.name}; VRAM={prop.total_memory/1024**3:.2f} GB")
    print(f"Fixed mappings/pair:            {len(VALID_NULL_IDS)}")
    print(f"Timed repeats/mapping/method:   {TIMED_REPEATS_PER_MAPPING}")
    print("Mapped submatrix excluded:       YES")
    print("Identical GPU Spearman backend:  YES")
    print("Score replay required:           YES")
    print("Whole-pipeline speedup claim:    NO")
    print("New scientific inference:        NO")
    print()
    print(f"Output: {OUT_JSON}")
    print(SEP)

if __name__ == "__main__":
    main()
