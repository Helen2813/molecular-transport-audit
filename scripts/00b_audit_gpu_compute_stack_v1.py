from __future__ import annotations

import json
import platform
import sys
import time
from pathlib import Path

SCRIPT_VERSION = "00b-audit-gpu-compute-stack-v1-no-cli"
OUT_DIR = Path(r"D:\paper4_tcbb_data\paper4_tcbb_gpu_audit_v1")


def main() -> None:
    print("=" * 104)
    print("Paper 4 / TCBB - GPU compute-stack audit")
    print("=" * 104)
    print(f"Script version: {SCRIPT_VERSION}")
    print(f"Python: {sys.version.split()[0]}")
    print(f"Platform: {platform.platform()}")

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    result = {
        "script_version": SCRIPT_VERSION,
        "python": sys.version,
        "platform": platform.platform(),
        "torch_importable": False,
        "cuda_available": False,
    }

    try:
        import torch
    except Exception as exc:
        print()
        print("PyTorch import: FAIL")
        print(f"Reason: {exc}")
        result["torch_import_error"] = repr(exc)
    else:
        result["torch_importable"] = True
        result["torch_version"] = torch.__version__
        result["torch_cuda_build"] = torch.version.cuda
        result["cuda_available"] = bool(torch.cuda.is_available())

        print()
        print(f"PyTorch: {torch.__version__}")
        print(f"PyTorch CUDA build: {torch.version.cuda}")
        print(f"CUDA available: {torch.cuda.is_available()}")

        if torch.cuda.is_available():
            n = torch.cuda.device_count()
            result["cuda_device_count"] = n
            devices = []

            for i in range(n):
                prop = torch.cuda.get_device_properties(i)
                total_gb = prop.total_memory / (1024 ** 3)
                devices.append(
                    {
                        "index": i,
                        "name": prop.name,
                        "total_memory_gb": total_gb,
                        "major": prop.major,
                        "minor": prop.minor,
                    }
                )
                print(
                    f"Device {i}: {prop.name}, VRAM={total_gb:.2f} GB, "
                    f"capability={prop.major}.{prop.minor}"
                )

            result["devices"] = devices

            # Small deterministic FP64 smoke test only.
            torch.manual_seed(20260916)
            device = torch.device("cuda:0")
            a = torch.randn((2048, 1024), dtype=torch.float64, device=device)
            b = torch.randn((1024, 2048), dtype=torch.float64, device=device)
            torch.cuda.synchronize()
            t0 = time.perf_counter()
            c = a @ b
            torch.cuda.synchronize()
            elapsed = time.perf_counter() - t0
            checksum = float(c[:16, :16].sum().cpu())

            print(f"FP64 smoke matmul: PASS ({elapsed:.3f} s)")
            print(f"Deterministic checksum: {checksum:.12g}")

            result["fp64_smoke_test"] = {
                "status": "PASS",
                "seconds": elapsed,
                "checksum": checksum,
            }

            del a, b, c
            torch.cuda.empty_cache()

    out = OUT_DIR / "gpu_compute_stack_audit_v1.json"
    out.write_text(json.dumps(result, indent=2), encoding="utf-8")

    print()
    print("=" * 104)
    if result["cuda_available"]:
        print("GPU AUDIT: CUDA READY")
    elif result["torch_importable"]:
        print("GPU AUDIT: PYTORCH PRESENT, CUDA NOT AVAILABLE")
    else:
        print("GPU AUDIT: PYTORCH NOT INSTALLED IN THIS PYTHON ENVIRONMENT")
    print(f"Output: {out}")
    print("=" * 104)


if __name__ == "__main__":
    main()
