$ErrorActionPreference = "Stop"

$SCRIPT_VERSION = "00c-install-pytorch-cuda-v1"

Write-Host "===================================================================================================="
Write-Host "Paper 4 / TCBB - install CUDA-enabled PyTorch into the ACTIVE Python virtual environment"
Write-Host "Script version: $SCRIPT_VERSION"
Write-Host "===================================================================================================="

if (-not (Get-Command python -ErrorAction SilentlyContinue)) {
    throw "python was not found on PATH. Activate the project .venv first."
}

$pythonExe = (Get-Command python).Source
Write-Host "Python executable: $pythonExe"

Write-Host ""
Write-Host "[1/4] Checking NVIDIA driver / CUDA compatibility ..."
if (-not (Get-Command nvidia-smi -ErrorAction SilentlyContinue)) {
    throw "nvidia-smi was not found. Install/update the NVIDIA driver first."
}

$smi = (& nvidia-smi | Out-String)
Write-Host $smi

$match = [regex]::Match($smi, 'CUDA Version:\s*([0-9]+)\.([0-9]+)')
if (-not $match.Success) {
    throw "Could not parse the maximum CUDA version from nvidia-smi."
}

$cudaMajor = [int]$match.Groups[1].Value
$cudaMinor = [int]$match.Groups[2].Value
$cudaAsNumber = [double]("$cudaMajor.$cudaMinor")

Write-Host "Maximum CUDA version reported by driver: $cudaAsNumber"

# Prefer the current stable PyTorch 2.14 / CUDA 13.0 wheel when the driver supports it.
# Otherwise use the well-supported CUDA 12.8 wheel, which is appropriate for recent NVIDIA GPUs.
if ($cudaAsNumber -ge 13.0) {
    $torchVersion = "2.14.0"
    $indexUrl = "https://download.pytorch.org/whl/cu130"
    $buildLabel = "CUDA 13.0"
}
elseif ($cudaAsNumber -ge 12.8) {
    $torchVersion = "2.11.0"
    $indexUrl = "https://download.pytorch.org/whl/cu128"
    $buildLabel = "CUDA 12.8"
}
else {
    throw "The NVIDIA driver reports CUDA $cudaAsNumber. Update the NVIDIA driver so it supports at least CUDA 12.8 before installing the GPU build."
}

Write-Host ""
Write-Host "[2/4] Selected PyTorch build:"
Write-Host "  torch version: $torchVersion"
Write-Host "  runtime:       $buildLabel"
Write-Host "  index:         $indexUrl"

Write-Host ""
Write-Host "[3/4] Installing torch into the active .venv ..."
& python -m pip install --upgrade "torch==$torchVersion" --index-url $indexUrl
if ($LASTEXITCODE -ne 0) {
    throw "PyTorch installation failed."
}

Write-Host ""
Write-Host "[4/4] Verifying CUDA from this exact Python environment ..."
$verify = @'
import json
import sys
import torch

print("python:", sys.executable)
print("torch:", torch.__version__)
print("torch CUDA build:", torch.version.cuda)
print("cuda available:", torch.cuda.is_available())

if not torch.cuda.is_available():
    raise SystemExit("CUDA is not available to PyTorch after installation.")

print("device count:", torch.cuda.device_count())
for i in range(torch.cuda.device_count()):
    p = torch.cuda.get_device_properties(i)
    print(f"device {i}: {p.name}, VRAM={p.total_memory / 1024**3:.2f} GB, capability={p.major}.{p.minor}")
'@

$tmp = Join-Path $env:TEMP "paper4_verify_torch_cuda_v1.py"
Set-Content -Path $tmp -Value $verify -Encoding UTF8
& python $tmp
$verifyExit = $LASTEXITCODE
Remove-Item $tmp -Force -ErrorAction SilentlyContinue

if ($verifyExit -ne 0) {
    throw "PyTorch installed, but CUDA verification failed."
}

Write-Host ""
Write-Host "===================================================================================================="
Write-Host "PYTORCH CUDA INSTALL: PASS"
Write-Host "===================================================================================================="
Write-Host "Next run:"
Write-Host "  python scripts\00b_audit_gpu_compute_stack_v1.py"
Write-Host "===================================================================================================="
