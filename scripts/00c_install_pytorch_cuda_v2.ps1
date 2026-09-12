$ErrorActionPreference = "Stop"

$SCRIPT_VERSION = "00c-install-pytorch-cuda-v2"

Write-Host "===================================================================================================="
Write-Host "Paper 4 / TCBB - install CUDA-enabled PyTorch into the ACTIVE Python virtual environment"
Write-Host "Script version: $SCRIPT_VERSION"
Write-Host "===================================================================================================="

if (-not (Get-Command python -ErrorAction SilentlyContinue)) {
    throw "python was not found on PATH. Activate the project .venv first."
}

$pythonExe = (Get-Command python).Source
Write-Host "Python executable: $pythonExe"

if ($pythonExe -notmatch '\\.venv\\Scripts\\python\.exe$') {
    Write-Warning "The active Python does not look like the project .venv. Continuing, but verify this is intentional."
}

Write-Host ""
Write-Host "[1/5] Checking NVIDIA driver / CUDA compatibility ..."

if (-not (Get-Command nvidia-smi -ErrorAction SilentlyContinue)) {
    throw "nvidia-smi was not found. Install/update the NVIDIA driver first."
}

$smi = (& nvidia-smi | Out-String)
Write-Host $smi

# Newer NVIDIA Windows drivers print "CUDA UMD Version", while older releases
# commonly print "CUDA Version". Accept both formats.
$patterns = @(
    'CUDA UMD Version:\s*([0-9]+(?:\.[0-9]+)?)',
    'CUDA Version:\s*([0-9]+(?:\.[0-9]+)?)'
)

$cudaVersionString = $null
foreach ($pattern in $patterns) {
    $match = [regex]::Match($smi, $pattern)
    if ($match.Success) {
        $cudaVersionString = $match.Groups[1].Value
        break
    }
}

if (-not $cudaVersionString) {
    throw "Could not parse CUDA compatibility from nvidia-smi. Expected 'CUDA UMD Version:' or 'CUDA Version:'."
}

$cudaParts = $cudaVersionString.Split('.')
$cudaMajor = [int]$cudaParts[0]
$cudaMinor = if ($cudaParts.Count -gt 1) { [int]$cudaParts[1] } else { 0 }
$cudaNumeric = $cudaMajor + ($cudaMinor / 10.0)

Write-Host "Maximum CUDA compatibility reported by driver: $cudaVersionString"

# Prefer CUDA 13.0 wheels on a CUDA 13.x driver.
# cu128 is retained as a safe fallback because newer NVIDIA drivers are backward compatible.
$candidateIndexes = @()

if ($cudaMajor -ge 13) {
    $candidateIndexes += "https://download.pytorch.org/whl/cu130"
    $candidateIndexes += "https://download.pytorch.org/whl/cu128"
}
elseif (($cudaMajor -eq 12 -and $cudaMinor -ge 8) -or $cudaMajor -gt 12) {
    $candidateIndexes += "https://download.pytorch.org/whl/cu128"
}
else {
    throw "The NVIDIA driver reports CUDA compatibility $cudaVersionString. Update the driver to support at least CUDA 12.8."
}

Write-Host ""
Write-Host "[2/5] Upgrading pip in the active .venv ..."
& python -m pip install --upgrade pip
if ($LASTEXITCODE -ne 0) {
    throw "pip upgrade failed."
}

Write-Host ""
Write-Host "[3/5] Installing CUDA-enabled torch ..."
$installed = $false
$usedIndex = $null

foreach ($indexUrl in $candidateIndexes) {
    Write-Host ""
    Write-Host "Trying PyTorch wheel index: $indexUrl"

    & python -m pip install --upgrade torch --index-url $indexUrl

    if ($LASTEXITCODE -eq 0) {
        $installed = $true
        $usedIndex = $indexUrl
        break
    }

    Write-Warning "Installation from $indexUrl failed. Trying the next compatible CUDA wheel index."
}

if (-not $installed) {
    throw "CUDA-enabled PyTorch installation failed from all compatible wheel indexes."
}

Write-Host ""
Write-Host "[4/5] Verifying torch/CUDA from this exact Python environment ..."

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
    print(
        f"device {i}: {p.name}, "
        f"VRAM={p.total_memory / 1024**3:.2f} GB, "
        f"capability={p.major}.{p.minor}"
    )

# Tiny FP64 correctness smoke test.
torch.manual_seed(20260916)
device = torch.device("cuda:0")
a = torch.randn((1024, 512), dtype=torch.float64, device=device)
b = torch.randn((512, 1024), dtype=torch.float64, device=device)
torch.cuda.synchronize()
c = a @ b
torch.cuda.synchronize()

checksum = float(c[:16, :16].sum().cpu())
print("FP64 smoke test: PASS")
print("FP64 checksum:", f"{checksum:.12g}")
'@

$tmp = Join-Path $env:TEMP "paper4_verify_torch_cuda_v2.py"
Set-Content -Path $tmp -Value $verify -Encoding UTF8

& python $tmp
$verifyExit = $LASTEXITCODE

Remove-Item $tmp -Force -ErrorAction SilentlyContinue

if ($verifyExit -ne 0) {
    throw "PyTorch installed, but CUDA verification failed."
}

Write-Host ""
Write-Host "[5/5] Recording installed package state ..."

$envOut = "D:\paper4_tcbb_data\paper4_tcbb_gpu_environment_v2"
New-Item -ItemType Directory -Force -Path $envOut | Out-Null

& python -m pip freeze | Set-Content -Encoding UTF8 (Join-Path $envOut "pip_freeze_after_torch_v2.txt")

$summary = @"
script_version=$SCRIPT_VERSION
python=$pythonExe
driver_cuda_compatibility=$cudaVersionString
pytorch_index=$usedIndex
"@

Set-Content -Path (Join-Path $envOut "gpu_install_summary_v2.txt") -Value $summary -Encoding UTF8

Write-Host ""
Write-Host "===================================================================================================="
Write-Host "PYTORCH CUDA INSTALL: PASS"
Write-Host "===================================================================================================="
Write-Host "Driver CUDA compatibility: $cudaVersionString"
Write-Host "PyTorch wheel index used:   $usedIndex"
Write-Host "Environment snapshot:       $envOut"
Write-Host ""
Write-Host "Next run:"
Write-Host "  python scripts\00b_audit_gpu_compute_stack_v1.py"
Write-Host "===================================================================================================="
