$ErrorActionPreference = "Stop"

$SCRIPT_VERSION = "04c-download-tcga-xena-pam50-v1"

$OutDir = "D:\paper4_tcbb_data\paper4_tcbb_tcga_pam50_xena_v1"
$RawFile = Join-Path $OutDir "BRCA_clinicalMatrix.gz"
$SafeFile = Join-Path $OutDir "tcga_pam50_safe_manifest_v1.tsv"
$AuditFile = Join-Path $OutDir "tcga_pam50_extraction_audit_v1.txt"
$HashFile = Join-Path $OutDir "SHA256SUMS.txt"

$Url = "https://tcga.xenahubs.net/download/TCGA.BRCA.sampleMap/BRCA_clinicalMatrix.gz"

Write-Host "===================================================================================================="
Write-Host "Paper 4 / TCBB - download and safely extract TCGA-BRCA PAM50Call_RNAseq from UCSC Xena"
Write-Host "Script version: $SCRIPT_VERSION"
Write-Host "===================================================================================================="
Write-Host ""
Write-Host "Scientific firewall:"
Write-Host "  Raw Xena phenotype file downloaded:                 YES"
Write-Host "  Downstream safe manifest retains:                   sample_id + PAM50Call_RNAseq ONLY"
Write-Host "  Survival/treatment/outcome fields retained:         NO"
Write-Host "  Preservation statistics calculated:                NO"
Write-Host ""

New-Item -ItemType Directory -Force -Path $OutDir | Out-Null

if (-not (Get-Command curl.exe -ErrorAction SilentlyContinue)) {
    throw "curl.exe not found."
}
if (-not (Get-Command python -ErrorAction SilentlyContinue)) {
    throw "python not found. Activate the project .venv first."
}

Write-Host "[1/4] Downloading UCSC Xena TCGA-BRCA phenotype matrix ..."
Write-Host "URL:  $Url"
Write-Host "FILE: $RawFile"

if (Test-Path $RawFile) {
    $existing = (Get-Item $RawFile).Length
    Write-Host ("Existing file: {0:N2} MB; attempting resume." -f ($existing / 1MB))
}

& curl.exe -L --fail --retry 8 --retry-all-errors --retry-delay 5 `
    --connect-timeout 30 -C - -o $RawFile $Url

if ($LASTEXITCODE -eq 33) {
    Write-Warning "Server refused resume; restarting file."
    Remove-Item $RawFile -Force -ErrorAction SilentlyContinue
    & curl.exe -L --fail --retry 8 --retry-all-errors --retry-delay 5 `
        --connect-timeout 30 -o $RawFile $Url
}

if ($LASTEXITCODE -ne 0) {
    throw "Download failed."
}

Write-Host ("Downloaded: {0:N2} MB" -f ((Get-Item $RawFile).Length / 1MB))

Write-Host ""
Write-Host "[2/4] Extracting PAM50Call_RNAseq ONLY into safe manifest ..."

$py = @'
from __future__ import annotations

import csv
import gzip
import os
import re
import sys
from collections import Counter

raw_path, safe_path, audit_path = sys.argv[1:4]

PAM_FIELD = "PAM50Call_RNAseq"
VALID_CANONICAL = {"LumA", "LumB", "Basal", "Her2", "Normal"}

def clean(x):
    if x is None:
        return ""
    s = str(x).strip().strip('"')
    if s.lower() in {"", "na", "nan", "n/a", "not available", "[not available]", "unknown"}:
        return ""
    return s

def canonicalize(x):
    s = clean(x)
    low = re.sub(r"[\s_\-]+", "", s.lower())
    mapping = {
        "luma": "LumA",
        "luminala": "LumA",
        "lumb": "LumB",
        "luminalb": "LumB",
        "basal": "Basal",
        "basallike": "Basal",
        "her2": "Her2",
        "her2enriched": "Her2",
        "her2e": "Her2",
        "normal": "Normal",
        "normallike": "Normal",
    }
    return mapping.get(low, s)

with gzip.open(raw_path, "rt", encoding="utf-8-sig", errors="replace", newline="") as f:
    rows = list(csv.reader(f, delimiter="\t"))

if not rows:
    raise RuntimeError("Downloaded Xena file is empty.")

manifest = []

# Case A: conventional row-oriented table with PAM50 in the header.
header = rows[0]
if PAM_FIELD in header:
    pam_i = header.index(PAM_FIELD)

    id_candidates = [
        "sample", "sampleID", "sample_id", "SAMPLE_ID",
        "_PATIENT", "patient", "PATIENT_ID"
    ]
    id_i = None
    for candidate in id_candidates:
        if candidate in header:
            id_i = header.index(candidate)
            break
    if id_i is None:
        id_i = 0

    for r in rows[1:]:
        if len(r) <= max(id_i, pam_i):
            continue
        sid = clean(r[id_i])
        pam = canonicalize(r[pam_i])
        if sid:
            manifest.append((sid, pam))

# Case B: Xena phenotype matrix with phenotype rows and sample columns.
else:
    pam_row = None
    for r in rows[1:]:
        if r and clean(r[0]) == PAM_FIELD:
            pam_row = r
            break

    if pam_row is None:
        raise RuntimeError(
            f"Could not find '{PAM_FIELD}' either as a column or phenotype row. "
            f"First header fields: {header[:20]}"
        )

    sample_ids = [clean(x) for x in header[1:]]
    values = [canonicalize(x) for x in pam_row[1:]]

    if len(sample_ids) != len(values):
        raise RuntimeError(
            f"Wide-matrix sample/PAM50 length mismatch: {len(sample_ids)} vs {len(values)}"
        )

    manifest = list(zip(sample_ids, values))

# Drop blank IDs; preserve missing PAM50 as blank.
manifest = [(sid, pam) for sid, pam in manifest if sid]

# Require unique sample IDs in the phenotype file.
ids = [sid for sid, _ in manifest]
if len(ids) != len(set(ids)):
    dup = [x for x, n in Counter(ids).items() if n > 1]
    raise RuntimeError(f"Duplicate sample IDs in extracted Xena PAM50 data: {dup[:20]}")

counts = Counter(pam for _, pam in manifest if pam)
canonical_counts = {k: counts.get(k, 0) for k in ["LumA", "LumB", "Basal", "Her2", "Normal"]}
other_counts = {k: v for k, v in counts.items() if k not in VALID_CANONICAL}

with open(safe_path, "w", encoding="utf-8", newline="") as out:
    w = csv.writer(out, delimiter="\t")
    w.writerow(["xena_sample_id", "pam50_subtype"])
    w.writerows(manifest)

with open(audit_path, "w", encoding="utf-8") as out:
    out.write(f"source_file={raw_path}\n")
    out.write(f"field={PAM_FIELD}\n")
    out.write(f"rows_with_sample_id={len(manifest)}\n")
    out.write(f"nonmissing_pam50={sum(canonical_counts.values()) + sum(other_counts.values())}\n")
    out.write(f"canonical_counts={canonical_counts}\n")
    out.write(f"other_nonempty_labels={other_counts}\n")

print("Extracted rows with sample ID:", len(manifest))
print("Nonmissing PAM50:", sum(canonical_counts.values()) + sum(other_counts.values()))
print("Canonical PAM50 counts:", canonical_counts)
print("Other nonempty labels:", other_counts)
'@

$tmpPy = Join-Path $env:TEMP "extract_tcga_pam50_xena_v1.py"
Set-Content -Path $tmpPy -Value $py -Encoding UTF8

& python $tmpPy $RawFile $SafeFile $AuditFile
$extractExit = $LASTEXITCODE
Remove-Item $tmpPy -Force -ErrorAction SilentlyContinue

if ($extractExit -ne 0) {
    throw "PAM50 extraction failed."
}

Write-Host ""
Write-Host "[3/4] Verifying gzip integrity ..."
$gzipCheck = @'
import gzip, sys
p = sys.argv[1]
n = 0
with gzip.open(p, "rb") as f:
    while True:
        b = f.read(8 * 1024 * 1024)
        if not b:
            break
        n += len(b)
print("GZIP_OK uncompressed_bytes=", n)
'@

$tmpCheck = Join-Path $env:TEMP "check_tcga_pam50_gzip_v1.py"
Set-Content -Path $tmpCheck -Value $gzipCheck -Encoding UTF8
& python $tmpCheck $RawFile
$checkExit = $LASTEXITCODE
Remove-Item $tmpCheck -Force -ErrorAction SilentlyContinue

if ($checkExit -ne 0) {
    throw "Gzip integrity check failed."
}

Write-Host ""
Write-Host "[4/4] Writing SHA-256 manifest ..."
Remove-Item $HashFile -Force -ErrorAction SilentlyContinue

foreach ($f in @($RawFile, $SafeFile, $AuditFile)) {
    $hash = (Get-FileHash -Algorithm SHA256 $f).Hash.ToLowerInvariant()
    "$hash  $([System.IO.Path]::GetFileName($f))" | Add-Content -Encoding UTF8 $HashFile
}

Write-Host ""
Write-Host "===================================================================================================="
Write-Host "04c TCGA PAM50 DOWNLOAD / SAFE EXTRACTION: PASS"
Write-Host "===================================================================================================="
Write-Host "Raw phenotype file: $RawFile"
Write-Host "Safe PAM50 manifest: $SafeFile"
Write-Host "Extraction audit:    $AuditFile"
Write-Host "SHA-256:             $HashFile"
Write-Host ""
Write-Host "IMPORTANT: downstream subtype analysis should read ONLY the safe PAM50 manifest, not the raw Xena clinical matrix."
Write-Host "===================================================================================================="
