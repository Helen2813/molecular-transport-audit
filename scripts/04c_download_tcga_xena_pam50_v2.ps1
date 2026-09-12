$ErrorActionPreference = "Stop"

$SCRIPT_VERSION = "04c-download-tcga-xena-pam50-v2"

$OutDir = "D:\paper4_tcbb_data\paper4_tcbb_tcga_pam50_xena_v2"
$RawFile = Join-Path $OutDir "BRCA_clinicalMatrix.txt"
$SafeFile = Join-Path $OutDir "tcga_pam50_safe_manifest_v2.tsv"
$AuditFile = Join-Path $OutDir "tcga_pam50_extraction_audit_v2.txt"
$HashFile = Join-Path $OutDir "SHA256SUMS.txt"

# Current direct S3 object first; legacy Xena endpoint second as fallback.
$Urls = @(
    "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.BRCA.sampleMap%2FBRCA_clinicalMatrix",
    "https://tcga.xenahubs.net/download/TCGA.BRCA.sampleMap/BRCA_clinicalMatrix"
)

Write-Host "===================================================================================================="
Write-Host "Paper 4 / TCBB - download and safely extract TCGA-BRCA PAM50Call_RNAseq from UCSC Xena"
Write-Host "Script version: $SCRIPT_VERSION"
Write-Host "===================================================================================================="
Write-Host ""
Write-Host "v2 correction:"
Write-Host "  - uses the direct current Xena S3 object first"
Write-Host "  - treats BRCA_clinicalMatrix as PLAIN TEXT, not gzip"
Write-Host "  - writes to a new v2 directory; the failed v1 download is not reused"
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

$downloaded = $false
$usedUrl = $null

foreach ($Url in $Urls) {
    Write-Host ""
    Write-Host "Trying:"
    Write-Host "  $Url"

    Remove-Item $RawFile -Force -ErrorAction SilentlyContinue

    & curl.exe -L --fail --retry 3 --retry-all-errors --retry-delay 3 `
        --connect-timeout 30 `
        -A "Mozilla/5.0" `
        -o $RawFile `
        $Url

    if ($LASTEXITCODE -eq 0 -and (Test-Path $RawFile) -and ((Get-Item $RawFile).Length -gt 1000) {
        $downloaded = $true
        $usedUrl = $Url
        break
    }

    Write-Warning "This endpoint failed or returned an implausibly small file."
    Remove-Item $RawFile -Force -ErrorAction SilentlyContinue
}

if (-not $downloaded) {
    throw "All Xena download endpoints failed."
}

Write-Host ("Downloaded: {0:N2} MB" -f ((Get-Item $RawFile).Length / 1MB))
Write-Host "Endpoint used: $usedUrl"

Write-Host ""
Write-Host "[2/4] Extracting PAM50Call_RNAseq ONLY ..."

$py = @'
from __future__ import annotations

import csv
import os
import re
import sys
from collections import Counter

raw_path, safe_path, audit_path = sys.argv[1:4]
PAM_FIELD = "PAM50Call_RNAseq"
VALID = {"LumA", "LumB", "Basal", "Her2", "Normal"}

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

with open(raw_path, "rt", encoding="utf-8-sig", errors="replace", newline="") as f:
    reader = csv.reader(f, delimiter="\t")
    rows = list(reader)

if not rows:
    raise RuntimeError("Downloaded Xena file is empty.")

header = rows[0]
manifest = []
layout = None

# Xena BRCA clinical matrix is normally row-oriented:
# sampleID <tab> phenotype1 <tab> ... <tab> PAM50Call_RNAseq ...
if PAM_FIELD in header:
    layout = "row_oriented"
    pam_i = header.index(PAM_FIELD)

    id_candidates = [
        "sampleID", "sample", "sample_id", "SAMPLE_ID",
        "_PATIENT", "patient", "PATIENT_ID",
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

# Defensive support for a transposed Xena matrix if format ever changes.
else:
    pam_row = None
    for r in rows[1:]:
        if r and clean(r[0]) == PAM_FIELD:
            pam_row = r
            break

    if pam_row is None:
        raise RuntimeError(
            f"Could not find '{PAM_FIELD}'. "
            f"Header has {len(header)} fields; first fields: {header[:30]}"
        )

    layout = "phenotype_rows"
    sample_ids = [clean(x) for x in header[1:]]
    values = [canonicalize(x) for x in pam_row[1:]]

    if len(sample_ids) != len(values):
        raise RuntimeError(
            f"Sample/PAM50 length mismatch: {len(sample_ids)} vs {len(values)}"
        )
    manifest = list(zip(sample_ids, values))

manifest = [(sid, pam) for sid, pam in manifest if sid]
ids = [sid for sid, _ in manifest]

if len(ids) != len(set(ids)):
    dup = [x for x, n in Counter(ids).items() if n > 1]
    raise RuntimeError(f"Duplicate Xena sample IDs: {dup[:20]}")

counts = Counter(pam for _, pam in manifest if pam)
canonical_counts = {k: counts.get(k, 0) for k in ["LumA", "LumB", "Basal", "Her2", "Normal"]}
other_counts = {k: v for k, v in counts.items() if k not in VALID}

with open(safe_path, "w", encoding="utf-8", newline="") as out:
    w = csv.writer(out, delimiter="\t")
    w.writerow(["xena_sample_id", "pam50_subtype"])
    w.writerows(manifest)

with open(audit_path, "w", encoding="utf-8") as out:
    out.write(f"source_file={raw_path}\n")
    out.write(f"layout={layout}\n")
    out.write(f"field={PAM_FIELD}\n")
    out.write(f"rows_with_sample_id={len(manifest)}\n")
    out.write(f"nonmissing_pam50={sum(counts.values())}\n")
    out.write(f"canonical_counts={canonical_counts}\n")
    out.write(f"other_nonempty_labels={other_counts}\n")

print("Xena layout:", layout)
print("Rows with sample ID:", len(manifest))
print("Nonmissing PAM50:", sum(counts.values()))
print("Canonical PAM50 counts:", canonical_counts)
print("Other nonempty labels:", other_counts)
'@

$tmpPy = Join-Path $env:TEMP "extract_tcga_pam50_xena_v2.py"
Set-Content -Path $tmpPy -Value $py -Encoding UTF8

& python $tmpPy $RawFile $SafeFile $AuditFile
$extractExit = $LASTEXITCODE
Remove-Item $tmpPy -Force -ErrorAction SilentlyContinue

if ($extractExit -ne 0) {
    throw "PAM50 extraction failed."
}

Write-Host ""
Write-Host "[3/4] Sanity-checking plain-text file ..."
$check = @'
import os, sys
p = sys.argv[1]
with open(p, "rb") as f:
    head = f.read(4096)
if not head:
    raise SystemExit("File is empty")
if b"PAM50Call_RNAseq" not in head and os.path.getsize(p) < 5000:
    raise SystemExit("File does not look like a valid Xena clinical matrix")
print("TEXT_FILE_OK bytes=", os.path.getsize(p))
'@

$tmpCheck = Join-Path $env:TEMP "check_tcga_pam50_text_v2.py"
Set-Content -Path $tmpCheck -Value $check -Encoding UTF8
& python $tmpCheck $RawFile
$checkExit = $LASTEXITCODE
Remove-Item $tmpCheck -Force -ErrorAction SilentlyContinue

if ($checkExit -ne 0) {
    throw "Plain-text Xena sanity check failed."
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
Write-Host "04c v2 TCGA PAM50 DOWNLOAD / SAFE EXTRACTION: PASS"
Write-Host "===================================================================================================="
Write-Host "Endpoint used:        $usedUrl"
Write-Host "Raw phenotype file:  $RawFile"
Write-Host "Safe PAM50 manifest: $SafeFile"
Write-Host "Extraction audit:    $AuditFile"
Write-Host "SHA-256:             $HashFile"
Write-Host ""
Write-Host "IMPORTANT: downstream subtype analysis should read ONLY the safe PAM50 manifest."
Write-Host "===================================================================================================="
