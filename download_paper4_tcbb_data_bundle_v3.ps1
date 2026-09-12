param(
    [string]$OutDir = (Join-Path (Get-Location) "paper4_tcbb_data_bundle_v3"),
    [switch]$FullBackup,
    [switch]$IncludeTranscriptLevel
)

$ErrorActionPreference = "Stop"
$SCRIPT_VERSION = "download-paper4-tcbb-data-bundle-v3"

Write-Host "======================================================================="
Write-Host "Paper 4 / TCBB - offline dataset bundle"
Write-Host "Script version: $SCRIPT_VERSION"
Write-Host "Output: $OutDir"
Write-Host "======================================================================="

New-Item -ItemType Directory -Force -Path $OutDir | Out-Null
$ScanBDir = Join-Path $OutDir "SCANB_GSE96058"
$StudyDir = Join-Path $OutDir "cBioPortal_study_archives"
New-Item -ItemType Directory -Force -Path $ScanBDir | Out-Null
New-Item -ItemType Directory -Force -Path $StudyDir | Out-Null

if (-not (Get-Command curl.exe -ErrorAction SilentlyContinue)) {
    throw "curl.exe was not found. Windows 10/11 normally includes it."
}

function Download-Resume {
    param(
        [Parameter(Mandatory=$true)][string]$Url,
        [Parameter(Mandatory=$true)][string]$Dest
    )

    Write-Host ""
    Write-Host "-----------------------------------------------------------------------"
    Write-Host "URL : $Url"
    Write-Host "FILE: $Dest"

    if (Test-Path $Dest) {
        $existing = (Get-Item $Dest).Length
        Write-Host ("Existing file: {0:N2} MB; attempting resume." -f ($existing/1MB))
    }

    & curl.exe -L --fail --retry 8 --retry-all-errors --retry-delay 5 `
        --connect-timeout 30 -C - -o $Dest $Url

    # curl 33 = server refused range resume; retry cleanly.
    if ($LASTEXITCODE -eq 33) {
        Write-Warning "Server refused resume; restarting this file."
        Remove-Item $Dest -Force -ErrorAction SilentlyContinue
        & curl.exe -L --fail --retry 8 --retry-all-errors --retry-delay 5 `
            --connect-timeout 30 -o $Dest $Url
    }

    if ($LASTEXITCODE -ne 0) {
        throw "Download failed (curl exit $LASTEXITCODE): $Url"
    }

    $f = Get-Item $Dest
    if ($f.Length -le 0) { throw "Downloaded file is empty: $Dest" }
    Write-Host ("Saved: {0:N2} MB" -f ($f.Length/1MB))
}

# ==========================================================================
# REQUIRED NEW COHORT: SCAN-B / GSE96058
# ==========================================================================

$geoBase = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE96nnn/GSE96058"

Download-Resume `
    -Url "$geoBase/suppl/GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv.gz" `
    -Dest (Join-Path $ScanBDir "GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv.gz")

# Both platform-specific series matrices: tiny, but important for technical-only manifest.
Download-Resume `
    -Url "$geoBase/matrix/GSE96058-GPL11154_series_matrix.txt.gz" `
    -Dest (Join-Path $ScanBDir "GSE96058-GPL11154_series_matrix.txt.gz")

Download-Resume `
    -Url "$geoBase/matrix/GSE96058-GPL18573_series_matrix.txt.gz" `
    -Dest (Join-Path $ScanBDir "GSE96058-GPL18573_series_matrix.txt.gz")

# Metadata backup. Useful if series-matrix parsing misses a sample-level field.
Download-Resume `
    -Url "$geoBase/soft/GSE96058_family.soft.gz" `
    -Dest (Join-Path $ScanBDir "GSE96058_family.soft.gz")

# Small genome annotation distributed with the study; not required for the
# primary gene-level analysis, but cheap insurance for offline work.
Download-Resume `
    -Url "$geoBase/suppl/GSE96058_UCSC_hg38_knownGenes_22sep2014.gtf.gz" `
    -Dest (Join-Path $ScanBDir "GSE96058_UCSC_hg38_knownGenes_22sep2014.gtf.gz")

# Transcript-level matrix is deliberately optional (~820.8 MB); gene-level
# data are sufficient for the planned MTA analysis.
if ($IncludeTranscriptLevel) {
    Download-Resume `
        -Url "$geoBase/suppl/GSE96058_transcript_expression_3273_samples_and_136_replicates.csv.gz" `
        -Dest (Join-Path $ScanBDir "GSE96058_transcript_expression_3273_samples_and_136_replicates.csv.gz")
}

# ==========================================================================
# OPTIONAL FULL OFFLINE BACKUP
# ==========================================================================
# These archives are NOT necessary to re-download if you already have the
# correct continuous matrices elsewhere. -FullBackup is only for insurance
# before losing broadband access.
# ==========================================================================

if ($FullBackup) {
    $cbio = "https://datahub.assets.cbioportal.org"

    Download-Resume `
        -Url "$cbio/brca_tcga_pan_can_atlas_2018.tar.gz" `
        -Dest (Join-Path $StudyDir "brca_tcga_pan_can_atlas_2018.tar.gz")

    Download-Resume `
        -Url "$cbio/brca_metabric.tar.gz" `
        -Dest (Join-Path $StudyDir "brca_metabric.tar.gz")

    # Optional cross-tissue stress controls. They are Supplement-only candidates,
    # not required for the main TCBB design.
    Download-Resume `
        -Url "$cbio/laml_tcga_pan_can_atlas_2018.tar.gz" `
        -Dest (Join-Path $StudyDir "laml_tcga_pan_can_atlas_2018.tar.gz")

    Download-Resume `
        -Url "$cbio/luad_tcga_pan_can_atlas_2018.tar.gz" `
        -Dest (Join-Path $StudyDir "luad_tcga_pan_can_atlas_2018.tar.gz")
}

# ==========================================================================
# INTEGRITY CHECKS
# ==========================================================================

$python = Get-Command python -ErrorAction SilentlyContinue
if ($python) {
    $py = @'
import csv, gzip, os, sys

root = sys.argv[1]
expr = os.path.join(root, "SCANB_GSE96058",
                    "GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv.gz")

def read_to_eof(path):
    n = 0
    with gzip.open(path, "rb") as f:
        while True:
            chunk = f.read(8 * 1024 * 1024)
            if not chunk:
                break
            n += len(chunk)
    print(f"GZIP_OK\t{os.path.basename(path)}\tuncompressed_bytes={n}")

for fn in os.listdir(os.path.join(root, "SCANB_GSE96058")):
    if fn.endswith(".gz") and fn != os.path.basename(expr):
        read_to_eof(os.path.join(root, "SCANB_GSE96058", fn))

rows = 0
with gzip.open(expr, "rt", newline="", encoding="utf-8-sig") as f:
    r = csv.reader(f)
    header = next(r)
    for _ in r:
        rows += 1
print(f"SCANB_EXPRESSION_OK\trows={rows}\tcolumns={len(header)}")
'@
    $tmp = Join-Path $OutDir "_integrity_check_v3.py"
    Set-Content -Path $tmp -Value $py -Encoding UTF8
    & python $tmp $OutDir
    if ($LASTEXITCODE -ne 0) {
        throw "SCAN-B gzip/CSV integrity check failed."
    }
    Remove-Item $tmp -Force
} else {
    Write-Warning "Python not found: gzip CRC / row-count validation skipped."
}

# Check tar archives without extracting if tar.exe exists.
if ($FullBackup -and (Get-Command tar.exe -ErrorAction SilentlyContinue)) {
    Get-ChildItem $StudyDir -Filter "*.tar.gz" | ForEach-Object {
        Write-Host "Checking archive: $($_.Name)"
        & tar.exe -tzf $_.FullName *> $null
        if ($LASTEXITCODE -ne 0) { throw "Corrupt tar.gz archive: $($_.FullName)" }
    }
}

# SHA-256 manifest for every downloaded file.
$hashFile = Join-Path $OutDir "SHA256SUMS.txt"
Remove-Item $hashFile -ErrorAction SilentlyContinue
Get-ChildItem $OutDir -Recurse -File |
    Where-Object { $_.Name -notin @("SHA256SUMS.txt", "_integrity_check_v3.py") } |
    Sort-Object FullName |
    ForEach-Object {
        $hash = (Get-FileHash -Algorithm SHA256 $_.FullName).Hash.ToLowerInvariant()
        $rel = $_.FullName.Substring($OutDir.Length).TrimStart('\')
        "$hash  $rel" | Add-Content -Encoding UTF8 $hashFile
    }

$readme = @"
Paper 4 / TCBB offline dataset bundle
Script version: $SCRIPT_VERSION

REQUIRED PRIMARY FILES
- SCAN-B gene-level expression:
  GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv.gz
- SCAN-B series matrices:
  GSE96058-GPL11154_series_matrix.txt.gz
  GSE96058-GPL18573_series_matrix.txt.gz

EXPECTED EXISTING cBioPortal CONTINUOUS FILES
- TCGA-BRCA PanCancer Atlas:
  data_mrna_seq_v2_rsem.txt
- METABRIC:
  data_mrna_illumina_microarray.txt

DO NOT SUBSTITUTE Z-SCORE MATRICES.

OPTIONAL FULL BACKUP ARCHIVES
- brca_tcga_pan_can_atlas_2018.tar.gz
- brca_metabric.tar.gz
- laml_tcga_pan_can_atlas_2018.tar.gz
- luad_tcga_pan_can_atlas_2018.tar.gz

NOT DOWNLOADED AUTOMATICALLY
- MSigDB Hallmark GMT: requires free MSigDB login/registration.
- TARGET-OS, GSE21257, GSE39055: reserved for Paper 6; do not use in Paper 4.

The transcript-level SCAN-B matrix is not needed for the planned analysis.
"@
Set-Content -Path (Join-Path $OutDir "README_DATA_BUNDLE.txt") -Value $readme -Encoding UTF8

Write-Host ""
Write-Host "======================================================================="
Write-Host "DONE"
Write-Host "Bundle:  $OutDir"
Write-Host "Hashes:  $hashFile"
Write-Host "======================================================================="
