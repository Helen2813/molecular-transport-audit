from __future__ import annotations

from pathlib import Path

SCRIPT_VERSION = "04g-probe-pretarget-audit-inputs-v1-no-cli"
DATA_ROOT = Path(r"D:\paper4_tcbb_data")

PATTERNS = {
    "TCGA continuous expression": [
        "**/data_mrna_seq_v2_rsem.txt",
    ],
    "SCAN-B gene expression": [
        "**/GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv.gz",
    ],
    "METABRIC continuous expression": [
        "**/data_mrna_illumina_microarray.txt",
    ],
    "SCAN-B pairing/primary manifest": [
        "**/*scanb*primary*profiles*frozen*.tsv",
        "**/*scanb*pair*.tsv",
    ],
    "Frozen TCGA source universe": [
        "**/*source*universe*.tsv",
        "**/*tcga*10000*.tsv",
    ],
    "Frozen source module membership": [
        "**/*module*membership*.tsv",
        "**/*module*membership*.csv",
    ],
    "Frozen source program/loadings": [
        "**/*source*program*.tsv",
        "**/*loading*.tsv",
        "**/*program*summary*.tsv",
    ],
    "Frozen source edge vectors": [
        "**/*edge*.tsv",
        "**/*edge*.csv",
        "**/*edge*.npy",
    ],
    "SCAN-B hg38 GTF": [
        "**/*.gtf",
        "**/*.gtf.gz",
    ],
}


def main() -> None:
    print("=" * 118)
    print("Paper 4 / TCBB - probe file inputs needed for the final pre-target audits")
    print("=" * 118)
    print(f"Script version: {SCRIPT_VERSION}")
    print()
    print("Scientific guard:")
    print("  Expression values read:                  NO")
    print("  Correlations/PCA/preservation computed:  NO")
    print("  Operation: filesystem discovery only")
    print("=" * 118)

    if not DATA_ROOT.exists():
        raise FileNotFoundError(f"Data root does not exist: {DATA_ROOT}")

    for label, patterns in PATTERNS.items():
        found = []
        seen = set()
        for pattern in patterns:
            for p in DATA_ROOT.glob(pattern):
                if p.is_file():
                    key = str(p.resolve()).lower()
                    if key not in seen:
                        seen.add(key)
                        found.append(p)

        found.sort(key=lambda p: (len(str(p)), str(p).lower()))

        print()
        print(f"{label}: {len(found)} candidate(s)")
        for p in found[:20]:
            try:
                size = p.stat().st_size
                if size >= 1024**3:
                    size_txt = f"{size / 1024**3:.2f} GiB"
                elif size >= 1024**2:
                    size_txt = f"{size / 1024**2:.2f} MiB"
                elif size >= 1024:
                    size_txt = f"{size / 1024:.2f} KiB"
                else:
                    size_txt = f"{size} B"
            except OSError:
                size_txt = "size unavailable"
            print(f"  {p}  [{size_txt}]")
        if len(found) > 20:
            print(f"  ... {len(found) - 20} more")

    print()
    print("=" * 118)
    print("04g PRE-TARGET INPUT PROBE: COMPLETE")
    print("=" * 118)
    print("No scientific data values were read.")
    print("Paste this output back before the executable identity/source-stability audit is built.")
    print("=" * 118)


if __name__ == "__main__":
    main()
