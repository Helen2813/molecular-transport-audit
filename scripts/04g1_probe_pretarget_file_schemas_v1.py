from __future__ import annotations

import csv
import gzip
from pathlib import Path

import pandas as pd


SCRIPT_VERSION = "04g1-probe-pretarget-file-schemas-v1-no-cli"

DATA_ROOT = Path(r"D:\paper4_tcbb_data")

UNIVERSE = (
    DATA_ROOT
    / "paper4_tcbb_tcga_source_universe_v1"
    / "tcga_source_gene_universe_frozen_v1.tsv"
)
MEMBERSHIP = (
    DATA_ROOT
    / "paper4_tcbb_tcga_source_modules_v1"
    / "tcga_source_module_membership_frozen_v1.tsv"
)
PROGRAM_DIR = DATA_ROOT / "paper4_tcbb_frozen_source_programs_v1"
WEIGHTS = PROGRAM_DIR / "tcga_frozen_source_program_weights_v1.tsv"
SUMMARY = PROGRAM_DIR / "tcga_frozen_source_program_summary_v1.tsv"

PAIR_DIR = DATA_ROOT / "paper4_tcbb_scanb_pairing_contract_v1"
TECH_PAIRS = PAIR_DIR / "scanb_technical_replicate_pairs_v1.tsv"
PRIMARY = PAIR_DIR / "scanb_primary_profiles_frozen_v1.tsv"

TCGA_EXPR = (
    DATA_ROOT
    / "paper4_tcbb_input_audit_v1"
    / "staged_continuous_inputs"
    / "TCGA_BRCA_PanCanAtlas2018"
    / "data_mrna_seq_v2_rsem.txt"
)
MET_EXPR = (
    DATA_ROOT
    / "paper4_tcbb_input_audit_v1"
    / "staged_continuous_inputs"
    / "METABRIC"
    / "data_mrna_illumina_microarray.txt"
)
SCANB_EXPR = (
    DATA_ROOT
    / "SCANB_GSE96058"
    / "GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv.gz"
)

GTF = (
    DATA_ROOT
    / "SCANB_GSE96058"
    / "GSE96058_UCSC_hg38_knownGenes_22sep2014.gtf.gz"
)


def require(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Required file does not exist: {path}")


def human_size(n: int) -> str:
    if n >= 1024**3:
        return f"{n / 1024**3:.2f} GiB"
    if n >= 1024**2:
        return f"{n / 1024**2:.2f} MiB"
    if n >= 1024:
        return f"{n / 1024:.2f} KiB"
    return f"{n} B"


def print_tsv_schema(label: str, path: Path, nrows: int = 3) -> None:
    require(path)
    df = pd.read_csv(path, sep="\t", dtype=str, nrows=nrows).fillna("")
    print()
    print(f"{label}")
    print(f"  path: {path}")
    print(f"  columns ({len(df.columns)}): {list(df.columns)}")
    print(f"  first {len(df)} row(s):")
    for rec in df.to_dict(orient="records"):
        print(f"    {rec}")


def read_delimited_header(path: Path, delimiter: str, gz: bool = False) -> list[str]:
    opener = gzip.open if gz else open
    with opener(path, "rt", encoding="utf-8-sig", errors="replace", newline="") as f:
        reader = csv.reader(f, delimiter=delimiter)
        return next(reader)


def print_expr_header(label: str, path: Path, delimiter: str, gz: bool = False) -> None:
    require(path)
    header = read_delimited_header(path, delimiter=delimiter, gz=gz)
    print()
    print(f"{label}")
    print(f"  path: {path}")
    print(f"  total columns: {len(header)}")
    print(f"  first 8 columns: {header[:8]}")
    print(f"  last 5 columns:  {header[-5:]}")


def main() -> None:
    print("=" * 126)
    print("Paper 4 / TCBB - detailed schema probe for final pre-target executable audits")
    print("=" * 126)
    print(f"Script version: {SCRIPT_VERSION}")
    print()
    print("Scientific guard:")
    print("  Expression VALUES read:                         NO")
    print("  Target correlations/PCA/preservation computed: NO")
    print("  Source preservation computed:                  NO")
    print("  Operation: file listing + headers + metadata rows only")
    print("=" * 126)

    # Exact schemas needed by the executable audit.
    print_tsv_schema("Frozen TCGA source universe", UNIVERSE)
    print_tsv_schema("Frozen source module membership", MEMBERSHIP)
    print_tsv_schema("Frozen source program weights/loadings", WEIGHTS)
    print_tsv_schema("Frozen source program summary", SUMMARY)
    print_tsv_schema("SCAN-B technical replicate pairs", TECH_PAIRS, nrows=5)
    print_tsv_schema("SCAN-B frozen primary manifest", PRIMARY, nrows=5)

    print_expr_header("TCGA continuous expression header", TCGA_EXPR, "\t")
    print_expr_header("METABRIC continuous expression header", MET_EXPR, "\t")
    print_expr_header("SCAN-B expression header", SCANB_EXPR, ",", gz=True)

    print()
    print("Frozen-program directory contents:")
    require(PROGRAM_DIR)
    files = sorted([p for p in PROGRAM_DIR.rglob("*") if p.is_file()])
    for p in files:
        print(f"  {p}  [{human_size(p.stat().st_size)}]")

    print()
    print("Searching for exact frozen source matrix / edge-like artifacts ...")
    patterns = [
        "**/tcga_source_datExpr_1082x10000_v1.rds",
        "**/*edge*",
        "**/*edges*",
        "**/*.rds",
        "**/*.npz",
        "**/*.npy",
        "**/*.parquet",
        "**/*.feather",
    ]
    found = []
    seen = set()
    for pat in patterns:
        for p in DATA_ROOT.glob(pat):
            if p.is_file():
                key = str(p.resolve()).lower()
                if key not in seen:
                    seen.add(key)
                    found.append(p)
    found.sort(key=lambda p: str(p).lower())
    if not found:
        print("  NONE")
    else:
        for p in found:
            print(f"  {p}  [{human_size(p.stat().st_size)}]")

    print()
    print("GTF structural check:")
    require(GTF)
    first_data = None
    with gzip.open(GTF, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            first_data = line.rstrip("\n")
            break
    if first_data is None:
        raise RuntimeError("No data line found in GTF.")
    fields = first_data.split("\t")
    print(f"  path: {GTF}")
    print(f"  field count: {len(fields)}")
    print(f"  first data feature: {fields[:8]}")
    print(f"  attributes preview: {fields[8][:500] if len(fields) > 8 else ''}")

    print()
    print("=" * 126)
    print("04g1 PRE-TARGET FILE-SCHEMA PROBE: COMPLETE")
    print("=" * 126)
    print("No expression values and no target preservation statistics were read/calculated.")
    print("Paste this output back; it is sufficient to build the executable identity/source-stability audit.")
    print("=" * 126)


if __name__ == "__main__":
    main()
