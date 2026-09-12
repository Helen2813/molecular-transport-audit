from __future__ import annotations

import csv
import gzip
import hashlib
import json
import math
import shutil
import tarfile
from pathlib import Path

import numpy as np
import pandas as pd


SCRIPT_VERSION = "01a-stage-and-audit-tcbb-inputs-v1-no-cli"

DATA_ROOT = Path(r"D:\paper4_tcbb_data")

SCANB_DIR = DATA_ROOT / "SCANB_GSE96058"
ARCHIVE_DIR = DATA_ROOT / "cBioPortal_study_archives"
OUT_DIR = DATA_ROOT / "paper4_tcbb_input_audit_v1"
STAGED_DIR = OUT_DIR / "staged_continuous_inputs"

SCANB_EXPR = (
    SCANB_DIR
    / "GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv.gz"
)
SCANB_SERIES = [
    SCANB_DIR / "GSE96058-GPL11154_series_matrix.txt.gz",
    SCANB_DIR / "GSE96058-GPL18573_series_matrix.txt.gz",
]

TCGA_ARCHIVE = ARCHIVE_DIR / "brca_tcga_pan_can_atlas_2018.tar.gz"
METABRIC_ARCHIVE = ARCHIVE_DIR / "brca_metabric.tar.gz"

TCGA_EXPR_NAME = "data_mrna_seq_v2_rsem.txt"
TCGA_META_NAME = "meta_mrna_seq_v2_rsem.txt"

METABRIC_EXPR_NAME = "data_mrna_illumina_microarray.txt"
METABRIC_META_NAME = "meta_mrna_illumina_microarray.txt"

CHUNK_ROWS = 256
GLOBAL_VALUE_SAMPLE_CAP = 1_000_000


def banner() -> None:
    print("=" * 112)
    print("Paper 4 / TCBB - stage and audit external-input matrices BEFORE preservation analysis")
    print("=" * 112)
    print(f"Script version: {SCRIPT_VERSION}")
    print(f"Data root:      {DATA_ROOT}")
    print(f"Output dir:     {OUT_DIR}")
    print()
    print("Safety / scientific contract:")
    print("  Target clinical outcomes loaded:                         NO")
    print("  Clinical patient/sample files opened:                    NO")
    print("  Target preservation correlations calculated:             NO")
    print("  Target edge/loading preservation statistics calculated:  NO")
    print("  Target-driven module selection/tuning performed:          NO")
    print("  Allowed input QC: schema/range/missingness/duplicates/platform/replicate titles: YES")
    print("=" * 112)


def require(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Required file does not exist: {path}")


def sha256_file(path: Path, chunk_size: int = 8 * 1024 * 1024) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        while True:
            block = f.read(chunk_size)
            if not block:
                break
            h.update(block)
    return h.hexdigest()


def extract_member_by_basename(archive: Path, basename: str, dest: Path) -> None:
    require(archive)
    dest.parent.mkdir(parents=True, exist_ok=True)

    if dest.exists() and dest.stat().st_size > 0:
        print(f"STAGED EXISTS: {dest}")
        return

    with tarfile.open(archive, "r:gz") as tf:
        matches = [
            m for m in tf.getmembers()
            if m.isfile() and Path(m.name).name == basename
        ]
        if len(matches) != 1:
            names = [m.name for m in matches]
            raise RuntimeError(
                f"Expected exactly one '{basename}' in {archive.name}; "
                f"found {len(matches)}: {names}"
            )
        member = matches[0]
        src = tf.extractfile(member)
        if src is None:
            raise RuntimeError(f"Could not read {member.name} from {archive}")
        with dest.open("wb") as out:
            shutil.copyfileobj(src, out, length=8 * 1024 * 1024)

    print(f"STAGED: {archive.name} :: {basename} -> {dest}")


def parse_cbio_meta(path: Path) -> dict[str, str]:
    require(path)
    out: dict[str, str] = {}
    with path.open("rt", encoding="utf-8-sig", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#") or ":" not in line:
                continue
            k, v = line.split(":", 1)
            out[k.strip()] = v.strip()
    return out


def _quantile_dict(values: np.ndarray) -> dict[str, float | None]:
    values = np.asarray(values, dtype=np.float64)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return {k: None for k in ("min", "q01", "q25", "median", "q75", "q99", "max")}
    q = np.quantile(values, [0, 0.01, 0.25, 0.5, 0.75, 0.99, 1.0])
    return {
        "min": float(q[0]),
        "q01": float(q[1]),
        "q25": float(q[2]),
        "median": float(q[3]),
        "q75": float(q[4]),
        "q99": float(q[5]),
        "max": float(q[6]),
    }


def audit_matrix(
    *,
    label: str,
    path: Path,
    sep: str,
    id_columns: int,
    compression: str | None = None,
) -> dict:
    require(path)

    header_df = pd.read_csv(
        path,
        sep=sep,
        compression=compression,
        nrows=0,
        low_memory=False,
    )
    columns = list(header_df.columns)
    if len(columns) <= id_columns:
        raise RuntimeError(f"{label}: no numeric/sample columns detected.")

    id_cols = columns[:id_columns]
    sample_cols = columns[id_columns:]
    duplicate_sample_columns = int(pd.Index(sample_cols).duplicated().sum())

    n_rows = 0
    total_values = 0
    missing_values = 0
    negative_values = 0
    zero_values = 0
    global_min = math.inf
    global_max = -math.inf

    first_id_seen: set[str] = set()
    duplicate_first_ids = 0
    blank_first_ids = 0

    gene_means: list[np.ndarray] = []
    gene_vars: list[np.ndarray] = []
    near_zero_var_genes = 0

    sampled_values: list[np.ndarray] = []
    sampled_value_count = 0

    reader = pd.read_csv(
        path,
        sep=sep,
        compression=compression,
        chunksize=CHUNK_ROWS,
        low_memory=False,
    )

    for chunk_idx, df in enumerate(reader, start=1):
        n_rows += len(df)

        first_ids = df[id_cols[0]].astype("string")
        for raw in first_ids:
            value = "" if pd.isna(raw) else str(raw).strip()
            if not value:
                blank_first_ids += 1
            if value in first_id_seen:
                duplicate_first_ids += 1
            else:
                first_id_seen.add(value)

        numeric = df[sample_cols].apply(pd.to_numeric, errors="coerce")
        arr = numeric.to_numpy(dtype=np.float64, copy=False)

        finite_mask = np.isfinite(arr)
        finite_n = int(finite_mask.sum())
        this_total = int(arr.size)
        total_values += this_total
        missing_values += this_total - finite_n

        if finite_n:
            finite = arr[finite_mask]
            negative_values += int((finite < 0).sum())
            zero_values += int((finite == 0).sum())
            global_min = min(global_min, float(finite.min()))
            global_max = max(global_max, float(finite.max()))

            if sampled_value_count < GLOBAL_VALUE_SAMPLE_CAP:
                remaining = GLOBAL_VALUE_SAMPLE_CAP - sampled_value_count
                take_target = min(5000, remaining)
                stride = max(1, finite.size // max(1, take_target))
                take = finite[::stride][:take_target].astype(np.float64, copy=True)
                sampled_values.append(take)
                sampled_value_count += int(take.size)

        with np.errstate(invalid="ignore", divide="ignore"):
            row_means = np.nanmean(arr, axis=1)
            row_vars = np.nanvar(arr, axis=1, ddof=1)

        gene_means.append(row_means.astype(np.float64, copy=False))
        gene_vars.append(row_vars.astype(np.float64, copy=False))
        near_zero_var_genes += int(np.nansum(row_vars <= 1e-12))

        if chunk_idx % 25 == 0:
            print(f"  {label}: audited {n_rows:,} rows ...")

    all_gene_means = np.concatenate(gene_means) if gene_means else np.array([], dtype=float)
    all_gene_vars = np.concatenate(gene_vars) if gene_vars else np.array([], dtype=float)
    global_sample = (
        np.concatenate(sampled_values) if sampled_values else np.array([], dtype=float)
    )

    finite_total = total_values - missing_values
    result = {
        "dataset": label,
        "path": str(path),
        "sha256": sha256_file(path),
        "rows_genes": int(n_rows),
        "id_columns": id_cols,
        "sample_columns": int(len(sample_cols)),
        "duplicate_sample_columns": duplicate_sample_columns,
        "first_id_column": id_cols[0],
        "duplicate_first_ids": int(duplicate_first_ids),
        "blank_first_ids": int(blank_first_ids),
        "numeric_values_total": int(total_values),
        "numeric_values_finite": int(finite_total),
        "missing_fraction": float(missing_values / total_values) if total_values else None,
        "negative_fraction_of_finite": (
            float(negative_values / finite_total) if finite_total else None
        ),
        "zero_fraction_of_finite": (
            float(zero_values / finite_total) if finite_total else None
        ),
        "global_min": None if math.isinf(global_min) else float(global_min),
        "global_max": None if math.isinf(global_max) else float(global_max),
        "sampled_value_quantiles": _quantile_dict(global_sample),
        "gene_mean_quantiles": _quantile_dict(all_gene_means),
        "gene_variance_quantiles": _quantile_dict(all_gene_vars),
        "near_zero_variance_genes": int(near_zero_var_genes),
    }
    return result


def _parse_geo_header_values(payload: str) -> list[str]:
    return next(csv.reader([payload], delimiter="\t", quotechar='"'))


def parse_scanb_technical_inventory(series_paths: list[Path], out_tsv: Path) -> dict:
    """
    Read ONLY GEO sample title/accession/platform header rows.
    Deliberately ignores !Sample_characteristics_ch1 and all clinical fields.
    """
    rows: list[dict[str, str]] = []

    for path in series_paths:
        require(path)
        fields: dict[str, list[str]] = {}
        with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
            for line in f:
                if line.startswith("!Sample_title"):
                    fields["title"] = _parse_geo_header_values(line.split("\t", 1)[1])
                elif line.startswith("!Sample_geo_accession"):
                    fields["geo_accession"] = _parse_geo_header_values(line.split("\t", 1)[1])
                elif line.startswith("!Sample_platform_id"):
                    fields["platform_id"] = _parse_geo_header_values(line.split("\t", 1)[1])

                if line.startswith("!series_matrix_table_begin"):
                    break

        required_keys = {"title", "geo_accession", "platform_id"}
        missing = required_keys.difference(fields)
        if missing:
            raise RuntimeError(f"{path.name}: missing safe technical header fields: {sorted(missing)}")

        lengths = {k: len(v) for k, v in fields.items()}
        if len(set(lengths.values())) != 1:
            raise RuntimeError(f"{path.name}: inconsistent sample header lengths: {lengths}")

        n = lengths["geo_accession"]
        for i in range(n):
            rows.append(
                {
                    "geo_accession": fields["geo_accession"][i].strip('"'),
                    "sample_title": fields["title"][i].strip('"'),
                    "platform_id": fields["platform_id"][i].strip('"'),
                    "source_series_matrix": path.name,
                }
            )

    df = pd.DataFrame(rows)
    df.to_csv(out_tsv, sep="\t", index=False)

    duplicate_geo = int(df["geo_accession"].duplicated().sum())
    title_lower = df["sample_title"].str.lower()
    apparent_replicates = int(
        title_lower.str.contains("repl|replicate|repeat|resequenc", regex=True, na=False).sum()
    )

    return {
        "rows": int(len(df)),
        "unique_geo_accessions": int(df["geo_accession"].nunique()),
        "duplicate_geo_accessions": duplicate_geo,
        "platform_counts": {
            str(k): int(v)
            for k, v in df["platform_id"].value_counts(dropna=False).to_dict().items()
        },
        "titles_with_apparent_replicate_keyword": apparent_replicates,
        "output_tsv": str(out_tsv),
        "clinical_characteristics_read": False,
    }


def main() -> None:
    banner()

    for path in [SCANB_EXPR, *SCANB_SERIES, TCGA_ARCHIVE, METABRIC_ARCHIVE]:
        require(path)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    STAGED_DIR.mkdir(parents=True, exist_ok=True)

    tcga_dir = STAGED_DIR / "TCGA_BRCA_PanCanAtlas2018"
    metabric_dir = STAGED_DIR / "METABRIC"

    tcga_expr = tcga_dir / TCGA_EXPR_NAME
    tcga_meta = tcga_dir / TCGA_META_NAME
    metabric_expr = metabric_dir / METABRIC_EXPR_NAME
    metabric_meta = metabric_dir / METABRIC_META_NAME

    print("\n[1/4] Staging only required continuous-expression files and expression metadata ...")
    extract_member_by_basename(TCGA_ARCHIVE, TCGA_EXPR_NAME, tcga_expr)
    extract_member_by_basename(TCGA_ARCHIVE, TCGA_META_NAME, tcga_meta)
    extract_member_by_basename(METABRIC_ARCHIVE, METABRIC_EXPR_NAME, metabric_expr)
    extract_member_by_basename(METABRIC_ARCHIVE, METABRIC_META_NAME, metabric_meta)

    tcga_meta_dict = parse_cbio_meta(tcga_meta)
    metabric_meta_dict = parse_cbio_meta(metabric_meta)

    if tcga_meta_dict.get("datatype", "").upper() != "CONTINUOUS":
        raise RuntimeError(
            f"TCGA expression metadata is not CONTINUOUS: {tcga_meta_dict}"
        )
    if metabric_meta_dict.get("datatype", "").upper() != "CONTINUOUS":
        raise RuntimeError(
            f"METABRIC expression metadata is not CONTINUOUS: {metabric_meta_dict}"
        )

    print("\n[2/4] Auditing matrix schema and expression-scale QC ONLY ...")

    tcga_qc = audit_matrix(
        label="TCGA_BRCA_PanCanAtlas2018",
        path=tcga_expr,
        sep="\t",
        id_columns=2,
        compression=None,
    )
    tcga_qc["cbio_meta"] = tcga_meta_dict

    metabric_qc = audit_matrix(
        label="METABRIC",
        path=metabric_expr,
        sep="\t",
        id_columns=2,
        compression=None,
    )
    metabric_qc["cbio_meta"] = metabric_meta_dict

    scanb_qc = audit_matrix(
        label="SCANB_GSE96058",
        path=SCANB_EXPR,
        sep=",",
        id_columns=1,
        compression="gzip",
    )
    scanb_qc["published_file_role"] = (
        "gene-level transformed expression; primary target input candidate"
    )

    print("\n[3/4] Building SCAN-B technical-only inventory (title/accession/platform only) ...")
    technical_tsv = OUT_DIR / "scanb_technical_inventory_v1.tsv"
    scanb_technical = parse_scanb_technical_inventory(SCANB_SERIES, technical_tsv)

    print("\n[4/4] Writing audit outputs ...")
    all_qc = {
        "script_version": SCRIPT_VERSION,
        "scientific_guard": {
            "target_outcomes_loaded": False,
            "clinical_files_opened": False,
            "target_preservation_statistics_calculated": False,
            "target_driven_tuning": False,
        },
        "TCGA_BRCA": tcga_qc,
        "SCANB_GSE96058": scanb_qc,
        "METABRIC": metabric_qc,
        "SCANB_technical_inventory": scanb_technical,
    }

    json_path = OUT_DIR / "input_qc_summary_v1.json"
    json_path.write_text(json.dumps(all_qc, indent=2), encoding="utf-8")

    table_rows = []
    for item in (tcga_qc, scanb_qc, metabric_qc):
        table_rows.append(
            {
                "dataset": item["dataset"],
                "rows_genes": item["rows_genes"],
                "sample_columns": item["sample_columns"],
                "duplicate_sample_columns": item["duplicate_sample_columns"],
                "duplicate_first_ids": item["duplicate_first_ids"],
                "missing_fraction": item["missing_fraction"],
                "negative_fraction_of_finite": item["negative_fraction_of_finite"],
                "zero_fraction_of_finite": item["zero_fraction_of_finite"],
                "global_min": item["global_min"],
                "global_max": item["global_max"],
                "near_zero_variance_genes": item["near_zero_variance_genes"],
                "sha256": item["sha256"],
                "path": item["path"],
            }
        )

    summary_tsv = OUT_DIR / "input_qc_matrix_summary_v1.tsv"
    pd.DataFrame(table_rows).to_csv(summary_tsv, sep="\t", index=False)

    staged_manifest = OUT_DIR / "staged_input_manifest_v1.tsv"
    staged_rows = []
    for dataset, path in [
        ("TCGA_BRCA_expression", tcga_expr),
        ("TCGA_BRCA_expression_meta", tcga_meta),
        ("METABRIC_expression", metabric_expr),
        ("METABRIC_expression_meta", metabric_meta),
        ("SCANB_expression", SCANB_EXPR),
        ("SCANB_series_GPL11154", SCANB_SERIES[0]),
        ("SCANB_series_GPL18573", SCANB_SERIES[1]),
    ]:
        staged_rows.append(
            {
                "dataset_role": dataset,
                "path": str(path),
                "bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            }
        )
    pd.DataFrame(staged_rows).to_csv(staged_manifest, sep="\t", index=False)

    print("\n" + "=" * 112)
    print("01a INPUT AUDIT: PASS")
    print("=" * 112)
    print(f"TCGA-BRCA: genes={tcga_qc['rows_genes']:,} samples={tcga_qc['sample_columns']:,} "
          f"min={tcga_qc['global_min']} max={tcga_qc['global_max']} "
          f"frac_neg={tcga_qc['negative_fraction_of_finite']:.6g}")
    print(f"SCAN-B:    genes={scanb_qc['rows_genes']:,} samples={scanb_qc['sample_columns']:,} "
          f"min={scanb_qc['global_min']} max={scanb_qc['global_max']} "
          f"frac_neg={scanb_qc['negative_fraction_of_finite']:.6g}")
    print(f"METABRIC:  genes={metabric_qc['rows_genes']:,} samples={metabric_qc['sample_columns']:,} "
          f"min={metabric_qc['global_min']} max={metabric_qc['global_max']} "
          f"frac_neg={metabric_qc['negative_fraction_of_finite']:.6g}")
    print()
    print(f"SCAN-B technical inventory: {scanb_technical['rows']} profiles")
    print(f"  platform counts: {scanb_technical['platform_counts']}")
    print(f"  titles containing replicate-like keyword: "
          f"{scanb_technical['titles_with_apparent_replicate_keyword']}")
    print()
    print("Outputs:")
    print(f"  {json_path}")
    print(f"  {summary_tsv}")
    print(f"  {technical_tsv}")
    print(f"  {staged_manifest}")
    print()
    print("NO preservation statistics were computed.")
    print("Next scientific step is to review these QC outputs and fill expression/input-related TO_FREEZE fields.")
    print("=" * 112)


if __name__ == "__main__":
    main()
