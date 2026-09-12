from __future__ import annotations

import csv
import gzip
import io
import json
import os
import re
import tarfile
from collections import Counter
from pathlib import Path
from typing import Iterable

import pandas as pd


SCRIPT_VERSION = "04b-audit-breast-subtype-metadata-v1-no-cli"

DATA_ROOT = Path(r"D:\paper4_tcbb_data")
ARCHIVE_ROOT = DATA_ROOT / "cBioPortal_study_archives"

TCGA_ARCHIVE = ARCHIVE_ROOT / "brca_tcga_pan_can_atlas_2018.tar.gz"
METABRIC_ARCHIVE = ARCHIVE_ROOT / "brca_metabric.tar.gz"

TCGA_EXPR = (
    DATA_ROOT
    / "paper4_tcbb_input_audit_v1"
    / "staged_continuous_inputs"
    / "TCGA_BRCA_PanCanAtlas2018"
    / "data_mrna_seq_v2_rsem.txt"
)
METABRIC_EXPR = (
    DATA_ROOT
    / "paper4_tcbb_input_audit_v1"
    / "staged_continuous_inputs"
    / "METABRIC"
    / "data_mrna_illumina_microarray.txt"
)

SCANB_SERIES = [
    DATA_ROOT / "SCANB_GSE96058" / "GSE96058-GPL11154_series_matrix.txt.gz",
    DATA_ROOT / "SCANB_GSE96058" / "GSE96058-GPL18573_series_matrix.txt.gz",
]
SCANB_PRIMARY = (
    DATA_ROOT
    / "paper4_tcbb_scanb_pairing_contract_v1"
    / "scanb_primary_profiles_frozen_v1.tsv"
)

OUT_DIR = DATA_ROOT / "paper4_tcbb_breast_subtype_audit_v1"

# Search likely thesis locations first; if absent, scan Desktop for clinical.tsv/csv.
LIKELY_THESIS_FILES = [
    Path(r"C:\Users\olegk\Desktop\Thesis_v3\data\drags\clinical.tsv"),
    Path(r"C:\Users\olegk\Desktop\Thesis_v3\data\drugs\clinical.tsv"),
    Path(r"C:\Users\olegk\Desktop\thesis_v3\data\drags\clinical.tsv"),
    Path(r"C:\Users\olegk\Desktop\thesis_v3\data\drugs\clinical.tsv"),
    Path(r"C:\Users\olegk\Desktop\Thesis\data\drags\clinical.tsv"),
    Path(r"C:\Users\olegk\Desktop\Thesis\data\drugs\clinical.tsv"),
]
DESKTOP_ROOT = Path(r"C:\Users\olegk\Desktop")

# Explicit PAM50 / intrinsic-subtype signals.
PAM50_PATTERNS = [
    re.compile(r"pam[\s_\-]*50", re.I),
    re.compile(r"intrinsic[\s_\-]*(subtype|cluster|class)", re.I),
    re.compile(r"(subtype|cluster|class)[\s_\-]*intrinsic", re.I),
    re.compile(r"molecular[\s_\-]*subtype", re.I),
]

# Useful backup clinical groupings. These are NOT PAM50 and must never be
# silently substituted for PAM50.
RECEPTOR_PATTERNS = [
    re.compile(r"\bestrogen\b|\ber[\s_\-]*(status|positive|negative)?\b", re.I),
    re.compile(r"\bprogesterone\b|\bpr[\s_\-]*(status|positive|negative)?\b", re.I),
    re.compile(r"\bher[\s_\-]*2\b|\bher2\b", re.I),
    re.compile(r"\btnbc\b|triple[\s_\-]*negative", re.I),
]


def require(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Required file does not exist: {path}")


def norm_value(x: object) -> str:
    if x is None:
        return ""
    s = str(x).strip().strip('"').strip()
    if s.lower() in {"", "na", "nan", "n/a", "not available", "unknown", "[not available]"}:
        return ""
    return s


def is_pam50_label(text: str) -> bool:
    return any(p.search(text or "") for p in PAM50_PATTERNS)


def is_receptor_label(text: str) -> bool:
    return any(p.search(text or "") for p in RECEPTOR_PATTERNS)


def canonicalize_pam50(value: str) -> str:
    """
    Conservative normalization only for reporting counts.
    Does NOT infer PAM50 from receptor status.
    """
    s = norm_value(value)
    if not s:
        return ""
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
        "claudinlow": "Claudin-low",
    }
    return mapping.get(low, s)


def count_levels(values: Iterable[str], canonicalize: bool = False) -> dict:
    vals = []
    for x in values:
        v = norm_value(x)
        if not v:
            continue
        vals.append(canonicalize_pam50(v) if canonicalize else v)
    return dict(Counter(vals).most_common())


def find_thesis_candidates() -> list[Path]:
    found: list[Path] = []
    seen: set[str] = set()

    for p in LIKELY_THESIS_FILES:
        if p.exists():
            key = str(p.resolve()).lower()
            if key not in seen:
                seen.add(key)
                found.append(p)

    if found or not DESKTOP_ROOT.exists():
        return found

    # Conservative fallback search. Prune obvious huge/non-data trees.
    prune_names = {
        ".venv", "venv", "node_modules", ".git", "__pycache__",
        "dist", "build", "target",
    }
    wanted_names = {"clinical.tsv", "clinical.csv"}

    for root, dirs, files in os.walk(DESKTOP_ROOT):
        dirs[:] = [d for d in dirs if d.lower() not in prune_names]
        for f in files:
            fl = f.lower()
            if fl in wanted_names or ("clinical" in fl and fl.endswith((".tsv", ".csv", ".txt"))):
                p = Path(root) / f
                key = str(p.resolve()).lower()
                if key not in seen:
                    seen.add(key)
                    found.append(p)
                    if len(found) >= 50:
                        return found
    return found


def expression_sample_headers(path: Path, sep: str, id_cols: int) -> list[str]:
    df = pd.read_csv(path, sep=sep, nrows=0, low_memory=False)
    return list(df.columns)[id_cols:]


def parse_cbio_clinical_member(
    archive: Path,
    member_basename: str,
    source_label: str,
    expression_samples: list[str] | None = None,
) -> tuple[list[dict], list[pd.DataFrame]]:
    """
    Reads cBioPortal clinical files safely:
      1) header/display/description schema
      2) only ID + candidate PAM50/subtype/receptor columns are retained
    Outcome/treatment columns are never persisted or returned.
    """
    results: list[dict] = []
    safe_tables: list[pd.DataFrame] = []

    with tarfile.open(archive, "r:gz") as tf:
        members = [
            m for m in tf.getmembers()
            if m.isfile() and Path(m.name).name == member_basename
        ]
        if not members:
            return results, safe_tables

        for member in members:
            fh = tf.extractfile(member)
            if fh is None:
                continue

            text = io.TextIOWrapper(fh, encoding="utf-8-sig", errors="replace")
            comment_rows: list[list[str]] = []
            header: list[str] | None = None
            data_lines: list[str] = []

            for line in text:
                stripped = line.rstrip("\r\n")
                if stripped.startswith("#") and header is None:
                    comment_rows.append(stripped[1:].split("\t"))
                    continue
                if header is None:
                    header = stripped.split("\t")
                    continue
                data_lines.append(stripped)

            if not header:
                continue

            # cBioPortal convention usually stores:
            # #Display Name, #Description, #Datatype, #Priority before actual header.
            display = comment_rows[0] if len(comment_rows) >= 1 else [""] * len(header)
            desc = comment_rows[1] if len(comment_rows) >= 2 else [""] * len(header)

            # Pad metadata rows if malformed/short.
            if len(display) < len(header):
                display += [""] * (len(header) - len(display))
            if len(desc) < len(header):
                desc += [""] * (len(header) - len(desc))

            pam_cols = []
            receptor_cols = []
            for i, col in enumerate(header):
                combined = " | ".join([col, display[i] if i < len(display) else "", desc[i] if i < len(desc) else ""])
                if is_pam50_label(combined):
                    pam_cols.append(i)
                elif is_receptor_label(combined):
                    receptor_cols.append(i)

            id_names = {"SAMPLE_ID", "PATIENT_ID", "CASE_ID"}
            id_cols = [i for i, c in enumerate(header) if c.upper() in id_names]
            selected = sorted(set(id_cols + pam_cols + receptor_cols))

            if not (pam_cols or receptor_cols):
                results.append(
                    {
                        "source": source_label,
                        "member": member.name,
                        "pam50_candidate_columns": [],
                        "receptor_candidate_columns": [],
                        "rows": len(data_lines),
                    }
                )
                continue

            rows = []
            for line in data_lines:
                parts = line.split("\t")
                if len(parts) < len(header):
                    parts += [""] * (len(header) - len(parts))
                rows.append({header[i]: parts[i] for i in selected})

            safe = pd.DataFrame(rows)
            safe_tables.append(safe)

            col_reports = []
            for idx in pam_cols:
                col = header[idx]
                vals = safe[col].tolist() if col in safe.columns else []
                col_reports.append(
                    {
                        "column": col,
                        "display_name": display[idx] if idx < len(display) else "",
                        "description": desc[idx] if idx < len(desc) else "",
                        "nonmissing": sum(bool(norm_value(v)) for v in vals),
                        "levels": count_levels(vals, canonicalize=True),
                    }
                )

            receptor_reports = []
            for idx in receptor_cols:
                col = header[idx]
                vals = safe[col].tolist() if col in safe.columns else []
                receptor_reports.append(
                    {
                        "column": col,
                        "display_name": display[idx] if idx < len(display) else "",
                        "description": desc[idx] if idx < len(desc) else "",
                        "nonmissing": sum(bool(norm_value(v)) for v in vals),
                        "levels": count_levels(vals, canonicalize=False),
                    }
                )

            match_info = {}
            if expression_samples:
                expr_set = set(expression_samples)
                patient_set = {"-".join(x.split("-")[:3]) for x in expression_samples if x.startswith("TCGA-")}
                for id_col in [c for c in safe.columns if c.upper() in id_names]:
                    vals = set(norm_value(x) for x in safe[id_col].tolist() if norm_value(x))
                    exact = len(vals & expr_set)
                    patient = len(vals & patient_set)
                    match_info[id_col] = {
                        "exact_expression_id_matches": exact,
                        "patient_id_matches": patient,
                    }

            results.append(
                {
                    "source": source_label,
                    "member": member.name,
                    "rows": len(data_lines),
                    "pam50_candidate_columns": col_reports,
                    "receptor_candidate_columns": receptor_reports,
                    "expression_id_match": match_info,
                }
            )

    return results, safe_tables


def parse_scanb_series(series_paths: list[Path]) -> tuple[pd.DataFrame, list[dict]]:
    rows = []
    reports = []

    for path in series_paths:
        require(path)

        sample_titles = None
        sample_geo = None
        sample_platform = None
        subtype_rows: list[tuple[str, list[str]]] = []

        with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
            for line in f:
                if line.startswith("!Sample_title"):
                    sample_titles = next(csv.reader([line.split("\t", 1)[1]], delimiter="\t"))
                elif line.startswith("!Sample_geo_accession"):
                    sample_geo = next(csv.reader([line.split("\t", 1)[1]], delimiter="\t"))
                elif line.startswith("!Sample_platform_id"):
                    sample_platform = next(csv.reader([line.split("\t", 1)[1]], delimiter="\t"))
                elif line.startswith("!Sample_characteristics_ch1"):
                    payload = next(csv.reader([line.split("\t", 1)[1]], delimiter="\t"))
                    # Inspect only the characteristic key before ":".
                    keys = []
                    for v in payload:
                        vv = str(v).strip('"')
                        key = vv.split(":", 1)[0].strip() if ":" in vv else ""
                        keys.append(key)
                    nonblank_keys = [k for k in keys if k]
                    if nonblank_keys:
                        dominant = Counter(nonblank_keys).most_common(1)[0]
                        if dominant[1] >= max(1, int(0.5 * len(nonblank_keys))):
                            key = dominant[0]
                            if is_pam50_label(key):
                                subtype_rows.append((key, payload))
                elif line.startswith("!series_matrix_table_begin"):
                    break

        if sample_titles is None or sample_geo is None or sample_platform is None:
            raise RuntimeError(f"Missing SCAN-B sample headers in {path}")

        n = len(sample_geo)
        if len(sample_titles) != n or len(sample_platform) != n:
            raise RuntimeError(f"SCAN-B header length mismatch in {path}")

        if not subtype_rows:
            reports.append(
                {
                    "source": "SCANB_GSE96058",
                    "series_matrix": path.name,
                    "samples": n,
                    "pam50_rows_found": 0,
                }
            )
            continue

        for key, payload in subtype_rows:
            if len(payload) != n:
                raise RuntimeError(f"Subtype row length mismatch in {path}: {key}")

            vals = []
            for v in payload:
                vv = str(v).strip('"')
                val = vv.split(":", 1)[1].strip() if ":" in vv else vv.strip()
                vals.append(canonicalize_pam50(val))

            for i in range(n):
                rows.append(
                    {
                        "sample_title": str(sample_titles[i]).strip('"'),
                        "geo_accession": str(sample_geo[i]).strip('"'),
                        "platform_id": str(sample_platform[i]).strip('"'),
                        "pam50_source_key": key,
                        "pam50_subtype": vals[i],
                        "source_series_matrix": path.name,
                    }
                )

            reports.append(
                {
                    "source": "SCANB_GSE96058",
                    "series_matrix": path.name,
                    "samples": n,
                    "pam50_rows_found": 1,
                    "key": key,
                    "nonmissing": sum(bool(norm_value(v)) for v in vals),
                    "levels": count_levels(vals, canonicalize=True),
                }
            )

    return pd.DataFrame(rows), reports


def audit_thesis_file(path: Path, tcga_samples: list[str]) -> dict | None:
    sep = "\t" if path.suffix.lower() in {".tsv", ".txt"} else ","
    try:
        header = pd.read_csv(path, sep=sep, nrows=0, low_memory=False)
    except Exception as exc:
        return {
            "path": str(path),
            "read_error": repr(exc),
        }

    cols = list(header.columns)
    pam_cols = [c for c in cols if is_pam50_label(c)]
    receptor_cols = [c for c in cols if (not is_pam50_label(c)) and is_receptor_label(c)]

    # Potential ID columns are safe to inspect.
    id_regex = re.compile(r"(sample|patient|case|submitter|barcode|bcr|id$)", re.I)
    id_cols = [c for c in cols if id_regex.search(c)]
    id_cols = id_cols[:10]

    if not pam_cols and not receptor_cols:
        return {
            "path": str(path),
            "rows": None,
            "pam50_candidate_columns": [],
            "receptor_candidate_columns": [],
            "header_columns": len(cols),
        }

    usecols = list(dict.fromkeys(id_cols + pam_cols + receptor_cols))
    try:
        df = pd.read_csv(path, sep=sep, usecols=usecols, dtype=str, low_memory=False)
    except Exception as exc:
        return {
            "path": str(path),
            "read_error": repr(exc),
            "pam50_candidate_columns": pam_cols,
            "receptor_candidate_columns": receptor_cols,
        }

    tcga_expr_set = set(tcga_samples)
    tcga_patient_set = {"-".join(x.split("-")[:3]) for x in tcga_samples if x.startswith("TCGA-")}

    pam_reports = []
    for c in pam_cols:
        pam_reports.append(
            {
                "column": c,
                "nonmissing": int(df[c].map(lambda x: bool(norm_value(x))).sum()),
                "levels": count_levels(df[c].tolist(), canonicalize=True),
            }
        )

    receptor_reports = []
    for c in receptor_cols:
        receptor_reports.append(
            {
                "column": c,
                "nonmissing": int(df[c].map(lambda x: bool(norm_value(x))).sum()),
                "levels": count_levels(df[c].tolist(), canonicalize=False),
            }
        )

    id_matches = {}
    for c in id_cols:
        vals = set(norm_value(x) for x in df[c].tolist() if norm_value(x))
        id_matches[c] = {
            "exact_expression_id_matches": len(vals & tcga_expr_set),
            "patient_id_matches": len(vals & tcga_patient_set),
        }

    return {
        "path": str(path),
        "rows": int(len(df)),
        "pam50_candidate_columns": pam_reports,
        "receptor_candidate_columns": receptor_reports,
        "expression_id_match": id_matches,
    }


def best_local_pam50_source(report: dict) -> dict:
    candidates = []

    for section in ["TCGA_cBioPortal", "METABRIC_cBioPortal", "Thesis_candidates"]:
        items = report.get(section, [])
        for item in items:
            for c in item.get("pam50_candidate_columns", []):
                candidates.append(
                    {
                        "dataset_group": section,
                        "source": item.get("source", item.get("path", "")),
                        "member_or_path": item.get("member", item.get("path", "")),
                        "column": c.get("column", c.get("key", "")),
                        "nonmissing": int(c.get("nonmissing", 0)),
                        "levels": c.get("levels", {}),
                    }
                )

    for item in report.get("SCANB", []):
        if item.get("pam50_rows_found", 0):
            candidates.append(
                {
                    "dataset_group": "SCANB",
                    "source": item.get("source", ""),
                    "member_or_path": item.get("series_matrix", ""),
                    "column": item.get("key", ""),
                    "nonmissing": int(item.get("nonmissing", 0)),
                    "levels": item.get("levels", {}),
                }
            )

    candidates = sorted(candidates, key=lambda x: x["nonmissing"], reverse=True)
    return {"candidates_ranked_by_nonmissing": candidates}


def main() -> None:
    print("=" * 128)
    print("Paper 4 / TCBB - SAFE breast-cancer subtype metadata audit before subtype-sensitivity contract")
    print("=" * 128)
    print(f"Script version: {SCRIPT_VERSION}")
    print()
    print("Scientific firewall:")
    print("  Survival/response/recurrence outcomes loaded:         NO")
    print("  Drug/treatment variables loaded:                      NO")
    print("  Target preservation statistics calculated:           NO")
    print("  PAM50/intrinsic-subtype metadata permitted:           YES")
    print("  ER/PR/HER2/TNBC inspected only as BACKUP grouping:    YES")
    print("  Receptor groups will NOT be substituted for PAM50:    YES")
    print("=" * 128)

    for p in [TCGA_ARCHIVE, METABRIC_ARCHIVE, TCGA_EXPR, METABRIC_EXPR, *SCANB_SERIES, SCANB_PRIMARY]:
        require(p)

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    tcga_samples = expression_sample_headers(TCGA_EXPR, sep="\t", id_cols=2)
    metabric_samples = expression_sample_headers(METABRIC_EXPR, sep="\t", id_cols=2)

    print("\n[1/4] Auditing TCGA-BRCA and METABRIC cBioPortal clinical schemas ...")
    tcga_reports = []
    tcga_safe_tables = []
    for base in ["data_clinical_sample.txt", "data_clinical_patient.txt"]:
        reps, tabs = parse_cbio_clinical_member(
            TCGA_ARCHIVE,
            base,
            "TCGA_BRCA_PanCanAtlas2018",
            expression_samples=tcga_samples,
        )
        tcga_reports.extend(reps)
        tcga_safe_tables.extend(tabs)

    metabric_reports = []
    metabric_safe_tables = []
    for base in ["data_clinical_sample.txt", "data_clinical_patient.txt"]:
        reps, tabs = parse_cbio_clinical_member(
            METABRIC_ARCHIVE,
            base,
            "METABRIC",
            expression_samples=metabric_samples,
        )
        metabric_reports.extend(reps)
        metabric_safe_tables.extend(tabs)

    for label, reps in [("TCGA", tcga_reports), ("METABRIC", metabric_reports)]:
        print(f"  {label}:")
        for r in reps:
            pam_cols = r.get("pam50_candidate_columns", [])
            rec_cols = r.get("receptor_candidate_columns", [])
            print(
                f"    {Path(r.get('member', 'unknown')).name}: "
                f"PAM50 candidates={len(pam_cols)}, receptor candidates={len(rec_cols)}"
            )
            for c in pam_cols:
                print(
                    f"      PAM50 {c['column']}: nonmissing={c['nonmissing']}, "
                    f"levels={c['levels']}"
                )

    print("\n[2/4] Auditing SCAN-B series matrices for PAM50 rows only ...")
    scanb_manifest, scanb_reports = parse_scanb_series(SCANB_SERIES)
    for r in scanb_reports:
        if r.get("pam50_rows_found", 0):
            print(
                f"  {r['series_matrix']}: key='{r['key']}', "
                f"nonmissing={r['nonmissing']}, levels={r['levels']}"
            )
        else:
            print(f"  {r['series_matrix']}: no PAM50 row found")

    if not scanb_manifest.empty:
        scanb_primary = pd.read_csv(SCANB_PRIMARY, sep="\t", dtype=str)
        primary_titles = set(scanb_primary["primary_title"].astype(str))
        scanb_manifest["is_frozen_primary_profile"] = (
            scanb_manifest["sample_title"].isin(primary_titles)
        ).astype(int)

        primary_scanb = scanb_manifest[
            scanb_manifest["is_frozen_primary_profile"] == 1
        ].copy()

        print(
            f"  SCAN-B primary profiles with PAM50 record: "
            f"{int(primary_scanb['pam50_subtype'].map(lambda x: bool(norm_value(x))).sum())}"
        )
        print(
            f"  SCAN-B frozen-primary PAM50 levels: "
            f"{count_levels(primary_scanb['pam50_subtype'].tolist(), canonicalize=True)}"
        )

        scanb_out = OUT_DIR / "scanb_pam50_safe_manifest_audit_v1.tsv"
        primary_scanb[
            [
                "sample_title",
                "geo_accession",
                "platform_id",
                "pam50_subtype",
                "source_series_matrix",
            ]
        ].to_csv(scanb_out, sep="\t", index=False)
    else:
        scanb_out = None

    print("\n[3/4] Auditing likely thesis clinical files ...")
    thesis_files = find_thesis_candidates()
    thesis_reports = []
    if not thesis_files:
        print("  No local thesis clinical candidate file found under expected paths/Desktop.")
    else:
        for p in thesis_files:
            rep = audit_thesis_file(p, tcga_samples)
            if rep is not None:
                thesis_reports.append(rep)
                print(f"  {p}")
                for c in rep.get("pam50_candidate_columns", []):
                    if isinstance(c, dict):
                        print(
                            f"    PAM50 candidate {c['column']}: "
                            f"nonmissing={c['nonmissing']}, levels={c['levels']}"
                        )
                if not rep.get("pam50_candidate_columns"):
                    print("    no PAM50/intrinsic-subtype column detected")

    print("\n[4/4] Writing SAFE subtype audit outputs ...")

    report = {
        "script_version": SCRIPT_VERSION,
        "status": "SUBTYPE_METADATA_AUDIT_ONLY_NOT_YET_FROZEN_FOR_ANALYSIS",
        "scientific_firewall": {
            "survival_response_recurrence_loaded": False,
            "drug_treatment_variables_loaded": False,
            "target_preservation_statistics_calculated": False,
            "pam50_intrinsic_subtype_permitted": True,
            "receptor_status_inspected_as_backup_only": True,
            "receptor_status_not_equated_to_pam50": True,
        },
        "TCGA_cBioPortal": tcga_reports,
        "METABRIC_cBioPortal": metabric_reports,
        "SCANB": scanb_reports,
        "Thesis_candidates": thesis_reports,
    }
    report["ranked_local_pam50_candidates"] = best_local_pam50_source(report)

    json_out = OUT_DIR / "breast_subtype_metadata_audit_v1.json"
    json_out.write_text(json.dumps(report, indent=2), encoding="utf-8")

    # Persist only subtype/receptor-safe tables from cBioPortal.
    for i, df in enumerate(tcga_safe_tables, start=1):
        df.to_csv(
            OUT_DIR / f"tcga_safe_subtype_candidate_table_{i}_v1.tsv",
            sep="\t",
            index=False,
        )
    for i, df in enumerate(metabric_safe_tables, start=1):
        df.to_csv(
            OUT_DIR / f"metabric_safe_subtype_candidate_table_{i}_v1.tsv",
            sep="\t",
            index=False,
        )

    ranked = report["ranked_local_pam50_candidates"]["candidates_ranked_by_nonmissing"]

    print("\n" + "=" * 128)
    print("04b BREAST SUBTYPE METADATA AUDIT: PASS")
    print("=" * 128)
    if ranked:
        print("Top local PAM50/intrinsic-subtype candidates:")
        for x in ranked[:10]:
            print(
                f"  {x['dataset_group']}: {x['column']} "
                f"nonmissing={x['nonmissing']} levels={x['levels']}"
            )
    else:
        print("No explicit local PAM50/intrinsic-subtype field was found.")

    print()
    print("No survival/response/treatment field was persisted.")
    print("No preservation statistic was calculated.")
    print("Next step: choose canonical subtype source(s), then freeze 04c subtype sensitivity contract.")
    print()
    print(f"Output: {json_out}")
    if scanb_out is not None:
        print(f"SCAN-B safe PAM50 manifest audit: {scanb_out}")
    print("=" * 128)


if __name__ == "__main__":
    main()
