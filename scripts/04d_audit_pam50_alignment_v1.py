from __future__ import annotations

import json
import re
from collections import Counter
from pathlib import Path

import pandas as pd


SCRIPT_VERSION = "04d-audit-pam50-alignment-v1-no-cli"

DATA_ROOT = Path(r"D:\paper4_tcbb_data")

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

TCGA_XENA_SAFE = (
    DATA_ROOT
    / "paper4_tcbb_tcga_pam50_xena_v3"
    / "tcga_pam50_safe_manifest_v3.tsv"
)

SUBTYPE_AUDIT_DIR = DATA_ROOT / "paper4_tcbb_breast_subtype_audit_v1"
SCANB_SAFE = SUBTYPE_AUDIT_DIR / "scanb_pam50_safe_manifest_audit_v1.tsv"
METABRIC_SAFE_GLOB = "metabric_safe_subtype_candidate_table_*_v1.tsv"

OUT_DIR = DATA_ROOT / "paper4_tcbb_pam50_alignment_audit_v1"

CANONICAL = ["LumA", "LumB", "Basal", "Her2", "Normal"]
CANONICAL_SET = set(CANONICAL)


def require(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Required file does not exist: {path}")


def clean(x: object) -> str:
    if x is None:
        return ""
    s = str(x).strip().strip('"')
    if s.lower() in {
        "", "na", "nan", "n/a", "unknown", "not available", "[not available]"
    }:
        return ""
    return s


def canon_subtype(x: object) -> str:
    s = clean(x)
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
        "nc": "NC",
    }
    return mapping.get(low, s)


TCGA_SAMPLE_RE = re.compile(r"^(TCGA-[A-Za-z0-9]{2}-[A-Za-z0-9]{4}-[A-Za-z0-9]{2})")


def canonical_tcga_sample_id(x: object) -> str:
    s = clean(x).upper()
    m = TCGA_SAMPLE_RE.match(s)
    return m.group(1) if m else s


def expression_samples(path: Path) -> list[str]:
    hdr = pd.read_csv(path, sep="\t", nrows=0, low_memory=False)
    cols = list(hdr.columns)
    if len(cols) < 3:
        raise RuntimeError(f"Expression file has too few columns: {path}")
    # Both frozen cBioPortal expression matrices use two gene-ID columns.
    return [clean(x) for x in cols[2:]]


def counts(values) -> dict[str, int]:
    c = Counter(canon_subtype(v) for v in values if canon_subtype(v))
    order = CANONICAL + sorted(k for k in c if k not in CANONICAL_SET)
    return {k: int(c[k]) for k in order if c[k] > 0}


def assert_no_conflicting_labels(df: pd.DataFrame, id_col: str, subtype_col: str, label: str) -> None:
    tmp = df[[id_col, subtype_col]].copy()
    tmp = tmp[(tmp[id_col] != "") & (tmp[subtype_col] != "")]
    conflicts = (
        tmp.groupby(id_col)[subtype_col]
        .nunique()
        .loc[lambda x: x > 1]
    )
    if len(conflicts):
        examples = conflicts.index.tolist()[:20]
        raise RuntimeError(f"{label}: conflicting subtype labels for IDs: {examples}")


def load_tcga_alignment(tcga_expr_ids: list[str]) -> tuple[pd.DataFrame, dict]:
    df = pd.read_csv(TCGA_XENA_SAFE, sep="\t", dtype=str).fillna("")
    required = {"xena_sample_id", "pam50_subtype"}
    if not required.issubset(df.columns):
        raise RuntimeError(f"TCGA safe manifest missing columns: {required - set(df.columns)}")

    df["tcga_sample_id"] = df["xena_sample_id"].map(canonical_tcga_sample_id)
    df["pam50_subtype"] = df["pam50_subtype"].map(canon_subtype)
    assert_no_conflicting_labels(df, "tcga_sample_id", "pam50_subtype", "TCGA Xena")

    subtype_by_sample = (
        df[df["pam50_subtype"] != ""]
        .drop_duplicates(["tcga_sample_id", "pam50_subtype"])
        .set_index("tcga_sample_id")["pam50_subtype"]
        .to_dict()
    )

    expr = pd.DataFrame({"expression_sample_id": tcga_expr_ids})
    expr["tcga_sample_id"] = expr["expression_sample_id"].map(canonical_tcga_sample_id)

    if expr["tcga_sample_id"].duplicated().any():
        dups = expr.loc[expr["tcga_sample_id"].duplicated(False), "tcga_sample_id"].unique().tolist()
        raise RuntimeError(f"TCGA frozen expression has duplicate canonical sample IDs: {dups[:20]}")

    expr["pam50_subtype"] = expr["tcga_sample_id"].map(subtype_by_sample).fillna("")
    expr["has_pam50"] = (expr["pam50_subtype"] != "").astype(int)

    expr_types = Counter()
    for sid in expr["tcga_sample_id"]:
        parts = sid.split("-")
        if len(parts) >= 4:
            expr_types[parts[3][:2]] += 1

    xena_ids = set(df.loc[df["pam50_subtype"] != "", "tcga_sample_id"])
    expr_ids = set(expr["tcga_sample_id"])

    report = {
        "frozen_expression_profiles": int(len(expr)),
        "frozen_expression_sample_type_codes": dict(expr_types),
        "pam50_matched_profiles": int(expr["has_pam50"].sum()),
        "pam50_missing_profiles": int((expr["has_pam50"] == 0).sum()),
        "pam50_coverage_fraction": float(expr["has_pam50"].mean()),
        "matched_counts": counts(expr["pam50_subtype"]),
        "xena_nonmissing_pam50_total": int((df["pam50_subtype"] != "").sum()),
        "xena_labeled_samples_not_in_frozen_expression": int(len(xena_ids - expr_ids)),
        "xena_outside_frozen_counts": counts(
            df.loc[
                (df["pam50_subtype"] != "")
                & (~df["tcga_sample_id"].isin(expr_ids)),
                "pam50_subtype",
            ]
        ),
    }
    return expr, report


def load_scanb_alignment() -> tuple[pd.DataFrame, dict]:
    df = pd.read_csv(SCANB_SAFE, sep="\t", dtype=str).fillna("")
    required = {"sample_title", "geo_accession", "platform_id", "pam50_subtype"}
    if not required.issubset(df.columns):
        raise RuntimeError(f"SCAN-B safe manifest missing columns: {required - set(df.columns)}")

    df["pam50_subtype"] = df["pam50_subtype"].map(canon_subtype)

    if df["sample_title"].duplicated().any():
        dups = df.loc[df["sample_title"].duplicated(False), "sample_title"].unique().tolist()
        raise RuntimeError(f"SCAN-B primary subtype manifest has duplicate sample titles: {dups[:20]}")

    report = {
        "frozen_primary_profiles": int(len(df)),
        "pam50_matched_profiles": int((df["pam50_subtype"] != "").sum()),
        "pam50_missing_profiles": int((df["pam50_subtype"] == "").sum()),
        "pam50_coverage_fraction": float((df["pam50_subtype"] != "").mean()),
        "pam50_counts": counts(df["pam50_subtype"]),
        "platform_counts": {
            str(k): int(v)
            for k, v in Counter(df["platform_id"]).items()
        },
    }
    return df, report


def find_metabric_table() -> tuple[Path, pd.DataFrame, str]:
    paths = sorted(SUBTYPE_AUDIT_DIR.glob(METABRIC_SAFE_GLOB))
    if not paths:
        raise FileNotFoundError(
            f"No METABRIC safe subtype candidate tables matching {METABRIC_SAFE_GLOB}"
        )

    matches = []
    for p in paths:
        df = pd.read_csv(p, sep="\t", dtype=str).fillna("")
        subtype_cols = [
            c for c in df.columns
            if c.upper() == "CLAUDIN_SUBTYPE"
            or "PAM50" in c.upper()
            or ("SUBTYPE" in c.upper() and "CLAUDIN" in c.upper())
        ]
        if subtype_cols:
            matches.append((p, df, subtype_cols[0]))

    if len(matches) != 1:
        detail = [(str(p), c) for p, _, c in matches]
        raise RuntimeError(
            "Expected exactly one safe METABRIC table containing CLAUDIN/PAM50 subtype; "
            f"found {len(matches)}: {detail}"
        )
    return matches[0]


def load_metabric_alignment(met_expr_ids: list[str]) -> tuple[pd.DataFrame, dict]:
    p, df, subtype_col = find_metabric_table()

    id_candidates = [
        c for c in df.columns
        if c.upper() in {"PATIENT_ID", "SAMPLE_ID", "CASE_ID"}
    ]
    if not id_candidates:
        raise RuntimeError(f"METABRIC safe table has no ID column: {p}")

    # Prefer exact ID column that overlaps the expression columns the most.
    expr_set = set(met_expr_ids)
    scored = []
    for c in id_candidates:
        vals = set(clean(x) for x in df[c] if clean(x))
        scored.append((len(vals & expr_set), c))
    scored.sort(reverse=True)
    best_overlap, id_col = scored[0]

    df = df.copy()
    df[id_col] = df[id_col].map(clean)
    df["subtype_raw"] = df[subtype_col].map(clean)
    df["subtype_canonical"] = df[subtype_col].map(canon_subtype)

    assert_no_conflicting_labels(df, id_col, "subtype_canonical", "METABRIC")

    subtype_by_id = (
        df[df["subtype_canonical"] != ""]
        .drop_duplicates([id_col, "subtype_canonical"])
        .set_index(id_col)["subtype_canonical"]
        .to_dict()
    )

    expr = pd.DataFrame({"expression_sample_id": met_expr_ids})
    expr["subtype_raw"] = expr["expression_sample_id"].map(
        df.drop_duplicates(id_col).set_index(id_col)["subtype_raw"].to_dict()
    ).fillna("")
    expr["subtype_canonical"] = expr["expression_sample_id"].map(subtype_by_id).fillna("")
    expr["pam50_compatible"] = expr["subtype_canonical"].isin(CANONICAL_SET).astype(int)

    report = {
        "safe_source_table": str(p),
        "safe_source_id_column": id_col,
        "safe_source_subtype_column": subtype_col,
        "frozen_expression_profiles": int(len(expr)),
        "exact_id_overlap_before_subtype_filter": int(best_overlap),
        "subtype_matched_profiles": int((expr["subtype_canonical"] != "").sum()),
        "all_subtype_counts": counts(expr["subtype_canonical"]),
        "pam50_compatible_profiles": int(expr["pam50_compatible"].sum()),
        "pam50_compatible_counts": counts(
            expr.loc[expr["pam50_compatible"] == 1, "subtype_canonical"]
        ),
        "non_pam50_extension_counts": counts(
            expr.loc[expr["pam50_compatible"] == 0, "subtype_canonical"]
        ),
    }
    return expr, report


def main() -> None:
    print("=" * 132)
    print("Paper 4 / TCBB - audit exact PAM50/subtype alignment to the frozen breast-cancer expression cohorts")
    print("=" * 132)
    print(f"Script version: {SCRIPT_VERSION}")
    print()
    print("Scientific firewall:")
    print("  Survival/response/recurrence outcomes loaded:         NO")
    print("  Drug/treatment variables loaded:                      NO")
    print("  Target gene-gene correlations calculated:             NO")
    print("  Target PCA/loadings calculated:                       NO")
    print("  Preservation statistics calculated:                   NO")
    print("  Operation: subtype ID alignment and sample counts only")
    print("=" * 132)

    for p in [TCGA_EXPR, METABRIC_EXPR, TCGA_XENA_SAFE, SCANB_SAFE]:
        require(p)

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print("\n[1/3] TCGA-BRCA Xena PAM50 -> frozen 1,082-source expression profiles ...")
    tcga_ids = expression_samples(TCGA_EXPR)
    tcga_manifest, tcga_report = load_tcga_alignment(tcga_ids)
    print(f"  frozen expression profiles: {tcga_report['frozen_expression_profiles']}")
    print(f"  sample type codes:          {tcga_report['frozen_expression_sample_type_codes']}")
    print(
        f"  PAM50 matched:              {tcga_report['pam50_matched_profiles']}/"
        f"{tcga_report['frozen_expression_profiles']} "
        f"({tcga_report['pam50_coverage_fraction']:.3f})"
    )
    print(f"  matched counts:             {tcga_report['matched_counts']}")
    print(
        f"  labeled Xena samples outside frozen expression: "
        f"{tcga_report['xena_labeled_samples_not_in_frozen_expression']}"
    )
    print(f"  outside-frozen counts:      {tcga_report['xena_outside_frozen_counts']}")

    print("\n[2/3] SCAN-B PAM50 -> frozen 3,273 primary profiles ...")
    scanb_manifest, scanb_report = load_scanb_alignment()
    print(
        f"  PAM50 matched:              {scanb_report['pam50_matched_profiles']}/"
        f"{scanb_report['frozen_primary_profiles']} "
        f"({scanb_report['pam50_coverage_fraction']:.3f})"
    )
    print(f"  counts:                     {scanb_report['pam50_counts']}")
    print(f"  platform counts:            {scanb_report['platform_counts']}")

    print("\n[3/3] METABRIC extended subtype -> frozen 1,980 expression profiles ...")
    met_ids = expression_samples(METABRIC_EXPR)
    met_manifest, met_report = load_metabric_alignment(met_ids)
    print(f"  safe subtype source:        {met_report['safe_source_table']}")
    print(f"  source subtype field:       {met_report['safe_source_subtype_column']}")
    print(
        f"  subtype matched:            {met_report['subtype_matched_profiles']}/"
        f"{met_report['frozen_expression_profiles']}"
    )
    print(f"  all subtype counts:         {met_report['all_subtype_counts']}")
    print(
        f"  PAM50-compatible:           {met_report['pam50_compatible_profiles']}/"
        f"{met_report['frozen_expression_profiles']}"
    )
    print(f"  compatible counts:          {met_report['pam50_compatible_counts']}")
    print(f"  extension counts:           {met_report['non_pam50_extension_counts']}")

    # Write SAFE alignment manifests only.
    tcga_out = OUT_DIR / "tcga_frozen_expression_pam50_alignment_v1.tsv"
    scanb_out = OUT_DIR / "scanb_frozen_primary_pam50_alignment_v1.tsv"
    met_out = OUT_DIR / "metabric_frozen_expression_subtype_alignment_v1.tsv"

    tcga_manifest.to_csv(tcga_out, sep="\t", index=False)
    scanb_manifest.to_csv(scanb_out, sep="\t", index=False)
    met_manifest.to_csv(met_out, sep="\t", index=False)

    summary = {
        "script_version": SCRIPT_VERSION,
        "status": "ALIGNMENT_AUDIT_ONLY_NOT_YET_SUBTYPE_ANALYSIS_CONTRACT",
        "scientific_firewall": {
            "survival_response_recurrence_loaded": False,
            "drug_treatment_loaded": False,
            "target_correlations_calculated": False,
            "target_pca_calculated": False,
            "preservation_statistics_calculated": False,
        },
        "TCGA_BRCA": tcga_report,
        "SCANB_GSE96058": scanb_report,
        "METABRIC": met_report,
        "important_interpretation_guard": {
            "tcga_raw_xena_counts_not_used_directly": True,
            "reason": (
                "The raw Xena phenotype matrix contains samples outside the frozen "
                "1,082 primary-tumour expression cohort. Only the exact aligned frozen "
                "source profiles may enter subtype sensitivity analyses."
            ),
            "metabric_claudin_subtype_not_equated_to_pure_pam50": True,
            "metabric_pam50_compatible_labels": CANONICAL,
            "metabric_extension_labels_kept_separate": ["Claudin-low", "NC"],
        },
    }

    json_out = OUT_DIR / "pam50_alignment_audit_v1.json"
    json_out.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print("\n" + "=" * 132)
    print("04d PAM50 / SUBTYPE ALIGNMENT AUDIT: PASS")
    print("=" * 132)
    print("No preservation statistic was calculated.")
    print("No subtype minimum-n rule has been selected yet.")
    print("Next: freeze the subtype-specific and subtype-residualized sensitivity contract using THESE aligned counts only.")
    print()
    print("Outputs:")
    print(f"  {json_out}")
    print(f"  {tcga_out}")
    print(f"  {scanb_out}")
    print(f"  {met_out}")
    print("=" * 132)


if __name__ == "__main__":
    main()
