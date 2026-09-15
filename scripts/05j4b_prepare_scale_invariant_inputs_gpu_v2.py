from __future__ import annotations

import gc
import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd
import torch


SCRIPT_VERSION = "05j4b-prepare-scale-invariant-inputs-gpu-v2-no-cli"

PROJECT_ROOT = Path(r"C:\Users\olegk\Desktop\molecular-transport-audit")
DATA_ROOT = Path(r"D:\paper4_tcbb_data")

SOURCE_05D = PROJECT_ROOT / "scripts" / "05d_run_controlled_mapping_degradation_gpu_v1.py"
SOURCE_05J = PROJECT_ROOT / "scripts" / "05j_run_structural_corruption_headtohead_gpu_v1.py"

CONTRACT = (
    DATA_ROOT
    / "paper4_tcbb_scale_invariant_corruption_audit_contract_v3"
    / "scale_invariant_corruption_audit_contract_v3.json"
)
HEADTOHEAD_REPLICATES = (
    DATA_ROOT
    / "paper4_tcbb_structural_corruption_headtohead_v1"
    / "structural_headtohead_replicates_v1.tsv"
)
CLASSIFICATION = (
    DATA_ROOT
    / "paper4_tcbb_final_mapping_specificity_null_v2"
    / "primary_pooled_final_classification_v2.tsv"
)
SPECIFICITY_ROOT = DATA_ROOT / "paper4_tcbb_final_mapping_specificity_null_v2"
DIRECT_ROOT = DATA_ROOT / "paper4_tcbb_primary_pooled_direct_preservation_v2"

SCANB_FEATURES = (
    DATA_ROOT / "paper4_tcbb_target_marginal_metrics_v2"
    / "scanb_target_marginal_matching_features_v2.tsv"
)
METABRIC_FEATURES = (
    DATA_ROOT / "paper4_tcbb_target_marginal_metrics_v2"
    / "metabric_target_marginal_matching_features_v2.tsv"
)
SCANB_EVAL = (
    DATA_ROOT / "paper4_tcbb_exact_variation_evaluability_correction_v1"
    / "scanb_target_exact_variation_evaluability_v1.tsv"
)
METABRIC_EVAL = (
    DATA_ROOT / "paper4_tcbb_exact_variation_evaluability_correction_v1"
    / "metabric_target_exact_variation_evaluability_v1.tsv"
)
WEIGHTS = (
    DATA_ROOT / "paper4_tcbb_frozen_source_programs_v1"
    / "tcga_frozen_source_program_weights_v1.tsv"
)

OUT_DIR = DATA_ROOT / "paper4_tcbb_scale_invariant_corruption_inputs_v3"
COMBO_DIR = OUT_DIR / "per_combo"

FULL_VALID = 1000
FIRST_REPLAY = 100
NEW_PARTIAL_FRACTIONS = [0.05, 0.15]
NEW_PARTIAL_REPLICATES = 100
MAX_ATTEMPTS = 10000
MAX_PARTIAL_ATTEMPTS = 1000

REPLAY_TOL = 5e-10
COMPARATOR_TOL = 5e-8


def require(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Required file does not exist: {path}")


def boolish(x: object) -> bool:
    if isinstance(x, bool):
        return x
    return str(x).strip().lower() in {"1", "true", "yes", "y"}


def canon_symbol(x: object) -> str:
    if x is None:
        return ""
    return " ".join(str(x).strip().split()).upper()


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot import {path}")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def pearson_np_against_fixed(fixed: np.ndarray, values: np.ndarray) -> float:
    x = np.asarray(fixed, dtype=np.float64)
    y = np.asarray(values, dtype=np.float64)
    xc = x - x.mean()
    yc = y - y.mean()
    den = np.linalg.norm(xc) * np.linalg.norm(yc)
    if den <= 0:
        return np.nan
    return float(np.dot(xc, yc) / den)


def build_mapping(
    exact05d,
    target_name: str,
    program_id: str,
    fraction: float,
    attempt_id: int,
    slot_genes: list[str],
    full_program_genes: set[str],
    features: pd.DataFrame,
    target_gene_index: dict[str, int],
    correct_target_idx: np.ndarray,
) -> tuple[np.ndarray, dict]:
    m = len(slot_genes)
    if fraction == 0.0:
        return correct_target_idx.copy(), {
            "panel_valid": True,
            "fraction_distance_0": np.nan,
            "fraction_distance_le1": np.nan,
            "maximum_distance": np.nan,
            "mean_distance": np.nan,
        }

    k = max(1, min(m, exact05d.round_half_up(fraction * m)))
    rng = exact05d.rng_for_attempt(target_name, program_id, fraction, attempt_id)
    corrupted = np.sort(rng.choice(m, size=k, replace=False))
    problem = exact05d.prepare_problem(
        slot_genes=slot_genes,
        corrupted_slot_idx=corrupted,
        full_program_genes=full_program_genes,
        features=features,
        target_gene_index=target_gene_index,
    )
    mapping_rows, distances, success, reason = exact05d.generate_mapping(problem, rng)
    if not success:
        raise RuntimeError(f"Mapping generation failed: {reason}")
    q = exact05d.quality(distances)
    if not q["panel_valid"]:
        raise RuntimeError("Mapping is not panel-valid.")

    candidate_rows = problem["candidate"].iloc[mapping_rows]
    replacement_target_idx = candidate_rows["target_corr_index"].to_numpy(dtype=np.int64)

    mapped = correct_target_idx.copy()
    mapped[corrupted] = replacement_target_idx

    if len(np.unique(mapped)) != len(mapped):
        raise RuntimeError("Duplicate mapped target genes.")

    return mapped, q


def score_mapping(
    mapped: np.ndarray,
    target_corr: torch.Tensor,
    ii: np.ndarray,
    jj: np.ndarray,
    ii_t: torch.Tensor,
    jj_t: torch.Tensor,
    source_subset: np.ndarray,
    source_rank: np.ndarray,
    source_rank_norm: float,
    tri_i_t: torch.Tensor,
    tri_j_t: torch.Tensor,
    source_all_centered: torch.Tensor,
    source_all_norm: torch.Tensor,
    source_kim_centered: torch.Tensor,
    source_kim_norm: torch.Tensor,
    exact05d,
    exact05j,
    device: torch.device,
) -> dict:
    mapped_t = torch.as_tensor(mapped, dtype=torch.long, device=device)
    sub = target_corr.index_select(0, mapped_t).index_select(1, mapped_t)

    subset_edges = sub[ii_t, jj_t].detach().cpu().numpy()

    mta = exact05d.spearman_against_fixed(
        source_rank,
        source_rank_norm,
        subset_edges,
    )
    pearson_subset = pearson_np_against_fixed(source_subset, subset_edges)

    target_all_edges = sub[tri_i_t, tri_j_t]
    corcor = exact05j.pearson_against_fixed(
        source_all_centered,
        source_all_norm,
        target_all_edges,
    )

    target_adj12 = torch.pow((1.0 + sub) / 2.0, 12.0)
    target_kim = torch.sum(target_adj12, dim=0)
    corkim = exact05j.pearson_against_fixed(
        source_kim_centered,
        source_kim_norm,
        target_kim,
    )

    del mapped_t, sub, target_all_edges, target_adj12, target_kim

    return {
        "mta_rho_edge": float(mta),
        "pearson_subset_matched": float(pearson_subset),
        "netrep_cor_cor": float(corcor),
        "wgcna_cor_kIM": float(corkim),
    }


def main() -> None:
    print("=" * 164)
    print("Paper 4 / TCBB - prepare scale-invariant audit v3 inputs")
    print("=" * 164)
    print(f"Script version: {SCRIPT_VERSION}")
    print()
    print("Execution:")
    print("  Existing 05j mappings reused at .10,.25,.50,.75: YES")
    print("  New exact-generator fractions .05,.15:            100 valid each")
    print("  Full-corruption exact-generator null:              1,000 valid")
    print("  First 100 full-null replay guard:                  REQUIRED")
    print("  Same-edge Pearson diagnostic:                      YES")
    print("  Primary MTA method changed:                        NO")
    print("=" * 164)

    for p in [
        SOURCE_05D,
        SOURCE_05J,
        CONTRACT,
        HEADTOHEAD_REPLICATES,
        CLASSIFICATION,
        WEIGHTS,
    ]:
        require(p)

    contract = json.loads(CONTRACT.read_text(encoding="utf-8"))
    if contract.get("status") != "FROZEN_POST_05J3_BEFORE_SCALE_INVARIANT_AUDIT_V3_RESULTS":
        raise RuntimeError("05j4b v3 contract has unexpected status.")

    exact05d = load_module(SOURCE_05D, "paper4_05d_exact_v3")
    exact05j = load_module(SOURCE_05J, "paper4_05j_exact_v3")

    if not torch.cuda.is_available():
        raise RuntimeError("CUDA is unavailable.")
    device = torch.device("cuda:0")
    prop = torch.cuda.get_device_properties(0)
    print(f"CUDA device: {prop.name}; VRAM={prop.total_memory/1024**3:.2f} GB")

    existing = pd.read_csv(HEADTOHEAD_REPLICATES, sep="\t", low_memory=False)
    if len(existing) != 11022:
        raise RuntimeError(f"Expected 11,022 05j rows, found {len(existing)}")

    classification = pd.read_csv(CLASSIFICATION, sep="\t", low_memory=False)
    classification = classification.loc[
        classification["primary_assessable"].map(boolish)
    ].copy()
    if len(classification) != 22:
        raise RuntimeError(f"Expected 22 assessable pairs, found {len(classification)}")

    weights = pd.read_csv(WEIGHTS, sep="\t", dtype=str).fillna("")
    weights["Hugo_Symbol"] = weights["Hugo_Symbol"].map(canon_symbol)

    source = exact05j.load_tcga_source()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    COMBO_DIR.mkdir(parents=True, exist_ok=True)

    all_perf = []
    all_fullnull = []

    for target_name in ["SCANB_GSE96058", "METABRIC"]:
        print("\n" + "=" * 164)
        print(f"TARGET: {target_name}")
        print("=" * 164)

        if target_name == "SCANB_GSE96058":
            corrected = exact05j.load_eval_symbols(SCANB_EVAL)
            features = exact05j.load_features(SCANB_FEATURES, corrected)
            target = exact05j.load_scanb_target(corrected)
        else:
            corrected = exact05j.load_eval_symbols(METABRIC_EVAL)
            features = exact05j.load_features(METABRIC_FEATURES, corrected)
            target = exact05j.load_metabric_target(corrected)

        print("  computing full target Pearson correlation on GPU ...")
        target_corr = exact05j.corr_gpu(target["x"], device)
        del target["x"]
        torch.cuda.empty_cache()
        gc.collect()

        rows_target = classification.loc[classification["target"] == target_name].copy()

        for combo_i, class_row in enumerate(rows_target.to_dict(orient="records"), start=1):
            program_id = str(class_row["program_id"])
            print(f"\n  [{combo_i:02d}/{len(rows_target):02d}] {program_id}")

            perf_path = COMBO_DIR / f"{target_name}__{program_id}__performance_grid_v3.tsv"
            full_path = COMBO_DIR / f"{target_name}__{program_id}__fullnull_1000_v3.tsv"
            meta_path = COMBO_DIR / f"{target_name}__{program_id}__scale_input_meta_v3.json"

            if perf_path.exists() and full_path.exists() and meta_path.exists():
                meta = json.loads(meta_path.read_text(encoding="utf-8"))
                if meta.get("status") == "COMPLETE":
                    print("      completed combo found — reusing")
                    all_perf.append(pd.read_csv(perf_path, sep="\t", low_memory=False))
                    all_fullnull.append(pd.read_csv(full_path, sep="\t", low_memory=False))
                    continue

            gene_file = (
                DIRECT_ROOT / "per_program" / target_name
                / f"{program_id}_evaluable_genes_v2.tsv"
            )
            require(gene_file)
            eg = pd.read_csv(gene_file, sep="\t", dtype=str).fillna("")
            slot_genes = [canon_symbol(x) for x in eg["Hugo_Symbol"]]
            m = len(slot_genes)

            full_program_genes = {
                canon_symbol(x)
                for x in weights.loc[
                    weights["program_id"] == program_id,
                    "Hugo_Symbol",
                ]
            }

            correct_target_idx = np.array(
                [target["gene_index"][g] for g in slot_genes],
                dtype=np.int64,
            )

            null_path = (
                SPECIFICITY_ROOT / "per_program" / target_name
                / f"{program_id}_mapping_specificity_null_v2.npz"
            )
            require(null_path)
            nz = np.load(null_path)
            ii = np.asarray(nz["edge_slot_i_0based"], dtype=np.int64)
            jj = np.asarray(nz["edge_slot_j_0based"], dtype=np.int64)
            source_subset = np.asarray(nz["source_edge_subset"], dtype=np.float64)
            source_rank, source_rank_norm = exact05d.centered_rank(source_subset)

            ii_t = torch.as_tensor(ii, dtype=torch.long, device=device)
            jj_t = torch.as_tensor(jj, dtype=torch.long, device=device)

            source_cols = np.array(
                [source["gene_index"][g] for g in slot_genes],
                dtype=np.int64,
            )
            source_corr = exact05j.corr_gpu(source["x"][:, source_cols], device)

            tri_i_np, tri_j_np = np.triu_indices(m, k=1)
            tri_i_t = torch.as_tensor(tri_i_np, dtype=torch.long, device=device)
            tri_j_t = torch.as_tensor(tri_j_np, dtype=torch.long, device=device)

            source_all_edges = source_corr[tri_i_t, tri_j_t]
            source_all_centered, source_all_norm = exact05j.fixed_pearson_vector(source_all_edges)

            source_adj12 = torch.pow((1.0 + source_corr) / 2.0, 12.0)
            source_kim = torch.sum(source_adj12, dim=0)
            source_kim_centered, source_kim_norm = exact05j.fixed_pearson_vector(source_kim)

            # -----------------------------------------------------------------
            # Performance grid: reuse exact 05j .10/.25/.50/.75 plus add .05/.15.
            # -----------------------------------------------------------------
            saved = existing.loc[
                (existing["target"] == target_name)
                & (existing["program_id"] == program_id)
            ].copy()
            saved = saved.sort_values(["fraction", "replicate_id"]).reset_index(drop=True)

            perf_rows = []

            # Fraction zero from existing exact 05j.
            zero_row = saved.loc[saved["fraction"] == 0.0]
            if len(zero_row) != 1:
                raise RuntimeError("Existing 05j fraction-zero row-count drift.")
            z = zero_row.iloc[0]
            mapped0 = correct_target_idx.copy()
            sc0 = score_mapping(
                mapped0, target_corr, ii, jj, ii_t, jj_t,
                source_subset, source_rank, source_rank_norm,
                tri_i_t, tri_j_t, source_all_centered, source_all_norm,
                source_kim_centered, source_kim_norm,
                exact05d, exact05j, device,
            )
            if abs(sc0["mta_rho_edge"] - float(z["mta_rho_edge_saved"])) > REPLAY_TOL:
                raise RuntimeError("Fraction-zero MTA replay mismatch.")
            if abs(sc0["netrep_cor_cor"] - float(z["netrep_cor_cor"])) > COMPARATOR_TOL:
                raise RuntimeError("Fraction-zero NetRep replay mismatch.")
            if abs(sc0["wgcna_cor_kIM"] - float(z["wgcna_cor_kIM"])) > COMPARATOR_TOL:
                raise RuntimeError("Fraction-zero WGCNA replay mismatch.")

            perf_rows.append(
                {
                    "target": target_name,
                    "program_id": program_id,
                    "fraction": 0.0,
                    "replicate_id": 0,
                    "attempt_id": 0,
                    **sc0,
                }
            )

            # Existing partial fractions.
            for fraction in [0.10, 0.25, 0.50, 0.75]:
                old = saved.loc[saved["fraction"] == fraction].sort_values("replicate_id")
                if len(old) != 100:
                    raise RuntimeError(f"Existing fraction {fraction} row-count drift.")

                for row in old.to_dict(orient="records"):
                    mapped, _ = build_mapping(
                        exact05d, target_name, program_id, fraction,
                        int(row["attempt_id"]), slot_genes, full_program_genes,
                        features, target["gene_index"], correct_target_idx,
                    )
                    sc = score_mapping(
                        mapped, target_corr, ii, jj, ii_t, jj_t,
                        source_subset, source_rank, source_rank_norm,
                        tri_i_t, tri_j_t, source_all_centered, source_all_norm,
                        source_kim_centered, source_kim_norm,
                        exact05d, exact05j, device,
                    )

                    if abs(sc["mta_rho_edge"] - float(row["mta_rho_edge_saved"])) > REPLAY_TOL:
                        raise RuntimeError("Existing partial MTA replay mismatch.")
                    if abs(sc["netrep_cor_cor"] - float(row["netrep_cor_cor"])) > COMPARATOR_TOL:
                        raise RuntimeError("Existing partial NetRep replay mismatch.")
                    if abs(sc["wgcna_cor_kIM"] - float(row["wgcna_cor_kIM"])) > COMPARATOR_TOL:
                        raise RuntimeError("Existing partial WGCNA replay mismatch.")

                    perf_rows.append(
                        {
                            "target": target_name,
                            "program_id": program_id,
                            "fraction": fraction,
                            "replicate_id": int(row["replicate_id"]),
                            "attempt_id": int(row["attempt_id"]),
                            **sc,
                        }
                    )

            # New .05/.15 exact-generator fractions.
            for fraction in NEW_PARTIAL_FRACTIONS:
                valid = 0
                attempts = 0
                while valid < NEW_PARTIAL_REPLICATES and attempts < MAX_PARTIAL_ATTEMPTS:
                    attempts += 1
                    try:
                        mapped, _ = build_mapping(
                            exact05d, target_name, program_id, fraction, attempts,
                            slot_genes, full_program_genes, features,
                            target["gene_index"], correct_target_idx,
                        )
                    except RuntimeError:
                        continue

                    sc = score_mapping(
                        mapped, target_corr, ii, jj, ii_t, jj_t,
                        source_subset, source_rank, source_rank_norm,
                        tri_i_t, tri_j_t, source_all_centered, source_all_norm,
                        source_kim_centered, source_kim_norm,
                        exact05d, exact05j, device,
                    )
                    valid += 1
                    perf_rows.append(
                        {
                            "target": target_name,
                            "program_id": program_id,
                            "fraction": fraction,
                            "replicate_id": valid,
                            "attempt_id": attempts,
                            **sc,
                        }
                    )

                if valid != NEW_PARTIAL_REPLICATES:
                    raise RuntimeError(
                        f"{target_name} {program_id} f={fraction}: "
                        f"{valid}/{NEW_PARTIAL_REPLICATES} valid after {attempts} attempts"
                    )
                print(f"      new fraction {fraction:.2f}: valid=100, attempts={attempts}")

            perf_df = pd.DataFrame(perf_rows)
            perf_df = perf_df.sort_values(["fraction", "replicate_id"]).reset_index(drop=True)
            if len(perf_df) != 601:
                raise RuntimeError(f"Performance-grid rows={len(perf_df)}, expected 601")
            perf_df.to_csv(perf_path, sep="\t", index=False)

            # -----------------------------------------------------------------
            # Extended full-corruption null: first 100 must replay 05j exactly.
            # -----------------------------------------------------------------
            saved_full = saved.loc[saved["fraction"] == 1.0].sort_values("replicate_id").reset_index(drop=True)
            if len(saved_full) != FIRST_REPLAY:
                raise RuntimeError("Existing full-null rows are not exactly 100.")

            full_rows = []
            attempts = 0
            valid = 0

            while valid < FULL_VALID and attempts < MAX_ATTEMPTS:
                attempts += 1

                try:
                    mapped, _ = build_mapping(
                        exact05d, target_name, program_id, 1.0, attempts,
                        slot_genes, full_program_genes, features,
                        target["gene_index"], correct_target_idx,
                    )
                except RuntimeError:
                    continue

                sc = score_mapping(
                    mapped, target_corr, ii, jj, ii_t, jj_t,
                    source_subset, source_rank, source_rank_norm,
                    tri_i_t, tri_j_t, source_all_centered, source_all_norm,
                    source_kim_centered, source_kim_norm,
                    exact05d, exact05j, device,
                )

                valid += 1

                if valid <= FIRST_REPLAY:
                    exp = saved_full.iloc[valid - 1]
                    if attempts != int(exp["attempt_id"]):
                        raise RuntimeError(
                            f"First-100 attempt replay mismatch: valid={valid}, "
                            f"attempt={attempts}, expected={int(exp['attempt_id'])}"
                        )
                    if abs(sc["mta_rho_edge"] - float(exp["mta_rho_edge_saved"])) > REPLAY_TOL:
                        raise RuntimeError("First-100 MTA full-null replay mismatch.")
                    if abs(sc["netrep_cor_cor"] - float(exp["netrep_cor_cor"])) > COMPARATOR_TOL:
                        raise RuntimeError("First-100 NetRep full-null replay mismatch.")
                    if abs(sc["wgcna_cor_kIM"] - float(exp["wgcna_cor_kIM"])) > COMPARATOR_TOL:
                        raise RuntimeError("First-100 WGCNA full-null replay mismatch.")

                full_rows.append(
                    {
                        "target": target_name,
                        "program_id": program_id,
                        "valid_null_id": valid,
                        "attempt_id": attempts,
                        **sc,
                    }
                )

                if valid % 100 == 0:
                    print(
                        f"      full null valid={valid:4d}/{FULL_VALID}; attempts={attempts}"
                    )

            if valid != FULL_VALID:
                raise RuntimeError(
                    f"{target_name} {program_id}: only {valid}/{FULL_VALID} full-null valid "
                    f"after {attempts} attempts"
                )

            full_df = pd.DataFrame(full_rows)
            full_df.to_csv(full_path, sep="\t", index=False)

            meta = {
                "script_version": SCRIPT_VERSION,
                "status": "COMPLETE",
                "target": target_name,
                "program_id": program_id,
                "performance_grid_rows": int(len(perf_df)),
                "new_fractions": NEW_PARTIAL_FRACTIONS,
                "full_null_valid_rows": int(len(full_df)),
                "full_null_attempts": int(attempts),
                "first_100_existing_replay": "PASS",
                "performance_file": str(perf_path),
                "full_null_file": str(full_path),
            }
            meta_path.write_text(json.dumps(meta, indent=2), encoding="utf-8")

            all_perf.append(perf_df)
            all_fullnull.append(full_df)

            del (
                source_corr,
                tri_i_t,
                tri_j_t,
                source_all_edges,
                source_all_centered,
                source_all_norm,
                source_adj12,
                source_kim,
                source_kim_centered,
                source_kim_norm,
                ii_t,
                jj_t,
            )
            torch.cuda.empty_cache()
            gc.collect()

        del target_corr, target, features
        torch.cuda.empty_cache()
        gc.collect()

    perf_all = pd.concat(all_perf, ignore_index=True)
    full_all = pd.concat(all_fullnull, ignore_index=True)

    perf_out = OUT_DIR / "scale_invariant_performance_grid_v3.tsv"
    full_out = OUT_DIR / "scale_invariant_fullnull_1000_v3.tsv"

    perf_all.to_csv(perf_out, sep="\t", index=False)
    full_all.to_csv(full_out, sep="\t", index=False)

    master = {
        "script_version": SCRIPT_VERSION,
        "status": "SCALE_INVARIANT_INPUTS_V3_COMPLETE",
        "target_program_pairs": 22,
        "performance_grid_rows": int(len(perf_all)),
        "performance_rows_per_pair": 601,
        "full_null_rows": int(len(full_all)),
        "full_null_valid_per_pair": FULL_VALID,
        "new_partial_fractions": NEW_PARTIAL_FRACTIONS,
        "performance_grid_file": str(perf_out),
        "full_null_file": str(full_out),
        "primary_mta_method_changed": False,
    }
    master_path = OUT_DIR / "scale_invariant_inputs_v3.json"
    master_path.write_text(json.dumps(master, indent=2), encoding="utf-8")

    print("\n" + "=" * 164)
    print("05j4b SCALE-INVARIANT INPUT PREPARATION v3: COMPLETE")
    print("=" * 164)
    print(f"Performance-grid rows:             {len(perf_all):,}")
    print(f"Extended full-null rows:           {len(full_all):,}")
    print("First 100 full-null replay guard:  PASS")
    print("Primary MTA method changed:        NO")
    print(f"Performance: {perf_out}")
    print(f"Full null:   {full_out}")
    print(f"Master:      {master_path}")
    print("=" * 164)


if __name__ == "__main__":
    main()
