# Legacy Scoring Contract Audit

> This report extracts scoring-, mapping-, weighting-, coverage-, and standardization-related logic from authoritative Paper 4 scripts.

It is descriptive only. No reusable scoring contract is defined until these legacy semantics have been reviewed.

- Tool version: `0.1.0`
- Generated: `2026-08-07T23:19:44.102402+00:00`

---

## `scripts/21_finalize_canine_transfer_programs_FIXED_V2.py`

- SHA-256: `8da88bbb0fef845dc2c964c24eee21f4382c90e86b98f55fb0b66b3c7af5205c`

### Relevant functions

#### `bootstrap_mean_ci`

- Lines: 177-185
- Signals: `mean`

```python
def bootstrap_mean_ci(values, reps=BOOTSTRAP_REPS, seed=RANDOM_SEED):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if len(values) == 0:
        return np.nan, np.nan
    rng = np.random.default_rng(seed)
    sampled = rng.choice(values, size=(reps, len(values)), replace=True)
    means = sampled.mean(axis=1)
    return float(np.quantile(means, 0.025)), float(np.quantile(means, 0.975))
```

#### `corrected_repeated_kfold_ttest`

- Lines: 188-218
- Signals: `mean`

```python
def corrected_repeated_kfold_ttest(fold_deltas):
    values = np.asarray(fold_deltas, dtype=float)
    values = values[np.isfinite(values)]
    result = {
        "corrected_t": np.nan,
        "corrected_t_df": np.nan,
        "corrected_t_p_two_sided": np.nan,
        "corrected_standard_error": np.nan,
        "correction_factor": np.nan,
    }
    if len(values) < 2:
        return result
    variance = np.var(values, ddof=1)
    test_train_ratio = 1.0 / (N_SPLITS - 1)
    correction = (1.0 / len(values)) + test_train_ratio
    standard_error = np.sqrt(correction * variance)
    if not np.isfinite(standard_error) or standard_error == 0:
        return result
    t_stat = float(np.mean(values) / standard_error)
    df = len(values) - 1
    p_value = float(2 * stats.t.sf(abs(t_stat), df=df))
    result.update(
        {
            "corrected_t": t_stat,
            "corrected_t_df": df,
            "corrected_t_p_two_sided": p_value,
            "corrected_standard_error": float(standard_error),
            "correction_factor": float(correction),
        }
    )
    return result
```

#### `repeat_level_scores`

- Lines: 221-238
- Signals: ``

```python
def repeat_level_scores(cv_folds):
    use = cv_folds.copy()
    use["c_index"] = pd.to_numeric(use["c_index"], errors="coerce")
    use = use[np.isfinite(use["c_index"])].copy()
    return (
        use.groupby(
            ["endpoint", "module_label", "model", "repeat"],
            dropna=False,
        )
        .agg(
            n_valid_folds=("c_index", "count"),
            repeat_mean_c_index=("c_index", "mean"),
            repeat_median_c_index=("c_index", "median"),
            repeat_min_c_index=("c_index", "min"),
            repeat_max_c_index=("c_index", "max"),
        )
        .reset_index()
    )
```

#### `summarize_paired_values`

- Lines: 241-262
- Signals: `mean, std`

```python
def summarize_paired_values(repeat_deltas, fold_deltas, seed):
    repeat_deltas = np.asarray(repeat_deltas, dtype=float)
    repeat_deltas = repeat_deltas[np.isfinite(repeat_deltas)]
    fold_deltas = np.asarray(fold_deltas, dtype=float)
    fold_deltas = fold_deltas[np.isfinite(fold_deltas)]
    ci_low, ci_high = bootstrap_mean_ci(repeat_deltas, seed=seed)
    corrected = corrected_repeated_kfold_ttest(fold_deltas)
    return {
        "n_repeats": len(repeat_deltas),
        "n_fold_pairs": len(fold_deltas),
        "mean_delta_c_index": float(np.mean(repeat_deltas)) if len(repeat_deltas) else np.nan,
        "median_delta_c_index": float(np.median(repeat_deltas)) if len(repeat_deltas) else np.nan,
        "std_delta_c_index": float(np.std(repeat_deltas, ddof=1)) if len(repeat_deltas) > 1 else np.nan,
        "bootstrap_repeat_mean_ci_low": ci_low,
        "bootstrap_repeat_mean_ci_high": ci_high,
        "fraction_repeat_deltas_positive": float((repeat_deltas > 0).mean()) if len(repeat_deltas) else np.nan,
        "wilcoxon_repeat_p_two_sided": safe_wilcoxon(repeat_deltas),
        "wilcoxon_repeat_p_greater": safe_wilcoxon(repeat_deltas, alternative="greater"),
        "sign_repeat_p_two_sided": safe_sign_test(repeat_deltas),
        "sign_repeat_p_greater": safe_sign_test(repeat_deltas, alternative="greater"),
        **corrected,
    }
```

#### `build_model_contrasts`

- Lines: 265-391
- Signals: `mean`

```python
def build_model_contrasts(cv_folds, repeat_scores):
    rows = []
    seed_counter = 0
    endpoints = sorted(cv_folds["endpoint"].dropna().astype(str).unique())
    modules = sorted(cv_folds["module_label"].dropna().astype(str).unique())

    for endpoint in endpoints:
        for module_label in modules:
            fold_part = cv_folds[
                cv_folds["endpoint"].astype(str).eq(endpoint)
                & cv_folds["module_label"].astype(str).eq(module_label)
            ].copy()
            repeat_part = repeat_scores[
                repeat_scores["endpoint"].astype(str).eq(endpoint)
                & repeat_scores["module_label"].astype(str).eq(module_label)
            ].copy()

            for spec in MODEL_CONTRASTS:
                model_a = spec["model_a"]
                model_b = spec["model_b"]

                fold_a = fold_part[fold_part["model"].eq(model_a)][
                    ["repeat", "fold", "c_index"]
                ].rename(columns={"c_index": "c_index_a"})
                fold_b = fold_part[fold_part["model"].eq(model_b)][
                    ["repeat", "fold", "c_index"]
                ].rename(columns={"c_index": "c_index_b"})
                fold_pairs = fold_a.merge(fold_b, on=["repeat", "fold"], how="inner")
                fold_pairs["delta"] = fold_pairs["c_index_a"] - fold_pairs["c_index_b"]

                repeat_a = repeat_part[repeat_part["model"].eq(model_a)][
                    ["repeat", "repeat_mean_c_index"]
                ].rename(columns={"repeat_mean_c_index": "repeat_c_index_a"})
                repeat_b = repeat_part[repeat_part["model"].eq(model_b)][
                    ["repeat", "repeat_mean_c_index"]
                ].rename(columns={"repeat_mean_c_index": "repeat_c_index_b"})
                repeat_pairs = repeat_a.merge(repeat_b, on="repeat", how="inner")
                repeat_pairs["delta"] = (
                    repeat_pairs["repeat_c_index_a"]
                    - repeat_pairs["repeat_c_index_b"]
                )

                if repeat_pairs.empty:
                    continue

                stats_row = summarize_paired_values(
                    repeat_deltas=repeat_pairs["delta"],
                    fold_deltas=fold_pairs["delta"],
                    seed=RANDOM_SEED + seed_counter,
                )
                seed_counter += 1
                rows.append(
                    {
                        "endpoint": endpoint,
                        "module_label": module_label,
                        "contrast": spec["contrast"],
                        "model_a": model_a,
                        "model_b": model_b,
                        "mean_repeat_c_index_a": float(repeat_pairs["repeat_c_index_a"].mean()),
                        "mean_repeat_c_index_b": float(repeat_pairs["repeat_c_index_b"].mean()),
                        **stats_row,
                    }
                )

            residual = repeat_part[
                repeat_part["model"].eq("residual_to_disjoint_proliferation")
            ].copy()
            residual_folds = fold_part[
                fold_part["model"].eq("residual_to_disjoint_proliferation")
            ].copy()
            if not residual.empty:
                repeat_delta = residual["repeat_mean_c_index"] - 0.5
                fold_delta = pd.to_numeric(
                    residual_folds["c_index"], errors="coerce"
                ) - 0.5
                stats_row = summarize_paired_values(
                    repeat_deltas=repeat_delta,
                    fold_deltas=fold_delta,
                    seed=RANDOM_SEED + seed_counter,
                )
                seed_counter += 1
                rows.append(
                    {
                        "endpoint": endpoint,
                        "module_label": module_label,
                        "contrast": "disjoint_residual_above_chance",
                        "model_a": "residual_to_disjoint_proliferation",
                        "model_b": "chance_0.50",
                        "mean_repeat_c_index_a": float(residual["repeat_mean_c_index"].mean()),
                        "mean_repeat_c_index_b": 0.5,
                        **stats_row,
                    }
                )

    contrasts = pd.DataFrame(rows)
    if contrasts.empty:
        return contrasts

    for p_col in [
        "wilcoxon_repeat_p_two_sided",
        "wilcoxon_repeat_p_greater",
        "corrected_t_p_two_sided",
    ]:
        contrasts[f"{p_col}_bh_q"] = np.nan
        for endpoint in contrasts["endpoint"].dropna().unique():
            mask = contrasts["endpoint"].eq(endpoint)
            contrasts.loc[mask, f"{p_col}_bh_q"] = bh_adjust(
                contrasts.loc[mask, p_col]
            )

    contrasts["descriptive_support"] = np.select(
        [
            (contrasts["mean_delta_c_index"] >= 0.02)
            & (contrasts["fraction_repeat_deltas_positive"] >= 0.75)
            & (contrasts["bootstrap_repeat_mean_ci_low"] > 0),
            (contrasts["mean_delta_c_index"] >= 0.01)
            & (contrasts["fraction_repeat_deltas_positive"] >= 0.60),
            (contrasts["mean_delta_c_index"] > 0),
        ],
        [
            "strong_descriptive_support",
            "moderate_descriptive_support",
            "weak_positive_descriptive_support",
        ],
        default="no_consistent_positive_support",
    )
    return contrasts
```

#### `detect_module_columns`

- Lines: 410-417
- Signals: ``

```python
def detect_module_columns(module_membership):
    module_candidates = ["module_label", "module", "module_id", "cluster", "cluster_id"]
    gene_candidates = ["gene", "gene_id", "expression_gene", "gene_column"]
    module_col = next((c for c in module_candidates if c in module_membership.columns), None)
    gene_col = next((c for c in gene_candidates if c in module_membership.columns), None)
    if module_col is None or gene_col is None:
        raise ValueError("Could not identify module and gene columns in module membership.")
    return module_col, gene_col
```

#### `get_module_score_mapping`

- Lines: 435-443
- Signals: ``

```python
def get_module_score_mapping(module_scores):
    mapping = {}
    for col in module_scores.columns:
        value = str(col)
        if value.startswith("module_") and value.endswith("_score"):
            mapping[value[len("module_") : -len("_score")]] = col
        else:
            mapping[value] = col
    return mapping
```

#### `build_module_gene_map`

- Lines: 457-478
- Signals: ``

```python
def build_module_gene_map(module_membership, expression, module_scores):
    membership = filter_full_cohort_membership(module_membership)
    module_col, gene_col = detect_module_columns(membership)
    score_mapping = get_module_score_mapping(module_scores)
    valid_labels = set(score_mapping)
    membership[module_col] = membership[module_col].astype(str)
    membership = membership[membership[module_col].isin(valid_labels)].copy()

    expression_columns = set(expression.columns.astype(str))
    symbol_to_columns = {}
    for col in expression.columns:
        symbol_to_columns.setdefault(clean_gene_symbol(col).upper(), []).append(col)

    module_gene_map = {}
    for module_label, part in membership.groupby(module_col):
        genes = []
        for value in part[gene_col].dropna().astype(str):
            resolved = resolve_expression_gene(value, expression_columns, symbol_to_columns)
            if resolved:
                genes.append(resolved)
        module_gene_map[str(module_label)] = list(dict.fromkeys(genes))
    return module_gene_map, score_mapping
```

#### `detect_ortholog_columns`

- Lines: 481-496
- Signals: ``

```python
def detect_ortholog_columns(ortholog):
    dog_candidates = ["gene", "canine_gene", "dog_gene", "expression_gene"]
    human_candidates = ["human_gene_symbol", "human_symbol", "human_gene"]
    status_candidates = ["ortholog_qc_status", "mapping_status", "ortholog_status"]
    dog_col = next((c for c in dog_candidates if c in ortholog.columns), None)
    human_col = next((c for c in human_candidates if c in ortholog.columns), None)
    status_col = next((c for c in status_candidates if c in ortholog.columns), None)
    if dog_col is None or human_col is None:
        raise ValueError("Could not identify dog and human ortholog columns.")
    if status_col is None:
        raise ValueError(
            "The ortholog QC status column is missing. Run "
            "scripts/17_ortholog_mapping_qc_transfer_sets.py and use "
            "GSE238110_RNA_master_candidate_evidence_table_with_ortholog_qc.csv."
        )
    return dog_col, human_col, status_col
```

#### `standardize_expression`

- Lines: 499-507
- Signals: `mean, std, zscore`

```python
def standardize_expression(expression, genes):
    x = expression[genes].apply(pd.to_numeric, errors="coerce")
    x = x.replace([np.inf, -np.inf], np.nan)
    x = x.fillna(x.median(axis=0))
    means = x.mean(axis=0)
    stds = x.std(axis=0).replace(0, np.nan)
    z = (x - means) / stds
    valid = z.columns[z.notna().all(axis=0)]
    return z[valid]
```

#### `safe_corr`

- Lines: 510-514
- Signals: `std, correlation`

```python
def safe_corr(a, b):
    frame = pd.concat([pd.Series(a, name="a"), pd.Series(b, name="b")], axis=1).dropna()
    if frame.shape[0] < 5 or frame["a"].std() == 0 or frame["b"].std() == 0:
        return np.nan
    return float(frame["a"].corr(frame["b"]))
```

#### `endpoint_module_row`

- Lines: 517-532
- Signals: ``

```python
def endpoint_module_row(module_associations, module_label, endpoint):
    if module_associations.empty:
        return pd.Series(dtype=object)
    use = module_associations.copy()
    use["module_label"] = use["module_label"].astype(str)
    use["endpoint"] = use["endpoint"].astype(str).str.upper()
    part = use[
        use["module_label"].eq(str(module_label))
        & use["endpoint"].eq(str(endpoint).upper())
    ].copy()
    if part.empty:
        return pd.Series(dtype=object)
    if "p" in part.columns:
        part["p"] = pd.to_numeric(part["p"], errors="coerce")
        part = part.sort_values("p", na_position="last")
    return part.iloc[0]
```

#### `build_weights_and_manifest`

- Lines: 554-778
- Signals: `zscore, loadings`

```python
def build_weights_and_manifest(
    expression,
    module_gene_map,
    module_scores,
    module_associations,
    priority,
    ortholog,
    audit,
    cv_summary,
    decision,
):
    dog_col, human_col, status_col = detect_ortholog_columns(ortholog)
    ortholog = ortholog.copy()
    ortholog[dog_col] = ortholog[dog_col].astype(str)
    ortholog["dog_symbol_key"] = ortholog[dog_col].map(clean_gene_symbol).str.upper()
    ortholog[human_col] = ortholog[human_col].astype(str).replace({"nan": "", "None": ""})

    print("")
    print("Ortholog mapping input:")
    print(f"  Dog gene column: {dog_col}")
    print(f"  Human gene column: {human_col}")
    print(f"  QC status column: {status_col}")
    print("  QC status counts:")
    print(ortholog[status_col].fillna("<missing>").value_counts(dropna=False).to_string())

    score_mapping = get_module_score_mapping(module_scores)
    weight_rows = []
    manifest_rows = []

    for module_label, genes in sorted(module_gene_map.items()):
        genes = [g for g in genes if g in expression.columns]
        if len(genes) < 3:
            continue

        z = standardize_expression(expression, genes)
        genes_used = list(z.columns)
        if len(genes_used) < 3:
            continue

        pca = PCA(n_components=1, random_state=RANDOM_SEED)
        score = pd.Series(pca.fit_transform(z).ravel(), index=z.index)
        raw_loadings = pd.Series(pca.components_[0], index=genes_used)

        score_col = score_mapping.get(module_label)
        orientation_corr = np.nan
        if score_col in module_scores.columns:
            reference = module_scores[score_col].reindex(score.index)
            orientation_corr = safe_corr(score, reference)
            if np.isfinite(orientation_corr) and orientation_corr < 0:
                score = -score
                raw_loadings = -raw_loadings

        dfi_row = endpoint_module_row(module_associations, module_label, "DFI")
        os_row = endpoint_module_row(module_associations, module_label, "OS")
        dfi_coef = get_numeric(dfi_row, ["coef", "module_coef"])
        os_coef = get_numeric(os_row, ["coef", "module_coef"])
        if not np.isfinite(dfi_coef):
            dfi_hr = get_numeric(dfi_row, ["hr", "hazard_ratio", "exp(coef)"])
            if np.isfinite(dfi_hr) and dfi_hr > 0:
                dfi_coef = float(np.log(dfi_hr))
        if not np.isfinite(os_coef):
            os_hr = get_numeric(os_row, ["hr", "hazard_ratio", "exp(coef)"])
            if np.isfinite(os_hr) and os_hr > 0:
                os_coef = float(np.log(os_hr))
        risk_reference_coef = dfi_coef if np.isfinite(dfi_coef) else os_coef
        risk_multiplier = 1.0 if not np.isfinite(risk_reference_coef) or risk_reference_coef >= 0 else -1.0
        risk_loadings = raw_loadings * risk_multiplier

        module_weight_rows = []
        for canine_gene, loading in risk_loadings.items():
            dog_symbol = clean_gene_symbol(canine_gene).upper()
            mappings = ortholog[ortholog["dog_symbol_key"].eq(dog_symbol)].copy()
            if mappings.empty:
                module_weight_rows.append(
                    {
                        "module_label": module_label,
                        "canine_gene": canine_gene,
                        "canine_gene_symbol": dog_symbol,
                        "human_gene_symbol": "",
                        "ortholog_qc_status": "not_found_in_ortholog_table",
                        "raw_pca_loading": float(raw_loadings.loc[canine_gene]),
                        "risk_oriented_loading": float(loading),
                    }
                )
                continue
            for _, mapping in mappings.iterrows():
                human_symbol = str(mapping[human_col]).strip().upper()
                module_weight_rows.append(
                    {
                        "module_label": module_label,
                        "canine_gene": canine_gene,
                        "canine_gene_symbol": dog_symbol,
                        "human_gene_symbol": human_symbol,
                        "ortholog_qc_status": str(mapping[status_col]),
                        "raw_pca_loading": float(raw_loadings.loc[canine_gene]),
                        "risk_oriented_loading": float(loading),
                    }
                )

        module_weights = pd.DataFrame(module_weight_rows)
        module_weights["absolute_risk_loading"] = module_weights["risk_oriented_loading"].abs()
        module_weights["is_strict_mapping"] = module_weights["ortholog_qc_status"].eq(
            "strict_symbol_concordant_one_to_one"
        ) & module_weights["human_gene_symbol"].ne("")
        module_weights["is_broad_mapping"] = (
            module_weights["human_gene_symbol"].ne("")
            & ~module_weights["ortholog_qc_status"].eq("not_transferable_or_unmapped")
            & ~module_weights["ortholog_qc_status"].eq("not_found_in_ortholog_table")
        )
        weight_rows.append(module_weights)

        strict = module_weights[module_weights["is_strict_mapping"]].copy()
        broad = module_weights[module_weights["is_broad_mapping"]].copy()

        priority_row = pd.Series(dtype=object)
        if not priority.empty and "module_label" in priority.columns:
            part = priority[priority["module_label"].astype(str).eq(module_label)]
            if not part.empty:
                priority_row = part.iloc[0]

        audit_row = pd.Series(dtype=object)
        if not audit.empty:
            part = audit[audit["module_label"].astype(str).eq(module_label)]
            if not part.empty:
                audit_row = part.iloc[0]

        manifest_rows.append(
            {
                "module_label": module_label,
                "validation_tier": assign_validation_tier(module_label),
                "multiplicity_family": (
                    "primary_confirmatory"
                    if module_label in PRIMARY_CLEAN_MODULES + PRIMARY_PROLIFERATION_AXIS_MODULES
                    else "secondary_prespecified"
                    if module_label in SECONDARY_SENSITIVITY_MODULES
                    else "exploratory"
                ),
                "provisional_program_label": PROVISIONAL_PROGRAM_LABELS.get(
                    module_label, f"exploratory_program_{module_label}"
                ),
                "program_label_requires_enrichment_confirmation": True,
                "canine_primary_endpoint": "DFI",
                "canine_secondary_endpoint": "OS_concordance",
                "positive_score_interpretation": "higher_score_higher_canine_DFI_risk",
                "n_canine_genes_used_for_pca": len(genes_used),
                "canine_pc1_explained_variance": float(pca.explained_variance_ratio_[0]),
                "canine_pca_orientation_correlation": orientation_corr,
                "dfi_full_cohort_coef": dfi_coef,
                "os_full_cohort_coef": os_coef,
                "risk_orientation_multiplier": risk_multiplier,
                "n_strict_human_genes": strict["human_gene_symbol"].replace("", np.nan).nunique(),
                "n_broad_human_genes": broad["human_gene_symbol"].replace("", np.nan).nunique(),
                "strict_transfer_eligible": strict["human_gene_symbol"].replace("", np.nan).nunique() >= MIN_STRICT_HUMAN_GENES,
                "module_transfer_qc_tier": priority_row.get("module_transfer_qc_tier", ""),
                "transfer_priority_score": priority_row.get("transfer_priority_score", np.nan),
                "fraction_strict_symbol_concordant": priority_row.get("fraction_strict_symbol_concordant", np.nan),
                "fraction_broad_transferable": priority_row.get("fraction_broad_transferable", np.nan),
                "raw_module_proliferation_correlation": audit_row.get("raw_module_proliferation_correlation", np.nan),
                "orthogonal_variance_fraction_1_minus_r2": audit_row.get("orthogonal_variance_fraction_1_minus_r2", np.nan),
                "n_overlap_symbols_with_proliferation": audit_row.get("n_overlap_symbols", np.nan),
                "primary_human_score": "strict_one_to_one_signed_mean_z",
                "secondary_human_score": "strict_one_to_one_canine_pca_weighted_z",
                "sensitivity_human_scores": "broad_mapped_mean_z;human_cohort_pc1;residual_to_disjoint_proliferation",
                "frozen_after_canine_script": "20_proliferation_overlap_crossfit_sensitivity.py",
            }
        )

    all_weights = pd.concat(weight_rows, axis=0, ignore_index=True)
    manifest = pd.DataFrame(manifest_rows)

    for endpoint in ["DFI", "OS"]:
        if not cv_summary.empty:
            part = cv_summary.copy()
            part["endpoint"] = part["endpoint"].astype(str).str.upper()
            part = part[part["endpoint"].eq(endpoint)].copy()
            for model in [
                "module_only",
                "module_plus_disjoint_proliferation",
                "residual_to_disjoint_proliferation",
                "residual_to_disjoint_proliferation_and_weight",
            ]:
                model_part = part[part["model"].eq(model)][
                    ["module_label", "mean_c_index", "std_c_index", "fraction_above_0_50"]
                ].copy()
                model_part = model_part.rename(
                    columns={
                        "mean_c_index": f"{endpoint.lower()}_{model}_mean_c_index",
                        "std_c_index": f"{endpoint.lower()}_{model}_std_c_index",
                        "fraction_above_0_50": f"{endpoint.lower()}_{model}_fraction_above_0_50",
                    }
                )
                manifest = manifest.merge(model_part, on="module_label", how="left")

    if not decision.empty:
        dfi_decision = decision[decision["endpoint"].astype(str).str.upper().eq("DFI")][
            ["module_label", "recommended_role"]
        ].rename(columns={"recommended_role": "script20_dfi_recommended_role"})
        os_decision = decision[decision["endpoint"].astype(str).str.upper().eq("OS")][
            ["module_label", "recommended_role"]
        ].rename(columns={"recommended_role": "script20_os_recommended_role"})
        manifest = manifest.merge(dfi_decision, on="module_label", how="left")
        manifest = manifest.merge(os_decision, on="module_label", how="left")

    manifest["manual_freeze_reason"] = np.select(
        [
            manifest["module_label"].eq("M34"),
            manifest["module_label"].eq("M11"),
            manifest["module_label"].eq("M24"),
            manifest["module_label"].eq("M40"),
            manifest["module_label"].isin(SECONDARY_SENSITIVITY_MODULES),
        ],
        [
            "strongest clean cross-endpoint non-proliferation program",
            "high transfer readiness and independent DFI signal",
            "compact high-readiness program with cross-endpoint signal",
            "proliferation-dominant axis with reproducible cross-fitted residual component",
            "prespecified sensitivity program based on performance or transferability",
        ],
        default="exploratory program retained without confirmatory status",
    )

    manifest = manifest.sort_values(
        ["multiplicity_family", "validation_tier", "module_label"]
    ).reset_index(drop=True)
    return all_weights, manifest
```

#### `deduplicate_human_weights`

- Lines: 781-797
- Signals: `weights`

```python
def deduplicate_human_weights(weights, mapping_flag):
    use = weights[weights[mapping_flag]].copy()
    use = use[use["human_gene_symbol"].astype(str).str.len() > 0].copy()
    if use.empty:
        return use
    use["human_symbol_duplicate_count"] = use.groupby(
        ["module_label", "human_gene_symbol"]
    )["human_gene_symbol"].transform("size")
    use = use.sort_values(
        ["module_label", "human_gene_symbol", "absolute_risk_loading"],
        ascending=[True, True, False],
    )
    use = use.drop_duplicates(["module_label", "human_gene_symbol"], keep="first")
    use["normalized_abs_sum_weight"] = use.groupby("module_label")[
        "risk_oriented_loading"
    ].transform(lambda x: x / x.abs().sum() if x.abs().sum() > 0 else np.nan)
    return use.reset_index(drop=True)
```

#### `write_gmt`

- Lines: 800-808
- Signals: `weights`

```python
def write_gmt(weights, path, suffix):
    lines = []
    for module_label, part in weights.groupby("module_label"):
        genes = part["human_gene_symbol"].dropna().astype(str).drop_duplicates().tolist()
        if not genes:
            continue
        name = f"CANINE_{module_label}_{suffix}"
        lines.append("\t".join([name, "frozen_canine_ortholog_transfer"] + genes))
    path.write_text("\n".join(lines), encoding="utf-8")
```

#### `build_scoring_specification`

- Lines: 811-865
- Signals: `zscore, weights, loadings, direction, coverage`

```python
def build_scoring_specification(manifest):
    rows = []
    for _, row in manifest.iterrows():
        rows.extend(
            [
                {
                    "module_label": row["module_label"],
                    "validation_tier": row["validation_tier"],
                    "score_name": "strict_one_to_one_unweighted_mean_z",
                    "analysis_role": "primary",
                    "gene_mapping": "strict_symbol_concordant_one_to_one",
                    "within_cohort_preprocessing": "z-score each available gene across human samples",
                    "score_formula": "arithmetic mean of signed gene z-scores using the sign of each frozen risk-oriented canine PCA loading",
                    "minimum_gene_rule": "at least 3 genes and at least 50% of frozen strict genes",
                    "outcome_use": "no outcome information used to construct the score",
                },
                {
                    "module_label": row["module_label"],
                    "validation_tier": row["validation_tier"],
                    "score_name": "strict_one_to_one_canine_pca_weighted_z",
                    "analysis_role": "secondary_zero_shot",
                    "gene_mapping": "strict_symbol_concordant_one_to_one",
                    "within_cohort_preprocessing": "z-score each available gene across human samples",
                    "score_formula": "sum of gene z-scores multiplied by frozen risk-oriented canine PCA loadings normalized by absolute weight sum",
                    "minimum_gene_rule": "at least 3 genes and at least 50% of frozen strict genes",
                    "outcome_use": "no outcome information used to construct the score",
                },
                {
                    "module_label": row["module_label"],
                    "validation_tier": row["validation_tier"],
                    "score_name": "broad_mapped_unweighted_mean_z",
                    "analysis_role": "mapping_sensitivity",
                    "gene_mapping": "broad transferable mapping including review-status mappings",
                    "within_cohort_preprocessing": "z-score each available gene across human samples",
                    "score_formula": "arithmetic mean of risk-oriented mapped gene z-scores",
                    "minimum_gene_rule": "at least 3 genes and at least 50% of frozen broad genes",
                    "outcome_use": "no outcome information used to construct the score",
                },
            ]
        )
        if row["module_label"] == "M40":
            rows.append(
                {
                    "module_label": row["module_label"],
                    "validation_tier": row["validation_tier"],
                    "score_name": "residual_to_disjoint_proliferation",
                    "analysis_role": "mechanistic_sensitivity",
                    "gene_mapping": "strict module genes plus strict proliferation genes after removing overlap",
                    "within_cohort_preprocessing": "construct both scores without outcomes; residualize M40 score on disjoint proliferation score",
                    "score_formula": "standardized residual from linear regression of M40 score on disjoint proliferation score",
                    "minimum_gene_rule": "same minimum-gene rules for both component scores",
                    "outcome_use": "no outcome information used during residualization",
                }
            )
    return pd.DataFrame(rows)
```

#### `write_readme`

- Lines: 876-909
- Signals: `weights, loadings, direction`

```python
def write_readme(manifest):
    primary = manifest[manifest["multiplicity_family"].eq("primary_confirmatory")]
    secondary = manifest[manifest["multiplicity_family"].eq("secondary_prespecified")]
    lines = [
        "Frozen canine-to-human osteosarcoma transfer programs",
        "",
        "Canine primary endpoint: DFI.",
        "Canine OS is a concordance/sensitivity endpoint.",
        "Program definitions, gene membership, score orientation, and validation tiers are frozen after script 20.",
        "No human outcome may be used to change module membership, gene weights, score direction, or tier assignment.",
        "",
        "Primary confirmatory programs:",
    ]
    for _, row in primary.iterrows():
        lines.append(
            f"  {row['module_label']}: {row['provisional_program_label']} | {row['validation_tier']}"
        )
    lines.append("")
    lines.append("Secondary prespecified programs:")
    for _, row in secondary.iterrows():
        lines.append(
            f"  {row['module_label']}: {row['provisional_program_label']} | {row['validation_tier']}"
        )
    lines.extend(
        [
            "",
            "Primary human score: strict one-to-one signed mean z-score using frozen canine loading signs.",
            "Secondary zero-shot score: frozen canine PCA-weighted score.",
            "M40 residualized score is a mechanistic sensitivity analysis, not a replacement for raw M40 or proliferation scores.",
            "TARGET-OS and GSE21257 must be treated as external datasets; cohort-specific preprocessing may not use outcomes.",
            "External validation, not canine repeated-CV p-values, determines translational support.",
        ]
    )
    OUTPUT_README.write_text("\n".join(lines), encoding="utf-8")
```

#### `main`

- Lines: 912-1113
- Signals: `weights, direction, coverage`

```python
def main():
    print("=" * 80)
    print("Finalize canine transfer programs")
    print("=" * 80)
    print(f"Script version: {SCRIPT_VERSION}")
    print(f"Project root: {PROJECT_ROOT}")
    print(f"Processed directory: {PROCESSED_DIR}")
    print(f"Results directory: {RESULTS_DIR}")
    print("")
    print("Design:")
    print("  Summarize paired repeated-CV model contrasts without treating 100 folds as independent evidence.")
    print("  Freeze confirmatory and sensitivity module tiers before human outcome analysis.")
    print("  Export strict and broad human ortholog gene sets and frozen canine PCA weights.")
    print("")

    expression = read_required_csv(PROCESSED_DIR / EXPRESSION_FILE, index_col=0)
    module_membership = read_required_csv(RESULTS_DIR / MODULE_MEMBERSHIP_FILE)
    module_scores = read_required_csv(RESULTS_DIR / MODULE_SCORE_FILE, index_col=0)
    module_associations = read_required_csv(RESULTS_DIR / MODULE_ASSOCIATION_FILE)
    priority = read_required_csv(RESULTS_DIR / MODULE_PRIORITY_FILE)
    ortholog = read_required_csv(RESULTS_DIR / ORTHOLOG_FILE)
    audit = read_required_csv(RESULTS_DIR / OVERLAP_AUDIT_FILE)
    cv_folds = read_required_csv(RESULTS_DIR / CV_FOLD_FILE)
    cv_summary = read_required_csv(RESULTS_DIR / CV_SUMMARY_FILE)
    decision = read_required_csv(RESULTS_DIR / DECISION_FILE)

    common_samples = expression.index.intersection(module_scores.index)
    expression = expression.loc[common_samples].copy()
    module_scores = module_scores.loc[common_samples].copy()

    print("")
    print("Matched data:")
    print(f"  Expression matrix: {expression.shape}")
    print(f"  Module score matrix: {module_scores.shape}")
    print(f"  Repeated-CV fold rows: {cv_folds.shape[0]}")

    repeat_scores = repeat_level_scores(cv_folds)
    contrasts = build_model_contrasts(cv_folds, repeat_scores)
    contrast_summary = build_contrast_summary(contrasts)

    module_gene_map, _ = build_module_gene_map(
        module_membership=module_membership,
        expression=expression,
        module_scores=module_scores,
    )
    all_weights, manifest = build_weights_and_manifest(
        expression=expression,
        module_gene_map=module_gene_map,
        module_scores=module_scores,
        module_associations=module_associations,
        priority=priority,
        ortholog=ortholog,
        audit=audit,
        cv_summary=cv_summary,
        decision=decision,
    )

    strict_weights = deduplicate_human_weights(all_weights, "is_strict_mapping")
    broad_weights = deduplicate_human_weights(all_weights, "is_broad_mapping")

    if strict_weights.empty:
        raise RuntimeError(
            "No strict one-to-one human ortholog weights were created. "
            "This usually means the pre-QC ortholog table was used instead of "
            "GSE238110_RNA_master_candidate_evidence_table_with_ortholog_qc.csv."
        )

    primary_modules = set(PRIMARY_CLEAN_MODULES + PRIMARY_PROLIFERATION_AXIS_MODULES)
    strict_primary_counts = (
        strict_weights[strict_weights["module_label"].isin(primary_modules)]
        .groupby("module_label")["human_gene_symbol"]
        .nunique()
        .reindex(sorted(primary_modules), fill_value=0)
    )
    insufficient_primary = strict_primary_counts[
        strict_primary_counts < MIN_STRICT_HUMAN_GENES
    ]
    if not insufficient_primary.empty:
        raise RuntimeError(
            "Primary frozen modules have too few strict human orthologs: "
            + ", ".join(
                f"{module}={int(count)}"
                for module, count in insufficient_primary.items()
            )
        )

    print("")
    print("Strict human ortholog counts for primary frozen modules:")
    print(strict_primary_counts.to_string())

    scoring_spec = build_scoring_specification(manifest)

    repeat_scores.to_csv(OUTPUT_REPEAT_SCORES, index=False)
    contrasts.to_csv(OUTPUT_CONTRASTS, index=False)
    contrast_summary.to_csv(OUTPUT_CONTRAST_SUMMARY, index=False)
    manifest.to_csv(OUTPUT_MANIFEST, index=False)
    strict_weights.to_csv(OUTPUT_STRICT_WEIGHTS, index=False)
    broad_weights.to_csv(OUTPUT_BROAD_WEIGHTS, index=False)
    scoring_spec.to_csv(OUTPUT_SCORING_SPEC, index=False)
    write_gmt(strict_weights, OUTPUT_STRICT_GMT, "STRICT")
    write_gmt(broad_weights, OUTPUT_BROAD_GMT, "BROAD")
    write_readme(manifest)

    freeze = {
        "frozen_after_script": "20_proliferation_overlap_crossfit_sensitivity.py",
        "primary_clean_modules": PRIMARY_CLEAN_MODULES,
        "primary_proliferation_axis_modules": PRIMARY_PROLIFERATION_AXIS_MODULES,
        "secondary_sensitivity_modules": SECONDARY_SENSITIVITY_MODULES,
        "primary_canine_endpoint": "DFI",
        "secondary_canine_endpoint": "OS_concordance",
        "primary_human_score": "strict_one_to_one_signed_mean_z",
        "secondary_human_score": "strict_one_to_one_canine_pca_weighted_z",
        "m40_residual_role": "mechanistic_sensitivity",
        "ortholog_qc_input_file": ORTHOLOG_FILE,
        "ortholog_qc_input_sha256": sha256_file(RESULTS_DIR / ORTHOLOG_FILE),
        "files": {},
    }

    frozen_paths = [
        OUTPUT_REPEAT_SCORES,
        OUTPUT_CONTRASTS,
        OUTPUT_CONTRAST_SUMMARY,
        OUTPUT_MANIFEST,
        OUTPUT_STRICT_WEIGHTS,
        OUTPUT_BROAD_WEIGHTS,
        OUTPUT_STRICT_GMT,
        OUTPUT_BROAD_GMT,
        OUTPUT_SCORING_SPEC,
        OUTPUT_README,
    ]
    for path in frozen_paths:
        freeze["files"][path.name] = {
            "sha256": sha256_file(path),
            "size_bytes": path.stat().st_size,
        }
    OUTPUT_FREEZE_JSON.write_text(json.dumps(freeze, indent=2), encoding="utf-8")

    print("")
    print("=" * 80)
    print("Frozen transfer program manifest")
    print("=" * 80)
    manifest_cols = [
        "module_label",
        "validation_tier",
        "multiplicity_family",
        "provisional_program_label",
        "n_canine_genes_used_for_pca",
        "n_strict_human_genes",
        "n_broad_human_genes",
        "raw_module_proliferation_correlation",
        "dfi_residual_to_disjoint_proliferation_mean_c_index",
        "os_residual_to_disjoint_proliferation_mean_c_index",
        "manual_freeze_reason",
    ]
    manifest_cols = [c for c in manifest_cols if c in manifest.columns]
    print(manifest[manifest_cols].to_string(index=False))

    print("")
    print("=" * 80)
    print("Key repeated-CV contrasts for frozen programs")
    print("=" * 80)
    key = contrast_summary[
        contrast_summary["contrast"].isin(
            [
                "incremental_module_beyond_disjoint_proliferation",
                "disjoint_residual_above_chance",
                "weight_adjusted_residual_vs_unadjusted_residual",
            ]
        )
    ].copy()
    key_cols = [
        "endpoint",
        "module_label",
        "contrast",
        "mean_repeat_c_index_a",
        "mean_repeat_c_index_b",
        "mean_delta_c_index",
        "bootstrap_repeat_mean_ci_low",
        "bootstrap_repeat_mean_ci_high",
        "fraction_repeat_deltas_positive",
        "wilcoxon_repeat_p_greater",
        "corrected_t_p_two_sided",
        "descriptive_support",
    ]
    key_cols = [c for c in key_cols if c in key.columns]
    print(key[key_cols].to_string(index=False))

    print("")
    print("=" * 80)
    print("Interpretation guardrails")
    print("=" * 80)
    print("Repeated cross-validation splits overlap; contrast p-values are descriptive stability diagnostics, not external validation.")
    print("The frozen primary clean programs are M34, M11, and M24.")
    print("M40 is frozen as a separate proliferation-dominant deviation axis.")
    print("M28, M38, M25, and M17 are prespecified sensitivity programs.")
    print("No human outcome may be used to revise module membership, weights, score direction, or validation tier.")

    print("")
    print("Saved:")
    for path in frozen_paths + [OUTPUT_FREEZE_JSON]:
        print(path)
    print("Done.")
```

### Relevant standalone lines

| Line | Signals | Code |
|---:|---|---|
| 184 | mean | `means = sampled.mean(axis=1)` |
| 206 | mean | `t_stat = float(np.mean(values) / standard_error)` |
| 251 | mean | `"mean_delta_c_index": float(np.mean(repeat_deltas)) if len(repeat_deltas) else np.nan,` |
| 253 | std | `"std_delta_c_index": float(np.std(repeat_deltas, ddof=1)) if len(repeat_deltas) > 1 else np.nan,` |
| 256 | mean | `"fraction_repeat_deltas_positive": float((repeat_deltas > 0).mean()) if len(repeat_deltas) else np.nan,` |
| 323 | mean | `"mean_repeat_c_index_a": float(repeat_pairs["repeat_c_index_a"].mean()),` |
| 324 | mean | `"mean_repeat_c_index_b": float(repeat_pairs["repeat_c_index_b"].mean()),` |
| 353 | mean | `"mean_repeat_c_index_a": float(residual["repeat_mean_c_index"].mean()),` |
| 499 | zscore | `def standardize_expression(expression, genes):` |
| 503 | mean | `means = x.mean(axis=0)` |
| 504 | std | `stds = x.std(axis=0).replace(0, np.nan)` |
| 512 | std | `if frame.shape[0] < 5 or frame["a"].std() == 0 or frame["b"].std() == 0:` |
| 514 | correlation | `return float(frame["a"].corr(frame["b"]))` |
| 588 | zscore | `z = standardize_expression(expression, genes)` |
| 623 | loadings | `for canine_gene, loading in risk_loadings.items():` |
| 635 | loadings | `"risk_oriented_loading": float(loading),` |
| 649 | loadings | `"risk_oriented_loading": float(loading),` |
| 781 | weights | `def deduplicate_human_weights(weights, mapping_flag):` |
| 782 | weights | `use = weights[weights[mapping_flag]].copy()` |
| 800 | weights | `def write_gmt(weights, path, suffix):` |
| 802 | weights | `for module_label, part in weights.groupby("module_label"):` |
| 823 | loadings, direction | `"score_formula": "arithmetic mean of signed gene z-scores using the sign of each frozen risk-oriented canine PCA loading",` |
| 834 | weights, loadings | `"score_formula": "sum of gene z-scores multiplied by frozen risk-oriented canine PCA loadings normalized by absolute weight sum",` |
| 858 | coverage | `"gene_mapping": "strict module genes plus strict proliferation genes after removing overlap",` |
| 860 | zscore | `"score_formula": "standardized residual from linear regression of M40 score on disjoint proliferation score",` |
| 885 | weights, direction | `"No human outcome may be used to change module membership, gene weights, score direction, or tier assignment.",` |
| 902 | loadings, direction | `"Primary human score: strict one-to-one signed mean z-score using frozen canine loading signs.",` |
| 903 | weights | `"Secondary zero-shot score: frozen canine PCA-weighted score.",` |
| 924 | weights | `print("  Export strict and broad human ortholog gene sets and frozen canine PCA weights.")` |
| 974 | weights | `"No strict one-to-one human ortholog weights were created. "` |
| 1103 | coverage | `print("Repeated cross-validation splits overlap; contrast p-values are descriptive stability diagnostics, not external validation.")` |
| 1107 | weights, direction | `print("No human outcome may be used to revise module membership, weights, score direction, or validation tier.")` |

### Referenced result-table schemas

#### `results/tables/GSE238110_RNA_full_cohort_module_associations.csv`

- Rows: 26
- SHA-256: `87e205e3c786d1028ee18922af5dd313c4cbe8f2565914b2e48ace1eb28f656f`
- Columns: `endpoint`, `module_label`, `n_genes`, `n`, `events`, `coef`, `hr`, `p`, `c_index`, `module_pc1_explained_variance`, `genes`, `gene_symbols`

#### `results/tables/GSE238110_RNA_full_cohort_transferable_module_priority.csv`

- Rows: 13
- SHA-256: `2d8ef5c813510f3ab3f1d89a736c25690c92cf1dec56ce443e4ed2bac77ded9f`
- Columns: `analysis`, `module_label`, `endpoint`, `fold`, `n_module_genes`, `n_strict_symbol_concordant`, `n_strict_one_to_one`, `n_broad_transferable`, `n_manual_review`, `fraction_strict_symbol_concordant`, `fraction_strict_one_to_one`, `fraction_broad_transferable`, `n_high_or_medium_rna_evidence`, `strict_human_symbols`, `broad_human_symbols`, `module_transfer_qc_tier`, `dfi_full_module_p`, `dfi_full_module_c_index`, `dfi_full_module_gene_symbols`, `os_full_module_p`, `os_full_module_c_index`, `os_full_module_gene_symbols`, `transfer_priority_score`

#### `results/tables/GSE238110_RNA_master_candidate_evidence_table_with_ortholog_qc.csv`

- Rows: 5013
- SHA-256: `26673d73d3e2b3769ae0eaf776d5f9fe239c594bf8c0f71e01b6ff50d1017586`
- Columns: `gene`, `gene_symbol_clean`, `univariate_candidate_rank`, `univariate_candidate_source`, `in_univariate_candidate_set`, `dfi_univ_rank`, `dfi_univ_coef`, `dfi_univ_hr_per_sd`, `dfi_univ_p`, `dfi_univ_q`, `dfi_univ_c_index`, `os_univ_rank`, `os_univ_coef`, `os_univ_hr_per_sd`, `os_univ_p`, `os_univ_q`, `os_univ_c_index`, `dfi_full_conditional_selected`, `dfi_full_conditional_rank`, `dfi_full_conditional_coef`, `dfi_full_conditional_hr_per_sd`, `dfi_full_conditional_p`, `dfi_full_conditional_c_index`, `os_full_conditional_selected`, `os_full_conditional_rank`, `os_full_conditional_coef`, `os_full_conditional_hr_per_sd`, `os_full_conditional_p`, `os_full_conditional_c_index`, `dfi_true_mb_n_configs`, `dfi_true_mb_algorithms`, `dfi_true_mb_alphas`, `dfi_true_mb_mean_rank`, `os_true_mb_n_configs`, `os_true_mb_algorithms`, `os_true_mb_alphas`, `os_true_mb_mean_rank`, `nested_dfi_conditional_selected_folds`, `nested_dfi_conditional_selection_frequency`, `nested_dfi_conditional_mean_rank`, `nested_dfi_conditional_plus_clinical_selected_folds`, `nested_dfi_conditional_plus_clinical_selection_frequency`, `nested_dfi_conditional_plus_clinical_mean_rank`, `nested_dfi_elasticnet_selected_folds`, `nested_dfi_elasticnet_selection_frequency`, `nested_dfi_elasticnet_mean_rank`, `nested_dfi_gsmb_selected_folds`, `nested_dfi_gsmb_selection_frequency`, `nested_dfi_gsmb_mean_rank`, `nested_dfi_iamb_selected_folds`, `nested_dfi_iamb_selection_frequency`, `nested_dfi_iamb_mean_rank`, `nested_dfi_univtop10_selected_folds`, `nested_dfi_univtop10_selection_frequency`, `nested_dfi_univtop10_mean_rank`, `nested_os_conditional_selected_folds`, `nested_os_conditional_selection_frequency`, `nested_os_conditional_mean_rank`, `nested_os_conditional_plus_clinical_selected_folds`, `nested_os_conditional_plus_clinical_selection_frequency`, `nested_os_conditional_plus_clinical_mean_rank`, `nested_os_elasticnet_selected_folds`, `nested_os_elasticnet_selection_frequency`, `nested_os_elasticnet_mean_rank`, `nested_os_gsmb_selected_folds`, `nested_os_gsmb_selection_frequency`, `nested_os_gsmb_mean_rank`, `nested_os_iamb_selected_folds`, `nested_os_iamb_selection_frequency`, `nested_os_iamb_mean_rank`, `nested_os_univtop10_selected_folds`, `nested_os_univtop10_selection_frequency`, `nested_os_univtop10_mean_rank`, `nested_any_max_selection_frequency`, `nested_any_total_selected_folds`, `dfi_os_both_univariate_available`, `dfi_os_same_direction`, `dfi_os_both_nominal_p05`, `dfi_os_both_fdr_q10`, `combined_univ_rank_score`, `rna_evidence_priority_score`, `rna_evidence_tier`, `gene_symbol_clean_upper`, `external_gene_name_upper`, `dog_ensembl_gene_id`, `dog_external_gene_name`, `dog_gene_biotype`, `human_ensembl_gene_id`, `human_gene_symbol`, `dog_human_orthology_type`, `dog_human_orthology_confidence`, `dog_human_orthology_confidence_numeric`, `dog_human_perc_id`, `human_dog_perc_id`, `has_human_homolog`, `is_one_to_one_ortholog`, `ortholog_mapping_status`, `dog_symbol_upper`, `human_symbol_upper`, `has_human_gene_symbol`, `human_symbol_same_as_dog_symbol`, `dog_symbol_problematic`, `human_symbol_problematic`, `ortholog_confidence_high`, `broad_transferable_ortholog`, `strict_transferable_ortholog`, `strict_symbol_concordant_transferable`, `needs_manual_ortholog_review`, `ortholog_qc_status`, `primary_human_validation_gene`, `sensitivity_human_validation_gene`

#### `results/tables/GSE238110_RNA_module_gene_membership.csv`

- Rows: 3227
- SHA-256: `62467463f491f25a9151eda96c65e115b206e0f3f62bc92d33b18f9939273888`
- Columns: `analysis`, `module_label`, `gene`, `gene_symbol_clean`, `endpoint`, `fold`

#### `results/tables/GSE238110_frozen_canine_transfer_program_manifest.csv`

- Rows: 13
- SHA-256: `fcd34b985ea99986f39dfc43ce8e4018db2fde9ecc8a225bb40ac93006d5372c`
- Columns: `module_label`, `validation_tier`, `multiplicity_family`, `provisional_program_label`, `program_label_requires_enrichment_confirmation`, `canine_primary_endpoint`, `canine_secondary_endpoint`, `positive_score_interpretation`, `n_canine_genes_used_for_pca`, `canine_pc1_explained_variance`, `canine_pca_orientation_correlation`, `dfi_full_cohort_coef`, `os_full_cohort_coef`, `risk_orientation_multiplier`, `n_strict_human_genes`, `n_broad_human_genes`, `strict_transfer_eligible`, `module_transfer_qc_tier`, `transfer_priority_score`, `fraction_strict_symbol_concordant`, `fraction_broad_transferable`, `raw_module_proliferation_correlation`, `orthogonal_variance_fraction_1_minus_r2`, `n_overlap_symbols_with_proliferation`, `primary_human_score`, `secondary_human_score`, `sensitivity_human_scores`, `frozen_after_canine_script`, `dfi_module_only_mean_c_index`, `dfi_module_only_std_c_index`, `dfi_module_only_fraction_above_0_50`, `dfi_module_plus_disjoint_proliferation_mean_c_index`, `dfi_module_plus_disjoint_proliferation_std_c_index`, `dfi_module_plus_disjoint_proliferation_fraction_above_0_50`, `dfi_residual_to_disjoint_proliferation_mean_c_index`, `dfi_residual_to_disjoint_proliferation_std_c_index`, `dfi_residual_to_disjoint_proliferation_fraction_above_0_50`, `dfi_residual_to_disjoint_proliferation_and_weight_mean_c_index`, `dfi_residual_to_disjoint_proliferation_and_weight_std_c_index`, `dfi_residual_to_disjoint_proliferation_and_weight_fraction_above_0_50`, `os_module_only_mean_c_index`, `os_module_only_std_c_index`, `os_module_only_fraction_above_0_50`, `os_module_plus_disjoint_proliferation_mean_c_index`, `os_module_plus_disjoint_proliferation_std_c_index`, `os_module_plus_disjoint_proliferation_fraction_above_0_50`, `os_residual_to_disjoint_proliferation_mean_c_index`, `os_residual_to_disjoint_proliferation_std_c_index`, `os_residual_to_disjoint_proliferation_fraction_above_0_50`, `os_residual_to_disjoint_proliferation_and_weight_mean_c_index`, `os_residual_to_disjoint_proliferation_and_weight_std_c_index`, `os_residual_to_disjoint_proliferation_and_weight_fraction_above_0_50`, `script20_dfi_recommended_role`, `script20_os_recommended_role`, `manual_freeze_reason`

#### `results/tables/GSE238110_frozen_transfer_gene_weights_broad.csv`

- Rows: 389
- SHA-256: `ee5f661677b173507e04856f848b6757ede2750dd69b9bb0d88a9f6a2a77bd4c`
- Columns: `module_label`, `canine_gene`, `canine_gene_symbol`, `human_gene_symbol`, `ortholog_qc_status`, `raw_pca_loading`, `risk_oriented_loading`, `absolute_risk_loading`, `is_strict_mapping`, `is_broad_mapping`, `human_symbol_duplicate_count`, `normalized_abs_sum_weight`

#### `results/tables/GSE238110_frozen_transfer_gene_weights_strict.csv`

- Rows: 321
- SHA-256: `4f065aa4c4edf117a0c74015840d2b4b2347929f172cd517e1818ba0f6163b91`
- Columns: `module_label`, `canine_gene`, `canine_gene_symbol`, `human_gene_symbol`, `ortholog_qc_status`, `raw_pca_loading`, `risk_oriented_loading`, `absolute_risk_loading`, `is_strict_mapping`, `is_broad_mapping`, `human_symbol_duplicate_count`, `normalized_abs_sum_weight`

#### `results/tables/GSE238110_frozen_transfer_scoring_specification.csv`

- Rows: 40
- SHA-256: `562aab92caa949691f43b0285cf2585c98008e78c62d2c15477a6d156d7f47a7`
- Columns: `module_label`, `validation_tier`, `score_name`, `analysis_role`, `gene_mapping`, `within_cohort_preprocessing`, `score_formula`, `minimum_gene_rule`, `outcome_use`

#### `results/tables/GSE238110_full_cohort_module_scores_for_proliferation_adjustment.csv`

- Rows: 186
- SHA-256: `9c8bba4a7ad3643170d8bedb07714c7bf8b1c17e0ef2366c7b05e7d381197cde`
- Columns: ``, `module_M11_score`, `module_M17_score`, `module_M20_score`, `module_M23_score`, `module_M24_score`, `module_M25_score`, `module_M27_score`, `module_M28_score`, `module_M30_score`, `module_M34_score`, `module_M38_score`, `module_M39_score`, `module_M40_score`

#### `results/tables/GSE238110_module_proliferation_overlap_audit.csv`

- Rows: 13
- SHA-256: `fb7cce6c16e4f5f4a0adb99227c9c9efd99ea4d7464faccec0e54d2188b02619`
- Columns: `module_label`, `n_module_genes_used`, `n_proliferation_genes_total`, `n_overlap_symbols`, `n_disjoint_proliferation_genes`, `fraction_module_symbols_in_proliferation`, `fraction_proliferation_symbols_in_module`, `raw_module_proliferation_correlation`, `orthogonal_variance_fraction_1_minus_r2`, `original_residual_sd_before_standardization`, `original_residual_post_correlation`, `overlap_symbols`, `disjoint_proliferation_genes`, `module_transfer_qc_tier`, `transfer_priority_score`, `n_module_genes`, `fraction_strict_symbol_concordant`, `fraction_broad_transferable`, `strict_human_symbols`

#### `results/tables/GSE238110_proliferation_independence_decision_table.csv`

- Rows: 26
- SHA-256: `576d5ed2921e57ee3fa7431b6e6be6cb647552a6a94c4f307b88586eb6fc8328`
- Columns: `module_label`, `n_module_genes_used`, `n_proliferation_genes_total`, `n_overlap_symbols`, `n_disjoint_proliferation_genes`, `fraction_module_symbols_in_proliferation`, `fraction_proliferation_symbols_in_module`, `raw_module_proliferation_correlation`, `orthogonal_variance_fraction_1_minus_r2`, `original_residual_sd_before_standardization`, `original_residual_post_correlation`, `overlap_symbols`, `disjoint_proliferation_genes`, `module_transfer_qc_tier`, `transfer_priority_score`, `n_module_genes`, `fraction_strict_symbol_concordant`, `fraction_broad_transferable`, `strict_human_symbols`, `endpoint`, `full_disjoint_residual_p`, `full_disjoint_residual_q`, `full_disjoint_residual_c_index`, `module_disjoint_proliferation_correlation`, `residual_disjoint_sd_before_standardization`, `n_valid_folds`, `cv_disjoint_residual_mean_c_index`, `cv_disjoint_residual_std_c_index`, `cv_disjoint_residual_median_c_index`, `cv_disjoint_residual_fraction_above_0_50`, `cv_disjoint_residual_fraction_above_0_55`, `mean_train_module_original_proliferation_correlation`, `mean_train_module_disjoint_proliferation_correlation`, `passes_full_cohort_disjoint_residual_fdr10`, `passes_cv_disjoint_residual_mean_055`, `passes_cv_majority_above_chance`, `is_extreme_proliferation_overlap`, `recommended_role`

#### `results/tables/GSE238110_repeated_cv_program_model_contrast_summary.csv`

- Rows: 112
- SHA-256: `2aec5db69a38b9dab0e8769334c977a86f0f505fcd16e6402b2bbe9f1e9832de`
- Columns: `endpoint`, `module_label`, `contrast`, `model_a`, `model_b`, `mean_repeat_c_index_a`, `mean_repeat_c_index_b`, `n_repeats`, `n_fold_pairs`, `mean_delta_c_index`, `median_delta_c_index`, `std_delta_c_index`, `bootstrap_repeat_mean_ci_low`, `bootstrap_repeat_mean_ci_high`, `fraction_repeat_deltas_positive`, `wilcoxon_repeat_p_two_sided`, `wilcoxon_repeat_p_greater`, `sign_repeat_p_two_sided`, `sign_repeat_p_greater`, `corrected_t`, `corrected_t_df`, `corrected_t_p_two_sided`, `corrected_standard_error`, `correction_factor`, `wilcoxon_repeat_p_two_sided_bh_q`, `wilcoxon_repeat_p_greater_bh_q`, `corrected_t_p_two_sided_bh_q`, `descriptive_support`

#### `results/tables/GSE238110_repeated_cv_program_model_contrasts.csv`

- Rows: 182
- SHA-256: `c4dd9f810ae1df0f3662754cb05349023bfb71a9ad1a5f3497401d8c0bbd5382`
- Columns: `endpoint`, `module_label`, `contrast`, `model_a`, `model_b`, `mean_repeat_c_index_a`, `mean_repeat_c_index_b`, `n_repeats`, `n_fold_pairs`, `mean_delta_c_index`, `median_delta_c_index`, `std_delta_c_index`, `bootstrap_repeat_mean_ci_low`, `bootstrap_repeat_mean_ci_high`, `fraction_repeat_deltas_positive`, `wilcoxon_repeat_p_two_sided`, `wilcoxon_repeat_p_greater`, `sign_repeat_p_two_sided`, `sign_repeat_p_greater`, `corrected_t`, `corrected_t_df`, `corrected_t_p_two_sided`, `corrected_standard_error`, `correction_factor`, `wilcoxon_repeat_p_two_sided_bh_q`, `wilcoxon_repeat_p_greater_bh_q`, `corrected_t_p_two_sided_bh_q`, `descriptive_support`

#### `results/tables/GSE238110_repeated_cv_program_scores_by_repeat.csv`

- Rows: 4160
- SHA-256: `90e1ba5ecdef36d1fc813ead372b6fad70bfe6e7f57dc7578e353939b3a09c75`
- Columns: `endpoint`, `module_label`, `model`, `repeat`, `n_valid_folds`, `repeat_mean_c_index`, `repeat_median_c_index`, `repeat_min_c_index`, `repeat_max_c_index`

#### `results/tables/GSE238110_repeated_cv_proliferation_sensitivity_fold_results.csv`

- Rows: 20800
- SHA-256: `f079774988dbaa2bdadf9b408f90b0867d4b413e60d96fc48a82e774520cfee3`
- Columns: `endpoint`, `repeat`, `fold`, `module_label`, `model`, `c_index`, `error`, `n_train`, `n_test`, `train_events`, `test_events`, `n_module_genes_used`, `n_original_proliferation_genes_used`, `n_disjoint_proliferation_genes_used`, `module_pc1_explained_variance`, `original_proliferation_pc1_explained_variance`, `disjoint_proliferation_pc1_explained_variance`, `train_module_original_proliferation_correlation`, `train_module_disjoint_proliferation_correlation`, `residual_sd_before_standardization`

#### `results/tables/GSE238110_repeated_cv_proliferation_sensitivity_summary.csv`

- Rows: 208
- SHA-256: `9e49c7841bc96ac50e3ec76bc224754672a9d190c621bd6d92849c8075cacf18`
- Columns: `endpoint`, `module_label`, `model`, `n_valid_folds`, `mean_c_index`, `std_c_index`, `median_c_index`, `q25_c_index`, `q75_c_index`, `min_c_index`, `max_c_index`, `fraction_above_0_50`, `fraction_above_0_55`, `mean_train_module_original_proliferation_correlation`, `mean_train_module_disjoint_proliferation_correlation`, `mean_n_module_genes_used`, `mean_n_disjoint_proliferation_genes_used`, `mean_c_index__disjoint_proliferation_only`, `mean_c_index__module_only`, `mean_c_index__module_plus_disjoint_proliferation`, `mean_c_index__module_plus_original_proliferation`, `mean_c_index__original_proliferation_only`, `mean_c_index__residual_to_disjoint_proliferation`, `mean_c_index__residual_to_disjoint_proliferation_and_weight`, `mean_c_index__residual_to_original_proliferation`, `delta_disjoint_residual_vs_original_residual`, `delta_joint_module_vs_disjoint_proliferation_only`

---

## `scripts/22_prepare_human_osteosarcoma_cohorts.py`

- SHA-256: `2b4b9734204b14dd6494ba168d08f1b52dc02d8d7ffd56637d75b8083fafffdd`

### Relevant functions

#### `verify_frozen_inputs`

- Lines: 146-195
- Signals: `weights`

```python
def verify_frozen_inputs() -> dict[str, Any]:
    if not FREEZE_JSON_FILE.exists():
        raise FileNotFoundError(
            f"Frozen-program manifest is missing: {FREEZE_JSON_FILE}. "
            "Run the corrected script 21 first."
        )

    freeze = json.loads(FREEZE_JSON_FILE.read_text(encoding="utf-8"))
    files = freeze.get("files", {})
    required = [
        FROZEN_MANIFEST_FILE,
        STRICT_WEIGHTS_FILE,
        BROAD_WEIGHTS_FILE,
        SCORING_SPEC_FILE,
    ]

    print("")
    print("Frozen input integrity check:")
    for path in required:
        if not path.exists():
            raise FileNotFoundError(f"Frozen input is missing: {path}")
        expected = files.get(path.name, {}).get("sha256")
        observed = sha256_file(path)
        if expected and expected != observed:
            raise RuntimeError(
                f"Frozen input hash mismatch for {path.name}. "
                "Do not continue after modifying a frozen file."
            )
        status = "verified" if expected else "hash_not_recorded_but_file_present"
        print(f"  {path.name}: {status}")

    strict = pd.read_csv(STRICT_WEIGHTS_FILE)
    if strict.empty:
        raise RuntimeError("The strict frozen weight table is empty.")

    primary_counts = (
        strict[strict["module_label"].isin(PRIMARY_MODULES)]
        .groupby("module_label")["human_gene_symbol"]
        .nunique()
    )
    missing_primary = [module for module in PRIMARY_MODULES if primary_counts.get(module, 0) < 3]
    if missing_primary:
        raise RuntimeError(
            "Primary frozen modules have insufficient strict genes: "
            + ", ".join(missing_primary)
        )

    print("  Primary strict gene counts:")
    print(primary_counts.to_string())
    return freeze
```

#### `standardize_target_clinical`

- Lines: 517-610
- Signals: `zscore`

```python
def standardize_target_clinical(cases: list[dict[str, Any]]) -> pd.DataFrame:
    rows = []
    for case in cases:
        demographic = case.get("demographic") or {}
        if isinstance(demographic, list):
            demographic = demographic[0] if demographic else {}
        diagnoses = case.get("diagnoses") or []
        follow_ups = case.get("follow_ups") or []

        vital_values = [demographic.get("vital_status")]
        vital_values.extend(item.get("vital_status") for item in follow_ups if isinstance(item, dict))
        vital_status = str(first_nonempty(vital_values, default="Unknown"))
        os_event = 1.0 if vital_status.lower() == "dead" else 0.0 if vital_status.lower() == "alive" else np.nan

        days_to_death = numeric_or_nan(demographic.get("days_to_death"))
        follow_up_times = [
            numeric_or_nan(item.get("days_to_follow_up"))
            for item in follow_ups
            if isinstance(item, dict)
        ]
        diagnosis_follow_up_times = []
        diagnosis_age_values = []
        diagnosis_metastasis_values = []
        for diagnosis in diagnoses:
            if not isinstance(diagnosis, dict):
                continue
            diagnosis_follow_up_times.extend(
                [
                    numeric_or_nan(diagnosis.get("days_to_last_follow_up")),
                    numeric_or_nan(diagnosis.get("days_to_last_known_disease_status")),
                ]
            )
            diagnosis_age_values.append(numeric_or_nan(diagnosis.get("age_at_diagnosis")))
            diagnosis_metastasis_values.extend(
                [
                    diagnosis.get("metastasis_at_diagnosis"),
                    diagnosis.get("ajcc_clinical_m"),
                    diagnosis.get("ajcc_pathologic_m"),
                ]
            )

        candidate_follow_up = [
            value
            for value in follow_up_times + diagnosis_follow_up_times
            if np.isfinite(value)
        ]
        days_to_last_follow_up = max(candidate_follow_up) if candidate_follow_up else np.nan
        if os_event == 1 and np.isfinite(days_to_death):
            os_time = days_to_death
        else:
            os_time = days_to_last_follow_up

        age_at_diagnosis_days = next(
            (value for value in diagnosis_age_values if np.isfinite(value)),
            np.nan,
        )
        age_at_diagnosis_years = age_at_diagnosis_days / 365.25 if np.isfinite(age_at_diagnosis_days) else np.nan

        metastasis_text = ";".join(
            sorted(
                {
                    str(value).strip()
                    for value in diagnosis_metastasis_values
                    if value is not None and str(value).strip()
                }
            )
        )

        rows.append(
            {
                "case_submitter_id": case.get("submitter_id", ""),
                "case_id": case.get("case_id", ""),
                "primary_site": case.get("primary_site", ""),
                "disease_type": case.get("disease_type", ""),
                "sex": demographic.get("gender", ""),
                "race": demographic.get("race", ""),
                "ethnicity": demographic.get("ethnicity", ""),
                "vital_status": vital_status,
                "os_time_days": os_time,
                "os_event": os_event,
                "days_to_death": days_to_death,
                "days_to_last_follow_up": days_to_last_follow_up,
                "age_at_diagnosis_days": age_at_diagnosis_days,
                "age_at_diagnosis_years": age_at_diagnosis_years,
                "metastasis_fields_raw": metastasis_text,
                "endpoint_preparation_note": "OS prepared from GDC vital status and death/follow-up fields; no outcome model fitted in script 22",
            }
        )

    clinical = pd.DataFrame(rows)
    clinical = clinical[clinical["case_submitter_id"].astype(str).str.len() > 0].copy()
    clinical = clinical.drop_duplicates("case_submitter_id", keep="first")
    clinical = clinical.set_index("case_submitter_id")
    return clinical
```

#### `prepare_gse21257`

- Lines: 734-818
- Signals: `loadings`

```python
def prepare_gse21257() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    GEOparse = import_geoparse()
    print("")
    print("Downloading/loading GSE21257 from NCBI GEO:")
    gse = GEOparse.get_GEO(geo=GSE_ACCESSION, destdir=str(GSE21257_RAW_DIR), silent=False)

    expression_probe_by_sample = gse.pivot_samples("VALUE")
    if expression_probe_by_sample.empty:
        raise RuntimeError("GSE21257 expression matrix is empty.")
    expression_probe_by_sample.index = expression_probe_by_sample.index.astype(str)

    phenotype = gse.phenotype_data.copy()
    phenotype.index = phenotype.index.astype(str)
    phenotype.to_csv(OUTPUT_GSE_PHENOTYPE_RAW)

    if not gse.gpls:
        raise RuntimeError("GSE21257 platform annotation was not loaded by GEOparse.")
    gpl_name = sorted(gse.gpls.keys())[0]
    platform = gse.gpls[gpl_name].table.copy()
    probe_col, symbol_col = detect_platform_columns(platform)

    annotation = platform[[probe_col, symbol_col]].copy()
    annotation[probe_col] = annotation[probe_col].astype(str)
    annotation["gene_symbol"] = annotation[symbol_col].map(normalize_unambiguous_gene_symbol)
    annotation = annotation[annotation["gene_symbol"].ne("")].copy()
    annotation = annotation.drop_duplicates(probe_col, keep="first")

    common_probes = expression_probe_by_sample.index.intersection(annotation[probe_col])
    if len(common_probes) == 0:
        raise RuntimeError("No GSE21257 probes matched the platform annotation.")

    expr = expression_probe_by_sample.loc[common_probes].apply(pd.to_numeric, errors="coerce")
    ann = annotation.set_index(probe_col).loc[common_probes]
    probe_variance = expr.var(axis=1)
    probe_map = pd.DataFrame(
        {
            "probe_id": common_probes,
            "gene_symbol": ann["gene_symbol"].values,
            "probe_variance": probe_variance.loc[common_probes].values,
        }
    )
    probe_map = probe_map.sort_values(
        ["gene_symbol", "probe_variance", "probe_id"],
        ascending=[True, False, True],
    )
    selected_map = probe_map.drop_duplicates("gene_symbol", keep="first")
    selected_map.to_csv(OUTPUT_GSE_PROBE_MAP, index=False)

    selected_expr = expr.loc[selected_map["probe_id"].tolist()].copy()
    selected_expr.index = selected_map.set_index("probe_id").loc[selected_expr.index, "gene_symbol"]
    expression = selected_expr.T
    expression.index.name = "geo_sample_id"
    expression = expression.loc[:, ~expression.columns.duplicated()].copy()
    expression = expression.loc[:, expression.var(axis=0) > 0]

    phenotype_text = phenotype.astype(str).agg(" | ".join, axis=1)
    clinical = phenotype.copy()
    clinical["metadata_text_combined"] = phenotype_text
    clinical["metastasis_within_5y"] = phenotype_text.map(classify_gse21257_metastasis)
    survival_parsed = phenotype_text.map(parse_gse21257_survival)
    clinical["os_time_months"] = survival_parsed.map(lambda value: value[0])
    clinical["os_time_days_approx"] = clinical["os_time_months"] * 30.4375
    clinical["os_event"] = survival_parsed.map(lambda value: value[1])
    clinical["age_months"] = phenotype_text.map(parse_gse21257_age_months)
    clinical["age_years"] = clinical["age_months"] / 12.0
    clinical["group_parsed"] = phenotype_text.map(lambda value: extract_gse_characteristic(value, "group"))
    clinical["status_parsed"] = phenotype_text.map(lambda value: extract_gse_characteristic(value, "status"))
    clinical["metastasis_endpoint_note"] = (
        "Pre-chemotherapy biopsy group parsed from GEO sample metadata; "
        "no classifier or association test fitted in script 22"
    )
    clinical["survival_endpoint_note"] = (
        "OS time/event parsed from GEO status text for sensitivity analysis; "
        "no survival model fitted in script 22"
    )
    clinical.index.name = "geo_sample_id"

    common_samples = expression.index.intersection(clinical.index)
    expression = expression.loc[common_samples].copy()
    clinical = clinical.loc[common_samples].copy()

    counts = clinical["metastasis_within_5y"].value_counts(dropna=False)
    print("GSE21257 metastasis label counts:")
    print(counts.to_string())
    return expression, clinical, selected_map
```

#### `zscore_columns`

- Lines: 821-827
- Signals: `mean, std`

```python
def zscore_columns(expression: pd.DataFrame) -> pd.DataFrame:
    x = expression.apply(pd.to_numeric, errors="coerce")
    x = x.replace([np.inf, -np.inf], np.nan)
    x = x.fillna(x.median(axis=0))
    stds = x.std(axis=0).replace(0, np.nan)
    z = (x - x.mean(axis=0)) / stds
    return z.loc[:, z.notna().all(axis=0)]
```

#### `zscore_series`

- Lines: 830-835
- Signals: `mean, std`

```python
def zscore_series(values: pd.Series) -> pd.Series:
    values = pd.to_numeric(values, errors="coerce")
    std = values.std()
    if not np.isfinite(std) or std == 0:
        return pd.Series(np.nan, index=values.index)
    return (values - values.mean()) / std
```

#### `safe_corr`

- Lines: 838-842
- Signals: `std, correlation`

```python
def safe_corr(a: pd.Series, b: pd.Series) -> float:
    frame = pd.concat([a.rename("a"), b.rename("b")], axis=1).dropna()
    if frame.shape[0] < 5 or frame["a"].std() == 0 or frame["b"].std() == 0:
        return np.nan
    return float(frame["a"].corr(frame["b"]))
```

#### `human_cohort_pc1`

- Lines: 845-854
- Signals: `correlation`

```python
def human_cohort_pc1(z: pd.DataFrame, reference: pd.Series | None = None) -> pd.Series:
    if z.shape[1] < 2:
        return pd.Series(np.nan, index=z.index)
    pca = PCA(n_components=1, random_state=RANDOM_SEED)
    score = pd.Series(pca.fit_transform(z).ravel(), index=z.index)
    if reference is not None:
        corr = safe_corr(score, reference)
        if np.isfinite(corr) and corr < 0:
            score = -score
    return zscore_series(score)
```

#### `compute_module_scores`

- Lines: 857-930
- Signals: `mean, weights, direction, coverage`

```python
def compute_module_scores(
    expression: pd.DataFrame,
    strict_weights: pd.DataFrame,
    broad_weights: pd.DataFrame,
    manifest: pd.DataFrame,
    cohort_name: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    expression = expression.copy()
    expression.columns = expression.columns.astype(str).str.upper()
    expression = expression.loc[:, ~expression.columns.duplicated()].copy()

    scores = pd.DataFrame(index=expression.index)
    coverage_rows = []

    for mapping_name, weights in [("strict", strict_weights), ("broad", broad_weights)]:
        for module_label, part in weights.groupby("module_label"):
            part = part.copy()
            part["human_gene_symbol"] = part["human_gene_symbol"].astype(str).str.upper()
            part = part.drop_duplicates("human_gene_symbol", keep="first")
            requested = part["human_gene_symbol"].tolist()
            available = [gene for gene in requested if gene in expression.columns]
            n_requested = len(requested)
            n_available = len(available)
            fraction = n_available / n_requested if n_requested else 0.0
            passed = n_available >= MIN_SCORE_GENES and fraction >= MIN_SCORE_FRACTION

            coverage_rows.append(
                {
                    "cohort": cohort_name,
                    "module_label": module_label,
                    "mapping": mapping_name,
                    "n_frozen_genes": n_requested,
                    "n_available_genes": n_available,
                    "coverage_fraction": fraction,
                    "minimum_rule_passed": passed,
                    "available_genes": ";".join(available),
                    "missing_genes": ";".join([gene for gene in requested if gene not in available]),
                }
            )
            if not passed:
                continue

            z = zscore_columns(expression[available])
            available = list(z.columns)
            if len(available) < MIN_SCORE_GENES:
                continue
            weight_indexed = part.set_index("human_gene_symbol").loc[available]
            signs = np.sign(
                pd.to_numeric(weight_indexed["risk_oriented_loading"], errors="coerce")
            ).replace(0, 1)
            signed_mean = z.mul(signs, axis=1).mean(axis=1)
            signed_mean = zscore_series(signed_mean)

            raw_weights = pd.to_numeric(
                weight_indexed["risk_oriented_loading"], errors="coerce"
            ).fillna(0.0)
            if raw_weights.abs().sum() > 0:
                normalized_weights = raw_weights / raw_weights.abs().sum()
                weighted = z.mul(normalized_weights, axis=1).sum(axis=1)
                weighted = zscore_series(weighted)
            else:
                weighted = pd.Series(np.nan, index=z.index)

            pc1 = human_cohort_pc1(z, reference=signed_mean)
            prefix = f"{module_label}__{mapping_name}"
            scores[f"{prefix}__signed_mean_z"] = signed_mean
            scores[f"{prefix}__canine_pca_weighted_z"] = weighted
            scores[f"{prefix}__human_pc1_z"] = pc1

    coverage = pd.DataFrame(coverage_rows)
    tier_map = manifest.set_index("module_label")["validation_tier"].to_dict()
    scores.insert(0, "cohort", cohort_name)
    coverage["validation_tier"] = coverage["module_label"].map(tier_map)
    return scores, coverage
```

#### `detect_ortholog_columns`

- Lines: 933-945
- Signals: ``

```python
def detect_ortholog_columns(table: pd.DataFrame) -> tuple[str, str, str]:
    dog_col = next(
        (column for column in ["gene", "gene_symbol_clean", "canine_gene_symbol"] if column in table.columns),
        None,
    )
    human_col = next(
        (column for column in ["human_gene_symbol", "human_symbol"] if column in table.columns),
        None,
    )
    status_col = "ortholog_qc_status" if "ortholog_qc_status" in table.columns else None
    if dog_col is None or human_col is None or status_col is None:
        raise ValueError("Ortholog QC table lacks required dog, human, or QC-status columns.")
    return dog_col, human_col, status_col
```

#### `build_strict_human_proliferation_mapping`

- Lines: 956-995
- Signals: ``

```python
def build_strict_human_proliferation_mapping(
    proliferation_table: pd.DataFrame,
    ortholog_qc: pd.DataFrame,
) -> pd.DataFrame:
    proliferation_gene_col = next(
        (
            column
            for column in ["gene", "expression_column", "gene_symbol", "canine_gene"]
            if column in proliferation_table.columns
        ),
        None,
    )
    if proliferation_gene_col is None:
        raise ValueError(
            f"Could not detect the proliferation gene column: {list(proliferation_table.columns)}"
        )

    dog_col, human_col, status_col = detect_ortholog_columns(ortholog_qc)
    mapping = ortholog_qc.copy()
    mapping["dog_symbol_key"] = mapping[dog_col].map(clean_canine_symbol)
    mapping[human_col] = mapping[human_col].astype(str).str.strip().str.upper()
    mapping = mapping[
        mapping[status_col].eq("strict_symbol_concordant_one_to_one")
        & mapping[human_col].ne("")
    ].copy()

    genes = proliferation_table[[proliferation_gene_col]].copy()
    genes["canine_proliferation_gene"] = genes[proliferation_gene_col].astype(str)
    genes["dog_symbol_key"] = genes[proliferation_gene_col].map(clean_canine_symbol)
    merged = genes.merge(
        mapping[["dog_symbol_key", human_col, status_col]],
        on="dog_symbol_key",
        how="left",
    )
    merged = merged.rename(columns={human_col: "human_gene_symbol"})
    merged = merged[merged["human_gene_symbol"].notna()].copy()
    merged["human_gene_symbol"] = merged["human_gene_symbol"].astype(str).str.upper()
    merged = merged.drop_duplicates("human_gene_symbol", keep="first")
    merged.to_csv(OUTPUT_PROLIFERATION_MAPPING, index=False)
    return merged
```

#### `compute_proliferation_scores`

- Lines: 998-1074
- Signals: `mean`

```python
def compute_proliferation_scores(
    expression: pd.DataFrame,
    proliferation_mapping: pd.DataFrame,
    strict_weights: pd.DataFrame,
    cohort_name: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    expression = expression.copy()
    expression.columns = expression.columns.astype(str).str.upper()
    expression = expression.loc[:, ~expression.columns.duplicated()].copy()

    frozen_genes = proliferation_mapping["human_gene_symbol"].dropna().astype(str).str.upper().drop_duplicates().tolist()
    available = [gene for gene in frozen_genes if gene in expression.columns]
    fraction = len(available) / len(frozen_genes) if frozen_genes else 0.0
    passed = len(available) >= MIN_PROLIFERATION_GENES and fraction >= MIN_PROLIFERATION_FRACTION

    coverage_rows = [
        {
            "cohort": cohort_name,
            "score": "strict_human_meta_proliferation_pc1",
            "n_frozen_genes": len(frozen_genes),
            "n_available_genes": len(available),
            "coverage_fraction": fraction,
            "minimum_rule_passed": passed,
            "available_genes": ";".join(available),
            "missing_genes": ";".join([gene for gene in frozen_genes if gene not in available]),
        }
    ]

    scores = pd.DataFrame(index=expression.index)
    if not passed:
        return scores, pd.DataFrame(coverage_rows)

    z = zscore_columns(expression[available])
    anchor_genes = [gene for gene in PROLIFERATION_ANCHOR_SYMBOLS if gene in z.columns]
    reference = z[anchor_genes].mean(axis=1) if len(anchor_genes) >= 3 else z.mean(axis=1)
    proliferation_pc1 = human_cohort_pc1(z, reference=reference)
    scores["strict_human_meta_proliferation_pc1_z"] = proliferation_pc1

    m40_genes = (
        strict_weights[strict_weights["module_label"].eq("M40")]["human_gene_symbol"]
        .dropna()
        .astype(str)
        .str.upper()
        .drop_duplicates()
        .tolist()
    )
    disjoint_frozen = [gene for gene in frozen_genes if gene not in set(m40_genes)]
    disjoint_available = [gene for gene in disjoint_frozen if gene in expression.columns]
    disjoint_fraction = len(disjoint_available) / len(disjoint_frozen) if disjoint_frozen else 0.0
    disjoint_passed = (
        len(disjoint_available) >= MIN_PROLIFERATION_GENES
        and disjoint_fraction >= MIN_PROLIFERATION_FRACTION
    )
    coverage_rows.append(
        {
            "cohort": cohort_name,
            "score": "M40_disjoint_strict_human_meta_proliferation_pc1",
            "n_frozen_genes": len(disjoint_frozen),
            "n_available_genes": len(disjoint_available),
            "coverage_fraction": disjoint_fraction,
            "minimum_rule_passed": disjoint_passed,
            "available_genes": ";".join(disjoint_available),
            "missing_genes": ";".join([gene for gene in disjoint_frozen if gene not in disjoint_available]),
        }
    )
    if disjoint_passed:
        z_disjoint = zscore_columns(expression[disjoint_available])
        anchor_disjoint = [gene for gene in PROLIFERATION_ANCHOR_SYMBOLS if gene in z_disjoint.columns]
        reference_disjoint = (
            z_disjoint[anchor_disjoint].mean(axis=1)
            if len(anchor_disjoint) >= 3
            else z_disjoint.mean(axis=1)
        )
        disjoint_pc1 = human_cohort_pc1(z_disjoint, reference=reference_disjoint)
        scores["M40_disjoint_strict_human_meta_proliferation_pc1_z"] = disjoint_pc1

    return scores, pd.DataFrame(coverage_rows)
```

#### `residualize_outcome_blind`

- Lines: 1077-1085
- Signals: `std, matrix_product`

```python
def residualize_outcome_blind(score: pd.Series, covariate: pd.Series) -> pd.Series:
    frame = pd.concat([score.rename("score"), covariate.rename("covariate")], axis=1).dropna()
    residual = pd.Series(np.nan, index=score.index)
    if frame.shape[0] < 10 or frame["covariate"].std() == 0:
        return residual
    x = np.column_stack([np.ones(frame.shape[0]), frame["covariate"].values])
    beta, _, _, _ = np.linalg.lstsq(x, frame["score"].values, rcond=None)
    residual.loc[frame.index] = frame["score"].values - x @ beta
    return zscore_series(residual)
```

#### `add_m40_residual_scores`

- Lines: 1088-1106
- Signals: ``

```python
def add_m40_residual_scores(scores: pd.DataFrame) -> pd.DataFrame:
    out = scores.copy()
    proliferation_col = "M40_disjoint_strict_human_meta_proliferation_pc1_z"
    if proliferation_col not in out.columns:
        return out
    for module_score_col in [
        "M40__strict__signed_mean_z",
        "M40__strict__canine_pca_weighted_z",
    ]:
        if module_score_col not in out.columns:
            continue
        residual_col = module_score_col.replace(
            "_z", "__residual_to_disjoint_proliferation_z"
        )
        out[residual_col] = residualize_outcome_blind(
            out[module_score_col],
            out[proliferation_col],
        )
    return out
```

#### `merge_score_components`

- Lines: 1109-1119
- Signals: ``

```python
def merge_score_components(
    module_scores: pd.DataFrame,
    proliferation_scores: pd.DataFrame,
) -> pd.DataFrame:
    cohort = module_scores["cohort"] if "cohort" in module_scores.columns else None
    module_numeric = module_scores.drop(columns=["cohort"], errors="ignore")
    out = module_numeric.join(proliferation_scores, how="outer")
    out = add_m40_residual_scores(out)
    if cohort is not None:
        out.insert(0, "cohort", cohort.reindex(out.index))
    return out
```

#### `write_readme`

- Lines: 1122-1124
- Signals: `weights, direction`

```python
def write_readme() -> None:
    text = f"""Human osteosarcoma cohort preparation\n\nScript version: {SCRIPT_VERSION}\n\nThis script prepares TARGET-OS and GSE21257 without fitting any outcome model.\nFrozen canine module membership, score direction, PCA weights, and validation tiers are not changed.\n\nPrimary external scores:\n- strict one-to-one signed mean z-score for M34, M11, M24, and M40\n- TARGET-OS: overall-survival metadata prepared from public GDC clinical fields\n- GSE21257: metastasis-within-five-years label parsed from GEO metadata\n\nSecondary/sensitivity scores:\n- strict canine PCA-weighted score\n- broad mapped score\n- human-cohort PC1 oriented without outcomes\n- M40 residual to a disjoint strict human proliferation PC1\n\nNo survival or metastasis association is tested in this script.\n"""
    OUTPUT_README.write_text(text, encoding="utf-8")
```

#### `print_coverage_summary`

- Lines: 1193-1214
- Signals: `coverage`

```python
def print_coverage_summary(coverage: pd.DataFrame, cohort: str) -> None:
    print("")
    print("=" * 80)
    print(f"Frozen score coverage: {cohort}")
    print("=" * 80)
    display = coverage.copy()
    keep = [
        "module_label",
        "mapping",
        "validation_tier",
        "n_frozen_genes",
        "n_available_genes",
        "coverage_fraction",
        "minimum_rule_passed",
    ]
    keep = [column for column in keep if column in display.columns]
    if keep:
        print(
            display[keep]
            .sort_values([column for column in ["validation_tier", "module_label", "mapping"] if column in keep])
            .to_string(index=False)
        )
```

#### `main`

- Lines: 1217-1427
- Signals: `zscore, weights`

```python
def main() -> None:
    print("=" * 80)
    print("Prepare human osteosarcoma validation cohorts")
    print("=" * 80)
    print(f"Script version: {SCRIPT_VERSION}")
    print(f"Project root: {PROJECT_ROOT}")
    print(f"Human processed directory: {HUMAN_PROCESSED_DIR}")
    print("")
    print("Design:")
    print("  Verify frozen canine transfer assets and hashes.")
    print("  Acquire public TARGET-OS RNA-seq and clinical metadata from the GDC API.")
    print("  Acquire GSE21257 expression and phenotype metadata from NCBI GEO.")
    print("  Harmonize expression to human gene symbols using outcome-blind rules.")
    print("  Construct frozen module scores without fitting any outcome model.")
    print("")

    freeze = verify_frozen_inputs()
    manifest = read_required_csv(FROZEN_MANIFEST_FILE)
    strict_weights = read_required_csv(STRICT_WEIGHTS_FILE)
    broad_weights = read_required_csv(BROAD_WEIGHTS_FILE)
    ortholog_qc = read_required_csv(ORTHOLOG_QC_FILE)
    proliferation_table = read_required_csv(PROLIFERATION_GENE_FILE)

    strict_weights["human_gene_symbol"] = strict_weights["human_gene_symbol"].astype(str).str.upper()
    broad_weights["human_gene_symbol"] = broad_weights["human_gene_symbol"].astype(str).str.upper()

    proliferation_mapping = build_strict_human_proliferation_mapping(
        proliferation_table,
        ortholog_qc,
    )
    print(f"Strict human proliferation genes: {proliferation_mapping.shape[0]}")

    target_file_manifest = query_target_os_expression_files()
    print("")
    print("TARGET-OS GDC expression manifest:")
    print(f"  Cases/files selected: {target_file_manifest.shape[0]}")
    print(
        target_file_manifest["sample_types"]
        .replace("", "<missing>")
        .value_counts()
        .head(10)
        .to_string()
    )
    target_downloaded = download_target_expression_files(target_file_manifest)
    target_expression, target_sample_map = build_target_expression(target_downloaded)
    target_cases = query_target_cases()
    target_clinical = standardize_target_clinical(target_cases)

    target_common = target_expression.index.intersection(target_clinical.index)
    target_expression = target_expression.loc[target_common].copy()
    target_clinical = target_clinical.loc[target_common].copy()
    target_sample_map = target_sample_map[
        target_sample_map["case_submitter_id"].isin(target_common)
    ].copy()

    target_module_scores, target_coverage = compute_module_scores(
        target_expression,
        strict_weights,
        broad_weights,
        manifest,
        "TARGET_OS",
    )
    target_proliferation_scores, target_proliferation_coverage = compute_proliferation_scores(
        target_expression,
        proliferation_mapping,
        strict_weights,
        "TARGET_OS",
    )
    target_scores = merge_score_components(
        target_module_scores,
        target_proliferation_scores,
    )
    target_coverage_all = pd.concat(
        [target_coverage, target_proliferation_coverage],
        axis=0,
        ignore_index=True,
        sort=False,
    )

    target_expression.to_csv(OUTPUT_TARGET_EXPRESSION)
    target_clinical.to_csv(OUTPUT_TARGET_CLINICAL)
    target_sample_map.to_csv(OUTPUT_TARGET_SAMPLE_MAP, index=False)
    target_scores.to_csv(OUTPUT_TARGET_SCORES)
    target_coverage_all.to_csv(OUTPUT_TARGET_COVERAGE, index=False)

    gse_expression, gse_clinical, _ = prepare_gse21257()
    gse_module_scores, gse_coverage = compute_module_scores(
        gse_expression,
        strict_weights,
        broad_weights,
        manifest,
        GSE_ACCESSION,
    )
    gse_proliferation_scores, gse_proliferation_coverage = compute_proliferation_scores(
        gse_expression,
        proliferation_mapping,
        strict_weights,
        GSE_ACCESSION,
    )
    gse_scores = merge_score_components(
        gse_module_scores,
        gse_proliferation_scores,
    )
    gse_coverage_all = pd.concat(
        [gse_coverage, gse_proliferation_coverage],
        axis=0,
        ignore_index=True,
        sort=False,
    )

    gse_expression.to_csv(OUTPUT_GSE_EXPRESSION)
    gse_clinical.to_csv(OUTPUT_GSE_CLINICAL)
    gse_scores.to_csv(OUTPUT_GSE_SCORES)
    gse_coverage_all.to_csv(OUTPUT_GSE_COVERAGE, index=False)

    summary = pd.DataFrame(
        [
            {
                "cohort": "TARGET_OS",
                "n_expression_samples": target_expression.shape[0],
                "n_expression_genes": target_expression.shape[1],
                "n_clinical_rows": target_clinical.shape[0],
                "n_os_complete": int(
                    target_clinical[["os_time_days", "os_event"]].dropna().shape[0]
                ),
                "n_metastasis_labels": np.nan,
                "n_frozen_score_columns": target_scores.shape[1] - int("cohort" in target_scores.columns),
            },
            {
                "cohort": GSE_ACCESSION,
                "n_expression_samples": gse_expression.shape[0],
                "n_expression_genes": gse_expression.shape[1],
                "n_clinical_rows": gse_clinical.shape[0],
                "n_os_complete": int(
                    gse_clinical[["os_time_months", "os_event"]].dropna().shape[0]
                ),
                "n_metastasis_labels": int(gse_clinical["metastasis_within_5y"].notna().sum()),
                "n_frozen_score_columns": gse_scores.shape[1] - int("cohort" in gse_scores.columns),
            },
        ]
    )
    summary.to_csv(OUTPUT_PREPARATION_SUMMARY, index=False)
    write_readme()
    create_preparation_manifest(
        freeze=freeze,
        target_expression=target_expression,
        target_clinical=target_clinical,
        target_scores=target_scores,
        gse_expression=gse_expression,
        gse_clinical=gse_clinical,
        gse_scores=gse_scores,
    )

    print("")
    print("=" * 80)
    print("Human cohort preparation summary")
    print("=" * 80)
    print(summary.to_string(index=False))
    print_coverage_summary(target_coverage, "TARGET_OS")
    print_coverage_summary(gse_coverage, GSE_ACCESSION)

    print("")
    print("=" * 80)
    print("Endpoint preparation audit")
    print("=" * 80)
    print("TARGET-OS OS fields:")
    print(
        target_clinical[["os_time_days", "os_event", "vital_status"]]
        .agg(["count"])
        .to_string()
    )
    print("")
    print("GSE21257 metastasis labels:")
    print(gse_clinical["metastasis_within_5y"].value_counts(dropna=False).to_string())
    print("")
    print("GSE21257 OS fields:")
    print(
        gse_clinical[["os_time_months", "os_event"]]
        .agg(["count"])
        .to_string()
    )

    print("")
    print("=" * 80)
    print("Interpretation guardrails")
    print("=" * 80)
    print("No human outcome was used to select genes, orient scores, tune weights, or revise validation tiers.")
    print("GSE21257 probe collapsing used the highest-variance probe per unambiguous gene symbol without outcome labels.")
    print("TARGET-OS expression used open GDC STAR-Counts files and was harmonized at the gene-symbol level.")
    print("M40 residual scores are mechanistic sensitivity scores based on a disjoint human proliferation PC1.")
    print("Outcome association, multiplicity control, and external performance estimation are deferred to script 23.")

    print("")
    print("Saved:")
    for path in [
        OUTPUT_TARGET_EXPRESSION,
        OUTPUT_TARGET_CLINICAL,
        OUTPUT_TARGET_SAMPLE_MAP,
        OUTPUT_TARGET_SCORES,
        OUTPUT_TARGET_COVERAGE,
        OUTPUT_GSE_EXPRESSION,
        OUTPUT_GSE_CLINICAL,
        OUTPUT_GSE_SCORES,
        OUTPUT_GSE_COVERAGE,
        OUTPUT_PROLIFERATION_MAPPING,
        OUTPUT_PREPARATION_SUMMARY,
        OUTPUT_PREPARATION_MANIFEST,
        OUTPUT_README,
    ]:
        print(path)
    print("Done.")
```

### Relevant standalone lines

| Line | Signals | Code |
|---:|---|---|
| 179 | weights | `raise RuntimeError("The strict frozen weight table is empty.")` |
| 517 | zscore | `def standardize_target_clinical(cases: list[dict[str, Any]]) -> pd.DataFrame:` |
| 737 | loadings | `print("Downloading/loading GSE21257 from NCBI GEO:")` |
| 825 | std | `stds = x.std(axis=0).replace(0, np.nan)` |
| 826 | mean | `z = (x - x.mean(axis=0)) / stds` |
| 832 | std | `std = values.std()` |
| 835 | mean | `return (values - values.mean()) / std` |
| 840 | std | `if frame.shape[0] < 5 or frame["a"].std() == 0 or frame["b"].std() == 0:` |
| 842 | correlation | `return float(frame["a"].corr(frame["b"]))` |
| 851 | correlation | `corr = safe_corr(score, reference)` |
| 852 | correlation | `if np.isfinite(corr) and corr < 0:` |
| 871 | weights | `for mapping_name, weights in [("strict", strict_weights), ("broad", broad_weights)]:` |
| 872 | weights | `for module_label, part in weights.groupby("module_label"):` |
| 904 | direction | `signs = np.sign(` |
| 907 | mean | `signed_mean = z.mul(signs, axis=1).mean(axis=1)` |
| 915 | weights | `weighted = z.mul(normalized_weights, axis=1).sum(axis=1)` |
| 916 | weights | `weighted = zscore_series(weighted)` |
| 918 | weights | `weighted = pd.Series(np.nan, index=z.index)` |
| 923 | weights | `scores[f"{prefix}__canine_pca_weighted_z"] = weighted` |
| 926 | coverage | `coverage = pd.DataFrame(coverage_rows)` |
| 929 | coverage | `coverage["validation_tier"] = coverage["module_label"].map(tier_map)` |
| 930 | coverage | `return scores, coverage` |
| 1032 | mean | `reference = z[anchor_genes].mean(axis=1) if len(anchor_genes) >= 3 else z.mean(axis=1)` |
| 1067 | mean | `z_disjoint[anchor_disjoint].mean(axis=1)` |
| 1069 | mean | `else z_disjoint.mean(axis=1)` |
| 1080 | std | `if frame.shape[0] < 10 or frame["covariate"].std() == 0:` |
| 1084 | matrix_product | `residual.loc[frame.index] = frame["score"].values - x @ beta` |
| 1123 | weights, direction | `text = f"""Human osteosarcoma cohort preparation\n\nScript version: {SCRIPT_VERSION}\n\nThis script prepares TARGET-OS and GSE21257 without fitting any outcome model.\nFrozen canine module membership, score direction, PCA weights, and validation tiers are not changed.\n\nPrimary external scores:\n- strict one-to-one signed mean z-score for M34, M11, M24, and M40\n- TARGET-OS: overall-survival metadata prepared from public GDC clinical fields\n- GSE21257: metastasis-within-five-years label parsed from GEO metadata\n\nSecondary/sensitivity scores:\n- strict canine PCA-weighted score\n- broad mapped score\n- human-cohort PC1 oriented without outcomes\n- M40 residual to a disjoint strict human proliferation PC1\n\nNo survival or metastasis association is tested in this script.\n"""` |
| 1193 | coverage | `def print_coverage_summary(coverage: pd.DataFrame, cohort: str) -> None:` |
| 1196 | coverage | `print(f"Frozen score coverage: {cohort}")` |
| 1198 | coverage | `display = coverage.copy()` |
| 1263 | zscore | `target_clinical = standardize_target_clinical(target_cases)` |
| 1403 | weights | `print("No human outcome was used to select genes, orient scores, tune weights, or revise validation tiers.")` |

### Referenced result-table schemas

#### `results/tables/GSE21257_GEO_phenotype_raw.csv`

- Rows: 53
- SHA-256: `2ddad092cab45a8fb676eb09af041f953abebd30ad6a6903fecf30a62fd3ef00`
- Columns: ``, `title`, `geo_accession`, `status`, `submission_date`, `last_update_date`, `type`, `channel_count`, `source_name_ch1`, `organism_ch1`, `taxid_ch1`, `characteristics_ch1.0.age`, `characteristics_ch1.1.gender`, `characteristics_ch1.2.histological subtype`, `characteristics_ch1.3.tumor location`, `characteristics_ch1.4.huvos grade`, `characteristics_ch1.5.group`, `characteristics_ch1.6.status`, `molecule_ch1`, `extract_protocol_ch1`, `label_ch1`, `label_protocol_ch1`, `hyb_protocol`, `scan_protocol`, `description`, `data_processing`, `platform_id`, `contact_name`, `contact_laboratory`, `contact_institute`, `contact_address`, `contact_city`, `contact_zip/postal_code`, `contact_country`, `supplementary_file`, `series_id`, `data_row_count`

#### `results/tables/GSE21257_frozen_transfer_score_coverage.csv`

- Rows: 28
- SHA-256: `94e9d7e7df4c359a5cd8fb43a6a782df8acf0b53a2ed75265786a0c7921dd845`
- Columns: `cohort`, `module_label`, `mapping`, `n_frozen_genes`, `n_available_genes`, `coverage_fraction`, `minimum_rule_passed`, `available_genes`, `missing_genes`, `validation_tier`, `score`

#### `results/tables/GSE21257_probe_to_gene_symbol_selected.csv`

- Rows: 24996
- SHA-256: `cdd7e04e9e3f2ca4a9bebf626e0b52f7edd3039b5341d2de8bd30fa22132aea6`
- Columns: `probe_id`, `gene_symbol`, `probe_variance`

#### `results/tables/GSE238110_RNA_master_candidate_evidence_table_with_ortholog_qc.csv`

- Rows: 5013
- SHA-256: `26673d73d3e2b3769ae0eaf776d5f9fe239c594bf8c0f71e01b6ff50d1017586`
- Columns: `gene`, `gene_symbol_clean`, `univariate_candidate_rank`, `univariate_candidate_source`, `in_univariate_candidate_set`, `dfi_univ_rank`, `dfi_univ_coef`, `dfi_univ_hr_per_sd`, `dfi_univ_p`, `dfi_univ_q`, `dfi_univ_c_index`, `os_univ_rank`, `os_univ_coef`, `os_univ_hr_per_sd`, `os_univ_p`, `os_univ_q`, `os_univ_c_index`, `dfi_full_conditional_selected`, `dfi_full_conditional_rank`, `dfi_full_conditional_coef`, `dfi_full_conditional_hr_per_sd`, `dfi_full_conditional_p`, `dfi_full_conditional_c_index`, `os_full_conditional_selected`, `os_full_conditional_rank`, `os_full_conditional_coef`, `os_full_conditional_hr_per_sd`, `os_full_conditional_p`, `os_full_conditional_c_index`, `dfi_true_mb_n_configs`, `dfi_true_mb_algorithms`, `dfi_true_mb_alphas`, `dfi_true_mb_mean_rank`, `os_true_mb_n_configs`, `os_true_mb_algorithms`, `os_true_mb_alphas`, `os_true_mb_mean_rank`, `nested_dfi_conditional_selected_folds`, `nested_dfi_conditional_selection_frequency`, `nested_dfi_conditional_mean_rank`, `nested_dfi_conditional_plus_clinical_selected_folds`, `nested_dfi_conditional_plus_clinical_selection_frequency`, `nested_dfi_conditional_plus_clinical_mean_rank`, `nested_dfi_elasticnet_selected_folds`, `nested_dfi_elasticnet_selection_frequency`, `nested_dfi_elasticnet_mean_rank`, `nested_dfi_gsmb_selected_folds`, `nested_dfi_gsmb_selection_frequency`, `nested_dfi_gsmb_mean_rank`, `nested_dfi_iamb_selected_folds`, `nested_dfi_iamb_selection_frequency`, `nested_dfi_iamb_mean_rank`, `nested_dfi_univtop10_selected_folds`, `nested_dfi_univtop10_selection_frequency`, `nested_dfi_univtop10_mean_rank`, `nested_os_conditional_selected_folds`, `nested_os_conditional_selection_frequency`, `nested_os_conditional_mean_rank`, `nested_os_conditional_plus_clinical_selected_folds`, `nested_os_conditional_plus_clinical_selection_frequency`, `nested_os_conditional_plus_clinical_mean_rank`, `nested_os_elasticnet_selected_folds`, `nested_os_elasticnet_selection_frequency`, `nested_os_elasticnet_mean_rank`, `nested_os_gsmb_selected_folds`, `nested_os_gsmb_selection_frequency`, `nested_os_gsmb_mean_rank`, `nested_os_iamb_selected_folds`, `nested_os_iamb_selection_frequency`, `nested_os_iamb_mean_rank`, `nested_os_univtop10_selected_folds`, `nested_os_univtop10_selection_frequency`, `nested_os_univtop10_mean_rank`, `nested_any_max_selection_frequency`, `nested_any_total_selected_folds`, `dfi_os_both_univariate_available`, `dfi_os_same_direction`, `dfi_os_both_nominal_p05`, `dfi_os_both_fdr_q10`, `combined_univ_rank_score`, `rna_evidence_priority_score`, `rna_evidence_tier`, `gene_symbol_clean_upper`, `external_gene_name_upper`, `dog_ensembl_gene_id`, `dog_external_gene_name`, `dog_gene_biotype`, `human_ensembl_gene_id`, `human_gene_symbol`, `dog_human_orthology_type`, `dog_human_orthology_confidence`, `dog_human_orthology_confidence_numeric`, `dog_human_perc_id`, `human_dog_perc_id`, `has_human_homolog`, `is_one_to_one_ortholog`, `ortholog_mapping_status`, `dog_symbol_upper`, `human_symbol_upper`, `has_human_gene_symbol`, `human_symbol_same_as_dog_symbol`, `dog_symbol_problematic`, `human_symbol_problematic`, `ortholog_confidence_high`, `broad_transferable_ortholog`, `strict_transferable_ortholog`, `strict_symbol_concordant_transferable`, `needs_manual_ortholog_review`, `ortholog_qc_status`, `primary_human_validation_gene`, `sensitivity_human_validation_gene`

#### `results/tables/GSE238110_frozen_canine_transfer_program_manifest.csv`

- Rows: 13
- SHA-256: `fcd34b985ea99986f39dfc43ce8e4018db2fde9ecc8a225bb40ac93006d5372c`
- Columns: `module_label`, `validation_tier`, `multiplicity_family`, `provisional_program_label`, `program_label_requires_enrichment_confirmation`, `canine_primary_endpoint`, `canine_secondary_endpoint`, `positive_score_interpretation`, `n_canine_genes_used_for_pca`, `canine_pc1_explained_variance`, `canine_pca_orientation_correlation`, `dfi_full_cohort_coef`, `os_full_cohort_coef`, `risk_orientation_multiplier`, `n_strict_human_genes`, `n_broad_human_genes`, `strict_transfer_eligible`, `module_transfer_qc_tier`, `transfer_priority_score`, `fraction_strict_symbol_concordant`, `fraction_broad_transferable`, `raw_module_proliferation_correlation`, `orthogonal_variance_fraction_1_minus_r2`, `n_overlap_symbols_with_proliferation`, `primary_human_score`, `secondary_human_score`, `sensitivity_human_scores`, `frozen_after_canine_script`, `dfi_module_only_mean_c_index`, `dfi_module_only_std_c_index`, `dfi_module_only_fraction_above_0_50`, `dfi_module_plus_disjoint_proliferation_mean_c_index`, `dfi_module_plus_disjoint_proliferation_std_c_index`, `dfi_module_plus_disjoint_proliferation_fraction_above_0_50`, `dfi_residual_to_disjoint_proliferation_mean_c_index`, `dfi_residual_to_disjoint_proliferation_std_c_index`, `dfi_residual_to_disjoint_proliferation_fraction_above_0_50`, `dfi_residual_to_disjoint_proliferation_and_weight_mean_c_index`, `dfi_residual_to_disjoint_proliferation_and_weight_std_c_index`, `dfi_residual_to_disjoint_proliferation_and_weight_fraction_above_0_50`, `os_module_only_mean_c_index`, `os_module_only_std_c_index`, `os_module_only_fraction_above_0_50`, `os_module_plus_disjoint_proliferation_mean_c_index`, `os_module_plus_disjoint_proliferation_std_c_index`, `os_module_plus_disjoint_proliferation_fraction_above_0_50`, `os_residual_to_disjoint_proliferation_mean_c_index`, `os_residual_to_disjoint_proliferation_std_c_index`, `os_residual_to_disjoint_proliferation_fraction_above_0_50`, `os_residual_to_disjoint_proliferation_and_weight_mean_c_index`, `os_residual_to_disjoint_proliferation_and_weight_std_c_index`, `os_residual_to_disjoint_proliferation_and_weight_fraction_above_0_50`, `script20_dfi_recommended_role`, `script20_os_recommended_role`, `manual_freeze_reason`

#### `results/tables/GSE238110_frozen_transfer_gene_weights_broad.csv`

- Rows: 389
- SHA-256: `ee5f661677b173507e04856f848b6757ede2750dd69b9bb0d88a9f6a2a77bd4c`
- Columns: `module_label`, `canine_gene`, `canine_gene_symbol`, `human_gene_symbol`, `ortholog_qc_status`, `raw_pca_loading`, `risk_oriented_loading`, `absolute_risk_loading`, `is_strict_mapping`, `is_broad_mapping`, `human_symbol_duplicate_count`, `normalized_abs_sum_weight`

#### `results/tables/GSE238110_frozen_transfer_gene_weights_strict.csv`

- Rows: 321
- SHA-256: `4f065aa4c4edf117a0c74015840d2b4b2347929f172cd517e1818ba0f6163b91`
- Columns: `module_label`, `canine_gene`, `canine_gene_symbol`, `human_gene_symbol`, `ortholog_qc_status`, `raw_pca_loading`, `risk_oriented_loading`, `absolute_risk_loading`, `is_strict_mapping`, `is_broad_mapping`, `human_symbol_duplicate_count`, `normalized_abs_sum_weight`

#### `results/tables/GSE238110_frozen_transfer_scoring_specification.csv`

- Rows: 40
- SHA-256: `562aab92caa949691f43b0285cf2585c98008e78c62d2c15477a6d156d7f47a7`
- Columns: `module_label`, `validation_tier`, `score_name`, `analysis_role`, `gene_mapping`, `within_cohort_preprocessing`, `score_formula`, `minimum_gene_rule`, `outcome_use`

#### `results/tables/GSE238110_meta_proliferation_gene_set.csv`

- Rows: 251
- SHA-256: `d21df2d8c2877e0584e4c410e25b23802c47bb3c4642827dadc62047a3702f8a`
- Columns: `gene`, `gene_symbol_clean`, `is_anchor_gene`, `correlation_with_anchor_score`

#### `results/tables/TARGET_OS_GDC_expression_file_manifest.csv`

- Rows: 88
- SHA-256: `f90276eb869b8c135be0058d1da7f5cb2f9e898a7d7b8a389cd235d7a18b4db0`
- Columns: `file_id`, `file_name`, `md5sum`, `file_size`, `created_datetime`, `updated_datetime`, `case_id`, `case_submitter_id`, `sample_id`, `sample_submitter_id`, `sample_types`, `is_primary_tumor`

#### `results/tables/TARGET_OS_frozen_transfer_score_coverage.csv`

- Rows: 28
- SHA-256: `6881e1b3c5be414b3b177be802b417d9583407088a312b587dce7b1c9a7af35f`
- Columns: `cohort`, `module_label`, `mapping`, `n_frozen_genes`, `n_available_genes`, `coverage_fraction`, `minimum_rule_passed`, `available_genes`, `missing_genes`, `validation_tier`, `score`

#### `results/tables/frozen_strict_human_proliferation_mapping.csv`

- Rows: 111
- SHA-256: `72e539d40b3eeb604986e224a345f5218a42f5fcd70ea48b27d000c178413de8`
- Columns: `gene`, `canine_proliferation_gene`, `dog_symbol_key`, `human_gene_symbol`, `ortholog_qc_status`

#### `results/tables/human_validation_cohort_preparation_summary.csv`

- Rows: 2
- SHA-256: `7165fbb009e27c39d16169032a29c7856d7e3e32a02ced14e58af3da497ee8e0`
- Columns: `cohort`, `n_expression_samples`, `n_expression_genes`, `n_clinical_rows`, `n_os_complete`, `n_metastasis_labels`, `n_frozen_score_columns`

---

## `scripts/23_external_human_validation.py`

- SHA-256: `677c29509c5165782d24323abad7c24fc571ba3bc35000e801e6b137103df21e`

### Relevant functions

#### `standardize_series`

- Lines: 166-171
- Signals: `mean, std, zscore`

```python
def standardize_series(series: pd.Series) -> pd.Series:
    values = safe_numeric(series)
    sd = values.std()
    if not np.isfinite(sd) or sd == 0:
        return pd.Series(np.nan, index=series.index)
    return (values - values.mean()) / sd
```

#### `clean_scores`

- Lines: 174-179
- Signals: ``

```python
def clean_scores(scores: pd.DataFrame) -> pd.DataFrame:
    out = scores.copy()
    out = out.drop(columns=["cohort"], errors="ignore")
    for col in out.columns:
        out[col] = safe_numeric(out[col])
    return out
```

#### `fit_cox_score`

- Lines: 218-322
- Signals: `std, zscore, rank`

```python
def fit_cox_score(
    clinical: pd.DataFrame,
    score: pd.Series,
    time_col: str,
    event_col: str,
    cohort: str,
    endpoint: str,
    module_label: str,
    score_name: str,
    analysis_tier: str,
    covariates: pd.DataFrame | None = None,
    proliferation_score: pd.Series | None = None,
    rng: np.random.Generator | None = None,
) -> dict[str, Any]:
    frame = pd.DataFrame({
        "time": safe_numeric(clinical[time_col]),
        "event": safe_numeric(clinical[event_col]),
        "score": standardize_series(score.reindex(clinical.index)),
    })
    model_covariates = []
    if proliferation_score is not None:
        frame["proliferation"] = standardize_series(proliferation_score.reindex(clinical.index))
        model_covariates.append("proliferation")
    if covariates is not None and not covariates.empty:
        for col in covariates.columns:
            frame[col] = safe_numeric(covariates[col].reindex(clinical.index))
            model_covariates.append(col)

    frame = frame.replace([np.inf, -np.inf], np.nan).dropna()
    row: dict[str, Any] = {
        "cohort": cohort,
        "endpoint": endpoint,
        "module_label": module_label,
        "score_name": score_name,
        "analysis_tier": analysis_tier,
        "n": int(frame.shape[0]),
        "events": int(frame["event"].sum()) if not frame.empty else 0,
        "n_covariates": len(model_covariates),
        "covariates": ";".join(model_covariates),
        "score_coef": np.nan,
        "score_hr_per_sd": np.nan,
        "score_ci_low": np.nan,
        "score_ci_high": np.nan,
        "score_p": np.nan,
        "model_c_index": np.nan,
        "fixed_score_c_index": np.nan,
        "fixed_score_c_index_ci_low": np.nan,
        "fixed_score_c_index_ci_high": np.nan,
        "bootstrap_valid": 0,
        "ph_test_p": np.nan,
        "error": "",
    }
    if frame.shape[0] < MIN_ANALYSIS_N:
        row["error"] = "too_few_complete_samples"
        return row
    if frame["event"].sum() < MIN_SURVIVAL_EVENTS:
        row["error"] = "too_few_events"
        return row
    if frame["score"].std() == 0:
        row["error"] = "zero_variance_score"
        return row

    try:
        fixed_c = float(concordance_index(frame["time"], -frame["score"], frame["event"]))
        row["fixed_score_c_index"] = fixed_c
        if rng is not None:
            low, high, n_valid = bootstrap_c_index(
                frame["time"].to_numpy(float),
                frame["event"].to_numpy(int),
                frame["score"].to_numpy(float),
                N_BOOTSTRAP,
                rng,
            )
            row["fixed_score_c_index_ci_low"] = low
            row["fixed_score_c_index_ci_high"] = high
            row["bootstrap_valid"] = n_valid
    except Exception:
        pass

    cph = CoxPHFitter(penalizer=COX_PENALIZER)
    try:
        fit_cols = ["time", "event", "score"] + model_covariates
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", ConvergenceWarning)
            cph.fit(
                frame[fit_cols],
                duration_col="time",
                event_col="event",
                fit_options={"max_steps": 500},
            )
        s = cph.summary.loc["score"]
        row["score_coef"] = float(s["coef"])
        row["score_hr_per_sd"] = float(s["exp(coef)"])
        row["score_ci_low"] = float(s["exp(coef) lower 95%"])
        row["score_ci_high"] = float(s["exp(coef) upper 95%"])
        row["score_p"] = float(s["p"])
        row["model_c_index"] = float(cph.concordance_index_)
        try:
            ph = proportional_hazard_test(cph, frame[fit_cols], time_transform="rank")
            row["ph_test_p"] = float(ph.summary.loc["score", "p"])
        except Exception:
            pass
    except Exception as exc:
        row["error"] = str(exc)[:500]
    return row
```

#### `fit_logistic_score`

- Lines: 351-449
- Signals: `std, zscore`

```python
def fit_logistic_score(
    clinical: pd.DataFrame,
    score: pd.Series,
    outcome_col: str,
    cohort: str,
    endpoint: str,
    module_label: str,
    score_name: str,
    analysis_tier: str,
    covariates: pd.DataFrame | None = None,
    proliferation_score: pd.Series | None = None,
    rng: np.random.Generator | None = None,
) -> dict[str, Any]:
    frame = pd.DataFrame({
        "outcome": safe_numeric(clinical[outcome_col]),
        "score": standardize_series(score.reindex(clinical.index)),
    })
    model_covariates = []
    if proliferation_score is not None:
        frame["proliferation"] = standardize_series(proliferation_score.reindex(clinical.index))
        model_covariates.append("proliferation")
    if covariates is not None and not covariates.empty:
        for col in covariates.columns:
            frame[col] = safe_numeric(covariates[col].reindex(clinical.index))
            model_covariates.append(col)
    frame = frame.replace([np.inf, -np.inf], np.nan).dropna()
    frame = frame[frame["outcome"].isin([0, 1])].copy()

    row: dict[str, Any] = {
        "cohort": cohort,
        "endpoint": endpoint,
        "module_label": module_label,
        "score_name": score_name,
        "analysis_tier": analysis_tier,
        "n": int(frame.shape[0]),
        "positives": int((frame["outcome"] == 1).sum()),
        "negatives": int((frame["outcome"] == 0).sum()),
        "n_covariates": len(model_covariates),
        "covariates": ";".join(model_covariates),
        "auc": np.nan,
        "auc_ci_low": np.nan,
        "auc_ci_high": np.nan,
        "average_precision": np.nan,
        "mann_whitney_p_two_sided": np.nan,
        "logistic_coef": np.nan,
        "logistic_or_per_sd": np.nan,
        "logistic_ci_low": np.nan,
        "logistic_ci_high": np.nan,
        "logistic_p": np.nan,
        "bootstrap_valid": 0,
        "error": "",
    }
    if frame.shape[0] < MIN_ANALYSIS_N:
        row["error"] = "too_few_complete_samples"
        return row
    if row["positives"] < MIN_CLASS_POS or row["negatives"] < MIN_CLASS_NEG:
        row["error"] = "insufficient_class_counts"
        return row
    if frame["score"].std() == 0:
        row["error"] = "zero_variance_score"
        return row

    y = frame["outcome"].to_numpy(int)
    s = frame["score"].to_numpy(float)
    try:
        row["auc"] = float(roc_auc_score(y, s))
        row["average_precision"] = float(average_precision_score(y, s))
        if rng is not None:
            low, high, n_valid = stratified_bootstrap_auc(y, s, N_BOOTSTRAP, rng)
            row["auc_ci_low"] = low
            row["auc_ci_high"] = high
            row["bootstrap_valid"] = n_valid
        pos_scores = frame.loc[frame["outcome"] == 1, "score"]
        neg_scores = frame.loc[frame["outcome"] == 0, "score"]
        row["mann_whitney_p_two_sided"] = float(
            stats.mannwhitneyu(pos_scores, neg_scores, alternative="two-sided").pvalue
        )
    except Exception as exc:
        row["error"] = str(exc)[:500]
        return row

    if HAS_STATSMODELS:
        try:
            x_cols = ["score"] + model_covariates
            x = sm.add_constant(frame[x_cols], has_constant="add")
            model = sm.Logit(frame["outcome"], x).fit(disp=False, maxiter=500)
            coef = float(model.params["score"])
            se = float(model.bse["score"])
            row["logistic_coef"] = coef
            row["logistic_or_per_sd"] = float(np.exp(coef))
            row["logistic_ci_low"] = float(np.exp(coef - 1.96 * se))
            row["logistic_ci_high"] = float(np.exp(coef + 1.96 * se))
            row["logistic_p"] = float(model.pvalues["score"])
        except Exception as exc:
            if not row["error"]:
                row["error"] = f"logistic_fit_failed:{str(exc)[:400]}"
    else:
        row["error"] = "statsmodels_not_installed_logistic_effect_not_fitted"
    return row
```

#### `build_target_adjustment_covariates`

- Lines: 467-481
- Signals: `zscore`

```python
def build_target_adjustment_covariates(clinical: pd.DataFrame) -> pd.DataFrame:
    out = pd.DataFrame(index=clinical.index)
    if "age_at_diagnosis_years" in clinical.columns:
        age = safe_numeric(clinical["age_at_diagnosis_years"])
        if age.notna().sum() >= 50 and age.nunique(dropna=True) >= 5:
            out["age_z"] = standardize_series(age)
    if "sex" in clinical.columns:
        sex = clinical["sex"].astype(str).str.strip().str.lower()
        if sex.isin(["male", "female"]).sum() >= 50 and sex[sex.isin(["male", "female"])].nunique() == 2:
            out["sex_male"] = sex.map({"female": 0.0, "male": 1.0})
    if "metastasis_fields_raw" in clinical.columns:
        metastatic = parse_metastasis_at_diagnosis(clinical["metastasis_fields_raw"])
        if metastatic.notna().sum() >= 40 and metastatic.nunique(dropna=True) == 2:
            out["metastatic_at_diagnosis"] = metastatic
    return out
```

#### `build_gse_adjustment_covariates`

- Lines: 484-490
- Signals: `zscore`

```python
def build_gse_adjustment_covariates(clinical: pd.DataFrame) -> pd.DataFrame:
    out = pd.DataFrame(index=clinical.index)
    if "age_years" in clinical.columns:
        age = safe_numeric(clinical["age_years"])
        if age.notna().sum() >= 40 and age.nunique(dropna=True) >= 5:
            out["age_z"] = standardize_series(age)
    return out
```

#### `score_column`

- Lines: 499-500
- Signals: ``

```python
def score_column(module: str, suffix: str) -> str:
    return f"{module}{suffix}"
```

#### `run_secondary_programs`

- Lines: 627-664
- Signals: ``

```python
def run_secondary_programs(
    target_clinical: pd.DataFrame,
    target_scores: pd.DataFrame,
    gse_clinical: pd.DataFrame,
    gse_scores: pd.DataFrame,
    rng: np.random.Generator,
) -> pd.DataFrame:
    rows = []
    for module in SECONDARY_MODULES:
        col = score_column(module, PRIMARY_SCORE_SUFFIX)
        if col in target_scores.columns:
            rows.append(fit_cox_score(
                target_clinical, target_scores[col], "os_time_days", "os_event",
                "TARGET_OS", "overall_survival", module, col,
                "secondary_prespecified", rng=rng,
            ))
        if col in gse_scores.columns:
            rows.append(fit_logistic_score(
                gse_clinical, gse_scores[col], "metastasis_within_5y",
                "GSE21257", "metastasis_within_5y", module, col,
                "secondary_prespecified", rng=rng,
            ))
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    out["p_for_fdr"] = np.where(
        out["endpoint"].eq("overall_survival"),
        safe_numeric(out.get("score_p", pd.Series(index=out.index, dtype=float))),
        safe_numeric(out.get("logistic_p", pd.Series(index=out.index, dtype=float))).where(
            safe_numeric(out.get("logistic_p", pd.Series(index=out.index, dtype=float))).notna(),
            safe_numeric(out.get("mann_whitney_p_two_sided", pd.Series(index=out.index, dtype=float))),
        ),
    )
    out["q_within_secondary_endpoint"] = np.nan
    for endpoint in out["endpoint"].dropna().unique():
        mask = out["endpoint"].eq(endpoint)
        out.loc[mask, "q_within_secondary_endpoint"] = bh_adjust(out.loc[mask, "p_for_fdr"]).values
    return out
```

#### `run_score_variant_sensitivity`

- Lines: 667-721
- Signals: ``

```python
def run_score_variant_sensitivity(
    target_clinical: pd.DataFrame,
    target_scores: pd.DataFrame,
    gse_clinical: pd.DataFrame,
    gse_scores: pd.DataFrame,
    rng: np.random.Generator,
) -> pd.DataFrame:
    rows = []
    variants = [
        (STRICT_WEIGHTED_SUFFIX, "strict_canine_pca_weighted"),
        (STRICT_PC1_SUFFIX, "strict_human_pc1"),
        (BROAD_SIGNED_SUFFIX, "broad_signed_mean"),
        (BROAD_WEIGHTED_SUFFIX, "broad_canine_pca_weighted"),
    ]
    for module in PRIMARY_MODULES:
        for suffix, variant_label in variants:
            col = score_column(module, suffix)
            if col in target_scores.columns:
                row = fit_cox_score(
                    target_clinical, target_scores[col], "os_time_days", "os_event",
                    "TARGET_OS", "overall_survival", module, col,
                    "score_variant_sensitivity", rng=rng,
                )
                row["score_variant"] = variant_label
                rows.append(row)
            if col in gse_scores.columns:
                row = fit_logistic_score(
                    gse_clinical, gse_scores[col], "metastasis_within_5y",
                    "GSE21257", "metastasis_within_5y", module, col,
                    "score_variant_sensitivity", rng=rng,
                )
                row["score_variant"] = variant_label
                rows.append(row)

    for col, variant_label in [
        (M40_RESIDUAL_SIGNED, "m40_residual_signed"),
        (M40_RESIDUAL_WEIGHTED, "m40_residual_weighted"),
    ]:
        if col in target_scores.columns:
            row = fit_cox_score(
                target_clinical, target_scores[col], "os_time_days", "os_event",
                "TARGET_OS", "overall_survival", "M40", col,
                "mechanistic_sensitivity", rng=rng,
            )
            row["score_variant"] = variant_label
            rows.append(row)
        if col in gse_scores.columns:
            row = fit_logistic_score(
                gse_clinical, gse_scores[col], "metastasis_within_5y",
                "GSE21257", "metastasis_within_5y", "M40", col,
                "mechanistic_sensitivity", rng=rng,
            )
            row["score_variant"] = variant_label
            rows.append(row)
    return pd.DataFrame(rows)
```

#### `strict_ortholog_universe`

- Lines: 798-814
- Signals: ``

```python
def strict_ortholog_universe(ortholog_qc: pd.DataFrame) -> list[str]:
    human_col = "human_gene_symbol" if "human_gene_symbol" in ortholog_qc.columns else "human_symbol"
    status_col = "ortholog_qc_status"
    if human_col not in ortholog_qc.columns or status_col not in ortholog_qc.columns:
        raise ValueError("Ortholog QC table lacks human symbol or QC status.")
    genes = (
        ortholog_qc.loc[
            ortholog_qc[status_col].eq("strict_symbol_concordant_one_to_one"),
            human_col,
        ]
        .dropna()
        .astype(str)
        .str.upper()
        .drop_duplicates()
        .tolist()
    )
    return genes
```

#### `zscore_expression`

- Lines: 817-826
- Signals: `mean, std`

```python
def zscore_expression(expression: pd.DataFrame) -> pd.DataFrame:
    x = expression.copy()
    x.columns = x.columns.astype(str).str.upper()
    x = x.loc[:, ~x.columns.duplicated()].copy()
    x = x.apply(pd.to_numeric, errors="coerce")
    x = x.replace([np.inf, -np.inf], np.nan)
    x = x.fillna(x.median(axis=0))
    sd = x.std(axis=0).replace(0, np.nan)
    z = (x - x.mean(axis=0)) / sd
    return z.loc[:, z.notna().all(axis=0)]
```

#### `assign_expression_bins`

- Lines: 829-855
- Signals: `mean, rank`

```python
def assign_expression_bins(expression: pd.DataFrame, n_bins: int = 5) -> pd.DataFrame:
    x = expression.apply(pd.to_numeric, errors="coerce")
    x = x.replace([np.inf, -np.inf], np.nan)
    x = x.fillna(x.median(axis=0))
    mean_expression = x.mean(axis=0)
    variance = x.var(axis=0)
    info = pd.DataFrame({"mean_expression": mean_expression, "variance": variance})
    try:
        info["mean_bin"] = pd.qcut(
            info["mean_expression"].rank(method="first"),
            n_bins,
            labels=False,
            duplicates="drop",
        )
    except Exception:
        info["mean_bin"] = 0
    try:
        info["var_bin"] = pd.qcut(
            info["variance"].rank(method="first"),
            n_bins,
            labels=False,
            duplicates="drop",
        )
    except Exception:
        info["var_bin"] = 0
    info["bin_key"] = info["mean_bin"].astype(str) + "_" + info["var_bin"].astype(str)
    return info
```

#### `random_control_for_module`

- Lines: 886-965
- Signals: `mean, zscore, direction`

```python
def random_control_for_module(
    cohort: str,
    endpoint: str,
    module: str,
    expression: pd.DataFrame,
    clinical: pd.DataFrame,
    strict_weights: pd.DataFrame,
    strict_universe: list[str],
    observed_metric: float,
    rng: np.random.Generator,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    raw = expression.copy()
    raw.columns = raw.columns.astype(str).str.upper()
    raw = raw.loc[:, ~raw.columns.duplicated()].copy()
    raw = raw.apply(pd.to_numeric, errors="coerce")
    raw = raw.replace([np.inf, -np.inf], np.nan)
    raw = raw.fillna(raw.median(axis=0))
    z = zscore_expression(raw)
    universe = [gene for gene in strict_universe if gene in z.columns]
    part = strict_weights[strict_weights["module_label"].eq(module)].copy()
    part["human_gene_symbol"] = part["human_gene_symbol"].astype(str).str.upper()
    part = part.drop_duplicates("human_gene_symbol", keep="first")
    part = part[part["human_gene_symbol"].isin(z.columns)].copy()
    module_genes = part["human_gene_symbol"].tolist()
    signs = np.sign(safe_numeric(part["risk_oriented_loading"])).replace(0, 1).to_numpy(float)
    bin_info = assign_expression_bins(raw[universe])
    rows = []
    metrics = []

    for repeat in range(1, N_RANDOM_GENE_SETS + 1):
        genes = sample_matched_random_genes(module_genes, universe, bin_info, rng)
        if len(genes) != len(module_genes) or len(genes) < 3:
            continue
        signs_perm = signs.copy()
        rng.shuffle(signs_perm)
        random_score = z[genes].mul(signs_perm, axis=1).mean(axis=1)
        random_score = standardize_series(random_score)
        metric = np.nan
        if endpoint == "overall_survival":
            frame = pd.DataFrame({
                "time": safe_numeric(clinical["os_time_days"]),
                "event": safe_numeric(clinical["os_event"]),
                "score": random_score.reindex(clinical.index),
            }).dropna()
            if frame.shape[0] >= MIN_ANALYSIS_N and frame["event"].sum() >= MIN_SURVIVAL_EVENTS:
                metric = float(concordance_index(frame["time"], -frame["score"], frame["event"]))
        elif endpoint == "metastasis_within_5y":
            frame = pd.DataFrame({
                "outcome": safe_numeric(clinical["metastasis_within_5y"]),
                "score": random_score.reindex(clinical.index),
            }).dropna()
            frame = frame[frame["outcome"].isin([0, 1])]
            if frame["outcome"].nunique() == 2:
                metric = float(roc_auc_score(frame["outcome"], frame["score"]))
        if np.isfinite(metric):
            metrics.append(metric)
            rows.append({
                "cohort": cohort,
                "endpoint": endpoint,
                "module_label": module,
                "repeat": repeat,
                "panel_size": len(genes),
                "random_metric": metric,
            })

    summary = {
        "cohort": cohort,
        "endpoint": endpoint,
        "module_label": module,
        "observed_metric": observed_metric,
        "n_module_genes_available": len(module_genes),
        "n_random_valid": len(metrics),
        "random_mean": float(np.mean(metrics)) if metrics else np.nan,
        "random_median": float(np.median(metrics)) if metrics else np.nan,
        "random_q90": float(np.quantile(metrics, 0.90)) if metrics else np.nan,
        "random_q95": float(np.quantile(metrics, 0.95)) if metrics else np.nan,
        "observed_percentile": float(np.mean(np.asarray(metrics) <= observed_metric)) if metrics and np.isfinite(observed_metric) else np.nan,
        "empirical_p_greater_equal": float((1 + np.sum(np.asarray(metrics) >= observed_metric)) / (1 + len(metrics))) if metrics and np.isfinite(observed_metric) else np.nan,
    }
    return summary, rows
```

#### `write_readme`

- Lines: 1082-1116
- Signals: `weights, direction`

```python
def write_readme() -> None:
    text = f"""Human external validation of frozen canine osteosarcoma programs

Script version: {SCRIPT_VERSION}

Primary external score:
- strict one-to-one ortholog signed-mean z-score
- fixed canine risk direction
- no human outcome used for gene selection, weighting, score orientation, or validation-tier revision

Primary external settings:
1. TARGET-OS overall survival: continuous fixed score, Cox HR per SD, and fixed-score Harrell C-index
2. GSE21257 metastasis within five years: continuous fixed score, logistic OR per SD, ROC-AUC, and PR-AUC

Multiplicity:
- BH correction within each primary setting across M34, M11, M24, and M40
- additional global BH correction across all eight primary tests

Sensitivity analyses:
- strict canine-PCA weighted score
- broad mapped score
- human-cohort PC1
- available clinical adjustment
- proliferation adjustment
- M40 residual to disjoint proliferation
- GSE21257 overall-survival association
- expression-matched random gene-set controls

Interpretation:
- External association is not proof of causality or clinical utility.
- GSE21257 is small; metastasis and OS results require replication.
- TARGET-OS has limited sample size and public clinical covariates.
- Random-gene-set controls are descriptive specificity diagnostics.
"""
    OUTPUT_README.write_text(text, encoding="utf-8")
```

#### `save_manifest`

- Lines: 1119-1146
- Signals: `weights, direction`

```python
def save_manifest(output_paths: list[Path], input_manifest: dict[str, Any]) -> None:
    manifest = {
        "script_version": SCRIPT_VERSION,
        "created_utc": utc_now_iso(),
        "preparation_manifest_sha256": sha256_file(PREPARATION_MANIFEST_FILE),
        "freeze_manifest_sha256": sha256_file(FREEZE_JSON_FILE),
        "primary_modules": PRIMARY_MODULES,
        "secondary_modules": SECONDARY_MODULES,
        "primary_score_suffix": PRIMARY_SCORE_SUFFIX,
        "n_bootstrap": N_BOOTSTRAP,
        "n_random_gene_sets": N_RANDOM_GENE_SETS,
        "random_seed": RANDOM_SEED,
        "multiplicity": {
            "within_setting": "Benjamini-Hochberg across four primary modules",
            "global": "Benjamini-Hochberg across eight primary tests",
        },
        "guardrail": "No human outcome revised any frozen gene set, score direction, weight, or validation tier.",
        "input_preparation_manifest": input_manifest,
        "files": {},
    }
    for path in output_paths:
        if path.exists():
            manifest["files"][path.name] = {
                "path": str(path),
                "sha256": sha256_file(path),
                "size_bytes": path.stat().st_size,
            }
    OUTPUT_MANIFEST.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
```

#### `main`

- Lines: 1186-1316
- Signals: `weights, direction`

```python
def main() -> None:
    print("=" * 80)
    print("External human validation of frozen canine osteosarcoma programs")
    print("=" * 80)
    print(f"Script version: {SCRIPT_VERSION}")
    print(f"Project root: {PROJECT_ROOT}")
    print(f"Human processed directory: {HUMAN_DIR}")
    print(f"Results directory: {RESULTS_DIR}")
    print("")
    print("Design:")
    print("  Test frozen strict ortholog scores in independent human osteosarcoma cohorts.")
    print("  Use TARGET-OS overall survival and GSE21257 metastasis within five years.")
    print("  Preserve canine score direction and validation tiers.")
    print("  Apply endpoint-specific and global multiplicity control.")
    print("  Run clinical, proliferation, score-variant, and random-gene-set sensitivity analyses.")

    preparation_manifest = verify_preparation_assets()
    if not FREEZE_JSON_FILE.exists():
        raise FileNotFoundError(f"Freeze file not found: {FREEZE_JSON_FILE}")

    manifest = read_required_csv(FROZEN_MANIFEST_FILE, index_col=None)
    strict_weights = read_required_csv(STRICT_WEIGHTS_FILE, index_col=None)
    broad_weights = read_required_csv(BROAD_WEIGHTS_FILE, index_col=None)
    ortholog_qc = read_required_csv(ORTHOLOG_QC_FILE, index_col=None)

    target_expression = read_required_csv(TARGET_EXPRESSION_FILE, index_col=0)
    target_clinical = read_required_csv(TARGET_CLINICAL_FILE, index_col=0)
    target_scores = clean_scores(read_required_csv(TARGET_SCORES_FILE, index_col=0))
    gse_expression = read_required_csv(GSE_EXPRESSION_FILE, index_col=0)
    gse_clinical = read_required_csv(GSE_CLINICAL_FILE, index_col=0)
    gse_scores = clean_scores(read_required_csv(GSE_SCORES_FILE, index_col=0))

    target_clinical, target_scores, target_expression = match_cohort(
        target_clinical, target_scores, target_expression
    )
    gse_clinical, gse_scores, gse_expression = match_cohort(
        gse_clinical, gse_scores, gse_expression
    )

    print("")
    print("Matched human validation data:")
    print(f"  TARGET-OS expression: {target_expression.shape}")
    print(f"  TARGET-OS clinical: {target_clinical.shape}")
    print(f"  TARGET-OS scores: {target_scores.shape}")
    target_complete = target_clinical[["os_time_days", "os_event"]].dropna()
    print(f"  TARGET-OS OS complete: {target_complete.shape[0]}")
    print(f"  TARGET-OS OS events: {int(safe_numeric(target_complete['os_event']).sum())}")
    print(f"  TARGET-OS censored: {int((safe_numeric(target_complete['os_event']) == 0).sum())}")
    print(f"  GSE21257 expression: {gse_expression.shape}")
    print(f"  GSE21257 clinical: {gse_clinical.shape}")
    print(f"  GSE21257 scores: {gse_scores.shape}")
    print(f"  GSE21257 metastasis positive: {int((safe_numeric(gse_clinical['metastasis_within_5y']) == 1).sum())}")
    print(f"  GSE21257 metastasis negative: {int((safe_numeric(gse_clinical['metastasis_within_5y']) == 0).sum())}")
    gse_os_complete = gse_clinical[["os_time_months", "os_event"]].dropna()
    print(f"  GSE21257 OS complete: {gse_os_complete.shape[0]}")
    print(f"  GSE21257 OS events: {int(safe_numeric(gse_os_complete['os_event']).sum())}")

    rng = np.random.default_rng(RANDOM_SEED)

    target_primary = run_primary_survival(
        target_clinical, target_scores, "TARGET_OS", "os_time_days", "os_event", rng
    )
    gse_primary = run_primary_metastasis(gse_clinical, gse_scores, rng)
    target_primary, gse_primary, multiplicity = add_primary_multiplicity(
        target_primary, gse_primary
    )

    secondary = run_secondary_programs(
        target_clinical, target_scores, gse_clinical, gse_scores, rng
    )
    score_variants = run_score_variant_sensitivity(
        target_clinical, target_scores, gse_clinical, gse_scores, rng
    )
    adjusted = run_adjusted_models(
        target_clinical, target_scores, gse_clinical, gse_scores
    )
    gse_os = run_gse_os_sensitivity(gse_clinical, gse_scores, rng)

    random_summary, random_distribution = run_random_controls(
        target_expression, target_clinical, target_primary,
        gse_expression, gse_clinical, gse_primary,
        strict_weights, ortholog_qc, rng,
    )

    synthesis = build_cross_cohort_synthesis(
        target_primary, gse_primary, gse_os, random_summary
    )

    target_primary.to_csv(OUTPUT_TARGET_PRIMARY, index=False)
    gse_primary.to_csv(OUTPUT_GSE_MET_PRIMARY, index=False)
    gse_os.to_csv(OUTPUT_GSE_OS_SENSITIVITY, index=False)
    secondary.to_csv(OUTPUT_SECONDARY, index=False)
    score_variants.to_csv(OUTPUT_SCORE_VARIANTS, index=False)
    adjusted.to_csv(OUTPUT_ADJUSTED, index=False)
    multiplicity.to_csv(OUTPUT_MULTIPLICITY, index=False)
    random_summary.to_csv(OUTPUT_RANDOM_CONTROLS, index=False)
    random_distribution.to_csv(OUTPUT_RANDOM_CONTROL_DISTRIBUTION, index=False)
    synthesis.to_csv(OUTPUT_SYNTHESIS, index=False)
    write_readme()

    output_paths = [
        OUTPUT_TARGET_PRIMARY,
        OUTPUT_GSE_MET_PRIMARY,
        OUTPUT_GSE_OS_SENSITIVITY,
        OUTPUT_SECONDARY,
        OUTPUT_SCORE_VARIANTS,
        OUTPUT_ADJUSTED,
        OUTPUT_MULTIPLICITY,
        OUTPUT_RANDOM_CONTROLS,
        OUTPUT_RANDOM_CONTROL_DISTRIBUTION,
        OUTPUT_SYNTHESIS,
        OUTPUT_README,
    ]
    save_manifest(output_paths, preparation_manifest)

    print_primary_results(target_primary, gse_primary, synthesis)

    print("")
    print("=" * 80)
    print("Interpretation guardrails")
    print("=" * 80)
    print("Primary inference uses only strict one-to-one signed-mean scores with frozen canine risk direction.")
    print("TARGET-OS and GSE21257 are treated as distinct primary external settings; global FDR across eight tests is also reported.")
    print("GSE21257 OS, broad mappings, human PC1, weighted scores, clinical adjustment, and proliferation adjustment are sensitivity analyses.")
    print("Random gene-set percentiles are descriptive specificity controls, not independent external cohorts.")
    print("External association does not establish causality, treatment response, or clinical readiness.")
    print("")
    print("Saved:")
    for path in output_paths + [OUTPUT_MANIFEST]:
        print(path)
    print("Done.")
```

### Relevant standalone lines

| Line | Signals | Code |
|---:|---|---|
| 166 | zscore | `def standardize_series(series: pd.Series) -> pd.Series:` |
| 168 | std | `sd = values.std()` |
| 171 | mean | `return (values - values.mean()) / sd` |
| 235 | zscore | `"score": standardize_series(score.reindex(clinical.index)),` |
| 239 | zscore | `frame["proliferation"] = standardize_series(proliferation_score.reindex(clinical.index))` |
| 276 | std | `if frame["score"].std() == 0:` |
| 316 | rank | `ph = proportional_hazard_test(cph, frame[fit_cols], time_transform="rank")` |
| 366 | zscore | `"score": standardize_series(score.reindex(clinical.index)),` |
| 370 | zscore | `frame["proliferation"] = standardize_series(proliferation_score.reindex(clinical.index))` |
| 409 | std | `if frame["score"].std() == 0:` |
| 472 | zscore | `out["age_z"] = standardize_series(age)` |
| 489 | zscore | `out["age_z"] = standardize_series(age)` |
| 824 | std | `sd = x.std(axis=0).replace(0, np.nan)` |
| 825 | mean | `z = (x - x.mean(axis=0)) / sd` |
| 833 | mean | `mean_expression = x.mean(axis=0)` |
| 838 | rank | `info["mean_expression"].rank(method="first"),` |
| 847 | rank | `info["variance"].rank(method="first"),` |
| 910 | direction | `signs = np.sign(safe_numeric(part["risk_oriented_loading"])).replace(0, 1).to_numpy(float)` |
| 921 | mean | `random_score = z[genes].mul(signs_perm, axis=1).mean(axis=1)` |
| 922 | zscore | `random_score = standardize_series(random_score)` |
| 958 | mean | `"random_mean": float(np.mean(metrics)) if metrics else np.nan,` |
| 962 | mean | `"observed_percentile": float(np.mean(np.asarray(metrics) <= observed_metric)) if metrics and np.isfinite(observed_metric) else np.nan,` |
| 1088 | direction | `- strict one-to-one ortholog signed-mean z-score` |
| 1089 | direction | `- fixed canine risk direction` |
| 1101 | weights | `- strict canine-PCA weighted score` |
| 1135 | weights, direction | `"guardrail": "No human outcome revised any frozen gene set, score direction, weight, or validation tier.",` |
| 1198 | direction | `print("  Preserve canine score direction and validation tiers.")` |
| 1307 | direction | `print("Primary inference uses only strict one-to-one signed-mean scores with frozen canine risk direction.")` |
| 1309 | weights | `print("GSE21257 OS, broad mappings, human PC1, weighted scores, clinical adjustment, and proliferation adjustment are sensitivity analyses.")` |

### Referenced result-table schemas

#### `results/tables/GSE21257_OS_frozen_program_sensitivity.csv`

- Rows: 7
- SHA-256: `f9052c4e43dce2daa21b0132328e0636e1b0959d9309257a124870cf523a93a8`
- Columns: `cohort`, `endpoint`, `module_label`, `score_name`, `analysis_tier`, `n`, `events`, `n_covariates`, `covariates`, `score_coef`, `score_hr_per_sd`, `score_ci_low`, `score_ci_high`, `score_p`, `model_c_index`, `fixed_score_c_index`, `fixed_score_c_index_ci_low`, `fixed_score_c_index_ci_high`, `bootstrap_valid`, `ph_test_p`, `error`, `q_within_gse_os_sensitivity`

#### `results/tables/GSE21257_frozen_transfer_score_coverage.csv`

- Rows: 28
- SHA-256: `94e9d7e7df4c359a5cd8fb43a6a782df8acf0b53a2ed75265786a0c7921dd845`
- Columns: `cohort`, `module_label`, `mapping`, `n_frozen_genes`, `n_available_genes`, `coverage_fraction`, `minimum_rule_passed`, `available_genes`, `missing_genes`, `validation_tier`, `score`

#### `results/tables/GSE21257_metastasis_primary_frozen_program_validation.csv`

- Rows: 4
- SHA-256: `66f8f57539d9fd9d666dd4fd27ffa017d0884fd694b68e0ec7fa5243842a2ac8`
- Columns: `cohort`, `endpoint`, `module_label`, `score_name`, `analysis_tier`, `n`, `positives`, `negatives`, `n_covariates`, `covariates`, `auc`, `auc_ci_low`, `auc_ci_high`, `average_precision`, `mann_whitney_p_two_sided`, `logistic_coef`, `logistic_or_per_sd`, `logistic_ci_low`, `logistic_ci_high`, `logistic_p`, `bootstrap_valid`, `error`, `primary_p`, `q_within_endpoint`, `direction_consistent`, `q_global_eight_tests`, `external_support_class`

#### `results/tables/GSE238110_RNA_master_candidate_evidence_table_with_ortholog_qc.csv`

- Rows: 5013
- SHA-256: `26673d73d3e2b3769ae0eaf776d5f9fe239c594bf8c0f71e01b6ff50d1017586`
- Columns: `gene`, `gene_symbol_clean`, `univariate_candidate_rank`, `univariate_candidate_source`, `in_univariate_candidate_set`, `dfi_univ_rank`, `dfi_univ_coef`, `dfi_univ_hr_per_sd`, `dfi_univ_p`, `dfi_univ_q`, `dfi_univ_c_index`, `os_univ_rank`, `os_univ_coef`, `os_univ_hr_per_sd`, `os_univ_p`, `os_univ_q`, `os_univ_c_index`, `dfi_full_conditional_selected`, `dfi_full_conditional_rank`, `dfi_full_conditional_coef`, `dfi_full_conditional_hr_per_sd`, `dfi_full_conditional_p`, `dfi_full_conditional_c_index`, `os_full_conditional_selected`, `os_full_conditional_rank`, `os_full_conditional_coef`, `os_full_conditional_hr_per_sd`, `os_full_conditional_p`, `os_full_conditional_c_index`, `dfi_true_mb_n_configs`, `dfi_true_mb_algorithms`, `dfi_true_mb_alphas`, `dfi_true_mb_mean_rank`, `os_true_mb_n_configs`, `os_true_mb_algorithms`, `os_true_mb_alphas`, `os_true_mb_mean_rank`, `nested_dfi_conditional_selected_folds`, `nested_dfi_conditional_selection_frequency`, `nested_dfi_conditional_mean_rank`, `nested_dfi_conditional_plus_clinical_selected_folds`, `nested_dfi_conditional_plus_clinical_selection_frequency`, `nested_dfi_conditional_plus_clinical_mean_rank`, `nested_dfi_elasticnet_selected_folds`, `nested_dfi_elasticnet_selection_frequency`, `nested_dfi_elasticnet_mean_rank`, `nested_dfi_gsmb_selected_folds`, `nested_dfi_gsmb_selection_frequency`, `nested_dfi_gsmb_mean_rank`, `nested_dfi_iamb_selected_folds`, `nested_dfi_iamb_selection_frequency`, `nested_dfi_iamb_mean_rank`, `nested_dfi_univtop10_selected_folds`, `nested_dfi_univtop10_selection_frequency`, `nested_dfi_univtop10_mean_rank`, `nested_os_conditional_selected_folds`, `nested_os_conditional_selection_frequency`, `nested_os_conditional_mean_rank`, `nested_os_conditional_plus_clinical_selected_folds`, `nested_os_conditional_plus_clinical_selection_frequency`, `nested_os_conditional_plus_clinical_mean_rank`, `nested_os_elasticnet_selected_folds`, `nested_os_elasticnet_selection_frequency`, `nested_os_elasticnet_mean_rank`, `nested_os_gsmb_selected_folds`, `nested_os_gsmb_selection_frequency`, `nested_os_gsmb_mean_rank`, `nested_os_iamb_selected_folds`, `nested_os_iamb_selection_frequency`, `nested_os_iamb_mean_rank`, `nested_os_univtop10_selected_folds`, `nested_os_univtop10_selection_frequency`, `nested_os_univtop10_mean_rank`, `nested_any_max_selection_frequency`, `nested_any_total_selected_folds`, `dfi_os_both_univariate_available`, `dfi_os_same_direction`, `dfi_os_both_nominal_p05`, `dfi_os_both_fdr_q10`, `combined_univ_rank_score`, `rna_evidence_priority_score`, `rna_evidence_tier`, `gene_symbol_clean_upper`, `external_gene_name_upper`, `dog_ensembl_gene_id`, `dog_external_gene_name`, `dog_gene_biotype`, `human_ensembl_gene_id`, `human_gene_symbol`, `dog_human_orthology_type`, `dog_human_orthology_confidence`, `dog_human_orthology_confidence_numeric`, `dog_human_perc_id`, `human_dog_perc_id`, `has_human_homolog`, `is_one_to_one_ortholog`, `ortholog_mapping_status`, `dog_symbol_upper`, `human_symbol_upper`, `has_human_gene_symbol`, `human_symbol_same_as_dog_symbol`, `dog_symbol_problematic`, `human_symbol_problematic`, `ortholog_confidence_high`, `broad_transferable_ortholog`, `strict_transferable_ortholog`, `strict_symbol_concordant_transferable`, `needs_manual_ortholog_review`, `ortholog_qc_status`, `primary_human_validation_gene`, `sensitivity_human_validation_gene`

#### `results/tables/GSE238110_frozen_canine_transfer_program_manifest.csv`

- Rows: 13
- SHA-256: `fcd34b985ea99986f39dfc43ce8e4018db2fde9ecc8a225bb40ac93006d5372c`
- Columns: `module_label`, `validation_tier`, `multiplicity_family`, `provisional_program_label`, `program_label_requires_enrichment_confirmation`, `canine_primary_endpoint`, `canine_secondary_endpoint`, `positive_score_interpretation`, `n_canine_genes_used_for_pca`, `canine_pc1_explained_variance`, `canine_pca_orientation_correlation`, `dfi_full_cohort_coef`, `os_full_cohort_coef`, `risk_orientation_multiplier`, `n_strict_human_genes`, `n_broad_human_genes`, `strict_transfer_eligible`, `module_transfer_qc_tier`, `transfer_priority_score`, `fraction_strict_symbol_concordant`, `fraction_broad_transferable`, `raw_module_proliferation_correlation`, `orthogonal_variance_fraction_1_minus_r2`, `n_overlap_symbols_with_proliferation`, `primary_human_score`, `secondary_human_score`, `sensitivity_human_scores`, `frozen_after_canine_script`, `dfi_module_only_mean_c_index`, `dfi_module_only_std_c_index`, `dfi_module_only_fraction_above_0_50`, `dfi_module_plus_disjoint_proliferation_mean_c_index`, `dfi_module_plus_disjoint_proliferation_std_c_index`, `dfi_module_plus_disjoint_proliferation_fraction_above_0_50`, `dfi_residual_to_disjoint_proliferation_mean_c_index`, `dfi_residual_to_disjoint_proliferation_std_c_index`, `dfi_residual_to_disjoint_proliferation_fraction_above_0_50`, `dfi_residual_to_disjoint_proliferation_and_weight_mean_c_index`, `dfi_residual_to_disjoint_proliferation_and_weight_std_c_index`, `dfi_residual_to_disjoint_proliferation_and_weight_fraction_above_0_50`, `os_module_only_mean_c_index`, `os_module_only_std_c_index`, `os_module_only_fraction_above_0_50`, `os_module_plus_disjoint_proliferation_mean_c_index`, `os_module_plus_disjoint_proliferation_std_c_index`, `os_module_plus_disjoint_proliferation_fraction_above_0_50`, `os_residual_to_disjoint_proliferation_mean_c_index`, `os_residual_to_disjoint_proliferation_std_c_index`, `os_residual_to_disjoint_proliferation_fraction_above_0_50`, `os_residual_to_disjoint_proliferation_and_weight_mean_c_index`, `os_residual_to_disjoint_proliferation_and_weight_std_c_index`, `os_residual_to_disjoint_proliferation_and_weight_fraction_above_0_50`, `script20_dfi_recommended_role`, `script20_os_recommended_role`, `manual_freeze_reason`

#### `results/tables/GSE238110_frozen_transfer_gene_weights_broad.csv`

- Rows: 389
- SHA-256: `ee5f661677b173507e04856f848b6757ede2750dd69b9bb0d88a9f6a2a77bd4c`
- Columns: `module_label`, `canine_gene`, `canine_gene_symbol`, `human_gene_symbol`, `ortholog_qc_status`, `raw_pca_loading`, `risk_oriented_loading`, `absolute_risk_loading`, `is_strict_mapping`, `is_broad_mapping`, `human_symbol_duplicate_count`, `normalized_abs_sum_weight`

#### `results/tables/GSE238110_frozen_transfer_gene_weights_strict.csv`

- Rows: 321
- SHA-256: `4f065aa4c4edf117a0c74015840d2b4b2347929f172cd517e1818ba0f6163b91`
- Columns: `module_label`, `canine_gene`, `canine_gene_symbol`, `human_gene_symbol`, `ortholog_qc_status`, `raw_pca_loading`, `risk_oriented_loading`, `absolute_risk_loading`, `is_strict_mapping`, `is_broad_mapping`, `human_symbol_duplicate_count`, `normalized_abs_sum_weight`

#### `results/tables/TARGET_OS_frozen_transfer_score_coverage.csv`

- Rows: 28
- SHA-256: `6881e1b3c5be414b3b177be802b417d9583407088a312b587dce7b1c9a7af35f`
- Columns: `cohort`, `module_label`, `mapping`, `n_frozen_genes`, `n_available_genes`, `coverage_fraction`, `minimum_rule_passed`, `available_genes`, `missing_genes`, `validation_tier`, `score`

#### `results/tables/TARGET_OS_primary_frozen_program_validation.csv`

- Rows: 4
- SHA-256: `1bc7d5881bae8f9bd63088c37aad57c6294d47f97d71baa06e7f9aa12fec04fe`
- Columns: `cohort`, `endpoint`, `module_label`, `score_name`, `analysis_tier`, `n`, `events`, `n_covariates`, `covariates`, `score_coef`, `score_hr_per_sd`, `score_ci_low`, `score_ci_high`, `score_p`, `model_c_index`, `fixed_score_c_index`, `fixed_score_c_index_ci_low`, `fixed_score_c_index_ci_high`, `bootstrap_valid`, `ph_test_p`, `error`, `primary_p`, `q_within_endpoint`, `direction_consistent`, `q_global_eight_tests`, `external_support_class`

#### `results/tables/human_external_validation_adjusted_models.csv`

- Rows: 14
- SHA-256: `88e2f4df8a6aea72f899abbfa446f820cf1e79eb6943f3fb6ee73acc824d5e18`
- Columns: `cohort`, `endpoint`, `module_label`, `score_name`, `analysis_tier`, `n`, `events`, `n_covariates`, `covariates`, `score_coef`, `score_hr_per_sd`, `score_ci_low`, `score_ci_high`, `score_p`, `model_c_index`, `fixed_score_c_index`, `fixed_score_c_index_ci_low`, `fixed_score_c_index_ci_high`, `bootstrap_valid`, `ph_test_p`, `error`, `adjustment_type`, `positives`, `negatives`, `auc`, `auc_ci_low`, `auc_ci_high`, `average_precision`, `mann_whitney_p_two_sided`, `logistic_coef`, `logistic_or_per_sd`, `logistic_ci_low`, `logistic_ci_high`, `logistic_p`

#### `results/tables/human_external_validation_cross_cohort_synthesis.csv`

- Rows: 4
- SHA-256: `f108ce1d1d94920e75664d464e08b24072d67aec6bd087c0f665b2af499450ae`
- Columns: `module_label`, `target_os_score_hr_per_sd`, `target_os_score_ci_low`, `target_os_score_ci_high`, `target_os_primary_p`, `target_os_q_within_endpoint`, `target_os_q_global_eight_tests`, `target_os_fixed_score_c_index`, `target_os_fixed_score_c_index_ci_low`, `target_os_fixed_score_c_index_ci_high`, `target_os_external_support_class`, `gse_met_auc`, `gse_met_auc_ci_low`, `gse_met_auc_ci_high`, `gse_met_average_precision`, `gse_met_logistic_or_per_sd`, `gse_met_logistic_ci_low`, `gse_met_logistic_ci_high`, `gse_met_primary_p`, `gse_met_q_within_endpoint`, `gse_met_q_global_eight_tests`, `gse_met_external_support_class`, `gse_os_score_hr_per_sd`, `gse_os_score_ci_low`, `gse_os_score_ci_high`, `gse_os_score_p`, `gse_os_q_within_gse_os_sensitivity`, `gse_os_fixed_score_c_index`, `target_random_percentile`, `target_random_empirical_p`, `gse_random_percentile`, `gse_random_empirical_p`, `direction_consistent_target_and_metastasis`, `cross_cohort_support_summary`

#### `results/tables/human_external_validation_primary_multiplicity.csv`

- Rows: 8
- SHA-256: `875038bf4d5a1857697ee7751454eecaabb3232f04bd6afd6faf19682431f448`
- Columns: `cohort`, `endpoint`, `module_label`, `score_name`, `primary_p`, `q_global_eight_tests`, `direction_consistent`

#### `results/tables/human_external_validation_random_gene_set_controls.csv`

- Rows: 8
- SHA-256: `80246ff94c05c2f3026447e1e7ec99342d6e9ff9795e2aa02e5437502ef5d658`
- Columns: `cohort`, `endpoint`, `module_label`, `observed_metric`, `n_module_genes_available`, `n_random_valid`, `random_mean`, `random_median`, `random_q90`, `random_q95`, `observed_percentile`, `empirical_p_greater_equal`

#### `results/tables/human_external_validation_random_gene_set_distribution.csv`

- Rows: 8000
- SHA-256: `5d8a51ff231f0c864c6bf18ab331e6e42ae94e819b0e1aad6a241a36a3df59d1`
- Columns: `cohort`, `endpoint`, `module_label`, `repeat`, `panel_size`, `random_metric`

#### `results/tables/human_external_validation_score_variant_sensitivity.csv`

- Rows: 36
- SHA-256: `460459d8c8b41c3c62fc528fd020774fe80ce0bf5d71c609b8c9cec82e8849ab`
- Columns: `cohort`, `endpoint`, `module_label`, `score_name`, `analysis_tier`, `n`, `events`, `n_covariates`, `covariates`, `score_coef`, `score_hr_per_sd`, `score_ci_low`, `score_ci_high`, `score_p`, `model_c_index`, `fixed_score_c_index`, `fixed_score_c_index_ci_low`, `fixed_score_c_index_ci_high`, `bootstrap_valid`, `ph_test_p`, `error`, `score_variant`, `positives`, `negatives`, `auc`, `auc_ci_low`, `auc_ci_high`, `average_precision`, `mann_whitney_p_two_sided`, `logistic_coef`, `logistic_or_per_sd`, `logistic_ci_low`, `logistic_ci_high`, `logistic_p`

#### `results/tables/human_external_validation_secondary_programs.csv`

- Rows: 6
- SHA-256: `88149c6b9c8563841db5d287272cdeece39d10733fc0b311c1dbe743d60a3ffb`
- Columns: `cohort`, `endpoint`, `module_label`, `score_name`, `analysis_tier`, `n`, `events`, `n_covariates`, `covariates`, `score_coef`, `score_hr_per_sd`, `score_ci_low`, `score_ci_high`, `score_p`, `model_c_index`, `fixed_score_c_index`, `fixed_score_c_index_ci_low`, `fixed_score_c_index_ci_high`, `bootstrap_valid`, `ph_test_p`, `error`, `positives`, `negatives`, `auc`, `auc_ci_low`, `auc_ci_high`, `average_precision`, `mann_whitney_p_two_sided`, `logistic_coef`, `logistic_or_per_sd`, `logistic_ci_low`, `logistic_ci_high`, `logistic_p`, `p_for_fdr`, `q_within_secondary_endpoint`

---

## `scripts/25_prepare_gse39055_third_human_cohort.py`

- SHA-256: `f5fd842d65adc2afd6828e055bf942efe764e14152aac5ff7f807d1d850a35f9`

### Relevant functions

#### `prepare_gse39055`

- Lines: 238-358
- Signals: `loadings`

```python
def prepare_gse39055() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    GEOparse = import_geoparse()

    print("")
    print("Downloading/loading GSE39055 from NCBI GEO:")
    gse = GEOparse.get_GEO(
        geo=GSE_ACCESSION,
        destdir=str(GSE_RAW_DIR),
        silent=False,
    )

    expression_probe_by_sample = gse.pivot_samples("VALUE")
    if expression_probe_by_sample.empty:
        raise RuntimeError("GSE39055 expression matrix is empty.")
    expression_probe_by_sample.index = expression_probe_by_sample.index.astype(str)

    phenotype = gse.phenotype_data.copy()
    phenotype.index = phenotype.index.astype(str)
    phenotype.to_csv(OUTPUT_PHENOTYPE_RAW)

    if not gse.gpls:
        raise RuntimeError("GSE39055 platform annotation was not loaded.")
    gpl_name = sorted(gse.gpls.keys())[0]
    platform = gse.gpls[gpl_name].table.copy()
    probe_col, symbol_col = detect_platform_columns(platform)

    annotation = platform[[probe_col, symbol_col]].copy()
    annotation[probe_col] = annotation[probe_col].astype(str)
    annotation["gene_symbol"] = annotation[symbol_col].map(
        normalize_unambiguous_gene_symbol
    )
    annotation = annotation[annotation["gene_symbol"].ne("")].copy()
    annotation = annotation.drop_duplicates(probe_col, keep="first")

    common_probes = expression_probe_by_sample.index.intersection(
        annotation[probe_col]
    )
    if len(common_probes) == 0:
        raise RuntimeError("No GSE39055 probes matched GPL14951 annotation.")

    expr = expression_probe_by_sample.loc[common_probes].apply(
        pd.to_numeric, errors="coerce"
    )
    ann = annotation.set_index(probe_col).loc[common_probes]
    probe_variance = expr.var(axis=1)

    probe_map = pd.DataFrame(
        {
            "probe_id": common_probes,
            "gene_symbol": ann["gene_symbol"].values,
            "probe_variance": probe_variance.loc[common_probes].values,
        }
    )
    probe_map = probe_map.sort_values(
        ["gene_symbol", "probe_variance", "probe_id"],
        ascending=[True, False, True],
    )
    selected_map = probe_map.drop_duplicates("gene_symbol", keep="first")
    selected_map.to_csv(OUTPUT_PROBE_MAP, index=False)

    selected_expr = expr.loc[selected_map["probe_id"].tolist()].copy()
    selected_expr.index = selected_map.set_index("probe_id").loc[
        selected_expr.index, "gene_symbol"
    ]
    expression = selected_expr.T
    expression.index.name = "geo_sample_id"
    expression = expression.loc[:, ~expression.columns.duplicated()].copy()
    expression = expression.loc[:, expression.var(axis=0) > 0]

    clinical_rows = []
    for sample_id, gsm in gse.gsms.items():
        characteristics = parse_characteristics(gsm.metadata)
        necrosis_value, good_necrosis = parse_percent_necrosis(
            characteristics.get("percent necrosis", "")
        )
        recurrence_event = binary_yes_no(
            characteristics.get("recurrence", "")
        )
        death_event = binary_yes_no(characteristics.get("death", ""))
        rfs_time = numeric_from_text(
            characteristics.get(
                "time until first recurrence or latest follow-up (months)",
                "",
            )
        )
        age_years = numeric_from_text(characteristics.get("age", ""))

        clinical_rows.append(
            {
                "geo_sample_id": sample_id,
                "title": str(gsm.metadata.get("title", [""])[0]),
                "age_years": age_years,
                "sex": characteristics.get("gender", ""),
                "chemotherapy": characteristics.get("chemotherapy", ""),
                "percent_necrosis_raw": characteristics.get(
                    "percent necrosis", ""
                ),
                "percent_necrosis_numeric": necrosis_value,
                "good_necrosis_response_ge90": good_necrosis,
                "recurrence_event": recurrence_event,
                "death_event_descriptive": death_event,
                "rfs_time_months": rfs_time,
                "tissue": characteristics.get("tissue", ""),
                "biopsy_resection_pair": characteristics.get(
                    "biopsy/resection pair", ""
                ),
                "metadata_text_combined": metadata_as_text(gsm.metadata),
                "endpoint_note": (
                    "RFS time is the GEO field 'time until first recurrence or "
                    "latest follow-up (months)'; event is recurrence Y/N. "
                    "Death status is descriptive because no death time is supplied."
                ),
            }
        )

    clinical = pd.DataFrame(clinical_rows).set_index("geo_sample_id")
    common_samples = expression.index.intersection(clinical.index)
    expression = expression.loc[common_samples].copy()
    clinical = clinical.loc[common_samples].copy()

    return expression, clinical, selected_map
```

#### `zscore_columns`

- Lines: 361-367
- Signals: `mean, std`

```python
def zscore_columns(expression: pd.DataFrame) -> pd.DataFrame:
    x = expression.apply(pd.to_numeric, errors="coerce")
    x = x.replace([np.inf, -np.inf], np.nan)
    x = x.fillna(x.median(axis=0))
    stds = x.std(axis=0).replace(0, np.nan)
    z = (x - x.mean(axis=0)) / stds
    return z.loc[:, z.notna().all(axis=0)]
```

#### `zscore_series`

- Lines: 370-375
- Signals: `mean, std`

```python
def zscore_series(values: pd.Series) -> pd.Series:
    values = pd.to_numeric(values, errors="coerce")
    std = values.std()
    if not np.isfinite(std) or std == 0:
        return pd.Series(np.nan, index=values.index)
    return (values - values.mean()) / std
```

#### `safe_corr`

- Lines: 378-382
- Signals: `std, correlation`

```python
def safe_corr(a: pd.Series, b: pd.Series) -> float:
    frame = pd.concat([a.rename("a"), b.rename("b")], axis=1).dropna()
    if frame.shape[0] < 5 or frame["a"].std() == 0 or frame["b"].std() == 0:
        return np.nan
    return float(frame["a"].corr(frame["b"]))
```

#### `human_cohort_pc1`

- Lines: 385-397
- Signals: `correlation`

```python
def human_cohort_pc1(
    z: pd.DataFrame,
    reference: pd.Series | None = None,
) -> pd.Series:
    if z.shape[1] < 2:
        return pd.Series(np.nan, index=z.index)
    pca = PCA(n_components=1, random_state=RANDOM_SEED)
    score = pd.Series(pca.fit_transform(z).ravel(), index=z.index)
    if reference is not None:
        corr = safe_corr(score, reference)
        if np.isfinite(corr) and corr < 0:
            score = -score
    return zscore_series(score)
```

#### `compute_module_scores`

- Lines: 400-484
- Signals: `mean, weights, direction, coverage`

```python
def compute_module_scores(
    expression: pd.DataFrame,
    strict_weights: pd.DataFrame,
    broad_weights: pd.DataFrame,
    manifest: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    expression = expression.copy()
    expression.columns = expression.columns.astype(str).str.upper()
    expression = expression.loc[:, ~expression.columns.duplicated()].copy()

    scores = pd.DataFrame(index=expression.index)
    coverage_rows = []

    for mapping_name, weights in [
        ("strict", strict_weights),
        ("broad", broad_weights),
    ]:
        for module_label, part in weights.groupby("module_label"):
            part = part.copy()
            part["human_gene_symbol"] = (
                part["human_gene_symbol"].astype(str).str.upper()
            )
            part = part.drop_duplicates("human_gene_symbol", keep="first")
            requested = part["human_gene_symbol"].tolist()
            available = [gene for gene in requested if gene in expression.columns]

            n_requested = len(requested)
            n_available = len(available)
            fraction = n_available / n_requested if n_requested else 0.0
            passed = (
                n_available >= MIN_SCORE_GENES
                and fraction >= MIN_SCORE_FRACTION
            )

            coverage_rows.append(
                {
                    "cohort": GSE_ACCESSION,
                    "module_label": module_label,
                    "mapping": mapping_name,
                    "n_frozen_genes": n_requested,
                    "n_available_genes": n_available,
                    "coverage_fraction": fraction,
                    "minimum_rule_passed": passed,
                    "available_genes": ";".join(available),
                    "missing_genes": ";".join(
                        gene for gene in requested if gene not in available
                    ),
                }
            )
            if not passed:
                continue

            z = zscore_columns(expression[available])
            available = list(z.columns)
            if len(available) < MIN_SCORE_GENES:
                continue

            weight_indexed = part.set_index("human_gene_symbol").loc[available]
            raw_loadings = pd.to_numeric(
                weight_indexed["risk_oriented_loading"],
                errors="coerce",
            ).fillna(0.0)
            signs = np.sign(raw_loadings).replace(0, 1)

            signed_mean = z.mul(signs, axis=1).mean(axis=1)
            signed_mean = zscore_series(signed_mean)

            if raw_loadings.abs().sum() > 0:
                normalized = raw_loadings / raw_loadings.abs().sum()
                weighted = z.mul(normalized, axis=1).sum(axis=1)
                weighted = zscore_series(weighted)
            else:
                weighted = pd.Series(np.nan, index=z.index)

            pc1 = human_cohort_pc1(z, reference=signed_mean)
            prefix = f"{module_label}__{mapping_name}"
            scores[f"{prefix}__signed_mean_z"] = signed_mean
            scores[f"{prefix}__canine_pca_weighted_z"] = weighted
            scores[f"{prefix}__human_pc1_z"] = pc1

    tier_map = manifest.set_index("module_label")["validation_tier"].to_dict()
    coverage = pd.DataFrame(coverage_rows)
    coverage["validation_tier"] = coverage["module_label"].map(tier_map)
    scores.insert(0, "cohort", GSE_ACCESSION)
    return scores, coverage
```

#### `compute_proliferation_scores`

- Lines: 487-581
- Signals: `mean`

```python
def compute_proliferation_scores(
    expression: pd.DataFrame,
    proliferation_mapping: pd.DataFrame,
    strict_weights: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    expression = expression.copy()
    expression.columns = expression.columns.astype(str).str.upper()
    expression = expression.loc[:, ~expression.columns.duplicated()].copy()

    frozen_genes = (
        proliferation_mapping["human_gene_symbol"]
        .dropna()
        .astype(str)
        .str.upper()
        .drop_duplicates()
        .tolist()
    )
    available = [gene for gene in frozen_genes if gene in expression.columns]
    fraction = len(available) / len(frozen_genes) if frozen_genes else 0.0
    passed = (
        len(available) >= MIN_PROLIFERATION_GENES
        and fraction >= MIN_PROLIFERATION_FRACTION
    )

    coverage_rows = [
        {
            "cohort": GSE_ACCESSION,
            "score": "strict_human_meta_proliferation_pc1",
            "n_frozen_genes": len(frozen_genes),
            "n_available_genes": len(available),
            "coverage_fraction": fraction,
            "minimum_rule_passed": passed,
        }
    ]

    scores = pd.DataFrame(index=expression.index)
    if not passed:
        return scores, pd.DataFrame(coverage_rows)

    z = zscore_columns(expression[available])
    anchors = [gene for gene in PROLIFERATION_ANCHOR_SYMBOLS if gene in z.columns]
    reference = z[anchors].mean(axis=1) if len(anchors) >= 3 else z.mean(axis=1)
    proliferation_pc1 = human_cohort_pc1(z, reference=reference)
    scores["strict_human_meta_proliferation_pc1_z"] = proliferation_pc1

    m40_genes = (
        strict_weights[strict_weights["module_label"].eq("M40")]
        ["human_gene_symbol"]
        .dropna()
        .astype(str)
        .str.upper()
        .drop_duplicates()
        .tolist()
    )
    disjoint_frozen = [
        gene for gene in frozen_genes if gene not in set(m40_genes)
    ]
    disjoint_available = [
        gene for gene in disjoint_frozen if gene in expression.columns
    ]
    disjoint_fraction = (
        len(disjoint_available) / len(disjoint_frozen)
        if disjoint_frozen else 0.0
    )
    disjoint_passed = (
        len(disjoint_available) >= MIN_PROLIFERATION_GENES
        and disjoint_fraction >= MIN_PROLIFERATION_FRACTION
    )
    coverage_rows.append(
        {
            "cohort": GSE_ACCESSION,
            "score": "M40_disjoint_strict_human_meta_proliferation_pc1",
            "n_frozen_genes": len(disjoint_frozen),
            "n_available_genes": len(disjoint_available),
            "coverage_fraction": disjoint_fraction,
            "minimum_rule_passed": disjoint_passed,
        }
    )

    if disjoint_passed:
        z_disjoint = zscore_columns(expression[disjoint_available])
        anchors_disjoint = [
            gene for gene in PROLIFERATION_ANCHOR_SYMBOLS
            if gene in z_disjoint.columns
        ]
        reference_disjoint = (
            z_disjoint[anchors_disjoint].mean(axis=1)
            if len(anchors_disjoint) >= 3
            else z_disjoint.mean(axis=1)
        )
        scores["M40_disjoint_strict_human_meta_proliferation_pc1_z"] = (
            human_cohort_pc1(z_disjoint, reference=reference_disjoint)
        )

    return scores, pd.DataFrame(coverage_rows)
```

#### `residualize_outcome_blind`

- Lines: 584-605
- Signals: `std, matrix_product`

```python
def residualize_outcome_blind(
    score: pd.Series,
    covariate: pd.Series,
) -> pd.Series:
    frame = pd.concat(
        [score.rename("score"), covariate.rename("covariate")],
        axis=1,
    ).dropna()
    residual = pd.Series(np.nan, index=score.index)
    if frame.shape[0] < 10 or frame["covariate"].std() == 0:
        return residual

    x = np.column_stack(
        [np.ones(frame.shape[0]), frame["covariate"].values]
    )
    beta, _, _, _ = np.linalg.lstsq(
        x,
        frame["score"].values,
        rcond=None,
    )
    residual.loc[frame.index] = frame["score"].values - x @ beta
    return zscore_series(residual)
```

#### `merge_score_components`

- Lines: 608-638
- Signals: ``

```python
def merge_score_components(
    module_scores: pd.DataFrame,
    proliferation_scores: pd.DataFrame,
) -> pd.DataFrame:
    cohort = module_scores["cohort"]
    out = module_scores.drop(columns=["cohort"]).join(
        proliferation_scores,
        how="outer",
    )

    proliferation_col = (
        "M40_disjoint_strict_human_meta_proliferation_pc1_z"
    )
    if proliferation_col in out.columns:
        for module_col in [
            "M40__strict__signed_mean_z",
            "M40__strict__canine_pca_weighted_z",
        ]:
            if module_col not in out.columns:
                continue
            residual_col = module_col.replace(
                "_z",
                "__residual_to_disjoint_proliferation_z",
            )
            out[residual_col] = residualize_outcome_blind(
                out[module_col],
                out[proliferation_col],
            )

    out.insert(0, "cohort", cohort.reindex(out.index))
    return out
```

#### `create_manifest`

- Lines: 641-690
- Signals: `weights, direction`

```python
def create_manifest(
    freeze: dict[str, Any],
    expression: pd.DataFrame,
    clinical: pd.DataFrame,
    scores: pd.DataFrame,
) -> None:
    output_paths = [
        OUTPUT_EXPRESSION,
        OUTPUT_CLINICAL,
        OUTPUT_SCORES,
        OUTPUT_COVERAGE,
        OUTPUT_PROBE_MAP,
        OUTPUT_PHENOTYPE_RAW,
        OUTPUT_PREPARATION_SUMMARY,
        OUTPUT_README,
    ]
    manifest = {
        "script_version": SCRIPT_VERSION,
        "created_utc": utc_now_iso(),
        "frozen_program_freeze_sha256": sha256_file(FREEZE_JSON_FILE),
        "frozen_program_definition": freeze,
        "source": {
            "accession": GSE_ACCESSION,
            "source": "NCBI GEO via GEOparse",
            "platform": "GPL14951",
            "sample_type": "diagnostic osteosarcoma biopsy",
        },
        "cohort_dimensions": {
            "expression": list(expression.shape),
            "clinical": list(clinical.shape),
            "scores": list(scores.shape),
        },
        "outcome_guardrail": (
            "Recurrence and follow-up fields were parsed after expression "
            "harmonization rules were fixed. Outcomes were not used to select "
            "probes, genes, weights, score direction, or validation tier."
        ),
        "files": {},
    }
    for path in output_paths:
        if path.exists():
            manifest["files"][path.name] = {
                "sha256": sha256_file(path),
                "size_bytes": path.stat().st_size,
            }

    OUTPUT_PREPARATION_MANIFEST.write_text(
        json.dumps(manifest, indent=2),
        encoding="utf-8",
    )
```

#### `main`

- Lines: 719-891
- Signals: `weights, direction, coverage`

```python
def main() -> None:
    print("=" * 80)
    print("Prepare GSE39055 third human osteosarcoma cohort")
    print("=" * 80)
    print(f"Script version: {SCRIPT_VERSION}")
    print(f"Project root: {PROJECT_ROOT}")
    print(f"Human processed directory: {HUMAN_PROCESSED_DIR}")
    print("")
    print("Design:")
    print("  Verify frozen canine transfer assets and hashes.")
    print("  Download GSE39055 expression and phenotype metadata from NCBI GEO.")
    print("  Collapse probes to human gene symbols using outcome-blind variance rules.")
    print("  Parse recurrence-free survival metadata.")
    print("  Construct frozen module scores without fitting any outcome model.")
    print("")

    freeze = verify_frozen_inputs()
    manifest = read_required_csv(FROZEN_MANIFEST_FILE)
    strict_weights = read_required_csv(STRICT_WEIGHTS_FILE)
    broad_weights = read_required_csv(BROAD_WEIGHTS_FILE)
    proliferation_mapping = read_required_csv(PROLIFERATION_MAPPING_FILE)

    strict_weights["human_gene_symbol"] = (
        strict_weights["human_gene_symbol"].astype(str).str.upper()
    )
    broad_weights["human_gene_symbol"] = (
        broad_weights["human_gene_symbol"].astype(str).str.upper()
    )
    proliferation_mapping["human_gene_symbol"] = (
        proliferation_mapping["human_gene_symbol"].astype(str).str.upper()
    )

    expression, clinical, _ = prepare_gse39055()

    module_scores, module_coverage = compute_module_scores(
        expression=expression,
        strict_weights=strict_weights,
        broad_weights=broad_weights,
        manifest=manifest,
    )
    proliferation_scores, proliferation_coverage = compute_proliferation_scores(
        expression=expression,
        proliferation_mapping=proliferation_mapping,
        strict_weights=strict_weights,
    )
    scores = merge_score_components(
        module_scores=module_scores,
        proliferation_scores=proliferation_scores,
    )

    coverage = pd.concat(
        [module_coverage, proliferation_coverage],
        axis=0,
        ignore_index=True,
        sort=False,
    )

    expression.to_csv(OUTPUT_EXPRESSION)
    clinical.to_csv(OUTPUT_CLINICAL)
    scores.to_csv(OUTPUT_SCORES)
    coverage.to_csv(OUTPUT_COVERAGE, index=False)

    rfs_complete = clinical[["rfs_time_months", "recurrence_event"]].dropna()
    summary = pd.DataFrame(
        [
            {
                "cohort": GSE_ACCESSION,
                "n_expression_samples": expression.shape[0],
                "n_expression_genes": expression.shape[1],
                "n_clinical_rows": clinical.shape[0],
                "n_rfs_complete": rfs_complete.shape[0],
                "n_recurrence_events": int(rfs_complete["recurrence_event"].sum()),
                "n_censored_without_recurrence": int(
                    (rfs_complete["recurrence_event"] == 0).sum()
                ),
                "n_frozen_score_columns": (
                    scores.shape[1] - int("cohort" in scores.columns)
                ),
            }
        ]
    )
    summary.to_csv(OUTPUT_PREPARATION_SUMMARY, index=False)

    write_readme()
    create_manifest(
        freeze=freeze,
        expression=expression,
        clinical=clinical,
        scores=scores,
    )

    print("")
    print("=" * 80)
    print("GSE39055 preparation summary")
    print("=" * 80)
    print(summary.to_string(index=False))

    print("")
    print("=" * 80)
    print("Frozen score coverage: GSE39055")
    print("=" * 80)
    display_cols = [
        "module_label",
        "mapping",
        "validation_tier",
        "n_frozen_genes",
        "n_available_genes",
        "coverage_fraction",
        "minimum_rule_passed",
    ]
    display_cols = [
        column for column in display_cols if column in module_coverage.columns
    ]
    print(
        module_coverage[display_cols]
        .sort_values(["validation_tier", "module_label", "mapping"])
        .to_string(index=False)
    )

    print("")
    print("=" * 80)
    print("Endpoint preparation audit")
    print("=" * 80)
    print("Recurrence event counts:")
    print(clinical["recurrence_event"].value_counts(dropna=False).to_string())
    print("")
    print("RFS fields:")
    print(
        clinical[["rfs_time_months", "recurrence_event"]]
        .agg(["count", "min", "median", "max"])
        .to_string()
    )
    print("")
    print("Clinical covariate availability:")
    print(
        clinical[
            [
                "age_years",
                "sex",
                "percent_necrosis_numeric",
                "good_necrosis_response_ge90",
                "death_event_descriptive",
            ]
        ]
        .agg(["count"])
        .to_string()
    )

    print("")
    print("=" * 80)
    print("Interpretation guardrails")
    print("=" * 80)
    print("No GSE39055 outcome was used to select probes, genes, weights, score direction, or validation tier.")
    print("Probe collapsing used the highest-variance probe per unambiguous gene symbol.")
    print("The primary future endpoint is recurrence-free survival using recurrence Y/N and the provided recurrence/follow-up time.")
    print("Death status is descriptive because no separate time-to-death field is provided.")
    print("Outcome association testing and multiplicity control are deferred to script 26.")

    print("")
    print("Saved:")
    for path in [
        OUTPUT_EXPRESSION,
        OUTPUT_CLINICAL,
        OUTPUT_SCORES,
        OUTPUT_COVERAGE,
        OUTPUT_PROBE_MAP,
        OUTPUT_PHENOTYPE_RAW,
        OUTPUT_PREPARATION_SUMMARY,
        OUTPUT_PREPARATION_MANIFEST,
        OUTPUT_README,
    ]:
        print(path)
    print("Done.")
```

### Relevant standalone lines

| Line | Signals | Code |
|---:|---|---|
| 242 | loadings | `print("Downloading/loading GSE39055 from NCBI GEO:")` |
| 365 | std | `stds = x.std(axis=0).replace(0, np.nan)` |
| 366 | mean | `z = (x - x.mean(axis=0)) / stds` |
| 372 | std | `std = values.std()` |
| 375 | mean | `return (values - values.mean()) / std` |
| 380 | std | `if frame.shape[0] < 5 or frame["a"].std() == 0 or frame["b"].std() == 0:` |
| 382 | correlation | `return float(frame["a"].corr(frame["b"]))` |
| 394 | correlation | `corr = safe_corr(score, reference)` |
| 395 | correlation | `if np.isfinite(corr) and corr < 0:` |
| 413 | weights | `for mapping_name, weights in [` |
| 417 | weights | `for module_label, part in weights.groupby("module_label"):` |
| 462 | direction | `signs = np.sign(raw_loadings).replace(0, 1)` |
| 464 | mean | `signed_mean = z.mul(signs, axis=1).mean(axis=1)` |
| 469 | weights | `weighted = z.mul(normalized, axis=1).sum(axis=1)` |
| 470 | weights | `weighted = zscore_series(weighted)` |
| 472 | weights | `weighted = pd.Series(np.nan, index=z.index)` |
| 477 | weights | `scores[f"{prefix}__canine_pca_weighted_z"] = weighted` |
| 481 | coverage | `coverage = pd.DataFrame(coverage_rows)` |
| 482 | coverage | `coverage["validation_tier"] = coverage["module_label"].map(tier_map)` |
| 484 | coverage | `return scores, coverage` |
| 528 | mean | `reference = z[anchors].mean(axis=1) if len(anchors) >= 3 else z.mean(axis=1)` |
| 573 | mean | `z_disjoint[anchors_disjoint].mean(axis=1)` |
| 575 | mean | `else z_disjoint.mean(axis=1)` |
| 593 | std | `if frame.shape[0] < 10 or frame["covariate"].std() == 0:` |
| 604 | matrix_product | `residual.loc[frame.index] = frame["score"].values - x @ beta` |
| 676 | weights, direction | `"probes, genes, weights, score direction, or validation tier."` |
| 769 | coverage | `coverage = pd.concat(` |
| 779 | coverage | `coverage.to_csv(OUTPUT_COVERAGE, index=False)` |
| 818 | coverage | `print("Frozen score coverage: GSE39055")` |
| 871 | weights, direction | `print("No GSE39055 outcome was used to select probes, genes, weights, score direction, or validation tier.")` |

### Referenced result-table schemas

#### `results/tables/GSE238110_frozen_canine_transfer_program_manifest.csv`

- Rows: 13
- SHA-256: `fcd34b985ea99986f39dfc43ce8e4018db2fde9ecc8a225bb40ac93006d5372c`
- Columns: `module_label`, `validation_tier`, `multiplicity_family`, `provisional_program_label`, `program_label_requires_enrichment_confirmation`, `canine_primary_endpoint`, `canine_secondary_endpoint`, `positive_score_interpretation`, `n_canine_genes_used_for_pca`, `canine_pc1_explained_variance`, `canine_pca_orientation_correlation`, `dfi_full_cohort_coef`, `os_full_cohort_coef`, `risk_orientation_multiplier`, `n_strict_human_genes`, `n_broad_human_genes`, `strict_transfer_eligible`, `module_transfer_qc_tier`, `transfer_priority_score`, `fraction_strict_symbol_concordant`, `fraction_broad_transferable`, `raw_module_proliferation_correlation`, `orthogonal_variance_fraction_1_minus_r2`, `n_overlap_symbols_with_proliferation`, `primary_human_score`, `secondary_human_score`, `sensitivity_human_scores`, `frozen_after_canine_script`, `dfi_module_only_mean_c_index`, `dfi_module_only_std_c_index`, `dfi_module_only_fraction_above_0_50`, `dfi_module_plus_disjoint_proliferation_mean_c_index`, `dfi_module_plus_disjoint_proliferation_std_c_index`, `dfi_module_plus_disjoint_proliferation_fraction_above_0_50`, `dfi_residual_to_disjoint_proliferation_mean_c_index`, `dfi_residual_to_disjoint_proliferation_std_c_index`, `dfi_residual_to_disjoint_proliferation_fraction_above_0_50`, `dfi_residual_to_disjoint_proliferation_and_weight_mean_c_index`, `dfi_residual_to_disjoint_proliferation_and_weight_std_c_index`, `dfi_residual_to_disjoint_proliferation_and_weight_fraction_above_0_50`, `os_module_only_mean_c_index`, `os_module_only_std_c_index`, `os_module_only_fraction_above_0_50`, `os_module_plus_disjoint_proliferation_mean_c_index`, `os_module_plus_disjoint_proliferation_std_c_index`, `os_module_plus_disjoint_proliferation_fraction_above_0_50`, `os_residual_to_disjoint_proliferation_mean_c_index`, `os_residual_to_disjoint_proliferation_std_c_index`, `os_residual_to_disjoint_proliferation_fraction_above_0_50`, `os_residual_to_disjoint_proliferation_and_weight_mean_c_index`, `os_residual_to_disjoint_proliferation_and_weight_std_c_index`, `os_residual_to_disjoint_proliferation_and_weight_fraction_above_0_50`, `script20_dfi_recommended_role`, `script20_os_recommended_role`, `manual_freeze_reason`

#### `results/tables/GSE238110_frozen_transfer_gene_weights_broad.csv`

- Rows: 389
- SHA-256: `ee5f661677b173507e04856f848b6757ede2750dd69b9bb0d88a9f6a2a77bd4c`
- Columns: `module_label`, `canine_gene`, `canine_gene_symbol`, `human_gene_symbol`, `ortholog_qc_status`, `raw_pca_loading`, `risk_oriented_loading`, `absolute_risk_loading`, `is_strict_mapping`, `is_broad_mapping`, `human_symbol_duplicate_count`, `normalized_abs_sum_weight`

#### `results/tables/GSE238110_frozen_transfer_gene_weights_strict.csv`

- Rows: 321
- SHA-256: `4f065aa4c4edf117a0c74015840d2b4b2347929f172cd517e1818ba0f6163b91`
- Columns: `module_label`, `canine_gene`, `canine_gene_symbol`, `human_gene_symbol`, `ortholog_qc_status`, `raw_pca_loading`, `risk_oriented_loading`, `absolute_risk_loading`, `is_strict_mapping`, `is_broad_mapping`, `human_symbol_duplicate_count`, `normalized_abs_sum_weight`

#### `results/tables/GSE238110_frozen_transfer_scoring_specification.csv`

- Rows: 40
- SHA-256: `562aab92caa949691f43b0285cf2585c98008e78c62d2c15477a6d156d7f47a7`
- Columns: `module_label`, `validation_tier`, `score_name`, `analysis_role`, `gene_mapping`, `within_cohort_preprocessing`, `score_formula`, `minimum_gene_rule`, `outcome_use`

#### `results/tables/GSE39055_GEO_phenotype_raw.csv`

- Rows: 37
- SHA-256: `3e1e8714429eed0e3572b7b76e39f9b5353b05074f3ebf297d67f30a3ab40071`
- Columns: ``, `title`, `geo_accession`, `status`, `submission_date`, `last_update_date`, `type`, `channel_count`, `source_name_ch1`, `organism_ch1`, `taxid_ch1`, `characteristics_ch1.0.age`, `characteristics_ch1.1.gender`, `characteristics_ch1.2.chemotherapy`, `characteristics_ch1.3.percent necrosis`, `characteristics_ch1.4.recurrence`, `characteristics_ch1.5.death`, `characteristics_ch1.6.time until first recurrence or latest follow-up (months)`, `characteristics_ch1.7.tissue`, `characteristics_ch1.8.biopsy/resection pair`, `molecule_ch1`, `extract_protocol_ch1`, `label_ch1`, `label_protocol_ch1`, `hyb_protocol`, `scan_protocol`, `data_processing`, `platform_id`, `contact_name`, `contact_email`, `contact_laboratory`, `contact_department`, `contact_institute`, `contact_address`, `contact_city`, `contact_state`, `contact_zip/postal_code`, `contact_country`, `supplementary_file`, `series_id`, `data_row_count`, `relation`

#### `results/tables/GSE39055_frozen_transfer_score_coverage.csv`

- Rows: 28
- SHA-256: `1945407e793a0c3d3213a55bc0e8cf7526dbfb8b9d65e716818b5f5748ec50e7`
- Columns: `cohort`, `module_label`, `mapping`, `n_frozen_genes`, `n_available_genes`, `coverage_fraction`, `minimum_rule_passed`, `available_genes`, `missing_genes`, `validation_tier`, `score`

#### `results/tables/GSE39055_preparation_summary.csv`

- Rows: 1
- SHA-256: `0a10c02a2f458982a0734c7e884f4b6042975fc484108f47d71aa8b97e597f95`
- Columns: `cohort`, `n_expression_samples`, `n_expression_genes`, `n_clinical_rows`, `n_rfs_complete`, `n_recurrence_events`, `n_censored_without_recurrence`, `n_frozen_score_columns`

#### `results/tables/GSE39055_probe_to_gene_symbol_selected.csv`

- Rows: 20793
- SHA-256: `00a773be24c545aaed9f22fa2e35594e946925026e5cb37c64a572b53288df99`
- Columns: `probe_id`, `gene_symbol`, `probe_variance`

#### `results/tables/frozen_strict_human_proliferation_mapping.csv`

- Rows: 111
- SHA-256: `72e539d40b3eeb604986e224a345f5218a42f5fcd70ea48b27d000c178413de8`
- Columns: `gene`, `canine_proliferation_gene`, `dog_symbol_key`, `human_gene_symbol`, `ortholog_qc_status`

---

## `scripts/26_validate_gse39055_rfs.py`

- SHA-256: `774bbf4dffb6cfc7d27066917682a3d32f5a8091c7f1f95c3c9b30a5bde5a34f`

### Relevant functions

#### `zscore_series`

- Lines: 129-134
- Signals: `mean, std`

```python
def zscore_series(values: pd.Series) -> pd.Series:
    values = pd.to_numeric(values, errors="coerce")
    std = values.std()
    if not np.isfinite(std) or std == 0:
        return pd.Series(np.nan, index=values.index, dtype=float)
    return (values - values.mean()) / std
```

#### `fit_cox`

- Lines: 165-232
- Signals: `std, rank`

```python
def fit_cox(
    data: pd.DataFrame,
    score_col: str,
    covariates: list[str] | None = None,
    time_col: str = TIME_COL,
    event_col: str = EVENT_COL,
    penalizer: float = COX_PENALIZER,
) -> dict[str, Any]:
    covariates = covariates or []
    columns = [time_col, event_col, score_col] + covariates
    frame = data[columns].replace([np.inf, -np.inf], np.nan).dropna().copy()

    result: dict[str, Any] = {
        "n": frame.shape[0],
        "events": int(frame[event_col].sum()) if frame.shape[0] else 0,
        "hr_per_sd": np.nan,
        "ci_low": np.nan,
        "ci_high": np.nan,
        "p": np.nan,
        "coef": np.nan,
        "se_coef": np.nan,
        "c_index": np.nan,
        "ph_test_p": np.nan,
        "error": "",
    }

    if frame.shape[0] < 20:
        result["error"] = "too_few_samples"
        return result
    if frame[event_col].sum() < 5:
        result["error"] = "too_few_events"
        return result
    if frame[score_col].std() == 0:
        result["error"] = "zero_variance_score"
        return result

    model = CoxPHFitter(penalizer=penalizer)
    try:
        model.fit(
            frame,
            duration_col=time_col,
            event_col=event_col,
            fit_options={"max_steps": 500},
        )
        summary = model.summary.loc[score_col]
        result.update(
            {
                "hr_per_sd": float(summary["exp(coef)"]),
                "ci_low": float(summary["exp(coef) lower 95%"]),
                "ci_high": float(summary["exp(coef) upper 95%"]),
                "p": float(summary["p"]),
                "coef": float(summary["coef"]),
                "se_coef": float(summary["se(coef)"]),
                "c_index": float(model.concordance_index_),
            }
        )
        try:
            ph = proportional_hazard_test(
                model,
                frame,
                time_transform="rank",
            )
            result["ph_test_p"] = float(ph.summary.loc[score_col, "p"])
        except Exception:
            pass
    except Exception as exc:
        result["error"] = str(exc)[:500]
    return result
```

#### `bootstrap_probability_coef_positive`

- Lines: 313-344
- Signals: `mean`

```python
def bootstrap_probability_coef_positive(
    data: pd.DataFrame,
    score_col: str,
    n_bootstrap: int,
    seed: int,
) -> tuple[float, float, float, int]:
    frame = data[[TIME_COL, EVENT_COL, score_col]].dropna().copy()
    event_idx = np.where(frame[EVENT_COL].values == 1)[0]
    censor_idx = np.where(frame[EVENT_COL].values == 0)[0]
    rng = np.random.default_rng(seed)
    coefs: list[float] = []

    for _ in range(n_bootstrap):
        sampled_event = rng.choice(event_idx, size=len(event_idx), replace=True)
        sampled_censor = rng.choice(censor_idx, size=len(censor_idx), replace=True)
        sampled = np.concatenate([sampled_event, sampled_censor])
        rng.shuffle(sampled)
        part = frame.iloc[sampled].copy()
        fit = fit_cox(part, score_col)
        coef = fit["coef"]
        if np.isfinite(coef):
            coefs.append(float(coef))

    if not coefs:
        return np.nan, np.nan, np.nan, 0
    coef_values = np.asarray(coefs, dtype=float)
    return (
        float(np.exp(np.quantile(coef_values, 0.025))),
        float(np.exp(np.quantile(coef_values, 0.975))),
        float(np.mean(coef_values > 0)),
        len(coef_values),
    )
```

#### `score_variant_sensitivity`

- Lines: 583-611
- Signals: `coverage`

```python
def score_variant_sensitivity(
    frame: pd.DataFrame,
    coverage: pd.DataFrame,
) -> pd.DataFrame:
    rows = []
    variant_suffixes = [
        ("strict_signed_mean", "__strict__signed_mean_z"),
        ("strict_canine_pca_weighted", "__strict__canine_pca_weighted_z"),
        ("strict_human_pc1", "__strict__human_pc1_z"),
        ("broad_signed_mean", "__broad__signed_mean_z"),
        ("broad_canine_pca_weighted", "__broad__canine_pca_weighted_z"),
        ("broad_human_pc1", "__broad__human_pc1_z"),
    ]

    for module in PRIMARY_MODULES + SECONDARY_MODULES:
        for variant_label, suffix in variant_suffixes:
            score_col = f"{module}{suffix}"
            if score_col not in frame.columns:
                continue
            fit = fit_cox(frame, score_col)
            rows.append(
                {
                    "module_label": module,
                    "variant": variant_label,
                    "score_column": score_col,
                    **fit,
                }
            )
    return pd.DataFrame(rows)
```

#### `leave_one_out_stability`

- Lines: 614-653
- Signals: `mean`

```python
def leave_one_out_stability(frame: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for module in PRIMARY_MODULES:
        score_col = f"{module}{PRIMARY_SCORE_SUFFIX}"
        if score_col not in frame.columns:
            continue
        values = []
        for sample_id in frame.index:
            part = frame.drop(index=sample_id)
            fit = fit_cox(part, score_col)
            values.append(
                {
                    "sample_removed": sample_id,
                    "hr": fit["hr_per_sd"],
                    "p": fit["p"],
                    "c_index": safe_fixed_direction_c_index(
                        part[TIME_COL],
                        part[EVENT_COL],
                        part[score_col],
                    ),
                }
            )
        table = pd.DataFrame(values)
        rows.append(
            {
                "module_label": module,
                "n_loo_fits": table.shape[0],
                "hr_min": table["hr"].min(),
                "hr_max": table["hr"].max(),
                "hr_median": table["hr"].median(),
                "fraction_hr_above_1": float((table["hr"] > 1).mean()),
                "c_index_min": table["c_index"].min(),
                "c_index_max": table["c_index"].max(),
                "c_index_median": table["c_index"].median(),
                "fraction_p_below_0_05": float(
                    (table["p"] < 0.05).mean()
                ),
            }
        )
    return pd.DataFrame(rows)
```

#### `expression_bins`

- Lines: 656-679
- Signals: `mean, std`

```python
def expression_bins(expression: pd.DataFrame) -> pd.DataFrame:
    stats = pd.DataFrame(
        {
            "mean": expression.mean(axis=0),
            "sd": expression.std(axis=0),
        }
    ).replace([np.inf, -np.inf], np.nan).dropna()
    stats = stats[stats["sd"] > 0].copy()
    stats["mean_bin"] = pd.qcut(
        stats["mean"],
        q=10,
        labels=False,
        duplicates="drop",
    )
    stats["sd_bin"] = pd.qcut(
        stats["sd"],
        q=10,
        labels=False,
        duplicates="drop",
    )
    stats["bin_key"] = (
        stats["mean_bin"].astype(str) + "_" + stats["sd_bin"].astype(str)
    )
    return stats
```

#### `random_gene_set_controls`

- Lines: 719-838
- Signals: `mean, std, weights, loadings, direction`

```python
def random_gene_set_controls(
    frame: pd.DataFrame,
    expression: pd.DataFrame,
    strict_weights: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    expression = expression.copy()
    expression.columns = expression.columns.astype(str).str.upper()
    expression = expression.loc[frame.index]
    expression = expression.loc[:, ~expression.columns.duplicated()]
    expression = expression.apply(pd.to_numeric, errors="coerce")
    expression = expression.fillna(expression.median(axis=0))
    expression = expression.loc[:, expression.std(axis=0) > 0]
    z = (expression - expression.mean(axis=0)) / expression.std(axis=0)

    stats = expression_bins(expression)
    rng = np.random.default_rng(RANDOM_SEED)
    summary_rows = []
    distribution_rows = []

    for module_index, module in enumerate(PRIMARY_MODULES):
        score_col = f"{module}{PRIMARY_SCORE_SUFFIX}"
        if score_col not in frame.columns:
            continue

        weights = strict_weights[
            strict_weights["module_label"].eq(module)
        ].copy()
        weights["human_gene_symbol"] = (
            weights["human_gene_symbol"].astype(str).str.upper()
        )
        weights = weights.drop_duplicates("human_gene_symbol", keep="first")
        target_genes = [
            gene
            for gene in weights["human_gene_symbol"].tolist()
            if gene in z.columns
        ]
        loadings = (
            weights.set_index("human_gene_symbol")
            .reindex(target_genes)["risk_oriented_loading"]
            .fillna(0.0)
        )
        signs = np.sign(loadings.values)
        signs[signs == 0] = 1

        observed = safe_fixed_direction_c_index(
            frame[TIME_COL],
            frame[EVENT_COL],
            frame[score_col],
        )
        excluded = set(target_genes)
        random_values: list[float] = []

        for repeat in range(1, N_RANDOM_SETS + 1):
            genes = draw_expression_matched_gene_set(
                target_genes=target_genes,
                stats=stats,
                excluded=excluded,
                rng=rng,
            )
            if len(genes) != len(target_genes):
                continue
            permuted_signs = rng.permutation(signs)
            random_score = z[genes].mul(permuted_signs, axis=1).mean(axis=1)
            value = safe_fixed_direction_c_index(
                frame[TIME_COL],
                frame[EVENT_COL],
                random_score,
            )
            if np.isfinite(value):
                random_values.append(value)
                distribution_rows.append(
                    {
                        "module_label": module,
                        "repeat": repeat,
                        "random_c_index": value,
                        "n_genes": len(genes),
                    }
                )

        random_array = np.asarray(random_values, dtype=float)
        empirical_p = (
            (1.0 + np.sum(random_array >= observed))
            / (len(random_array) + 1.0)
            if len(random_array)
            else np.nan
        )
        summary_rows.append(
            {
                "cohort": "GSE39055",
                "endpoint": "recurrence_free_survival",
                "module_label": module,
                "observed_c_index": observed,
                "n_module_genes_available": len(target_genes),
                "n_random_valid": len(random_array),
                "random_mean": (
                    float(np.mean(random_array)) if len(random_array) else np.nan
                ),
                "random_median": (
                    float(np.median(random_array)) if len(random_array) else np.nan
                ),
                "random_q95": (
                    float(np.quantile(random_array, 0.95))
                    if len(random_array)
                    else np.nan
                ),
                "observed_percentile": (
                    float(np.mean(random_array <= observed))
                    if len(random_array)
                    else np.nan
                ),
                "empirical_p_greater_equal": empirical_p,
            }
        )

    summary = pd.DataFrame(summary_rows)
    if not summary.empty:
        summary["empirical_q_bh"] = bh_adjust(
            summary["empirical_p_greater_equal"]
        )
    return summary, pd.DataFrame(distribution_rows)
```

#### `create_manifest`

- Lines: 925-960
- Signals: `direction`

```python
def create_manifest(
    preparation_manifest: dict[str, Any],
    outputs: list[Path],
) -> None:
    payload = {
        "script_version": SCRIPT_VERSION,
        "created_utc": utc_now_iso(),
        "preparation_manifest_sha256": sha256_file(PREPARATION_MANIFEST_FILE),
        "preparation_script_version": preparation_manifest.get("script_version"),
        "analysis_design": {
            "primary_modules": PRIMARY_MODULES,
            "primary_score": "strict one-to-one signed-mean z-score",
            "primary_endpoint": "recurrence-free survival",
            "primary_time_rule": "exclude nonpositive recorded times",
            "zero_time_sensitivity": (
                f"replace nonpositive time with {EPSILON_TIME_MONTHS:.6f} months"
            ),
            "bootstrap_repetitions": N_BOOTSTRAP,
            "permutation_repetitions": N_PERMUTATIONS,
            "random_gene_set_repetitions": N_RANDOM_SETS,
            "multiplicity": "BH across the four frozen primary modules",
        },
        "files": {},
    }

    for path in outputs:
        if path.exists():
            payload["files"][path.name] = {
                "sha256": sha256_file(path),
                "size_bytes": path.stat().st_size,
            }

    OUTPUT_MANIFEST.write_text(
        json.dumps(payload, indent=2),
        encoding="utf-8",
    )
```

#### `write_readme`

- Lines: 963-992
- Signals: `weights, direction`

```python
def write_readme() -> None:
    text = f"""GSE39055 frozen-program recurrence-free-survival validation
Script version: {SCRIPT_VERSION}

Primary analysis
----------------
- Four frozen primary canine programs: M34, M11, M24, M40
- Strict one-to-one signed-mean human score
- Recurrence-free survival using the GEO recurrence/follow-up time and recurrence event
- Cox HR per SD, fixed-direction C-index, BH FDR across four modules

Zero-time handling
------------------
The primary Cox analysis excludes nonpositive recorded times.
A prespecified sensitivity analysis replaces zero time with one day
({EPSILON_TIME_MONTHS:.6f} months).

Adjustment hierarchy
--------------------
Age and sex are baseline sensitivity covariates.
Human proliferation PC1 is a mechanistic sensitivity adjustment.
Percent necrosis is post-treatment and is not treated as a primary baseline confounder.

Interpretation
--------------
GSE39055 is a small third human cohort. Bootstrap, leave-one-out,
score-variant, and expression-matched random-panel analyses are robustness diagnostics.
No result may be used to change frozen module membership, weights, direction, or tier.
"""
    OUTPUT_README.write_text(text, encoding="utf-8")
```

#### `main`

- Lines: 995-1208
- Signals: `weights, direction, coverage`

```python
def main() -> None:
    print("=" * 80)
    print("GSE39055 external RFS validation of frozen canine programs")
    print("=" * 80)
    print(f"Script version: {SCRIPT_VERSION}")
    print(f"Project root: {PROJECT_ROOT}")
    print(f"Results directory: {RESULTS_DIR}")
    print("")
    print("Design:")
    print("  Preserve frozen canine module definitions and score direction.")
    print("  Test four strict one-to-one signed-mean scores for recurrence-free survival.")
    print("  Apply BH multiplicity control across four primary modules.")
    print("  Audit zero-time handling, baseline/proliferation adjustment, PH, leave-one-out stability, and matched random panels.")
    print("")

    preparation_manifest = verify_preparation_inputs()

    expression = read_required_csv(EXPRESSION_FILE, index_col=0)
    clinical = read_required_csv(CLINICAL_FILE, index_col=0)
    scores = read_required_csv(SCORES_FILE, index_col=0)
    coverage = read_required_csv(COVERAGE_FILE)
    frozen_programs = read_required_csv(FROZEN_PROGRAM_MANIFEST_FILE)
    strict_weights = read_required_csv(STRICT_WEIGHTS_FILE)

    common = clinical.index.intersection(scores.index).intersection(expression.index)
    clinical = clinical.loc[common].copy()
    scores = scores.loc[common].copy()
    expression = expression.loc[common].copy()

    primary_frame, epsilon_frame, zero_time_audit = create_analysis_frames(
        clinical,
        scores,
    )
    zero_time_audit.to_csv(OUTPUT_ZERO_TIME_AUDIT)

    print("")
    print("Matched GSE39055 validation data:")
    print(f"  Expression: {expression.shape}")
    print(f"  Clinical: {clinical.shape}")
    print(f"  Scores: {scores.shape}")
    print(f"  Complete RFS with positive time: {primary_frame.shape[0]}")
    print(f"  Recurrence events in primary frame: {int(primary_frame[EVENT_COL].sum())}")
    print(f"  Nonpositive-time records excluded from primary: {zero_time_audit.shape[0]}")

    primary = primary_validation(primary_frame, epsilon_frame)
    primary.to_csv(OUTPUT_PRIMARY, index=False)
    primary[
        [
            "module_label",
            "primary_p",
            "q_within_gse39055",
            "permutation_c_index_p_two_sided",
            "permutation_c_index_q_bh",
        ]
    ].to_csv(OUTPUT_MULTIPLICITY, index=False)

    adjusted = adjusted_validation(primary_frame)
    adjusted.to_csv(OUTPUT_ADJUSTED, index=False)

    variants = score_variant_sensitivity(primary_frame, coverage)
    variants.to_csv(OUTPUT_VARIANTS, index=False)

    loo = leave_one_out_stability(primary_frame)
    loo.to_csv(OUTPUT_LOO, index=False)

    random_summary, random_distribution = random_gene_set_controls(
        primary_frame,
        expression,
        strict_weights,
    )
    random_summary.to_csv(OUTPUT_RANDOM, index=False)
    random_distribution.to_csv(OUTPUT_RANDOM_DISTRIBUTION, index=False)

    synthesis = create_cross_cohort_synthesis(primary, random_summary)
    synthesis.to_csv(OUTPUT_CROSS_COHORT, index=False)

    write_readme()
    create_manifest(
        preparation_manifest,
        [
            OUTPUT_PRIMARY,
            OUTPUT_MULTIPLICITY,
            OUTPUT_ADJUSTED,
            OUTPUT_VARIANTS,
            OUTPUT_LOO,
            OUTPUT_RANDOM,
            OUTPUT_RANDOM_DISTRIBUTION,
            OUTPUT_ZERO_TIME_AUDIT,
            OUTPUT_CROSS_COHORT,
            OUTPUT_README,
        ],
    )

    print("")
    print("=" * 80)
    print("GSE39055 primary frozen-program RFS validation")
    print("=" * 80)
    display_cols = [
        "module_label",
        "n",
        "events",
        "hr_per_sd",
        "ci_low",
        "ci_high",
        "primary_p",
        "q_within_gse39055",
        "fixed_score_c_index",
        "fixed_score_c_index_ci_low",
        "fixed_score_c_index_ci_high",
        "permutation_c_index_q_bh",
        "ph_test_p",
        "gse39055_support_class",
    ]
    print(primary[display_cols].to_string(index=False))

    print("")
    print("=" * 80)
    print("Zero-time sensitivity")
    print("=" * 80)
    print(
        primary[
            [
                "module_label",
                "hr_per_sd",
                "primary_p",
                "epsilon_zero_time_hr_per_sd",
                "epsilon_zero_time_p",
            ]
        ].to_string(index=False)
    )

    print("")
    print("=" * 80)
    print("GSE39055 adjustment robustness")
    print("=" * 80)
    adjusted_display = adjusted[
        adjusted["module_label"].isin(PRIMARY_MODULES)
    ].copy()
    print(
        adjusted_display[
            [
                "module_label",
                "score_variant",
                "adjustment",
                "n",
                "events",
                "hr_per_sd",
                "ci_low",
                "ci_high",
                "p",
                "c_index",
                "ph_test_p",
            ]
        ].to_string(index=False)
    )

    print("")
    print("=" * 80)
    print("GSE39055 leave-one-out stability")
    print("=" * 80)
    print(loo.to_string(index=False))

    print("")
    print("=" * 80)
    print("GSE39055 expression-matched random controls")
    print("=" * 80)
    print(random_summary.to_string(index=False))

    print("")
    print("=" * 80)
    print("Three-cohort external evidence synthesis")
    print("=" * 80)
    synthesis_cols = [
        "module_label",
        "target_hr_per_sd",
        "gse_met_auc",
        "gse39055_rfs_hr_per_sd",
        "gse39055_rfs_q",
        "gse39055_rfs_c_index",
        "gse39055_random_empirical_p",
        "three_cohort_evidence_grade",
    ]
    synthesis_cols = [
        column for column in synthesis_cols if column in synthesis.columns
    ]
    print(synthesis[synthesis_cols].to_string(index=False))

    print("")
    print("=" * 80)
    print("Interpretation guardrails")
    print("=" * 80)
    print("The strict one-to-one signed-mean score and canine risk direction remain frozen.")
    print("Primary Cox analyses exclude nonpositive recorded follow-up times; one-day replacement is only a sensitivity analysis.")
    print("Age and sex are baseline sensitivity covariates; percent necrosis is post-treatment and cannot be treated as a primary baseline confounder.")
    print("GSE39055 contains only 37 samples and 18 recurrences before zero-time exclusion; external effect sizes require cautious interpretation.")
    print("No GSE39055 result may alter module membership, score weights, risk direction, or validation tier.")

    print("")
    print("Saved:")
    for path in [
        OUTPUT_PRIMARY,
        OUTPUT_MULTIPLICITY,
        OUTPUT_ADJUSTED,
        OUTPUT_VARIANTS,
        OUTPUT_LOO,
        OUTPUT_RANDOM,
        OUTPUT_RANDOM_DISTRIBUTION,
        OUTPUT_ZERO_TIME_AUDIT,
        OUTPUT_CROSS_COHORT,
        OUTPUT_README,
        OUTPUT_MANIFEST,
    ]:
        print(path)
    print("Done.")
```

### Relevant standalone lines

| Line | Signals | Code |
|---:|---|---|
| 131 | std | `std = values.std()` |
| 134 | mean | `return (values - values.mean()) / std` |
| 197 | std | `if frame[score_col].std() == 0:` |
| 225 | rank | `time_transform="rank",` |
| 342 | mean | `float(np.mean(coef_values > 0)),` |
| 585 | coverage | `coverage: pd.DataFrame,` |
| 644 | mean | `"fraction_hr_above_1": float((table["hr"] > 1).mean()),` |
| 649 | mean | `(table["p"] < 0.05).mean()` |
| 659 | mean | `"mean": expression.mean(axis=0),` |
| 660 | std | `"sd": expression.std(axis=0),` |
| 730 | std | `expression = expression.loc[:, expression.std(axis=0) > 0]` |
| 731 | mean, std | `z = (expression - expression.mean(axis=0)) / expression.std(axis=0)` |
| 743 | weights | `weights = strict_weights[` |
| 746 | weights | `weights["human_gene_symbol"] = (` |
| 747 | weights | `weights["human_gene_symbol"].astype(str).str.upper()` |
| 749 | weights | `weights = weights.drop_duplicates("human_gene_symbol", keep="first")` |
| 752 | weights | `for gene in weights["human_gene_symbol"].tolist()` |
| 755 | loadings | `loadings = (` |
| 756 | weights | `weights.set_index("human_gene_symbol")` |
| 760 | loadings, direction | `signs = np.sign(loadings.values)` |
| 781 | mean | `random_score = z[genes].mul(permuted_signs, axis=1).mean(axis=1)` |
| 814 | mean | `float(np.mean(random_array)) if len(random_array) else np.nan` |
| 825 | mean | `float(np.mean(random_array <= observed))` |
| 936 | direction | `"primary_score": "strict one-to-one signed-mean z-score",` |
| 970 | direction | `- Strict one-to-one signed-mean human score` |
| 972 | direction | `- Cox HR per SD, fixed-direction C-index, BH FDR across four modules` |
| 990 | weights, direction | `No result may be used to change frozen module membership, weights, direction, or tier.` |
| 1004 | direction | `print("  Preserve frozen canine module definitions and score direction.")` |
| 1005 | direction | `print("  Test four strict one-to-one signed-mean scores for recurrence-free survival.")` |
| 1015 | coverage | `coverage = read_required_csv(COVERAGE_FILE)` |
| 1054 | coverage | `variants = score_variant_sensitivity(primary_frame, coverage)` |
| 1186 | direction | `print("The strict one-to-one signed-mean score and canine risk direction remain frozen.")` |
| 1190 | weights, direction | `print("No GSE39055 result may alter module membership, score weights, risk direction, or validation tier.")` |

### Referenced result-table schemas

#### `results/tables/GSE238110_frozen_canine_transfer_program_manifest.csv`

- Rows: 13
- SHA-256: `fcd34b985ea99986f39dfc43ce8e4018db2fde9ecc8a225bb40ac93006d5372c`
- Columns: `module_label`, `validation_tier`, `multiplicity_family`, `provisional_program_label`, `program_label_requires_enrichment_confirmation`, `canine_primary_endpoint`, `canine_secondary_endpoint`, `positive_score_interpretation`, `n_canine_genes_used_for_pca`, `canine_pc1_explained_variance`, `canine_pca_orientation_correlation`, `dfi_full_cohort_coef`, `os_full_cohort_coef`, `risk_orientation_multiplier`, `n_strict_human_genes`, `n_broad_human_genes`, `strict_transfer_eligible`, `module_transfer_qc_tier`, `transfer_priority_score`, `fraction_strict_symbol_concordant`, `fraction_broad_transferable`, `raw_module_proliferation_correlation`, `orthogonal_variance_fraction_1_minus_r2`, `n_overlap_symbols_with_proliferation`, `primary_human_score`, `secondary_human_score`, `sensitivity_human_scores`, `frozen_after_canine_script`, `dfi_module_only_mean_c_index`, `dfi_module_only_std_c_index`, `dfi_module_only_fraction_above_0_50`, `dfi_module_plus_disjoint_proliferation_mean_c_index`, `dfi_module_plus_disjoint_proliferation_std_c_index`, `dfi_module_plus_disjoint_proliferation_fraction_above_0_50`, `dfi_residual_to_disjoint_proliferation_mean_c_index`, `dfi_residual_to_disjoint_proliferation_std_c_index`, `dfi_residual_to_disjoint_proliferation_fraction_above_0_50`, `dfi_residual_to_disjoint_proliferation_and_weight_mean_c_index`, `dfi_residual_to_disjoint_proliferation_and_weight_std_c_index`, `dfi_residual_to_disjoint_proliferation_and_weight_fraction_above_0_50`, `os_module_only_mean_c_index`, `os_module_only_std_c_index`, `os_module_only_fraction_above_0_50`, `os_module_plus_disjoint_proliferation_mean_c_index`, `os_module_plus_disjoint_proliferation_std_c_index`, `os_module_plus_disjoint_proliferation_fraction_above_0_50`, `os_residual_to_disjoint_proliferation_mean_c_index`, `os_residual_to_disjoint_proliferation_std_c_index`, `os_residual_to_disjoint_proliferation_fraction_above_0_50`, `os_residual_to_disjoint_proliferation_and_weight_mean_c_index`, `os_residual_to_disjoint_proliferation_and_weight_std_c_index`, `os_residual_to_disjoint_proliferation_and_weight_fraction_above_0_50`, `script20_dfi_recommended_role`, `script20_os_recommended_role`, `manual_freeze_reason`

#### `results/tables/GSE238110_frozen_transfer_gene_weights_broad.csv`

- Rows: 389
- SHA-256: `ee5f661677b173507e04856f848b6757ede2750dd69b9bb0d88a9f6a2a77bd4c`
- Columns: `module_label`, `canine_gene`, `canine_gene_symbol`, `human_gene_symbol`, `ortholog_qc_status`, `raw_pca_loading`, `risk_oriented_loading`, `absolute_risk_loading`, `is_strict_mapping`, `is_broad_mapping`, `human_symbol_duplicate_count`, `normalized_abs_sum_weight`

#### `results/tables/GSE238110_frozen_transfer_gene_weights_strict.csv`

- Rows: 321
- SHA-256: `4f065aa4c4edf117a0c74015840d2b4b2347929f172cd517e1818ba0f6163b91`
- Columns: `module_label`, `canine_gene`, `canine_gene_symbol`, `human_gene_symbol`, `ortholog_qc_status`, `raw_pca_loading`, `risk_oriented_loading`, `absolute_risk_loading`, `is_strict_mapping`, `is_broad_mapping`, `human_symbol_duplicate_count`, `normalized_abs_sum_weight`

#### `results/tables/GSE39055_RFS_adjustment_robustness.csv`

- Rows: 30
- SHA-256: `746e15d82a3b529a9bbb0c435dcaf97ddb30a3edff72c24934bc8a9f6ce626f4`
- Columns: `module_label`, `score_variant`, `adjustment`, `covariates`, `n`, `events`, `hr_per_sd`, `ci_low`, `ci_high`, `p`, `coef`, `se_coef`, `c_index`, `ph_test_p`, `error`, `interpretation_note`

#### `results/tables/GSE39055_RFS_leave_one_out_stability.csv`

- Rows: 4
- SHA-256: `294015c8df103c6e280b02b5cefbcd72bbe1237a9ff9e156aeb1ffb198459d47`
- Columns: `module_label`, `n_loo_fits`, `hr_min`, `hr_max`, `hr_median`, `fraction_hr_above_1`, `c_index_min`, `c_index_max`, `c_index_median`, `fraction_p_below_0_05`

#### `results/tables/GSE39055_RFS_primary_frozen_program_validation.csv`

- Rows: 4
- SHA-256: `1945cca593030f57cbaf6eef3006391b4fc3ad2ff334d99e34cc81436533951e`
- Columns: `module_label`, `score_column`, `n`, `events`, `hr_per_sd`, `ci_low`, `ci_high`, `primary_p`, `fixed_score_c_index`, `fixed_score_c_index_ci_low`, `fixed_score_c_index_ci_high`, `c_index_bootstrap_valid`, `ph_test_p`, `bootstrap_hr_ci_low`, `bootstrap_hr_ci_high`, `bootstrap_probability_coef_positive`, `hr_bootstrap_valid`, `permutation_c_index_p_one_sided`, `permutation_c_index_p_two_sided`, `epsilon_zero_time_hr_per_sd`, `epsilon_zero_time_ci_low`, `epsilon_zero_time_ci_high`, `epsilon_zero_time_p`, `error`, `q_within_gse39055`, `permutation_c_index_q_bh`, `gse39055_support_class`

#### `results/tables/GSE39055_RFS_primary_multiplicity.csv`

- Rows: 4
- SHA-256: `66dc8e035dc52a60659348fe9c6fca999ec57e3d14e95d09ce8385d27318bc92`
- Columns: `module_label`, `primary_p`, `q_within_gse39055`, `permutation_c_index_p_two_sided`, `permutation_c_index_q_bh`

#### `results/tables/GSE39055_RFS_random_gene_set_controls.csv`

- Rows: 4
- SHA-256: `2031aa3de080aa1b0b785ed4e1a86ab01a90f8660236dd09ea3111307ffe9b72`
- Columns: `cohort`, `endpoint`, `module_label`, `observed_c_index`, `n_module_genes_available`, `n_random_valid`, `random_mean`, `random_median`, `random_q95`, `observed_percentile`, `empirical_p_greater_equal`, `empirical_q_bh`

#### `results/tables/GSE39055_RFS_random_gene_set_distribution.csv`

- Rows: 8000
- SHA-256: `855203a55f412afaf3eba1ee7ba808f75bc847cfaa70d58ea2bbc009579241d3`
- Columns: `module_label`, `repeat`, `random_c_index`, `n_genes`

#### `results/tables/GSE39055_RFS_score_variant_sensitivity.csv`

- Rows: 45
- SHA-256: `af0df2d0b1088ed3a255acaa76203f43caf732b8cbe9547ab0306076ad674469`
- Columns: `module_label`, `variant`, `score_column`, `n`, `events`, `hr_per_sd`, `ci_low`, `ci_high`, `p`, `coef`, `se_coef`, `c_index`, `ph_test_p`, `error`

#### `results/tables/GSE39055_frozen_transfer_score_coverage.csv`

- Rows: 28
- SHA-256: `1945407e793a0c3d3213a55bc0e8cf7526dbfb8b9d65e716818b5f5748ec50e7`
- Columns: `cohort`, `module_label`, `mapping`, `n_frozen_genes`, `n_available_genes`, `coverage_fraction`, `minimum_rule_passed`, `available_genes`, `missing_genes`, `validation_tier`, `score`

#### `results/tables/GSE39055_zero_time_endpoint_audit.csv`

- Rows: 1
- SHA-256: `5bce5fe93c71ddd7062f70986ce2623be62974cfa2360b463a3a1f75ebd8281b`
- Columns: `geo_sample_id`, `title`, `age_years`, `sex`, `chemotherapy`, `percent_necrosis_raw`, `percent_necrosis_numeric`, `good_necrosis_response_ge90`, `recurrence_event`, `death_event_descriptive`, `rfs_time_months`, `tissue`, `biopsy_resection_pair`, `metadata_text_combined`, `endpoint_note`, `M11__strict__signed_mean_z`, `M11__strict__canine_pca_weighted_z`, `M11__strict__human_pc1_z`, `M24__strict__signed_mean_z`, `M24__strict__canine_pca_weighted_z`, `M24__strict__human_pc1_z`, `M25__strict__signed_mean_z`, `M25__strict__canine_pca_weighted_z`, `M25__strict__human_pc1_z`, `M27__strict__signed_mean_z`, `M27__strict__canine_pca_weighted_z`, `M27__strict__human_pc1_z`, `M28__strict__signed_mean_z`, `M28__strict__canine_pca_weighted_z`, `M28__strict__human_pc1_z`, `M34__strict__signed_mean_z`, `M34__strict__canine_pca_weighted_z`, `M34__strict__human_pc1_z`, `M38__strict__signed_mean_z`, `M38__strict__canine_pca_weighted_z`, `M38__strict__human_pc1_z`, `M39__strict__signed_mean_z`, `M39__strict__canine_pca_weighted_z`, `M39__strict__human_pc1_z`, `M40__strict__signed_mean_z`, `M40__strict__canine_pca_weighted_z`, `M40__strict__human_pc1_z`, `M11__broad__signed_mean_z`, `M11__broad__canine_pca_weighted_z`, `M11__broad__human_pc1_z`, `M17__broad__signed_mean_z`, `M17__broad__canine_pca_weighted_z`, `M17__broad__human_pc1_z`, `M20__broad__signed_mean_z`, `M20__broad__canine_pca_weighted_z`, `M20__broad__human_pc1_z`, `M24__broad__signed_mean_z`, `M24__broad__canine_pca_weighted_z`, `M24__broad__human_pc1_z`, `M25__broad__signed_mean_z`, `M25__broad__canine_pca_weighted_z`, `M25__broad__human_pc1_z`, `M27__broad__signed_mean_z`, `M27__broad__canine_pca_weighted_z`, `M27__broad__human_pc1_z`, `M28__broad__signed_mean_z`, `M28__broad__canine_pca_weighted_z`, `M28__broad__human_pc1_z`, `M30__broad__signed_mean_z`, `M30__broad__canine_pca_weighted_z`, `M30__broad__human_pc1_z`, `M34__broad__signed_mean_z`, `M34__broad__canine_pca_weighted_z`, `M34__broad__human_pc1_z`, `M38__broad__signed_mean_z`, `M38__broad__canine_pca_weighted_z`, `M38__broad__human_pc1_z`, `M39__broad__signed_mean_z`, `M39__broad__canine_pca_weighted_z`, `M39__broad__human_pc1_z`, `M40__broad__signed_mean_z`, `M40__broad__canine_pca_weighted_z`, `M40__broad__human_pc1_z`, `strict_human_meta_proliferation_pc1_z`, `M40_disjoint_strict_human_meta_proliferation_pc1_z`, `M40__strict__signed_mean__residual_to_disjoint_proliferation_z`, `M40__strict__canine_pca_weighted__residual_to_disjoint_proliferation_z`

#### `results/tables/human_external_validation_robust_evidence_summary.csv`

- Rows: 4
- SHA-256: `3420c827b16e1cadad3877edd06c2d843534cf0a517024fa6a83d81e4427103c`
- Columns: `module_label`, `target_hr_per_sd`, `target_primary_p`, `target_q_within_endpoint`, `target_c_index`, `gse_met_auc`, `gse_q_within_endpoint`, `global_q_eight_tests`, `gse_robust_logistic_or_per_sd`, `gse_permutation_auc_q_bh`, `direction_consistent`, `target_loo_fraction_hr_above_1`, `gse_loo_fraction_auc_above_0_50`, `target_random_empirical_p`, `gse_random_empirical_p`, `robust_external_evidence_grade`, `interpretation`

#### `results/tables/human_external_validation_three_cohort_synthesis.csv`

- Rows: 4
- SHA-256: `d40b6c6446a24173c671ca9d9a32416478a0e620f762dd29c1a1fe2d60771fbd`
- Columns: `module_label`, `target_hr_per_sd`, `target_primary_p`, `target_q_within_endpoint`, `target_c_index`, `gse_met_auc`, `gse_q_within_endpoint`, `global_q_eight_tests`, `gse_robust_logistic_or_per_sd`, `gse_permutation_auc_q_bh`, `direction_consistent`, `target_loo_fraction_hr_above_1`, `gse_loo_fraction_auc_above_0_50`, `target_random_empirical_p`, `gse_random_empirical_p`, `robust_external_evidence_grade`, `interpretation`, `gse39055_rfs_hr_per_sd`, `gse39055_rfs_ci_low`, `gse39055_rfs_ci_high`, `gse39055_rfs_p`, `gse39055_rfs_q`, `gse39055_rfs_c_index`, `gse39055_rfs_ph_p`, `gse39055_support_class`, `gse39055_random_percentile`, `gse39055_random_empirical_p`, `gse39055_random_empirical_q`, `three_cohort_evidence_grade`

---

## `scripts/28_conservative_module_preservation_audit.py`

- SHA-256: `a2f9631f1ebc358b6d5db00daa6a1ee754eabcf0ac4e0f73326dbe7e27c33b88`

### Relevant functions

#### `zscore_matrix`

- Lines: 103-109
- Signals: `mean, std`

```python
def zscore_matrix(expression: pd.DataFrame) -> pd.DataFrame:
    x = expression.apply(pd.to_numeric, errors="coerce")
    x = x.replace([np.inf, -np.inf], np.nan)
    x = x.fillna(x.median(axis=0))
    std = x.std(axis=0).replace(0, np.nan)
    z = (x - x.mean(axis=0)) / std
    return z.loc[:, z.notna().all(axis=0)]
```

#### `safe_spearman`

- Lines: 112-118
- Signals: `std`

```python
def safe_spearman(a: np.ndarray, b: np.ndarray) -> float:
    mask = np.isfinite(a) & np.isfinite(b)
    if mask.sum() < 3:
        return np.nan
    if np.std(a[mask]) == 0 or np.std(b[mask]) == 0:
        return np.nan
    return float(stats.spearmanr(a[mask], b[mask]).statistic)
```

#### `safe_pearson`

- Lines: 121-127
- Signals: `std`

```python
def safe_pearson(a: np.ndarray, b: np.ndarray) -> float:
    mask = np.isfinite(a) & np.isfinite(b)
    if mask.sum() < 3:
        return np.nan
    if np.std(a[mask]) == 0 or np.std(b[mask]) == 0:
        return np.nan
    return float(np.corrcoef(a[mask], b[mask])[0, 1])
```

#### `orient_pc1_to_frozen_score`

- Lines: 135-157
- Signals: `mean, loadings, matrix_product, correlation`

```python
def orient_pc1_to_frozen_score(
    z: pd.DataFrame,
    frozen_signs: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, float]:
    values = z.values.astype(float)
    _, singular_values, vt = np.linalg.svd(values, full_matrices=False)
    loading = vt[0].copy()
    score = values @ loading
    frozen_score = (values * frozen_signs).mean(axis=1)

    corr = safe_pearson(score, frozen_score)
    if np.isfinite(corr) and corr < 0:
        loading = -loading
        score = -score
        corr = -corr

    denominator = float(np.sum(singular_values ** 2))
    variance_explained = (
        float(singular_values[0] ** 2 / denominator)
        if denominator > 0
        else np.nan
    )
    return loading, score, variance_explained
```

#### `align_module_matrices`

- Lines: 160-201
- Signals: `weights`

```python
def align_module_matrices(
    canine_expression: pd.DataFrame,
    human_expression: pd.DataFrame,
    weights: pd.DataFrame,
    module: str,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    part = weights[weights["module_label"].eq(module)].copy()
    part["canine_gene"] = part["canine_gene"].astype(str)
    part["human_gene_symbol"] = (
        part["human_gene_symbol"].astype(str).str.upper()
    )
    part = part.drop_duplicates("human_gene_symbol", keep="first")
    part = part.drop_duplicates("canine_gene", keep="first")

    canine_columns = {str(column): column for column in canine_expression.columns}

    human = human_expression.copy()
    human.columns = human.columns.astype(str).str.upper()
    human = human.loc[:, ~human.columns.duplicated()].copy()

    available = part[
        part["canine_gene"].isin(canine_columns)
        & part["human_gene_symbol"].isin(human.columns)
    ].copy()

    if available.shape[0] < 3:
        return (
            pd.DataFrame(index=canine_expression.index),
            pd.DataFrame(index=human.index),
            available,
        )

    canine = canine_expression[
        [canine_columns[gene] for gene in available["canine_gene"]]
    ].copy()
    canine.columns = available["human_gene_symbol"].tolist()

    human_matrix = human[
        available["human_gene_symbol"].tolist()
    ].copy()

    return canine, human_matrix, available
```

#### `correlation_edge_audit`

- Lines: 204-271
- Signals: `mean, direction`

```python
def correlation_edge_audit(
    canine_z: pd.DataFrame,
    human_z: pd.DataFrame,
    seed: int,
) -> dict[str, float]:
    canine_corr = np.corrcoef(canine_z.values, rowvar=False)
    human_corr = np.corrcoef(human_z.values, rowvar=False)

    canine_edges = upper_triangle(canine_corr)
    human_edges = upper_triangle(human_corr)

    observed_spearman = safe_spearman(canine_edges, human_edges)
    observed_pearson = safe_pearson(canine_edges, human_edges)

    valid = (
        np.isfinite(canine_edges)
        & np.isfinite(human_edges)
        & (canine_edges != 0)
        & (human_edges != 0)
    )
    sign_concordance = (
        float(
            np.mean(
                np.sign(canine_edges[valid])
                == np.sign(human_edges[valid])
            )
        )
        if valid.sum() > 0
        else np.nan
    )

    rng = np.random.default_rng(seed)
    null = np.empty(N_EDGE_PERMUTATIONS, dtype=float)
    n_genes = human_corr.shape[0]

    for index in range(N_EDGE_PERMUTATIONS):
        permutation = rng.permutation(n_genes)
        permuted_human = human_corr[np.ix_(permutation, permutation)]
        null[index] = safe_spearman(
            canine_edges,
            upper_triangle(permuted_human),
        )

    empirical_p_positive = (
        (1.0 + np.sum(null >= observed_spearman))
        / (N_EDGE_PERMUTATIONS + 1.0)
        if np.isfinite(observed_spearman)
        else np.nan
    )
    empirical_p_two_sided = (
        (
            1.0
            + np.sum(np.abs(null) >= abs(observed_spearman))
        )
        / (N_EDGE_PERMUTATIONS + 1.0)
        if np.isfinite(observed_spearman)
        else np.nan
    )

    return {
        "edge_spearman": observed_spearman,
        "edge_pearson": observed_pearson,
        "edge_sign_concordance": sign_concordance,
        "edge_permutation_p_positive": empirical_p_positive,
        "edge_permutation_p_two_sided": empirical_p_two_sided,
        "edge_null_mean": float(np.nanmean(null)),
        "edge_null_q95": float(np.nanquantile(null, 0.95)),
    }
```

#### `loading_audit`

- Lines: 274-332
- Signals: `mean, direction`

```python
def loading_audit(
    human_z: pd.DataFrame,
    frozen_loadings: np.ndarray,
    seed: int,
) -> dict[str, float]:
    frozen_signs = np.sign(frozen_loadings)
    frozen_signs[frozen_signs == 0] = 1.0

    human_loading, _, variance_explained = orient_pc1_to_frozen_score(
        human_z,
        frozen_signs,
    )

    observed_spearman = safe_spearman(
        human_loading,
        frozen_loadings,
    )
    observed_sign_concordance = float(
        np.mean(
            np.sign(human_loading)
            == np.sign(frozen_loadings)
        )
    )

    rng = np.random.default_rng(seed)
    null = np.empty(N_LOADING_PERMUTATIONS, dtype=float)

    for index in range(N_LOADING_PERMUTATIONS):
        permuted = rng.permutation(human_loading)
        null[index] = safe_spearman(
            permuted,
            frozen_loadings,
        )

    empirical_p_positive = (
        (1.0 + np.sum(null >= observed_spearman))
        / (N_LOADING_PERMUTATIONS + 1.0)
        if np.isfinite(observed_spearman)
        else np.nan
    )
    empirical_p_two_sided = (
        (
            1.0
            + np.sum(np.abs(null) >= abs(observed_spearman))
        )
        / (N_LOADING_PERMUTATIONS + 1.0)
        if np.isfinite(observed_spearman)
        else np.nan
    )

    return {
        "human_pc1_variance_explained": variance_explained,
        "loading_spearman": observed_spearman,
        "loading_sign_concordance": observed_sign_concordance,
        "loading_permutation_p_positive": empirical_p_positive,
        "loading_permutation_p_two_sided": empirical_p_two_sided,
        "loading_null_mean": float(np.nanmean(null)),
        "loading_null_q95": float(np.nanquantile(null, 0.95)),
    }
```

#### `unique_small_module_splits`

- Lines: 335-345
- Signals: ``

```python
def unique_small_module_splits(n_genes: int) -> list[tuple[np.ndarray, np.ndarray]]:
    first_size = n_genes // 2
    all_indices = np.arange(n_genes)
    splits: list[tuple[np.ndarray, np.ndarray]] = []

    for first_tuple in combinations(range(n_genes), first_size):
        first = np.asarray(first_tuple, dtype=int)
        second = np.setdiff1d(all_indices, first)
        splits.append((first, second))

    return splits
```

#### `random_large_module_splits`

- Lines: 348-368
- Signals: ``

```python
def random_large_module_splits(
    n_genes: int,
    seed: int,
) -> list[tuple[np.ndarray, np.ndarray]]:
    rng = np.random.default_rng(seed)
    first_size = n_genes // 2
    all_indices = np.arange(n_genes)
    splits: list[tuple[np.ndarray, np.ndarray]] = []

    for _ in range(N_SPLIT_REPEATS):
        first = np.sort(
            rng.choice(
                all_indices,
                size=first_size,
                replace=False,
            )
        )
        second = np.setdiff1d(all_indices, first)
        splits.append((first, second))

    return splits
```

#### `disjoint_split_half_reliability`

- Lines: 371-421
- Signals: `mean, direction`

```python
def disjoint_split_half_reliability(
    human_z: pd.DataFrame,
    frozen_loadings: np.ndarray,
    seed: int,
) -> tuple[dict[str, float], pd.DataFrame]:
    signs = np.sign(frozen_loadings)
    signs[signs == 0] = 1.0
    signed_values = human_z.values.astype(float) * signs

    n_genes = human_z.shape[1]
    if n_genes <= 12:
        splits = unique_small_module_splits(n_genes)
        method = "all_balanced_disjoint_splits"
    else:
        splits = random_large_module_splits(n_genes, seed)
        method = "random_balanced_disjoint_splits"

    rows: list[dict[str, Any]] = []

    for iteration, (first, second) in enumerate(splits, start=1):
        first_score = signed_values[:, first].mean(axis=1)
        second_score = signed_values[:, second].mean(axis=1)
        correlation = safe_pearson(first_score, second_score)

        rows.append(
            {
                "iteration": iteration,
                "method": method,
                "n_genes_first_half": len(first),
                "n_genes_second_half": len(second),
                "disjoint_half_score_correlation": correlation,
            }
        )

    table = pd.DataFrame(rows)
    values = pd.to_numeric(
        table["disjoint_half_score_correlation"],
        errors="coerce",
    ).dropna()

    summary = {
        "split_half_method": method,
        "n_split_half_valid": int(values.shape[0]),
        "split_half_mean": float(values.mean()),
        "split_half_median": float(values.median()),
        "split_half_q05": float(values.quantile(0.05)),
        "split_half_q95": float(values.quantile(0.95)),
        "split_half_fraction_positive": float((values > 0).mean()),
        "split_half_fraction_above_0_30": float((values > 0.30).mean()),
    }
    return summary, table
```

#### `main`

- Lines: 464-738
- Signals: `weights, loadings, direction, correlation`

```python
def main() -> None:
    print("=" * 80)
    print("Conservative cross-cohort module-preservation audit")
    print("=" * 80)
    print(f"Script version: {SCRIPT_VERSION}")
    print(f"Project root: {PROJECT_ROOT}")
    print("")
    print("Purpose:")
    print("  Replace mechanically inflated preservation metrics with direct tests.")
    print("  Compare canine and human within-module correlation matrices.")
    print("  Test PC1-loading concordance by gene-label permutation.")
    print("  Estimate reliability from non-overlapping gene halves.")
    print("  Preserve all frozen genes, weights, directions, and outcome results.")
    print("")

    weights = read_required_csv(STRICT_WEIGHTS_FILE)
    canine_expression = read_required_csv(
        CANINE_EXPRESSION_FILE,
        index_col=0,
    )
    human_expression = {
        "TARGET_OS": read_required_csv(
            TARGET_EXPRESSION_FILE,
            index_col=0,
        ),
        "GSE21257": read_required_csv(
            GSE21257_EXPRESSION_FILE,
            index_col=0,
        ),
        "GSE39055": read_required_csv(
            GSE39055_EXPRESSION_FILE,
            index_col=0,
        ),
    }

    old_preservation = (
        read_required_csv(SCRIPT27_PRESERVATION_FILE)
        if SCRIPT27_PRESERVATION_FILE.exists()
        else pd.DataFrame()
    )
    old_synthesis = (
        read_required_csv(SCRIPT27_SYNTHESIS_FILE)
        if SCRIPT27_SYNTHESIS_FILE.exists()
        else pd.DataFrame()
    )

    summary_rows: list[dict[str, Any]] = []
    split_tables: list[pd.DataFrame] = []

    job = 0
    total_jobs = len(HUMAN_COHORTS) * len(PRIMARY_MODULES)

    for cohort_index, cohort in enumerate(HUMAN_COHORTS):
        for module_index, module in enumerate(PRIMARY_MODULES):
            job += 1
            print(f"  Job {job}/{total_jobs}: {cohort} {module}")

            canine_matrix, human_matrix, mapping = align_module_matrices(
                canine_expression=canine_expression,
                human_expression=human_expression[cohort],
                weights=weights,
                module=module,
            )

            if mapping.shape[0] < 3:
                summary_rows.append(
                    {
                        "cohort": cohort,
                        "module_label": module,
                        "n_genes_shared": mapping.shape[0],
                        "error": "fewer_than_three_shared_genes",
                    }
                )
                continue

            canine_z = zscore_matrix(canine_matrix)
            human_z = zscore_matrix(human_matrix)

            shared = canine_z.columns.intersection(human_z.columns)
            canine_z = canine_z[shared]
            human_z = human_z[shared]
            mapping_indexed = (
                mapping.set_index("human_gene_symbol")
                .loc[shared]
            )
            frozen_loadings = pd.to_numeric(
                mapping_indexed["risk_oriented_loading"],
                errors="coerce",
            ).fillna(0.0).values

            seed = (
                RANDOM_SEED
                + cohort_index * 100000
                + module_index * 1000
            )

            edge = correlation_edge_audit(
                canine_z,
                human_z,
                seed=seed,
            )
            loading = loading_audit(
                human_z,
                frozen_loadings,
                seed=seed + 100,
            )
            split_summary, split_table = disjoint_split_half_reliability(
                human_z,
                frozen_loadings,
                seed=seed + 200,
            )

            split_table.insert(0, "module_label", module)
            split_table.insert(0, "cohort", cohort)
            split_tables.append(split_table)

            summary_rows.append(
                {
                    "cohort": cohort,
                    "module_label": module,
                    "n_canine_samples": canine_z.shape[0],
                    "n_human_samples": human_z.shape[0],
                    "n_genes_shared": len(shared),
                    **edge,
                    **loading,
                    **split_summary,
                    "error": "",
                }
            )

    audit = pd.DataFrame(summary_rows)
    audit["edge_permutation_q_positive"] = bh_adjust(
        audit["edge_permutation_p_positive"]
    )
    audit["loading_permutation_q_positive"] = bh_adjust(
        audit["loading_permutation_p_positive"]
    )
    audit["conservative_preservation_class"] = audit.apply(
        conservative_classification,
        axis=1,
    )

    if not old_synthesis.empty:
        old_columns = [
            "cohort",
            "module_label",
            "structural_support_class",
            "outcome_direction_concordant",
            "joint_preservation_outcome_class",
        ]
        old_columns = [
            column for column in old_columns if column in old_synthesis.columns
        ]
        audit = audit.merge(
            old_synthesis[old_columns],
            on=["cohort", "module_label"],
            how="left",
        )

    audit.to_csv(OUTPUT_EDGE, index=False)

    split_output = (
        pd.concat(split_tables, ignore_index=True)
        if split_tables
        else pd.DataFrame()
    )
    split_output.to_csv(OUTPUT_SPLIT_HALF, index=False)

    conservative_columns = [
        "cohort",
        "module_label",
        "n_genes_shared",
        "edge_spearman",
        "edge_permutation_p_positive",
        "edge_permutation_q_positive",
        "loading_spearman",
        "loading_permutation_p_positive",
        "loading_permutation_q_positive",
        "split_half_median",
        "split_half_q05",
        "split_half_fraction_positive",
        "conservative_preservation_class",
        "structural_support_class",
        "outcome_direction_concordant",
        "joint_preservation_outcome_class",
    ]
    conservative_columns = [
        column for column in conservative_columns if column in audit.columns
    ]
    conservative = audit[conservative_columns].copy()
    conservative.to_csv(OUTPUT_CONSERVATIVE, index=False)

    readme = f"""Conservative module-preservation audit
Script version: {SCRIPT_VERSION}

Why this audit was added
------------------------
The script 27 signed-mean/weighted-score correlation and overlapping-subset
reliability metrics can be high partly because the compared scores share most
or all genes. They are useful diagnostics but should not, by themselves,
define cross-cohort structural preservation.

Primary preservation evidence in this audit
-------------------------------------------
1. Spearman preservation of the full within-module gene-correlation matrix
   between canine DOG2 and each human cohort.
2. Concordance of human PC1 loadings with frozen canine risk-oriented loadings.
3. Correlation between scores formed from non-overlapping gene halves.

Permutation tests
-----------------
Gene-label permutations are used for correlation-matrix and loading
concordance. BH correction is applied across the 12 human cohort-module
comparisons for each test family.

Interpretation
--------------
Outcome association and representation preservation remain separate.
Frozen score direction is never reversed after viewing a human outcome.
"""
    OUTPUT_README.write_text(readme, encoding="utf-8")

    manifest = {
        "script_version": SCRIPT_VERSION,
        "created_utc": utc_now_iso(),
        "frozen_weights_sha256": sha256_file(STRICT_WEIGHTS_FILE),
        "edge_permutations": N_EDGE_PERMUTATIONS,
        "loading_permutations": N_LOADING_PERMUTATIONS,
        "large_module_split_repeats": N_SPLIT_REPEATS,
        "files": {},
    }
    for path in [
        OUTPUT_EDGE,
        OUTPUT_SPLIT_HALF,
        OUTPUT_CONSERVATIVE,
        OUTPUT_README,
    ]:
        if path.exists():
            manifest["files"][path.name] = {
                "sha256": sha256_file(path),
                "size_bytes": path.stat().st_size,
            }

    OUTPUT_MANIFEST.write_text(
        json.dumps(manifest, indent=2),
        encoding="utf-8",
    )

    print("")
    print("=" * 80)
    print("Conservative preservation classification")
    print("=" * 80)
    print(conservative.to_string(index=False))

    print("")
    print("=" * 80)
    print("Interpretation guardrails")
    print("=" * 80)
    print("Script 27 metrics remain descriptive and are not deleted.")
    print("High correlation between scores sharing the same genes is not treated as independent preservation evidence.")
    print("Non-overlapping split halves replace overlapping gene-subset reliability for the conservative classification.")
    print("Outcome-direction discordance is never repaired by flipping a frozen score.")
    print("Small seven-gene modules require cautious interpretation because edge and loading estimates are discrete and unstable.")

    print("")
    print("Saved:")
    for path in [
        OUTPUT_EDGE,
        OUTPUT_SPLIT_HALF,
        OUTPUT_CONSERVATIVE,
        OUTPUT_README,
        OUTPUT_MANIFEST,
    ]:
        print(path)
    print("Done.")
```

### Relevant standalone lines

| Line | Signals | Code |
|---:|---|---|
| 107 | std | `std = x.std(axis=0).replace(0, np.nan)` |
| 108 | mean | `z = (x - x.mean(axis=0)) / std` |
| 116 | std | `if np.std(a[mask]) == 0 or np.std(b[mask]) == 0:` |
| 125 | std | `if np.std(a[mask]) == 0 or np.std(b[mask]) == 0:` |
| 141 | loadings | `loading = vt[0].copy()` |
| 142 | loadings, matrix_product | `score = values @ loading` |
| 143 | mean | `frozen_score = (values * frozen_signs).mean(axis=1)` |
| 145 | correlation | `corr = safe_pearson(score, frozen_score)` |
| 146 | correlation | `if np.isfinite(corr) and corr < 0:` |
| 147 | loadings | `loading = -loading` |
| 149 | correlation | `corr = -corr` |
| 157 | loadings | `return loading, score, variance_explained` |
| 163 | weights | `weights: pd.DataFrame,` |
| 166 | weights | `part = weights[weights["module_label"].eq(module)].copy()` |
| 226 | mean | `np.mean(` |
| 227 | direction | `np.sign(canine_edges[valid])` |
| 228 | direction | `== np.sign(human_edges[valid])` |
| 269 | mean | `"edge_null_mean": float(np.nanmean(null)),` |
| 279 | direction | `frozen_signs = np.sign(frozen_loadings)` |
| 292 | mean | `np.mean(` |
| 293 | direction | `np.sign(human_loading)` |
| 294 | direction | `== np.sign(frozen_loadings)` |
| 330 | mean | `"loading_null_mean": float(np.nanmean(null)),` |
| 376 | direction | `signs = np.sign(frozen_loadings)` |
| 391 | mean | `first_score = signed_values[:, first].mean(axis=1)` |
| 392 | mean | `second_score = signed_values[:, second].mean(axis=1)` |
| 414 | mean | `"split_half_mean": float(values.mean()),` |
| 418 | mean | `"split_half_fraction_positive": float((values > 0).mean()),` |
| 419 | mean | `"split_half_fraction_above_0_30": float((values > 0.30).mean()),` |
| 474 | loadings | `print("  Test PC1-loading concordance by gene-label permutation.")` |
| 476 | weights | `print("  Preserve all frozen genes, weights, directions, and outcome results.")` |
| 479 | weights | `weights = read_required_csv(STRICT_WEIGHTS_FILE)` |
| 524 | weights | `weights=weights,` |
| 565 | loadings | `loading = loading_audit(` |
| 588 | loadings | `**loading,` |
| 661 | weights, direction | `The script 27 signed-mean/weighted-score correlation and overlapping-subset` |
| 668 | correlation | `1. Spearman preservation of the full within-module gene-correlation matrix` |
| 670 | loadings | `2. Concordance of human PC1 loadings with frozen canine risk-oriented loadings.` |
| 675 | loadings | `Gene-label permutations are used for correlation-matrix and loading` |
| 682 | direction | `Frozen score direction is never reversed after viewing a human outcome.` |
| 725 | direction | `print("Outcome-direction discordance is never repaired by flipping a frozen score.")` |
| 726 | loadings | `print("Small seven-gene modules require cautious interpretation because edge and loading estimates are discrete and unstable.")` |

### Referenced result-table schemas

#### `results/tables/GSE238110_frozen_transfer_gene_weights_strict.csv`

- Rows: 321
- SHA-256: `4f065aa4c4edf117a0c74015840d2b4b2347929f172cd517e1818ba0f6163b91`
- Columns: `module_label`, `canine_gene`, `canine_gene_symbol`, `human_gene_symbol`, `ortholog_qc_status`, `raw_pca_loading`, `risk_oriented_loading`, `absolute_risk_loading`, `is_strict_mapping`, `is_broad_mapping`, `human_symbol_duplicate_count`, `normalized_abs_sum_weight`

#### `results/tables/cross_cohort_conservative_preservation_classification.csv`

- Rows: 12
- SHA-256: `184fd8eade2b83039c776aa0a6dcae422897907469b571f3ae65b69be7f5f25e`
- Columns: `cohort`, `module_label`, `n_genes_shared`, `edge_spearman`, `edge_permutation_p_positive`, `edge_permutation_q_positive`, `loading_spearman`, `loading_permutation_p_positive`, `loading_permutation_q_positive`, `split_half_median`, `split_half_q05`, `split_half_fraction_positive`, `conservative_preservation_class`, `structural_support_class`, `outcome_direction_concordant`, `joint_preservation_outcome_class`

#### `results/tables/cross_cohort_disjoint_split_half_reliability.csv`

- Rows: 6210
- SHA-256: `79e0a786216c3f156c478a8e3ad3df7f220e99d68466298d3105a695a230084d`
- Columns: `cohort`, `module_label`, `iteration`, `method`, `n_genes_first_half`, `n_genes_second_half`, `disjoint_half_score_correlation`

#### `results/tables/cross_cohort_edge_preservation_audit.csv`

- Rows: 12
- SHA-256: `e29862c72b0f5f0d404e72dd1e1c4f81d1284dc07fabe0de058062f60220b861`
- Columns: `cohort`, `module_label`, `n_canine_samples`, `n_human_samples`, `n_genes_shared`, `edge_spearman`, `edge_pearson`, `edge_sign_concordance`, `edge_permutation_p_positive`, `edge_permutation_p_two_sided`, `edge_null_mean`, `edge_null_q95`, `human_pc1_variance_explained`, `loading_spearman`, `loading_sign_concordance`, `loading_permutation_p_positive`, `loading_permutation_p_two_sided`, `loading_null_mean`, `loading_null_q95`, `split_half_method`, `n_split_half_valid`, `split_half_mean`, `split_half_median`, `split_half_q05`, `split_half_q95`, `split_half_fraction_positive`, `split_half_fraction_above_0_30`, `error`, `edge_permutation_q_positive`, `loading_permutation_q_positive`, `conservative_preservation_class`, `structural_support_class`, `outcome_direction_concordant`, `joint_preservation_outcome_class`

#### `results/tables/cross_cohort_module_preservation_outcome_synthesis.csv`

- Rows: 12
- SHA-256: `0da091252528976a82d7c13d000a541752af577911b28b3b35fe851a3463b742`
- Columns: `cohort`, `module_label`, `n_samples`, `n_genes_available`, `pc1_variance_explained`, `frozen_loading_sign_concordance`, `frozen_loading_spearman`, `mean_signed_pairwise_correlation`, `median_signed_pairwise_correlation`, `mean_absolute_pairwise_correlation`, `signed_mean_vs_canine_weighted_correlation`, `signed_mean_vs_human_pc1_correlation`, `subset_method`, `n_subset_iterations`, `subset_reliability_mean`, `subset_reliability_median`, `subset_reliability_q05`, `subset_reliability_min`, `n_random_valid`, `signed_corr_random_mean`, `signed_corr_random_q95`, `signed_corr_random_percentile`, `signed_corr_empirical_p`, `pc1_variance_random_mean`, `pc1_variance_random_q95`, `pc1_variance_random_percentile`, `pc1_variance_empirical_p`, `error`, `signed_corr_empirical_q`, `pc1_variance_empirical_q`, `structural_support_class`, `target_os_hr_per_sd`, `target_os_p`, `target_os_q`, `target_os_fixed_c_index`, `gse21257_metastasis_auc`, `gse21257_metastasis_or_per_sd`, `gse21257_permutation_q`, `gse39055_rfs_hr_per_sd`, `gse39055_rfs_p`, `gse39055_rfs_q`, `gse39055_rfs_fixed_c_index`, `target_direction_concordant`, `gse21257_direction_concordant`, `gse39055_direction_concordant`, `n_human_settings_direction_concordant`, `outcome_direction_concordant`, `joint_preservation_outcome_class`

#### `results/tables/cross_cohort_module_representation_preservation.csv`

- Rows: 16
- SHA-256: `430428f52a2a66f8b87ee9258cbc331fba2456266fece84fa01b8cccbea0b983`
- Columns: `cohort`, `module_label`, `n_samples`, `n_genes_available`, `pc1_variance_explained`, `frozen_loading_sign_concordance`, `frozen_loading_spearman`, `mean_signed_pairwise_correlation`, `median_signed_pairwise_correlation`, `mean_absolute_pairwise_correlation`, `signed_mean_vs_canine_weighted_correlation`, `signed_mean_vs_human_pc1_correlation`, `subset_method`, `n_subset_iterations`, `subset_reliability_mean`, `subset_reliability_median`, `subset_reliability_q05`, `subset_reliability_min`, `n_random_valid`, `signed_corr_random_mean`, `signed_corr_random_q95`, `signed_corr_random_percentile`, `signed_corr_empirical_p`, `pc1_variance_random_mean`, `pc1_variance_random_q95`, `pc1_variance_random_percentile`, `pc1_variance_empirical_p`, `error`, `signed_corr_empirical_q`, `pc1_variance_empirical_q`

---

## `scripts/31_gse39055_assay_quality_diagnostic_v2.py`

- SHA-256: `758832ea7c193fa73be07048a6f97b9297f57b00518e322e1a3f6179cd8faf66`

### Relevant functions

#### `build_probe_quality`

- Lines: 340-373
- Signals: `mean, std`

```python
def build_probe_quality(
    values_probe_by_sample: pd.DataFrame,
    detection_probe_by_sample: pd.DataFrame,
    annotation: pd.DataFrame,
) -> pd.DataFrame:
    quality = pd.DataFrame(
        {
            "probe_id": values_probe_by_sample.index.astype(str),
            "expression_mean": values_probe_by_sample.mean(axis=1).values,
            "expression_median": values_probe_by_sample.median(axis=1).values,
            "expression_variance": values_probe_by_sample.var(axis=1).values,
            "expression_sd": values_probe_by_sample.std(axis=1).values,
            "detection_p_median": detection_probe_by_sample.median(axis=1).values,
            "detection_p_mean": detection_probe_by_sample.mean(axis=1).values,
            "detected_fraction_p_lt_0_01": (
                detection_probe_by_sample.lt(
                    DETECTION_THRESHOLD_STRICT
                ).mean(axis=1).values
            ),
            "detected_fraction_p_lt_0_05": (
                detection_probe_by_sample.lt(
                    DETECTION_THRESHOLD_RELAXED
                ).mean(axis=1).values
            ),
            "n_samples": values_probe_by_sample.shape[1],
        }
    )

    quality = quality.merge(annotation, on="probe_id", how="left")
    quality["gene_symbol"] = (
        quality["gene_symbol"].fillna("").astype(str).str.upper()
    )
    quality["unambiguous_gene_symbol"] = quality["gene_symbol"].ne("")
    return quality
```

#### `zscore_columns`

- Lines: 456-462
- Signals: `mean, std`

```python
def zscore_columns(expression: pd.DataFrame) -> pd.DataFrame:
    x = expression.apply(pd.to_numeric, errors="coerce")
    x = x.replace([np.inf, -np.inf], np.nan)
    x = x.fillna(x.median(axis=0))
    std = x.std(axis=0).replace(0, np.nan)
    z = (x - x.mean(axis=0)) / std
    return z.loc[:, z.notna().all(axis=0)]
```

#### `zscore_series`

- Lines: 465-470
- Signals: `mean, std`

```python
def zscore_series(values: pd.Series) -> pd.Series:
    values = pd.to_numeric(values, errors="coerce")
    std = values.std()
    if not np.isfinite(std) or std == 0:
        return pd.Series(np.nan, index=values.index, dtype=float)
    return (values - values.mean()) / std
```

#### `fit_cox_score`

- Lines: 548-610
- Signals: `std`

```python
def fit_cox_score(
    frame: pd.DataFrame,
    score_col: str,
) -> dict[str, Any]:
    data = frame[
        [TIME_COL, EVENT_COL, score_col]
    ].replace([np.inf, -np.inf], np.nan).dropna().copy()
    data = data[data[TIME_COL].gt(0)].copy()

    result: dict[str, Any] = {
        "n": data.shape[0],
        "events": (
            int(data[EVENT_COL].sum())
            if data.shape[0]
            else 0
        ),
        "hr_per_sd": np.nan,
        "ci_low": np.nan,
        "ci_high": np.nan,
        "p": np.nan,
        "fixed_direction_c_index": np.nan,
        "direction_concordant_with_canine": np.nan,
        "error": "",
    }

    if (
        data.shape[0] < 20
        or data[EVENT_COL].sum() < 5
        or data[score_col].std() == 0
    ):
        result["error"] = "insufficient_data"
        return result

    result["fixed_direction_c_index"] = fixed_direction_c_index(
        data[TIME_COL],
        data[EVENT_COL],
        data[score_col],
    )

    model = CoxPHFitter(penalizer=COX_PENALIZER)
    try:
        model.fit(
            data,
            duration_col=TIME_COL,
            event_col=EVENT_COL,
            fit_options={"max_steps": 500},
        )
        summary = model.summary.loc[score_col]
        result["hr_per_sd"] = float(summary["exp(coef)"])
        result["ci_low"] = float(
            summary["exp(coef) lower 95%"]
        )
        result["ci_high"] = float(
            summary["exp(coef) upper 95%"]
        )
        result["p"] = float(summary["p"])
        result["direction_concordant_with_canine"] = bool(
            result["hr_per_sd"] > 1.0
        )
    except Exception as exc:
        result["error"] = str(exc)[:500]

    return result
```

#### `compute_module_scores`

- Lines: 613-722
- Signals: `mean, weights, loadings, direction`

```python
def compute_module_scores(
    strict_weights: pd.DataFrame,
    values_probe_by_sample: pd.DataFrame,
    selection_strategies: dict[str, tuple[pd.DataFrame, float | None]],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    scores = pd.DataFrame(index=values_probe_by_sample.columns)
    coverage_rows: list[dict[str, Any]] = []

    weights = strict_weights.copy()
    weights["human_gene_symbol"] = (
        weights["human_gene_symbol"].astype(str).str.upper()
    )

    for strategy_name, (
        selection,
        min_detected_fraction,
    ) in selection_strategies.items():
        gene_expression, selected = selection_to_gene_expression(
            selection=selection,
            values_probe_by_sample=values_probe_by_sample,
            min_detected_fraction=min_detected_fraction,
        )
        selected_index = selected.set_index("gene_symbol")

        for module in PRIMARY_MODULES:
            module_weights = weights[
                weights["module_label"].eq(module)
            ].copy()
            module_weights = module_weights.drop_duplicates(
                "human_gene_symbol",
                keep="first",
            )
            requested_genes = module_weights[
                "human_gene_symbol"
            ].tolist()
            available_genes = [
                gene
                for gene in requested_genes
                if gene in gene_expression.columns
            ]

            n_requested = len(requested_genes)
            n_available = len(available_genes)
            coverage_fraction = (
                n_available / n_requested
                if n_requested
                else np.nan
            )

            detection_values = (
                selected_index.loc[
                    available_genes,
                    "detected_fraction_p_lt_0_01",
                ]
                if available_genes
                else pd.Series(dtype=float)
            )

            coverage_rows.append(
                {
                    "module_label": module,
                    "strategy": strategy_name,
                    "n_frozen_genes": n_requested,
                    "n_available_genes": n_available,
                    "coverage_fraction": coverage_fraction,
                    "median_gene_detected_fraction_p_lt_0_01": (
                        float(detection_values.median())
                        if not detection_values.empty
                        else np.nan
                    ),
                    "minimum_gene_detected_fraction_p_lt_0_01": (
                        float(detection_values.min())
                        if not detection_values.empty
                        else np.nan
                    ),
                    "available_genes": ";".join(available_genes),
                    "missing_or_filtered_genes": ";".join(
                        gene
                        for gene in requested_genes
                        if gene not in available_genes
                    ),
                }
            )

            if n_available < 3:
                continue

            z = zscore_columns(
                gene_expression[available_genes]
            )
            available_genes = list(z.columns)
            if len(available_genes) < 3:
                continue

            loadings = (
                module_weights
                .set_index("human_gene_symbol")
                .loc[available_genes, "risk_oriented_loading"]
            )
            signs = np.sign(
                pd.to_numeric(loadings, errors="coerce").fillna(0.0)
            )
            signs = signs.replace(0, 1)

            score = z.mul(signs, axis=1).mean(axis=1)
            scores[f"{module}__{strategy_name}"] = zscore_series(
                score
            )

    return scores, pd.DataFrame(coverage_rows)
```

#### `score_correlations_with_locked`

- Lines: 725-769
- Signals: `std, correlation`

```python
def score_correlations_with_locked(
    scores: pd.DataFrame,
    locked_scores: pd.DataFrame,
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []

    for module in PRIMARY_MODULES:
        locked_col = f"{module}__strict__signed_mean_z"
        if locked_col not in locked_scores.columns:
            continue

        for score_col in [
            column
            for column in scores.columns
            if column.startswith(f"{module}__")
        ]:
            frame = pd.concat(
                [
                    locked_scores[locked_col].rename("locked"),
                    scores[score_col].rename("diagnostic"),
                ],
                axis=1,
            ).dropna()

            correlation = (
                float(frame["locked"].corr(frame["diagnostic"]))
                if (
                    frame.shape[0] >= 5
                    and frame["locked"].std() > 0
                    and frame["diagnostic"].std() > 0
                )
                else np.nan
            )
            rows.append(
                {
                    "module_label": module,
                    "strategy": score_col.replace(
                        f"{module}__",
                        "",
                    ),
                    "score_correlation_with_locked": correlation,
                }
            )

    return pd.DataFrame(rows)
```

#### `run_score_rfs_sensitivity`

- Lines: 772-808
- Signals: ``

```python
def run_score_rfs_sensitivity(
    clinical: pd.DataFrame,
    scores: pd.DataFrame,
    correlations: pd.DataFrame,
) -> pd.DataFrame:
    frame = clinical.join(scores, how="inner")
    rows: list[dict[str, Any]] = []

    correlation_index = correlations.set_index(
        ["module_label", "strategy"]
    )

    for score_col in scores.columns:
        module, strategy = score_col.split("__", 1)
        fit = fit_cox_score(frame, score_col)

        correlation = np.nan
        if (module, strategy) in correlation_index.index:
            correlation = correlation_index.loc[
                (module, strategy),
                "score_correlation_with_locked",
            ]

        rows.append(
            {
                "module_label": module,
                "strategy": strategy,
                "score_column": score_col,
                "score_correlation_with_locked": correlation,
                **fit,
            }
        )

    result = pd.DataFrame(rows)
    return result.sort_values(
        ["module_label", "strategy"]
    ).reset_index(drop=True)
```

#### `build_sample_quality`

- Lines: 811-879
- Signals: `mean, weights`

```python
def build_sample_quality(
    detection_probe_by_sample: pd.DataFrame,
    annotation: pd.DataFrame,
    best_detected: pd.DataFrame,
    strict_weights: pd.DataFrame,
) -> pd.DataFrame:
    sample_qc = pd.DataFrame(
        index=detection_probe_by_sample.columns
    )
    sample_qc["n_all_probes"] = detection_probe_by_sample.shape[0]
    sample_qc["n_detected_p_lt_0_01"] = (
        detection_probe_by_sample.lt(
            DETECTION_THRESHOLD_STRICT
        ).sum(axis=0)
    )
    sample_qc["fraction_detected_p_lt_0_01"] = (
        detection_probe_by_sample.lt(
            DETECTION_THRESHOLD_STRICT
        ).mean(axis=0)
    )
    sample_qc["n_detected_p_lt_0_05"] = (
        detection_probe_by_sample.lt(
            DETECTION_THRESHOLD_RELAXED
        ).sum(axis=0)
    )
    sample_qc["fraction_detected_p_lt_0_05"] = (
        detection_probe_by_sample.lt(
            DETECTION_THRESHOLD_RELAXED
        ).mean(axis=0)
    )
    sample_qc["median_detection_p_all_probes"] = (
        detection_probe_by_sample.median(axis=0)
    )

    weights = strict_weights.copy()
    weights["human_gene_symbol"] = (
        weights["human_gene_symbol"].astype(str).str.upper()
    )
    best_index = best_detected.set_index("gene_symbol")

    for module in PRIMARY_MODULES:
        genes = (
            weights[weights["module_label"].eq(module)]
            ["human_gene_symbol"]
            .drop_duplicates()
            .tolist()
        )
        available = [
            gene for gene in genes if gene in best_index.index
        ]
        probes = best_index.loc[available, "probe_id"].tolist()
        probes = [
            probe
            for probe in probes
            if probe in detection_probe_by_sample.index
        ]
        if not probes:
            continue

        sample_qc[
            f"{module}_best_detected_probe_fraction_p_lt_0_01"
        ] = (
            detection_probe_by_sample.loc[probes]
            .lt(DETECTION_THRESHOLD_STRICT)
            .mean(axis=0)
        )

    sample_qc.index.name = "geo_sample_id"
    return sample_qc
```

#### `run_sample_quality_endpoint_diagnostics`

- Lines: 882-999
- Signals: `std`

```python
def run_sample_quality_endpoint_diagnostics(
    sample_qc: pd.DataFrame,
    clinical: pd.DataFrame,
    locked_scores: pd.DataFrame,
) -> pd.DataFrame:
    frame = clinical.join(sample_qc, how="inner").join(
        locked_scores.drop(columns=["cohort"], errors="ignore"),
        how="left",
    )

    predictor_cols = [
        column
        for column in sample_qc.columns
        if (
            "fraction_detected" in column
            or column.endswith("fraction_p_lt_0_01")
        )
    ]

    rows: list[dict[str, Any]] = []

    for predictor in predictor_cols:
        frame[f"{predictor}__z"] = zscore_series(
            frame[predictor]
        )
        fit = fit_cox_score(
            frame,
            f"{predictor}__z",
        )

        event_values = frame.loc[
            frame[EVENT_COL].eq(1),
            predictor,
        ].dropna()
        censor_values = frame.loc[
            frame[EVENT_COL].eq(0),
            predictor,
        ].dropna()

        mann_whitney_p = np.nan
        if (
            event_values.shape[0] >= 3
            and censor_values.shape[0] >= 3
        ):
            mann_whitney_p = float(
                stats.mannwhitneyu(
                    event_values,
                    censor_values,
                    alternative="two-sided",
                ).pvalue
            )

        rows.append(
            {
                "diagnostic_type": "sample_quality_vs_endpoint",
                "predictor": predictor,
                "module_label": "",
                "spearman_rho": np.nan,
                "spearman_p": np.nan,
                "event_group_median": (
                    float(event_values.median())
                    if not event_values.empty
                    else np.nan
                ),
                "censored_group_median": (
                    float(censor_values.median())
                    if not censor_values.empty
                    else np.nan
                ),
                "mann_whitney_p": mann_whitney_p,
                **fit,
            }
        )

    for module in PRIMARY_MODULES:
        score_col = f"{module}__strict__signed_mean_z"
        if score_col not in frame.columns:
            continue

        for predictor in predictor_cols:
            pair = frame[[score_col, predictor]].dropna()
            rho = np.nan
            p_value = np.nan
            if (
                pair.shape[0] >= 5
                and pair[score_col].std() > 0
                and pair[predictor].std() > 0
            ):
                result = stats.spearmanr(
                    pair[score_col],
                    pair[predictor],
                )
                rho = float(result.statistic)
                p_value = float(result.pvalue)

            rows.append(
                {
                    "diagnostic_type": "sample_quality_vs_locked_score",
                    "predictor": predictor,
                    "module_label": module,
                    "spearman_rho": rho,
                    "spearman_p": p_value,
                    "event_group_median": np.nan,
                    "censored_group_median": np.nan,
                    "mann_whitney_p": np.nan,
                    "n": pair.shape[0],
                    "events": np.nan,
                    "hr_per_sd": np.nan,
                    "ci_low": np.nan,
                    "ci_high": np.nan,
                    "p": np.nan,
                    "fixed_direction_c_index": np.nan,
                    "direction_concordant_with_canine": np.nan,
                    "error": "",
                }
            )

    return pd.DataFrame(rows)
```

#### `run_canonical_gene_audit`

- Lines: 1060-1188
- Signals: `weights, loadings, direction`

```python
def run_canonical_gene_audit(
    panel: pd.DataFrame,
    probe_quality: pd.DataFrame,
    values_probe_by_sample: pd.DataFrame,
    clinical: pd.DataFrame,
    selection_comparison: pd.DataFrame,
    strict_weights: pd.DataFrame,
) -> pd.DataFrame:
    comparison = selection_comparison.set_index("gene_symbol")
    weights = strict_weights.copy()
    weights["human_gene_symbol"] = (
        weights["human_gene_symbol"].astype(str).str.upper()
    )
    frozen_loading_map = (
        weights.groupby("human_gene_symbol")[
            "risk_oriented_loading"
        ]
        .first()
        .to_dict()
    )

    frame_base = clinical.copy()
    rows: list[dict[str, Any]] = []

    for _, panel_row in panel.iterrows():
        gene = panel_row["gene_symbol"]
        role = panel_row["panel_role"]

        probes = probe_quality[
            probe_quality["gene_symbol"].eq(gene)
        ].copy()

        for _, probe_row in probes.iterrows():
            probe_id = probe_row["probe_id"]
            if probe_id not in values_probe_by_sample.index:
                continue

            score = zscore_series(
                values_probe_by_sample.loc[probe_id]
            )
            frame = frame_base.join(
                score.rename("probe_expression_z"),
                how="inner",
            )
            fit = fit_cox_score(
                frame,
                "probe_expression_z",
            )

            max_probe = (
                comparison.loc[
                    gene,
                    "max_variance_probe_id",
                ]
                if gene in comparison.index
                else np.nan
            )
            best_probe = (
                comparison.loc[
                    gene,
                    "best_detected_probe_id",
                ]
                if gene in comparison.index
                else np.nan
            )

            rows.append(
                {
                    "gene_symbol": gene,
                    "panel_role": role,
                    "probe_id": probe_id,
                    "selected_by_max_variance": (
                        probe_id == max_probe
                    ),
                    "selected_by_best_detected": (
                        probe_id == best_probe
                    ),
                    "frozen_risk_oriented_loading": (
                        frozen_loading_map.get(gene, np.nan)
                    ),
                    "expression_variance": probe_row[
                        "expression_variance"
                    ],
                    "detected_fraction_p_lt_0_01": probe_row[
                        "detected_fraction_p_lt_0_01"
                    ],
                    "detected_fraction_p_lt_0_05": probe_row[
                        "detected_fraction_p_lt_0_05"
                    ],
                    "detection_p_median": probe_row[
                        "detection_p_median"
                    ],
                    **fit,
                }
            )

    result = pd.DataFrame(rows)

    expected_notes = {
        "MKI67": (
            "Proliferation marker; higher mRNA is a descriptive "
            "aggressiveness sanity check, not a definitive validation control."
        ),
        "TOP2A": (
            "Proliferation marker; higher mRNA is a descriptive "
            "aggressiveness sanity check."
        ),
        "BIRC5": (
            "Cell-cycle/survival marker; higher mRNA is a descriptive "
            "aggressiveness sanity check."
        ),
        "UBE2C": (
            "Cell-cycle marker; higher mRNA is a descriptive "
            "aggressiveness sanity check."
        ),
        "EZR": (
            "Ezrin evidence is largely protein/localization based; "
            "mRNA direction is not a definitive positive control."
        ),
    }
    if not result.empty:
        result["interpretation_note"] = result[
            "gene_symbol"
        ].map(expected_notes).fillna(
            "Frozen M34-loading diagnostic; no independent expected "
            "outcome direction is imposed."
        )

    return result
```

#### `build_module_summary`

- Lines: 1191-1325
- Signals: `mean, direction, coverage`

```python
def build_module_summary(
    coverage: pd.DataFrame,
    rfs: pd.DataFrame,
    locked_rfs: pd.DataFrame,
) -> pd.DataFrame:
    locked = locked_rfs.set_index("module_label")
    rows: list[dict[str, Any]] = []

    for module in PRIMARY_MODULES:
        module_rfs = rfs[
            rfs["module_label"].eq(module)
        ].copy()
        module_coverage = coverage[
            coverage["module_label"].eq(module)
        ].copy()

        locked_hr = (
            locked.loc[module, "hr_per_sd"]
            if module in locked.index
            else np.nan
        )
        locked_direction = (
            bool(locked_hr > 1)
            if np.isfinite(locked_hr)
            else np.nan
        )

        valid_directions = module_rfs[
            "direction_concordant_with_canine"
        ].dropna()
        valid_correlations = module_rfs[
            "score_correlation_with_locked"
        ].dropna()

        rows.append(
            {
                "module_label": module,
                "locked_hr_per_sd": locked_hr,
                "locked_direction_concordant_with_canine": locked_direction,
                "n_diagnostic_score_strategies": module_rfs.shape[0],
                "n_direction_concordant_strategies": int(
                    valid_directions.sum()
                ),
                "fraction_direction_concordant_strategies": (
                    float(valid_directions.mean())
                    if not valid_directions.empty
                    else np.nan
                ),
                "minimum_score_correlation_with_locked": (
                    float(valid_correlations.min())
                    if not valid_correlations.empty
                    else np.nan
                ),
                "median_score_correlation_with_locked": (
                    float(valid_correlations.median())
                    if not valid_correlations.empty
                    else np.nan
                ),
                "best_detected_50_coverage_fraction": (
                    module_coverage.loc[
                        module_coverage["strategy"].eq(
                            "best_detected_p01_ge50"
                        ),
                        "coverage_fraction",
                    ].iloc[0]
                    if (
                        module_coverage["strategy"].eq(
                            "best_detected_p01_ge50"
                        ).any()
                    )
                    else np.nan
                ),
                "best_detected_80_coverage_fraction": (
                    module_coverage.loc[
                        module_coverage["strategy"].eq(
                            "best_detected_p01_ge80"
                        ),
                        "coverage_fraction",
                    ].iloc[0]
                    if (
                        module_coverage["strategy"].eq(
                            "best_detected_p01_ge80"
                        ).any()
                    )
                    else np.nan
                ),
                "diagnostic_interpretation": "",
            }
        )

    summary = pd.DataFrame(rows)

    def interpretation(row: pd.Series) -> str:
        fraction = row[
            "fraction_direction_concordant_strategies"
        ]
        minimum_corr = row[
            "minimum_score_correlation_with_locked"
        ]
        coverage_50 = row[
            "best_detected_50_coverage_fraction"
        ]

        if (
            np.isfinite(fraction)
            and fraction == 1.0
            and np.isfinite(minimum_corr)
            and minimum_corr >= 0.90
            and np.isfinite(coverage_50)
            and coverage_50 >= 0.70
        ):
            return (
                "Direction and score representation are stable across "
                "detection-aware rules."
            )
        if (
            np.isfinite(fraction)
            and fraction == 0.0
            and np.isfinite(minimum_corr)
            and minimum_corr >= 0.90
        ):
            return (
                "Discordant direction persists across detection-aware "
                "rules despite high score concordance."
            )
        return (
            "Assay-rule sensitivity is material or coverage is limited; "
            "interpret GSE39055 cautiously."
        )

    summary["diagnostic_interpretation"] = summary.apply(
        interpretation,
        axis=1,
    )
    return summary
```

#### `interpretation`

- Lines: 1283-1319
- Signals: `direction, coverage`

```python
def interpretation(row: pd.Series) -> str:
        fraction = row[
            "fraction_direction_concordant_strategies"
        ]
        minimum_corr = row[
            "minimum_score_correlation_with_locked"
        ]
        coverage_50 = row[
            "best_detected_50_coverage_fraction"
        ]

        if (
            np.isfinite(fraction)
            and fraction == 1.0
            and np.isfinite(minimum_corr)
            and minimum_corr >= 0.90
            and np.isfinite(coverage_50)
            and coverage_50 >= 0.70
        ):
            return (
                "Direction and score representation are stable across "
                "detection-aware rules."
            )
        if (
            np.isfinite(fraction)
            and fraction == 0.0
            and np.isfinite(minimum_corr)
            and minimum_corr >= 0.90
        ):
            return (
                "Discordant direction persists across detection-aware "
                "rules despite high score concordance."
            )
        return (
            "Assay-rule sensitivity is material or coverage is limited; "
            "interpret GSE39055 cautiously."
        )
```

#### `write_readme`

- Lines: 1328-1364
- Signals: `weights, direction`

```python
def write_readme() -> None:
    text = f"""GSE39055 assay-quality diagnostic
Script version: {SCRIPT_VERSION}

Purpose
-------
This script audits the FFPE WG-DASL assay layer without changing the frozen
program definitions or the locked primary analyses.

Data used
---------
- GEO sample-level normalized VALUE
- GEO sample-level Detection PVal
- GPL14951 probe-to-gene annotation
- Frozen strict canine-to-human genes and risk-oriented signs

Outcome-blind probe rules
-------------------------
1. Highest-variance probe per unambiguous gene.
2. Best-detected probe per unambiguous gene.
3. Highest-variance probe filtered to Detection PVal < 0.01 in at least 50% of samples.
4. Best-detected probe filtered to Detection PVal < 0.01 in at least 50% of samples.
5. Best-detected probe filtered to Detection PVal < 0.01 in at least 80% of samples.

Interpretation restriction
--------------------------
RFS associations under alternative assay rules are diagnostic sensitivities.
They do not replace script 26, change frozen weights, reverse score direction,
or reopen the locked evidence hierarchy from script 29.

Canonical-gene restriction
--------------------------
MKI67, TOP2A, BIRC5, UBE2C, and EZR are descriptive assay-direction checks.
They are not treated as gold-standard positive controls, and no result is used
to select or orient a module.
"""
    OUTPUT_README.write_text(text, encoding="utf-8")
```

#### `create_manifest`

- Lines: 1367-1409
- Signals: `loadings`

```python
def create_manifest(
    input_paths: list[Path],
    output_paths: list[Path],
    detection_columns: dict[str, str],
) -> None:
    payload = {
        "script_version": SCRIPT_VERSION,
        "created_utc": utc_now_iso(),
        "detection_thresholds": {
            "strict_p": DETECTION_THRESHOLD_STRICT,
            "relaxed_p": DETECTION_THRESHOLD_RELAXED,
            "minimum_fraction_50": MIN_DETECTED_FRACTION_50,
            "minimum_fraction_80": MIN_DETECTED_FRACTION_80,
        },
        "detection_columns_by_sample": detection_columns,
        "guardrails": [
            "No outcome-guided probe selection.",
            "No change to frozen genes, loadings, directions, or tiers.",
            "No replacement of the locked script 26 primary analysis.",
            "Canonical genes are descriptive assay checks only.",
        ],
        "inputs": {},
        "outputs": {},
    }

    for path in input_paths:
        if path.exists():
            payload["inputs"][path.name] = {
                "sha256": sha256_file(path),
                "size_bytes": path.stat().st_size,
            }

    for path in output_paths:
        if path.exists():
            payload["outputs"][path.name] = {
                "sha256": sha256_file(path),
                "size_bytes": path.stat().st_size,
            }

    OUTPUT_MANIFEST.write_text(
        json.dumps(payload, indent=2),
        encoding="utf-8",
    )
```

#### `main`

- Lines: 1412-1771
- Signals: `mean, weights, loadings, direction, coverage`

```python
def main() -> None:
    print("=" * 80)
    print("GSE39055 FFPE DASL assay-quality diagnostic")
    print("=" * 80)
    print(f"Script version: {SCRIPT_VERSION}")
    print(f"Project root: {PROJECT_ROOT}")
    print("")
    print("Design:")
    print("  Read sample-level GEO Detection PVal values.")
    print("  Compare highest-variance and best-detected probe selection.")
    print("  Reconstruct frozen scores under outcome-blind detection filters.")
    print("  Audit sample quality, canonical genes, and RFS direction stability.")
    print("  Preserve the locked script 26 and script 29 conclusions.")
    print("")

    strict_weights = read_required_csv(STRICT_WEIGHTS_FILE)
    clinical = read_required_csv(CLINICAL_FILE, index_col=0)
    locked_scores = read_required_csv(
        LOCKED_SCORES_FILE,
        index_col=0,
    )
    locked_rfs = read_required_csv(LOCKED_RFS_FILE)

    strict_weights["human_gene_symbol"] = (
        strict_weights["human_gene_symbol"]
        .astype(str)
        .str.upper()
    )

    gse = load_gse39055()
    (
        values_probe_by_sample,
        detection_probe_by_sample,
        detection_columns,
    ) = extract_probe_matrices(gse)
    annotation = build_annotation(gse)

    common_samples = (
        values_probe_by_sample.columns
        .intersection(clinical.index)
        .intersection(locked_scores.index)
    )
    values_probe_by_sample = values_probe_by_sample[
        common_samples
    ]
    detection_probe_by_sample = detection_probe_by_sample[
        common_samples
    ]
    clinical = clinical.loc[common_samples].copy()
    locked_scores = locked_scores.loc[common_samples].copy()

    probe_quality = build_probe_quality(
        values_probe_by_sample=values_probe_by_sample,
        detection_probe_by_sample=detection_probe_by_sample,
        annotation=annotation,
    )
    probe_quality.to_csv(OUTPUT_PROBE_QC, index=False)

    (
        max_variance,
        best_detected,
        selection_comparison,
    ) = select_probes(probe_quality)
    selection_comparison = add_locked_probe_map_comparison(
        selection_comparison
    )
    selection_comparison.to_csv(
        OUTPUT_SELECTION,
        index=False,
    )

    selection_strategies = {
        "max_variance_unfiltered": (
            max_variance,
            None,
        ),
        "max_variance_p01_ge50": (
            max_variance,
            MIN_DETECTED_FRACTION_50,
        ),
        "best_detected_unfiltered": (
            best_detected,
            None,
        ),
        "best_detected_p01_ge50": (
            best_detected,
            MIN_DETECTED_FRACTION_50,
        ),
        "best_detected_p01_ge80": (
            best_detected,
            MIN_DETECTED_FRACTION_80,
        ),
    }

    scores, coverage = compute_module_scores(
        strict_weights=strict_weights,
        values_probe_by_sample=values_probe_by_sample,
        selection_strategies=selection_strategies,
    )
    coverage.to_csv(
        OUTPUT_MODULE_COVERAGE,
        index=False,
    )
    scores.to_csv(OUTPUT_SCORES)

    correlations = score_correlations_with_locked(
        scores=scores,
        locked_scores=locked_scores,
    )
    rfs = run_score_rfs_sensitivity(
        clinical=clinical,
        scores=scores,
        correlations=correlations,
    )
    rfs.to_csv(OUTPUT_RFS, index=False)

    sample_qc = build_sample_quality(
        detection_probe_by_sample=detection_probe_by_sample,
        annotation=annotation,
        best_detected=best_detected,
        strict_weights=strict_weights,
    )
    sample_qc.to_csv(OUTPUT_SAMPLE_QC)

    endpoint_qc = run_sample_quality_endpoint_diagnostics(
        sample_qc=sample_qc,
        clinical=clinical,
        locked_scores=locked_scores,
    )
    endpoint_qc.to_csv(
        OUTPUT_ENDPOINT_QC,
        index=False,
    )

    canonical_panel = build_canonical_panel(
        strict_weights
    )
    canonical = run_canonical_gene_audit(
        panel=canonical_panel,
        probe_quality=probe_quality,
        values_probe_by_sample=values_probe_by_sample,
        clinical=clinical,
        selection_comparison=selection_comparison,
        strict_weights=strict_weights,
    )
    canonical.to_csv(
        OUTPUT_CANONICAL,
        index=False,
    )

    summary = build_module_summary(
        coverage=coverage,
        rfs=rfs,
        locked_rfs=locked_rfs,
    )
    summary.to_csv(OUTPUT_SUMMARY, index=False)

    write_readme()

    input_paths = [
        SOFT_FILE,
        STRICT_WEIGHTS_FILE,
        CLINICAL_FILE,
        LOCKED_SCORES_FILE,
        LOCKED_PROBE_MAP_FILE,
        LOCKED_RFS_FILE,
    ]
    output_paths = [
        OUTPUT_PROBE_QC,
        OUTPUT_SELECTION,
        OUTPUT_MODULE_COVERAGE,
        OUTPUT_SCORES,
        OUTPUT_RFS,
        OUTPUT_CANONICAL,
        OUTPUT_SAMPLE_QC,
        OUTPUT_ENDPOINT_QC,
        OUTPUT_SUMMARY,
        OUTPUT_README,
    ]
    create_manifest(
        input_paths=input_paths,
        output_paths=output_paths,
        detection_columns=detection_columns,
    )

    print("")
    print("=" * 80)
    print("Detection P-value extraction")
    print("=" * 80)
    unique_detection_columns = sorted(
        set(detection_columns.values())
    )
    print(f"Samples parsed: {len(detection_columns)}")
    print(f"Detection columns: {unique_detection_columns}")
    print(
        f"Probe-by-sample matrix: "
        f"{values_probe_by_sample.shape[0]} probes x "
        f"{values_probe_by_sample.shape[1]} samples"
    )

    print("")
    print("=" * 80)
    print("Sample-level assay quality")
    print("=" * 80)
    print(
        sample_qc[
            [
                "n_detected_p_lt_0_01",
                "fraction_detected_p_lt_0_01",
                "n_detected_p_lt_0_05",
                "fraction_detected_p_lt_0_05",
                "median_detection_p_all_probes",
            ]
        ]
        .describe()
        .to_string()
    )

    print("")
    print("=" * 80)
    print("Selected-probe detectability comparison")
    print("=" * 80)
    comparison_summary = pd.DataFrame(
        {
            "metric": [
                "genes_compared",
                "same_probe_fraction",
                "locked_matches_recomputed_max_variance_fraction",
                "median_max_variance_detected_fraction_p01",
                "median_best_detected_detected_fraction_p01",
            ],
            "value": [
                selection_comparison.shape[0],
                selection_comparison[
                    "same_selected_probe"
                ].mean(),
                selection_comparison[
                    "locked_matches_recomputed_max_variance"
                ].dropna().mean(),
                selection_comparison[
                    "max_variance_detected_fraction_p_lt_0_01"
                ].median(),
                selection_comparison[
                    "best_detected_detected_fraction_p_lt_0_01"
                ].median(),
            ],
        }
    )
    print(comparison_summary.to_string(index=False))

    print("")
    print("=" * 80)
    print("Detection-aware frozen-module coverage")
    print("=" * 80)
    print(
        coverage[
            [
                "module_label",
                "strategy",
                "n_frozen_genes",
                "n_available_genes",
                "coverage_fraction",
                "median_gene_detected_fraction_p_lt_0_01",
                "minimum_gene_detected_fraction_p_lt_0_01",
            ]
        ].to_string(index=False)
    )

    print("")
    print("=" * 80)
    print("Detection-aware RFS direction sensitivity")
    print("=" * 80)
    print(
        rfs[
            [
                "module_label",
                "strategy",
                "score_correlation_with_locked",
                "n",
                "events",
                "hr_per_sd",
                "ci_low",
                "ci_high",
                "p",
                "fixed_direction_c_index",
                "direction_concordant_with_canine",
            ]
        ].to_string(index=False)
    )

    print("")
    print("=" * 80)
    print("Module-level assay diagnostic summary")
    print("=" * 80)
    print(summary.to_string(index=False))

    print("")
    print("=" * 80)
    print("Canonical and M34-loading probe audit")
    print("=" * 80)
    canonical_display = canonical[
        canonical["selected_by_max_variance"]
        | canonical["selected_by_best_detected"]
    ].copy()
    if canonical_display.empty:
        print("No selected canonical probes were available.")
    else:
        print(
            canonical_display[
                [
                    "gene_symbol",
                    "panel_role",
                    "probe_id",
                    "selected_by_max_variance",
                    "selected_by_best_detected",
                    "detected_fraction_p_lt_0_01",
                    "detection_p_median",
                    "hr_per_sd",
                    "p",
                    "fixed_direction_c_index",
                ]
            ].to_string(index=False)
        )

    print("")
    print("=" * 80)
    print("Assay-quality endpoint diagnostics")
    print("=" * 80)
    print(
        endpoint_qc[
            endpoint_qc["diagnostic_type"].eq(
                "sample_quality_vs_endpoint"
            )
        ][
            [
                "predictor",
                "event_group_median",
                "censored_group_median",
                "mann_whitney_p",
                "hr_per_sd",
                "p",
            ]
        ].to_string(index=False)
    )

    print("")
    print("=" * 80)
    print("Interpretation guardrails")
    print("=" * 80)
    print("This script is a diagnostic addendum and does not replace the locked script 26 analysis.")
    print("Probe selection and detection filters are outcome-blind.")
    print("Frozen genes, risk-oriented signs, weights, and validation tiers are unchanged.")
    print("Canonical-gene results are descriptive assay checks, not positive-control validation.")
    print("A stable discordant direction across detection-aware rules supports cohort/platform heterogeneity; an unstable direction indicates assay-rule sensitivity.")

    print("")
    print("Saved:")
    for path in output_paths + [OUTPUT_MANIFEST]:
        print(path)
    print("Done.")
```

### Relevant standalone lines

| Line | Signals | Code |
|---:|---|---|
| 348 | mean | `"expression_mean": values_probe_by_sample.mean(axis=1).values,` |
| 351 | std | `"expression_sd": values_probe_by_sample.std(axis=1).values,` |
| 353 | mean | `"detection_p_mean": detection_probe_by_sample.mean(axis=1).values,` |
| 357 | mean | `).mean(axis=1).values` |
| 362 | mean | `).mean(axis=1).values` |
| 460 | std | `std = x.std(axis=0).replace(0, np.nan)` |
| 461 | mean | `z = (x - x.mean(axis=0)) / std` |
| 467 | std | `std = values.std()` |
| 470 | mean | `return (values - values.mean()) / std` |
| 576 | std | `or data[score_col].std() == 0` |
| 621 | weights | `weights = strict_weights.copy()` |
| 622 | weights | `weights["human_gene_symbol"] = (` |
| 623 | weights | `weights["human_gene_symbol"].astype(str).str.upper()` |
| 638 | weights | `module_weights = weights[` |
| 639 | weights | `weights["module_label"].eq(module)` |
| 707 | loadings | `loadings = (` |
| 712 | direction | `signs = np.sign(` |
| 713 | loadings | `pd.to_numeric(loadings, errors="coerce").fillna(0.0)` |
| 717 | mean | `score = z.mul(signs, axis=1).mean(axis=1)` |
| 750 | correlation | `float(frame["locked"].corr(frame["diagnostic"]))` |
| 753 | std | `and frame["locked"].std() > 0` |
| 754 | std | `and frame["diagnostic"].std() > 0` |
| 829 | mean | `).mean(axis=0)` |
| 839 | mean | `).mean(axis=0)` |
| 845 | weights | `weights = strict_weights.copy()` |
| 846 | weights | `weights["human_gene_symbol"] = (` |
| 847 | weights | `weights["human_gene_symbol"].astype(str).str.upper()` |
| 853 | weights | `weights[weights["module_label"].eq(module)]` |
| 875 | mean | `.mean(axis=0)` |
| 967 | std | `and pair[score_col].std() > 0` |
| 968 | std | `and pair[predictor].std() > 0` |
| 1069 | weights | `weights = strict_weights.copy()` |
| 1070 | weights | `weights["human_gene_symbol"] = (` |
| 1071 | weights | `weights["human_gene_symbol"].astype(str).str.upper()` |
| 1074 | weights | `weights.groupby("human_gene_symbol")[` |
| 1177 | direction | `"mRNA direction is not a definitive positive control."` |
| 1184 | loadings | `"Frozen M34-loading diagnostic; no independent expected "` |
| 1185 | direction | `"outcome direction is imposed."` |
| 1192 | coverage | `coverage: pd.DataFrame,` |
| 1203 | coverage | `module_coverage = coverage[` |
| 1204 | coverage | `coverage["module_label"].eq(module)` |
| 1235 | mean | `float(valid_directions.mean())` |
| 1303 | direction | `"Direction and score representation are stable across "` |
| 1313 | direction | `"Discordant direction persists across detection-aware "` |
| 1317 | coverage | `"Assay-rule sensitivity is material or coverage is limited; "` |
| 1355 | weights, direction | `They do not replace script 26, change frozen weights, reverse score direction,` |
| 1360 | direction | `MKI67, TOP2A, BIRC5, UBE2C, and EZR are descriptive assay-direction checks.` |
| 1384 | loadings | `"No change to frozen genes, loadings, directions, or tiers.",` |
| 1423 | direction | `print("  Audit sample quality, canonical genes, and RFS direction stability.")` |
| 1506 | coverage | `scores, coverage = compute_module_scores(` |
| 1511 | coverage | `coverage.to_csv(` |
| 1563 | coverage | `coverage=coverage,` |
| 1647 | mean | `].mean(),` |
| 1650 | mean | `].dropna().mean(),` |
| 1664 | coverage | `print("Detection-aware frozen-module coverage")` |
| 1667 | coverage | `coverage[` |
| 1682 | direction | `print("Detection-aware RFS direction sensitivity")` |
| 1710 | loadings | `print("Canonical and M34-loading probe audit")` |
| 1763 | weights | `print("Frozen genes, risk-oriented signs, weights, and validation tiers are unchanged.")` |
| 1765 | direction | `print("A stable discordant direction across detection-aware rules supports cohort/platform heterogeneity; an unstable direction indicates assay-rule sensitivity.")` |

### Referenced result-table schemas

#### `results/tables/GSE238110_frozen_transfer_gene_weights_strict.csv`

- Rows: 321
- SHA-256: `4f065aa4c4edf117a0c74015840d2b4b2347929f172cd517e1818ba0f6163b91`
- Columns: `module_label`, `canine_gene`, `canine_gene_symbol`, `human_gene_symbol`, `ortholog_qc_status`, `raw_pca_loading`, `risk_oriented_loading`, `absolute_risk_loading`, `is_strict_mapping`, `is_broad_mapping`, `human_symbol_duplicate_count`, `normalized_abs_sum_weight`

#### `results/tables/GSE39055_RFS_primary_frozen_program_validation.csv`

- Rows: 4
- SHA-256: `1945cca593030f57cbaf6eef3006391b4fc3ad2ff334d99e34cc81436533951e`
- Columns: `module_label`, `score_column`, `n`, `events`, `hr_per_sd`, `ci_low`, `ci_high`, `primary_p`, `fixed_score_c_index`, `fixed_score_c_index_ci_low`, `fixed_score_c_index_ci_high`, `c_index_bootstrap_valid`, `ph_test_p`, `bootstrap_hr_ci_low`, `bootstrap_hr_ci_high`, `bootstrap_probability_coef_positive`, `hr_bootstrap_valid`, `permutation_c_index_p_one_sided`, `permutation_c_index_p_two_sided`, `epsilon_zero_time_hr_per_sd`, `epsilon_zero_time_ci_low`, `epsilon_zero_time_ci_high`, `epsilon_zero_time_p`, `error`, `q_within_gse39055`, `permutation_c_index_q_bh`, `gse39055_support_class`

#### `results/tables/GSE39055_assay_quality_endpoint_diagnostics.csv`

- Rows: 30
- SHA-256: `16d889f92e140adc9af0be4774b843bec778357cac5022674479ef18cc144e6b`
- Columns: `diagnostic_type`, `predictor`, `module_label`, `spearman_rho`, `spearman_p`, `event_group_median`, `censored_group_median`, `mann_whitney_p`, `n`, `events`, `hr_per_sd`, `ci_low`, `ci_high`, `p`, `fixed_direction_c_index`, `direction_concordant_with_canine`, `error`

#### `results/tables/GSE39055_assay_quality_module_summary.csv`

- Rows: 4
- SHA-256: `0fdc2fb790819356e2d460dc508b7aa1db31d31443bd73a4eca38d04eeb26a41`
- Columns: `module_label`, `locked_hr_per_sd`, `locked_direction_concordant_with_canine`, `n_diagnostic_score_strategies`, `n_direction_concordant_strategies`, `fraction_direction_concordant_strategies`, `minimum_score_correlation_with_locked`, `median_score_correlation_with_locked`, `best_detected_50_coverage_fraction`, `best_detected_80_coverage_fraction`, `diagnostic_interpretation`

#### `results/tables/GSE39055_canonical_gene_assay_audit.csv`

- Rows: 23
- SHA-256: `4bd5f7001d53d6bccb9fbe6d50c179d7d86ad51beb1d9352e572c50fcb08175c`
- Columns: `gene_symbol`, `panel_role`, `probe_id`, `selected_by_max_variance`, `selected_by_best_detected`, `frozen_risk_oriented_loading`, `expression_variance`, `detected_fraction_p_lt_0_01`, `detected_fraction_p_lt_0_05`, `detection_p_median`, `n`, `events`, `hr_per_sd`, `ci_low`, `ci_high`, `p`, `fixed_direction_c_index`, `direction_concordant_with_canine`, `error`, `interpretation_note`

#### `results/tables/GSE39055_detection_aware_RFS_sensitivity.csv`

- Rows: 16
- SHA-256: `4bf8dea35f0469e90272352765cd6476ef9fcf3bd669e95fd17f82deb2dbe478`
- Columns: `module_label`, `strategy`, `score_column`, `score_correlation_with_locked`, `n`, `events`, `hr_per_sd`, `ci_low`, `ci_high`, `p`, `fixed_direction_c_index`, `direction_concordant_with_canine`, `error`

#### `results/tables/GSE39055_detection_aware_module_coverage.csv`

- Rows: 20
- SHA-256: `b54c838913171e8c9c3d586100cb1593a24d0987dab0c8024092716249695ca0`
- Columns: `module_label`, `strategy`, `n_frozen_genes`, `n_available_genes`, `coverage_fraction`, `median_gene_detected_fraction_p_lt_0_01`, `minimum_gene_detected_fraction_p_lt_0_01`, `available_genes`, `missing_or_filtered_genes`

#### `results/tables/GSE39055_gene_probe_selection_comparison.csv`

- Rows: 20793
- SHA-256: `ade10b84086b292e10535a21e45e559513c1ed686f808bd50558b8566d4bea53`
- Columns: `max_variance_probe_id`, `max_variance_expression_mean`, `max_variance_expression_median`, `max_variance_expression_variance`, `max_variance_expression_sd`, `max_variance_detection_p_median`, `max_variance_detection_p_mean`, `max_variance_detected_fraction_p_lt_0_01`, `max_variance_detected_fraction_p_lt_0_05`, `max_variance_n_samples`, `gene_symbol`, `max_variance_unambiguous_gene_symbol`, `best_detected_probe_id`, `best_detected_expression_mean`, `best_detected_expression_median`, `best_detected_expression_variance`, `best_detected_expression_sd`, `best_detected_detection_p_median`, `best_detected_detection_p_mean`, `best_detected_detected_fraction_p_lt_0_01`, `best_detected_detected_fraction_p_lt_0_05`, `best_detected_n_samples`, `best_detected_unambiguous_gene_symbol`, `same_selected_probe`, `locked_probe_id`, `locked_matches_recomputed_max_variance`

#### `results/tables/GSE39055_probe_assay_quality.csv`

- Rows: 29377
- SHA-256: `045704a7132bc28a5e6cebbc888f38ba31190fae816e9f879c0c876b66232357`
- Columns: `probe_id`, `expression_mean`, `expression_median`, `expression_variance`, `expression_sd`, `detection_p_median`, `detection_p_mean`, `detected_fraction_p_lt_0_01`, `detected_fraction_p_lt_0_05`, `n_samples`, `gene_symbol`, `unambiguous_gene_symbol`

#### `results/tables/GSE39055_probe_to_gene_symbol_selected.csv`

- Rows: 20793
- SHA-256: `00a773be24c545aaed9f22fa2e35594e946925026e5cb37c64a572b53288df99`
- Columns: `probe_id`, `gene_symbol`, `probe_variance`

#### `results/tables/GSE39055_sample_assay_quality.csv`

- Rows: 37
- SHA-256: `6c89881ba9104a36810afa292ee8789427c3bda4bc6bcb94a75d386de08175a6`
- Columns: `geo_sample_id`, `n_all_probes`, `n_detected_p_lt_0_01`, `fraction_detected_p_lt_0_01`, `n_detected_p_lt_0_05`, `fraction_detected_p_lt_0_05`, `median_detection_p_all_probes`, `M34_best_detected_probe_fraction_p_lt_0_01`, `M11_best_detected_probe_fraction_p_lt_0_01`, `M24_best_detected_probe_fraction_p_lt_0_01`, `M40_best_detected_probe_fraction_p_lt_0_01`

---

## `scripts/42_score_ammons_single_cell_localization.py`

- SHA-256: `5f876e15aa10dca589e1a7878de26d34ca3e1a4116beeb43ef038b5b74cb641b`

### Relevant functions

#### `read_inputs`

- Lines: 170-224
- Signals: `weights`

```python
def read_inputs() -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    for path in [
        METADATA_FILE,
        EXPRESSION_FILE,
        PREFLIGHT_MANIFEST_FILE,
        STRICT_WEIGHTS_FILE,
    ]:
        if not path.exists():
            raise FileNotFoundError(f"Required file not found: {path}")

    metadata = pd.read_csv(
        METADATA_FILE,
        sep="\t",
        dtype=str,
        low_memory=False,
    )
    weights = pd.read_csv(STRICT_WEIGHTS_FILE)
    manifest = json.loads(
        PREFLIGHT_MANIFEST_FILE.read_text(encoding="utf-8")
    )

    if bool(manifest.get("outcome_loaded", True)):
        raise RuntimeError(
            "The single-cell preflight manifest does not confirm "
            "outcome_loaded=false."
        )

    required_metadata_columns = [
        CELL_ID_COLUMN,
        DOG_COLUMN,
        *ANNOTATION_LEVELS,
    ]
    missing = [
        column
        for column in required_metadata_columns
        if column not in metadata.columns
    ]
    if missing:
        raise ValueError(
            f"Required metadata columns are missing: {missing}"
        )

    weights["human_gene_symbol"] = (
        weights["human_gene_symbol"].astype(str).str.upper()
    )
    weights["risk_oriented_loading"] = pd.to_numeric(
        weights["risk_oriented_loading"],
        errors="coerce",
    )
    weights = weights[
        weights["module_label"].isin(ALL_MODULES)
        & weights["risk_oriented_loading"].notna()
    ].copy()

    return metadata, weights, manifest
```

#### `extract_selected_expression`

- Lines: 243-298
- Signals: `mean`

```python
def extract_selected_expression(
    requested_genes: set[str],
    expected_cells: int,
) -> tuple[list[str], np.ndarray, dict[str, int]]:
    values_by_gene: dict[str, list[np.ndarray]] = {}
    row_counts: dict[str, int] = {}

    with gzip.open(
        EXPRESSION_FILE,
        "rt",
        encoding="utf-8",
    ) as handle:
        handle.readline()

        for line_number, line in enumerate(handle, start=2):
            if not line.strip():
                continue

            gene_token, values_text = line.split("\t", 1)
            gene = clean_gene_symbol(gene_token)

            if gene not in requested_genes:
                continue

            values = np.fromstring(
                values_text,
                sep="\t",
                dtype=np.float32,
            )
            if values.size != expected_cells:
                raise RuntimeError(
                    f"Gene {gene} on line {line_number} has "
                    f"{values.size} values; expected {expected_cells}."
                )

            values_by_gene.setdefault(gene, []).append(values)
            row_counts[gene] = row_counts.get(gene, 0) + 1

    genes = sorted(values_by_gene)
    matrix = np.empty(
        (len(genes), expected_cells),
        dtype=np.float32,
    )

    for index, gene in enumerate(genes):
        arrays = values_by_gene[gene]
        if len(arrays) == 1:
            matrix[index] = arrays[0]
        else:
            matrix[index] = np.mean(
                np.stack(arrays, axis=0),
                axis=0,
                dtype=np.float32,
            )

    return genes, matrix, row_counts
```

#### `rank_gaussian_transform`

- Lines: 301-329
- Signals: `std`

```python
def rank_gaussian_transform(
    expression: np.ndarray,
) -> np.ndarray:
    transformed = np.empty_like(expression, dtype=np.float32)
    n_cells = expression.shape[1]

    for gene_index in range(expression.shape[0]):
        values = expression[gene_index].astype(float)
        finite = np.isfinite(values)

        if finite.sum() < 2 or np.nanstd(values[finite]) == 0:
            transformed[gene_index] = 0.0
            continue

        fill_value = float(np.nanmedian(values[finite]))
        values[~finite] = fill_value

        ranks = stats.rankdata(values, method="average")
        probabilities = (ranks - 0.5) / n_cells
        probabilities = np.clip(
            probabilities,
            RANK_GAUSSIAN_EPSILON,
            1.0 - RANK_GAUSSIAN_EPSILON,
        )
        transformed[gene_index] = stats.norm.ppf(
            probabilities
        ).astype(np.float32)

    return transformed
```

#### `build_score_weights`

- Lines: 332-404
- Signals: `weights, loadings`

```python
def build_score_weights(
    weights: pd.DataFrame,
    detected_genes: list[str],
) -> tuple[dict[str, dict[str, np.ndarray]], pd.DataFrame]:
    gene_index = {
        gene: index for index, gene in enumerate(detected_genes)
    }
    score_weights: dict[str, dict[str, np.ndarray]] = {}
    coverage_rows = []

    for module in ALL_MODULES:
        part = weights[
            weights["module_label"].eq(module)
        ].drop_duplicates(
            "human_gene_symbol",
            keep="first",
        )

        available = part[
            part["human_gene_symbol"].isin(gene_index)
        ].copy()
        indices = np.asarray(
            [
                gene_index[gene]
                for gene in available["human_gene_symbol"]
            ],
            dtype=int,
        )
        loadings = available[
            "risk_oriented_loading"
        ].to_numpy(dtype=float)

        positive_mask = loadings > 0
        negative_mask = loadings < 0

        score_weights[module] = {
            "indices": indices,
            "loadings": loadings,
            "positive_indices": indices[positive_mask],
            "positive_weights": loadings[positive_mask],
            "negative_indices": indices[negative_mask],
            "negative_weights": np.abs(loadings[negative_mask]),
        }

        coverage_rows.append(
            {
                "module_label": module,
                "n_frozen_genes": int(part.shape[0]),
                "n_detected_genes": int(available.shape[0]),
                "coverage_fraction": float(
                    available.shape[0] / part.shape[0]
                ),
                "n_positive_detected": int(positive_mask.sum()),
                "n_negative_detected": int(negative_mask.sum()),
                "signed_score_estimable": bool(
                    available.shape[0] >= 3
                ),
                "positive_component_estimable": bool(
                    positive_mask.sum() >= 2
                ),
                "negative_component_estimable": bool(
                    negative_mask.sum() >= 2
                ),
                "component_guardrail": (
                    "M34 positive component is descriptive because "
                    "only two positive-loading genes are available."
                    if module == "M34"
                    else ""
                ),
            }
        )

    return score_weights, pd.DataFrame(coverage_rows)
```

#### `weighted_component`

- Lines: 407-431
- Signals: `weights, matrix_product`

```python
def weighted_component(
    expression: np.ndarray,
    indices: np.ndarray,
    weights: np.ndarray,
) -> np.ndarray:
    if indices.size == 0:
        return np.full(
            expression.shape[1],
            np.nan,
            dtype=np.float32,
        )

    denominator = float(np.sum(np.abs(weights)))
    if denominator <= 0:
        return np.full(
            expression.shape[1],
            np.nan,
            dtype=np.float32,
        )

    return (
        weights.astype(np.float32)
        @ expression[indices]
        / denominator
    ).astype(np.float32)
```

#### `calculate_cell_scores`

- Lines: 434-468
- Signals: `loadings, direction`

```python
def calculate_cell_scores(
    transformed_expression: np.ndarray,
    score_weights: dict[str, dict[str, np.ndarray]],
) -> pd.DataFrame:
    score_columns: dict[str, np.ndarray] = {}

    for module, specification in score_weights.items():
        indices = specification["indices"]
        loadings = specification["loadings"]

        signed = weighted_component(
            transformed_expression,
            indices,
            loadings,
        )
        positive = weighted_component(
            transformed_expression,
            specification["positive_indices"],
            specification["positive_weights"],
        )
        negative_expression = weighted_component(
            transformed_expression,
            specification["negative_indices"],
            specification["negative_weights"],
        )

        score_columns[f"{module}_signed_risk_score"] = signed
        score_columns[
            f"{module}_positive_component_expression"
        ] = positive
        score_columns[
            f"{module}_negative_component_expression"
        ] = negative_expression

    return pd.DataFrame(score_columns)
```

#### `aggregate_scores`

- Lines: 504-537
- Signals: `mean`

```python
def aggregate_scores(
    cells: pd.DataFrame,
    group_columns: list[str],
    score_columns: list[str],
    minimum_cells: int,
) -> pd.DataFrame:
    grouped = cells.groupby(
        group_columns,
        dropna=False,
        observed=True,
    )

    counts = grouped.size().rename("n_cells")
    means = grouped[score_columns].mean()
    medians = grouped[score_columns].median()
    detection = grouped[
        [f"{module}_detected_gene_fraction" for module in ALL_MODULES]
    ].mean()

    means.columns = [f"{column}_mean" for column in means.columns]
    medians.columns = [
        f"{column}_median" for column in medians.columns
    ]
    detection.columns = [
        f"{column}_mean" for column in detection.columns
    ]

    result = pd.concat(
        [counts, means, medians, detection],
        axis=1,
    ).reset_index()
    return result[
        result["n_cells"].ge(minimum_cells)
    ].reset_index(drop=True)
```

#### `exact_sign_flip_test`

- Lines: 540-590
- Signals: `mean`

```python
def exact_sign_flip_test(
    differences: np.ndarray,
) -> tuple[float, float]:
    values = np.asarray(differences, dtype=float)
    values = values[np.isfinite(values)]

    if values.size == 0:
        return np.nan, np.nan

    observed = float(np.mean(values))
    n = values.size

    if n <= 20:
        statistics = []
        for signs in product([-1.0, 1.0], repeat=n):
            statistics.append(
                float(
                    np.mean(
                        values * np.asarray(signs, dtype=float)
                    )
                )
            )
        statistics_array = np.asarray(statistics)
        p = float(
            np.mean(
                np.abs(statistics_array)
                >= abs(observed) - 1e-12
            )
        )
    else:
        rng = np.random.default_rng(RANDOM_SEED)
        signs = rng.choice(
            [-1.0, 1.0],
            size=(100000, n),
        )
        statistics_array = np.mean(
            signs * values[None, :],
            axis=1,
        )
        p = float(
            (
                1
                + np.sum(
                    np.abs(statistics_array)
                    >= abs(observed)
                )
            )
            / (statistics_array.size + 1)
        )

    return observed, p
```

#### `bootstrap_mean_ci`

- Lines: 593-613
- Signals: `mean`

```python
def bootstrap_mean_ci(
    differences: np.ndarray,
) -> tuple[float, float]:
    values = np.asarray(differences, dtype=float)
    values = values[np.isfinite(values)]

    if values.size == 0:
        return np.nan, np.nan

    rng = np.random.default_rng(RANDOM_SEED)
    indices = rng.integers(
        0,
        values.size,
        size=(N_BOOTSTRAP, values.size),
    )
    means = values[indices].mean(axis=1)

    return (
        float(np.quantile(means, 0.025)),
        float(np.quantile(means, 0.975)),
    )
```

#### `targeted_tests`

- Lines: 673-749
- Signals: `loadings, direction`

```python
def targeted_tests(
    compartment_scores: pd.DataFrame,
    l1_scores: pd.DataFrame,
) -> pd.DataFrame:
    tests = []

    tests.append(
        paired_contrast(
            pseudobulk=compartment_scores,
            group_column="broad_compartment",
            group_a="osteoblast_lineage",
            group_b="immune_combined",
            score_column="M34_signed_risk_score_mean",
            contrast_name=(
                "M34 signed risk: osteoblast lineage versus immune"
            ),
            hypothesis_direction=(
                "Higher in osteoblast lineage"
            ),
        )
    )
    tests.append(
        paired_contrast(
            pseudobulk=compartment_scores,
            group_column="broad_compartment",
            group_a="immune_combined",
            group_b="osteoblast_lineage",
            score_column=(
                "M34_negative_component_expression_mean"
            ),
            contrast_name=(
                "M34 negative-loading expression: immune versus "
                "osteoblast lineage"
            ),
            hypothesis_direction="Higher in immune compartment",
        )
    )
    tests.append(
        paired_contrast(
            pseudobulk=l1_scores,
            group_column="celltype.l1",
            group_a="Osteoblast_cycling",
            group_b="Osteoblast",
            score_column="M40_signed_risk_score_mean",
            contrast_name=(
                "M40 signed risk: cycling versus non-cycling "
                "osteoblast"
            ),
            hypothesis_direction="Higher in cycling osteoblast",
        )
    )
    tests.append(
        paired_contrast(
            pseudobulk=l1_scores,
            group_column="celltype.l1",
            group_a="Osteoblast_cycling",
            group_b="Osteoblast",
            score_column=(
                "M40_positive_component_expression_mean"
            ),
            contrast_name=(
                "M40 positive-loading expression: cycling versus "
                "non-cycling osteoblast"
            ),
            hypothesis_direction="Higher in cycling osteoblast",
        )
    )

    result = pd.DataFrame(tests)
    result.loc[
        ~result["estimable"],
        "exact_sign_flip_p",
    ] = np.nan
    result["primary_bh_q"] = bh_adjust(
        result["exact_sign_flip_p"]
    )
    return result
```

#### `build_celltype_summary`

- Lines: 752-811
- Signals: `mean`

```python
def build_celltype_summary(
    dog_celltype: pd.DataFrame,
    score_columns: list[str],
) -> pd.DataFrame:
    rows = []

    for annotation_level in ANNOTATION_LEVELS:
        part = dog_celltype[
            dog_celltype["annotation_level"].eq(annotation_level)
        ]

        for cell_type, cell_part in part.groupby(
            "cell_type",
            observed=True,
        ):
            if cell_part[DOG_COLUMN].nunique() < MIN_DOGS_PER_CELLTYPE:
                continue

            for score_column in score_columns:
                values = cell_part[
                    f"{score_column}_mean"
                ].to_numpy(dtype=float)
                values = values[np.isfinite(values)]

                if values.size == 0:
                    continue

                rows.append(
                    {
                        "annotation_level": annotation_level,
                        "cell_type": cell_type,
                        "score_column": score_column,
                        "n_dogs": int(
                            cell_part[DOG_COLUMN].nunique()
                        ),
                        "total_cells": int(
                            cell_part["n_cells"].sum()
                        ),
                        "median_dog_score": float(
                            np.median(values)
                        ),
                        "mean_dog_score": float(
                            np.mean(values)
                        ),
                        "iqr_low": float(
                            np.quantile(values, 0.25)
                        ),
                        "iqr_high": float(
                            np.quantile(values, 0.75)
                        ),
                        "minimum_dog_score": float(
                            np.min(values)
                        ),
                        "maximum_dog_score": float(
                            np.max(values)
                        ),
                    }
                )

    return pd.DataFrame(rows)
```

#### `rank_stability`

- Lines: 814-917
- Signals: `mean`

```python
def rank_stability(
    dog_celltype: pd.DataFrame,
    score_columns: list[str],
) -> pd.DataFrame:
    rows = []

    for annotation_level in ANNOTATION_LEVELS:
        part = dog_celltype[
            dog_celltype["annotation_level"].eq(annotation_level)
        ].copy()

        eligible_types = (
            part.groupby("cell_type")[DOG_COLUMN]
            .nunique()
            .loc[lambda series: series >= MIN_DOGS_PER_CELLTYPE]
            .index.tolist()
        )
        part = part[part["cell_type"].isin(eligible_types)]

        dogs = sorted(part[DOG_COLUMN].dropna().unique())
        if len(dogs) < MIN_DOGS_PER_CELLTYPE:
            continue

        for score_column in score_columns:
            score_name = f"{score_column}_mean"
            rank_records = []

            full_means = (
                part.groupby("cell_type")[score_name]
                .mean()
                .sort_values(ascending=False)
            )
            full_ranks = pd.Series(
                np.arange(1, len(full_means) + 1),
                index=full_means.index,
            )

            for left_out_dog in dogs:
                loo = part[
                    part[DOG_COLUMN].ne(left_out_dog)
                ]
                loo_means = (
                    loo.groupby("cell_type")[score_name]
                    .mean()
                    .sort_values(ascending=False)
                )
                loo_ranks = pd.Series(
                    np.arange(1, len(loo_means) + 1),
                    index=loo_means.index,
                )

                common = full_ranks.index.intersection(
                    loo_ranks.index
                )
                for cell_type in common:
                    rank_records.append(
                        {
                            "cell_type": cell_type,
                            "left_out_dog": left_out_dog,
                            "full_rank": int(
                                full_ranks.loc[cell_type]
                            ),
                            "loo_rank": int(
                                loo_ranks.loc[cell_type]
                            ),
                        }
                    )

            rank_table = pd.DataFrame(rank_records)
            if rank_table.empty:
                continue

            for cell_type, cell_part in rank_table.groupby(
                "cell_type"
            ):
                rows.append(
                    {
                        "annotation_level": annotation_level,
                        "score_column": score_column,
                        "cell_type": cell_type,
                        "full_rank": int(
                            cell_part["full_rank"].iloc[0]
                        ),
                        "minimum_loo_rank": int(
                            cell_part["loo_rank"].min()
                        ),
                        "maximum_loo_rank": int(
                            cell_part["loo_rank"].max()
                        ),
                        "median_loo_rank": float(
                            cell_part["loo_rank"].median()
                        ),
                        "maximum_absolute_rank_change": int(
                            np.max(
                                np.abs(
                                    cell_part["loo_rank"]
                                    - cell_part["full_rank"]
                                )
                            )
                        ),
                    }
                )

    return pd.DataFrame(rows)
```

#### `gene_celltype_localization`

- Lines: 920-991
- Signals: `mean, weights`

```python
def gene_celltype_localization(
    transformed_expression: np.ndarray,
    genes: list[str],
    metadata: pd.DataFrame,
    weights: pd.DataFrame,
) -> pd.DataFrame:
    gene_to_index = {
        gene: index for index, gene in enumerate(genes)
    }
    rows = []

    metadata_l1 = metadata["celltype.l1"].astype(str)
    eligible_types = (
        metadata_l1.value_counts()
        .loc[lambda series: series >= 100]
        .index.tolist()
    )

    for module in PRIMARY_MODULES:
        module_weights = weights[
            weights["module_label"].eq(module)
        ].drop_duplicates("human_gene_symbol")

        for weight_row in module_weights.itertuples(index=False):
            gene = str(weight_row.human_gene_symbol)
            if gene not in gene_to_index:
                continue

            values = transformed_expression[
                gene_to_index[gene]
            ]

            type_means = {}
            for cell_type in eligible_types:
                mask = metadata_l1.eq(cell_type).to_numpy()
                type_means[cell_type] = float(
                    np.mean(values[mask])
                )

            if not type_means:
                continue

            ordered = sorted(
                type_means.items(),
                key=lambda item: item[1],
                reverse=True,
            )

            for cell_type, mean_value in ordered:
                rows.append(
                    {
                        "module_label": module,
                        "gene_symbol": gene,
                        "risk_oriented_loading": float(
                            weight_row.risk_oriented_loading
                        ),
                        "loading_sign": (
                            "positive"
                            if weight_row.risk_oriented_loading > 0
                            else "negative"
                        ),
                        "cell_type_l1": cell_type,
                        "mean_rank_gaussian_expression": mean_value,
                        "highest_expression_cell_type": ordered[0][0],
                        "highest_expression_value": ordered[0][1],
                        "specificity_range": (
                            ordered[0][1] - ordered[-1][1]
                        ),
                    }
                )

    return pd.DataFrame(rows)
```

#### `create_heatmap`

- Lines: 994-1075
- Signals: `direction`

```python
def create_heatmap(
    celltype_summary: pd.DataFrame,
) -> None:
    selected_scores = [
        "M34_signed_risk_score",
        "M34_negative_component_expression",
        "M40_signed_risk_score",
        "M40_positive_component_expression",
    ]

    part = celltype_summary[
        celltype_summary["annotation_level"].eq("celltype.l1")
        & celltype_summary["score_column"].isin(selected_scores)
        & celltype_summary["n_dogs"].ge(MIN_DOGS_PER_CELLTYPE)
    ].copy()

    if part.empty:
        return

    matrix = part.pivot(
        index="cell_type",
        columns="score_column",
        values="median_dog_score",
    ).reindex(columns=selected_scores)

    matrix = matrix.loc[
        matrix.notna().any(axis=1)
    ].copy()

    order = matrix[
        "M34_signed_risk_score"
    ].sort_values(ascending=False).index
    matrix = matrix.loc[order]

    figure_height = max(5.0, 0.35 * matrix.shape[0] + 1.5)
    fig, ax = plt.subplots(
        figsize=(8.5, figure_height)
    )
    image = ax.imshow(
        matrix.to_numpy(dtype=float),
        aspect="auto",
    )

    labels = [
        "M34 signed risk",
        "M34 negative component",
        "M40 signed risk",
        "M40 positive component",
    ]
    ax.set_xticks(np.arange(len(labels)))
    ax.set_xticklabels(labels, rotation=25, ha="right")
    ax.set_yticks(np.arange(matrix.shape[0]))
    ax.set_yticklabels(matrix.index.tolist())
    ax.set_title(
        "Ammons canine OS atlas: dog-level median cell-type localization"
    )

    for row_index in range(matrix.shape[0]):
        for column_index in range(matrix.shape[1]):
            value = matrix.iloc[row_index, column_index]
            text = "NA" if not np.isfinite(value) else f"{value:.2f}"
            ax.text(
                column_index,
                row_index,
                text,
                ha="center",
                va="center",
                fontsize=7,
            )

    fig.colorbar(image, ax=ax, shrink=0.8)
    fig.tight_layout()
    fig.savefig(
        HEATMAP_PNG,
        dpi=300,
        bbox_inches="tight",
    )
    fig.savefig(
        HEATMAP_PDF,
        bbox_inches="tight",
    )
    plt.close(fig)
```

#### `write_readme`

- Lines: 1128-1186
- Signals: `loadings, direction, coverage, rank`

```python
def write_readme(
    score_coverage: pd.DataFrame,
    targeted: pd.DataFrame,
) -> None:
    text = f"""Ammons canine osteosarcoma single-cell localization
Script version: {SCRIPT_VERSION}

Data
----
- 32,613 cells from the UCSC Cell Browser "All Cells" leaf.
- Six dog/sample identifiers are expected in orig.ident.
- Three annotation levels: celltype.l1, celltype.l2, celltype.l3.

Scoring
-------
Only frozen primary-program genes are extracted from the compressed matrix.
Each gene is rank-Gaussian transformed across all cells. Frozen risk-oriented
loadings are then applied without outcome information.

For each module:
- signed_risk_score uses signed frozen loadings normalized by their absolute sum;
- positive_component_expression uses positive loadings only;
- negative_component_expression uses the absolute magnitude of negative
  loadings and therefore represents expression of the protective/opposite-sign
  component.

Biological replication
----------------------
Primary summaries use dog x cell type or dog x broad compartment pseudobulks.
Individual cells are never treated as independent biological replicates.

Primary targeted contrasts
--------------------------
1. M34 signed risk: osteoblast lineage versus immune compartment.
2. M34 negative-loading expression: immune versus osteoblast lineage.
3. M40 signed risk: cycling versus non-cycling osteoblast.
4. M40 positive-loading expression: cycling versus non-cycling osteoblast.

The exact sign-flip test uses paired dog-level differences, and BH correction
is applied across the four targeted contrasts.

Important guardrails
--------------------
- M34 has only two detected positive-loading genes; that positive component is
  descriptive and is not a primary test.
- M11 and M24 contain no detected negative-loading genes and remain secondary.
- Single-cell localization is same-species biological annotation, not external
  outcome validation.
- A high M34 score in osteoblast-lineage cells and high negative-component
  expression in immune cells would support an immune-depletion/exclusion
  interpretation, not prove causal immune exclusion.

Score coverage:
{score_coverage.to_string(index=False)}

Targeted tests:
{targeted.to_string(index=False)}
"""
    OUTPUT_README.write_text(text, encoding="utf-8")
```

#### `main`

- Lines: 1189-1688
- Signals: `mean, weights, loadings, direction, coverage, rank`

```python
def main() -> None:
    print("=" * 80)
    print("Score frozen programs in Ammons canine OS single-cell atlas")
    print("=" * 80)
    print(f"Script version: {SCRIPT_VERSION}")
    print(f"Project root: {PROJECT_ROOT}")
    print("")
    print("Design:")
    print("  Extract only frozen program genes from the 563.7-MB matrix.")
    print("  Rank-Gaussian transform each gene across cells.")
    print("  Compute signed, positive, and negative component scores.")
    print("  Aggregate by dog x cell type and dog x broad compartment.")
    print("  Use exact paired dog-level tests for four targeted contrasts.")
    print("")

    metadata, weights, preflight_manifest = read_inputs()
    expression_cells, gene_column = expression_header()

    if len(expression_cells) != metadata.shape[0]:
        raise RuntimeError(
            "Expression and metadata cell counts differ."
        )

    metadata_indexed = (
        metadata.set_index(CELL_ID_COLUMN, drop=False)
        .reindex(expression_cells)
    )
    if metadata_indexed[CELL_ID_COLUMN].isna().any():
        raise RuntimeError(
            "Some expression cells are missing from metadata."
        )
    metadata = metadata_indexed.reset_index(drop=True)

    requested_genes = set(
        weights["human_gene_symbol"].dropna()
    )
    genes, raw_expression, duplicate_rows = (
        extract_selected_expression(
            requested_genes=requested_genes,
            expected_cells=len(expression_cells),
        )
    )

    if not genes:
        raise RuntimeError(
            "No frozen-program genes were extracted."
        )

    print(
        f"Extracted {len(genes)} unique frozen-program genes "
        f"across {len(expression_cells)} cells."
    )

    transformed_expression = rank_gaussian_transform(
        raw_expression
    )
    score_weights, score_coverage = build_score_weights(
        weights=weights,
        detected_genes=genes,
    )
    cell_score_values = calculate_cell_scores(
        transformed_expression=transformed_expression,
        score_weights=score_weights,
    )

    cell_scores = metadata[
        [
            CELL_ID_COLUMN,
            DOG_COLUMN,
            *ANNOTATION_LEVELS,
        ]
    ].copy()
    cell_scores["broad_compartment"] = (
        cell_scores["celltype.l1"].map(broad_compartment)
    )
    cell_scores["immune_combined"] = np.where(
        cell_scores["broad_compartment"].isin(
            ["myeloid_innate", "lymphoid"]
        ),
        "immune_combined",
        cell_scores["broad_compartment"],
    )

    for column in cell_score_values.columns:
        cell_scores[column] = cell_score_values[column].to_numpy()

    gene_to_index = {
        gene: index for index, gene in enumerate(genes)
    }
    for module in ALL_MODULES:
        module_genes = (
            weights[
                weights["module_label"].eq(module)
                & weights["human_gene_symbol"].isin(gene_to_index)
            ]["human_gene_symbol"]
            .drop_duplicates()
            .tolist()
        )
        if module_genes:
            indices = [
                gene_to_index[gene] for gene in module_genes
            ]
            detected_fraction = np.mean(
                raw_expression[indices] > 0,
                axis=0,
            )
        else:
            detected_fraction = np.full(
                raw_expression.shape[1],
                np.nan,
            )
        cell_scores[
            f"{module}_detected_gene_fraction"
        ] = detected_fraction.astype(np.float32)

    score_columns = [
        column
        for column in cell_score_values.columns
        if column.endswith(
            (
                "_signed_risk_score",
                "_positive_component_expression",
                "_negative_component_expression",
            )
        )
    ]

    dog_celltype_tables = []
    for annotation_level in ANNOTATION_LEVELS:
        aggregated = aggregate_scores(
            cells=cell_scores,
            group_columns=[DOG_COLUMN, annotation_level],
            score_columns=score_columns,
            minimum_cells=MIN_CELLS_PER_DOG_CELLTYPE,
        )
        aggregated = aggregated.rename(
            columns={annotation_level: "cell_type"}
        )
        aggregated.insert(
            1,
            "annotation_level",
            annotation_level,
        )
        dog_celltype_tables.append(aggregated)

    dog_celltype = pd.concat(
        dog_celltype_tables,
        ignore_index=True,
    )

    l1_scores = dog_celltype[
        dog_celltype["annotation_level"].eq("celltype.l1")
    ].rename(columns={"cell_type": "celltype.l1"})

    compartment_source = cell_scores.copy()
    compartment_source["broad_compartment"] = (
        compartment_source["immune_combined"]
    )
    compartment_scores = aggregate_scores(
        cells=compartment_source,
        group_columns=[DOG_COLUMN, "broad_compartment"],
        score_columns=score_columns,
        minimum_cells=MIN_CELLS_PER_DOG_COMPARTMENT,
    )

    targeted = targeted_tests(
        compartment_scores=compartment_scores,
        l1_scores=l1_scores,
    )
    celltype_summary = build_celltype_summary(
        dog_celltype=dog_celltype,
        score_columns=score_columns,
    )
    stability = rank_stability(
        dog_celltype=dog_celltype,
        score_columns=score_columns,
    )
    gene_localization = gene_celltype_localization(
        transformed_expression=transformed_expression,
        genes=genes,
        metadata=metadata,
        weights=weights,
    )

    np.savez_compressed(
        OUTPUT_SELECTED_EXPRESSION,
        genes=np.asarray(genes, dtype=object),
        cell_ids=np.asarray(expression_cells, dtype=object),
        rank_gaussian_expression=transformed_expression,
    )
    OUTPUT_SELECTED_EXPRESSION_METADATA.write_text(
        json.dumps(
            {
                "script_version": SCRIPT_VERSION,
                "gene_column": gene_column,
                "n_genes": len(genes),
                "n_cells": len(expression_cells),
                "duplicate_expression_rows": duplicate_rows,
                "expression_sha256": sha256_file(EXPRESSION_FILE),
                "metadata_sha256": sha256_file(METADATA_FILE),
            },
            indent=2,
        ),
        encoding="utf-8",
    )

    cell_scores.to_csv(
        OUTPUT_CELL_SCORES,
        index=False,
        compression="gzip",
    )
    dog_celltype.to_csv(OUTPUT_DOG_CELLTYPE, index=False)
    celltype_summary.to_csv(OUTPUT_CELLTYPE_SUMMARY, index=False)
    compartment_scores.to_csv(OUTPUT_COMPARTMENT, index=False)
    targeted.to_csv(OUTPUT_TARGETED, index=False)
    gene_localization.to_csv(
        OUTPUT_GENE_LOCALIZATION,
        index=False,
    )
    stability.to_csv(OUTPUT_RANK_STABILITY, index=False)
    score_coverage.to_csv(OUTPUT_SCORE_COVERAGE, index=False)

    create_heatmap(celltype_summary)

    figure_paths = [HEATMAP_PNG, HEATMAP_PDF]
    figure_paths.extend(
        create_paired_plot(
            pseudobulk=compartment_scores,
            group_column="broad_compartment",
            group_a="osteoblast_lineage",
            group_b="immune_combined",
            score_column="M34_signed_risk_score_mean",
            title=(
                "M34 signed risk: osteoblast lineage versus immune"
            ),
            filename_stem=(
                "Ammons_M34_signed_osteoblast_vs_immune"
            ),
        )
    )
    figure_paths.extend(
        create_paired_plot(
            pseudobulk=compartment_scores,
            group_column="broad_compartment",
            group_a="immune_combined",
            group_b="osteoblast_lineage",
            score_column=(
                "M34_negative_component_expression_mean"
            ),
            title=(
                "M34 negative-loading expression: immune versus "
                "osteoblast lineage"
            ),
            filename_stem=(
                "Ammons_M34_negative_immune_vs_osteoblast"
            ),
        )
    )
    figure_paths.extend(
        create_paired_plot(
            pseudobulk=l1_scores,
            group_column="celltype.l1",
            group_a="Osteoblast_cycling",
            group_b="Osteoblast",
            score_column="M40_signed_risk_score_mean",
            title=(
                "M40 signed risk: cycling versus non-cycling osteoblast"
            ),
            filename_stem=(
                "Ammons_M40_signed_cycling_vs_osteoblast"
            ),
        )
    )
    figure_paths.extend(
        create_paired_plot(
            pseudobulk=l1_scores,
            group_column="celltype.l1",
            group_a="Osteoblast_cycling",
            group_b="Osteoblast",
            score_column=(
                "M40_positive_component_expression_mean"
            ),
            title=(
                "M40 positive component: cycling versus "
                "non-cycling osteoblast"
            ),
            filename_stem=(
                "Ammons_M40_positive_cycling_vs_osteoblast"
            ),
        )
    )

    write_readme(
        score_coverage=score_coverage,
        targeted=targeted,
    )

    input_paths = [
        METADATA_FILE,
        EXPRESSION_FILE,
        PREFLIGHT_MANIFEST_FILE,
        STRICT_WEIGHTS_FILE,
    ]
    output_paths = [
        OUTPUT_SELECTED_EXPRESSION,
        OUTPUT_SELECTED_EXPRESSION_METADATA,
        OUTPUT_CELL_SCORES,
        OUTPUT_DOG_CELLTYPE,
        OUTPUT_CELLTYPE_SUMMARY,
        OUTPUT_COMPARTMENT,
        OUTPUT_TARGETED,
        OUTPUT_GENE_LOCALIZATION,
        OUTPUT_RANK_STABILITY,
        OUTPUT_SCORE_COVERAGE,
        OUTPUT_README,
        *[path for path in figure_paths if path.exists()],
    ]

    manifest = {
        "script_version": SCRIPT_VERSION,
        "created_utc": utc_now_iso(),
        "outcome_loaded": False,
        "biological_replicate": DOG_COLUMN,
        "cell_id_column": CELL_ID_COLUMN,
        "annotation_levels": ANNOTATION_LEVELS,
        "scoring": {
            "gene_transform": (
                "rank-Gaussian across all cells"
            ),
            "signed_score": (
                "sum(risk-oriented loading * transformed expression) "
                "/ sum(abs(loading))"
            ),
            "positive_component": (
                "positive-loading weighted mean"
            ),
            "negative_component": (
                "absolute negative-loading weighted mean expression"
            ),
        },
        "targeted_tests": targeted[
            [
                "contrast_name",
                "score_column",
                "n_paired_dogs",
                "mean_paired_difference",
                "exact_sign_flip_p",
                "primary_bh_q",
            ]
        ].to_dict(orient="records"),
        "guardrails": [
            "No clinical outcome or endpoint was loaded.",
            "Frozen program loadings and signs were not changed.",
            "Individual cells were not used as biological replicates.",
            "Primary inference used paired dog-level pseudobulk contrasts.",
            "M34 positive component is descriptive because only two positive-loading genes were detected.",
            "Single-cell localization does not establish causality or external outcome validation.",
        ],
        "inputs": {
            str(path): {
                "sha256": sha256_file(path),
                "size_bytes": path.stat().st_size,
            }
            for path in input_paths
        },
        "outputs": {
            str(path): {
                "sha256": sha256_file(path),
                "size_bytes": path.stat().st_size,
            }
            for path in output_paths
            if path.exists()
        },
    }
    OUTPUT_MANIFEST.write_text(
        json.dumps(manifest, indent=2),
        encoding="utf-8",
    )

    print("")
    print("=" * 80)
    print("Score gene coverage")
    print("=" * 80)
    print(score_coverage.to_string(index=False))

    print("")
    print("=" * 80)
    print("Primary targeted single-cell localization tests")
    print("=" * 80)
    print(
        targeted[
            [
                "contrast_name",
                "n_paired_dogs",
                "mean_paired_difference",
                "median_paired_difference",
                "bootstrap_ci_low",
                "bootstrap_ci_high",
                "exact_sign_flip_p",
                "primary_bh_q",
                "estimable",
            ]
        ].to_string(index=False)
    )

    print("")
    print("=" * 80)
    print("Top cell types by dog-level median score")
    print("=" * 80)
    selected_summary = celltype_summary[
        celltype_summary["annotation_level"].eq("celltype.l1")
        & celltype_summary["score_column"].isin(
            [
                "M34_signed_risk_score",
                "M34_negative_component_expression",
                "M40_signed_risk_score",
                "M40_positive_component_expression",
            ]
        )
    ]
    top = (
        selected_summary.sort_values(
            ["score_column", "median_dog_score"],
            ascending=[True, False],
        )
        .groupby("score_column", as_index=False)
        .head(10)
    )
    print(
        top[
            [
                "score_column",
                "cell_type",
                "n_dogs",
                "total_cells",
                "median_dog_score",
                "iqr_low",
                "iqr_high",
            ]
        ].to_string(index=False)
    )

    print("")
    print("=" * 80)
    print("Leave-one-dog-out rank stability")
    print("=" * 80)
    stable_top = (
        stability[
            stability["annotation_level"].eq("celltype.l1")
        ]
        .sort_values(
            ["score_column", "full_rank"],
        )
        .groupby("score_column", as_index=False)
        .head(10)
    )
    print(
        stable_top[
            [
                "score_column",
                "cell_type",
                "full_rank",
                "minimum_loo_rank",
                "maximum_loo_rank",
                "maximum_absolute_rank_change",
            ]
        ].to_string(index=False)
    )

    print("")
    print("=" * 80)
    print("Interpretation guardrails")
    print("=" * 80)
    print("No clinical outcome was loaded.")
    print("Dog, not cell, is the unit of primary inference.")
    print("M34 positive component is descriptive because it contains only two detected genes.")
    print("The M34 negative component represents expression of opposite-sign/protective genes.")
    print("Single-cell localization can support immune-depletion or proliferation localization, but cannot prove causality.")

    print("")
    print("Saved:")
    for path in [
        OUTPUT_SELECTED_EXPRESSION,
        OUTPUT_SELECTED_EXPRESSION_METADATA,
        OUTPUT_CELL_SCORES,
        OUTPUT_DOG_CELLTYPE,
        OUTPUT_CELLTYPE_SUMMARY,
        OUTPUT_COMPARTMENT,
        OUTPUT_TARGETED,
        OUTPUT_GENE_LOCALIZATION,
        OUTPUT_RANK_STABILITY,
        OUTPUT_SCORE_COVERAGE,
        OUTPUT_README,
        OUTPUT_MANIFEST,
        HEATMAP_PNG,
        HEATMAP_PDF,
    ]:
        print(path)
    print(f"Figures directory: {FIGURES_DIR}")
    print("Done.")
```

### Relevant standalone lines

| Line | Signals | Code |
|---:|---|---|
| 186 | weights | `weights = pd.read_csv(STRICT_WEIGHTS_FILE)` |
| 212 | weights | `weights["human_gene_symbol"] = (` |
| 213 | weights | `weights["human_gene_symbol"].astype(str).str.upper()` |
| 215 | weights | `weights["risk_oriented_loading"] = pd.to_numeric(` |
| 216 | weights | `weights["risk_oriented_loading"],` |
| 219 | weights | `weights = weights[` |
| 220 | weights | `weights["module_label"].isin(ALL_MODULES)` |
| 221 | weights | `& weights["risk_oriented_loading"].notna()` |
| 224 | weights | `return metadata, weights, manifest` |
| 292 | mean | `matrix[index] = np.mean(` |
| 311 | std | `if finite.sum() < 2 or np.nanstd(values[finite]) == 0:` |
| 333 | weights | `weights: pd.DataFrame,` |
| 343 | weights | `part = weights[` |
| 344 | weights | `weights["module_label"].eq(module)` |
| 360 | loadings | `loadings = available[` |
| 364 | loadings | `positive_mask = loadings > 0` |
| 365 | loadings | `negative_mask = loadings < 0` |
| 369 | loadings | `"loadings": loadings,` |
| 371 | loadings | `"positive_weights": loadings[positive_mask],` |
| 373 | loadings | `"negative_weights": np.abs(loadings[negative_mask]),` |
| 397 | loadings | `"only two positive-loading genes are available."` |
| 410 | weights | `weights: np.ndarray,` |
| 419 | weights | `denominator = float(np.sum(np.abs(weights)))` |
| 428 | weights | `weights.astype(np.float32)` |
| 429 | matrix_product | `@ expression[indices]` |
| 442 | loadings | `loadings = specification["loadings"]` |
| 444 | direction | `signed = weighted_component(` |
| 447 | loadings | `loadings,` |
| 460 | direction | `score_columns[f"{module}_signed_risk_score"] = signed` |
| 517 | mean | `means = grouped[score_columns].mean()` |
| 521 | mean | `].mean()` |
| 549 | mean | `observed = float(np.mean(values))` |
| 557 | mean | `np.mean(` |
| 564 | mean | `np.mean(` |
| 575 | mean | `statistics_array = np.mean(` |
| 608 | mean | `means = values[indices].mean(axis=1)` |
| 687 | direction | `"M34 signed risk: osteoblast lineage versus immune"` |
| 704 | loadings | `"M34 negative-loading expression: immune versus "` |
| 718 | direction | `"M40 signed risk: cycling versus non-cycling "` |
| 734 | loadings | `"M40 positive-loading expression: cycling versus "` |
| 794 | mean | `np.mean(values)` |
| 843 | mean | `.mean()` |
| 857 | mean | `.mean()` |
| 924 | weights | `weights: pd.DataFrame,` |
| 939 | weights | `module_weights = weights[` |
| 940 | weights | `weights["module_label"].eq(module)` |
| 956 | mean | `np.mean(values[mask])` |
| 1038 | direction | `"M34 signed risk",` |
| 1040 | direction | `"M40 signed risk",` |
| 1144 | rank | `Each gene is rank-Gaussian transformed across all cells. Frozen risk-oriented` |
| 1145 | loadings | `loadings are then applied without outcome information.` |
| 1148 | loadings, direction | `- signed_risk_score uses signed frozen loadings normalized by their absolute sum;` |
| 1149 | loadings | `- positive_component_expression uses positive loadings only;` |
| 1151 | loadings, direction | `loadings and therefore represents expression of the protective/opposite-sign` |
| 1161 | direction | `1. M34 signed risk: osteoblast lineage versus immune compartment.` |
| 1162 | loadings | `2. M34 negative-loading expression: immune versus osteoblast lineage.` |
| 1163 | direction | `3. M40 signed risk: cycling versus non-cycling osteoblast.` |
| 1164 | loadings | `4. M40 positive-loading expression: cycling versus non-cycling osteoblast.` |
| 1166 | direction | `The exact sign-flip test uses paired dog-level differences, and BH correction` |
| 1171 | loadings | `- M34 has only two detected positive-loading genes; that positive component is` |
| 1173 | loadings | `- M11 and M24 contain no detected negative-loading genes and remain secondary.` |
| 1180 | coverage | `Score coverage:` |
| 1198 | rank | `print("  Rank-Gaussian transform each gene across cells.")` |
| 1199 | direction | `print("  Compute signed, positive, and negative component scores.")` |
| 1204 | weights | `metadata, weights, preflight_manifest = read_inputs()` |
| 1223 | weights | `weights["human_gene_symbol"].dropna()` |
| 1246 | weights | `weights=weights,` |
| 1280 | weights | `weights[` |
| 1281 | weights | `weights["module_label"].eq(module)` |
| 1282 | weights | `& weights["human_gene_symbol"].isin(gene_to_index)` |
| 1291 | mean | `detected_fraction = np.mean(` |
| 1370 | weights | `weights=weights,` |
| 1422 | direction | `"M34 signed risk: osteoblast lineage versus immune"` |
| 1439 | loadings | `"M34 negative-loading expression: immune versus "` |
| 1455 | direction | `"M40 signed risk: cycling versus non-cycling osteoblast"` |
| 1516 | rank | `"rank-Gaussian across all cells"` |
| 1519 | loadings | `"sum(risk-oriented loading * transformed expression) "` |
| 1520 | loadings | `"/ sum(abs(loading))"` |
| 1523 | weights, loadings | `"positive-loading weighted mean"` |
| 1526 | weights, loadings | `"absolute negative-loading weighted mean expression"` |
| 1541 | loadings | `"Frozen program loadings and signs were not changed.",` |
| 1544 | loadings | `"M34 positive component is descriptive because only two positive-loading genes were detected.",` |
| 1570 | coverage | `print("Score gene coverage")` |
| 1633 | rank | `print("Leave-one-dog-out rank stability")` |
| 1665 | direction | `print("The M34 negative component represents expression of opposite-sign/protective genes.")` |

### Referenced result-table schemas

#### `results/tables/Ammons_scRNA_celltype_rank_stability.csv`

- Rows: 972
- SHA-256: `07f8337eeaa33e59bef505b9e931fc6645464313f084e03d35de5ce507d7c70e`
- Columns: `annotation_level`, `score_column`, `cell_type`, `full_rank`, `minimum_loo_rank`, `maximum_loo_rank`, `median_loo_rank`, `maximum_absolute_rank_change`

#### `results/tables/Ammons_scRNA_celltype_score_summary.csv`

- Rows: 810
- SHA-256: `bb9722efe5fd8461f5d4fb5837659d4fe9d8a79a470243c73a096fc2c6fddf1f`
- Columns: `annotation_level`, `cell_type`, `score_column`, `n_dogs`, `total_cells`, `median_dog_score`, `mean_dog_score`, `iqr_low`, `iqr_high`, `minimum_dog_score`, `maximum_dog_score`

#### `results/tables/Ammons_scRNA_dog_celltype_pseudobulk_scores.csv`

- Rows: 548
- SHA-256: `53f917ddd151c65d21da54b0668e68cc31b7bdb03b25fe13217384f31dbe6195`
- Columns: `orig.ident`, `annotation_level`, `cell_type`, `n_cells`, `M34_signed_risk_score_mean`, `M34_positive_component_expression_mean`, `M34_negative_component_expression_mean`, `M40_signed_risk_score_mean`, `M40_positive_component_expression_mean`, `M40_negative_component_expression_mean`, `M11_signed_risk_score_mean`, `M11_positive_component_expression_mean`, `M11_negative_component_expression_mean`, `M24_signed_risk_score_mean`, `M24_positive_component_expression_mean`, `M24_negative_component_expression_mean`, `M34_signed_risk_score_median`, `M34_positive_component_expression_median`, `M34_negative_component_expression_median`, `M40_signed_risk_score_median`, `M40_positive_component_expression_median`, `M40_negative_component_expression_median`, `M11_signed_risk_score_median`, `M11_positive_component_expression_median`, `M11_negative_component_expression_median`, `M24_signed_risk_score_median`, `M24_positive_component_expression_median`, `M24_negative_component_expression_median`, `M34_detected_gene_fraction_mean`, `M40_detected_gene_fraction_mean`, `M11_detected_gene_fraction_mean`, `M24_detected_gene_fraction_mean`

#### `results/tables/Ammons_scRNA_dog_compartment_scores.csv`

- Rows: 37
- SHA-256: `f9a7cc1ab8807db12f791543a4f0ec9b39da7066d670df01789676c9506a0041`
- Columns: `orig.ident`, `broad_compartment`, `n_cells`, `M34_signed_risk_score_mean`, `M34_positive_component_expression_mean`, `M34_negative_component_expression_mean`, `M40_signed_risk_score_mean`, `M40_positive_component_expression_mean`, `M40_negative_component_expression_mean`, `M11_signed_risk_score_mean`, `M11_positive_component_expression_mean`, `M11_negative_component_expression_mean`, `M24_signed_risk_score_mean`, `M24_positive_component_expression_mean`, `M24_negative_component_expression_mean`, `M34_signed_risk_score_median`, `M34_positive_component_expression_median`, `M34_negative_component_expression_median`, `M40_signed_risk_score_median`, `M40_positive_component_expression_median`, `M40_negative_component_expression_median`, `M11_signed_risk_score_median`, `M11_positive_component_expression_median`, `M11_negative_component_expression_median`, `M24_signed_risk_score_median`, `M24_positive_component_expression_median`, `M24_negative_component_expression_median`, `M34_detected_gene_fraction_mean`, `M40_detected_gene_fraction_mean`, `M11_detected_gene_fraction_mean`, `M24_detected_gene_fraction_mean`

#### `results/tables/Ammons_scRNA_gene_celltype_localization.csv`

- Rows: 4978
- SHA-256: `44ff18ea215ed9175077c86d361baa0126471f2cd1c8d4d7393ac67c74859aa1`
- Columns: `module_label`, `gene_symbol`, `risk_oriented_loading`, `loading_sign`, `cell_type_l1`, `mean_rank_gaussian_expression`, `highest_expression_cell_type`, `highest_expression_value`, `specificity_range`

#### `results/tables/Ammons_scRNA_score_gene_coverage.csv`

- Rows: 4
- SHA-256: `13ae84b614bc12dd08ef47485a5cf70ff56c74df6c43ca91d7806ee3a283cf49`
- Columns: `module_label`, `n_frozen_genes`, `n_detected_genes`, `coverage_fraction`, `n_positive_detected`, `n_negative_detected`, `signed_score_estimable`, `positive_component_estimable`, `negative_component_estimable`, `component_guardrail`

#### `results/tables/Ammons_scRNA_targeted_localization_tests.csv`

- Rows: 4
- SHA-256: `5f75eac6ccc4f12b710c02f4220e78ee13e5692a7c694097b4fbb0e9049f2203`
- Columns: `contrast_name`, `score_column`, `group_column`, `group_a`, `group_b`, `difference_definition`, `hypothesis_direction`, `n_paired_dogs`, `mean_paired_difference`, `median_paired_difference`, `bootstrap_ci_low`, `bootstrap_ci_high`, `exact_sign_flip_p`, `dog_level_differences`, `estimable`, `primary_bh_q`

#### `results/tables/GSE238110_frozen_transfer_gene_weights_strict.csv`

- Rows: 321
- SHA-256: `4f065aa4c4edf117a0c74015840d2b4b2347929f172cd517e1818ba0f6163b91`
- Columns: `module_label`, `canine_gene`, `canine_gene_symbol`, `human_gene_symbol`, `ortholog_qc_status`, `raw_pca_loading`, `risk_oriented_loading`, `absolute_risk_loading`, `is_strict_mapping`, `is_broad_mapping`, `human_symbol_duplicate_count`, `normalized_abs_sum_weight`

---

## `scripts/44_recompute_ammons_six_dog_localization.py`

- SHA-256: `dbbf957c450f64138171336048dbf3c727be3dc418f3e97c3664f1aa62f4d1a0`

### Relevant functions

#### `exact_sign_flip_test`

- Lines: 179-207
- Signals: `mean`

```python
def exact_sign_flip_test(
    differences: np.ndarray,
) -> tuple[float, float]:
    values = np.asarray(differences, dtype=float)
    values = values[np.isfinite(values)]

    if values.size == 0:
        return np.nan, np.nan

    observed = float(np.mean(values))
    statistics = []

    for signs in product([-1.0, 1.0], repeat=values.size):
        statistics.append(
            float(
                np.mean(
                    values * np.asarray(signs, dtype=float)
                )
            )
        )

    permutation_statistics = np.asarray(statistics, dtype=float)
    p = float(
        np.mean(
            np.abs(permutation_statistics)
            >= abs(observed) - 1e-12
        )
    )
    return observed, p
```

#### `bootstrap_mean_ci`

- Lines: 210-230
- Signals: `mean`

```python
def bootstrap_mean_ci(
    differences: np.ndarray,
) -> tuple[float, float]:
    values = np.asarray(differences, dtype=float)
    values = values[np.isfinite(values)]

    if values.size == 0:
        return np.nan, np.nan

    rng = np.random.default_rng(RANDOM_SEED)
    indices = rng.integers(
        0,
        values.size,
        size=(N_BOOTSTRAP, values.size),
    )
    means = values[indices].mean(axis=1)

    return (
        float(np.quantile(means, 0.025)),
        float(np.quantile(means, 0.975)),
    )
```

#### `estimable_score_columns`

- Lines: 233-252
- Signals: `coverage`

```python
def estimable_score_columns(
    coverage: pd.DataFrame,
) -> list[str]:
    columns: list[str] = []

    for row in coverage.itertuples(index=False):
        module = str(row.module_label)

        if bool(row.signed_score_estimable):
            columns.append(f"{module}_signed_risk_score")
        if bool(row.positive_component_estimable):
            columns.append(
                f"{module}_positive_component_expression"
            )
        if bool(row.negative_component_estimable):
            columns.append(
                f"{module}_negative_component_expression"
            )

    return sorted(set(columns))
```

#### `validate_six_dog_mapping`

- Lines: 255-308
- Signals: ``

```python
def validate_six_dog_mapping(
    metadata: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    required = [CELL_ID_COLUMN, LIBRARY_COLUMN, DOG_COLUMN]
    missing = [
        column for column in required if column not in metadata.columns
    ]
    if missing:
        raise ValueError(
            f"Required metadata columns are missing: {missing}"
        )

    mapping = (
        metadata[[LIBRARY_COLUMN, DOG_COLUMN]]
        .dropna()
        .astype(str)
        .drop_duplicates()
        .sort_values([DOG_COLUMN, LIBRARY_COLUMN])
        .reset_index(drop=True)
    )

    if mapping[LIBRARY_COLUMN].nunique() != 8:
        raise RuntimeError(
            "Expected eight Cell Ranger library identifiers."
        )
    if mapping[DOG_COLUMN].nunique() != 6:
        raise RuntimeError(
            "Expected the metadata name column to define six dogs."
        )

    libraries_per_dog = (
        mapping.groupby(DOG_COLUMN)[LIBRARY_COLUMN]
        .nunique()
        .sort_values(ascending=False)
    )

    if sorted(libraries_per_dog.tolist()) != [1, 1, 1, 1, 2, 2]:
        raise RuntimeError(
            "Expected two dogs with two libraries and four dogs "
            "with one library."
        )

    counts = (
        metadata[DOG_COLUMN]
        .astype(str)
        .value_counts()
        .rename_axis("dog_id")
        .reset_index(name="n_cells")
        .sort_values("dog_id")
    )
    counts["fraction_cells"] = counts["n_cells"] / metadata.shape[0]
    counts["n_libraries"] = counts["dog_id"].map(libraries_per_dog)

    return mapping.rename(columns={DOG_COLUMN: "dog_id"}), counts
```

#### `aggregate_scores`

- Lines: 369-399
- Signals: `mean`

```python
def aggregate_scores(
    cells: pd.DataFrame,
    group_columns: list[str],
    score_columns: list[str],
    minimum_cells: int,
) -> pd.DataFrame:
    grouped = cells.groupby(
        group_columns,
        dropna=False,
        observed=True,
    )

    counts = grouped.size().rename("n_cells")
    means = grouped[score_columns].mean()
    medians = grouped[score_columns].median()

    means.columns = [
        f"{column}_mean" for column in means.columns
    ]
    medians.columns = [
        f"{column}_median" for column in medians.columns
    ]

    result = pd.concat(
        [counts, means, medians],
        axis=1,
    ).reset_index()

    return result[
        result["n_cells"].ge(minimum_cells)
    ].reset_index(drop=True)
```

#### `targeted_tests`

- Lines: 483-549
- Signals: `loadings, direction`

```python
def targeted_tests(
    compartment_scores: pd.DataFrame,
    l1_scores: pd.DataFrame,
) -> pd.DataFrame:
    tests = [
        paired_contrast(
            pseudobulk=compartment_scores,
            group_column="analysis_compartment",
            group_a="osteoblast_lineage",
            group_b="immune_combined",
            score_column="M34_signed_risk_score_mean",
            contrast_name=(
                "M34 signed risk: osteoblast lineage versus immune"
            ),
            hypothesis_direction="Higher in osteoblast lineage",
        ),
        paired_contrast(
            pseudobulk=compartment_scores,
            group_column="analysis_compartment",
            group_a="immune_combined",
            group_b="osteoblast_lineage",
            score_column=(
                "M34_negative_component_expression_mean"
            ),
            contrast_name=(
                "M34 negative-loading expression: immune versus "
                "osteoblast lineage"
            ),
            hypothesis_direction="Higher in immune compartment",
        ),
        paired_contrast(
            pseudobulk=l1_scores,
            group_column="celltype.l1",
            group_a="Osteoblast_cycling",
            group_b="Osteoblast",
            score_column="M40_signed_risk_score_mean",
            contrast_name=(
                "M40 signed risk: cycling versus non-cycling "
                "osteoblast"
            ),
            hypothesis_direction="Higher in cycling osteoblast",
        ),
        paired_contrast(
            pseudobulk=l1_scores,
            group_column="celltype.l1",
            group_a="Osteoblast_cycling",
            group_b="Osteoblast",
            score_column=(
                "M40_positive_component_expression_mean"
            ),
            contrast_name=(
                "M40 positive-loading expression: cycling versus "
                "non-cycling osteoblast"
            ),
            hypothesis_direction="Higher in cycling osteoblast",
        ),
    ]

    result = pd.DataFrame(tests)
    result.loc[
        ~result["estimable"],
        "exact_sign_flip_p",
    ] = np.nan
    result["primary_bh_q"] = bh_adjust(
        result["exact_sign_flip_p"]
    )
    return result
```

#### `celltype_summary`

- Lines: 552-611
- Signals: `mean`

```python
def celltype_summary(
    dog_celltype: pd.DataFrame,
    score_columns: list[str],
) -> pd.DataFrame:
    rows = []

    for annotation_level in ANNOTATION_LEVELS:
        part = dog_celltype[
            dog_celltype["annotation_level"].eq(annotation_level)
        ]

        for cell_type, cell_part in part.groupby(
            "cell_type",
            observed=True,
        ):
            n_dogs = int(cell_part["dog_id"].nunique())
            if n_dogs < MIN_DOGS_PER_CELLTYPE:
                continue

            for score_column in score_columns:
                values = pd.to_numeric(
                    cell_part[f"{score_column}_mean"],
                    errors="coerce",
                ).to_numpy(dtype=float)
                values = values[np.isfinite(values)]

                if values.size == 0:
                    continue

                rows.append(
                    {
                        "annotation_level": annotation_level,
                        "cell_type": cell_type,
                        "score_column": score_column,
                        "n_dogs": n_dogs,
                        "total_cells": int(
                            cell_part["n_cells"].sum()
                        ),
                        "median_dog_score": float(
                            np.median(values)
                        ),
                        "mean_dog_score": float(
                            np.mean(values)
                        ),
                        "iqr_low": float(
                            np.quantile(values, 0.25)
                        ),
                        "iqr_high": float(
                            np.quantile(values, 0.75)
                        ),
                        "minimum_dog_score": float(
                            np.min(values)
                        ),
                        "maximum_dog_score": float(
                            np.max(values)
                        ),
                    }
                )

    return pd.DataFrame(rows)
```

#### `rank_stability`

- Lines: 614-716
- Signals: `mean`

```python
def rank_stability(
    dog_celltype: pd.DataFrame,
    score_columns: list[str],
) -> pd.DataFrame:
    rows = []

    for annotation_level in ANNOTATION_LEVELS:
        part = dog_celltype[
            dog_celltype["annotation_level"].eq(annotation_level)
        ].copy()

        eligible_types = (
            part.groupby("cell_type")["dog_id"]
            .nunique()
            .loc[lambda values: values >= MIN_DOGS_PER_CELLTYPE]
            .index.tolist()
        )
        part = part[part["cell_type"].isin(eligible_types)]

        dogs = sorted(part["dog_id"].dropna().unique())
        if len(dogs) < MIN_DOGS_PER_CELLTYPE:
            continue

        for score_column in score_columns:
            score_name = f"{score_column}_mean"

            full_means = (
                part.groupby("cell_type")[score_name]
                .mean()
                .sort_values(ascending=False)
            )
            full_ranks = pd.Series(
                np.arange(1, len(full_means) + 1),
                index=full_means.index,
            )

            records = []
            for left_out_dog in dogs:
                loo = part[part["dog_id"].ne(left_out_dog)]
                loo_means = (
                    loo.groupby("cell_type")[score_name]
                    .mean()
                    .sort_values(ascending=False)
                )
                loo_ranks = pd.Series(
                    np.arange(1, len(loo_means) + 1),
                    index=loo_means.index,
                )

                common = full_ranks.index.intersection(
                    loo_ranks.index
                )

                for cell_type in common:
                    records.append(
                        {
                            "cell_type": cell_type,
                            "left_out_dog": left_out_dog,
                            "full_rank": int(
                                full_ranks.loc[cell_type]
                            ),
                            "loo_rank": int(
                                loo_ranks.loc[cell_type]
                            ),
                        }
                    )

            record_table = pd.DataFrame(records)
            if record_table.empty:
                continue

            for cell_type, cell_part in record_table.groupby(
                "cell_type"
            ):
                rows.append(
                    {
                        "annotation_level": annotation_level,
                        "score_column": score_column,
                        "cell_type": cell_type,
                        "full_rank": int(
                            cell_part["full_rank"].iloc[0]
                        ),
                        "minimum_loo_rank": int(
                            cell_part["loo_rank"].min()
                        ),
                        "maximum_loo_rank": int(
                            cell_part["loo_rank"].max()
                        ),
                        "median_loo_rank": float(
                            cell_part["loo_rank"].median()
                        ),
                        "maximum_absolute_rank_change": int(
                            np.max(
                                np.abs(
                                    cell_part["loo_rank"]
                                    - cell_part["full_rank"]
                                )
                            )
                        ),
                    }
                )

    return pd.DataFrame(rows)
```

#### `build_locked_single_cell`

- Lines: 799-946
- Signals: `loadings, direction`

```python
def build_locked_single_cell(
    targeted: pd.DataFrame,
) -> pd.DataFrame:
    indexed = targeted.set_index("contrast_name")

    m34_signed_name = (
        "M34 signed risk: osteoblast lineage versus immune"
    )
    m34_negative_name = (
        "M34 negative-loading expression: immune versus "
        "osteoblast lineage"
    )
    m40_signed_name = (
        "M40 signed risk: cycling versus non-cycling osteoblast"
    )
    m40_positive_name = (
        "M40 positive-loading expression: cycling versus "
        "non-cycling osteoblast"
    )

    rows = [
        {
            "module_label": "M34",
            "single_cell_localization_class": (
                "immune_negative_component_with_"
                "osteoblast_high_signed_risk"
            ),
            "n_biological_dogs": 6,
            "primary_contrast_1": m34_signed_name,
            "primary_effect_1": indexed.loc[
                m34_signed_name,
                "mean_paired_difference",
            ],
            "primary_q_1": indexed.loc[
                m34_signed_name,
                "primary_bh_q",
            ],
            "primary_contrast_2": m34_negative_name,
            "primary_effect_2": indexed.loc[
                m34_negative_name,
                "mean_paired_difference",
            ],
            "primary_q_2": indexed.loc[
                m34_negative_name,
                "primary_bh_q",
            ],
            "locked_single_cell_interpretation": (
                "Across six biological dogs, the M34 signed "
                "risk-oriented score was consistently higher in "
                "osteoblast-lineage than immune compartments, while "
                "the negative-loading component was consistently "
                "higher in immune than osteoblast-lineage cells. "
                "Because 153 of 155 detected M34 genes have negative "
                "loadings, M34 is best interpreted as an inverse "
                "immune-lineage or immune-depletion/exclusion axis."
            ),
            "replicate_guardrail": (
                "Dogs 1 and 2 each had two technical replicate "
                "libraries, which were combined by metadata name "
                "before dog-level pseudobulk inference."
            ),
        },
        {
            "module_label": "M40",
            "single_cell_localization_class": (
                "pan_cycling_program_with_"
                "cycling_osteoblast_enrichment"
            ),
            "n_biological_dogs": 6,
            "primary_contrast_1": m40_signed_name,
            "primary_effect_1": indexed.loc[
                m40_signed_name,
                "mean_paired_difference",
            ],
            "primary_q_1": indexed.loc[
                m40_signed_name,
                "primary_bh_q",
            ],
            "primary_contrast_2": m40_positive_name,
            "primary_effect_2": indexed.loc[
                m40_positive_name,
                "mean_paired_difference",
            ],
            "primary_q_2": indexed.loc[
                m40_positive_name,
                "primary_bh_q",
            ],
            "locked_single_cell_interpretation": (
                "Across six biological dogs, M40 signed and "
                "positive-loading scores were consistently higher in "
                "cycling than non-cycling osteoblasts. Together with "
                "high scores in cycling T-cell and osteoclast "
                "populations, M40 is a broad cycling/proliferation "
                "axis with clear tumor-lineage enrichment rather than "
                "an osteoblast-specific program."
            ),
            "replicate_guardrail": (
                "Dogs 1 and 2 each had two technical replicate "
                "libraries, which were combined by metadata name "
                "before dog-level pseudobulk inference."
            ),
        },
        {
            "module_label": "M11",
            "single_cell_localization_class": (
                "secondary_positive_component_only"
            ),
            "n_biological_dogs": 6,
            "primary_contrast_1": "",
            "primary_effect_1": np.nan,
            "primary_q_1": np.nan,
            "primary_contrast_2": "",
            "primary_effect_2": np.nan,
            "primary_q_2": np.nan,
            "locked_single_cell_interpretation": (
                "M11 remains a secondary positive-component-only "
                "localization and does not alter its locked "
                "cross-species evidence grade."
            ),
            "replicate_guardrail": (
                "No primary single-cell localization test was "
                "specified for M11."
            ),
        },
        {
            "module_label": "M24",
            "single_cell_localization_class": (
                "secondary_positive_component_only"
            ),
            "n_biological_dogs": 6,
            "primary_contrast_1": "",
            "primary_effect_1": np.nan,
            "primary_q_1": np.nan,
            "primary_contrast_2": "",
            "primary_effect_2": np.nan,
            "primary_q_2": np.nan,
            "locked_single_cell_interpretation": (
                "M24 remains a secondary positive-component-only "
                "localization and cannot compensate for limited "
                "cross-species representation and outcome evidence."
            ),
            "replicate_guardrail": (
                "No primary single-cell localization test was "
                "specified for M24."
            ),
        },
    ]
    return pd.DataFrame(rows)
```

#### `write_sentences`

- Lines: 1001-1067
- Signals: `direction`

```python
def write_sentences(
    locked: pd.DataFrame,
    targeted: pd.DataFrame,
) -> None:
    indexed = locked.set_index("module_label")

    all_floor = bool(
        targeted["all_same_nonzero_sign"].all()
        and np.allclose(
            targeted["exact_sign_flip_p"].to_numpy(dtype=float),
            0.03125,
        )
    )

    lines = [
        "Locked six-dog single-cell localization results",
        "================================================",
        "",
        "Replicate correction",
        "--------------------",
        (
            "The Cell Browser dataset contained eight Cell Ranger "
            "libraries but six biological dogs. Dogs 1 and 2 each "
            "contributed two technical replicate libraries, which "
            "were combined using the metadata name column before "
            "biological inference."
        ),
        "",
        "M34",
        "---",
        indexed.loc[
            "M34",
            "locked_single_cell_interpretation",
        ],
        "",
        "M40",
        "---",
        indexed.loc[
            "M40",
            "locked_single_cell_interpretation",
        ],
        "",
        "Statistical interpretation",
        "--------------------------",
        (
            "All four prespecified contrasts retained the same "
            "direction in all six biological dogs. The exact "
            "two-sided sign-flip P value was 0.03125, the minimum "
            "attainable value for six paired dogs."
            if all_floor
            else
            "Dog-level exact paired results are reported in "
            "Ammons_scRNA_targeted_localization_tests_six_dogs.csv."
        ),
        "",
        "Supersession rule",
        "-----------------",
        (
            "The six-dog analysis supersedes the earlier "
            "eight-orig.ident P=0.0078125 analysis."
        ),
    ]

    OUTPUT_SENTENCES.write_text(
        "\n".join(lines),
        encoding="utf-8",
    )
```

#### `write_readme`

- Lines: 1070-1111
- Signals: `direction`

```python
def write_readme(
    mapping: pd.DataFrame,
    targeted: pd.DataFrame,
) -> None:
    text = f"""Ammons six-dog single-cell localization
Script version: {SCRIPT_VERSION}

Why this analysis is required
-----------------------------
The official study generated eight raw scRNA-seq samples from six dogs.
Dogs 1 and 2 each had two technical replicate libraries. The Cell Browser
metadata reflects this structure:
- 8 orig.ident library identifiers
- 6 name values

This script combines technical replicate libraries by name before any
biological inference. The earlier eight-orig.ident targeted tests are retained
only for audit comparison and are superseded by this six-dog analysis.

Statistical implication
-----------------------
For six paired biological dogs, the minimum attainable two-sided exact sign-flip
P value is 2 / 2^6 = 0.03125.

Mapping
-------
{mapping.to_string(index=False)}

Corrected targeted results
--------------------------
{targeted.to_string(index=False)}

Guardrails
----------
- Dog is the biological replicate.
- Individual cells and technical libraries are not independent replicates.
- The four targeted contrasts are biologically related and are not four
  independent replications.
- Single-cell localization supports biological interpretation, not causal or
  external prognostic validation.
"""
    OUTPUT_README.write_text(text, encoding="utf-8")
```

#### `main`

- Lines: 1114-1472
- Signals: `direction, coverage`

```python
def main() -> None:
    print("=" * 80)
    print("Recompute Ammons localization using six biological dogs")
    print("=" * 80)
    print(f"Script version: {SCRIPT_VERSION}")
    print(f"Project root: {PROJECT_ROOT}")
    print("")
    print("Design:")
    print("  Collapse eight Cell Ranger libraries to six dogs using metadata name.")
    print("  Combine technical replicate cells before pseudobulk aggregation.")
    print("  Recalculate all four prespecified paired contrasts.")
    print("  Replace eight-orig.ident P values with six-dog exact inference.")
    print("")

    metadata = read_required_csv(
        METADATA_FILE,
        sep="\t",
        dtype=str,
        low_memory=False,
    )
    cell_scores = read_required_csv(
        CELL_SCORES_FILE,
        compression="gzip",
        low_memory=False,
    )
    coverage = read_required_csv(SCORE_COVERAGE_FILE)
    original_targeted = read_required_csv(
        ORIGINAL_TARGETED_FILE
    )
    script42_manifest = read_required_json(
        SCRIPT42_MANIFEST_FILE
    )

    if bool(script42_manifest.get("outcome_loaded", True)):
        raise RuntimeError(
            "Script 42 manifest does not confirm outcome_loaded=false."
        )

    mapping, dog_counts = validate_six_dog_mapping(metadata)
    cells = merge_dog_ids(metadata, cell_scores)

    if cells["dog_id"].nunique() != 6:
        raise RuntimeError(
            "Merged cell-score table does not contain six dogs."
        )

    score_columns = estimable_score_columns(coverage)
    missing_scores = [
        column
        for column in score_columns
        if column not in cells.columns
    ]
    if missing_scores:
        raise ValueError(
            f"Expected score columns are missing: {missing_scores}"
        )

    dog_celltype_tables = []
    for annotation_level in ANNOTATION_LEVELS:
        aggregated = aggregate_scores(
            cells=cells,
            group_columns=["dog_id", annotation_level],
            score_columns=score_columns,
            minimum_cells=MIN_CELLS_PER_DOG_CELLTYPE,
        )
        aggregated = aggregated.rename(
            columns={annotation_level: "cell_type"}
        )
        aggregated.insert(
            1,
            "annotation_level",
            annotation_level,
        )
        dog_celltype_tables.append(aggregated)

    dog_celltype = pd.concat(
        dog_celltype_tables,
        ignore_index=True,
    )

    l1_scores = dog_celltype[
        dog_celltype["annotation_level"].eq("celltype.l1")
    ].rename(columns={"cell_type": "celltype.l1"})

    compartment_cells = cells.copy()
    if "immune_combined" not in compartment_cells.columns:
        raise ValueError(
            "The saved cell-score table is missing immune_combined."
        )
    compartment_cells["analysis_compartment"] = (
        compartment_cells["immune_combined"]
    )

    compartment_scores = aggregate_scores(
        cells=compartment_cells,
        group_columns=["dog_id", "analysis_compartment"],
        score_columns=score_columns,
        minimum_cells=MIN_CELLS_PER_DOG_COMPARTMENT,
    )

    targeted = targeted_tests(
        compartment_scores=compartment_scores,
        l1_scores=l1_scores,
    )
    comparison = compare_with_eight_sample(
        six_dog=targeted,
        eight_sample=original_targeted,
    )
    summary = celltype_summary(
        dog_celltype=dog_celltype,
        score_columns=score_columns,
    )
    stability = rank_stability(
        dog_celltype=dog_celltype,
        score_columns=score_columns,
    )
    locked = build_locked_single_cell(targeted)

    mapping.to_csv(OUTPUT_MAPPING, index=False)
    dog_counts.to_csv(OUTPUT_DOG_COUNTS, index=False)
    dog_celltype.to_csv(OUTPUT_DOG_CELLTYPE, index=False)
    compartment_scores.to_csv(
        OUTPUT_DOG_COMPARTMENT,
        index=False,
    )
    targeted.to_csv(OUTPUT_TARGETED, index=False)
    comparison.to_csv(OUTPUT_COMPARISON, index=False)
    summary.to_csv(OUTPUT_CELLTYPE_SUMMARY, index=False)
    stability.to_csv(OUTPUT_RANK_STABILITY, index=False)
    locked.to_csv(OUTPUT_LOCKED_SINGLE_CELL, index=False)

    if MULTIDIMENSIONAL_FILE.exists():
        multidimensional = pd.read_csv(
            MULTIDIMENSIONAL_FILE
        )
        updated = multidimensional.merge(
            locked[
                [
                    "module_label",
                    "single_cell_localization_class",
                    "n_biological_dogs",
                    "locked_single_cell_interpretation",
                    "replicate_guardrail",
                ]
            ],
            on="module_label",
            how="left",
        )
        updated.to_csv(
            OUTPUT_UPDATED_MASTER,
            index=False,
        )
    else:
        updated = pd.DataFrame()

    paired_plot(
        pseudobulk=compartment_scores,
        group_column="analysis_compartment",
        group_a="osteoblast_lineage",
        group_b="immune_combined",
        score_column="M34_signed_risk_score_mean",
        title=(
            "M34 signed risk: osteoblast lineage versus immune "
            "(six dogs)"
        ),
        png_path=M34_PAIRED_PNG,
        pdf_path=M34_PAIRED_PDF,
    )
    paired_plot(
        pseudobulk=l1_scores,
        group_column="celltype.l1",
        group_a="Osteoblast_cycling",
        group_b="Osteoblast",
        score_column="M40_signed_risk_score_mean",
        title=(
            "M40 signed risk: cycling versus non-cycling "
            "osteoblast (six dogs)"
        ),
        png_path=M40_PAIRED_PNG,
        pdf_path=M40_PAIRED_PDF,
    )

    write_sentences(locked, targeted)
    write_readme(mapping, targeted)

    input_paths = [
        METADATA_FILE,
        CELL_SCORES_FILE,
        SCORE_COVERAGE_FILE,
        ORIGINAL_TARGETED_FILE,
        SCRIPT42_MANIFEST_FILE,
    ]
    if MULTIDIMENSIONAL_FILE.exists():
        input_paths.append(MULTIDIMENSIONAL_FILE)

    output_paths = [
        OUTPUT_MAPPING,
        OUTPUT_DOG_COUNTS,
        OUTPUT_DOG_CELLTYPE,
        OUTPUT_DOG_COMPARTMENT,
        OUTPUT_TARGETED,
        OUTPUT_COMPARISON,
        OUTPUT_CELLTYPE_SUMMARY,
        OUTPUT_RANK_STABILITY,
        OUTPUT_LOCKED_SINGLE_CELL,
        OUTPUT_SENTENCES,
        OUTPUT_README,
        M34_PAIRED_PNG,
        M34_PAIRED_PDF,
        M40_PAIRED_PNG,
        M40_PAIRED_PDF,
    ]
    if not updated.empty:
        output_paths.append(OUTPUT_UPDATED_MASTER)

    manifest = {
        "script_version": SCRIPT_VERSION,
        "created_utc": utc_now_iso(),
        "outcome_loaded": False,
        "biological_replicate_column": DOG_COLUMN,
        "n_cell_ranger_libraries": 8,
        "n_biological_dogs": 6,
        "technical_replicate_structure": {
            "library_counts_per_dog": (
                mapping.groupby("dog_id")[LIBRARY_COLUMN]
                .nunique()
                .to_dict()
            ),
        },
        "supersedes": [
            "Ammons_scRNA_targeted_localization_tests.csv",
            "paper4_locked_single_cell_biological_localization.csv",
        ],
        "guardrails": [
            "Dogs 1 and 2 had two technical replicate libraries.",
            "Technical replicate cells were combined before dog-level pseudobulk inference.",
            "The six-dog analysis supersedes the eight-orig.ident exact tests.",
            "No clinical outcome or endpoint was loaded.",
            "Single-cell localization is biological annotation, not prognostic validation.",
        ],
        "inputs": {
            str(path): {
                "sha256": sha256_file(path),
                "size_bytes": path.stat().st_size,
            }
            for path in input_paths
        },
        "outputs": {
            str(path): {
                "sha256": sha256_file(path),
                "size_bytes": path.stat().st_size,
            }
            for path in output_paths
            if path.exists()
        },
    }
    OUTPUT_MANIFEST.write_text(
        json.dumps(manifest, indent=2),
        encoding="utf-8",
    )

    print("=" * 80)
    print("Confirmed library-to-dog mapping")
    print("=" * 80)
    print(mapping.to_string(index=False))

    print("")
    print("=" * 80)
    print("Corrected six-dog targeted tests")
    print("=" * 80)
    print(
        targeted[
            [
                "contrast_name",
                "n_paired_dogs",
                "n_positive_differences",
                "n_negative_differences",
                "mean_paired_difference",
                "median_paired_difference",
                "bootstrap_ci_low",
                "bootstrap_ci_high",
                "exact_sign_flip_p",
                "minimum_attainable_two_sided_p",
                "primary_bh_q",
            ]
        ].to_string(index=False)
    )

    print("")
    print("=" * 80)
    print("Eight-sample versus six-dog comparison")
    print("=" * 80)
    print(
        comparison[
            [
                "contrast_name",
                "eight_sample_mean_difference",
                "six_dog_mean_difference",
                "mean_effect_change",
                "eight_sample_exact_p",
                "six_dog_exact_p",
                "six_dog_bh_q",
                "all_same_nonzero_sign",
            ]
        ].to_string(index=False)
    )

    print("")
    print("=" * 80)
    print("Corrected locked single-cell localization")
    print("=" * 80)
    print(
        locked[
            [
                "module_label",
                "single_cell_localization_class",
                "n_biological_dogs",
                "primary_effect_1",
                "primary_q_1",
                "primary_effect_2",
                "primary_q_2",
            ]
        ].to_string(index=False)
    )

    print("")
    print("=" * 80)
    print("Interpretation guardrails")
    print("=" * 80)
    print("The metadata name column defines six biological dogs.")
    print("Dogs 1 and 2 each contributed two technical replicate libraries.")
    print("The earlier eight-orig.ident P=0.0078125 results are superseded.")
    print("For six paired dogs, the exact two-sided floor is P=0.03125.")
    print("Biological interpretation remains valid only if dog-level directions remain concordant.")

    print("")
    print("Saved:")
    for path in [
        OUTPUT_MAPPING,
        OUTPUT_DOG_COUNTS,
        OUTPUT_DOG_CELLTYPE,
        OUTPUT_DOG_COMPARTMENT,
        OUTPUT_TARGETED,
        OUTPUT_COMPARISON,
        OUTPUT_CELLTYPE_SUMMARY,
        OUTPUT_RANK_STABILITY,
        OUTPUT_LOCKED_SINGLE_CELL,
        OUTPUT_UPDATED_MASTER,
        OUTPUT_SENTENCES,
        OUTPUT_README,
        OUTPUT_MANIFEST,
        M34_PAIRED_PNG,
        M34_PAIRED_PDF,
        M40_PAIRED_PNG,
        M40_PAIRED_PDF,
    ]:
        if path.exists():
            print(path)
    print("Done.")
```

### Relevant standalone lines

| Line | Signals | Code |
|---:|---|---|
| 188 | mean | `observed = float(np.mean(values))` |
| 194 | mean | `np.mean(` |
| 202 | mean | `np.mean(` |
| 225 | mean | `means = values[indices].mean(axis=1)` |
| 234 | coverage | `coverage: pd.DataFrame,` |
| 238 | coverage | `for row in coverage.itertuples(index=False):` |
| 382 | mean | `means = grouped[score_columns].mean()` |
| 495 | direction | `"M34 signed risk: osteoblast lineage versus immune"` |
| 508 | loadings | `"M34 negative-loading expression: immune versus "` |
| 520 | direction | `"M40 signed risk: cycling versus non-cycling "` |
| 534 | loadings | `"M40 positive-loading expression: cycling versus "` |
| 594 | mean | `np.mean(values)` |
| 642 | mean | `.mean()` |
| 655 | mean | `.mean()` |
| 805 | direction | `"M34 signed risk: osteoblast lineage versus immune"` |
| 808 | loadings | `"M34 negative-loading expression: immune versus "` |
| 812 | direction | `"M40 signed risk: cycling versus non-cycling osteoblast"` |
| 815 | loadings | `"M40 positive-loading expression: cycling versus "` |
| 846 | direction | `"Across six biological dogs, the M34 signed "` |
| 849 | loadings | `"the negative-loading component was consistently "` |
| 852 | loadings | `"loadings, M34 is best interpreted as an inverse "` |
| 887 | direction | `"Across six biological dogs, M40 signed and "` |
| 888 | loadings | `"positive-loading scores were consistently higher in "` |
| 1047 | direction | `"direction in all six biological dogs. The exact "` |
| 1048 | direction | `"two-sided sign-flip P value was 0.03125, the minimum "` |
| 1091 | direction | `For six paired biological dogs, the minimum attainable two-sided exact sign-flip` |
| 1139 | coverage | `coverage = read_required_csv(SCORE_COVERAGE_FILE)` |
| 1160 | coverage | `score_columns = estimable_score_columns(coverage)` |
| 1276 | direction | `"M34 signed risk: osteoblast lineage versus immune "` |
| 1289 | direction | `"M40 signed risk: cycling versus non-cycling "` |

### Referenced result-table schemas

#### `results/tables/Ammons_scRNA_celltype_rank_stability_six_dogs.csv`

- Rows: 820
- SHA-256: `5b0af4bcf4c1e80c506d9b62a458458e7ddab6e1d04855b428557e72271fbd66`
- Columns: `annotation_level`, `score_column`, `cell_type`, `full_rank`, `minimum_loo_rank`, `maximum_loo_rank`, `median_loo_rank`, `maximum_absolute_rank_change`

#### `results/tables/Ammons_scRNA_celltype_score_summary_six_dogs.csv`

- Rows: 820
- SHA-256: `ea1675623fe6c9da7370ccbca8ad1ff0d9546957f0bcef25923a4ae7831b3ea9`
- Columns: `annotation_level`, `cell_type`, `score_column`, `n_dogs`, `total_cells`, `median_dog_score`, `mean_dog_score`, `iqr_low`, `iqr_high`, `minimum_dog_score`, `maximum_dog_score`

#### `results/tables/Ammons_scRNA_confirmed_library_to_dog_mapping.csv`

- Rows: 8
- SHA-256: `19432f4696ef8fd4e11f598f8fd4fe8b793d8edc9e396489e35ad0b4d326c08e`
- Columns: `orig.ident`, `dog_id`

#### `results/tables/Ammons_scRNA_score_gene_coverage.csv`

- Rows: 4
- SHA-256: `13ae84b614bc12dd08ef47485a5cf70ff56c74df6c43ca91d7806ee3a283cf49`
- Columns: `module_label`, `n_frozen_genes`, `n_detected_genes`, `coverage_fraction`, `n_positive_detected`, `n_negative_detected`, `signed_score_estimable`, `positive_component_estimable`, `negative_component_estimable`, `component_guardrail`

#### `results/tables/Ammons_scRNA_six_dog_cell_counts.csv`

- Rows: 6
- SHA-256: `92951d6b30a238b3856e009afe7c64427eba895959a003418ca0f78d08e9fe0c`
- Columns: `dog_id`, `n_cells`, `fraction_cells`, `n_libraries`

#### `results/tables/Ammons_scRNA_six_dog_celltype_pseudobulk_scores.csv`

- Rows: 430
- SHA-256: `976007babfc4f6f407a40d673c8a6fe4bb0317c8ddaf66b82970a67de174adf6`
- Columns: `dog_id`, `annotation_level`, `cell_type`, `n_cells`, `M11_positive_component_expression_mean`, `M11_signed_risk_score_mean`, `M24_positive_component_expression_mean`, `M24_signed_risk_score_mean`, `M34_negative_component_expression_mean`, `M34_positive_component_expression_mean`, `M34_signed_risk_score_mean`, `M40_negative_component_expression_mean`, `M40_positive_component_expression_mean`, `M40_signed_risk_score_mean`, `M11_positive_component_expression_median`, `M11_signed_risk_score_median`, `M24_positive_component_expression_median`, `M24_signed_risk_score_median`, `M34_negative_component_expression_median`, `M34_positive_component_expression_median`, `M34_signed_risk_score_median`, `M40_negative_component_expression_median`, `M40_positive_component_expression_median`, `M40_signed_risk_score_median`

#### `results/tables/Ammons_scRNA_six_dog_compartment_scores.csv`

- Rows: 29
- SHA-256: `a3dda4b5a0bc986c3404427c4dba5007a682534856d0122b30069606b1d165ba`
- Columns: `dog_id`, `analysis_compartment`, `n_cells`, `M11_positive_component_expression_mean`, `M11_signed_risk_score_mean`, `M24_positive_component_expression_mean`, `M24_signed_risk_score_mean`, `M34_negative_component_expression_mean`, `M34_positive_component_expression_mean`, `M34_signed_risk_score_mean`, `M40_negative_component_expression_mean`, `M40_positive_component_expression_mean`, `M40_signed_risk_score_mean`, `M11_positive_component_expression_median`, `M11_signed_risk_score_median`, `M24_positive_component_expression_median`, `M24_signed_risk_score_median`, `M34_negative_component_expression_median`, `M34_positive_component_expression_median`, `M34_signed_risk_score_median`, `M40_negative_component_expression_median`, `M40_positive_component_expression_median`, `M40_signed_risk_score_median`

#### `results/tables/Ammons_scRNA_targeted_eight_sample_vs_six_dog_comparison.csv`

- Rows: 4
- SHA-256: `c198453f93f0005a34c4603b80df5b65c4d2d798980e5631188e057f26013630`
- Columns: `contrast_name`, `n_orig_ident_units`, `eight_sample_mean_difference`, `eight_sample_median_difference`, `eight_sample_ci_low`, `eight_sample_ci_high`, `eight_sample_exact_p`, `eight_sample_bh_q`, `n_biological_dogs`, `six_dog_mean_difference`, `six_dog_median_difference`, `six_dog_ci_low`, `six_dog_ci_high`, `six_dog_exact_p`, `six_dog_bh_q`, `all_same_nonzero_sign`, `mean_effect_change`, `mean_effect_ratio`, `supersession_note`

#### `results/tables/Ammons_scRNA_targeted_localization_tests.csv`

- Rows: 4
- SHA-256: `5f75eac6ccc4f12b710c02f4220e78ee13e5692a7c694097b4fbb0e9049f2203`
- Columns: `contrast_name`, `score_column`, `group_column`, `group_a`, `group_b`, `difference_definition`, `hypothesis_direction`, `n_paired_dogs`, `mean_paired_difference`, `median_paired_difference`, `bootstrap_ci_low`, `bootstrap_ci_high`, `exact_sign_flip_p`, `dog_level_differences`, `estimable`, `primary_bh_q`

#### `results/tables/Ammons_scRNA_targeted_localization_tests_six_dogs.csv`

- Rows: 4
- SHA-256: `acc4f4b2486ddd7beffd15560bc76c1c993d7c56e4c21f6653bed0b8ca12c5fb`
- Columns: `contrast_name`, `score_column`, `group_column`, `group_a`, `group_b`, `difference_definition`, `hypothesis_direction`, `n_paired_dogs`, `n_positive_differences`, `n_negative_differences`, `n_zero_differences`, `all_same_nonzero_sign`, `mean_paired_difference`, `median_paired_difference`, `bootstrap_ci_low`, `bootstrap_ci_high`, `exact_sign_flip_p`, `minimum_attainable_two_sided_p`, `dog_level_differences`, `estimable`, `primary_bh_q`

#### `results/tables/paper4_locked_multidimensional_transport_evidence.csv`

- Rows: 4
- SHA-256: `d6f00375286e78e184a99c0774290a91cd502771595beb815907cb57f6bbcbf5`
- Columns: `module_label`, `n_strong_structure_settings`, `n_partial_structure_settings`, `n_limited_structure_settings`, `structure_evidence_class`, `variable_only_latent_class`, `targeted_mofa_representation_class`, `core_natural_coverage_fraction`, `core_ubiquitous_median_capture`, `core_non_ffpe_median_capture`, `detection_aware_ubiquitous_median_capture`, `no_ffpe_ubiquitous_median_capture`, `core_median_max_absolute_factor_cosine`, `n_projectwide_fdr_supported_human_settings`, `n_nominal_supported_human_settings`, `minimum_projectwide_q_12`, `outcome_evidence_class`, `gse39055_assay_stability_class`, `multidimensional_transport_class`, `locked_multidimensional_interpretation`, `structure_evidence_score`, `independent_latent_recurrence_score`, `outcome_transport_evidence_score`, `gse39055_assay_stability_score`

#### `results/tables/paper4_locked_multidimensional_transport_evidence_with_single_cell_six_dogs.csv`

- Rows: 4
- SHA-256: `a096f3381e2854bf49499056ddd69f59791dd544fa94e61aeed01d0b80dc6537`
- Columns: `module_label`, `n_strong_structure_settings`, `n_partial_structure_settings`, `n_limited_structure_settings`, `structure_evidence_class`, `variable_only_latent_class`, `targeted_mofa_representation_class`, `core_natural_coverage_fraction`, `core_ubiquitous_median_capture`, `core_non_ffpe_median_capture`, `detection_aware_ubiquitous_median_capture`, `no_ffpe_ubiquitous_median_capture`, `core_median_max_absolute_factor_cosine`, `n_projectwide_fdr_supported_human_settings`, `n_nominal_supported_human_settings`, `minimum_projectwide_q_12`, `outcome_evidence_class`, `gse39055_assay_stability_class`, `multidimensional_transport_class`, `locked_multidimensional_interpretation`, `structure_evidence_score`, `independent_latent_recurrence_score`, `outcome_transport_evidence_score`, `gse39055_assay_stability_score`, `single_cell_localization_class`, `n_biological_dogs`, `locked_single_cell_interpretation`, `replicate_guardrail`

#### `results/tables/paper4_locked_single_cell_biological_localization.csv`

- Rows: 4
- SHA-256: `9ac746cc177aad28067e358d538e845661f757d3c2b8f4f7ebbae4e242b694fa`
- Columns: `module_label`, `single_cell_localization_class`, `primary_contrast_1`, `primary_effect_1`, `primary_q_1`, `primary_contrast_2`, `primary_effect_2`, `primary_q_2`, `locked_single_cell_interpretation`, `replicate_guardrail`

#### `results/tables/paper4_locked_single_cell_biological_localization_six_dogs.csv`

- Rows: 4
- SHA-256: `511b1817ed542c1ba0d26d70b5842f38906a3f356d86de6210c64ca02dd20633`
- Columns: `module_label`, `single_cell_localization_class`, `n_biological_dogs`, `primary_contrast_1`, `primary_effect_1`, `primary_q_1`, `primary_contrast_2`, `primary_effect_2`, `primary_q_2`, `locked_single_cell_interpretation`, `replicate_guardrail`

---

## `scripts/46_gse239948_external_canine_representation_v2.py`

- SHA-256: `eade633a0d4a11419af829380a46e92dac3adf17a4d2d09c6691f6df28aca106`

### Relevant functions

#### `sample_columns`

- Lines: 235-262
- Signals: `mean`

```python
def sample_columns(table: pd.DataFrame) -> list[str]:
    columns = [str(column) for column in table.columns]
    matched = [
        column
        for column in columns
        if re.match(r"^CCB\d+", column, flags=re.IGNORECASE)
        or re.match(r"^GSM\d+", column, flags=re.IGNORECASE)
    ]

    if len(matched) >= 30:
        return matched

    numeric_fraction = {}
    for column in columns:
        numeric = pd.to_numeric(table[column], errors="coerce")
        numeric_fraction[column] = float(numeric.notna().mean())

    candidates = [
        column
        for column in columns
        if numeric_fraction[column] >= 0.95
    ]

    if len(candidates) < 30:
        raise RuntimeError(
            "Could not identify at least 30 numeric sample columns."
        )
    return candidates
```

#### `frozen_identifier_columns`

- Lines: 328-356
- Signals: `weights`

```python
def frozen_identifier_columns(weights: pd.DataFrame) -> list[str]:
    preferred = [
        "canine_gene_symbol",
        "canine_gene",
        "dog_gene_symbol",
        "dog_gene",
        "canine_ensembl_gene_id",
        "dog_ensembl_gene_id",
        "ensembl_gene_id",
    ]
    columns = []

    for column in preferred:
        if column in weights.columns:
            columns.append(column)

    for column in weights.columns:
        lower = str(column).lower()
        if (
            column not in columns
            and any(token in lower for token in ["canine", "dog"])
            and any(
                token in lower
                for token in ["gene", "symbol", "ensembl", "id"]
            )
        ):
            columns.append(str(column))

    return columns
```

#### `identifier_pair_audit`

- Lines: 359-478
- Signals: `weights`

```python
def identifier_pair_audit(
    table: pd.DataFrame,
    sample_cols: list[str],
    weights: pd.DataFrame,
) -> pd.DataFrame:
    external_columns = [
        str(column)
        for column in table.columns
        if str(column) not in sample_cols
    ]
    frozen_columns = frozen_identifier_columns(weights)

    if not external_columns:
        raise RuntimeError(
            "No non-sample columns were available for gene identifiers."
        )
    if not frozen_columns:
        raise RuntimeError(
            "No canine identifier column was found in the frozen weights."
        )

    primary_weights = weights[
        weights["module_label"].isin(PRIMARY_MODULES)
    ].copy()

    rows = []
    for frozen_column in frozen_columns:
        for scheme in ["symbol", "ensembl", "raw"]:
            frozen_keys = primary_weights[frozen_column].map(
                lambda value: normalize_identifier(value, scheme)
            )
            frozen_keys = frozen_keys[frozen_keys.ne("")]

            if frozen_keys.empty:
                continue

            frozen_set = set(frozen_keys)
            module_key_sets = {
                module: set(
                    primary_weights.loc[
                        primary_weights["module_label"].eq(module),
                        frozen_column,
                    ].map(
                        lambda value: normalize_identifier(value, scheme)
                    )
                )
                - {""}
                for module in PRIMARY_MODULES
            }

            for external_column in external_columns:
                external_keys = table[external_column].map(
                    lambda value: normalize_identifier(value, scheme)
                )
                external_set = set(external_keys) - {""}
                overlap_set = frozen_set.intersection(external_set)

                module_overlaps = {
                    module: len(
                        module_key_sets[module].intersection(external_set)
                    )
                    for module in PRIMARY_MODULES
                }
                modules_ge3 = sum(
                    count >= MIN_MODULE_GENES
                    for count in module_overlaps.values()
                )

                rows.append(
                    {
                        "external_column": external_column,
                        "frozen_identifier_column": frozen_column,
                        "identifier_scheme": scheme,
                        "n_external_unique_ids": len(external_set),
                        "n_frozen_unique_ids": len(frozen_set),
                        "n_unique_overlap": len(overlap_set),
                        "n_primary_modules_with_at_least_3_genes": modules_ge3,
                        **{
                            f"{module}_overlap": module_overlaps[module]
                            for module in PRIMARY_MODULES
                        },
                        "example_overlap": ";".join(
                            sorted(overlap_set)[:20]
                        ),
                    }
                )

    if not rows:
        raise RuntimeError(
            "No usable external/frozen identifier combinations were found."
        )

    audit = pd.DataFrame(rows)
    preference = {
        "canine_gene_symbol": 4,
        "canine_gene": 3,
        "dog_gene_symbol": 2,
        "dog_gene": 1,
    }
    scheme_preference = {
        "symbol": 3,
        "ensembl": 2,
        "raw": 1,
    }
    audit["frozen_column_preference"] = audit[
        "frozen_identifier_column"
    ].map(preference).fillna(0)
    audit["scheme_preference"] = audit[
        "identifier_scheme"
    ].map(scheme_preference).fillna(0)

    return audit.sort_values(
        [
            "n_primary_modules_with_at_least_3_genes",
            "n_unique_overlap",
            "frozen_column_preference",
            "scheme_preference",
        ],
        ascending=[False, False, False, False],
    ).reset_index(drop=True)
```

#### `remap_reference_expression`

- Lines: 552-612
- Signals: `weights`

```python
def remap_reference_expression(
    reference: pd.DataFrame,
    weights: pd.DataFrame,
    frozen_identifier_column: str,
    identifier_scheme: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    if "canine_gene" not in weights.columns:
        raise ValueError(
            "Frozen weights are missing the canine_gene column needed "
            "to map the DOG2 reference expression matrix."
        )

    mapping = weights[
        ["canine_gene", frozen_identifier_column]
    ].copy()
    mapping["canine_gene"] = mapping["canine_gene"].astype(str)
    mapping["analysis_gene_id"] = mapping[
        frozen_identifier_column
    ].map(
        lambda value: normalize_identifier(
            value,
            identifier_scheme,
        )
    )
    mapping = mapping[
        mapping["analysis_gene_id"].ne("")
    ].drop_duplicates(
        ["canine_gene", "analysis_gene_id"],
        keep="first",
    )

    exact_map = (
        mapping.drop_duplicates("canine_gene")
        .set_index("canine_gene")["analysis_gene_id"]
        .to_dict()
    )

    rename_map = {}
    mapping_rows = []
    for column in reference.columns.astype(str):
        analysis_id = exact_map.get(column, "")
        if not analysis_id and identifier_scheme == "symbol":
            analysis_id = normalize_identifier(column, "symbol")
        if not analysis_id:
            continue
        rename_map[column] = analysis_id
        mapping_rows.append(
            {
                "reference_expression_column": column,
                "analysis_gene_id": analysis_id,
            }
        )

    remapped = reference.rename(columns=rename_map)
    remapped = remapped.loc[:, list(rename_map.values())]
    remapped = remapped.loc[
        :,
        ~pd.Index(remapped.columns).duplicated(keep="first"),
    ]

    return remapped, pd.DataFrame(mapping_rows)
```

#### `choose_transform`

- Lines: 613-641
- Signals: `mean`

```python
def choose_transform(expression: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    values = expression.to_numpy(dtype=float)
    finite = values[np.isfinite(values)]

    if finite.size == 0:
        raise RuntimeError("The external expression matrix is empty.")

    minimum = float(np.min(finite))
    maximum = float(np.max(finite))
    median = float(np.median(finite))
    integer_like_fraction = float(
        np.mean(np.isclose(finite, np.round(finite)))
    )

    if minimum >= 0 and (maximum > 50 or median > 10):
        transformed = np.log2(expression.clip(lower=0) + 1.0)
        transform = "log2_x_plus_1"
    else:
        transformed = expression.copy()
        transform = "as_provided"

    diagnostics = {
        "minimum": minimum,
        "median": median,
        "maximum": maximum,
        "integer_like_fraction": integer_like_fraction,
        "transform": transform,
    }
    return transformed, diagnostics
```

#### `zscore_columns`

- Lines: 644-649
- Signals: `mean, std`

```python
def zscore_columns(expression: pd.DataFrame) -> pd.DataFrame:
    x = expression.apply(pd.to_numeric, errors="coerce")
    x = x.fillna(x.median(axis=0))
    std = x.std(axis=0).replace(0, np.nan)
    z = (x - x.mean(axis=0)) / std
    return z.loc[:, z.notna().all(axis=0)]
```

#### `loading_permutation_p`

- Lines: 687-710
- Signals: ``

```python
def loading_permutation_p(
    frozen_loadings: np.ndarray,
    external_loadings: np.ndarray,
    observed: float,
    seed: int,
) -> float:
    if not np.isfinite(observed):
        return np.nan
    rng = np.random.default_rng(seed)
    count = 0

    for _ in range(N_GENE_LABEL_PERMUTATIONS):
        permuted = rng.permutation(external_loadings)
        value = stats.spearmanr(
            frozen_loadings,
            permuted,
        ).statistic

        if abs(value) >= abs(observed):
            count += 1

    return float(
        (count + 1) / (N_GENE_LABEL_PERMUTATIONS + 1)
    )
```

#### `frozen_signed_score`

- Lines: 713-723
- Signals: `mean, loadings, direction`

```python
def frozen_signed_score(
    expression_z: pd.DataFrame,
    genes: list[str],
    loadings: pd.Series,
) -> pd.Series:
    signs = np.sign(loadings.reindex(genes)).replace(0, 1)
    score = expression_z[genes].mul(
        signs,
        axis=1,
    ).mean(axis=1)
    return score
```

#### `split_half_reliability`

- Lines: 726-782
- Signals: `loadings, correlation`

```python
def split_half_reliability(
    expression_z: pd.DataFrame,
    genes: list[str],
    loadings: pd.Series,
    seed: int,
) -> dict[str, float]:
    if len(genes) < 4:
        return {
            "split_half_median": np.nan,
            "split_half_q05": np.nan,
            "split_half_q95": np.nan,
            "split_half_valid_repeats": 0,
        }

    rng = np.random.default_rng(seed)
    correlations = []

    for _ in range(N_SPLIT_HALF_REPEATS):
        shuffled = np.asarray(genes, dtype=object)
        shuffled = rng.permutation(shuffled)
        midpoint = len(shuffled) // 2
        first = shuffled[:midpoint].tolist()
        second = shuffled[midpoint:].tolist()

        first_score = frozen_signed_score(
            expression_z,
            first,
            loadings,
        )
        second_score = frozen_signed_score(
            expression_z,
            second,
            loadings,
        )
        correlation = first_score.corr(second_score)

        if np.isfinite(correlation):
            correlations.append(float(correlation))

    if not correlations:
        return {
            "split_half_median": np.nan,
            "split_half_q05": np.nan,
            "split_half_q95": np.nan,
            "split_half_valid_repeats": 0,
        }

    return {
        "split_half_median": float(np.median(correlations)),
        "split_half_q05": float(
            np.quantile(correlations, 0.05)
        ),
        "split_half_q95": float(
            np.quantile(correlations, 0.95)
        ),
        "split_half_valid_repeats": len(correlations),
    }
```

#### `gene_leave_one_out`

- Lines: 785-819
- Signals: `loadings, correlation`

```python
def gene_leave_one_out(
    expression_z: pd.DataFrame,
    module: str,
    genes: list[str],
    loadings: pd.Series,
) -> pd.DataFrame:
    full_score = frozen_signed_score(
        expression_z,
        genes,
        loadings,
    )
    rows = []

    if len(genes) <= MIN_MODULE_GENES:
        return pd.DataFrame()

    for gene in genes:
        subset = [item for item in genes if item != gene]
        score = frozen_signed_score(
            expression_z,
            subset,
            loadings,
        )
        rows.append(
            {
                "module_label": module,
                "left_out_gene": gene,
                "n_genes_remaining": len(subset),
                "correlation_with_full_score": float(
                    full_score.corr(score)
                ),
            }
        )

    return pd.DataFrame(rows)
```

#### `analyze_module`

- Lines: 822-984
- Signals: `weights, correlation, coverage`

```python
def analyze_module(
    module: str,
    reference_z: pd.DataFrame,
    external_z: pd.DataFrame,
    weights: pd.DataFrame,
    seed: int,
) -> tuple[dict[str, Any], dict[str, Any], pd.DataFrame, pd.Series]:
    part = weights[
        weights["module_label"].astype(str).eq(module)
    ].copy()
    part["canine_gene_symbol"] = (
        part["canine_gene_symbol"].map(clean_symbol)
    )
    part = part.drop_duplicates(
        "canine_gene_symbol",
        keep="first",
    )
    part = part.set_index("canine_gene_symbol")

    genes = [
        gene
        for gene in part.index
        if gene in reference_z.columns
        and gene in external_z.columns
    ]

    coverage = {
        "module_label": module,
        "n_frozen_genes": int(part.shape[0]),
        "n_common_genes": len(genes),
        "coverage_fraction": (
            len(genes) / part.shape[0]
            if part.shape[0]
            else np.nan
        ),
        "common_genes": ";".join(genes),
        "missing_external_genes": ";".join(
            gene
            for gene in part.index
            if gene not in external_z.columns
        ),
    }

    if len(genes) < MIN_MODULE_GENES:
        result = {
            "module_label": module,
            "n_common_genes": len(genes),
            "edge_spearman": np.nan,
            "edge_permutation_p": np.nan,
            "loading_spearman": np.nan,
            "loading_permutation_p": np.nan,
            "external_pc1_variance_explained": np.nan,
            "pc1_orientation_correlation_with_frozen_score": np.nan,
            "split_half_median": np.nan,
            "split_half_q05": np.nan,
            "split_half_q95": np.nan,
            "split_half_valid_repeats": 0,
            "estimable": False,
            "nonestimable_reason": "fewer_than_three_common_genes",
        }
        return result, coverage, pd.DataFrame(), pd.Series(
            np.nan,
            index=external_z.index,
        )

    reference_correlation = np.corrcoef(
        reference_z[genes].to_numpy(dtype=float),
        rowvar=False,
    )
    external_correlation = np.corrcoef(
        external_z[genes].to_numpy(dtype=float),
        rowvar=False,
    )

    reference_edges = upper_triangle(reference_correlation)
    external_edges = upper_triangle(external_correlation)
    edge_spearman = float(
        stats.spearmanr(
            reference_edges,
            external_edges,
        ).statistic
    )
    edge_p = permutation_gene_label_p(
        reference_edges,
        external_correlation,
        edge_spearman,
        seed,
    )

    frozen_loadings = pd.to_numeric(
        part.loc[genes, "risk_oriented_loading"],
        errors="coerce",
    ).fillna(0.0)

    signed_score = frozen_signed_score(
        external_z,
        genes,
        frozen_loadings,
    )

    pca = PCA(n_components=1, random_state=RANDOM_SEED)
    pc_score = pd.Series(
        pca.fit_transform(
            external_z[genes].to_numpy(dtype=float)
        ).ravel(),
        index=external_z.index,
    )
    pc_loadings = pca.components_[0].copy()

    orientation_correlation = pc_score.corr(signed_score)
    if np.isfinite(orientation_correlation) and orientation_correlation < 0:
        pc_score = -pc_score
        pc_loadings = -pc_loadings
        orientation_correlation = -orientation_correlation

    loading_spearman = float(
        stats.spearmanr(
            frozen_loadings.to_numpy(dtype=float),
            pc_loadings,
        ).statistic
    )
    loading_p = loading_permutation_p(
        frozen_loadings.to_numpy(dtype=float),
        pc_loadings,
        loading_spearman,
        seed + 1,
    )

    reliability = split_half_reliability(
        external_z,
        genes,
        frozen_loadings,
        seed + 2,
    )

    result = {
        "module_label": module,
        "n_common_genes": len(genes),
        "edge_spearman": edge_spearman,
        "edge_permutation_p": edge_p,
        "loading_spearman": loading_spearman,
        "loading_permutation_p": loading_p,
        "external_pc1_variance_explained": float(
            pca.explained_variance_ratio_[0]
        ),
        "pc1_orientation_correlation_with_frozen_score": (
            float(orientation_correlation)
            if np.isfinite(orientation_correlation)
            else np.nan
        ),
        **reliability,
        "estimable": True,
        "nonestimable_reason": "",
    }

    loo = gene_leave_one_out(
        external_z,
        module,
        genes,
        frozen_loadings,
    )

    return result, coverage, loo, signed_score
```

#### `variability_bins`

- Lines: 987-1007
- Signals: `rank`

```python
def variability_bins(
    reference_z: pd.DataFrame,
    external_z: pd.DataFrame,
) -> pd.Series:
    common = reference_z.columns.intersection(
        external_z.columns
    )
    reference_variance = reference_z[common].var(axis=0)
    external_variance = external_z[common].var(axis=0)

    combined_rank = (
        reference_variance.rank(pct=True)
        + external_variance.rank(pct=True)
    ) / 2.0

    return pd.qcut(
        combined_rank.rank(method="average"),
        q=min(N_VARIABILITY_BINS, len(combined_rank)),
        labels=False,
        duplicates="drop",
    )
```

#### `random_panel_controls`

- Lines: 1010-1147
- Signals: `weights`

```python
def random_panel_controls(
    reference_z: pd.DataFrame,
    external_z: pd.DataFrame,
    weights: pd.DataFrame,
    observed: pd.DataFrame,
) -> pd.DataFrame:
    common = reference_z.columns.intersection(
        external_z.columns
    )
    bins = variability_bins(reference_z, external_z)

    frozen_union = set(
        weights.loc[
            weights["module_label"].isin(PRIMARY_MODULES),
            "canine_gene_symbol",
        ].map(clean_symbol)
    )
    candidate_genes = [
        gene for gene in common if gene not in frozen_union
    ]

    rng = np.random.default_rng(RANDOM_SEED)
    rows = []

    for module_index, module in enumerate(PRIMARY_MODULES):
        part = weights[
            weights["module_label"].astype(str).eq(module)
        ].copy()
        part["canine_gene_symbol"] = part[
            "canine_gene_symbol"
        ].map(clean_symbol)
        part = part.drop_duplicates("canine_gene_symbol")
        target_genes = [
            gene
            for gene in part["canine_gene_symbol"]
            if gene in common
        ]

        if len(target_genes) < MIN_MODULE_GENES:
            continue

        target_bins = bins.reindex(target_genes).to_numpy()
        null_values = []

        for repeat in range(N_RANDOM_PANELS):
            selected = []
            used = set()

            for target_bin in target_bins:
                pool = [
                    gene
                    for gene in candidate_genes
                    if gene not in used
                    and bins.get(gene, np.nan) == target_bin
                ]

                if not pool:
                    pool = [
                        gene
                        for gene in candidate_genes
                        if gene not in used
                    ]

                if not pool:
                    selected = []
                    break

                gene = str(rng.choice(pool))
                selected.append(gene)
                used.add(gene)

            if len(selected) != len(target_genes):
                continue

            reference_correlation = np.corrcoef(
                reference_z[selected].to_numpy(dtype=float),
                rowvar=False,
            )
            external_correlation = np.corrcoef(
                external_z[selected].to_numpy(dtype=float),
                rowvar=False,
            )
            edge = stats.spearmanr(
                upper_triangle(reference_correlation),
                upper_triangle(external_correlation),
            ).statistic

            if np.isfinite(edge):
                null_values.append(float(edge))

        observed_edge = float(
            observed.loc[
                observed["module_label"].eq(module),
                "edge_spearman",
            ].iloc[0]
        )

        null_array = np.asarray(null_values, dtype=float)
        empirical_p = (
            float(
                (
                    1
                    + np.sum(
                        np.abs(null_array) >= abs(observed_edge)
                    )
                )
                / (len(null_array) + 1)
            )
            if len(null_array) and np.isfinite(observed_edge)
            else np.nan
        )

        rows.append(
            {
                "module_label": module,
                "n_module_genes": len(target_genes),
                "n_random_panels": len(null_array),
                "observed_edge_spearman": observed_edge,
                "random_edge_median": (
                    float(np.median(null_array))
                    if len(null_array)
                    else np.nan
                ),
                "random_edge_q05": (
                    float(np.quantile(null_array, 0.05))
                    if len(null_array)
                    else np.nan
                ),
                "random_edge_q95": (
                    float(np.quantile(null_array, 0.95))
                    if len(null_array)
                    else np.nan
                ),
                "random_panel_empirical_p": empirical_p,
            }
        )

    return pd.DataFrame(rows)
```

#### `classify_results`

- Lines: 1150-1239
- Signals: `loadings`

```python
def classify_results(
    structure: pd.DataFrame,
    random_controls: pd.DataFrame,
) -> pd.DataFrame:
    result = structure.merge(
        random_controls[
            [
                "module_label",
                "random_panel_empirical_p",
            ]
        ],
        on="module_label",
        how="left",
    )

    direct_p = pd.concat(
        [
            result[
                ["module_label", "edge_permutation_p"]
            ].rename(
                columns={"edge_permutation_p": "p"}
            ).assign(test="edge"),
            result[
                ["module_label", "loading_permutation_p"]
            ].rename(
                columns={"loading_permutation_p": "p"}
            ).assign(test="loading"),
        ],
        ignore_index=True,
    )
    direct_p["q_bh_8"] = bh_adjust(direct_p["p"])

    result = result.merge(
        direct_p[direct_p["test"].eq("edge")][
            ["module_label", "q_bh_8"]
        ].rename(columns={"q_bh_8": "edge_q_bh_8"}),
        on="module_label",
        how="left",
    )
    result = result.merge(
        direct_p[direct_p["test"].eq("loading")][
            ["module_label", "q_bh_8"]
        ].rename(
            columns={"q_bh_8": "loading_q_bh_8"}
        ),
        on="module_label",
        how="left",
    )

    classes = []
    for row in result.itertuples(index=False):
        edge_supported = bool(
            np.isfinite(row.edge_q_bh_8)
            and row.edge_q_bh_8 < 0.05
        )
        loading_supported = bool(
            np.isfinite(row.loading_q_bh_8)
            and row.loading_q_bh_8 < 0.05
        )
        random_supported = bool(
            np.isfinite(row.random_panel_empirical_p)
            and row.random_panel_empirical_p < 0.05
        )
        reliable = bool(
            np.isfinite(row.split_half_median)
            and row.split_half_median >= 0.60
        )

        if edge_supported and loading_supported and reliable:
            label = "strong_external_canine_representation_preservation"
        elif (
            (edge_supported or loading_supported)
            and reliable
        ):
            label = "partial_external_canine_representation_preservation"
        elif (
            random_supported
            or edge_supported
            or loading_supported
        ):
            label = "limited_specific_external_canine_signal"
        else:
            label = "no_clear_external_canine_representation_preservation"

        classes.append(label)

    result[
        "external_canine_representation_class"
    ] = classes
    return result
```

#### `create_heatmap`

- Lines: 1242-1290
- Signals: `loadings`

```python
def create_heatmap(classification: pd.DataFrame) -> None:
    matrix = classification.set_index(
        "module_label"
    ).reindex(PRIMARY_MODULES)[
        [
            "edge_spearman",
            "loading_spearman",
            "split_half_median",
        ]
    ]

    fig, ax = plt.subplots(figsize=(6.6, 4.4))
    image = ax.imshow(
        matrix.to_numpy(dtype=float),
        aspect="auto",
        vmin=-1,
        vmax=1,
    )
    ax.set_xticks(np.arange(3))
    ax.set_xticklabels(
        [
            "Correlation-edge\npreservation",
            "Frozen-loading\nconcordance",
            "Split-half\nreliability",
        ]
    )
    ax.set_yticks(np.arange(len(PRIMARY_MODULES)))
    ax.set_yticklabels(PRIMARY_MODULES)
    ax.set_title(
        "Independent canine cohort representation preservation"
    )

    for row_index in range(matrix.shape[0]):
        for column_index in range(matrix.shape[1]):
            value = matrix.iloc[row_index, column_index]
            text = "NA" if not np.isfinite(value) else f"{value:.2f}"
            ax.text(
                column_index,
                row_index,
                text,
                ha="center",
                va="center",
            )

    fig.colorbar(image, ax=ax, shrink=0.8)
    fig.tight_layout()
    fig.savefig(HEATMAP_PNG, dpi=300, bbox_inches="tight")
    fig.savefig(HEATMAP_PDF, bbox_inches="tight")
    plt.close(fig)
```

#### `write_readme`

- Lines: 1293-1342
- Signals: `zscore, loadings, direction, correlation`

```python
def write_readme(
    raw_path: Path,
    transform_diagnostics: dict[str, Any],
) -> None:
    text = f"""GSE239948 independent canine representation validation
Script version: {SCRIPT_VERSION}

Dataset
-------
Accession: GSE239948
Samples: expected 43 fresh-frozen canine osteosarcoma tumors.
Processed source: {raw_path}

Transformation
--------------
{json.dumps(transform_diagnostics, indent=2)}

Purpose
-------
Evaluate whether the frozen DOG2 programs preserve their transcriptional
representation in an independent canine osteosarcoma RNA-seq cohort.

No outcome model is fitted because GSE239948 does not provide a standardized
survival endpoint in the processed matrix.

Primary representation metrics
------------------------------
1. Spearman concordance of within-module gene-gene correlation edges between
   DOG2 and GSE239948.
2. Concordance between frozen DOG2 risk-oriented loadings and GSE239948 PC1
   loadings. The external PC1 sign is oriented using the frozen signed score,
   without outcomes.
3. Repeated non-overlapping split-half reliability of the frozen signed score.
4. Gene leave-one-out stability.
5. Variability-matched random-panel specificity controls.

Multiplicity
------------
BH correction is applied jointly across four edge-preservation and four
loading-concordance tests.

Guardrails
----------
- Frozen genes, loadings, signs, and tiers are unchanged.
- No GSE239948 outcome is used.
- GSE239948 tumors underwent different therapies; therefore treatment-related
  distribution differences are possible.
- Representation preservation is not prognostic validation.
"""
    OUTPUT_README.write_text(text, encoding="utf-8")
```

#### `main`

- Lines: 1345-1716
- Signals: `mean, std, weights, loadings, direction, coverage`

```python
def main() -> None:
    print("=" * 80)
    print("GSE239948 independent canine representation validation")
    print("=" * 80)
    print(f"Script version: {SCRIPT_VERSION}")
    print(f"Project root: {PROJECT_ROOT}")
    print("")
    print("Design:")
    print("  Load the independent 43-sample canine osteosarcoma cohort.")
    print("  Match external genes to frozen canine symbols outcome-blind.")
    print("  Test edge preservation, loading concordance, and score reliability.")
    print("  Use gene-label permutations and variability-matched random panels.")
    print("")

    raw_path = locate_or_download_raw_file()

    reference = pd.read_csv(
        REFERENCE_EXPRESSION_FILE,
        index_col=0,
    )
    weights = pd.read_csv(STRICT_WEIGHTS_FILE)
    freeze = (
        json.loads(FREEZE_FILE.read_text(encoding="utf-8"))
        if FREEZE_FILE.exists()
        else {}
    )

    raw_table = read_raw_table(raw_path)
    samples = sample_columns(raw_table)

    identifier_audit = identifier_pair_audit(
        table=raw_table,
        sample_cols=samples,
        weights=weights,
    )
    (
        gene_column,
        frozen_identifier_column,
        identifier_scheme,
    ) = choose_identifier_pair(identifier_audit)

    identifier_audit["selected_identifier_pair"] = (
        identifier_audit["external_column"].eq(gene_column)
        & identifier_audit["frozen_identifier_column"].eq(
            frozen_identifier_column
        )
        & identifier_audit["identifier_scheme"].eq(
            identifier_scheme
        )
    )
    identifier_audit.to_csv(
        OUTPUT_IDENTIFIER_PAIR_AUDIT,
        index=False,
    )

    weights["canine_gene_symbol"] = weights[
        frozen_identifier_column
    ].map(
        lambda value: normalize_identifier(
            value,
            identifier_scheme,
        )
    )

    external_raw, gene_map = collapse_to_identifiers(
        table=raw_table,
        gene_column=gene_column,
        sample_cols=samples,
        identifier_scheme=identifier_scheme,
    )
    external, transform_diagnostics = choose_transform(
        external_raw
    )

    reference, reference_gene_map = remap_reference_expression(
        reference=reference,
        weights=weights,
        frozen_identifier_column=frozen_identifier_column,
        identifier_scheme=identifier_scheme,
    )

    external = external.loc[
        :,
        ~pd.Index(external.columns).duplicated(keep="first"),
    ]

    reference_z = zscore_columns(reference)
    external_z = zscore_columns(external)

    input_audit = identifier_audit.copy()
    input_audit["n_detected_sample_columns"] = len(samples)
    input_audit["external_matrix_rows"] = raw_table.shape[0]
    input_audit["external_matrix_columns"] = raw_table.shape[1]
    input_audit["selected_external_column"] = gene_column
    input_audit["selected_frozen_identifier_column"] = (
        frozen_identifier_column
    )
    input_audit["selected_identifier_scheme"] = identifier_scheme
    input_audit.to_csv(OUTPUT_INPUT_AUDIT, index=False)

    gene_map = gene_map.merge(
        reference_gene_map,
        on="analysis_gene_id",
        how="left",
    )
    gene_map.to_csv(OUTPUT_GENE_MAP, index=False)
    external.to_csv(OUTPUT_EXPRESSION)

    print("")
    print("Selected identifier mapping:")
    print(f"  External column: {gene_column}")
    print(f"  Frozen identifier column: {frozen_identifier_column}")
    print(f"  Identifier scheme: {identifier_scheme}")
    print(
        identifier_audit.head(10)[
            [
                "external_column",
                "frozen_identifier_column",
                "identifier_scheme",
                "n_unique_overlap",
                "n_primary_modules_with_at_least_3_genes",
                "M34_overlap",
                "M11_overlap",
                "M24_overlap",
                "M40_overlap",
                "selected_identifier_pair",
            ]
        ].to_string(index=False)
    )

    structure_rows = []
    coverage_rows = []
    reliability_rows = []
    loo_tables = []
    score_table = pd.DataFrame(index=external.index)

    for module_index, module in enumerate(PRIMARY_MODULES):
        result, coverage, loo, score = analyze_module(
            module=module,
            reference_z=reference_z,
            external_z=external_z,
            weights=weights,
            seed=RANDOM_SEED + module_index * 100,
        )
        structure_rows.append(result)
        coverage_rows.append(coverage)

        reliability_rows.append(
            {
                "module_label": module,
                "n_common_genes": result["n_common_genes"],
                "split_half_median": result[
                    "split_half_median"
                ],
                "split_half_q05": result["split_half_q05"],
                "split_half_q95": result["split_half_q95"],
                "split_half_valid_repeats": result[
                    "split_half_valid_repeats"
                ],
                "minimum_gene_loo_correlation": (
                    float(
                        loo["correlation_with_full_score"].min()
                    )
                    if not loo.empty
                    else np.nan
                ),
                "median_gene_loo_correlation": (
                    float(
                        loo["correlation_with_full_score"].median()
                    )
                    if not loo.empty
                    else np.nan
                ),
            }
        )

        if not loo.empty:
            loo_tables.append(loo)

        score_table[
            f"{module}__strict_signed_mean_z"
        ] = (
            (score - score.mean()) / score.std()
            if score.notna().sum() > 1
            else score
        )

    structure = pd.DataFrame(structure_rows)
    coverage = pd.DataFrame(coverage_rows)
    reliability = pd.DataFrame(reliability_rows)
    gene_loo = (
        pd.concat(loo_tables, ignore_index=True)
        if loo_tables
        else pd.DataFrame()
    )

    random_controls = random_panel_controls(
        reference_z=reference_z,
        external_z=external_z,
        weights=weights,
        observed=structure,
    )
    classification = classify_results(
        structure=structure,
        random_controls=random_controls,
    )

    score_table.to_csv(OUTPUT_SCORES)
    coverage.to_csv(OUTPUT_COVERAGE, index=False)
    structure.to_csv(OUTPUT_STRUCTURE, index=False)
    reliability.to_csv(OUTPUT_RELIABILITY, index=False)
    gene_loo.to_csv(OUTPUT_GENE_LOO, index=False)
    random_controls.to_csv(OUTPUT_RANDOM, index=False)
    classification.to_csv(
        OUTPUT_CLASSIFICATION,
        index=False,
    )

    create_heatmap(classification)
    write_readme(raw_path, transform_diagnostics)

    input_paths = [
        raw_path,
        REFERENCE_EXPRESSION_FILE,
        STRICT_WEIGHTS_FILE,
    ]
    if FREEZE_FILE.exists():
        input_paths.append(FREEZE_FILE)

    output_paths = [
        OUTPUT_INPUT_AUDIT,
        OUTPUT_IDENTIFIER_PAIR_AUDIT,
        OUTPUT_GENE_MAP,
        OUTPUT_EXPRESSION,
        OUTPUT_SCORES,
        OUTPUT_COVERAGE,
        OUTPUT_STRUCTURE,
        OUTPUT_RELIABILITY,
        OUTPUT_GENE_LOO,
        OUTPUT_RANDOM,
        OUTPUT_CLASSIFICATION,
        OUTPUT_README,
        HEATMAP_PNG,
        HEATMAP_PDF,
    ]

    manifest = {
        "script_version": SCRIPT_VERSION,
        "created_utc": utc_now_iso(),
        "accession": ACCESSION,
        "n_external_samples": int(external.shape[0]),
        "n_external_genes": int(external.shape[1]),
        "selected_gene_identifier_column": gene_column,
        "selected_frozen_identifier_column": frozen_identifier_column,
        "selected_identifier_scheme": identifier_scheme,
        "transform_diagnostics": transform_diagnostics,
        "frozen_definition": freeze,
        "gene_label_permutations": N_GENE_LABEL_PERMUTATIONS,
        "random_panels": N_RANDOM_PANELS,
        "split_half_repeats": N_SPLIT_HALF_REPEATS,
        "outcome_loaded": False,
        "guardrails": [
            "No external outcome was loaded.",
            "Frozen module genes, loadings, signs, and tiers were unchanged.",
            "External PC1 orientation used only the frozen signed score.",
            "Representation preservation is not prognostic validation.",
            "Therapy heterogeneity may contribute to expression-domain shift.",
        ],
        "inputs": {
            str(path): {
                "sha256": sha256_file(path),
                "size_bytes": path.stat().st_size,
            }
            for path in input_paths
        },
        "outputs": {
            str(path): {
                "sha256": sha256_file(path),
                "size_bytes": path.stat().st_size,
            }
            for path in output_paths
            if path.exists()
        },
    }
    OUTPUT_MANIFEST.write_text(
        json.dumps(manifest, indent=2),
        encoding="utf-8",
    )

    print("")
    print("=" * 80)
    print("External matrix preflight")
    print("=" * 80)
    print(
        pd.DataFrame(
            [
                {
                    "raw_rows": raw_table.shape[0],
                    "raw_columns": raw_table.shape[1],
                    "sample_columns": len(samples),
                    "selected_gene_column": gene_column,
                    "selected_frozen_identifier_column": frozen_identifier_column,
                    "identifier_scheme": identifier_scheme,
                    "processed_samples": external.shape[0],
                    "processed_genes": external.shape[1],
                    "transform": transform_diagnostics["transform"],
                }
            ]
        ).to_string(index=False)
    )
    print("")
    print(
        input_audit.head(10).to_string(index=False)
    )

    print("")
    print("=" * 80)
    print("Frozen-program external canine coverage")
    print("=" * 80)
    print(
        coverage[
            [
                "module_label",
                "n_frozen_genes",
                "n_common_genes",
                "coverage_fraction",
            ]
        ].to_string(index=False)
    )

    print("")
    print("=" * 80)
    print("External canine representation preservation")
    print("=" * 80)
    print(
        classification[
            [
                "module_label",
                "n_common_genes",
                "edge_spearman",
                "edge_q_bh_8",
                "loading_spearman",
                "loading_q_bh_8",
                "split_half_median",
                "random_panel_empirical_p",
                "external_canine_representation_class",
            ]
        ].to_string(index=False)
    )

    print("")
    print("=" * 80)
    print("Frozen-score reliability")
    print("=" * 80)
    print(reliability.to_string(index=False))

    print("")
    print("=" * 80)
    print("Interpretation guardrails")
    print("=" * 80)
    print("GSE239948 is an independent canine expression cohort, not an outcome-validation cohort.")
    print("No external survival or treatment-response endpoint was used.")
    print("Strong preservation requires edge, loading, and split-half support.")
    print("Therapy heterogeneity may contribute to external expression differences.")
    print("Human transfer tiers and project-wide multiplicity remain unchanged.")

    print("")
    print("Saved:")
    for path in output_paths + [OUTPUT_MANIFEST]:
        if path.exists():
            print(path)
    print("Done.")
```

### Relevant standalone lines

| Line | Signals | Code |
|---:|---|---|
| 250 | mean | `numeric_fraction[column] = float(numeric.notna().mean())` |
| 328 | weights | `def frozen_identifier_columns(weights: pd.DataFrame) -> list[str]:` |
| 341 | weights | `if column in weights.columns:` |
| 344 | weights | `for column in weights.columns:` |
| 362 | weights | `weights: pd.DataFrame,` |
| 369 | weights | `frozen_columns = frozen_identifier_columns(weights)` |
| 377 | weights | `"No canine identifier column was found in the frozen weights."` |
| 380 | weights | `primary_weights = weights[` |
| 381 | weights | `weights["module_label"].isin(PRIMARY_MODULES)` |
| 554 | weights | `weights: pd.DataFrame,` |
| 558 | weights | `if "canine_gene" not in weights.columns:` |
| 560 | weights | `"Frozen weights are missing the canine_gene column needed "` |
| 564 | weights | `mapping = weights[` |
| 624 | mean | `np.mean(np.isclose(finite, np.round(finite)))` |
| 647 | std | `std = x.std(axis=0).replace(0, np.nan)` |
| 648 | mean | `z = (x - x.mean(axis=0)) / std` |
| 716 | loadings | `loadings: pd.Series,` |
| 718 | loadings, direction | `signs = np.sign(loadings.reindex(genes)).replace(0, 1)` |
| 722 | mean | `).mean(axis=1)` |
| 729 | loadings | `loadings: pd.Series,` |
| 753 | loadings | `loadings,` |
| 758 | loadings | `loadings,` |
| 760 | correlation | `correlation = first_score.corr(second_score)` |
| 789 | loadings | `loadings: pd.Series,` |
| 794 | loadings | `loadings,` |
| 806 | loadings | `loadings,` |
| 814 | correlation | `full_score.corr(score)` |
| 826 | weights | `weights: pd.DataFrame,` |
| 829 | weights | `part = weights[` |
| 830 | weights | `weights["module_label"].astype(str).eq(module)` |
| 848 | coverage | `coverage = {` |
| 882 | coverage | `return result, coverage, pd.DataFrame(), pd.Series(` |
| 931 | correlation | `orientation_correlation = pc_score.corr(signed_score)` |
| 984 | coverage | `return result, coverage, loo, signed_score` |
| 998 | rank | `reference_variance.rank(pct=True)` |
| 999 | rank | `+ external_variance.rank(pct=True)` |
| 1003 | rank | `combined_rank.rank(method="average"),` |
| 1013 | weights | `weights: pd.DataFrame,` |
| 1022 | weights | `weights.loc[` |
| 1023 | weights | `weights["module_label"].isin(PRIMARY_MODULES),` |
| 1035 | weights | `part = weights[` |
| 1036 | weights | `weights["module_label"].astype(str).eq(module)` |
| 1176 | loadings | `).assign(test="loading"),` |
| 1190 | loadings | `direct_p[direct_p["test"].eq("loading")][` |
| 1264 | loadings | `"Frozen-loading\nconcordance",` |
| 1315 | zscore | `No outcome model is fitted because GSE239948 does not provide a standardized` |
| 1320 | correlation | `1. Spearman concordance of within-module gene-gene correlation edges between` |
| 1322 | loadings | `2. Concordance between frozen DOG2 risk-oriented loadings and GSE239948 PC1` |
| 1323 | loadings, direction | `loadings. The external PC1 sign is oriented using the frozen signed score,` |
| 1325 | direction | `3. Repeated non-overlapping split-half reliability of the frozen signed score.` |
| 1332 | loadings | `loading-concordance tests.` |
| 1336 | loadings | `- Frozen genes, loadings, signs, and tiers are unchanged.` |
| 1355 | loadings | `print("  Test edge preservation, loading concordance, and score reliability.")` |
| 1365 | weights | `weights = pd.read_csv(STRICT_WEIGHTS_FILE)` |
| 1378 | weights | `weights=weights,` |
| 1400 | weights | `weights["canine_gene_symbol"] = weights[` |
| 1421 | weights | `weights=weights,` |
| 1482 | coverage | `result, coverage, loo, score = analyze_module(` |
| 1486 | weights | `weights=weights,` |
| 1490 | coverage | `coverage_rows.append(coverage)` |
| 1527 | mean, std | `(score - score.mean()) / score.std()` |
| 1533 | coverage | `coverage = pd.DataFrame(coverage_rows)` |
| 1544 | weights | `weights=weights,` |
| 1553 | coverage | `coverage.to_csv(OUTPUT_COVERAGE, index=False)` |
| 1608 | loadings | `"Frozen module genes, loadings, signs, and tiers were unchanged.",` |
| 1609 | direction | `"External PC1 orientation used only the frozen signed score.",` |
| 1662 | coverage | `print("Frozen-program external canine coverage")` |
| 1665 | coverage | `coverage[` |
| 1707 | loadings | `print("Strong preservation requires edge, loading, and split-half support.")` |

### Referenced result-table schemas

#### `results/tables/GSE238110_frozen_transfer_gene_weights_strict.csv`

- Rows: 321
- SHA-256: `4f065aa4c4edf117a0c74015840d2b4b2347929f172cd517e1818ba0f6163b91`
- Columns: `module_label`, `canine_gene`, `canine_gene_symbol`, `human_gene_symbol`, `ortholog_qc_status`, `raw_pca_loading`, `risk_oriented_loading`, `absolute_risk_loading`, `is_strict_mapping`, `is_broad_mapping`, `human_symbol_duplicate_count`, `normalized_abs_sum_weight`

#### `results/tables/GSE239948_external_frozen_program_coverage.csv`

- Rows: 4
- SHA-256: `8a3f4d87bb9ebf1bf0d0a3c1c69f38e95d020abed5968a350f5c1db21578cf80`
- Columns: `module_label`, `n_frozen_genes`, `n_common_genes`, `coverage_fraction`, `common_genes`, `missing_external_genes`

#### `results/tables/GSE239948_external_gene_mapping.csv`

- Rows: 15324
- SHA-256: `56aff5474085f70f9ffb7b54024ba3ccbe40517bef5b48b09fc6e6b39fd034d3`
- Columns: `Name`, `analysis_gene_id`, `row_variance`, `reference_expression_column`

#### `results/tables/GSE239948_external_identifier_pair_audit.csv`

- Rows: 4
- SHA-256: `e8db70407048465652527a0f3e5d66637fbe141d811ce083b875b142396b2b88`
- Columns: `external_column`, `frozen_identifier_column`, `identifier_scheme`, `n_external_unique_ids`, `n_frozen_unique_ids`, `n_unique_overlap`, `n_primary_modules_with_at_least_3_genes`, `M34_overlap`, `M11_overlap`, `M24_overlap`, `M40_overlap`, `example_overlap`, `frozen_column_preference`, `scheme_preference`, `selected_identifier_pair`

#### `results/tables/GSE239948_external_input_audit.csv`

- Rows: 4
- SHA-256: `fc70f9ef0219df153f56d224439bda07b4838935bca2c407f2e92b95a144923f`
- Columns: `external_column`, `frozen_identifier_column`, `identifier_scheme`, `n_external_unique_ids`, `n_frozen_unique_ids`, `n_unique_overlap`, `n_primary_modules_with_at_least_3_genes`, `M34_overlap`, `M11_overlap`, `M24_overlap`, `M40_overlap`, `example_overlap`, `frozen_column_preference`, `scheme_preference`, `selected_identifier_pair`, `n_detected_sample_columns`, `external_matrix_rows`, `external_matrix_columns`, `selected_external_column`, `selected_frozen_identifier_column`, `selected_identifier_scheme`

#### `results/tables/GSE239948_external_module_gene_leave_one_out.csv`

- Rows: 272
- SHA-256: `af1c55371afd999ae4d05a517d9bc3cd90a0ee9dad9437b7e94d1bdaf895a37d`
- Columns: `module_label`, `left_out_gene`, `n_genes_remaining`, `correlation_with_full_score`

#### `results/tables/GSE239948_external_module_score_reliability.csv`

- Rows: 4
- SHA-256: `3a9412a26547cbb926eaccf38451373d82d6282884da5f6c0724be31825e183a`
- Columns: `module_label`, `n_common_genes`, `split_half_median`, `split_half_q05`, `split_half_q95`, `split_half_valid_repeats`, `minimum_gene_loo_correlation`, `median_gene_loo_correlation`

#### `results/tables/GSE239948_external_module_structure_preservation.csv`

- Rows: 4
- SHA-256: `3a23d206704be171b611cae997e9ebeebbbaf943b8d6a6b09e821f4ed1c9e305`
- Columns: `module_label`, `n_common_genes`, `edge_spearman`, `edge_permutation_p`, `loading_spearman`, `loading_permutation_p`, `external_pc1_variance_explained`, `pc1_orientation_correlation_with_frozen_score`, `split_half_median`, `split_half_q05`, `split_half_q95`, `split_half_valid_repeats`, `estimable`, `nonestimable_reason`

#### `results/tables/GSE239948_external_random_panel_controls.csv`

- Rows: 4
- SHA-256: `472b0676317ed64e2e56d02772e918375b2b96013ed4bb7219f5dac18c979195`
- Columns: `module_label`, `n_module_genes`, `n_random_panels`, `observed_edge_spearman`, `random_edge_median`, `random_edge_q05`, `random_edge_q95`, `random_panel_empirical_p`

#### `results/tables/GSE239948_external_representation_classification.csv`

- Rows: 4
- SHA-256: `b4c966b753d7fd8d4132445ec71a1365be3774aee8f2d745bff4cbed62aa812e`
- Columns: `module_label`, `n_common_genes`, `edge_spearman`, `edge_permutation_p`, `loading_spearman`, `loading_permutation_p`, `external_pc1_variance_explained`, `pc1_orientation_correlation_with_frozen_score`, `split_half_median`, `split_half_q05`, `split_half_q95`, `split_half_valid_repeats`, `estimable`, `nonestimable_reason`, `random_panel_empirical_p`, `edge_q_bh_8`, `loading_q_bh_8`, `external_canine_representation_class`

---

## `scripts/49_gse239948_blind_de_novo_rediscovery_FIXED_V2.py`

- SHA-256: `40540bd7f10f7eb1c8f9bbf4822c54c0ee74d6e12e657b1277c7efae0f4d2247`

### Relevant functions

#### `select_discovery_universe`

- Lines: 122-137
- Signals: `rank`

```python
def select_discovery_universe(expression: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    variances = expression.var(axis=0, ddof=1).sort_values(ascending=False)
    n_keep = min(DISCOVERY_TOP_VARIABLE_GENES, variances.shape[0])
    selected = variances.iloc[:n_keep].index.tolist()
    x = expression[selected].copy()

    ranks = variances.rank(method="average", pct=True)
    audit = pd.DataFrame(
        {
            "gene_symbol": variances.index,
            "variance": variances.values,
            "variance_percentile": ranks.reindex(variances.index).values,
            "selected_for_blind_discovery": variances.index.isin(selected),
        }
    )
    return x, audit
```

#### `rank_transform`

- Lines: 147-155
- Signals: `mean, std, rank`

```python
def rank_transform(expression: pd.DataFrame) -> np.ndarray:
    ranked = expression.rank(axis=0, method="average")
    values = ranked.to_numpy(dtype=float)
    values -= values.mean(axis=0, keepdims=True)
    std = values.std(axis=0, ddof=1, keepdims=True)
    valid = np.isfinite(std) & (std > 0)
    std[~valid] = 1.0
    values /= std
    return values
```

#### `spearman_gene_correlation`

- Lines: 166-181
- Signals: `matrix_product, correlation`

```python
def spearman_gene_correlation(expression: pd.DataFrame) -> np.ndarray:
    ranked = rank_transform(expression)
    n_samples = ranked.shape[0]
    if n_samples < 2:
        return np.eye(ranked.shape[1], dtype=float)

    corr = (ranked.T @ ranked) / float(n_samples - 1)
    corr = np.asarray(corr, dtype=float)
    corr = np.clip(corr, -1.0, 1.0)
    np.fill_diagonal(corr, 1.0)

    if not np.isfinite(corr).all():
        raise RuntimeError(
            "Non-finite Spearman correlation remained after variance filtering."
        )
    return corr
```

#### `cluster_expression`

- Lines: 184-219
- Signals: `correlation`

```python
def cluster_expression(expression: pd.DataFrame) -> tuple[pd.Series, float]:
    all_genes = pd.Index(expression.columns)
    labels_all = pd.Series(0, index=all_genes, dtype=int)

    clusterable = clusterable_gene_mask(expression)
    x = expression.loc[:, clusterable].copy()
    if x.shape[1] < MIN_DISCOVERED_MODULE_SIZE:
        return labels_all, np.nan

    corr = spearman_gene_correlation(x)
    distance = 1.0 - corr
    distance = np.clip(distance, 0.0, 2.0)
    np.fill_diagonal(distance, 0.0)

    condensed = squareform(distance, checks=False)
    if not np.isfinite(condensed).all():
        raise RuntimeError(
            "Non-finite clustering distances remained after removing "
            "zero-variance genes."
        )
    tree = linkage(condensed, method="average", optimal_ordering=False)

    r_cut = correlation_threshold(x.shape[0])
    distance_cut = 1.0 - r_cut
    raw_labels = fcluster(tree, t=distance_cut, criterion="distance")
    labels = pd.Series(raw_labels, index=x.columns, dtype=int)

    counts = labels.value_counts()
    valid_modules = counts[counts >= MIN_DISCOVERED_MODULE_SIZE].index
    labels = labels.where(labels.isin(valid_modules), 0)

    nonzero = sorted(v for v in labels.unique() if v != 0)
    remap = {old: new for new, old in enumerate(nonzero, start=1)}
    labels = labels.map(lambda value: remap.get(value, 0)).astype(int)
    labels_all.loc[labels.index] = labels
    return labels_all, r_cut
```

#### `summarize_modules`

- Lines: 222-238
- Signals: `mean, correlation`

```python
def summarize_modules(expression: pd.DataFrame, labels: pd.Series) -> pd.DataFrame:
    rows = []
    for module_id in sorted(v for v in labels.unique() if v != 0):
        genes = labels.index[labels.eq(module_id)].tolist()
        part = expression[genes]
        corr = spearman_gene_correlation(part)
        upper = corr[np.triu_indices(len(genes), k=1)]
        rows.append(
            {
                "discovered_module_id": int(module_id),
                "n_genes": len(genes),
                "mean_pairwise_spearman": float(np.nanmean(upper)) if upper.size else np.nan,
                "median_pairwise_spearman": float(np.nanmedian(upper)) if upper.size else np.nan,
                "genes": ";".join(genes),
            }
        )
    return pd.DataFrame(rows)
```

#### `module_sets_from_labels`

- Lines: 241-245
- Signals: ``

```python
def module_sets_from_labels(labels: pd.Series) -> dict[int, set[str]]:
    result = {}
    for module_id in sorted(v for v in labels.unique() if v != 0):
        result[int(module_id)] = set(labels.index[labels.eq(module_id)])
    return result
```

#### `best_set_match`

- Lines: 258-287
- Signals: `coverage`

```python
def best_set_match(query: set[str], candidates: dict[int, set[str]]) -> dict:
    best = {
        "module_id": np.nan,
        "overlap": 0,
        "jaccard": 0.0,
        "f1": 0.0,
        "query_recall": 0.0,
        "module_precision": 0.0,
        "module_size": 0,
    }
    if not query or not candidates:
        return best

    for module_id, genes in candidates.items():
        overlap = len(query & genes)
        score = f1_overlap(query, genes)
        if (
            score > best["f1"]
            or (np.isclose(score, best["f1"]) and overlap > best["overlap"])
        ):
            best = {
                "module_id": int(module_id),
                "overlap": int(overlap),
                "jaccard": float(jaccard(query, genes)),
                "f1": float(score),
                "query_recall": float(overlap / len(query)) if query else np.nan,
                "module_precision": float(overlap / len(genes)) if genes else np.nan,
                "module_size": int(len(genes)),
            }
    return best
```

#### `subsample_stability`

- Lines: 290-347
- Signals: `mean, coverage`

```python
def subsample_stability(
    expression: pd.DataFrame,
    full_labels: pd.Series,
    rng: np.random.Generator,
) -> pd.DataFrame:
    full_modules = module_sets_from_labels(full_labels)
    records = []
    n_samples = expression.shape[0]
    subset_n = max(10, int(round(n_samples * SUBSAMPLE_FRACTION)))

    for repeat in range(1, SUBSAMPLE_REPEATS + 1):
        chosen = rng.choice(n_samples, size=subset_n, replace=False)
        subset = expression.iloc[np.sort(chosen)]
        clusterable = clusterable_gene_mask(subset)
        n_clusterable = int(clusterable.sum())
        n_dropped = int((~clusterable).sum())

        labels, r_cut = cluster_expression(subset)
        repeat_modules = module_sets_from_labels(labels)

        for module_id, genes in full_modules.items():
            match = best_set_match(genes, repeat_modules)
            records.append(
                {
                    "repeat": repeat,
                    "full_discovered_module_id": module_id,
                    "subsample_n": subset_n,
                    "subsample_r_cut": r_cut,
                    "n_clusterable_genes": n_clusterable,
                    "n_zero_or_near_zero_variance_genes_dropped": n_dropped,
                    "best_subsample_module_id": match["module_id"],
                    "best_jaccard": match["jaccard"],
                    "best_f1": match["f1"],
                    "overlap_genes": match["overlap"],
                }
            )

    raw = pd.DataFrame(records)
    if raw.empty:
        return raw

    summary = (
        raw.groupby("full_discovered_module_id", as_index=False)
        .agg(
            n_repeats=("repeat", "count"),
            median_clusterable_genes=("n_clusterable_genes", "median"),
            max_zero_or_near_zero_variance_genes_dropped=(
                "n_zero_or_near_zero_variance_genes_dropped",
                "max",
            ),
            median_best_jaccard=("best_jaccard", "median"),
            q05_best_jaccard=("best_jaccard", lambda x: float(np.quantile(x, 0.05))),
            median_best_f1=("best_f1", "median"),
            fraction_jaccard_ge_0_50=("best_jaccard", lambda x: float(np.mean(np.asarray(x) >= 0.50))),
            fraction_f1_ge_0_50=("best_f1", lambda x: float(np.mean(np.asarray(x) >= 0.50))),
        )
    )
    return summary
```

#### `find_canine_gene_column`

- Lines: 350-354
- Signals: `weights`

```python
def find_canine_gene_column(weights: pd.DataFrame) -> str:
    for column in ["canine_gene_symbol", "canine_gene", "gene"]:
        if column in weights.columns:
            return column
    raise ValueError("No canine gene-symbol column found in frozen weights.")
```

#### `frozen_gene_sets`

- Lines: 357-367
- Signals: `weights`

```python
def frozen_gene_sets(weights: pd.DataFrame) -> dict[str, set[str]]:
    gene_col = find_canine_gene_column(weights)
    output = {}
    for module in PRIMARY_MODULES:
        part = weights[weights["module_label"].astype(str).eq(module)]
        output[module] = {
            clean_gene_symbol(value)
            for value in part[gene_col]
            if clean_gene_symbol(value)
        }
    return output
```

#### `make_variance_bins`

- Lines: 389-395
- Signals: `rank`

```python
def make_variance_bins(universe_audit: pd.DataFrame) -> pd.Series:
    part = universe_audit[universe_audit["selected_for_blind_discovery"]].copy()
    values = part.set_index("gene_symbol")["variance"]
    ranks = values.rank(method="first", pct=True)
    bins = np.minimum((ranks * VARIANCE_BINS).astype(int), VARIANCE_BINS - 1)
    bins = bins.clip(lower=0)
    return bins.astype(int)
```

#### `rediscovery_analysis`

- Lines: 428-548
- Signals: `mean, weights, coverage`

```python
def rediscovery_analysis(
    discovery_genes: set[str],
    discovered_modules: dict[int, set[str]],
    stability: pd.DataFrame,
    weights: pd.DataFrame,
    universe_audit: pd.DataFrame,
    rng: np.random.Generator,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    frozen = frozen_gene_sets(weights)
    all_primary_frozen = set().union(*frozen.values())
    variance_bins = make_variance_bins(universe_audit)

    stability_index = (
        stability.set_index("full_discovered_module_id")
        if not stability.empty
        else pd.DataFrame()
    )

    result_rows = []
    null_rows = []

    for module in PRIMARY_MODULES:
        frozen_all = frozen[module]
        frozen_in_universe = frozen_all & discovery_genes
        observed = best_set_match(frozen_in_universe, discovered_modules)

        null_f1 = []
        null_jaccard = []
        for iteration in range(1, RANDOM_PANELS + 1):
            panel = matched_random_panel(
                frozen_in_universe,
                variance_bins,
                forbidden=all_primary_frozen,
                rng=rng,
            )
            match = best_set_match(panel, discovered_modules)
            null_f1.append(match["f1"])
            null_jaccard.append(match["jaccard"])
            null_rows.append(
                {
                    "module_label": module,
                    "iteration": iteration,
                    "random_panel_size": len(panel),
                    "maximum_match_f1": match["f1"],
                    "maximum_match_jaccard": match["jaccard"],
                }
            )

        null_f1_arr = np.asarray(null_f1, dtype=float)
        null_j_arr = np.asarray(null_jaccard, dtype=float)
        empirical_p = (
            1.0 + np.sum(null_f1_arr >= observed["f1"])
        ) / (1.0 + len(null_f1_arr))

        best_id = observed["module_id"]
        stability_median = np.nan
        stability_q05 = np.nan
        stability_fraction = np.nan
        if (
            np.isfinite(best_id)
            and not stability.empty
            and int(best_id) in stability_index.index
        ):
            row = stability_index.loc[int(best_id)]
            stability_median = float(row["median_best_jaccard"])
            stability_q05 = float(row["q05_best_jaccard"])
            stability_fraction = float(row["fraction_jaccard_ge_0_50"])

        result_rows.append(
            {
                "module_label": module,
                "n_frozen_canine_genes": len(frozen_all),
                "n_frozen_genes_in_blind_discovery_universe": len(frozen_in_universe),
                "discovery_universe_coverage_fraction": (
                    len(frozen_in_universe) / len(frozen_all) if frozen_all else np.nan
                ),
                "best_discovered_module_id": best_id,
                "best_discovered_module_size": observed["module_size"],
                "overlap_genes": observed["overlap"],
                "frozen_gene_recall_within_discovery_universe": observed["query_recall"],
                "discovered_module_precision": observed["module_precision"],
                "best_match_jaccard": observed["jaccard"],
                "best_match_f1": observed["f1"],
                "best_module_subsample_stability_median_jaccard": stability_median,
                "best_module_subsample_stability_q05_jaccard": stability_q05,
                "best_module_fraction_subsamples_jaccard_ge_0_50": stability_fraction,
                "random_max_f1_mean": float(np.nanmean(null_f1_arr)),
                "random_max_f1_q95": float(np.nanquantile(null_f1_arr, 0.95)),
                "random_max_jaccard_q95": float(np.nanquantile(null_j_arr, 0.95)),
                "empirical_max_match_p": float(empirical_p),
            }
        )

    results = pd.DataFrame(result_rows)
    results["empirical_max_match_q_bh_4"] = bh_adjust(results["empirical_max_match_p"])

    classes = []
    for _, row in results.iterrows():
        n_in = int(row["n_frozen_genes_in_blind_discovery_universe"])
        q = row["empirical_max_match_q_bh_4"]
        recall = row["frozen_gene_recall_within_discovery_universe"]
        stability_median = row["best_module_subsample_stability_median_jaccard"]

        if n_in < 3:
            label = "insufficient_blind_discovery_coverage"
        elif (
            np.isfinite(q)
            and q < 0.05
            and recall >= 0.30
            and np.isfinite(stability_median)
            and stability_median >= 0.50
        ):
            label = "strong_blind_independent_rediscovery"
        elif np.isfinite(q) and q < 0.10 and recall >= 0.20:
            label = "partial_blind_independent_rediscovery"
        else:
            label = "no_clear_blind_independent_rediscovery"
        classes.append(label)

    results["blind_rediscovery_class"] = classes
    return results, pd.DataFrame(null_rows)
```

#### `create_readme`

- Lines: 551-588
- Signals: `correlation, coverage`

```python
def create_readme(r_cut: float) -> None:
    text = f"""GSE239948 blind de novo rediscovery audit
=========================================

Script version
--------------
{SCRIPT_VERSION}

Purpose
-------
Test whether de novo co-expression modules formed in GSE239948 recover frozen canine programs without using frozen membership during module discovery.

Discovery procedure
-------------------
- Only GSE239948 expression is used during feature selection and clustering.
- The top {DISCOVERY_TOP_VARIABLE_GENES} genes by GSE239948 variance define the discovery universe.
- Gene-gene Spearman correlation is clustered by average linkage.
- Genes with zero or near-zero variance inside a sample subsample are excluded from that subsample clustering and retained as unassigned genes.
- The full-cohort correlation cut corresponds to a one-sided nominal positive-correlation P value of {CORRELATION_GRAPH_NOMINAL_P}; the observed full-cohort r cut was {r_cut:.4f}.
- Clusters smaller than {MIN_DISCOVERED_MODULE_SIZE} genes are not treated as discovered modules.
- Stability is evaluated across {SUBSAMPLE_REPEATS} outcome-blind {int(SUBSAMPLE_FRACTION * 100)}% sample subsamples.

Frozen-program matching
-----------------------
Frozen M34, M11, M24, and M40 gene sets are loaded only after de novo modules have been constructed.
Each frozen program is compared with every discovered module, and the best F1 overlap is retained.
Expression-variance-matched random panels also take their maximum match across all discovered modules, so the empirical null includes the best-of-many-module search.
BH correction is applied across the four frozen primary programs.

Interpretation guardrails
-------------------------
- This is a representation rediscovery analysis, not an outcome validation analysis.
- Variable-gene filtering is outcome-blind but may reduce coverage of frozen programs.
- A significant best-match result supports independent rediscovery of a related co-expression structure; it does not imply identical module boundaries.
- The cohort has 43 samples, so module boundaries and small-module results require caution.
- No human or canine outcome is loaded or used.
"""
    OUTPUT_README.write_text(text, encoding="utf-8")
```

#### `main`

- Lines: 633-791
- Signals: `weights, coverage`

```python
def main() -> None:
    print("=" * 80)
    print("GSE239948 blind de novo consensus rediscovery of canine programs")
    print("=" * 80)
    print(f"Script version: {SCRIPT_VERSION}")
    print(f"Project root: {PROJECT_ROOT}")
    print("")
    print("Design:")
    print("  Build de novo modules using GSE239948 expression only.")
    print("  Do not load frozen module membership until discovery is complete.")
    print("  Audit module stability across outcome-blind sample subsamples.")
    print("  Test best frozen-program matches against variance-matched random panels.")
    print("  Correct the null for searching across all discovered modules.")
    print("")

    raw = read_required_csv(EXPRESSION_FILE, index_col=0)
    expression = prepare_expression(raw)
    discovery_expression, universe_audit = select_discovery_universe(expression)

    print("Blind discovery data:")
    print(f"  Samples: {discovery_expression.shape[0]}")
    print(f"  Expression genes available: {expression.shape[1]}")
    print(f"  Genes in blind discovery universe: {discovery_expression.shape[1]}")

    labels, r_cut = cluster_expression(discovery_expression)
    discovered_modules = module_sets_from_labels(labels)
    module_summary = summarize_modules(discovery_expression, labels)

    membership = pd.DataFrame(
        {
            "gene_symbol": labels.index,
            "discovered_module_id": labels.values,
        }
    )
    membership = membership.merge(
        universe_audit[["gene_symbol", "variance", "variance_percentile"]],
        on="gene_symbol",
        how="left",
    )

    print("")
    print("=" * 80)
    print("Blind de novo module discovery")
    print("=" * 80)
    print(f"Positive-correlation cut: r >= {r_cut:.4f}")
    print(f"Discovered modules >= {MIN_DISCOVERED_MODULE_SIZE} genes: {len(discovered_modules)}")
    if not module_summary.empty:
        display = module_summary.sort_values("n_genes", ascending=False).head(30)
        print(display[[
            "discovered_module_id",
            "n_genes",
            "mean_pairwise_spearman",
            "median_pairwise_spearman",
        ]].to_string(index=False))

    rng = np.random.default_rng(RANDOM_SEED)
    stability = subsample_stability(discovery_expression, labels, rng)

    print("")
    print("=" * 80)
    print("Subsample module stability")
    print("=" * 80)
    if stability.empty:
        print("No stable discovered modules were available for subsample analysis.")
    else:
        display = stability.sort_values(
            ["median_best_jaccard", "median_best_f1"],
            ascending=False,
        ).head(30)
        print(display.to_string(index=False))

    # Frozen membership is intentionally loaded only after blind discovery is complete.
    weights = read_required_csv(FROZEN_WEIGHTS_FILE)
    if SCRIPT47_LOCK_FILE.exists():
        read_required_csv(SCRIPT47_LOCK_FILE)
    if SCRIPT47_MANIFEST_FILE.exists():
        print(f"Loaded: {SCRIPT47_MANIFEST_FILE}")
        manifest47 = json.loads(SCRIPT47_MANIFEST_FILE.read_text(encoding="utf-8"))
        if manifest47.get("script_version") != "47-lock-gse239948-independent-canine-evidence-v2":
            raise RuntimeError(
                "Script 47 lock is not from the expected v2 independence audit."
            )

    rediscovery, random_controls = rediscovery_analysis(
        discovery_genes=set(discovery_expression.columns),
        discovered_modules=discovered_modules,
        stability=stability,
        weights=weights,
        universe_audit=universe_audit,
        rng=rng,
    )

    print("")
    print("=" * 80)
    print("Blind frozen-program rediscovery")
    print("=" * 80)
    columns = [
        "module_label",
        "n_frozen_canine_genes",
        "n_frozen_genes_in_blind_discovery_universe",
        "discovery_universe_coverage_fraction",
        "best_discovered_module_id",
        "best_discovered_module_size",
        "overlap_genes",
        "frozen_gene_recall_within_discovery_universe",
        "discovered_module_precision",
        "best_match_jaccard",
        "best_match_f1",
        "best_module_subsample_stability_median_jaccard",
        "random_max_f1_q95",
        "empirical_max_match_p",
        "empirical_max_match_q_bh_4",
        "blind_rediscovery_class",
    ]
    print(rediscovery[columns].to_string(index=False))

    print("")
    print("=" * 80)
    print("Interpretation guardrails")
    print("=" * 80)
    print("Frozen genes are not used to construct the de novo modules.")
    print("The random-panel null repeats the best-of-many discovered-module search.")
    print("Variable-gene filtering is outcome-blind and natural frozen-gene coverage is reported.")
    print("A recovered module need not have identical boundaries to the frozen canine program.")
    print("No outcome data are loaded or tested.")

    universe_audit.to_csv(OUTPUT_GENE_UNIVERSE, index=False)
    membership.to_csv(OUTPUT_MODULE_MEMBERSHIP, index=False)
    module_summary.to_csv(OUTPUT_MODULE_SUMMARY, index=False)
    stability.to_csv(OUTPUT_SUBSAMPLE_STABILITY, index=False)
    rediscovery.to_csv(OUTPUT_REDISCOVERY, index=False)
    random_controls.to_csv(OUTPUT_RANDOM_CONTROLS, index=False)
    create_readme(r_cut)

    output_paths = [
        OUTPUT_GENE_UNIVERSE,
        OUTPUT_MODULE_MEMBERSHIP,
        OUTPUT_MODULE_SUMMARY,
        OUTPUT_SUBSAMPLE_STABILITY,
        OUTPUT_REDISCOVERY,
        OUTPUT_RANDOM_CONTROLS,
        OUTPUT_README,
    ]
    create_manifest(
        input_paths=[
            EXPRESSION_FILE,
            FROZEN_WEIGHTS_FILE,
            SCRIPT47_LOCK_FILE,
            SCRIPT47_MANIFEST_FILE,
        ],
        output_paths=output_paths,
        r_cut=r_cut,
    )

    print("")
    print("Saved:")
    for path in output_paths + [OUTPUT_MANIFEST]:
        print(path)
    print("Done.")
```

### Relevant standalone lines

| Line | Signals | Code |
|---:|---|---|
| 128 | rank | `ranks = variances.rank(method="average", pct=True)` |
| 148 | rank | `ranked = expression.rank(axis=0, method="average")` |
| 150 | mean | `values -= values.mean(axis=0, keepdims=True)` |
| 151 | std | `std = values.std(axis=0, ddof=1, keepdims=True)` |
| 172 | matrix_product, correlation | `corr = (ranked.T @ ranked) / float(n_samples - 1)` |
| 173 | correlation | `corr = np.asarray(corr, dtype=float)` |
| 174 | correlation | `corr = np.clip(corr, -1.0, 1.0)` |
| 175 | correlation | `np.fill_diagonal(corr, 1.0)` |
| 177 | correlation | `if not np.isfinite(corr).all():` |
| 179 | correlation | `"Non-finite Spearman correlation remained after variance filtering."` |
| 181 | correlation | `return corr` |
| 193 | correlation | `corr = spearman_gene_correlation(x)` |
| 194 | correlation | `distance = 1.0 - corr` |
| 227 | correlation | `corr = spearman_gene_correlation(part)` |
| 228 | correlation | `upper = corr[np.triu_indices(len(genes), k=1)]` |
| 233 | mean | `"mean_pairwise_spearman": float(np.nanmean(upper)) if upper.size else np.nan,` |
| 261 | coverage | `"overlap": 0,` |
| 272 | coverage | `overlap = len(query & genes)` |
| 276 | coverage | `or (np.isclose(score, best["f1"]) and overlap > best["overlap"])` |
| 280 | coverage | `"overlap": int(overlap),` |
| 283 | coverage | `"query_recall": float(overlap / len(query)) if query else np.nan,` |
| 284 | coverage | `"module_precision": float(overlap / len(genes)) if genes else np.nan,` |
| 323 | coverage | `"overlap_genes": match["overlap"],` |
| 343 | mean | `fraction_jaccard_ge_0_50=("best_jaccard", lambda x: float(np.mean(np.asarray(x) >= 0.50))),` |
| 344 | mean | `fraction_f1_ge_0_50=("best_f1", lambda x: float(np.mean(np.asarray(x) >= 0.50))),` |
| 350 | weights | `def find_canine_gene_column(weights: pd.DataFrame) -> str:` |
| 352 | weights | `if column in weights.columns:` |
| 354 | weights | `raise ValueError("No canine gene-symbol column found in frozen weights.")` |
| 357 | weights | `def frozen_gene_sets(weights: pd.DataFrame) -> dict[str, set[str]]:` |
| 358 | weights | `gene_col = find_canine_gene_column(weights)` |
| 361 | weights | `part = weights[weights["module_label"].astype(str).eq(module)]` |
| 392 | rank | `ranks = values.rank(method="first", pct=True)` |
| 432 | weights | `weights: pd.DataFrame,` |
| 436 | weights | `frozen = frozen_gene_sets(weights)` |
| 506 | coverage | `"overlap_genes": observed["overlap"],` |
| 514 | mean | `"random_max_f1_mean": float(np.nanmean(null_f1_arr)),` |
| 567 | correlation | `- Gene-gene Spearman correlation is clustered by average linkage.` |
| 576 | coverage | `Each frozen program is compared with every discovered module, and the best F1 overlap is retained.` |
| 583 | coverage | `- Variable-gene filtering is outcome-blind but may reduce coverage of frozen programs.` |
| 705 | weights | `weights = read_required_csv(FROZEN_WEIGHTS_FILE)` |
| 720 | weights | `weights=weights,` |
| 755 | coverage | `print("Variable-gene filtering is outcome-blind and natural frozen-gene coverage is reported.")` |

### Referenced result-table schemas

#### `results/tables/GSE238110_frozen_transfer_gene_weights_strict.csv`

- Rows: 321
- SHA-256: `4f065aa4c4edf117a0c74015840d2b4b2347929f172cd517e1818ba0f6163b91`
- Columns: `module_label`, `canine_gene`, `canine_gene_symbol`, `human_gene_symbol`, `ortholog_qc_status`, `raw_pca_loading`, `risk_oriented_loading`, `absolute_risk_loading`, `is_strict_mapping`, `is_broad_mapping`, `human_symbol_duplicate_count`, `normalized_abs_sum_weight`

#### `results/tables/GSE239948_blind_discovered_module_membership.csv`

- Rows: 3000
- SHA-256: `6c9287c871f52643166227c7457b8bf201456e7f82889137d3da5c549a5542cd`
- Columns: `gene_symbol`, `discovered_module_id`, `variance`, `variance_percentile`

#### `results/tables/GSE239948_blind_discovered_module_subsample_stability.csv`

- Rows: 134
- SHA-256: `7f9a147427c9ffa6266cf54d2d433376feb28d768a7e6c9f8971a403692734e3`
- Columns: `full_discovered_module_id`, `n_repeats`, `median_clusterable_genes`, `max_zero_or_near_zero_variance_genes_dropped`, `median_best_jaccard`, `q05_best_jaccard`, `median_best_f1`, `fraction_jaccard_ge_0_50`, `fraction_f1_ge_0_50`

#### `results/tables/GSE239948_blind_discovered_module_summary.csv`

- Rows: 134
- SHA-256: `1260a33a86200176b78312270f8ee19e404897c27efb175e240883ecbc681bb7`
- Columns: `discovered_module_id`, `n_genes`, `mean_pairwise_spearman`, `median_pairwise_spearman`, `genes`

#### `results/tables/GSE239948_blind_discovery_gene_universe.csv`

- Rows: 14899
- SHA-256: `3820ab6e2cd6e7800b333c616da30c9c07ded02791b2984da7595d5c41e5cef8`
- Columns: `gene_symbol`, `variance`, `variance_percentile`, `selected_for_blind_discovery`

#### `results/tables/GSE239948_blind_frozen_program_random_controls.csv`

- Rows: 8000
- SHA-256: `aa25d096e9d1fa995055c8152ca476447055f0b4103907fd274a899c43805299`
- Columns: `module_label`, `iteration`, `random_panel_size`, `maximum_match_f1`, `maximum_match_jaccard`

#### `results/tables/GSE239948_blind_frozen_program_rediscovery.csv`

- Rows: 4
- SHA-256: `473895e65c3d9ec8367036e38e5be29209389f522b237afc5faa80c0aeca9446`
- Columns: `module_label`, `n_frozen_canine_genes`, `n_frozen_genes_in_blind_discovery_universe`, `discovery_universe_coverage_fraction`, `best_discovered_module_id`, `best_discovered_module_size`, `overlap_genes`, `frozen_gene_recall_within_discovery_universe`, `discovered_module_precision`, `best_match_jaccard`, `best_match_f1`, `best_module_subsample_stability_median_jaccard`, `best_module_subsample_stability_q05_jaccard`, `best_module_fraction_subsamples_jaccard_ge_0_50`, `random_max_f1_mean`, `random_max_f1_q95`, `random_max_jaccard_q95`, `empirical_max_match_p`, `empirical_max_match_q_bh_4`, `blind_rediscovery_class`

#### `results/tables/paper4_locked_independent_canine_representation.csv`

- Rows: 4
- SHA-256: `642ee4b0a1e5a04a5cdefb3d03d5329bde5ce920f77939473f151e51b3241d89`
- Columns: `module_label`, `edge_spearman`, `edge_permutation_p`, `loading_spearman`, `loading_permutation_p`, `external_pc1_variance_explained`, `pc1_orientation_correlation_with_frozen_score`, `split_half_median`, `split_half_q05`, `split_half_q95`, `split_half_valid_repeats`, `estimable`, `nonestimable_reason`, `random_panel_empirical_p`, `edge_q_bh_8`, `loading_q_bh_8`, `external_canine_representation_class`, `n_frozen_genes`, `n_common_genes`, `coverage_fraction`, `minimum_gene_loo_correlation`, `median_gene_loo_correlation`, `cohort_independence_evidence_class`, `independent_cohort_wording_allowed`, `locked_external_canine_interpretation`, `large_module_leave_one_out_guardrail`, `manuscript_guardrail`
