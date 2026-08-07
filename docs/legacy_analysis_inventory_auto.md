# Legacy Analysis Inventory — Automated Scan

> This file is generated automatically. Axis assignments are suggestions only and are not scientific classifications.

## Source snapshot

- Source repository: `paper4_sarcoma_dog`
- Remote: `https://github.com/Helen2813/paper4_sarcoma_dog.git`
- Branch: `main`
- HEAD: `2e82099920ccb169527543a7435d33977f8caf71`
- Source worktree dirty: `True`
- Inventory tool version: `0.1.0`
- Generated: `2026-08-07T20:31:27.709062+00:00`

## Scan summary

- Scripts scanned: **56**
- Scripts with detected randomization: **38**
- Scripts with detected version text: **38**
- Numeric script IDs represented by multiple files: **6**

### Automated axis suggestions

- `clinical_outcome_transport`: 27
- `measurement_robustness`: 12
- `molecular_representation_preservation`: 10
- `multiplicity_provenance`: 41
- `predictive_specificity`: 25
- `unsupervised_recurrence`: 14

## Script inventory

| # | Script | Version | SHA-256 | Git commit | Random | Automated axis suggestion |
|---:|---|---|---|---|:---:|---|
| 1 | `01_load_clean_explore_gse238110.py` |  | `2004f9f3ae` | `4f1e0d68e2` |  |  |
| 5 | `05_univariate_cox_screen.py` |  | `f1a76dd27d` | `395882ef49` |  | clinical_outcome_transport |
| 6 | `06_conditional_cox_mb_selection.py` |  | `2fbbaf1b4d` | `833d912b51` |  | predictive_specificity |
| 7 | `07_true_iamb_mb_ablation.py` |  | `ef397e3e0f` | `833d912b51` | yes | clinical_outcome_transport; predictive_specificity |
| 8 | `08_permutation_univariate_cox_control.py` |  | `4e8c830452` | `833d912b51` | yes | clinical_outcome_transport |
| 9 | `09_pca_clinical_structure_qc.py` |  | `99355af832` | `833d912b51` | yes | clinical_outcome_transport |
| 10 | `10_nested_cv_survival_benchmark.py` |  | `4b7070a945` | `833d912b51` | yes | predictive_specificity |
| 11 | `11_nested_cv_random_control_diagnostics.py` |  | `ed12587f8e` | `f8566b78e2` | yes | predictive_specificity |
| 12 | `12_build_rna_master_candidate_evidence_table.py` |  | `6d90c3f4a3` | `f8566b78e2` |  | predictive_specificity |
| 13 | `13_rna_module_signal_analysis.py` |  | `5519bfd75a` | `f8566b78e2` | yes | predictive_specificity |
| 14 | `14_rna_module_random_control_benchmark.py` |  | `bd77fd3f7e` | `f8566b78e2` | yes | predictive_specificity |
| 15 | `15_ortholog_mapping_dog_to_human.py` |  | `ed72d23f2f` | `f8566b78e2` |  |  |
| 16 | `16_statistical_tests_for_random_controls.py` |  | `ab138dda53` | `f8566b78e2` |  | predictive_specificity |
| 17 | `17_ortholog_mapping_qc_transfer_sets.py` |  | `884010edc9` | `f8566b78e2` |  |  |
| 18 | `18_meta_proliferation_adjustment.py` |  | `f001e890bd` | `f8566b78e2` | yes | clinical_outcome_transport; multiplicity_provenance |
| 19 | `19_residualized_proliferation_signal_analysis.py` |  | `b21760297a` | `f8566b78e2` |  | clinical_outcome_transport; multiplicity_provenance |
| 20 | `20_proliferation_overlap_crossfit_sensitivity.py` |  | `dfc9bafc77` | `364e8b8e9c` | yes | predictive_specificity |
| 21 | `21_finalize_canine_transfer_programs.py` |  | `ebcde60aa2` | `364e8b8e9c` | yes | multiplicity_provenance; clinical_outcome_transport |
| 21 | `21_finalize_canine_transfer_programs_FIXED_V2.py` | 21-fixed-ortholog-qc-v2 | `8da88bbb0f` | `364e8b8e9c` | yes | multiplicity_provenance; clinical_outcome_transport |
| 22 | `22_prepare_human_osteosarcoma_cohorts.py` | 22-human-cohort-preparation-v1 | `2b4b973420` | `364e8b8e9c` | yes | multiplicity_provenance; clinical_outcome_transport |
| 23 | `23_external_human_validation.py` | 23-human-external-validation-v1 | `677c29509c` | `364e8b8e9c` | yes | multiplicity_provenance; clinical_outcome_transport; predictive_specificity |
| 24 | `24_external_validation_robustness_audit.py` | 24-external-validation-robustness-audit-v1 | `92233907cc` | `364e8b8e9c` | yes | clinical_outcome_transport; multiplicity_provenance; predictive_specificity |
| 25 | `25_prepare_gse39055_third_human_cohort.py` | 25-gse39055-preparation-v1 | `f5fd842d65` | `364e8b8e9c` | yes | multiplicity_provenance; clinical_outcome_transport |
| 26 | `26_validate_gse39055_rfs.py` | 26-gse39055-rfs-validation-v1 | `774bbf4dff` | `364e8b8e9c` | yes | predictive_specificity; multiplicity_provenance; clinical_outcome_transport |
| 27 | `27_cross_cohort_module_preservation.py` | 27-cross-cohort-module-preservation-v1 | `052725018e` | `364e8b8e9c` | yes | clinical_outcome_transport; predictive_specificity; multiplicity_provenance |
| 28 | `28_conservative_module_preservation_audit.py` | 28-conservative-module-preservation-audit-v1 | `a2f9631f1e` | `364e8b8e9c` | yes | molecular_representation_preservation; multiplicity_provenance; clinical_outcome_transport |
| 29 | `29_lock_cross_species_evidence.py` | 29-lock-cross-species-evidence-v1 | `11fa0408b1` | `364e8b8e9c` |  | clinical_outcome_transport; multiplicity_provenance |
| 30 | `30_generate_paper4_locked_figures_tables.py` | 30-paper4-locked-figures-tables-v1 | `e9846e5721` | `364e8b8e9c` | yes | clinical_outcome_transport; molecular_representation_preservation; multiplicity_provenance |
| 31 | `31_gse39055_assay_quality_diagnostic.py` | 31-gse39055-assay-quality-diagnostic-v1 | `a34e101c96` | `3db7185a3d` |  | measurement_robustness; multiplicity_provenance; clinical_outcome_transport |
| 31 | `31_gse39055_assay_quality_diagnostic_v2.py` | 31-gse39055-assay-quality-diagnostic-v2 | `758832ea7c` | `3db7185a3d` |  | measurement_robustness; multiplicity_provenance; clinical_outcome_transport |
| 32 | `32_patkar_tme_subtype_convergence.py` | 32-patkar-tme-convergence-v1 | `9c762ee503` | `3db7185a3d` | yes | multiplicity_provenance |
| 33 | `33_projectwide_multiplicity_assay_aware_lock.py` | 33-projectwide-multiplicity-assay-aware-lock-v1 | `647c20f3aa` | `3db7185a3d` | yes | multiplicity_provenance; clinical_outcome_transport; measurement_robustness |
| 34 | `34_prepare_multistudy_factor_inputs.py` | 34-prepare-multistudy-factor-inputs-v1 | `1a08efdab1` | `3db7185a3d` |  | measurement_robustness; clinical_outcome_transport; multiplicity_provenance; unsupervised_recurrence |
| 35 | `35_run_mofapy2_multigroup.py` | 35-mofapy2-multigroup-factor-v1 | `ca7c18c911` | `3db7185a3d` | yes | unsupervised_recurrence; multiplicity_provenance |
| 36 | `36_align_frozen_modules_to_mofa.py` | 36-align-frozen-modules-to-mofa-v1 | `3d85958f29` | `3db7185a3d` | yes | predictive_specificity; measurement_robustness; unsupervised_recurrence; multiplicity_provenance |
| 36 | `36_align_frozen_modules_to_mofa_v2.py` | 36-align-frozen-modules-to-mofa-v2 | `74aa8325dc` | `3db7185a3d` | yes | predictive_specificity; measurement_robustness; unsupervised_recurrence; multiplicity_provenance |
| 37 | `37_prepare_variable_only_mofa_inputs.py` | 37-prepare-variable-only-mofa-inputs-v1 | `c4e404ef19` | `3db7185a3d` |  | unsupervised_recurrence; measurement_robustness; multiplicity_provenance; clinical_outcome_transport |
| 38 | `38_run_variable_only_mofapy2_multigroup.py` | 38-variable-only-mofapy2-multigroup-v1 | `89be9343bd` | `3db7185a3d` | yes | unsupervised_recurrence; multiplicity_provenance |
| 39 | `39_align_variable_only_mofa.py` | 39-align-variable-only-mofa-v1 | `c2cf9da917` | `3db7185a3d` | yes | predictive_specificity; measurement_robustness; multiplicity_provenance; unsupervised_recurrence |
| 40 | `40_lock_multidimensional_transport_evidence.py` | 40-lock-multidimensional-transport-evidence-v1 | `504e04e1c4` | `3db7185a3d` |  | clinical_outcome_transport; multiplicity_provenance; measurement_robustness; molecular_representation_preservation; unsupervised_recurrence |
| 41 | `41_download_preflight_ammons_single_cell.py` | 41-download-preflight-ammons-single-cell-v1 | `76b159b06b` | `3db7185a3d` |  | multiplicity_provenance |
| 41 | `41_download_preflight_ammons_single_cell_v2.py` | 41-download-preflight-ammons-single-cell-v2 | `96cbe107b2` | `2e82099920` |  | multiplicity_provenance |
| 42 | `42_score_ammons_single_cell_localization.py` | 42-score-ammons-single-cell-localization-v1 | `5f876e15aa` | `` | yes | multiplicity_provenance |
| 43 | `43_audit_lock_ammons_single_cell.py` | 43-audit-lock-ammons-single-cell-v1 | `770eecc518` | `` | yes | multiplicity_provenance |
| 44 | `44_recompute_ammons_six_dog_localization.py` | 44-recompute-ammons-six-dog-localization-v1 | `dbbf957c45` | `` | yes | multiplicity_provenance |
| 45 | `45_dog2_pathological_necrosis_response.py` | 45-dog2-pathological-necrosis-response-v1 | `2b1bad1324` | `` | yes | multiplicity_provenance; predictive_specificity |
| 46 | `46_gse239948_external_canine_representation.py` | 46-gse239948-external-canine-representation-v1 | `54800118be` | `` | yes | molecular_representation_preservation; multiplicity_provenance; predictive_specificity |
| 46 | `46_gse239948_external_canine_representation_v2.py` | 46-gse239948-external-canine-representation-v2 | `eade633a0d` | `` | yes | molecular_representation_preservation; multiplicity_provenance; predictive_specificity |
| 47 | `47_audit_and_lock_gse239948_external_canine_evidence.py` | 47-lock-gse239948-independent-canine-evidence-v2 | `6438f1a818` | `` |  | molecular_representation_preservation; multiplicity_provenance; predictive_specificity; clinical_outcome_transport |
| 48 | `48_gse39055_pathological_necrosis_response.py` | 48-gse39055-pathological-necrosis-response-v1 | `031ef4c592` | `` | yes | multiplicity_provenance; measurement_robustness |
| 49 | `49_gse239948_blind_de_novo_rediscovery.py` | 49-gse239948-blind-de-novo-rediscovery-v1 | `97ff6c046e` | `` | yes | multiplicity_provenance; predictive_specificity; unsupervised_recurrence |
| 49 | `49_gse239948_blind_de_novo_rediscovery_FIXED_V2.py` | 49-gse239948-blind-de-novo-rediscovery-v2 | `40540bd7f1` | `` | yes | multiplicity_provenance; predictive_specificity; unsupervised_recurrence |
| 50 | `50_wgcna_module_preservation_benchmark.py` | 50-wgcna-module-preservation-benchmark-v1 | `8cedf0d018` | `` | yes | molecular_representation_preservation; multiplicity_provenance; unsupervised_recurrence |
| 51 | `51_final_lock_external_canine_triangulation.py` | 51-final-lock-external-canine-triangulation-v1 | `5b4787b5c1` | `` | yes | molecular_representation_preservation; clinical_outcome_transport; multiplicity_provenance; unsupervised_recurrence; predictive_specificity |
| 52 | `52_generate_final_paper4_manuscript_assets.py` | 52-generate-final-paper4-manuscript-assets-v1 | `95398b1984` | `` |  | clinical_outcome_transport; multiplicity_provenance; molecular_representation_preservation; unsupervised_recurrence; measurement_robustness; predictive_specificity |
| 53 | `53_generate_paper4_manuscript_draft.py` | 53-generate-paper4-manuscript-draft-v1 | `676edd33af` | `` | yes | clinical_outcome_transport; molecular_representation_preservation; multiplicity_provenance; unsupervised_recurrence; measurement_robustness; predictive_specificity |

## Review workflow

`legacy_analysis_inventory_auto.csv` is regenerated on every scan.

`legacy_analysis_inventory_review.csv` preserves manually reviewed fields across rescans by matching `relative_path`.

Manual review must determine which scripts are authoritative, reusable, superseded, or study-specific.
