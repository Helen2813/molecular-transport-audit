# Osteosarcoma Locked Reference Registry

This directory contains exact copies of selected authoritative outputs from the frozen canine-human osteosarcoma study.

These artifacts are regression targets for the reusable Molecular Transport Audit implementation.

- Registry version: `0.1.0`
- Source repository HEAD at build: `aba92fa57af165b2ef433ae9e581a9278a04e7bf`
- Artifact count: **24**

## Artifacts

| ID | Axis | Regression class | Source SHA-256 |
|---|---|---|---|
| `frozen_canine_program_manifest` | program_freeze | exact_artifact | `fcd34b985ea9` |
| `frozen_strict_gene_weights` | program_freeze | exact_artifact | `4f065aa4c4ed` |
| `frozen_broad_gene_weights` | program_freeze | exact_artifact | `ee5f661677b1` |
| `frozen_scoring_specification` | program_freeze | exact_artifact | `562aab92caa9` |
| `projectwide_primary_multiplicity` | clinical_outcome_transport | deterministic_table | `a1bf689c3adc` |
| `human_primary_outcomes` | clinical_outcome_transport | exact_artifact | `48881f473c4e` |
| `human_structural_preservation_manifest` | molecular_representation_preservation | exact_artifact | `4abe64e72bfd` |
| `human_structural_preservation_classification` | molecular_representation_preservation | deterministic_table | `184fd8eade2b` |
| `gse39055_assay_manifest` | measurement_robustness | exact_artifact | `219dc54b1f4a` |
| `variable_only_mofa_model_manifest` | unsupervised_recurrence | latent_model_manifest | `e96930c1d68b` |
| `variable_only_mofa_alignment_manifest` | unsupervised_recurrence | latent_model_manifest | `be60f367df45` |
| `external_canine_direct_manifest` | molecular_representation_preservation | stochastic_manifest | `eae7ad5cd8f6` |
| `external_canine_direct_structure` | molecular_representation_preservation | stochastic_table | `3a23d206704b` |
| `external_canine_random_controls` | predictive_specificity | stochastic_table | `472b0676317e` |
| `external_canine_blind_manifest` | unsupervised_recurrence | stochastic_manifest | `b57dfff98c8d` |
| `external_canine_blind_rediscovery` | unsupervised_recurrence | stochastic_table | `473895e65c3d` |
| `external_canine_wgcna_manifest` | molecular_representation_preservation | wgcna_manifest | `505fbe4316ea` |
| `external_canine_wgcna_results` | molecular_representation_preservation | wgcna_table | `0b98b6fe9e17` |
| `single_cell_six_dog_manifest` | biological_localization | exact_artifact | `9bd7cec81c81` |
| `single_cell_six_dog_primary_tests` | biological_localization | deterministic_table | `acc4f4b2486d` |
| `external_canine_evidence_lock` | evidence_lock | exact_artifact | `78459cfdb608` |
| `final_analysis_lock` | evidence_lock | exact_artifact | `19e5d7332517` |
| `final_external_canine_triangulation` | evidence_lock | exact_artifact | `196dace12b76` |
| `final_module_interpretation` | evidence_lock | exact_artifact | `ccab4f6e11ce` |

## Duplicate-content audit

No duplicate reference artifacts detected.

## Important

A reference artifact is an immutable test oracle. It is not automatically a reusable software input.
