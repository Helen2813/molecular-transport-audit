# Molecular Transport Audit

A reproducible bioinformatics framework for auditing whether molecular programs preserve their biological representation and clinical associations across cohorts, species, tissues, and measurement platforms.

## Project status

This repository is under active development.

The project originates from a comparative canine-to-human osteosarcoma study in which molecular representation preservation and prognostic transport were evaluated as distinct properties rather than treated as a single biomarker-replication criterion.

The original osteosarcoma analysis remains the scientific source of truth. This repository generalizes the computational audit into reusable research software without changing the frozen scientific results of the original study.

## Motivation

Cross-cohort and cross-species biomarker studies often treat replication of a clinical association as evidence that the underlying molecular program has transferred successfully.

These are different questions.

A molecular program may:

- remain structurally recognizable while losing its clinical association;
- retain a clinical association despite weak structural representation;
- recur independently in unsupervised analyses;
- be sensitive to measurement or assay choices;
- fail to outperform appropriately matched random controls.

Molecular Transport Audit is designed to evaluate these properties separately.

## Audit dimensions

The framework is organized around five primary dimensions.

### 1. Predictive specificity

Tests whether a selected molecular representation provides evidence beyond matched random controls.

### 2. Molecular representation preservation

Evaluates whether the internal structure of a frozen molecular program is preserved in a target dataset using measures such as edge concordance, loading concordance, and reliability.

### 3. Unsupervised recurrence

Tests whether a frozen program can be recovered without forcing the target data to reproduce it, including latent-factor and blind module-rediscovery analyses.

### 4. Clinical outcome transport

Evaluates frozen scores and directions in target cohorts without outcome-driven re-selection, re-weighting, or re-orientation.

### 5. Measurement robustness

Evaluates whether conclusions are stable across predefined outcome-blind assay or measurement-processing rules.

Multiplicity control and provenance tracking are handled as cross-cutting components of the framework.

## Planned architecture

```text
React
  |
  v
FastAPI
  |
  v
Nextflow DSL2
  |
  +-- Python scientific core
  +-- R / WGCNA
  +-- MOFA2
  +-- Scanpy / AnnData
  |
  v
Versioned result bundle
  |
  +-- metrics
  +-- figures
  +-- manifests
  +-- provenance
  +-- report
The scientific core is intentionally independent of the web interface.
The same audit should eventually be runnable through:
a Nextflow workflow;
a Python API;
a React/FastAPI interface.
Scientific provenance
The original osteosarcoma analysis is maintained separately in:
paper4_sarcoma_dog
The original repository contains the historical analysis scripts, frozen outputs, dataset-specific preprocessing, and manuscript-specific analyses.
This repository contains the reusable implementation.
Reference outputs imported from the original study will include provenance metadata identifying:
source repository;
source script;
script version;
source Git commit;
input hashes where available;
output hashes;
scientific lock status.
The generalized implementation will be regression-tested against those frozen outputs.
Reproducibility principles
The project follows several constraints:
scientific conclusions from the source study are frozen;
target outcomes must not be used to redefine frozen molecular programs;
reusable software must reproduce locked reference analyses before replacing legacy implementations;
stochastic analyses require explicitly documented randomization and resampling policies;
scientific provenance must be machine-readable;
dataset-specific preprocessing must remain separate from reusable computational logic.
See SCIENCE_LOCK.md for the complete guardrails.
Planned bioinformatics stack
The project will use established bioinformatics and reproducibility tools where they are scientifically appropriate:
Python
R
Nextflow DSL2
WGCNA
MOFA2 / mofapy2
AnnData
Scanpy
Ensembl / biomaRt mappings
Docker
Conda or micromamba
SLURM execution profiles
pytest
React
TypeScript
FastAPI
Tools will not be introduced solely for software-stack breadth. Each dependency must serve a defined analysis or reproducibility role.
Scope
This project focuses on auditing frozen molecular programs.
It is not intended to perform multimodal representation learning.
ATAC-seq integration, DNA methylation, proteomics, CNV integration, spatial multi-omics, cross-modal neural representation learning, and new multimodal methods are intentionally reserved for separate future work.
Development roadmap
The implementation roadmap is maintained in PLAN.md.
The software architecture is described in ARCHITECTURE.md.
License
A software license will be selected before the first public research release.
Citation
Citation information will be added with the first versioned research release.
