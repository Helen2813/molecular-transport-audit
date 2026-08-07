```markdown
# Development Plan

## Objective

Build a reusable, regression-tested bioinformatics implementation of the molecular transportability audit developed in the osteosarcoma study.

The project has two simultaneous goals:

1. provide a rigorous computational artifact supporting the scientific paper;
2. provide a usable research workflow demonstrating modern bioinformatics software engineering.

Scientific scope is controlled by `SCIENCE_LOCK.md`.

---

# Phase 0 — Project boundary and architecture

## Goal

Establish the project before migrating scientific code.

## Deliverables

- [x] README
- [x] SCIENCE_LOCK
- [x] ARCHITECTURE
- [x] PLAN
- [ ] .gitignore
- [ ] initial repository structure
- [ ] initial versioning convention
- [ ] license decision deferred until public release

## Exit condition

No scientific implementation begins until the source study and reusable software are clearly separated.

---

# Phase 1 — Legacy analysis inventory

## Goal

Map the existing osteosarcoma analysis to reusable computational components.

Do not refactor code during this phase.

## Inventory fields

For every relevant source script record:

- script filename;
- script version;
- Git commit if available;
- scientific purpose;
- framework axis;
- input files;
- output files;
- randomization;
- external packages;
- dataset-specific logic;
- reusable logic;
- known fixes or superseded versions;
- authoritative locked output;
- proposed destination module;
- proposed regression strategy.

## Framework axes

Each reusable analysis should map to one or more of:

1. predictive specificity;
2. molecular representation preservation;
3. unsupervised recurrence;
4. clinical outcome transport;
5. measurement robustness;
6. multiplicity/provenance support.

## Deliverable

`docs/legacy_analysis_inventory.csv`

and a human-readable summary:

`docs/legacy_analysis_inventory.md`

## Exit condition

We know exactly which legacy computations will be migrated and which will remain study-specific.

---

# Phase 2 — Reference-result registry

## Goal

Create authoritative regression targets before changing implementations.

## Tasks

- identify final locked source outputs;
- reject superseded outputs;
- record exact generating script/version;
- record source Git commit where available;
- compute output hashes;
- identify source input hashes where available;
- define metric-specific regression rules;
- document stochastic procedures.

## Proposed structure

```text
reference_results/
`-- osteosarcoma_locked/
    |-- manifests/
    |-- metrics/
    |-- tables/
    `-- tolerances.yaml
Important rule
Reference results must exist before corresponding reusable functions are considered validated.
Exit condition
Every migrated scientific component has an authoritative regression target or an explicit documented reason why one cannot exist.

Phase 3 — Data contracts and scientific core
Goal
Extract reusable scientific logic without changing frozen statistical definitions.
Initial modules
transport_audit/
|-- schemas.py
|-- io.py
|-- validation.py
|-- scoring.py
|-- random_controls.py
|-- preservation.py
|-- recurrence.py
|-- outcome_transport.py
|-- assay_robustness.py
|-- multiplicity.py
|-- provenance.py
`-- reporting.py
Rules
no hard-coded dataset accessions in core;
no hard-coded project directories;
no target-outcome-driven optimization;
inputs and outputs must have explicit schemas;
stochastic functions must accept explicit RNG configuration;
public functions require concise documentation and type hints.
Exit condition
Core components can operate independently of the historical repository.

Phase 4 — Scientific regression testing
Goal
Demonstrate that the reusable implementation reproduces the locked study.
Test categories
Unit tests
Small deterministic function-level tests.
Scientific regression tests
Compare reusable outputs with frozen osteosarcoma reference outputs.
Stochastic regression tests
Reuse persisted resampling schedules where available or apply predefined stochastic-validation rules.
Invariance tests
Test properties such as latent-factor sign or ordering invariance where mathematically appropriate.
Critical requirement
Scientific regression status must be visible and machine-readable.
Exit condition
Migrated analyses reproduce the frozen reference results within predefined metric-specific rules.

Phase 5 — Nextflow DSL2 orchestration
Goal
Create a real bioinformatics workflow around the validated scientific core.
Planned processes
input validation;
feature mapping;
frozen scoring;
random controls;
structural preservation;
latent recurrence;
blind rediscovery;
outcome transport;
assay robustness;
multiplicity adjustment;
report generation;
provenance generation.
Profiles
local;
Docker;
SLURM.
Important rule
Nextflow orchestrates scientific programs but does not redefine their statistical logic.
Exit condition
The osteosarcoma worked example can be executed reproducibly through Nextflow.

Phase 6 — Bioinformatics adapters
Goal
Integrate established tools where scientifically justified.
Planned tools
WGCNA
For network/module preservation.
MOFA2 / mofapy2
For outcome-blind latent recurrence.
AnnData / Scanpy
For supported single-cell localization workflows.
Ensembl / biomaRt
For preparation of frozen cross-species feature maps.
Important rule
Tools are added because they implement required scientific functionality, not to increase the number of technologies listed by the project.
Exit condition
External-tool components are versioned, containerized where appropriate, and included in provenance.

Phase 7 — FastAPI application layer
Goal
Expose workflow execution without coupling the API to statistical implementation details.
Initial responsibilities
create project;
register datasets;
register molecular programs;
validate inputs;
create run;
start workflow;
query run status;
retrieve results;
retrieve manifest.
Exit condition
A complete audit can be launched and inspected through API calls.

Phase 8 — React interface
Goal
Provide a researcher-friendly interface over the same validated workflow.
Technology
React;
TypeScript;
Vite.
Initial screens
project setup;
dataset upload;
input QC;
molecular-program upload;
analysis selection;
run monitor;
results dashboard;
provenance viewer;
artifact export.
Development strategy
Frontend development may begin early against mocked API contracts.
Scientific calculations must never be implemented in the browser.
Exit condition
A researcher can complete the worked example without using the command line.

Phase 9 — Containers and HPC
Goal
Make execution portable.
Deliverables
Docker images;
Docker execution profile;
reproducible Python/R environments;
pinned important scientific dependencies;
SLURM profile;
documented local execution.
Exit condition
The same test workflow passes locally and in the reproducible container environment.

Phase 10 — Worked osteosarcoma example
Goal
Provide a small, understandable example tied directly to the source paper.
The example must use frozen/reference inputs or an appropriately reduced distributable dataset.
It must not require users to reconstruct the entire historical project directory.
Deliverables
example config;
documented inputs;
expected outputs;
reference manifest;
tutorial.
Exit condition
A new user can reproduce the documented example from the standalone repository.

Phase 11 — Documentation and research release
Deliverables
complete README;
installation guide;
Nextflow guide;
UI guide;
input-format documentation;
output-format documentation;
scientific-method documentation;
provenance documentation;
developer documentation;
CITATION.cff;
license;
versioned release.
Exit condition
The repository can be cited and used independently of the original osteosarcoma repository.

Optional Phase 12 — Independent audit-only demonstration
This phase is deliberately NOT part of the initial implementation requirement.
After the software and manuscript are assembled, we will reassess whether a second disease-domain demonstration is scientifically necessary.
If added, it should preferably test the audit framework on a previously defined external molecular signature rather than launch a new biomarker-discovery study.
A full second-domain discovery pipeline is outside the current plan.

Explicitly reserved for future multimodal work
The following are not Paper 1 development tasks:
multimodal representation learning;
ATAC integration;
DNA methylation modeling;
proteomics integration;
CNV integration;
spatial multi-omics;
modality-specific encoders;
contrastive multimodal learning;
shared/private latent neural representations;
a newly trained cross-modal AI model;
broad multi-disease method benchmarking;
method-development simulations with known ground truth.
These belong to a distinct future research project.

Project decision rule
Before adding a task, ask:
Does this task improve implementation, reproducibility, usability, or validation of the existing transportability audit?
If yes, it may belong here.
If the task changes the biological question or introduces a new inferential objective, stop and classify it before implementation.
