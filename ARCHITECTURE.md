# Software Architecture

## 1. Design goals

Molecular Transport Audit is designed as reusable research software rather than as a web application with analysis code embedded inside it.

The architecture separates:

1. scientific computation;
2. workflow orchestration;
3. application services;
4. user interface;
5. reproducibility and provenance.

The scientific core must remain executable and testable without React or FastAPI.

---

## 2. High-level architecture

```text
                    +--------------------+
                    |    React + TS      |
                    |    user interface  |
                    +---------+----------+
                              |
                              | HTTP / JSON
                              v
                    +--------------------+
                    |      FastAPI       |
                    | application layer  |
                    +---------+----------+
                              |
                              | validated run config
                              v
                    +--------------------+
                    |   Nextflow DSL2    |
                    | workflow engine    |
                    +---------+----------+
                              |
             +----------------+----------------+
             |                |                |
             v                v                v
      +-------------+   +-----------+   +-------------+
      |   Python    |   |     R     |   | external    |
      | core logic  |   | WGCNA etc |   | bioinfo     |
      +------+------+   +-----+-----+   | tools       |
             |                |         +------+------+
             +----------------+----------------+
                              |
                              v
                    +--------------------+
                    | Result artifacts   |
                    | metrics / figures  |
                    | manifests / report |
                    +--------------------+

3. Planned repository structure
molecular-transport-audit/
|
|-- core/
|   `-- transport_audit/
|       |-- io.py
|       |-- schemas.py
|       |-- validation.py
|       |-- scoring.py
|       |-- random_controls.py
|       |-- preservation.py
|       |-- recurrence.py
|       |-- outcome_transport.py
|       |-- assay_robustness.py
|       |-- multiplicity.py
|       |-- provenance.py
|       `-- reporting.py
|
|-- workflow/
|   |-- main.nf
|   |-- nextflow.config
|   |-- modules/
|   `-- profiles/
|
|-- r/
|   `-- analysis adapters
|
|-- backend/
|   `-- FastAPI application
|
|-- frontend/
|   `-- React + TypeScript + Vite application
|
|-- configs/
|-- examples/
|-- reference_results/
|-- tests/
|-- containers/
`-- docs/
The exact structure may evolve during implementation, but the separation of responsibilities should remain.

4. Scientific core
The core package is responsible for analysis logic only.
It must not:
start web servers;
manage browser state;
depend on React;
contain UI-specific formatting;
launch jobs through HTTP;
assume osteosarcoma-specific file paths;
contain hard-coded GEO accession logic.
The core should operate on validated data structures and explicit parameters.
Planned modules include:
scoring
Frozen molecular-program scoring.
random_controls
Matched random-panel generation and predictive-specificity testing.
preservation
Structural representation analyses such as edge concordance, loading concordance, and reliability.
recurrence
Outcome-blind latent recurrence and blind rediscovery utilities.
outcome_transport
Frozen clinical-outcome evaluation without target-driven program optimization.
assay_robustness
Evaluation across predefined outcome-blind measurement-processing rules.
multiplicity
Explicit correction within predefined inferential families.
provenance
Hashes, configuration metadata, software versions, randomization information, and result manifests.
reporting
Machine-readable and human-readable result summaries.

5. Workflow layer
Nextflow DSL2 will orchestrate the computational components.
The workflow layer is responsible for:
process dependencies;
intermediate artifacts;
execution isolation;
retries and failures;
local execution;
container execution;
future SLURM execution;
reproducible process invocation.
Scientific statistics should remain in Python or R modules rather than being implemented directly inside Nextflow scripts.
Planned workflow processes include:
VALIDATE_INPUTS
MAP_FEATURES
SCORE_PROGRAMS
RUN_RANDOM_CONTROLS
RUN_PRESERVATION
RUN_LATENT_RECURRENCE
RUN_BLIND_REDISCOVERY
RUN_OUTCOME_TRANSPORT
RUN_ASSAY_ROBUSTNESS
ADJUST_MULTIPLICITY
BUILD_REPORT
WRITE_MANIFEST
Not every run must execute every process.

6. Language interoperability
Python and R will communicate through explicit files and workflow artifacts rather than in-memory bridges where possible.
This avoids making scientific execution depend on interfaces such as rpy2 or reticulate.
Example:
Python preprocessing
        |
        v
versioned tabular artifact
        |
        v
R / WGCNA
        |
        v
structured result file
        |
        v
Python reporting
All intermediate formats must eventually have documented schemas.

7. Bioinformatics tools
Established tools will be used where they correspond to a real scientific component.
Planned examples include:
WGCNA for network/module preservation;
MOFA2 or mofapy2 for outcome-blind latent-factor analyses;
AnnData as a standardized container where appropriate;
Scanpy for single-cell analysis components;
Ensembl or biomaRt for preparation of cross-species mappings.
The workflow will distinguish two orthology modes.
Reproducibility mode
Uses the exact frozen mapping associated with a reference study.
New-project mode
Allows preparation of a new mapping, which must then be materialized, versioned, and hashed before downstream analysis.
Live mappings must not silently change an existing analysis.

8. API layer
FastAPI acts as an application layer, not as a statistical engine.
Responsibilities include:
project creation;
file registration;
validation requests;
run configuration;
workflow launch;
job-state monitoring;
result discovery;
artifact delivery.
A preliminary API may eventually expose endpoints such as:
POST /projects
POST /datasets
POST /programs
POST /runs

GET /runs/{id}
GET /runs/{id}/status
GET /runs/{id}/results
GET /runs/{id}/manifest
The API contract will be defined before frontend integration.

9. Frontend
The web interface will use:
React;
TypeScript;
Vite.
The frontend is intended to make a technically rigorous workflow usable by researchers who do not want to operate Nextflow directly.
Planned views include:
project setup;
data upload;
input QC;
molecular-program registration;
audit selection;
run monitoring;
transportability results;
provenance;
artifact export.
The frontend will initially work against mocked API responses while the scientific backend is being developed.

10. Execution modes
The same scientific workflow should eventually support:
Local
For development and small datasets.
Docker
For reproducible execution.
SLURM
For HPC environments.
The scientific code must not change across execution profiles.

11. Result bundle
A successful run should eventually produce a structured bundle similar to:
run/
|-- metrics/
|-- figures/
|-- tables/
|-- logs/
|-- report/
|-- provenance/
|   `-- manifest.json
`-- config/
    `-- resolved_config.yaml
The manifest should make it possible to determine exactly how the run was produced.

12. Reference-study relationship
The original osteosarcoma repository is a scientific reference implementation.
Molecular Transport Audit must not import code at runtime from that repository.
Instead:
reusable logic is migrated deliberately;
reference outputs are copied with provenance;
the new implementation is regression-tested against those outputs;
the new repository remains independently executable.
This avoids creating an undocumented runtime dependency on the historical research codebase.

13. Future boundary
The current architecture may eventually make multimodal extensions technically possible.
They are nevertheless outside the present scientific scope.
Software extensibility must not be confused with permission to expand the current paper's scientific claims.

