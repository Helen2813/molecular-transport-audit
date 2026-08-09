"""Pydantic models for API contracts."""

from __future__ import annotations

from typing import Literal

from pydantic import BaseModel, ConfigDict


class StrictModel(BaseModel):
    model_config = ConfigDict(
        extra="forbid",
    )


RunStatus = Literal[
    "queued",
    "running",
    "completed",
    "failed",
]

RunStage = Literal[
    "queued",
    "starting",
    "direct_preservation",
    "permutation_inference",
    "reliability",
    "random_controls",
    "summary_assembly",
    "verification",
    "completed",
    "failed",
]


class HealthResponse(StrictModel):
    status: Literal["ok"]
    service: str
    api_version: str
    summary_schema_version: str


class Coverage(StrictModel):
    n_frozen_genes: int | None
    n_common_genes: int | None
    fraction: float | None


class DirectPreservation(StrictModel):
    edge_spearman: float | None
    loading_spearman: float | None
    external_pc1_variance_explained: (
        float | None
    )
    pc1_orientation_correlation: (
        float | None
    )


class Inference(StrictModel):
    edge_permutation_p: float | None
    edge_q_bh_8: float | None
    edge_extreme_count: int | None
    loading_permutation_p: float | None
    loading_q_bh_8: float | None
    loading_extreme_count: int | None


class Reliability(StrictModel):
    split_half_median: float | None
    split_half_q05: float | None
    split_half_q95: float | None
    valid_repeats: int | None
    minimum_gene_loo_correlation: (
        float | None
    )
    median_gene_loo_correlation: (
        float | None
    )


class RandomControl(StrictModel):
    matching_mode: str
    n_random_panels: int | None
    null_median: float | None
    null_q05: float | None
    null_q95: float | None
    empirical_p: float | None
    fallback_draws: int | None


ClassificationShort = Literal[
    "strong",
    "partial",
    "limited",
    "no_clear",
]

ClassificationFull = Literal[
    (
        "strong_external_canine_"
        "representation_preservation"
    ),
    (
        "partial_external_canine_"
        "representation_preservation"
    ),
    (
        "limited_specific_"
        "external_canine_signal"
    ),
    (
        "no_clear_external_canine_"
        "representation_preservation"
    ),
]


class Classification(StrictModel):
    short: ClassificationShort
    full: ClassificationFull


class ModuleResult(StrictModel):
    module: str
    coverage: Coverage
    direct_preservation: DirectPreservation
    inference: Inference
    reliability: Reliability
    random_control: RandomControl
    classification: Classification


class AuditRun(StrictModel):
    analysis_id: str
    status: Literal["completed"]
    reference_cohort: str
    external_cohort: str
    outcome_loaded: bool


class AuditParameters(StrictModel):
    permutations: int
    split_half_repeats: int
    random_panels: int
    random_panel_matching_mode: str


class AuditDiagnostics(StrictModel):
    gene_leave_one_out_rows: int
    module_count: int


class AuditSummary(StrictModel):
    schema_version: Literal["0.1.0"]
    assembler_version: str
    generated_at_utc: str
    run: AuditRun
    parameters: AuditParameters
    diagnostics: AuditDiagnostics
    modules: list[ModuleResult]


class ManagedRunState(StrictModel):
    run_id: str
    analysis_id: str
    status: RunStatus
    stage: RunStage

    created_at_utc: str
    started_at_utc: str | None = None
    finished_at_utc: str | None = None

    reference_cohort: str
    external_cohort: str
    outcome_loaded: bool

    pipeline_version: str
    summary_schema_version: (
        str | None
    ) = None

    module_count: int | None = None

    validated: bool = False
    validation_status: (
        Literal["PASS"] | None
    ) = None

    worker_pid: int | None = None
    scientific_pid: int | None = None

    error: str | None = None


class RunDetail(StrictModel):
    run_id: str
    analysis_id: str

    status: RunStatus
    stage: RunStage

    validated: bool
    validation_status: (
        Literal["PASS"] | None
    )

    pipeline_version: str | None
    schema_version: str | None

    created_at_utc: str
    started_at_utc: str | None
    finished_at_utc: str | None
    generated_at_utc: str | None

    reference_cohort: str | None
    external_cohort: str | None
    outcome_loaded: bool | None

    module_count: int | None
    error: str | None


class RunListResponse(StrictModel):
    runs: list[RunDetail]


class RunCreateRequest(StrictModel):
    analysis: Literal[
        "gse239948_external_canine_audit"
    ]


class RunCreateResponse(StrictModel):
    run_id: str
    status: RunStatus
    stage: RunStage
    created_at_utc: str
