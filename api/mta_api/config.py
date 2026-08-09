"""Application configuration for the Molecular Transport Audit API."""

from __future__ import annotations

import os
from pathlib import Path


API_VERSION = "0.2.0"
SERVICE_NAME = "molecular-transport-audit-api"
SUMMARY_SCHEMA_VERSION = "0.1.0"

SUPPORTED_ANALYSIS_ID = (
    "gse239948_external_canine_audit"
)

REPOSITORY_ROOT = (
    Path(__file__).resolve().parents[2]
)

REPORTS_ROOT = Path(
    os.environ.get(
        "MTA_REPORTS_ROOT",
        str(REPOSITORY_ROOT / "reports"),
    )
).expanduser().resolve()

MANAGED_RUNS_ROOT = (
    REPORTS_ROOT
    / "api_runs"
)

LEGACY_RUN_ROOT = (
    REPORTS_ROOT
    / "unified_audit_validation"
)

SCIENTIFIC_OUTPUT_ROOT = (
    REPOSITORY_ROOT
    / "reports"
    / "unified_audit_validation"
)

SCIENTIFIC_RUNNER = (
    REPOSITORY_ROOT
    / "tools"
    / "run_unified_audit_validation.py"
)
