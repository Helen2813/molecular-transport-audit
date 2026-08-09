"""Health routes."""

from __future__ import annotations

from fastapi import APIRouter

from mta_api.config import (
    API_VERSION,
    SERVICE_NAME,
    SUMMARY_SCHEMA_VERSION,
)
from mta_api.models import HealthResponse


router = APIRouter(
    tags=["health"],
)


@router.get(
    "/health",
    response_model=HealthResponse,
)
def health() -> HealthResponse:
    return HealthResponse(
        status="ok",
        service=SERVICE_NAME,
        api_version=API_VERSION,
        summary_schema_version=(
            SUMMARY_SCHEMA_VERSION
        ),
    )
