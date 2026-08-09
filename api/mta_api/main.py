"""FastAPI application entry point."""

from __future__ import annotations

from fastapi import FastAPI

from mta_api.config import API_VERSION
from mta_api.routes.health import (
    router as health_router,
)
from mta_api.routes.runs import (
    router as runs_router,
)


app = FastAPI(
    title="Molecular Transport Audit API",
    description=(
        "Validated API access to molecular "
        "transport audit results."
    ),
    version=API_VERSION,
)

app.include_router(
    health_router
)
app.include_router(
    runs_router
)
