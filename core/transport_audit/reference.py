"""Utilities for accessing locked scientific reference artifacts."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import yaml


def repository_root() -> Path:
    """Return the Molecular Transport Audit repository root."""

    return Path(__file__).resolve().parents[2]


def reference_root() -> Path:
    """Return the locked osteosarcoma reference directory."""

    return (
        repository_root()
        / "reference_results"
        / "osteosarcoma_locked"
    )


def load_reference_registry() -> dict[str, Any]:
    """Load the locked osteosarcoma reference registry."""

    path = (
        reference_root()
        / "registry.json"
    )

    if not path.exists():
        raise FileNotFoundError(
            f"Reference registry not found: {path}"
        )

    with path.open(
        "r",
        encoding="utf-8",
    ) as handle:
        payload = json.load(handle)

    if not isinstance(payload, dict):
        raise ValueError(
            "Reference registry must contain a JSON object."
        )

    return payload


def get_reference_record(
    artifact_id: str,
) -> dict[str, Any]:
    """Return one reference-registry record by artifact ID."""

    registry = load_reference_registry()

    matches = [
        artifact
        for artifact in registry.get(
            "artifacts",
            [],
        )
        if artifact.get("id") == artifact_id
    ]

    if not matches:
        raise KeyError(
            f"Reference artifact not found: {artifact_id}"
        )

    if len(matches) > 1:
        raise RuntimeError(
            "Reference registry contains duplicate artifact IDs: "
            f"{artifact_id}"
        )

    return matches[0]


def get_reference_artifact_path(
    artifact_id: str,
) -> Path:
    """Resolve the local immutable copy of a reference artifact."""

    record = get_reference_record(
        artifact_id
    )

    relative_path = record[
        "registry_copy"
    ]["path"]

    path = (
        repository_root()
        / relative_path
    )

    if not path.exists():
        raise FileNotFoundError(
            f"Reference artifact copy not found: {path}"
        )

    return path


def load_tolerance_policy() -> dict[str, Any]:
    """Load metric-specific scientific regression tolerances."""

    path = (
        reference_root()
        / "tolerances.yaml"
    )

    if not path.exists():
        raise FileNotFoundError(
            f"Tolerance policy not found: {path}"
        )

    with path.open(
        "r",
        encoding="utf-8",
    ) as handle:
        payload = yaml.safe_load(handle)

    if not isinstance(payload, dict):
        raise ValueError(
            "Tolerance policy must contain a YAML mapping."
        )

    return payload
