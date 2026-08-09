"""Smoke-test the Molecular Transport Audit API."""

from __future__ import annotations

import sys
from pathlib import Path

from fastapi.testclient import TestClient


TOOL_VERSION = "0.1.0"

ROOT = Path(__file__).resolve().parents[1]
API_DIR = ROOT / "api"

if str(API_DIR) not in sys.path:
    sys.path.insert(
        0,
        str(API_DIR),
    )

from mta_api.main import app


RUN_ID = "unified_audit_validation"

EXPECTED_CLASSIFICATION = {
    "M34": "strong",
    "M11": "no_clear",
    "M24": "no_clear",
    "M40": "strong",
}


def require(
    condition: bool,
    message: str,
) -> None:
    if not condition:
        raise RuntimeError(message)


def main() -> None:
    print("=" * 80)
    print(
        "Molecular Transport Audit "
        "- API smoke test"
    )
    print("=" * 80)
    print(
        f"Tool version: {TOOL_VERSION}"
    )

    with TestClient(app) as client:
        health_response = client.get(
            "/health"
        )

        require(
            health_response.status_code == 200,
            "GET /health failed.",
        )

        health = health_response.json()

        require(
            health["status"] == "ok",
            "Health status is not ok.",
        )

        require(
            health[
                "summary_schema_version"
            ]
            == "0.1.0",
            "Unexpected summary schema.",
        )

        print(
            "GET /health: PASS"
        )

        runs_response = client.get(
            "/runs"
        )

        require(
            runs_response.status_code == 200,
            "GET /runs failed.",
        )

        runs_payload = (
            runs_response.json()
        )

        run_ids = {
            item["run_id"]
            for item in runs_payload["runs"]
        }

        require(
            RUN_ID in run_ids,
            (
                "Validated unified audit "
                "was not discovered."
            ),
        )

        print(
            "GET /runs: PASS"
        )

        detail_response = client.get(
            f"/runs/{RUN_ID}"
        )

        require(
            detail_response.status_code
            == 200,
            "GET /runs/{run_id} failed.",
        )

        detail = detail_response.json()

        require(
            detail["validated"] is True,
            "Run is not marked validated.",
        )

        require(
            detail["validation_status"]
            == "PASS",
            "Run validation is not PASS.",
        )

        require(
            detail["analysis_id"]
            == (
                "gse239948_external_"
                "canine_audit"
            ),
            "Unexpected analysis ID.",
        )

        require(
            detail["outcome_loaded"]
            is False,
            (
                "API reports that outcome "
                "data were loaded."
            ),
        )

        print(
            "GET /runs/{run_id}: PASS"
        )

        summary_response = client.get(
            f"/runs/{RUN_ID}/summary"
        )

        require(
            summary_response.status_code
            == 200,
            (
                "GET /runs/{run_id}/summary "
                "failed."
            ),
        )

        summary = summary_response.json()

        require(
            summary["schema_version"]
            == "0.1.0",
            "Unexpected summary schema.",
        )

        require(
            summary["run"][
                "outcome_loaded"
            ]
            is False,
            "Summary outcome guard failed.",
        )

        require(
            summary["parameters"][
                "random_panel_matching_mode"
            ]
            == (
                "prestandardization_"
                "variance"
            ),
            (
                "Unexpected random-panel "
                "matching mode."
            ),
        )

        observed_classes = {
            item["module"]:
            item["classification"]["short"]
            for item in summary["modules"]
        }

        require(
            observed_classes
            == EXPECTED_CLASSIFICATION,
            (
                "Unexpected module "
                "classification."
            ),
        )

        print(
            "GET /runs/{run_id}/summary: "
            "PASS"
        )

        missing_response = client.get(
            "/runs/does-not-exist"
        )

        require(
            missing_response.status_code
            == 404,
            "Missing run did not return 404.",
        )

        print(
            "404 contract: PASS"
        )

    print("")
    print("=" * 80)
    print(
        "Molecular Transport Audit API: "
        "PASS"
    )
    print("=" * 80)
    print(
        "Validated run: "
        f"{RUN_ID}"
    )
    print(
        "Summary schema: v0.1.0"
    )
    print(
        "Outcome loaded: False"
    )
    print(
        "M34: strong"
    )
    print(
        "M11: no_clear"
    )
    print(
        "M24: no_clear"
    )
    print(
        "M40: strong"
    )
    print("=" * 80)


if __name__ == "__main__":
    main()
