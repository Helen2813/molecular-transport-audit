"""Validate the persistent API run manager end to end."""

from __future__ import annotations

import sys
import time
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


EXPECTED_CLASSES = {
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
        raise RuntimeError(
            message
        )


def main() -> None:
    print("=" * 80)
    print(
        "Molecular Transport Audit "
        "- API run manager validation"
    )
    print("=" * 80)
    print(
        f"Tool version: {TOOL_VERSION}"
    )

    with TestClient(app) as client:
        response = client.post(
            "/runs",
            json={
                "analysis":
                (
                    "gse239948_external_"
                    "canine_audit"
                )
            },
        )

        require(
            response.status_code == 202,
            (
                "POST /runs failed:\n"
                f"{response.text}"
            ),
        )

        created = response.json()

        run_id = created[
            "run_id"
        ]

        print("")
        print(
            "POST /runs: PASS"
        )
        print(
            f"Run ID: {run_id}"
        )

        early_summary = client.get(
            f"/runs/{run_id}/summary"
        )

        require(
            early_summary.status_code
            == 409,
            (
                "Unverified summary was "
                "unexpectedly exposed."
            ),
        )

        print(
            "Pre-verification summary "
            "guard: PASS"
        )

        previous_stage = None

        start = time.monotonic()

        while True:
            detail_response = (
                client.get(
                    f"/runs/{run_id}"
                )
            )

            require(
                detail_response.status_code
                == 200,
                (
                    "GET /runs/{run_id} "
                    "failed."
                ),
            )

            detail = (
                detail_response.json()
            )

            stage = detail[
                "stage"
            ]

            if stage != previous_stage:
                print(
                    "Stage: "
                    f"{stage}"
                )
                previous_stage = stage

            status_value = detail[
                "status"
            ]

            if (
                status_value
                == "completed"
            ):
                break

            if status_value == "failed":
                log_response = client.get(
                    (
                        f"/runs/{run_id}/log"
                        "?tail=200"
                    )
                )

                raise RuntimeError(
                    "Managed run failed:\n"
                    f"{detail.get('error')}\n\n"
                    f"{log_response.text}"
                )

            if (
                time.monotonic()
                - start
                > 3600
            ):
                raise RuntimeError(
                    "Managed run validation "
                    "timed out."
                )

            time.sleep(
                2.0
            )

        require(
            detail[
                "validated"
            ]
            is True,
            "Completed run is not validated.",
        )

        require(
            detail[
                "validation_status"
            ]
            == "PASS",
            "Validation status is not PASS.",
        )

        print(
            "Managed scientific run: PASS"
        )

        summary_response = client.get(
            f"/runs/{run_id}/summary"
        )

        require(
            summary_response.status_code
            == 200,
            (
                "Validated summary was not "
                "available."
            ),
        )

        summary = (
            summary_response.json()
        )

        require(
            summary[
                "schema_version"
            ]
            == "0.1.0",
            "Unexpected summary schema.",
        )

        require(
            summary[
                "run"
            ][
                "outcome_loaded"
            ]
            is False,
            (
                "Outcome-loaded guard "
                "failed."
            ),
        )

        observed_classes = {
            item["module"]:
            item[
                "classification"
            ][
                "short"
            ]
            for item
            in summary[
                "modules"
            ]
        }

        require(
            observed_classes
            == EXPECTED_CLASSES,
            (
                "Unexpected module "
                "classification."
            ),
        )

        log_response = client.get(
            f"/runs/{run_id}/log"
            "?tail=50"
        )

        require(
            log_response.status_code
            == 200,
            "Run log endpoint failed.",
        )

        print(
            "GET /runs/{run_id}/summary: "
            "PASS"
        )
        print(
            "GET /runs/{run_id}/log: PASS"
        )

        list_response = client.get(
            "/runs"
        )

        require(
            list_response.status_code
            == 200,
            "GET /runs failed.",
        )

        listed_ids = {
            item["run_id"]
            for item
            in list_response.json()[
                "runs"
            ]
        }

        require(
            run_id in listed_ids,
            "Managed run not found "
            "in registry.",
        )

        print(
            "Persistent run registry: PASS"
        )

    print("")
    print("=" * 80)
    print(
        "API run manager validation: PASS"
    )
    print("=" * 80)
    print(
        f"Run ID: {run_id}"
    )
    print(
        "Status: completed"
    )
    print(
        "Validation: PASS"
    )
    print(
        "Summary schema: v0.1.0"
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
