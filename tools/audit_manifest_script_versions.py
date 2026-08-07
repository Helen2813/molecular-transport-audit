from __future__ import annotations

import ast
import json
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT.parent / "paper4_sarcoma_dog"

CHECKS = [
    {
        "script": (
            SOURCE
            / "scripts"
            / "47_audit_and_lock_gse239948_external_canine_evidence.py"
        ),
        "manifest": (
            SOURCE
            / "results"
            / "tables"
            / "paper4_external_canine_evidence_manifest.json"
        ),
    },
    {
        "script": (
            SOURCE
            / "scripts"
            / "48_gse39055_pathological_necrosis_response.py"
        ),
        "manifest": (
            SOURCE
            / "results"
            / "tables"
            / "GSE39055_necrosis_analysis_manifest.json"
        ),
    },
]


VERSION_NAMES = {
    "SCRIPT_VERSION",
    "VERSION",
    "ANALYSIS_VERSION",
}


def extract_python_version(
    path: Path,
) -> str:
    text = path.read_text(
        encoding="utf-8",
        errors="replace",
    )

    try:
        tree = ast.parse(text)
    except SyntaxError:
        return ""

    for node in tree.body:
        if not isinstance(
            node,
            (ast.Assign, ast.AnnAssign),
        ):
            continue

        if isinstance(node, ast.Assign):
            targets = node.targets
            value = node.value
        else:
            targets = [node.target]
            value = node.value

        for target in targets:
            if (
                isinstance(target, ast.Name)
                and target.id in VERSION_NAMES
                and isinstance(value, ast.Constant)
                and isinstance(value.value, str)
            ):
                return value.value

    match = re.search(
        r"Script version:\s*([A-Za-z0-9._-]+)",
        text,
        re.IGNORECASE,
    )

    if match:
        return match.group(1)

    return ""


def find_manifest_versions(
    payload: object,
    prefix: str = "",
) -> list[tuple[str, str]]:
    found: list[tuple[str, str]] = []

    if isinstance(payload, dict):
        for key, value in payload.items():
            path = (
                f"{prefix}.{key}"
                if prefix
                else key
            )

            if (
                isinstance(value, str)
                and "version" in key.lower()
            ):
                found.append(
                    (path, value)
                )

            found.extend(
                find_manifest_versions(
                    value,
                    path,
                )
            )

    elif isinstance(payload, list):
        for index, value in enumerate(payload):
            path = f"{prefix}[{index}]"

            found.extend(
                find_manifest_versions(
                    value,
                    path,
                )
            )

    return found


def main() -> None:
    print("=" * 80)
    print(
        "Manifest / source script version audit"
    )
    print("=" * 80)

    for check in CHECKS:
        script = check["script"]
        manifest = check["manifest"]

        print()
        print(f"Script:   {script.name}")

        if not script.exists():
            print("  SCRIPT MISSING")
            continue

        version = extract_python_version(
            script
        )

        print(
            f"  Declared source version: "
            f"{version or '<none>'}"
        )

        print(
            f"Manifest: {manifest.name}"
        )

        if not manifest.exists():
            print("  MANIFEST MISSING")
            continue

        payload = json.loads(
            manifest.read_text(
                encoding="utf-8",
                errors="replace",
            )
        )

        versions = find_manifest_versions(
            payload
        )

        if not versions:
            print(
                "  No version fields detected."
            )

        else:
            for key, value in versions:
                print(
                    f"  {key}: {value}"
                )

    print()
    print("=" * 80)


if __name__ == "__main__":
    main()
