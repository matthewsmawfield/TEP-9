"""Lightweight checks for repository layout and shared utilities."""

from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent


def test_core_pipeline_scripts_exist():
    steps = ROOT / "scripts" / "steps"
    assert steps.is_dir()
    for rel, _ in __import__(
        "scripts.run_all", fromlist=["CORE_STEPS"]
    ).CORE_STEPS:
        assert (steps / rel).is_file(), f"missing step script: {rel}"


def test_extract_step_number():
    from scripts.utils.step_logger import extract_step_number

    assert extract_step_number("step_033_iri_trajectory_profile") == "033"
    assert extract_step_number("step_040a_extract_3d_vectors") == "040A"


def test_required_results_schema_when_present():
    """If primary products exist, they must be non-empty JSON objects."""
    catalog = ROOT / "results" / "step003_archival_flyby_catalog.json"
    if not catalog.is_file():
        pytest.skip("archival catalog not built in this checkout")
    import json

    data = json.loads(catalog.read_text(encoding="utf-8"))
    assert isinstance(data, dict)
    assert "flybys" in data
    assert len(data["flybys"]) >= 1
