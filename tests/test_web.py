"""Tests for the FastAPI web interface (uses the ASGI test client)."""
from __future__ import annotations

import time
from pathlib import Path

from starlette.testclient import TestClient  # requires httpx

from pipeline import orchestrator
from pipeline.web.app import dapp, _STATE, _LOCK


def test_health():
    with TestClient(dapp) as c:
        r = c.get("/api/health")
        assert r.status_code == 200
        assert r.json() == {"status": "ok"}


def test_index_renders_form():
    with TestClient(dapp) as c:
        r = c.get("/")
        assert r.status_code == 200
        assert "New docking run" in r.text
        assert "receptor" in r.text
        assert "Run docking" in r.text


def test_health_status_endpoint():
    with TestClient(dapp) as c:
        r = c.get("/api/status")
        assert r.status_code == 200
        assert r.json()["running"] is False


def test_submit_missing_fields_returns_errors():
    # Empty form: both receptor and ligand missing -> friendly error re-render.
    with TestClient(dapp) as c:
        r = c.post("/submit", data={"box_size": "20,20,20"})
        assert r.status_code == 200
        assert "Receptor is required" in r.text
        assert "Ligand is required" in r.text


def test_submit_bad_box_returns_errors():
    with TestClient(dapp) as c:
        r = c.post(
            "/submit",
            data={
                "receptor": "demo/receptor.pdb",
                "ligand": "C1C2C(C(=O)CCSC(=O)NCCC1)SCCC2C(=O)O",
                "box_center": "a,b,c",
                "box_size": "20,20,20",
            },
        )
        assert r.status_code == 200
        assert "box" in r.text and "Invalid" in r.text  # parsed via validation


def test_submit_box_conflict_returns_errors():
    with TestClient(dapp) as c:
        r = c.post(
            "/submit",
            data={
                "receptor": "demo/receptor.pdb",
                "ligand": "C1C2C(C(=O)CCSC(=O)NCCC1)SCCC2C(=O)O",
                "box_center": "1,2,3",
                "box_source": "demo/pocket_btn.pdb",
                "box_size": "20,20,20",
            },
        )
        assert r.status_code == 200
        assert "EITHER" in r.text


def test_results_empty_when_unknown_dir():
    with TestClient(dapp) as c:
        r = c.get("/results", params={"out": "C:/does/not/exist"})
        assert r.status_code == 200  # graceful "no results" fallback, not a 500
        assert "No results yet" in r.text


def test_results_error_page_renders():
    with _LOCK:
        import copy

        old = copy.deepcopy(_STATE)
        _STATE["error"] = "SomeError: test failure"
        _STATE["results"] = None
    try:
        with TestClient(dapp) as c:
            r = c.get("/results")
            assert r.status_code == 200
            assert "Run failed" in r.text
            assert "test failure" in r.text
    finally:
        with _LOCK:
            _STATE.clear()
            _STATE.update(old)


def test_submit_valid_launches_run(monkeypatch):
    """A valid /submit must kick off a background job that lands in _STATE."""

    def fake_run(cfg, verbose=False):
        # Simulate the pipeline writing results.json under cfg.resolved_out().
        out = cfg.resolved_out()
        out.mkdir(parents=True, exist_ok=True)
        (out / "results.json").write_text(
            '{"docking": {"best_pose": {"affinity_kcal_mol": "-4.2"}}}',
            encoding="utf-8",
        )
        return {"docking": {"best_pose": {"affinity_kcal_mol": "-4.2"}}}

    monkeypatch.setattr(orchestrator, "run", fake_run)

    with _LOCK:
        old = dict(_STATE)
        _STATE.update({"running": False, "results": None, "out_dir": None,
                       "error": None})
    try:
        with TestClient(dapp) as c:
            r = c.post(
                "/submit",
                data={
                    "receptor": "demo/receptor.pdb",
                    "ligand": "demo/ligand.sdf",
                    "box_source": "demo/pocket_btn.pdb",
                    "box_size": "10,10,10",
                    "exhaustiveness": "1",
                    "num_modes": "1",
                    "seed": "42",
                },
            )
            assert r.status_code == 200
            assert "Running your docking" in r.text

            # Wait for the background thread to finish writing state.
            deadline = time.time() + 20
            while time.time() < deadline:
                with _LOCK:
                    done = _STATE["running"] is False and _STATE["results"] is not None
                if done:
                    break
                time.sleep(0.2)
            with _LOCK:
                assert _STATE["error"] is None
                assert _STATE["results"]["docking"]["best_pose"]["affinity_kcal_mol"] == "-4.2"
                assert _STATE["out_dir"] is not None
                assert (Path(_STATE["out_dir"]) / "results.json").exists()
    finally:
        with _LOCK:
            _STATE.clear()
            _STATE.update(old)