"""Tests for the FastAPI web interface (uses the ASGI test client, no socket needed)."""
from __future__ import annotations

from starlette.testclient import TestClient  # requires httpx

from pipeline.web.app import dapp


def test_health():
    with TestClient(dapp) as c:
        r = c.get("/api/health")
        assert r.status_code == 200
        assert r.json() == {"status": "ok"}


def test_index_renders():
    with TestClient(dapp) as c:
        r = c.get("/")
        assert r.status_code == 200
        assert "Molecular Docking Pipeline" in r.text


def test_results_empty_when_unknown_dir():
    with TestClient(dapp) as c:
        r = c.get("/results", params={"out": "C:/does/not/exist"})
        assert r.status_code == 200  # graceful "no results" fallback, not a 500


def test_run_rejects_invalid_receptor_with_friendly_error():
    # A receptor that is neither an existing file nor a PDB ID should yield a
    # friendly error page (400) rather than a crash.
    with TestClient(dapp) as c:
        r = c.post(
            "/run",
            data={
                "receptor": "a path that is not a file and is not an id!",
                "ligand": "OC(=O)c1ccccc1C(=O)O",
                "box_size": "20,20,20",
            },
        )
        assert r.status_code == 400
        assert "error" in r.text  # rendered through the friendly results template


def test_run_box_conflict_rejected():
    with TestClient(dapp) as c:
        r = c.post(
            "/run",
            data={
                "receptor": "demo/receptor.pdb",
                "ligand": "OC(=O)c1ccccc1C(=O)O",
                "box_center": "1,2,3",
                "box_source": "demo/pocket_btn.pdb",
                "box_size": "20,20,20",
            },
        )
        assert r.status_code == 400