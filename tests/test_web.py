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