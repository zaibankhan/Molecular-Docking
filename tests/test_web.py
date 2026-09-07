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


def test_index_renders_landing_options():
    # The landing page shows the tool options, not the docking form.
    with TestClient(dapp) as c:
        r = c.get("/")
        assert r.status_code == 200
        assert "Choose a tool" in r.text
        assert "New docking run" not in r.text
        assert "/dock" in r.text
        assert "/blast" in r.text
        assert "/msa" in r.text
        assert "/history" in r.text


def test_dock_page_renders_form():
    # The docking form lives on the dedicated /dock page.
    with TestClient(dapp) as c:
        r = c.get("/dock")
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


HBA = ("MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKK"
       "VADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHAS"
       "LDKFLASVSTVLTSKYR")


def test_blast_saves_history_and_links_to_msa(monkeypatch, tmp_path):
    import pipeline.web.app as webapp
    monkeypatch.setattr(webapp, "_SEQLAB_HISTORY", tmp_path / "hist.json")
    with TestClient(dapp) as c:
        r = c.post("/blast", data={"query": HBA, "db_type": "protein"})
        assert r.status_code == 200
        assert "hit(s)" in r.text
        assert "Align query + top hits" in r.text  # BLAST -> MSA link

        recs = webapp._load_seqlab_history()
        assert len(recs) == 1
        assert recs[0]["kind"] == "blast"
        assert recs[0]["hits"]
        assert recs[0]["hits"][0]["subject_id"] == "HBA_HUMAN"

        h = c.get("/history")
        assert h.status_code == 200
        assert "Sequence laboratory" in h.text
        assert "BLAST" in h.text

        v = c.get("/seqlab", params={"ts": recs[0]["ts"]})
        assert v.status_code == 200
        assert "Align query + top hits" in v.text  # link survives re-viewing


def test_msa_saves_history_and_links_to_blast(monkeypatch, tmp_path):
    import pipeline.web.app as webapp
    monkeypatch.setattr(webapp, "_SEQLAB_HISTORY", tmp_path / "hist.json")
    fasta = ">A seq one\nACDEFGHIKLMNPQRSTVWY\n>B seq two\nACDEFGHIKLMNPQRSTVWYC\n"
    with TestClient(dapp) as c:
        r = c.post("/msa", data={"sequences": fasta})
        assert r.status_code == 200
        assert "columns" in r.text
        assert "MSA &rarr; BLAST" in r.text  # MSA -> BLAST link

        recs = webapp._load_seqlab_history()
        assert len(recs) == 1
        assert recs[0]["kind"] == "msa"
        assert recs[0]["columns"] > 0
        assert len(recs[0]["sequences"]) == 2

        h = c.get("/history")
        assert h.status_code == 200
        assert "MSA" in h.text

        v = c.get("/seqlab", params={"ts": recs[0]["ts"]})
        assert v.status_code == 200
        assert "Run BLAST" in v.text  # link survives re-viewing


def test_seqlab_view_unknown_ts_redirects_to_history(monkeypatch, tmp_path):
    import pipeline.web.app as webapp
    monkeypatch.setattr(webapp, "_SEQLAB_HISTORY", tmp_path / "hist.json")
    with TestClient(dapp) as c:
        r = c.get("/seqlab", params={"ts": "19990101_000000"},
                  follow_redirects=False)
        assert r.status_code == 307


def test_blast_page_has_no_database_options():
    with TestClient(dapp) as c:
        page = c.get("/blast").text
        assert 'name="db_text"' not in page
        assert 'name="db_type"' not in page
        assert 'name="db_file"' not in page
        assert 'name="query"' in page
        assert 'name="query_file"' in page


def test_blast_auto_db_and_sensitivity_fallback(monkeypatch, tmp_path):
    import pipeline.web.app as webapp
    monkeypatch.setattr(webapp, "_SEQLAB_HISTORY", tmp_path / "hist.json")
    # A query that shares no k-mer word with the built-in protein database
    # must still produce a result table (sensitivity fallback).
    with TestClient(dapp) as c:
        r = c.post("/blast", data={"query": "WWWWWWWWWWWWWWWWWW"})
        assert r.status_code == 200
        assert "sensitivity" in r.text
        assert "hit(s)" in r.text