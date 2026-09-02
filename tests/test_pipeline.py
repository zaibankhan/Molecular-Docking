"""Unit tests for the docking pipeline.

Run with:  python -m pytest tests/ -q
"""
from __future__ import annotations

from pathlib import Path

from pipeline import box, docking, analysis, receptor_prep, report


class TestDockingParse:
    def test_parse_log_poses(self):
        log = (
            "AutoDock Vina v1.2.5\n"
            "    1  -8.50  0.000  0.000\n"
            "    2  -7.10  2.123  1.500\n"
        )
        poses, version = docking.parse_log(log, expected_modes=2)
        assert version == "AutoDock Vina v1.2.5"
        assert len(poses) == 2
        assert poses[0].affinity == -8.50
        assert poses[1].index == 2
        assert poses[1].rmsd_lb == 2.123
        assert poses[1].rmsd_ub == 1.500


class TestPaddedColumns:
    def test_right_aligned_numeric(self):
        rec = receptor_prep._padded(
            [
                ("ATOM", 6, "l"),
                ("1", 5, "r"),
                (" ", 1, "l"),
                ("O", 4, "l"),
                (" ", 1, "l"),
                ("ALA", 4, "l"),
                ("A", 1, "l"),
                ("13", 4, "r"),
                (" ", 1, "l"),
                (" ", 3, "l"),
                ("  22.637", 8, "l"),
                ("   5.768", 8, "l"),
                ("  11.762", 8, "l"),
                ("1.00", 6, "l"),
                ("0.00", 6, "l"),
                (" ", 4, "l"),
                (" -0.32", 6, "l"),
                (" ", 2, "l"),
            ]
        ) + (" " + "N")[-2:]
        # Vina reads serial from columns 7-11, right-aligned.
        assert rec[6:11] == "    1"
        # Atom type at columns 79-80.
        assert rec[78:80] == " N"
        # Coordinates at 31-54.
        assert rec[30:38] == "  22.637"


class TestBox:
    def test_explicit_center_override(self):
        class FakeMol:
            pass

        b = box.default_box(FakeMol(), explicit_center=(1.0, 2.0, 3.0), explicit_size=(10, 20, 30))
        assert b.center == (1.0, 2.0, 3.0)
        assert b.size == (10.0, 20.0, 30.0)


class TestAffinityBand:
    def test_bands(self):
        assert "very strong" in report.affinity_band(-9.5)
        assert "strong" in report.affinity_band(-8.0)
        assert "moderate" in report.affinity_band(-6.0)
        assert "weak" in report.affinity_band(-4.0)
        assert "non-specific" in report.affinity_band(-2.0)


class TestAdmetLite:
    def test_lipinski_pass(self):
        # Phthalic acid: small, drug-like.
        res = analysis.admet_lite("O=C(O)c1ccccc1C(=O)O")
        assert res["lipinski_violations"] == 0
        assert res["lipinski_pass"] is True

    def test_invalid_smiles(self):
        res = analysis.admet_lite("not-a-smiles")
        assert "error" in res