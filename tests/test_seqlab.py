"""Tests for the local BLAST-style search and multiple sequence alignment."""
from __future__ import annotations

from pipeline import seqlab


def test_parse_single_sequence():
    entries = seqlab.parse_fasta("MVLSPADKTNVKAAWGK")
    assert len(entries) == 1
    assert entries[0]["seq"] == "MVLSPADKTNVKAAWGK"
    assert entries[0]["id"] == "query"


def test_parse_fasta_multiple():
    text = ">one desc\nACGT\n>two\nACGA\n"
    entries = seqlab.parse_fasta(text)
    assert len(entries) == 2
    assert entries[0]["id"] == "one"
    assert entries[1]["id"] == "two"
    assert entries[1]["seq"] == "ACGA"


def test_blast_finds_identical_hit():
    db = [{"id": "x", "description": "", "seq": "ACDEFGHIK"}]
    hits = seqlab.blast_search("ACDEFGHIK", db)
    assert len(hits) == 1
    assert hits[0].identity_pct == 100.0
    assert hits[0].e_value < 1e-3


def test_blast_ranks_best_hit_first():
    db = [
        {"id": "far", "description": "", "seq": "IIIIIIIIIIIIIIIIIIII"},
        {"id": "near", "description": "", "seq": "MVLSPADKTNVKAAWGK"},
    ]
    hits = seqlab.blast_search("MVLSPADKTNVKAAWGK", db)
    assert hits
    assert hits[0].subject_id == "near"


def test_blast_nucleotide():
    db = [{"id": "r1", "description": "", "seq": "ATGGCGATCGAATTCCCC"}]
    hits = seqlab.blast_search("ATGGCGATCGAATTCCCC", db)
    assert hits
    assert hits[0].identity_pct == 100.0


def test_blast_short_query_finds_segment_in_long_subject():
    # NCBI-BLAST-like behaviour: a short query must find its conserved segment
    # inside a much longer subject (local, not global, alignment). This was the
    # reported bug where short queries returned "no hits".
    db = [{"id": "x", "description": "", "seq": "MVLSPADKTNVKAAWGK" * 8}]
    hits = seqlab.blast_search("MVLSPADKTNVKAAWGK", db)
    assert hits, "short query must hit its segment in the long subject"
    top = hits[0]
    assert top.identity_pct == 100.0
    assert top.coverage >= 90.0
    assert top.e_value < 1e-6


def test_blast_query_coverage_is_query_based():
    # A short matching fragment of a much longer query shows NCBI-style
    # "query cover" (fraction of the query that aligns), not subject coverage.
    db = [{"id": "x", "description": "", "seq": "MVLSPADKTNVKAAWGK"}]
    hits = seqlab.blast_search("GGGGGGGGGGMVLSPADKTNVKAAWGKGGGGGGGGGG", db)
    assert hits
    assert 40.0 < hits[0].coverage < 70.0


def test_msa_aligns_identical_sequences():
    result = seqlab.align_multiple([
        {"id": "a", "seq": "ACDEFG"},
        {"id": "b", "seq": "ACDEFG"},
    ])
    assert result.columns == 6
    assert ["a", "b"] == [s["id"] for s in result.sequences]
    assert all(s["seq"] == "ACDEFG" for s in result.sequences)
    assert all(c == 100.0 for c in result.conservation)


def test_msa_three_sequences_guide_tree():
    result = seqlab.align_multiple([
        {"id": "a", "seq": "MVLSPADKTNVKAAWGK"},
        {"id": "b", "seq": "MVHLTPEEKSAVTALWK"},
        {"id": "c", "seq": "MGLSDGEWQLVLNVWGK"},
    ])
    assert len(result.sequences) == 3
    assert result.columns > 0
    # Every output row must be the same aligned width.
    widths = {len(s["seq"]) for s in result.sequences}
    assert widths == {result.columns}


def test_msa_single_sequence():
    result = seqlab.align_multiple([{"id": "a", "seq": "ACDEFG"}])
    assert len(result.sequences) == 1
    assert result.columns == 6
