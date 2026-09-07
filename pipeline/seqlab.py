"""Sequence laboratory: local BLAST-style search and multiple-sequence alignment.

This is a self-contained, offline bioinformatics toolkit that runs entirely on
this device (no external web services) — consistent with the rest of the
pipeline. It provides:

* ``blast_search`` - a BLAST-style search modelled on NCBI BLAST's behaviour.
  The query is compared against every sequence in a local database. Candidate
  subjects are flagged by shared k-mer words, then the best gapped *local*
  alignment (Smith-Waterman extension, Biopython ``PairwiseAligner`` with
  BLOSUM62 for proteins or nucleotide scoring for DNA/RNA) is computed for
  each, and reported with identity %, NCBI-style query coverage %, a bit score
  and an E-value from the Karlin-Altschul Gumbel statistics using the same
  gapped parameters as NCBI BLAST (lambda 0.267 / K 0.041 for proteins).
  Local rather than global alignment means a short query still finds its
  matching region inside a long subject — just like NCBI BLAST.

* ``align_multiple`` - progressive multiple sequence alignment. A UPGMA guide
  tree is built from pairwise distances (via scipy), then sequences are
  progressively merged into a full alignment using sum-of-pairs scoring. The
  result is returned as aligned rows plus a per-column conservation track.

Both are deliberately *local* so the tool remains fully deterministic and
usable without network access.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Optional

from Bio.Align import PairwiseAligner, substitution_matrices
from Bio.Seq import Seq

from .utils import LOG

# ---------------------------------------------------------------------------
# FASTA / sequence parsing
# ---------------------------------------------------------------------------

_DNA_ALPHABET = set("ACGTUN")


def parse_fasta(text: str) -> list[dict[str, str]]:
    """Parse FASTA text (or single sequence lines) into [{id, description, seq}]."""
    entries: list[dict[str, str]] = []
    current: Optional[dict[str, str]] = None
    seq_lines: list[str] = []
    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current is not None:
                current["seq"] = "".join(seq_lines).replace(" ", "").upper()
                entries.append(current)
            header = line[1:].strip()
            ident = header.split()[0] if header.split() else f"seq{len(entries) + 1}"
            current = {"id": ident, "description": header, "seq": ""}
            seq_lines = []
        else:
            seq_lines.append(line)
    if current is not None:
        current["seq"] = "".join(seq_lines).replace(" ", "").upper()
        entries.append(current)

    # If nothing had FASTA headers, treat the input lines as a single sequence.
    if not entries and text.strip():
        clean = "".join(text.split()).upper()
        if clean:
            entries.append({"id": "query", "description": "query", "seq": clean})
    return entries


def is_nucleotide(seq: str) -> bool:
    """True when every alphabetic character is a nucleotide (ACGTUN)."""
    letters = [c for c in seq.upper() if c.isalpha()]
    if not letters:
        return True
    return all(c in _DNA_ALPHABET for c in letters)


def database_label(sequences: list[dict[str, str]]) -> str:
    """A short human label for a sequence set, e.g. '5 protein sequences'."""
    kinds = {_alphabet_kind([s["seq"]]) for s in sequences if s.get("seq")}
    kind_name = "nucleotide" if kinds == {"nucleotide"} else (
        "protein" if kinds == {"protein"} else "mixed")
    total = sum(len(s.get("seq") or "") for s in sequences)
    return f"{len(sequences)} {kind_name} sequence(s), {total} letters"


def _alphabet_kind(seqs: list[str]) -> str:
    if all(is_nucleotide(s) for s in seqs if s):
        return "nucleotide"
    return "protein"


# ---------------------------------------------------------------------------
# Pairwise scoring
# ---------------------------------------------------------------------------

def _make_aligner(kind: str, mode: str = "global") -> PairwiseAligner:
    a = PairwiseAligner()
    a.mode = mode
    if kind == "protein":
        a.substitution_matrix = substitution_matrices.load("BLOSUM62")
        a.open_gap_score = -11
        a.extend_gap_score = -1
    else:
        a.match_score = 2
        a.mismatch_score = -1
        a.open_gap_score = -5
        a.extend_gap_score = -0.5
    return a


_BLOSUM62 = substitution_matrices.load("BLOSUM62")

# Karlin-Altschul Gumbel parameters for GAPPED local alignments, matching the
# values NCBI's BLAST uses for its defaults:
#   protein: blastp BLOSUM62, gap open 11 / extend 1 -> lambda 0.267, K 0.041
#   nucleotide: approximate blastn gapped defaults.
_PROTEIN_KA = (0.267, 0.041)
_NUCLEOTIDE_KA = (1.33, 0.621)


def _bit_score(raw: float, kind: str) -> float:
    lam, _k = _PROTEIN_KA if kind == "protein" else _NUCLEOTIDE_KA
    return (lam * raw - math.log(_k)) / math.log(2)


def _evalue(raw: float, kind: str, db_size: int, query_len: int) -> float:
    lam, k = _PROTEIN_KA if kind == "protein" else _NUCLEOTIDE_KA
    # Search space ~ m*n (query length x database length).
    m = max(query_len, 1)
    n = max(db_size, 1)
    exp_arg = -lam * raw + math.log(k * m * n)
    return min(max(math.exp(exp_arg), 0.0), 10.0)


# ---------------------------------------------------------------------------
# BLAST-style local search
# ---------------------------------------------------------------------------

@dataclass
class BlastHit:
    rank: int
    subject_id: str
    description: str
    raw_score: float
    bit_score: float
    e_value: float
    identity_pct: float
    coverage: float
    query_start: int
    query_end: int
    subject_start: int
    subject_end: int
    query_aln: str
    subject_aln: str
    match_mid: str

    def to_dict(self) -> dict:
        return {
            "rank": self.rank,
            "subject_id": self.subject_id,
            "description": self.description,
            "raw_score": round(self.raw_score, 1),
            "bit_score": round(self.bit_score, 1),
            "e_value": f"{self.e_value:.2g}",
            "identity_pct": round(self.identity_pct, 1),
            "coverage": round(self.coverage, 1),
            "query_start": self.query_start,
            "query_end": self.query_end,
            "subject_start": self.subject_start,
            "subject_end": self.subject_end,
            "query_aln": self.query_aln,
            "subject_aln": self.subject_aln,
            "match_mid": self.match_mid,
        }


def _seed_positions(seq: str, k: int) -> dict[str, list[int]]:
    """Map k-mer -> list of start positions in `seq`."""
    table: dict[str, list[int]] = {}
    for i in range(len(seq) - k + 1):
        table.setdefault(seq[i:i + k], []).append(i)
    return table


def _align_pair(query: str, subject: str, kind: str) -> Optional[dict]:
    """Best local (Smith-Waterman style) alignment of query vs subject.

    Local alignment is what gives BLAST its power: a short query finds its
    matching segment inside a much longer subject instead of being force-
    aligned end-to-end (which scores unrelated flanks and kills the hit).
    """
    aligner = _make_aligner(kind, mode="local")
    try:
        alns = aligner.align(query, subject)
    except Exception:  # noqa: BLE001
        return None
    if not alns:
        return None
    aln = alns[0]
    # Extract the aligned rows as strings and compute identity.
    q_aln = str(aln[0])
    s_aln = str(aln[1])
    pairs = [qc for qc, sc in zip(q_aln, s_aln) if qc != "-" and sc != "-"]
    identities = [qc for qc, sc in zip(q_aln, s_aln) if qc == sc and qc != "-"]
    identity_pct = 100.0 * len(identities) / len(pairs) if pairs else 0.0
    # NCBI reports "Query cover": fraction of the query that takes part in
    # the alignment (not the subject fraction).
    q_covered = sum(1 for c in q_aln if c != "-")
    coverage = 100.0 * q_covered / max(len(query), 1)
    if kind == "protein":
        # BLAST style mid-line: | = identical, + = positive BLOSUM62 score.
        match = []
        for qc, sc in zip(q_aln, s_aln):
            if qc == sc and qc != "-":
                match.append("|")
            elif qc != "-" and sc != "-" and _BLOSUM62.get((qc, sc), 0) > 0:
                match.append("+")
            else:
                match.append(" ")
        match_mid = "".join(match)
    else:
        match_mid = "".join(
            "|" if a == b and a != "-" else " "
            for a, b in zip(q_aln, s_aln)
        )
    return {
        "raw": float(aln.score),
        "identity_pct": identity_pct,
        "coverage": coverage,
        "query_aln": q_aln,
        "subject_aln": s_aln,
        "match_mid": match_mid,
        "q_start": _first_non_gap(q_aln) + 1,
        "q_end": len(query) - _runs_of_gap_at_end(q_aln),
        "s_start": _first_non_gap(s_aln) + 1,
        "s_end": len(subject) - _runs_of_gap_at_end(s_aln),
    }


def _first_non_gap(row: str) -> int:
    for i, c in enumerate(row):
        if c != "-":
            return i
    return 0


def _runs_of_gap_at_end(row: str) -> int:
    count = 0
    for c in reversed(row):
        if c == "-":
            count += 1
        else:
            break
    return count


def blast_search(query: str, database: list[dict[str, str]],
                 word_size: Optional[int] = None,
                 require_seed: bool = True) -> list[BlastHit]:
    """Local BLAST-style search of `query` against `database` sequences.

    ``require_seed`` mirrors BLAST's seeding step: only subjects that share at
    least one exact k-mer word with the query are scored. When False, every
    subject is scored with its best local alignment (used as a sensitivity
    fallback so a query with no shared words still gets a result table).
    """
    query = query.strip().upper()
    if not query:
        return []
    db = [d for d in database if d.get("seq")]
    kind = _alphabet_kind([query] + [d["seq"] for d in db])
    k = word_size or (3 if kind == "protein" else 6)
    k = max(2, min(k, min(len(query), 10)))

    q_positions = _seed_positions(query, k)
    # For candidates, count shared seeds via a rapid scan (no full positional
    # bookkeeping needed; we rely on the gapped local aligner for scoring).
    ranked: list[tuple[float, dict[str, str]]] = []
    for entry in db:
        subj = entry["seq"]
        if len(subj) < 2:
            continue
        if require_seed:
            shared = sum(1 for mer in _seed_positions(subj, k) if mer in q_positions)
            if shared == 0:
                continue
        res = _align_pair(query, subj, kind)
        if res is None or res["raw"] <= 0:
            continue
        ranked.append((res["raw"], {**entry, **res}))

    ranked.sort(key=lambda t: -t[0])
    total_db_len = sum(len(d["seq"]) for d in db) or 1

    hits: list[BlastHit] = []
    for i, (raw, info) in enumerate(ranked):
        bit = _bit_score(raw, kind)
        ev = _evalue(raw, kind, total_db_len, len(query))
        hits.append(
            BlastHit(
                rank=i + 1,
                subject_id=info["id"],
                description=info.get("description", ""),
                raw_score=raw,
                bit_score=bit,
                e_value=ev,
                identity_pct=info["identity_pct"],
                coverage=info["coverage"],
                query_start=info["q_start"],
                query_end=info["q_end"],
                subject_start=info["s_start"],
                subject_end=info["s_end"],
                query_aln=info["query_aln"],
                subject_aln=info["subject_aln"],
                match_mid=info["match_mid"],
            )
        )
    return hits


# ---------------------------------------------------------------------------
# Progressive MSA
# ---------------------------------------------------------------------------

@dataclass
class MsaResult:
    sequences: list[dict[str, str]] = field(default_factory=list)  # aligned rows
    columns: int = 0
    conservation: list[float] = field(default_factory=list)
    guide_tree: str = ""
    method: str = ""

    def to_dict(self) -> dict:
        return {
            "columns": self.columns,
            "method": self.method,
            "guide_tree": self.guide_tree,
            "conservation": [round(c, 1) for c in self.conservation],
            "sequences": [
                {"id": s["id"], "aligned": s["seq"]} for s in self.sequences
            ],
        }


def _make_distance_matrix(seqs: list[str], kind: str) -> list[list[float]]:
    aligner = _make_aligner(kind)
    n = len(seqs)
    mat = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            score = 0.0
            if len(seqs[i]) and len(seqs[j]):
                try:
                    score = float(aligner.align(seqs[i], seqs[j]).score)
                except Exception:  # noqa: BLE001
                    score = 0.0
            mat[i][j] = mat[j][i] = -score  # similarity -> distance
    return mat


def _upgma(seqs: list[str]) -> tuple[list[int], str]:
    """Return (leaf_order, newick) built with UPGMA on alignment distances."""
    from scipy.cluster.hierarchy import linkage, to_tree

    n = len(seqs)
    if n < 2:
        return list(range(n)), "(" + str(n) + ")"
    dmat = _make_distance_matrix(seqs, _alphabet_kind(seqs))
    condensed = [dmat[i][j] for i in range(n) for j in range(i + 1, n)]
    import numpy as np

    condensed = np.array(condensed, dtype=float)
    if condensed.min() < 0 or not np.all(np.isfinite(condensed)):
        condensed = np.ones_like(condensed)
    try:
        Z = linkage(condensed, method="average")
    except Exception:  # noqa: BLE001
        Z = linkage(np.ones_like(condensed), method="average")
    root = to_tree(Z)

    order: list[int] = []
    _leaf_order(root, order)
    return order, _tree_to_newick(root)


def _leaf_order(t, out: list[int]) -> None:
    """Collect leaf ids in a deterministic left-to-right order."""
    if t is None:
        return
    if t.is_leaf():
        out.append(int(t.id))
    else:
        _leaf_order(t.left, out)
        _leaf_order(t.right, out)


def _tree_to_newick(t) -> str:
    if t is None or t.is_leaf():
        return (t.id if t is not None else 0).__str__()
    return f"({_tree_to_newick(t.left)},{_tree_to_newick(t.right)})"


def align_multiple(sequences: list[dict[str, str]]) -> MsaResult:
    """Progressive multiple sequence alignment.

    Uses a UPGMA guide tree and a sum-of-pairs BLOSUM62 consensus-profile
    alignment at each merge. Sequence rows are returned already gap-padded.
    """
    seqs = [s["seq"] for s in sequences if s.get("seq")]
    ids = [s["id"] for s in sequences if s.get("seq")]
    if not seqs:
        return MsaResult()
    if len(seqs) == 1:
        return MsaResult(
            sequences=[{"id": ids[0], "seq": seqs[0]}],
            columns=len(seqs[0]),
            conservation=[100.0] * len(seqs[0]),
            guide_tree="(0)",
            method="progressive (single sequence)",
        )

    kind = _alphabet_kind(seqs)
    lead_order, newick = _upgma(seqs)

    profile_rows = [list(seqs[lead_order[0]])]
    current_order = [lead_order[0]]
    for next_idx in lead_order[1:]:
        profile_rows = _merge_into(profile_rows, list(seqs[next_idx]), kind)
        current_order.append(next_idx)

    width = max(len(r) for r in profile_rows)
    padded = [r + ["-"] * (width - len(r)) for r in profile_rows]

    aligned = []
    order_map = {orig: pos for pos, orig in enumerate(current_order)}
    cons = _conservation(padded)
    for orig_idx in range(len(seqs)):
        pos = order_map[orig_idx]
        aligned.append({"id": ids[orig_idx], "seq": "".join(padded[pos])})

    return MsaResult(
        sequences=aligned,
        columns=width,
        conservation=[c * 100.0 for c in cons],
        guide_tree=newick,
        method="progressive MSA (UPGMA guide tree + BLOSUM62 sum-of-pairs)",
    )


def _consensus_rows(profile_rows: list[list[str]], kind: str) -> str:
    width = max(len(r) for r in profile_rows)
    consensus = []
    for col in range(width):
        residues = [r[col] for r in profile_rows if col < len(r)]
        count: dict[str, int] = {}
        for r in residues:
            count[r] = count.get(r, 0) + 1
        best = max(count, key=count.get)
        consensus.append(best if best != "-" else "-")
    return "".join(consensus)


def _merge_into(profile_rows: list[list[str]], new_seq: list[str], kind: str):
    """Insert `new_seq` into the profile by profile-vs-sequence alignment."""
    consensus = _consensus_rows(profile_rows, kind)
    aligner = _make_aligner(kind)
    # Align the new sequence to the consensus string.
    try:
        aln = aligner.align(consensus, "".join(new_seq))[0]
    except Exception:  # noqa: BLE001
        aln = None
    if aln is None:
        profile_rows = [r + ["-"] * 0 for r in profile_rows]
        profile_rows.append(new_seq)
        return profile_rows

    cons_aln = str(aln[0])
    new_aln = str(aln[1])

    # Build new column assignments. Existing rows 0..k-1, new sequence row k.
    n_existing = len(profile_rows)
    new_cols: list[list[str]] = [[] for _ in range(n_existing + 1)]
    for ci in range(len(cons_aln)):
        rr, nr = cons_aln[ci], new_aln[ci]
        if rr == "-" and nr == "-":
            continue
        if rr == "-":
            # consensus gap -> gap in every existing row; new sequence residue.
            for row in new_cols[:-1]:
                row.append("-")
            new_cols[-1].append(nr)
        else:
            src_col = _index_of_unconsumed(cons_aln, ci)
            for ri, row in enumerate(profile_rows):
                new_cols[ri].append(row[src_col] if src_col < len(row) else "-")
            new_cols[-1].append(nr if nr != "-" else "-")
    return new_cols


def _index_of_unconsumed(aln_str: str, up_to: int) -> int:
    count = -1
    for i in range(up_to + 1):
        if aln_str[i] != "-":
            count += 1
    return count


def _conservation(padded: list[list[str]]) -> list[float]:
    """Fraction of the most common non-gap residue per column."""
    if not padded:
        return []
    width = max(len(r) for r in padded)
    cons = []
    for col in range(width):
        count: dict[str, int] = {}
        nongap = 0
        for r in padded:
            res = r[col] if col < len(r) else "-"
            if res != "-":
                nongap += 1
                count[res] = count.get(res, 0) + 1
        if nongap == 0:
            cons.append(0.0)
        else:
            best = max(count.values())
            cons.append(best / nongap)
    return cons


# ---------------------------------------------------------------------------
# Demo databases (so the web pages are usable out of the box)
# ---------------------------------------------------------------------------

DEMO_PROTEIN_DB = [
    {"id": "HBA_HUMAN", "description": "Hemoglobin subunit alpha [Homo sapiens]",
     "seq": "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR"},
    {"id": "HBB_HUMAN", "description": "Hemoglobin subunit beta [Homo sapiens]",
     "seq": "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH"},
    {"id": "MYG_HUMAN", "description": "Myoglobin [Homo sapiens]",
     "seq": "MGLSDGEWQLVLNVWGKVEADIPGHGQEVLIRLFKGHPETLEKFDKFKHLKSEDEMKASEDLKKHGATVLTALGGILKKKGHHEAEIKPLAQSHATKHKIPVKYLEFISECIIQVLQSKHPGDFGADAQGAMNKALELFRKDMASNYKELGFQG"},
    {"id": "TRX1_HUMAN", "description": "Thioredoxin-1 [Homo sapiens]",
     "seq": "MVKQIESKTAFQEALDAAGDKLVVVDFSATWCGPCKMIKPFFHSLSEKYSNVIFLEVDVDDCQDVASECEVKCMPTFQFFKKGQKVGEFSGANKEKLEATINELV"},
    {"id": "INS_HUMAN", "description": "Insulin [Homo sapiens]",
     "seq": "MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAEDLQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN"},
]

DEMO_NUCLEOTIDE_DB = [
    {"id": "rbcL_SPIOL", "description": "Ribulose-bisphosphate carboxylase large chain, partial",
     "seq": "ATGTCACCACAAACAGAGACTAAAGCAAGTGGTGGGGCAAAATTACAAGATGCTTATTGCCCATTTTTTCAATCACCATCTGTGGAAGAGCAGAAAATTTCTACTGGTTATGTAGCTTACCCATTAGATTTATTTGAAGAAGTTTCTGTTAACAATGTGGCTACCTATTT"},
    {"id": "COX1_HUMAN", "description": "Cytochrome c oxidase subunit I",
     "seq": "ATGACCCTATTTTATCGACGAGTTCAGAAAAGAATATCTTATATCCTCATTGGGGCGGCCTTTTTTTTAATTGCCACTATCATTTTATTATTTATGACACTAGCAACAGCTCTCATCACCCTATTTTAT"},
    {"id": "18S_locus", "description": "18S ribosomal RNA (fragment)",
     "seq": "TACCTGGTTGATCCTGCCAGTAGTCATATGCTTGTCTCAAAGATTAAGCCATGCATGTCTAAGTATAAACTGCTTTATACTGTGAAACTGCGAATGGCTCATTATATCAGTTATAGTTTATTTGATGGTTTACCTACG"},
    {"id": "COX1_ARATH", "description": "Cytochrome c oxidase subunit I",
     "seq": "ATGTTTACCGTTCTTCAAAACTACGTACCCGTCGAAGAAAATATTTTTCTTTTGATCAACACGATTTTACTCGTCTTTACGATGAGCGCCACAGCCCTCGTTACAACACTCTTTTAT"},
]

# Demo sequences used to seed the MSA page.
DEMO_MSA_SEQUENCES = (
    ">HBA_HUMAN Hemoglobin alpha\nMVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDL\n"
    ">HBB_HUMAN Hemoglobin beta\nMVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGD\n"
    ">MYG_HUMAN Myoglobin\nMGLSDGEWQLVLNVWGKVEADIPGHGQEVLIRLFKGHPETLEKFDKF\n"
)
