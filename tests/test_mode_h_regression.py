"""End-to-end regression test for issue #19 (silently dropped reads).

Symptom before the fix: ``vacmap -mode H`` on the bundled testdata wrote a
SAM with **zero** alignment records.  The DP did produce three blocks, but
the second and third start at q>0, so their CIGARs open with a soft clip;
``nm_from_cigar()`` was handed the inner query slice instead of the full
read, walked ``q_pos`` past the end of that slice, raised IndexError, and
the output worker swallowed the exception with a bare ``continue``.

Expected (and documented in testdata/README) output for one read:

    read_1  16     chr1:6814  6830S10821M1D5383M6811S
    read_1  2048   chr1:23017 23015S6814M16S
    read_1  2048   chr1:1     6811M23034S

The middle record is the leading-soft-clip case that used to kill the read.
"""
import os
import pathlib
import subprocess
import sys

import pytest

REPO = pathlib.Path(__file__).resolve().parents[1]

EXPECTED_CIGARS = {
    "6830S10821M1D5383M6811S",
    "23015S6814M16S",
    "6811M23034S",
}


def _run(mode, out_path):
    env = dict(os.environ, PYTHONPATH=str(REPO / "src"))
    cmd = [
        sys.executable,
        str(REPO / "src" / "vacmap" / "vacmap"),
        "-ref", "testdata/reference.fasta",
        "-read", "testdata/read.fasta",
        "-mode", mode,
        "-t", "1",
        "-o", str(out_path),
        "--MD",
    ]
    subprocess.run(cmd, check=True, env=env, cwd=str(REPO),
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


@pytest.mark.slow
@pytest.mark.parametrize("mode", ["H", "S"])
def test_all_alignments_are_emitted(tmp_path, mode):
    out = tmp_path / ("out_%s.sam" % mode)
    _run(mode, out)

    recs = [l.split("\t") for l in out.read_text().splitlines()
            if not l.startswith("@")]

    # The bug: this was 0 even though the DP found three blocks.
    assert len(recs) == 3, "expected 3 records, got %d" % len(recs)
    assert {r[5] for r in recs} == EXPECTED_CIGARS
    assert all(r[0] == "read_1" for r in recs)

    # The alignment that used to disappear: it starts at q>0 (leading soft clip).
    assert any(r[5] == "23015S6814M16S" for r in recs), \
        "leading-soft-clip alignment missing: %r" % [r[5] for r in recs]

    # NM must stay consistent with the CIGAR (I+D bound) and with MD, and the
    # MD:Z:/cs:Z: tags must not be empty (an M-only CIGAR used to return '')
    import re
    op_re = re.compile(r"(\d+)([MIDSH=X])")
    for r in recs:
        tags = {t.split(":")[0]: t.split(":", 2)[2] for t in r[11:]}
        nm = int(tags["NM"])
        assert tags["MD"] != "", "empty MD:Z for %s" % r[5]
        assert tags["cs"] != "", "empty cs:Z for %s" % r[5]
        indels = sum(int(n) for n, op in op_re.findall(r[5]) if op in "ID")
        # mismatches are the bare letters; bases after '^' are deletions
        md, i, mismatches = tags["MD"], 0, 0
        while i < len(md):
            if md[i] == "^":
                i += 1
                while i < len(md) and md[i].isalpha():
                    i += 1
            elif md[i].isdigit():
                while i < len(md) and md[i].isdigit():
                    i += 1
            else:
                mismatches += 1
                i += 1
        assert nm == mismatches + indels, \
            "NM=%d != MD mismatches %d + I/D %d" % (nm, mismatches, indels)