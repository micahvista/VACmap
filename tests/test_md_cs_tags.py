"""MD:Z: / cs:Z: must be populated for CIGARs built from 'M' operations.

VACmap emits 'M'-only CIGARs unless ``--eqx`` is given, but the MD/CS
builders only handled '=', 'X', 'D' and 'I' and returned ('', '') on anything
else.  Every alignment therefore shipped ``MD:Z:`` and ``cs:Z:`` with an empty
value (issue #19 reporter ran with ``--MD`` and got nothing usable).

Convention (unchanged, see get_bam_dict_str call sites): the query argument is
the aligned slice, i.e. the CIGAR may open with nS but the sequence does not
contain the clipped bases - soft clips are skipped without advancing readloc.
"""
import importlib

import pytest

MODULES = [
    "mammap_clrnano",
    "mammap_sensitive",
    "mammap_asm",
    "mammap_ccs",
    "mammap_noprefercloser",
]

# (cigar, target, query, md, cs-short) - query is the ALIGNED SLICE
CASES = [
    ("6M", "ACGTAC", "ACGTAC", "6", ":6"),
    # leading soft clip + one mismatch in the aligned block
    ("2S6M", "ACGTAC", "ACGAAC", "3T2", ":3*ta:2"),
    # insertion, deletion and two mismatches
    ("3M1I2M1D2M", "ACGTACGT", "ACGTTACG", "5^C0G0T", ":3+t:2-c*gc*tg"),
]


@pytest.mark.parametrize("module", MODULES)
@pytest.mark.parametrize("cigar,target,query,md,cs", CASES)
def test_md_cs_populated_for_M_ops(module, cigar, target, query, md, cs):
    mergecigar_md_cs = importlib.import_module(
        "vacmap." + module).mergecigar_md_cs
    cigar_out, md_out, cs_out, _ = mergecigar_md_cs(cigar, target, query, True)
    assert md_out == md, "%s %s -> MD %r (want %r)" % (module, cigar, md_out, md)
    assert cs_out == cs, "%s %s -> cs %r (want %r)" % (module, cigar, cs_out, cs)
    # the CIGAR handed to the caller must survive untouched
    assert cigar_out == cigar


@pytest.mark.parametrize("module", MODULES)
def test_long_cs_form(module):
    mergecigar_md_cs = importlib.import_module(
        "vacmap." + module).mergecigar_md_cs
    _, md_out, cs_out, _ = mergecigar_md_cs("2S6M", "ACGTAC", "ACGAAC", False)
    assert md_out == "3T2"
    assert cs_out == "=ACG*ta=AC"