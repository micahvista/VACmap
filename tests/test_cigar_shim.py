"""The vendored Cigar helper must behave exactly like the PyPI `cigar` package.

`cigar` cannot be installed on a current machine (its sdist bootstraps
setuptools 2.0.1 from a dead URL), so it was replaced by `vacmap._cigar`.  These
tests pin the replacement: fixed expectations that always run, plus a
randomised differential test against the upstream package whenever it happens
to be importable.
"""
import random

import pytest

from vacmap._cigar import Cigar as ShimCigar


def test_known_query_span():
    assert len(ShimCigar("100M")) == 100
    assert len(ShimCigar("20H20M20S")) == 40          # hard clips are free
    assert len(ShimCigar("10M20S10M")) == 40
    assert len(ShimCigar("5M2D3N4M")) == 9            # D/N consume reference only
    assert len(ShimCigar("4=1X2I3=")) == 10
    assert len(ShimCigar("0M")) == 0


def test_known_items_and_reference_span():
    assert list(ShimCigar("2S3M1I4M2D5M").items()) == [
        (2, "S"), (3, "M"), (1, "I"), (4, "M"), (2, "D"), (5, "M")]
    assert ShimCigar("3M2D4M").reference_length() == 9
    assert ShimCigar("2S3M1I4M2D5M").reference_length() == 14
    assert str(ShimCigar("3M2D4M")) == "3M2D4M"
    assert repr(ShimCigar("3M2D4M")) == "Cigar('3M2D4M')"


@pytest.mark.parametrize("cigar", [
    "6811M23034S",
    "6830S10821M1D5383M6811S",
    "23015S6814M16S",
    "1S1H1S5H1S5M10H",
])
def test_real_cigars_from_the_issue(cigar):
    """CIGAR strings produced by v47.0007 on the issue-19 data set."""
    query_span = sum(int(n) for n, op in ShimCigar(cigar).items()
                     if op in ("M", "I", "S", "=", "X"))
    assert len(ShimCigar(cigar)) == query_span


def test_star_behaves_like_upstream():
    """Upstream yields (0, None) then raises, which PEP 479 turns into RuntimeError."""
    with pytest.raises(RuntimeError):
        list(ShimCigar("*").items())


def test_random_differential_against_upstream():
    upstream = pytest.importorskip(
        "cigar", reason="upstream `cigar` not installed; differential test skipped")
    rng = random.Random(20261109)
    ops = "MIDNSHP=X"
    cigars = ["".join("%d%s" % (rng.randint(0, 500), rng.choice(ops))
                      for _ in range(rng.randint(1, 8))) for _ in range(500)]
    cigars += ["100M", "20H20M20S", "10M20S10M", "5=3X2I4D", "0M", "1M"]
    for cigar in cigars:
        assert list(ShimCigar(cigar).items()) == list(upstream.Cigar(cigar).items()), cigar
        assert len(ShimCigar(cigar)) == len(upstream.Cigar(cigar)), cigar
        assert (ShimCigar(cigar).reference_length()
                == upstream.Cigar(cigar).reference_length()), cigar
        assert str(ShimCigar(cigar)) == str(upstream.Cigar(cigar)), cigar