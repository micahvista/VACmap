"""Regression tests for the 'few chains' early exit of decode_hit (mode R).

decode_hit returned a tuple whose 4th element was the local `factor`, but on
the `len(one_mapinfo) <= 2` path `factor` is never assigned: the name lookup
raised NameError, and the worker's bare `except` swallowed it, so every read
with <=2 raw chains was silently dropped from the output.

The 4th element is the redo ratio; the value hit2work_1 uses when there is
nothing to redo is -1., so the early exit must return that.
"""

import numpy as np
import pytest

from vacmap import mammap_noprefercloser as M


class _NoHitsIndex:
    """Stub of vacmap_index's index object: reports zero k-mer hits."""

    def map(self, seq, check_num=None, mid_occ=-1):
        return []


def _call_decode_hit(monkeypatch, chains):
    monkeypatch.setattr(M, 'get_reversed_chain_numpy_rough', lambda arr, L: (False, chains))
    return M.decode_hit(_NoHitsIndex(), {}, 'ACGT' * 25, 100, {}, 15, {})


def test_decode_hit_no_chain_returns_five_tuple(monkeypatch):
    out = _call_decode_hit(monkeypatch, [])
    assert len(out) == 5, "decode_hit must always return 5 values"
    mapq, scores, path, factor, path_list = out
    assert mapq == 0
    assert scores == 0.0
    assert path == []
    assert factor == -1.0, "redo ratio must be hit2work_1's 'no redo' sentinel, not NameError"
    assert path_list == []


def test_decode_hit_two_chains_returns_five_tuple(monkeypatch):
    out = _call_decode_hit(monkeypatch, [1, 2])
    assert len(out) == 5
    assert out[3] == -1.0
    assert out[4] == []