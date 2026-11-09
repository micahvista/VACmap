"""Regression tests for the in-process reference-index builder.

Before `vacmap/index_builder.py`, the CLI shelled out to `minimap2 -d`; on any
machine without that binary the run died with a bare ``FileNotFoundError``
(observed on the author's box: `which minimap2` empty, testdata run failed).
`vacmap_index.Aligner(..., fn_idx_out=...)` builds the very same file, so the
external tool is no longer required.
"""

import hashlib
import os

import pytest

from vacmap import index_builder

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REF = os.path.join(ROOT, "testdata", "reference.fasta")
SHIPPED = REF + ".w10_k15.mmi"
# md5 of the index committed in testdata/ (built by `minimap2 -d ... -w 10 -k 15`)
SHIPPED_MD5 = "d2e95bb25c50136cbfb49b9e05639eda"


def _md5(path):
    with open(path, "rb") as handle:
        return hashlib.md5(handle.read()).hexdigest()


def test_index_path_matches_cli_convention(tmp_path):
    ref = str(tmp_path / "ref.fa")
    assert index_builder.index_path(ref, 10, 15) == ref + ".w10_k15.mmi"


def test_build_index_is_byte_identical_to_shipped(tmp_path):
    if not os.path.isfile(SHIPPED):
        pytest.skip("bundled testdata index is missing")
    out = str(tmp_path / "built.mmi")
    ok, why = index_builder.build_index(REF, out, w=10, k=15)
    assert ok, why
    assert _md5(out) == _md5(SHIPPED)
    assert _md5(out) == SHIPPED_MD5


def test_build_index_needs_no_minimap2_binary(tmp_path, monkeypatch):
    def boom(*args, **kwargs):
        raise AssertionError("minimap2 fallback must not be reached")

    monkeypatch.setattr(index_builder, "_build_with_minimap2", boom)
    out = str(tmp_path / "built.mmi")
    ok, why = index_builder.build_index(REF, out, w=10, k=15)
    assert ok, why
    assert os.path.getsize(out) > 0


def test_missing_reference_reports_failure_without_partial_file(tmp_path):
    ref = str(tmp_path / "absent.fa")
    out = ref + ".w10_k15.mmi"
    ok, why = index_builder.build_index(ref, out, w=10, k=15)
    assert not ok
    assert "reference not found" in why and "absent.fa" in why
    assert not os.path.exists(out)