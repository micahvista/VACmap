"""Tests for nm_from_cigar in output_functions."""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from vacmap.output_functions import nm_from_cigar


def test_perfect_match():
    assert nm_from_cigar('5M', 'ACGTA', 'ACGTA') == 0


def test_single_mismatch():
    # position 2: G vs T
    assert nm_from_cigar('5M', 'ACGTA', 'ACTTA') == 1


def test_insertion():
    # 3M 2I 3M: 2 inserted bases → NM=2
    # query: ACG + TT (inserted) + GTA  = ACGTTGTA (8 bases)
    # ref:   ACG + GTA = ACGGTA (6 bases)
    assert nm_from_cigar('3M2I3M', 'ACGTTGTA', 'ACGGTA') == 2


def test_deletion():
    # 3M 1D 3M: 1 deleted base → NM=1
    # query: ACG + GTA = ACGGTA (6 bases)
    # ref:   ACG + X (deleted) + GTA = ACGXGTA (7 bases)
    assert nm_from_cigar('3M1D3M', 'ACGGTA', 'ACGXGTA') == 1


def test_soft_clip_not_counted():
    # 3S 5M 2S: soft clips do not contribute to NM
    # query: AAA (soft) + GGGCC (aligned) + TT (soft) = AAAGGGCCTT
    # ref for aligned region: GGGCC
    assert nm_from_cigar('3S5M2S', 'AAAGGGCCTT', 'GGGCC') == 0


def test_soft_clip_with_mismatch():
    # 3S 5M 2S: 1 mismatch in aligned region (position 3 of aligned: X vs C)
    # query aligned portion: GGGXC (pos 3,4,5,6,7 of full read)
    assert nm_from_cigar('3S5M2S', 'AAAGGGXCTT', 'GGGCC') == 1


def test_hard_clip_ignored():
    # 2H 5M 2H: hard clips — no sequence stored, q_pos unchanged
    assert nm_from_cigar('2H5M2H', 'ACGTA', 'ACGTA') == 0


def test_case_insensitive():
    assert nm_from_cigar('5M', 'acgta', 'ACGTA') == 0


def test_mixed_ops():
    # 3M 1I 2M 1D 3M
    # query: ACG (3M) + T (1I) + AG (2M) + GTA (3M) = ACGTAGGTA  (9 query bases)
    # ref:   ACG (3M) + TG (2M) + X (1D) + GTA (3M) = ACGTGXGTA  (9 ref bases)
    # mismatches: 2M: A vs T = 1 mismatch; I: +1; D: +1 → NM=3
    query = 'ACGTAGGTA'
    ref   = 'ACGTGXGTA'
    assert nm_from_cigar('3M1I2M1D3M', query, ref) == 3
