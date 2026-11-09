"""Query/reference span bookkeeping for CIGAR strings.

VACmap used to import the PyPI package `cigar` (0.1.3) for exactly one
operation: the number of query bases a CIGAR string consumes, i.e.
``len(Cigar(cigar_string))``.  That package publishes an sdist only, and its
``setup.py`` bootstraps setuptools 2.0.1 from a URL that no longer serves that
file, so ``pip install cigar`` fails on any current machine (reproduced on the
GitHub runner: ``Failed to build 'cigar' when getting requirements to build
wheel`` -> ``tarfile.ReadError``).

The parsing and the two length helpers below are reproduced verbatim from the
upstream implementation, including its quirks (``Cigar("*").items()`` yields
``(0, None)`` and then raises, which under PEP 479 surfaces as RuntimeError),
so that replacing the dependency cannot change any code path that previously
worked.
"""
from itertools import groupby


class Cigar(object):
    read_consuming_ops = ("M", "I", "S", "=", "X")
    ref_consuming_ops = ("M", "D", "N", "=", "X")

    def __init__(self, cigar_string):
        self.cigar = cigar_string

    def items(self):
        if self.cigar == "*":
            yield (0, None)
            raise StopIteration
        cig_iter = groupby(self.cigar, lambda c: c.isdigit())
        for _group, digits in cig_iter:
            yield int("".join(digits)), "".join(next(cig_iter)[1])

    def __len__(self):
        """Sum of M/I/S/=/X ops, i.e. the number of query bases consumed."""
        return sum(l for l, op in self.items()
                   if op in Cigar.read_consuming_ops)

    def reference_length(self):
        """Sum of M/D/N/=/X ops, i.e. the number of reference bases consumed."""
        return sum(l for l, op in self.items()
                   if op in Cigar.ref_consuming_ops)

    def __str__(self):
        return self.cigar

    def __repr__(self):
        return "Cigar(%r)" % self.cigar