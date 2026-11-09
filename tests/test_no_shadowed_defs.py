"""Guard against re-introduced shadowed top-level definitions.

Python keeps only the LAST binding of a top-level name, so a second `def foo`
silently shadows the first.  The mammap_* modules used to carry 6-12 copies of
the same function (547 dead copies / 82,277 lines); patching an early copy had
no effect at runtime, which is what made the #19 / MD-CS bugs so hard to see.
This test fails the moment a duplicate comes back.
"""
import ast
import glob
import os

import pytest

SRC = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                   "src", "vacmap")
MODULES = sorted(glob.glob(os.path.join(SRC, "mammap_*.py")))


def top_level_defs(path):
    tree = ast.parse(open(path, encoding="utf-8").read())
    return [n.name for n in tree.body
            if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef))]


def test_modules_found():
    assert len(MODULES) == 5, MODULES


@pytest.mark.parametrize("path", MODULES, ids=os.path.basename)
def test_no_duplicate_top_level_defs(path):
    names = top_level_defs(path)
    seen, dups = set(), []
    for n in names:
        if n in seen:
            dups.append(n)
        seen.add(n)
    assert not dups, (
        f"{os.path.basename(path)}: {len(dups)} shadowed definition(s) — only the "
        f"last one runs: {sorted(set(dups))}")


@pytest.mark.parametrize("path", MODULES, ids=os.path.basename)
def test_module_body_parses_and_has_defs(path):
    # cheap smoke test: the file still parses and is not accidentally emptied
    assert len(top_level_defs(path)) > 30