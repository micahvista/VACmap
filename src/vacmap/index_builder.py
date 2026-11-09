"""Reference index construction for the vacmap CLI.

The `.mmi` file consumed by `vacmap_index.Aligner` is a plain minimap2 index,
so it can be produced without the external `minimap2` binary: passing
`fn_idx_out=` to `Aligner` builds and writes it in one call.  Verified byte
identical to `minimap2 -d out.mmi ref.fa -w 10 -k 15 -I 1000G`
(md5 d2e95bb25c50136cbfb49b9e05639eda on testdata/reference.fasta).

Memory: the index is built in RAM, same as `minimap2 -d`; the .mmi is roughly
2-3x the reference size for a vertebrate genome (7.3 GB for GRCh38 with
w=10/k=15), so budget accordingly.  The external binary remains available as a
fallback and is required only if the in-process build fails.
"""

import logging
import os
import shutil
import subprocess

__all__ = ["build_index", "index_path"]


def index_path(refpath, w, k):
    """Canonical name the CLI uses for the index of `refpath`."""
    return f"{refpath}.w{w}_k{k}.mmi"


def _build_in_process(refpath, index_name, w, k):
    import vacmap_index

    vacmap_index.Aligner(refpath, fn_idx_out=index_name, w=int(w), k=int(k))
    return os.path.isfile(index_name)


def _build_with_minimap2(refpath, index_name, w, k):
    exe = shutil.which("minimap2")
    if exe is None:
        return False, "minimap2 not found on PATH"
    cmd = [exe, "-d", index_name, refpath, "-w", str(w), "-k", str(k), "-I", "1000G"]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        tail = (proc.stderr or proc.stdout or "").strip().splitlines()[-3:]
        return False, f"minimap2 exited {proc.returncode}: {' / '.join(tail)}"
    return os.path.isfile(index_name), "ok"


def build_index(refpath, index_name=None, w=10, k=15):
    """Build `index_name` from `refpath`.

    Returns (ok, message).  On failure no partial index is left behind by the
    in-process path; a partial file from a failed minimap2 run is removed too.
    """
    if index_name is None:
        index_name = index_path(refpath, w, k)

    if not os.path.isfile(refpath):
        # Aligner() creates nothing for a missing file and does not raise
        return False, f"reference not found: {refpath}"

    logging.info("Building reference index %s (w=%s, k=%s)", index_name, w, k)
    try:
        if _build_in_process(refpath, index_name, w, k):
            return True, "built with vacmap_index"
        reason = "vacmap_index produced no file"
    except Exception as error:  # noqa: BLE001 - fall back to the external binary
        reason = f"vacmap_index build failed: {error}"

    logging.warning("%s; falling back to the minimap2 binary", reason)
    ok, message = _build_with_minimap2(refpath, index_name, w, k)
    if not ok and os.path.isfile(index_name):
        try:
            os.remove(index_name)
        except OSError:
            pass
    if not ok:
        logging.error(
            "Failed to build reference index for %s: %s; %s", refpath, reason, message
        )
        return False, f"{reason}; {message}"
    return True, message