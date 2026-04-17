# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What is VACmap

VACmap is a long-read sequence aligner specifically designed for complex structural variation (SV) discovery. It uses non-linear chaining to detect translocations, inversions, and other complex SVs that linear aligners miss. The underlying index (`vacmap_index`) is a C extension wrapping a minimap2-like k-mer/minimizer index.

## Installation

Do not use provided installation instructions. Use docker image instead.

```bash
docker pull hydragenetics/vacmap:1.0.0
```

To test DeepSomatic on VACmap-produced PACBIO alignments, use docker image:

```bash
docker pull google/deepsomatic:1.9.0
```

## Running VACmap

Use docker image above to run VACmap commands.
Bind host directories by running docker with `-v /Users/andgu885/Documents/GitHub/VACmap/dir_to_bin:/dir`

```bash
# Validation with test data
vacmap -ref testdata/reference.fasta -read testdata/read.fasta -mode S -t 4

# Typical ONT/CLR run
vacmap -ref ref.fasta -read reads.fasta -mode H -t 8 -o out.sorted.bam

# HiFi/CCS run
vacmap -ref ref.fasta -read reads.fasta -mode L -t 8 -o out.sorted.bam

# Assembly-to-reference
vacmap -ref ref.fasta -read asm.fasta -mode asm -t 8 --H --fakecigar -workdir /tmp/asm_workdir -o out.sorted.bam
```

Output is auto-sorted when the `-o` path ends in `.sorted.bam`.

VACmap-aligned PACBIO (HiFi) reads are stored in `test_bam/COLO829_aligned.bam`.
Unaligned PACBIO (HiFi) reads are stored in `test_bam/COLO829_unaligned.bam`.

Run DeepSomatic using flag `--model-type PACBIO`.

## DeepSomatic compatibility

VACmap output is compatible with DeepSomatic without extra flags. A default read group (`ID: 1, SM: sample`) is injected automatically when `--rg-id` is not provided. The `PL` field is intentionally omitted from the default — pass `--rg-pl PACBIO` (or `ONT`) to include it. Omitting `--rg-pl` will trigger a `MISSING_PLATFORM_VALUE` warning from Picard ValidateSamFile.


## Linting

```bash
flake8 . --count --select=E9,F63,F7,F82 --show-source --statistics
```

CI runs flake8 on push/PR to `main` (`.github/workflows/python-app.yml`). There are no automated tests beyond the lint check.

## Architecture

### Entry point

`src/vacmap/vacmap` — a Python script (no `.py` extension) installed as the `vacmap` CLI. It parses arguments via `argparse`, builds a parameter dict `pdict`, then dispatches to the appropriate mode module.

### Mode dispatch

Based on `-mode`, the entry point dynamically imports one of:

| Mode | Module | Use case |
|------|--------|----------|
| `H`  | `mammap_clrnano.py` | ONT / PacBio CLR (high error rate) |
| `L`  | `mammap_ccs.py`     | PacBio HiFi (low error rate) |
| `S`  | `mammap_sensitive.py` | Sensitive mode for small variants (<100 bp) |
| `R`  | `mammap_noprefercloser.py` | Fixed variation penalty; sensitive to translocations/gene conversion |
| `asm`| `mammap_asm.py`     | Full genome / assembly alignment |

Each mode module implements the core alignment logic independently (seed lookup → non-linear chaining → extension → SAM output). They share helpers from `output_functions.py`.

### Key dependencies

- **`vacmap_index`** (pip package `vacmap-index`): C extension providing the minimizer index (`mp.Aligner`). Wraps minimap2-style seeding.
- **`numba`**: JIT-compiled inner loops for chaining (`@njit` decorators).
- **`pysam`**: SAM/BAM I/O.
- **`edlib`**: Edit-distance-based alignment extension.
- **`mappy`**: Used only in `index.py` (standalone index builder) and `vacsim/`.

### Pre-built index

`index.py` is a standalone script to pre-build a `.mmi` index file:
```bash
python index.py ref.fasta ref.mmi
```
Pass the `.mmi` path as `-ref` to skip index rebuild at runtime.

### vacsim

`vacsim/vacsim.py` is a read simulator (separate from the aligner) used to generate synthetic test data.

## Key implementation details

- The non-linear chaining step is what differentiates VACmap from minimap2 — it allows anchors to form chains that represent rearranged or inverted segments.
- `asm` mode writes intermediate files to `-workdir`; other modes are stateless (stream stdin → stdout/BAM).
- Multiprocessing: reads are chunked and dispatched via `multiprocessing.Pool`; results are collected and written in order.
- Output to `.sorted.bam` triggers an internal `pysam.sort` call after alignment.
