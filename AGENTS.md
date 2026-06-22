# TP-BioInfo — Agent Guide

## Quick start
- Python 3.10+, BioPython ≥1.83. No build system — standalone scripts.
- External binaries required: `blastp`, `makeblastdb` (NCBI BLAST+), `clustalo`, `pepstats`, `patmatmotifs`, `prosextract` (EMBOSS).
- Setup: `python -m venv .venv && pip install -r requirements.txt`

## Project structure
```
HTT_mRNA.gbk           # Input: human HTT mRNA (NM_002111), GenBank format
HTT_ORFs.fas           # Output of Ex1: 6 translated reading frames
blast_local_report.xml # Output of Ex2: local BLAST report (XML)
src/Ex1.py .. Ex5.py   # 5 exercise scripts (run directly)
tests/test_ex1.py      # Only test file
data/orthologs/        # Cached SwissProt FASTA downloads (Ex3)
results/               # MSA outputs (Ex3)
```

## Running
```
python src/Ex1.py                          # GenBank → 6 ORFs FASTA
python src/Ex2.py                          # BLAST local (needs Swissport/)
python src/Ex3.py                          # MSA via Clustal Omega
python src/Ex4.py "Mus musculus"           # Filter BLAST report by pattern
python src/Ex4.py "Huntingtin" --fetch     # + download matching sequences
python src/Ex5.py                          # EMBOSS: pepstats + patmatmotifs (PROSITE)
python src/Ex5.py --prosite-dir <dir>      # idem with custom PROSITE data dir
python -m unittest discover -s tests -p "test_*.py"
```

## Key facts
- **Correct frame = Frame_2_Fwd** — CDS starts at mRNA position 146 (146 mod 3 = 2). Verified by `test_detecta_frame_correcto_HTT`.
- **6 frame naming:** `Frame_{1,2,3}_{Fwd,Rev}` (3 forward + 3 reverse complement).
- **BLAST DB path:** `./Swissport/swissprot` — created via `makeblastdb -in swissprot -dbtype prot`. The `Swissport/` dir is gitignored (large files).
- **Ex4 NCBI fetches** use `Entrez.email = "brios@frba.utn.edu.ar"` for Entrez access.
- **Ex3 caches** ortholog FASTA in `data/orthologs/` (avoids re-download).
- **Proxy SSL fix** (UTN networks): set `trusted-host` in `pip.ini` — see README.md.
- All docs in Spanish, code identifiers in English.

## Testing quirks
- Single test file `test_ex1.py` using `unittest`. Run with `python -m unittest`.
- Tests create temp files for output; no fixtures or external services needed.
- No tests for Ex2–Ex4.

## Conventions
- Relative paths resolved from `Path(__file__).resolve().parent.parent` (repo root).
- No type checking, no linter config — just raw Python.
- Outputs written to repo root (`HTT_ORFs.fas`, `blast_local_report.xml`) or `results/`.
