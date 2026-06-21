#!/usr/bin/env python3
"""
Ex5.py — Ejercicio 5: EMBOSS

Pipeline:
    1. Convierte el mRNA (GenBank) a FASTA de nucleótidos con seqret.
    2. Calcula los ORFs / proteínas posibles con getorf.
    3. Busca dominios/motivos funcionales PROSITE en esas proteínas
       con patmatmotifs.

Input:  archivo de mRNA en formato GenBank (.gbk) -- por defecto HTT_mRNA.gbk
Output: results/<basename>_domains.patmatmotifs (dominios encontrados)

Uso:
    python3 Ex5.py [archivo_mrna.gbk]

Si no se pasa argumento, usa HTT_mRNA.gbk en la raíz del proyecto.
"""

import subprocess
import sys
from pathlib import Path

# ----------------------------------------------------------------------------
# 0. Configuración
# ----------------------------------------------------------------------------

PROJECT_ROOT = Path(__file__).resolve().parent.parent
RESULTS_DIR = PROJECT_ROOT / "results"

# Directorio donde se preparó la base PROSITE con prosextract
# (ajustar si se usó otra ruta al correr `sudo prosextract -prositedir ...`)
PROSITE_DATA_DIR = Path("/usr/share/EMBOSS/data/PROSITE")

# ORF mínimo en nucleótidos a reportar (150 nt = 50 aa). Filtra ruido.
MIN_ORF_SIZE = 150

REQUIRED_PROGRAMS = ["seqret", "getorf", "patmatmotifs"]


# ----------------------------------------------------------------------------
# Utilidades
# ----------------------------------------------------------------------------

def check_program(name: str) -> None:
    """Verifica que un programa EMBOSS esté disponible en el PATH."""
    found = subprocess.run(
        ["which", name], capture_output=True, text=True
    )
    if found.returncode != 0:
        sys.exit(
            f"ERROR: no se encontró el programa EMBOSS '{name}' en el PATH.\n"
            f"       Verificá la instalación con: embossversion"
        )


def run(cmd: list[str], step_name: str) -> None:
    """Corre un comando externo, mostrando qué se ejecuta y abortando con
    mensaje claro si falla."""
    print(f"      $ {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        sys.exit(
            f"ERROR en el paso '{step_name}':\n"
            f"{result.stderr.strip()}"
        )


def split_fasta(multi_fasta: Path, out_dir: Path) -> list[Path]:
    """Separa un multi-FASTA en archivos individuales (una secuencia c/u).
    patmatmotifs procesa de a una secuencia por corrida de forma confiable,
    así que conviene iterar en vez de pasarle el multi-FASTA entero."""
    seq_files: list[Path] = []
    current_lines: list[str] = []
    current_path: Path | None = None
    idx = 0

    def flush():
        if current_path is not None:
            current_path.write_text("".join(current_lines))
            seq_files.append(current_path)

    for line in multi_fasta.read_text().splitlines(keepends=True):
        if line.startswith(">"):
            flush()
            idx += 1
            current_path = out_dir / f"seq_{idx}.fas"
            current_lines = [line]
        else:
            current_lines.append(line)
    flush()

    return seq_files


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------

def main() -> None:
    input_gbk = Path(sys.argv[1]) if len(sys.argv) > 1 else PROJECT_ROOT / "HTT_mRNA.gbk"

    RESULTS_DIR.mkdir(exist_ok=True)
    basename = input_gbk.stem

    mrna_fasta = RESULTS_DIR / f"{basename}_mrna.fasta"
    orfs_fasta = RESULTS_DIR / f"{basename}_ORFs_filtered.fas"
    domains_out = RESULTS_DIR / f"{basename}_domains.patmatmotifs"

    print("=" * 70)
    print(" Ejercicio 5 - EMBOSS: ORFs + dominios PROSITE")
    print("=" * 70)
    print(f"Input (mRNA GenBank): {input_gbk}")
    print(f"Resultados en:        {RESULTS_DIR}")
    print()

    # -- 1. Chequeos previos --------------------------------------------------
    if not input_gbk.is_file():
        sys.exit(f"ERROR: no se encontró el archivo de entrada '{input_gbk}'")

    for prog in REQUIRED_PROGRAMS:
        check_program(prog)

    if not (PROSITE_DATA_DIR / "prosite.lines").is_file():
        sys.exit(
            f"ERROR: no se encontró la base PROSITE indexada en '{PROSITE_DATA_DIR}'.\n"
            f"       Corré primero: sudo prosextract -prositedir ~/PROSITE"
        )

    # -- 2. GenBank -> FASTA de nucleótidos (seqret) --------------------------
    print("[1/3] Convirtiendo GenBank a FASTA con seqret...")
    run(
        [
            "seqret",
            "-sequence", f"genbank::{input_gbk}",
            "-outseq", f"fasta::{mrna_fasta}",
            "-auto",
        ],
        "seqret",
    )
    print(f"      -> {mrna_fasta}\n")

    # -- 3. Cálculo de ORFs y traducción a proteína (getorf) ------------------
    # -find 0    : ORFs clásicos delimitados por START (ATG) y STOP, que es
    #              lo que interesa para "posibles secuencias de proteína".
    # -minsize   : tamaño mínimo del ORF en nucleótidos, filtra ruido.
    # -reverse Y : también busca ORFs en la cadena reversa-complementaria.
    print(f"[2/3] Calculando ORFs con getorf (tamaño mínimo: {MIN_ORF_SIZE} nt)...")
    run(
        [
            "getorf",
            "-sequence", str(mrna_fasta),
            "-outseq", str(orfs_fasta),
            "-find", "0",
            "-minsize", str(MIN_ORF_SIZE),
            "-reverse", "Y",
            "-auto",
        ],
        "getorf",
    )

    n_orfs = orfs_fasta.read_text().count(">")
    print(f"      -> {orfs_fasta}  ({n_orfs} ORFs encontrados)\n")

    # -- 4. Búsqueda de dominios PROSITE (patmatmotifs) -----------------------
    print("[3/3] Buscando dominios PROSITE con patmatmotifs...")

    import tempfile
    with tempfile.TemporaryDirectory() as tmp_dir_str:
        tmp_dir = Path(tmp_dir_str)
        seq_files = split_fasta(orfs_fasta, tmp_dir)

        with domains_out.open("w") as out_f:
            for seq_file in seq_files:
                out_file = seq_file.with_suffix(".out")
                run(
                    [
                        "patmatmotifs",
                        "-sequence", str(seq_file),
                        "-outfile", str(out_file),
                        "-auto",
                    ],
                    "patmatmotifs",
                )
                out_f.write(out_file.read_text())
                out_f.write("\n")

    n_hits = domains_out.read_text().count("Motif = ")

    print(f"      -> {domains_out}\n")
    print("=" * 70)
    print(f" Listo. {n_orfs} ORFs analizados, {n_hits} dominios PROSITE encontrados.")
    print(f" Resultado final: {domains_out}")
    print("=" * 70)


if __name__ == "__main__":
    main()