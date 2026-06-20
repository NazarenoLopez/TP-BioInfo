"""
Ex4.py — Parseo del reporte BLAST con filtro por patrón.

Recibe un patrón de búsqueda (case-insensitive) y lista todos los hits del
reporte BLAST cuya descripción lo contenga. Con --fetch descarga además las
secuencias completas desde NCBI y las escribe en un archivo FASTA.

Uso:
    python src/Ex4.py <pattern>
    python src/Ex4.py <pattern> --fetch

Ejemplos:
    python src/Ex4.py "Mus musculus"
    python src/Ex4.py "Huntingtin" --fetch
    python src/Ex4.py "Homo sapiens"
"""

import os
import sys

from Bio import Entrez, SeqIO
from Bio.Blast import NCBIXML

# --- Configuración ---
Entrez.email = "brios@frba.utn.edu.ar"

script_dir = os.path.dirname(os.path.abspath(__file__))
proyecto_dir = os.path.dirname(script_dir)
BLAST_REPORT = os.path.join(proyecto_dir, "blast_local_report.xml")
OUTPUT_FASTA = os.path.join(proyecto_dir, "matching_hits.fas")


def parse_accession(hit_def: str) -> str:
    """
    Extrae el accession (sin número de versión) del campo hit_def del XML.

    Ejemplo:
        "P42858.2 RecName: Full=Huntingtin; ..."  →  "P42858"
    """
    token = hit_def.split()[0]   # "P42858.2"
    return token.split(".")[0]   # "P42858"


def buscar_hits(pattern: str, fetch: bool = False) -> list:
    """
    Parsea el reporte BLAST y devuelve los hits cuya descripción contiene el patrón.

    Args:
        pattern: cadena de búsqueda (case-insensitive).
        fetch:   si es True, descarga las secuencias completas desde NCBI.

    Returns:
        Lista de dicts con la info de cada hit que coincide.
    """
    if not os.path.exists(BLAST_REPORT):
        print(f"Error: No se encuentra el reporte BLAST en '{BLAST_REPORT}'")
        sys.exit(1)

    pattern_lower = pattern.lower()
    matches = []

    with open(BLAST_REPORT) as f:
        blast_records = list(NCBIXML.parse(f))

    print(f"Buscando patrón: '{pattern}' en {len(blast_records)} frame(s)\n")
    print("=" * 80)

    for record in blast_records:
        frame_id = record.query

        for hit in record.alignments:
            if pattern_lower not in hit.hit_def.lower():
                continue

            hsp = hit.hsps[0]  # mejor HSP del hit
            accession = parse_accession(hit.hit_def)

            matches.append({
                "frame":        frame_id,
                "accession":    accession,
                "description":  hit.hit_def,
                "e_value":      hsp.expect,
                "bit_score":    hsp.bits,
                "identities":   hsp.identities,
                "align_length": hsp.align_length,
            })

            print(f"  Frame     : {frame_id}")
            print(f"  Accesión  : {accession}")
            # Truncar descripción larga para legibilidad
            # desc = hit.hit_def if len(hit.hit_def) <= 90 else hit.hit_def[:87] + "..."
            print(f"  Descripción: {hit.hit_def}")
            print(f"  Descripción: {hit.hit_def}")
            print(f"  E-value   : {hsp.expect:.2e}")
            print(f"  Bit score : {hsp.bits:.1f}")
            print(f"  Identidad : {hsp.identities}/{hsp.align_length}")
            print()

    print("=" * 80)
    print(f"Total de hits encontrados: {len(matches)}")

    if fetch and matches:
        _fetch_sequences(matches)
    elif fetch and not matches:
        print("No hay hits que descargar.")

    return matches


def _fetch_sequences(matches: list):
    """
    Descarga las secuencias completas de NCBI (db=protein) y las guarda en FASTA.

    Elimina duplicados de accession manteniendo el orden.
    """
    # unique accessions en orden de aparición
    seen = set()
    accessions = []
    for m in matches:
        if m["accession"] not in seen:
            seen.add(m["accession"])
            accessions.append(m["accession"])

    print(f"\nDescargando {len(accessions)} secuencia(s) de NCBI: {accessions}")

    try:
        handle = Entrez.efetch(
            db="protein",
            id=",".join(accessions),
            rettype="fasta",
            retmode="text",
        )
        records = list(SeqIO.parse(handle, "fasta"))
        handle.close()

        SeqIO.write(records, OUTPUT_FASTA, "fasta")
        print(f"Secuencias guardadas en: '{OUTPUT_FASTA}' ({len(records)} registro(s))")

    except Exception as e:
        print(f"Error al descargar secuencias de NCBI: {e}")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)

    patron = sys.argv[1]
    descargar = "--fetch" in sys.argv

    buscar_hits(patron, fetch=descargar)
