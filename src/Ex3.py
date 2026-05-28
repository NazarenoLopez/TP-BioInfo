"""
Ejercicio 3 - Multiple Sequence Alignment (MSA) de la proteina Huntingtina.

Descarga ortologos de SwissProt (UniProt curado), los combina con la
proteina humana extraida del CDS del .gbk usado en el Ej1, y corre
Clustal Omega para producir el alineamiento multiple.

Salidas:
    results/HTT_combined.fasta  -> input combinado (query + ortologos)
    results/HTT_msa.aln         -> alineamiento en formato Clustal (legible)
    results/HTT_msa.fasta       -> alineamiento en formato FASTA (para analisis)
"""

import subprocess
import sys
import urllib.request
from pathlib import Path

from Bio import SeqIO


script_dir = Path(__file__).resolve().parent
proyecto_dir = script_dir.parent
gbk_input = proyecto_dir / "HTT_mRNA.gbk"
orthologs_dir = proyecto_dir / "data" / "orthologs"
results_dir = proyecto_dir / "results"
combined_fasta = results_dir / "HTT_combined.fasta"
clustal_aln = results_dir / "HTT_msa.aln"
fasta_aligned = results_dir / "HTT_msa.fasta"

# Ortologos de SwissProt para Huntingtina (verificados manualmente).
# Cada hit aparece en el reporte BLAST del Ej2 con E-value 0.
ORTOLOGOS = [
    ("P42859", "HD_MOUSE", "Mus musculus"),
    ("P51111", "HD_RAT", "Rattus norvegicus"),
    ("P51112", "HD_TAKRU", "Takifugu rubripes"),
]

UNIPROT_URL = "https://rest.uniprot.org/uniprotkb/{}.fasta"


def descargar_ortologos() -> list[Path]:
    """Descarga FASTA de SwissProt (uno por accession). Cachea localmente."""
    orthologs_dir.mkdir(parents=True, exist_ok=True)
    paths = []
    for accession, sp_id, organismo in ORTOLOGOS:
        destino = orthologs_dir / f"{accession}.fasta"
        if destino.exists() and destino.stat().st_size > 0:
            print(f"[OK]   Cache: {destino.name}  ({sp_id}, {organismo})")
        else:
            url = UNIPROT_URL.format(accession)
            print(f"[INFO] Descargando {accession} ({organismo}) desde UniProt...")
            urllib.request.urlretrieve(url, destino)
            print(f"[OK]   Guardado: {destino}")
        paths.append(destino)
    return paths


def extraer_query_humana(gbk_path: Path) -> str:
    """Extrae la proteina anotada del primer CDS del archivo GenBank."""
    record = SeqIO.read(str(gbk_path), "genbank")
    for feature in record.features:
        if feature.type == "CDS":
            translation = feature.qualifiers.get("translation", [None])[0]
            if translation:
                header = (
                    ">sp|P42858|HD_HUMAN Huntingtin "
                    f"OS=Homo sapiens GN=HTT (CDS de {gbk_path.name})"
                )
                return f"{header}\n{translation}\n"
    raise RuntimeError(f"No se encontro un CDS con 'translation' en {gbk_path}")


def combinar_secuencias(query_fasta: str, ortologos_paths: list[Path]) -> int:
    """Concatena la query humana + ortologos en un unico FASTA."""
    results_dir.mkdir(parents=True, exist_ok=True)
    n = 1  # la query humana
    with open(combined_fasta, "w", encoding="ascii", newline="\n") as out:
        out.write(query_fasta)
        for path in ortologos_paths:
            contenido = path.read_text(encoding="ascii")
            if not contenido.endswith("\n"):
                contenido += "\n"
            out.write(contenido)
            n += contenido.count(">")
    return n


def correr_clustalo() -> None:
    """Genera dos formatos de salida: Clustal (legible) y FASTA (para parseo)."""
    base_cmd = ["clustalo", "-i", str(combined_fasta), "--force"]

    print(f"[INFO] Generando alineamiento Clustal...")
    subprocess.run(
        base_cmd + ["-o", str(clustal_aln), "--outfmt", "clu", "-v"],
        check=True,
    )

    print(f"[INFO] Generando alineamiento FASTA...")
    subprocess.run(
        base_cmd + ["-o", str(fasta_aligned), "--outfmt", "fa"],
        check=True,
    )


def main() -> int:
    if not gbk_input.exists():
        print(f"[ERROR] No se encuentra {gbk_input}", file=sys.stderr)
        return 1

    print(f"[INFO] Extrayendo query humana del CDS de {gbk_input.name}...")
    query_fasta = extraer_query_humana(gbk_input)

    print(f"[INFO] Descargando/recuperando ortologos...")
    paths = descargar_ortologos()

    print(f"[INFO] Combinando secuencias...")
    n = combinar_secuencias(query_fasta, paths)
    print(f"[OK]   Archivo combinado: {n} secuencias en {combined_fasta.name}")

    print(f"[INFO] Corriendo Clustal Omega...")
    correr_clustalo()

    print(f"[OK]   MSA completado:")
    print(f"       - {clustal_aln.relative_to(proyecto_dir)}")
    print(f"       - {fasta_aligned.relative_to(proyecto_dir)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
