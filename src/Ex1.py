from pathlib import Path

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def procesar_secuencia_y_detectar_orf(archivo_entrada, archivo_salida):
    # 1. Cargar el registro GenBank
    record = SeqIO.read(archivo_entrada, "genbank")
    secuencia_dna = record.seq

    # 2. Obtener la proteína real desde la anotación CDS para comparar
    proteina_real = ""
    for feature in record.features:
        if feature.type == "CDS":
            proteina_real = str(feature.qualifiers.get('translation', [''])[0])
            break

    marcos_traducidos = []
    marco_correcto_id = "No detectado"

    # 3. Generar y comparar los 6 marcos: directos (1, 2, 3) y reversos (4, 5, 6)
    for i in range(3):
        # Sentido directo
        dna_fwd = secuencia_dna[i:]
        trim_fwd = dna_fwd[: len(dna_fwd) - (len(dna_fwd) % 3)]
        traduccion_fwd = trim_fwd.translate()
        id_fwd = f"Frame_{i+1}_Fwd"
        marcos_traducidos.append(
            SeqRecord(traduccion_fwd, id=id_fwd, description="")
        )

        # Sentido reverso
        dna_rev = secuencia_dna.reverse_complement()[i:]
        trim_rev = dna_rev[: len(dna_rev) - (len(dna_rev) % 3)]
        traduccion_rev = trim_rev.translate()
        id_rev = f"Frame_{i+1}_Rev"
        marcos_traducidos.append(
            SeqRecord(traduccion_rev, id=id_rev, description="")
        )

        # ¿Alguna traducción contiene la proteína anotada?
        if proteina_real and proteina_real in str(traduccion_fwd):
            marco_correcto_id = id_fwd
        elif proteina_real and proteina_real in str(traduccion_rev):
            marco_correcto_id = id_rev

    # 4. Guardar el archivo FASTA
    SeqIO.write(marcos_traducidos, archivo_salida, "fasta")

    # 5. Informar resultados
    print(f"--- Proceso Completado ---")
    print(f"Archivo generado: {archivo_salida}")
    print(f"El marco de lectura correcto es: {marco_correcto_id}")
    if marco_correcto_id != "No detectado":
        print("Razón: Coincide con la secuencia de proteína anotada en el CDS del GenBank.")

    return marco_correcto_id


if __name__ == "__main__":
    proyecto_dir = Path(__file__).resolve().parent.parent
    procesar_secuencia_y_detectar_orf(
        str(proyecto_dir / "HTT_mRNA.gbk"),
        str(proyecto_dir / "HTT_ORFs.fas"),
    )
