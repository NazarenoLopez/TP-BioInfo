from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def procesar_secuencia_y_detectar_orf(archivo_entrada, archivo_salida):
    # 1. Cargar el registro GenBank
    record = SeqIO.read(archivo_entrada, "genbank")
    secuencia_dna = record.seq
    
    # 2. Obtener la proteína real desde la anotación CDS del archivo para comparar
    proteina_real = ""
    for feature in record.features:
        if feature.type == "CDS":
            proteina_real = str(feature.qualifiers.get('translation', [''])[0])
            break

    marcos_traducidos = []
    marco_correcto_id = "No detectado"

    # 3. Generar y comparar los 6 marcos
    # Directos (1, 2, 3) y Reversos (4, 5, 6)
    for i in range(3):
        # Sentido Directo
        dna_fwd = secuencia_dna[i:]
        trim_fwd = dna_fwd[:len(dna_fwd) - (len(dna_fwd) % 3)]
        prot_fwd = str(trim_fwd.translate())
        
        id_fwd = f"Frame_{i+1}_Fwd"
        marcos_traducidos.append(SeqRecord(trim_fwd.translate(), id=id_fwd, description=""))
        
        # Sentido Reverso
        dna_rev = secuencia_dna.reverse_complement()[i:]
        trim_rev = dna_rev[:len(dna_rev) - (len(dna_rev) % 3)]
        prot_rev = str(trim_rev.translate())
        
        id_rev = f"Frame_{i+1}_Rev"
        marcos_traducidos.append(SeqRecord(trim_rev.translate(), id=id_rev, description=""))

        # Lógica de detección: ¿Esta traducción contiene la proteína real?
        if proteina_real and proteina_real in prot_fwd:
            marco_correcto_id = id_fwd
        elif proteina_real and proteina_real in prot_rev:
            marco_correcto_id = id_rev

    # 4. Guardar el archivo FASTA
    SeqIO.write(marcos_traducidos, archivo_salida, "fasta")
    
    # 5. Informar resultados
    print(f"--- Proceso Completado ---")
    print(f"Archivo generado: {archivo_salida}")
    print(f"El marco de lectura correcto es: {marco_correcto_id}")
    if marco_correcto_id != "No detectado":
        print(f"Razón: Coincide con la secuencia de proteína anotada en el CDS del GenBank.")

# Ejecución
procesar_secuencia_y_detectar_orf("HTT_mRNA.gbk", "HTT_ORFs.fas")