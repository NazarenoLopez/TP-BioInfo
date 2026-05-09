from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastpCommandline
import os

# 1. Definición de rutas relativas
# Obtener el directorio del script actual
script_dir = os.path.dirname(os.path.abspath(__file__))
proyecto_dir = os.path.dirname(script_dir)  # Subir un nivel al directorio raíz

# El archivo HTT_ORFs está en la raíz del proyecto
archivo_query = os.path.join(proyecto_dir, "HTT_ORFs.fas")
# La base de datos está dentro de la subcarpeta Swissport
ruta_db = os.path.join(proyecto_dir, "Swissport", "swissprot")
archivo_output = "blast_local_report.xml"

def ejecutar_blast_local():
    print(f"Iniciando BLAST local contra la base de datos en: {ruta_db}...")
    
    # Verificamos que el archivo input existe
    if not os.path.exists(archivo_query):
        print(f"Error: No se encuentra el archivo {archivo_query}")
        return

    # 2. Configuración del comando blastp
    # outfmt=5 genera un XML compatible con los parsers de BioPython
    blastp_cline = NcbiblastpCommandline(
        query=archivo_query, 
        db=ruta_db, 
        outfmt=5, 
        out=archivo_output
    )

    # 3. Ejecución del comando
    try:
        stdout, stderr = blastp_cline()
        print(f"Éxito: El reporte se ha generado en '{archivo_output}'")
    except Exception as e:
        print(f"Ocurrió un error al ejecutar BLAST local: {e}")

if __name__ == "__main__":
    ejecutar_blast_local()