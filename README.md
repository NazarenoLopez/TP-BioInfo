# TP-BioInfo — Trabajo Práctico de Bioinformática

Trabajo práctico de la materia Bioinformática (UTN). Análisis del gen **HTT**
(Huntingtina, NM_002111), asociado a la enfermedad de Huntington.

## Estructura

```
src/
  Ex1.py                  # GenBank → 6 ORFs FASTA + detección de frame correcto
  Ex2.py                  # BLAST local contra SwissProt
  Ex3.py                  # MSA con Clustal Omega contra ortólogos
  Ex4.py                  # Parseo del reporte BLAST con filtro por patrón
data/
  orthologs/              # Ortólogos descargados de UniProt (cacheados)
results/                  # Outputs generados (MSA, alineamientos)
tests/                    # Unit tests
docs/                     # Documentación de las entregas
HTT_mRNA.gbk              # mRNA humano (NM_002111) en formato GenBank
HTT_ORFs.fas              # 6 marcos de lectura traducidos (output de Ex1)
blast_local_report.xml    # Reporte BLAST local (output de Ex2)
```

## Requisitos

- Python 3.10+
- BLAST+ (`blastp`, `makeblastdb`)
- Clustal Omega (`clustalo`)

En Ubuntu / WSL:

```bash
sudo apt install -y ncbi-blast+ clustalo python3-venv
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

En Windows (PowerShell):

```powershell
python -m venv .venv
Set-ExecutionPolicy -ExecutionPolicy RemoteSigned -Scope CurrentUser
.\.venv\Scripts\Activate.ps1
pip install -r requirements.txt
```

> **Red con proxy SSL (UTN u otras redes corporativas):** si `pip install` falla con
> `CERTIFICATE_VERIFY_FAILED`, creá el archivo `%APPDATA%\pip\pip.ini` con el
> siguiente contenido para que pip confíe en los hosts de PyPI:
>
> ```ini
> [global]
> trusted-host = pypi.org
>     files.pythonhosted.org
>     pypi.python.org
> ```

## Ejecución

```bash
# Ej1 — generar 6 ORFs y detectar el frame correcto
python src/Ex1.py

# Ej2 — BLAST local (requiere la base de datos SwissProt en ./Swissport/)
python src/Ex2.py

# Ej3 — MSA con ortólogos descargados de SwissProt
python src/Ex3.py

# Ej4 — parsear el reporte BLAST filtrando por patrón
#   Busca en blast_local_report.xml todos los hits cuya descripción contenga
#   el patrón indicado (búsqueda case-insensitive) e imprime accesión,
#   descripción, E-value, bit score e identidad de cada coincidencia.
python src/Ex4.py "Mus musculus"
python src/Ex4.py "Huntingtin"

#   Con --fetch descarga además las secuencias completas desde NCBI (requiere
#   conexión a internet) y las guarda en matching_hits.fas
python src/Ex4.py "Huntingtin" --fetch
```

## Tests

```bash
python -m unittest discover -s tests -p "test_*.py"
```
