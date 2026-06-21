# Ejercicio 3 — Multiple Sequence Alignment (MSA)

## 1. Objetivo

Realizar un alineamiento múltiple entre la secuencia proteica de **Huntingtina (HTT)**
humana y los ortólogos detectados en el reporte BLAST del Ejercicio 2, para evaluar
la conservación evolutiva de la proteína y caracterizar sus regiones variables.

## 2. Selección de secuencias

A partir de los hits significativos del BLAST local contra SwissProt (Ej2), se
seleccionaron tres ortólogos representativos de distintas distancias evolutivas
respecto del humano. Todos provienen de SwissProt (UniProt curado manualmente):

| Especie | UniProt ID | Identificador | Largo (aa) |
|---|---|---|---|
| *Homo sapiens* (query) | P42858 | HD_HUMAN | 3 144 |
| *Mus musculus* (ratón) | P42859 | HD_MOUSE | 3 119 |
| *Rattus norvegicus* (rata) | P51111 | HD_RAT | 3 110 |
| *Takifugu rubripes* (pez globo japonés) | P51112 | HD_TAKRU | 3 148 |

La query humana se obtuvo extrayendo la traducción anotada del **CDS** del archivo
`HTT_mRNA.gbk` ya utilizado en el Ej1 (la misma proteína cuyo BLAST en el Ej2
confirmó el frame correcto). Los ortólogos se descargaron automáticamente desde
la API REST de UniProt (`https://rest.uniprot.org/uniprotkb/{ID}.fasta`).

## 3. Metodología y herramientas

Se desarrolló el script `src/Ex3.py` en Python (Biopython) que automatiza:

1. **Extracción de la query humana** desde el CDS del GenBank (`HTT_mRNA.gbk`).
2. **Descarga de ortólogos** desde la API REST de UniProt (con caché local en
   `data/orthologs/` para evitar re-descargas).
3. **Combinación** de las 4 secuencias en `results/HTT_combined.fasta`.
4. **Alineamiento múltiple** con **Clustal Omega 1.2.4** (instalado vía
   `apt install clustalo`), invocado mediante `subprocess`.
5. **Generación de dos formatos de salida**:
   - `results/HTT_msa.aln` — formato Clustal (legible, con marcadores de conservación).
   - `results/HTT_msa.fasta` — formato FASTA alineado (apto para análisis posteriores).

### Comando de Clustal Omega

```bash
clustalo -i results/HTT_combined.fasta -o results/HTT_msa.aln --outfmt clu --force -v
```

### Instrucciones de ejecución

```bash
# Pre-requisitos: clustalo instalado y biopython en el venv
sudo apt install -y clustalo
pip install -r requirements.txt

# Correr el ejercicio
python src/Ex3.py
```

## 4. Resultados

### 4.1. Estadísticas generales del alineamiento

| Métrica | Valor |
|---|---|
| Número de secuencias | 4 |
| Longitud del alineamiento | 3 226 columnas |
| Columnas 100 % conservadas (sin gaps) | **2 113 (65,5 %)** |
| Columnas con al menos un gap | 171 (5,3 %) |

### 4.2. Identidad por pares respecto a HD_HUMAN

| Comparación | Identidad |
|---|---|
| HD_HUMAN vs HD_MOUSE  | **91,1 %** |
| HD_HUMAN vs HD_RAT    | **90,8 %** |
| HD_HUMAN vs HD_TAKRU  | **72,5 %** |

## 5. Interpretación

### 5.1. Conservación general

La proteína Huntingtina presenta una **muy alta conservación** entre vertebrados.
Los porcentajes de identidad por pares siguen el árbol filogenético esperado:
- Humano ↔ ratón / rata: ≈ 91 %, coherente con la divergencia mamífero–mamífero
  (~90 millones de años).
- Humano ↔ fugu (vertebrado no-mamífero): ≈ 73 %, coherente con la divergencia
  mucho más antigua (~450 millones de años).

Que el **65 % de las columnas estén perfectamente conservadas** (sin sustituciones
ni gaps) en una proteína de más de 3 100 aa es un indicador inequívoco de que
la mayor parte de la cadena cumple funciones esenciales sometidas a fuerte
**presión selectiva purificadora**.

### 5.2. Región Poly-Q: la "región de baja complejidad" predicha por Mount

El primer bloque del alineamiento (Clustal `HTT_msa.aln`) muestra de forma
inmediata el rasgo más característico de HTT — la **expansión de poliglutaminas
(Poly-Q)** cerca del extremo N-terminal:

```
HD_HUMAN  MATLEKLMKAFESLKSF QQQQQQQQQQQQQQQQQQQQQQQ PPPPPPPPPPP...
HD_MOUSE  MATLEKLMKAFESLKSF QQQQQQQP--------------- PPQAPPPPPPP...
HD_RAT    -------MKAFESLKSF QQQQQQQQ--------------- PPPQPPPPPPP...
HD_TAKRU  MATMEKLMKAFESLKSF QQQQG------------------ -----------...
```

| Especie | Glutaminas consecutivas (Poly-Q) |
|---|---|
| Humano   | 23 Q |
| Ratón    | 7 Q |
| Rata     | 8 Q |
| Fugu     | 4 Q |

Este resultado es justamente el fenómeno descripto por **David Mount (cap. 4)**
sobre las **regiones de baja complejidad**: aunque el resto de la proteína se
conserva fuertemente (los **flanqueos** `MATLEKLMKAFESLKSF` y `PPPPPP…` aparecen
prácticamente idénticos en mamíferos), la **longitud del tracto repetitivo es
altamente variable entre especies**.

Este rasgo es de relevancia clínica directa: en humanos, la **expansión patológica
del tracto Poly-Q por encima de ~36-40 repeticiones** es justamente la mutación
responsable de la enfermedad de Huntington. El alineamiento muestra que el
"baseline" funcional es notoriamente menor (4-23 Q según la especie), confirmando
que la patología no está en la presencia del tracto sino en su elongación
descontrolada.

### 5.3. Región C-terminal y dominios HEAT

El resto del alineamiento (no mostrado por brevedad) presenta **bloques largos
de identidad casi total** intercalados con regiones más variables — patrón
consistente con la presencia de los repetidos **HEAT** (Huntingtin–Elongation
factor 3–PP2A–TOR1) que conforman la estructura modular alfa-solenoide de HTT.
Estos dominios son fundamentales para sus funciones de scaffold en transporte
vesicular y desarrollo neuronal, lo cual justifica biológicamente la fuerte
presión selectiva observada.

## 6. Conclusión

El alineamiento múltiple confirma que **HTT es una proteína evolutivamente muy
conservada** (>90 % identidad entre mamíferos, ~73 % con un vertebrado lejano),
en consistencia con su rol esencial en el desarrollo neuronal. La **única
región notoriamente variable** entre especies es la **expansión Poly-Q** del
N-terminal, que en humanos es además el sustrato de la mutación causante de
la enfermedad de Huntington. El resultado del MSA respalda y profundiza la
interpretación biológica iniciada con el reporte BLAST del Ej2.

## 7. Archivos generados

```
data/orthologs/P42859.fasta        # HD_MOUSE (descargado de UniProt)
data/orthologs/P51111.fasta        # HD_RAT
data/orthologs/P51112.fasta        # HD_TAKRU
results/HTT_combined.fasta         # Input combinado (4 secuencias)
results/HTT_msa.aln                # MSA en formato Clustal (legible)
results/HTT_msa.fasta              # MSA en formato FASTA (para parseo)
src/Ex3.py                         # Script del ejercicio
```
