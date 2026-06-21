# TP de Bioinformática:

## Setup inicial
Para poder desarrollar el trabajo práctico deberán tener instalados el sistema operativo Linux (opcional), el lenguaje de programación de elección (Perl, Python o Java entre otros) con las librerías correspondientes (de BioPerl o de BioPython, …) y la suite de programas de Blast.

### Linux
Pueden realizar una partición del disco o bien instalar una máquina virtual.
Las distribuciones Linux que suelen utilizarse en bioinformática son Debian o Ubuntu. Existen distribuciones como Bio-Linux cuyos paquetes pueden ser instalados sobre Debian o Ubuntu, o la máquina virtual de DNALinux, ambas ya tienen Perl preinstalado, EMBOSS y otras herramientas que les serán de utilidad ya incluidas.

### Perl (lenguaje recomendado, aunque otros lenguajes también se pueden)
Los programas Perl son llamados scripts y tienen extensión *.pl (o bien *.cgi si son aplicaciones web). Perl es una lenguaje interpretado o de scripting (aunque también hay compiladores) cuya estructura deriva del C y toma cosas de la programación shell. Instalación y documentación en:
http://www.perl.org

### BioPerl
BioPerl es proyecto comunitario open source de módulos Perl integrados para trabajar con secuencias y anotaciones, acceder a bases de datos remotas, parsear el output de programas como BLAST, FASTA, etc. Es prácticamente esencial para todo bioinformático.
BioPerl-core tiene los módulos principales, BioPerl-run es una colección de módulos que facilitan la ejecución de programas locales como EMBOSS suite. Instalar BioPerl desde
http://www.bioperl.org


### BLAST local (BLAST+ is a new suite of BLAST tools that utilizes the NCBI)
La descarga de los programas blast la hacen desde:
https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html#downloadblastdata
Existen distintas maneras de correr Blast de manera local.
Se puede ejecutar con algún comando específico of bioperl o bien desde línea de comando:

```
> blastp -d swissprot -i demo.fasta -o myblast.report (con standalone Blast)
> blastall -p blastp -d swissprot -i demo.fasta -o myblast.report (con Blast+)
```

Desde BioPerl:

Utilizando el objeto Bio::Tools::Run::StandAloneBlast.pm (con Blast instalado local)
Utilizando el objeto Bio::Tools::Run::RemoteBlast.pm (usando el Blast del NCBI)
BLAST help: http://www.ncbi.nlm.nih.gov/books/NBK52637/2


# Trabajo Práctico (parte 1)
El presente trabajo práctico tiene por objetivo adquirir las primeras habilidades en el campo de la Bioinformática. Se incluyen cuatro ejercicios donde deberán desarrollar pequeños scripts para resolver problemas específicos. Los mismos pueden ser desarrollados utilizando cualquiera de los
lenguajes de programación bioinformática de código abierto como BioPerl, BioJava y BioRuby, que son ampliamente utilizados en la investigación bioinformática y de biología computacional, aunque se sugiere la utilización de BioPerl para facilitar la resolución de los ejercicios. Las herramientas
computacionales escritas en estos lenguajes proporcionan múltiples  uncionalidades para crear soluciones personalizadas y realizar análisis de datos biológicos. Un quinto ejercicio esta relacionado con la comprensión de la información en bases de datos de biología molecular.
Para comenzar el trabajo deben entrar en la base de datos Online Mendelian Inheritance in Man (OMIM) donde encontrarán el catálogo online genes humanos asociados a trastornos genéticos más importante de la actualidad. En grupo decidan sobre que enfermedad hereditaria quieren investigar
y luego a partir de la información en OMIM seleccionen un gen asociado a esta patología para comenzar con el ejercicio 1. La secuencia de este mismo gen/proteína debe utilizarse luego en el
ejercicio 5 (en la parte 2 del TP).
Cuando entreguen el TP completo (ambas partes) cada grupo tendrá 10 minutos para exponer como
realizó el trabajo práctico y comentar sobre su investigación (no tanto sobre el código realizado).
Por favor preparen una presentación. La correcta exposición del trabajo realizado por los miembros
del grupo también entra en la evaluación.

### Ejercicio 1 – PROCESAMIENTO DE SECUENCIAS. 
Escribir un script que lea una o más secuencias (de nucleótidos) de un archivo que contenga la información en formato GenBank de un mRNA de referencia (NM_xxxx) de su gen (o genes) de interés, las traduzca a sus secuencias de amino ácidos posibles (tener en cuenta los 6 Reading Frames posibles) y escriba los resultados en un archivo en formato FASTA. Ustedes deben generarse su archivo GenBank de secuencias input, por ejemplo realizando una consulta de los mRNA del gen INS (que está asociado a la Diabetes) en la base de datos de NCBI-Gene y obtener uno o más resultados en formato GenBank en un archivo de texto. Si no desean seguir trabajando con las seis secuencias
de aa posibles, pueden utilizar alguna función o programa que les permita saber cual el es marco de lectura correcto (ORF, Open Reading Frames) y seguir con esa secuencia. Informar el cual de los 6 Reading frames
está codificada la secuencia.

_NOTA: Ver aclaración de este ejercicio al final del documento._

- Input: Archivo de secuencias Genbank (ej. NMxxxx.gbk con una o más secuencias).
- Output: Archivo de secuencias Fasta de cada ORF (ej. Xxxxx.fas con una o más secuencias de aminoácidos). Indicar cuál es el marco de lectura correcto.

Deben entregar el script Ex1.pm (si lo hacen con BioPerl, sino será otra extensión) y el input file que utilicen
con una breve descripción de lo que hicieron y como se debe ejecutar para probarlo.

### Ejercicio 2.a - BLAST. 
Escribir un script que realice un BLAST de una o varias secuencias (si son varias se realiza
un Blast por cada secuencia input) y escriba el resultado (blast output) en un archivo. Nota: Pueden ejecutar
BLAST de manera remota o bien localmente (¡si hacen ambos tienen más puntos!), para esto deben instalarse
BLAST localmente del FTP del NCBI, luego bajarse la base de datos:
ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/swissprot.gz
y descomprimirla en un dir por ej. ncbi-blast-2.3.0+/data/, luego usar el comando ncbi-blast-
2.3.0+/bin/makeblastdb sobre el archivo swissprot (el original ya está en formato FASTA) para darle formato
de BLAST DB. Dependiendo de la versión de Blast suite que tengan instalado puede que en vez de
makeblastdb deban utilizar el comando formatdb.
- Input: Secuencia Fasta (ej. Xxxx.fas con una o más secuencias de aminoácidos obtenidas en Ej.1).
- Output: Reporte Blast (ej. blast.out, si deciden hacer múltiples pueden generar un único o varios
archivos).

Deben entregar el script Ex2.pm y su input file con una breve descripción de lo que hicieron, con una
interpretación de los resultados del Blast, y mencionar como se debe ejecutar para probarlo.

### Ejercicio 2.b – Interpretación del resultado del Blast. 

Dar una explicación del resultado blast obtenido en
términos de las secuencias encontradas y dar una explicación sobre que significan los valores estadísticos
asociados a las secuencias encontradas (el capítulo 4 del libro de David Mount puede ayudarlos).

### Ejercicio 3 – Multiple Sequence Alignment (MSA). 
Descargar las secuencias (en formato fasta) de 3 o más
organismos distintos pertenecientes a otras especies que hayan salido en los resultados del Blast y realizar
un alineamiento múltiple con tu secuencia de consulta más estas otras encontradas. Si no pueden hacerlo
localmente pueden utilizar algún programa de MSA online. Intenten realizar una interpretación del resultado
del alineamiento múltiple. Entregar información del MSA.

### Aclaración para el Ejercicio 1:
Para bajar una secuencia de nuestro gen elegido que funcione para en el Ejercicio 1 y para los demás, deberán bajarse alguno de los RNA mensajeros maduros (transcripto) de su gen de interés. Es decir, una secuencia de
mRNA que ya haya sido procesada y no tenga intrones, esta es la secuencia que deben bajarse en formato Genbank y hacer la traducción a su secuencia de aminoácidos. Por ejemplo, para el gen de la insulina humano (INS Homo sapiens):
1. Hacer una búsqueda en la base de datos de Genes e ir a las secuencias de Referencia y seleccionar algunos de los mRNA (NMxxxx)
2. Seleccionar uno de los transcriptos del gen en formato GenBank (en lo posible la isoforma 1)
3. ORF (Open Reading Frame)
Una vez que tienen la secuencia bajada tengan en cuenta que ustedes  esconocen cuál es el marco de lectura correcto de los 6 posibles. Por lo tanto deberán calcular los 6 marcos de lectura posibles, evaluar todos ellos en el ejercicio 2, y así darse cuenta cuál de los 6 es el real. Existen funciones en BioPerl para hacer esto, o mismo pueden usar el programa OrfFinder para ayudarse.


# PARTE 2 DEL TP

### Ejercicio 4 – BLAST OUTPUT. 
Escribir un script para analizar (parsear) un reporte de salida de blast que identifique los hits que en su descripción aparezca un Pattern determinado que le damos como parámetro de entrada. El pattern puede ser una palabra. Punto extra: pueden a su vez parsear cuál es el ACCESSION del hit identificado (donde hay una coincidencia del Pattern) y con el módulo Bio::DB::GenBank obtener la secuencia completa del hit en formato FASTA y escribirla a un archivo, es decir, levantar las secuencias originales completas de los hits seleccionados.
− Input: Reporte Blast (blast.out del ej. 2) y un Pattern (por ej. “Mus Musculus”).
− Output: Lista de los hits que coincidan con el pattern (por ej. solo los hits de Ratones).
Deben entregar el script Ex4.pm y su input file con una breve descripción.

### Ejercicio 5 - EMBOSS. 
Instalar EMBOSS. Escribir un script que llame a uno o más programas EMBOSS para hacer algún análisis sobre la secuencia de nucleótidos del mRNA y/o sobre la secuencia de amino ácidos de la proteína investigada.
Por ejemplo, pueden correr un programa que calcule los ORFs y obtenga las secuencias de proteínas posibles o algún otro programa EMBOSS que les resulte de interés.
Luego deben bájense la base de datos PROSITE de dominios/motivos funcionales conocidos (archivo prosite.dat), y por medio del llamado a otro programa EMBOSS realizar un análisis de dominios de las secuencias de aminoácidos obtenidas y escribir los resultados en un archivo de salida.
− Input: Archivo de secuencias Fasta (por ej. Xxxxx.fas con una o más secuencias de aa).
− Output: Archivo de resultados de dominios funcionales encontrados en las secuencias de aa.

### Ejercicio 6. Trabajo con Bases de Datos Biológicas.
a) A partir del gen o proteína de interés para ustedes dar su link a NCBI-Gene como una entrada de Entrez,
por ej.: http://www.ncbi.nlm.nih.gov/gene/3630
Expliquen brevemente lo que hace la proteína y por qué la eligieron.
b) ¿Cuántos genes / proteínas homólogas se conocen en otros organismos? Utilicen la información que está en la base de datos de HomoloGene y en las bases de datos Ensembl. Describan los resultados en ambas bases de datos, y en qué se diferencian. Mencionen sobre qué tan común creen son estos genes o proteínas y a qué grupos taxonómicos pertenecen (sólo en las bacterias, en los vertebrados, etc.)
c) ¿Cuántos transcriptos y cuántas formas alternativas de splicing son conocidos para este gen / proteína?
¿Cuáles de estos splicing alternativos se expresan? ¿Tienen funciones alternativas? Buscar evidencia de
esto en las bases de datos de NCBI y en los transcriptos de Ensembl ¿Cómo el número de splicings alternativos diferente entre las dos bases de datos y cuál piensan que es más precisa y por qué?
d) ¿Con cuántas otras proteínas interactúa el producto génico de su gen? ¿Existe un patrón o relación entre las interacciones? Mencione las interacciones interesantes o inusuales. Usted encontrará las interacciones de su gene/proteína tanto en la base de datos NCBI Gene como en la base de datos UniProt. Compare las
dos tablas entre sí. ¿Hay proteínas que interactúan únicas para cada tabla?
e) Expliquen brevemente de qué componente celular forma parte su proteína (pista: se puede estudiar la información de Gene Ontology - GO), ¿A qué procesos biológicos pertenece (pista idem)? y ¿En qué función molecular trabaja esta proteína? Los términos ontológicos de genes los pueden encontrar tanto en NCBI
Gene y en la base de datos UniProt como haciendo una búsqueda en AmiGO.
f) Discutan brevemente en qué estructura o vías metabólicas específicas (pathways) estaría participando su gen / proteína? (Reactome, KEGG son algunas bases de datos de pathways).
g) Entrar en la base de datos de variantes genéticas dbSNP e intentar interpretar o encontrar info sobre alguna variante (reference SNP - rsXXXX) asociada con la patología investigada en su gen de interés. ¿Qué variante es? ¿Hay información sobre la frecuencia que tiene esta variante en la población? ¿Qué grupo étnico parece ser el más afectado?

NOTA: Para hacer este ejercicio les pueden servir algunas otras bases de datos como estas (entre otras):

http://www.genecards.org
https://www.ncbi.nlm.nih.gov/snp/ (para obtener información de la variante en la población)
http://www.ncbi.nlm.nih.gov/clinvar/ (para obtener información clínica del gen y sus variantes)
https://ghr.nlm.nih.gov

### Ejercicio 7. Armar una presentación donde expliquen la enfermedad que investigaron, lo que hicieron y los resultados que fueron obteniendo en los ejercicios del TP.
Los integrantes de cada grupo tendrán un máximo de 10 minutos (~15 diapositivas) para exponer el trabajo
práctico. Comentar sobre sus investigaciones: Intro - Métodos (código implementado, muy brevemente) -
Resultados. La correcta presentación del trabajo realizado es también parte importante de la evaluación
(solo enviar el pdf de la presentación).