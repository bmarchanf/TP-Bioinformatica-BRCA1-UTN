# Trabajo Práctico de Bioinformática — Cáncer de Mama y Ovario Hereditario (`BRCA1`)

**Universidad Tecnológica Nacional — Facultad Regional Buenos Aires (UTN.BA)**  
**Asignatura:** Bioinformática — Curso K5561 (2° Cuatrimestre 2026)  
**Docente:** Dr. Patricio Yankilevich  

## Integrantes del Equipo

| Apellido y Nombre | Correo Electrónico Institucional |
| :--- | :--- |
| **Enriquez, Sylvina** | `sylvina64@frba.utn.edu.ar` |
| **Kesel, Daniel** | `dkesel@frba.utn.edu.ar` |
| **Marchan Figuera, Brisabed** | `bmarchanfiguera@frba.utn.edu.ar` |
| **Vasquez Landaeta, Carlos** | `cvasquez@frba.utn.edu.ar` |

---

## Descripción del Proyecto

En el presente repositorio se desarrolla el procesamiento y análisis bioinformático del gen supresor de tumores humano **`BRCA1`** (*BRCA1 DNA repair associated*, OMIM `*113705`, NCBI Gene ID `672`, locus `17q21.31`), cuyas variantes germinales patogénicas se asocian al **Síndrome de susceptibilidad hereditaria al cáncer de mama y ovario tipo 1** (OMIM `#604370`).

### Herramientas Utilizadas
* **Python 3** y biblioteca **Biopython** (`Bio.SeqIO`, `Bio.Seq`, `Bio.Blast.NCBIWWW`, `Bio.Blast.NCBIXML`, `Bio.AlignIO`).
* **NCBI BLAST+ (`blastp`)** contra la base de datos curada **UniProtKB/SwissProt**.
* **Clustal Omega (EMBL-EBI)** para Alineamiento Múltiple de Secuencias (MSA).
* **PROSITE** para detección de dominios y motivos funcionales (`PS00518 ZF_RING_1`, `PS50172 BRCT`, `PS00015 NLS_BP`).

---

## Estructura del Repositorio

* **`Ex1.py`**: Procesamiento del transcripto maduro de referencia `NM_007294.4` (`NM_007294.gbk`, 7088 pb), traducción en los **6 marcos de lectura posibles (`+1, +2, +3, -1, -2, -3`)** y extracción del marco abierto de lectura (ORF) codificante real en el **Marco +3** (1863 aminoácidos, proteína `NP_009225.1` / `P38398.2`).
* **`Ex2.py`**: Ejecución de búsqueda de homología mediante **BLASTp** (remota vía `NCBIWWW` y local vía `blastp`) sobre la base de datos **SwissProt** y generación del cuadro resumen de métricas estadísticas (`blast_summary_table.txt` y `blast_brca1.out`).
* **`Ex3.py`**: Descarga de secuencias ortólogas de `BRCA1` en 8 especies de mamíferos (*Homo sapiens*, *Pan troglodytes*, *Gorilla gorilla*, *Pongo pygmaeus*, *Macaca mulatta*, *Canis lupus familiaris*, *Bos taurus* y *Mus musculus*), alineamiento múltiple con **Clustal Omega** y cálculo de la matriz de identidad porcentual pareada (`msa_identity_report.txt`).
* **`Ex4.py`**: Analizador sintáctico (*Blast Output Parser*) por patrón taxonómico y escaneo de motivos funcionales de la base de datos **PROSITE** sobre la proteína `BRCA1` (`prosite_brca1_motifs.txt`).

### Archivos de Entrada y Salida Incluidos
* `NM_007294.gbk` / `NM_007294.fasta`: Secuencia de referencia del ARNm de `BRCA1` en formatos GenBank y FASTA.
* `brca1_6frames.fasta`: Traducción completa de los 6 marcos de lectura.
* `brca1_orf.fasta`: Secuencia de aminoácidos del ORF principal (1863 aa).
* `blast_brca1.out` / `blast_brca1.xml` / `blast_summary_table.txt`: Resultados completos y tabla resumen de BLASTp.
* `msa_input.fasta` / `msa_output.aln` / `msa_output.fasta` / `msa_identity_report.txt`: Archivos del alineamiento múltiple de secuencias y matriz de identidad.
* `prosite_brca1_motifs.txt` / `blast_match_Pan_troglodytes.fasta`: Salidas del parser de BLAST y búsqueda de motivos PROSITE.

---

## Instrucciones de Ejecución

```bash
# 1. Instalar dependencias
pip install biopython

# 2. Ejercicio 1: Traducción de 6 marcos de lectura y obtención del ORF
python3 Ex1.py

# 3. Ejercicio 2: Búsqueda BLASTp y tabla de métricas estadísticas
python3 Ex2.py

# 4. Ejercicio 3: Alineamiento Múltiple de Secuencias (MSA) y matriz de identidad
python3 Ex3.py

# 5. Ejercicio 4: Parser de BLAST por patrón y motivos funcionales PROSITE
python3 Ex4.py
```
