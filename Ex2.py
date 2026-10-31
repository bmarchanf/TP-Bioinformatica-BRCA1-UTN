#!/usr/bin/env python3
"""
Trabajo Práctico de Bioinformática - Curso K5561 (UTN.BA)
Ejercicio 2: Alineamiento de Secuencias con BLAST (Remoto y Local)
Gen de estudio: BRCA1 (Homo sapiens, 1863 aminoácidos)
"""

import os
import re
import subprocess
from Bio import SeqIO
from Bio.Blast import NCBIWWW


def ejecutar_blast_remoto(fasta_file="brca1_orf.fasta", xml_out="blast_brca1.xml", db="swissprot", hits=15):
    record = SeqIO.read(fasta_file, "fasta")
    print(f"[BLAST Remoto] Consultando NCBI BLASTp ({db}) para {record.id} ({len(record.seq)} aa)...")
    result_handle = NCBIWWW.qblast("blastp", db, str(record.seq), hitlist_size=hits)
    with open(xml_out, "w", encoding="utf-8") as out_f:
        out_f.write(result_handle.read())
    result_handle.close()
    print(f"[BLAST Remoto] Resultado almacenado en '{xml_out}'.")


def ejecutar_blast_local(fasta_file="brca1_orf.fasta", out_file="blast_brca1.out", db="swissprot", hits=20):
    print(f"[BLAST Local] Ejecutando blastp contra la base de datos local '{db}'...")
    cmd = ["blastp", "-query", fasta_file, "-db", db, "-out", out_file, "-max_target_seqs", str(hits)]
    subprocess.run(cmd, check=True)


def resumir_resultados_blast(blast_out_file="blast_brca1.out", summary_file="blast_summary_table.txt"):
    with open(blast_out_file, "r", encoding="utf-8") as f:
        text = f.read()

    bloques = text.split("\n>")[1:13]
    lineas = [
        "=" * 102,
        "EJERCICIO 2: RESUMEN DE HITS SIGNIFICATIVOS - BLASTP (Base de datos: SwissProt)",
        "=" * 102,
        f"{'#':<3} | {'Acceso / Proteína y Organismo':<46} | {'Bit Score':<9} | {'Raw Score':<9} | {'E-value':<8} | {'Identidad':<10} | {'Positivos':<9}",
        "-" * 102,
    ]

    for idx, bloque in enumerate(bloques, start=1):
        encabezado = " ".join(bloque.splitlines()[:2]).replace("  ", " ")
        m_score = re.search(r"Score = ([\d.]+) bits \((\d+)\),\s*Expect = ([\de.-]+)", bloque)
        m_id = re.search(
            r"Identities = \d+/\d+ \((\d+%)\),\s*Positives = \d+/\d+ \((\d+%)\),\s*Gaps = (\d+/\d+ \(\d+%\))",
            bloque,
        )
        if m_score and m_id:
            lineas.append(
                f"{idx:<3} | {encabezado[:46]:<46} | {m_score.group(1):<9} | {m_score.group(2):<9} | "
                f"{m_score.group(3):<8} | {m_id.group(1):<10} | {m_id.group(2):<9}"
            )

    lineas.append("-" * 102)
    salida = "\n".join(lineas) + "\n"
    print(salida)
    with open(summary_file, "w", encoding="utf-8") as f:
        f.write(salida)


if __name__ == "__main__":
    if not os.path.exists("blast_brca1.out"):
        ejecutar_blast_remoto()
    resumir_resultados_blast()
