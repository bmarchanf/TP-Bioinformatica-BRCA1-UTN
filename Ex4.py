#!/usr/bin/env python3
"""
Trabajo Práctico de Bioinformática - Curso K5561 (UTN.BA)
Ejercicio 4: Parser de Salida BLAST por Patrón y Búsqueda de Motivos Funcionales (PROSITE)
Gen de estudio: BRCA1 (Homo sapiens, P38398 / NP_009225.1)
"""

import os
import re
from Bio import SeqIO


def filtrar_blast_por_patron(blast_out_file="blast_brca1.out", patron="Pan troglodytes"):
    with open(blast_out_file, "r", encoding="utf-8") as f:
        text = f.read()

    out_fasta = f"blast_match_{patron.replace(' ', '_')}.fasta"
    coincidencias = []

    print("=" * 88)
    print(f"EJERCICIO 4A: PARSER DE SALIDA BLAST (Patrón buscado: '{patron}')")
    print("=" * 88)

    for bloque in text.split("\n>")[1:]:
        encabezado = " ".join(bloque.splitlines()[:2]).strip()
        if patron.lower() in encabezado.lower():
            m_score = re.search(r"Score = ([\d.]+) bits \((\d+)\),\s*Expect = ([\de.-]+)", bloque)
            sbjct_partes = re.findall(r"Sbjct\s+\d+\s+([A-Z-]+)\s+\d+", bloque)
            seq_completa = "".join(sbjct_partes).replace("-", "")
            bits, raw, ev = m_score.groups() if m_score else ("N/A", "N/A", "N/A")
            coincidencias.append((encabezado, raw, bits, ev, seq_completa))
            print(f"[MATCH] {encabezado[:75]}")
            print(f"        Tipo: Proteína (BLASTp) | Score: {raw} ({bits} bits) | E-value: {ev} | Longitud: {len(seq_completa)} aa")

    with open(out_fasta, "w", encoding="utf-8") as f:
        for titulo, _, _, _, seq in coincidencias:
            f.write(f">{titulo}\n{seq}\n")
    print(f"[ARCHIVO] Secuencias exportadas en '{out_fasta}' ({len(coincidencias)} hits).\n")


def escanear_motivos_prosite(fasta_file="brca1_orf.fasta", reporte_motivos="prosite_brca1_motifs.txt"):
    rec = SeqIO.read(fasta_file, "fasta")
    seq = str(rec.seq)

    motivos_prosite = [
        ("PS00518", "ZF_RING_1 (Zinc Finger RING-type signature)", r"C.H.[LIVMFY]C.{2}C[LIVMYA]",
         "Dominio E3 ubiquitina ligasa N-terminal; media la heterodimerización con BARD1."),
        ("PS50172", "BRCT (Dominio BRCA1 C-Terminal de unión a fosfopéptidos)", r"VNKRMSMVVSGLTPEEFMLVYKFARKHHITLTNLITEETTHVVMKTDAEFVCERTLKYFL",
         "Dominio C-terminal de reconocimiento de proteínas fosforiladas por ATM/ATR (ABRAXAS1, BRIP1, CtIP)."),
        ("PS00015", "NLS_BP (Señal de Localización Nuclear)", r"K[RK]K[RK][RK]P",
         "Señal de importación nuclear activa de BRCA1 hacia el nucleoplasma."),
        ("PS00001", "ASN_GLYCOSYLATION (Sitio de N-glicosilación)", r"N[^P][ST][^P]",
         "Sitio consenso de modificación postraduccional por N-glicosilación."),
    ]

    lineas = [
        "=" * 88,
        "EJERCICIO 4B: ANÁLISIS DE DOMINIOS Y MOTIVOS PROSITE EN BRCA1",
        "=" * 88,
        f"Secuencia analizada: {rec.id} ({len(seq)} aminoácidos)\n",
    ]

    for codigo, nombre, regex, funcion in motivos_prosite:
        matches = list(re.finditer(regex, seq))
        lineas.append(f"Motivo PROSITE: {codigo} - {nombre}")
        lineas.append(f"Función Biológica: {funcion}")
        lineas.append(f"Ocurrencias detectadas: {len(matches)}")
        for m in matches[:4]:
            lineas.append(f"  -> Posición {m.start() + 1}..{m.end()}: '{m.group(0)}'")
        lineas.append("-" * 88)

    salida = "\n".join(lineas) + "\n"
    print(salida)
    with open(reporte_motivos, "w", encoding="utf-8") as f:
        f.write(salida)


if __name__ == "__main__":
    if os.path.exists("blast_brca1.out"):
        filtrar_blast_por_patron("blast_brca1.out", "Pan troglodytes")
    escanear_motivos_prosite()
