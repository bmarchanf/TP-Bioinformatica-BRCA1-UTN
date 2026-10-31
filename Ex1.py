#!/usr/bin/env python3
"""
Trabajo Práctico de Bioinformática - Curso K5561 (UTN.BA)
Ejercicio 1: Procesamiento de Secuencias y Análisis de los 6 Marcos de Lectura (ORF)
Gen de estudio: BRCA1 (Homo sapiens, RefSeq mRNA NM_007294.4)
"""

import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def analizar_marcos_de_lectura(gbk_file="NM_007294.gbk", out_frames="brca1_6frames.fasta", out_orf="brca1_orf.fasta"):
    record = SeqIO.read(gbk_file, "genbank")
    mrna_seq = record.seq

    print("=" * 78)
    print("EJERCICIO 1: PROCESAMIENTO DE SECUENCIAS CON BIOPYTHON (BRCA1)")
    print("=" * 78)
    print(f"Accession  : {record.id}")
    print(f"Descripción: {record.description}")
    print(f"Longitud   : {len(mrna_seq)} pb\n")

    print("-" * 78)
    print(f"{'Marco':^8} | {'Hebra':^9} | {'Codones':^10} | {'Codones Stop (*)':^18} | {'ORF Máximo (aa)':^18}")
    print("-" * 78)

    registros_marcos = []
    mejor_orf = ("", 0, "")

    hebras = [("+", "Directa", mrna_seq), ("-", "Reversa", mrna_seq.reverse_complement())]

    for signo, nombre_hebra, sec_hebra in hebras:
        for offset in range(3):
            marco_id = f"{signo}{offset + 1}"
            sub_seq = sec_hebra[offset:]
            # Ajustar a múltiplo exacto de 3 nucleótidos
            sub_seq = sub_seq[:len(sub_seq) - (len(sub_seq) % 3)]
            traduccion = sub_seq.translate(to_stop=False)

            num_stops = str(traduccion).count("*")
            orfs = re.findall(r"M[^*]*", str(traduccion))
            orf_max = max(orfs, key=len) if orfs else ""

            print(f"{marco_id:^8} | {nombre_hebra:^9} | {len(traduccion):^10} | {num_stops:^18} | {len(orf_max):^18}")

            registros_marcos.append(
                SeqRecord(
                    traduccion,
                    id=f"{record.id}_marco_{marco_id}",
                    description=f"Traduccion completa marco {marco_id} ({len(traduccion)} aa, {num_stops} stops)"
                )
            )

            if len(orf_max) > mejor_orf[1]:
                mejor_orf = (marco_id, len(orf_max), orf_max)

    print("-" * 78)
    print(f"\n[RESULTADO] Marco de lectura codificante identificado: Marco {mejor_orf[0]}")
    print(f"[RESULTADO] Longitud de la proteína BRCA1 (ORF principal): {mejor_orf[1]} aminoácidos")
    print(f"[RESULTADO] Primeros 60 aa: {mejor_orf[2][:60]}...")

    # Guardar los 6 marcos traducidos y el ORF principal en formato FASTA
    SeqIO.write(registros_marcos, out_frames, "fasta")
    orf_record = SeqRecord(
        Seq(mejor_orf[2]),
        id=f"{record.id}_ORF_{mejor_orf[0]}",
        description=f"BRCA1 isoform 1 [Homo sapiens] | {mejor_orf[1]} aa | Marco {mejor_orf[0]}"
    )
    SeqIO.write([orf_record], out_orf, "fasta")
    print(f"[ARCHIVOS] Generados '{out_frames}' y '{out_orf}' exitosamente.\n")


if __name__ == "__main__":
    analizar_marcos_de_lectura()
