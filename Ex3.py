#!/usr/bin/env python3
"""
Trabajo Práctico de Bioinformática - Curso K5561 (UTN.BA)
Ejercicio 3: Alineamiento Múltiple de Secuencias (MSA) con Clustal Omega y Biopython
Gen de estudio: BRCA1 (Ortólogos en Mamíferos)
"""

import os
import time
import urllib.request
import urllib.parse
from Bio import AlignIO


ORTOGOS_BRCA1 = [
    ("Homo_sapiens", "P38398", "Humano (Consulta)"),
    ("Pan_troglodytes", "Q9GKK8", "Chimpancé"),
    ("Gorilla_gorilla", "Q6J6I8", "Gorila"),
    ("Pongo_pygmaeus", "Q6J6J0", "Orangután"),
    ("Macaca_mulatta", "Q6J6I9", "Macaco Rhesus"),
    ("Canis_lupus_familiaris", "Q95153", "Perro"),
    ("Bos_taurus", "Q864U1", "Bovino"),
    ("Mus_musculus", "P48754", "Ratón"),
]


def descargar_ortologos(fasta_in="msa_input.fasta"):
    bloques = []
    for especie, acc, _ in ORTOGOS_BRCA1:
        url = f"https://rest.uniprot.org/uniprotkb/{acc}.fasta"
        req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
        with urllib.request.urlopen(req) as resp:
            lineas = resp.read().decode("utf-8").strip().splitlines()
        secuencia = "".join(lineas[1:])
        bloques.append(f">{especie}_BRCA1_{acc}\n{secuencia}")
    with open(fasta_in, "w", encoding="utf-8") as f:
        f.write("\n".join(bloques) + "\n")


def analizar_msa(fasta_aln="msa_output.fasta", reporte_txt="msa_identity_report.txt"):
    aln = AlignIO.read(fasta_aln, "fasta")
    n = len(aln)
    ncols = aln.get_alignment_length()

    matriz = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i == j:
                matriz[i][j] = 100.0
            else:
                iguales = sum(1 for a, b in zip(aln[i].seq, aln[j].seq) if a == b and a != "-")
                validos = sum(1 for a, b in zip(aln[i].seq, aln[j].seq) if a != "-" or b != "-")
                matriz[i][j] = 100.0 * iguales / validos if validos else 0.0

    invariantes = sum(
        1 for c in range(ncols)
        if "-" not in [aln[i].seq[c] for i in range(n)] and len({aln[i].seq[c] for i in range(n)}) == 1
    )

    lineas = [
        "=" * 92,
        "EJERCICIO 3: ALINEAMIENTO MÚLTIPLE DE SECUENCIAS (MSA) - BRCA1 (CLUSTAL OMEGA)",
        "=" * 92,
        f"Secuencias analizadas          : {n} especies de mamíferos",
        f"Longitud total del alineamiento: {ncols} columnas",
        f"Posiciones 100% invariantes    : {invariantes} columnas ({100.0 * invariantes / ncols:.2f}%)",
        "-" * 92,
        "MATRIZ DE IDENTIDAD PORCENTUAL PAREADA (%):",
        "-" * 92,
    ]

    encabezado = f"{'Especie / Accesión UniProt':<34} | " + " | ".join(f"[{k+1:>1}]" for k in range(n))
    lineas.append(encabezado)
    lineas.append("-" * 92)

    for i, rec in enumerate(aln):
        nombre = f"[{i+1}] {rec.id[:28]}"
        valores = " | ".join(f"{matriz[i][j]:5.1f}" for j in range(n))
        lineas.append(f"{nombre:<34} | {valores}")

    lineas.append("-" * 92)
    texto_salida = "\n".join(lineas) + "\n"
    print(texto_salida)
    with open(reporte_txt, "w", encoding="utf-8") as f:
        f.write(texto_salida)


if __name__ == "__main__":
    if not os.path.exists("msa_input.fasta"):
        descargar_ortologos()
    analizar_msa()
