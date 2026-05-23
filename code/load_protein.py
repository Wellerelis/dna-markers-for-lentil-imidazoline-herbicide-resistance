#!/usr/bin/env python3
#load_protein.py — загрузка FASTA и перевод ДНК в белок

import sys

genetic_code = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G"
}


def load_fasta(filename):
    with open(filename, "r") as file:
        lines = file.readlines()

    header = ""
    sequence = ""

    for line in lines:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            header = line[1:]
        else:
            sequence = sequence + line

    return header, sequence


def check_sequence(dna_seq):
    errors = []

    for i, letter in enumerate(dna_seq):
        if letter not in "ATGC":
            errors.append(f"позиция {i}: неправильная буква '{letter}'")

    if len(dna_seq) % 3 != 0:
        errors.append(f"размер {len(dna_seq)} не кратен 3")

    if not dna_seq.startswith("ATG"):
        errors.append("не начинается с ATG")

    last_codon = dna_seq[-3:]
    if last_codon not in ["TAA", "TAG", "TGA"]:
        errors.append(f"не заканчивается стоп-кодоном, последний кодон: {last_codon}")

    return errors


def translate(dna_seq):
    protein = ""
    for i in range(0, len(dna_seq), 3):
        codon = dna_seq[i:i+3].upper()
        if len(codon) < 3:
            break
        protein = protein + genetic_code.get(codon, "X")
    return protein


def protein_stats(protein):
    stop_count = protein.count("*")
    at_end = stop_count == 1 and protein[-1] == "*"

    if protein.endswith("*"):
        length_no_stop = len(protein) - 1
    else:
        length_no_stop = len(protein)

    return {
        "length": len(protein),
        "first_aa": protein[0] if protein else "",
        "stop_codons": stop_count,
        "stop_at_end": at_end,
        "actual_length": length_no_stop,
    }


def main():
    input_file = "data/input/ahas_lentil.fasta"
    output_file = "data/output/protein_sequence.txt"

    print("анализ последовательности AHAS чечевицы")

    try:
        header, dna_sequence = load_fasta(input_file)
        print(f"последовательность: {header}")
        print(f"размер ДНК: {len(dna_sequence)} пар оснований")

        errors = check_sequence(dna_sequence)
        if errors:
            print("ошибки в последовательности:")
            for error in errors:
                print(f"  {error}")
            sys.exit(1)

        protein = translate(dna_sequence)
        stats = protein_stats(protein)

        print(f"размер белка: {stats['actual_length']} аминокислот")
        print(f"первая аминокислота: {stats['first_aa']}")
        print(f"стоп-кодонов: {stats['stop_codons']}, в конце: {stats['stop_at_end']}")

        print("последовательность белка:")
        for i in range(0, len(protein), 60):
            print(f"  {protein[i:i+60]}")

        with open(output_file, "w") as f:
            f.write(f">Lcu.2RBY.4g013440.1 (переведённый белок)\n")
            f.write(f"{protein}\n")

        print(f"сохранено в {output_file}")

    except FileNotFoundError:
        print(f"файл не найден: {input_file}")
        sys.exit(1)


if __name__ == "__main__":
    main()
