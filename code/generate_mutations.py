#!/usr/bin/env python3
#generate_mutations.py — генерация мутаций в участке связывания AHAS

import csv
import json
import sys
import os

aminoacids = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L",
              "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]

properties = {
    "A": {"hydrophobic": 1, "charge": 0, "size": 1, "polar": 0},
    "C": {"hydrophobic": 1, "charge": 0, "size": 1, "polar": 0},
    "D": {"hydrophobic": 0, "charge": -1, "size": 2, "polar": 1},
    "E": {"hydrophobic": 0, "charge": -1, "size": 2, "polar": 1},
    "F": {"hydrophobic": 1, "charge": 0, "size": 2, "polar": 0},
    "G": {"hydrophobic": 0, "charge": 0, "size": 0, "polar": 0},
    "H": {"hydrophobic": 0, "charge": 0.5, "size": 2, "polar": 1},
    "I": {"hydrophobic": 1, "charge": 0, "size": 2, "polar": 0},
    "K": {"hydrophobic": 0, "charge": 1, "size": 2, "polar": 1},
    "L": {"hydrophobic": 1, "charge": 0, "size": 2, "polar": 0},
    "M": {"hydrophobic": 1, "charge": 0, "size": 2, "polar": 0},
    "N": {"hydrophobic": 0, "charge": 0, "size": 1, "polar": 1},
    "P": {"hydrophobic": 0.5, "charge": 0, "size": 1, "polar": 0},
    "Q": {"hydrophobic": 0, "charge": 0, "size": 2, "polar": 1},
    "R": {"hydrophobic": 0, "charge": 1, "size": 2, "polar": 1},
    "S": {"hydrophobic": 0, "charge": 0, "size": 1, "polar": 1},
    "T": {"hydrophobic": 0, "charge": 0, "size": 1, "polar": 1},
    "V": {"hydrophobic": 1, "charge": 0, "size": 2, "polar": 0},
    "W": {"hydrophobic": 1, "charge": 0, "size": 3, "polar": 0},
    "Y": {"hydrophobic": 0.5, "charge": 0, "size": 3, "polar": 1},
}

BINDING_POSITIONS = [84, 85, 86, 88, 145, 146, 147, 151, 161, 163, 164, 165, 167, 168, 171]

FIELDS = ["id", "type", "position", "wt_aa", "mut_aa", "position_2", "wt_aa_2", "mut_aa_2", "description"]


def is_viable(wt_aa, mut_aa):
    if wt_aa == mut_aa:
        return False
    wt = properties[wt_aa]
    mt = properties[mut_aa]
    if wt_aa == "P" and mut_aa != "A":
        return False
    if mut_aa == "P" and wt_aa not in ["S", "T", "A", "G"]:
        return False
    if abs(wt["charge"] - mt["charge"]) > 1.5:
        return False
    if abs(wt["size"] - mt["size"]) > 2:
        return False
    return True


def generate_single(protein_seq, target_positions):
    mutations = []
    for pos in target_positions:
        idx = pos - 1
        if idx < 0 or idx >= len(protein_seq):
            continue
        wt_aa = protein_seq[idx]
        for mut_aa in aminoacids:
            if is_viable(wt_aa, mut_aa):
                mutations.append({
                    "id": f"SNP_{pos}_{wt_aa}{mut_aa}",
                    "type": "single",
                    "position": pos,
                    "wt_aa": wt_aa,
                    "mut_aa": mut_aa,
                    "position_2": "",
                    "wt_aa_2": "",
                    "mut_aa_2": "",
                    "description": f"{wt_aa}{pos}{mut_aa}",
                })
    return mutations


def generate_double(single_mutations, max_combos=60):
    multi = []
    count = 0
    for i in range(len(single_mutations)):
        for j in range(i + 1, len(single_mutations)):
            if count >= max_combos:
                break
            m1 = single_mutations[i]
            m2 = single_mutations[j]
            if m1["position"] == m2["position"]:
                continue
            multi.append({
                "id": f"DBL_{m1['position']}_{m1['mut_aa']}_{m2['position']}_{m2['mut_aa']}",
                "type": "double",
                "position": m1["position"],
                "wt_aa": m1["wt_aa"],
                "mut_aa": m1["mut_aa"],
                "position_2": m2["position"],
                "wt_aa_2": m2["wt_aa"],
                "mut_aa_2": m2["mut_aa"],
                "description": f"{m1['wt_aa']}{m1['position']}{m1['mut_aa']}+{m2['wt_aa']}{m2['position']}{m2['mut_aa']}",
            })
            count += 1
        if count >= max_combos:
            break
    return multi


def load_protein(filename):
    protein = ""
    with open(filename) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                continue
            protein += line
    return protein


def main():
    protein_file = "data/output/protein_sequence.txt"
    if not os.path.exists(protein_file):
        print(f"файл не найден: {protein_file}")
        sys.exit(1)

    protein = load_protein(protein_file)
    print(f"белок загружен: {len(protein)} аминокислот")

    single = generate_single(protein, BINDING_POSITIONS)
    double = generate_double(single, max_combos=60)
    all_mut = single + double

    print(f"одиночных: {len(single)}, двойных: {len(double)}, итого: {len(all_mut)}")

    os.makedirs("data/output", exist_ok=True)
    out = "data/output/mutations_library.csv"
    with open(out, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(all_mut)

    print(f"сохранено в {out}")

    with open("data/output/mutation_meta.json", "w") as f:
        json.dump({"binding_positions": BINDING_POSITIONS, "n_single": len(single),
                   "n_double": len(double), "protein_length": len(protein)}, f, indent=2)


if __name__ == "__main__":
    main()
