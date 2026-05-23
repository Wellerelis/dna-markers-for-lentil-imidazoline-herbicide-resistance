#!/usr/bin/env python3
#design_primers.py — дизайн qPCR праймеров для топ-10 кандидатов

import csv
import json
import os

#полная ДНК последовательность (564 bp)
DNA = "ATGATTAACCGAATCGCTGGCGTTTTTGCGCGTCGCGGCTTTAATATCGAGAGCCTGGTCGTCGGCCTGAATCAGGATCGCGCACTGTTTAGCATCGTTGTTTATGGCGCTGATAACATCATCCGCCAGCTGCAGAAGCTCGTCAATGTGCTCAAAGTATCGCTGAGCCCATTCTGCGATATTTTGCTAATTTGCGTTTGCCGTCAGCGTATGAACAGTAACCTGCGTGTGAAGGATCTGAGTAACGAGCCGCAGGTGGAGCGTGAGCTTATGCTGCTGAAGGTTCATGCTGACCCGCAGAATCAGTCAGAAGTCCCGGCGCTCCTGTTCCTCGTCTGCTGTATCCAGACGTGCGTTCTGTTCCACTGCCTGCTGTTGTTTATCTATCTGTTCGATTTCTTTGTCGTCGGCAGTTTGCAAATCAAATGGCTCGTTGATACCTTTCGCGCAAAAATCGTCGATATATCGGAGGATTATGTGACTACGGAAGTCATTGGTGATCCGGGCAAAATGTTTGCTGTCCAGCGTATATTTAAGTTCGGCATCAAAGAGATTGCACGTACGGGC"

codon_table = {
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
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

#предпочтительные кодоны для каждой аминокислоты (E. coli / plant)
preferred_codon = {
    "A": "GCT", "C": "TGC", "D": "GAT", "E": "GAA", "F": "TTT",
    "G": "GGT", "H": "CAT", "I": "ATT", "K": "AAA", "L": "CTG",
    "M": "ATG", "N": "AAT", "P": "CCT", "Q": "CAA", "R": "CGT",
    "S": "TCT", "T": "ACT", "V": "GTT", "W": "TGG", "Y": "TAT",
}

complement = str.maketrans("ATGCatgc", "TACGtacg")


def rev_comp(seq):
    return seq.translate(complement)[::-1]


def gc_content(seq):
    gc = sum(1 for b in seq.upper() if b in "GC")
    return round(100 * gc / len(seq), 1)


def tm_basic(seq):
    #формула Уоллеса для коротких праймеров
    at = seq.upper().count("A") + seq.upper().count("T")
    gc = seq.upper().count("G") + seq.upper().count("C")
    return 2 * at + 4 * gc


def design_primer_pair(pos_1based, wt_aa, mut_aa, flank=10):
    #позиция в ДНК (0-based начало кодона)
    codon_start = (pos_1based - 1) * 3
    codon_end = codon_start + 3

    new_codon = preferred_codon[mut_aa]

    #мутантная ДНК
    mut_dna = DNA[:codon_start] + new_codon + DNA[codon_end:]

    #форвард праймер: flank нуклеотидов до мутации + мутантный кодон + flank после
    fwd_start = max(0, codon_start - flank)
    fwd_end = min(len(mut_dna), codon_end + flank)
    fwd = mut_dna[fwd_start:fwd_end]

    #ревёрс праймер: комплементарный к участку за мутацией
    rev_start = max(0, codon_start - flank)
    rev_end = min(len(mut_dna), codon_end + flank)
    rev = rev_comp(mut_dna[rev_start:rev_end])

    return {
        "forward": fwd,
        "reverse": rev,
        "new_codon": new_codon,
        "fwd_tm": tm_basic(fwd),
        "rev_tm": tm_basic(rev),
        "fwd_gc": gc_content(fwd),
        "rev_gc": gc_content(rev),
        "fwd_len": len(fwd),
        "rev_len": len(rev),
    }


def main():
    top = list(csv.DictReader(open("data/output/top_candidates.csv")))

    results = []
    for r in top:
        pos = int(r["position"])
        wt = r["wt_aa"]
        mut = r["mut_aa"]
        primers = design_primer_pair(pos, wt, mut)

        results.append({
            "id": r["id"],
            "description": r["description"],
            "position": pos,
            "wt_aa": wt,
            "mut_aa": mut,
            "dg_predicted": r["dg_predicted"],
            "ddg": r["ddg"],
            "new_codon": primers["new_codon"],
            "forward_primer": primers["forward"],
            "reverse_primer": primers["reverse"],
            "fwd_len": primers["fwd_len"],
            "rev_len": primers["rev_len"],
            "fwd_tm": primers["fwd_tm"],
            "rev_tm": primers["rev_tm"],
            "fwd_gc_pct": primers["fwd_gc"],
            "rev_gc_pct": primers["rev_gc"],
        })

    out = "data/output/primers_designed.csv"
    with open(out, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=results[0].keys())
        writer.writeheader()
        writer.writerows(results)

    print(f"праймеры для {len(results)} кандидатов сохранены в {out}")
    print()
    for r in results:
        print(f"{r['description']}: FWD {r['fwd_tm']}°C/{r['fwd_gc_pct']}%GC  REV {r['rev_tm']}°C/{r['rev_gc_pct']}%GC")


if __name__ == "__main__":
    main()
