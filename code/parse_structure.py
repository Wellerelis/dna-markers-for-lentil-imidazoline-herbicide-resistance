#!/usr/bin/env python3
#parse_structure.py — разбор CIF файлов AlphaFold, извлечение pLDDT и геометрии участка

import os
import sys
import json
import numpy as np

def parse_cif_atoms(filepath):
    atoms = []
    in_atom_block = False
    fields = []

    with open(filepath) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith("_atom_site."):
                fields.append(line.split(".")[1])
                in_atom_block = True
                continue
            if in_atom_block and line.startswith("ATOM"):
                parts = line.split()
                rec = dict(zip(fields, parts))
                atoms.append(rec)
            elif in_atom_block and line.startswith("#"):
                in_atom_block = False

    return atoms


def residue_plddt(atoms):
    #средний pLDDT (B_iso) по каждому остатку
    sums = {}
    counts = {}
    names = {}
    for a in atoms:
        seq = int(a["label_seq_id"])
        b = float(a["B_iso_or_equiv"])
        sums[seq] = sums.get(seq, 0) + b
        counts[seq] = counts.get(seq, 0) + 1
        names[seq] = a["label_comp_id"]
    result = {}
    for seq in sums:
        result[seq] = {
            "plddt": sums[seq] / counts[seq],
            "resname": names[seq],
        }
    return result


def find_binding_region(res_plddt, threshold=78.0, min_run=5):
    #находит непрерывные участки с pLDDT выше порога
    high = sorted([r for r, v in res_plddt.items() if v["plddt"] >= threshold])
    if not high:
        return []
    regions = []
    run = [high[0]]
    for r in high[1:]:
        if r == run[-1] + 1:
            run.append(r)
        else:
            if len(run) >= min_run:
                regions.append(run[:])
            run = [r]
    if len(run) >= min_run:
        regions.append(run)
    return regions


def region_center(atoms, res_list):
    xs, ys, zs = [], [], []
    res_set = set(res_list)
    for a in atoms:
        if int(a["label_seq_id"]) in res_set and a["label_atom_id"] == "CA":
            xs.append(float(a["Cartn_x"]))
            ys.append(float(a["Cartn_y"]))
            zs.append(float(a["Cartn_z"]))
    if not xs:
        return None
    return [np.mean(xs), np.mean(ys), np.mean(zs)]


def main():
    cif_dir = "data/input"
    out_dir = "data/output"
    os.makedirs(out_dir, exist_ok=True)


    summary = {}
    cif_files = [f for f in os.listdir(cif_dir) if f.endswith('.cif')]
    for cif_file in cif_files:
        path = os.path.join(cif_dir, cif_file)
        atoms = parse_cif_atoms(path)
        rplddt = residue_plddt(atoms)

        mean_plddt = np.mean([v["plddt"] for v in rplddt.values()])
        regions = find_binding_region(rplddt, threshold=78.0, min_run=4)

        #выбираем регион с наибольшим средним pLDDT
        best_region = None
        best_score = 0
        for reg in regions:
            score = np.mean([rplddt[r]["plddt"] for r in reg])
            if score > best_score:
                best_score = score
                best_region = reg

        center = region_center(atoms, best_region) if best_region else None

        summary[cif_file] = {
            "mean_plddt": round(mean_plddt, 2),
            "n_atoms": len(atoms),
            "n_residues": len(rplddt),
            "high_plddt_regions": [{"residues": r, "mean_plddt": round(np.mean([rplddt[x]["plddt"] for x in r]), 2)} for r in regions],
            "best_region": best_region,
            "best_region_plddt": round(best_score, 2),
            "best_region_center": [round(c, 2) for c in center] if center else None,
        }

        print(f"{cif_file}: mean_pLDDT={mean_plddt:.2f}, лучший регион={best_region}, центр={center}")

    with open(os.path.join(out_dir, "structure_summary.json"), "w") as f:
        json.dump(summary, f, indent=2)

    #определяем лучшую модель и участок связывания для докинга
    best_model = max(summary, key=lambda k: summary[k]["best_region_plddt"])
    binding_site = {
        "model": best_model,
        "residues": summary[best_model]["best_region"],
        "center": summary[best_model]["best_region_center"],
        "plddt": summary[best_model]["best_region_plddt"],
    }
    with open(os.path.join(out_dir, "binding_site.json"), "w") as f:
        json.dump(binding_site, f, indent=2)

    print(f"\nлучшая модель: {best_model}, участок связывания: {binding_site['residues'][:5]}...{binding_site['residues'][-5:]}")
    print(f"центр участка: {binding_site['center']}")


if __name__ == "__main__":
    main()
