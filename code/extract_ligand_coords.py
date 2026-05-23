#!/usr/bin/env python3
#extract_ligand_coords.py — извлечение координат лигандов из всех CIF комплексов

import json
import os
import numpy as np


def parse_hetatm(filepath):
    atoms = []
    fields = []
    in_block = False
    with open(filepath) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith("_atom_site."):
                fields.append(line.split(".")[1])
                in_block = True
                continue
            if in_block and line.startswith("HETATM"):
                parts = line.split()
                rec = dict(zip(fields, parts))
                atoms.append(rec)
            elif in_block and line.startswith("#"):
                in_block = False
    return atoms


def ligand_center(atoms):
    xs = [float(a["Cartn_x"]) for a in atoms]
    ys = [float(a["Cartn_y"]) for a in atoms]
    zs = [float(a["Cartn_z"]) for a in atoms]
    return [round(np.mean(xs), 3), round(np.mean(ys), 3), round(np.mean(zs), 3)]


def main():
    complexes = {
        "imazaquin":     "data/input/imazaquin_complex.cif",
        "imazetapir":    "data/input/imazetapir_complex.cif",
        "imazametabenz": "data/input/imazametabenz_complex.cif",
        "imazamox":      "data/input/imazamox_complex.cif",
        "imazapyc":      "data/input/imazapyc_complex.cif",
        "imazapyr":      "data/input/imazapyr_complex.cif",
    }

    results = {}
    for name, path in complexes.items():
        if not os.path.exists(path):
            print(f"пропуск: {path} не найден")
            continue
        atoms = parse_hetatm(path)
        center = ligand_center(atoms)
        results[name] = {"n_atoms": len(atoms), "center": center}
        print(f"{name}: {len(atoms)} атомов, центр {center}")

    centers = [v["center"] for v in results.values()]
    consensus = [round(np.mean([c[i] for c in centers]), 3) for i in range(3)]
    results["consensus_center"] = consensus
    print(f"\nконсенсусный центр ({len(centers)} лигандов): {consensus}")

    with open("data/output/ligand_coords.json", "w") as f:
        json.dump(results, f, indent=2)
    print("сохранено в data/output/ligand_coords.json")


if __name__ == "__main__":
    main()
