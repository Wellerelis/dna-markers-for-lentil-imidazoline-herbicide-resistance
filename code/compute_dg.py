#!/usr/bin/env python3
#compute_dg.py — расчёт ΔG связывания на основе геометрии и свойств мутаций
#этап 4: у нас нет AutoDock, поэтому используем физико-химическую модель
#на основе: (1) расстояния мутации от центра лиганда, (2) изменения свойств аминокислоты

import csv
import json
import os
import math
import numpy as np

R = 1.987e-3   #ккал/(моль*К)
T = 298.15


#гидрофобность по Kyte-Doolittle
hydro_kd = {
    "A": 1.8, "C": 2.5, "D": -3.5, "E": -3.5, "F": 2.8,
    "G": -0.4, "H": -3.2, "I": 4.5, "K": -3.9, "L": 3.8,
    "M": 1.9, "N": -3.5, "P": -1.6, "Q": -3.5, "R": -4.5,
    "S": -0.8, "T": -0.7, "V": 4.2, "W": -0.9, "Y": -1.3,
}

#объём боковой цепи (Å³, приблизительно)
volume = {
    "A": 67, "C": 86, "D": 91, "E": 109, "F": 135,
    "G": 48, "H": 118, "I": 124, "K": 135, "L": 124,
    "M": 124, "N": 96, "P": 90, "Q": 114, "R": 148,
    "S": 73, "T": 93, "V": 105, "W": 163, "Y": 141,
}

charge = {
    "D": -1, "E": -1, "K": 1, "R": 1, "H": 0.5,
}


def aa_charge(aa):
    return charge.get(aa, 0)


def parse_cif_ca(filepath):
    #координаты Cα каждого остатка
    atoms = {}
    fields = []
    in_block = False
    with open(filepath) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith("_atom_site."):
                fields.append(line.split(".")[1])
                in_block = True
                continue
            if in_block and line.startswith("ATOM"):
                parts = line.split()
                rec = dict(zip(fields, parts))
                if rec.get("label_atom_id") == "CA":
                    seq = int(rec["label_seq_id"])
                    atoms[seq] = [float(rec["Cartn_x"]), float(rec["Cartn_y"]), float(rec["Cartn_z"])]
            elif in_block and line.startswith("#"):
                in_block = False
    return atoms


def dist(a, b):
    return math.sqrt(sum((x-y)**2 for x, y in zip(a, b)))


def estimate_dg(wt_aa, mut_aa, ca_coord, ligand_center):
    #базовая ΔG дикого типа: примерно -6.5 ккал/моль (типичная для имидазолинов)
    dg_wt = -6.5

    d = dist(ca_coord, ligand_center)

    #влияние расстояния: чем дальше, тем меньше эффект
    distance_weight = math.exp(-d / 8.0)

    #изменение гидрофобности
    delta_hydro = hydro_kd[mut_aa] - hydro_kd[wt_aa]

    #изменение объёма
    delta_vol = (volume[mut_aa] - volume[wt_aa]) / 100.0

    #изменение заряда
    delta_charge = abs(aa_charge(mut_aa) - aa_charge(wt_aa))

    #вклад в ΔΔG:
    #гидрофобная замена вблизи лиганда снижает ΔG (улучшает связывание)
    #объёмная замена вблизи лиганда повышает ΔG (стерический конфликт)
    #зарядовая замена нестабилизирует
    ddg = distance_weight * (
        -0.4 * delta_hydro
        + 0.8 * abs(delta_vol)
        + 1.2 * delta_charge
    )

    dg = dg_wt + ddg
    return round(dg, 3), round(ddg, 3)


def main():
    with open("data/output/ligand_coords.json") as f:
        ligand_data = json.load(f)
    ligand_center = ligand_data["consensus_center"]


    # Используем все .cif файлы из input
    import glob
    cif_files = glob.glob("data/input/*.cif")
    ca_coords_all = {}
    for cif_file in cif_files:
        ca_coords = parse_cif_ca(cif_file)
        ca_coords_all[os.path.basename(cif_file)] = ca_coords

    # Для простоты: объединяем все координаты (если позиции совпадают, берём из первого файла)
    ca_coords = {}
    for coords in ca_coords_all.values():
        for k, v in coords.items():
            if k not in ca_coords:
                ca_coords[k] = v

    mutations = []
    with open("data/output/mutations_library.csv") as f:
        mutations = list(csv.DictReader(f))

    results = []
    for mut in mutations:
        pos = int(mut["position"])
        wt_aa = mut["wt_aa"]
        mut_aa = mut["mut_aa"]
        mut_type = mut["type"]

        if pos not in ca_coords:
            continue

        ca = ca_coords[pos]
        dg, ddg = estimate_dg(wt_aa, mut_aa, ca, ligand_center)

        pos2_val = ""
        wt2_val  = ""
        mut2_val = ""
        dist2_val = ""

        #для двойных — усредняем вклады и сохраняем данные второй замены
        if mut_type == "double" and mut.get("position_2"):
            pos2 = int(mut["position_2"])
            wt2 = mut["wt_aa_2"]
            mut2 = mut["mut_aa_2"]
            pos2_val = pos2
            wt2_val  = wt2
            mut2_val = mut2
            if pos2 in ca_coords:
                ca2 = ca_coords[pos2]
                dist2_val = round(dist(ca2, ligand_center), 2)
                dg2, ddg2 = estimate_dg(wt2, mut2, ca2, ligand_center)
                dg = round((dg + dg2) / 2 + ddg2 * 0.3, 3)
                ddg = round(dg - (-6.5), 3)

        row = {
            "id": mut["id"],
            "type": mut_type,
            "position": pos,
            "wt_aa": wt_aa,
            "mut_aa": mut_aa,
            "position_2": pos2_val,
            "wt_aa_2": wt2_val,
            "mut_aa_2": mut2_val,
            "description": mut["description"],
            "dg_kcal_mol": dg,
            "ddg_kcal_mol": ddg,
            "dg_wt": -6.5,
            "distance_to_ligand": round(dist(ca, ligand_center), 2),
            "distance_to_ligand_2": dist2_val,
        }

        results.append(row)

    out = "data/output/docking_results.csv"
    with open(out, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=results[0].keys())
        writer.writeheader()
        writer.writerows(results)

    print(f"рассчитано {len(results)} мутаций, сохранено в {out}")

    dgs = [r["dg_kcal_mol"] for r in results]
    print(f"ΔG: min={min(dgs):.2f}, max={max(dgs):.2f}, mean={sum(dgs)/len(dgs):.2f} ккал/моль")


if __name__ == "__main__":
    main()
