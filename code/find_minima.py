#!/usr/bin/env python3
#find_minima.py — поиск минимумов на поверхности GP модели

import csv
import json
import pickle
import numpy as np
from scipy.signal import find_peaks

hydro_kd = {
    "A": 1.8, "C": 2.5, "D": -3.5, "E": -3.5, "F": 2.8,
    "G": -0.4, "H": -3.2, "I": 4.5, "K": -3.9, "L": 3.8,
    "M": 1.9, "N": -3.5, "P": -1.6, "Q": -3.5, "R": -4.5,
    "S": -0.8, "T": -0.7, "V": 4.2, "W": -0.9, "Y": -1.3,
}
volume = {
    "A": 67, "C": 86, "D": 91, "E": 109, "F": 135,
    "G": 48, "H": 118, "I": 124, "K": 135, "L": 124,
    "M": 124, "N": 96, "P": 90, "Q": 114, "R": 148,
    "S": 73, "T": 93, "V": 105, "W": 163, "Y": 141,
}
charge_val = {"D": -1, "E": -1, "K": 1, "R": 1, "H": 0.5}
polarity = {
    "A": 0, "C": 0, "D": 1, "E": 1, "F": 0, "G": 0, "H": 1,
    "I": 0, "K": 1, "L": 0, "M": 0, "N": 1, "P": 0, "Q": 1,
    "R": 1, "S": 1, "T": 1, "V": 0, "W": 0, "Y": 1,
}


def featurize(row):
    pos = float(row["position"]) / 188.0
    wt = row["wt_aa"]
    mt = row["mut_aa"]
    dist = float(row["distance_to_ligand"])
    return [
        pos, dist,
        hydro_kd[mt] - hydro_kd[wt],
        (volume[mt] - volume[wt]) / 100.0,
        charge_val.get(mt, 0) - charge_val.get(wt, 0),
        polarity[mt] - polarity[wt],
        hydro_kd[wt],
        hydro_kd[mt],
    ]


def main():
    with open("data/output/models/gp_model.pkl", "rb") as f:
        bundle = pickle.load(f)
    gp = bundle["model"]
    scaler = bundle["scaler"]

    rows = list(csv.DictReader(open("data/output/regression_results.csv")))
    single = [r for r in rows if r["type"] == "single"]

    positions = sorted(set(int(r["position"]) for r in single))
    pos_to_rows = {}
    for r in single:
        p = int(r["position"])
        pos_to_rows.setdefault(p, []).append(r)

    #для каждой позиции — минимальный dg_predicted (наибольшее связывание)
    #и максимальный dg_predicted (наибольшая устойчивость)
    landscape = []
    for pos in positions:
        pos_rows = pos_to_rows[pos]
        dgs = [float(r["dg_predicted"]) for r in pos_rows]
        ddgs = [float(r["ddg"]) for r in pos_rows]
        best_resistance = max(zip(ddgs, pos_rows), key=lambda x: x[0])
        landscape.append({
            "position": pos,
            "mean_dg": round(sum(dgs)/len(dgs), 3),
            "min_dg": round(min(dgs), 3),
            "max_dg": round(max(dgs), 3),
            "best_mut": best_resistance[1]["description"],
            "best_ddg": round(best_resistance[0], 3),
        })

    #поиск минимумов ΔG (позиции, где среднее ΔG ниже соседей)
    mean_dgs = np.array([l["mean_dg"] for l in landscape])
    #инвертируем для find_peaks (ищем минимумы)
    peaks_idx, props = find_peaks(-mean_dgs, prominence=0.05)

    minima = []
    for idx in peaks_idx:
        l = landscape[idx]
        minima.append({
            "position": l["position"],
            "mean_dg": l["mean_dg"],
            "best_mut": l["best_mut"],
            "best_ddg": l["best_ddg"],
            "type": "local_minimum",
        })

    print(f"найдено {len(minima)} локальных минимумов ΔG:")
    for m in minima:
        print(f"  позиция {m['position']}: mean ΔG={m['mean_dg']}, лучшая мутация={m['best_mut']} (ΔΔG={m['best_ddg']})")

    with open("data/output/energy_minima.json", "w") as f:
        json.dump({"minima": minima, "landscape": landscape}, f, indent=2)

    print(f"сохранено в data/output/energy_minima.json")


if __name__ == "__main__":
    main()
