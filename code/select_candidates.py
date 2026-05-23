#!/usr/bin/env python3
#select_candidates.py — выбор топ кандидатов и сохранение итоговой таблицы

import csv
import json
import os


def main():
    rows = []
    with open("data/output/regression_results.csv") as f:
        rows = list(csv.DictReader(f))

    for r in rows:
        r["dg_predicted"] = float(r["dg_predicted"])
        r["dg_std"] = float(r["dg_std"])
        r["ddg"] = float(r["ddg"])
        r["distance_to_ligand"] = float(r["distance_to_ligand"])

    #сортируем одиночные по dg_predicted (чем ниже, тем сильнее связывание)
    #для устойчивости нужно слабое связывание (высокий ΔG), то есть нас интересуют
    #мутации с максимальным ΔΔG > 0 (ослабление связывания с гербицидом)
    single = [r for r in rows if r["type"] == "single"]

    #ранжируем по ddg убывающий (наибольшее ослабление связывания = наибольшая устойчивость)
    single_sorted = sorted(single, key=lambda x: -x["ddg"])

    top10 = single_sorted[:10]

    print("топ-10 кандидатов по устойчивости (наибольший ΔΔG):")
    print(f"{'#':<3} {'позиция':<8} {'мутация':<12} {'ΔG':<8} {'ΔΔG':<8} {'расст.':<8}")
    for i, r in enumerate(top10, 1):
        print(f"{i:<3} {r['position']:<8} {r['description']:<12} {r['dg_predicted']:<8.3f} {r['ddg']:<8.3f} {r['distance_to_ligand']:<8.2f}")

    out = "data/output/top_candidates.csv"
    with open(out, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=top10[0].keys())
        writer.writeheader()
        writer.writerows(top10)

    print(f"\nсохранено в {out}")

    #также сохраняем полную таблицу для отчёта
    all_sorted = sorted(single, key=lambda x: -x["ddg"])
    with open("data/output/all_single_ranked.csv", "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=all_sorted[0].keys())
        writer.writeheader()
        writer.writerows(all_sorted)


if __name__ == "__main__":
    main()
