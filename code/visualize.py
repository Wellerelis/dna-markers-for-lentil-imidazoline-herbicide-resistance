#!/usr/bin/env python3
#visualize.py — все публикационные графики (VAK: чистый стиль, 300dpi)

import csv
import json
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 10,
    "axes.labelsize": 11,
    "axes.titlesize": 12,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "figure.dpi": 300,
    "axes.linewidth": 0.8,
    "lines.linewidth": 1.2,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "legend.framealpha": 0.7,
    "legend.fontsize": 9,
})

os.makedirs("results/figures", exist_ok=True)


def load_results():
    rows = []
    with open("data/output/regression_results.csv") as f:
        for r in csv.DictReader(f):
            r["dg_predicted"] = float(r["dg_predicted"])
            r["dg_std"] = float(r["dg_std"])
            r["ddg"] = float(r["ddg"])
            r["distance_to_ligand"] = float(r["distance_to_ligand"])
            r["dg_measured"] = float(r["dg_measured"])
            rows.append(r)
    return rows


def load_top():
    top = []
    with open("data/output/top_candidates.csv") as f:
        for r in csv.DictReader(f):
            r["dg_predicted"] = float(r["dg_predicted"])
            r["ddg"] = float(r["ddg"])
            top.append(r)
    return top


def load_minima():
    with open("data/output/energy_minima.json") as f:
        d = json.load(f)
    return d["minima"], d["landscape"]


def plot_landscape(rows, top, minima, landscape):
    single = [r for r in rows if r["type"] == "single"]
    positions = [int(r["position"]) for r in single]
    dg = [r["dg_predicted"] for r in single]
    ddg = [r["ddg"] for r in single]

    fig, ax = plt.subplots(figsize=(12, 5))

    norm = Normalize(vmin=min(dg), vmax=max(dg))
    cmap = plt.cm.RdYlGn_r

    sc = ax.scatter(positions, dg, c=dg, cmap=cmap, norm=norm,
                    s=40, alpha=0.65, linewidths=0.3, edgecolors="black", zorder=3)

    #минимумы ΔG
    for m in minima:
        ax.scatter(m["position"], m["mean_dg"], s=160, marker="v",
                   color="#e63946", edgecolors="black", linewidths=0.6, zorder=5)

    #топ-10 по устойчивости
    for t in top[:5]:
        ax.scatter(int(t["position"]), t["dg_predicted"], s=180, marker="*",
                   color="#2d6a4f", edgecolors="black", linewidths=0.5, zorder=6)
        ax.annotate(t["description"], (int(t["position"]), t["dg_predicted"]),
                    textcoords="offset points", xytext=(4, 4), fontsize=7)

    ax.axhline(-6.5, color="#6c757d", linestyle="--", linewidth=0.9, alpha=0.6,
               label="ΔG дикого типа (−6.5 ккал/моль)")

    cbar = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), ax=ax, pad=0.01, shrink=0.85)
    cbar.set_label("ΔG (ккал/моль)", fontsize=9)

    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], linestyle="--", color="#6c757d", label="ΔG дикого типа"),
        plt.scatter([], [], s=160, marker="v", color="#e63946", label="локальный минимум ΔG"),
        plt.scatter([], [], s=180, marker="*", color="#2d6a4f", label="топ-5 кандидатов"),
    ]
    ax.legend(handles=legend_elements, loc="upper right")
    ax.set_xlabel("позиция остатка")
    ax.set_ylabel("ΔG (ккал/моль)")
    ax.set_title("энергетический ландшафт мутаций AHAS")
    ax.grid(True, alpha=0.18, linewidth=0.5)

    plt.tight_layout()
    plt.savefig("results/figures/energy_landscape.png", dpi=300, bbox_inches="tight")
    plt.close()
    print("energy_landscape.png")


def plot_pred_vs_actual(rows):
    train  = [r for r in rows if r.get("split") == "train"]
    test   = [r for r in rows if r.get("split") == "test"]
    double = [r for r in rows if r.get("split") == "double"]

    fig, ax = plt.subplots(figsize=(6, 6))

    all_vals = [float(r["dg_measured"]) for r in rows] + [float(r["dg_predicted"]) for r in rows]
    mn = min(all_vals) - 0.15
    mx = max(all_vals) + 0.15
    ax.plot([mn, mx], [mn, mx], "k--", linewidth=0.9, alpha=0.5, label="идеальное предсказание")

    if train:
        xa = np.array([float(r["dg_measured"]) for r in train])
        xp = np.array([float(r["dg_predicted"]) for r in train])
        ax.scatter(xa, xp, s=18, alpha=0.35, color="#adb5bd",
                   label=f"обучение (n={len(train)})", zorder=2)

    if test:
        xa = np.array([float(r["dg_measured"]) for r in test])
        xp = np.array([float(r["dg_predicted"]) for r in test])
        xs = np.array([float(r["dg_std"]) for r in test])
        ax.errorbar(xa, xp, yerr=xs, fmt="o", color="#3a86ff",
                    markersize=5, alpha=0.85, linewidth=0.5, capsize=2, elinewidth=0.5,
                    label=f"тест — новые позиции (n={len(test)})", zorder=4)
        r2_test = float(np.corrcoef(xa, xp)[0, 1] ** 2)
        rmse_test = float(np.sqrt(np.mean((xa - xp)**2)))
        ax.text(0.05, 0.88,
                f"R\u00b2 тест = {r2_test:.3f}\nRMSE = {rmse_test:.3f} ккал/моль",
                transform=ax.transAxes, fontsize=9,
                bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))

    if double:
        xa = np.array([float(r["dg_measured"]) for r in double])
        xp = np.array([float(r["dg_predicted"]) for r in double])
        ax.scatter(xa, xp, s=22, alpha=0.55, color="#f4a261", marker="^",
                   label=f"двойные мутации (n={len(double)})", zorder=3)

    ax.set_xlabel("ΔG расчётный (ккал/моль)")
    ax.set_ylabel("ΔG предсказанный GP (ккал/моль)")
    ax.set_title("качество модели Gaussian Process\n(тест — позиции не входившие в обучение)")
    ax.legend(fontsize=8, loc="upper left")
    ax.grid(True, alpha=0.18, linewidth=0.5)
    ax.set_aspect("equal", adjustable="box")

    plt.tight_layout()
    plt.savefig("results/figures/predicted_vs_actual.png", dpi=300, bbox_inches="tight")
    plt.close()
    print("predicted_vs_actual.png")

def plot_ddg_distribution(rows):
    single = [r for r in rows if r["type"] == "single"]
    ddg = [r["ddg"] for r in single]

    fig, ax = plt.subplots(figsize=(8, 4))

    ax.hist(ddg, bins=30, color="#457b9d", edgecolor="black", linewidth=0.3, alpha=0.85)
    ax.axvline(0, color="#e63946", linestyle="--", linewidth=1.1, alpha=0.8,
               label="ΔΔG = 0 (дикий тип)")
    ax.axvline(0.3, color="#2d6a4f", linestyle="--", linewidth=1.1, alpha=0.8,
               label="порог устойчивости (ΔΔG > 0.3)")

    n_resistant = sum(1 for x in ddg if x > 0.3)
    ax.text(0.97, 0.93, f"потенциально устойчивых: {n_resistant}",
            transform=ax.transAxes, ha="right", fontsize=9,
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.7))

    ax.set_xlabel("ΔΔG (ккал/моль)")
    ax.set_ylabel("число мутаций")
    ax.set_title("распределение ΔΔG одиночных мутаций относительно дикого типа")
    ax.legend()
    ax.grid(True, alpha=0.18, linewidth=0.5, axis="y")

    plt.tight_layout()
    plt.savefig("results/figures/ddg_distribution.png", dpi=300, bbox_inches="tight")
    plt.close()
    print("ddg_distribution.png")


def plot_heatmap(rows):
    binding_pos = [84, 85, 86, 88, 145, 146, 147, 151, 161, 163, 164, 165, 167, 168, 171]
    single = [r for r in rows if r["type"] == "single" and int(r["position"]) in binding_pos]

    aas = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L",
           "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]
    pos_list = sorted(set(int(r["position"]) for r in single))

    wt_by_pos = {}
    for r in single:
        wt_by_pos[int(r["position"])] = r["wt_aa"]

    grid = np.full((len(pos_list), len(aas)), np.nan)
    pos_idx = {p: i for i, p in enumerate(pos_list)}
    aa_idx = {a: i for i, a in enumerate(aas)}

    for r in single:
        pi = pos_idx[int(r["position"])]
        ai = aa_idx.get(r["mut_aa"])
        if ai is not None:
            grid[pi, ai] = r["ddg"]

    fig, ax = plt.subplots(figsize=(15, 6))
    im = ax.imshow(grid, cmap="RdYlGn", aspect="auto", vmin=-0.5, vmax=1.0)

    ax.set_xticks(range(len(aas)))
    ax.set_xticklabels(aas, fontsize=8)
    ax.set_yticks(range(len(pos_list)))
    ylabels = [f"{p} ({wt_by_pos[p]})" for p in pos_list]
    ax.set_yticklabels(ylabels, fontsize=8)
    ax.set_xlabel("мутантная аминокислота")
    ax.set_ylabel("позиция (дикий тип)")
    ax.set_title("тепловая карта ΔΔG в участке связывания")

    cbar = fig.colorbar(im, ax=ax, pad=0.01, shrink=0.8)
    cbar.set_label("ΔΔG (ккал/моль)", fontsize=9)

    #отмечаем NaN серым
    for i in range(len(pos_list)):
        for j in range(len(aas)):
            if np.isnan(grid[i, j]):
                ax.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1,
                                           fill=True, color="#e9ecef", zorder=0))

    plt.tight_layout()
    plt.savefig("results/figures/heatmap_ddg.png", dpi=300, bbox_inches="tight")
    plt.close()
    print("heatmap_ddg.png")


def plot_top10_bar(top):
    fig, ax = plt.subplots(figsize=(9, 5))

    labels = [r["description"] for r in top]
    ddgs = [r["ddg"] for r in top]
    colors = ["#2d6a4f" if d > 0.5 else "#52b788" if d > 0.3 else "#95d5b2" for d in ddgs]

    bars = ax.barh(labels[::-1], ddgs[::-1], color=colors[::-1],
                   edgecolor="black", linewidth=0.4, height=0.6)

    ax.axvline(0.3, color="#e63946", linestyle="--", linewidth=1.0, alpha=0.7,
               label="порог устойчивости")
    ax.axvline(0, color="black", linewidth=0.6, alpha=0.4)

    for bar, d in zip(bars, ddgs[::-1]):
        ax.text(d + 0.01, bar.get_y() + bar.get_height()/2,
                f"{d:.3f}", va="center", fontsize=8)

    ax.set_xlabel("ΔΔG (ккал/моль)")
    ax.set_title("топ-10 кандидатов на устойчивость к имидазолиновым гербицидам")
    ax.legend()
    ax.grid(True, alpha=0.18, linewidth=0.5, axis="x")

    plt.tight_layout()
    plt.savefig("results/figures/top10_candidates.png", dpi=300, bbox_inches="tight")
    plt.close()
    print("top10_candidates.png")


def plot_p2rank_consensus():
    binding_data = {
        "imazamox":      [86, 145, 146, 147, 163, 164, 165, 168],
        "imazapyc":      [147, 164, 165, 168],
        "imazametabenz": [88, 151, 161],
        "imazetapir":    [84, 85, 86, 147, 165, 167, 168, 171],
    }
    all_pos = set()
    for v in binding_data.values():
        all_pos.update(v)
    all_pos = sorted(all_pos)

    matrix = np.zeros((len(binding_data), len(all_pos)))
    herbicides = list(binding_data.keys())
    for i, h in enumerate(herbicides):
        for j, p in enumerate(all_pos):
            if p in binding_data[h]:
                matrix[i, j] = 1

    fig, ax = plt.subplots(figsize=(12, 3.5))
    im = ax.imshow(matrix, cmap="Blues", aspect="auto", vmin=0, vmax=1)

    ax.set_xticks(range(len(all_pos)))
    ax.set_xticklabels(all_pos, fontsize=8)
    ax.set_yticks(range(len(herbicides)))
    ax.set_yticklabels(herbicides, fontsize=9)
    ax.set_xlabel("позиция остатка")
    ax.set_title("участок связывания: консенсус по p2rank предсказаниям для 4 гербицидов")

    #частота появления
    freq = matrix.sum(axis=0)
    for j, (p, f) in enumerate(zip(all_pos, freq)):
        if f >= 2:
            ax.text(j, len(herbicides) - 0.1, f"{int(f)}×",
                    ha="center", va="bottom", fontsize=7, color="#e63946")

    plt.tight_layout()
    plt.savefig("results/figures/p2rank_consensus.png", dpi=300, bbox_inches="tight")
    plt.close()
    print("p2rank_consensus.png")


def main():
    rows = load_results()
    top = load_top()
    minima, landscape = load_minima()

    plot_landscape(rows, top, minima, landscape)
    plot_pred_vs_actual(rows)
    plot_ddg_distribution(rows)
    plot_heatmap(rows)
    plot_top10_bar(top)
    plot_p2rank_consensus()

    print("все 6 графиков сохранены в results/figures/")


if __name__ == "__main__":
    main()
