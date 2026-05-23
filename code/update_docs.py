#!/usr/bin/env python3
#update_docs.py — обновляем статусные документы по итогам этапов 2-5

import json
import csv
import os


def main():
    with open("data/output/regression_metrics.json") as f:
        metrics = json.load(f)

    top = []
    with open("data/output/top_candidates.csv") as f:
        top = list(csv.DictReader(f))

    status = f"""# Текущий статус проекта AHAS

Дата обновления: 23 май 2026

## Что готово

Этапы 1–5 завершены.

Этап 1: последовательность белка (188 аминокислот), скрипты load_protein.py и generate_mutations.py.

Этап 2: участок связывания определён из p2rank предсказаний для 6 имидазолиновых гербицидов. Консенсусные позиции: 84, 85, 86, 88, 145, 146, 147, 151, 161, 163, 164, 165, 167, 168, 171.

Этап 3: сгенерировано 307 мутаций (247 одиночных + 60 двойных). ΔG рассчитан по физико-химической модели на основе координат из AlphaFold (rank_4) и консенсусного центра трёх лигандов (имазаквин, имазетапир, имазаметабенз).

Этап 4: Gaussian Process регрессия обучена на одиночных мутациях. R² = {metrics['r2_test']}, RMSE = {metrics['rmse']} ккал/моль.

Этап 5: выбраны топ-10 кандидатов по максимальному ΔΔG.

## Топ-10 кандидатов

"""
    for i, r in enumerate(top, 1):
        status += f"{i}. {r['description']} — ΔG = {float(r['dg_predicted']):.3f}, ΔΔG = {float(r['ddg']):.3f} ккал/моль\n"

    status += """
## Файлы

```
data/input/rank_1..5.cif                   AlphaFold модели белка
data/input/imazaquin_complex.cif           комплекс с имазаквином
data/input/imazetapir_complex.cif          комплекс с имазетапиром
data/input/imazametabenz_complex.cif       комплекс с имазаметабензом
data/input/prediction_*.json               p2rank предсказания (6 гербицидов)
data/output/protein_sequence.txt           белок 188 аминокислот
data/output/binding_site.json             участок связывания
data/output/ligand_coords.json            координаты центров лигандов
data/output/mutations_library.csv         307 мутаций
data/output/docking_results.csv           ΔG для всех мутаций
data/output/regression_results.csv        предсказания GP модели
data/output/regression_metrics.json       метрики модели
data/output/top_candidates.csv            топ-10
data/output/all_single_ranked.csv         все одиночные по рангу
data/output/models/gp_model.pkl           обученная GP модель
results/figures/energy_landscape.png      ландшафт ΔG
results/figures/predicted_vs_actual.png   предсказано vs расчёт
results/figures/ddg_distribution.png      распределение ΔΔG
results/figures/heatmap_ddg.png           тепловая карта
```

## Следующие шаги

Этап 6: поиск минимумов на поверхности GP модели.
Этап 7: дизайн праймеров qPCR для топ-10 кандидатов.
Этап 8: финальный отчёт и презентация.
"""

    with open("TEXTUSHII_STATUS.md", "w") as f:
        f.write(status)
    print("TEXTUSHII_STATUS.md обновлён")


if __name__ == "__main__":
    main()
