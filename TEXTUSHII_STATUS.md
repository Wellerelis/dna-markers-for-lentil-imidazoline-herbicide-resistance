# Текущий статус проекта AHAS

Дата обновления: 23 май 2026

## Что готово

Этапы 1–5 завершены.

Этап 1: последовательность белка (188 аминокислот), скрипты load_protein.py и generate_mutations.py.

Этап 2: участок связывания определён из p2rank предсказаний для 6 имидазолиновых гербицидов. Консенсусные позиции: 84, 85, 86, 88, 145, 146, 147, 151, 161, 163, 164, 165, 167, 168, 171.

Этап 3: сгенерировано 307 мутаций (247 одиночных + 60 двойных). ΔG рассчитан по физико-химической модели на основе координат из AlphaFold (rank_4) и консенсусного центра трёх лигандов (имазаквин, имазетапир, имазаметабенз).

Этап 4: Gaussian Process регрессия обучена на одиночных мутациях. R² = 0.9985, RMSE = 0.0082 ккал/моль.

Этап 5: выбраны топ-10 кандидатов по максимальному ΔΔG.

## Топ-10 кандидатов

1. I151R — ΔG = -5.778, ΔΔG = 0.723 ккал/моль
2. V163R — ΔG = -5.804, ΔΔG = 0.696 ккал/моль
3. I151D — ΔG = -5.825, ΔΔG = 0.676 ккал/моль
4. I151K — ΔG = -5.826, ΔΔG = 0.673 ккал/моль
5. I151E — ΔG = -5.846, ΔΔG = 0.655 ккал/моль
6. V163K — ΔG = -5.852, ΔΔG = 0.648 ккал/моль
7. T161R — ΔG = -5.862, ΔΔG = 0.640 ккал/моль
8. V163D — ΔG = -5.892, ΔΔG = 0.608 ккал/моль
9. V163E — ΔG = -5.903, ΔΔG = 0.597 ккал/моль
10. I164R — ΔG = -5.920, ΔΔG = 0.579 ккал/моль

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
