# Извлечение данных и расчеты ЯМР анилита

Редакция 2, 26 ноября 2026 г. В этой папке сохраняются исходные выгрузки и их аудит. Порошковая подгонка и пересчет затухания эха выполняются отдельными скриптами в соседних папках `powder_model/` и `relaxation_analysis/`.

## Исходные файлы

Для численного извлечения нужны:

- `upload/Anilite.OPJ`
- `upload/Anilite14072016.zi`
- `upload/CuS.zi`
- `upload/Attachments_irek.mukhamedshin@gmail.com_2016-07-04_09-39-32.zip`

Расширения `.zi` менять не требуется: это ZIP-архивы. Все десять пользовательских приложений входят в общий пакет; остальные документы задают экспериментальный и литературный контекст.

## Зависимости и начальное извлечение

Нужны Python 3.10+, NumPy, SciPy и Matplotlib. Для выполненного анализа использованы NumPy 2.3.5, SciPy 1.17.0, Matplotlib 3.10.8. Из корневой папки распакованного пакета:

```bash
python data_analysis/analyze_measurements.py
python data_analysis/summarize_manuscript.py
```

`analyze_measurements.py` безопасно извлекает три архива, читает 30 колонок Origin, сохраняет таблицы, метаданные, флаги и диагностические рисунки. Макросы/формулы Origin не исполняются. `summarize_manuscript.py` воспроизводит описательные суммы, частотные максимумы и аналитические интервалы νQ. Эти два скрипта сами не подгоняют порошковый спектр и огибающую T2.

## Пересчет затухания эха

После извлечения `opj_T2a.csv`:

```bash
python relaxation_analysis/analyze_t2.py
```

Принято `TwoTau` в микросекундах и `t=2τ`. Основание — экспериментальный контекст автора и дипломные масштабы времени. Амплитуда равна `sqrt(IntRe**2+IntIm**2)`; фаза показана отдельно. Основной результат обобщенной экспоненты — эффективное время около 0.955 мс и β≈1.255. Подмножества и другие определения амплитуды служат проверкой чувствительности, не доверительными интервалами.

Итоги находятся в `relaxation_analysis/export/`: `t2_summary.json`, `t2_fit_sensitivity.csv`, `t2_envelope_points.csv`, `figure_t2_reanalysis.png/.pdf`.

## Порошковая модель

Автор подтвердил одинаковый неориентированный порошок для ЯКР/ЯМР. Подгонка использует две серии интегралов 20/500 мс. Окна±0.12 МГц около 75.941 и 81.346 МГц исключают узкие линии меди датчика из целевой функции. После исключения конфликтующих строк 81.40 МГц и аппаратных окон остается по 245 точек.

Главные результаты:

- `powder_model/MODEL_ANALYSIS_RU.txt` — определения и физические ограничения.
- `powder_model/FINAL_MODEL_COMPARISON.csv` — итоговая таблица семейств и ветвей решений.
- `powder_model/central_exact_validation.json` — прямой пересчет на плотной ориентационной сетке.
- `powder_model/powder_model_comparison.png/.pdf` — основной рисунок.

Сопоставимые решения центральных переходов дают существенно разные η первой слабой группы. Результаты условны относительно модели; четыре группы s1–s4 не отождествлены с пятью атомными позициями Cu1–Cu5.

Для полного основного пересчета, из корневой папки пакета:

```bash
python powder_model/joint_fit.py
python powder_model/refine_exact_joint.py
python powder_model/central_control.py
python powder_model/joint_literature_sensitivity.py
python powder_model/finalize_controls.py
python powder_model/extent_examples.py
```

`joint_fit.py` создает `joint_results.json`; `refine_exact_joint.py` уточняет его решения и создает `exact_joint_refinements.json`; `central_control.py` использует это уточнение; `joint_literature_sensitivity.py` проверяет другой наборЯКР. `finalize_controls.py` читает сохраненные результаты этих расчетов, прямо пересчитывает центральные переходы и формирует итоговую таблицу/рисунок. Его можно запустить отдельно для повторного построения итогов по уже включенным в пакет JSON.

Библиотеки `*.npz` создаются автоматически при отсутствии. Они являются промежуточным кэшем и в пакет не включаются. Полный поиск требует существенно больше времени, чем чтение CSV.

Исторические диагностические семейства можно воспроизвести дополнительно:

```bash
python powder_model/powder_fit.py --set thesis2009 --extended
python powder_model/powder_fit.py --set publication2013
python powder_model/validate_and_plot.py
python powder_model/finalize_results.py
```

Последний скрипт также использует `joint_results.json`. Файлы `summary_*.json`, `fit_*.json`, `model_comparison.csv` и `joint_candidate_parameters.csv` относятся к промежуточным этапам и могут отличаться от последующего уточнения. Актуальный итог — `FINAL_MODEL_COMPARISON.csv`, с прописным `FINAL`.

## Выгрузки и аудит

- `opj_Spec05122016.csv`: исходная интегральная таблица, включая конфликтующие записи.
- `opj_Spec05122016_audit.csv`: статусы строк.
- `opj_spec20ms.csv`, `opj_spec500ms.csv`: обработанные спектры, по 8651 точке.
- `opj_T2a.csv`:35 сохраненных точек эха; исходные числа времени не изменены, масштабмкс применяется в `relaxation_analysis`.
- `july_integrals_audited.csv`: июльская сводная таблица с флагами повторений и смены масштаба.
- `july_raw_trace_audit.csv`: проверка 600 временных трасс.
- `opj_graph_metadata.json`, `opj_inventory.json`, `source_archive_manifest.json`: связи графиков, подписи и структура источников.
- `manuscript_window_ratios.csv`, `manuscript_window_points.csv`: описательные отношения интенсивностей и использованные точки.
- `manuscript_peak_summary.csv`, `manuscript_nqr_bounds.csv`, `manuscript_numeric_summary.json`: числовая сводка исходных массивов и ограничений первой стадии анализа.

`figure_opj_manuscript.png/.pdf` показывает исходные данные до модельного анализа; текущие рисунки статьи находятся в `powder_model` и `relaxation_analysis/export`.

`NOTES_DATA_RU.txt` — текущие пояснения редакции 2. `forensics_report.md` — исходный аудит до новых материалов и авторских уточнений, сохраненный для прослеживаемости; его прежние отложения расчета η/T2 не являются текущим итогом.

Новые источники и проверенные страницы:

- `new_evidence/thesis/THESIS_EVIDENCE_RU.txt`
- `new_evidence/papers/NQR_PUBLICATIONS_AND_COIL_RU.txt`

## Сборка рукописи

После наличия итоговых рисунков нужны python-docx 1.2.0 и установленный Pandoc:

```bash
python draft/build_document.py
```

Сборка читает `draft/manuscript_v2.md` и `draft/powder_results_v2.md`, создает `draft/manuscript_v2_complete.md` и записывает `output/Anilite_NMR_draft_ru.docx`.

Для кратностей опубликованной структуры отдельно выполняется:

```bash
python research_literature/reproduce_anilite_sites.py
```

CSV сохраняют вычислительные знаки; они не задают экспериментальную точность. Результаты разных моделей и определений амплитуды не следует смешивать в один набор физических параметров.

## Сведения о формате Origin

Разбор data list соответствует описанию Origin 7.0552 Juliusz Gonera: https://github.com/jgonera/openopj/blob/master/docs/opj_format.markdown . Связи графиков и аннотаций сверены с liborigin: https://sources.debian.org/src/labplot/2.5.0-2/liborigin/OriginAnyParser.cpp . Считыватель предназначен для данного проекта и не является универсальной заменой Origin.
