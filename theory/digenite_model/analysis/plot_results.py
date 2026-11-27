#!/usr/bin/env python3
"""Plot saved calculated results; never generate or refit scientific data here."""
from pathlib import Path
import csv
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = Path(__file__).resolve().parents[1]
REPO = BASE.parents[1]
FIGURES = BASE / "figures"
FIGURES.mkdir(exist_ok=True)
plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 11,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.titleweight": "bold", "axes.titlepad": 13,
    "savefig.dpi": 180, "pdf.fonttype": 42, "svg.fonttype": "none",
})
BLUE, RED, ORANGE, DARK = "#1665A7", "#B23A48", "#D18A21", "#202A35"

def rows(path):
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))

def save(fig, name):
    for ext in ("png", "pdf", "svg"):
        path = FIGURES / f"{name}.{ext}"
        fig.savefig(path, bbox_inches="tight", facecolor="white")
        if ext == 'svg':
            path.write_text('\n'.join(line.rstrip() for line in path.read_text().splitlines()) + '\n')
    plt.close(fig)

def frequency_comparison():
    continuous = json.loads((BASE / "results/continuous_charge_optimum.json").read_text())
    best = continuous["global_best"]
    observed = np.array(continuous["observations_MHz"])
    uniform = rows(BASE / "results/assignments_best6.csv")[0]
    ani = rows(REPO / "theory/results/efg/frequency_scaling_best_by_scenario.csv")[0]
    x = np.arange(1, 5)
    fig, ax = plt.subplots(figsize=(9.7, 6.1))
    ax.plot(x, [float(uniform[f"fitted_s{i}_MHz"]) for i in x],
            "^--", color=ORANGE, lw=1.6, ms=7,
            label="Дигенит, одинаковые заряды Cu; RMSE 8.71 МГц")
    ax.plot(x, best["actual_selected_predictions_MHz"], "s--", color=RED,
            lw=1.8, ms=7,
            label="Дигенит, оптимум непрерывной семьи зарядов; 5.43 МГц")
    ax.plot(x, [float(ani[f"predicted_s{i}_MHz"]) for i in x], "D--",
            color=BLUE, lw=1.8, ms=6,
            label="Анилит, прежняя базовая модель; 0.58 МГц")
    ax.plot(x, observed, "o-", color=DARK, lw=2.2, ms=7,
            label="Измеренные частоты ⁶³Cu при 4.2 K")
    for xi, yi in zip(x, observed):
        ax.annotate(f"{yi:.2f}", (xi, yi), xytext=(0, -19),
                    textcoords="offset points", ha="center", color=DARK, fontsize=10)
    ax.set(xticks=x, xticklabels=[f"s{i}" for i in x], ylim=(0, 26.4),
           ylabel="Частота ⁶³Cu, МГц",
           xlabel="Спектральные группы (назначения к Cu оптимизируются отдельно)",
           title="Относительные частоты: один общий масштаб на модель")
    ax.grid(axis="y", color="#DCE2E7", lw=.7)
    ax.legend(loc="upper left", frameon=False, fontsize=9.4)
    fig.text(.125, -.015,
             "Фиксированные исторические координаты; решётка точечных зарядов. "
             "Валентный ГЭП не вычислен.\n"
             "Сравниваются четыре группы; пятая позиция каждой модели может не входить в этот диапазон.",
             fontsize=9, color="#4C5865")
    save(fig, "frequency_model_comparison")

def geometry_sensitivity():
    data = rows(BASE / "results/sensitivities_uniform_charges.csv")
    linear = [r for r in data if r["variation"] == "uniform_linear_scale"
              and r["evaluation"] == "analytic_lambda_minus3" and r["site"] == "DgA"]
    sites = [f"Dg{s}" for s in "ABCDE"]
    maximum = [max(abs(float(r["relative_frequency_change"])) for r in data
                   if r["variation"] == "inversion_paired_axial_displacement_A"
                   and r["site"] == s) * 100 for s in sites]
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.7), gridspec_kw={"wspace": .35})
    x = [(float(r["value"]) - 1) * 100 for r in linear]
    y = [float(r["relative_frequency_change"]) * 100 for r in linear]
    axes[0].plot(x, y, "o-", color=BLUE, lw=2)
    axes[0].axhline(0, color="#A0AAB3", lw=.7)
    axes[0].axvline(0, color="#A0AAB3", lw=.7)
    axes[0].set(xlabel="Изменение всех линейных размеров, %",
                ylabel="Изменение частоты, %",
                title="Однородная деформация: λ⁻³")
    bars = axes[1].bar(sites, maximum, color=[BLUE, BLUE, BLUE, BLUE, ORANGE], width=.65)
    axes[1].bar_label(bars, labels=[f"{v:.2f}%" for v in maximum], padding=4, fontsize=10)
    axes[1].set(ylim=(0, 36), ylabel="Наибольшее |изменение частоты|, %",
                title="Парные осевые смещения ±0.01 Å")
    for ax in axes:
        ax.grid(axis="y", color="#E1E5E9", lw=.7)
        ax.set_axisbelow(True)
    fig.suptitle("Чувствительность геометрии — это не рассчитанная ν(T)",
                 fontsize=14, fontweight="bold", y=1.04)
    fig.text(.125, -.08,
             "Справа: максимум по отдельным смещениям каждой Cu-пары B–E; "
             "S фиксированы. Во всех этих вариантах η = 0.\n"
             "Физический температурный интервал и вероятности смещений не назначались.",
             fontsize=9, color="#4C5865")
    save(fig, "geometry_sensitivity")

def nmr_axial_control():
    folder = REPO / "theory/digenite_nmr_control"
    fig, axes = plt.subplots(2, 2, figsize=(11.2, 6.6), sharex="col",
                             gridspec_kw={"height_ratios": [3, 1], "hspace": .08, "wspace": .2})
    for col, series in enumerate(("500ms", "20ms")):
        axial = np.genfromtxt(folder / f"axial_exact1280_{series}.csv", delimiter=",", names=True)
        free = np.genfromtxt(folder / f"archived_free_eta_seed67_exact160_{series}.csv", delimiter=",", names=True)
        x = axial["frequency_MHz"]
        assert np.allclose(x, free["frequency_MHz"])
        used = axial["used_in_fit"] > .5
        scale = np.max(axial["observed"][used])
        upper, lower = axes[0, col], axes[1, col]
        upper.plot(x[used], axial["observed"][used]/scale, ".", ms=3.5,
                   color="#8F989F", label="Эксперимент")
        upper.plot(x, free["predicted"]/scale, color=BLUE, lw=1.6, label="Свободные η")
        upper.plot(x, axial["predicted"]/scale, color=RED, lw=1.5, label="Все η = 0")
        for source, color in ((free, BLUE), (axial, RED)):
            residual = np.where(used, source["residual"]/scale, np.nan)
            lower.plot(x, residual, color=color, lw=1.1)
        lower.axhline(0, color=DARK, lw=.7)
        for ax in (upper, lower):
            for centre in (75.941, 81.346):
                ax.axvspan(centre-.12, centre+.12, color="#DBDFE3", alpha=.6, lw=0)
            ax.set_xlim(72, 85)
        upper.set(ylim=(-.03, 1.075), title=f"Период повторения {series.replace('ms', ' мс')}")
        limit = max(np.max(np.abs(source["residual"][used]/scale)) for source in (free, axial)) * 1.08
        lower.set(xlabel="Частота, МГц", ylim=(-limit, limit), yticks=[-.2, 0, .2])
    axes[0, 0].set_ylabel("Интенсивность / максимум серии")
    axes[1, 0].set_ylabel("Остаток")
    axes[0, 0].legend(frameon=False, loc="upper left", fontsize=9)
    fig.suptitle("Совместный контроль аксиальности ГЭП: 5 K, 6.720566 T", fontsize=14, fontweight="bold", y=.99)
    fig.text(.125, -.025,
             "Измеренные ЯКР-частоты фиксированы; это тест η, а не полного спектра Donnay. "
             "Серые полосы исключены из подгонки.\n"
             "SSE: 0.539970 → 0.820887 (+52.02%); RMS: 0.033196 → 0.040930 (+23.30%). "
             "Критерий не является калиброванным χ².",
             fontsize=9, color="#4C5865")
    save(fig, "nmr_axial_control")

if __name__ == "__main__":
    frequency_comparison()
    geometry_sensitivity()
    nmr_axial_control()
    print(f"Saved three figures in PNG, PDF and SVG: {FIGURES}")
