"""Reproducible reanalysis of the stored Anilite T2a echo integrals.

TwoTau is interpreted as microseconds, informed by the pulse timing and
relaxation scales documented in the author's 2009 thesis. The original OPJ
graph label 'mS' is retained in the metadata as a conflicting annotation.
The fitted time is an effective echo-envelope time at t = 2*tau.
No spin-lattice or crystallographic-site assignment is made by this script.
"""

from pathlib import Path
import csv
import json
import numpy as np
from scipy.optimize import least_squares
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
SOURCE = ROOT.parent / 'data_analysis/export/opj_T2a.csv'
OUT = ROOT / 'export'
OUT.mkdir(parents=True, exist_ok=True)
a = np.genfromtxt(SOURCE, delimiter=',', names=True)
t = a['TwoTau']
mag = np.hypot(a['IntRe'], a['IntIm'])
phase = np.degrees(np.arctan2(a['IntIm'], a['IntRe']))

def fit_signal(y, mask, model):
    tx, yy = t[mask], y[mask]
    if model == 'single_exponential':
        evaluate = lambda p, x: p[0] * np.exp(-x / p[1])
        initial = [1500.0, 950.0]
        bounds = ([0, 1], [10000, 100000])
    else:
        evaluate = lambda p, x: p[0] * np.exp(-(x / p[1]) ** p[2])
        initial = [1340.0, 955.0, 1.25]
        bounds = ([0, 1, .1], [10000, 100000, 4])
    result = least_squares(lambda p: (evaluate(p, tx) - yy) / 1000,
                           initial, bounds=bounds)
    residual = evaluate(result.x, tx) - yy
    return {
        'n_points': int(mask.sum()),
        'amplitude': float(result.x[0]),
        'effective_time_us': float(result.x[1]),
        'beta': 1.0 if model == 'single_exponential' else float(result.x[2]),
        'rms_residual_arbitrary_units': float(np.sqrt(np.mean(residual**2))),
        'success': bool(result.success),
    }

masks = {
    'all_35': t > 0,
    'first_sweep_19': a['origin_row'] <= 19,
    'second_sweep_16': a['origin_row'] >= 20,
    'all_t_ge_100_us': t >= 100,
    'first_sweep_t_ge_100_us': (a['origin_row'] <= 19) & (t >= 100),
    'all_t_ge_200_us': t >= 200,
}
fits = []
for channel, y in [('complex_magnitude', mag), ('IntIm', a['IntIm']),
                   ('IntQ', a['IntQ'])]:
    for subset, mask in masks.items():
        for model in ['single_exponential', 'generalized_exponential']:
            fits.append({'channel': channel, 'subset': subset, 'model': model,
                         **fit_signal(y, mask, model)})

with (OUT / 't2_fit_sensitivity.csv').open('w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=list(fits[0]))
    writer.writeheader(); writer.writerows(fits)
with (OUT / 't2_envelope_points.csv').open('w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['origin_row', 'two_tau_us', 'IntRe', 'IntIm', 'IntQ',
                     'complex_magnitude', 'phase_degrees', 'sweep'])
    for i in range(len(a)):
        writer.writerow([int(a['origin_row'][i]), t[i], a['IntRe'][i],
                         a['IntIm'][i], a['IntQ'][i], mag[i], phase[i],
                         1 if a['origin_row'][i] <= 19 else 2])

main = next(f for f in fits if f['channel'] == 'complex_magnitude'
            and f['subset'] == 'all_35' and f['model'] == 'generalized_exponential')
single = next(f for f in fits if f['channel'] == 'complex_magnitude'
              and f['subset'] == 'all_35' and f['model'] == 'single_exponential')
summary = {
    'source': 'data_analysis/export/opj_T2a.csv; Anilite.OPJ, sheet T2a',
    'temperature_K': 5,
    'frequency_MHz': 51,
    'nominal_field_kOe': 45.7,
    'time_definition': 't = TwoTau = 2*tau; microseconds used in this analysis',
    'unit_evidence': 'Author 2009 thesis: pulse timing in microseconds and '
                     'NQR inverse T2 in ms^-1, Fig. 11, p.29 of version _5.',
    'original_graph_unit_label': 'mS (inconsistent with the adopted scale)',
    'echo_amplitude_definition': 'sqrt(IntRe**2 + IntIm**2); this differs from stored IntQ',
    'main_model': 'A * exp(-(t/T2e)**beta), zero offset, unweighted amplitude residuals',
    'main_fit': main,
    'single_exponential_fit': single,
    'phase_range_degrees': [float(phase.min()), float(phase.max())],
    'interpretation': 'T2e is the effective 1/e envelope time. beta is descriptive. '
                      'Sensitivity runs are not confidence intervals. '
                      'No microscopic relaxation mechanism is inferred.',
    'nqr_comparison': {'source': 'Evlampiev thesis 2009, Fig.11, p.29 (_5)',
                       'approx_inverse_T2_ms_minus1_lowT': [6, 9],
                       'approx_T2_us_lowT': [1000/9, 1000/6],
                       'numeric_precision': 'Read as a scale from a figure, not tabulated measurements'},
}
(OUT / 't2_summary.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2))

plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10,
                     'axes.spines.top': False, 'axes.spines.right': False})
fig, axs = plt.subplots(2, 1, figsize=(7.2, 6.2), sharex=True,
                        gridspec_kw={'height_ratios': [2.15, 1]}, layout='constrained')
colors = ['#0B6E8C', '#B95728']
for sweep, color, marker in [(1, colors[0], 'o'), (2, colors[1], 's')]:
    m = (a['origin_row'] <= 19) if sweep == 1 else (a['origin_row'] >= 20)
    axs[0].scatter(t[m]/1000, mag[m]/1000, color=color, marker=marker,
                   s=25, label=f'Data block {sweep}', zorder=3)
    axs[1].scatter(t[m]/1000, phase[m], color=color, marker=marker, s=23)
grid = np.linspace(0, 1300, 500)
axs[0].plot(grid/1000, main['amplitude'] / 1000 *
            np.exp(-(grid/main['effective_time_us'])**main['beta']),
            color='#202020', lw=1.5, label='Generalized exponential')
axs[0].plot(grid/1000, single['amplitude'] / 1000 *
            np.exp(-grid/single['effective_time_us']),
            color='#777777', lw=1.1, ls='--', label='Single exponential')
axs[0].set_ylabel(r'$|\mathrm{IntRe}+i\,\mathrm{IntIm}|$ ($10^3$ a.u.)')
axs[0].text(.97, .96, '5 K, 51 MHz, nominal 4.57 T\n'
            + rf'$T_{{2e}}={main["effective_time_us"]/1000:.3f}$ ms, '
            + rf'$\beta={main["beta"]:.2f}$',
            ha='right', va='top', transform=axs[0].transAxes, fontsize=10)
axs[0].legend(loc='lower left', frameon=False, fontsize=9)
axs[0].set_ylim(0, 1.6)
axs[1].set_ylabel('Integral phase (deg)')
axs[1].set_xlabel(r'Echo time $t=2\tau$ (ms; stored TwoTau interpreted as $\mu$s)')
axs[1].set_ylim(47, 70)
for ax in axs:
    ax.grid(axis='y', color='#DDDDDD', linewidth=.6)
    ax.set_xlim(0, 1.32)
for ext in ['png', 'pdf']:
    fig.savefig(OUT / f'figure_t2_reanalysis.{ext}', dpi=220)
plt.close(fig)
print(json.dumps(summary, ensure_ascii=False, indent=2))
