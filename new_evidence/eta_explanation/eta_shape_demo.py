"""Illustrate EFG-asymmetry sensitivity at fixed NQR frequency.

This is a forward-model demonstration, not an additional experimental fit.
Uses the exact spin-3/2 Hamiltonian and transverse-RF matrix elements.
"""
from pathlib import Path
import json
import numpy as np
from scipy.ndimage import gaussian_filter1d
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
spin = 1.5
m = np.arange(spin, -spin - 1, -1)
jp = np.zeros((4, 4), complex)
for j in range(1, 4):
    jp[j-1, j] = np.sqrt(spin*(spin+1) - m[j]*(m[j]+1))
ix, iy, iz = (jp+jp.T)/2, (jp-jp.T)/(2j), np.diag(m)
operators = np.array([ix, iy, iz])
q0 = 3*iz@iz - spin*(spin+1)*np.eye(4)
qe = ix@ix - iy@iy
u = (np.arange(256)+.5)/256
phi = (np.arange(512)+.5)/512*np.pi/2
u, phi = np.meshgrid(u, phi, indexing='ij')
directions = np.column_stack([
    np.sqrt(1-u.ravel()**2)*np.cos(phi.ravel()),
    np.sqrt(1-u.ravel()**2)*np.sin(phi.ravel()),
    u.ravel(),
])
x = np.arange(-2.5, 1.5001, .0025)
step = x[1]-x[0]
nu_z, nu_nqr, sigma = 76., 16.18, .03
curves, metadata = [], []
for eta in (0., .3, .6, .9):
    nu_q = nu_nqr/np.sqrt(1+eta**2/3)
    ham = -nu_z*np.einsum('na,aij->nij', directions, operators) + nu_q/6*(q0+eta*qe)
    values, vectors = np.linalg.eigh(ham)
    elements = np.einsum('nki,akl,nlj->naij', vectors.conj(), operators, vectors, optimize=True)
    parallel = np.einsum('na,naij->nij', directions, elements)
    weight = (np.sum(abs(elements[:, :, 1, 2])**2, axis=1)-abs(parallel[:, 1, 2])**2)/2
    offsets = values[:, 2]-values[:, 1]-nu_z
    position = (offsets-x[0])/step
    lower = np.floor(position).astype(int)
    fraction = position-lower
    profile = np.zeros_like(x)
    keep = (lower>=0)&(lower<len(x)-1)
    np.add.at(profile, lower[keep], weight[keep]*(1-fraction[keep]))
    np.add.at(profile, lower[keep]+1, weight[keep]*fraction[keep])
    profile = gaussian_filter1d(profile, sigma/step)
    profile /= np.trapezoid(profile, x)
    curves.append(profile)
    metadata.append({'eta':eta, 'nu_Q_MHz':float(nu_q),
                     'central_transition_extent_MHz':[float(offsets.min()),float(offsets.max())]})

plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.spines.top':False,'axes.spines.right':False})
fig, ax = plt.subplots(figsize=(7.5, 4.6), constrained_layout=True)
for row, profile, color in zip(metadata, curves, ['#264653','#2a9d8f','#d49a20','#c44e52']):
    ax.plot(x, profile, lw=1.9, color=color,
            label=rf'$\eta={row["eta"]:.1f}$; $\nu_Q={row["nu_Q_MHz"]:.2f}$ MHz')
ax.set(xlim=(-1.9, 1.0), xlabel=r'Frequency offset from $\nu_Z$ (MHz)',
       ylabel='Unit-area central-transition intensity')
ax.set_title(r'Different powder shapes at the same NQR frequency', loc='left', fontweight='bold', pad=27)
ax.text(0, 1.025, r'$^{63}$Cu, $I=3/2$; $\nu_{\mathrm{NQR}}=16.18$ MHz, $\nu_Z=76.00$ MHz, $\sigma=0.03$ MHz',
        transform=ax.transAxes, fontsize=10)
ax.legend(frameon=False, fontsize=9.6, loc='upper right')
ax.grid(alpha=.15)
fig.savefig(ROOT/'eta_shape_demo.png',dpi=200)
fig.savefig(ROOT/'eta_shape_demo.pdf')
np.savetxt(ROOT/'eta_shape_demo.csv',np.column_stack([x,*curves]),delimiter=',',
           header='offset_MHz,eta_0p0,eta_0p3,eta_0p6,eta_0p9',comments='')
report={'purpose':'Forward-model illustration only; not additional experimental data or fit.',
        'nu_NQR_MHz':nu_nqr,'nu_Z_MHz':nu_z,'gaussian_sigma_MHz':sigma,
        'orientation_grid_in_octant':[256,512],'models':metadata}
(ROOT/'eta_shape_demo.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report))
