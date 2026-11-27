"""Reproduce the analytic steps used in efg_derivation.tex.

All checks are algebraic identities. They are independent of the archived
numerical Ewald implementation and do not fit experimental data.
Run with Python 3 and SymPy; outputs are written alongside this script.
"""
from pathlib import Path
import csv
import json
import sympy as sp

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "results"
OUT.mkdir(exist_ok=True)
steps = []
checks = []


def save(name, expression, meaning):
    steps.append({"name": name, "meaning": meaning,
                  "sympy": str(expression), "latex": sp.latex(expression)})


def verify(name, expression):
    if isinstance(expression, sp.MatrixBase):
        residual = expression.applyfunc(sp.simplify)
        passed = residual == sp.zeros(*residual.shape)
    else:
        residual = sp.simplify(expression)
        passed = residual == 0
    checks.append({"check": name, "passed": bool(passed),
                   "residual": str(residual)})
    assert passed, (name, residual)


# 1. Unscreened Coulomb potential: exact Cartesian differentiation.
x, y, z = sp.symbols("x y z", real=True)
rv = sp.Matrix([x, y, z])
radius = sp.sqrt(x*x + y*y + z*z)
coulomb_gradient = sp.Matrix([sp.diff(1/radius, a) for a in rv])
coulomb_hessian = sp.hessian(1/radius, rv)
kernel = (3*rv*rv.T - radius**2*sp.eye(3))/radius**5
save("coulomb_potential", 1/radius, "Unit-charge potential, excluding r=0.")
save("coulomb_gradient", coulomb_gradient, "First Cartesian derivatives.")
save("coulomb_hessian", coulomb_hessian, "Second Cartesian derivatives.")
save("coulomb_tensor_kernel", kernel, "Compact Cartesian Hessian.")
verify("coulomb_hessian_equals_tensor_kernel", coulomb_hessian-kernel)
verify("coulomb_kernel_is_traceless_away_from_source", sp.trace(kernel))

# 2. The real-space screened Hessian and the removed self Gaussian.
r, alpha = sp.symbols("r alpha", positive=True)
f = sp.erfc(alpha*r)/r
fp, fpp = sp.diff(f, r), sp.diff(f, r, 2)
aa = sp.simplify((fpp-fp/r)/r**2)
bb = sp.simplify(fp/r)
aa_explicit = (3*sp.erfc(alpha*r)/r**5
               + 6*alpha*sp.exp(-alpha**2*r**2)/(sp.sqrt(sp.pi)*r**4)
               + 4*alpha**3*sp.exp(-alpha**2*r**2)/(sp.sqrt(sp.pi)*r**2))
bb_explicit = (-sp.erfc(alpha*r)/r**3
               - 2*alpha*sp.exp(-alpha**2*r**2)/(sp.sqrt(sp.pi)*r**2))
for name, expr, meaning in [
    ("screened_potential", f, "Real-space part erfc(alpha*r)/r."),
    ("screened_first_radial_derivative", fp, "f'(r)."),
    ("screened_second_radial_derivative", fpp, "f''(r)."),
    ("screened_dyadic_coefficient_A", aa_explicit, "Coefficient of r_a*r_b."),
    ("screened_identity_coefficient_B", bb_explicit, "Coefficient of delta_ab."),
    ("screened_trace", sp.simplify(aa*r*r+3*bb), "Trace before reciprocal and self terms.")]:
    save(name, expr, meaning)
verify("screened_A_matches_radial_derivatives", aa-aa_explicit)
verify("screened_B_matches_radial_derivatives", bb-bb_explicit)
verify("screened_trace_equals_gaussian_source", aa*r*r+3*bb-4*alpha**3/sp.sqrt(sp.pi)*sp.exp(-alpha**2*r**2))
self_series = sp.series(sp.erf(alpha*r)/r, r, 0, 8)
self_hessian_diagonal = 2*self_series.removeO().coeff(r, 2)
save("smooth_self_potential_series", self_series,
     "Taylor expansion of the smooth self potential near the observation nucleus.")
save("smooth_self_hessian_diagonal", self_hessian_diagonal,
     "Hessian of the self Gaussian already present in the reciprocal sum.")
save("self_removal_hessian_diagonal", -self_hessian_diagonal,
     "Positive correction required to REMOVE that self contribution.")
verify("self_hessian_coefficient", self_hessian_diagonal+4*alpha**3/(3*sp.sqrt(sp.pi)))

# 3. Differentiation of the reciprocal potential.
kx, ky, kz, phase = sp.symbols("k_x k_y k_z phase", real=True)
kv = sp.Matrix([kx, ky, kz])
reciprocal_cosine = sp.cos(kv.dot(rv)+phase)
reciprocal_hessian = sp.hessian(reciprocal_cosine, rv)
save("reciprocal_cosine_hessian", reciprocal_hessian,
     "The minus sign in the reciprocal EFG follows from differentiating cosine twice.")
verify("reciprocal_hessian_sign", reciprocal_hessian+kv*kv.T*reciprocal_cosine)

# 4. EFG principal values and the mirror constraint.
vzz, eta = sp.symbols("V_ZZ eta", real=True)
vxx, vyy = (eta-1)*vzz/2, -(eta+1)*vzz/2
save("principal_values_from_Vzz_eta", sp.diag(vxx, vyy, vzz),
     "Traceless principal tensor reconstructed from V_ZZ and eta.")
verify("principal_values_trace", vxx+vyy+vzz)
verify("principal_values_asymmetry", vxx-vyy-eta*vzz)
vaa, vbb, vcc, vab, vac, vbc = sp.symbols("V_aa V_bb V_cc V_ab V_ac V_bc", real=True)
general_v = sp.Matrix([[vaa,vab,vac],[vab,vbb,vbc],[vac,vbc,vcc]])
mirror = sp.diag(1,-1,1)
mirror_difference = mirror*general_v*mirror.T-general_v
save("mirror_symmetry_difference", mirror_difference,
     "Reflection perpendicular to b requires V_ab=V_bc=0; it does not fix eta.")
verify("mirror_allowed_tensor", mirror_difference.subs({vab:0,vbc:0}))
tetra = [sp.Matrix(a) for a in [(1,1,1),(1,-1,-1),(-1,1,-1),(-1,-1,1)]]
tetra_sum = sp.zeros(3)
for a in tetra:
    rr = sp.sqrt(a.dot(a))
    tetra_sum += (3*a*a.T-rr**2*sp.eye(3))/rr**5
save("regular_tetrahedron_tensor", sp.simplify(tetra_sum),
     "Four identical point charges on an ideal tetrahedron give a zero tensor; eta is undefined.")
verify("ideal_tetrahedron_EFG_vanishes", tetra_sum)

# 5. Exact spin-3/2 quadrupolar Hamiltonian and its characteristic polynomial.
spin = sp.Rational(3,2)
ms = [spin-j for j in range(4)]
ip = sp.zeros(4)
for j in range(1,4):
    ip[j-1,j] = sp.sqrt(spin*(spin+1)-ms[j]*(ms[j]+1))
im = ip.T
ix, iy, iz = (ip+im)/2, (ip-im)/(2*sp.I), sp.diag(*ms)
nuq, nuz = sp.symbols("nu_Q nu_Z", real=True)
# charpoly uses an assumption-free generator; keep the comparison identical.
lam = sp.Symbol("lambda")
qmat = sp.simplify(nuq/6*(3*iz**2-spin*(spin+1)*sp.eye(4)+eta*(ix**2-iy**2)))
e2 = nuq**2*(3+eta**2)/12
save("spin_basis_m", sp.Matrix(ms), "Basis order |3/2>,|1/2>,|-1/2>,|-3/2>.")
save("I_plus", ip, "Dimensionless raising operator.")
save("I_X", ix, "Dimensionless spin X operator.")
save("I_Y", iy, "Dimensionless spin Y operator.")
save("I_Z", iz, "Dimensionless spin Z operator.")
save("H_Q_over_h", qmat, "Quadrupolar Hamiltonian in frequency units.")
save("H_Q_block", qmat.extract([0,2],[0,2]), "One of two equivalent 2x2 blocks.")
save("H_Q_characteristic_polynomial", sp.factor(qmat.charpoly(lam).as_expr()),
     "Each eigenvalue occurs twice.")
save("NQR_eigenvalues", sp.Matrix([-sp.sqrt(e2),-sp.sqrt(e2),sp.sqrt(e2),sp.sqrt(e2)]),
     "Ordered zero-field eigenvalues in frequency units.")
save("NQR_frequency", sp.Abs(nuq)*sp.sqrt(1+eta**2/3),
     "Energy separation divided by h.")
verify("spin_commutator_XY", ix*iy-iy*ix-sp.I*iz)
verify("quadrupole_H_squared", qmat*qmat-e2*sp.eye(4))
verify("quadrupole_characteristic_polynomial", qmat.charpoly(lam).as_expr()-(lam**2-e2)**2)

# 6. Axial high-field perturbation theory, with all intermediate state shifts.
c = sp.symbols("c", real=True)
sq = sp.sqrt(1-c*c)
efg_axis_spin = c*iz+sq*ix
qrot = sp.simplify(nuq/6*(3*efg_axis_spin**2-spin*(spin+1)*sp.eye(4)))
first = [sp.simplify(qrot[i,i]) for i in range(4)]
second = [sp.simplify(sum(qrot[i,j]*qrot[j,i]/(nuz*(ms[j]-ms[i]))
                          for j in range(4) if j!=i)) for i in range(4)]
ct_first = sp.simplify(first[2]-first[1])
ct_second = sp.factor(second[2]-second[1])
ct_expected = 3*nuq**2/(16*nuz)*(1-c*c)*(1-9*c*c)
save("axial_rotated_H_Q_over_h", qrot, "EFG Z axis at cos(theta)=c relative to B0.")
save("first_order_energy_shifts", sp.Matrix(first), "State-by-state first-order shifts, in Hz if nu is in Hz.")
save("second_order_energy_shifts", sp.Matrix(second), "State-by-state nondegenerate perturbation sums.")
save("CT_first_order_shift", ct_first, "Central first-order quadrupole shift vanishes.")
save("CT_second_order_shift", ct_second, "Central-transition second-order shift.")
save("CT_minimum_second_order", -nuq**2/(3*nuz), "Minimum at cos(theta)^2=5/9.")
save("CT_maximum_second_order", 3*nuq**2/(16*nuz), "Maximum at cos(theta)^2=0.")
verify("CT_first_order_zero", ct_first)
verify("CT_second_order_formula", ct_second-ct_expected)
verify("CT_minimum", ct_second.subs(c,sp.sqrt(5)/3)+nuq**2/(3*nuz))
verify("CT_maximum", ct_second.subs(c,0)-3*nuq**2/(16*nuz))
low_efg = {"nu_Q_MHz":1.5,"nu_Z_MHz":76.,
           "CT_second_order_min_kHz":float((-nuq**2/(3*nuz)).subs({nuq:1.5,nuz:76})*1000),
           "CT_second_order_max_kHz":float((3*nuq**2/(16*nuz)).subs({nuq:1.5,nuz:76})*1000),
           "satellite_first_order_outer_bounds_MHz":[74.5,77.5],
           "status":"Illustration at eta=0; not a measured or predicted Cu3 frequency."}

# 7. One-scale least-squares calibration is derived, not guessed.
scale = sp.symbols("a", real=True)
ps = sp.symbols("P_1:5", real=True)
fs = sp.symbols("f_1:5", real=True)
loss = sum((scale*p-f)**2 for p,f in zip(ps,fs))
scale_solution = sum(p*f for p,f in zip(ps,fs))/sum(p*p for p in ps)
save("common_scale_objective", loss,
     "Conditional empirical calibration for a fixed four-site permutation.")
save("common_scale_derivative", sp.diff(loss,scale), "First derivative with respect to the single scale.")
save("common_scale_solution", scale_solution, "Least-squares scale; positive for positive input frequencies.")
verify("common_scale_stationary_point", sp.diff(loss,scale).subs(scale,scale_solution))

report = {"status":"All symbolic identities passed; no experimental data were fitted.",
          "sympy_version":sp.__version__, "checks":checks,
          "steps":steps, "low_EFG_illustration":low_efg}
(OUT/"symbolic_algebra.json").write_text(json.dumps(report,indent=2)+"\n")
lines = ["ANILITE EFG/NQR DERIVATION: REPRODUCIBLE SYMBOLIC STEPS", ""]
for item in steps:
    lines.extend([item["name"],item["meaning"],item["sympy"],""])
lines.extend(["IDENTITY CHECKS",json.dumps(checks,indent=2),"LOW-EFG ILLUSTRATION",json.dumps(low_efg,indent=2)])
(OUT/"symbolic_algebra.txt").write_text("\n".join(lines)+"\n")
with (OUT/"symbolic_checks.csv").open("w",newline="") as stream:
    writer=csv.DictWriter(stream,fieldnames=["check","passed","residual"])
    writer.writeheader();writer.writerows(checks)
print(json.dumps({"checks_passed":len(checks),"saved_steps":len(steps),
                  "outputs":[str((OUT/name).relative_to(ROOT.parent)) for name in
                             ["symbolic_algebra.json","symbolic_algebra.txt","symbolic_checks.csv"]],
                  "low_EFG_illustration":low_efg},indent=2))
