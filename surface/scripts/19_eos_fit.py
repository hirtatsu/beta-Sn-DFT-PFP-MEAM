"""Fit γ(a) for (001) slab EOS scan and compare with BFGS/OptC5."""
import os
import json
import numpy as np
import matplotlib.pyplot as plt

HA_EV = 27.211386245988
EV_A2_TO_MJ_M2 = 16021.766

E_bulk_Ha = -284.157741623577
E_bulk_per_atom_eV = (E_bulk_Ha / 4.0) * HA_EV

# EOS scan points (from log.txt parsing)
data = [
    {'a': 5.82, 'E_Ha': -2273.156004576403, 'fmax': 7.26e-4, 'md': 221},
    {'a': 5.90, 'E_Ha': -2273.167681298297, 'fmax': 2.97e-4, 'md': 155},
    {'a': 5.97, 'E_Ha': -2273.172231377641, 'fmax': 1.75e-4, 'md': 77},
    {'a': 6.05, 'E_Ha': -2273.168561403050, 'fmax': 2.68e-4, 'md': 172},
    {'a': 6.12, 'E_Ha': -2273.162176950991, 'fmax': 2.97e-4, 'md': 114},
]

n_slab = 32

print(f"{'a (Å)':>7} {'E_slab (eV)':>13} {'dE (eV)':>8} {'A (Å²)':>8} "
      f"{'γ (mJ/m²)':>11} {'fmax':>10}")
print('-' * 66)

a_list, g_list, E_list = [], [], []
for d in data:
    a = d['a']
    E_slab_eV = d['E_Ha'] * HA_EV
    dE = E_slab_eV - n_slab * E_bulk_per_atom_eV
    A = a * a
    gamma = dE / (2 * A) * EV_A2_TO_MJ_M2
    print(f"{a:>7.2f} {E_slab_eV:>13.3f} {dE:>8.3f} {A:>8.3f} {gamma:>11.2f} {d['fmax']:>10.2e}")
    a_list.append(a); g_list.append(gamma); E_list.append(E_slab_eV)

a_arr = np.array(a_list)
g_arr = np.array(g_list)
E_arr = np.array(E_list)

# Fit γ vs a with 2nd-order polynomial
# Use all 5 points; also try 4 points (excluding a=5.82 which is poorly converged)
print()
print('--- quadratic fit γ(a) = p0 + p1*a + p2*a² ---')
p_all = np.polyfit(a_arr, g_arr, 2)
a_opt_all = -p_all[1] / (2 * p_all[0])
g_min_all = np.polyval(p_all, a_opt_all)
print(f"5 points : a_opt = {a_opt_all:.4f} Å,  γ_min = {g_min_all:.1f} mJ/m²")

mask4 = a_arr != 5.82
p4 = np.polyfit(a_arr[mask4], g_arr[mask4], 2)
a_opt_4 = -p4[1] / (2 * p4[0])
g_min_4 = np.polyval(p4, a_opt_4)
print(f"4 points : a_opt = {a_opt_4:.4f} Å,  γ_min = {g_min_4:.1f} mJ/m²")

# Also fit E(a) and locate minimum (cell relaxation proper)
print()
print('--- quadratic fit E(a) ---')
pE_all = np.polyfit(a_arr, E_arr, 2)
a_opt_E_all = -pE_all[1] / (2 * pE_all[0])
E_min_all = np.polyval(pE_all, a_opt_E_all)
# γ at fitted minimum of E, using fitted a_opt
g_at_E_min = (E_min_all - n_slab * E_bulk_per_atom_eV) / (2 * a_opt_E_all**2) * EV_A2_TO_MJ_M2
print(f"5 points : a_opt = {a_opt_E_all:.4f} Å,  E_min = {E_min_all:.3f} eV,  γ = {g_at_E_min:.1f} mJ/m²")

pE_4 = np.polyfit(a_arr[mask4], E_arr[mask4], 2)
a_opt_E_4 = -pE_4[1] / (2 * pE_4[0])
E_min_4 = np.polyval(pE_4, a_opt_E_4)
g_at_E_min_4 = (E_min_4 - n_slab * E_bulk_per_atom_eV) / (2 * a_opt_E_4**2) * EV_A2_TO_MJ_M2
print(f"4 points : a_opt = {a_opt_E_4:.4f} Å,  E_min = {E_min_4:.3f} eV,  γ = {g_at_E_min_4:.1f} mJ/m²")

# Reference BFGS (atom-only @ a=5.97) and OptC5 (broken symmetry)
bfgs_gamma = 550.5
optc5_gamma = 501.0

# Figure
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.5))
a_fit = np.linspace(5.78, 6.15, 100)
ax1.plot(a_arr, g_arr, 'o', color='#3366cc', ms=9, label='EOS data (atom-only BFGS at fixed a)')
ax1.plot(a_fit, np.polyval(p_all, a_fit), '-', color='#3366cc', alpha=0.5,
         label=f'2nd fit (5 pts): min at {a_opt_all:.3f} Å, γ={g_min_all:.1f}')
ax1.plot(a_fit, np.polyval(p4, a_fit), '--', color='#dc3912', alpha=0.7,
         label=f'2nd fit (4 pts, excl. a=5.82): min at {a_opt_4:.3f} Å, γ={g_min_4:.1f}')
ax1.axvline(5.97, color='k', ls=':', alpha=0.5, label='DFT bulk a=5.97')
ax1.axhline(bfgs_gamma, color='green', ls=':', alpha=0.5, label=f'BFGS γ={bfgs_gamma:.1f}')
ax1.axhline(optc5_gamma, color='purple', ls=':', alpha=0.5, label=f'OptC5 γ={optc5_gamma:.1f}')
ax1.set_xlabel('in-plane lattice a = b (Å)')
ax1.set_ylabel('γ (mJ/m²)')
ax1.set_title('β-Sn (001) γ vs a (EOS scan)')
ax1.legend(fontsize=8)
ax1.grid(alpha=0.3)

ax2.plot(a_arr, E_arr, 'o', color='#3366cc', ms=9)
ax2.plot(a_fit, np.polyval(pE_all, a_fit), '-', color='#3366cc', alpha=0.5,
         label=f'E(a) fit 5pts: min at {a_opt_E_all:.3f} Å')
ax2.plot(a_fit, np.polyval(pE_4, a_fit), '--', color='#dc3912', alpha=0.7,
         label=f'E(a) fit 4pts: min at {a_opt_E_4:.3f} Å')
ax2.axvline(5.97, color='k', ls=':', alpha=0.5, label='DFT bulk a=5.97')
ax2.set_xlabel('in-plane lattice a = b (Å)')
ax2.set_ylabel('E_slab (eV)')
ax2.set_title('β-Sn (001) E(a) and fitted minimum')
ax2.legend(fontsize=8)
ax2.grid(alpha=0.3)

plt.tight_layout()
outp = os.path.expanduser('~/projects/beta-sn-elasticity/surface_gb/figures/slab_001_EOS.png')
plt.savefig(outp, dpi=150, bbox_inches='tight')
print(f"\nsaved -> {outp}")
