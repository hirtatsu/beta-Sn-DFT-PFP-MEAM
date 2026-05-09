"""
05_parse_openmx_gamma.py (run locally on WSL)
Parses OpenMX SP (or Relax) logs on intel via ssh, computes surface energies,
and prints a comparison table against PFP.
"""
import json
import re
import subprocess
import sys

HARTREE_TO_EV = 27.211386245988
EV_A2_TO_MJ_M2 = 16021.766

SSH = ['ssh', 'intel']

# PFP reference gammas at matching thickness (from PFP convergence run)
PFP_GAMMA = {'100_L12': 834.5, '001_L8': 736.2, '110_L8': 737.5, '101_L12': 833.1}


def fetch(path):
    r = subprocess.run(SSH + [f'cat {path}'], capture_output=True, text=True)
    if r.returncode != 0:
        return None
    return r.stdout


def last_utot_hartree(log_text):
    vals = re.findall(r'Utot\s*=\s*(-?\d+\.\d+)', log_text)
    return float(vals[-1]) if vals else None


def get_cell_area(manifest, sysname):
    for s in manifest['slabs']:
        if s['sysname'] == sysname or s['sysname'].replace('_SP', '') == sysname.replace('_SP', '').replace('_Relax', ''):
            return s['area_A2'], s['n_atoms']
    raise KeyError(sysname)


def main(mode='SP'):
    # remote dirs
    rdir = {'SP': '~/work_dir/openmx_betaSn_surface',
            'Relax': '~/work_dir/openmx_betaSn_surface_relax'}[mode]

    man_text = fetch(f'~/work_dir/openmx_betaSn_surface/manifest.json')
    manifest = json.loads(man_text)

    # bulk reference (always from SP dir)
    bulk_log = fetch('~/work_dir/openmx_betaSn_surface/bulk_pfp_SP.log')
    if bulk_log is None:
        print("bulk SP log not found")
        return
    E_bulk_Ha = last_utot_hartree(bulk_log)
    if E_bulk_Ha is None:
        print("bulk SP log does not yet contain Utot")
        return
    E_bulk_eV = E_bulk_Ha * HARTREE_TO_EV
    E_bulk_per_atom = E_bulk_eV / 4.0
    print(f"Bulk SP (PFP lattice): Utot = {E_bulk_Ha:.6f} Ha  = {E_bulk_eV:.4f} eV")
    print(f"E_bulk_per_atom = {E_bulk_per_atom:.6f} eV\n")

    slabs = [('100', 12), ('001', 8), ('110', 8), ('101', 12)]
    rows = []
    suffix = '_SP' if mode == 'SP' else '_Relax'
    for tag, L in slabs:
        sysname = f'slab_{tag}_L{L}{suffix}'
        log_path = f'{rdir}/{sysname}.log'
        txt = fetch(log_path)
        if txt is None:
            print(f'{sysname}: log not found (not yet started or running)')
            rows.append({'sysname': sysname, 'status': 'pending'})
            continue
        E_slab_Ha = last_utot_hartree(txt)
        if E_slab_Ha is None:
            print(f'{sysname}: no Utot yet (still running?)')
            rows.append({'sysname': sysname, 'status': 'running'})
            continue
        # match manifest entry
        sp_name = f'slab_{tag}_L{L}_SP'
        entry = next(s for s in manifest['slabs'] if s['sysname'] == sp_name)
        A = entry['area_A2']
        n = entry['n_atoms']
        E_slab_eV = E_slab_Ha * HARTREE_TO_EV
        gamma = (E_slab_eV - n * E_bulk_per_atom) / (2.0 * A) * EV_A2_TO_MJ_M2
        rows.append({
            'sysname': sysname, 'tag': tag, 'layers': L, 'n_atoms': n,
            'area_A2': A, 'E_slab_Ha': E_slab_Ha, 'E_slab_eV': E_slab_eV,
            'gamma_mJ_m2': gamma,
            'pfp_gamma': PFP_GAMMA[f'{tag}_L{L}'],
            'diff': gamma - PFP_GAMMA[f'{tag}_L{L}'],
            'diff_pct': 100 * (gamma - PFP_GAMMA[f'{tag}_L{L}']) / PFP_GAMMA[f'{tag}_L{L}'],
        })

    print(f"\n{'='*78}")
    print(f"Surface energy comparison ({mode}): OpenMX/PBE  vs  PFP/PBE+D3")
    print(f"{'='*78}")
    print(f"{'face':>6} {'L':>3} {'nat':>4}  {'γ_DFT (mJ/m²)':>14}  {'γ_PFP':>8}  {'Δ':>8}  {'Δ %':>7}")
    print('-' * 78)
    for r in rows:
        if 'gamma_mJ_m2' not in r:
            print(f"{r['sysname']:>40}  {r.get('status','?')}")
            continue
        print(f"{r['tag']:>6} {r['layers']:>3} {r['n_atoms']:>4}  "
              f"{r['gamma_mJ_m2']:>14.1f}  {r['pfp_gamma']:>8.1f}  "
              f"{r['diff']:>+8.1f}  {r['diff_pct']:>+7.1f}")
    print('=' * 78)

    # save (write next to script; original was hardcoded under <PROJECT_ROOT>/surface_gb/)
    with open(f'gamma_{mode}.json', 'w') as f:
        json.dump({'bulk_E_per_atom_eV': E_bulk_per_atom, 'rows': rows}, f, indent=2)


if __name__ == '__main__':
    mode = sys.argv[1] if len(sys.argv) > 1 else 'SP'
    main(mode)
