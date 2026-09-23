"""Figure for the off-plane exploration note: degenerate GVM loci A/D/C of the
(f,s,f) interaction as functions of the azimuth phi (three panels, one measure each)."""
import csv, numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

COL = {'0': '#2a78d6', '45': '#eb6834', '90': '#1baf7a'}      # fixed categorical order A, D, C
LBL = {'0': 'A (0°)', '45': 'D (45°)', '90': 'C (90°)'}

rows = list(csv.DictReader(open('offplane_degenerate_loci.csv')))
epm = list(csv.DictReader(open('offplane_epm_locus.csv')))


def series(target, octant):
    """octant 1: theta < 90 (continuation of the manuscript's theta <= 90 cut); 2: mirror."""
    out = []
    if target == '45':
        for r in epm:
            phi = float(r['phi_deg']); th = float(r['theta_deg'])
            if octant == 1:
                out.append((phi, th, float(r['pump_nm']), float(r['abs_deff']), float(r['rho_max_mrad']), float(r['fom'])))
            else:
                out.append((phi, 180 - th, float(r['pump_nm']), float(r['abs_deff_mirror']), float(r['rho_max_mrad']), float(r['fom_mirror'])))
    else:
        for r in rows:
            if r['assign'] != 'fsf' or r['target'] != target:
                continue
            th = float(r['theta_deg'])
            if (octant == 1) == (th < 90):
                out.append((float(r['phi_deg']), th, float(r['pump_nm']), float(r['abs_deff']), float(r['rho_max_mrad']), float(r['fom'])))
    out.sort()
    return np.array(out)


fig, ax = plt.subplots(1, 4, figsize=(13, 3.3))
for t in ('0', '45', '90'):
    for oc, ls in ((1, '-'), (2, '--')):
        d = series(t, oc)
        lab = LBL[t] + (' , θ<90°' if oc == 1 else ' , mirror octant')
        ax[0].plot(d[:, 0], d[:, 1], ls, color=COL[t], lw=2, label=lab)
        ax[1].plot(d[:, 0], d[:, 2], ls, color=COL[t], lw=2)
        ax[2].plot(d[:, 0], d[:, 3], ls, color=COL[t], lw=2)
        ax[3].plot(d[:, 0], d[:, 5] * 1e3, ls, color=COL[t], lw=2)
for a, yl, title in zip(ax, ('cut angle θ (deg)', 'pump wavelength (nm)', '|d_eff| (pm/V)', 'd_eff/ρ (pm/V per rad)'),
                        ('polar angle of the locus', 'degenerate pump', 'nonlinearity (octant-split)', 'walk-off figure of merit')):
    a.set_xlabel('azimuth φ from the x–z plane (deg)'); a.set_ylabel(yl); a.set_title(title, fontsize=10)
    a.set_xlim(0, 90); a.grid(alpha=0.25); a.spines[['top', 'right']].set_visible(False)
ax[0].legend(fontsize=7, frameon=False, loc='lower left')
# in-plane reference marks
for t, y in (('0', 26.1), ('45', 27.7), ('90', 27.1)):
    ax[3].axhline(y, color=COL[t], lw=0.8, alpha=0.4)
ax[3].text(60, 29.5, 'thin lines: in-plane values', fontsize=7, color='#555')
fig.suptitle('BIBO (f,s,f) = x–z oeo interaction continued off the principal plane: degenerate GVM loci vs azimuth', fontsize=10)
fig.tight_layout()
fig.savefig('offplane_loci.png', dpi=160)
print('saved offplane_loci.png')
