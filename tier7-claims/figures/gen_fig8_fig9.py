"""Generate Fig 8 (8q per-patch K* vs random) and Fig 9 (dF scaling across n)."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({
    'font.family': 'serif',
    'font.size': 9,
    'axes.labelsize': 10,
    'xtick.labelsize': 8,
    'ytick.labelsize': 8,
    'legend.fontsize': 8,
    'figure.dpi': 300,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.15,
})

COL_WIDTH = 3.375
DBL_WIDTH = 6.75

# ── Figure 8: 8-qubit per-patch K* vs random (two-panel) ──────────────

fig8, (ax8a, ax8b) = plt.subplots(2, 1, figsize=(COL_WIDTH + 0.4, 5.2))

# Panel (a): Per-patch fidelity comparison
patches = ['A\n(q0-3)', 'B\n(q1-4)', 'C\n(q2-5)', 'D\n(q3-6)', 'E\n(q4-7)']
f_kstar  = [0.603, 0.582, 0.758, 0.709, 0.689]
f_random = [0.085, 0.117, 0.255, 0.116, 0.098]

# Bootstrap 95% CIs for K*
ci_kstar_lo = [0.505, 0.455, 0.595, 0.583, 0.593]
ci_kstar_hi = [0.620, 0.601, 0.772, 0.729, 0.749]
# Bootstrap 95% CIs for dF
dF = [fk - fr for fk, fr in zip(f_kstar, f_random)]
dF_ci_lo = [0.298, 0.314, 0.356, 0.398, 0.377]
dF_ci_hi = [0.535, 0.499, 0.646, 0.622, 0.639]

x = np.arange(len(patches))
w = 0.35

err_kstar = [[fk - lo for fk, lo in zip(f_kstar, ci_kstar_lo)],
             [hi - fk for fk, hi in zip(f_kstar, ci_kstar_hi)]]

bars_k = ax8a.bar(x - w/2, f_kstar, w, label='$K^*$ structured',
                   color='#2166ac', edgecolor='black', linewidth=0.5,
                   yerr=err_kstar, capsize=3, error_kw={'linewidth': 0.8})
bars_r = ax8a.bar(x + w/2, f_random, w, label='Random Pauli',
                   color='#d6604d', edgecolor='black', linewidth=0.5)

# dF annotations
for i in range(5):
    y_top = f_kstar[i] + (ci_kstar_hi[i] - f_kstar[i]) + 0.03
    ax8a.text(x[i], y_top,
              f'$\\Delta F = +{dF[i]:.2f}$',
              ha='center', va='bottom', fontsize=6.5, fontweight='bold')

ax8a.set_ylabel('Reconstruction fidelity $F$')
ax8a.set_xticks(x)
ax8a.set_xticklabels(patches, fontsize=7)
ax8a.set_ylim(0, 1.05)
ax8a.axhline(y=1/16, color='gray', linewidth=0.5, linestyle='--', alpha=0.5)
ax8a.text(4.35, 1/16 + 0.015, '$1/d$', fontsize=6, color='gray',
          ha='left', va='bottom')
ax8a.legend(loc='upper left', frameon=True, edgecolor='gray', fontsize=7)
ax8a.set_title('(a) 8-qubit $W_8$ state, per-patch (Rigetti Ankaa-3)', fontsize=8)

# Panel (b): Boundary consistency
boundaries = ['A-B', 'B-C', 'C-D', 'D-E']
dtr_w8 =      [0.194, 0.156, 0.157, 0.214]
dtr_product = [0.032, 0.041, 0.053, 0.035]
dtr_w8_ci =      [[0.195, 0.179, 0.186, 0.183],  # lo
                   [0.343, 0.338, 0.359, 0.349]]  # hi
dtr_prod_ci =    [[0.031, 0.039, 0.043, 0.038],  # lo
                   [0.072, 0.081, 0.090, 0.093]]  # hi

xb = np.arange(len(boundaries))
wb = 0.35

err_w8 = [[max(0, d - lo) for d, lo in zip(dtr_w8, dtr_w8_ci[0])],
          [max(0, hi - d) for d, hi in zip(dtr_w8, dtr_w8_ci[1])]]
err_pr = [[max(0, d - lo) for d, lo in zip(dtr_product, dtr_prod_ci[0])],
          [max(0, hi - d) for d, hi in zip(dtr_product, dtr_prod_ci[1])]]

ax8b.bar(xb - wb/2, dtr_w8, wb, label='$W_8$ (entangled)',
         color='#b2182b', edgecolor='black', linewidth=0.5,
         yerr=err_w8, capsize=3, error_kw={'linewidth': 0.8})
ax8b.bar(xb + wb/2, dtr_product, wb, label='$|{+}\\rangle^{\\otimes 8}$ (product)',
         color='#4393c3', edgecolor='black', linewidth=0.5,
         yerr=err_pr, capsize=3, error_kw={'linewidth': 0.8})

ax8b.set_ylabel('Boundary $D_{\\mathrm{tr}}$')
ax8b.set_xticks(xb)
ax8b.set_xticklabels(boundaries, fontsize=7)
ax8b.set_ylim(0, 0.45)
ax8b.legend(loc='upper right', frameon=True, edgecolor='gray', fontsize=6.5)
ax8b.set_title('(b) Patch-boundary consistency', fontsize=8)

fig8.tight_layout(h_pad=1.5)
fig8.savefig('fig8_8qubit_results.pdf')
fig8.savefig('fig8_8qubit_results.png')
print("Figure 8: fig8_8qubit_results.pdf")


# ── Figure 9: dF scaling across qubit counts ────────────────────────

fig9, ax9 = plt.subplots(figsize=(DBL_WIDTH, 3.2))

# Data points: n, dF, platform, marker
# n=2 Bell (IBM): dF = +0.052  (single run)
# n=3 W3 (Rigetti): dF = +0.009  (overcomplete)
# n=4 W4 (IBM, 4 runs): dF = +0.33 +/- 0.07
# n=4 W4 (Rigetti): dF = +0.248
# n=8 W8 compositional (Rigetti, mean of 5 patches): dF = +0.534

# IBM points
n_ibm = [2, 4]
dF_ibm = [0.052, 0.330]
dF_ibm_err = [None, 0.07]
err_ibm = [0 if e is None else e for e in dF_ibm_err]

ax9.errorbar(n_ibm, dF_ibm, yerr=err_ibm, fmt='s', color='#2166ac',
             markersize=9, markeredgecolor='black', markeredgewidth=0.5,
             capsize=4, linewidth=1.2, label='IBM ibm\\_fez', zorder=5)

# Rigetti points
n_rig = [3, 4, 8]
dF_rig = [0.009, 0.248, 0.534]
# For n=8, use std across 5 patches as uncertainty
dF_rig_err = [None, None, 0.051]
err_rig = [0 if e is None else e for e in dF_rig_err]

ax9.errorbar(n_rig, dF_rig, yerr=err_rig, fmt='o', color='#b2182b',
             markersize=9, markeredgecolor='black', markeredgewidth=0.5,
             capsize=4, linewidth=1.2, label='Rigetti Ankaa-3', zorder=5)

# Connect the scaling trend (IBM n=2 -> n=4, Rigetti n=3 -> n=4 -> n=8)
ax9.plot([2, 4], [0.052, 0.330], '--', color='#2166ac', linewidth=0.8, alpha=0.5, zorder=2)
ax9.plot([3, 4, 8], [0.009, 0.248, 0.534], '--', color='#b2182b', linewidth=0.8, alpha=0.5, zorder=2)

# Annotations — more breathing room at double width
ax9.annotate('Bell ($n{=}2$)', xy=(2, 0.052), xytext=(2.5, 0.14),
             fontsize=8, color='#2166ac',
             arrowprops=dict(arrowstyle='->', color='#2166ac', lw=0.8))
ax9.annotate('overcomplete ($M/4^n = 42\\%$)', xy=(3, 0.009), xytext=(3.6, -0.06),
             fontsize=7.5, color='#b2182b',
             arrowprops=dict(arrowstyle='->', color='#b2182b', lw=0.8))
ax9.annotate('compositional (5 patches)', xy=(8, 0.534), xytext=(6.3, 0.63),
             fontsize=8, color='#b2182b', fontweight='bold',
             arrowprops=dict(arrowstyle='->', color='#b2182b', lw=0.8))

# Sweet spot annotation
ax9.axvline(x=4, color='gray', linewidth=0.5, linestyle=':', alpha=0.5)
ax9.text(4.15, 0.60, '$n{=}4$ sweet spot ($M/4^n = 54\\%$)',
         fontsize=7.5, color='gray', va='top')

ax9.axhline(y=0, color='black', linewidth=0.3)
ax9.set_xlabel('Number of qubits $n$')
ax9.set_ylabel('$K^*$ advantage $\\Delta F$')
ax9.set_xticks([2, 3, 4, 5, 6, 7, 8])
ax9.set_xlim(1.5, 8.8)
ax9.set_ylim(-0.10, 0.72)
ax9.legend(loc='upper left', frameon=True, edgecolor='gray', fontsize=8)

fig9.tight_layout()
fig9.savefig('fig9_scaling.pdf')
fig9.savefig('fig9_scaling.png')
print("Figure 9: fig9_scaling.pdf")

print("\nFigures 8-9 generated.")
