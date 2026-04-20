"""Generate all companion manuscript figures."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

# K* Verification style defaults
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

COL_WIDTH = 3.375  # single-column width (inches)
DBL_WIDTH = 6.75   # double-column width (inches)

# ── Figure 1: Head-to-head fidelity comparison ──────────────────────

fig1, ax1 = plt.subplots(figsize=(COL_WIDTH + 0.8, 2.8))

states = [r'$|{+}\rangle^{\otimes 4}$', r'Bell$^*$', r'$W_4$', r'$\mathrm{GHZ}_4$']
# 4-run paper values (Table I): W4 mean +/- std over 4 runs
f_struct = [0.996, 0.518, 0.872, 0.496]
f_rand   = [0.980, 0.465, 0.54,  0.498]
# Use exact JSON advantage values for annotations (not computed from rounded fidelities)
# Product: 0.01589, Bell: 0.05245, W4 4-run mean: 0.32970, GHZ: -0.00250
delta_f  = [+0.016, +0.052, +0.330, -0.002]
err_struct = [None, None, 0.021, None]   # W4: 4-run std; others single-run
err_rand   = [None, None, 0.08,  None]

x = np.arange(len(states))
w = 0.32

# Convert None errors to NaN so matplotlib skips the cap entirely
err_s = [e if e is not None else float('nan') for e in err_struct]
err_r = [e if e is not None else float('nan') for e in err_rand]

bars1 = ax1.bar(x - w/2, f_struct, w, label='Structured ($K^*$)',
                color='#2166ac', edgecolor='black', linewidth=0.5,
                yerr=err_s, capsize=3, error_kw={'linewidth': 0.8})
bars2 = ax1.bar(x + w/2, f_rand, w, label='Random',
                color='#d6604d', edgecolor='black', linewidth=0.5,
                yerr=err_r, capsize=3, error_kw={'linewidth': 0.8})

# Delta F annotations
for i, (d, fs) in enumerate(zip(delta_f, f_struct)):
    sign = '+' if d >= 0 else ''
    weight = 'bold' if abs(d) > 0.05 else 'normal'
    es = err_s[i] if not np.isnan(err_s[i]) else 0
    er = err_r[i] if not np.isnan(err_r[i]) else 0
    y_top = max(f_struct[i] + es, f_rand[i] + er) + 0.035
    ax1.text(x[i], y_top, f'$\\Delta F={sign}{d:.3f}$',
             ha='center', va='bottom', fontsize=7, fontweight=weight)

ax1.set_ylabel('Reconstruction fidelity $F$')
ax1.set_xticks(x)
ax1.set_xticklabels(states)
ax1.set_ylim(0, 1.25)
ax1.legend(loc='upper center', frameon=True, edgecolor='gray',
           ncol=2, bbox_to_anchor=(0.65, 1.0), columnspacing=0.8)
ax1.axhline(y=0.5, color='gray', linewidth=0.5, linestyle='--', alpha=0.5)
ax1.set_xlabel(r'Prepared state (IBM $\mathtt{ibm\_fez}$, $n{=}4$; $^*$Bell: $n{=}2$)')

fig1.tight_layout()
fig1.savefig('fig1_fidelity.pdf')
fig1.savefig('fig1_fidelity.png')
print("Figure 1: fig1_fidelity.pdf")

# ── Figure 2: Krawtchouk spectral correspondence ────────────────────
#
# Panel (a): Hamming-shell to Pauli-weight mapping
#   Left column: lattice shells |m|^2 = 0..5, with shell multiplicities r_4(k)
#   Right column: five Pauli weight classes w=0..4, with operator counts M_w
#   Arrows connect shells to the weight classes they populate (via mod-2 reduction)
#
# Panel (b): Gram eigenvalue spectrum
#   Five bars showing lambda_w = 16 * c_w / C(4,w), colored by weight class
#   Annotations: saturated vs subsampled, operator budget M_w / A_w

from matplotlib.patches import FancyArrowPatch
import matplotlib.gridspec as gridspec

fig2 = plt.figure(figsize=(COL_WIDTH + 0.6, 6.8))
gs = gridspec.GridSpec(2, 1, height_ratios=[1.15, 1], hspace=0.35)
ax2a = fig2.add_subplot(gs[0])
ax2b = fig2.add_subplot(gs[1])

# ── Panel (a): Shell → Weight mapping (horizontal layout) ──

ax2a.set_xlim(-0.3, 5.3)
ax2a.set_ylim(-0.6, 6.7)
ax2a.axis('off')
ax2a.set_title('(a) Lattice shells $\\to$ Pauli weight classes',
               fontsize=9, pad=6)

# Shell data: |m|^2 = k, multiplicity r_4(k)
shells = [
    (0, 1,   {0: 1}),
    (1, 8,   {1: 8}),
    (2, 24,  {2: 24}),
    (3, 32,  {1: 4, 3: 24, 2: 4}),
    (4, 24,  {0: 8, 2: 16}),
    (5, 48,  {1: 0, 2: 8, 3: 24, 4: 16}),
]

# Weight class data
weight_labels = ['$w{=}0$', '$w{=}1$', '$w{=}2$', '$w{=}3$', '$w{=}4$']
c_w = [9, 56, 24, 32, 16]
M_w = [1, 12, 54, 54, 16]
A_w = [1, 12, 54, 108, 81]
w_colors = ['#4a1486', '#2166ac', '#4393c3', '#f4a582', '#b2182b']

# Draw shells (left column)
shell_x = 0.6
for i, (k, r, _) in enumerate(shells):
    y = 5.5 - i * 0.9
    ax2a.add_patch(plt.Rectangle((shell_x - 0.5, y - 0.22), 1.0, 0.44,
                                  facecolor='#f0f0f0', edgecolor='#666666',
                                  linewidth=0.8, zorder=2))
    ax2a.text(shell_x, y, f'$|m|^2\\!={k}$',
              ha='center', va='center', fontsize=7.5, zorder=3)
    ax2a.text(shell_x, y - 0.42, f'$r={r}$',
              ha='center', va='top', fontsize=6, color='#555555')

# Draw weight classes (right column)
wt_x = 4.5
wt_ys = [5.5, 4.6, 3.7, 2.8, 1.9]
for i, (label, mw, aw) in enumerate(zip(weight_labels, M_w, A_w)):
    y = wt_ys[i]
    ax2a.add_patch(plt.Rectangle((wt_x - 0.6, y - 0.22), 1.2, 0.44,
                                  facecolor=w_colors[i], edgecolor='black',
                                  linewidth=0.8, alpha=0.25, zorder=2))
    ax2a.text(wt_x, y, label,
              ha='center', va='center', fontsize=8, fontweight='bold', zorder=3)
    sat_mark = '$\\checkmark$' if mw == aw else f'{mw}/{aw}'
    ax2a.text(wt_x, y - 0.42, f'$M_w$={sat_mark}',
              ha='center', va='top', fontsize=6, color='#333333')

# Arrows: shell → weight class
shell_ys = [5.5, 4.6, 3.7, 2.8, 1.9, 1.0]

def draw_arrow(sy_idx, wy_idx, rad=0.0, lw=0.6, alpha=0.5):
    sy = shell_ys[sy_idx]
    wy = wt_ys[wy_idx]
    arr = FancyArrowPatch((shell_x + 0.5, sy),
                          (wt_x - 0.6, wy),
                          arrowstyle='->', color=w_colors[wy_idx],
                          linewidth=lw, alpha=alpha,
                          connectionstyle=f'arc3,rad={rad}',
                          shrinkA=2, shrinkB=2, zorder=1)
    ax2a.add_patch(arr)

draw_arrow(0, 0, rad=0.0, lw=1.4, alpha=0.7)
draw_arrow(1, 1, rad=0.0, lw=1.4, alpha=0.7)
draw_arrow(2, 2, rad=0.0, lw=1.4, alpha=0.7)
draw_arrow(3, 1, rad=0.12, lw=0.7, alpha=0.4)
draw_arrow(3, 2, rad=0.08, lw=0.7, alpha=0.4)
draw_arrow(3, 3, rad=0.0, lw=1.4, alpha=0.7)
draw_arrow(4, 0, rad=-0.18, lw=0.7, alpha=0.4)
draw_arrow(4, 2, rad=0.08, lw=0.7, alpha=0.4)
draw_arrow(5, 2, rad=0.12, lw=0.7, alpha=0.4)
draw_arrow(5, 3, rad=0.08, lw=0.7, alpha=0.4)
draw_arrow(5, 4, rad=0.0, lw=1.4, alpha=0.7)

# Column headers — evenly space: header, subtitle, first box top (5.72)
ax2a.text(shell_x, 6.35, 'Lattice shells', ha='center', va='bottom',
          fontsize=8, fontweight='bold', color='#444444')
ax2a.text(shell_x, 6.03, '$m \\in \\mathbb{Z}^4$', ha='center', va='center',
          fontsize=7, color='#666666')
ax2a.text(wt_x, 6.35, 'Pauli weights', ha='center', va='bottom',
          fontsize=8, fontweight='bold', color='#444444')
ax2a.text(wt_x, 6.03, '$w(\\bar{m})$', ha='center', va='center',
          fontsize=7, color='#666666')

# Mod-2 arrow label
ax2a.text(2.55, 6.35, '$\\mathrm{mod}\\;2$', ha='center', va='bottom',
          fontsize=7.5, color='#888888', fontstyle='italic')
ax2a.annotate('', xy=(3.6, 6.3), xytext=(1.5, 6.3),
              arrowprops=dict(arrowstyle='->', color='#888888', lw=1.0))

# Cutoff annotation (below the r=48 label of the |m|^2=5 shell)
ax2a.plot([shell_x - 0.65, shell_x + 0.65], [0.25, 0.25],
          color='red', linewidth=1.0, linestyle='--', alpha=0.7)
ax2a.text(shell_x, 0.02, '$K^*{=}5$ cutoff', ha='center', va='top',
          fontsize=7, color='red', fontstyle='italic')

# ── Panel (b): Gram eigenvalue spectrum ──

weights = [0, 1, 2, 3, 4]
lambda_w = [144, 224, 64, 128, 256]
binom_w = [1, 4, 6, 4, 1]
pct_trace = [6.6, 40.9, 17.5, 23.4, 11.7]
saturated = [True, True, True, False, False]

ax2b.set_title('(b) Gram eigenvalues $\\lambda_w$', fontsize=9, pad=6)

bars = ax2b.bar(weights, lambda_w, color=w_colors, edgecolor='black',
                linewidth=0.6, alpha=0.8, width=0.6)

for i, (lam, bn, pct, sat) in enumerate(zip(lambda_w, binom_w, pct_trace, saturated)):
    ax2b.text(i, lam + 6, f'{lam}',
              ha='center', va='bottom', fontsize=7.5, fontweight='bold')
    mult_label = f'$\\times\\binom{{4}}{{{i}}}$={bn}'
    ax2b.text(i, lam * 0.5, mult_label,
              ha='center', va='center', fontsize=6, color='white',
              fontweight='bold')
    # Mass percentage: left of label for rightmost bar, right otherwise
    if i == 4:  # 256 at right edge — place % to the left
        ax2b.text(i - 0.22, lam + 6, f'({pct:.1f}%)',
                  ha='right', va='bottom', fontsize=5.5, color='#666666')
    else:
        pct_offset = 0.15 if lam < 100 else 0.22
        ax2b.text(i + pct_offset, lam + 6, f'({pct:.1f}%)',
                  ha='left', va='bottom', fontsize=5.5, color='#666666')
    marker = '$\\bullet$' if sat else '$\\circ$'
    ax2b.text(i, -16, marker, ha='center', va='top', fontsize=7)

ax2b.set_xlabel('Pauli weight $w$')
ax2b.set_ylabel('Eigenvalue $\\lambda_w$')
ax2b.set_xticks(weights)
ax2b.set_xlim(-0.5, 4.5)
ax2b.set_ylim(-10, 290)

ax2b.text(0.02, 0.97, '$\\bullet$ saturated  $\\circ$ subsampled',
          fontsize=6.5, ha='left', va='top', transform=ax2b.transAxes,
          bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                    edgecolor='gray', alpha=0.8))
fig2.savefig('fig2_spectral_correspondence.pdf')
fig2.savefig('fig2_spectral_correspondence.png')
print("Figure 2: fig2_spectral_correspondence.pdf")

# ── Figure 3: Compositional architecture ────────────────────────────

fig3, ax3 = plt.subplots(figsize=(COL_WIDTH, 1.6))
ax3.set_xlim(-0.5, 8.5)
ax3.set_ylim(-1.2, 3.5)
ax3.set_aspect('equal')
ax3.axis('off')

# Draw 8 qubits
for i in range(8):
    circle = plt.Circle((i, 0), 0.3, fill=True, facecolor='white',
                         edgecolor='black', linewidth=1.0, zorder=3)
    ax3.add_patch(circle)
    ax3.text(i, 0, f'$q_{i}$', ha='center', va='center', fontsize=7, zorder=4)

# Draw connections between adjacent qubits
for i in range(7):
    ax3.plot([i + 0.3, i + 0.7], [0, 0], 'k-', linewidth=0.8, zorder=1)

# Draw 4-qubit patches with color groups
color_map = {
    0: '#2166ac',  # blue
    1: '#b2182b',  # red
    2: '#e08214',  # orange (colorblind-safe, replaces green)
    3: '#762a83',  # purple
}
patches = [
    (0, 3, 0, 'Patch 1'),
    (1, 4, 1, 'Patch 2'),
    (2, 5, 2, 'Patch 3'),
    (3, 6, 3, 'Patch 4'),
    (4, 7, 0, 'Patch 5'),
]

for start, end, color_group, label in patches:
    y_offset = 1.0 + color_group * 0.55
    color = color_map[color_group]
    # Bracket
    ax3.plot([start - 0.35, start - 0.35], [0.4, y_offset],
             color=color, linewidth=1.0, alpha=0.7)
    ax3.plot([end + 0.35, end + 0.35], [0.4, y_offset],
             color=color, linewidth=1.0, alpha=0.7)
    ax3.plot([start - 0.35, end + 0.35], [y_offset, y_offset],
             color=color, linewidth=1.5, alpha=0.7)
    ax3.text((start + end) / 2, y_offset + 0.15, label,
             ha='center', va='bottom', fontsize=6, color=color, fontweight='bold')

# Legend for color groups
for cg in range(4):
    ax3.plot([], [], color=color_map[cg], linewidth=2,
             label=f'Color {cg+1}')
ax3.legend(fontsize=5.5, loc='lower right', ncol=4,
           frameon=True, edgecolor='gray', borderpad=0.3,
           columnspacing=0.8)

fig3.tight_layout()
fig3.text(0.5, 0.01, '116 circuits (4 colors $\\times$ 29 bases), any chain length',
          ha='center', va='bottom', fontsize=6.5, style='italic')
fig3.savefig('fig3_compositional.pdf')
fig3.savefig('fig3_compositional.png')
print("Figure 3: fig3_compositional.pdf")

# ── Figure 4: SOTA strategy comparison ──────────────────────────────

fig4, ax4 = plt.subplots(figsize=(DBL_WIDTH, 3.4))

# Data from Table V (Table tab:sota in manuscript)
strategies = [
    '$K^*$', 'D-opt', 'E-opt', 'A-opt', 'WA-rand',
    'DR-shadow', 'Adapt(3R)', 'Adapt(5R)', 'Unif-rand'
]
# Per state-class mean fidelities
f_product = [0.991, 0.987, 0.992, 0.993, 0.888, 0.952, 0.120, 0.034, 0.711]
f_W       = [0.932, 0.811, 0.833, 0.833, 0.414, 0.567, 0.096, 0.092, 0.291]
f_GHZ     = [0.502, 0.225, 0.225, 0.229, 0.352, 0.114, 0.139, 0.075, 0.215]
f_Haar    = [0.397, 0.398, 0.317, 0.316, 0.242, 0.066, 0.062, 0.062, 0.063]
f_Mixed   = [0.513, 0.502, 0.466, 0.462, 0.391, 0.227, 0.224, 0.224, 0.227]
f_overall = [0.501, 0.479, 0.429, 0.428, 0.347, 0.198, 0.140, 0.133, 0.179]

# Standard deviations across 20 trials (SEM = std/sqrt(20))
n_trials_sota = 20
s_product = np.array([0.001, 0.002, 0.002, 0.002, 0.210, 0.015, 0.237, 0.093, 0.290])
s_W       = np.array([0.012, 0.015, 0.018, 0.017, 0.224, 0.167, 0.052, 0.115, 0.259])
s_GHZ     = np.array([0.017, 0.020, 0.013, 0.017, 0.154, 0.063, 0.202, 0.105, 0.184])
s_Haar    = np.array([0.018, 0.018, 0.016, 0.016, 0.060, 0.062, 0.015, 0.015, 0.060])
s_Mixed   = np.array([0.015, 0.015, 0.014, 0.013, 0.050, 0.020, 0.022, 0.020, 0.050])
sem_product = s_product / np.sqrt(n_trials_sota)
sem_W       = s_W / np.sqrt(n_trials_sota)
sem_GHZ     = s_GHZ / np.sqrt(n_trials_sota)
sem_Haar    = s_Haar / np.sqrt(n_trials_sota)
sem_Mixed   = s_Mixed / np.sqrt(n_trials_sota)

xs = np.arange(len(strategies))
bw = 0.16

cat_colors = {
    'Product': '#a6cee3',
    'W':       '#2166ac',
    'GHZ':     '#b2182b',
    'Haar':    '#f4a582',
    'Mixed':   '#762a83',
}

offsets = [-2, -1, 0, 1, 2]
data_lists = [f_product, f_W, f_GHZ, f_Haar, f_Mixed]
sem_lists  = [sem_product, sem_W, sem_GHZ, sem_Haar, sem_Mixed]
labels = ['Product', 'W', 'GHZ', 'Haar', 'Mixed']

for off, data, sem, lab in zip(offsets, data_lists, sem_lists, labels):
    ax4.bar(xs + off * bw, data, bw, label=lab,
            color=cat_colors[lab], edgecolor='black', linewidth=0.3,
            yerr=sem, capsize=1.5, error_kw={'linewidth': 0.5, 'capthick': 0.5})

# Overall markers
ax4.scatter(xs, f_overall, marker='D', s=18, color='black', zorder=5,
            label='Overall')

ax4.set_xticks(xs)
ax4.set_xticklabels(strategies, rotation=30, ha='right', fontsize=8)
ax4.set_ylabel('Mean fidelity $\\bar{F}$')
ax4.set_ylim(0, 1.15)
ax4.legend(fontsize=7, ncol=6, loc='upper right', frameon=True,
           edgecolor='gray', borderpad=0.4, columnspacing=0.8,
           handlelength=1.4)
ax4.axhline(y=f_overall[0], color='#2166ac', linewidth=0.5,
            linestyle=':', alpha=0.6)
ax4.text(6.3, f_overall[0] + 0.0234,
         f'$K^*$ overall = {f_overall[0]:.3f}',
         ha='center', va='bottom', fontsize=7.2, color='#2166ac')

fig4.tight_layout()
fig4.savefig('fig4_sota_comparison.pdf')
fig4.savefig('fig4_sota_comparison.png')
print("Figure 4: fig4_sota_comparison.pdf")

# ── Figure 5: Dimensional dependence of Delta(S-AR) ────────────────

fig5, ax5 = plt.subplots(figsize=(COL_WIDTH, 2.2))

# Source: external_simulations_output.txt lines 137-140
n_vals = [2, 3, 4, 5]
delta_SAR = [0.573, 0.331, 0.100, 0.042]

ax5.plot(n_vals, delta_SAR, 'o-', color='#2166ac', linewidth=1.5,
         markersize=7, markeredgecolor='black', markeredgewidth=0.5,
         zorder=3, label=r'$\Delta(S{-}AR)$, MLE')

# Highlight n=4 sweet spot
ax5.axvline(x=4, color='#b2182b', linewidth=0.8, linestyle='--', alpha=0.6)
ax5.annotate('$2n = 2^{n-1}$\n(dim. rigidity)',
             xy=(4, delta_SAR[2]), xytext=(4.4, 0.35),
             fontsize=7, color='#b2182b',
             arrowprops=dict(arrowstyle='->', color='#b2182b',
                             lw=0.8),
             ha='left', va='center')

# Allocation fraction = D(AR-UR)/D(S-UR) from external_simulations_output.txt
# n=2: (0.414-0.479)/(0.987-0.479) = -13% (allocation hurts at n=2)
# n=3: (0.538-0.284)/(0.868-0.284) = 43%
# n=4: (0.733-0.302)/(0.833-0.302) = 81%  (hardware-confirmed: 76±12%)
# n=5: no S/AR/UR breakdown available
alloc_frac = [-13, 43, 81, None]
label_offsets = [(0.15, 0.04), (-0.15, 0.0), (-0.15, 0.0), (0.15, 0.0)]
label_ha =     ['left',         'right',       'right',       'left']
for ni, di, af, (ox, oy), ha in zip(n_vals, delta_SAR, alloc_frac, label_offsets, label_ha):
    if af is not None:
        ax5.annotate(f'${af:+d}$%', xy=(ni, di),
                     xytext=(ni + ox, di + oy),
                     fontsize=6.5, color='#555555', ha=ha, va='center')

ax5.set_xlabel('Number of qubits $n$')
ax5.set_ylabel(r'Residual $\Delta(S{-}AR)$')
ax5.set_xticks(n_vals)
ax5.set_xlim(1.5, 5.5)
ax5.set_ylim(-0.02, 0.65)

ax5.spines['top'].set_visible(False)
ax5.spines['right'].set_visible(False)
fig5.tight_layout(pad=1.2)
fig5.savefig('fig5_dimensional_dependence.pdf')
fig5.savefig('fig5_dimensional_dependence.png')
print("Figure 5: fig5_dimensional_dependence.pdf")

# ── Figure 6: Circuit efficiency ─────────────────────────────────────

fig6, (ax6a, ax6b) = plt.subplots(2, 1, figsize=(COL_WIDTH + 0.4, 5.6))

# Top panel: bases required at n=4, W state (IBM hardware)
# Source: w_repeat_results.json (K*, rand), oq_grouped_results_20260316.json (Rigetti)
methods = ['$K^*$', 'Random', 'Full\ntomography']
bases   = [29, 50, 81]
# W-state fidelities only; other states differ (see Table I)
fid_ibm = [0.872, 0.54, None]
fid_rig = [0.816, 0.568, None]
bar_colors = ['#2166ac', '#d6604d', '#999999']

bars6 = ax6a.bar(range(3), bases, color=bar_colors, edgecolor='black',
                  linewidth=0.5, width=0.6)

# Annotate fidelity on K* and Random bars (W state, IBM hardware)
for i in range(2):
    ax6a.text(i, bases[i] + 1.5,
              f'$F = {fid_ibm[i]:.2f}$\n($W_4$, IBM)',
              ha='center', va='bottom', fontsize=6, color='black')


ax6a.set_xticks(range(3))
ax6a.set_xticklabels(methods, fontsize=7)
ax6a.set_ylabel('Bases (circuits)')
ax6a.set_ylim(0, 100)
ax6a.set_title('$n = 4$ qubits', fontsize=8)

# Bottom panel: circuit scaling with processor size
# All lines count measurement BASES (= circuit executions), not operators.
# K*: 137 operators fit in 29 bases; Random: 50 bases; Full tomo: 3^n bases.
n_qubits_chain = np.array([4, 8, 12, 16, 20, 24, 28, 32])
kstar_circuits = np.empty_like(n_qubits_chain, dtype=float)
for idx, nq in enumerate(n_qubits_chain):
    n_patches = nq - 3
    n_colors = min(n_patches, 4)
    kstar_circuits[idx] = n_colors * 29  # 29 bases per color group

rand_circuits = np.array([min(nq - 3, 4) * 50
                           for nq in n_qubits_chain], dtype=float)
# Direct full tomography (impractical): 3^n
full_tomo = 3**n_qubits_chain

ax6b.semilogy(n_qubits_chain, kstar_circuits, 's-', color='#2166ac',
              linewidth=1.5, markersize=5, label='$K^*$ (compositional)',
              markeredgecolor='black', markeredgewidth=0.5)
ax6b.semilogy(n_qubits_chain, rand_circuits, 'o--', color='#d6604d',
              linewidth=1.2, markersize=5, label='Random (compositional)',
              markeredgecolor='black', markeredgewidth=0.5)
ax6b.semilogy(n_qubits_chain, full_tomo, '^:', color='#999999',
              linewidth=1.0, markersize=5, label='Full tomography',
              markeredgecolor='black', markeredgewidth=0.5)

ax6b.set_xlabel('Processor qubits')
ax6b.set_ylabel('Total circuits')
ax6b.set_xticks(n_qubits_chain)
ax6b.set_xticklabels(n_qubits_chain, fontsize=6)
ax6b.legend(fontsize=6, loc='upper left')
# (sub-title moved to caption)
# Mark the O(1) region
ax6b.axhline(y=116, color='#2166ac', linewidth=0.4, linestyle=':', alpha=0.5)
ax6b.text(18, 250, '$K^*$: 116', fontsize=7, color='#2166ac', ha='center',
          bbox=dict(facecolor='white', edgecolor='none', pad=1.5, alpha=0.8))

fig6.tight_layout(h_pad=2.0)
fig6.savefig('fig6_circuit_efficiency.pdf')
fig6.savefig('fig6_circuit_efficiency.png')
print("Figure 6: fig6_circuit_efficiency.pdf")

# ── Figure 7: GHZ eigenvalue-mass mismatch and subspace recovery ────

fig7, (ax7a, ax7b) = plt.subplots(2, 1, figsize=(COL_WIDTH + 0.4, 4.8))

# Left panel: information distribution vs K* allocation by weight
weights_ghz = [0, 1, 2, 3, 4]
# GHZ nonzero expectations by weight: w0=1(IIII), w1=0, w2=6(ZZ), w3=0, w4=9(phase+ZZZZ)
# But "information" = number of nonzero expectations (excluding identity for info fraction)
ghz_info_pct = [0, 0, 40, 0, 60]  # 6/15=40%, 9/15=60% of nontrivial info

# K* eigenvalue mass allocation (from fig2 data)
kstar_alloc_pct = [6.6, 40.9, 17.5, 23.4, 11.7]

x7 = np.arange(len(weights_ghz))
bw7 = 0.35

bars_info = ax7a.bar(x7 - bw7/2, ghz_info_pct, bw7,
                      label='GHZ information', color='#b2182b',
                      edgecolor='black', linewidth=0.5)
bars_alloc = ax7a.bar(x7 + bw7/2, kstar_alloc_pct, bw7,
                       label='$K^*$ eigenvalue mass', color='#2166ac',
                       edgecolor='black', linewidth=0.5)


ax7a.set_xlabel('Pauli weight $w$')
ax7a.set_ylabel('Fraction (\\%)')
ax7a.set_xticks(x7)
ax7a.set_ylim(0, 72)
ax7a.legend(fontsize=6.5, loc='upper left')
# (sub-title moved to caption)

# Bottom panel: fidelity recovery from same K* hardware data
# Source: data/ghz_dfe_results.json (derived from 214922.json via ghz_dfe_reanalysis.py)
# 0.496: full-state MLE of GHZ hardware run (f_mle_full)
# 0.926: subspace projection into {|0000>,|1111>} (f_subspace_mle)
# 0.618: direct fidelity estimation lower bound (f_dfe_informative)
estimators = ['Full-state\nMLE', 'Subspace\nestimation', 'DFE\n(lower bound)', 'Auto-switch\n(fingerprint)']
fidelities = [0.496, 0.926, 0.618, 0.887]
bar_colors_r = ['#d6604d', '#2166ac', '#92c5de', '#e08214']

bars_r = ax7b.bar(range(4), fidelities, color=bar_colors_r,
                   edgecolor='black', linewidth=0.5, width=0.6)

# Annotate values
for i, f in enumerate(fidelities):
    ax7b.text(i, f + 0.02, f'$F = {f:.3f}$',
              ha='center', va='bottom', fontsize=7)

ax7b.set_ylabel('GHZ fidelity $F$')
ax7b.set_xticks(range(4))
ax7b.set_xticklabels(estimators, fontsize=7)
ax7b.set_ylim(0, 1.12)
ax7b.axhline(y=0.5, color='gray', linewidth=0.5, linestyle='--', alpha=0.5)
# (sub-title moved to caption)


fig7.tight_layout(h_pad=2.0)
fig7.savefig('fig7_ghz_mismatch.pdf')
fig7.savefig('fig7_ghz_mismatch.png')
print("Figure 7: fig7_ghz_mismatch.pdf")

print("\nAll figures generated.")
