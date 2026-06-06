"""
Example usage of the pymplot function API.
Requires:  pip install -e .  (from the repo root)
Run from anywhere:  python3 example/test_api.py
Output files are written to /tmp/pymplot_test_*.pdf
"""

import numpy as np
import matplotlib.pyplot as plt
from pymplot import (showmatrix, showcontour, showslice, showslicon,
                     showwiggle, showcolorbar, showgraph, showvolume)

out = '/tmp/pymplot_test_{}.pdf'

# ── synthetic data ────────────────────────────────────────────────────────────

np.random.seed(42)

x = np.linspace(0, 2 * np.pi, 120)
y = np.linspace(0, np.pi, 80)
XX, YY = np.meshgrid(x, y)
data2d = np.sin(XX) * np.cos(YY) + 0.1 * np.random.randn(80, 120)
data2d = data2d.astype('float32')

i1, i2, i3 = np.mgrid[0:40, 0:50, 0:30]
data3d = np.sin(i1 / 6.0) * np.cos(i2 / 8.0) * np.exp(-i3 / 20.0)
data3d = data3d.astype('float32')

t = np.linspace(0, 4 * np.pi, 200)

# ── subplot stress test ───────────────────────────────────────────────────────

print('\n--- 2x2 subplot stress test ---')
fig, axes = plt.subplots(2, 2, figsize=(12, 9))

showmatrix(data2d, colormap='seismic', legend=True, unit='P-wave velocity (m/s)',   ax=axes[0, 0])
showcontour(data2d, overlay=True,                     ax=axes[0, 1])
showwiggle(data2d[:, :30], fill=1,                    ax=axes[1, 0])
showvolume(data3d,                                    ax=axes[1, 1])

for ax, label in zip(axes.flat, 'abcd'):
    ax.text(0.02, 1.2, f'({label})',
            transform=ax.transAxes,
            fontsize=14, fontweight='bold',
            va='top', ha='center')

plt.tight_layout()
plt.savefig(out.format('subplot_2x2'), dpi=150, bbox_inches='tight')
print(f'  → {out.format("subplot_2x2")}')
plt.show()
# plt.savefig(out.format('subplot_2x2'), dpi=150, bbox_inches='tight')

exit(0)

# ── individual function tests ─────────────────────────────────────────────────

print('\n--- showmatrix ---')
showmatrix(data2d)
showmatrix(data2d, colormap='jet', legend=True, outfile=out.format('matrix'))
showmatrix(data2d, background=np.abs(data2d), outfile=out.format('matrix_bg'))

print('\n--- showcontour ---')
showcontour(data2d)
showcontour(data2d, outfile=out.format('contour'))
showcontour(data2d, overlay=True, outfile=out.format('contour_overlay'))

print('\n--- showslice ---')
showslice(data3d)
showslice(data3d, slice1=20, slice2=25, slice3=15,
          outfile=out.format('slice'))
showslice(data3d, background=np.abs(data3d), outfile=out.format('slice_bg'))

print('\n--- showslicon ---')
showslicon(data3d)
showslicon(data3d, outfile=out.format('slicon'))

print('\n--- showwiggle ---')
showwiggle(data2d[:, :40])
showwiggle(data2d[:, :40], fill=1, outfile=out.format('wiggle'))

print('\n--- showcolorbar ---')
showcolorbar(colormap='jet', cmin=-1, cmax=1, unit='Amplitude',
             outfile=out.format('colorbar_right'))
showcolorbar(lloc='bottom', colormap='gray',
             outfile=out.format('colorbar_bottom'))

print('\n--- showgraph ---')
showgraph(np.sin(t))
showgraph(np.column_stack([t, np.sin(t)]), outfile=out.format('graph'))

print('\nDone. Output PDFs written to /tmp/pymplot_test_*.pdf')
