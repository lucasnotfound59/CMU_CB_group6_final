#!/usr/bin/env python3
"""Generate a pipeline diagram for the B-factor / docking experiment."""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

fig, ax = plt.subplots(figsize=(14, 9))
ax.set_xlim(0, 14)
ax.set_ylim(0, 9)
ax.axis('off')

# Colors
C_TITLE = '#1a1a2e'
C_STEP = '#16213e'
C_BOX = '#e8f0fe'
C_BORDER = '#3a7ca5'
C_ARROW = '#d35400'
C_KEY = '#27ae60'
C_RESULT = '#8e44ad'
C_BG = '#fafafa'

fig.patch.set_facecolor(C_BG)
ax.set_facecolor(C_BG)

# Title
ax.text(7, 8.5, 'Experimental Pipeline', fontsize=18, fontweight='bold',
        ha='center', va='center', color=C_TITLE,
        fontfamily='sans-serif')
ax.text(7, 8.1, 'Can B-factors identify flexible binding-site residues\nthat reduce rigid docking accuracy?',
        fontsize=10, ha='center', va='center', color='#555', style='italic')

# --- Helper functions ---
def draw_box(x, y, w, h, text, subtitle=None, color=C_BOX, border=C_BORDER, fontsize=9):
    box = FancyBboxPatch((x - w/2, y - h/2), w, h,
                         boxstyle="round,pad=0.1", facecolor=color,
                         edgecolor=border, linewidth=1.5)
    ax.add_patch(box)
    if subtitle:
        ax.text(x, y + 0.15, text, fontsize=fontsize, fontweight='bold',
                ha='center', va='center', color=C_STEP)
        ax.text(x, y - 0.2, subtitle, fontsize=7.5, ha='center', va='center', color='#555')
    else:
        ax.text(x, y, text, fontsize=fontsize, fontweight='bold',
                ha='center', va='center', color=C_STEP)

def draw_arrow(x1, y1, x2, y2, text='', color=C_ARROW):
    ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle='->', color=color, lw=2))
    if text:
        mx, my = (x1+x2)/2, (y1+y2)/2
        ax.text(mx + 0.15, my, text, fontsize=7, color=color, style='italic')

# === ROW 1: Input structures ===
draw_box(2.5, 7.0, 2.8, 0.7, 'PDB Structures', 'Same target (e.g. CDK2)', fontsize=9)
draw_box(6.5, 7.0, 2.8, 0.7, 'Holo Ligands', 'Co-crystallized drugs', fontsize=9)
draw_box(10.5, 7.0, 2.8, 0.7, 'B-factor Extraction', 'Per-residue thermal\ndisplacement', fontsize=8.5)

# === ROW 2: Cross-docking ===
draw_box(4.5, 5.5, 3.5, 0.8, 'Cross-Docking', 'Receptor_i  ×  Ligand_j   (i ≠ j)', fontsize=10,
         color='#fff3e0', border='#e65100')
ax.text(4.5, 5.1, 'AutoDock Vina', fontsize=7.5, ha='center', color='#e65100', style='italic')

# === ROW 3: Two outputs ===
draw_box(2.5, 3.8, 2.8, 0.7, 'Docking RMSD', 'Predicted vs crystal pose\nRMSD > 2Å → Failure', fontsize=8.5,
         color='#fce4ec', border='#c62828')
draw_box(6.5, 3.8, 2.8, 0.7, 'B-factor Profile', 'Mean / Max / Relative\nper binding-site residue', fontsize=8.5,
         color='#e8f5e9', border='#2e7d32')

# === ROW 4: Correlation ===
draw_box(4.5, 2.3, 3.8, 0.8, 'Correlation Analysis', 'Spearman ρ:  B-factor  vs  RMSD', fontsize=10,
         color='#f3e5f5', border='#6a1b9a')

# === ROW 5: Conclusion ===
draw_box(4.5, 0.8, 4.5, 0.8,
         'Conclusion',
         fontsize=10, color='#e0f2f1', border='#00695c')

# Arrows
# Row 1 → Row 2
draw_arrow(2.5, 6.6, 4.0, 5.9, 'receptor')
draw_arrow(6.5, 6.6, 5.0, 5.9, 'ligand')
draw_arrow(10.5, 6.6, 10.5, 5.9)

# Row 2 → Row 3
draw_arrow(3.8, 5.1, 2.8, 4.2, 'RMSD')
draw_arrow(5.2, 5.1, 6.2, 4.2, 'B-factors')

# Row 3 → Row 4
draw_arrow(2.5, 3.4, 4.0, 2.7)
draw_arrow(6.5, 3.4, 5.0, 2.7)

# Row 4 → Row 5
draw_arrow(4.5, 1.9, 4.5, 1.2)

# Key insight box
ax.text(10.5, 3.8, 'Key Logic', fontsize=9, fontweight='bold', ha='center', color=C_KEY)
ax.text(10.5, 3.3, 'Self-docking:  RMSD ≈ 0  ✓', fontsize=8, ha='center', color='#555',
        family='monospace')
ax.text(10.5, 2.9, 'Cross-docking:  RMSD = ?  →', fontsize=8, ha='center', color='#555',
        family='monospace')
ax.text(10.5, 2.5, 'If RMSD ↑ when B-factor ↑', fontsize=8, ha='center', color=C_KEY,
        fontweight='bold')
ax.text(10.5, 2.1, '→ B-factors predict\n   docking failure', fontsize=8, ha='center', color=C_KEY)

# Conclusion text
ax.text(4.5, 0.8, 'B-factors partially predict which flexible residues\ncause rigid docking to fail — but backbone–sidechain\ndecoupling (Najmanovich 2000) limits accuracy',
        fontsize=7.5, ha='center', va='center', color='#00695c')

plt.tight_layout()
plt.savefig('/Users/xinlu/Desktop/CMU_CB/experiment/pipeline_diagram.png', dpi=200,
            bbox_inches='tight', facecolor=C_BG)
print('Saved: experiment/pipeline_diagram.png')
