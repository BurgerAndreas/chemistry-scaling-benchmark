import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.patches import FancyBboxPatch
import numpy as np

# Create figure and axis
fig, ax = plt.subplots(figsize=(14, 10))

# Define positions along the diagonal (x, y coordinates)
# Lower left = low cost, low accuracy; Upper right = high cost, high accuracy
positions = {
    'Classical Force Fields': (0.5, 0.5),
    'Tight Binding': (1.5, 1.5),
    'Hartree-Fock': (2.8, 2.8),
    # DFT methods (will be in a box)
    'DFT_box': (4.5, 4.5),
    'LDA': (3.8, 3.8),
    'GGA': (4.2, 4.2),
    'Meta-GGA': (4.6, 4.6),
    'Hybrid': (5.0, 5.0),
    'range-sep Hybrid': (5.4, 5.4),
    'Double Hybrid': (5.8, 5.8),
    # ab initio methods (will be in a box)
    'ab_initio_box': (7.5, 7.5),
    'MP2': (6.8, 6.8),
    'CCSD': (7.2, 7.2),
    'CCSD(T)': (7.6, 7.6),
    'CCSD(Q)': (8.0, 8.0),
    'FCI': (8.4, 8.4),
}

# Plot main methods (outside boxes)
main_methods = ['Classical Force Fields', 'Tight Binding', 'Hartree-Fock']
for method in main_methods:
    x, y = positions[method]
    ax.plot(x, y, 'o', markersize=12, color='#2E86AB', zorder=3)
    ax.text(x, y, method, fontsize=11, ha='center', va='center',
            bbox=dict(boxstyle='round,pad=0.5', facecolor='white', edgecolor='#2E86AB', linewidth=1.5),
            zorder=4)

# DFT dotted box
dft_methods = ['LDA', 'GGA', 'Meta-GGA', 'Hybrid', 'range-sep Hybrid', 'Double Hybrid']
dft_x_coords = [positions[m][0] for m in dft_methods]
dft_y_coords = [positions[m][1] for m in dft_methods]

dft_box = FancyBboxPatch(
    (min(dft_x_coords) - 0.4, min(dft_y_coords) - 0.3),
    max(dft_x_coords) - min(dft_x_coords) + 0.8,
    max(dft_y_coords) - min(dft_y_coords) + 0.6,
    boxstyle="round,pad=0.1",
    linestyle='--',
    linewidth=2,
    edgecolor='#A23B72',
    facecolor='#F18F01',
    alpha=0.15,
    zorder=1
)
ax.add_patch(dft_box)

# Add DFT label
ax.text(positions['DFT_box'][0], max(dft_y_coords) + 0.5,
        'Density Functional Theory (DFT)',
        fontsize=12, ha='center', va='bottom', fontweight='bold',
        color='#A23B72', zorder=4)

# Plot DFT methods
for method in dft_methods:
    x, y = positions[method]
    ax.plot(x, y, 's', markersize=8, color='#A23B72', zorder=3)
    ax.text(x, y + 0.15, method, fontsize=9, ha='center', va='bottom',
            color='#A23B72', zorder=4)

# ab initio dotted box
ab_initio_methods = ['MP2', 'CCSD', 'CCSD(T)', 'CCSD(Q)', 'FCI']
ab_x_coords = [positions[m][0] for m in ab_initio_methods]
ab_y_coords = [positions[m][1] for m in ab_initio_methods]

ab_box = FancyBboxPatch(
    (min(ab_x_coords) - 0.4, min(ab_y_coords) - 0.3),
    max(ab_x_coords) - min(ab_x_coords) + 0.8,
    max(ab_y_coords) - min(ab_y_coords) + 0.6,
    boxstyle="round,pad=0.1",
    linestyle='--',
    linewidth=2,
    edgecolor='#C73E1D',
    facecolor='#6A994E',
    alpha=0.15,
    zorder=1
)
ax.add_patch(ab_box)

# Add ab initio label
ax.text(positions['ab_initio_box'][0], max(ab_y_coords) + 0.5,
        'ab initio (Wavefunction Methods)',
        fontsize=12, ha='center', va='bottom', fontweight='bold',
        color='#C73E1D', zorder=4)

# Plot ab initio methods
for method in ab_initio_methods:
    x, y = positions[method]
    ax.plot(x, y, '^', markersize=8, color='#C73E1D', zorder=3)
    ax.text(x, y + 0.15, method, fontsize=9, ha='center', va='bottom',
            color='#C73E1D', zorder=4)

# Draw diagonal arrow showing the trend
arrow_start = (0.2, 0.2)
arrow_end = (8.7, 8.7)
ax.annotate('', xy=arrow_end, xytext=arrow_start,
            arrowprops=dict(arrowstyle='->', lw=2, color='gray', alpha=0.4))

# Set labels and title
ax.set_xlabel('Computational Cost →', fontsize=14, fontweight='bold')
ax.set_ylabel('Accuracy →', fontsize=14, fontweight='bold')
ax.set_title('Accuracy vs Computational Cost in Computational Chemistry Methods',
             fontsize=16, fontweight='bold', pad=20)

# Set limits and remove ticks
ax.set_xlim(0, 9)
ax.set_ylim(0, 9)
ax.set_xticks([])
ax.set_yticks([])

# Add grid for better readability
ax.grid(True, alpha=0.2, linestyle=':', linewidth=0.5)

# Set aspect ratio to be equal
ax.set_aspect('equal')

# Adjust layout
plt.tight_layout()

# Save figure
plt.savefig('chemistry_accuracy_cost_tradeoff.png', dpi=300, bbox_inches='tight')
print("Figure saved as 'chemistry_accuracy_cost_tradeoff.png'")

# Show the plot
plt.show()
