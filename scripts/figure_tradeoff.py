import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import numpy as np
import seaborn as sns


BOX_LABEL_FONT_SIZE = 20
BOX_CONTENT_FONT_SIZE = 18
OUTSIDE_BOX_FONT_SIZE = 20

# Set seaborn poster theme
sns.set_theme(style="white", context="poster")

# Create figure and axis
fig, ax = plt.subplots(figsize=(10, 10))

# Define positions along the diagonal (x, y coordinates)
# Lower left = low cost, low accuracy; Upper right = high cost, high accuracy
positions = {
    "Force Fields": (0.8, 0.8),
    "Tight Binding": (1.3, 1.3),
    "Hartree-Fock": (1.8, 1.8),
    # DFT methods (will be in a box)
    "DFT_box": (3.2, 3.2),
    "LDA": (2.5, 2.5),
    "GGA": (2.9, 2.9),
    "Meta-GGA": (3.3, 3.3),
    "Hybrid": (3.7, 3.7),
    "Range-Sep Hybrid": (4.1, 4.1),
    "Double Hybrid": (4.5, 4.5),
    # ab initio methods (will be in a box)
    "ab_initio_box": (6.2, 6.2),
    "MP2": (5.5, 5.5),
    "CCSD": (5.9, 5.9),
    "CCSD(T)": (6.3, 6.3),
    "CCSD(Q)": (6.7, 6.7),
    "FCI": (7.1, 7.1),
}

# Draw diagonal arrow showing the trend (draw first so it's behind boxes)
TOP = 7.7
arrow_start = (0.1, 0.1)
arrow_end = (TOP, TOP)
ax.annotate(
    "",
    xy=arrow_end,
    xytext=arrow_start,
    arrowprops=dict(arrowstyle="->", lw=2, color="#B8B8B8", alpha=0.5),
    zorder=0,
)

# Plot main methods (outside boxes)
main_methods = ["Force Fields", "Tight Binding", "Hartree-Fock"]
for method in main_methods:
    x, y = positions[method]
    x += 0.1
    y += 0.1
    ax.plot(
        x,
        y,
        "s",
        markersize=8,
        color="#2E86AB",
        zorder=3,
    )
    ax.text(
        x,
        y,
        method,
        # ha="right",
        # va="bottom",
        ha="center",
        va="center",
        fontsize=OUTSIDE_BOX_FONT_SIZE,
        bbox=dict(
            boxstyle="round,pad=0.2",
            facecolor="#D4E5EE",
            edgecolor="none",
            # facecolor="white",
            # edgecolor="#2E86AB",
            linewidth=1.5,
            alpha=1.0,
        ),
        color="#2E86AB",
        zorder=4,
    )

# DFT dotted box
dft_methods = ["LDA", "GGA", "Meta-GGA", "Hybrid", "Range-Sep Hybrid", "Double Hybrid"]
dft_x_coords = [positions[m][0] for m in dft_methods]
dft_y_coords = [positions[m][1] for m in dft_methods]
# dft_y_coords = [y - 0.01 for y in dft_y_coords]

dft_box = FancyBboxPatch(
    xy=(min(dft_x_coords) - 0.4, min(dft_y_coords) - 0.1),
    width=max(dft_x_coords) - min(dft_x_coords) + 0.8,
    height=max(dft_y_coords) - min(dft_y_coords) + 0.4,
    boxstyle="round,pad=0.1",
    # linestyle='--',
    linestyle="",
    linewidth=2,
    edgecolor="#A23B72",
    # facecolor="#F18F01",
    facecolor="#A23B72",
    alpha=0.15,
    zorder=1,
)
ax.add_patch(dft_box)

# Add DFT label
ax.text(
    positions["DFT_box"][0] + 0.3,
    max(dft_y_coords) + 0.45,
    "Density Functional \nTheory (DFT)",
    ha="center",
    va="bottom",
    fontweight="bold",
    color="#A23B72",
    zorder=4,
    fontsize=BOX_LABEL_FONT_SIZE,
)

# Plot DFT methods
for method in dft_methods:
    x, y = positions[method]
    ax.plot(x, y, "s", markersize=8, color="#A23B72", zorder=3)
    ax.text(
        x + 0.09,
        y + 0.05,
        method,
        ha="right",
        va="bottom",
        color="#A23B72",
        zorder=4,
        fontsize=BOX_CONTENT_FONT_SIZE,
    )

SHIFT_AB_INITIO = 0.3
# ab initio dotted box
ab_initio_methods = ["MP2", "CCSD", "CCSD(T)", "CCSD(Q)", "FCI"]
ab_x_coords = [positions[m][0] - SHIFT_AB_INITIO for m in ab_initio_methods]
ab_y_coords = [positions[m][1] - SHIFT_AB_INITIO for m in ab_initio_methods]

ab_box = FancyBboxPatch(
    (min(ab_x_coords) - 0.35, min(ab_y_coords) - 0.1),
    max(ab_x_coords) - min(ab_x_coords) + 0.6,
    max(ab_y_coords) - min(ab_y_coords) + 0.4,
    boxstyle="round,pad=0.1",
    # linestyle='--',
    linestyle="",
    linewidth=2,
    edgecolor="#C73E1D",
    facecolor="#C73E1D",
    alpha=0.15,
    zorder=1,
)
ax.add_patch(ab_box)

# Add ab initio label
ax.text(
    positions["ab_initio_box"][0] + 0.05 - SHIFT_AB_INITIO,
    max(ab_y_coords) + 0.45,
    "Ab Initio \n(Wavefunction)",
    ha="center",
    va="bottom",
    fontweight="bold",
    color="#C73E1D",
    zorder=4,
    fontsize=BOX_LABEL_FONT_SIZE,
)

# Plot ab initio methods
for method in ab_initio_methods:
    x, y = positions[method]
    x -= SHIFT_AB_INITIO
    y -= SHIFT_AB_INITIO
    ax.plot(x, y, "s", markersize=8, color="#C73E1D", zorder=3)
    ax.text(
        x + 0.09,
        y + 0.05,
        method,
        ha="right",
        va="bottom",
        color="#C73E1D",
        zorder=4,
        fontsize=BOX_CONTENT_FONT_SIZE,
    )

# Set labels and title
ax.set_ylabel("Compute →")
ax.set_xlabel("Accuracy →")
# ax.set_title('Accuracy vs Computational Cost in Computational Chemistry Methods', pad=20)

# Set limits and remove ticks
ax.set_xlim(0, TOP + 0.1)
ax.set_ylim(0, TOP + 0.1)
ax.set_xticks([])
ax.set_yticks([])

# Add grid for better readability
# ax.grid(True, alpha=0.2, linestyle=':', linewidth=0.5)

# Remove top and right spines
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

# Add arrow tips to the axes using seaborn's despine with offset
sns.despine(ax=ax, left=False, bottom=False, offset=0, trim=False)
ax.spines["left"].set_position(("data", 0))
ax.spines["bottom"].set_position(("data", 0))
ax.spines["left"].set_visible(True)
ax.spines["bottom"].set_visible(True)

# Set aspect ratio to be equal
ax.set_aspect("equal")

# Adjust layout
plt.tight_layout(pad=0.0)

# Save figure
fname = "results/accuracy_cost_tradeoff.png"
plt.savefig(fname, dpi=300, bbox_inches="tight")
print(f"Figure saved as '{fname}'")

# Show the plot
plt.show()
