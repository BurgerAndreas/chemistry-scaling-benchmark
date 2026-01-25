#!/usr/bin/env python3
"""
Plot benchmarking results showing time and memory scaling across methods.
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib.ticker import LogLocator

# Set seaborn poster theme
sns.set_theme(style="whitegrid", context="poster")

# Read data
results_file = "results/results_600.csv"
df = pd.read_csv(results_file)

# Exclude rows with a non-empty error_message
df = df[df['error_message'].isna() | (df['error_message'] == '')]

# Filter only successful calculations
df_success = df[df['success'] == True].copy()

# Exclude all rows where the time is the max time (600)
df_success = df_success[df_success['time_seconds'] < 600]

print(f"Total calculations: {len(df)}")
print(f"Successful: {len(df_success)}")

# Create method labels
df_success['method_label'] = df_success.apply(
    lambda row: f"{row['method_name']}" if pd.isna(row['basis_set'])
    else f"{row['method_name']}/{row['basis_set']}",
    axis=1
)

# Average over repetitions and different molecules for the same natoms
df_plot_data = df_success.groupby(['method_label', 'natoms']).agg(
    time_seconds=('time_seconds', 'mean'),
    peak_memory_mb=('peak_memory_mb', 'mean'),
    method_name=('method_name', 'first'),
    basis_set=('basis_set', 'first'),
    method_category=('method_category', 'first')
).reset_index()

# Define method order by computational cost (ladder to heaven)
# basis_set = '6-31G(d)'
for basis_set in ['STO-3G', '6-31G(d)']:
    method_order = [
        'UFF',
        'GFN2-xTB',
        # f'RHF/{basis_set}',
        # f'SVWN/{basis_set}',      # LDA
        # f'PBE/{basis_set}',       # GGA
        # f'TPSS/{basis_set}',      # Meta-GGA
        # f'B3LYP/{basis_set}',     # Hybrid
        f'wB97X/{basis_set}',     # Range-separated hybrid
        f'MP2/{basis_set}',
        f'CCSD/{basis_set}',
        # 'CCSD(T)/{basis_set}',
        f'FCI/{basis_set}',
    ]

    # Filter to only include STO-3G or None (NaN) basis sets
    df_combined = df_plot_data[df_plot_data['basis_set'].isna() | (df_plot_data['basis_set'] == basis_set)]

    # Create combined plot
    fig, axes = plt.subplots(2, 1, figsize=(16, 20))

    # Get methods present in data, ordered by cost
    available_methods = df_combined['method_label'].unique()
    methods_combined = [m for m in method_order if m in available_methods]
    # Add any methods not in the predefined order (shouldn't happen but safe)
    # methods_combined += [m for m in available_methods if m not in method_order]
    palette = sns.color_palette("husl", n_colors=len(methods_combined))

    # Plot 1: Time scaling
    ax1 = axes[0]
    for i, method in enumerate(methods_combined):
        data = df_combined[df_combined['method_label'] == method].sort_values('natoms')
        ax1.scatter(data['natoms'], data['time_seconds'],
                label=method, s=150, alpha=0.7, color=palette[i])
        if len(data) > 1:
            ax1.plot(data['natoms'], data['time_seconds'],
                    alpha=0.4, linewidth=2, color=palette[i])

    ax1.set_xlabel('Number of Atoms')
    ax1.set_ylabel('Time (seconds)')
    ax1.set_title('Computational Time Scaling (Averaged)', fontweight='bold', pad=20)
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=14)
    ax1.grid(True, alpha=0.3, which='both')
    ax1.xaxis.set_minor_locator(LogLocator(base=10, subs=np.arange(2, 10)))
    ax1.yaxis.set_minor_locator(LogLocator(base=10, subs=np.arange(2, 10)))
    ax1.tick_params(axis='both', which='minor', length=4)
    ax1.tick_params(axis='both', which='major', length=8)

    # Plot 2: Memory scaling
    ax2 = axes[1]
    for i, method in enumerate(methods_combined):
        data = df_combined[df_combined['method_label'] == method].sort_values('natoms')
        ax2.scatter(data['natoms'], data['peak_memory_mb'],
                label=method, s=150, alpha=0.7, color=palette[i])
        if len(data) > 1:
            ax2.plot(data['natoms'], data['peak_memory_mb'],
                    alpha=0.4, linewidth=2, color=palette[i])

    ax2.set_xlabel('Number of Atoms')
    ax2.set_ylabel('Peak Memory (MB)')
    ax2.set_title('Memory Usage Scaling (Averaged)', fontweight='bold', pad=20)
    ax2.set_xscale('log')
    ax2.set_yscale('log')
    ax2.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=14)
    ax2.grid(True, alpha=0.3, which='both')
    ax2.xaxis.set_minor_locator(LogLocator(base=10, subs=np.arange(2, 10)))
    ax2.yaxis.set_minor_locator(LogLocator(base=10, subs=np.arange(2, 10)))
    ax2.tick_params(axis='both', which='minor', length=4)
    ax2.tick_params(axis='both', which='major', length=8)

    plt.tight_layout()
    plt.savefig(f'results/scaling_{basis_set}.png', dpi=150, bbox_inches='tight')
    print(f"\nSaved: results/scaling_{basis_set}.png")


############################################################################################
# Plot by method category
for category in df_plot_data['method_category'].unique():
    df_cat = df_plot_data[df_plot_data['method_category'] == category]
    methods_cat = sorted(df_cat['method_label'].unique())
    palette_cat = sns.color_palette("husl", n_colors=len(methods_cat))

    fig, axes = plt.subplots(1, 2, figsize=(20, 8))

    # Time plot
    for i, method in enumerate(methods_cat):
        data = df_cat[df_cat['method_label'] == method].sort_values('natoms')
        axes[0].scatter(data['natoms'], data['time_seconds'],
                       label=method, s=150, alpha=0.7, color=palette_cat[i])
        if len(data) > 1:
            axes[0].plot(data['natoms'], data['time_seconds'],
                        alpha=0.4, linewidth=2, color=palette_cat[i])

    axes[0].set_xlabel('Number of Atoms')
    axes[0].set_ylabel('Time (seconds)')
    axes[0].set_title(f'{category.replace("_", " ").title()} - Time (Averaged)', fontweight='bold')
    axes[0].set_xscale('log')
    axes[0].set_yscale('log')
    axes[0].legend()
    axes[0].grid(True, alpha=0.3, which='both')
    axes[0].xaxis.set_minor_locator(LogLocator(base=10, subs=np.arange(2, 10)))
    axes[0].yaxis.set_minor_locator(LogLocator(base=10, subs=np.arange(2, 10)))
    axes[0].tick_params(axis='both', which='minor', length=4)
    axes[0].tick_params(axis='both', which='major', length=8)

    # Memory plot
    for i, method in enumerate(methods_cat):
        data = df_cat[df_cat['method_label'] == method].sort_values('natoms')
        axes[1].scatter(data['natoms'], data['peak_memory_mb'],
                       label=method, s=150, alpha=0.7, color=palette_cat[i])
        if len(data) > 1:
            axes[1].plot(data['natoms'], data['peak_memory_mb'],
                        alpha=0.4, linewidth=2, color=palette_cat[i])

    axes[1].set_xlabel('Number of Atoms')
    axes[1].set_ylabel('Peak Memory (MB)')
    axes[1].set_title(f'{category.replace("_", " ").title()} - Memory (Averaged)', fontweight='bold')
    axes[1].set_xscale('log')
    axes[1].set_yscale('log')
    axes[1].legend()
    axes[1].grid(True, alpha=0.3, which='both')
    axes[1].xaxis.set_minor_locator(LogLocator(base=10, subs=np.arange(2, 10)))
    axes[1].yaxis.set_minor_locator(LogLocator(base=10, subs=np.arange(2, 10)))
    axes[1].tick_params(axis='both', which='minor', length=4)
    axes[1].tick_params(axis='both', which='major', length=8)

    plt.tight_layout()
    plt.savefig(f'results/scaling_{category}.png', dpi=150, bbox_inches='tight')
    print(f"Saved: results/scaling_{category}.png")

print("\n" + "="*60)
print("SUMMARY (Combined Plot - Ordered by Cost)")
print("="*60)
for method in methods_combined:
    data = df_combined[df_combined['method_label'] == method]
    if not data.empty:
        print(f"\n{method}:")
        print(f"  Averaged calculations: {len(data)}")
        print(f"  Atoms: {data['natoms'].min()}-{data['natoms'].max()}")
        print(f"  Time: {data['time_seconds'].min():.2f}-{data['time_seconds'].max():.2f}s")
        print(f"  Memory: {data['peak_memory_mb'].min():.0f}-{data['peak_memory_mb'].max():.0f} MB")
