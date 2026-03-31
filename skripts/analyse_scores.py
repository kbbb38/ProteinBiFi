import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
import os

# Read floats from file
def read_floats(filename):
    with open(filename, 'r') as f:
        data = [float(line.strip()) for line in f if line.strip()]
    return np.array(data)

# ================= CONFIGURATION =================
# Add second file to overlay (leave as None for single plot)
input_file_1 = '/storage/mi/malek01/ForschungsPraktikumRKI/data/human/tanimoto_scores/gt_FDR1_bs1.txt'
input_file_2 = '/storage/mi/malek01/ForschungsPraktikumRKI/data/human/tanimoto_scores/all_bs1.txt'  # ← Set to 'your_file_2.txt' to overlay two densities
output_file = '/storage/mi/malek01/ForschungsPraktikumRKI/data/human/tanimoto_scores/img/FDR1_bs1_allvsgt.png'
# =================================================

# Read data
data1 = read_floats(input_file_1)
data2 = read_floats(input_file_2) if input_file_2 else None

# Server-friendly backend
plt.switch_backend('Agg')

# Create figure
fig, ax = plt.subplots(figsize=(10, 6), dpi=150)

# ggplot2-style background
ax.set_facecolor('#EBEBEB')
ax.grid(True, color='white', linestyle='-', linewidth=0.8, alpha=0.9)
ax.set_axisbelow(True)

# Color palette (ggplot2 inspired)
colors = ['#0072B2', '#D55E00']  # Blue and Orange for contrast
labels = ['Ground Truth scores', 'Best scores per library spectra']

# Plot first density
kde1 = gaussian_kde(data1)
x_range1 = np.linspace(max(0, data1.min()), data1.max(), 1000)
density1 = kde1(x_range1)
ax.fill_between(x_range1, density1, color=colors[0], alpha=0.35, label=labels[0])
ax.plot(x_range1, density1, color=colors[0], linewidth=2.2)

# Plot second density (if provided)
if data2 is not None:
    kde2 = gaussian_kde(data2)
    x_range2 = np.linspace(max(0, data2.min()), data2.max(), 1000)
    density2 = kde2(x_range2)
    ax.fill_between(x_range2, density2, color=colors[1], alpha=0.35, label=labels[1])
    ax.plot(x_range2, density2, color=colors[1], linewidth=2.2)

# Labels
ax.set_xlabel('Tanimoto Score', fontsize=13, fontweight='normal', color='#333333')
ax.set_ylabel('Density', fontsize=13, fontweight='normal', color='#333333')
ax.set_title('Tanimoto Score Distribution', fontsize=15, fontweight='bold', pad=20, color='#222222')

# Style spines
for spine in ax.spines.values():
    spine.set_edgecolor('#AAAAAA')
    spine.set_linewidth(0.8)

# Ticks
ax.tick_params(axis='both', colors='#555555', labelsize=11)

# Force x-axis to start at 0
ax.set_xlim(left=0)

# Legend
legend = ax.legend(fontsize=11, frameon=True, facecolor='white', 
                   edgecolor='#CCCCCC', shadow=False)
legend.get_frame().set_linewidth(0.5)

# Save
plt.tight_layout()
plt.savefig(output_file, dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
print(f"✅ Image saved: {os.path.abspath(output_file)}")
plt.close(fig)

# Print stats
print(f"\n📊 Data Summary:")
print(f"   Dataset 1: {len(data1):,} points | Mean: {data1.mean():.3f} | Std: {data1.std():.3f}")
if data2 is not None:
    print(f"   Dataset 2: {len(data2):,} points | Mean: {data2.mean():.3f} | Std: {data2.std():.3f}")