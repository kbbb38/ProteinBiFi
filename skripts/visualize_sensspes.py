import matplotlib.pyplot as plt

# ============== INPUT DATA ==============
# Insert your two lists of values here
list1 = [0.9168, 0.9168, 0.9167, 0.9165, 0.9012, 0.5277, 0.2352, 0.0797, 0.0145, 0.0007]
list2 = [0.0000, 0.0001, 0.0004, 0.0019, 0.0459, 0.5844, 0.8394, 0.9486, 0.9912, 0.9998]

# X-axis values (labels)
x_values = [0.001, 0.005, 0.01, 0.02, 0.04, 0.08, 0.10, 0.12, 0.15, 0.20]  # or ['Jan', 'Feb', 'Mar', ...] for text labels

# ============== PLOT CONFIGURATION ==============
# Create the figure and axis
fig, ax = plt.subplots(figsize=(10, 6))

# Plot both lines
ax.plot(x_values, list1, label='Sensitivity', marker='o', linewidth=2, color='blue')
ax.plot(x_values, list2, label='Specificity', marker='s', linewidth=2, color='red')

# ============== LABELS AND TITLE ==============
# Name the axes
ax.set_xlabel('X-Axis Label', fontsize=12, fontweight='bold')
ax.set_ylabel('Y-Axis Label', fontsize=12, fontweight='bold')

# Add title
ax.set_title('Sensitivity vs Specificity, bin size = 0.2, FDR = 5%', fontsize=14, fontweight='bold', pad=15)

# Add legend
ax.legend(loc='best', fontsize=10, frameon=True)

# ============== OPTIONAL ENHANCEMENTS ==============
# Add grid for better readability
ax.grid(True, alpha=0.3, linestyle='--')

# Set x-ticks to match your x_values
ax.set_xticks(x_values)

# Adjust layout to prevent label cutoff
plt.tight_layout()

# Display the plot
plt.show()

# ============== SAVE OPTION (Optional) ==============
# To save the plot instead of displaying:
plt.savefig('/storage/mi/malek01/ForschungsPraktikumRKI/data/human/sensandspe/img/fdr5_bs02.png', dpi=300, bbox_inches='tight')