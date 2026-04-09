import os
import csv

# Define your reference files
ground_truth_file = "/storage/mi/malek01/ForschungsPraktikumRKI/data/human/sensandspe/gt_FDR5.tsv"
true_negatives_file = "/storage/mi/malek01/ForschungsPraktikumRKI/data/human/sensandspe/tn_FDR5.tsv"
folder_path = "/storage/mi/malek01/ForschungsPraktikumRKI/data/human/sensandspe/fdr5"  # Change this if your files are in a different folder

def load_peptides(filename):
    """Reads a one-column tsv into a set of strings."""
    peptides = set()
    with open(filename, 'r') as f:
        for line in f:
            if line.strip():
                peptides.add(line.split('\t')[0].strip())
    return peptides

def calculate_metrics(predicted, actual_pos, actual_neg):
    tp = len(predicted.intersection(actual_pos))
    fp = len(predicted.intersection(actual_neg))
    fn = len(actual_pos) - tp
    tn = len(actual_neg) - fp

    sensitivity = tp / len(actual_pos) if len(actual_pos) > 0 else 0
    specificity = tn / len(actual_neg) if len(actual_neg) > 0 else 0

    return {
        'TP': tp,
        'FP': fp,
        'TN': tn,
        'FN': fn,
        'Sensitivity': sensitivity,
        'Specificity': specificity
    }

# Load reference data
actual_positives = load_peptides(ground_truth_file)
actual_negatives = load_peptides(true_negatives_file)

# Find all TSV files in the folder, excluding reference files
exclude_files = {ground_truth_file, true_negatives_file}
prediction_files = [
    f for f in os.listdir(folder_path)
    if f.endswith('.tsv') and f not in exclude_files
]

# Store results
results = []

print(f"Analyzing {len(prediction_files)} files...\n")
print(f"{'File':<40} {'Sensitivity':<12} {'Specificity':<12}")
print("-" * 64)

for filename in sorted(prediction_files):
    filepath = os.path.join(folder_path, filename)
    predicted = load_peptides(filepath)
    metrics = calculate_metrics(predicted, actual_positives, actual_negatives)

    results.append({
        'File': filename,
        'TP': metrics['TP'],
        'FP': metrics['FP'],
        'TN': metrics['TN'],
        'FN': metrics['FN'],
        'Sensitivity': metrics['Sensitivity'],
        'Specificity': metrics['Specificity']
    })

    print(f"{filename:<40} {metrics['Sensitivity']:<12.4f} {metrics['Specificity']:<12.4f}")

# Save results to CSV
output_csv = "results_summary.csv"
with open(output_csv, 'w', newline='') as csvfile:
    fieldnames = ['File', 'TP', 'FP', 'TN', 'FN', 'Sensitivity', 'Specificity']
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(results)

print("-" * 64)
print(f"\nResults saved to {output_csv}")