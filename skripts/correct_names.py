import pandas as pd
import re

# Read the TSV file
df = pd.read_csv("/storage/mi/malek01/ForschungsPraktikumRKI/data/fragpipe/human5_1/ground_truth.tsv", sep="\t")

# Convert SpecId to Thermo format
def convert_to_thermo_id(spec_id):
    # Extract scan number (usually second or third number in the ID)
    # Example: OrbiElite04554.2345.2345.2 -> scan 2345
    numbers = re.findall(r'\d+', spec_id)
    if len(numbers) >= 2:
        scan_num = int(numbers[1])  # Remove leading zeros by converting to int
        return f"controllerType=0 controllerNumber=1 scan={scan_num}"
    return spec_id  # Fallback if format doesn't match

# Apply conversion
df["Spectrum"] = df["Spectrum"].apply(convert_to_thermo_id)

# Save back to TSV
df.to_csv("/storage/mi/malek01/ForschungsPraktikumRKI/data/fragpipe/human5_1/ground_truth.tsv", sep="\t", index=False)

print(f"Converted {len(df)} entries")
print(f"Example: {df['Spectrum'].iloc[0]}")