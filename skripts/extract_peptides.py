input_file = "/storage/mi/malek01/ForschungsPraktikumRKI/data/fragpipe/human5_1/ground_truth.tsv"
output_file = "/storage/mi/malek01/ForschungsPraktikumRKI/data/human/sensandspe/gt_FDR5.tsv"

try:
    with open(input_file, 'r', encoding='utf-8') as f_in:
        lines = f_in.readlines()
    
    print(f"Total lines read: {len(lines)}")
    if len(lines) > 0:
        print(f"First line (header): '{lines[0].strip()}'")
    if len(lines) > 1:
        print(f"Second line (data): '{lines[1].strip()}'")
    
    if len(lines) <= 1:
        print("Error: File appears to have only a header or is empty.")
    else:
        with open(output_file, 'w', encoding='utf-8') as f_out:
            for line in lines[1:]:  # Skip header
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    f_out.write(parts[1] + '\n')
        print(f"Done! Peptides saved to {output_file}")

except FileNotFoundError:
    print(f"Error: File '{input_file}' not found in current directory.")
    print("Current working directory:", __import__('os').getcwd())
except Exception as e:
    print(f"Unexpected error: {e}")