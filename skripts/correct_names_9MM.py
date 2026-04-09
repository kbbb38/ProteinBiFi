#!/usr/bin/env python3
"""
Process a TSV file by replacing Spectrum IDs with their corresponding 
spectrum titles from an mzML file, using the numeric index for lookup.
"""

import csv
import xml.etree.ElementTree as ET
from pathlib import Path
import sys


MZML_NS = '{http://psi.hupo.org/ms/mzml}'


def extract_index_from_spectrum_id(spectrum_id: str) -> str | None:
    """Extract numeric index from spectrum ID like '9MM_FASP.00372.00372.4'."""
    parts = spectrum_id.strip().split('.')
    if len(parts) >= 3:
        index_str = parts[2].lstrip('0')
        return index_str if index_str else '0'
    return None


def build_index_to_title_map(mzml_path: str) -> dict[str, str]:
    """Parse mzML file and build index → spectrum title mapping."""
    index_to_title = {}
    print(f"  Parsing mzML file: {mzml_path}")
    
    tree = ET.parse(mzml_path)
    root = tree.getroot()
    
    for spectrum in root.iter(f'{MZML_NS}spectrum'):
        index = spectrum.get('index')
        if index is None:
            continue
        for cvparam in spectrum.findall(f'{MZML_NS}cvParam'):
            if cvparam.get('name') == 'spectrum title':
                title = cvparam.get('value')
                if title:
                    index_to_title[index] = title
                break
    
    print(f"  Found {len(index_to_title)} spectra with titles")
    return index_to_title


def process_tsv_file(
    input_tsv: str, 
    mzml_file: str, 
    output_tsv: str,
    spectrum_col: str = 'Spectrum',
    peptide_col: str = 'Peptide'
) -> None:
    """Main processing function with robust TSV handling."""
    
    # ⚠️ Safety check: don't overwrite input file
    if Path(input_tsv).resolve() == Path(output_tsv).resolve():
        print(f"ERROR: Input and output files are the same: {input_tsv}")
        print("Please specify a different output file to avoid data loss!")
        sys.exit(1)
    
    # Step 1: Build mzML index map
    print(f"Building spectrum index map from mzML file...")
    index_to_title = build_index_to_title_map(mzml_file)
    
    if not index_to_title:
        print("WARNING: No spectrum titles found in mzML file.")
    
    # Step 2: Read and process TSV with robust encoding/delimiter handling
    print(f"Processing TSV file: {input_tsv}")
    
    rows_processed = 0
    rows_updated = 0
    
    # Use utf-8-sig to handle BOM, explicitly set delimiter to tab
    with open(input_tsv, 'r', newline='', encoding='utf-8-sig') as infile:
        # Read first line to debug header
        first_line = infile.readline()
        print(f"  First line of TSV (raw): {repr(first_line[:100])}")
        
        # Reset file pointer
        infile.seek(0)
        
        # Explicitly use tab delimiter, skip auto-detection
        reader = csv.DictReader(infile, delimiter='\t')
        
        if reader.fieldnames is None:
            print("ERROR: Could not read TSV header. File may be empty or malformed.")
            print(f"  First few bytes of file: {repr(first_line[:200])}")
            sys.exit(1)
        
        print(f"  Detected columns: {reader.fieldnames}")
        
        if spectrum_col not in reader.fieldnames:
            # Try to find a close match (case-insensitive)
            candidates = [c for c in reader.fieldnames if c.lower() == spectrum_col.lower()]
            if candidates:
                print(f"  Warning: '{spectrum_col}' not found, but found similar: {candidates}")
                spectrum_col = candidates[0]
            else:
                raise ValueError(
                    f"Column '{spectrum_col}' not found in TSV. "
                    f"Available columns: {reader.fieldnames}"
                )
        
        # Read all rows into memory first (safer than streaming when overwriting)
        rows = list(reader)
    
    # Process rows
    for row in rows:
        original_spectrum = row[spectrum_col].strip()
        index = extract_index_from_spectrum_id(original_spectrum)
        
        if index and index in index_to_title:
            row[spectrum_col] = index_to_title[index]
            rows_updated += 1
        elif index:
            print(f"  Warning: Index '{index}' (from '{original_spectrum}') not found in mzML")
        
        rows_processed += 1
    
    # Write output
    with open(output_tsv, 'w', newline='', encoding='utf-8') as outfile:
        writer = csv.DictWriter(outfile, fieldnames=reader.fieldnames, delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)
    
    print(f"✓ Output written to: {output_tsv}")
    print(f"  Processed: {rows_processed} rows, Updated: {rows_updated} spectrum IDs")


def main():
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Map TSV Spectrum IDs to mzML spectrum titles via index lookup'
    )
    parser.add_argument('input_tsv', help='Input TSV file path')
    parser.add_argument('mzml_file', help='Reference mzML file path')
    parser.add_argument('output_tsv', help='Output TSV file path (must differ from input!)')
    parser.add_argument('--spectrum-col', default='Spectrum', 
                       help='Name of Spectrum column (default: Spectrum)')
    parser.add_argument('--peptide-col', default='Peptide',
                       help='Name of Peptide column (default: Peptide)')
    
    args = parser.parse_args()
    
    # Validate input files
    for path in [args.input_tsv, args.mzml_file]:
        if not Path(path).exists():
            print(f"ERROR: File not found: {path}")
            sys.exit(1)
    
    # Ensure output directory exists
    Path(args.output_tsv).parent.mkdir(parents=True, exist_ok=True)
    
    process_tsv_file(
        input_tsv=args.input_tsv,
        mzml_file=args.mzml_file,
        output_tsv=args.output_tsv,
        spectrum_col=args.spectrum_col,
        peptide_col=args.peptide_col
    )
    return 0


if __name__ == '__main__':
    sys.exit(main())