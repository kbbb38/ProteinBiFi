# ProteinBiFi (Protein Bit Filter)

A high-performance C++ tool for filtering and matching mass spectrometry spectra using bit set representations and Tanimoto similarity scoring.

## Overview

ProteinBiFi is designed to efficiently process spectral libraries and experimental mass spectrometry data. It uses a novel bit set representation for spectra, enabling fast similarity calculations and filtering based on Tanimoto scores.

## Features

- **Fast Spectrum Matching**: Uses bit set representations for efficient spectrum comparison
- **Tanimoto Similarity Scoring**: Calculates similarity between experimental and library spectra
- **Cutoff Analysis**: Analyzes ground truth matches to determine optimal filtering thresholds
- **Multi-format Support**: Reads both MGF (experimental) and MSP (library) file formats
- **Parallel Processing**: Optimized with multi-threading support
- **Ground Truth Validation**: Supports validation against known ground truth data

## Requirements

- CMake ≥ 3.10
- C++20 compatible compiler (GCC, Clang, MSVC)
- Python 3.x (for analysis scripts, optional)

### Python Dependencies (for analysis scripts)

```bash
pip install numpy matplotlib scipy
```

## Building

```bash
# Clone the repository
git clone <repository-url>
cd ProteinBiFi

# Create build directory
mkdir build && cd build

# Configure with CMake
cmake ..

# Build (Release mode with optimizations)
make -j$(nproc)

# Optional: Build with AddressSanitizer for debugging
cmake -DENABLE_ASAN=ON ..
make -j$(nproc)
```

## Usage

### Basic Usage

```bash
./ProteinBiFi -l <library_path> -e <experimental_path> -o <output_dir> -r <resolution> -g <ground_truth_path>
```

### Command Line Options

| Option | Short | Description | Required |
|--------|-------|-------------|----------|
| `--library` | `-l` | Path to spectral library file (MSP format) or directory | Yes |
| `--experimental` | `-e` | Path to experimental file (MGF format) or directory | Yes |
| `--out` | `-o` | Output directory for results | Yes |
| `--resolution` | `-r` | Resolution for bit sets (float value) | Yes |
| `--ground-truth` | `-g` | Path to ground truth file for validation | Yes |

### Example

```bash
./ProteinBiFi \
    -l /path/to/library.msp \
    -e /path/to/experimental.mgf \
    -o ./results \
    -r 1.0 \
    -g /path/to/ground_truth.txt
```

## Output

The program generates the following outputs in the specified output directory:

- **Tanimoto Scores**: Similarity scores between matched spectra
- **Filtered Spectra**: Spectra that pass the cutoff threshold
- **Cutoff Analysis Results**: Statistical analysis of ground truth matches

## Analysis Scripts

The `skripts/` directory contains Python scripts for post-processing and visualization:

### analyse_scores.py
Visualizes Tanimoto score distributions with density plots.

```python
# Configure input/output paths in the script
python skripts/analyse_scores.py
```

### calculate_metrics.py
Computes performance metrics from the results.

```python
python skripts/calculate_metrics.py
```

### extract_peptides.py
Extracts peptide information from the results.

```python
python skripts/extract_peptides.py
```

### visualize_sensspes.py
Creates sensitivity/specificity visualizations.

```python
python skripts/visualize_sensspes.py
```

### correct_names.py / correct_names_9MM.py
Utilities for correcting spectrum names in output files.

```python
python skripts/correct_names.py
```

## Project Structure

```
ProteinBiFi/
├── ProteinBiFi.cpp          # Main entry point
├── CMakeLists.txt           # CMake build configuration
├── include/
│   └── CLI11.hpp            # CLI11 command-line parsing library
├── src/
│   ├── SpectrumBitSet.h/cpp # Core bit set spectrum representation
│   ├── ExperimentalSpectra.h/cpp  # Experimental spectrum handling
│   ├── LibrarySpectra.h/cpp       # Library spectrum handling
│   ├── CutoffAnalysis.h/cpp       # Cutoff threshold analysis
│   ├── Config.h             # Configuration structures
│   └── Constants.h          # Application constants
├── skripts/                 # Python analysis and visualization scripts
└── build/                   # Build directory (generated)
```

## Algorithm

1. **Loading**: Experimental (MGF) and library (MSP) spectra are loaded into memory
2. **Bit Set Conversion**: Spectra are converted to compact bit set representations
3. **Cutoff Analysis**: Ground truth matches are analyzed to determine optimal thresholds
4. **Matching**: Tanimoto scores are calculated between all spectrum pairs
5. **Filtering**: Spectra are filtered based on the determined cutoff values
6. **Output**: Results are written to the specified output directory

## Performance

- Optimized with `-Ofast` compilation flag
- Multi-threaded execution (`-pthread`)
- Efficient bit set operations for fast similarity calculations
- Memory-efficient spectrum storage

## Debugging

Enable AddressSanitizer for memory error detection:

```bash
cd build
cmake -DENABLE_ASAN=ON ..
make -j$(nproc)
```

## License

This project uses [CLI11](https://github.com/CLIUtils/CLI11) for command-line parsing (BSD-3-Clause license).

## Citation

If you use ProteinBiFi in your research, please cite appropriately.

## Contact

For issues, questions, or contributions, please open an issue on the project repository.
