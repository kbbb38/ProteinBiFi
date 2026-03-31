# ProteinBiFi - Protein Bit Filter

ProteinBiFi is a high-performance C++ tool for filtering and analyzing mass spectrometry data using bit set representations of spectra. It efficiently matches experimental spectra against spectral libraries and applies cutoff analysis for quality filtering.

## Features

- **Fast Spectrum Matching**: Uses bit set representations for efficient spectrum comparison
- **Cutoff Analysis**: Statistical analysis to determine optimal filtering thresholds
- **Tanimoto Scoring**: Calculates similarity scores between spectra
- **Multi-threaded**: Built with pthread support for parallel processing
- **Flexible Input/Output**: Supports various file formats for experimental and library spectra

## Requirements

- CMake ≥ 3.10
- C++20 compatible compiler (GCC, Clang, etc.)
- Python 3 (for analysis scripts)
  - NumPy
  - Matplotlib
  - SciPy
  - Pandas

## Building

```bash
mkdir build
cd build
cmake ..
make
```

### Optional: Enable AddressSanitizer

For debugging memory issues:

```bash
mkdir build
cd build
cmake -DENABLE_ASAN=ON ..
make
```

## Usage

```bash
./ProteinBiFi -l <library> -e <experimental> -o <output> -r <resolution> [options]
```

### Required Arguments

| Argument | Short | Description |
|----------|-------|-------------|
| `--library` | `-l` | Path to spectral library file |
| `--experimental` | `-e` | Path to experimental file or directory |
| `--out` | `-o` | Output directory |
| `--resolution` | `-r` | Resolution for the bit sets (float) |

### Optional Arguments

| Argument | Short | Description |
|----------|-------|-------------|
| `--filter_experimental` | `-f` | Filter experimental spectra instead of spectra in the search library |

### Example

```bash
./ProteinBiFi -l library.mgf -e experimental.mgf -o results/ -r 0.01
```

## Project Structure

```
/workspace
├── CMakeLists.txt          # Build configuration
├── ProteinBiFi.cpp         # Main application entry point
├── include/
│   └── CLI11.hpp           # Command-line parsing library
├── src/
│   ├── Config.h            # Configuration structures
│   ├── Constants.h         # Application constants
│   ├── SpectrumBitSet.h    # Bit set representation for spectra
│   ├── SpectrumBitSet.cpp
│   ├── ExperimentalSpectra.h # Experimental spectrum handling
│   ├── ExperimentalSpectra.cpp
│   ├── LibrarySpectra.h    # Library spectrum handling
│   ├── LibrarySpectra.cpp
│   ├── CutoffAnalysis.h    # Statistical cutoff analysis
│   └── CutoffAnalysis.cpp
└── skripts/
    ├── analyse_scores.py   # Tanimoto score distribution analysis
    └── correct_names.py    # Spectrum ID format conversion
```

## Scripts

### analyse_scores.py

Analyzes and visualizes Tanimoto score distributions from ProteinBiFi output.

**Requirements:**
```bash
pip install numpy matplotlib scipy
```

**Configuration:**
Edit the script to set input/output paths:
```python
input_file_1 = 'path/to/scores.txt'
input_file_2 = None  # Optional second dataset for comparison
output_file = 'path/to/output.png'
```

**Usage:**
```bash
python skripts/analyse_scores.py
```

### correct_names.py

Converts spectrum IDs to Thermo Scientific format for compatibility with ground truth files.

**Requirements:**
```bash
pip install pandas
```

**Usage:**
Edit the script to set the input TSV file path, then run:
```bash
python skripts/correct_names.py
```

## Output

ProteinBiFi generates:
- Filtered spectrum files in the specified output directory
- Match statistics and scoring information
- Cutoff analysis results

## Performance

The tool is optimized for performance with:
- `-Ofast` optimization level
- Multi-threading support via pthread
- Efficient bit set operations for spectrum comparison

Typical processing times depend on:
- Number of spectra
- Resolution setting
- Available CPU cores

## License

This project is part of a research practicum at ForschungsPraktikumRKI.

## Acknowledgments

- Uses [CLI11](https://github.com/CLIUtils/CLI11) for command-line parsing
- Developed for mass spectrometry data analysis in proteomics research
