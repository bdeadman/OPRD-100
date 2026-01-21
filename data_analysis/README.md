# OPRD-100 Data Analysis

This directory contains clean (using claude Sonnet 4.5), Jupyter notebooks and associated code for analyzing the OPRD-100 dataset.

## Directory Structure

```
data_analysis/
├── notebooks/           # Jupyter notebooks for analysis
│   ├── dataset_analysis.ipynb           # OPRD-100 dataset analysis
│   └── reaxys_comparison_clean.ipynb    # OPRD-100 vs Reaxys comparison
├── src/                 # Python source files
│   ├── PlottingFunctions.py             # Visualization utilities
│   ├── ReaxysComparison.py              # Reaxys data extraction
│   ├── Stereochemistry.py               # Stereochemistry analysis
│   ├── StereoChem.py                    # Stereochemistry utilities
│   ├── reaction_networks.py             # Reaction network generation
│   ├── molecule_presence_checker.py     # PubChem CAS number lookup
│   └── CHEM21_solvent_dict.json         # Solvent classification data
├── output/              # Generated outputs
│   ├── analysis_figures/                # Figures from dataset_analysis.ipynb
│   ├── reaxys_comparison_figures/       # Figures from reaxys_comparison_clean.ipynb
│   └── reaxys_comparison_data/          # CSV/Excel files from Reaxys comparison
├── data/                # Input data files
│   ├── OPRD-100.json                    # Main dataset file
│   └── Reaxys_data/
│       └── reaxys_100_articles_data.xml # Reaxys XML data
└── README.md            # This file
```

## Notebooks

### 1. dataset_analysis.ipynb

Comprehensive analysis of the OPRD-100 dataset covering:
- **Data Locations**: Where reaction data appears in articles (schemes, tables, experimental sections)
- **Yield Types**: Distribution of different yield metrics (isolated, assay, ee, dr, etc.)
- **Yield Distributions**: Statistical analysis of yields by data location
- **Stereochemistry**: Analysis of stereogenic centers gained/lost/preserved
- **Solvent Usage**: Frequency, volume, and CHEM21 environmental classification
- **Reaction Conditions**: Temperature and time distributions
- **Catalyst Usage**: Frequency of different metal catalysts

**Outputs**: 
- Figures saved to `output/analysis_figures/`
- All plots generated from `data/OPRD-100.json`

### 2. reaxys_comparison_clean.ipynb

Comparison between OPRD-100 and Reaxys datasets from the same 100 articles:
- **Data Location Comparison**: Comparison of where data is reported
- **Yield Type Comparison**: Types and distributions across datasets
- **Reaction Networks**: Molecule (node) and reaction (edge) overlap analysis
- **Stereochemistry Conflicts**: Detailed classification of stereochemical differences
- **Molecule Identity**: Classification as intermediates vs reagents (using CAS lookup)
- **Example Networks**: Visualization of reaction networks

**Outputs**:
- Figures saved to `output/reaxys_comparison_figures/`
- Data files (CSV/Excel) saved to `output/reaxys_comparison_data/`

## Installation

### Option 1: Using Conda (Recommended)

The easiest way to install all dependencies, especially rdkit and cairosvg which have system dependencies:

```bash
# Create the conda environment with Python 3.12
conda env create -f environment.yml

# Activate the environment
conda activate oprd100
```
## Running the Notebooks

1. **Activate Environment** (if using conda):
   ```bash
   conda activate oprd100
   ```

2. **Start Jupyter**:
   ```bash
   cd notebooks/
   jupyter notebook
   ```

3. **Run Notebooks**:
   - Open either notebook
   - Run cells sequentially (the notebooks handle path setup automatically)
   - Outputs will be generated in the `output/` directory

## Notes

- The notebooks use relative paths that work from the `notebooks/` directory
- Python modules are automatically imported from `../src/`
- All original analysis logic is preserved without refactoring
- Notebooks include comprehensive markdown documentation for each section

## Data Files

- **OPRD-100.json**: Main dataset containing 3878 reactions from 100 OPRD articles
- **reaxys_100_articles_data.xml**: Reaxys data extracted from the same 100 articles
- **CHEM21_solvent_dict.json**: Solvent classification data for environmental assessment
