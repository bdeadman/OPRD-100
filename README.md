# OPRD-100: Ground-Truth Reaction Dataset and Validation Similarity Metrics

OPRD-100 is a ground-truth dataset of reaction data curated from 100 papers in Organic Process Research & Development (OPRD). This repository also provides code and example notebooks to compute and visualize similarity metrics for validating automated data extraction on the same literature.

## Key features
- Curated ground-truth reaction data from 100 OPRD papers.
- Re-extraction validation subset and scoring pipeline.
- Multiple similarity metrics for comparing structured reaction data.
- Example notebooks with code and plots illustrating the metrics.

## Repository layout
- data/
  - OPRD-100/ — curated ground-truth data (JSON)
  - OPRD-100-raw/ — raw extractions used to form the ground truth
  - OPRD-100-validation/ — validation subset for re-extraction (30 locations)
- src/ — metric implementations and utilities
- notebooks/ — examples and visualizations of metrics and matching
- results — stores visualisation of reaction similarity metrics (initialised with human-extraced results as an example)
- README.md — this file

Note: Some folders may be created when you first run notebooks or scripts.

## Dataset overview
- Scope: 100 OPRD papers; reactions extracted from Schemes, Tables, and Experimentals.
- Format: Reactions processed to JSON for programmatic comparison.
- Intended use: Benchmarking and validating automated reaction data extraction.

## Validation methodology (summary)
A validation subset (OPRD-100-validation) of 30 data locations (10 Schemes, 10 Tables, 10 Experimentals) was randomly selected from OPRD-100-raw. Curators re-extracted reactions from these locations using the same protocol, avoiding articles they originally processed. For schemes and tables, all reactions were extracted; for experimentals, only the first procedure was considered (to handle multistep protocols common in Org. Process Res. Dev.).

Both original (OPRD-100) and validation JSON data were compared with the following metrics. Reaction pairing used the Hungarian algorithm (SciPy’s linear_sum_assignment) based on the average of metrics 2–8, termed total_similarity. If reaction counts differed, the lowest-scoring unmatched reactions were discarded from analysis.

## Similarity metrics
Unless stated otherwise, scores are binary (1/0) or Jaccard-style similarities computed on the exact reported text strings.

1) Reaction entry count
  - Equality of the number of reactions per location (for diagnostics; not in total_similarity).

2) SMILES similarity
  - String-based comparison of reported SMILES (exact match or set overlap where applicable).

3) Reaction step number equality
  - Binary equality of step counts.

4) Yield data similarity
  - Comparison of reported yields and their associations with products.

5) Reagent similarity
  - Product of two sub-scores:
    - Name match score
    - Amount-of-substance score

6) Solvent similarity
  - Product of two sub-scores:
    - Name match score
    - Amount-of-substance score

7) Reaction time similarity
  - 1 minus the normalized difference between reported times.

8) Temperature similarity
  - Binary comparison of reported min/max temperatures.

total_similarity = average of metrics 2–8 for a candidate reaction pair.

## Pairing strategy
- Optimal pairing via Hungarian algorithm using total_similarity as the cost/score.
- When counts differ, lowest-scoring unmatched reactions are excluded (affected 6 reactions in total during human validation).


## Compute Scoring metrics for your dataset!

  - Use the ``example_scoring.ipynb`` notebook in the notebooks folder, change the filepath to your own and simply run all cells.


## Reproducibility tips
- Keep raw, ground-truth, and validation JSONs immutable; write derived artifacts to a results/ folder.

## Contributing
- Open an issue to propose features, metrics, or bug fixes.

## Citation
If you use OPRD-100 or the validation metrics in a publication, please cite this repository. A BibTeX entry can be added here once the persistent identifier (e.g., DOI) is available.

## Contact
- For questions or feedback, please open an issue in this repository.