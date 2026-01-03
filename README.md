# OPRD-100: Ground-Truth Reaction Dataset and Validation Similarity Metrics

OPRD-100 is a ground-truth dataset of reaction data curated from 100 papers in Organic Process Research & Development (OPRD). This repository also provides code and example notebooks to compute and visualize similarity metrics for validating automated data extraction on the same literature.

## Key features
- Curated ground-truth reaction data from 100 OPRD papers.
- Re-extraction validation subset and scoring pipeline.
- Multiple similarity metrics for comparing structured reaction data.
- Example notebooks with code and plots illustrating the metrics.

# Compute Scoring metrics for your dataset!

There are two ways to evaluate your extraction method:

### Option 1: Use the Jupyter Notebook (Manual)
- Use the `example_scoring.ipynb` notebook in the notebooks folder
- Change the filepath to your own dataset
- Run all cells to generate scores and visualizations

### Option 2: Submit to the Leaderboard (Automated)
Submit your extracted data via pull request for automated scoring and public leaderboard inclusion:

1. **Extract reactions** from OPRD-100 papers using your method
2. **Format your data** following the [submission template](data/submissions/README.md)
3. **Create a pull request** with your submission file in `data/submissions/`
4. **Automated scoring** runs and adds your results to the leaderboard below

See detailed instructions in [data/submissions/README.md](data/submissions/README.md).

## 🏆 Leaderboard

Submissions are automatically scored against the OPRD-100 validation dataset. Scores represent similarity metrics between extracted and ground-truth data.

| Rank | Submitter | Combined | Experimental | Table | Scheme | Reactions | Date | Details |
|------|-----------|----------|--------------|-------|--------|-----------|------|----------|
| 1 | Test User | 0.882 | 0.817 | 0.892 | 0.847 | 189 | 2026-01-03 | [PR #0](../../pull/0) |

**Metrics explanation:**
- **Combined**: Average similarity across all reaction types (higher is better, max 1.0)
- **Experimental/Table/Scheme**: Average similarity for each data location type
- **Reactions**: Total number of reactions scored
## Reproducibility tips
- Keep raw, ground-truth, and validation JSONs immutable; write derived artifacts to a results/ folder.

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




## Contributing
- Open an issue to propose features, metrics, or bug fixes.

## Citation
If you use OPRD-100 or the validation metrics in a publication, please cite this repository. A BibTeX entry can be added here once the persistent identifier (e.g., DOI) is available.

## Contact
- For questions or feedback, please open an issue in this repository.