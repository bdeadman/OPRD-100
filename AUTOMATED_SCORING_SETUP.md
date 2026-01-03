# Automated Scoring System - Setup Overview

This document provides an overview of the automated scoring system for OPRD-100 submissions.

## 🎯 What This Does

When someone submits their extracted reaction data via a pull request:
1. **Automated scoring** runs against the OPRD-100 ground-truth dataset
2. **Plots are generated** and saved in `results/submitter_name_timestamp/`
3. **Leaderboard updates** automatically in README.md
4. **PR comments** show the submitter their scores immediately

## 📂 File Structure

```
OPRD-100/
├── .github/workflows/
│   └── score_submission.yml          # GitHub Actions workflow
├── data/submissions/
│   ├── README.md                     # Submission instructions
│   └── example_submission.json       # Example format (ignored by CI)
├── scripts/
│   ├── run_scoring.py               # Run scoring analysis
│   ├── extract_scores.py            # Extract key metrics
│   ├── update_leaderboard.py        # Update README & leaderboard.json
│   └── README.md                    # Scripts documentation
├── leaderboard.json                 # Leaderboard data storage
└── README.md                        # Updated with leaderboard section
```

## 🔄 Workflow Process

### 1. User Submits Data
User creates a PR with `data/submissions/your_name.json` containing:
- Submitter information
- Method description
- Extracted reaction data in OPRD-100 format

### 2. GitHub Actions Triggers
When PR touches `data/submissions/**`, the workflow:
- Checks out the PR branch
- Sets up Python environment
- Installs dependencies

### 3. Scoring Runs
`run_scoring.py`:
- Loads submission data
- Initializes DataComparer with submission
- Computes similarity scores vs ground truth
- Generates all plots from example_scoring.ipynb
- Saves results to `results/submitter_timestamp/`

### 4. Metrics Extraction
`extract_scores.py`:
- Reads the CSV results
- Calculates mean, median, std for all categories
- Saves to `scores.json`

### 5. Leaderboard Update
`update_leaderboard.py`:
- Updates `leaderboard.json` with new scores
- Sorts by combined mean score
- Updates README.md table

### 6. Commit & Comment
Workflow:
- Commits results back to PR branch
- Posts score summary as PR comment
- Links to full leaderboard

## 🧪 Testing Locally

Before pushing to GitHub, test the system:

```bash
# Create test submission
python scripts/run_scoring.py \
  --submission-file data/submissions/example_submission.json \
  --output-dir results/test_run

# Extract scores
python scripts/extract_scores.py \
  --results-dir results/test_run \
  --output results/test_run/scores.json

# Update leaderboard
python scripts/update_leaderboard.py \
  --submitter "Test" \
  --scores-file results/test_run/scores.json \
  --pr-number 0
```

## 🔧 Configuration

### Modify Workflow Triggers
Edit `.github/workflows/score_submission.yml`:
```yaml
on:
  pull_request:
    paths:
      - 'data/submissions/**'  # Only trigger on submission files
```

### Adjust Scoring Metrics
Modify `scripts/run_scoring.py` to change:
- Which plots are generated
- Which metrics are computed
- Plot styling and parameters

### Customize Leaderboard Display
Edit `scripts/update_leaderboard.py` to change:
- Table columns
- Sorting criteria
- Formatting

## 📝 Submission Format

Submissions must be JSON files with:
```json
{
  "submitter_name": "string (required)",
  "method_description": "string (required)",
  "repository_url": "string (optional)",
  "reactions": [/* OPRD-100 format */] (required)
}
```

See `data/submissions/README.md` for full schema.

## 🔒 Permissions

The workflow needs:
- `contents: write` - to commit results back to PR
- `pull-requests: write` - to comment on PR

These are configured in the workflow file.

## 🐛 Troubleshooting

### Workflow doesn't trigger
- Check that files are in `data/submissions/`
- Ensure JSON is valid
- Check GitHub Actions permissions

### Scoring fails
- Verify submission JSON matches schema
- Check that required fields are present
- Ensure reaction data follows OPRD-100 format

### Leaderboard doesn't update
- Check for merge conflicts in README.md
- Verify `leaderboard.json` is valid JSON
- Check script permissions

## 🚀 Next Steps

1. **Enable GitHub Pages** (optional) for web leaderboard:
   - Settings → Pages → Source: main branch
   - Create `docs/index.html` with interactive table

2. **Add badges** to README showing stats:
   - Number of submissions
   - Top score
   - Last updated

3. **Email notifications** when new top score achieved

4. **Enhanced visualizations** comparing submissions

## 📚 Resources

- [GitHub Actions Documentation](https://docs.github.com/en/actions)
- [OPRD-100 Schema](data/submissions/README.md)
- [Example Scoring Notebook](notebooks/example_scoring.ipynb)
