# Submission Instructions

Submit your automatically extracted reaction data to benchmark against the OPRD-100 ground-truth dataset!

## 📋 File Format

Create a JSON file in this directory with the following structure:

```json
{
  "submitter_name": "YourName_or_TeamName",
  "method_description": "Brief description of your extraction method (e.g., LLM-based extraction using GPT-4, Rule-based NLP pipeline, etc.)",
  "repository_url": "https://github.com/your/repo",
  "paper_url": "https://doi.org/your.paper (optional)",
  "contact_email": "your.email@example.com (optional)",
  "reactions": [
    {
      "DOI": "10.1021/acs.oprd.example",
      "data_location": "scheme_1",
      "location": {
        "Type": "scheme",
        "Num": "1"
      },
      "reactions": [
        {
          "step_num": 1,
          "reaction_smiles": "CC(=O)O.CCO>>CC(=O)OCC",
          "reagents": [
            {
              "reagent_name": "acetic acid",
              "reagent_smiles": "CC(=O)O",
              "amount_of_substance": {
                "value": 1.0,
                "unit": "mol",
                "reported_value": "1.0 mol"
              }
            }
          ],
          "solvents": [
            {
              "solvent_name": "ethanol",
              "solvent_smiles": "CCO",
              "amount_of_substance": {
                "value": 100,
                "unit": "mL",
                "reported_value": "100 mL"
              }
            }
          ],
          "time": {
            "value": 24,
            "unit": "h",
            "reported_value": "24 h"
          },
          "temperature": {
            "min": 80,
            "max": 80,
            "unit": "°C",
            "reported_value": "80 °C"
          },
          "yield_data": [
            {
              "product_smiles": "CC(=O)OCC",
              "yield_percent": 85.0,
              "reported_value": "85%"
            }
          ]
        }
      ]
    }
  ]
}
```

## 📝 Field Descriptions

### Required Fields
- `submitter_name`: Your name or team identifier (will be displayed on leaderboard)
- `method_description`: Brief description of your extraction approach
- `reactions`: Array of reaction data following the OPRD-100 schema

### Optional Fields
- `repository_url`: Link to your code repository
- `paper_url`: Link to associated paper or preprint
- `contact_email`: Your contact email

### Reaction Data Format
Follow the OPRD-100 schema structure. See the [validation_reactions.json](../validation_reactions.json) file for complete examples.

Key reaction fields:
- `DOI`: Paper DOI (must match OPRD-100 papers)
- `data_location`: Location identifier (e.g., "scheme_1", "table_2", "experimental_1")
- `location`: Dict with Type and Num
- `reactions`: Array of reaction objects with:
  - `step_num`: Reaction step number
  - `reaction_smiles`: Reaction SMILES string
  - `reagents`, `solvents`, `time`, `temperature`, `yield_data`: Structured extraction data

## 🚀 How to Submit

1. **Fork this repository**
2. **Extract reactions** from OPRD-100 papers using your method
3. **Create your submission file** as `data/submissions/your_name.json`
4. **Create a pull request** with:
   - Title: `[SUBMISSION] Your Name/Team`
   - Description: Brief overview of your method
5. **Automated evaluation** will:
   - Run scoring against the ground-truth dataset
   - Generate plots in `results/your_name_timestamp/`
   - Update the leaderboard in README.md
   - Comment results on your PR

## ✅ Validation Tips

Before submitting:
- Ensure your JSON is valid (use a JSON validator)
- Match the OPRD-100 schema structure
- Extract from the same papers/locations as the ground-truth dataset
- Include all required fields

## 📊 Scoring Metrics

Your submission will be evaluated using:
- **Reaction SMILES similarity**: Exact string matching
- **Reagent similarity**: Name and amount matching
- **Solvent similarity**: Name and amount matching  
- **Time/Temperature similarity**: Value comparison
- **Yield similarity**: Product and percentage matching

See [example_scoring.ipynb](../notebooks/example_scoring.ipynb) for detailed metric descriptions.

## 🤝 Questions?

Open an issue if you need help with the submission format or have questions about the evaluation process.
