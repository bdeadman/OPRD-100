"""
Script to extract key metrics from scoring results.
"""
import argparse
import json
import pandas as pd
import numpy as np
import os


def extract_scores(results_dir, output_file):
    """Extract key metrics from scoring results."""
    
    print(f"Extracting scores from: {results_dir}")
    
    # Load results CSV files
    exp_results = pd.read_csv(os.path.join(results_dir, "experimental_results.csv"))
    table_results = pd.read_csv(os.path.join(results_dir, "table_results.csv"))
    scheme_results = pd.read_csv(os.path.join(results_dir, "scheme_results.csv"))
    combined_results = pd.read_csv(os.path.join(results_dir, "combined_results.csv"))
    
    # Extract key statistics
    scores = {
        # Overall combined scores
        "combined_mean": float(combined_results['total_similarity'].mean()),
        "combined_median": float(combined_results['total_similarity'].median()),
        "combined_std": float(combined_results['total_similarity'].std()),
        
        # Experimental scores
        "exp_mean": float(exp_results['total_similarity'].mean()) if len(exp_results) > 0 else 0.0,
        "exp_median": float(exp_results['total_similarity'].median()) if len(exp_results) > 0 else 0.0,
        "exp_std": float(exp_results['total_similarity'].std()) if len(exp_results) > 0 else 0.0,
        
        # Table scores
        "table_mean": float(table_results['total_similarity'].mean()) if len(table_results) > 0 else 0.0,
        "table_median": float(table_results['total_similarity'].median()) if len(table_results) > 0 else 0.0,
        "table_std": float(table_results['total_similarity'].std()) if len(table_results) > 0 else 0.0,
        
        # Scheme scores
        "scheme_mean": float(scheme_results['total_similarity'].mean()) if len(scheme_results) > 0 else 0.0,
        "scheme_median": float(scheme_results['total_similarity'].median()) if len(scheme_results) > 0 else 0.0,
        "scheme_std": float(scheme_results['total_similarity'].std()) if len(scheme_results) > 0 else 0.0,
        
        # Individual metric means (combined)
        "reagent_name_mean": float(combined_results['reagent_name_similarity'].mean()),
        "reagent_amount_mean": float(combined_results['reagent_amount_similarity'].mean()),
        "reaction_smiles_mean": float(combined_results['reaction_smiles_similarity'].mean()),
        "solvent_mean": float(combined_results['solvent_similarity'].mean()),
        "time_mean": float(combined_results['time_similarity'].mean()),
        "temperature_mean": float(combined_results['temperature_similarity'].mean()),
        "yield_mean": float(combined_results['yield_similarity'].mean()),
        "reaction_steps_mean": float(combined_results['reaction_steps_similarity'].mean()),
        
        # Counts
        "total_reactions": len(combined_results),
        "num_experimental": len(exp_results),
        "num_table": len(table_results),
        "num_scheme": len(scheme_results)
    }
    
    # Save scores
    with open(output_file, 'w') as f:
        json.dump(scores, f, indent=2)
    
    print(f"Scores saved to: {output_file}")
    print(f"Combined mean similarity: {scores['combined_mean']:.3f}")
    print(f"Combined median similarity: {scores['combined_median']:.3f}")
    
    return scores


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Extract scores from results')
    parser.add_argument('--results-dir', required=True, help='Directory containing results CSV files')
    parser.add_argument('--output', required=True, help='Output JSON file path')
    args = parser.parse_args()
    
    extract_scores(args.results_dir, args.output)
