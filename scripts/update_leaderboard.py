"""
Script to update README.md leaderboard and leaderboard.json with new submission scores.
"""
import argparse
import json
import re
from datetime import datetime
import os


def update_leaderboard(submitter, scores_file, pr_number):
    """Update README.md leaderboard and leaderboard.json"""
    
    print(f"Updating leaderboard for {submitter}")
    
    with open(scores_file, 'r') as f:
        scores = json.load(f)
    
    # Update leaderboard.json
    leaderboard_path = 'leaderboard.json'
    try:
        with open(leaderboard_path, 'r') as f:
            leaderboard = json.load(f)
    except FileNotFoundError:
        leaderboard = {'entries': []}
    
    # Check if this submitter already exists and update or add new entry
    existing_entry = None
    for i, entry in enumerate(leaderboard['entries']):
        if entry['submitter'] == submitter:
            existing_entry = i
            break
    
    entry = {
        'submitter': submitter,
        'date': datetime.now().strftime('%Y-%m-%d'),
        'pr_number': pr_number,
        'scores': scores
    }
    
    if existing_entry is not None:
        # Update existing entry
        leaderboard['entries'][existing_entry] = entry
        print(f"Updated existing entry for {submitter}")
    else:
        # Add new entry
        leaderboard['entries'].append(entry)
        print(f"Added new entry for {submitter}")
    
    # Sort by combined_mean score (descending)
    leaderboard['entries'].sort(key=lambda x: x['scores']['combined_mean'], reverse=True)
    
    # Add rank to each entry
    for i, entry in enumerate(leaderboard['entries'], 1):
        entry['rank'] = i
    
    with open(leaderboard_path, 'w') as f:
        json.dump(leaderboard, f, indent=2)
    
    print(f"Leaderboard saved to {leaderboard_path}")
    
    # Update README.md
    update_readme_table(leaderboard)
    print("README.md updated")


def update_readme_table(leaderboard):
    """Insert/update leaderboard table in README.md"""
    
    readme_path = 'README.md'
    
    with open(readme_path, 'r') as f:
        readme = f.read()
    
    # Generate table
    table = "## 🏆 Leaderboard\n\n"
    table += "Submissions are automatically scored against the OPRD-100 validation dataset. Scores represent similarity metrics between extracted and ground-truth data.\n\n"
    table += "| Rank | Submitter | Combined | Experimental | Table | Scheme | Reactions | Date | Details |\n"
    table += "|------|-----------|----------|--------------|-------|--------|-----------|------|----------|\n"
    
    for entry in leaderboard['entries']:
        s = entry['scores']
        rank = entry['rank']
        submitter = entry['submitter']
        combined = f"{s['combined_mean']:.3f}"
        exp = f"{s['exp_mean']:.3f}" if s.get('num_experimental', 0) > 0 else "N/A"
        table_score = f"{s['table_mean']:.3f}" if s.get('num_table', 0) > 0 else "N/A"
        scheme = f"{s['scheme_mean']:.3f}" if s.get('num_scheme', 0) > 0 else "N/A"
        num_reactions = s.get('total_reactions', 0)
        date = entry['date']
        pr_num = entry['pr_number']
        
        table += f"| {rank} | {submitter} | {combined} | {exp} | {table_score} | {scheme} | {num_reactions} | {date} | [PR #{pr_num}](../../pull/{pr_num}) |\n"
    
    table += "\n**Metrics explanation:**\n"
    table += "- **Combined**: Average similarity across all reaction types (higher is better, max 1.0)\n"
    table += "- **Experimental/Table/Scheme**: Average similarity for each data location type\n"
    table += "- **Reactions**: Total number of reactions scored\n\n"
    
    # Replace or insert leaderboard section
    leaderboard_pattern = r'## 🏆 Leaderboard.*?(?=\n## |\Z)'
    if re.search(leaderboard_pattern, readme, re.DOTALL):
        readme = re.sub(leaderboard_pattern, table.rstrip(), readme, flags=re.DOTALL)
    else:
        # Insert before "Contributing" section if it exists, otherwise at the end
        contributing_pattern = r'\n## Contributing'
        if re.search(contributing_pattern, readme):
            readme = re.sub(contributing_pattern, f"\n{table}\n## Contributing", readme)
        else:
            readme += f"\n\n{table}"
    
    with open(readme_path, 'w') as f:
        f.write(readme)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Update leaderboard with new scores')
    parser.add_argument('--submitter', required=True, help='Submitter name')
    parser.add_argument('--scores-file', required=True, help='Path to scores JSON file')
    parser.add_argument('--pr-number', required=True, type=int, help='Pull request number')
    args = parser.parse_args()
    
    update_leaderboard(args.submitter, args.scores_file, args.pr_number)
