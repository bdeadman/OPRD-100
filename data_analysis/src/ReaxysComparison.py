import json
from lxml import etree
from rdkit import Chem
import re
import numpy as np
import pandas as pd

"""
Python file to compare reaxys XML files to JSON processed manually extracted reaction data.

** Logic: **
1. Load the JSON file containing manually extracted reaction data.
2. Parse the XML file from Reaxys.
3. Extract the reaxys data into the same JSON format as the manually extracted data.
4. Compare the reactions in the JSON file with those in the XML file:
    4.1 Extract relevant data from both sources (e.g., reactants, products, conditions).

5. Identify discrepancies in reaction data, such as missing reactions, different reactants/products, etc.
6. Output the results of the comparison, highlighting any differences found.
"""


class ReaxysDataExtractor:

    def __init__(self, xml_file_path):
        self.xml_file_path = xml_file_path

    def extract_smiles(self, xml_obj, tag):
        """
        Extract SMILES from the XML object based on the provided tag.
        """
        smiles = []
        for element in xml_obj.findall(tag):
            mol_block = element.text
            if mol_block:
                mol = Chem.MolFromMolBlock(mol_block)
                
                if mol:
                    smiles.append(Chem.MolToSmiles(mol))
        smiles = ".".join(smiles)
        return smiles

    def extract_reaction_data(self):
        """
        Extract reaction data from the Reaxys XML file and store it in a structured format.
        """
        reactions = []
        # Parse the XML file
        tree = etree.parse(self.xml_file_path)
        root = tree.getroot()
        # Iterate through each <reaction> element
        for reaction in root.findall('.//reaction'):
            reaction_data = {
                'Reaction': None,
                'Reference': None,
                'Reagents': None,
                'Solvents': None,
                'Yields': None,
                'Yield Other': None,
                'Temperatures': None,
                'Times': None,
                'Location': None,
                'Num Steps': None,
                'Num Stages': None,
                'Reaxys ID': None
            }
            # Create the Reaction SMILES
            reactants_smiles = self.extract_smiles(reaction, './/RY.RCT')
            products_smiles = self.extract_smiles(reaction, './/RY.PRO')
            reaction_smiles = reactants_smiles + ">>" + products_smiles
            reaction_data['Reaction'] = reaction_smiles
            # Extract the reference
            doi_element = reaction.find(".//CIT.DOI/hi")
            if doi_element is not None:
                # Remove the prefix '10.1021/' from the DOI
                stripped_doi = re.sub(r'10.1021/', '', doi_element.text.strip())
                reaction_data['Reference'] = stripped_doi
            # Extract reagents
            rgts = reaction.findall(".//RXD.RGT")
            if rgts:
                reaction_data['Reagents'] = {elem.text for elem in rgts}
            # Extract Temperature
            temps = reaction.findall(".//RXD.T")
            if temps:
                reaction_data['Temperatures'] = {elem.text for elem in temps}
            # Extract Solvent
            solvents = reaction.findall(".//RXD.SOL")
            if solvents:
                reaction_data['Solvents'] = {elem.text for elem in solvents}
            # Extract Times
            times = reaction.findall(".//RXD.TIM")
            if times:
                reaction_data['Times'] = {elem.text for elem in times}
            # Extract Yields
            yields = reaction.findall(".//RXD.YD")
            if yields:
                reaction_data['Yields'] = {elem.text for elem in yields}
            # Extract Yield Other
            yield_other = reaction.findall(".//RXD.YDO")
            if yield_other:
                reaction_data['Yield Other'] = {elem.text for elem in yield_other}
            # Extract Location
            location_element = reaction.findall(".//RXD.RXDES")
            if location_element:
                reaction_data['Location'] = {elem.text for elem in location_element}
            # Extract Num Steps
            num_steps_element = reaction.findall(".//RXD.STP")
            if num_steps_element:
                reaction_data['Num Steps'] = max(int(elem.text) for elem in num_steps_element)
            # Extract Num Stages
            num_stages_element = reaction.findall(".//RXD.SNR")
            if num_stages_element:
                reaction_data['Num Stages'] = max(int(elem.text) for elem in num_stages_element)
            # Extract Reaxys ID
            reaxys_id = reaction.findall(".//RX.ID")
            if reaxys_id:
                reaction_data['Reaxys ID'] = [elem.text for elem in reaxys_id][0]
            # Append the reaction data to the list
            reactions.append(reaction_data)
        return reactions
    
def load_json(json_file_path):
    """
    Load the JSON file containing manually extracted reaction data.
    """
    with open(json_file_path, 'r') as file:
        return json.load(file)


def CompareNumReactions(manual: pd.DataFrame, reaxys: pd.DataFrame):
    """
    Takes in two JSON files, one with manually extracted reactions and one with Reaxys XML data.
    Compares the number of reaction entries between the unique DOI refereneces in both files.
    Outputs the number of reactions in each file and the difference.
    """

    grouped_1 = manual.groupby('Reference').size()
    grouped_1.name = 'Count_manual'

    grouped_2 = reaxys.groupby('Reference').size()
    grouped_2.name = 'Count_reaxys'
    # merge the two grouped dataframes
    merged = pd.merge(grouped_1, grouped_2, left_index=True, right_index=True, how='inner')
    merged['Count_manual-reaxys'] = merged['Count_manual'] - merged['Count_reaxys']
    return merged

