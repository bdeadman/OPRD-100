
import json
import re
from collections import defaultdict
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdChemReactions
from rdkit.Chem import rdChemReactions
from rdkit.Chem import rdCIPLabeler 
from rxnmapper import RXNMapper


class ReactionParser:

    """
    Parse reactions in JSON file and count stereochemistry changes.
    """
    def __init__(self, json_file_path):
        self.json_file_path = json_file_path
        self.data = self.load_json(json_file_path)
        self.total_num_centres = 0
    
    def map_reaction(self, reaction_smiles)->None:
        """
        Map reaction SMILES using RXNMapper.
        """
        # Check if the reaction SMILES is valid
        reaction = rdChemReactions.ReactionFromSmarts(reaction_smiles)
        if not reaction:
            raise ValueError("Invalid reaction SMILES")
        # Map the reaction using RXNMapper
        rxn_mapper = RXNMapper()
        mapped_reaction = rxn_mapper.get_attention_guided_atom_maps([reaction_smiles])
        if not mapped_reaction:
            raise ValueError("Failed to map the reaction")
        # Add the mapped reaction to the data dictionary
        mapped_reaction= mapped_reaction[0].get("mapped_rxn")
        return mapped_reaction

    def get_mapped_rxn_smiles(self, reaction_data):
            # Parse the JSON string to get the reaction SMILES
            try:
                return self.map_reaction(reaction_data.get("Reaction"))
            except KeyError:
                raise KeyError("Key 'Reaction' not found in the JSON entry")
            except Exception as e:
                raise Exception(f"An error occurred while getting reaction SMILES: {e}")

    def classify_centre(self, centre):
        """
        Classify a stereocentre based on its properties.
        """
        smiles, root, controlling, specified, descriptor, type_, (cip_type, cip_label) = centre
        if cip_label:
            return (root, cip_label)
        if specified == "Unspecified" and descriptor == "Fused ring":
            return (root, "unspecified Fused-ring")
        if specified == "Unspecified" and type_ == "Atom_Tetrahedral":
            return (root, "Unspecified Tetrahedral Atom")
        if specified == "Unspecified" and type_ == "Bond_Double":
            return (root, "Unspecified Double Bond")
         

    def get_stereochemistry_tuples(self, mol):
        """
        Get stereochemistry tuples from a molecule.
        """
        try:
            checker = CheckStereochemistry(mol)
        except Exception as e:
            raise Exception(f"An error occurred while getting stereochemistry tuples: {e}")
        if checker.stereocentres is None:
            return set()
        # clean up the tuples to remove redundant info
        centres = set()
        for centre in checker.stereocentres:
            centre_class_tuple = self.classify_centre(centre)
            centres.add(centre_class_tuple)
        return centres
    
    def get_reaction_obj(self, smiles):
        # Parse the reaction SMILES and create a reaction object
        reaction = Chem.rdChemReactions.ReactionFromSmarts(smiles)
        if not reaction:
            raise ValueError("Invalid reaction SMILES")
        return reaction


    def split_reaction(self, reaction):
        # Split the reaction SMILES into reactants and products
        reactants = reaction.GetReactants()
        products = reaction.GetProducts()
        if not reactants or not products:
            raise ValueError("No reactants or products found in the reaction")
        return reactants, products
    
    def parse_reaction(self, reaction_data):
        """
        Parse a reaction and get the stereochemistry changes.
        """
        try:
            unmapped_smiles = reaction_data.get("Reaction")
            if len(unmapped_smiles.strip()) > 512:
                return unmapped_smiles,set(),set()
            rxn_smiles = self.get_mapped_rxn_smiles(reaction_data)
            reaction = self.get_reaction_obj(rxn_smiles)
            reactants, products = self.split_reaction(reaction)
            reactant_centres = set()
            for reactant in reactants:
                reactant_centres.update(self.get_stereochemistry_tuples(reactant))
            product_centres = set()
            for product in products:
                product_centres.update(self.get_stereochemistry_tuples(product))
            return rxn_smiles,reactant_centres, product_centres
        except Exception as e:
            raise Exception(f"An error occurred while parsing the reaction: {e}")
    
    def _extract_centre_class(self, centre):
        """
        Extract the root and class from a stereocentre tuple.
        """
        root, cls = centre
        return cls

    def compare_stereocentres(self, reactant_centres, product_centres):
        """
        Compare stereocentres between reactants and products.
        """
        gained = product_centres - reactant_centres
        lost = reactant_centres - product_centres
        unchanged = reactant_centres & product_centres
        # clean up the tuples to remove redundant info
        gained_counts = defaultdict(int)
        for centre in gained:
            cls = self._extract_centre_class(centre)
            gained_counts[cls] += 1
        lost_counts = defaultdict(int)
        for centre in lost:
            cls = self._extract_centre_class(centre)
            lost_counts[cls] += 1
        unchanged_counts = defaultdict(int)
        for centre in unchanged:
            cls = self._extract_centre_class(centre)
            unchanged_counts[cls] += 1
        return dict(gained_counts), dict(lost_counts), dict(unchanged_counts)
    
    def load_json(self, json_file_path):
        # Open the JSON file and load it
        try:
            with open(json_file_path, "r") as file:
                data = json.load(file)
                return data
        except FileNotFoundError:
            raise FileNotFoundError(f"File not found: {json_file_path}")
        except json.JSONDecodeError as e:
            raise ValueError(f"Invalid JSON file: {e}")
        except Exception as e:
            raise Exception(f"An error occurred: {e}")

    def get_stereochemistry_data(self):
        """
        debugging function to get stereochemistry data for all reactions.
        """
        results = []
        for entry in self.data:
            try:
                rxn_smiles, reactant_centres, product_centres = self.parse_reaction(entry)
                gained, lost, unchanged = self.compare_stereocentres(reactant_centres, product_centres)
                result = {
                    "Reaction": entry.get("Reaction"),
                    "Mapped Reaction": rxn_smiles,
                    "Gained": gained,
                    "Lost": lost,
                    "Unchanged": unchanged
                }
                results.append(result)
            except Exception as e:
                print(f"An error occurred while processing an entry: {e}")
        return results

    def canonical_smiles(self, smi: str) -> str | None:
        """Return RDKit-canonical SMILES (or None)."""
        m = Chem.MolFromSmiles(smi)
        return Chem.MolToSmiles(m, canonical=True) if m else None

    def smiles_to_inchikey(self, smi: str) -> str | None:
        """Return InChIKey for a SMILES (or None)."""
        m = Chem.MolFromSmiles(smi)
        return Chem.MolToInchiKey(m) if m else None

    def parse_reaction_smiles(self, rxn: str) -> tuple[list[str], list[str]]:
        """Split reaction SMILES into canonicalized reactant/product SMILES lists."""
        reactants, products = rxn.split(">>")
        rsmi = [s.strip() for s in reactants.split(".") if s.strip()]
        psmi = [s.strip() for s in products.split(".") if s.strip()]
        rsmi = [self.canonical_smiles(s) for s in rsmi]
        psmi = [self.canonical_smiles(s) for s in psmi]
        rsmi = [s for s in rsmi if s]
        psmi = [s for s in psmi if s]
        return rsmi, psmi
    
    def rxn_key(self, rxn):
        """
        Create a unique key for a reaction from InChIKeys of reactants/products
        computed directly from canonical SMILES (no InChI/SMILES round-trips).
        """
        rsmi, psmi = self.parse_reaction_smiles(rxn)
        rkeys = [self.smiles_to_inchikey(s) for s in rsmi]
        pkeys = [self.smiles_to_inchikey(s) for s in psmi]
        rkeys = [k for k in rkeys if k]
        pkeys = [k for k in pkeys if k]
        if not rkeys or not pkeys:
            return None
        reactants_key = ".".join(sorted(rkeys))
        products_key = ".".join(sorted(pkeys))
        return f"{reactants_key}>>{products_key}"

    def check_stereo_yield(self, reaction_data):
        """
        Check if a reaction has stereochemistry yield data.
        """
        stereo_yield_types = ["Selectivity","er","er R:S","er S:R","ee","ee (R)","ee (S)","Purity","dr","E:Z","Z:E"]
        try:
            steps = reaction_data.get("Steps")
            last_step = steps[-1]
            yield_info = last_step.get("Yield")
            if yield_info is None:
                return False
            has_stereo_yield = False
            for type_ in stereo_yield_types:
                data = yield_info.get(type_)
                if isinstance(data, list):
                    if any(d is not None for d in data):
                        return True
                elif data is not None:
                    return True
            # if any(yield_info[type_] is not None for type_ in stereo_yield_types):
            #     return True
        except Exception as e:
            print(f"An error occurred while checking stereochemistry yield: {e}")

    def count_stereochemistry_changes(self, unique_only=False, yield_check=False):
        """
        Count stereochemistry changes for all reactions in the JSON file.
        if unique_only is true, return count for only unique reactions
        """
        gained = defaultdict(int)
        lost = defaultdict(int)
        unchanged = defaultdict(int)
        gained_has_stereo_yield = defaultdict(int)
        lost_has_stereo_yield = defaultdict(int)
        unchanged_has_stereo_yield = defaultdict(int)
        seen_keys = set()
        for entry in self.data:
            has_stereo_yield = self.check_stereo_yield(entry)
            try:
                rxn_smiles, reactant_centres, product_centres = self.parse_reaction(entry)
                if unique_only:
                    key = self.rxn_key(rxn_smiles)
                    if not key or key in seen_keys:
                        continue
                    seen_keys.add(key)
                g, l, u = self.compare_stereocentres(reactant_centres, product_centres)
                for cls, count in g.items():
                    if yield_check and has_stereo_yield:
                        gained_has_stereo_yield[cls] += count
                    gained[cls] += count
                for cls, count in l.items():
                    if yield_check and has_stereo_yield:
                        lost_has_stereo_yield[cls] += count
                    lost[cls] += count
                for cls, count in u.items():
                    if yield_check and has_stereo_yield:
                        unchanged_has_stereo_yield[cls] += count
                    unchanged[cls] += count
            except Exception as e:
                print(f"An error occurred while processing an entry: {rxn_smiles}")
        if yield_check:
            return (dict(gained), dict(lost), dict(unchanged),
                    dict(gained_has_stereo_yield),
                    dict(lost_has_stereo_yield),
                    dict(unchanged_has_stereo_yield))
        return dict(gained), dict(lost), dict(unchanged),

class CheckStereochemistry:
    
    def __init__(self, mol):
        # self.mol = mol
        # self.smiles = Chem.MolToSmiles(mol)
        # try:
        #     rdCIPLabeler.AssignCIPLabels(self.mol)
        # except Exception as e:
        #     print(f"Failed to assign CIP labels to smiles {self.smiles}: {e}")

         # Work on a copy
        self.original_mol = mol

        # 1. Convert to a clean, fully sanitized molecule (ReactionFromSmarts can yield query atoms)
        try:
            tmp_smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
            self.mol = Chem.MolFromSmiles(tmp_smiles)  # fresh standard molecule
        except Exception:
            # Fallback: take the original if conversion fails
            self.mol = Chem.Mol(mol)

        # 2. Ensure property cache & implicit valences are computed
        try:
            self.mol.UpdatePropertyCache(strict=False)
            Chem.SanitizeMol(self.mol)
        except Exception:
            # Try a lighter sanitization if full sanitization fails
            try:
                self.mol.UpdatePropertyCache(strict=False)
            except Exception:
                pass

        # 3. Assign chiral tags from 3D / structure (needed before CIP sometimes)
        try:
            Chem.AssignAtomChiralTagsFromStructure(self.mol, replaceExistingTags=True)
        except Exception:
            pass

        self.smiles = Chem.MolToSmiles(self.mol, isomericSmiles=True)

        # 4. Attempt CIP assignment with defensive fallbacks
        try:
            rdCIPLabeler.AssignCIPLabels(self.mol)
        except Exception as e:
            # Retry after forcing another re-build
            try:
                rebuilt = Chem.MolFromSmiles(self.smiles)
                rebuilt.UpdatePropertyCache(strict=False)
                Chem.AssignAtomChiralTagsFromStructure(rebuilt, replaceExistingTags=True)
                rdCIPLabeler.AssignCIPLabels(rebuilt)
                self.mol = rebuilt
            except Exception as e2:
                print(f"Failed to assign CIP labels to smiles {self.smiles}: {e2}")

        self.rxn_smiles = None
        self.reaction = None
        self.reactants = None
        self.products = None
        self.stereocentres = set()
        self.get_rdkit_stereocentres()
        
        #self.get_sulphoxide_centres()
        #self.get_phosphine_centres()
        #self.get_phosphate_centres()



    def _map_or_index(self, atom_idx:int):
        """
        Return atom map number if non-zero else the original atom index.
        """
        if atom_idx < 0 or atom_idx >= self.mol.GetNumAtoms():
            return atom_idx
        a = self.mol.GetAtomWithIdx(atom_idx)
        amap = a.GetAtomMapNum()
        return amap if amap != 0 else atom_idx

    def get_rdkit_stereocentres(self):
        """
        Populate self.stereocentres with tuples:
        (smiles, root_atom_or_bond_tuple, controlling_atoms_tuple, specified, descriptor, type, ('CIP', cip_label))
        Root and controlling atom indices are replaced by atom map numbers when present (non‑zero),
        otherwise the original indices are retained.
        """
        self.stereocentres = set()
        si = Chem.FindPotentialStereo(self.mol)
        for element in si:
            specified = str(element.specified)
            descriptor = str(element.descriptor)
            type_ = str(element.type)
            cip_label = None

            # Raw controlling atoms (may contain large sentinels like 1000+ we ignore)
            raw_ca = [int(i) for i in element.controllingAtoms if 0 <= int(i) < self.mol.GetNumAtoms()]
            mapped_controlling = tuple(self._map_or_index(i) for i in raw_ca)

            # Root centre
            root = element.centeredOn  # for atoms OR bond index (for double bonds)
            if type_ == "Bond_Double":
                bond = self.mol.GetBonds()[root]
                if bond.GetBondTypeAsDouble() == 2.0:
                    a1 = bond.GetBeginAtomIdx()
                    a2 = bond.GetEndAtomIdx()
                    root_mapped = (self._map_or_index(a1), self._map_or_index(a2))
                    # CIP E/Z from bond property if specified
                    if specified == "Specified" and bond.HasProp("_CIPCode"):
                        cip_label = bond.GetProp("_CIPCode")
                else:
                    # Skip non true double
                    continue
            else:
                root_mapped = self._map_or_index(root)
                # CIP R/S for tetrahedral if specified
                if specified == "Specified" and type_ == "Atom_Tetrahedral":
                    atom = self.mol.GetAtomWithIdx(root)
                    if atom.HasProp("_CIPCode"):
                        cip_label = atom.GetProp("_CIPCode")

            # Handle fused ring heuristic for unspecified sp3 bridgeheads
            if specified == "Unspecified":
                if type_ == "Atom_Tetrahedral":
                    atom = self.mol.GetAtomWithIdx(root)
                    if (atom.GetAtomicNum() == 6 and atom.IsInRing() and
                        atom.GetHybridization() == Chem.rdchem.HybridizationType.SP3):
                        ring_bond_count = sum(1 for b in atom.GetBonds() if b.IsInRing())
                        if ring_bond_count == 3:
                            descriptor = "Fused ring"

            centre_tuple = (
                self.smiles,
                root_mapped,
                mapped_controlling,
                specified,
                descriptor,
                type_,
                ("CIP", cip_label)
            )
            self.stereocentres.add(centre_tuple)