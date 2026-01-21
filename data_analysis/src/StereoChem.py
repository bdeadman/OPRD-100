# This script parses the reaction SMILES of a given reaction, maps it with rxnmapper, finds the sterocentres in the reactants, and then finds any new stereochemistry generated in the products.
# We check for the presence of E/Z isomerism, R/S isomerism on C, and any S, P or N isomerism.


#from rxnmapper import RXNMapper
import json
import re
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdChemReactions
from rdkit.Chem import rdChemReactions
from rdkit.Chem import rdCIPLabeler 
from rxnmapper import RXNMapper
from collections import Counter

class ReactionParser:

    def __init__(self, json_file_path):
        self.json_file_path = json_file_path
        self.data = self.load_json(json_file_path)
        # Sanity check of molecule count
        self.molecule_count = 0
        # initialise empty arrays of length 8 to store the stereocentre counts
        self.stereo_union = np.zeros(8)
        self.stereo_produces = np.zeros(8)
        # initialise an array of column headings
        self.column_headings = ["E", "Z", "R", "S", "Sulfoxide", "Phosphate", "Phosphine", "Fused ring"]
        # self.apply_functions()
        # # Create a dataframe from the stereo_union and stereo_produces arrays
        # self.df = pd.DataFrame(np.vstack((self.stereo_union, self.stereo_produces)), columns=self.column_headings)
        self.cip_change_counter = Counter()
        self.reaction_cip_changes = []      # list of (reaction_index, [(root, class_label), ...])
        self.reaction_diffs = []            # optional full diff objects per reaction
        self.cip_change_counts_df = pd.DataFrame(columns=["class_label", "count"])


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
            
    def simplify_centres_tuple(self, centres):
        """
        remove SMILES from the stereocentre tuples for easier comparison
        """
        simplified = set()
        for c in centres:
            try:
                # Expected 7-tuple with CIP pair last
                if (isinstance(c, (tuple, list)) and len(c) == 7 and
                        isinstance(c[-1], (tuple, list)) and len(c[-1]) == 2 and c[-1][0] == "CIP"):
                    smiles, root_mapped, mapped_controlling, specified, descriptor, type_, cip_pair = c
                    simplified.add((root_mapped, (mapped_controlling, specified, descriptor, type_, cip_pair[1])))
                else:
                    # Fallback: try to interpret 6‑tuple (no CIP pair)
                    if isinstance(c, (tuple, list)) and len(c) == 6:
                        smiles, root_mapped, mapped_controlling, specified, descriptor, type_ = c
                        simplified.add((root_mapped, (mapped_controlling, specified, descriptor, type_, None)))
            except Exception:
                continue
        return simplified
    
    def diff_simplified_stereocentres(self, reactant_centres, product_centres):
        """
        Compare two sets produced by simplify_centres_tuple(), ignoring controlling atom differences.
        
        Each element is:
            (root, (controlling, specified, descriptor, type_, cip_label))
        
        Returns list of dicts with:
            root, status (changed|unchanged|gained|lost), changes (excluding controlling),
            old, new (each with controlling, specified, descriptor, type, CIP)
        """
        def to_index(ed):
            idx = {}
            for root, features in ed:
                idx.setdefault(root, []).append(features)
            return idx

        # Helper to strip controlling for comparison
        def feature_key(f):
            # f: (controlling, specified, descriptor, type, CIP)
            return (f[1], f[2], f[3], f[4])  # exclude controlling

        r_idx = to_index(reactant_centres)
        p_idx = to_index(product_centres)

        all_roots = set(r_idx.keys()) | set(p_idx.keys())
        results = []

        compare_fields = ["specified", "descriptor", "type", "CIP"]

        for root in sorted(all_roots, key=lambda x: (isinstance(x, tuple), x)):
            r_feats = r_idx.get(root, [])
            p_feats = p_idx.get(root, [])

            if not r_feats and p_feats:
                for nf in p_feats:
                    results.append({
                        "root": root,
                        "status": "gained",
                        "changes": {},
                        "old": None,
                        "new": {
                            "controlling": nf[0],
                            "specified": nf[1],
                            "descriptor": nf[2],
                            "type": nf[3],
                            "CIP": nf[4]
                        }
                    })
                continue

            if r_feats and not p_feats:
                for of in r_feats:
                    results.append({
                        "root": root,
                        "status": "lost",
                        "changes": {},
                        "old": {
                            "controlling": of[0],
                            "specified": of[1],
                            "descriptor": of[2],
                            "type": of[3],
                            "CIP": of[4]
                        },
                        "new": None
                    })
                continue

            p_remaining = p_feats.copy()

            # First match identical (ignoring controlling atoms) to mark unchanged
            matched_p = []
            for of in r_feats:
                key_of = feature_key(of)
                match_idx = next((i for i, pf in enumerate(p_remaining) if feature_key(pf) == key_of), None)
                if match_idx is not None:
                    pf = p_remaining.pop(match_idx)
                    results.append({
                        "root": root,
                        "status": "unchanged",
                        "changes": {},
                        "old": {
                            "controlling": of[0],
                            "specified": of[1],
                            "descriptor": of[2],
                            "type": of[3],
                            "CIP": of[4]
                        },
                        "new": {
                            "controlling": pf[0],
                            "specified": pf[1],
                            "descriptor": pf[2],
                            "type": pf[3],
                            "CIP": pf[4]
                        }
                    })
                    matched_p.append(of)

            r_unmatched = [f for f in r_feats if f not in matched_p]

            # Pair remaining reactant features to remaining product features heuristically
            for of in r_unmatched:
                if not p_remaining:
                    # Orphan => lost (should not happen here because earlier branch covers pure losses)
                    results.append({
                        "root": root,
                        "status": "lost",
                        "changes": {},
                        "old": {
                            "controlling": of[0],
                            "specified": of[1],
                            "descriptor": of[2],
                            "type": of[3],
                            "CIP": of[4]
                        },
                        "new": None
                    })
                    continue
                nf = p_remaining.pop(0)
                changes = {}
                for i, field in enumerate(compare_fields, start=1):
                    old_val = of[i]
                    new_val = nf[i]
                    if old_val != new_val:
                        changes[field] = f"{old_val}->{new_val}"
                results.append({
                    "root": root,
                    "status": "changed" if changes else "unchanged",
                    "changes": changes,
                    "old": {
                        "controlling": of[0],
                        "specified": of[1],
                        "descriptor": of[2],
                        "type": of[3],
                        "CIP": of[4]
                    },
                    "new": {
                        "controlling": nf[0],
                        "specified": nf[1],
                        "descriptor": nf[2],
                        "type": nf[3],
                        "CIP": nf[4]
                    }
                })

            # Any product features still left are gains
            for nf in p_remaining:
                results.append({
                    "root": root,
                    "status": "gained",
                    "changes": {},
                    "old": None,
                    "new": {
                        "controlling": nf[0],
                        "specified": nf[1],
                        "descriptor": nf[2],
                        "type": nf[3],
                        "CIP": nf[4]
                    }
                })

        return results

    def classify_cip_changes(self, diffs):
        """
        From diff_simplified_stereocentres output, extract concise CIP change classes.
        Returns list of tuples: (root, class_string)
        class_string examples:
            R->S, S->R, None->R, R->None, E->Z, Z->E, unchanged, gained(R), lost(R)
        """
        classes = []
        for d in diffs:
            root = d["root"]
            status = d["status"]
            old_cip = d["old"]["CIP"] if d["old"] else None
            new_cip = d["new"]["CIP"] if d["new"] else None

            if status == "unchanged":
                classes.append((root, "unchanged"))
            elif status == "gained":
                classes.append((root, f"gained({new_cip})"))
            elif status == "lost":
                classes.append((root, f"lost({old_cip})"))
            else:
                if old_cip != new_cip:
                    classes.append((root, f"{old_cip}->{new_cip}"))
                else:
                    classes.append((root, "changed(non-CIP)"))
        return classes

    def apply_functions(self):
            """
            Process all reactions:
              - Map reaction
              - Collect stereocentres (reactants/products)
              - Count overall & produced stereochemistry
              - Compute detailed diffs & CIP change classes
              - Aggregate CIP change class counts into a dataframe
            """
            for rxn_idx, reaction_data in enumerate(self.data):
                try:
                    smiles = self.get_mapped_rxn_smiles(reaction_data)
                    rxn = self.get_reaction_obj(smiles)
                except Exception as e:
                    # Skip problematic reaction
                    continue

                reactants, products = self.split_reaction(rxn)

                reactant_centres = set()
                product_centres = set()

                for reactant in reactants:
                    self.molecule_count += 1
                    sc_reactants = CheckStereochemistry(reactant)
                    reactant_centres.update(sc_reactants.stereocentres)

                for product in products:
                    self.molecule_count += 1
                    sc_products = CheckStereochemistry(product)
                    product_centres.update(sc_products.stereocentres)

                # Union counts (all centres present in either side)
                all_centres = reactant_centres.union(product_centres)
                union_counts = self.count_stereochemistry(all_centres)
                self.stereo_union += union_counts

                # Produced = centres appearing only in products
                produced_centres = product_centres - reactant_centres
                produced_counts = self.count_stereochemistry(produced_centres)
                self.stereo_produces += produced_counts

                # Detailed diff & CIP change classes
                react_simple = self.simplify_centres_tuple(reactant_centres)
                prod_simple = self.simplify_centres_tuple(product_centres)
                diffs = self.diff_simplified_stereocentres(react_simple, prod_simple)
                cip_classes = self.classify_cip_changes(diffs)

                # Store per-reaction info
                self.reaction_cip_changes.append((rxn_idx, cip_classes))
                self.reaction_diffs.append((rxn_idx, diffs))

                # Update aggregate class counts
                for _, cls in cip_classes:
                    self.cip_change_counter[cls] += 1

            # Build / refresh summary dataframe
            if self.cip_change_counter:
                self.cip_change_counts_df = pd.DataFrame(
                    sorted(self.cip_change_counter.items(), key=lambda x: (-x[1], x[0])),
                    columns=["class_label", "count"]
                )
            else:
                self.cip_change_counts_df = pd.DataFrame(columns=["class_label", "count"])

            return self.cip_change_counts_df

    # def apply_functions(self):
    #         for reaction_data in self.data:
    #             # Get the reaction SMILES
    #             smiles = self.get_mapped_rxn_smiles(reaction_data)
    #             # Parse the reaction SMILES and create a reaction object
    #             rxn = self.get_reaction_obj(smiles)
    #             # Split the reaction into reactants and products
    #             reactants, products = self.split_reaction(rxn)
    #             # Check for stereochemistry in the reactants and products
    #             # iterate through reatcants
    #             reactant_centres = set()
    #             product_centres = set()
    #             for reactant in reactants:
    #                 self.molecule_count += 1
    #                 # Check for stereochemistry in the reactant
    #                 sc_reactants = CheckStereochemistry(reactant)
    #                 reactant_centres.update(sc_reactants.stereocentres)
    #             # Check for stereochemistry in the product
    #             for product in products:
    #                 self.molecule_count += 1
    #                 # Check for stereochemistry in the product
    #                 sc_products = CheckStereochemistry(product)
    #                 product_centres.update(sc_products.stereocentres)
    #             # get the union of the reactant and product centres
    #             all_centres = reactant_centres.union(product_centres)
    #             # Count the stereochemistry in the reaction
    #             counts = self.count_stereochemistry(all_centres)
    #             # Update the counts
    #             self.stereo_union += counts
    #             # Get the stereocentres produced in the reaction
    #             produced_centres = product_centres - reactant_centres
    #             # Count the stereochemistry in the products
    #             counts = self.count_stereochemistry(produced_centres)
    #             # Update the counts
    #             self.stereo_produces += counts
    
    # def count_stereochemistry(self, centres):
    #         """
    #         Count the number of stereocentres in the reaction.
    #         """
    #         count_E = 0
    #         count_Z = 0
    #         count_R = 0
    #         count_S = 0
    #         count_sulfoxide = 0
    #         count_phosphate = 0
    #         count_phosphine = 0
    #         count_FR = 0
    #         # Iterate through the centres and count the stereochemistry
    #         for centre in centres:
    #             atom_number, controlling_atoms, specified, descriptor, type = centre
    #             if descriptor == "Bond_Trans":
    #                 count_Z += 1
    #             elif descriptor == "Bond_Cis":
    #                 count_E += 1
    #             if descriptor == "Tet_CCW":
    #                 count_S += 1
    #             elif descriptor == "Tet_CW":
    #                 count_R += 1
    #             elif descriptor == "Fused ring":
    #                 count_FR += 1
    #             elif type == "Sulfoxide":
    #                 count_sulfoxide += 1
    #             elif type == "Phosphate":
    #                 count_phosphate += 1
    #             elif type == "Phosphine":
    #                 count_phosphine += 1
    #         # Return the counts as an array
    #         return np.array([count_E, count_Z, count_R, count_S, count_sulfoxide, count_phosphate, count_phosphine, count_FR])

    def count_stereochemistry(self, centres):
            """
            Count stereocentres using CIP labels where available.
            Order: [E, Z, R, S, Sulfoxide, Phosphate, Phosphine, Fused ring]
            Falls back to legacy descriptor names if CIP not assigned.
            """
            count_E = count_Z = count_R = count_S = 0
            count_sulfoxide = count_phosphate = count_phosphine = 0
            count_FR = 0

            for centre in centres:
                # Centre tuple may have shapes:
                # (smiles, atom_number, controlling_atoms, specified, descriptor, type, ('CIP', cip_label))
                # or legacy shorter variants.
                cip_label = None
                descriptor = None
                type_name = None

                # Try to detect CIP pair
                if isinstance(centre, (tuple, list)) and len(centre) >= 6:
                    # With smiles + CIP
                    if len(centre) >= 7 and isinstance(centre[-1], (tuple, list)) and len(centre[-1]) == 2 and centre[-1][0] == "CIP":
                        cip_label = centre[-1][1]
                        descriptor = centre[4]
                        type_name = centre[5]
                    else:
                        # Possibly without CIP pair
                        descriptor = centre[4]
                        type_name = centre[5] if len(centre) > 5 else None
                else:
                    # Unexpected shape; skip
                    continue

                # Normalise CIP label (RDKit gives 'R','S','E','Z')
                if cip_label:
                    if cip_label == "R":
                        count_R += 1
                    elif cip_label == "S":
                        count_S += 1
                    elif cip_label == "E":
                        count_E += 1
                    elif cip_label == "Z":
                        count_Z += 1
                else:
                    # Fallback to legacy descriptor strings if no CIP
                    if descriptor == "Bond_Trans":
                        # Note: Bond_Trans corresponds to E
                        count_E += 1
                    elif descriptor == "Bond_Cis":
                        # Bond_Cis corresponds to Z
                        count_Z += 1
                    elif descriptor == "Tet_CW":
                        # CW generally maps to R
                        count_R += 1
                    elif descriptor == "Tet_CCW":
                        # CCW generally maps to S
                        count_S += 1
                    elif descriptor == "Fused ring":
                        count_FR += 1

                # Hetero stereochemistry types
                if type_name == "Sulfoxide":
                    count_sulfoxide += 1
                elif type_name == "Phosphate":
                    count_phosphate += 1
                elif type_name == "Phosphine":
                    count_phosphine += 1
                elif descriptor == "Fused ring":
                    # Already counted above, but ensure captured if descriptor path only
                    count_FR += 0  # no double count

            return np.array([
                count_E,
                count_Z,
                count_R,
                count_S,
                count_sulfoxide,
                count_phosphate,
                count_phosphine,
                count_FR
            ])

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
        # if use_map_numbers:
        #     self.substitute_atom_map_numbers()


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


    def check_substructure(self, mol, substructure_smarts, decomposition_rxn_smarts, num_diff_groups=2):
        """
        Check for a substructure with chiral centre, and if found, run the decomposition reaction
        to check if the atom is chiral. If it is, return the atom index.
        """
        substructure_match = mol.GetSubstructMatches(Chem.MolFromSmarts(substructure_smarts))
        r_groups = set()
        matching_atoms = []
        controlling_atoms = []
        if substructure_match:
            for match in substructure_match:
                ca = []
                # Create the reaction object from the decomposition SMARTS
                rxn = rdChemReactions.ReactionFromSmarts(decomposition_rxn_smarts)
                products = rxn.RunReactants((mol,))
                print(Chem.MolToSmiles(mol))
                # Process the products
                for product in products:
                    product = product[0]  # Get the first product
                    # Restore atom map numbers
                    for atom in product.GetAtoms():
                        if atom.HasProp('react_atom_idx'):
                            orig_atom = mol.GetAtomWithIdx(atom.GetIntProp('react_atom_idx'))
                            atom.SetAtomMapNum(orig_atom.GetAtomMapNum())
                    # Get the controlling atoms - atoms surrounding the centre
                    controlling_atom = product.GetAtoms()[1] # Skip the H dummy atom
                    controlling_atom_map_num = controlling_atom.GetAtomMapNum()
                    # get the atom id of the controlling atom in the original molecule
                    atom_map_numbers = [atom.GetAtomMapNum() for atom in mol.GetAtoms()]
                    controlling_atom_idx = atom_map_numbers.index(controlling_atom_map_num)
                    controlling_atoms.append(controlling_atom_idx)
                    
                    # controlling_atoms.append(controlling_atom.GetAtomMapNum())
                    # print(f"Controlling atom: {controlling_atom.GetAtomMapNum()}")
                    # print(Chem.MolToSmiles(product))
                    # Remove atom map numbers
                    for atom in product.GetAtoms():
                        atom.SetAtomMapNum(0)
                    # Get the InChI key
                    inchi_key = Chem.MolToInchiKey(product, )
                    r_groups.add(inchi_key)
                    # print(Chem.MolToSmiles(product))
                #ca = tuple(ca)
                # If multiple unique R groups are found, add the atom index
                if len(r_groups) == num_diff_groups:
                    matching_atoms.append(match[0])
                    # print(ca)
                    # controlling_atoms.append(ca)

        return matching_atoms, tuple(controlling_atoms)
    
    # def get_rdkit_stereocentres(self):
    #     """
    #     Return a set of tuples contianing stereocentre information.
    #     uses rdkit atom indices which need to be corrected to the atom map numbers later.
    #     """
    #     si = Chem.FindPotentialStereo(self.mol)
    #     for element in si:
    #         specified = str(element.specified)
    #         descriptor = str(element.descriptor)
    #         atom_number = element.centeredOn
    #         controlling_atoms = np.array(element.controllingAtoms)
    #         type = str(element.type)
    #         cip_label = None  # Default to None
            
    #         if type == "Bond_Double":
    #             centred_on = self.mol.GetBonds()[element.centeredOn]
    #             controlling_atoms = controlling_atoms[controlling_atoms < 1000]
    #             # Get the atoms that are involved in the double bond
    #             if centred_on.GetBondTypeAsDouble() == 2.0:
    #                 atom1 = centred_on.GetBeginAtomIdx()
    #                 atom2 = centred_on.GetEndAtomIdx()
    #                 atom_number = (atom1, atom2)
    #                 # CIP for double bond (E/Z) if specified
    #                 if specified == "Specified":
    #                     db = self.mol.GetBondBetweenAtoms(atom1, atom2)
    #                     if db and db.HasProp("_CIPCode"):
    #                         cip_label = db.GetProp("_CIPCode")

    #         if specified == "Unspecified":
    #             # Check for fused ring unspecified stereogenic centres
    #             atom = self.mol.GetAtoms()[element.centeredOn]
    #             if atom.GetAtomicNum() == 6 and atom.IsInRing() and atom.GetHybridization() == Chem.rdchem.HybridizationType.SP3:
    #                 in_ring_count = 0
    #                 for bond in atom.GetBonds():
    #                     if bond.IsInRing():
    #                         in_ring_count += 1
    #                 if in_ring_count == 3:
    #                     descriptor = "Fused ring"
    #         else:
    #             # CIP for tetrahedral (R/S) if specified
    #             if type == "Atom_Tetrahedral":
    #                 a = self.mol.GetAtomWithIdx(element.centeredOn)
    #                 if a.HasProp("_CIPCode"):
    #                     cip_label = a.GetProp("_CIPCode")

    #         controlling_atoms = tuple(map(int, controlling_atoms))
    #         centre_tuple = (self.smiles, atom_number, controlling_atoms, specified, descriptor, type, ("CIP", cip_label))
    #         self.stereocentres.add(centre_tuple)

    def get_sulphoxide_centres(self):
        """
        Check the molecule for sulfoxide centres, update the stereocentres set.
        """
        sulfoxide_SMARTS = "[S;$([S;D3]([*])([*])=O)]"
        sulfoxide_decomp = "[*;$([*]S([*])=O):1]-[S;$(S([*])([*])=O)]>>[#1]-[*:1]"
        sulfoxide_centres, ca = self.check_substructure(self.mol, sulfoxide_SMARTS, sulfoxide_decomp, num_diff_groups=2)
        if sulfoxide_centres:
            for atom_number in sulfoxide_centres:
                specified = "Unspecified"
                descriptor = "Sulfoxide"
                type = "Sulfoxide"
                centre_tuple = (self.smiles, atom_number, ca, specified, descriptor, type, ("CIP", None))  # add CIP=None
                if centre_tuple not in self.stereocentres:
                    # get all atom numbers in sterocentres
                    atom_numbers = [x[0] for x in self.stereocentres]
                    if atom_number in atom_numbers:
                        self.stereocentres = set([x for x in self.stereocentres if x[0] != atom_number])
                    self.stereocentres.add(centre_tuple)

    def get_phosphine_centres(self):
        phosphine_SMARTS = "[P;$([P;D3]([*])([*])[*])]"
        phosphine_decomp = "[*;$([*][P;D3]([*])[*]):1]-[P;$([P;D3]([*])([*])[*])]>>[#1]-[*:1]"
        p_chiral_centres, ca = self.check_substructure(self.mol, phosphine_SMARTS, phosphine_decomp, num_diff_groups=3)
        if p_chiral_centres:
            for atom_number in p_chiral_centres:
                specified = "Unspecified"
                descriptor = "Phosphine"
                type = "Phosphine"
                centre_tuple = (self.smiles, atom_number, ca, specified, descriptor, type, ("CIP", None))  # add CIP=None
                if centre_tuple not in self.stereocentres:
                    atom_numbers = [x[0] for x in self.stereocentres]
                    if atom_number in atom_numbers:
                        self.stereocentres = set([x for x in self.stereocentres if x[0] != atom_number])
                    self.stereocentres.add(centre_tuple)
    
    def get_phosphate_centres(self):
        phosphate_SMARTS = "[P;$([P;D4]([O][*])([O][*])([O][*])=[O])]"
        phosphate_decomp = "[*;$([*][P;D4]([*])([*])=[O]):1]-[P;$([P;D4]([*])([*])([*])=[O])]>>[#1]-[*:1]"
        p_chiral_centres, ca = self.check_substructure(self.mol, phosphate_SMARTS, phosphate_decomp, num_diff_groups=3)
        if p_chiral_centres:
            for atom_number in p_chiral_centres:
                specified = "Unspecified"
                descriptor = "Phosphate"
                type = "Phosphate"
                centre_tuple = (self.smiles, atom_number, ca, specified, descriptor, type, ("CIP", None))  # add CIP=None
                if centre_tuple not in self.stereocentres:
                    atom_numbers = [x[0] for x in self.stereocentres]
                    if atom_number in atom_numbers:
                        self.stereocentres = set([x for x in self.stereocentres if x[0] != atom_number])
                    self.stereocentres.add(centre_tuple)
       
    
    def substitute_atom_map_numbers(self):
        """
        Substitute the atom map numbers in the stereocentres set with the atom indices in the mapped molecule.
        """
        updated_stereocentres = set()
        for item in self.stereocentres:
            # Split out optional CIP pair
            cip_pair = ("CIP", None)
            core = item
            if isinstance(item, (tuple, list)) and len(item) >= 7 and isinstance(item[-1], (tuple, list)) and len(item[-1]) == 2 and item[-1][0] == "CIP":
                cip_pair = tuple(item[-1])
                core = item[:-1]

            # Unpack with/without SMILES
            if len(core) == 6:
                smiles, atom_number, ca, specified, descriptor, type = core
            elif len(core) == 5:
                atom_number, ca, specified, descriptor, type = core
                smiles = self.smiles
            else:
                # Unexpected shape; keep as-is
                updated_stereocentres.add(tuple(item))
                continue

            # update the controlling atoms to the atom map numbers
            controlling_atoms = []
            for atom in ca:
                atom_obj = self.mol.GetAtomWithIdx(int(atom))
                atom_map_number = atom_obj.GetAtomMapNum()
                controlling_atoms.append(atom_map_number if atom_map_number != 0 else -1)
            controlling_atoms = tuple(controlling_atoms)

            if isinstance(atom_number, tuple):
                atom1 = self.mol.GetAtomWithIdx(atom_number[0])
                atom2 = self.mol.GetAtomWithIdx(atom_number[1])
                atom_map_number1 = atom1.GetAtomMapNum() or -1
                atom_map_number2 = atom2.GetAtomMapNum() or -1
                updated_stereocentres.add((smiles, (atom_map_number1, atom_map_number2), controlling_atoms, specified, descriptor, type, cip_pair))
            else:
                atom_obj = self.mol.GetAtomWithIdx(atom_number)
                atom_map_number = atom_obj.GetAtomMapNum()
                if atom_map_number != 0:
                    updated_stereocentres.add((smiles, atom_map_number, ca, specified, descriptor, type, cip_pair))
                else:
                    updated_stereocentres.add((smiles, -1, controlling_atoms, specified, descriptor, type, cip_pair))
        self.stereocentres = updated_stereocentres