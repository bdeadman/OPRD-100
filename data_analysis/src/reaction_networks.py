"""
Python file to generate reaction network maps for comparisons between reaction datasets.

This file is constructed to operate on two datasets, which both have data from the same reference reaction set (i.e. a research article).

The code generates reactions networks first by collecting common molecules found in both datasets and storing them as nodes. The reactions between these molecules are then stored as edges in a network graph. The resulting graph can be used to visualize the relationships between the reactions and molecules in both datasets - with edge and node colours representing the dataset they belong to.
"""

import matplotlib.pyplot as plt
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
import cairosvg
import json
import networkx as nx
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
from io import BytesIO
from PIL import Image, ImageChops
import numpy as np
import pandas as pd
from ReaxysComparison import ReaxysDataExtractor
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import rdDepictor
from rdkit import DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem import rdChemReactions as Reactions

class ReactionNetworkGenerator:
    """
    Class to generate reaction networks from two datasets.
    It can load data from JSON, CSV, and XML formats, extract reactions based on DOI references,
    and build a reaction network graph.
    """

    def __init__(
            self, 
            dataset1_path, 
            dataset2_path, 
            doi, 
            dataset1_name="dataset_1", 
            dataset2_name="dataset_2",
            yields_1=None,
            yields_2=None,
            mol_labels=None
    ):
        """
        Initialize with paths to two datasets.
        Args:
            dataset1_path (str): Path to the first dataset file, 
            file must have Reaction SMILES stored under "Reaction" tag or column header.
            dataset2_path (str): Path to the second dataset file, 
            file must have Reaction SMILES stored under "Reaction" tag or column header.
            doi (str): DOI reference to filter reactions.
            dataset1_name (str): Name for the first dataset.
            dataset2_name (str): Name for the second dataset.
            yields_1 (dict, optional): {InChi Key: Yield} for reactions in the first dataset.
            yields_2 (dict, optional): {InChi Key: Yield} for reactions in the second dataset.
            mol_labels (dict, optional): {InChi Key: Label} for molecules in the datasets.

        Attributes:
            dataset1_path (str): Path to the first dataset.
            dataset2_path (str): Path to the second dataset.
            dataset1_name (str): Name of the first dataset.
            dataset2_name (str): Name of the second dataset.
            doi (str): DOI reference for filtering reactions.
            reactions_1 (list): Extracted reactions from the first dataset based on DOI.
            reactions_2 (list): Extracted reactions from the second dataset based on DOI.

        """
        self.dataset1_path = dataset1_path
        self.dataset2_path = dataset2_path
        self.dataset1_name = dataset1_name
        self.dataset2_name = dataset2_name
        self.doi = doi
        self.yields_1 = {}
        self.yields_2 = {}
        self.mol_labels = self.smiles_dict_to_inchi_dict(mol_labels) if mol_labels else {}
        # Load datasets
        self.load_data()
        # Extract reactions based on DOI reference
        self.reactions_1 = self.load_reaction_data(self.dataset1, self.doi)
        self.reactions_2 = self.load_reaction_data(self.dataset2, self.doi)
        self.reaction_keys_1 = {self.rxn_key(r): r for r in self.reactions_1}
        self.reaction_keys_2 = {self.rxn_key(r): r for r in self.reactions_2}
        self.num_invalid_reactions = 0
        # Generate node and edge lists for the reaction network
        self.nodes, self.edges = self.build_directed_graph()
        # Build the directed graph from nodes and edges
        self.G = None  # Will be set later when building the graph
        

        
    def smiles_dict_to_inchi_dict(self, smiles_dict):
        """
        Convert a dictionary of SMILES strings to InChI keys.
        """
        inchi_dict = {}
        for smiles, label in smiles_dict.items():
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                inchi = Chem.MolToInchiKey(mol)
                inchi_dict[inchi] = label
        return inchi_dict

    def load_data(self):
        """
        Load data from the specified dataset paths.
        """
        self.dataset1 = self.load_dataset(self.dataset1_path)
        self.dataset2 = self.load_dataset(self.dataset2_path)

    def load_dataset(self, file_path):
        """
        Load a dataset from a file.
        """
        if file_path.endswith(".json"):
            return self.load_json_data(file_path)
        elif file_path.endswith(".csv"):
            return self.load_csv_data(file_path)
        elif file_path.endswith(".xlsx"):
            return self.load_excel_data(file_path)
        elif file_path.endswith(".xml"):
            return file_path  # Return the file path for XML, as it will be processed later
        else:
            raise ValueError(f"Error on file import {file_path}.")
        
    def load_json_data(self, file_path):
        """Load JSON data from a file."""
        with open(file_path, 'r') as file:
            data = json.load(file)
        return data

    def load_csv_data(self, file_path):
        """Load CSV data from a file."""
        df = pd.read_csv(file_path)
        return df

    def load_excel_data(self, file_path):
        """Load Excel data from a file."""
        df = pd.read_excel(file_path)
        return df

    def load_reaction_data(self, dataset, doi=None):
        """
        Extract reactions from the dataset based on DOI reference.
        """
        if isinstance(dataset, list):
            return self.extract_reactions_from_json(dataset, doi)
        elif isinstance(dataset, pd.DataFrame):
            return self.extract_reactions_from_csv(dataset, doi)
        elif isinstance(dataset, str) and dataset.endswith(".xml"):
            return self.extract_reactions_from_xml(dataset, doi)

    def extract_reactions_from_json(self, data, doi=None):
        """Extract reactions from JSON data based on DOI reference."""
        reactions = []
        for entry in data:
            if doi and entry.get("Reference") == doi:
                reaction = entry.get("Reaction")
                if reaction:
                    reactions.append(reaction)
        return reactions

    def extract_reactions_from_csv(self, df, doi=None):
        """Extract reactions from a DataFrame based on DOI reference."""
        reactions = []
        has_yield = []
        for index, row in df.iterrows():
            if doi and row.get("Reference") == doi:
                reaction = row.get("Reaction")
                if reaction:
                    reactions.append(reaction)
                # rxn_data = row.get("Reaction data")
                # if rxn_data:
                #     # Split into steps
                #     steps = re.findall(r"\d+: ([^;]+)", rxn_data)
                #     # Get the last step and get yield
                #     yield_data = steps[-1].split(" / ")[-1]
                #     if yield_data != "None":
                #         has_yield.append(True)
                #     else:
                #         has_yield.append(False)
        return reactions

    def extract_reactions_from_xml(self, file_path, doi=None):
        """Extract reactions from an XML file based on DOI reference."""
        rde = ReaxysDataExtractor(file_path)
        df = pd.DataFrame(rde.extract_reaction_data()) 
        if doi:
            df = df[df['Reference'] == doi]
        reactions = df["Reaction"].tolist()
        
        return reactions

    def compute_unique_molecule_tanimoto(
        self,
        dataset_a=None,                # defaults to self.dataset1_name
        dataset_b=None,                # defaults to self.dataset2_name
        radius=2,
        nBits=2048,
        threshold=0.8,
        remove_stereo=False,
        return_matrix=True,
        return_pairs=False,
    ):
        """
        Compute the Tanimoto similarity between molecules unique to dataset A vs dataset B.

        Parameters
        - dataset_a: str or None. Dataset name; defaults to self.dataset1_name.
        - dataset_b: str or None. Dataset name; defaults to self.dataset2_name.
        - radius: int. Morgan fingerprint radius.
        - nBits: int. Bit vector length.
        - threshold: float in [0,1]. Count pairs with similarity > threshold.
        - remove_stereo: bool. If True, stereochemistry is removed before fingerprinting.
        - return_matrix: bool. If True, returns a pandas.DataFrame similarity matrix (A x B).
        - return_pairs: bool. If True, also returns a list of (smiles_a, smiles_b, similarity) above threshold.

        Returns
        - result: dict with keys:
            'count': int, number of cross-dataset pairs with similarity > threshold
            'matrix': pandas.DataFrame or None, similarity matrix (rows=unique A, cols=unique B) if return_matrix=True
            'pairs': list or None, list of tuples (smi_a, smi_b, sim) if return_pairs=True
            'rows': list of SMILES (row labels)
            'cols': list of SMILES (column labels)
        Notes
        - Uses self.get_common_and_unique_nodes(), which requires self.G to be built first.
        """
        if getattr(self, "G", None) is None:
            raise RuntimeError("Graph not built. Call build_graph_from_nodes_edges() before computing similarities.")

        dataset_a = dataset_a or self.dataset1_name
        dataset_b = dataset_b or self.dataset2_name

        node_dict = self.get_common_and_unique_nodes()
        key_a = f"unique_nodes_{dataset_a}"
        key_b = f"unique_nodes_{dataset_b}"
        if key_a not in node_dict or key_b not in node_dict:
            raise ValueError(f"Unique node keys not found. Got keys: {list(node_dict.keys())}")

        smiles_a = sorted({s for s in node_dict[key_a] if isinstance(s, str) and s})
        smiles_b = sorted({s for s in node_dict[key_b] if isinstance(s, str) and s})

        # Short-circuit if empty
        if not smiles_a or not smiles_b:
            import pandas as pd
            mat = pd.DataFrame(index=smiles_a, columns=smiles_b, dtype=float) if return_matrix else None
            return {"count": 0, "matrix": mat, "pairs": [] if return_pairs else None, "rows": smiles_a, "cols": smiles_b}

        # Build RDKit molecules (optionally remove stereo)
        from rdkit import Chem
        from rdkit.Chem import AllChem
        from rdkit import DataStructs
        import numpy as np
        import pandas as pd

        def to_mol(smi):
            m = Chem.MolFromSmiles(smi)
            if not m:
                return None
            if remove_stereo:
                Chem.RemoveStereochemistry(m)
            return m

        mols_a = [to_mol(s) for s in smiles_a]
        mols_b = [to_mol(s) for s in smiles_b]

        # Filter out None entries consistently in labels and mols
        keep_a = [(s, m) for s, m in zip(smiles_a, mols_a) if m is not None]
        keep_b = [(s, m) for s, m in zip(smiles_b, mols_b) if m is not None]
        if not keep_a or not keep_b:
            mat = pd.DataFrame(index=[s for s, _ in keep_a], columns=[s for s, _ in keep_b], dtype=float) if return_matrix else None
            return {"count": 0, "matrix": mat, "pairs": [] if return_pairs else None,
                    "rows": [s for s, _ in keep_a], "cols": [s for s, _ in keep_b]}

        smiles_a, mols_a = zip(*keep_a)
        smiles_b, mols_b = zip(*keep_b)

        # Fingerprints
        fps_a = [AllChem.GetMorganFingerprintAsBitVect(m, radius, nBits=nBits) for m in mols_a]
        fps_b = [AllChem.GetMorganFingerprintAsBitVect(m, radius, nBits=nBits) for m in mols_b]

        # Similarity matrix (A x B); use BulkTanimotoSimilarity for speed
        nA, nB = len(fps_a), len(fps_b)
        sims = np.zeros((nA, nB), dtype=float)
        above_pairs = [] if return_pairs else None
        count = 0
        for i, fpa in enumerate(fps_a):
            row = DataStructs.BulkTanimotoSimilarity(fpa, fps_b)
            sims[i, :] = row
            # Count pairs strictly greater than threshold
            cnt_i = sum(1 for val in row if val > threshold)
            count += cnt_i
            if return_pairs:
                for j, val in enumerate(row):
                    if val > threshold:
                        above_pairs.append((smiles_a[i], smiles_b[j], float(val)))

        mat = None
        if return_matrix:
            mat = pd.DataFrame(sims, index=list(smiles_a), columns=list(smiles_b))

        return {
            "count": int(count),
            "matrix": mat,
            "pairs": above_pairs,
            "rows": list(smiles_a),
            "cols": list(smiles_b),
        }



    def canonical_smiles(self, smi: str) -> str | None:
        """Return RDKit-canonical SMILES (or None)."""
        m = Chem.MolFromSmiles(smi)
        return Chem.MolToSmiles(m, canonical=True) if m else None


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

    def build_directed_graph(self):
        """
        Build a directed graph from the reactions in both datasets.
        Nodes are InChIKeys (id) with canonical SMILES stored as 'smiles'.
        Edges are built using reaction keys derived from InChIKeys.
        """

        def collect_ik_to_smiles_map(reactions: list[str]) -> dict[str, str]:
            """Collect IK -> canonical SMILES from a list of reaction SMILES."""
            ik_to_smi: dict[str, str] = {}
            for rxn in reactions:
                try:
                    rsmi, psmi = self.parse_reaction_smiles(rxn)
                except Exception:
                    continue
                for s in rsmi + psmi:
                    ik = self.smiles_to_inchikey(s)
                    if ik and ik not in ik_to_smi:
                        ik_to_smi[ik] = s
            return ik_to_smi

        # Build IK->SMILES maps for each dataset
        ik_map_1 = collect_ik_to_smiles_map(self.reactions_1)
        ik_map_2 = collect_ik_to_smiles_map(self.reactions_2)

        set_ik_1 = set(ik_map_1.keys())
        set_ik_2 = set(ik_map_2.keys())
        all_iks = set_ik_1 | set_ik_2

        # Create molecule nodes keyed by IK with canonical SMILES
        nodes: dict[str, dict] = {}
        for ik in all_iks:
            in_1 = ik in set_ik_1
            in_2 = ik in set_ik_2
            if in_1 and in_2:
                dataset = "both"
                color = "green"
            elif in_1:
                dataset = self.dataset1_name
                color = "blue"
            else:
                dataset = self.dataset2_name
                color = "red"

            # Prefer dataset1 SMILES when present, else dataset2
            smiles = ik_map_1.get(ik, ik_map_2.get(ik))
            nodes[ik] = {
                "id": ik,
                "shape": "o",
                "color": color,
                "type": "molecule",
                "dataset": dataset,
                "smiles": smiles,         # canonical SMILES
                "label": self.mol_labels.get(ik, None),
            }

        # Build reaction maps using IK-based keys
        rxn_map_1 = {self.rxn_key(r): r for r in self.reactions_1 if self.rxn_key(r) is not None}
        rxn_map_2 = {self.rxn_key(r): r for r in self.reactions_2 if self.rxn_key(r) is not None}
        all_rxn_keys = set(rxn_map_1) | set(rxn_map_2)

        edges = []
        for rk in all_rxn_keys:
            in_1 = rk in rxn_map_1
            in_2 = rk in rxn_map_2
            dataset = "both" if in_1 and in_2 else (self.dataset1_name if in_1 else self.dataset2_name)
            edge_color = "green" if dataset == "both" else ("blue" if dataset == self.dataset1_name else "red")

            rxn = rxn_map_1.get(rk, rxn_map_2.get(rk))
            try:
                rsmi, psmi = self.parse_reaction_smiles(rxn)
            except Exception:
                self.num_invalid_reactions += 1
                continue

            reactant_iks = [self.smiles_to_inchikey(s) for s in rsmi]
            product_iks = [self.smiles_to_inchikey(s) for s in psmi]
            reactant_iks = [k for k in reactant_iks if k]
            product_iks = [k for k in product_iks if k]

            intermediate_reactant = None
            intermediate_product = None

            # Intermediate nodes for multiple reactants/products
            if len(reactant_iks) > 1:
                intR_id = f"{rk}_intR"
                if intR_id not in nodes:
                    nodes[intR_id] = {
                        "id": intR_id,
                        "shape": "s",
                        "color": edge_color,
                        "type": "int_reactant",
                        "dataset": dataset,
                    }
                intermediate_reactant = intR_id
                for r in reactant_iks:
                    edges.append({
                        "source": r,
                        "target": intR_id,
                        "color": edge_color,
                        "dataset": dataset,
                        "reaction_key": rk,
                    })

            if len(product_iks) > 1:
                intP_id = f"{rk}_intP"
                if intP_id not in nodes:
                    nodes[intP_id] = {
                        "id": intP_id,
                        "shape": "^",
                        "color": edge_color,
                        "type": "int_product",
                        "dataset": dataset,
                    }
                intermediate_product = intP_id
                for p in product_iks:
                    edges.append({
                        "source": intP_id,
                        "target": p,
                        "color": edge_color,
                        "dataset": dataset,
                        "reaction_key": rk,
                    })

            # Connect accordingly
            if intermediate_reactant and intermediate_product:
                edges.append({
                    "source": intermediate_reactant,
                    "target": intermediate_product,
                    "color": edge_color,
                    "dataset": dataset,
                    "reaction_key": rk,
                })
            elif intermediate_reactant and not intermediate_product and product_iks:
                edges.append({
                    "source": intermediate_reactant,
                    "target": product_iks[0],
                    "color": edge_color,
                    "dataset": dataset,
                    "reaction_key": rk,
                })
            elif not intermediate_reactant and intermediate_product and reactant_iks:
                edges.append({
                    "source": reactant_iks[0],
                    "target": intermediate_product,
                    "color": edge_color,
                    "dataset": dataset,
                    "reaction_key": rk,
                })
            else:
                try:
                    r, p = reactant_iks[0], product_iks[0]
                    edges.append({
                        "source": r,
                        "target": p,
                        "color": edge_color,
                        "dataset": dataset,
                        "reaction_key": rk,
                    })
                except Exception:
                    self.num_invalid_reactions += 1

        return list(nodes.values()), edges


    def smiles_to_inchikey(self, smiles):
        """Convert SMILES to InChIKey, or return None if invalid."""
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return None
        inchikey = Chem.MolToInchiKey(mol)
        return inchikey

    def smiles_to_node_id(self, smiles, nodes):
        """Return node ID (InChIKey) for a given SMILES, or None if not found."""
        can = self.canonical_smiles(smiles)
        target_ik = self.smiles_to_inchikey(can) if can else None
        if not target_ik:
            return None
        for node in nodes:
            if node.get("id") == target_ik:
                return node["id"]
        return None


    def get_unique_molecule_locations(self, dataset="dataset1"):
        """
        Find Location values (from the dataset file) for molecule nodes unique to a given dataset.
        dataset: "dataset1", "dataset2", or the actual dataset name (self.dataset1_name/self.dataset2_name)
        Returns: pandas.DataFrame with columns:
            ['DOI','Dataset','InChIKey','SMILES','Locations','Location_Count']
        Notes:
            - Uses rows for this instance's DOI (Reference == self.doi) when present.
            - Requires columns: Reaction, Location (and optionally Reference) in the chosen dataset.
        """
        # Resolve dataset selection
        if dataset in ("dataset1", self.dataset1_name):
            df = self.dataset1
            dataset_name = self.dataset1_name
        elif dataset in ("dataset2", self.dataset2_name):
            df = self.dataset2
            dataset_name = self.dataset2_name
        else:
            raise ValueError(f"Unknown dataset selector: {dataset}")

        # Must be a table with Reaction/Location
        if not isinstance(df, pd.DataFrame):
            raise ValueError(f"Selected dataset '{dataset_name}' is not a tabular file (CSV/XLSX).")
        if "Reaction" not in df.columns:
            raise ValueError(f"Selected dataset '{dataset_name}' missing 'Reaction' column.")
        if "Location" not in df.columns:
            raise ValueError(f"Selected dataset '{dataset_name}' missing 'Location' column.")

        # Restrict to this DOI if Reference column exists
        if "Reference" in df.columns:
            df_slice = df[df["Reference"] == self.doi].copy()
        else:
            df_slice = df.copy()

        # Determine unique molecule keys for the chosen dataset
        def iter_nodes_from_G_or_cache():
            if getattr(self, "G", None) is not None:
                for n, d in self.G.nodes(data=True):
                    yield n, d
            else:
                for nd in self.nodes:
                    yield nd["id"], nd

        unique_keys = set()
        key_to_meta = {}
        for nid, nd in iter_nodes_from_G_or_cache():
            if nd.get("type") == "molecule" and nd.get("dataset") == dataset_name:
                ikey = nid                      # InChIKey used as node id
                unique_keys.add(ikey)
                key_to_meta[ikey] = {
                    "smiles": nd.get("smiles"),
                }

        # Nothing to do
        if not unique_keys or df_slice.empty:
            return pd.DataFrame(columns=["DOI","Dataset","InChIKey","SMILES","Locations","Location_Count"])

        # Accumulate locations for each unique InChIKey
        from collections import defaultdict
        key_to_locations = defaultdict(set)

        for _, row in df_slice.iterrows():
            rxn = row.get("Reaction")
            if not isinstance(rxn, str) or ">>" not in rxn:
                continue
            loc = row.get("Location", None)
            # Parse reaction molecules
            try:
                reactants, products = rxn.split(">>")
            except Exception:
                continue
            parts = []
            parts.extend([s.strip() for s in reactants.split(".") if isinstance(s, str) and s.strip()])
            parts.extend([s.strip() for s in products.split(".") if isinstance(s, str) and s.strip()])

            # Convert to InChIKey and check uniqueness membership
            for smi in parts:
                ikey = self.smiles_to_inchikey(smi)
                if not ikey:
                    continue
                if ikey in unique_keys and isinstance(loc, str) and loc.strip():
                    key_to_locations[ikey].add(loc.strip())

        # Build result dataframe
        rows = []
        for ikey in sorted(unique_keys):
            meta = key_to_meta.get(ikey, {})
            locs = sorted(key_to_locations.get(ikey, set()))
            rows.append({
                "DOI": self.doi,
                "Dataset": dataset_name,
                "InChIKey": ikey,
                "SMILES": meta.get("smiles"),
                "Locations": ", ".join(locs),
                "Location_Count": len(locs),
            })

        return pd.DataFrame(rows)

    def get_unique_edge_locations(self, dataset="dataset1"):
        """
        Find Location values for reactions (edges) unique to a given dataset.
        dataset: "dataset1", "dataset2", or the actual dataset name (self.dataset1_name/self.dataset2_name)
        Returns: pandas.DataFrame with columns:
            ['DOI','Dataset','ReactionKey','Reaction','Locations','Location_Count']
        Notes:
            - Uses rows for this instance's DOI (Reference == self.doi) when present.
            - Requires columns: Reaction, Location (and optionally Reference) in the chosen dataset.
            - ReactionKey is computed from normalized (InChI-based) reactants/products.
        """
        # Resolve dataset selection
        if dataset in ("dataset1", self.dataset1_name):
            df = self.dataset1
            dataset_name = self.dataset1_name
            keys_this = set(self.reaction_keys_1.keys())
            keys_other = set(self.reaction_keys_2.keys())
            key_to_rxn = self.reaction_keys_1
        elif dataset in ("dataset2", self.dataset2_name):
            df = self.dataset2
            dataset_name = self.dataset2_name
            keys_this = set(self.reaction_keys_2.keys())
            keys_other = set(self.reaction_keys_1.keys())
            key_to_rxn = self.reaction_keys_2
        else:
            raise ValueError(f"Unknown dataset selector: {dataset}")

        # Validate table + columns
        if not isinstance(df, pd.DataFrame):
            raise ValueError(f"Selected dataset '{dataset_name}' is not a tabular file (CSV/XLSX).")
        for col in ("Reaction", "Location"):
            if col not in df.columns:
                raise ValueError(f"Selected dataset '{dataset_name}' missing '{col}' column.")

        # Filter to this DOI if present
        df_slice = df[df["Reference"] == self.doi].copy() if "Reference" in df.columns else df.copy()

        # Unique reaction keys for the chosen dataset
        unique_keys = keys_this - keys_other
        if not unique_keys or df_slice.empty:
            return pd.DataFrame(columns=["DOI","Dataset","ReactionKey","Reaction","Locations","Location_Count"])

        # Accumulate locations per unique reaction key
        from collections import defaultdict
        key_to_locations = defaultdict(set)

        for _, row in df_slice.iterrows():
            rxn = row.get("Reaction")
            if not isinstance(rxn, str) or ">>" not in rxn:
                continue
            k = self.rxn_key(rxn)
            if not k or k not in unique_keys:
                continue
            loc = row.get("Location")
            if isinstance(loc, str) and loc.strip():
                key_to_locations[k].add(loc.strip())

        # Build result dataframe
        out_rows = []
        for k in sorted(unique_keys):
            rxn_str = key_to_rxn.get(k, None)
            locs = sorted(key_to_locations.get(k, set()))
            out_rows.append({
                "DOI": self.doi,
                "Dataset": dataset_name,
                "ReactionKey": k,
                "Reaction": rxn_str,
                "Locations": ", ".join(locs),
                "Location_Count": len(locs),
            })
        return pd.DataFrame(out_rows)


    def get_inchi(self, smi):
        mol = Chem.MolFromSmiles(smi)
        return Chem.MolToInchi(mol) if mol else None

    def get_inchi_key(self, inchi):
        mol = Chem.MolFromInchi(inchi)
        return Chem.MolToInchiKey(mol) if mol else None

    def get_smiles_from_inchi(self, inchi):
        mol = Chem.MolFromInchi(inchi)
        #mol = self.imine_to_amide(mol)
        return Chem.MolToSmiles(mol) if mol else None
    

    def build_graph_from_nodes_edges(
        self,
        from_smiles=None,
        to_smiles=None,
        hide_datasets=None,
        include_connected=True
    ):
        """
        Construct a NetworkX DiGraph from lists of node and edge attribute dicts.
        Optionally restrict to the subgraph between two SMILES.
        Optionally hide edges with a dataset in hide_datasets (a list of dataset names).
        Each node dict must have at least 'id'.
        Each edge dict must have 'source' and 'target'.

        Returns:
            G (networkx.DiGraph): The constructed graph.
        """
        if hide_datasets is None:
            hide_datasets = []

        G = nx.DiGraph()
        # Add nodes with attributes
        for node in self.nodes:
            node_id = node['id']
            dataset = node.get('dataset', None)
            if dataset in hide_datasets:
                continue
            G.add_node(node_id, **{k: v for k, v in node.items() if k != 'id'})

        # Add edges with attributes, skip if dataset is hidden
        for edge in self.edges:
            src = edge['source']
            tgt = edge['target']
            dataset = edge.get('dataset', None)
            if dataset in hide_datasets:
                continue
            attrs = {k: v for k, v in edge.items() if k not in ['source', 'target']}
            G.add_edge(src, tgt, **attrs)

        # If restricting to subgraph between two nodes
        if from_smiles and to_smiles:
            source_id = self.smiles_to_node_id(from_smiles, self.nodes)
            target_id = self.smiles_to_node_id(to_smiles, self.nodes)
            if not source_id:
                raise ValueError(
                    f"Could not find node for the given SMILES: {from_smiles}"
                )
            if not target_id:
                raise ValueError(
                    f"Could not find node for the given SMILES: {to_smiles}"
                )
            else:
                G = self.find_subgraph_between_nodes(G, source_id, target_id, include_connected=include_connected)
                if G.number_of_nodes() == 0:
                    print("No path found between those nodes. No network generated.")
                    return None
        self.G = G


    def find_subgraph_between_nodes(self, G, source_id, target_id, include_connected=True):
        """
        Return a subgraph containing all nodes/edges on any path from source_id to target_id,
        AND also any directly connected molecule ('o'/circle) nodes to triangle ('^') or square ('s') nodes in the subgraph.
        """
        paths = list(nx.all_simple_paths(G, source=source_id, target=target_id))
        nodes_on_paths = set()
        edges_on_paths = set()
        for path in paths:
            nodes_on_paths.update(path)
            edges_on_paths.update(zip(path[:-1], path[1:]))

        # Start with just those nodes/edges
        subG = G.subgraph(nodes_on_paths).copy()
        # Remove edges not in edges_on_paths
        for u, v in list(subG.edges()):
            if (u, v) not in edges_on_paths:
                subG.remove_edge(u, v)

        if include_connected:
            # --- Add immediate molecule neighbors of triangles and squares
            extra_nodes = set()
            extra_edges = set()
            for n in subG.nodes:
                shape = G.nodes[n].get('shape', None)
                if shape in ('^', 's'):
                    # Check all neighbors (predecessors and successors)
                    for neighbor in set(G.predecessors(n)).union(G.successors(n)):
                        if neighbor not in subG.nodes:
                            neighbor_shape = G.nodes[neighbor].get('shape', None)
                            if neighbor_shape == 'o':
                                extra_nodes.add(neighbor)
                                # Add all edges between n and neighbor (in either direction)
                                if G.has_edge(n, neighbor):
                                    extra_edges.add((n, neighbor))
                                if G.has_edge(neighbor, n):
                                    extra_edges.add((neighbor, n))
            # Add the extra nodes and edges
            for n in extra_nodes:
                subG.add_node(n, **G.nodes[n])
            for u, v in extra_edges:
                # Copy attributes from the original edge
                if G.has_edge(u, v):
                    subG.add_edge(u, v, **G.edges[u, v])
        return subG
    
    

    def draw_reaction_network(
        self,
        max_per_layer=10,
        img_dpi=300,
        figsize=(8, 12),

        # Layout
        base_ys={0: 0.8, 1: 0.7, 2: 0.5, 3: 0.3, 4: 0.2},
        sublayer_dy=0.08,
        y_offset_per_node=0.02,

        # Node drawing
        node_size=200,
        node_alpha=0.8,

        # Edge drawing
        edge_alpha=0.3,
        edge_width=2,

        # Molecule images
        img_y_top=0.9,
        img_y_bottom=0.1,
        img_y_offset=1,
        img_connection_offset=0.04,
        img_zoom=1.0,

        # Legend
        legend_fontsize=10,
        legend_bbox_to_anchor=(1.15, 1.0),
    ):
        """
        Draws a reaction network with tuneable parameters for layout, node/edge style, molecule images, and legend.
        """

        # Assign some space in the figure for each layer
        layer_nodes = {i: [] for i in base_ys}
        node_layer = {}
        for n in self.G.nodes:
            layer = self.get_node_layer(self.G, n)
            node_layer[n] = layer
            layer_nodes[layer].append(n)

        # Split layers into sublayers if they have too many nodes.
        pos = {}
        for layer, nodes in layer_nodes.items():
            if not nodes:
                continue
            chunks = list(self.chunk_list(nodes, max_per_layer))
            for ci, chunk in enumerate(chunks):
                base_y = base_ys[layer] - ci * sublayer_dy  # Adjust dy for sublayers
                n_chunk = len(chunk)
                x_coords = [i / (n_chunk - 1) if n_chunk > 1 else 0.5 for i in range(n_chunk)]
                offsets = self.alternating_constant_offsets(n_chunk, y_offset_per_node)
                for xi, n in enumerate(chunk):
                    y = base_y + offsets[xi]
                    pos[n] = (x_coords[xi], y)
        for n in self.G.nodes:
            if n not in pos:
                pos[n] = (0.5, 0.0)

        fig, ax = plt.subplots(figsize=figsize, dpi=img_dpi)
        ax.axis('off')

        # Draw edges grouped by color
        edge_color_map = {}
        for u, v, d in self.G.edges(data=True):
            color = d.get('color', 'black')
            edge_color_map.setdefault(color, []).append((u, v))
        for color, edgelist in edge_color_map.items():
            nx.draw_networkx_edges(
                self.G, pos, ax=ax, edgelist=edgelist,
                edge_color=color, alpha=edge_alpha, arrows=True, width=edge_width
            )

        # Draw nodes by shape
        shape_map = {'o': [], 's': [], '^': []}
        color_map = {'o': [], 's': [], '^': []}
        for n, d in self.G.nodes(data=True):
            shape = d.get('shape', 'o')
            color = d.get('color', 'black')
            shape_map.setdefault(shape, []).append(n)
            color_map.setdefault(shape, []).append(color)
        for shape, nodelist in shape_map.items():
            if nodelist:
                nx.draw_networkx_nodes(
                    self.G, pos, nodelist=nodelist, node_shape=shape,
                    node_color=color_map[shape], node_size=node_size,
                    ax=ax, alpha=node_alpha
                )

        # Draw molecule images at the top and bottom
        # Starting materials (top)
        starting_material_nodes = [n for n in self.G.nodes if node_layer[n] == 0]
        for i, n in enumerate(starting_material_nodes):
            smiles = self.G.nodes[n].get('smiles')
            if not smiles:
                continue
            mol = Chem.MolFromSmiles(smiles)
            img_size = self.get_img_size_for_mol(mol)
            try:
                png_bytes = self.mol_to_png_bytes(smiles, img_size=img_size)
            except Exception as e:
                print(f"Error drawing molecule for node {n} ({smiles}): {e}")
                continue
            node_x, node_y = pos[n]
            img_x = node_x
            img_y = img_y_top
            self.add_svg_image_to_axes(ax, png_bytes, xy=(img_x, img_y), zoom=img_zoom)
            ax.plot([node_x, img_x], [node_y, img_y], linestyle=(0, (3, 7)),
                color='gray', linewidth=1.5, zorder=9)

        # Products (bottom)
        product_nodes = [n for n in self.G.nodes if node_layer[n] == 4]
        for i, n in enumerate(product_nodes):
            smiles = self.G.nodes[n].get('smiles')
            if not smiles:
                continue
            mol = Chem.MolFromSmiles(smiles)
            img_size = self.get_img_size_for_mol(mol)
            try:
                png_bytes = self.mol_to_png_bytes(smiles, img_size=img_size)
            except Exception as e:
                print(f"Error drawing molecule for node {n} ({smiles}): {e}")
                continue
            node_x, node_y = pos[n]
            img_x = node_x
            img_y = img_y_bottom
            self.add_svg_image_to_axes(ax, png_bytes, xy=(img_x, img_y), zoom=img_zoom)
            ax.plot([node_x, img_x], [node_y, img_y], linestyle=(0, (3, 7)),
                color='gray', linewidth=1.5, zorder=9)

        # Legend
        legend_elements = [
            plt.Line2D([0], [1], color='green', lw=3, label='Both Datasets'),
            plt.Line2D([0], [1], color='blue', lw=3, label='BNN Dataset'),
            plt.Line2D([0], [1], color='red', lw=3, label='Reaxys Dataset'),
            plt.Line2D([0], [0], marker='s', color='w', markerfacecolor='black', markersize=6, label='Connecting Multiple Reactants', linestyle='None'),
            plt.Line2D([0], [0], marker='^', color='w', markerfacecolor='black', markersize=6, label='Connecting Multiple Products', linestyle='None'),
            plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='black', markersize=6, label='Molecule Node', linestyle='None'),
        ]
        ax.legend(
            handles=legend_elements, loc='upper right',
            fontsize=legend_fontsize, bbox_to_anchor=legend_bbox_to_anchor
        )

        plt.show()

    def get_img_size_for_mol(self, mol, min_size=120, max_size=1000):
        n_atoms = mol.GetNumAtoms()
        atoms_per_px = 40
        if n_atoms > 20:
            atoms_per_px = 50
        size = min_size + atoms_per_px * n_atoms
        return min(size, max_size)

    # def mol_to_png_bytes(
    #         self,
    #         smi,
    #         img_size=150
    #     ):
    #     mol = Chem.MolFromSmiles(smi)
    #     if mol is None:
    #         raise ValueError("Invalid SMILES: %s" % smi)
    #     Chem.rdDepictor.Compute2DCoords(mol)
    #     Chem.rdDepictor.StraightenDepiction(mol)
    #     drawer = rdMolDraw2D.MolDraw2DSVG(img_size, img_size)
    #     #Draw.SetACS1996Mode(drawer.drawOptions(), Draw.MeanBondLength(mol))
    #     rdMolDraw2D.PrepareAndDrawMolecule(drawer, mol)
    #     drawer.FinishDrawing()
    #     svg = drawer.GetDrawingText()
    #     # Convert SVG to high-res PNG by specifying output_width/output_height
    #     png_bytes = cairosvg.svg2png(bytestring=svg.encode('utf-8'),
    #                                 output_width=img_size, output_height=img_size)
    #     return png_bytes

    def mol_to_png_bytes(
        self,
        smi,
        img_size=150,
        bond_line_width=3.0,          # base line width in px
        scale_line_width=True,        # scale with image size
        ):
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            raise ValueError("Invalid SMILES: %s" % smi)
        Chem.rdDepictor.Compute2DCoords(mol)
        Chem.rdDepictor.StraightenDepiction(mol)
        drawer = rdMolDraw2D.MolDraw2DSVG(img_size, img_size)

        # Thicken bonds (RDKit option). Scale with image size so larger images have proportionally thicker lines.
        lw = float(bond_line_width) * (img_size / 250.0) if scale_line_width else float(bond_line_width)
        opts = drawer.drawOptions()
        if hasattr(opts, "bondLineWidth"):
            opts.bondLineWidth = lw
        elif hasattr(opts, "lineWidth"):  # fallback for older RDKit
            opts.lineWidth = lw

        rdMolDraw2D.PrepareAndDrawMolecule(drawer, mol)
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText()

        # Safety: ensure all strokes are thick (covers wedge/double bonds too)
        try:
            import re as _re
            svg = _re.sub(r"stroke-width:\s*[\d.]+px", f"stroke-width:{lw:.2f}px", svg)
        except Exception:
            pass

        # Convert SVG to PNG
        png_bytes = cairosvg.svg2png(
            bytestring=svg.encode('utf-8'),
            output_width=img_size,
            output_height=img_size
        )
        return png_bytes

    def trim_whitespace(self, im):
        bg = Image.new(im.mode, im.size, im.getpixel((0,0)))
        diff = ImageChops.difference(im, bg)
        diff = ImageChops.add(diff, diff)
        bbox = diff.getbbox()
        if bbox:
            return im.crop(bbox)
        
    def add_svg_image_to_axes(self, ax, png_bytes, xy, zoom=0.5):
        image = Image.open(BytesIO(png_bytes))
        trimmed = self.trim_whitespace(image) or image
        arr = np.array(trimmed)
        imagebox = OffsetImage(arr, zoom=zoom, dpi_cor=False)
        ab = AnnotationBbox(imagebox, xy, frameon=False)
        ax.add_artist(ab)
        return ab
    

    def chunk_list(self, lst, n):
        """Yield successive n-sized chunks from lst."""
        for i in range(0, len(lst), n):
            yield lst[i:i + n]

    def alternating_constant_offsets(self, n, amplitude):
        # Alternates: [-A, +A, -A, +A, ...]
        if n == 1:
            return [0]
        return [(-amplitude if i % 2 == 0 else amplitude) for i in range(n)]
        
    def get_node_layer(self, G, n):
        """
        Classify nodes based on number of in and out connections and shape
        """
        node = G.nodes[n]
        shape = node.get('shape', 'o')
        if shape == 's':
            return 2  # intermediate reactant (square)
        elif shape == '^':
            return 3  # intermediate product (triangle)
        else:  # molecule
            in_deg = G.in_degree(n)
            out_deg = G.out_degree(n)
            if in_deg == 0 and out_deg > 0:
                return 0  # only reactant (top)
            elif in_deg > 0 and out_deg == 0:
                return 4  # only product (bottom)
            elif in_deg > 0 and out_deg > 0:
                return 1  # both reactant and product (middle)
            
    def draw_reaction_network_circular(
        self,
        top_smiles=None,   # List of SMILES to place at the top
        bottom_smiles=None, # List of SMILES to place at the bottom
        img_dpi=300,
        figsize=(8, 12),
        label_fontsize=10,
        # Layout
        base_ys={0: 0.8, 4: 0.2},
        svg_bond_line_width_top=2.0,
        svg_bond_line_width_bottom=2.0,
        img_y_top=0.9,
        img_y_bottom=0.1,
        node_size=200,
        node_alpha=0.8,
        edge_alpha=0.3,
        edge_width=2,
        img_zoom=1.0,
        legend_fontsize=10,
        legend_bbox_to_anchor=(1.15, 1.0),
        circle_radius=0.25,
        circle_center=(0.5, 0.5),
        inner_node_scatter_radius=0.18,
        seed=42,
        show_legend=True,
        ylim=(-0.02, 1.02),
        savefig=False,
        filename="reaction_network.png",
    ):
        """
        Draws a reaction network with selected molecules at top/bottom (by SMILES),
        other molecules ('o') on the circle,
        and all intermediates ('s', '^') inside the circle.
        """

        G = self.G
        if top_smiles is None:
            top_smiles = []
        if bottom_smiles is None:
            bottom_smiles = []

        # Convert SMILES to node IDs (InChI) for top and bottom
        def smiles_to_id(smiles):
            inchi = self.get_inchi_key(self.get_inchi(smiles))
            for n in G.nodes:
                if G.nodes[n].get("id", n) == inchi or G.nodes[n].get("smiles", "") == smiles:
                    return n
            return None

        top_ids = [smiles_to_id(s) for s in top_smiles if smiles_to_id(s) is not None]
        bottom_ids = [smiles_to_id(s) for s in bottom_smiles if smiles_to_id(s) is not None]

        # 1. Categorize nodes
        all_circle_nodes = []
        inner_nodes = []
        for n in G.nodes:
            shape = G.nodes[n].get("shape", "o")
            if shape == 'o':
                if n in top_ids or n in bottom_ids:
                    continue
                all_circle_nodes.append(n)
            elif shape in ('s', '^'):
                inner_nodes.append(n)

        # 2. Position nodes:
        pos = {}

        # Top (starting material)
        n_top = len(top_ids)
        for i, n in enumerate(top_ids):
            pos[n] = ((i+1)/(n_top+1), base_ys[0])

        # Bottom (product)
        n_bottom = len(bottom_ids)
        for i, n in enumerate(bottom_ids):
            pos[n] = ((i+1)/(n_bottom+1), base_ys[4])

        # Circle (all other non-terminal molecules, shape 'o')
        n_circle = len(all_circle_nodes)
        if n_circle > 0:
            angles = np.linspace(0, 2*np.pi, n_circle, endpoint=False)
            for idx, n in enumerate(all_circle_nodes):
                x = circle_center[0] + circle_radius * np.cos(angles[idx])
                y = circle_center[1] + circle_radius * np.sin(angles[idx])
                pos[n] = (x, y)

        # Inner (intermediate/connecting: 's', '^')
        n_inner = len(inner_nodes)
        import random
        rng = random.Random(seed)
        for n in inner_nodes:
            r = inner_node_scatter_radius * np.sqrt(rng.random())
            theta = rng.uniform(0, 2*np.pi)
            x = circle_center[0] + r * np.cos(theta)
            y = circle_center[1] + r * np.sin(theta)
            pos[n] = (x, y)

        # 3. Plotting
        fig, ax = plt.subplots(figsize=figsize, dpi=img_dpi)
        ax.axis('off')
        # Ensure images at y=img_y_top/img_y_bottom are visible even without connector lines
        #ax.set_xlim(-0.02, 1.02)
        ax.set_ylim(ylim)

        # 4. Draw edges by color
        edge_color_map = {}
        for u, v, d in G.edges(data=True):
            color = d.get('color', 'black')
            edge_color_map.setdefault(color, []).append((u, v))
        for color, edgelist in edge_color_map.items():
            nx.draw_networkx_edges(
                G, pos, ax=ax, edgelist=edgelist,
                edge_color=color, alpha=edge_alpha, arrows=True, width=edge_width
            )

        # 5. Draw nodes by shape
        shape_map = {'o': [], 's': [], '^': []}
        color_map = {'o': [], 's': [], '^': []}
        for n, d in G.nodes(data=True):
            shape = d.get('shape', 'o')
            color = d.get('color', 'black')
            if n in pos:  # Only plot placed nodes
                shape_map.setdefault(shape, []).append(n)
                color_map.setdefault(shape, []).append(color)
        for shape, nodelist in shape_map.items():
            if nodelist:
                nx.draw_networkx_nodes(
                    G, pos, nodelist=nodelist, node_shape=shape,
                    node_color=color_map[shape], node_size=node_size,
                    ax=ax, alpha=node_alpha
                )

        labels = {}
        for n, d in G.nodes(data=True):
            if "label" in d and d["label"]:
                labels[n] = d["label"]
        nx.draw_networkx_labels(
            G, pos, labels=labels, font_color="white", font_size=label_fontsize, ax=ax
        )
        # 6. Draw molecule images at top and bottom only for specified molecules
        for n in top_ids:
            smiles = G.nodes[n].get('smiles')
            if not smiles:
                continue
            mol = Chem.MolFromSmiles(smiles)
            img_size = self.get_img_size_for_mol(mol)
            try:
                png_bytes = self.mol_to_png_bytes(smiles, img_size=img_size, bond_line_width=svg_bond_line_width_top)
            except Exception as e:
                print(f"Error drawing molecule for node {n} ({smiles}): {e}")
                continue
            node_x, node_y = pos[n]
            img_x = node_x
            img_y = img_y_top
            self.add_svg_image_to_axes(ax, png_bytes, xy=(img_x, img_y), zoom=img_zoom)

            # ax.plot([node_x, img_x], [node_y, img_y], linestyle=(0, (3, 7)),
            #     color='gray', linewidth=1.5, alpha=0.8, zorder=9)

        for n in bottom_ids:
            smiles = G.nodes[n].get('smiles')
            if not smiles:
                continue
            mol = Chem.MolFromSmiles(smiles)
            img_size = self.get_img_size_for_mol(mol)
            try:
                png_bytes = self.mol_to_png_bytes(smiles, img_size=img_size, bond_line_width=svg_bond_line_width_bottom)
            except Exception as e:
                print(f"Error drawing molecule for node {n} ({smiles}): {e}")
                continue
            node_x, node_y = pos[n]
            img_x = node_x
            img_y = img_y_bottom
            self.add_svg_image_to_axes(ax, png_bytes, xy=(img_x, img_y), zoom=img_zoom)
            # ax.plot(
            #     [node_x, img_x], [node_y+0.2, img_y-0.2], linestyle=(0, (3, 7)),
            #     color='gray', linewidth=1.5, alpha=0.8, zorder=9)

        # 7. Legend
        if show_legend:
            legend_elements = [
                plt.Line2D([0], [1], color='green', lw=3, label='Both Datasets'),
                plt.Line2D([0], [1], color='blue', lw=3, label='BNN Dataset'),
                plt.Line2D([0], [1], color='red', lw=3, label='Reaxys Dataset'),
                plt.Line2D([0], [0], marker='s', color='w', markerfacecolor='grey', markersize=10, label='Connecting Multiple Reactants', linestyle='None'),
                plt.Line2D([0], [0], marker='^', color='w', markerfacecolor='grey', markersize=10, label='Connecting Multiple Products', linestyle='None'),
                plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='grey', markersize=10, label='Molecule Node', linestyle='None'),
            ]
            ax.legend(
                handles=legend_elements, loc='upper right',
                fontsize=legend_fontsize, bbox_to_anchor=legend_bbox_to_anchor
            )
        if savefig:
            plt.savefig(filename, bbox_inches='tight')
        else:
            plt.show()

    def draw_reaction_network_circular_clean(
        self,
        top_smiles=None,
        bottom_smiles=None,
        auto_terminals=True,
        figsize=(9, 10),
        img_dpi=300,
        # ring layout
        circle_center=(0.5, 0.5),
        circle_radius=0.33,
        inner_radius=0.22,
        oval_ratio=1.0,  # 0<ratio<=1; y-radius = ratio * x-radius (1.0 => circle, <1 => flatter oval)
    ring_shuffle=False,
    ring_shuffle_seed=42,
        # terminal bands
        top_y=0.92,
        bottom_y=0.08,
        # style
        node_size_molecule=230,
        node_size_intermediate=170,
        node_alpha=0.9,
        edge_alpha=0.35,
        edge_width=2.0,
    # optional edge border (outline) controls
    edge_border=False,
    edge_border_color="black",
    edge_border_width=None,  # if None, defaults to edge_width + 1.0
        # dataset colors (applied to nodes and edges)
        color_dataset1="blue",
        color_dataset2="red",
        color_both="green",
    # per-dataset node alphas (None -> fallback to node_alpha)
    color_dataset1_alpha=None,
    color_dataset2_alpha=None,
    color_both_alpha=None,
    edge_alpha_dataset1=0.35, 
    edge_alpha_dataset2=0.35, 
    edge_alpha_both=0.6,
    # node border controls
    node_edgecolor="black",
    node_edgewidth=1.0,
        label_mode="all",  # 'none' | 'terminals' | 'all'
        label_fontsize=9,
    label_fontcolor="black",
        # layout refinement for inner nodes
        use_spring_for_inner=True,
        spring_iterations=200,
        spring_k=None,
        clamp_inner_to_radius=True,
        seed=42,
        # output
        show_legend=True,
        savefig=False,
        filename="reaction_network_clean.png",
    ):
        """
        Draw a clean circular network without molecule images.
        - Molecule nodes ('o') are placed on a ring.
        - Intermediate nodes ('s', '^') are placed inside the ring and relaxed
          via a constrained spring layout to reduce overlaps.
        - "From" (only reactant: in_deg=0) molecules are lifted above the ring,
          and "To" (only product: out_deg=0) molecules are dropped below the ring.

        Options:
        - top_smiles / bottom_smiles: optional explicit terminals by SMILES.
        - auto_terminals: when True, infer terminals from graph degrees.
                - label_mode: 'none' (no labels), 'terminals' (only top/bottom), 'all' (every labeled node).
                - color_dataset1, color_dataset2, color_both: RGB/hex names to color nodes and edges for
                    dataset1, dataset2, and shared ('both') respectively.
                - oval_ratio: scales the vertical radius of the ring and inner region. Use values <1 (e.g., 0.6)
                    to flatten the ring into an oval and improve edge visibility.
                - edge_border / edge_border_color / edge_border_width: enable an outline behind all edges with a
                    uniform color and thickness (drawn first), then draw colored edges on top.
        """

        # Ensure the graph exists
        if getattr(self, "G", None) is None:
            # Build full graph with defaults
            self.build_graph_from_nodes_edges()
        G = self.G

        # Dataset -> color mapping (used for both nodes and edges)
        dataset_color = {
            self.dataset1_name: color_dataset1,
            self.dataset2_name: color_dataset2,
            "both": color_both,
        }
        # Node alpha per dataset (fallback to global node_alpha if None)
        dataset_alpha = {
            self.dataset1_name: (node_alpha if color_dataset1_alpha is None else float(color_dataset1_alpha)),
            self.dataset2_name: (node_alpha if color_dataset2_alpha is None else float(color_dataset2_alpha)),
            "both": (node_alpha if color_both_alpha is None else float(color_both_alpha)),
        }

        # Helpers
        def smiles_to_id(smiles: str):
            return self.smiles_to_node_id(smiles, self.nodes)

        def even_x_positions(n: int, cx: float, width: float):
            if n <= 0:
                return []
            # distribute within [cx - width/2, cx + width/2]
            left = cx - width / 2.0
            span = width
            return [left + (i + 1) * span / (n + 1) for i in range(n)]

        # Determine terminal nodes (top/bottom)
        top_ids, bottom_ids = set(), set()
        # If explicit lists provided, map to ids
        if top_smiles:
            for s in top_smiles:
                nid = smiles_to_id(s)
                if nid:
                    top_ids.add(nid)
        if bottom_smiles:
            for s in bottom_smiles:
                nid = smiles_to_id(s)
                if nid:
                    bottom_ids.add(nid)

        # Auto-detect terminals if asked (and allow union with explicit ones)
        if auto_terminals:
            for n in G.nodes:
                d = G.nodes[n]
                if d.get("shape", "o") != "o":
                    continue
                in_deg = G.in_degree(n)
                out_deg = G.out_degree(n)
                if in_deg == 0 and out_deg > 0:
                    top_ids.add(n)
                elif in_deg > 0 and out_deg == 0:
                    bottom_ids.add(n)

        # Categorize nodes
        circle_nodes = []
        inner_nodes = []
        for n in G.nodes:
            shape = G.nodes[n].get("shape", "o")
            if shape == "o":
                if n not in top_ids and n not in bottom_ids:
                    circle_nodes.append(n)
            elif shape in ("s", "^"):
                inner_nodes.append(n)

        # Initial positions
        pos = {}

        # Place ring molecules on an ellipse (oval)
        n_circle = len(circle_nodes)
        if n_circle > 0:
            # Determine placement order: sorted or shuffled
            ring_nodes = list(circle_nodes)
            if ring_shuffle:
                import random as _random
                _rng = _random.Random(ring_shuffle_seed)
                _rng.shuffle(ring_nodes)
            else:
                ring_nodes = sorted(ring_nodes)
            angles = np.linspace(0, 2 * np.pi, n_circle, endpoint=False)
            for i, n in enumerate(ring_nodes):
                x = circle_center[0] + circle_radius * np.cos(angles[i])
                y = circle_center[1] + (circle_radius * float(oval_ratio)) * np.sin(angles[i])
                pos[n] = (x, y)

        # Place terminals above/below (evenly across the ring width)
        top_ids_sorted = sorted(top_ids)
        bottom_ids_sorted = sorted(bottom_ids)
        top_xs = even_x_positions(len(top_ids_sorted), circle_center[0], 2 * circle_radius)
        bot_xs = even_x_positions(len(bottom_ids_sorted), circle_center[0], 2 * circle_radius)
        for x, n in zip(top_xs, top_ids_sorted):
            pos[n] = (x, top_y)
        for x, n in zip(bot_xs, bottom_ids_sorted):
            pos[n] = (x, bottom_y)

        # Place inner nodes randomly within inner ellipse
        rng = np.random.default_rng(seed)
        for n in inner_nodes:
            # random point inside circle
            r = inner_radius * np.sqrt(rng.random())
            theta = rng.random() * 2 * np.pi
            x = circle_center[0] + r * np.cos(theta)
            y = circle_center[1] + (r * float(oval_ratio)) * np.sin(theta)
            pos[n] = (x, y)

        # Optionally relax inner nodes with spring layout, keeping ring + terminals fixed
        if use_spring_for_inner and (inner_nodes):
            fixed_nodes = set(circle_nodes) | set(top_ids) | set(bottom_ids)
            # Ensure every graph node has a starting position (fallback at center)
            for n in G.nodes:
                if n not in pos:
                    pos[n] = circle_center
            pos = nx.spring_layout(
                G,
                pos=pos,
                fixed=list(fixed_nodes),
                iterations=int(spring_iterations),
                k=spring_k,
                scale=1.0,
                center=circle_center,
                seed=seed,
            )
            # Keep inner nodes inside the inner ellipse if requested
            if clamp_inner_to_radius:
                cx, cy = circle_center
                rx = float(inner_radius)
                ry = float(inner_radius) * float(oval_ratio)
                for n in inner_nodes:
                    x, y = pos.get(n, circle_center)
                    dx, dy = x - cx, y - cy
                    # Check ellipse equation: (dx/rx)^2 + (dy/ry)^2 <= 1
                    denom = (dx / rx) ** 2 + (dy / ry) ** 2 if rx > 0 and ry > 0 else 0.0
                    if denom > 1.0 and denom > 0:
                        t = denom ** -0.5  # project onto ellipse boundary along the ray
                        pos[n] = (cx + dx * t, cy + dy * t)

        # Draw
        fig, ax = plt.subplots(figsize=figsize, dpi=img_dpi)
        ax.axis("off")
        ax.set_xlim(0.0, 1.0)
        ax.set_ylim(0.0, 1.0)
        ax.set_aspect("equal", adjustable="box")

        # Group edges by dataset for independent alpha control
        edge_groups = {self.dataset1_name: [], self.dataset2_name: [], "both": []}
        for u, v, d in G.edges(data=True):
            ds = d.get("dataset")
            if ds in edge_groups:
                edge_groups[ds].append((u, v))
        all_edges = edge_groups[self.dataset1_name] + edge_groups[self.dataset2_name] + edge_groups["both"]

        # Optional edge border (drawn first as a thicker underlay)
        if edge_border and all_edges:
            bw = float(edge_border_width) if edge_border_width is not None else (float(edge_width) + 1.0)
            nx.draw_networkx_edges(
                G,
                pos,
                ax=ax,
                edgelist=all_edges,
                edge_color=edge_border_color,
                alpha=edge_alpha,
                arrows=True,
                width=bw,
            )
        # Per-dataset alpha mapping (fallback to global edge_alpha)
        alpha_map = {
            self.dataset1_name: (edge_alpha if edge_alpha_dataset1 is None else float(edge_alpha_dataset1)),
            self.dataset2_name: (edge_alpha if edge_alpha_dataset2 is None else float(edge_alpha_dataset2)),
            "both": (edge_alpha if edge_alpha_both is None else float(edge_alpha_both)),
        }
        # Draw each dataset group with its own alpha and color
        if edge_groups[self.dataset1_name]:
            nx.draw_networkx_edges(
                G, pos, ax=ax, edgelist=edge_groups[self.dataset1_name],
                edge_color=dataset_color[self.dataset1_name], alpha=alpha_map[self.dataset1_name],
                arrows=True, width=edge_width,
            )
        if edge_groups[self.dataset2_name]:
            nx.draw_networkx_edges(
                G, pos, ax=ax, edgelist=edge_groups[self.dataset2_name],
                edge_color=dataset_color[self.dataset2_name], alpha=alpha_map[self.dataset2_name],
                arrows=True, width=edge_width,
            )
        if edge_groups["both"]:
            nx.draw_networkx_edges(
                G, pos, ax=ax, edgelist=edge_groups["both"],
                edge_color=dataset_color["both"], alpha=alpha_map["both"],
                arrows=True, width=edge_width,
            )

        # Nodes by shape
        mol_nodes = [n for n in G.nodes if G.nodes[n].get("shape", "o") == "o" and n in pos]
        int_nodes = [n for n in G.nodes if G.nodes[n].get("shape", "o") in ("s", "^") and n in pos]
        # Build per-node RGBA colors with dataset-specific alpha
        from matplotlib import colors as mcolors
        mol_colors = [
            mcolors.to_rgba(
                dataset_color.get(G.nodes[n].get("dataset"), "grey"),
                dataset_alpha.get(G.nodes[n].get("dataset"), node_alpha),
            )
            for n in mol_nodes
        ]
        nx.draw_networkx_nodes(
            G,
            pos,
            nodelist=mol_nodes,
            node_shape="o",
            node_color=mol_colors,
            node_size=node_size_molecule,
            edgecolors=node_edgecolor,
            linewidths=float(node_edgewidth),
            ax=ax,
        )
        # Draw squares and triangles separately to honor shapes
        sq_nodes = [n for n in int_nodes if G.nodes[n].get("shape") == "s"]
        tr_nodes = [n for n in int_nodes if G.nodes[n].get("shape") == "^"]
        if sq_nodes:
            sq_colors = [
                mcolors.to_rgba(
                    dataset_color.get(G.nodes[n].get("dataset"), "grey"),
                    dataset_alpha.get(G.nodes[n].get("dataset"), node_alpha),
                )
                for n in sq_nodes
            ]
            nx.draw_networkx_nodes(
                G,
                pos,
                nodelist=sq_nodes,
                node_shape="s",
                node_color=sq_colors,
                node_size=node_size_intermediate,
                edgecolors=node_edgecolor,
                linewidths=float(node_edgewidth),
                ax=ax,
            )
        if tr_nodes:
            tr_colors = [
                mcolors.to_rgba(
                    dataset_color.get(G.nodes[n].get("dataset"), "grey"),
                    dataset_alpha.get(G.nodes[n].get("dataset"), node_alpha),
                )
                for n in tr_nodes
            ]
            nx.draw_networkx_nodes(
                G,
                pos,
                nodelist=tr_nodes,
                node_shape="^",
                node_color=tr_colors,
                node_size=node_size_intermediate,
                edgecolors=node_edgecolor,
                linewidths=float(node_edgewidth),
                ax=ax,
            )

        # Labels
        if label_mode != "none":
            labels = {}
            if label_mode == "all":
                for n, d in G.nodes(data=True):
                    if n in pos and d.get("label"):
                        labels[n] = d["label"]
            elif label_mode == "terminals":
                for n in list(top_ids) + list(bottom_ids):
                    d = G.nodes[n]
                    if d.get("label"):
                        labels[n] = d["label"]
            if labels:
                nx.draw_networkx_labels(
                    G,
                    pos,
                    labels=labels,
                    font_color=label_fontcolor,
                    font_size=label_fontsize,
                    ax=ax,
                )

        if show_legend:
            legend_elements = [
                plt.Line2D([0], [1], color=color_both, lw=3, label="Both Datasets"),
                plt.Line2D([0], [1], color=color_dataset1, lw=3, label=self.dataset1_name),
                plt.Line2D([0], [1], color=color_dataset2, lw=3, label=self.dataset2_name),
                plt.Line2D([0], [0], marker="s", color="w", markerfacecolor="grey", markersize=9, label="Connecting Reactants", linestyle="None"),
                plt.Line2D([0], [0], marker="^", color="w", markerfacecolor="grey", markersize=9, label="Connecting Products", linestyle="None"),
                plt.Line2D([0], [0], marker="o", color="w", markerfacecolor="grey", markersize=9, label="Molecule", linestyle="None"),
            ]
            ax.legend(handles=legend_elements, loc="upper right", fontsize=10)

        if savefig:
            plt.savefig(filename, bbox_inches="tight")
        else:
            plt.show()

    def count_common_and_different_nodes_edges(self):
        """
        Counts the number of nodes and edges in common and unique to each dataset in self.G,
        grouped by both dataset and node shape.
        Returns:
            dict with counts for:
                - nodes_by_dataset_and_shape: {dataset: {shape: count}}
                - common_edges
                - unique_edges_dataset1
                - unique_edges_dataset2
        """
        # Prepare node counts by dataset and shape
        node_counts = {}
        for n, d in self.G.nodes(data=True):
            dataset = d.get('dataset')
            shape = d.get('shape', 'o')
            if dataset not in node_counts:
                node_counts[dataset] = {}
            if shape not in node_counts[dataset]:
                node_counts[dataset][shape] = 0
            node_counts[dataset][shape] += 1

        # Edges: use (source, target) tuples
        edges_1 = {(u, v) for u, v, d in self.G.edges(data=True) if d.get('dataset') == self.dataset1_name}
        edges_2 = {(u, v) for u, v, d in self.G.edges(data=True) if d.get('dataset') == self.dataset2_name}
        edges_both = {(u, v) for u, v, d in self.G.edges(data=True) if d.get('dataset') == 'both'}

        # Common and unique edges
        common_edges = len(edges_both)
        unique_edges_1 = len(edges_1 - edges_2 - edges_both)
        unique_edges_2 = len(edges_2 - edges_1 - edges_both)

        return {
            "nodes_by_dataset_and_shape": node_counts,
            "common_edges": common_edges,
            "unique_edges_{}".format(self.dataset1_name): unique_edges_1,
            "unique_edges_{}".format(self.dataset2_name): unique_edges_2,
        }
    
    def get_common_and_unique_nodes(self):
        """
        Returns:
            dict with:
                - common_nodes: set of SMILES in both datasets
                - unique_nodes_dataset1: set of SMILES unique to dataset 1
                - unique_nodes_dataset2: set of SMILES unique to dataset 2
        """
        nodes_1 = {n for n, d in self.G.nodes(data=True) if d.get('dataset') == self.dataset1_name}
        nodes_2 = {n for n, d in self.G.nodes(data=True) if d.get('dataset') == self.dataset2_name}
        nodes_both = {n for n, d in self.G.nodes(data=True) if d.get('dataset') == 'both'}

        unique_nodes_1 = nodes_1 - nodes_2 - nodes_both
        unique_nodes_2 = nodes_2 - nodes_1 - nodes_both

        # Helper to get SMILES from node
        def get_smiles(n):
            return self.G.nodes[n].get('smiles', None)

        return {
            "common_nodes": {get_smiles(n) for n in nodes_both if get_smiles(n)},
            "unique_nodes_{}".format(self.dataset1_name): {get_smiles(n) for n in unique_nodes_1 if get_smiles(n)},
            "unique_nodes_{}".format(self.dataset2_name): {get_smiles(n) for n in unique_nodes_2 if get_smiles(n)},
        }

    def get_common_and_unique_edges(self, return_smiles=False):
        """
        Returns:
            dict with:
                - common_edges: set of (source, target) in both datasets
                - unique_edges_dataset1: set of (source, target) unique to dataset 1
                - unique_edges_dataset2: set of (source, target) unique to dataset 2
            If return_smiles=True, returns a set of reaction SMILES for each group.
        """
        edges_1 = {(u, v, d.get('reaction_key')) for u, v, d in self.G.edges(data=True) if d.get('dataset') == self.dataset1_name}
        edges_2 = {(u, v, d.get('reaction_key')) for u, v, d in self.G.edges(data=True) if d.get('dataset') == self.dataset2_name}
        edges_both = {(u, v, d.get('reaction_key')) for u, v, d in self.G.edges(data=True) if d.get('dataset') == 'both'}

        unique_edges_1 = edges_1 - edges_2 - edges_both
        unique_edges_2 = edges_2 - edges_1 - edges_both

        def rxn_smiles_set(edge_set):
            smiles = set()
            for _, _, rxn_key in edge_set:
                if rxn_key in self.reaction_keys_1:
                    smiles.add(self.reaction_keys_1[rxn_key])
                elif rxn_key in self.reaction_keys_2:
                    smiles.add(self.reaction_keys_2[rxn_key])
            return smiles

        if return_smiles:
            return {
                "common_edges": rxn_smiles_set(edges_both),
                f"unique_edges_{self.dataset1_name}": rxn_smiles_set(unique_edges_1),
                f"unique_edges_{self.dataset2_name}": rxn_smiles_set(unique_edges_2),
            }
        else:
            return {
                "common_edges": {(u, v) for u, v, _ in edges_both},
                f"unique_edges_{self.dataset1_name}": {(u, v) for u, v, _ in unique_edges_1},
                f"unique_edges_{self.dataset2_name}": {(u, v) for u, v, _ in unique_edges_2},
            }

    def is_unique_molecule_between_shared(
        self,
        node_id: str,
        require_direct: bool = False,
        enforce_edge_dataset_match: bool = False,
    ) -> bool:
        """
        Check whether a given node is a unique molecule that lies "between" two shared molecule nodes.

        Definition used:
        - The node must be a molecule (shape == 'o').
        - The node's dataset must be unique to one dataset (i.e., dataset == self.dataset1_name or self.dataset2_name), not 'both'.
        - There must exist at least one upstream shared molecule and at least one downstream shared molecule.
          Upstream means a shared molecule reachable via either:
            * a direct edge: shared('o') -> node_id, or
            * a 2-step path via a single intermediate node: shared('o') -> ('s'/'^') -> node_id.
          Downstream is defined analogously in the forward direction from node_id to a shared('o').

        Parameters
        - node_id: The node to test.
        - require_direct: If True, only consider direct molecule-to-molecule edges ('o' -> 'o' and 'o' <- 'o').
                          If False (default), also consider one intermediate hop via 's'/'^'.
        - enforce_edge_dataset_match: If True, require that the edges along the qualifying upstream and downstream paths
                                      belong to the same dataset as the node (or 'both').

        Returns
        - bool: True if the criteria are satisfied; False otherwise.
        """
        # Ensure graph exists
        if getattr(self, "G", None) is None:
            self.build_graph_from_nodes_edges()
        G = self.G

        if node_id not in G:
            return False

        nd = G.nodes[node_id]
        if nd.get("shape", "o") != "o":
            return False
        node_ds = nd.get("dataset")
        if node_ds not in (self.dataset1_name, self.dataset2_name):
            # Not unique to one dataset
            return False

        def edges_match_dataset(edges: list[tuple[str, str]]) -> bool:
            if not enforce_edge_dataset_match:
                return True
            for (u, v) in edges:
                eds = G.edges[u, v].get("dataset")
                if eds not in (node_ds, "both"):
                    return False
            return True

        def has_shared_upstream() -> bool:
            # Direct predecessors
            for p in G.predecessors(node_id):
                pd = G.nodes[p]
                if pd.get("shape") == "o" and pd.get("dataset") == "both":
                    if edges_match_dataset([(p, node_id)]):
                        return True
                if not require_direct and pd.get("shape") in ("s", "^"):
                    for pp in G.predecessors(p):
                        ppd = G.nodes[pp]
                        if ppd.get("shape") == "o" and ppd.get("dataset") == "both":
                            if edges_match_dataset([(pp, p), (p, node_id)]):
                                return True
            return False

        def has_shared_downstream() -> bool:
            # Direct successors
            for s in G.successors(node_id):
                sd = G.nodes[s]
                if sd.get("shape") == "o" and sd.get("dataset") == "both":
                    if edges_match_dataset([(node_id, s)]):
                        return True
                if not require_direct and sd.get("shape") in ("s", "^"):
                    for ss in G.successors(s):
                        ssd = G.nodes[ss]
                        if ssd.get("shape") == "o" and ssd.get("dataset") == "both":
                            if edges_match_dataset([(node_id, s), (s, ss)]):
                                return True
            return False

        return has_shared_upstream() and has_shared_downstream()

    def get_unique_molecules_between_shared(
        self,
        require_direct: bool = False,
        enforce_edge_dataset_match: bool = False,
        return_labels: bool = False,
        return_dataset_counts = False,
    ):
        """
        Convenience: list all molecule node IDs that are unique to one dataset and lie
        between shared molecule nodes (per is_unique_molecule_between_shared).

        Parameters
        - require_direct: See is_unique_molecule_between_shared.
        - enforce_edge_dataset_match: See is_unique_molecule_between_shared.
        - return_labels: If True, returns list of (node_id, label) tuples; else just node_id list.
        - return_dataset_counts: If True, returns a dictionary of interemdiate molecule node counts for each dataset
        Returns
        - list of node ids or (id, label) pairs.
        """
        if getattr(self, "G", None) is None:
            self.build_graph_from_nodes_edges()
        G = self.G

        result = []
        dataset_counts = {
            self.dataset1_name: 0,
            self.dataset2_name: 0
        }
        for n, d in G.nodes(data=True):
            if d.get("shape", "o") != "o":
                continue
            if d.get("dataset") not in (self.dataset1_name, self.dataset2_name):
                continue
            if self.is_unique_molecule_between_shared(n, require_direct=require_direct, enforce_edge_dataset_match=enforce_edge_dataset_match):
                if return_labels:
                    result.append((n, d.get("dataset"), d.get("smiles")))
                if return_dataset_counts:
                    label = d.get("dataset")
                    if label == self.dataset1_name:
                        dataset_counts[self.dataset1_name] += 1
                    if label == self.dataset2_name:
                        dataset_counts[self.dataset2_name] += 1
                else:
                    result.append(n)
        if return_dataset_counts:
            return dataset_counts, result
        else:
            return result


# def build_reaction_network_data(reactions_1, reactions_2, dataset1_name="dataset_1", dataset2_name="dataset_2"):
#     """
#     Builds node and edge lists for a reaction network with correct shapes and colors.

#     Returns:
#         nodes: list of dicts (id, shape, color, type, dataset)
#         edges: list of dicts (source, target, color, dataset, reaction_key)
#     """
#     def get_inchi(smi):
#         mol = Chem.MolFromSmiles(smi)
#         return Chem.MolToInchiKey(mol) if mol else None

#     # 1. Collect all molecules and their dataset memberships
#     all_smiles_1 = [smi for rxn in reactions_1 for part in rxn.split(">>") for smi in part.split(".") if smi.strip()]
#     all_smiles_2 = [smi for rxn in reactions_2 for part in rxn.split(">>") for smi in part.split(".") if smi.strip()]
#     inchi_set_1 = {get_inchi(smi) for smi in all_smiles_1}
#     inchi_set_2 = {get_inchi(smi) for smi in all_smiles_2}
#     all_inchis = inchi_set_1 | inchi_set_2

#     # 2. Make node records for molecules
#     nodes = {}
#     for inchi in all_inchis:
#         if not inchi: continue
#         if inchi in inchi_set_1 and inchi in inchi_set_2:
#             dataset = "both"
#             color = "green"
#         elif inchi in inchi_set_1:
#             dataset = dataset1_name
#             color = "blue"
#         elif inchi in inchi_set_2:
#             dataset = dataset2_name
#             color = "red"
#         else:
#             dataset = "other"
#             color = "black"
#         nodes[inchi] = {
#             "id": inchi,
#             "shape": "o",
#             "color": color,
#             "type": "molecule",
#             "dataset": dataset
#         }

#     # 3. Get reaction keys and which dataset they're in
#     def rxn_key(rxn):
#         reactants, products = rxn.split(">>")
#         rsmi = [s.strip() for s in reactants.split(".") if s.strip()]
#         psmi = [s.strip() for s in products.split(".") if s.strip()]
#         rkeys = [get_inchi(s) for s in rsmi]
#         pkeys = [get_inchi(s) for s in psmi]
#         reactants_key = '.'.join(sorted(rkeys))
#         products_key = '.'.join(sorted(pkeys))
#         return reactants_key + '>>' + products_key

#     rxn_map_1 = {rxn_key(r): r for r in reactions_1}
#     rxn_map_2 = {rxn_key(r): r for r in reactions_2}
#     all_rxn_keys = set(rxn_map_1) | set(rxn_map_2)

#     # 4. Build nodes and edges for reactions/intermediates
#     edges = []
#     for rxn_key_val in all_rxn_keys:
#         in_1 = rxn_key_val in rxn_map_1
#         in_2 = rxn_key_val in rxn_map_2
#         dataset = "both" if in_1 and in_2 else (dataset1_name if in_1 else dataset2_name)
#         edge_color = "green" if dataset == "both" else ("blue" if dataset == dataset1_name else "red")

#         rxn = rxn_map_1.get(rxn_key_val, rxn_map_2.get(rxn_key_val))
#         reactants, products = rxn.split(">>")
#         reactant_smis = [s.strip() for s in reactants.split(".") if s.strip()]
#         product_smis = [s.strip() for s in products.split(".") if s.strip()]
#         reactant_inchis = [get_inchi(s) for s in reactant_smis]
#         product_inchis = [get_inchi(s) for s in product_smis]

#         intermediate_reactant = None
#         intermediate_product = None

#         # --- Intermediate reactant node if needed
#         if len(reactant_inchis) > 1:
#             intR_id = f"{rxn_key_val}_intR"
#             if intR_id not in nodes:
#                 nodes[intR_id] = {
#                     "id": intR_id,
#                     "shape": "s",
#                     "color": edge_color,
#                     "type": "int_reactant",
#                     "dataset": dataset
#                 }
#             intermediate_reactant = intR_id
#             for r in reactant_inchis:
#                 if r:
#                     edges.append({
#                         "source": r,
#                         "target": intR_id,
#                         "color": edge_color,
#                         "dataset": dataset,
#                         "reaction_key": rxn_key_val
#                     })
#         # --- Intermediate product node if needed
#         if len(product_inchis) > 1:
#             intP_id = f"{rxn_key_val}_intP"
#             if intP_id not in nodes:
#                 nodes[intP_id] = {
#                     "id": intP_id,
#                     "shape": "^",
#                     "color": edge_color,
#                     "type": "int_product",
#                     "dataset": dataset
#                 }
#             intermediate_product = intP_id
#             for p in product_inchis:
#                 if p:
#                     edges.append({
#                         "source": intP_id,
#                         "target": p,
#                         "color": edge_color,
#                         "dataset": dataset,
#                         "reaction_key": rxn_key_val
#                     })
#         # --- Connect intermediates or direct
#         if intermediate_reactant and intermediate_product:
#             edges.append({
#                 "source": intermediate_reactant,
#                 "target": intermediate_product,
#                 "color": edge_color,
#                 "dataset": dataset,
#                 "reaction_key": rxn_key_val
#             })
#         elif intermediate_reactant and not intermediate_product:
#             p = product_inchis[0]
#             if p:
#                 edges.append({
#                     "source": intermediate_reactant,
#                     "target": p,
#                     "color": edge_color,
#                     "dataset": dataset,
#                     "reaction_key": rxn_key_val
#                 })
#         elif not intermediate_reactant and intermediate_product:
#             r = reactant_inchis[0]
#             if r:
#                 edges.append({
#                     "source": r,
#                     "target": intermediate_product,
#                     "color": edge_color,
#                     "dataset": dataset,
#                     "reaction_key": rxn_key_val
#                 })
#         else:
#             r, p = reactant_inchis[0], product_inchis[0]
#             if r and p:
#                 edges.append({
#                     "source": r,
#                     "target": p,
#                     "color": edge_color,
#                     "dataset": dataset,
#                     "reaction_key": rxn_key_val
#                 })
#     return list(nodes.values()), edges


# def build_graph_from_nodes_edges(nodes, edges):
#     """
#     Construct a NetworkX DiGraph from lists of node and edge attribute dicts.
#     Each node dict must have at least 'id'.
#     Each edge dict must have 'source' and 'target'.

#     Returns:
#         G (networkx.DiGraph): The constructed graph.
#     """
#     G = nx.DiGraph()
#     # Add nodes with attributes
#     for node in nodes:
#         node_id = node['id']
#         G.add_node(node_id, **{k: v for k, v in node.items() if k != 'id'})
#     # Add edges with attributes
#     for edge in edges:
#         src = edge['source']
#         tgt = edge['target']
#         attrs = {k: v for k, v in edge.items() if k not in ['source', 'target']}
#         G.add_edge(src, tgt, **attrs)
#     return G


# def draw_reaction_network_from_G(G, figsize=(8, 8), layout='spring'):
#     """
#     Draws a reaction network from a networkx.DiGraph G.
#     Assumes node attributes: shape, color, type, dataset
#     Assumes edge attributes: color, dataset, reaction_key

#     Returns the matplotlib figure object.
#     """
#     # Layout
#     if layout == 'kamada':
#         pos = nx.kamada_kawai_layout(G)
#     elif layout == 'circular':
#         pos = nx.circular_layout(G)
#     else:
#         pos = nx.spring_layout(G)

#     fig, ax = plt.subplots(figsize=figsize)
#     ax.axis('off')

    

#     # Group nodes by shape for drawing
#     shape_map = {'o': [], 's': [], '^': []}
#     color_map = {'o': [], 's': [], '^': []}

#     for n, d in G.nodes(data=True):
#         shape = d.get('shape', 'o')
#         color = d.get('color', 'black')
#         if shape not in shape_map:
#             shape_map[shape] = []
#             color_map[shape] = []
#         shape_map[shape].append(n)
#         color_map[shape].append(color)

#     for shape, nodelist in shape_map.items():
#         if nodelist:
#             nx.draw_networkx_nodes(
#                 G,
#                 pos,
#                 nodelist=nodelist,
#                 node_shape=shape,
#                 node_color=color_map[shape],
#                 node_size=200,
#                 ax=ax,
#                 alpha=0.8
#             )

#     # --- Group and draw edges by color ---
#     edge_color_map = {}
#     for u, v, d in G.edges(data=True):
#         color = d.get('color', 'black')
#         edge_color_map.setdefault(color, []).append((u, v))
#     for color, edgelist in edge_color_map.items():
#         nx.draw_networkx_edges(G, pos, ax=ax, edgelist=edgelist, edge_color=color, alpha=0.15, arrows=True, width=2)

#     # # Draw edges colored by dataset
#     # edge_colors = [edata["color"] for _, _, edata in G.edges(data=True)]
#     # nx.draw_networkx_edges(G, pos, ax=ax, edge_color=edge_colors, alpha=0.5, arrows=True)
#     # Optionally draw node labels (can be customized)
#     #nx.draw_networkx_labels(G, pos, ax=ax, font_size=8)

#     # Legend for node colors (dataset)
#     legend_elements = [
#         plt.Line2D([0], [1], color='green', lw=3, label='Both Datasets'),
#         plt.Line2D([0], [1], color='blue', lw=3, label='BNN Dataset'),
#         plt.Line2D([0], [1], color='red', lw=3, label='Reaxys Dataset'),
#         plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='black', markersize=10, label='Molecule'),
#         plt.Line2D([0], [0], marker='s', color='w', markerfacecolor='black', markersize=10, label='Multiple Reactants'),
#         plt.Line2D([0], [0], marker='^', color='w', markerfacecolor='black', markersize=10, label='Multiple Products')
#     ]
#     ax.legend(handles=legend_elements, loc='upper right', fontsize=10)

#     plt.tight_layout()
#     plt.show()

# def get_node_layer(G, n):
#     node = G.nodes[n]
#     shape = node.get('shape', 'o')
#     if shape == 's':
#         return 2  # intermediate reactant (square)
#     elif shape == '^':
#         return 3  # intermediate product (triangle)
#     else:  # molecule
#         in_deg = G.in_degree(n)
#         out_deg = G.out_degree(n)
#         if in_deg == 0 and out_deg > 0:
#             return 0  # only reactant (top)
#         elif in_deg > 0 and out_deg == 0:
#             return 4  # only product (bottom)
#         elif in_deg > 0 and out_deg > 0:
#             return 1  # both reactant and product (middle)
#         else:
#             return 1  # fallback

# def chunk_list(lst, n):
#     """Yield successive n-sized chunks from lst."""
#     for i in range(0, len(lst), n):
#         yield lst[i:i + n]

# def draw_reaction_network_layered_zigzag(G, figsize=(12, 8), max_per_layer=10, dy=0.08, y_amplitude=0.04):
#     # Assign nodes to layers
#     base_ys = {0: 1.0, 1: 0.8, 2: 0.6, 3: 0.4, 4: 0.2}
#     layer_nodes = {i: [] for i in base_ys}
#     for n in G.nodes:
#         layer = get_node_layer(G, n)
#         layer_nodes[layer].append(n)

#     # Calculate positions: split layers into sublayers if needed
#     pos = {}
#     for layer, nodes in layer_nodes.items():
#         if not nodes:
#             continue
#         chunks = list(chunk_list(nodes, max_per_layer))
#         for ci, chunk in enumerate(chunks):
#             base_y = base_ys[layer] - ci * dy
#             n_chunk = len(chunk)
#             x_coords = [i / (n_chunk - 1) if n_chunk > 1 else 0.5 for i in range(n_chunk)]
#             for xi, n in enumerate(chunk):
#                 if n_chunk > 1:
#                     sign = (-1) ** xi
#                     y = base_y + sign * y_amplitude * ((xi + 1) // 2)
#                 else:
#                     y = base_y
#                 pos[n] = (x_coords[xi], y)

#     # Fallback for any missing nodes
#     for n in G.nodes:
#         if n not in pos:
#             pos[n] = (0.5, 0.0)

#     fig, ax = plt.subplots(figsize=figsize)
#     ax.axis('off')

#     # Draw edges grouped by color
#     edge_color_map = {}
#     for u, v, d in G.edges(data=True):
#         color = d.get('color', 'black')
#         edge_color_map.setdefault(color, []).append((u, v))
#     for color, edgelist in edge_color_map.items():
#         nx.draw_networkx_edges(G, pos, ax=ax, edgelist=edgelist, edge_color=color, alpha=0.1, arrows=True, width=2)

#     # Draw nodes by shape
#     shape_map = {'o': [], 's': [], '^': []}
#     color_map = {'o': [], 's': [], '^': []}
#     for n, d in G.nodes(data=True):
#         shape = d.get('shape', 'o')
#         color = d.get('color', 'black')
#         shape_map.setdefault(shape, []).append(n)
#         color_map.setdefault(shape, []).append(color)
#     for shape, nodelist in shape_map.items():
#         if nodelist:
#             nx.draw_networkx_nodes(
#                 G, pos, nodelist=nodelist, node_shape=shape,
#                 node_color=color_map[shape], node_size=200, ax=ax, alpha=0.8
#             )

#     # Optionally draw molecule labels
#     # nx.draw_networkx_labels(G, pos, ax=ax, font_size=7)

#     # Legend: lines for datasets, markers for intermediates
#     legend_elements = [
#         plt.Line2D([0], [1], color='green', lw=3, label='Both Datasets'),
#         plt.Line2D([0], [1], color='blue', lw=3, label='BNN Dataset'),
#         plt.Line2D([0], [1], color='red', lw=3, label='Reaxys Dataset'),
#         plt.Line2D([0], [0], marker='s', color='w', markerfacecolor='grey', markersize=10, label='Connecting Reactants', linestyle='None'),
#         plt.Line2D([0], [0], marker='^', color='w', markerfacecolor='grey', markersize=10, label='Connecting Products', linestyle='None'),
#         plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='grey', markersize=10, label='Molecule', linestyle='None'),
#     ]
#     ax.legend(handles=legend_elements, loc='upper right', fontsize=10)

#     plt.tight_layout()
#     plt.show()


# def create_molecule_dictionary(reactions_1, reactions_2, name_1="dataset_1", name_2="dataset_2"):
#     """
#     Creates a dictionary of molecules from two sets of reactions, which can be accessed by the molecule inChi key.
#     The keys are the InChI keys which references a dictionary including the molecule object and the dataset it belongs to, or "both" if it is found in both datasets.
#     """
#     reaction_smiles_1 = [i for row in reactions_1 for part in row.split(">>") for i in part.split(".")]
#     reactions_1_inchis = {Chem.MolToInchiKey(Chem.MolFromSmiles(smi.strip())) for smi in reaction_smiles_1}
#     reaction_smiles_2 = [i for row in reactions_2 for part in row.split(">>") for i in part.split(".")]
#     reactions_2_inchis = {Chem.MolToInchiKey(Chem.MolFromSmiles(smi.strip())) for smi in reaction_smiles_2}
#     # iterate through the reactions and create a dictionary of inChi keys
#     inchi_dict = {}
#     for reaction in reactions_1 + reactions_2:
#         # Split the reaction into reactants and products
#         reactants, products = reaction.split('>>')
#         # Get the individual molecules from reactants and products
#         smiles = reactants.split('.') + products.split('.')
#         for smi in smiles:
#             # Create a molecule object from the SMILES string
#             mol = Chem.MolFromSmiles(smi.strip())
#             if mol:
#                 inchi_key = Chem.MolToInchiKey(mol)
#                 if inchi_key not in inchi_dict:
#                     inchi_dict[inchi_key] = {
#                         "mol": mol,
#                         "dataset": None,
#                     }
#                     # Add the dataset to the dictionary
#                     if inchi_key in reactions_1_inchis and inchi_key in reactions_2_inchis:
#                         inchi_dict[inchi_key]["dataset"] = "both"
#                     elif inchi_key in reactions_1_inchis and inchi_key not in reactions_2_inchis:
#                         inchi_dict[inchi_key]["dataset"] = name_1
#                     elif inchi_key in reactions_2_inchis and inchi_key not in reactions_1_inchis:
#                         inchi_dict[inchi_key]["dataset"] = name_2
#                 else:
#                     # If the molecule is already in the dictionary, we do not need to add it again
#                     continue      
#     return inchi_dict

# def create_reaction_key(reaction):
#     """
#     Create a unique key for a reaction based on its reactants and products.
#     """
#     reactants, products = reaction.split('>>')
#     # Get the individual molecules from reactants and products
#     reactant_smiles = [s.strip() for s in reactants.split('.')]
#     product_smiles = [s.strip() for s in products.split('.')]
#     reactant_inchis = [Chem.MolToInchiKey(Chem.MolFromSmiles(smi.strip())) for smi in reactant_smiles]
#     product_inchis = [Chem.MolToInchiKey(Chem.MolFromSmiles(smi.strip())) for smi in product_smiles]
#     # Create a unique key for the reaction using the InChI keys of reactants and products
#     reactants_key = '.'.join(sorted(reactant_inchis)) # sort to ensure consistent ordering
#     products_key = '.'.join(sorted(product_inchis)) # sort to ensure consistent ordering
#     reaction_key = reactants_key + '>>' + products_key
#     return reaction_key

# def create_reaction_network(reactions_1, reactions_2):
#     """
#     Creates a reaction network graph from two sets of reactions.
#     The nodes are the molecules (circles) and the edges are the reactions between them.
#     Intermediate nodes (squares for reactants, triangles for products) are added if there are multiple reactants or products.
#     """
#     G = nx.DiGraph()
#     inchi_dict = create_molecule_dictionary(reactions_1, reactions_2)
    
#     # Add molecule nodes
#     for index, (inchi_key, data) in enumerate(inchi_dict.items()):
#         G.add_node(
#             inchi_key,
#             mol=data["mol"],
#             dataset=data["dataset"],
#             index=index,
#             shape='o',   # circle for molecules
#             type='molecule'
#         )

#     # Build all reaction keys (unique by InChI)
#     all_reactions = reactions_1 + reactions_2
#     reaction_keys = {create_reaction_key(r): r for r in all_reactions}
#     # Map for dataset type by reaction key
#     reaction_datasets = {}
#     rkeys_1 = {create_reaction_key(r) for r in reactions_1}
#     rkeys_2 = {create_reaction_key(r) for r in reactions_2}
#     for k in reaction_keys:
#         if k in rkeys_1 and k in rkeys_2:
#             reaction_datasets[k] = 'both'
#         elif k in rkeys_1:
#             reaction_datasets[k] = 'dataset_1'
#         elif k in rkeys_2:
#             reaction_datasets[k] = 'dataset_2'
#         else:
#             reaction_datasets[k] = 'other'

#     node_counter = len(G.nodes)
#     for reaction_key, reaction in reaction_keys.items():
#         reactants, products = reaction.split('>>')
#         reactant_smiles = [s.strip() for s in reactants.split('.') if s.strip()]
#         product_smiles = [s.strip() for s in products.split('.') if s.strip()]
#         reactant_inchis = [Chem.MolToInchiKey(Chem.MolFromSmiles(smi)) for smi in reactant_smiles]
#         product_inchis = [Chem.MolToInchiKey(Chem.MolFromSmiles(smi)) for smi in product_smiles]
#         dataset = reaction_datasets[reaction_key]

#         # Determine intermediate nodes
#         intermediate_reactant = None
#         intermediate_product = None

#         # Add intermediate reactant node if multiple reactants
#         if len(reactant_inchis) > 1:
#             intermediate_reactant = f"{reaction_key}_intR"
#             G.add_node(
#                 intermediate_reactant,
#                 mol=None,
#                 dataset=dataset,
#                 index=node_counter,
#                 shape='s',   # square
#                 type='int_reactant'
#             )
#             node_counter += 1
#             for r in reactant_inchis:
#                 if r in G:
#                     G.add_edge(r, intermediate_reactant, reaction=reaction_key, dataset=dataset)

#         # Add intermediate product node if multiple products
#         if len(product_inchis) > 1:
#             intermediate_product = f"{reaction_key}_intP"
#             G.add_node(
#                 intermediate_product,
#                 mol=None,
#                 dataset=dataset,
#                 index=node_counter,
#                 shape='^',   # triangle
#                 type='int_product'
#             )
#             node_counter += 1
#             for p in product_inchis:
#                 if p in G:
#                     G.add_edge(intermediate_product, p, reaction=reaction_key, dataset=dataset)

#         # Connect intermediates or direct edges
#         if intermediate_reactant and intermediate_product:
#             # Multiple reactants & products: intR → intP
#             G.add_edge(intermediate_reactant, intermediate_product, reaction=reaction_key, dataset=dataset)
#         elif intermediate_reactant and not intermediate_product:
#             # Multiple reactants, single product: intR → product
#             p = product_inchis[0]
#             if p in G:
#                 G.add_edge(intermediate_reactant, p, reaction=reaction_key, dataset=dataset)
#         elif not intermediate_reactant and intermediate_product:
#             # Single reactant, multiple products: reactant → intP
#             r = reactant_inchis[0]
#             if r in G:
#                 G.add_edge(r, intermediate_product, reaction=reaction_key, dataset=dataset)
#         else:
#             # Single reactant, single product: reactant → product
#             r, p = reactant_inchis[0], product_inchis[0]
#             if r in G and p in G:
#                 G.add_edge(r, p, reaction=reaction_key, dataset=dataset)

#     return G




# def hide_edges(G, hide_both=True, hide_dataset_1=False, hide_dataset_2=False):
#     """
#     Hides edges based on the dataset they belong to.
#     """
#     edges_to_remove = []
#     for u, v, data in G.edges(data=True):
#         if (hide_both and G.nodes[u]['dataset'] == 'both' and G.nodes[v]['dataset'] == 'both') or \
#            (hide_dataset_1 and G.nodes[u]['dataset'] == 'dataset_1' and G.nodes[v]['dataset'] == 'dataset_1') or \
#            (hide_dataset_2 and G.nodes[u]['dataset'] == 'dataset_2' and G.nodes[v]['dataset'] == 'dataset_2'):
#             edges_to_remove.append((u, v))
#     G.remove_edges_from(edges_to_remove)
#     return G

# def hide_nodes(G, hide_both=True, hide_dataset_1=False, hide_dataset_2=False):
#     """
#     Hides nodes based on the dataset they belong to.
#     """
#     nodes_to_remove = []
#     for node, data in G.nodes(data=True):
#         if (hide_both and data['dataset'] == 'both') or \
#            (hide_dataset_1 and data['dataset'] == 'dataset_1') or \
#            (hide_dataset_2 and data['dataset'] == 'dataset_2'):
#             nodes_to_remove.append(node)
#     G.remove_nodes_from(nodes_to_remove)
#     return G

# def show_png(mol_bytes):
#     """Convert RDKit Cairo PNG bytes to a numpy array for imshow."""
#     bio = BytesIO(mol_bytes)
#     img = Image.open(bio)
#     return np.array(img)

# def highlight_mol(smi,color):
    
#     mol = Chem.MolFromSmiles(smi)
    
#     if color == 'darkred':
#         rgba = (0.55, 0.0, 0.0, 0.3)       
#     elif color == 'red':
#         rgba = (1.0, 0.0, 0.0, 0.2)    
#     elif color == 'orange':
#         rgba = (1.0, 0.65, 0.0, 0.2)        
#     elif color == 'yellow':
#         rgba = (1.0, 1.0, 0.0, 0.4)        
#     elif color == 'green':
#         rgba = (0.0, 0.50, 0.0, 0.15)        
#     elif color == 'lightblue':
#         rgba = (0.68, 0.85, 0.90, 0.5)   
#     elif color == 'blue':
#         rgba = (0.0, 0.0, 1.0, 0.1)    
#     else: # no color
#         rgba = (1,1,1,1) 
             
#     atoms = []
#     for a in mol.GetAtoms():
#         atoms.append(a.GetIdx())
    
#     bonds = []
#     for bond in mol.GetBonds():
#         aid1 = atoms[bond.GetBeginAtomIdx()]
#         aid2 = atoms[bond.GetEndAtomIdx()]
#         bonds.append(mol.GetBondBetweenAtoms(aid1,aid2).GetIdx())

#     drawer = rdMolDraw2D.MolDraw2DCairo(150,150)
#     Draw.SetACS1996Mode(drawer.drawOptions(),Draw.MeanBondLength(mol))
#     drawer.drawOptions().fillHighlights=True
#     drawer.drawOptions().setHighlightColour((rgba))
#     drawer.drawOptions().highlightBondWidthMultiplier=15
#     # drawer.drawOptions().legendFontSize=60
#     drawer.drawOptions().clearBackground = False
#     rdMolDraw2D.PrepareAndDrawMolecule(drawer, mol, highlightAtoms=atoms, highlightBonds=bonds)
    
#     mol_png = drawer.GetDrawingText()
#     return mol_png

# def draw_reaction_network(G, figsize=(8, 8), layout='spring', draw_mols=False):
#     """
#     Draws the reaction network graph with different node shapes and colors for dataset.
#     Returns the matplotlib figure object.
#     """
#     # Layout
#     if layout == 'kamada':
#         pos = nx.kamada_kawai_layout(G)
#     elif layout == 'circular':
#         pos = nx.circular_layout(G)
#     else:
#         pos = nx.spring_layout(G)

#     fig, ax = plt.subplots(figsize=figsize)
#     ax.axis('off')
#     nx.draw_networkx_edges(G, pos, ax=ax, alpha=0.5, arrows=True)

#     if draw_mols:
#         # Transform from data coordinates (scaled between xlim and ylim) to display coordinates
#         tr_figure = ax.transData.transform
#         # Transform from display to figure coordinates
#         tr_axes = fig.transFigure.inverted().transform

#         # Select the size of the image (relative to the X axis)
#         struct_size = (ax.get_xlim()[1] - ax.get_xlim()[0]) * 0.1
#         struct_center = struct_size / 2.0
#         # Draw molecule nodes with different colours based on the dataset
#         for inChi,values in G.nodes.items():
#             mol = values['mol']
#             smi = Chem.MolToSmiles(mol)
#             dataset = values['dataset']
#             if dataset == 'both':
#                 colour = 'green'
#             elif dataset == 'dataset_1':
#                 colour = 'blue'
#             elif dataset == 'dataset_2':
#                 colour = 'red'
#             else:
#                 colour = 'yellow'
#             xf, yf = tr_figure(pos[inChi])
#             xa, ya = tr_axes((xf, yf))
#             # get overlapped axes and plot structure
#             a = plt.axes([xa - struct_center, ya - struct_center, struct_size, struct_size])
#             a.imshow(show_png(highlight_mol(smi, colour)))
#             a.axis("off")
#     else:
#         # Group nodes by shape
#         shape_nodes = {'o': [], 's': [], '^': []}
#         shape_colors = {'o': [], 's': [], '^': []}
#         for node, data in G.nodes(data=True):
#             shape = data.get('shape', 'o')
#             dataset = data.get('dataset', 'other')
#             if dataset == 'both':
#                 color = 'green'
#             elif dataset == 'dataset_1':
#                 color = 'blue'
#             elif dataset == 'dataset_2':
#                 color = 'red'
#             else:
#                 color = 'black'
#             if shape in shape_nodes:
#                 shape_nodes[shape].append(node)
#                 shape_colors[shape].append(color)
#             else:
#                 shape_nodes[shape] = [node]
#                 shape_colors[shape] = [color]

#         # Draw each shape separately
#         for shape, nodelist in shape_nodes.items():
#             if nodelist:
#                 nx.draw_networkx_nodes(
#                     G,
#                     pos,
#                     nodelist=nodelist,
#                     node_shape=shape,
#                     node_color=shape_colors[shape],
#                     node_size=200,
#                     ax=ax,
#                     alpha=0.7
#                 )
#         #nx.draw_networkx_labels(G, pos, ax=ax, font_size=8)

#     # Legend for dataset colors (not shape)
#     legend_elements = [
#         plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='green', markersize=10, label='Both Datasets'),
#         plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='blue', markersize=10, label='Dataset 1'),
#         plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='red', markersize=10, label='Dataset 2'),
#         plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='black', markersize=10, label='Other'),
#     ]
#     ax.legend(handles=legend_elements, loc='upper right', fontsize=10)

#     plt.tight_layout()
#     plt.show()
