"""
Utility: Check molecule name presence in PubChem and CIRpy.

Features:
- Accept a single name or a list of names.
- Prefer python libs (pubchempy, cirpy) when available.
- Fallback to HTTP APIs when libraries are not installed.
- Simple, dependency-light (only stdlib for fallbacks).

Returns, per name: { 'pubchem': bool, 'cirpy': bool, 'exists_any': bool }
Also supports:
- name_to_smiles(name) -> list[str]
- smiles_to_names(smiles) -> list[str]
- batch resolvers for both directions
"""
from __future__ import annotations

from typing import Iterable, Union, Dict, Any, List
from functools import lru_cache
from urllib import request, parse, error
import requests  # for CAS Common Chemistry API

# RDKit (optional but recommended)
try:
    from rdkit import Chem  # type: ignore
    try:
        from rdkit.Chem import inchi as rd_inchi  # type: ignore
    except Exception:  # pragma: no cover
        rd_inchi = None  # type: ignore
except Exception:  # pragma: no cover
    Chem = None  # type: ignore
    rd_inchi = None  # type: ignore

# Optional libs
try:  # type: ignore
    import pubchempy as pcp  # noqa: F401
except Exception:  # pragma: no cover
    pcp = None  # type: ignore

try:  # type: ignore
    import cirpy  # noqa: F401
except Exception:  # pragma: no cover
    cirpy = None  # type: ignore


class MoleculePresenceChecker:
    """
    Check if molecule names exist in PubChem and CIRpy.

    Example:
        checker = MoleculePresenceChecker()
        result = checker.check(["aspirin", "not_a_real_chemical_123"])
        # {
        #   'aspirin': {'pubchem': True, 'cirpy': True, 'exists_any': True},
        #   'not_a_real_chemical_123': {'pubchem': False, 'cirpy': False, 'exists_any': False}
        # }
    """

    def __init__(self, timeout: float = 12.0, use_http_fallback: bool = True) -> None:
        self.timeout = float(timeout)
        self.use_http_fallback = bool(use_http_fallback)
        
    # ------------------------ Helpers (internal) ------------------
    def _as_list(self, value: Any) -> List[str]:
        if value is None:
            return []
        if isinstance(value, (list, tuple, set)):
            return [str(v) for v in value if v is not None]
        return [str(value)]

    # ------------------------ PubChem checks ------------------------
    @lru_cache(maxsize=4096)
    def check_pubchem(self, name: str) -> bool:
        name = (name or "").strip()
        if not name:
            return False
        # Prefer pubchempy, if available
        if pcp is not None:
            try:
                # Fast existence check: get at most one compound by name
                comps = pcp.get_compounds(name, "name")  # type: ignore[attr-defined]
                return bool(comps)
            except Exception:
                # Fallback to HTTP if allowed
                if not self.use_http_fallback:
                    return False
        # HTTP fallback via PUG REST
        if self.use_http_fallback:
            try:
                qname = parse.quote(name)
                url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{qname}/cids/TXT"
                txt = self._http_get_text(url, timeout=self.timeout)
                # Non-empty response with one or more CIDs means it exists
                return bool(txt and txt.strip())
            except Exception:
                return False
        return False

    # ------------------------ CIRpy checks -------------------------
    @lru_cache(maxsize=4096)
    def check_cirpy(self, name: str) -> bool:
        name = (name or "").strip()
        if not name:
            return False
        if cirpy is not None:
            try:
                # Try resolving to SMILES first, then InChIKey
                smi = cirpy.resolve(name, "smiles")  # type: ignore[attr-defined]
                if smi:
                    return True
                ik = cirpy.resolve(name, "inchikey")  # type: ignore[attr-defined]
                if ik:
                    return True
            except Exception:
                if not self.use_http_fallback:
                    return False
        # HTTP fallback to CACTUS (analogous to cirpy backend)
        if self.use_http_fallback:
            try:
                qname = parse.quote(name)
                url = f"https://cactus.nci.nih.gov/chemical/structure/{qname}/smiles"
                txt = self._http_get_text(url, timeout=self.timeout)
                return bool(txt and txt.strip() and "Not Found" not in txt)
            except Exception:
                return False
        return False

    # ------------------------ Name -> SMILES -----------------------
    @lru_cache(maxsize=4096)
    def name_to_smiles(self, name: str, prefer_isomeric: bool = True) -> List[str]:
        """
        Resolve a chemical name to one or more SMILES strings.
        Order is not guaranteed; duplicates removed.

        Strategy:
        1) cirpy.resolve(name, 'smiles') if available
        2) pubchempy.get_compounds(name, 'name') -> canonical_smiles
        3) HTTP fallbacks:
           - CACTUS: /chemical/structure/{name}/smiles
           - PubChem PUG REST properties (CanonicalSMILES/IsomericSMILES)
        """
        name = (name or "").strip()
        if not name:
            return []

        smiles_set = set()
        # cirpy first
        if cirpy is not None:
            try:
                smi = cirpy.resolve(name, "smiles")  # type: ignore[attr-defined]
                for s in self._as_list(smi):
                    if s:
                        smiles_set.add(s.strip())
            except Exception:
                pass

        # pubchempy compounds
        if pcp is not None:
            try:
                comps = pcp.get_compounds(name, "name")  # type: ignore[attr-defined]
                for c in comps or []:
                    try:
                        cs = getattr(c, "isomeric_smiles", None) if prefer_isomeric else getattr(c, "canonical_smiles", None)
                        if not cs:
                            cs = getattr(c, "canonical_smiles", None)
                        if cs:
                            smiles_set.add(cs.strip())
                    except Exception:
                        continue
            except Exception:
                pass

        # HTTP fallbacks
        if self.use_http_fallback:
            # CACTUS
            try:
                qname = parse.quote(name)
                url = f"https://cactus.nci.nih.gov/chemical/structure/{qname}/smiles"
                txt = self._http_get_text(url, timeout=self.timeout)
                if txt and txt.strip() and "Not Found" not in txt:
                    smiles_set.add(txt.strip())
            except Exception:
                pass
            # PubChem PUG REST properties
            try:
                prop = "IsomericSMILES" if prefer_isomeric else "CanonicalSMILES"
                qname = parse.quote(name)
                url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{qname}/property/{prop}/TXT"
                txt = self._http_get_text(url, timeout=self.timeout)
                if txt and txt.strip():
                    for line in txt.splitlines():
                        s = line.strip()
                        if s:
                            smiles_set.add(s)
            except Exception:
                pass

        return sorted(smiles_set)

    # ------------------------ SMILES -> Names ----------------------
    @lru_cache(maxsize=4096)
    def smiles_to_names(self, smiles: str, max_names: int | None = 100) -> List[str]:
        """
        Resolve a SMILES string to one or more names (synonyms).
        Returns a de-duplicated list (empty if no mapping).
        """
        smiles = (smiles or "").strip()
        if not smiles:
            return []

        names: List[str] = []
        seen = set()

        # pubchempy first (rich synonyms)
        if pcp is not None:
            try:
                comps = pcp.get_compounds(smiles, "smiles")  # type: ignore[attr-defined]
                for c in comps or []:
                    # Try attached synonyms first
                    syns = []
                    try:
                        syns = getattr(c, "synonyms", None) or []
                    except Exception:
                        syns = []
                    # If not present, fetch via CID
                    if not syns:
                        try:
                            cid = getattr(c, "cid", None)
                            if cid:
                                # pcp.get_synonyms returns list of dicts; extract synonyms
                                data = pcp.get_synonyms(int(cid))  # type: ignore[attr-defined]
                                for entry in data or []:
                                    for s in entry.get("Synonym", []) or []:
                                        if s and s not in seen:
                                            names.append(s)
                                            seen.add(s)
                        except Exception:
                            pass
                    else:
                        for s in syns:
                            if s and s not in seen:
                                names.append(s)
                                seen.add(s)
            except Exception:
                pass

        # cirpy (may support 'names' via CACTUS backend)
        if cirpy is not None:
            try:
                cir_names = cirpy.resolve(smiles, "names")  # type: ignore[attr-defined]
                for s in self._as_list(cir_names):
                    s2 = s.strip()
                    if s2 and s2 not in seen:
                        names.append(s2)
                        seen.add(s2)
            except Exception:
                pass

        # HTTP fallbacks
        if self.use_http_fallback:
            # CACTUS names
            try:
                qsmiles = parse.quote(smiles)
                url = f"https://cactus.nci.nih.gov/chemical/structure/{qsmiles}/names"
                txt = self._http_get_text(url, timeout=self.timeout)
                if txt and txt.strip() and "Not Found" not in txt:
                    for line in txt.splitlines():
                        s = line.strip()
                        if s and s not in seen:
                            names.append(s)
                            seen.add(s)
            except Exception:
                pass
            # PubChem synonyms by SMILES
            try:
                qsmiles = parse.quote(smiles)
                url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{qsmiles}/synonyms/TXT"
                txt = self._http_get_text(url, timeout=self.timeout)
                if txt and txt.strip():
                    for line in txt.splitlines():
                        s = line.strip()
                        if s and s not in seen:
                            names.append(s)
                            seen.add(s)
            except Exception:
                pass

        if max_names is not None:
            return names[: int(max_names)]
        return names

    # ------------------------ Batch resolvers ----------------------
    def names_to_smiles(self, names: Union[str, Iterable[str]]) -> Dict[str, List[str]]:
        """Batch: {name -> [smiles...]}."""
        items = [names] if isinstance(names, str) else list(names)
        return {n: self.name_to_smiles(n) for n in items}

    def smiles_list_to_names(self, smiles_list: Union[str, Iterable[str]], max_names: int | None = 100) -> Dict[str, List[str]]:
        """Batch: {smiles -> [names...]}."""
        items = [smiles_list] if isinstance(smiles_list, str) else list(smiles_list)
        return {s: self.smiles_to_names(s, max_names=max_names) for s in items}

    # ------------------------ Batch API ----------------------------
    def check(self, names: Union[str, Iterable[str]]) -> Dict[str, Dict[str, bool]]:
        """
        Check one or many names across PubChem and CIRpy.

        Returns a dict mapping each input name to:
          {'pubchem': bool, 'cirpy': bool, 'exists_any': bool}
        """
        if isinstance(names, str):
            items = [names]
        else:
            items = list(names)

        out: Dict[str, Dict[str, bool]] = {}
        for name in items:
            exists_pubchem = self.check_pubchem(name)
            exists_cirpy = self.check_cirpy(name)
            out[name] = {
                "pubchem": bool(exists_pubchem),
                "cirpy": bool(exists_cirpy),
                "exists_any": bool(exists_pubchem or exists_cirpy),
                "Has_CAS": self.has_cas(name),
            }
        return out

    # ------------------------ Helpers ------------------------------
    def _http_get_text(self, url: str, timeout: float) -> str | None:
        req = request.Request(url, method="GET")
        try:
            with request.urlopen(req, timeout=timeout) as resp:
                # 2xx indicates OK; read small text body
                if 200 <= resp.status < 300:
                    data = resp.read()
                    return data.decode("utf-8", errors="replace")
                return None
        except error.HTTPError as e:
            # 404/400 -> not found
            if e.code in (400, 404):
                return None
            # Other errors -> treat as not found / transient failure
            return None

    # ------------------------ Convenience --------------------------
    def exists(self, name: str) -> bool:
        """Return True if name found in PubChem or CIRpy; False otherwise."""
        res = self.check(name).get(name, {})
        return bool(res.get("exists_any", False))

    # ------------------------ CAS via CAS Common Chemistry --------
    @lru_cache(maxsize=4096)
    def _cas_registry_from_smiles(self, smiles: str) -> str | None:
        """Return CAS RN for a SMILES using the CAS Common Chemistry API, or None.

        Follows the approach requested: SMILES -> RDKit mol -> InChI -> CAS /search API.
        """
        smiles = (smiles or "").strip()
        if not smiles or Chem is None:
            return None
        try:
            mol = Chem.MolFromSmiles(smiles)
        except Exception:
            mol = None
        if not mol:
            return None

        inchi_str: str | None = None
        try:
            if rd_inchi is not None and hasattr(rd_inchi, "MolToInchi"):
                val = rd_inchi.MolToInchi(mol)  # may return str or (str, ...)
                inchi_str = val[0] if isinstance(val, (list, tuple)) else val
            else:
                # Fallback to older API
                val = getattr(Chem, "rdinchi", None)
                if val is not None and hasattr(val, "MolToInchi"):
                    inchi_str = val.MolToInchi(mol)[0]
        except Exception:
            inchi_str = None
        if not inchi_str:
            return None

        try:
            base_url = "https://commonchemistry.cas.org/api"
            endpoint = "/search"
            params = {"q": inchi_str, "size": "1"}
            resp = requests.get(base_url + endpoint, params=params, timeout=self.timeout)
            if resp.status_code == 200:
                data = resp.json()
                try:
                    if int(data.get("count", 0)) > 0:
                        return data["results"][0].get("rn")
                except Exception:
                    return None
            return None
        except Exception:
            return None

    @lru_cache(maxsize=4096)
    def has_cas_for_smiles(self, smiles: str) -> bool:
        rn = self._cas_registry_from_smiles(smiles)
        return bool(rn)

    @lru_cache(maxsize=4096)
    def has_cas(self, name: str, smiles_limit: int = 3) -> bool:
        """Return True if any SMILES resolved for this name has a CAS RN via CAS API."""
        # Resolve name to SMILES first (uses cirpy/pubchempy/http fallbacks)
        smiles_list = self.name_to_smiles(name)
        for smi in smiles_list[: int(smiles_limit)]:
            if self.has_cas_for_smiles(smi):
                return True
        return False

