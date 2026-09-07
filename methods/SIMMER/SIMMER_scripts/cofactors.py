"""
Standard metabolic cofactor list for reaction-SMILES stripping.

Matching is done on the RDKit "skeleton" InChIKey block (first 14 characters),
which is invariant to protonation/charge state and stereochemistry. This lets
a single reference SMILES (e.g. neutral ATP) match ATP written with any mix
of [O-]/[OH] or isomeric annotations, which is how the same cofactor shows up
differently across Rhea/MetaCyc/drug reaction SMILES.
"""

from rdkit import Chem
from rdkit.Chem import inchi

# Reference SMILES for common cofactors / cosubstrates and small
# inorganic leaving groups that are not the "reaction of interest" chemistry.
COFACTOR_SMILES = {
    # Water / simple inorganics
    "water": "O",
    "hydron": "[H+]",
    "hydroxide": "[OH-]",
    "co2": "O=C=O",
    "bicarbonate": "OC(=O)[O-]",
    "oxygen": "O=O",
    "hydrogen_peroxide": "OO",
    "superoxide": "[O-][O]",
    "ammonia": "N",
    "ammonium": "[NH4+]",
    "hydrogen": "[H][H]",
    "phosphate": "OP(=O)(O)O",
    "diphosphate": "OP(=O)(O)OP(=O)(O)O",
    "sulfate": "OS(=O)(=O)O",
    "sulfite": "OS(=O)O",
    "chloride": "[Cl-]",
    "nitric_oxide": "[N]=O",
    "nitrite": "[O-]N=O",
    "nitrate": "[O-][N+](=O)[O-]",
    # Adenosine phosphates
    "amp": "Nc1ncnc2c1ncn2C1OC(COP(=O)(O)O)C(O)C1O",
    "adp": "Nc1ncnc2c1ncn2C1OC(COP(=O)(O)OP(=O)(O)O)C(O)C1O",
    "atp": "Nc1ncnc2c1ncn2C1OC(COP(=O)(O)OP(=O)(O)OP(=O)(O)O)C(O)C1O",
    # Guanosine phosphates
    "gmp": "Nc1nc2c(ncn2C2OC(COP(=O)(O)O)C(O)C2O)c(=O)[nH]1",
    "gdp": "Nc1nc2c(ncn2C2OC(COP(=O)(O)OP(=O)(O)O)C(O)C2O)c(=O)[nH]1",
    "gtp": "Nc1nc2c(ncn2C2OC(COP(=O)(O)OP(=O)(O)OP(=O)(O)O)C(O)C2O)c(=O)[nH]1",
    # Uridine / cytidine phosphates
    "udp": "O=c1ccn(C2OC(COP(=O)(O)OP(=O)(O)O)C(O)C2O)c(=O)[nH]1",
    "utp": "O=c1ccn(C2OC(COP(=O)(O)OP(=O)(O)OP(=O)(O)O)C(O)C2O)c(=O)[nH]1",
    "cdp": "Nc1ccn(C2OC(COP(=O)(O)OP(=O)(O)O)C(O)C2O)c(=O)n1",
    "ctp": "Nc1ccn(C2OC(COP(=O)(O)OP(=O)(O)OP(=O)(O)O)C(O)C2O)c(=O)n1",
    # Redox cofactors
    "nad+": "NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1",
    "nadh": "NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1",
    "nadp+": "NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)c1",
    "nadph": "NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)C=CC1",
    "fad": "Cc1cc2nc3c(=O)[nH]c(=O)nc-3n(CC(O)C(O)C(O)COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)c2cc1C",
    "fadh2": "Cc1cc2c(cc1C)N(CC(O)C(O)C(O)COP(=O)(O)OP(=O)(O)OCC1OC(n3cnc4c(N)ncnc43)C(O)C1O)c1[nH]c(=O)[nH]c(=O)c1N2",
    "fmn": "Cc1cc2nc3c(=O)[nH]c(=O)nc-3n(CC(O)C(O)C(O)COP(=O)(O)O)c2cc1C",
    "fmnh2": "Cc1cc2c(cc1C)Nc1c(=O)[nH]c(=O)[nH]c1N2CC(O)C(O)C(O)COP(=O)(O)O",
    "ubiquinone": "O=C1C(OC)=C(OC)C(=O)C(C)=C1CC=C(C)C",
    "ubiquinol": "Oc1c(OC)c(OC)c(CC=C(C)C)cc1C",
    "glutathione_reduced": "NC(CCC(=O)NC(CS)C(=O)NCC(=O)O)C(=O)O",
    "glutathione_oxidized": "NC(CCC(=O)NC(CSSCC(NC(=O)CCC(N)C(=O)O)C(=O)NCC(=O)O)C(=O)NCC(=O)O)C(=O)O",
    # Coenzyme A / acyl carriers
    "coa": "CC(C)(COP(=O)(O)OP(=O)(O)OCC1OC(n2cnc3c(N)ncnc32)C(OP(=O)(O)O)C1O)C(O)C(=O)NCCC(=O)NCCS",
    "acetyl_coa": "CC(=O)SCCNC(=O)CCNC(=O)C(O)C(C)(C)COP(=O)(O)OP(=O)(O)OCC1OC(n2cnc3c(N)ncnc32)C(OP(=O)(O)O)C1O",
    # One-carbon / methyl transfer
    "sam": "C[S+](CCC(N)C(=O)O)CC1OC(n2cnc3c(N)ncnc32)C(O)C1O",
    "sah": "NC(CCSCC1OC(n2cnc3c(N)ncnc32)C(O)C1O)C(=O)O",
    "thf": "Nc1nc2NCC(CNc3ccc(cc3)C(=O)NC(CCC(=O)O)C(=O)O)Nc2c(=O)[nH]1",
    "methyl_thf": "Nc1nc2N(C)C(CNc3ccc(cc3)C(=O)NC(CCC(=O)O)C(=O)O)CNc2c(=O)[nH]1",
    "biotin": "O=C1NC2C(NC1=O)SCC2CCCCC(=O)O",
    # Pyridoxal phosphate
    "plp": "Cc1ncc(COP(=O)(O)O)c(C=O)c1O",
    "pmp": "Cc1ncc(COP(=O)(O)O)c(CN)c1O",
    # Thiamine pyrophosphate
    "tpp": "Cc1ncc(C[n+]2csc(CCOP(=O)(O)OP(=O)(O)O)c2C)c(N)n1",
    # Metal-ion redox/cofactor centers common in oxidoreductase (class 1) reactions
    "fe0": "[Fe]",
    "fe2": "[Fe+2]",
    "fe3": "[Fe+3]",
    "sulfide": "[SH-]",
    "cu1": "[Cu+]",
    "cu2": "[Cu+2]",
    "zn2": "[Zn+2]",
    "mg2": "[Mg+2]",
    "mn2": "[Mn+2]",
    "ni2": "[Ni+2]",
    "co2_ion": "[Co+2]",
    "mo": "[Mo]",
}


def _skeleton_key(smiles: str) -> str | None:
    """First block of the InChIKey: connectivity only, ignores charge/tautomer/stereo."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    try:
        key = inchi.MolToInchiKey(mol)
    except Exception:
        return None
    if not key or key.startswith("MOSTPRZ"):  # RDKit's inchi failure sentinel
        return None
    return key.split("-")[0]


def build_cofactor_keyset() -> set[str]:
    keys = set()
    for name, smi in COFACTOR_SMILES.items():
        key = _skeleton_key(smi)
        if key is None:
            raise ValueError(f"Could not parse cofactor reference SMILES for {name!r}: {smi}")
        keys.add(key)
    return keys


COFACTOR_KEYS = build_cofactor_keyset()


def strip_cofactor_components(side_smiles: str) -> tuple[str, list[str]]:
    """
    Given one side of a reaction ('.'-joined component SMILES), drop any
    component whose skeleton InChIKey matches a known cofactor.
    Returns (remaining_side_smiles, dropped_component_smiles).
    """
    components = side_smiles.split(".")
    kept, dropped = [], []
    for comp in components:
        key = _skeleton_key(comp)
        if key is not None and key in COFACTOR_KEYS:
            dropped.append(comp)
        else:
            kept.append(comp)
    return ".".join(kept), dropped
