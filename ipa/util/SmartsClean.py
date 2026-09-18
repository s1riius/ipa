# ------------------------------------------------------------
# RDKit functions to exclude bad RDKit protonation states
# ------------------------------------------------------------

from typing import Iterable,List,Tuple,Optional
from rdkit import Chem

BAD_TAUTOMER_SMARTS = {
    #amide/enamide iminol forms: R-C(OH)=N-R'
    "amide_iminol": "[CX3]([OX2H])=[NX2]",
    #enol-amino motif: C=C(OH)-N
    "enol_amino": "[CX3]=[CX3]([OX2H])-[NX3;H0,H1,H2]",
    #protonated amide nitrogen: R-C(=O)-[NH2+]-R / R-C(=O)-[N+](H)(H)R
    "protonated_amide": "[CX3](=[OX1])-[NX4+]",
    #protonated aniline-like nitrogen attached to aromatic carbon
    "protonated_aniline": "[NX4+;H1,H2,H3]-c"
}

def clean_rdkit_tautomer_smiles(smiles_list: Iterable[str]) -> Tuple[List[str], List[Tuple[str, str]]]:
    """Remove chemically unwanted RDKit-generated tautomer/protonation states.
       Returns:
           kept: list of accepted SMILES
           rejected: list of (SMILES, reason)"""
    patterns = {name: Chem.MolFromSmarts(smarts)
                for name,smarts in BAD_TAUTOMER_SMARTS.items()}
    kept: List[str] = []
    rejected: List[Tuple[str,str]] = []
    for smi in smiles_list:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            rejected.append((smi,"invalid_smiles"))
            continue
        reason: Optional[str] = None
        for name,patt in patterns.items():
            if patt is not None and mol.HasSubstructMatch(patt):
                reason = name
                break
        if reason is None:
            kept.append(smi)
        else:
            rejected.append((smi,reason))
    return kept,rejected
