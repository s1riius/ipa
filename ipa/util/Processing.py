# ------------------------------------------------------------
# Util for processing active residue data
# ------------------------------------------------------------

from typing import List,Union

#---------------------------------------------------------------------
# Exception-handling
#---------------------------------------------------------------------
class PocketOptError(Exception):
    """Custom exception for pocket optimisation preparation errors."""

def process_active_residues_base(activeres: str) -> List[List[Union[int,str]]]:
    #function that processes and parses the input residues into a list that can be used to extract what matters.
    import operator
    listactiveres: List[List[Union[int,str]]] = []
    data = activeres.replace(",",";").split(";")
    for entry in data:
        entry = entry.strip()
        if not entry: continue
        chainidentifier = ""
        if entry[0].isalpha():
            chainidentifier = entry[0]
        elif entry[-1].isalpha():
            chainidentifier = entry[-1]
        numeric_part = entry.replace(chainidentifier,"")
        if not chainidentifier:
            raise PocketOptError(f"Residue '{entry}' is invalid. Use chain-prefixed format like A45")
        if "-" in numeric_part:
            limits = numeric_part.split("-")
            if len(limits) != 2:
                raise PocketOptError(f"Invalid residue range: {entry}")
            llim = int(limits[0])
            ulim = int(limits[1]) + 1
            for ires in range(llim,ulim):
                listactiveres.append([ires,chainidentifier.upper()])
        else:
            listactiveres.append([int(numeric_part),chainidentifier.upper()])
    listactiveres = sorted(listactiveres,key = operator.itemgetter(1,0))
    if not listactiveres:
        raise PocketOptError("No valid residues parsed from --residues")
    return listactiveres
