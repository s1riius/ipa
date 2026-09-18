# ------------------------------------------------------------
# CCD list generation
# ------------------------------------------------------------

#Reads mon_lib_list.cif (in the same folder as this script) and writes
#a Python file containing:
#_ALL_CCDs_MONOMER_LIBRARY = { ... }

from pathlib import Path
import gemmi

def main() -> None:
    here = Path(__file__).resolve().parent
    cif_path = here / "monomers_library/list/mon_lib_list.cif"
    out_path = here / "monomer_library_ccds.py"
    doc = gemmi.cif.read(str(cif_path))
    ccds = set()
    #columns that commonly contain component IDs in monomer/CCD list files
    candidate_tags = ("_chem_comp.id",
                      "_chem_comp_synonym.comp_id",
                      "_chem_comp_deriv.comp_id",
                      "_chem_comp_deriv.source_comp_id")
    for block in doc:
        for tag in candidate_tags:
            col = block.find_loop(tag)
            if col is None:
                continue
            for v in col:
                s = str(v).strip()
                if s and s not in (".","?"):
                    ccds.add(s)
    ccds_sorted = sorted(ccds)
    #write a stable, readable Python file
    lines = []
    lines.append("# Auto-generated from mon_lib_list.cif — do not edit by hand.\n")
    lines.append("_ALL_CCDs_MONOMER_LIBRARY = {\n")
    for i in range(0,len(ccds_sorted),50):
        chunk = ccds_sorted[i:i + 50]
        joined = ", ".join(repr(c) for c in chunk)
        lines.append(f"    {joined},\n")
    lines.append("}\n")
    out_path.write_text("".join(lines),encoding = "utf-8")
    print(f"Wrote {len(ccds_sorted)} CCD IDs to: {out_path}")

if __name__ == "__main__":
    main()
