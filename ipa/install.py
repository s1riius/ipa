# ------------------------------------------------------------
# Installation script
# ------------------------------------------------------------

import subprocess

#list of (source file, output executable)
compilation_tasks = [
    ("gfn2xtb_in_pocket_opt.cpp","gfn2xtb_in_pocket_opt.exe"),
    ("gfn2xtb_in_pocket_general_opt.cpp","gfn2xtb_in_pocket_general_opt.exe"),
    ("gfn2xtb_density.cpp","gfn2xtb_density.exe"),
    ("EDA_GFN2_3mol_solv.cpp","EDA_GFN2_3mol_solv.exe"),
    ("align.cpp","align.exe"),
    ("rmsd.cpp","rmsd.exe"),
    ("CCD_Protonate.cpp","CCD_Protonate.exe"),
    ("Protein_Protonate.cpp","Protein_Protonate.exe"),
    ("VSEPR_Protonate.cpp","VSEPR_Protonate.exe"),
    ("SDFMOL2_Protonate.cpp","SDFMOL2_Protonate.exe")
]
print("Add the installation of RDKit, Gemmi, and mrcfile")
compiler = "g++"
flags = ["-O3"]
for source,output in compilation_tasks:
    print(f"\nCompiling {source} -> {output}")
    result = subprocess.run([compiler,*flags,source,"-o",output],
                            stdout = subprocess.PIPE,
                            stderr = subprocess.PIPE,
                            text = True)
    if result.returncode != 0:
        print(f"Failed to compile {source}:\n{result.stderr}")
    else:
        print(f"{output} successfully compiled")
print("\nAll compilation tasks done.")