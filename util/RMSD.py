# ------------------------------------------------------------
# Geometry-related utilitarians
# ------------------------------------------------------------

import subprocess
from pathlib import Path
#import util.Formatting as FORMAT
import util.Geometry as GEOM

#FMCM: this should go to geometry
def CleanSDF(sdffile):
    #function that cleans the SDF file
    rfile = open(sdffile,"r")
    rfilelines = rfile.readlines()
    rfile.close()
    wfile = open(sdffile,"w")
    for line in rfilelines:
        wfile.write(line)
        if line.startswith("M  END"):
            wfile.write("$$$$\n")
            break
    wfile.close()

def RMSDcalculation(structure1, structure2):
    #function to get the RMSD between two molecules
    REPO_ROOT = Path(__file__).resolve().parent.parent
    exefile = REPO_ROOT / "rmsd.exe"
    #before getting the RMSD, we need to clean the SDF files
    CleanSDF(structure1)
    CleanSDF(structure2)
    atoms1,_ = GEOM.ReadSDF(structure1.replace(".sdf",""))
    atoms2,_ = GEOM.ReadSDF(structure2.replace(".sdf",""))
    RMSD = 10000.0
    if len(atoms1) == len(atoms2):
        proc = subprocess.Popen([exefile,structure1,structure2],stdout = subprocess.PIPE,stderr = subprocess.PIPE)
        out,_ = proc.communicate()
        outputfile = out.decode("utf-8").split("\n")
        RMSD = float(outputfile[1].replace("RMSD:",""))
    return RMSD
