# ------------------------------------------------------------
# Interface to QM executables
# ------------------------------------------------------------

import subprocess
from pathlib import Path
import subprocess
import util.Formatting as FORMAT

# -----------------------------------------------------------------------------
# Resolve executables relative to repository root (NOT current working directory)
# -----------------------------------------------------------------------------
REPO_ROOT = Path(__file__).resolve().parent.parent

#main IP executable
GFN2_IP_EXE = REPO_ROOT / "gfn2xtb_in_pocket_opt.exe"
if not GFN2_IP_EXE.exists():
    raise FileNotFoundError(f"gfn2xtb executable not found at: {GFN2_IP_EXE}\n"
                            "Did you compile it (see README.md) or run install.py?")

#general IP executable
GEN_GFN2_IP_EXE = REPO_ROOT / "gfn2xtb_in_pocket_general_opt.exe"
if not GEN_GFN2_IP_EXE.exists():
    raise FileNotFoundError(f"generalised gfn2xtb executable not found at: {GFN2_IP_EXE}\n"
                            "Did you compile it (see README.md) or run install.py?")

#density executable
GFN2_DENSITY_EXE = REPO_ROOT / "gfn2xtb_density.exe"
if not GFN2_DENSITY_EXE.exists():
    raise FileNotFoundError(f"Density executable not found at: {GFN2_DENSITY_EXE}\n"
                            "Did you compile it (see README.md) or run install.py?")

#EDDA executable
EDDA_EXE = REPO_ROOT / "EDA_GFN2_3mol_solv.exe"
if not EDDA_EXE.exists():
    raise FileNotFoundError(f"EDDA executable not found at: {EDDA_EXE}\n"
                            "Did you compile it (see README.md) or run install.py?")

def RunIP(geometryfile, charge, alpha, kappa, ignH, ip, Telec, solvation, solventname, geomfile, energythresh, gradientthresh):
    """Run the In-Pocket executable (GFN2-xTB based)."""
    proc = subprocess.Popen([str(GFN2_IP_EXE),geometryfile,charge,alpha,kappa,ignH,ip,Telec,solvation,solventname,geomfile,energythresh,gradientthresh],
                             stdout = subprocess.PIPE,stderr = subprocess.PIPE)
    out,err = proc.communicate()
    return out,err

def RunGenIP(xyz_in, skip_file, xyz_out, ignH, constrained, solv_flag, solventname, charge, elec_temp, alpha, kappa, ethresh, gthresh):
    """Run the generalised In-Pocket executable (GFN2-xTB based)."""
    cmd = [str(GEN_GFN2_IP_EXE.resolve()),xyz_in,charge,alpha,kappa,ignH,constrained,elec_temp,solv_flag,solventname,xyz_out,ethresh,gthresh,skip_file]
    completed = subprocess.run(cmd,check = True,capture_output = True,text = True)
    return completed

def GetDensity(geometryfile, charge, Telec, solvation, solventname, densfile):
    """Run the density calculation executable."""
    proc = subprocess.Popen([str(GFN2_DENSITY_EXE),geometryfile,charge,Telec,solvation,solventname,densfile],stdout = subprocess.PIPE,stderr = subprocess.PIPE)
    proc.communicate()
    return None

def RunEDDA(geomAB, chargeAB, geomA, chargeA, geomB, chargeB):
    """Run EDDA calculation."""
    proc = subprocess.Popen([str(EDDA_EXE),geomAB,chargeAB,geomA,chargeA,geomB,chargeB],stdout = subprocess.PIPE,stderr = subprocess.PIPE)
    out,err = proc.communicate()
    return out,err

def ParseOutput(output):
    #function that parses the output of a calculation
    outputfile = output.decode("utf-8").split("\n")
    etotal = 0.0
    fetchatomicprops = False
    atmlist = []
    charges = []
    polarisabilities = []
    totalpol = 0.0
    rmsd = 0.0
    converged = True
    for iline in range(len(outputfile)):
        line = outputfile[iline]
        #check whether array data is available
        if line.startswith(">atom;charge;pol"): 
            fetchatomicprops = True
            continue
        elif line.startswith("<atom;charge;pol"): 
            fetchatomicprops = False
            continue
        #fetch data
        if line.startswith("Optimized g energy:"): etotal = float(line.replace("Optimized g energy:","").replace("Hartree","").replace(" ",""))
        elif line.startswith("RMSD with initial structure:"): rmsd = float(line.replace("RMSD with initial structure:","").replace(" ",""))
        elif "Total Polarizability" in line: totalpol = float(line.replace("Total Polarizability","").replace(" ",""))
        elif fetchatomicprops:
            cleanline = FORMAT.CleanTheLine(line)
            data = cleanline.split(" ")
            atmlist.append(int(data[0]))
            charges.append(float(data[2]))
            polarisabilities.append(float(data[3]))
        elif "no convergence in" in line: converged = False
    return etotal,rmsd,totalpol,atmlist,charges,polarisabilities,converged

def Write2File(etotal, rmsd, totalpol, atmlist, charges, polarisabilities, converged, outfile):
    wfileres = open(outfile,"w")
    wfileres.write("Energy: " + str(etotal) + " Eh\n")
    wfileres.write("RMSD with starting structure: " + str(rmsd) + " Angstrom\n")
    wfileres.write("Optimisation Converged? " + str(converged) + "\n")
    wfileres.write("Atomic properties: charge    polarisability\n")
    totalcharge = 0.0
    for idatm in range(len(atmlist)):
        wfileres.write("    " + str(atmlist[idatm]) + "   " + str(charges[idatm]) + "   " + str(polarisabilities[idatm]) + "\n")
        totalcharge += charges[idatm]
    wfileres.write("Total Charge: " + str(totalcharge) + "\n")
    wfileres.write("Total Polarisability: " + str(totalpol) + " a0^3\n")
    wfileres.close()
    return 0

def ParseEDDAOutput(output, atomsA, atomsB):
    #function parsing the EDDA output
    ElectrostaticE = 0.0
    PolarisationE = 0.0
    ChargeTransferE = 0.0
    ExchangeE = 0.0
    RepulsionE = 0.0
    DispersionE = 0.0
    SolvationE = 0.0
    fetch_ES = False
    fetch_POL = False
    fetch_CT = False
    fetch_EXC = False
    fetch_REP = False
    fetch_DISP = False
    fetch_SOLV = False
    ESmat = []
    POLmat = []
    CTmat = []
    EXCmat = []
    REPmat = []
    DISPmat = []
    SOLVmat = []
    INTmat = []
    outputfile = output.decode("utf-8").split("\n")
    for iline in range(len(outputfile)):
        line = outputfile[iline]
        if line.startswith("electrostatic energies:"):
            fetch_ES = True
            continue
        elif line.startswith("polarization energies:"):
            fetch_ES = False
            fetch_POL = True
            continue
        elif line.startswith("charge transfer energies:"):
            fetch_POL = False
            fetch_CT = True
            continue
        elif line.startswith("intermolecular exchange:"):
            fetch_CT = False
            fetch_EXC = True
            continue
        elif line.startswith("overlap repulsion:"):
            fetch_EXC = False
            fetch_REP = True
            continue
        elif line.startswith("dispersion:"):
            fetch_REP = False
            fetch_DISP = True
            continue
        elif line.startswith("solvation:"):
            fetch_DISP = False
            fetch_SOLV = True
            continue
        elif line.startswith("Ees        <<"): 
            fetch_SOLV = False
            ElectrostaticE = float(line.replace("Ees        <<","").strip().rstrip())
        elif line.startswith("Epol       <<"): 
            PolarisationE = float(line.replace("Epol       <<","").strip().rstrip())
            continue
        elif line.startswith("Ect        <<"): 
            ChargeTransferE = float(line.replace("Ect        <<","").strip().rstrip())
            continue
        elif line.startswith("Einterexc  <<"): 
            ExchangeE = float(line.replace("Einterexc  <<","").strip().rstrip())
            continue
        elif line.startswith("Erep       <<"): 
            RepulsionE = float(line.replace("Erep       <<","").strip().rstrip())
            continue
        elif line.startswith("Edisp      <<"): 
            DispersionE = float(line.replace("Edisp      <<","").strip().rstrip())
            continue
        elif line.startswith("Esolv      <<"): 
            SolvationE = float(line.replace("Esolv      <<","").strip().rstrip())
            continue
        if fetch_ES:
            data = line.strip().rstrip().replace("  "," ").split(" ")
            aux = []
            for idata in data: aux.append(float(idata))
            ESmat.append(aux)
        elif fetch_POL:
            data = line.strip().rstrip().replace("  "," ").split(" ")
            aux = []
            for idata in data: aux.append(float(idata))
            POLmat.append(aux)
        elif fetch_CT:
            data = line.strip().rstrip().replace("  "," ").split(" ")
            aux = []
            for idata in data: aux.append(float(idata))
            CTmat.append(aux)
        elif fetch_EXC:
            data = line.strip().rstrip().replace("  "," ").split(" ")
            aux = []
            for idata in data: aux.append(float(idata))
            EXCmat.append(aux)
        elif fetch_REP:
            data = line.strip().rstrip().replace("  "," ").split(" ")
            aux = []
            for idata in data: aux.append(float(idata))
            REPmat.append(aux)
        elif fetch_DISP:
            data = line.strip().rstrip().replace("  "," ").split(" ")
            aux = []
            for idata in data: aux.append(float(idata))
            DISPmat.append(aux)
        elif fetch_SOLV:
            data = line.strip().rstrip().replace("  "," ").split(" ")
            aux = []
            for idata in data: aux.append(float(idata))
            SOLVmat.append(aux)
    #now sum the interaction matrices
    for iel in range(len(ESmat)):
        aux = []
        for jel in range(len(ESmat[iel])):
            lint = ESmat[iel][jel] + POLmat[iel][jel] + CTmat[iel][jel] + EXCmat[iel][jel] + REPmat[iel][jel] + DISPmat[iel][jel] + SOLVmat[iel][jel]
            aux.append(lint)
        INTmat.append(aux)
    #now we need the local indices
    natomsA = len(atomsA)
    natomsB = len(atomsB)
    LocalRepulsionA = []
    LocalRepulsionB = []
    LocalInteractionA = []
    LocalInteractionB = []
    for iAtm in range(natomsA):
        rep_energy = 0.0
        int_energy = 0.0
        for iBtm in range(natomsB):
            rep_energy += REPmat[iAtm][natomsA + iBtm]
            int_energy += INTmat[iAtm][natomsA + iBtm]
        LocalRepulsionA.append(rep_energy)
        LocalInteractionA.append(int_energy)
    for iBtm in range(natomsB):
        rep_energy = 0.0
        int_energy = 0.0
        for iAtm in range(natomsA):
            rep_energy += REPmat[natomsA + iBtm][iAtm]
            int_energy += INTmat[natomsA + iBtm][iAtm]
        LocalRepulsionB.append(rep_energy)
        LocalInteractionB.append(int_energy)
    InteractionE = ElectrostaticE + PolarisationE + ChargeTransferE + ExchangeE + RepulsionE + DispersionE + SolvationE
    return RepulsionE,InteractionE,LocalRepulsionA,LocalRepulsionB,LocalInteractionA,LocalInteractionB

