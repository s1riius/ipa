# ------------------------------------------------------------
# List of implicit solvents available
# ------------------------------------------------------------

def SolventList(isolventname):
    osolventname = ""
    if (isolventname == "water") or (isolventname == "h2o") or (isolventname == "o"): osolventname = "water"
    elif (isolventname == "acetone") or (isolventname == "cc(o)c"): osolventname = "acetone"
    elif (isolventname == "acetonitrile") or (isolventname == "ch3cn") or (isolventname == "ccn"): osolventname = "acetonitrile"
    elif (isolventname == "aniline") or (isolventname == "phnh2") or (isolventname == "nc1ccccc1") or (isolventname == "c1ccc(cc1)n"): osolventname = "aniline"
    elif (isolventname == "benzaldehyde") or (isolventname == "phcho") or (isolventname == "occ1ccccc1") or (isolventname == "c1ccc(cc1)co"): osolventname = "benzaldehyde"
    elif (isolventname == "benzene") or (isolventname == "c6h6") or (isolventname == "phh") or (isolventname == "c1ccccc1"): osolventname = "benzene"
    elif (isolventname == "dichloromethane") or (isolventname == "ch2cl2") or (isolventname == "c(cl)cl") or (isolventname == "c(cl)(cl)"): osolventname = "dichloromethane"
    elif (isolventname == "chloroform") or (isolventname == "chcl3") or (isolventname == "c(cl)(cl)cl") or (isolventname == "c(cl)(cl)(cl)"): osolventname = "chloroform"
    elif (isolventname == "carbon disulfide") or (isolventname == "carbondisulfide") or (isolventname == "cs2") or (isolventname == "scs"): osolventname = "carbon disulfide"
    elif (isolventname == "dioxane") or (isolventname == "o1ccocc1"): osolventname = "dioxane"
    elif (isolventname == "dmf") or (isolventname == "dimethylformamide") or (isolventname == "cn(c)co"): osolventname = "dmf"
    elif (isolventname == "dmso") or (isolventname == "dimethylsulfoxide") or (isolventname == "cs(o)c") or (isolventname == "cs(c)o"): osolventname = "dmso"
    elif (isolventname == "ethanol") or (isolventname == "etoh") or (isolventname == "ch3ch2oh") or (isolventname == "cco"): osolventname = "ethanol"
    elif (isolventname == "diethyl ether") or (isolventname == "etoet") or (isolventname == "ccocc") or (isolventname == "ch3ch2och2ch3"): osolventname = "diethyl ether"
    elif (isolventname == "ethyl acetate") or (isolventname == "acoet") or (isolventname == "etoac") or (isolventname == "ccoc(o)c"): osolventname = "ethyl acetate"
    elif (isolventname == "furan") or (isolventname == "furane"): osolventname = "furane"
    elif (isolventname == "hexadecane") or (isolventname == "c16"): osolventname = "hexadecane"
    elif (isolventname == "hexane") or (isolventname == "cccccc") or (isolventname == "c6"): osolventname = "hexane"
    elif (isolventname == "methanol") or (isolventname == "meoh") or (isolventname == "co") or (isolventname == "ch3oh"): osolventname = "methanol"
    elif (isolventname == "nitromethane") or (isolventname == "meno2") or (isolventname == "ch3no2") or (isolventname == "cn(o)o"): osolventname = "nitromethane"
    elif (isolventname == "octanol") or (isolventname == "cccccccc") or (isolventname == "c8"): osolventname = "octanol"
    elif (isolventname == "phenol") or (isolventname == "phoh") or (isolventname == "oc1ccccc1") or (isolventname == "c1ccc(cc1)o"): osolventname = "phenol"
    elif (isolventname == "thf") or (isolventname == "tetrahydrofuran"): osolventname = "thf"
    elif (isolventname == "toluene") or (isolventname == "phme") or (isolventname == "cc1ccccc1"): osolventname = "toluene"
    elif (isolventname == "octanol wet") or (isolventname == "octanolwet") or (isolventname == "wet octanol") or (isolventname == "wetoctanol"): osolventname = "octanol wet"
    elif (isolventname == "gas"): osolventname = "gas"
    return osolventname
