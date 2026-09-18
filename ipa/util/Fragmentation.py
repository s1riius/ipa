# ------------------------------------------------------------
# Utils associated with system fragmentation
# ------------------------------------------------------------

import operator

def ProcessActiveResidues(activeres):
    #function processing a list into active residues
    listactiveres = []
    if activeres == "all": 
        for ires in range(1,protein_resnumbr[len(protein_resnumbr) - 1] + 1): listactiveres.append([protein_resnumbr[iatm][1],protein_resnumbr[iatm][2]])
    else:
        nares = activeres.replace(",",";")
        data = nares.split(";")
        for idt in range(len(data)):
            chainidentifier = ""
            if data[idt][0].isalpha(): chainidentifier = data[idt][0]
            elif data[idt][len(data[idt]) - 1].isalpha(): chainidentifier = data[idt][len(data[idt]) - 1]
            process = data[idt].replace(chainidentifier,"")
            if "-" in process:
                limits = process.split("-")
                llim = int(limits[0])
                ulim = int(limits[1]) + 1
                for ires in range(llim,ulim): listactiveres.append([ires,chainidentifier.upper()])
            else: listactiveres.append([int(process),chainidentifier.upper()])
    #list of active residues might be unstructured; fix the issue here
    listactiveres = sorted(listactiveres,key = operator.itemgetter(1,0))
    return listactiveres

def ExtractSubSystem(geom_protein, atp_protein, atoms_protein, ires_protein, res_protein, chn_protein, ter_protein, restp_protein, geom_other, atp_other, atoms_other, ires_other, res_other, chn_other, ter_other, restp_other, listactiveres, connectors, cap, totalnumberofresidues, flipHIS, flipASN, flipGLN, connectresidues):
    #function that extracts part of a system of interest
    geom_sel = []
    atp_sel = []
    atoms_sel = []
    res_sel = []
    ires_sel = []
    chn_sel = []
    ter_sel = []
    restp_sel = []
    icount = 0
    rHIS = []
    rGLN = []
    rASN = []
    rall = []
    for ires in range(totalnumberofresidues):
        protein_res = ires_protein[ires][1]
        other_res = ires_other[ires][1]
        if other_res == -1000:
            chainres = ires_protein[ires][2]
            found = 0
            #look for protein residue, so cap or complete with missing residues
            for jres in range(len(listactiveres)):
                if (protein_res == listactiveres[jres][0]) and (chainres == listactiveres[jres][1]):
                    found = 1
                    break
            if (found == 1):
                geom_aux = geom_protein[ires]
                atp_aux = atp_protein[ires]
                atoms_aux = atoms_protein[ires]
                res_aux = res_protein[ires]
                ires_aux = ires_protein[ires]
                chn_aux = chn_protein[ires]
                ter_aux = ter_protein[ires]
                restp_aux = restp_protein[ires]
                if cap and (ires > 0) and (ter_protein[ires] != "N"):
                    previousres = [ires_protein[ires - 1][1],chainres]
                    if (not (previousres in listactiveres)) and (not (previousres in connectors)):
                        for idatm in range(len(atoms_protein[ires - 1])):
                            if (atp_protein[ires - 1][idatm] == "C") or (atp_protein[ires - 1][idatm] == "O"):
                                geom_aux.append(geom_protein[ires - 1][idatm])
                                atp_aux.append(atp_protein[ires - 1][idatm])
                                atoms_aux.append(atoms_protein[ires - 1][idatm])
                geom_sel.append(geom_aux)
                atp_sel.append(atp_aux)
                atoms_sel.append(atoms_aux)
                res_sel.append(res_aux)
                ires_sel.append(ires_aux)
                chn_sel.append(chn_aux)
                ter_sel.append(ter_aux)
                restp_sel.append(restp_aux)
                if (res_aux == "HIS") and flipHIS: 
                    rHIS.append(icount)
                    rall.append(icount)
                elif (res_aux == "ASN") and flipASN: 
                    rASN.append(icount)
                    rall.append(icount)
                elif (res_aux == "GLN") and flipGLN: 
                    rGLN.append(icount)
                    rall.append(icount)
                icount += 1
            elif (found == 0) and connectresidues:
                for jres in range(len(connectors)):
                    if (protein_res == connectors[jres][0]) and (chainres == connectors[jres][1]):
                        found = 2
                        break
                if (found == 2):
                    #add main chain only
                    geom_aux = []
                    atp_aux = []
                    atoms_aux = []
                    res_aux = res_protein[ires]
                    ires_aux = ires_protein[ires]
                    chn_aux = chn_protein[ires]
                    ter_aux = ter_protein[ires]
                    restp_aux = restp_protein[ires]
                    for idatm in range(len(atoms_protein[ires])):
                        if (atp_protein[ires][idatm] == "C") or (atp_protein[ires][idatm] == "O") or (atp_protein[ires][idatm] == "N") or (atp_protein[ires][idatm] == "CA"):
                            geom_aux.append(geom_protein[ires][idatm])
                            atp_aux.append(atp_protein[ires][idatm])
                            atoms_aux.append(atoms_protein[ires][idatm])
                    geom_sel.append(geom_aux)
                    atp_sel.append(atp_aux)
                    atoms_sel.append(atoms_aux)
                    res_sel.append(res_aux)
                    ires_sel.append(ires_aux)
                    chn_sel.append(chn_aux)
                    ter_sel.append(ter_aux)
                    restp_sel.append(restp_aux)
                    icount += 1
        else: 
            #look for other residue, so just add, no cap
            chainres = ires_other[ires][2]
            found = 0
            for jres in range(len(listactiveres)):
                if (other_res == listactiveres[jres][0]) and (chainres == listactiveres[jres][1]):
                    found = 1
                    break
            if (found == 1):
                geom_sel.append(geom_other[ires])
                atp_sel.append(atp_other[ires])
                atoms_sel.append(atoms_other[ires])
                res_sel.append(res_other[ires])
                ires_sel.append(ires_other[ires])
                chn_sel.append(chn_other[ires])
                ter_sel.append(ter_other[ires])
                restp_sel.append(restp_other[ires])
                icount += 1
    return geom_sel,atp_sel,atoms_sel,res_sel,ires_sel,chn_sel,ter_sel,restp_sel,rHIS,rGLN,rASN,rall
