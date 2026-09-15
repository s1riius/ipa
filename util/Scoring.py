# ------------------------------------------------------------
# Util for QM short-contact scoring function
# ------------------------------------------------------------

def ShortContactScore(atomlist, atomrep, atomint, geometry):
    #function that takes an atomic list, the calculated atomic repulsion indices, and the calculated interaction indices, and returns the H-coupled lists, as well as the respective short-contact scores
    heavyatm_list = []
    heavyatm_rep = []
    heavyatm_int = []
    index_conversion = []
    for idatm in range(len(atomlist)):
        element = atomlist[idatm]
        #protons skipped for now, done in second loop to ensure that all heavy atoms are allocated
        if element == "H": continue
        else:
            #add to list
            heavyatm_list.append(element)
            heavyatm_rep.append(atomrep[idatm])
            heavyatm_int.append(atomint[idatm])
            index_conversion.append([idatm,len(heavyatm_list) - 1])
    for idatm in range(len(atomlist)):
        element = atomlist[idatm]
        if element == "H":
            #find the closest heavy atom
            for idbtm in range(len(atomlist)):
                if idbtm == idatm: continue
                r2dist = 0.0
                for idcoord in range(3):
                    aux = geometry[idatm][idcoord] - geometry[idbtm][idcoord]
                    r2dist += aux*aux
                if r2dist < 1.96:                                              #1.4*1.4
                    #found connected atom
                    for idctm in range(len(index_conversion)):
                        if index_conversion[idctm][0] == idbtm:
                            heavyatm_rep[index_conversion[idctm][1]] += atomrep[idatm]
                            heavyatm_int[index_conversion[idctm][1]] += atomint[idatm]
                            break
                    break
    shortcontact = []
    for idheavy in range(len(heavyatm_list)):
        denom = heavyatm_int[idheavy]
        # Guard against divide-by-zero
        if abs(denom) < 1.0e-12:
            aux = 0.0
        else:
            aux = heavyatm_rep[idheavy]/denom
        shortcontact.append(aux)
    return heavyatm_list,heavyatm_rep,heavyatm_int,shortcontact
