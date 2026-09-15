# ------------------------------------------------------------
# Functions related to energy renormalisation due to changes in protonation state
# ------------------------------------------------------------

def AcidBaseCorrection(number_of_protons):
    #function that give the acid-base corrections based on the number of protons
    return -0.16518534373*number_of_protons
