/*ULYSSES, a semi-empirical package
    Copyright (C) 2023- Filipe Menezes (filipe.menezes@helmholtz-munich.de)
                        Federico Ballabio (federico.ballabio@helmholtz-munich.de)
                        Grzegorz Popowicz (grzegorz.popowicz@helmholtz-munich.de)

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; either
    version 2.1 of the License, or (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public
    License along with this library; if not, write to the Free Software 
    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA*/

#ifndef _Biological_Data_Package_
#define _Biological_Data_Package_
#include <vector>
#include <fstream>
#include <iostream>

//description:
//functions containing some definitions for biological systems

bool ValidAminoAcids(std::string & residue) {
  //list of valid amino acids
  bool valid = false;
  //base amino acids
  if ((residue == "ALA")||(residue == "ARG")||(residue == "ASN")||(residue == "ASP")||(residue == "CYS")) {valid = true;}
  else if ((residue == "GLN")||(residue == "GLU")||(residue == "GLY")||(residue == "HIS")||(residue == "HID")) {valid = true;}
  else if ((residue == "HSD")||(residue == "HIE")||(residue == "HSE")||(residue == "HIP")||(residue == "ILE")) {valid = true;}
  else if ((residue == "LEU")||(residue == "LYS")||(residue == "MET")||(residue == "PHE")||(residue == "PRO")) {valid = true;}
  else if ((residue == "SER")||(residue == "THR")||(residue == "TRP")||(residue == "TYR")||(residue == "VAL")) {valid = true;}
  //rare amino acids
  else if ((residue == "SEC")||(residue == "MSE")) {valid = true;}
  //PTMs
  else if ((residue == "SEP")||(residue == "TPO")||(residue == "PTR")||(residue == "NEP")) {valid = true;}                           //phosphorylation
  else if ((residue == "ALY")||(residue == "AME")||(residue == "FME")) {valid = true;}                                               //acylation
  else if ((residue == "MLY")||(residue == "MLZ")||(residue == "M3L")||(residue == "AGM")||(residue == "MMO")) {valid = true;}       //methylation
  else if ((residue == "MGN")||(residue == "HIC")||(residue == "MHS")||(residue == "A1I9G")||(residue == "SMC")) {valid = true;}     //methylation
  else if ((residue == "TRX")||(residue == "TOQ")||(residue == "4HT")||(residue == "4PQ")||(residue == "HRP")) {valid = true;}       //hydroxylation
  else if ((residue == "HTR")||(residue == "UX8")||(residue == "HYP")||(residue == "BHD")||(residue == "BH2")) {valid = true;}       //hydroxylation
  else if ((residue == "CSO")||(residue == "TYQ")||(residue == "LYZ")||(residue == "ILX")) {valid = true;}                           //hydroxylation
  else if (residue == "KCX") {valid = true;}                                                                                         //carboxylation
  else if ((residue == "PCA")||(residue == "7YO")) {valid = true;}                                                                   //cyclisation
  else if ((residue == "SCH")||(residue == "2CO")||(residue == "CSD")||(residue == "0CS")) {valid = true;}                           //other
  else if ((residue == "CSX")||(residue == "CSS")||(residue == "TSY")||(residue == "TQZ")) {valid = true;}                           //other
  return valid;
}
bool ValidFragment(std::string & residue) {
  //list of valid fragments
  bool valid = false;
  if ((residue == "ACE")||(residue == "NME")||(residue == "NHE")) {valid = true;}
  return valid;
}
bool ValidProteinAtom(std::string & atmtype) {
  //function listing all valid protein atoms
  bool valid = false;
  if (atmtype == "H") {valid = true;}
  else if (atmtype == "H1") {valid = true;}
  else if (atmtype == "H2") {valid = true;}
  else if (atmtype == "H3") {valid = true;}
  else if (atmtype == "HA") {valid = true;}
  else if (atmtype == "HA2") {valid = true;}
  else if (atmtype == "HA3") {valid = true;}
  else if (atmtype == "HB") {valid = true;}
  else if (atmtype == "HB1") {valid = true;}
  else if (atmtype == "HB2") {valid = true;}
  else if (atmtype == "HB3") {valid = true;}
  else if (atmtype == "HG") {valid = true;}
  else if (atmtype == "HG1") {valid = true;}
  else if (atmtype == "HG11") {valid = true;}
  else if (atmtype == "HG12") {valid = true;}
  else if (atmtype == "HG13") {valid = true;}
  else if (atmtype == "HG2") {valid = true;}
  else if (atmtype == "HG21") {valid = true;}
  else if (atmtype == "HG22") {valid = true;}
  else if (atmtype == "HG23") {valid = true;}
  else if (atmtype == "HG3") {valid = true;}
  else if (atmtype == "HD") {valid = true;}
  else if (atmtype == "HD1") {valid = true;}
  else if (atmtype == "HD11") {valid = true;}
  else if (atmtype == "HD12") {valid = true;}
  else if (atmtype == "HD13") {valid = true;}
  else if (atmtype == "HD2") {valid = true;}
  else if (atmtype == "HD21") {valid = true;}
  else if (atmtype == "HD22") {valid = true;}
  else if (atmtype == "HD23") {valid = true;}
  else if (atmtype == "HD3") {valid = true;}
  else if (atmtype == "HE") {valid = true;}
  else if (atmtype == "HE1") {valid = true;}
  else if (atmtype == "HE2") {valid = true;}
  else if (atmtype == "HE21") {valid = true;}
  else if (atmtype == "HE22") {valid = true;}
  else if (atmtype == "HE3") {valid = true;}
  else if (atmtype == "HH") {valid = true;}
  else if (atmtype == "HH1") {valid = true;}
  else if (atmtype == "HH11") {valid = true;}
  else if (atmtype == "HH12") {valid = true;}
  else if (atmtype == "HH2") {valid = true;}
  else if (atmtype == "HH21") {valid = true;}
  else if (atmtype == "HH22") {valid = true;}
  else if (atmtype == "HZ") {valid = true;}
  else if (atmtype == "HZ1") {valid = true;}
  else if (atmtype == "HZ2") {valid = true;}
  else if (atmtype == "HZ3") {valid = true;}
  else if (atmtype == "HN") {valid = true;}
  else if (atmtype == "HXT") {valid = true;}
  else if (atmtype == "C") {valid = true;}
  else if (atmtype == "CA") {valid = true;}
  else if (atmtype == "CB") {valid = true;}
  else if (atmtype == "CG") {valid = true;}
  else if (atmtype == "CG1") {valid = true;}
  else if (atmtype == "CG2") {valid = true;}
  else if (atmtype == "CD") {valid = true;}
  else if (atmtype == "CD1") {valid = true;}
  else if (atmtype == "CD2") {valid = true;}
  else if (atmtype == "CE") {valid = true;}
  else if (atmtype == "CE1") {valid = true;}
  else if (atmtype == "CE2") {valid = true;}
  else if (atmtype == "CE3") {valid = true;}
  else if (atmtype == "CH") {valid = true;}
  else if (atmtype == "CH1") {valid = true;}
  else if (atmtype == "CH2") {valid = true;}
  else if (atmtype == "CZ") {valid = true;}
  else if (atmtype == "CZ1") {valid = true;}
  else if (atmtype == "CZ2") {valid = true;}
  else if (atmtype == "CZ3") {valid = true;}
  else if (atmtype == "N") {valid = true;}
  else if (atmtype == "ND") {valid = true;}
  else if (atmtype == "ND1") {valid = true;}
  else if (atmtype == "ND2") {valid = true;}
  else if (atmtype == "NE") {valid = true;}
  else if (atmtype == "NE1") {valid = true;}
  else if (atmtype == "NE2") {valid = true;}
  else if (atmtype == "NH") {valid = true;}
  else if (atmtype == "NH1") {valid = true;}
  else if (atmtype == "NH2") {valid = true;}
  else if (atmtype == "NZ") {valid = true;}
  else if (atmtype == "O") {valid = true;}
  else if (atmtype == "OG") {valid = true;}
  else if (atmtype == "OG1") {valid = true;}
  else if (atmtype == "OD") {valid = true;}
  else if (atmtype == "OD1") {valid = true;}
  else if (atmtype == "OD2") {valid = true;}
  else if (atmtype == "OE") {valid = true;}
  else if (atmtype == "OE1") {valid = true;}
  else if (atmtype == "OE2") {valid = true;}
  else if (atmtype == "OH") {valid = true;}
  else if (atmtype == "OXT") {valid = true;}
  else if (atmtype == "SG") {valid = true;}
  else if (atmtype == "SD") {valid = true;}
  return valid;
}
std::string Convert2Element(std::string & proteintype) {
  //function listing all valid protein atoms
  std::string element = "";
  if ((proteintype == "H")||(proteintype == "H1")||(proteintype == "H2")||(proteintype == "H3")) {element = "H";}
  else if ((proteintype == "HA")||(proteintype == "HA2")||(proteintype == "HA3")) {element = "H";}
  else if ((proteintype == "HB")||(proteintype == "HB1")||(proteintype == "HB2")||(proteintype == "HB3")) {element = "H";}
  else if ((proteintype == "HG")||(proteintype == "HG1")||(proteintype == "HG2")||(proteintype == "HG3")) {element = "H";}
  else if ((proteintype == "HG11")||(proteintype == "HG12")||(proteintype == "HG13")) {element = "H";}
  else if ((proteintype == "HG21")||(proteintype == "HG22")||(proteintype == "HG23")) {element = "H";}
  else if ((proteintype == "HD")||(proteintype == "HD1")||(proteintype == "HD2")||(proteintype == "HD3")) {element = "H";}
  else if ((proteintype == "HD11")||(proteintype == "HD12")||(proteintype == "HD13")) {element = "H";}
  else if ((proteintype == "HD21")||(proteintype == "HD22")||(proteintype == "HD23")) {element = "H";}
  else if ((proteintype == "HE")||(proteintype == "HE1")||(proteintype == "HE3")) {element = "H";}
  else if (proteintype == "HOB") {element = "H";}
  else if ((proteintype == "HE2")||(proteintype == "HE21")||(proteintype == "HE22")) {element = "H";}
  else if ((proteintype == "HH")||(proteintype == "HH1")||(proteintype == "HH11")||(proteintype == "HH12")) {element = "H";}
  else if ((proteintype == "HH2")||(proteintype == "HH21")||(proteintype == "HH22")) {element = "H";}
  else if ((proteintype == "HZ")||(proteintype == "HZ1")||(proteintype == "HZ2")||(proteintype == "HZ3")) {element = "H";}
  else if ((proteintype == "HN")||(proteintype == "HD")||(proteintype == "HXT")) {element = "H";}
  else if ((proteintype == "C")||(proteintype == "CA")||(proteintype == "CB")||(proteintype == "CB1")||(proteintype == "CB2")) {element = "C";}
  else if ((proteintype == "CG")||(proteintype == "CG1")||(proteintype == "CG2")) {element = "C";}
  else if ((proteintype == "CD")||(proteintype == "CD1")||(proteintype == "CD2")) {element = "C";}
  else if ((proteintype == "CE")||(proteintype == "CE1")||(proteintype == "CE2")||(proteintype == "CE3")) {element = "C";}
  else if ((proteintype == "CH")||(proteintype == "CH1")||(proteintype == "CH2")||(proteintype == "CH3")) {element = "C";}
  else if ((proteintype == "CN")||(proteintype == "CF")||(proteintype == "CT1")||(proteintype == "CT2")||(proteintype == "CS")) {element = "C";}  
  else if (proteintype == "CX") {element = "C";}
  else if ((proteintype == "CM")||(proteintype == "CM1")||(proteintype == "CM2")||(proteintype == "CM3")) {element = "C";}  
  else if ((proteintype == "CZ")||(proteintype == "CZ1")||(proteintype == "CZ2")||(proteintype == "CZ3")) {element = "C";}
  else if ((proteintype == "N")||(proteintype == "ND")||(proteintype == "ND1")||(proteintype == "ND2")) {element = "N";}
  else if ((proteintype == "NE")||(proteintype == "NE1")||(proteintype == "NE2")||(proteintype == "N5")) {element = "N";}
  else if ((proteintype == "NH")||(proteintype == "NH1")||(proteintype == "NH2")||(proteintype == "NZ")) {element = "N";}
  else if (proteintype == "O") {element = "O";}
  else if ((proteintype == "OB")||(proteintype == "OZ")) {element = "O";}
  else if ((proteintype == "O1")||(proteintype == "O6")||(proteintype == "O7")||(proteintype == "O9")||(proteintype == "OT")) {element = "O";}
  else if ((proteintype == "OG")||(proteintype == "OG1")||(proteintype == "OD")||(proteintype == "OD1")||(proteintype == "OD2")||(proteintype == "OD3")) {element = "O";}
  else if ((proteintype == "OE")||(proteintype == "OE1")||(proteintype == "OE2")||(proteintype == "OE3")) {element = "O";}
  else if ((proteintype == "OH")||(proteintype == "OH2")) {element = "O";}
  else if ((proteintype == "OX")||(proteintype == "OXT")||(proteintype == "OAD")) {element = "O";}
  else if ((proteintype == "OQ1")||(proteintype == "OQ2")) {element = "O";}
  else if ((proteintype == "O1P")||(proteintype == "O2P")||(proteintype == "O3P")) {element = "O";}
  else if (proteintype == "P") {element = "P";}
  else if ((proteintype == "SG")||(proteintype == "SD")||(proteintype == "S1")||(proteintype == "S2")||(proteintype == "S3")||(proteintype == "S4")) {element = "S";}
  else if (proteintype == "SE") {element = "Se";}
  else if (proteintype == "ZN") {element = "Zn";}
  return element;
}
bool ValidFragmentAtom(std::string & atmtype) {
  //function listing all valid fragment atoms
  bool valid = false;
  if (atmtype == "H") {valid = true;}
  else if (atmtype == "HH31") {valid = true;}
  else if (atmtype == "HH32") {valid = true;}
  else if (atmtype == "HH33") {valid = true;}
  else if (atmtype == "HN1") {valid = true;}
  else if (atmtype == "HN2") {valid = true;}
  else if (atmtype == "C") {valid = true;}
  else if (atmtype == "CH3") {valid = true;}
  else if (atmtype == "N") {valid = true;}
  else if (atmtype == "O") {valid = true;}
  return valid;
}

#endif //_Biological_Data_Package_
