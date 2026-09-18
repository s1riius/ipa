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

#ifndef _File_Format_Package_
#define _File_Format_Package_
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <cctype>
#include "BiologicalData.hpp"
#include "../atoms/AtomPackage.hpp"

//description:
//functions that read or write specific file-types associated to molecules

std::string GetFileExtension(std::string _file) {
  //function that gets the file extension
  std::string extension = "";
  int fetch = -1;
  size_t ncharacters = _file.size();
  for (size_t idel = ncharacters - 1; idel >= 0; --idel) {
    if (std::string(1,_file[idel]) == ".") {
      fetch = idel;
      break;
    }
  }
  for (size_t idel = fetch + 1; idel < ncharacters; ++idel) {
    extension += _file[idel];
  }
  return extension;
}
void ReadXYZFormat(std::string geomfile, size_t & natoms, matrixE & geometry, std::vector<size_t> & atoms) {
  //function to read geometries and atomic lists from a given xyz file
  std::ifstream gfile(geomfile,std::ios::in);
  if (!gfile.is_open()) {throw std::string("ERROR: FileFormats.hpp: ReadXYZFormat(): xyz file could not be open");}
  std::string read;
  std::getline(gfile,read);
  std::istringstream (read) >> natoms;
  std::getline(gfile,read);           //skip comment
  geometry.resize(natoms,3);
  atoms.resize(natoms);
  for (size_t idatm = 0; idatm < natoms; ++idatm) {
    gfile >> read;
    atoms[idatm] = Symbol2AtomNr(read);
    for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
      gfile >> geometry(idatm + 1,idcoord + 1);
    }
  }
  gfile.close();
}
void WriteXYZFormat(std::string name, std::vector<size_t> & atoms, matrixE & geometry, int counter = -1, int prc = 7) {
    //function to write geometries into a given xyz file
    std::string filename = name;
    double precision = pow(10.0,-prc);
    double number;
    int natoms = atoms.size();
    if (counter >= 0) {filename += szt2str(counter);}
    filename += ".xyz";
    std::ofstream gfile(filename,std::ios::out);
    if (!gfile.is_open()) {throw std::string("ERROR: FileFormats.hpp: WriteXYZFormat(): could not open file for writing geometry");}
    gfile << std::fixed;
    gfile << std::setprecision(prc);
    gfile << natoms << "\n";
    gfile << "\n";
    for (size_t idx = 0; idx < natoms; ++idx) {
      filename = AtomNr2Symbol(atoms[idx]);
      gfile << filename << "    ";
      if (filename.length() == 1) {gfile << " ";}
      for (size_t idx2 = 0; idx2 < 3; ++idx2) {
        number = geometry(idx + 1,idx2 + 1);
        if (fabs(number) < precision) {number = 0.0;}
        if (number >= -precision) {gfile << " ";}
        gfile << number << "    ";
      }
      gfile << "\n";
    }
    gfile.close();
}
void ReadMOL2Format(std::string mol2file, size_t & Natoms, matrixE & geometry, std::vector<size_t> & atoms, int maxnelements = 100000000) {
  //function to read geometries and atomic lists from a given mol2 file
  std::ifstream gfile(mol2file,std::ios::in);
  if (!gfile.is_open()) {throw std::string("ERROR: FileFormats.hpp: ReadMOL2Format(): mol2 file could not be open");}
  std::string read;
  std::string discard;
  std::string build;
  std::string newelement;
  char element;
  bool proteinnaming = false;
  //move to where geometry begins
  for (size_t idline = 0; idline < maxnelements; ++idline) {
    std::getline(gfile,read);
    if (read == "@<TRIPOS>ATOM") {break;}
  }
  //now read geometry
  Natoms = 1;
  for (size_t idchar = 0; idchar < maxnelements; ++idchar) {
    gfile >> read;
    if (read == "@<TRIPOS>BOND") {break;}
    geometry.resize(Natoms,3);
    atoms.resize(Natoms);
    discard = read;
    gfile >> read;
    //process atomic number; mol2 adds stuff after the elements, take this out
    build = "";
    for (size_t idel = 0; idel < read.size(); ++idel) {
      element = static_cast<unsigned char>(read[idel]);
      if (isalpha(element)) {build += read[idel];}
      else {break;}
    }
    proteinnaming = false;
    if (ValidProteinAtom(build)) {proteinnaming = true;}
    if (proteinnaming) {newelement = Convert2Element(build);}
    else {newelement = build;}
    atoms[Natoms - 1] = Symbol2AtomNr(newelement);
    for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
      geometry(Natoms,idcoord + 1) = 0.0;
      gfile >> geometry(Natoms,idcoord + 1);
    }
    for (size_t idrepeat = 0; idrepeat < 4; ++idrepeat) {
      gfile >> read;
      discard = read;
    }
    ++Natoms;
  }
  --Natoms;
  gfile.close();
}
void FullReadMOL2Format(std::string mol2file, int & Natoms, matrixE & geometry, std::vector<size_t> & atoms, matrix<int> & Connectivity, matrixE & BondOrders, int maxneigh = 8, int maxnelements = 100000000) {
  //function to read geometries, atomic lists, connectivities and bond orders from a given mol2 file
  std::ifstream gfile(mol2file,std::ios::in);
  if (!gfile.is_open()) {throw std::string("ERROR: FileFormats.hpp: FullReadMOL2Format(): mol2 file could not be open");}
  std::string read;
  std::string discard;
  std::string build;
  std::string newelement;
  int atmA;
  int atmB;
  int counter = 0;
  int nbonds = 0;
  int natms = 0;
  char element;
  double bo;
  bool proteinnaming = false;
  //move to where geometry begins
  for (size_t idline = 0; idline < maxnelements; ++idline) {
    std::getline(gfile,read);
    if (read == "@<TRIPOS>MOLECULE") {
      std::getline(gfile,read);       //skip comment
      std::getline(gfile,read);
      //this line contains number of atoms and of bonds
      build = "";
      for (size_t idchars = 0; idchars < read.size(); ++idchars) {
        if (isblank(read[idchars])) {
          if (build.size() > 0) {
            std::stringstream convert;
            convert << build;
            if (counter == 0) {
              convert >> natms;
              ++counter;
            }
            else if (counter == 1) {
              convert >> nbonds;
              ++counter;
            }
          }
          build = "";
          continue;
        }
        build += read[idchars];
      }
    }
    if (read == "@<TRIPOS>ATOM") {break;}
  }
  //now read geometry
  Natoms = 1;
  for (size_t idchar = 0; idchar < maxnelements; ++idchar) {
    gfile >> read;
    if (read == "@<TRIPOS>BOND") {break;}
    geometry.resize(Natoms,3);
    atoms.resize(Natoms);
    discard = read;
    gfile >> read;
    //process atomic number; mol2 adds stuff after the elements, take this out
    build = "";
    for (size_t idel = 0; idel < read.size(); ++idel) {
      element = static_cast<unsigned char>(read[idel]);
      if (isalpha(element)) {build += read[idel];}
      else {break;}
    }
    proteinnaming = false;
    if (ValidProteinAtom(build)) {proteinnaming = true;}
    if (proteinnaming) {newelement = Convert2Element(build);}
    else {newelement = build;}
    atoms[Natoms - 1] = Symbol2AtomNr(newelement);
    for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
      geometry(Natoms,idcoord + 1) = 0.0;
      gfile >> geometry(Natoms,idcoord + 1);
    }
    for (size_t idrepeat = 0; idrepeat < 4; ++idrepeat) {
      gfile >> read;
      discard = read;
    }
    ++Natoms;
  }
  --Natoms;
  if (Natoms != natms) {std::cout << "WARNING: FileFormats.hpp: FullReadMOL2Format(): number of atoms declared and counted do not match" << std::endl;}
  //prepare connectivity and bond orders
  std::vector<int> neighbours(Natoms,0);
  Connectivity.resize(Natoms,maxneigh);           //assume a maximum of 8 neighbours per atom
  BondOrders.resize(Natoms,maxneigh);             //assume a maximum of 8 neighbours per atom
  //now read them
  for (size_t idchar = 0; idchar < maxnelements; ++idchar) {
    gfile >> read;
    if (read == "@<TRIPOS>SUBSTRUCTURE") {break;}
    discard = read;
    std::stringstream convert;
    convert << discard;
    convert >> counter;
    gfile >> atmA;
    gfile >> atmB;
    gfile >> read;
    if (read == "am") {bo = 1.0;}
    else if (read == "ar") {bo = 1.5;}
    else if (read == "1") {bo = 1.0;}
    else if (read == "2") {bo = 2.0;}
    else if (read == "3") {bo = 3.0;}
    ++neighbours[atmA - 1];
    ++neighbours[atmB - 1];
    Connectivity(atmA,neighbours[atmA - 1]) = atmB;
    Connectivity(atmB,neighbours[atmB - 1]) = atmA;
    BondOrders(atmA,neighbours[atmA - 1]) = bo;
    BondOrders(atmB,neighbours[atmB - 1]) = bo;
    if (counter == nbonds) {break;}
  }
  gfile.close();
}
void ReadSDFFormat(std::string geomfile, size_t & natoms, matrixE & geometry, std::vector<size_t> & atoms, int & charge, int maxnelements = 100000000) {
  //function to read geometries and atomic lists from a given sdf file
  std::ifstream gfile(geomfile,std::ios::in);
  if (!gfile.is_open()) {throw std::string("ERROR: FileFormats.hpp: ReadSDFFormat(): sdf file could not be open");}
  int sizewords = 120;          //molecule with 30 atoms
  int increment = sizewords;
  int start = -1;
  int countbetweenatoms = 0;
  int auxcount = 0;
  int firstatom = -1;
  bool e0alpha;
  bool e1alpha;
  std::vector<std::string> read(increment);
  //read the whole file to string
  for (size_t idchar = 0; idchar < maxnelements; ++idchar) {
    if (idchar + 1 == sizewords) {
      sizewords += increment;
      read.resize(sizewords);
    }
    gfile >> read[idchar];
    if (read[idchar] == "$$$$") {
      sizewords = idchar + 1;
      break;
    }
  }
  gfile.close();
  //determine where to start reading the geometry
  for (size_t idel = 0; idel < sizewords; ++idel) {
    if ((read[idel] == "V2000")||(read[idel] == "v2000")) {
      start = idel + 1;
      break;
    }
  }
  //determine what is the frequency with which atoms are declared
  for (size_t idel = start; idel < sizewords; ++idel) {
    if (read[idel].size() <= 2) {
      e0alpha = isalpha(read[idel][0]);
      e1alpha = true;
      if (read[idel].size() == 2) {e1alpha = isalpha(read[idel][1]);}
      if ((e0alpha)&&(e1alpha)) {
        ++auxcount;
        if (auxcount == 1) {firstatom = idel;}
      }
    }
    if (auxcount == 1) {++countbetweenatoms;}
    else if (auxcount == 2) {break;}
  }
  //count the number of atoms; use the fact that after coordinates there are bond orders
  natoms = 0;
  for (size_t idel = firstatom; idel < sizewords; idel += countbetweenatoms) {
    if (!isalpha(read[idel][0])) {break;}
    ++natoms;
  }
  geometry.resize(natoms,3);
  atoms.resize(natoms);
  firstatom = 0;
  //now write geometry and atom list
  for (size_t idel = start; idel < sizewords; ++idel) {
    if (isalpha(read[idel][0])) {
      atoms[firstatom] = Symbol2AtomNr(read[idel]);
      geometry(firstatom + 1,1) = stod(read[idel - 3]);
      geometry(firstatom + 1,2) = stod(read[idel - 2]);
      geometry(firstatom + 1,3) = stod(read[idel - 1]);
      ++firstatom;
    }
    if (firstatom == natoms) {
      //store last position to check whether charge was declared
      firstatom = idel + 1;
      break;
    }
  }
  //check for potentially declared charge
  for (size_t idel = firstatom; idel < sizewords; ++idel) {
    if (read[idel] == "<Formal_Charge>") {
      //charge was indeed declared
      charge = stoi(read[idel + 1]);
    }
  }
}
void FullReadSDFFormat(std::string geomfile, int & natoms, matrixE & geometry, std::vector<size_t> & atoms, int & charge, matrix<int> & Connectivity, matrixE & BondOrders, int maxneigh = 8, int maxnelements = 100000000) {
  //function to read geometries and atomic lists from a given sdf file
  std::ifstream gfile(geomfile,std::ios::in);
  if (!gfile.is_open()) {throw std::string("ERROR: FileFormats.hpp: FullReadSDFFormat(): sdf file could not be open");}
  int sizewords = 120;          //molecule with 30 atoms
  int increment = sizewords;
  int start = -1;
  int countbetweenatoms = 0;
  int auxcount = 0;
  int firstatom = -1;
  int atmA;
  int atmB;
  int counter;
  bool found_charge = false;
  bool e0alpha;
  bool e1alpha;
  double bond;
  std::vector<std::string> read(increment);
  std::string auxstrng;
  std::vector<int> neighbours;
  //read the whole file to string
  for (size_t idchar = 0; idchar < maxnelements; ++idchar) {
    if (idchar + 1 == sizewords) {
      sizewords += increment;
      read.resize(sizewords);
    }
    gfile >> read[idchar];
    if (read[idchar] == "$$$$") {
      sizewords = idchar + 1;
      break;
    }
  }
  gfile.close();
  //determine where to start reading the geometry
  for (size_t idel = 0; idel < sizewords; ++idel) {
    if ((read[idel] == "V2000")||(read[idel] == "v2000")) {
      start = idel + 1;
      break;
    }
  }
  //determine what is the frequency with which atoms are declared
  for (size_t idel = start; idel < sizewords; ++idel) {
    if (read[idel].size() <= 2) {
      e0alpha = isalpha(read[idel][0]);
      e1alpha = true;
      if (read[idel].size() == 2) {e1alpha = isalpha(read[idel][1]);}
      if ((e0alpha)&&(e1alpha)) {
        ++auxcount;
        if (auxcount == 1) {firstatom = idel;}
      }
    }
    if (auxcount == 1) {++countbetweenatoms;}
    else if (auxcount == 2) {break;}
  }
  //count the number of atoms; use the fact that after coordinates there are bond orders
  natoms = 0;
  for (size_t idel = firstatom; idel < sizewords; idel += countbetweenatoms) {
    if (!isalpha(read[idel][0])) {break;}
    ++natoms;
  }
  geometry.resize(natoms,3);
  atoms.resize(natoms);
  neighbours.resize(natoms);
  Connectivity.resize(natoms,maxneigh);
  BondOrders.resize(natoms,maxneigh);
  firstatom = 0;
  for (size_t idatm = 0; idatm < natoms; ++idatm) {
    neighbours[idatm] = 0;
    for (size_t idneigh = 0; idneigh < maxneigh; ++idneigh) {
      Connectivity(idatm + 1,idneigh + 1) = 0;
      BondOrders(idatm + 1,idneigh + 1) = 0;
    }
  }
  //now write geometry and atom list
  for (size_t idel = start; idel < sizewords; ++idel) {
    if (isalpha(read[idel][0])) {
      atoms[firstatom] = Symbol2AtomNr(read[idel]);
      geometry(firstatom + 1,1) = stod(read[idel - 3]);
      geometry(firstatom + 1,2) = stod(read[idel - 2]);
      geometry(firstatom + 1,3) = stod(read[idel - 1]);
      ++firstatom;
    }
    if (firstatom == natoms) {
      //store last position to check whether charge was declared
      firstatom = idel + 1 + countbetweenatoms - 4;
      break;
    }
  }
  //check how many zeroes are being given between bond declaration
  auxcount = 0;
  for (size_t idaux = 1; idaux < 301; ++idaux) {
    if ((read[firstatom + idaux + 3] == "M")&&(read[firstatom + idaux + 4] == "END")) {break;}
    if ((stoi(read[firstatom + idaux + 3]) != 0)&&(stoi(read[firstatom + idaux + 4]) != 0)&&(stoi(read[firstatom + idaux + 2]) != 0)) {break;}
    ++auxcount;
  }
  //check for bonds
  for (size_t idel = firstatom; idel < sizewords; ++idel) {
    if (((read[idel] == "M")&&(read[idel + 1] == "END"))||((read[idel] == "M")&&(read[idel + 1] == "CHG"))) {
      firstatom = idel + 2;
      if ((read[idel] == "M")&&(read[idel + 1] == "CHG")) {found_charge = true;}
      break;
    }
    countbetweenatoms = 0;
    atmA = stoi(read[idel]);
    atmB = stoi(read[idel + 1]);
    bond = stod(read[idel + 2]);
    if (natoms > 99) {
      //in this case it is possible that numbers are merged
      if (atmA > 999) {
        //then it did happen
        bond = atmB;
        countbetweenatoms = 1;
        if (read[idel].size() == 4) {
          auxstrng = read[idel][0];
          atmA = stoi(auxstrng);
          atmB = stoi(read[idel]) - atmA*10000;
        }
        else if (read[idel].size() == 5) {
          std::stringstream convert;
          convert << read[idel];
          convert >> atmA;
          atmA /= 1000;
          atmB = stoi(read[idel]) - atmA*1000;
        }
        else if (read[idel].size() == 6) {
          auxstrng = read[idel][0] + read[idel][1] + read[idel][2];
          atmA = stoi(auxstrng);
          atmB = stoi(read[idel]) - atmA*1000000;
        }
      }
    }
    if (bond == 4.0) {bond = 1.5;}
    else if (bond >= 5.0) {bond = 1.0;}
    ++neighbours[atmA - 1];
    ++neighbours[atmB - 1];
    Connectivity(atmA,neighbours[atmA - 1]) = atmB;
    Connectivity(atmB,neighbours[atmB - 1]) = atmA;
    BondOrders(atmA,neighbours[atmA - 1]) = bond;
    BondOrders(atmB,neighbours[atmB - 1]) = bond;
    idel += auxcount + 2 - countbetweenatoms;
  }
  //check for potentially declared charge
  if (found_charge) {
    counter = stoi(read[firstatom]);
    ++firstatom;
    charge = 0;
    for (size_t idatm = 0; idatm < counter; ++idatm) {
      charge += stoi(read[firstatom + 1]);
      firstatom += 2;
    }
  }
  else {
    for (size_t idel = firstatom; idel < sizewords; ++idel) {
      if (read[idel] == "<Formal_Charge>") {
        //charge was indeed declared
        charge = stoi(read[idel + 1]);
      }
    }
  }
}
void WriteSDFFormat(std::string filename, std::vector<size_t> & atoms, matrixE & geometry, matrix<int> & connectivity, matrix<int> & borders, std::string moleculename, int counter = -1) {
  //function that writes SDF files
  double precision = pow(10.0,-4);
  double number;
  int natoms = atoms.size();
  int nbonds = 0;
  int maxneigh = borders.cols();
  int idBtm;
  std::string spaces;
  //set filename to full
  if (counter >= 0) {filename += szt2str(counter);}
  filename += ".sdf";
  std::ofstream gfile(filename,std::ios::out);
  if (!gfile.is_open()) {throw std::string("ERROR: FileFormats.hpp: WriteSDFFormat(): could not open file for writing geometry");}
  //set options
  gfile << std::fixed;
  gfile << std::setprecision(4);
  //count bonds
  for (size_t idAtm = 0; idAtm < natoms; ++idAtm) {
    for (size_t idneigh = 0; idneigh < maxneigh; ++idneigh) {
      if (borders(idAtm + 1,idneigh + 1) == 0) {break;}
      idBtm = connectivity(idAtm + 1,idneigh + 1) - 1;
      nbonds += (idAtm < idBtm);
    }
  }
  //header
  gfile << moleculename << "\n";
  gfile << "Written by ULYSSES \n\n";
  //number of atoms
  if (natoms > 99) {spaces = "";}
  else if (natoms > 9) {spaces = " ";}
  else {spaces = "  ";}
  gfile << spaces << natoms;
  //number of bonds
  if (nbonds > 99) {spaces = "";}
  else if (nbonds > 9) {spaces = " ";}
  else {spaces = "  ";}
  gfile << spaces << nbonds << "  0  0  0  0  0  0  0  0999 V2000\n";
  //geometry
  for (size_t idAtm = 0; idAtm < natoms; ++idAtm) {
    for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
      number = geometry(idAtm + 1,idcoord + 1);
      if (fabs(number) < precision) {number = 0.0;}
      gfile << std::fixed << std::setw(10) << std::setprecision(4) << number;
    }
    gfile << " " << AtomNr2Symbol(atoms[idAtm]);
    if (AtomNr2Symbol(atoms[idAtm]).size() == 1) {gfile << " ";}
    gfile << "  0  0  0  0  0\n";
  }
  //bonds
  for (size_t idAtm = 0; idAtm < natoms; ++idAtm) {
    for (size_t idneigh = 0; idneigh < maxneigh; ++idneigh) {
      if (borders(idAtm + 1,idneigh + 1) == 0) {break;}
      idBtm = connectivity(idAtm + 1,idneigh + 1) - 1;
      if (idAtm < idBtm) {
        //write bond
        if (idAtm + 1 > 99) {spaces = "";}
        else if (idAtm + 1 > 9) {spaces = " ";}
        else {spaces = "  ";}
        gfile << spaces << idAtm + 1;
        if (idBtm + 1 > 99) {spaces = "";}
        else if (idBtm + 1 > 9) {spaces = " ";}
        else {spaces = "  ";}
        gfile << spaces << idBtm + 1;
        gfile << "  " << borders(idAtm + 1,idneigh + 1) << "  0  0  0  0\n";
      }
    }
  }
  gfile << "M  END\n$$$$\n";
  gfile.close();
}
std::string TrimSpaces(const std::string & str) {
  size_t start = 0;
  while ((start < str.size())&&(std::isspace(str[start]))) {++start;}
  size_t end = str.size();
  while ((end > start)&&(std::isspace(str[end - 1]))) {--end;}
  return str.substr(start,end - start);
}
void ReadPDBFormat(std::string geomfile, size_t & Natoms, matrixE & geometry, std::vector<size_t> & atoms) {
  //function that reads PDB files for geometries and atoms
  std::ifstream gfile(geomfile,std::ios::in);
  if (!gfile.is_open()) {throw std::string("ERROR: FileFormats.hpp: ReadPDBFormat(): pdb file could not be open");}
  std::string line;
  std::string residue;
  std::vector<std::string> atomtype;
  std::vector<std::string> atomname;
  std::vector<std::string> atomsymbol;
  Natoms = 0;
  //get data from PDBfile
  while (std::getline(gfile,line)) {
    // Check for ATOM or HETATM records
    if ((line.substr(0,4) == "ATOM")||(line.substr(0,6) == "HETATM")) {
      ++Natoms;
      geometry.resize(Natoms,3);
      residue = line.substr(17,3);
      if (line.substr(0,4) == "ATOM") {atomtype.push_back("p");}
      else if (ValidAminoAcids(residue)) {atomtype.push_back("p");}
      else {atomtype.push_back("l");}
      atomname.push_back(TrimSpaces(line.substr(12,4)));                      //atom name (columns 13-16)
      atomsymbol.push_back(TrimSpaces(line.substr(76,2)));                    //element symbol (columns 77-78)      
      geometry(Natoms,1) = std::atof(line.substr(30,8).c_str());              //x coordinate (columns 31-38)
      geometry(Natoms,2) = std::atof(line.substr(38,8).c_str());              //y coordinate (columns 39-46)
      geometry(Natoms,3) = std::atof(line.substr(46,8).c_str());              //z coordinate (columns 47-54)
    }
  }
  gfile.close();
  //convert atom names to atomic numbers
  atoms.resize(Natoms);
  for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
    //ligand atom type
    if (atomtype[idAtm] == "l") {
      if ((atomname[idAtm] == "PA")||(atomname[idAtm] == "Pa")||(atomname[idAtm] == "pa")) {atomname[idAtm] = "P";}
      else if ((atomname[idAtm] == "PB")||(atomname[idAtm] == "Pb")||(atomname[idAtm] == "pb")) {atomname[idAtm] = "P";}
      else if ((atomname[idAtm] == "PC")||(atomname[idAtm] == "Pc")||(atomname[idAtm] == "pc")) {atomname[idAtm] = "P";}
      atoms[idAtm] = Symbol2AtomNr(AtomName2Element(atomname[idAtm]));
    }
    //protein
    else {
      atoms[idAtm] = Symbol2AtomNr(Convert2Element(atomname[idAtm]));
      if ((atoms[idAtm] == 0)&&(atomsymbol[idAtm].size() != 0)) {
        size_t atmnr_from_symbol = Symbol2AtomNr(AtomName2Element(atomsymbol[idAtm]));
        if (atmnr_from_symbol != 0) {atoms[idAtm] = atmnr_from_symbol;}
      }
    }
  }
}
void FullReadPDBFormat(std::string geomfile, int & Natoms, matrixE & geometry, std::vector<size_t> & atoms, int & nheteroatoms, matrixE & hetgeom, std::vector<size_t> & hetatoms, std::vector<std::string> & atomnames, matrix<int> & chains, std::vector<std::string> & residue, std::vector<int> & resnumber, std::vector<int> & res_start, bool & hasH, int & nres) {
  //function that reads PDB files for geometries and atoms
  std::ifstream gfile(geomfile,std::ios::in);
  if (!gfile.is_open()) {throw std::string("ERROR: FileFormats.hpp: FullReadPDBFormat(): pdb file could not be open");}
  std::vector<std::string> atomtype;
  std::vector<std::string> atomsymbol;
  std::vector<std::string> chainid;
  std::string line;
  double coords[3];
  int ihetcnt = 0;
  int ires = 1;
  int nchains = 1;
  int iter = 0;
  Natoms = 0;
  nres = 0;
  nheteroatoms = 0;
  hasH = false;
  //get data from PDBfile
  while (std::getline(gfile,line)) {
    // Check for ATOM or HETATM records
    if ((line.substr(0,4) == "ATOM")||(line.substr(0,6) == "HETATM")) {
      ++Natoms;
      geometry.resize(Natoms,3);
      atomnames.push_back(TrimSpaces(line.substr(12,4)));                      //atom name (columns 13-16)
      residue.push_back(TrimSpaces(line.substr(17,3)));
      resnumber.push_back(std::atoi(line.substr(22,4).c_str()));
      chainid.push_back(line.substr(21,1).c_str());
      if (iter > 0) {nres += (chainid[Natoms] != chainid[Natoms - 1]);}
      if (atomnames[Natoms - 1] == "CA") {++nres;}                            //count the number of residues from the Calpha
      coords[0] = std::atof(line.substr(30,8).c_str());
      coords[1] = std::atof(line.substr(38,8).c_str());
      coords[2] = std::atof(line.substr(46,8).c_str());
      atomsymbol.push_back(TrimSpaces(line.substr(76,2)));
      geometry(Natoms,1) = coords[0];                                         //x coordinate (columns 31-38)
      geometry(Natoms,2) = coords[1];                                         //y coordinate (columns 39-46)
      geometry(Natoms,3) = coords[2];                                         //z coordinate (columns 47-54)
      if (line.substr(0,4) == "ATOM") {atomtype.push_back("p");}
      else {
        ++nheteroatoms;
        atomtype.push_back("l");
        hetgeom.resize(nheteroatoms,3);
        hetgeom(nheteroatoms,1) = coords[0];
        hetgeom(nheteroatoms,1) = coords[1];
        hetgeom(nheteroatoms,2) = coords[2];
      }
      ++iter;
    }
  }
  gfile.close();
  //convert atom names to atomic numbers
  atoms.resize(Natoms);
  hetatoms.resize(nheteroatoms);
  chains.resize(nchains,2);
  chains(1,1) = 1;
  if (nres == 0) {nres = 1;}
  res_start.resize(nres);
  res_start[0] = 0;
  iter = 1;
  for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
    //ligand atom type
    if (atomtype[idAtm] == "l") {
      if ((atomnames[idAtm] == "PA")||(atomnames[idAtm] == "Pa")||(atomnames[idAtm] == "pa")) {atomnames[idAtm] = "P1";}
      else if ((atomnames[idAtm] == "PB")||(atomnames[idAtm] == "Pb")||(atomnames[idAtm] == "pb")) {atomnames[idAtm] = "P2";}
      else if ((atomnames[idAtm] == "PC")||(atomnames[idAtm] == "Pc")||(atomnames[idAtm] == "pc")) {atomnames[idAtm] = "P3";}
      if (atomsymbol[idAtm].size() != 0) {atoms[idAtm] = Symbol2AtomNr(AtomName2Element(atomsymbol[idAtm]));}
      else {atoms[idAtm] = Symbol2AtomNr(AtomName2Element(atomnames[idAtm]));}
      hetatoms[ihetcnt] = atoms[idAtm];
      ++ihetcnt;
    }
    //protein
    else {
      atoms[idAtm] = Symbol2AtomNr(Convert2Element(atomnames[idAtm]));
      if ((atoms[idAtm] == 0)&&(atomsymbol[idAtm].size() != 0)) {
        size_t atmnr_from_symbol = Symbol2AtomNr(AtomName2Element(atomsymbol[idAtm]));
        if (atmnr_from_symbol != 0) {atoms[idAtm] = atmnr_from_symbol;}
      }
    }
    if (atoms[idAtm] == 1) {hasH = true;}
    if ((idAtm > 0)&&(resnumber[idAtm] != resnumber[idAtm - 1])) {
      res_start[ires] = idAtm;
      ++ires;
    }
    if ((idAtm > 0)&&(chainid[idAtm] != chainid[idAtm - 1])) {
      chains(iter,2) = idAtm;
      ++iter;
      chains(iter,1) = idAtm + 1;
    }
  }
}

#endif //_File_Format_Package_
