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

#ifndef _PATTERNS_
#define _PATTERNS_
#include "Molecule.hpp"
#include "other/VSEPR.hpp"
#include "other/AddProtons.hpp"

//descrition:
//The Patterns class for assigning VSEPR codes to structures

class Patterns {
  //this is the implementation of the patterns method
  //patterns seeks to give VSEPR classifications to atoms, which is then used for further action
  //the available classification is limited to organic molecules
  //possible codes below
  //    1,0  -> terminal atoms (1,-1 is for proton)
  //    2,0  -> linear arrangement (AX2E0)
  //    3,0  -> trigonal planar arrangement (AX3E0)
  //    2,1  -> bent arrangement (AX2E1)
  //    4,0  -> tetrahedral arrangement (AX4E0)
  //    3,1  -> trigonal pyramidal arrangement (AX3E1)
  //    2,2  -> bent arrangement (AX2E2)
  //    5,0  -> trigonal bipyramidal arrangement (AX5E0)
  //    4,1  -> seesaw arrangement (AX4E1)
  //    3,2  -> t-shaped arrangement (AX3E2)
  //    2,3  -> linear arrangement (AX2E3)
  //    6,0  -> octahedral arrangement (AX6E0)
  //    5,1  -> square pyramidal arrangement (AX5E1)
  //    4,2  -> square planar arrangement (AX4E2)
  //    7,0  -> pentagonal bipyramidal arrangement (AX7E0)
  //    6,1  -> pentagonal pyramidal pyramidal arrangement (AX6E1)
  //    5,2  -> pentagonal planar arrangement (AX5E2)
  //    8,0  -> square antiprimatic arrangement (AX8E0)
protected:
  int atomA;
  int atomB;
  int atomC;
  int atomD;
  int atomE;
  int atomF;
  int atomG;
  int nC;
  int nCar;
  int nN;
  int nO;
  int molecularcharge;
  bool heavy_metal;
  size_t Natoms;
  size_t maxneigh;                                                     //maxneigh is the maximum number of neighbours allowed; this is limited by the maximum allowed by VSEPR
  double scalefactor;                                                  //scalefactor is the scale of summed covalent radii to consider a bond
  double radian2degree;
  Molecule TheSystem;
  Molecule Metals;
  std::vector<double> torsions;
  std::vector<double> UnpairedElectrons;
  std::vector<double> maxangle;
  std::vector<int> formalcharges;
  std::vector<int> neighbours;
  std::vector<int> resnumber;
  std::vector<size_t> atoms;
  std::vector<std::string> atom_types;
  std::vector<std::string> residues;
  matrixE UnpairedElectronsNeighbours;
  matrixE modVSEPR;
  matrixE geom;
  matrixE therings;
  matrix<int> connectivity;
  matrix<int> connectivity_interm;
  matrix<int> borders;
public:
  Patterns(Molecule & _mol) {
    torsions.resize(4);
    maxangle.resize(3);
    heavy_metal = false;
    maxneigh = 8;
    scalefactor = 1.25;
    this->setSystem(_mol);
    radian2degree = 180.0/pi;
    molecularcharge = 0;
  }
  ~Patterns() {}
  //getters
  matrixE VSEPRcodes() {return modVSEPR;}
  size_t MaxNeighbours() {return maxneigh;}
  double ScaleFactor() {return scalefactor;}
  std::vector<double> UnpairedElectronVector() {return UnpairedElectrons;}
  std::vector<int> AtomicCharges() {return formalcharges;}
  int MolecularCharge() {
    int mcharge = molecularcharge;
    for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
      mcharge += formalcharges[idAtm];
    }
    return mcharge;
  }
  Molecule GetMolecule() {return TheSystem;}
  //setters
  void setVSEPRcodes(matrixE & modVSEPR_) {modVSEPR = modVSEPR_;}
  void setSystem(Molecule & _mol) {
    geom = _mol.Geometry();
    atoms = _mol.Atoms();
    Natoms = atoms.size();
    TheSystem.set2System(geom,atoms,0,1,"C1");
    connectivity.resize(Natoms,maxneigh);
    //get the connectivity matrix
    TheSystem.ConnectivityMatrix(connectivity,scalefactor);
    //resize other arrays
    modVSEPR.resize(Natoms,2);
    UnpairedElectronsNeighbours.resize(Natoms,maxneigh);
    UnpairedElectrons.resize(Natoms);
    formalcharges.resize(Natoms);
  }
  void setMaxNeighbours(size_t maxn) {maxneigh = maxn;}
  void setScalingFacotr(double scal) {scalefactor = scal;}
  //writing
  void WriteSDF(std::string sdffile, std::string moleculename, int counter = -1) {
    //function that writes an SDF file
    WriteSDFFormat(sdffile,atoms,geom,connectivity,borders,moleculename,counter);
  }
  //treating metals
  void Check4Metals() {
    //function checking whether a molecule contains metals
    heavy_metal = false;
    for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
      if ((MetalHead(atoms[idAtm]))&&(atoms[idAtm] != 5)) {
        heavy_metal = true;
        break;
      }
    }
  }
  void IsolateMetals() {
    //function to remove heavy metals from the main molecular system
    matrixE new_coord(Natoms,3);
    matrixE metal_coord;
    Molecule newsystem;
    std::vector<size_t> new_atoms(Natoms);
    std::vector<size_t> metal_atoms;
    int icnt = 0;
    int imetal = 0;
    for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
      if ((!MetalHead(atoms[idAtm]))||(atoms[idAtm] == 5)) {
        new_coord(icnt + 1,1) = geom(idAtm + 1,1);
        new_coord(icnt + 1,2) = geom(idAtm + 1,2);
        new_coord(icnt + 1,3) = geom(idAtm + 1,3);
        new_atoms[icnt] = atoms[idAtm];
        ++icnt;
      }
      else {
        ++imetal;
        metal_atoms.push_back(atoms[idAtm]);
        metal_coord.resize(imetal,3);
        metal_coord(imetal,1) = geom(idAtm + 1,1);
        metal_coord(imetal,2) = geom(idAtm + 1,2);
        metal_coord(imetal,3) = geom(idAtm + 1,3);
      }
    }
    //set molecules with and without metals
    new_atoms.resize(icnt);
    new_coord.resize(icnt,3);
    newsystem.set2System(new_coord,new_atoms,0,1,"C1");
    Metals.set2System(metal_coord,metal_atoms,0,1,"C1");
    this->setSystem(newsystem);
  }
  void AddBackMetal() {
    //adding back the metals
    Molecule newsystem = TheSystem;
    newsystem.AddMolecule(Metals);
    this->setSystem(newsystem);
  }
  void AddCompleteStructure() {
    //function that makes the molecule complete again
    Molecule newsystem = TheSystem;
    this->setSystem(newsystem);
  }
  void ExpandConnectivityAndBOs() {
    //function that fixes the bond order and connectivity matrices after expanding a structure
    //start with connectivity
    for (size_t idAtm = 0; idAtm < connectivity_interm.rows(); ++idAtm) {
      for (size_t idneigh = 0; idneigh < maxneigh; ++idneigh) {
        if (connectivity_interm(idAtm + 1,idneigh + 1) == 0) {break;}
        connectivity(idAtm + 1,idneigh + 1) = connectivity_interm(idAtm + 1,idneigh + 1);
      }
    }
    //now the bond orders
    borders.resize(Natoms,maxneigh);
    for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
      for (size_t idneigh = 0; idneigh < maxneigh; ++idneigh) {
        if (connectivity(idAtm + 1,idneigh + 1) == 0) {break;}
        if ((connectivity(idAtm + 1,idneigh + 1) != 0)&&(borders(idAtm + 1,idneigh + 1) == 0)) {borders(idAtm + 1,idneigh + 1) = 1;}
      }
    }
  }
  void ReorderMoshPit(double factor = 1.25) {
    //function that looks for clashes between protons and metals
    Molecule newsystem;
    std::vector<size_t> new_atoms(Natoms);
    std::vector<int> new_formalcharges(Natoms);
    std::vector<size_t> atoms2exclude;
    matrixE new_coord(Natoms,3);
    int icnt = 0;
    bool found;
    double distance;
    double rABcov;
    //determine which atoms we want to exclude from the system because they overlap with the metal
    for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
      if ((MetalHead(atoms[idAtm]))&&(!(atoms[idAtm] == 5))) {
        for (size_t idBtm = 0; idBtm < Natoms; ++idBtm) {
          if (atoms[idBtm] == 1) {
            distance = Distance(idAtm + 1,idBtm + 1,geom);
            rABcov = factor*(AtmRadii(atoms[idAtm]) + AtmRadii(atoms[idBtm]));
            if (distance < rABcov) {atoms2exclude.push_back(idBtm);}
          }
        }
      }
    }
    atoms2exclude.push_back(-1);
    //keep only unique atoms
    for (size_t idx = 0; idx < atoms2exclude.size(); ++idx) {
      for (size_t jdx = idx + 1; jdx < atoms2exclude.size(); ++jdx) {
        if (atoms2exclude[jdx] == atoms2exclude[idx]) {
          //remove jdx
          for (size_t kdx = jdx + 1; kdx < atoms2exclude.size(); ++kdx) {
            atoms2exclude[kdx - 1] = atoms2exclude[kdx];
          }
        }
      }
    }
    //reduce vector size
    for (size_t idx = 0; idx < atoms2exclude.size(); ++idx) {
      if (atoms2exclude[idx] == -1) {
        atoms2exclude.resize(idx);
        break;
      }
    }
    //now keep only what we want
    for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
      found = false;
      for (size_t idx = 0; idx < atoms2exclude.size(); ++idx) {
        if (idAtm == atoms2exclude[idx]) {
          found = true;
          break;
        }
      }
      if (!found) {
        new_atoms[icnt] = atoms[idAtm];
        new_formalcharges[icnt] = formalcharges[idAtm];
        new_coord(icnt + 1,1) = geom(idAtm + 1,1);
        new_coord(icnt + 1,2) = geom(idAtm + 1,2);
        new_coord(icnt + 1,3) = geom(idAtm + 1,3);
        ++icnt;
      }
      else {
        //find the neighbour to the proton to be removed
        for (size_t idneigh = 0; idneigh < maxneigh; ++idneigh) {
          atomB = connectivity(idAtm + 1,idneigh + 1) - 1;
          if (atomB < 0) {break;}
          else if ((MetalHead(atoms[atomB]))&&(!(atoms[atomB] == 5))) {continue;}
          else if (atoms[atomB] == 1) {continue;}
          new_formalcharges[atomB] -= 1;
          break;
        }
      }
    }
    new_coord.resize(icnt,3);
    new_atoms.resize(icnt);
    newsystem.set2System(new_coord,new_atoms,0,1,"C1");
    formalcharges = new_formalcharges;
    this->setSystem(newsystem);
    formalcharges.resize(Natoms);
  }
  void ClassifyMetals() {
    //function assigning VSEPR codes to metal atoms
    int nneighbours;
    int Oneighbours;
    int icnt;
    int ang2check = 28;
    int idxmin;
    std::vector<int> anglecounter(4);
    std::vector<double> dist2idealangle(4);
    std::vector<double> idealangles(4);
    std::vector<double> allangles(ang2check);
    idealangles[0] = 90.0;
    idealangles[1] = 109.5;
    idealangles[2] = 120.0;
    idealangles[3] = 180.0;
    for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
      if ((MetalHead(atoms[idAtm]))&&(!(atoms[idAtm] == 5))) {
        //zeroing stuff
        nneighbours = 0;
        Oneighbours = 0;
        icnt = 0;
        for (size_t idang = 0; idang < ang2check; ++idang) {
          allangles[idang] = 0;
        }
        for (size_t idneigh = 0; idneigh < maxneigh; ++idneigh) {
          atomB = connectivity(idAtm + 1,idneigh + 1);
          if (atomB == 0) {break;}
          ++nneighbours;
          Oneighbours += (atoms[atomB - 1] == 8);
          for (size_t idneigh2 = idneigh + 1; idneigh2 < maxneigh; ++idneigh2) {
            atomC = connectivity(idAtm + 1,idneigh2 + 1);
            if (atomC == 0) {break;}
            allangles[icnt] = radian2degree*Angle(atomB,idAtm + 1,atomC,geom);
            ++icnt;
          }
        }
        //count the angles
        anglecounter[0] = 0;         //90 degrees
        anglecounter[1] = 0;         //109.5 degrees
        anglecounter[2] = 0;         //120 degrees
        anglecounter[3] = 0;         //180 degrees
        for (size_t idang = 0; idang < ang2check; ++idang) {
          //get distances to ideal angles
          if (allangles[idang] < 1.0) {break;}
          for (size_t jdang = 0; jdang < 4; ++jdang) {
            dist2idealangle[jdang] = fabs(idealangles[jdang] - allangles[idang]);
          }
          idxmin = 0;
          for (size_t idx = 1; idx < 4; ++idx) {
            if (dist2idealangle[idx] < dist2idealangle[idxmin]) {idxmin = idx;}
          }
          ++anglecounter[idxmin];
        }
        if (nneighbours == 1) {SetTerminal(idAtm);}
        else if (nneighbours == 2) {
          if (anglecounter[3] != 0) {SetSP_20(idAtm);}                     //linear molecule; not distinguishing from more complex linearity
          else if (anglecounter[2] != 0) {SetSP2_21(idAtm);}               //if the angles are around 120, then sp2-like
          else {SetSP3_22(idAtm);}                                         //else, something like water
        }
        else if (nneighbours == 3) {
          if (anglecounter[0] != 0) {SetTShape(idAtm);}                    //if there are angles of 90 degrees, then T-shaped
          else if (anglecounter[2] != 0) {SetPlanar(idAtm);}               //if the angles are around 120, then 3,0
          else {SetSP3_31(idAtm);}                                         //else, something like trigonal pyramidal
        }
        else if (nneighbours == 4) {
          if (anglecounter[1] != 0) {SetSP3_40(idAtm);}                    //if angles are around 109.5 then tetrahedral
          else {SetSquarePlanar(idAtm);}                                   //otherwise then square planar; we're ignoring seesaw
        }
        else if (nneighbours == 5) {
          if (anglecounter[2] != 0) {SetTrigonalBipyramidal(idAtm);}       //if there is at least one angle around 120 degrees
          else {SetSquarePyramidal(idAtm);}                                //if not, then has to be this
        }
        else if (nneighbours == 6) {SetOctahedral(idAtm);}
        formalcharges[idAtm] = MetalCharge(atoms[idAtm],nneighbours,Oneighbours);
      }
    }
  }
  //proton-related functions
  void RemoveProtons() {
    //function that removes protons
    int nprotons = 0;
    for (size_t idatm = 0; idatm < Natoms; ++idatm) {
      if (atoms[idatm] == 1) {++nprotons;}
    }
    if (nprotons != 0) {
      Molecule newsystem = TheSystem;
      geom = newsystem.Geometry();
      atoms = newsystem.Atoms();
      Natoms = atoms.size();
      matrixE newgeom;
      std::vector<size_t> newatoms;
      int non_protons = 0;
      for (size_t idatm = 0; idatm < Natoms; ++idatm) {
        if (atoms[idatm] != 1) {
          ++non_protons;
          newgeom.resize(non_protons,3);
          newatoms.resize(non_protons);
          newatoms[non_protons - 1] = atoms[idatm];
          newgeom(non_protons,1) = geom(idatm + 1,1);
          newgeom(non_protons,2) = geom(idatm + 1,2);
          newgeom(non_protons,3) = geom(idatm + 1,3);
        }
      }
      newsystem.set2System(newgeom,newatoms,0,1,"C1");
      this->setSystem(newsystem);
    }
  }
  void ClassifyProtons() {
    //function assigning VSEPR codes to protons
    for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
      if (atoms[idAtm] == 1) {SetProton(idAtm);}
    }
  }
  //VSEPR from bond orders
  int VSEPRfromBOs(matrixE & BondOrders) {
    //function assigning VSEPR codes from bond order data
    int nsingle = 0;
    int ndouble = 0;
    int narom = 0;
    int ntriple = 0;
    int nneigh = 0;
    int charge = 0;
    int n_Csp2_neigh = 0;
    int atomA;
    int atomB;
    int atomC;
    int atomD;
    double dihedral = 0.0;
    for (size_t idatm = 0; idatm < Natoms; ++idatm) {
      nsingle = 0;
      narom = 0;
      ndouble = 0;
      ntriple = 0;
      nneigh = 0;
      n_Csp2_neigh = 0;
      for (size_t idneigh = 0; idneigh < maxneigh; ++idneigh) {
        if (BondOrders(idatm + 1,idneigh + 1) == 0) {break;}
        ndouble += (BondOrders(idatm + 1,idneigh + 1) == 2.0);
        narom += (BondOrders(idatm + 1,idneigh + 1) == 1.5);
        ntriple += (BondOrders(idatm + 1,idneigh + 1) == 3.0);
        nsingle += (BondOrders(idatm + 1,idneigh + 1) == 1.0);
        ++nneigh;
        if ((atoms[connectivity(idatm + 1,idneigh + 1) - 1] == 6)&&(modVSEPR(connectivity(idatm + 1,idneigh + 1),1) == 3.0)&&(modVSEPR(connectivity(idatm + 1,idneigh + 1),2) == 0.0)) {++n_Csp2_neigh;}
      }
      //clearer way is to do this by atom
      if (atoms[idatm] == 1) {
        //proton
        if (nsingle > 1) {std::cout << "ERROR: Patterns.hpp: Patterns: VSEPRfromBOs(): hydrogen atom with more than 1 neighbour." << std::endl;}
        else {SetProton(idatm);}
      }
      else if ((atoms[idatm] == 5)||(atoms[idatm] == 13)) {
        //boron and aluminium
        if (nneigh == 3) {SetPlanar(idatm);}
        else if (nneigh == 4) {
          SetSP3_40(idatm);
          --charge;
        }
      }
      else if (atoms[idatm] == 6) {
        //carbon
        if ((ndouble == 0)&&(narom == 0)&&(ntriple == 0)) {SetSP3_40(idatm);}           //no insaturation
        else if ((ndouble == 1)||(narom > 1)) {SetPlanar(idatm);}                       //sp2 systems
        else if ((ndouble == 2)||(ntriple == 1)) {SetSP_20(idatm);}                     //sp systems
      }
      else if (atoms[idatm] == 7) {
        //nitrogen
        if (nneigh == 4) {
          SetSP3_40(idatm);
          ++charge;
        }
        else if ((ndouble == 0)&&(narom == 0)&&(ntriple == 0)) {                        //no insaturation
          if (n_Csp2_neigh == 0) {SetSP3_31(idatm);}
          else {SetPlanar(idatm);}
        }
        else if ((ndouble == 1)||(narom > 1)) {SetSP2_21(idatm);}                       //sp2 systems
        else if ((ndouble == 2)||(ntriple == 1)) {SetSP_20(idatm);}                     //sp systems
      }
      else if (atoms[idatm] == 8) {
        //oxygen
        if (nneigh == 3) {
          SetSP3_31(idatm);
          ++charge;
        }
        else if ((ndouble == 0)&&(narom == 0)&&(ntriple == 0)) {SetSP3_22(idatm);}      //no insaturation
        else if (ndouble == 1) {SetTerminal(idatm);}                                    //sp2 systems
        else if (narom > 1) {SetSP2_21(idatm);}                                         //sp2 systems
      }
      else if (atoms[idatm] == 9) {
        //fluorine
        if ((nneigh == 1)&&(nsingle == 1)) {SetTerminal(idatm);}
      }
      else if (atoms[idatm] == 15) {
        //phosphorous
        if (nneigh == 3) {SetSP3_31(idatm);}
        if (nneigh == 4) {SetSP3_40(idatm);}
      }
      else if ((atoms[idatm] == 16)||(atoms[idatm] == 34)) {
        //sulphur and selenium
        if (nneigh == 4) {SetSP3_40(idatm);}
        else if (nneigh == 3) {
          SetSP3_31(idatm);
          ++charge;
        }
        else if ((ndouble == 0)&&(narom == 0)&&(ntriple == 0)) {SetSP3_22(idatm);}      //no insaturation
        else if (ndouble == 1) {SetTerminal(idatm);}                                    //sp2 systems
        else if (narom > 1) {SetSP2_21(idatm);}                                         //sp2 systems
      }
      else if ((atoms[idatm] == 17)||(atoms[idatm] == 35)||(atoms[idatm] == 53)) {
        //chlorine, bromine, iodine
        if ((nneigh == 1)&&(nsingle == 1)) {SetTerminal(idatm);}
      }
      else if ((atoms[idatm] == 12)||(atoms[idatm] == 20)||(atoms[idatm] == 38)||(atoms[idatm] == 25)||(atoms[idatm] == 26)||(atoms[idatm] == 27)||(atoms[idatm] == 28)||(atoms[idatm] == 29)||(atoms[idatm] == 30)) {
        //magnesium, calcium, strontium, manganese, iron, cobalt, nickel, copper, zinc
        if (nneigh == 4) {
          //distinguish between square planar and tetrahedral
          //his is done by checking the dihedral around neighbours
          atomA = connectivity(idatm + 1,1);
          atomB = connectivity(idatm + 1,2);
          atomC = connectivity(idatm + 1,3);
          atomD = connectivity(idatm + 1,4);
          dihedral = Torsion(size_t(atomA),size_t(atomB),size_t(atomC),size_t(atomD),geom)*rad2degree;
          if (dihedral < 10.0) {SetSquarePlanar(idatm);}
          else {SetSP3_40(idatm);}
        }
      }
    }
    //second run needed for some refinements
    for (size_t idatm = 0; idatm < Natoms; ++idatm) {
      n_Csp2_neigh = 0;
      for (size_t idneigh = 0; idneigh < maxneigh; ++idneigh) {
        if (BondOrders(idatm + 1,idneigh + 1) == 0) {break;}
        if ((atoms[connectivity(idatm + 1,idneigh + 1) - 1] == 6)&&(modVSEPR(connectivity(idatm + 1,idneigh + 1),1) == 3.0)&&(modVSEPR(connectivity(idatm + 1,idneigh + 1),2) == 0.0)) {++n_Csp2_neigh;}
      }
      //still doing this by atom
      //if (atoms[idatm] == 1) {}
      //else if ((atoms[idatm] == 5)||(atoms[idatm] == 13)) {}
      //else if (atoms[idatm] == 6) {}
      if (atoms[idatm] == 7) {
        if ((modVSEPR(idatm + 1,1) == 3.0)&&(modVSEPR(idatm + 1,2) == 1.0)) {           //no insaturation
          if (n_Csp2_neigh != 0) {SetPlanar(idatm);}
        }
      }
      //else if (atoms[idatm] == 8) {}
      //else if (atoms[idatm] == 9) {}
      //else if (atoms[idatm] == 15) {}
      //else if ((atoms[idatm] == 16)||(atoms[idatm] == 34)) {}
      //else if ((atoms[idatm] == 17)||(atoms[idatm] == 35)||(atoms[idatm] == 53)) {}
    }
    return charge;
  }
  //reset bond orders
  void GetNewBOs(matrixE & BondOrders) {
    //function that takes a bond order matrix, and converts to the integer-only bond order matrix
    int new_bo;
    int Natom_local = BondOrders.rows();
    for (size_t idatm = 0; idatm < Natom_local; ++idatm) {
      for (size_t idneigh = 0; idneigh < maxneigh; ++idneigh) {
        new_bo = int(BondOrders(idatm + 1,idneigh + 1));
        if (new_bo == 0) {break;}
        if (BondOrders(idatm + 1,idneigh + 1) == 1.5) {new_bo = 4;}
        borders(idatm + 1,idneigh + 1) = new_bo;
      }
    }
  }
  //protonation of SDF and mol2 formats
  void ProtonateFromSDFMOL2(std::string molfile) {
    //function that takes an SDF or MOL2 file, gets bond orders, and protonates accordingly
    std::string extension = GetFileExtension(molfile);
    int charge = 0;
    matrixE BondOrders;
    int Natoms_ = 0;
    if (extension == "mol2") {FullReadMOL2Format(molfile,Natoms_,geom,atoms,connectivity,BondOrders);}
    else if (extension == "sdf") {FullReadSDFFormat(molfile,Natoms_,geom,atoms,charge,connectivity,BondOrders);}
    Natoms = Natoms_;
    charge = this->VSEPRfromBOs(BondOrders);
    AddProtons(TheSystem,modVSEPR,connectivity);
    AddCompleteStructure();
    ExpandConnectivityAndBOs();
    atoms = TheSystem.Atoms();
    Natoms = atoms.size();
    geom = TheSystem.Geometry();
    this->GetNewBOs(BondOrders);
  }
  void WriteProtonatedSDF(std::string filename, std::string moleculename) {
    //function that writes an SDF file for a protonated molecular system
    WriteSDFFormat(filename,atoms,geom,connectivity,borders,moleculename);
  }
  //protonation of proteins
  void ProtonateProteinFromPDB(std::string pdbfile) {
    //function that protonates a peptide (protein) based on PDB
    this->ReadMoleculeFromPDB(pdbfile);
    Residue2VSEPR();
    AddProtons(TheSystem,modVSEPR,connectivity);
    AddCompleteStructure();
  }
  void Residue2VSEPR() {
    //function that classifies amino acid atoms according to VSEPR
    bool found;
    modVSEPR.resize(Natoms,2);
    formalcharges.resize(Natoms);
    for (size_t idatm = 0; idatm < Natoms; ++idatm) {
      formalcharges[idatm] = 0;
      if (atom_types[idatm] == "C") {SetPlanar(idatm);}
      else if (atom_types[idatm] == "CA") {SetSP3_40(idatm);}
      else if ((atom_types[idatm] == "CB")||(atom_types[idatm] == "CB1")||(atom_types[idatm] == "CB2")) {SetSP3_40(idatm);}
      else if (atom_types[idatm] == "CG") {
        if ((residues[idatm] == "ARG")||(residues[idatm] == "MMO")||(residues[idatm] == "AGM")||(residues[idatm] == "GLN")||(residues[idatm] == "MGN")||(residues[idatm] == "GLH")||(residues[idatm] == "GLU")||(residues[idatm] == "LEU")||(residues[idatm] == "LYS")||(residues[idatm] == "LYZ")||(residues[idatm] == "KCX")||(residues[idatm] == "MLZ")||(residues[idatm] == "MLY")||(residues[idatm] == "M3L")||(residues[idatm] == "ALY")||(residues[idatm] == "MET")||(residues[idatm] == "MSE")||(residues[idatm] == "AME")||(residues[idatm] == "FME")||(residues[idatm] == "PRO")||(residues[idatm] == "PCA")||(residues[idatm] == "HYP")||(residues[idatm] == "7YO")) {SetSP3_40(idatm);}
        else if ((residues[idatm] == "ASN")||(residues[idatm] == "ASH")||(residues[idatm] == "ASP")||(residues[idatm] == "BHD")||(residues[idatm] == "BH2")||(residues[idatm] == "HIS")||(residues[idatm] == "HID")||(residues[idatm] == "HSD")||(residues[idatm] == "HIE")||(residues[idatm] == "HSE")||(residues[idatm] == "HIP")||(residues[idatm] == "NEP")||(residues[idatm] == "HIC")||(residues[idatm] == "MHS")||(residues[idatm] == "A1I9G")||(residues[idatm] == "PHE")||(residues[idatm] == "TYR")||(residues[idatm] == "TYQ")||(residues[idatm] == "PTR")) {
          SetPlanar(idatm);
          formalcharges[idatm] += (residues[idatm] == "ASP") + (residues[idatm] == "ASH") + (residues[idatm] == "BHD") + (residues[idatm] == "BH2");
        }
        else if ((residues[idatm] == "TRP")||(residues[idatm] == "TRX")||(residues[idatm] == "TOQ")||(residues[idatm] == "4HT")||(residues[idatm] == "4PQ")||(residues[idatm] == "HRP")||(residues[idatm] == "HTR")||(residues[idatm] == "UX8")) {
          SetPlanar(idatm);
        }        
      }
      else if ((atom_types[idatm] == "CG1")||(atom_types[idatm] == "CG2")) {SetSP3_40(idatm);}
      else if (atom_types[idatm] == "CD") {
        if ((residues[idatm] == "GLN")||(residues[idatm] == "MGN")||(residues[idatm] == "GLH")||(residues[idatm] == "GLU")||(residues[idatm] == "PCA")||(residues[idatm] == "7YO")) {
          SetPlanar(idatm);
          formalcharges[idatm] += (residues[idatm] == "GLU") + (residues[idatm] == "GLH");
        }
        else {SetSP3_40(idatm);}
      }
      else if (atom_types[idatm] == "CE") {SetSP3_40(idatm);}
      else if ((atom_types[idatm] == "CD1")||(atom_types[idatm] == "CD2")) {
        if ((residues[idatm] == "ILE")||(residues[idatm] == "LEU")||(residues[idatm] == "ILX")) {SetSP3_40(idatm);}
        else if ((residues[idatm] == "HIS")||(residues[idatm] == "HID")||(residues[idatm] == "HSD")||(residues[idatm] == "HIE")||(residues[idatm] == "HSE")||(residues[idatm] == "HIP")||(residues[idatm] == "NEP")||(residues[idatm] == "HIC")||(residues[idatm] == "MHS")||(residues[idatm] == "A1I9G")||(residues[idatm] == "PHE")||(residues[idatm] == "TYR")||(residues[idatm] == "TYQ")||(residues[idatm] == "PTR")) {SetPlanar(idatm);}
        else if ((residues[idatm] == "TRP")||(residues[idatm] == "TRX")||(residues[idatm] == "TOQ")||(residues[idatm] == "4HT")||(residues[idatm] == "4PQ")||(residues[idatm] == "HRP")||(residues[idatm] == "HTR")||(residues[idatm] == "UX8")) {SetPlanar(idatm);}          
      }
      else if ((atom_types[idatm] == "CE1")||(atom_types[idatm] == "CE2")||(atom_types[idatm] == "CE3")) {
        if (residues[idatm] == "AGM") {SetSP3_40(idatm);}
        else {SetPlanar(idatm);}
      }
      else if ((atom_types[idatm] == "CZ")||(atom_types[idatm] == "CZ2")||(atom_types[idatm] == "CZ3")) {
        if (residues[idatm] == "HIC") {SetSP3_40(idatm);}
        else {
          SetPlanar(idatm);
          formalcharges[idatm] += (atom_types[idatm] == "CZ")*((residues[idatm] == "ARG") + (residues[idatm] == "MMO") + (residues[idatm] == "AGM"));
        }
      }
      else if ((atom_types[idatm] == "CH")||(atom_types[idatm] == "CX")) {
        SetPlanar(idatm);
        formalcharges[idatm] += (atom_types[idatm] == "CX");
      }
      else if (atom_types[idatm] == "CH1") {SetSP3_40(idatm);}      
      else if (atom_types[idatm] == "CH2") {
        if ((residues[idatm] == "TRP")||(residues[idatm] == "TRX")||(residues[idatm] == "TOQ")||(residues[idatm] == "4HT")||(residues[idatm] == "4PQ")||(residues[idatm] == "HRP")||(residues[idatm] == "HTR")||(residues[idatm] == "UX8")) {SetPlanar(idatm);}
        else if (residues[idatm] == "MLY") {SetSP3_40(idatm);}
      }
      else if ((atom_types[idatm] == "CH3")||(atom_types[idatm] == "CS")||(atom_types[idatm] == "CT2")) {SetSP3_40(idatm);}
      else if ((atom_types[idatm] == "CF")||(atom_types[idatm] == "CM")||(atom_types[idatm] == "CM1")||(atom_types[idatm] == "CM2")||(atom_types[idatm] == "CM3")) {SetSP3_40(idatm);}                  
      else if (atom_types[idatm] == "CN") {
        if (residues[idatm] == "FME") {SetPlanar(idatm);}
        else if (residues[idatm] == "MMO") {SetSP3_40(idatm);}
      }
      else if (atom_types[idatm] == "CT1") {SetPlanar(idatm);}
      else if (atom_types[idatm] == "N") {
        if ((residues[idatm] == "AME")||(residues[idatm] == "FME")||(residues[idatm] == "PCA")||(residues[idatm] == "7YO")) {SetPlanar(idatm);}
        else if (residues[idatm] == "MMO") {
          SetSP3_40(idatm);
          ++formalcharges[idatm];
        }
        else {
          //check whether there is a C neighbour
          found = false;
          for (size_t idneigh = 0; idneigh < maxneigh; ++idneigh) {
            if (connectivity(idatm + 1,idneigh + 1) == 0) {break;}
            if (atom_types[connectivity(idatm + 1,idneigh + 1) - 1] == "C") {
              found = true;
              break;
            }
          }
          if (found) {SetPlanar(idatm);}
          else {
            SetSP3_40(idatm);
            ++formalcharges[idatm];
          }
        }
      }
      //                            HIS and variants
      else if (atom_types[idatm] == "ND1") {
        if ((residues[idatm] == "HIS")||(residues[idatm] == "HIC")||(residues[idatm] == "MHS")||(residues[idatm] == "A1I9G")) {
          if (connectivity(idatm + 1,3) != 0) {
            SetPlanar(idatm);
            found = false;
            for (size_t idbtm = idatm + 1; idbtm < Natoms; ++idbtm) {
              if (resnumber[idbtm] != resnumber[idatm]) {break;}
              if (atom_types[idbtm] == "NE2") {
                found = true;
                if ((modVSEPR(idbtm + 1,1) == 0)&&(modVSEPR(idbtm + 1,2) == 0)) {
                  if (connectivity(idbtm + 1,3) == 0) {SetSP2_21(idbtm);}
                  else {
                    SetPlanar(idbtm);
                    ++formalcharges[idatm];
                  }
                }
                else if ((modVSEPR(idbtm + 1,1) == 3)&&(modVSEPR(idbtm + 1,2) == 0)) {
                  if (connectivity(idbtm + 1,3) == 0) {SetSP2_21(idbtm);}
                  else {++formalcharges[idatm];}
                }
                break;
              }
            }
            if (!found) {
              for (int idbtm = idatm - 1; idbtm >= 0; --idbtm) {
                if (resnumber[idbtm] != resnumber[idatm]) {break;}
                if (atom_types[idbtm] == "NE2") {
                  if ((modVSEPR(idbtm + 1,1) == 0)&&(modVSEPR(idbtm + 1,2) == 0)) {
                    if (connectivity(idbtm + 1,3) == 0) {SetSP2_21(idbtm);}
                    else {
                      SetPlanar(idbtm);
                      ++formalcharges[idatm];
                    }
                  }
                  else if ((modVSEPR(idbtm + 1,1) == 3)&&(modVSEPR(idbtm + 1,2) == 0)) {
                    if (connectivity(idbtm + 1,3) == 0) {SetSP2_21(idbtm);}
                    else {++formalcharges[idatm];}
                  }
                  break;
                }
              }
            }
          }
          else {SetSP2_21(idatm);}
        }
        else {
          if (residues[idatm] == "HID") {SetPlanar(idatm);}
          else if (residues[idatm] == "HSD") {SetPlanar(idatm);}
          else if ((residues[idatm] == "HIE")||(residues[idatm] == "NEP")||(residues[idatm] == "HSE")) {SetSP2_21(idatm);}
          else if (residues[idatm] == "HIP") {
            SetPlanar(idatm);
            ++formalcharges[idatm];
          }
        }
      }
      //                            ASN
      else if (atom_types[idatm] == "ND2") {SetPlanar(idatm);}
      //                             ARG                          TRP                           GLN,HIS + variants
      else if ((atom_types[idatm] == "NE")||(atom_types[idatm] == "NE1")||(atom_types[idatm] == "NE2")) {
        if ((modVSEPR(idatm + 1,1) == 0)&&(modVSEPR(idatm + 1,2) == 0)) {
          if ((residues[idatm] == "ARG")||(residues[idatm] == "AGM")||(residues[idatm] == "MMO")||(residues[idatm] == "GLN")||(residues[idatm] == "MGN")||
              (residues[idatm] == "HIS")||(residues[idatm] == "HIE")||(residues[idatm] == "HSE")||(residues[idatm] == "HIP")||(residues[idatm] == "NEP")||
              (residues[idatm] == "HIC")||(residues[idatm] == "MHS")||(residues[idatm] == "A1I9G")||
              (residues[idatm] == "TRP")||(residues[idatm] == "TRX")||(residues[idatm] == "TOQ")||(residues[idatm] == "4HT")||
              (residues[idatm] == "4PQ")||(residues[idatm] == "HRP")||(residues[idatm] == "HTR")||(residues[idatm] == "UX8")) {SetPlanar(idatm);}
          else if ((residues[idatm] == "HID")||(residues[idatm] == "HSD")) {SetSP2_21(idatm);}
          else {SetSP3_40(idatm);}
        }
      }
      else if (atom_types[idatm] == "NZ") {
        if ((residues[idatm] == "LYS")||(residues[idatm] == "LYZ")||(residues[idatm] == "MLZ")||(residues[idatm] == "MLY")||(residues[idatm] == "M3L")) {
          ++formalcharges[idatm];
          SetSP3_40(idatm);
        }
        else if ((residues[idatm] == "ALY")||(residues[idatm] == "KCX")) {SetPlanar(idatm);}
      }
      else if ((atom_types[idatm] == "NH1")||(atom_types[idatm] == "NH2")||(atom_types[idatm] == "N5")) {SetPlanar(idatm);}
      else if ((residues[idatm] == "HOH")||(residues[idatm] == "SOL")||(residues[idatm] == "WAT")) {SetSP3_22(idatm);}
      else if (atom_types[idatm] == "O") {SetTerminal(idatm);}
      else if (atom_types[idatm] == "OD") {
        if ((residues[idatm] == "CSO")||(residues[idatm] == "2CO")) {SetSP3_22(idatm);}
        else if ((residues[idatm] == "CSX")) {SetTerminal(idatm);}  
      }
      else if ((atom_types[idatm] == "OG")||(atom_types[idatm] == "OG1")||(atom_types[idatm] == "O9")||(atom_types[idatm] == "OB")||(atom_types[idatm] == "OD3")) {SetSP3_22(idatm);}
      else if ((atom_types[idatm] == "OD1")||(atom_types[idatm] == "OD2")) {
        if ((residues[idatm] == "ASH")&&(atom_types[idatm] == "OD2")) {SetSP2_21(idatm);}
        else {
          if ((residues[idatm] == "HYP")||(residues[idatm] == "ILX")) {SetSP3_22(idatm);}
          else {
            if (connectivity(idatm + 1,2) != 0) {
              //oxygen atom with two neighbours
              SetSP2_21(idatm);
            }
            else {
              SetTerminal(idatm);
              formalcharges[idatm] -= (residues[idatm] == "ASP") + (residues[idatm] == "ASH") + (residues[idatm] == "BHD") + (residues[idatm] == "BH2") + (residues[idatm] == "CSD");
            }
          }
        }
      }
      else if ((atom_types[idatm] == "O1")||(atom_types[idatm] == "OE")) {
        if (residues[idatm] == "2CO") {SetSP3_22(idatm);}
        else {SetTerminal(idatm);}
      }
      else if ((atom_types[idatm] == "OE1")||(atom_types[idatm] == "OE2")) {
        if ((residues[idatm] == "GLH")&&(atom_types[idatm] == "OE2")) {SetSP2_21(idatm);}
        else {
          if (connectivity(idatm + 1,2) != 0) {
            //oxygen atom with two neighbours
            SetSP2_21(idatm);
          }
          else {
            SetTerminal(idatm);
            formalcharges[idatm] -= (residues[idatm] == "GLU") + (residues[idatm] == "GLH");
          }
        }
      }
      else if ((atom_types[idatm] == "OQ1")||(atom_types[idatm] == "OQ2")) {
        if (connectivity(idatm + 1,2) != 0) {
          //oxygen atom with two neighbours
          SetSP2_21(idatm);
        }
        else {
          SetTerminal(idatm);
          formalcharges[idatm] -= 1;
        }
      }      
      else if ((atom_types[idatm] == "OH")||(atom_types[idatm] == "OZ")) {
        if ((residues[idatm] == "TYR")||(residues[idatm] == "TYQ")) {SetSP2_21(idatm);}
        else if (residues[idatm] == "PTR") {SetSP2_21(idatm);} // phosphotyrosine bridge oxygen (CZ-O-P)
        else if (residues[idatm] == "ALY") {SetTerminal(idatm);}
        else if ((residues[idatm] == "HTR")||(residues[idatm] == "LYZ")) {SetSP3_22(idatm);}
      }
      else if ((atom_types[idatm] == "OH2")||(atom_types[idatm] == "O6")||(atom_types[idatm] == "O7")||(atom_types[idatm] == "OE3")||(atom_types[idatm] == "OAD")||(atom_types[idatm] == "OX")) {
        SetSP2_21(idatm);
      }      
      else if (atom_types[idatm] == "OXT") {
        SetTerminal(idatm);
        formalcharges[idatm] = -1;
      }
      else if ((atom_types[idatm] == "O1P")||(atom_types[idatm] == "O2P")||(atom_types[idatm] == "O3P")) {
        SetTerminal(idatm);
        formalcharges[idatm] = -1;
      }
      else if ((atom_types[idatm] == "O1")||(atom_types[idatm] == "OT")) {
        SetTerminal(idatm);
      }      
      else if (atom_types[idatm] == "P") {
        SetSP3_40(idatm);
        formalcharges[idatm] = 1;
      }
      else if (atom_types[idatm] == "SG") {
        if ((residues[idatm] == "CSD")||(residues[idatm] == "0CS")||(residues[idatm] == "CSX")) {
          formalcharges[idatm] += (residues[idatm] == "CSD");
          SetSP3_31(idatm);
        }
        else {SetSP3_22(idatm);}
      }
      else if ((atom_types[idatm] == "SD")||(atom_types[idatm] == "SE")||(atom_types[idatm] == "S1")||(atom_types[idatm] == "S2")||(atom_types[idatm] == "S3")||(atom_types[idatm] == "S4")) {SetSP3_22(idatm);}
    }
    
    ///for (size_t idx = 0; idx < residues.size(); ++idx) {
    ///  std::cout << residues[idx] << "  " << modVSEPR(idx + 1,1) << "," << modVSEPR(idx + 1,2) << " " << formalcharges[idx] << std::endl;
    ///}
  }
  //Chemical Definition File functions
  void ProtonateFromCDF(std::string pdbfile, std::string cdffile) {
    //function that protonates a ligand based on a CDF file
    this->ReadMoleculeFromPDB(pdbfile);
    this->ReadCCD(cdffile);
    this->FixConnectivity();
    this->CCD2VSEPR();
    connectivity_interm = connectivity;
    AddProtons(TheSystem,modVSEPR,connectivity);
    AddCompleteStructure();
    ExpandConnectivityAndBOs();
  }
  void ReadMoleculeFromPDB(std::string pdbfile) {
    //function to get molecular data from PDB
    std::string extension = GetFileExtension(pdbfile);
    if (extension != "pdb") {throw std::string("ERROR: Patterns.hpp: ReadMoleculeFromPDB(): expected PDB file");}
    int discard;
    int discardtoo;
    int totalatoms;
    bool throwaway;
    matrixE garbage;
    matrix<int> also2discard;
    std::vector<size_t> garbaggio;
    std::vector<int> res_start;
    FullReadPDBFormat(pdbfile,totalatoms,geom,atoms,discard,garbage,garbaggio,atom_types,also2discard,residues,resnumber,res_start,throwaway,discardtoo);
    Natoms = totalatoms;
    TheSystem.set2System(geom,atoms,0,1,"C1");
  }
  void ReadCCD(std::string cdffile) {
    std::ifstream gfile(cdffile,std::ios::in);
    if (!gfile.is_open()) {throw std::string("ERROR: Patterns.hpp: ReadCCD(): ccd file could not be open");}
    std::vector<std::string> read(4);
    std::vector<double> coords(3);
    std::vector<std::string> atom1;
    std::vector<std::string> atom2;
    std::vector<std::string> bo;
    bool found;
    int icnt = 1;
    int npairs;
    int idxA;
    int idxB;
    Natoms = atom_types.size();
    connectivity.resize(Natoms,maxneigh);
    borders.resize(Natoms,maxneigh);
    neighbours.resize(Natoms);
    //get the important data from CDF
    for (size_t idx = 0; idx < 1000000000000; ++idx) {
      gfile >> read[0];
      if (read[0] == "_chem_comp.pdbx_formal_charge") {
        gfile >> read[0];
        molecularcharge = std::atoi(read[0].c_str());
      }
      if (read[0] == "_chem_comp_bond.comp_id") {
        gfile >> read[0];
        if (read[0] == "_chem_comp_bond.atom_id_1") {
          gfile >> read[0];
          if (read[0] == "_chem_comp_bond.atom_id_2") {
            gfile >> read[0];
            if (read[0] == "_chem_comp_bond.value_order") {
              gfile >> read[0];
              gfile >> read[0];
              gfile >> read[0];
              //discard
              for (size_t jdx = 0; jdx < 1000000000000; ++jdx) {
                gfile >> read[0];
                if ((read[0] == "#")||(read[0] == "loop_")) {break;}
                gfile >> read[1];
                gfile >> read[2];
                gfile >> read[3];
                gfile >> read[0];
                gfile >> read[0];
                gfile >> read[0];
                if (read[1][0] == '"') {
                  read[0] = read[1].substr(1,read[1].size() - 2);
                  read[1] = read[0];
                }
                if (read[2][0] == '"') {
                  read[0] = read[2].substr(1,read[2].size() - 2);
                  read[2] = read[0];
                }
                atom1.push_back(read[1]);
                atom2.push_back(read[2]);
                bo.push_back(read[3]);
              }
              break;
            }
          }
        }
      }
    }
    gfile.close();
    npairs = atom1.size();
    //get connectivity and list of atomic numbers
    for (size_t idatm = 0; idatm < Natoms; ++idatm) {
      neighbours[idatm] = 0;
      for (size_t idneigh = 0; idneigh < maxneigh; ++idneigh) {
        connectivity(idatm + 1,idneigh + 1) = 0;
      }
    }
    for (size_t ipair = 0; ipair < npairs; ++ipair) {
      idxA = -1;
      idxB = -1;
      //correct
      if (atom1[ipair] == "PA") {atom1[ipair] = "P1";}
      else if (atom1[ipair] == "PB") {atom1[ipair] = "P2";}
      else if (atom1[ipair] == "PC") {atom1[ipair] = "P3";}
      if (atom2[ipair] == "PA") {atom2[ipair] = "P1";}
      else if (atom2[ipair] == "PB") {atom2[ipair] = "P2";}
      else if (atom2[ipair] == "PC") {atom2[ipair] = "P3";}
      for (size_t idatm = 0; idatm < Natoms; ++idatm) {
        if (atom_types[idatm] == atom1[ipair]) {idxA = idatm;}
        else if (atom_types[idatm] == atom2[ipair]) {idxB = idatm;}
        if ((idxA >= 0)&&(idxB >= 0)) {break;}
      }
      if ((idxA >= 0)&&(idxB >= 0)) {
        ++neighbours[idxA];
        ++neighbours[idxB];
        connectivity(idxA + 1,neighbours[idxA]) = idxB + 1;
        connectivity(idxB + 1,neighbours[idxB]) = idxA + 1;
        borders(idxA + 1,neighbours[idxA]) = (bo[ipair] == "SING") + (bo[ipair] == "sing") + 2*(bo[ipair] == "DOUB") + 2*(bo[ipair] == "doub") + 3*(bo[ipair] == "TRIP") + 3*(bo[ipair] == "trip");
        borders(idxB + 1,neighbours[idxB]) = (bo[ipair] == "SING") + (bo[ipair] == "sing") + 2*(bo[ipair] == "DOUB") + 2*(bo[ipair] == "doub") + 3*(bo[ipair] == "TRIP") + 3*(bo[ipair] == "trip");
      }
      //there are no implicit double or triple bonds
      else if ((idxA >= 0)&&(idxB < 0)) {++neighbours[idxA];}
      else if ((idxA < 0)&&(idxB >= 0)) {++neighbours[idxB];}
    }
  }
  void FixConnectivity() {
    //the connectivity from CCD file might contain holes due to missing atoms
    int neighbors_search;
    for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
      neighbors_search = neighbours[idAtm];
      for (size_t idneigh = 0; idneigh < neighbors_search; ++idneigh) {
        if (connectivity(idAtm + 1,idneigh + 1) == 0) {
          //shift
          for (size_t idneigh2 = idneigh + 1; idneigh2 < maxneigh; ++idneigh2) {
            connectivity(idAtm + 1,idneigh2) = connectivity(idAtm + 1,idneigh2 + 1);
          }
          --idneigh;
          --neighbors_search;
        }
      }
    }
  }
  void CCD2VSEPR() {
    //function to convert CCD read data to VSEPR codes
    int maxbo;
    modVSEPR.resize(Natoms,2);
    for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
      maxbo = 1;
      for (size_t idneigh = 0; idneigh < neighbours[idAtm]; ++idneigh) {
        if ((borders(idAtm + 1,idneigh + 1) == 2)&&(maxbo < 2)) {maxbo = 2;}
        else if (borders(idAtm + 1,idneigh + 1) == 3) {maxbo = 3;}
      }
      if (maxbo == 1) {
        if (neighbours[idAtm] == 4) {
          if ((atoms[idAtm] == 5)||(atoms[idAtm] == 6)||(atoms[idAtm] == 14)) {SetSP3_40(idAtm);}
          else if ((atoms[idAtm] == 7)|| (atoms[idAtm] == 15)) {SetSP3_40(idAtm);}
        }
        else if (neighbours[idAtm] == 3) {
          if (atoms[idAtm] == 5) {SetPlanar(idAtm);}
          if ((atoms[idAtm] == 6)|| (atoms[idAtm] == 14)) {SetPlanar(idAtm);}
          else if ((atoms[idAtm] == 7)|| (atoms[idAtm] == 15)) {SetUNDEFINED(idAtm);}
        }
        else if (neighbours[idAtm] == 2) {
          if ((atoms[idAtm] == 8)|| (atoms[idAtm] == 16)) {SetSP3_22(idAtm);}
        }
      }
      else if (maxbo == 2) {
        if ((atoms[idAtm] == 6)|| (atoms[idAtm] == 14)) {SetPlanar(idAtm);}
        else if ((atoms[idAtm] == 7)|| (atoms[idAtm] == 15)) {SetSP2_21(idAtm);}
        else if ((atoms[idAtm] == 8)|| (atoms[idAtm] == 16)) {SetTerminal(idAtm);}
      }
      else if (maxbo == 3) {SetSP_20(idAtm);}
    }
    for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
      if ((modVSEPR(idAtm + 1,1) == 0.0)&&((modVSEPR(idAtm + 1,2) == 0.0))) {
        if (atoms[idAtm] == 7) {
          maxbo = 0;        //number of sp2 neighbours
          for (size_t idneigh = 0; idneigh < neighbours[idAtm]; ++idneigh) {
            atomB = connectivity(idAtm + 1,idneigh + 1);
            if (atomB > 0) {if (((modVSEPR(atomB,1) == 3.0)&&(modVSEPR(atomB,2) == 0.0))||((modVSEPR(atomB,1) == 2.0)&&(modVSEPR(atomB,2) == 1.0))) {++maxbo;}}
          }
          if (maxbo > 0) {SetPlanar(idAtm);}
          else {SetSP3_31(idAtm);}
        }
      }
    }
  }
  //functions to determine atomic VSEPR codes
  void SetProton(int apos) {
    //set atom as proton classifier
    modVSEPR(apos + 1,1) = 1.0;
    modVSEPR(apos + 1,2) = -1.0;
  }
  void SetTerminal(int apos) {
    //set atom as terminal
    modVSEPR(apos + 1,1) = 1.0;
    modVSEPR(apos + 1,2) = 0.0;
  }
  void SetSP3_40(int apos) {
    //set atom as sp3 with 4 neighbours
    modVSEPR(apos + 1,1) = 4.0;
    modVSEPR(apos + 1,2) = 0.0;
  }
  void SetSP3_31(int apos) {
    //set atom as sp3 with 3 neighbours
    modVSEPR(apos + 1,1) = 3.0;
    modVSEPR(apos + 1,2) = 1.0;
  }
  void SetSP3_22(int apos) {
    //set atom as sp3 with 2 neighbours
    modVSEPR(apos + 1,1) = 2.0;
    modVSEPR(apos + 1,2) = 2.0;
  }
  void SetPlanar(int apos) {
    //set atom as planar
    modVSEPR(apos + 1,1) = 3.0;
    modVSEPR(apos + 1,2) = 0.0;
  }
  void SetSP2_21(int apos) {
    //set atom as sp2 with 2 neighbours
    modVSEPR(apos + 1,1) = 2.0;
    modVSEPR(apos + 1,2) = 1.0;
  }
  void SetSP_20(int apos) {
    //set atom as sp
    modVSEPR(apos + 1,1) = 2.0;
    modVSEPR(apos + 1,2) = 0.0;
  }
  void SetUNDEFINED(int apos) {
    //set atom as undefined
    modVSEPR(apos + 1,1) = 0.0;
    modVSEPR(apos + 1,2) = 0.0;
  }
  //for metals
  void SetOctahedral(int apos) {
    //set atom as octahedral
    modVSEPR(apos + 1,1) = 6.0;
    modVSEPR(apos + 1,2) = 0.0;
  }
  void SetSquarePyramidal(int apos) {
    //set atom as octahedral
    modVSEPR(apos + 1,1) = 5.0;
    modVSEPR(apos + 1,2) = 1.0;
  }
  void SetTrigonalBipyramidal(int apos) {
    //set atom as octahedral
    modVSEPR(apos + 1,1) = 5.0;
    modVSEPR(apos + 1,2) = 0.0;
  }
  void SetSquarePlanar(int apos) {
    //set atom as octahedral
    modVSEPR(apos + 1,1) = 4.0;
    modVSEPR(apos + 1,2) = 2.0;
  }
  void SetTShape(int apos) {
    //set atom as octahedral
    modVSEPR(apos + 1,1) = 3.0;
    modVSEPR(apos + 1,2) = 2.0;
  }
  void MultipleTerminals(int atomB, int atomC, int atomD) {
    //set several atoms to terminal?
    if ((atoms[atomB] == 8)&&(modVSEPR(atomB + 1,1) == 2.0)&&(modVSEPR(atomB + 1,2) == 2.0)&&(connectivity(atomB + 1,2) == 0)) {SetTerminal(atomB);}
    if ((atoms[atomC] == 8)&&(modVSEPR(atomC + 1,1) == 2.0)&&(modVSEPR(atomC + 1,2) == 2.0)&&(connectivity(atomC + 1,2) == 0)) {SetTerminal(atomC);}
    if ((atoms[atomD] == 8)&&(modVSEPR(atomD + 1,1) == 2.0)&&(modVSEPR(atomD + 1,2) == 2.0)&&(connectivity(atomD + 1,2) == 0)) {SetTerminal(atomD);}
  }
  void LoopTerminals(int idAtm, int neighA) {
    for (size_t idneigh = 0; idneigh < neighA; ++idneigh) {
      atomC = connectivity(idAtm + 1,idneigh + 1) - 1;
      if ((atoms[atomC] == 8)&&(modVSEPR(atomC + 1,1) == 2.0)&&(modVSEPR(atomC + 1,2) == 2.0)&&(connectivity(atomC + 1,2) == 0)) {
        SetTerminal(atomC);
        break;
      }
    }
  }
  //bond-length related functions
  double MaximalBondLength(int atom1, int atom2, int border) {
    //function determining the max limits of bond lengths with bond order
    int atomA = atom1;
    int atomB = atom2;
    if (atomA > atomB) {
      atomA = atom2;
      atomB = atom1;
    }
    double blength = 10.0;
    if ((atomA == 6)&&(atomB == 6)) {
      if (border == 2) {blength = 1.425;}
    }
    else if ((atomA == 6)&&(atomB == 7)) {
      if (border == 2) {blength = 1.375;}
    }
    else if ((atomA == 6)&&(atomB == 8)) {
      if (border == 2) {blength = 1.350;}
    }
    return blength;
  }
  double MinimalBondLength(int atom1, int atom2, int border) {
    //function determining the min limits of bond lengths with bond order
    int atomA = atom1;
    int atomB = atom2;
    if (atomA > atomB) {
      atomA = atom2;
      atomB = atom1;
    }
    double blength = 10.0;
    if ((atomA == 6)&&(atomB == 6)) {
      if (border == 1) {blength = 1.475;}
      if (border == 2) {blength = 1.250;}
    }
    else if ((atomA == 6)&&(atomB == 7)) {
      if (border == 1) {blength = 1.450;}
      if (border == 2) {blength = 1.175;}
    }
    else if ((atomA == 6)&&(atomB == 8)) {
      if (border == 1) {blength = 1.400;}
      if (border == 2) {blength = 1.150;}
    }
    else if ((atomA == 6)&&(atomB == 16)) {
      if (border == 1) {blength = 1.800;}
    }
    else if ((atomA == 7)&&(atomB == 8)) {
      if (border == 1) {blength = 1.350;}
    }
    return blength;
  }
  //angle checking functions
  int CheckOOP(int atomA, int atomB, int atomC, int atomD, double thresholdC, double thresholdN) {
    //function checking whether an atom is out-of-plane with its neighbours
    int nsp2 = 0;
    torsions[0] = fabs(radian2degree*Torsion(atomA,atomB,atomC,atomD,geom));
    torsions[1] = fabs(radian2degree*Torsion(atomA,atomB,atomD,atomC,geom));
    torsions[2] = fabs(radian2degree*Torsion(atomA,atomC,atomB,atomD,geom));
    torsions[3] = fabs(radian2degree*Torsion(atomA,atomC,atomD,atomB,geom));
    for (size_t idtor = 0; idtor < 4; ++idtor) {
      if (torsions[idtor] < 180.0) {nsp2 += (torsions[idtor] < thresholdC*(atoms[atomA - 1] == 6) + thresholdN*(atoms[atomA - 1] == 7));}
      else {nsp2 += (fabs(torsions[idtor] - 360.0) < thresholdC*(atoms[atomA - 1] == 6) + thresholdN*(atoms[atomA - 1] == 7));}
    }
    return nsp2;
  }
  //checkers
  bool IsAtomInRing(matrixE & therings, int ringsize, int iring, int atomA) {
    //function checking whether an atom is in a ring
    bool isinring = false;
    int atomB;
    for (int idatm = 0; idatm < therings.rows(); ++idatm) {
      if (idatm > ringsize) {break;}
      atomB = therings(idatm + 1,iring + 1) - 1;
      if (atomB == atomA) {
        isinring = true;
        break;
      }
    }
    return isinring;
  }
  bool ReduceUnpairedElectrons(int idAtm, double & sum) {
    UnpairedElectrons[idAtm] = 0.0;
    sum -= 0.5;
    return true;
  }
  void RemoveUnpairedElectron(int idAtm, bool checkoxygens) {
    UnpairedElectrons[idAtm] = 0.0;
    for (size_t idneigh = 0; idneigh < maxneigh; ++idneigh) {
      atomB = connectivity(idAtm + 1,idneigh + 1) - 1;
      if (atomB < 0) {break;}
      UnpairedElectronsNeighbours(idAtm + 1,idneigh + 1) = 0.0;
      if (checkoxygens) {
        if ((atoms[atomB] == 8)&&(modVSEPR(atomB + 1,1) == 1.0)&&(modVSEPR(atomB + 1,2) == 0.0)) {
          UnpairedElectrons[atomB] = 0.0;
          UnpairedElectronsNeighbours(atomB + 1,1) = 0.0;
        }
      }
      for (size_t idneigh2 = 0; idneigh2 < maxneigh; ++idneigh2) {
        atomC = connectivity(atomB + 1,idneigh2 + 1) - 1;
        if (atomC == idAtm) {
          UnpairedElectronsNeighbours(atomB + 1,idneigh2 + 1) = 0.0;
          break;
        }
      }
    }
  }
  //printing functions
  void PrintElectronPairing() {
    for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
      std::cout << idAtm + 1 << ": " << UnpairedElectrons[idAtm] << "    ";
      for (size_t idneigh = 0; idneigh < 5; ++idneigh) {
        std::cout << UnpairedElectronsNeighbours(idAtm + 1,idneigh + 1) << " ";
      }
      std::cout << formalcharges[idAtm] << std::endl;
    }
  }
  //counting functions
  void AngleBasedHybridisationNeighbours_3center(std::vector<double> & maxangle, int & nsp2, int & nsp3, int & lessthansp3, double threshang) {
    //function getting the number of neighbours of a certain hybridization based on angular arrangements
    nsp2 = 0;
    nsp3 = 0;
    lessthansp3 = 0;
    for (size_t idx = 0; idx < 3; ++idx) {
      nsp2 += ((maxangle[idx] > 120.0 - threshang)&&(maxangle[idx] < 120.0 + 3.0*threshang));
      nsp3 += ((maxangle[idx] > 109.5 - 2.0*threshang)&&(109.5 < 109.5 + threshang));
      lessthansp3 += (maxangle[idx] > 109.5 - 2.0*threshang);
    }
  }
  int Matchcounter(matrixE & therings, int idring, int jdring, int niring, int njring) {
    //function counting how many atoms are common between two rings
    int match = 0;
    for (size_t idAtm = 0; idAtm < niring; ++idAtm) {
      for (size_t idBtm = 0; idBtm < njring; ++idBtm) {
        if (therings(idAtm + 1,idring + 1) == therings(idBtm + 1,jdring + 1)) {
          ++match;
          break;
        }
      }
    }
    return match;
  }
  int CountSP2Neighbours(int atomA, int nneighA) {
    //function that counts sp2 neighbours
    int nsp2_neigh = 0;            //count the number of sp2 atoms in the neighbourhood
    for (size_t idneigh = 0; idneigh < nneighA; ++idneigh) {
      atomB = connectivity(atomA + 1,idneigh + 1) - 1;
      if ((modVSEPR(atomB + 1,1) == 3.0)&&(modVSEPR(atomB + 1,2) == 0.0)) {++nsp2_neigh;}
      else if ((modVSEPR(atomB + 1,1) == 2.0)&&(modVSEPR(atomB + 1,2) == 1.0)) {++nsp2_neigh;}
      else if ((modVSEPR(atomB + 1,1) == 1.0)&&(modVSEPR(atomB + 1,2) == 0.0)) {++nsp2_neigh;}
    }
    return nsp2_neigh;
  }
  //other setters
  void set3Neighbours(int atomA, int & atomB, int & atomC, int & atomD) {
    atomB = connectivity(atomA + 1,1) - 1;
    atomC = connectivity(atomA + 1,2) - 1;
    atomD = connectivity(atomA + 1,3) - 1;
  }
  //copy functions
  void CopyPath(std::vector<int> & connectionpath, matrix<int> & paths, int ifirstatm, int pathlength, int path2save, int atomsize) {
    connectionpath.resize(pathlength);
    connectionpath[0] = ifirstatm + 1;
    for (size_t iatm = 0; iatm < atomsize; ++iatm) {
      if (paths(path2save + 1,iatm + 1) == 0) {break;}
      connectionpath[iatm + 1] = paths(path2save + 1,iatm + 1);
    }
  }
  //ring functions
  void DeleteRing(matrixE & therings, int idring, int ring_size) {
    //function that deletes a ring
    for (size_t idAtm = 0; idAtm < ring_size; ++idAtm) {
      therings(idAtm + 1,idring + 1) = 0.0;
    }
  }
  int RingSize(int iring) {
    //determine the ring size
    int ringsize = therings.rows() - 1;
    for (size_t idAtm = 0; idAtm < therings.rows(); ++idAtm) {
      if (therings(idAtm + 1,iring + 1) == 0) {
        ringsize = idAtm - 1;
        break;
      }
    }
    return ringsize;
  }
  void RingFinder() {
    //function to identify rings
    int nneigh;
    int totalpaths = 0;
    int npaths = 3*Natoms + 20;
    int currentpos;
    int firstadd;
    int iBtm;
    int iCtm;
    int prevatm;
    int openpaths = 0;
    int iring = 0;
    int nrings = 0;
    int largestring = 0;
    int startpos;
    int thelastatm;
    int thelastpos;
    bool writeatm;
    std::vector<int> active(npaths,0);
    std::vector<int> hasring(npaths,0);
    matrix<int> paths(npaths,Natoms + 20);
    matrixE therings_temp(1,1);
    for (size_t idneigh = 0; idneigh < neighbours[0] ; ++idneigh) {
      paths(idneigh + 1,1) = 1;
      paths(idneigh + 1,2) = connectivity(1,idneigh + 1);
      active[idneigh] = 1;
      ++totalpaths;
    }
    currentpos = 2;
nextatom:
    if (totalpaths > npaths - 20) {
      npaths = 2*totalpaths;
      active.resize(npaths);
      hasring.resize(npaths);
      paths.resize(npaths,Natoms + 20);
    }
    //add new paths
    for (size_t ipath = 0; ipath < npaths; ++ipath) {
      if (active[ipath]) {
        iBtm = paths(ipath + 1,currentpos);
        prevatm = paths(ipath + 1,currentpos - 1);
        if (iBtm == 0) {
          active[ipath] = 0;
          continue;
        }
        if (paths(ipath + 1,currentpos + 1) != 0) {continue;}
        firstadd = 0;
        for (size_t idneigh = 0; idneigh < maxneigh; ++idneigh) {
          iCtm = connectivity(iBtm,idneigh + 1);
          if (iCtm == prevatm) {continue;}
          else if (iCtm - 1 < 0) {break;}
          if (!firstadd) {
            firstadd = 1;
            paths(ipath + 1,currentpos + 1) = iCtm;
          }
          else {
            active[totalpaths] = 1;
            ++totalpaths;
            //copy path
            for (size_t idx = 0; idx < currentpos + 1; ++idx) {
              paths(totalpaths,idx + 1) = paths(ipath + 1,idx + 1);
            }
            paths(totalpaths,currentpos + 1) = iCtm;
          }
        }
      }
    }
    //close paths?
    openpaths = 0;
    for (size_t ipath = 0; ipath < npaths; ++ipath) {
      if (active[ipath]) {
        if (neighbours[paths(ipath + 1,currentpos + 1) - 1] == 1) {
          active[ipath] = 0;
          continue;
        }
        ++openpaths;
        for (size_t jpath = 0; jpath < currentpos; ++jpath) {
          if (paths(ipath + 1,jpath + 1) == paths(ipath + 1,currentpos + 1)) {
            active[ipath] = 0;
            hasring[ipath] = 1;
            --openpaths;
            break;
          }
        }
      }
    }
    if (openpaths) {
      ++currentpos;
      goto nextatom;
    }
    for (size_t ipath = 0; ipath < npaths; ++ipath) {
      nrings += hasring[ipath];
    }
    therings_temp.resize(Natoms,nrings);
    active.resize(nrings);
    for (size_t ipath = 0; ipath < npaths; ++ipath) {
      if (hasring[ipath]) {
        for (size_t jpath = 1; jpath < npaths; ++jpath) {
          if (paths(ipath + 1,jpath + 1) == 0) {
            thelastatm = paths(ipath + 1,jpath);
            thelastpos = jpath;
            break;
          }
        }
        if (thelastpos < 3) {
          hasring[ipath] = 0;
          continue;
        }
        startpos = -1;
        for (size_t jpath = 0; jpath < thelastpos - 1; ++jpath) {
          if (paths(ipath + 1,jpath + 1) == thelastatm) {
            startpos = jpath;
            if (thelastpos - jpath > largestring) {largestring = thelastpos - jpath;}
            break;
          }
        }
        if (startpos > -1) {
          ++iring;
          iBtm = 0;
          active[iring - 1] = thelastpos - startpos - 1;
          for (size_t jpath = startpos; jpath < thelastpos - 1; ++jpath) {
            ++iBtm;
            therings_temp(iBtm,iring) = paths(ipath + 1,jpath + 1);
          }
        }
      }
    }
    for (size_t idring = 0; idring < nrings; ++idring) {
      if ((active[idring] == 0)&&(therings_temp(1,idring + 1) == 0)) {continue;}
      for (size_t jdring = idring + 1; jdring < nrings; ++jdring) {
        if ((active[jdring] == 0)&&(therings_temp(1,jdring + 1) == 0)) {continue;}
        if (active[idring] == active[jdring]) {
          //then the atom lists must match exactly
          iBtm = Matchcounter(therings_temp,idring,jdring,active[idring],active[jdring]);
          if (iBtm == active[idring]) {
            //remove ring j, it is repeated
            DeleteRing(therings_temp,jdring,active[jdring]);
            active[jdring] = 0;
          }
        }
        else {
          //see whether they are subrings
          if (active[idring] < active[jdring]) {
            iBtm = Matchcounter(therings_temp,idring,jdring,active[idring],active[jdring]);
            if (iBtm == active[idring]) {
              DeleteRing(therings_temp,jdring,active[jdring]);
              active[jdring] = 0;
            }
          }
          else {
            iBtm = Matchcounter(therings_temp,jdring,idring,active[jdring],active[idring]);
            if (iBtm == active[jdring]) {
              DeleteRing(therings_temp,idring,active[idring]);
              active[idring] = 0;
              break;
            }
          }
        }
      }
    }
    therings.resize(largestring,nrings);
    iring = 0;
    for (size_t idring = 0; idring < nrings; ++idring) {
      if (therings_temp(1,idring + 1) > 0.0) {
        //copy
        for (size_t idAtm = 0; idAtm < active[idring]; ++idAtm) {
          therings(idAtm + 1,iring + 1) = therings_temp(idAtm + 1,idring + 1);
        }
        ++iring;
      }
    }
    therings.resize(largestring,iring);
  }
  void CountElectronsAndCharges() {
    //function that counts electrons and charges based on a certain VSEPR classification
    for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
      formalcharges[idAtm] = 0;
      UnpairedElectrons[idAtm] = 0.0;
      nO = 0;         //terminal oxygens
      nN = 0;         //planar nitrogens
      nC = 0;         //sp3 carbons
      for (size_t idneigh = 0; idneigh < maxneigh; ++idneigh) {
        atomB = connectivity(idAtm + 1,idneigh + 1) - 1;
        if (atomB < 0) {break;}
        if ((atoms[atomB] == 8)&&(modVSEPR(atomB + 1,1) == 1.0)&&(modVSEPR(atomB + 1,2) == 0.0)) {++nO;}
        else if ((atoms[atomB] == 7)&&(modVSEPR(atomB + 1,1) == 3.0)&&(modVSEPR(atomB + 1,2) == 0.0)) {++nN;}
        else if ((atoms[atomB] == 6)&&(modVSEPR(atomB + 1,1) == 4.0)&&(modVSEPR(atomB + 1,2) == 0.0)) {++nC;}
        if ((atoms[atomB] == 6)||(atoms[atomB] == 14)) {
          if ((modVSEPR(atomB + 1,1) == 3.0)&&(modVSEPR(atomB + 1,2) == 0.0)) {UnpairedElectronsNeighbours(idAtm + 1,idneigh + 1) = 0.5;}
          else if ((modVSEPR(atomB + 1,1) == 2.0)&&(modVSEPR(atomB + 1,2) == 0.0)) {UnpairedElectronsNeighbours(idAtm + 1,idneigh + 1) = 1.0;}
        }
        else if ((atoms[atomB] == 7)||(atoms[atomB] == 15)) {
          if ((modVSEPR(atomB + 1,1) == 2.0)&&(modVSEPR(atomB + 1,2) == 1.0)) {UnpairedElectronsNeighbours(idAtm + 1,idneigh + 1) = 0.5;}
          else if ((modVSEPR(atomB + 1,1) == 2.0)&&(modVSEPR(atomB + 1,2) == 0.0)) {UnpairedElectronsNeighbours(idAtm + 1,idneigh + 1) = 1.0;}
        }
        else if ((atoms[atomB] == 8)||(atoms[atomB] == 16)||(atoms[atomB] == 34)) {
          if ((modVSEPR(atomB + 1,1) == 1.0)&&(modVSEPR(atomB + 1,2) == 0.0)) {UnpairedElectronsNeighbours(idAtm + 1,idneigh + 1) = 0.5;}
        }
      }
      if ((atoms[idAtm] == 6)||(atoms[idAtm] == 14)) {
        if ((modVSEPR(idAtm + 1,1) == 3.0)&&(modVSEPR(idAtm + 1,2) == 0.0)) {
          UnpairedElectrons[idAtm] = 0.5;
          //carboxylate
          if (nO == 2) {formalcharges[idAtm] = -1;}
          //guanidino
          else if ((nO == 0)&&(nN == 3)) {formalcharges[idAtm] = 1;}
          //iminium
          else if ((nO == 0)&&(nN == 2)&&(nC == 1)) {formalcharges[idAtm] = 1;}
        }
        else if ((modVSEPR(idAtm + 1,1) == 2.0)&&(modVSEPR(idAtm + 1,2) == 0.0)) {UnpairedElectrons[idAtm] = 1.0;}
      }
      else if ((atoms[idAtm] == 7)||(atoms[idAtm] == 15)) {
        if ((modVSEPR(idAtm + 1,1) == 2.0)&&(modVSEPR(idAtm + 1,2) == 1.0)) {UnpairedElectrons[idAtm] = 0.5;}
        else if ((modVSEPR(idAtm + 1,1) == 2.0)&&(modVSEPR(idAtm + 1,2) == 0.0)) {UnpairedElectrons[idAtm] = 1.0;}
        else if ((atoms[idAtm] == 7)&&(modVSEPR(idAtm + 1,1) == 4.0)&&(modVSEPR(idAtm + 1,2) == 0.0)) {formalcharges[idAtm] = 1;}
        else if ((atoms[idAtm] == 15)&&(modVSEPR(idAtm + 1,1) == 4.0)&&(modVSEPR(idAtm + 1,2) == 0.0)) {formalcharges[idAtm] = - (nO - 1);}
      }
      else if ((atoms[idAtm] == 8)||(atoms[idAtm] == 16)||(atoms[idAtm] == 34)) {
        if ((modVSEPR(idAtm + 1,1) == 1.0)&&(modVSEPR(idAtm + 1,2) == 0.0)) {UnpairedElectrons[idAtm] = 0.5;}
        else if ((atoms[idAtm] == 16)&&(modVSEPR(idAtm + 1,1) == 4.0)&&(modVSEPR(idAtm + 1,2) == 0.0)) {formalcharges[idAtm] = - (nO - 2);}
      }
    }
    //fixes
    for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
      nO = 0;         //terminal oxygens
      for (size_t idneigh = 0; idneigh < maxneigh; ++idneigh) {
        atomB = connectivity(idAtm + 1,idneigh + 1) - 1;
        if (atomB < 0) {break;}
        if ((atoms[atomB] == 8)&&(modVSEPR(atomB + 1,1) == 1.0)&&(modVSEPR(atomB + 1,2) == 0.0)) {++nO;}
      }
      //carboxylates are to ignore
      if ((atoms[idAtm] == 6)&&(nO == 2)) {RemoveUnpairedElectron(idAtm,true);}
      //nitros too
      if ((atoms[idAtm] == 7)&&(nO == 2)&&(modVSEPR(idAtm + 1,1) == 3.0)&&(modVSEPR(idAtm + 1,2) == 0.0)) {RemoveUnpairedElectron(idAtm,true);}
      //so are phosphates
      if (((atoms[idAtm] == 15)||(atoms[idAtm] == 16))&&(nO > 0)) {RemoveUnpairedElectron(idAtm,true);}
    }
  }
  int sp2Connected(int idAtm, int idBtm, std::vector<int> & connectionpath) {
    //function that determines whether two atoms are connected in double bonds
    int idCtm;
    int idDtm;
    int npaths = Natoms + 30;
    int openpaths = 0;
    int themagicalpath = -1;
    int thelastatm;
    int firstpath;
    int counter = 0;
    bool repeated;
    matrix<int> paths(npaths,Natoms + 10);
    //allocate paths
    for (size_t idneigh = 0; idneigh < maxneigh; ++idneigh) {
      idCtm = connectivity(idAtm + 1,idneigh + 1) - 1;
      if (idCtm < 0) {break;}
      if (((modVSEPR(idCtm + 1,1) == 3.0)&&(modVSEPR(idCtm + 1,2) == 0.0))||((modVSEPR(idCtm + 1,1) == 2.0)&&(modVSEPR(idCtm + 1,2) == 1.0))) {
        paths(openpaths + 1,1) = idCtm + 1;
        ++openpaths;
      }
      if (idCtm == idBtm) {themagicalpath = openpaths - 1;}
    }
redoiteration:
    ++counter;
    openpaths = 0;
    //create paths as long as there are sp2 atoms around
    for (size_t ipath = 0; ipath < paths.rows(); ++ipath) {
      if (paths(ipath + 1,1) != 0) {
        ++openpaths;
        thelastatm = -1;
        for (size_t iatm = 1; iatm < Natoms + 10; ++iatm) {
          if (paths(ipath + 1,iatm + 1) == 0) {
            thelastatm = iatm;
            idCtm = paths(ipath + 1,iatm);
            break;
          }
        }
        if (thelastatm > 0) {
          firstpath = 0;
          for (size_t idneigh = 0; idneigh < maxneigh; ++idneigh) {
            idDtm = connectivity(idCtm,idneigh + 1) - 1;
            if (idDtm < 0) {break;}
            else if (idDtm == idBtm) {
              themagicalpath = ipath;
              paths(ipath + 1,thelastatm + 1) = idDtm + 1;
              break;
            }
            if (((modVSEPR(idDtm + 1,1) == 3.0)&&(modVSEPR(idDtm + 1,2) == 0.0))||((modVSEPR(idDtm + 1,1) == 2.0)&&(modVSEPR(idDtm + 1,2) == 1.0))) {
              if (firstpath == 0) {
                firstpath = 1;
                paths(ipath + 1,thelastatm + 1) = idDtm + 1;
              }
              else {
                for (size_t jpath = 0; jpath < paths.rows(); ++jpath) {
                  if (paths(jpath + 1,1) == 0) {
                    //start by "zeroing"
                    repeated = false;
                    for (size_t icol = 0; icol < paths.cols(); ++icol) {
                      paths(jpath + 1,icol + 1) = paths(ipath + 1,icol + 1);
                      if (paths(jpath + 1,icol + 1) == idDtm + 1) {
                        repeated = true;
                        break;
                      }
                    }
                    if (repeated) {paths(jpath + 1,1) = 0;}
                    else {paths(jpath + 1,thelastatm + 1) = idDtm + 1;}
                    break;
                  }
                }
              }
            }
          }
          idDtm = paths(ipath + 1,thelastatm + 1);
          for (size_t icol = 0; icol < thelastatm; ++icol) {
            if (paths(ipath + 1,icol + 1) == idDtm) {
              paths(ipath + 1,1) = 0;
              break;
            }
          }
        }
      }
      if (themagicalpath >= 0) {break;}
    }
    if ((openpaths > 0)&&(themagicalpath < 0)&&(counter < 100)) {goto redoiteration;}
    if (themagicalpath != -1) {
      int pathlength = 1;
      int npl = 0;
      for (size_t iatm = 0; iatm < Natoms + 10; ++iatm) {
        if (paths(themagicalpath + 1,iatm + 1) == 0) {break;}
        ++pathlength;
        npl += (modVSEPR(paths(themagicalpath + 1,iatm + 1),1) == 3.0)*(modVSEPR(paths(themagicalpath + 1,iatm + 1),2) == 0.0)*(atoms[paths(themagicalpath + 1,iatm + 1) - 1] == 7);
      }
      if ((pathlength - npl)%2 == 0) {CopyPath(connectionpath,paths,idAtm,pathlength,themagicalpath,Natoms + 10);}
      else if (npl != 0) {
        CopyPath(connectionpath,paths,idAtm,pathlength,themagicalpath,Natoms + 10);
        themagicalpath = -2;
      }
      else if (pathlength != 0) {
        CopyPath(connectionpath,paths,idAtm,pathlength,themagicalpath,Natoms + 10);
        themagicalpath = -3;
      }
      else {themagicalpath = -1;}
    }
    return themagicalpath;
  }
  void CoupleElectrons(std::vector<int> & flatcenter) {
    //function that couples electrons to leave only the radicals
    int namides;
    int nconnectivity;
    int isringatom;
    int matchingshell;
    int nrings = therings.cols();
    int maxnatm = therings.rows();
    int theringend;
    int firstrun = 0;
    int nflatatoms;
    int total_nsp2_double_count;
    double sum_unpaired_electrons = 0.0;
    std::vector<int> doublebondcountsINRING(8);
    std::vector<int> doublebondcountsOFFRING(8);
    std::vector<int> connectionpath;
    //go over rings first
doringloop:
    for (size_t iring = 0; iring < nrings; ++iring) {
      theringend = therings.rows();
      nflatatoms = 0;
      //get the size of the ring and the number of flat atoms in it
      for (size_t idAtm = 0; idAtm < therings.rows(); ++idAtm) {
        atomB = therings(idAtm + 1,iring + 1);
        if (atomB == 0) {
          theringend = idAtm;
          break;
        }
        nflatatoms += (flatcenter[atomB - 1] == 1);
        if ((flatcenter[atomB - 1] == 0)&&(connectivity(atomB,3) == 0)) {++nflatatoms;}
      }
      if (firstrun == 0) {
        if (nflatatoms != theringend) {continue;}
        else if (theringend > 6) {continue;}
      }
      else {if (nflatatoms == theringend) {continue;}}
      std::cout << theringend << std::endl;
      doublebondcountsINRING.resize(theringend);
      doublebondcountsOFFRING.resize(theringend);
      total_nsp2_double_count = 0;
      for (size_t idAtm = 0; idAtm < theringend; ++idAtm) {
        atomA = therings(idAtm + 1,iring + 1) - 1;
        doublebondcountsINRING[idAtm] = 0;
        doublebondcountsOFFRING[idAtm] = 0;
        if (UnpairedElectrons[atomA] == 0.0) {
          doublebondcountsINRING[idAtm] = -1;
          doublebondcountsOFFRING[idAtm] = -1;
          //neighbours do not matter
          for (size_t idneigh = 0; idneigh < maxneigh; ++idneigh) {
            atomB = connectivity(atomA + 1,idneigh + 1) - 1;
            if (atomB < 0) {break;}
            UnpairedElectronsNeighbours(atomA + 1,idneigh + 1) = 0.0;
            for (size_t idneigh2 = 0; idneigh2 < maxneigh; ++idneigh2) {
              atomC = connectivity(atomB + 1,idneigh + 1) - 1;
              if (atomC == atomA) {
                UnpairedElectronsNeighbours(atomB + 1,idneigh + 1) = 0.0;
                break;
              }
            }
          }
        }
        else {
          for (size_t idneigh = 0; idneigh < maxneigh; ++idneigh) {
            atomB = connectivity(atomA + 1,idneigh + 1) - 1;
            if (atomB < 0) {break;}
            if (UnpairedElectrons[atomB] == 0.0) {continue;}
            isringatom = 0;
            for (size_t idBtm = 0; idBtm < theringend; ++idBtm) {
              if (therings(idBtm + 1,iring + 1) - 1 == atomB) {
                isringatom = 1;
                break;
              }
            }
            if (isringatom == 1) {++doublebondcountsINRING[idAtm];}
            else {++doublebondcountsOFFRING[idAtm];}
          }
        }
        if (doublebondcountsINRING[idAtm] == 1) {
          //block this double bond
          for (size_t idneigh = 0; idneigh < maxneigh; ++idneigh) {
            atomB = connectivity(atomA + 1,idneigh + 1) - 1;
            if (atomB < 0) {break;}
            if (UnpairedElectrons[atomB] == 0.0) {continue;}
            isringatom = 0;
            for (size_t idBtm = 0; idBtm < theringend; ++idBtm) {
              if (therings(idBtm + 1,iring + 1) - 1 == atomB) {
                isringatom = 1;
                break;
              }
            }
            if (isringatom == 1) {
              //found the neighbour
              UnpairedElectrons[atomA] = 0.0;
              UnpairedElectrons[atomB] = 0.0;
              for (size_t idneigh2 = 0; idneigh2 < maxneigh; ++idneigh2) {
                atomC = connectivity(atomA + 1,idneigh2 + 1) - 1;
                atomD = connectivity(atomB + 1,idneigh2 + 1) - 1;
                if ((atomC < 0)&&(atomD < 0)) {break;}
                if (atomC >= 0) {UnpairedElectronsNeighbours(atomA + 1,idneigh2 + 1) = 0.0;}
                if (atomD >= 0) {UnpairedElectronsNeighbours(atomB + 1,idneigh2 + 1) = 0.0;}
                for (size_t idneigh3 = 0; idneigh3 < maxneigh; ++idneigh3) {
                  atomE = connectivity(atomC + 1,idneigh3 + 1) - 1;
                  if (atomE == atomA) {
                    UnpairedElectronsNeighbours(atomC + 1,idneigh3 + 1) = 0.0;
                    break;
                  }
                }
                for (size_t idneigh3 = 0; idneigh3 < maxneigh; ++idneigh3) {
                  atomE = connectivity(atomD + 1,idneigh3 + 1) - 1;
                  if (atomE == atomB) {
                    UnpairedElectronsNeighbours(atomD + 1,idneigh3 + 1) = 0.0;
                    break;
                  }
                }
              }
            }
          }
          --iring;
          std::cout << "  " << atomA + 1 << " " << doublebondcountsINRING[idAtm] << " " << doublebondcountsOFFRING[idAtm] << std::endl;
          break;
        }
        else if (doublebondcountsINRING[idAtm] == 2) {++total_nsp2_double_count;}
        std::cout << "  " << atomA + 1 << " " << doublebondcountsINRING[idAtm] << " " << doublebondcountsOFFRING[idAtm] << std::endl;
      }
      if (total_nsp2_double_count == theringend) {
        if (theringend%2 == 0) {
          //aromatic ring, couple all
          for (size_t idAtm = 0; idAtm < theringend; ++idAtm) {
            RemoveUnpairedElectron(therings(idAtm + 1,iring + 1) - 1,false);
          }
        }
        else {
          matchingshell = 0;             //whether a nitrogen was found that can be protonated
          for (size_t idAtm = 0; idAtm < theringend; ++idAtm) {
            atomA = therings(idAtm + 1,iring + 1) - 1;
            if (connectivity(atomA + 1,3) != 0) {
              for (size_t idneigh = 0; idneigh < maxneigh; ++idneigh) {
                atomB = connectivity(atomA + 1,idneigh + 1) - 1;
                if (atomB < 0) {break;}
                if ((atoms[atomB] == 7)&&(modVSEPR(atomB + 1,1) == 2.0)&&(modVSEPR(atomB + 1,2) == 1.0)) {
                  for (size_t idBtm = 0; idBtm < theringend; ++idBtm) {
                    if (therings(idBtm + 1,iring + 1) - 1 == atomB) {
                      //found a suitable nitrogen to protonate, get rid of it from radical counting
                      matchingshell = 1;
                      SetPlanar(atomB);
                      RemoveUnpairedElectron(atomB,false);
                      --iring;
                      break;
                    }
                  }
                }
                if (matchingshell == 1) {break;}
              }
            }
            if (matchingshell == 1) {break;}
          }
          if (matchingshell == 0) {
            //we kind of need to decide here what to do
          }
        }
      }
    }
    if (firstrun == 0) {
      firstrun = 1;
      goto doringloop;
    }
    //look for specific cases
    for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
      if (UnpairedElectrons[idAtm] != 0.0) {
        //phenols and nitroso
        if ((atoms[idAtm] == 8)&&(connectivity(idAtm + 1,2) == 0)) {
          if ((modVSEPR(idAtm + 1,1) == 1.0)&&(modVSEPR(idAtm + 1,2) == 0.0)) {
            atomB = connectivity(idAtm + 1,1) - 1;
            if ((atoms[atomB] == 6)&&(modVSEPR(atomB + 1,1) == 3.0)&&(modVSEPR(atomB + 1,2) == 0.0)) {
              total_nsp2_double_count = 0;                                     //number of sp2 neighbours
              for (size_t idneigh = 0; idneigh < maxneigh; ++idneigh) {
                atomC = connectivity(atomB + 1,idneigh + 1) - 1;
                if (atomC < 0) {break;}
                if ((atoms[atomC] == 6)&&(modVSEPR(atomC + 1,1) == 3.0)&&(modVSEPR(atomC + 1,2) == 0.0)) {++total_nsp2_double_count;}
                else if ((atoms[atomC] == 7)&&(modVSEPR(atomC + 1,1) == 2.0)&&(modVSEPR(atomC + 1,2) == 1.0)) {++total_nsp2_double_count;}
              }
              if (total_nsp2_double_count == 2) {
                SetSP3_22(idAtm);
                UnpairedElectrons[idAtm] = 0.0;
              }
            }
          }
          else {UnpairedElectrons[idAtm] = 0.0;}
        }
        //indoles that were processed after the 6 member ring
        if ((atoms[idAtm] == 7)&&(modVSEPR(idAtm + 1,1) == 2.0)&&(modVSEPR(idAtm + 1,2) == 1.0)) {
          //check if the nitrogen is in ring
          matchingshell = 0;
          for (size_t iring = 0; iring < nrings; ++iring) {
            for (size_t idBtm = 0; idBtm < therings.rows(); ++idBtm) {
              atomB = therings(idBtm + 1,iring + 1) - 1;
              if (atomB < 0) {break;}
              else if (atomB == idAtm) {
                SetPlanar(idAtm);
                UnpairedElectrons[idAtm] = 0.0;
                matchingshell = 1;
                break;
              }
            }
            if (matchingshell == 1) {break;}
          }
          if (matchingshell == 0) {
            //we kind of need to decide what to do here
          }
        }
        //find isolated insaturations
        matchingshell = 0;     //isolated insaturations
        isringatom = 0;        //planar nitrogens
        for (size_t idneigh = 0; idneigh < maxneigh; ++idneigh) {
          atomB = connectivity(idAtm + 1,idneigh + 1) - 1;
          if (atomB < 0) {break;}
          if (UnpairedElectrons[idAtm] == UnpairedElectrons[atomB]) {++matchingshell;}
          if ((atoms[atomB] == 7)&&(modVSEPR(atomB + 1,1) == 3.0)&&(modVSEPR(atomB + 1,2) == 0.0)) {++isringatom;}
        }
        if ((matchingshell == 0)&&(isringatom > 1)) {
          //give it charge
          formalcharges[idAtm] = 1;
          UnpairedElectrons[idAtm] = 0.0;
          for (size_t idneigh = 0; idneigh < maxneigh; ++idneigh) {
            atomB = connectivity(idAtm + 1,idneigh + 1) - 1;
            if (atomB < 0) {break;}
            if ((atoms[atomB] == 7)&&(modVSEPR(atomB + 1,1) == 3.0)&&(modVSEPR(atomB + 1,2) == 0.0)) {
              for (size_t idneigh2 = 0; idneigh2 < maxneigh; ++idneigh2) {
                atomC = connectivity(atomB + 1,idneigh2 + 1) - 1;
                if (idAtm == atomC) {
                  UnpairedElectronsNeighbours(atomB + 1,idneigh2 + 1) = 0.0;
                  break;
                }
              }
            }
          }
        }
        if (matchingshell == 1) {
          for (size_t idneigh = 0; idneigh < maxneigh; ++idneigh) {
            atomB = connectivity(idAtm + 1,idneigh + 1) - 1;
            if (atomB < 0) {break;}
            for (size_t idneigh2 = 0; idneigh2 < maxneigh; ++idneigh2) {
              atomC = connectivity(atomB + 1,idneigh2 + 1) - 1;
              if (atomC == idAtm) {
                UnpairedElectronsNeighbours(atomB + 1,idneigh2 + 1) = 0.0;
                break;
              }
            }
            if (UnpairedElectrons[idAtm] == UnpairedElectrons[atomB]) {
              UnpairedElectronsNeighbours(idAtm + 1,idneigh + 1) = 0.0;
              RemoveUnpairedElectron(atomB,false);
            }
          }
          UnpairedElectrons[idAtm] = 0.0;
        }
      }
      sum_unpaired_electrons += UnpairedElectrons[idAtm];
    }
    if (sum_unpaired_electrons > 0.0) {
      std::cout << "still undone" << std::endl;
      for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
        if (UnpairedElectrons[idAtm] > 0.0) {
          matchingshell = 0;
          theringend = 0;        //number of neighbours
          sum_unpaired_electrons = 0.0;
          atomE = 0;             //neighboring chalcogens
          nflatatoms = 0;        //sp2 nitrogens
          atomD = 0;             //planar nitrogens nearby
          namides = 0;
          nconnectivity = 0;
          for (size_t idneigh = 0; idneigh < maxneigh; ++idneigh) {
            atomB = connectivity(idAtm + 1,idneigh + 1) - 1;
            if (atomB < 0) {break;}
            sum_unpaired_electrons += UnpairedElectrons[atomB];
            ++theringend;
            atomE += (atoms[atomB] == 8) + (atoms[atomB] == 16);
            nflatatoms += (atoms[atomB] == 8)*(modVSEPR(atomB + 1,1) == 2.0)*(modVSEPR(atomB + 1,2) == 2.0)*(connectivity(atomB + 1,2) == 0);
            atomD += (atoms[atomB] == 7)*(modVSEPR(atomB + 1,1) == 3.0)*(modVSEPR(atomB + 1,2) == 0.0);
            if ((atoms[atomB] == 7)&&(modVSEPR(atomB + 1,1) == 3.0)&&(modVSEPR(atomB + 1,2) == 0.0)) {nconnectivity = connectivity(atomB + 1,3);}
            namides += (flatcenter[atomB] == 3);
          }
          if (((theringend == 1)||((theringend == 2)&&(atomE == 1)))&&(sum_unpaired_electrons == 0.0)) {
            std::cout << "  1 " << idAtm + 1 << " " << atoms[idAtm] << std::endl;
            //nothing around, make it sp3
            matchingshell = 1;
            UnpairedElectrons[idAtm] = 0.0;
            if ((atoms[idAtm] == 6)||(atoms[idAtm] == 14)) {SetSP3_40(idAtm);}
            else if ((atoms[idAtm] == 7)||(atoms[idAtm] == 15)) {SetSP3_31(idAtm);}
            else if ((atoms[idAtm] == 8)||(atoms[idAtm] == 16)) {SetSP3_22(idAtm);}
          }
          else if ((atomD > 0)&&((namides != atomD)||(namides > 1))) {
            std::cout << "  2" << std::endl;
            if (nflatatoms == 1) {
              for (size_t idneigh = 0; idneigh < maxneigh; ++idneigh) {
                atomB = connectivity(idAtm + 1,idneigh + 1) - 1;
                if (atomB < 0) {break;}
                if ((atoms[atomB] == 8)&&(modVSEPR(atomB + 1,1) == 2.0)&&(modVSEPR(atomB + 1,2) == 2.0)&&(connectivity(atomB + 1,2) == 0)) {
                  SetTerminal(atomB);
                  UnpairedElectrons[idAtm] = 0.0;
                  matchingshell = 1;
                  break;
                }
              }
            }
            else if (nconnectivity == 0) {
              //transform nitrogen to 2,1
              namides = 0;
              for (size_t idneigh = 0; idneigh < maxneigh; ++idneigh) {
                atomB = connectivity(idAtm + 1,idneigh + 1) - 1;
                if (atomB < 0) {break;}
                if (flatcenter[atomB] == 3) {continue;}
                if ((atoms[atomB] == 7)&&(namides == 0)&&(modVSEPR(atomB + 1,1) == 3.0)&&(modVSEPR(atomB + 1,2) == 0.0)&&(connectivity(atomB + 1,3) == 0)) {
                  SetSP2_21(atomB);
                  UnpairedElectrons[idAtm] = 0.0;
                  ++namides;
                  matchingshell = 1;
                }
              }
            }
          }
          if (matchingshell == 0) {
            std::cout << "  3" << std::endl;
            for (size_t idBtm = 0; idBtm < Natoms; ++idBtm) {
              if (idBtm == idAtm) {continue;}
              if ((UnpairedElectrons[idBtm] > 0.0)||((atoms[idBtm] == 8)&&(connectivity(idBtm + 1,2) == 0)&&(modVSEPR(idBtm + 1,1) == 2.0)&&(modVSEPR(idBtm + 1,2) == 2.0))) {
                matchingshell = sp2Connected(idAtm,idBtm,connectionpath);
                std::cout << "correction? " << matchingshell << std::endl;
                if (matchingshell >= 0) {
                  matchingshell = connectionpath.size();
                  for (size_t iconn = 0; iconn < matchingshell; ++iconn) {
                    atomB = connectionpath[iconn] - 1;
                    UnpairedElectrons[atomB] = 0.0;
                  }
                  if ((atoms[idBtm] == 8)&&(connectivity(idBtm + 1,2) == 0)&&(modVSEPR(idBtm + 1,1) == 2.0)&&(modVSEPR(idBtm + 1,2) == 2.0)) {SetTerminal(idBtm);}
                  break;
                }
                else if (matchingshell == -2) {
                  //find the nitrogen and "planarise it"
                  for (size_t iconn = 0; iconn < connectionpath.size(); ++iconn) {
                    atomB = connectionpath[iconn] - 1;
                    UnpairedElectrons[atomB] = 0.0;
                    if ((atoms[atomB] == 7)&&(modVSEPR(atomB + 1,1) == 2.0)&&(modVSEPR(atomB + 1,2) == 1.0)) {SetPlanar(atomB);}
                  }
                }
                else if (matchingshell == -3) {
                  //the connection is there, but the number of sp2 atoms is uneven, so one goes to sp3.
                  if ((atoms[idAtm] == 8)&&(connectivity(idAtm + 1,2) == 0)) {SetCarbonyl(idAtm,idBtm,nrings);}
                  else if ((atoms[idBtm] == 8)&&(connectivity(idBtm + 1,2) == 0)) {SetCarbonyl(idBtm,idAtm,nrings);}
                }
              }
            }
          }
          //if there is still no solution, then force sp3
          if ((theringend == 2)&&(sum_unpaired_electrons == 0.0)&&(UnpairedElectrons[idAtm] > 0.0)) {
            UnpairedElectrons[idAtm] = 0.0;
            if ((atoms[idAtm] == 6)||(atoms[idAtm] == 14)) {SetSP3_40(idAtm);}
            else if ((atoms[idAtm] == 7)||(atoms[idAtm] == 15)) {SetSP3_31(idAtm);}
            else if ((atoms[idAtm] == 8)||(atoms[idAtm] == 16)) {SetSP3_22(idAtm);}
          }
        }
      }
    }
  }
  void FixNplanar(int atomA) {
    //function fixing isolated 3,0 nitrogens and nitroso
    int nsp2_neigh = 0;            //count the number of sp2 atoms in the neighbourhood
    double distance;
    for (size_t idneigh = 0; idneigh < neighbours[atomA]; ++idneigh) {
      atomB = connectivity(atomA + 1,idneigh + 1) - 1;
      if ((modVSEPR(atomB + 1,1) == 3.0)&&(modVSEPR(atomB + 1,2) == 0.0)) {++nsp2_neigh;}
      else if ((modVSEPR(atomB + 1,1) == 2.0)&&(modVSEPR(atomB + 1,2) == 1.0)) {++nsp2_neigh;}
      else if ((modVSEPR(atomB + 1,1) == 1.0)&&(modVSEPR(atomB + 1,2) == 0.0)) {++nsp2_neigh;}
    }
    if (nsp2_neigh == 0) {SetSP3_31(atomA);}
    if (neighbours[atomA] == 2) {
      set3Neighbours(atomA,atomB,atomC,atomD);
      if (((atoms[atomB] == 8)||(atoms[atomB] == 16))&&(neighbours[atomB] == 1)) {
        distance = Distance(atomA + 1,atomB + 1,geom);
        if (distance < MinimalBondLength(atoms[atomA],atoms[atomB],1)) {
          SetSP2_21(atomA);
          SetTerminal(atomB);
        }
        else {
          SetPlanar(atomA);
          SetSP3_22(atomB);
        }
      }
    }
  }
  void FixN31Next2Sp2(int atomA, int neighA) {
    //function that checks and fixes sp3 nitrogens bound to sp2 atoms
    int nsp2_neigh = 0;            //count the number of sp2 atoms in the neighbourhood
    int atomB;
    for (size_t idneigh = 0; idneigh < neighA; ++idneigh) {
      atomB = connectivity(atomA + 1,idneigh + 1) - 1;
      if ((atoms[atomB] == 6)&&(modVSEPR(atomB + 1,1) == 3.0)&&(modVSEPR(atomB + 1,2) == 0.0)) {++nsp2_neigh;}
      else if ((atoms[atomB] == 7)&&(modVSEPR(atomB + 1,1) == 2.0)&&(modVSEPR(atomB + 1,2) == 1.0)) {++nsp2_neigh;}
      else if ((atoms[atomB] == 8)&&(modVSEPR(atomB + 1,1) == 1.0)&&(modVSEPR(atomB + 1,2) == 0.0)) {++nsp2_neigh;}
    }
    if (nsp2_neigh != 0) {SetPlanar(atomA);}
  }
  void AminolCounter(int atomC, int atomD, int atomE,int & counter) {
    if ((atoms[atomC] == 8)&&(modVSEPR(atomC + 1,1) == 2.0)&&(modVSEPR(atomC + 1,2) == 2.0)&&(neighbours[atomC] == 1)) {
      if ((atoms[atomD] == 7)&&(modVSEPR(atomD + 1,1) == 3.0)&&(modVSEPR(atomD + 1,2) == 0.0)) {++counter;}
      else if ((atoms[atomE] == 7)&&(modVSEPR(atomE + 1,1) == 3.0)&&(modVSEPR(atomE + 1,2) == 0.0)) {++counter;}
    }
  }
  void ChangeAminol(int atomC, int atomD, int atomE) {
    if ((atoms[atomC] == 8)&&(modVSEPR(atomC + 1,1) == 2.0)&&(modVSEPR(atomC + 1,2) == 2.0)&&(neighbours[atomC] == 1)) {
      if ((atoms[atomD] == 7)&&(modVSEPR(atomD + 1,1) == 3.0)&&(modVSEPR(atomD + 1,2) == 0.0)) {SetTerminal(atomC);}
      else if ((atoms[atomE] == 7)&&(modVSEPR(atomE + 1,1) == 3.0)&&(modVSEPR(atomE + 1,2) == 0.0)) {SetTerminal(atomC);}
    }
  }
  bool IsCarbonyl(int atomA, int nneighbours) {
    //function to check whether a carbon is carbonyl group
    bool iscarbonyl = false;
    if (atoms[atomA] == 6) {
      int nTermO = 0;
      int atomB = -1;
      for (size_t idneigh = 0; idneigh < nneighbours; ++idneigh) {
        atomB = connectivity(atomA + 1,idneigh + 1) - 1;
        if (atomB < 0) {break;}
        if ((atoms[atomB] == 8)&&(modVSEPR(atomB + 1,1) == 1.0)&&(modVSEPR(atomB + 1,2) == 0.0)) {++nTermO;}
        else if ((atoms[atomB] == 8)&&(modVSEPR(atomB + 1,1) == 2.0)&&(modVSEPR(atomB + 1,2) == 0.0)&&(connectivity(atomB + 1,2) == 0)) {++nTermO;}
      }
      if (nTermO > 0) {iscarbonyl = true;}
    }
    return iscarbonyl;
  }
  void MissedCarbonyl(int idAtm, int atomB, int atomC, int atomD) {
    SetPlanar(idAtm);
    if ((atoms[atomB] == 8)&&(neighbours[atomB] == 1)) {SetTerminal(atomB);}
    if ((atoms[atomC] == 8)&&(neighbours[atomC] == 1)) {SetTerminal(atomC);}
    if ((atoms[atomD] == 8)&&(neighbours[atomD] == 1)) {SetTerminal(atomD);}
  }
  void SetCarbonyl(int idAtm, int idBtm, int nrings) {
    //function to set a carbonyl group
    int isringatom;
    int nflatatoms;
    atomC = connectivity(idAtm + 1,1) - 1;
    if ((atoms[atomC] == 6)&&(modVSEPR(atomC + 1,1) == 3.0)&&(modVSEPR(atomC + 1,2) == 0.0)) {
      //carbonyl group
      SetTerminal(idAtm);
      UnpairedElectrons[idAtm] = 0.0;
      for (size_t iring = 0; iring < nrings; ++iring) {
        isringatom = 0;
        nflatatoms = 0;
        for (size_t idCtm = 0; idCtm < therings.rows(); ++idCtm) {
          atomD = therings(idCtm + 1,iring + 1) - 1;
          if (atomD < 0) {
            nflatatoms = idCtm;
            break;
          }
          isringatom += (atomD == idBtm) + (atomD == atomC);
        }
        if (isringatom == 2) {break;}
      }
      if ((isringatom == 2)&&(nflatatoms%2 != 0)) {
        //carbonyl in a ring with odd number of atoms
        SetPlanar(idBtm);
        UnpairedElectrons[idBtm] = 0.0;
      }
    }
  }
  void VSEPRComponents(double threshang = 6.3, int distance_filter = 1, int verbose = 1) {
    //function that applies some modified form of VSEPR to determine the type of geometry around each atom
    //the return matrix has dimensions of the total number of atoms times two; for each atom we give the subscript for X (number of neighbours) and for E (which is obtained from the geometry)
    //possible outcomes for each atom:
    //if second position outputs -1, then terminal atom; else if second position outputs -2, then atom not identified
    //threshang is the error allowed in angles for assigning geometry
    //distance_filter decides whether to use bond length checking
    int nneigh;
    int nneighB;
    int nsp;
    int nsp2;
    int nsp3;
    int undefO;
    int CisArom;
    int DisArom;
    int theringend;
    int torsionatm[7];
    int replacement[4];
    int nflattorsions;
    int lessthansp3;
    int conndoublecount;
    int atomcount;
    double sum_unpaired_electrons = 0.0;
    double minangle[3];
    double auxvar;
    std::vector<int> connectionpath;
    std::vector<int> friendlyneighbours(maxneigh,0);
    std::vector<int> angles(5);
    neighbours.resize(Natoms);
    std::vector<int> ForceHybridisation(Natoms,0);
    std::vector<int> flatcenter(Natoms,0);                //1 - yes; -1 - no; 0 - no idea; 3 - amide
    std::vector<int> isringatom(Natoms,0);
    matrix<int> ConnectedDoubles(1,1);
    matrix<int> islinear(maxneigh,maxneigh);
    matrix<int> istrigonal(maxneigh,maxneigh);
    matrix<int> istetrahedral(maxneigh,maxneigh);
    matrix<int> isorthogonal(maxneigh,maxneigh);
    matrix<int> ispentagonal(maxneigh,maxneigh);
    matrixE theta(maxneigh,maxneigh);
    //get the number of neighbours
    for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
      nneigh = 0;
      for (size_t idneigh = 0; idneigh < maxneigh; ++idneigh) {
        if (connectivity(idAtm + 1,idneigh + 1) == 0) {break;}
        ++nneigh;
      }
      neighbours[idAtm] = nneigh;
    }
    //first round of classifications
    for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
      atomA = atoms[idAtm];
      if (atomA == 1) {
        SetProton(idAtm);
        continue;
      }
      nneigh = neighbours[idAtm];
      if (nneigh == 1) {
        if ((atomA == 9)||(atomA == 17)||(atomA == 35)||(atomA == 53)) {SetTerminal(idAtm);}
        continue;                           //other terminal, do later
      }
      for (size_t idneigh = 0; idneigh < maxneigh; ++idneigh) {
        friendlyneighbours[idneigh] = int(connectivity(idAtm + 1,idneigh + 1)) - 1;
      }
      theta.zero();
      //now get all the angles and classify them
      GetAngles(angles,friendlyneighbours,geom,islinear,istrigonal,istetrahedral,isorthogonal,ispentagonal,theta,idAtm,nneigh,threshang);
      maxangle[0] = theta(2,1);
      maxangle[1] = theta(3,1);
      maxangle[2] = theta(3,2);
      if ((atomA == 5)||(atomA == 6)) {
        if (nneigh == 2) {
          if ((theta(2,1) > 109.5 - 2.0*threshang)&&(theta(2,1) < 109.5 + threshang)) {SetSP3_40(idAtm);}
          else if ((theta(2,1) > 120.0 - threshang)&&(theta(2,1) < 120.0 + 3.0*threshang)) {SetPlanar(idAtm);}
          else if ((theta(2,1) > 180.0 - threshang)&&(theta(2,1) < 180.0 + threshang)) {SetSP_20(idAtm);}
        }
        else if (nneigh == 3) {
          AngleBasedHybridisationNeighbours_3center(maxangle,nsp2,nsp3,lessthansp3,threshang);
          if (nsp3 > nsp2) {SetSP3_40(idAtm);}
          else {SetPlanar(idAtm);}
        }
        else if (nneigh == 4) {SetSP3_40(idAtm);}
      }
      else if ((atomA == 7)||(atomA == 15)) {
        if (nneigh == 2) {
          if ((theta(2,1) > 109.5 - 2.0*threshang)&&(theta(2,1) < 109.5 + threshang)) {SetSP3_31(idAtm);}
          else if ((theta(2,1) > 120.0 - threshang)&&(theta(2,1) < 120.0 + 3.0*threshang)) {SetPlanar(idAtm);}
          else if ((theta(2,1) > 180.0 - threshang)&&(theta(2,1) < 180.0 + threshang)) {SetSP_20(idAtm);}
        }
        else if (nneigh == 3) {
          AngleBasedHybridisationNeighbours_3center(maxangle,nsp2,nsp3,lessthansp3,threshang);
          if (nsp3 > nsp2) {SetSP3_31(idAtm);}
          else {
            if (lessthansp3 > 0) {SetSP3_31(idAtm);}
            else {SetPlanar(idAtm);}
          }
        }
        else if (nneigh == 4) {SetSP3_40(idAtm);}
      }
      else if ((atomA == 8)||(atomA == 16)) {
        if (nneigh == 2) {SetSP3_22(idAtm);}
        else if ((atomA == 16)&&(nneigh == 4)) {SetSP3_40(idAtm);}
        else if ((atomA == 16)&&(nneigh == 3)) {SetSP3_31(idAtm);}
        else {continue;}
      }
    }
    //check ring atoms and fix them
    RingFinder();
    for (size_t iring = 0; iring < therings.cols(); ++iring) {
      theringend = RingSize(iring);
      if (theringend < 11) {
        for (size_t idAtm = 0; idAtm < theringend; ++idAtm) {
          isringatom[therings(idAtm + 1,iring + 1) - 1] = 1;
        }
      }
      for (int idAtm = 0; idAtm < therings.rows(); ++idAtm) {
        if (idAtm > theringend) {break;}
        //get indices for torsions
        torsionatm[0] = idAtm - 3;
        torsionatm[1] = idAtm - 2;
        torsionatm[2] = idAtm - 1;
        if (torsionatm[2] < 0) {
          torsionatm[0] = theringend - 2;
          torsionatm[1] = theringend - 1;
          torsionatm[2] = theringend;
        }
        if (torsionatm[1] < 0) {
          torsionatm[0] = theringend - 1;
          torsionatm[1] = theringend;
        }
        if (torsionatm[0] < 0) {torsionatm[0] = theringend;}
        torsionatm[3] = idAtm;
        torsionatm[4] = idAtm + 1;
        torsionatm[5] = idAtm + 2;
        torsionatm[6] = idAtm + 3;
        if (torsionatm[4] > theringend) {
          torsionatm[4] = 0;
          torsionatm[5] = 1;
          torsionatm[6] = 2;
        }
        if (torsionatm[5] > theringend) {
          torsionatm[5] = 0;
          torsionatm[6] = 1;
        }
        if (torsionatm[6] > theringend) {torsionatm[6] = 0;}
        //get torsions
        for (size_t idx = 0; idx < 4; ++idx) {
          replacement[idx] = therings(torsionatm[idx + 3] + 1,iring + 1);
        }
        if ((theringend == 2)||(theringend == 3)) {
          atomA = therings(torsionatm[3] + 1,iring + 1);
          atomB = therings(torsionatm[4] + 1,iring + 1);
          atomC = therings(torsionatm[5] + 1,iring + 1);
          atomD = therings(torsionatm[6] + 1,iring + 1);
          atomE = therings(torsionatm[2] + 1,iring + 1);
          atomF = therings(torsionatm[1] + 1,iring + 1);
          //3 and 4 member rings are always flat, find the replacement atoms
          for (size_t idneigh = 0; idneigh < maxneigh; ++idneigh) {
            if ((connectivity(atomE,idneigh + 1) != 0)&&(connectivity(atomE,idneigh + 1) != atomA)&&(connectivity(atomE,idneigh + 1) != atomF)) {replacement[0] = connectivity(atomE,idneigh + 1);}
            if ((connectivity(atomA,idneigh + 1) != 0)&&(connectivity(atomA,idneigh + 1) != atomB)&&(connectivity(atomA,idneigh + 1) != atomE)) {replacement[1] = connectivity(atomA,idneigh + 1);}
            if ((connectivity(atomB,idneigh + 1) != 0)&&(connectivity(atomB,idneigh + 1) != atomC)&&(connectivity(atomB,idneigh + 1) != atomA)) {replacement[2] = connectivity(atomB,idneigh + 1);}
            if ((connectivity(atomC,idneigh + 1) != 0)&&(connectivity(atomC,idneigh + 1) != atomD)&&(connectivity(atomC,idneigh + 1) != atomB)) {replacement[3] = connectivity(atomC,idneigh + 1);}
          }
        }
        for (size_t idx = 0; idx < 4; ++idx) {
          torsions[idx] = radian2degree*Torsion(therings(torsionatm[idx] + 1,iring + 1),therings(torsionatm[idx + 1] + 1,iring + 1),therings(torsionatm[idx + 2] + 1,iring + 1),replacement[idx],geom);
        }
        nflattorsions = 0;
        for (size_t idx = 0; idx < 4; ++idx) {
          if (fabs(torsions[idx]) < 2.0*threshang) {++nflattorsions;}
          else if (fabs(torsions[idx] - 180.0) < 2.0*threshang) {++nflattorsions;}
          else if (fabs(torsions[idx] - 360.0) < 2.0*threshang) {++nflattorsions;}
        }
        if ((theringend == 2)||(theringend == 3)) {
          //only check the one with non-ring atom
          if ((nflattorsions == 1)&&(fabs(torsions[1]) < 2.0*threshang)&&(replacement[1] != therings(torsionatm[4] + 1,iring + 1))) {nflattorsions = 2;}
          else if ((nflattorsions == 1)&&(fabs(torsions[1] - 180.0) < 2.0*threshang)&&(replacement[1] != therings(torsionatm[4] + 1,iring + 1))) {nflattorsions = 2;}
          else if ((nflattorsions == 1)&&(fabs(torsions[1] - 360.0) < 2.0*threshang)&&(replacement[1] != therings(torsionatm[4] + 1,iring + 1))) {nflattorsions = 2;}
          else {nflattorsions = 0;}
        }
        if (nflattorsions > 1) {
          if (neighbours[therings(torsionatm[3] + 1,iring + 1) - 1] == 4) {SetSP3_40(therings(torsionatm[3] + 1,iring + 1) - 1);}
          //geometry around atom is flat, so make it sp2 (ring)
          else {
            if ((atoms[therings(torsionatm[3] + 1,iring + 1) - 1] == 5)||(atoms[therings(torsionatm[3] + 1,iring + 1) - 1] == 6)||(atoms[therings(torsionatm[3] + 1,iring + 1) - 1] == 14)) {SetPlanar(therings(torsionatm[3] + 1,iring + 1) - 1);}
            else if ((atoms[therings(torsionatm[3] + 1,iring + 1) - 1] == 7)||(atoms[therings(torsionatm[3] + 1,iring + 1) - 1] == 15)) {
              if (neighbours[therings(torsionatm[3] + 1,iring + 1) - 1] == 2) {SetSP2_21(therings(torsionatm[3] + 1,iring + 1) - 1);}
              else if (neighbours[therings(torsionatm[3] + 1,iring + 1) - 1] == 3) {SetPlanar(therings(torsionatm[3] + 1,iring + 1) - 1);}
              else if (neighbours[therings(torsionatm[3] + 1,iring + 1) - 1] == 4) {SetSP3_31(therings(torsionatm[3] + 1,iring + 1) - 1);}
            }
            else if ((atoms[therings(torsionatm[3] + 1,iring + 1) - 1] == 8)||(atoms[therings(torsionatm[3] + 1,iring + 1) - 1] == 16)) {SetSP3_22(therings(torsionatm[3] + 1,iring + 1) - 1);}
          }
        }
        else {
          //geometry around atom is not flat, so make it sp3 (ring)
          //but only do this if the atoms was not previously classified
          if (((modVSEPR(therings(torsionatm[3] + 1,iring + 1),1) == 0.0)&&(modVSEPR(therings(torsionatm[3] + 1,iring + 1),2) == 0.0))||(theringend == 2)) {
            if ((atoms[therings(torsionatm[3] + 1,iring + 1) - 1] == 5)||(atoms[therings(torsionatm[3] + 1,iring + 1) - 1] == 6)||(atoms[therings(torsionatm[3] + 1,iring + 1) - 1] == 14)) {SetSP3_40(therings(torsionatm[3] + 1,iring + 1) - 1);}
            else if ((atoms[therings(torsionatm[3] + 1,iring + 1) - 1] == 7)||(atoms[therings(torsionatm[3] + 1,iring + 1) - 1] == 15)) {
              if (neighbours[therings(torsionatm[3] + 1,iring + 1) - 1] == 4) {SetSP3_40(therings(torsionatm[3] + 1,iring + 1) - 1);}
              else {SetSP3_31(therings(torsionatm[3] + 1,iring + 1) - 1);}
            }
            else if ((atoms[therings(torsionatm[3] + 1,iring + 1) - 1] == 8)||(atoms[therings(torsionatm[3] + 1,iring + 1) - 1] == 16)) {SetSP3_22(therings(torsionatm[3] + 1,iring + 1) - 1);}
          }
        }
      }
      lessthansp3 = 0;
      nO = 0;            //number of lactams
      for (int idAtm = 0; idAtm < therings.rows(); ++idAtm) {
        if (idAtm > theringend) {break;}
        atomA = therings(idAtm + 1,iring + 1) - 1;
        if (((modVSEPR(atomA + 1,1) == 3.0)&&(modVSEPR(atomA + 1,2) == 0.0))||((modVSEPR(atomA + 1,1) == 2.0)&&(modVSEPR(atomA + 1,2) == 1.0))) {++lessthansp3;}
        if ((atoms[atomA] == 6)&&(neighbours[atomA] == 3)) {
          set3Neighbours(atomA,atomB,atomC,atomD);
          if (atoms[atomB] == 7) {nO += IsAtomInRing(therings,theringend,iring,atomB)*((atoms[atomC] == 8)*(neighbours[atomC] == 1) + (atoms[atomD] == 8)*(neighbours[atomD] == 1));}
          else if (atoms[atomC] == 7) {nO += IsAtomInRing(therings,theringend,iring,atomC)*((atoms[atomB] == 8)*(neighbours[atomB] == 1) + (atoms[atomD] == 8)*(neighbours[atomD] == 1));}
          else if (atoms[atomD] == 7) {nO += IsAtomInRing(therings,theringend,iring,atomD)*((atoms[atomB] == 8)*(neighbours[atomB] == 1) + (atoms[atomC] == 8)*(neighbours[atomC] == 1));}
        }
      }
      if ((lessthansp3 == theringend + 1)&&(nO == 0)&&(theringend + 1 == 6)) {
        //aromatic ring, mark it
        for (int idAtm = 0; idAtm < therings.rows(); ++idAtm) {
          if (idAtm > theringend) {break;}
          atomA = therings(idAtm + 1,iring + 1) - 1;
          isringatom[atomA] = 2;
        }
      }
    }
    //terminal atoms
    for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
      nneigh = neighbours[idAtm];
      if (nneigh == 1) {
        atomB = connectivity(idAtm + 1,1) - 1;
        nneighB = neighbours[atomB];
        nC = 0;
        nN = 0;
        nsp2 = 0;
        for (size_t idneigh = 0; idneigh < nneighB; ++ idneigh) {
          atomC = connectivity(atomB + 1,idneigh + 1) - 1;
          if (atomC == idAtm) {continue;}
          if (((modVSEPR(atomC + 1,1) == 3.0)&&(modVSEPR(atomC + 1,2) == 0.0))||((modVSEPR(atomC + 1,1) == 2.0)&&(modVSEPR(atomC + 1,2) == 1.0))) {++nsp2;}
          if (atoms[atomC] == 6) {++nC;}
          else if (atoms[atomC] == 7) {++nN;}
        }
        if (((modVSEPR(atomB + 1,1) == 4.0)&&(modVSEPR(atomB + 1,2) == 0.0))||((modVSEPR(atomB + 1,1) == 3.0)&&(modVSEPR(atomB + 1,2) == 1.0))||((modVSEPR(atomB + 1,1) == 2.0)&&(modVSEPR(atomB + 1,2) == 2.0))) {
          if (atoms[idAtm] == 5) {SetPlanar(idAtm);}
          else if ((atoms[idAtm] == 6)||(atoms[idAtm] == 14)) {SetSP3_40(idAtm);}
          else if ((atoms[idAtm] == 7)||(atoms[idAtm] == 15)) {SetSP3_31(idAtm);}
          else if ((atoms[idAtm] == 8)||(atoms[idAtm] == 16)) {
            if ((atoms[atomB] == 15)||(atoms[atomB] == 16)) {SetTerminal(idAtm);}
            else {SetSP3_22(idAtm);}
          }
        }
        else if (((modVSEPR(atomB + 1,1) == 3.0)&&(modVSEPR(atomB + 1,2) == 0.0))||((modVSEPR(atomB + 1,1) == 2.0)&&(modVSEPR(atomB + 1,2) == 1.0))) {
          if (atoms[idAtm] == 5) {SetSP2_21(idAtm);}
          else if ((atoms[idAtm] == 6)||(atoms[idAtm] == 14)) {
            if ((nC == 2)&&(nsp2 == 2)) {SetSP3_40(idAtm);}
            else if ((nC == 1)&&(nN == 1)&&(nsp2 == 2)) {SetSP3_40(idAtm);}
            else if ((nN == 2)&&(nsp2 == 2)) {SetSP3_40(idAtm);}
            else if (atoms[atomB] == 7) {SetSP3_40(idAtm);}
            else {SetPlanar(idAtm);}
          }
          else if ((atoms[idAtm] == 7)||(atoms[idAtm] == 15)) {SetPlanar(idAtm);}
          else if ((atoms[idAtm] == 8)||(atoms[idAtm] == 16)) {SetTerminal(idAtm);}
        }
        else if ((modVSEPR(atomB + 1,1) == 2.0)&&(modVSEPR(atomB + 1,2) == 0.0)) {
          if (atoms[idAtm] == 5) {SetSP_20(idAtm);}
          else if ((atoms[idAtm] == 6)||(atoms[idAtm] == 14)) {SetSP_20(idAtm);}
          else if ((atoms[idAtm] == 7)||(atoms[idAtm] == 15)) {SetSP_20(idAtm);}
          else if ((atoms[idAtm] == 8)||(atoms[idAtm] == 16)) {SetUNDEFINED(idAtm);}
        }
      }
    }
    //check the actual angles around atoms with 3 neighbours
    for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
      if (neighbours[idAtm] == 3) {
        set3Neighbours(idAtm,atomB,atomC,atomD);
        nsp2 = CheckOOP(idAtm + 1,atomB + 1,atomC + 1,atomD + 1,18.5,17.5);
        if (nsp2 == 4) {
          flatcenter[idAtm] = 1;
          nsp2 = (atoms[atomB] == 8) + (atoms[atomC] == 8) + (atoms[atomD] == 8);
          nsp3 = ((modVSEPR(atomB + 1,1) == 4.0)*(modVSEPR(atomB + 1,2) == 0.0) + (modVSEPR(atomB + 1,1) == 3.0)*(modVSEPR(atomB + 1,2) == 1.0)) + ((modVSEPR(atomC + 1,1) == 4.0)*(modVSEPR(atomC + 1,2) == 0.0) + (modVSEPR(atomC + 1,1) == 3.0)*(modVSEPR(atomC + 1,2) == 1.0)) + ((modVSEPR(atomD + 1,1) == 4.0)*(modVSEPR(atomD + 1,2) == 0.0) + (modVSEPR(atomD + 1,1) == 3.0)*(modVSEPR(atomD + 1,2) == 1.0));
          if (atoms[idAtm] == 6) {
            if (nsp2 > 1) {MultipleTerminals(atomB,atomC,atomD);}                    //lactone
            else if ((nsp2 == 1)&&(nsp3 == 2)) {
              SetPlanar(idAtm);
              MultipleTerminals(atomB,atomC,atomD);
            }
          }
          else if (atoms[idAtm] == 7) {
            if (nsp2 > 1) {MultipleTerminals(atomB,atomC,atomD);}                    //nitro or nitrate
          }
          nsp2 = 4;
        }
        else if (nsp2 < 4) {flatcenter[idAtm] = -1;}
        if ((nsp2 == 4)&&(modVSEPR(idAtm + 1,1) == 3.0)&&(modVSEPR(idAtm + 1,2) == 1.0)) {SetPlanar(idAtm);}
        else if ((nsp2 == 4)&&(modVSEPR(idAtm + 1,1) == 4.0)&&(modVSEPR(idAtm + 1,2) == 0.0)) {SetPlanar(idAtm);}
        else if ((nsp2 < 4)&&(modVSEPR(idAtm + 1,1) == 2.0)&&(modVSEPR(idAtm + 1,2) == 1.0)) {
          if ((atoms[idAtm] == 6)||(atoms[idAtm] == 14)) {SetSP3_40(idAtm);}
          else if ((atoms[idAtm] == 7)||(atoms[idAtm] == 15)) {SetSP3_31(idAtm);}
          else if ((atoms[idAtm] == 8)||(atoms[idAtm] == 16)) {SetSP3_22(idAtm);}
        }
        else if ((nsp2 < 4)&&(modVSEPR(idAtm + 1,1) == 3.0)&&(modVSEPR(idAtm + 1,2) == 0.0)) {
          if ((atoms[idAtm] == 6)||(atoms[idAtm] == 14)) {SetSP3_40(idAtm);}
          else if ((atoms[idAtm] == 7)||(atoms[idAtm] == 15)) {SetSP3_31(idAtm);}
          else if ((atoms[idAtm] == 8)||(atoms[idAtm] == 16)) {SetSP3_22(idAtm);}
        }
        if ((nsp2 == 4)&&(modVSEPR(idAtm + 1,1) == 3.0)&&(modVSEPR(idAtm + 1,2) == 0.0)) {ForceHybridisation[idAtm] = 1;}
        else if ((nsp2 < 4)&&(modVSEPR(idAtm + 1,1) == 4.0)&&(modVSEPR(idAtm + 1,2) == 0.0)) {ForceHybridisation[idAtm] = 1;}
        else if ((nsp2 < 4)&&(modVSEPR(idAtm + 1,1) == 3.0)&&(modVSEPR(idAtm + 1,2) == 1.0)) {ForceHybridisation[idAtm] = 1;}
      }
      else if (neighbours[idAtm] == 4) {flatcenter[idAtm] = -1;}
    }
    //check for 6 member rings with one sp3 carbon, 4 sp2 carbons, and 1 2,1 nitrogen
    for (size_t iring = 0; iring < therings.cols(); ++iring) {
      theringend = RingSize(iring) + 1;
      if (theringend != 6) {continue;}
      nsp2 = 0;
      nsp3 = 0;
      nN = 0;
      for (int idAtm = 0; idAtm < therings.rows(); ++idAtm) {
        if (idAtm == theringend) {break;}
        if ((modVSEPR(therings(idAtm + 1,iring + 1),1) == 3.0)&&(modVSEPR(therings(idAtm + 1,iring + 1),2) == 0.0)) {++nsp2;}
        else if ((modVSEPR(therings(idAtm + 1,iring + 1),1) == 2.0)&&(modVSEPR(therings(idAtm + 1,iring + 1),2) == 1.0)) {++nsp2;}
        else if ((modVSEPR(therings(idAtm + 1,iring + 1),1) == 4.0)&&(modVSEPR(therings(idAtm + 1,iring + 1),2) == 0.0)) {++nsp3;}
        else if ((modVSEPR(therings(idAtm + 1,iring + 1),1) == 3.0)&&(modVSEPR(therings(idAtm + 1,iring + 1),2) == 1.0)) {++nsp3;}
        if (atoms[therings(idAtm + 1,iring + 1) - 1] == 7) {++nN;}
      }
      if ((nsp2 == 5)&&(nsp3 == 1)&&(nN == 1)) {
        //find and neutralise the nitrogen
        for (int idAtm = 0; idAtm < therings.rows(); ++idAtm) {
          if (atoms[therings(idAtm + 1,iring + 1) - 1] == 7) {
            SetPlanar(therings(idAtm + 1,iring + 1) - 1);
            break;
          }
        }
      }
    }
    //terminal carbonyls and unassigned atoms
    for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
      if ((modVSEPR(idAtm + 1,1) == 0.0)&&(modVSEPR(idAtm + 1,2) == 0.0)) {
        if ((atoms[idAtm] == 7)&&(neighbours[idAtm] == 2)) {
          nsp2 = 0;
          for (size_t idneigh = 0; idneigh < neighbours[idAtm]; ++ idneigh) {
            atomB = connectivity(idAtm + 1,idneigh + 1) - 1;
            if (((modVSEPR(atomB + 1,1) == 3.0)&&(modVSEPR(atomB + 1,2) == 0.0))||((modVSEPR(atomB + 1,1) == 2.0)&&(modVSEPR(atomB + 1,2) == 1.0))) {++nsp2;}
          }
          if (nsp2 > 0) {SetPlanar(idAtm);}
          else {SetSP3_31(idAtm);}
        }
      }
      else {
        //check for terminal carbonyls/carboxyls and related groups
        undefO = 0;
        for (size_t ineigh = 0; ineigh < maxneigh; ++ineigh) {
          atomB = connectivity(idAtm + 1,ineigh + 1) - 1;
          if ((atomB >= 0)&&(atoms[atomB] == 8)) {
            if ((modVSEPR(atomB + 1,1) == 0)&&(modVSEPR(atomB + 1,2) == 0)) {++undefO;}
          }
        }
        if (undefO > 0) {
          for (size_t ineigh = 0; ineigh < maxneigh; ++ineigh) {
            atomB = connectivity(idAtm + 1,ineigh + 1) - 1;
            if ((atomB >= 0)&&(atoms[atomB] != 8)) {
              nneigh = neighbours[atomB];
              if (nneigh == 1) {
                if (atoms[atomB] == 6) {SetSP3_40(atomB);}
                else if (atoms[atomB] == 7) {SetPlanar(atomB);}
              }
            }
          }
        }
      }
    }
    //fix amides and other weird things
    for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
      if ((atoms[idAtm] == 6)&&(neighbours[idAtm] == 3)) {
        nC = 0;
        nN = 0;
        nO = 0;
        nsp2 = 0;        //only sp2 oxygens
        //count some properties
        for (size_t idneigh = 0; idneigh < neighbours[idAtm]; ++idneigh) {
          atomB = connectivity(idAtm + 1,idneigh + 1) - 1;
          if (atoms[atomB] == 6) {++nC;}
          else if (atoms[atomB] == 7) {++nN;}
          else if (atoms[atomB] == 8) {
            ++nO;
            if ((modVSEPR(atomB + 1,1) == 1.0)&&(modVSEPR(atomB + 1,2) == 0.0)) {++nsp2;}
          }
        }
        //amides
        if ((nN > 0)&&(nO > 0)&&(nsp2 > 0)) {
          //search for the nitrogen
          for (size_t idneigh = 0; idneigh < neighbours[idAtm]; ++idneigh) {
            atomB = connectivity(idAtm + 1,idneigh + 1) - 1;
            //this is supposed to be amide, so fix
            if (atoms[atomB] == 7) {
              if ((modVSEPR(atomB + 1,1) == 2.0)&&(modVSEPR(atomB + 1,2) == 1.0)) {SetPlanar(atomB);}
              flatcenter[atomB] = 3;
            }
            //this is supposed to be terminal from amide, so fix
            else if (atoms[atomB] == 6) {
              if ((modVSEPR(atomB + 1,1) == 3.0)&&(modVSEPR(atomB + 1,2) == 0.0)) {
                if (neighbours[atomB] == 1) {SetSP3_40(atomB);}
                else {
                  nsp2 = 0;        //number of sp2 atoms
                  for (size_t idneigh2 = 0; idneigh2 < neighbours[atomB]; ++idneigh2) {
                    atomC = connectivity(atomB + 1,idneigh2 + 1) - 1;
                    if (atomC < 0) {break;}
                    else if (atomC == idAtm) {continue;}
                    else if ((modVSEPR(atomC + 1,1) == 3.0)&&(modVSEPR(atomC + 1,2) == 0.0)&&(atoms[atomC] == 6)) {++nsp2;}
                    else if ((modVSEPR(atomC + 1,1) == 2.0)&&(modVSEPR(atomC + 1,2) == 1.0)&&(atoms[atomC] == 7)) {++nsp2;}
                    else if ((atoms[atomC] == 8)&&(neighbours[atomC] == 1)) {++nsp2;}                   //just in case
                  }
                  if (nsp2 == 0) {SetSP3_40(atomB);}
                }
              }
            }
          }
        }
        else if ((nN > 0)&&(nO > 0)&&(nsp2 == 0)) {
          if (isringatom[idAtm] == 0) {
            for (size_t idneigh = 0; idneigh < neighbours[idAtm]; ++idneigh) {
              atomB = connectivity(idAtm + 1,idneigh + 1) - 1;
              if ((atoms[atomB] == 8)&&(connectivity(atomB + 1,2) == 0)) {nsp2 = 1;}
            }
            if (nsp2 == 1) {
              //if we are here, then lets check angles
              set3Neighbours(idAtm,atomB,atomC,atomD);
              nsp2 = CheckOOP(idAtm + 1,atomB + 1,atomC + 1,atomD + 1,18.5,16);
              if (nsp2 == 4) {
                //missed amide
                for (size_t idneigh = 0; idneigh < neighbours[idAtm]; ++idneigh) {
                  atomB = connectivity(idAtm + 1,idneigh + 1) - 1;
                  if ((atoms[atomB] == 6)&&(neighbours[atomB] == 3)) {
                    set3Neighbours(atomB,atomC,atomD,atomE);
                    nsp2 = CheckOOP(atomB + 1,atomC + 1,atomD + 1,atomE + 1,15.0,16.0);
                    if (nsp2 != 4) {SetSP3_40(atomB);}
                  }
                  else if (atoms[atomB] == 7) {SetPlanar(atomB);}
                  else if ((atoms[atomB] == 8)&&(connectivity(atomB + 1,2) == 0)) {SetTerminal(atomB);}
                }
              }
            }
          }
        }
        //carbonyls
        else if ((nN == 0)&&(nO > 0)&&(nC > 0)&&(nsp2 > 0)) {
          //find the carbons and check whether they are terminal
          for (size_t idneigh = 0; idneigh < neighbours[idAtm]; ++idneigh) {
            atomB = connectivity(idAtm + 1,idneigh + 1) - 1;
            if ((atoms[atomB] == 6)&&(neighbours[atomB] == 1)&&(modVSEPR(atomB + 1,1) != 4.0)) {SetSP3_40(atomB);}
          }
        }
      }
      //isolated sp2 carbons shall not exist (radical)
      if ((atoms[idAtm] == 6)&&(modVSEPR(idAtm + 1,1) == 3.0)&&(modVSEPR(idAtm + 1,2) == 0.0)) {
        nsp2 = 0;            //count the number of sp2 atoms in the neighbourhood
        for (size_t idneigh = 0; idneigh < neighbours[idAtm]; ++idneigh) {
          atomB = connectivity(idAtm + 1,idneigh + 1) - 1;
          if ((modVSEPR(atomB + 1,1) == 3.0)&&(modVSEPR(atomB + 1,2) == 0.0)) {++nsp2;}
          else if ((modVSEPR(atomB + 1,1) == 2.0)&&(modVSEPR(atomB + 1,2) == 1.0)) {++nsp2;}
          else if ((modVSEPR(atomB + 1,1) == 1.0)&&(modVSEPR(atomB + 1,2) == 0.0)) {++nsp2;}
        }
        if (nsp2 == 0) {         //nsp2 controls whether any modification took place
          //check whether in position 3 there is sp2 carbon
          if (ForceHybridisation[idAtm] == 1) {
            for (size_t idneigh = 0; idneigh < neighbours[idAtm]; ++idneigh) {
              atomB = connectivity(idAtm + 1,idneigh + 1) - 1;
              if (neighbours[atomB] == 4) {continue;}
              if ((atoms[atomB] == 6)&&(modVSEPR(atomB + 1,1) == 4.0)&&(modVSEPR(atomB + 1,2) == 0.0)) {
                for (size_t idneigh2 = 0; idneigh2 < neighbours[atomB]; ++idneigh2) {
                  atomC = connectivity(atomB + 1,idneigh2 + 1) - 1;
                  if (atomC == idAtm) {continue;}
                  if ((modVSEPR(atomC + 1,1) == 3.0)&&(modVSEPR(atomC + 1,2) == 0.0)&&(atoms[atomC] == 6)) {      //change atomB
                    nsp2 = 1;
                    SetPlanar(atomB);
                    break;
                  }
                }
              }
              if (nsp2 == 1) {break;}
            }
          }
          if (nsp2 == 0) {SetSP3_40(idAtm);}
        }
      }
      //similar for sp carbons
      if ((atoms[idAtm] == 6)&&(modVSEPR(idAtm + 1,1) == 2.0)&&(modVSEPR(idAtm + 1,2) == 0.0)) {
        set3Neighbours(idAtm,atomB,atomC,atomD);
        if (atomC < 0) {continue;}
        nsp = (modVSEPR(atomB + 1,1) == 2.0)*(modVSEPR(atomB + 1,2) == 0.0) + (modVSEPR(atomC + 1,1) == 2.0)*(modVSEPR(atomC + 1,2) == 0.0);
        if (nsp == 0) {
          if (atoms[idAtm] == 6) {SetSP3_40(idAtm);}
          else if (atoms[idAtm] == 7) {SetSP3_31(idAtm);}
        }
        else if (nsp == 2) {
          maxangle[0] = 0.0;             //rAB2
          maxangle[1] = 0.0;             //rAC2
          for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
            maxangle[2] = geom(idAtm + 1,idcoord + 1) - geom(atomB + 1,idcoord + 1);
            maxangle[0] += maxangle[2]*maxangle[2];
            maxangle[2] = geom(idAtm + 1,idcoord + 1) - geom(atomC + 1,idcoord + 1);
            maxangle[1] += maxangle[2]*maxangle[2];
          }
          if (fabs(maxangle[0] - maxangle[1]) > 0.5) {
            //no allene, must be changed
            if ((neighbours[atomB] == 1)||(neighbours[atomB] == 3)) {
              if (atoms[atomB] == 6) {SetSP3_40(atomB);}
              else if (atoms[atomB] == 7) {SetSP3_31(atomB);}
            }
            else if ((neighbours[atomC] == 1)||(neighbours[atomC] == 3)) {
              if (atoms[atomC] == 6) {SetSP3_40(atomC);}
              else if (atoms[atomC] == 7) {SetSP3_31(atomC);}
            }
          }
        }
      }
      //other special cases
      if ((atoms[idAtm] == 7)&&(modVSEPR(idAtm + 1,1) == 3.0)&&(modVSEPR(idAtm + 1,2) == 0.0)&&(neighbours[idAtm] == 2)) {
        set3Neighbours(idAtm,atomB,atomC,atomD);
        if ((atoms[atomB] == 7)||(atoms[atomC] == 7)) {
          if (atoms[atomC] == 7) {
            //swap the atoms
            set3Neighbours(idAtm,atomC,atomB,atomD);
          }
          //diazo?
          if ((modVSEPR(atomB + 1,1) == 3.0)&&(modVSEPR(atomB + 1,2) == 0.0)&&(neighbours[atomB] == 2)) {
            //check the neighbours of idAtm; is atomC aromatic?
            CisArom = 0;
            if ((atoms[atomC] == 6)&&(neighbours[atomC] == 3)&&(modVSEPR(atomC + 1,1) == 3.0)&&(modVSEPR(atomC + 1,2) == 0.0)) {
              for (size_t idneigh = 0; idneigh < neighbours[atomC]; ++idneigh) {
                atomE = connectivity(atomC + 1,idneigh + 1) - 1;
                if (atomE == idAtm) {continue;}
                if ((modVSEPR(atomE + 1,1) == 3.0)&&(modVSEPR(atomE + 1,2) == 0.0)) {++CisArom;}
              }
            }
            DisArom = 0;
            //look for the neighbours of atomB
            for (size_t idneigh = 0; idneigh < neighbours[atomB]; ++idneigh) {
              atomD = connectivity(atomB + 1,idneigh + 1) - 1;
              if (atomD != idAtm) {break;}
            }
            if ((atoms[atomD] == 6)&&(neighbours[atomD] == 3)&&(modVSEPR(atomD + 1,1) == 3.0)&&(modVSEPR(atomD + 1,2) == 0.0)) {
              for (size_t idneigh = 0; idneigh < neighbours[atomD]; ++idneigh) {
                atomE = connectivity(atomD + 1,idneigh + 1) - 1;
                if (atomE == atomB) {continue;}
                if ((modVSEPR(atomE + 1,1) == 3.0)&&(modVSEPR(atomE + 1,2) == 0.0)) {++DisArom;}
              }
            }
            if ((CisArom == 2)&&(DisArom == 2)) {
              SetSP2_21(idAtm);
              SetSP2_21(atomB);
            }
          }
        }
      }
      //isolated planar nitrogens cannot occur
      if ((atoms[idAtm] == 7)&&(modVSEPR(idAtm + 1,1) == 3.0)&&(modVSEPR(idAtm + 1,2) == 0.0)) {FixNplanar(idAtm);}
      //tetrahedral nitrogens next to sp2 must be checked
      if ((atoms[idAtm] == 7)&&(modVSEPR(idAtm + 1,1) == 3.0)&&(modVSEPR(idAtm + 1,2) == 1.0)) {FixN31Next2Sp2(idAtm,neighbours[idAtm]);}
    }
    if (distance_filter > 0) {
      for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
        if (isringatom[idAtm] == 2) {continue;}
        if ((atoms[idAtm] == 7)&&(modVSEPR(idAtm + 1,1) == 3.0)&&(modVSEPR(idAtm + 1,2) == 0.0)) {FixNplanar(idAtm);}
        if (neighbours[idAtm] == 3) {
          set3Neighbours(idAtm,atomB,atomC,atomD);
          nsp2 = CheckOOP(idAtm + 1,atomB + 1,atomC + 1,atomD + 1,21.0,16);
          lessthansp3 = (modVSEPR(idAtm + 1,1) == 3.0)*(modVSEPR(idAtm + 1,2) == 0.0) + (modVSEPR(idAtm + 1,1) == 2.0)*(modVSEPR(idAtm + 1,2) == 1.0);
          //only check if the atom is not sp2
          if ((nsp2 == 4)&&(lessthansp3 == 0)) {
            //number of nitrogens
            nN = (atoms[atomB] == 7) + (atoms[atomC] == 7) + (atoms[atomD] == 7);
            //number of terminal oxygens
            nO = (atoms[atomB] == 8)*(connectivity(atomB + 1,2) == 0) + (atoms[atomC] == 8)*(connectivity(atomC + 1,2) == 0) + (atoms[atomD] == 8)*(connectivity(atomD + 1,2) == 0);
            //check neighbours, if there is something suitable, modify too.
            if (nN == 3) {
              SetPlanar(idAtm);
              SetPlanar(atomB);
              SetPlanar(atomC);
              SetPlanar(atomD);
            }
            else if (nO >= 1) {SetPlanar(idAtm);}
            else {
              nsp3 = 0;
              for (size_t idneigh = 0; idneigh < 3; ++idneigh) {
                atomB = connectivity(idAtm + 1,idneigh + 1) - 1;
                nsp3 += (modVSEPR(atomB + 1,1) == 4.0)*(modVSEPR(atomB + 1,2) == 0.0) + (modVSEPR(atomB + 1,1) == 3.0)*(modVSEPR(atomB + 1,2) == 1.0);
              }
              if (nsp3 == 3) {
                nsp3 = 0;
                for (size_t idneigh = 0; idneigh < 3; ++idneigh) {
                  atomB = connectivity(idAtm + 1,idneigh + 1) - 1;
                  set3Neighbours(atomB,atomC,atomD,atomE);
                  if (neighbours[atomB] == 3) {
                    nsp2 = CheckOOP(atomB + 1,atomC + 1,atomD + 1,atomE + 1,21.0,16);
                    if (nsp2 == 4) {
                      SetPlanar(idAtm);
                      SetPlanar(atomB);
                      break;
                    }
                  }
                  else if (neighbours[atomB] == 2) {
                    auxvar = radian2degree*Angle(atomC + 1,atomB + 1,atomD + 1,geom);
                    if (auxvar > 109.5 + threshang) {
                      SetPlanar(idAtm);
                      SetPlanar(atomB);
                      break;
                    }
                  }
                  else if (neighbours[atomB] == 1) {++nsp3;}
                }
              }
            }
          }
        }
        nsp2 = 0;
        nsp3 = 0;
        nN = 0;                 //number of sp2 nitrogens
        nO = 0;                 //number of sp2 oxygens
        undefO = 0;             //number of terminal oxygens
        for (size_t idneigh = 0; idneigh < neighbours[idAtm]; ++idneigh) {
          atomB = connectivity(idAtm + 1,idneigh + 1) - 1;
          if (((modVSEPR(atomB + 1,1) == 3.0)&&(modVSEPR(atomB + 1,2) == 0.0))||((modVSEPR(atomB + 1,1) == 2.0)&&(modVSEPR(atomB + 1,2) == 1.0))) {
            ++nsp2;
            if (atoms[atomB] == 7) {++nN;}
          }
          else if ((modVSEPR(atomB + 1,1) == 1.0)&&(modVSEPR(atomB + 1,2) == 0.0)) {
            ++nsp2;
            if (atoms[atomB] == 8) {++nO;}
          }
          else if (((modVSEPR(atomB + 1,1) == 4.0)&&(modVSEPR(atomB + 1,2) == 0.0))||((modVSEPR(atomB + 1,1) == 3.0)&&(modVSEPR(atomB + 1,2) == 1.0))||((modVSEPR(atomB + 1,1) == 2.0)&&(modVSEPR(atomB + 1,2) == 2.0))) {++nsp3;}
          if ((atoms[atomB] == 8)&&(connectivity(atomB + 1,2) == 0)) {++undefO;}
        }
        if ((nsp2 == 1)&&(((atoms[idAtm] == 6)&&(modVSEPR(idAtm + 1,1) == 3.0)&&(modVSEPR(idAtm + 1,2) == 0.0))||((atoms[idAtm] == 7)&&(modVSEPR(idAtm + 1,1) == 2.0)&&(modVSEPR(idAtm + 1,2) == 1.0)))) {
          //sp2 atom with a single sp2 neighbour
          if ((undefO == 1)&&(nO == 0)&&(nN >= 1)) {
            //missed amide
            for (size_t idneigh = 0; idneigh < neighbours[idAtm]; ++idneigh) {
              atomB = connectivity(idAtm + 1,idneigh + 1) - 1;
              if ((atoms[atomB] == 8)&&(connectivity(atomB + 1,2) == 0)) {SetTerminal(atomB);}
              if ((atoms[atomB] == 7)&&(modVSEPR(atomB + 1,1) == 2.0)&&(modVSEPR(atomB + 1,2) == 1.0)) {SetSP2_21(atomB);}
            }
          }
          else {
            for (size_t idneigh = 0; idneigh < neighbours[idAtm]; ++idneigh) {
              atomB = connectivity(idAtm + 1,idneigh + 1) - 1;
              if (((modVSEPR(atomB + 1,1) == 3.0)&&(modVSEPR(atomB + 1,2) == 0.0))||((modVSEPR(atomB + 1,1) == 2.0)&&(modVSEPR(atomB + 1,2) == 1.0))||(((atoms[atomB] == 8)||(atoms[atomB] == 16))&&(modVSEPR(atomB + 1,1) == 1.0)&&(modVSEPR(atomB + 1,2) == 0.0))) {
                auxvar = Distance(idAtm + 1,atomB + 1,geom) + 0.5*(atoms[idAtm] == 6)*(atoms[atomB] == 7)*(neighbours[atomB] == 1);
                if (auxvar >= MinimalBondLength(atoms[idAtm],atoms[atomB],1)) {
                  //make it sp3
                  if (atoms[idAtm] == 6) {SetSP3_40(idAtm);}
                  else if ((atoms[idAtm] == 7)||(atoms[idAtm] == 15)) {SetSP3_31(idAtm);}
                  else if ((atoms[idAtm] == 8)||(atoms[idAtm] == 16)) {SetSP3_22(idAtm);}
                }
                break;
              }
            }
          }
        }
        else if ((nsp2 == 2)&&(((modVSEPR(idAtm + 1,1) == 3.0)&&(modVSEPR(idAtm + 1,2) == 0.0))||((modVSEPR(idAtm + 1,1) == 2.0)&&(modVSEPR(idAtm + 1,2) == 1.0)))) {
          //sp2 atom with a two sp2 neighbours
          nsp2 = 0;
          nsp3 = 0;
          for (size_t idneigh = 0; idneigh < neighbours[idAtm]; ++idneigh) {
            atomB = connectivity(idAtm + 1,idneigh + 1) - 1;
            if (((modVSEPR(atomB + 1,1) == 3.0)&&(modVSEPR(atomB + 1,2) == 0.0))||((modVSEPR(atomB + 1,1) == 2.0)&&(modVSEPR(atomB + 1,2) == 1.0))||((modVSEPR(atomB + 1,1) == 1.0)&&(modVSEPR(atomB + 1,2) == 0.0))) {
              auxvar = Distance(idAtm + 1,atomB + 1,geom);
              if (auxvar >= MinimalBondLength(atoms[idAtm],atoms[atomB],1)) {++nsp3;}
            }
          }
          if (nsp3 == 2) {
            //make it sp3?
            if (neighbours[idAtm] == 3) {
              if (neighbours[idAtm] == 3) {
                set3Neighbours(idAtm,atomB,atomC,atomD);
                nsp2 = CheckOOP(idAtm + 1,atomB + 1,atomC + 1,atomD + 1,10.0,10.0);
              }
              else {nsp2 = 2;}
            }
            if (nsp2 != 4) {
              if (atoms[idAtm] == 6) {SetSP3_40(idAtm);}
              else if (atoms[idAtm] == 7) {SetSP3_31(idAtm);}
              else if (atoms[idAtm] == 8) {SetSP3_22(idAtm);}
              else if (atoms[idAtm] == 15) {SetSP3_31(idAtm);}
              else if (atoms[idAtm] == 16) {SetSP3_22(idAtm);}
            }
          }
        }
      }
    }
    //count all electrons
    CountElectronsAndCharges();
    if (verbose > 0) {PrintElectronPairing();}
    //and couple them
    CoupleElectrons(flatcenter);
    for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
      sum_unpaired_electrons += UnpairedElectrons[idAtm];
    }
    if (sum_unpaired_electrons > 0.0) {
      std::cout << "frustration" << std::endl;
      bool found = false;
      for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
        if (UnpairedElectrons[idAtm] > 0.0) {
          std::cout << "  " << idAtm + 1 << std::endl;
          //look for a suitable neighbour
          for (size_t idneigh = 0; idneigh < maxneigh; ++idneigh) {
            atomC = connectivity(idAtm + 1,idneigh + 1) - 1;
            if (atomC < 0) {break;}
            if ((modVSEPR(atomC + 1,1) == 3.0)&&(modVSEPR(atomC + 1,2) == 0.0)&&(atoms[atomC] == 7)&&(connectivity(atomC + 1,3) == 0)) {
              //it's conversion time
              UnpairedElectrons[idAtm] = 0.0;
              SetSP2_21(atomC);
            }
          }
          //find a terminal sp2 atom
          found = false;
          for (size_t idBtm = 0; idBtm < Natoms; ++idBtm) {
            if (idAtm == idBtm) {continue;}
            if (((modVSEPR(idBtm + 1,1) == 3.0)&&(modVSEPR(idBtm + 1,2) == 0.0)&&(atoms[idBtm] == 6))||((modVSEPR(idBtm + 1,1) == 2.0)&&(modVSEPR(idBtm + 1,2) == 1.0)&&(atoms[idBtm] == 7))) {
              nsp2 = 0;
              for (size_t idneigh = 0; idneigh < maxneigh; ++idneigh) {
                atomC = connectivity(idBtm + 1,idneigh + 1) - 1;
                if (atomC < 0) {break;}
                if (((modVSEPR(atomC + 1,1) == 3.0)&&(modVSEPR(atomC + 1,2) == 0.0)&&(atoms[atomC] == 6))||((modVSEPR(atomC + 1,1) == 2.0)&&(modVSEPR(atomC + 1,2) == 1.0)&&(atoms[atomC] == 7))) {++nsp2;}
              }
              if (nsp2 == 1) {
                //terminal sp2, check whether they are connected
                lessthansp3 = sp2Connected(idAtm,idBtm,connectionpath);
                if ((connectionpath.size() > 0)&&(connectionpath.size()%2 != 0)) {
                  //change hybridization of terminal atom
                  if (atoms[idBtm] == 6) {SetSP3_40(idBtm);}
                  else if (atoms[idBtm] == 7) {SetPlanar(idBtm);}
                  found = ReduceUnpairedElectrons(idAtm,sum_unpaired_electrons);
                  break;
                }
                else if (connectionpath.size() == 2) {
                  //no problem
                  found = ReduceUnpairedElectrons(idAtm,sum_unpaired_electrons);
                  break;
                }
              }
            }
            else if ((modVSEPR(idBtm + 1,1) == 3.0)&&(modVSEPR(idBtm + 1,2) == 0.0)&&(atoms[idBtm] == 7)&&(connectivity(idBtm + 1,3) == 0)) {
              nsp2 = 0;
              for (size_t idneigh = 0; idneigh < maxneigh; ++idneigh) {
                atomC = connectivity(idBtm + 1,idneigh + 1) - 1;
                if (atomC < 0) {break;}
                if (((modVSEPR(atomC + 1,1) == 3.0)&&(modVSEPR(atomC + 1,2) == 0.0)&&(atoms[atomC] == 6))||((modVSEPR(atomC + 1,1) == 2.0)&&(modVSEPR(atomC + 1,2) == 1.0)&&(atoms[atomC] == 7))) {++nsp2;}
              }
              if (nsp2 == 1) {
                //convert nitrogen to sp2?
                lessthansp3 = sp2Connected(idAtm,idBtm,connectionpath);
                if ((connectionpath.size() > 0)&&(connectionpath.size()%2 != 0)) {
                  //change hybridization of nitrogen
                  SetSP2_21(idBtm);
                  found = ReduceUnpairedElectrons(idAtm,sum_unpaired_electrons);
                  break;
                }
              }
            }
          }
          if (!found) {
            if (atoms[idAtm] == 6) {SetSP3_40(idAtm);}
            else if (atoms[idAtm] == 7) {SetSP3_31(idAtm);}
            ReduceUnpairedElectrons(idAtm,sum_unpaired_electrons);
          }
        }
        if (sum_unpaired_electrons < 0.1) {break;}
      }
    }
    //final checks
    for (size_t iring = 0; iring < therings.cols(); ++iring) {
      theringend = RingSize(iring) + 1;
      if (theringend == 5) {
        nC = 0;         //sp3 carbons
        nCar = 0;       //number of aromatic carbons
        nsp3 = 0;       //sp2 carbons bound to 2,2 oxygen that has one neighbour and bound to the nitrogen
        nN = 0;         //number of planar nitrogens
        nsp2 = 0;       //number of nitrogens bound to one aromatic carbon
        for (size_t idAtm = 0; idAtm < theringend; ++idAtm) {
          atomB = therings(idAtm + 1,iring + 1) - 1;
          if ((atoms[atomB] == 6)&&((modVSEPR(atomB + 1,1) == 4.0)&&(modVSEPR(atomB + 1,2) == 0.0))) {++nC;}
          if ((atoms[atomB] == 6)&&((modVSEPR(atomB + 1,1) == 3.0)&&(modVSEPR(atomB + 1,2) == 0.0))) {
            ++nCar;
            if (neighbours[atomB] == 3) {
              set3Neighbours(atomB,atomC,atomD,atomE);
              AminolCounter(atomC,atomD,atomE,nsp3);
              AminolCounter(atomD,atomC,atomE,nsp3);
              AminolCounter(atomE,atomC,atomD,nsp3);
            }
          }
          else if ((atoms[atomB] == 7)&&(neighbours[atomB] == 2)&&(modVSEPR(atomB + 1,1) == 2.0)&&(modVSEPR(atomB + 1,2) == 1.0)) {
            set3Neighbours(atomB,atomC,atomD,atomE);
            if ((atoms[atomD] == 6)&&(modVSEPR(atomD + 1,1) == 3.0)&&(modVSEPR(atomD + 1,2) == 0.0)) {++nsp2;}
            else if ((atoms[atomC] == 6)&&(modVSEPR(atomC + 1,1) == 3.0)&&(modVSEPR(atomC + 1,2) == 0.0)) {++nsp2;}
          }
          else if ((atoms[atomB] == 7)&&(modVSEPR(atomB + 1,1) == 3.0)&&(modVSEPR(atomB + 1,2) == 0.0)) {++nN;}
        }
        if ((nC == 1)&&(nCar == 3)&&(nsp2 == 1)) {
          //fix indole
          for (size_t idAtm = 0; idAtm < theringend; ++idAtm) {
            atomB = therings(idAtm + 1,iring + 1) - 1;
            if ((atoms[atomB] == 6)&&((modVSEPR(atomB + 1,1) == 4.0)&&(modVSEPR(atomB + 1,2) == 0.0))) {SetPlanar(atomB);}
            else if ((atoms[atomB] == 7)&&(neighbours[atomB] == 2)&&(modVSEPR(atomB + 1,1) == 2.0)&&(modVSEPR(atomB + 1,2) == 1.0)) {SetPlanar(atomB);}
          }
        }
        else if ((nsp3 == 2)&&(nN == 1)) {
          //change the oxygens
          for (size_t idAtm = 0; idAtm < theringend; ++idAtm) {
            atomB = therings(idAtm + 1,iring + 1) - 1;
            if ((atoms[atomB] == 6)&&((modVSEPR(atomB + 1,1) == 3.0)&&(modVSEPR(atomB + 1,2) == 0.0))) {
              if (neighbours[atomB] == 3) {
                set3Neighbours(atomB,atomC,atomD,atomE);
                ChangeAminol(atomC,atomD,atomE);
                ChangeAminol(atomD,atomC,atomE);
                ChangeAminol(atomE,atomC,atomD);
              }
            }
          }
        }
      }
      else if (theringend == 6) {
        nsp3 = 0;       //sp2 carbons bound to 2,2 oxygen that has one neighbour and bound to the nitrogen
        nN = 0;         //number of planar nitrogens
        for (size_t idAtm = 0; idAtm < theringend; ++idAtm) {
          atomB = therings(idAtm + 1,iring + 1) - 1;
          if ((atoms[atomB] == 6)&&((modVSEPR(atomB + 1,1) == 3.0)&&(modVSEPR(atomB + 1,2) == 0.0))) {
            if (neighbours[atomB] == 3) {
              set3Neighbours(atomB,atomC,atomD,atomE);
              AminolCounter(atomC,atomD,atomE,nsp3);
              AminolCounter(atomD,atomC,atomE,nsp3);
              AminolCounter(atomE,atomC,atomD,nsp3);
            }
          }
          else if ((atoms[atomB] == 7)&&(modVSEPR(atomB + 1,1) == 3.0)&&(modVSEPR(atomB + 1,2) == 0.0)) {++nN;}
        }
        if ((nsp3 == 2)&&(nN == 2)) {
          //change the oxygens
          for (size_t idAtm = 0; idAtm < theringend; ++idAtm) {
            atomB = therings(idAtm + 1,iring + 1) - 1;
            if ((atoms[atomB] == 6)&&((modVSEPR(atomB + 1,1) == 3.0)&&(modVSEPR(atomB + 1,2) == 0.0))) {
              if (neighbours[atomB] == 3) {
                set3Neighbours(atomB,atomC,atomD,atomE);
                ChangeAminol(atomC,atomD,atomE);
                ChangeAminol(atomD,atomC,atomE);
                ChangeAminol(atomE,atomC,atomD);
              }
            }
          }
        }
      }
    }
    for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
      //thiazole check
      if ((atoms[idAtm] == 6)&&(neighbours[idAtm] <= 3)) {
        if ((modVSEPR(idAtm + 1,1) == 2.0)&&((modVSEPR(idAtm + 1,2) == 0.0))) {continue;}
        nN = 0;             //number planar nitrogens
        nO = 0;             //number chalcogens
        undefO = 0;         //number carbonyls
        nflattorsions = 0;  //number of non-planar nitrogens
        nneighB = 0;        //number of sp2 neighbours
        //count some properties
        for (size_t idneigh = 0; idneigh < neighbours[idAtm]; ++idneigh) {
          atomB = connectivity(idAtm + 1,idneigh + 1) - 1;
          if ((atoms[atomB] == 6)&&(modVSEPR(atomB + 1,1) == 3.0)&&(modVSEPR(atomB + 1,2) == 0.0)) {++nneighB;}
          else if ((atoms[atomB] == 7)&&(modVSEPR(atomB + 1,1) == 2.0)&&(modVSEPR(atomB + 1,2) == 1.0)) {++nneighB;}
          else if (atoms[atomB] == 7) {
            if ((modVSEPR(atomB + 1,1) == 3.0)&&(modVSEPR(atomB + 1,2) == 0.0)) {
              ++nN;
              ++nneighB;
              atomF = atomB;
            }
            else {++nflattorsions;}
          }
          else if ((atoms[atomB] == 8)&&(modVSEPR(atomB + 1,1) == 2.0)&&(modVSEPR(atomB + 1,2) == 2.0)) {
            ++nO;
            atomG = atomB;
          }
          else if ((atoms[atomB] == 16)&&(modVSEPR(atomB + 1,1) == 2.0)&&(modVSEPR(atomB + 1,2) == 2.0)) {
            ++nO;
            atomG = atomB;
          }
          else if ((atoms[atomB] == 8)&&(modVSEPR(atomB + 1,1) == 1.0)&&(modVSEPR(atomB + 1,2) == 0.0)) {++undefO;}
          else if ((atoms[atomB] == 16)&&(modVSEPR(atomB + 1,1) == 1.0)&&(modVSEPR(atomB + 1,2) == 0.0)) {++undefO;}
        }
        if ((formalcharges[idAtm] == 1)&&(nO > 0)&&(undefO == 0)) {
          //make double bond
          formalcharges[idAtm] = 0;
          SetPlanar(idAtm);
          for (size_t idneigh = 0; idneigh < neighbours[idAtm]; ++idneigh) {
            atomB = connectivity(idAtm + 1,idneigh + 1) - 1;
            if ((atoms[atomB] == 8)&&(modVSEPR(atomB + 1,1) == 2.0)&&(modVSEPR(atomB + 1,2) == 2.0)) {
              SetTerminal(atomB);
              break;
            }
          }
        }
        if ((nN == 1)&&(nO == 1)&&(undefO == 0)) {
          //potential thiazole
          theringend = -1;
          for (size_t iring = 0; iring < therings.cols(); ++iring) {
            theringend = therings.rows() - 1;
            DisArom = 0;                          //check whether the 3 atoms of interest are in the ring
            for (size_t idCtm = 0; idCtm < theringend; ++idCtm) {
              atomE = therings(idCtm + 1,iring + 1) - 1;
              if (atomE < 0) {
                theringend = idCtm;
                break;
              }
              DisArom += (atomE == idAtm) + (atomE == atomF) + (atomE == atomG);
            }
            if (neighbours[idAtm] == 3) {
              set3Neighbours(idAtm,atomC,atomD,atomE);
              nsp2 = CheckOOP(idAtm + 1,atomC + 1,atomD + 1,atomE + 1,21.0,16);
            }
            else {nsp2 = 4;}
            if ((DisArom == 3)&&(theringend == 5)&&(nsp2 == 4)) {
              //thiazole found
              SetPlanar(idAtm);
              for (size_t idneigh = 0; idneigh < neighbours[idAtm]; ++idneigh) {
                atomB = connectivity(idAtm + 1,idneigh + 1) - 1;
                if ((atoms[atomB] == 7)&&(atomB == atomF)) {
                  if (neighbours[atomB] == 3) {formalcharges[idAtm] = 1;}
                  else {SetSP2_21(atomB);}
                  break;
                }
              }
              break;
            }
          }
        }
        else if ((undefO > 0)&&(nflattorsions > 0)) {
          //fix amide
          for (size_t idneigh = 0; idneigh < neighbours[idAtm]; ++idneigh) {
            atomB = connectivity(idAtm + 1,idneigh + 1) - 1;
            if ((atoms[atomB] == 7)&&(((modVSEPR(atomB + 1,1) == 2.0)&&(modVSEPR(atomB + 1,2) == 1.0))||((modVSEPR(atomB + 1,1) == 3.0)&&(modVSEPR(atomB + 1,2) == 1.0))||((modVSEPR(atomB + 1,1) == 4.0)&&(modVSEPR(atomB + 1,2) == 0.0)))) {
              SetPlanar(atomB);
              break;
            }
          }
        }
      }
      //check special nitrogens
      if ((atoms[idAtm] == 7)&&(modVSEPR(idAtm + 1,1) == 3.0)&&(modVSEPR(idAtm + 1,2) == 0.0)) {FixNplanar(idAtm);}
      if ((atoms[idAtm] == 7)&&(modVSEPR(idAtm + 1,1) == 3.0)&&(modVSEPR(idAtm + 1,2) == 1.0)) {
        FixN31Next2Sp2(idAtm,neighbours[idAtm]);
        if (neighbours[idAtm] == 1) {
          atomB = connectivity(idAtm + 1,1) - 1;
          if ((atoms[atomB] == 16)&&(neighbours[atomB] == 4)) {
            set3Neighbours(atomB,atomC,atomD,atomE);
            atomF = connectivity(atomB + 1,4) - 1;
            nN = (atoms[atomC] == 7)*(neighbours[atomC] <= 2) + (atoms[atomD] == 7)*(neighbours[atomD] <= 2) + (atoms[atomE] == 7)*(neighbours[atomE] <= 2) + (atoms[atomF] == 7)*(neighbours[atomF] <= 2);
            nO = (atoms[atomC] == 8)*(neighbours[atomC] == 1) + (atoms[atomD] == 8)*(neighbours[atomD] == 1) + (atoms[atomE] == 8)*(neighbours[atomE] == 1) + (atoms[atomF] == 8)*(neighbours[atomF] == 1);
            if (((nN == 2)&&(nO == 0))||((nN == 1)&&(nO == 1))) {SetSP2_21(idAtm);}
          }
        }
      }
      //if at this point something is undefined, just make it sp3
      if ((modVSEPR(idAtm + 1,1) == 0.0)&&(modVSEPR(idAtm + 1,2) == 0.0)) {
        if (atoms[idAtm] == 6) {
          if ((neighbours[idAtm] == 2)&&(radian2degree*fabs(Angle(connectivity(idAtm + 1,1),idAtm + 1,connectivity(idAtm + 1,2),geom)) > 160.0)) {
            if (((atoms[connectivity(idAtm + 1,1) - 1] == 6)&&(neighbours[connectivity(idAtm + 1,1) - 1] < 3))||((atoms[connectivity(idAtm + 1,1) - 1] == 7)&&(neighbours[connectivity(idAtm + 1,1) - 1] < 2))) {
              SetSP_20(idAtm);
              SetSP_20(connectivity(idAtm + 1,1) - 1);
            }
            else if (((atoms[connectivity(idAtm + 1,2) - 1] == 6)&&(neighbours[connectivity(idAtm + 1,2) - 1] < 3))||((atoms[connectivity(idAtm + 1,2) - 1] == 7)&&(neighbours[connectivity(idAtm + 1,2) - 1] < 2))) {
              SetSP_20(idAtm);
              SetSP_20(connectivity(idAtm + 1,2) - 1);
            }
            else {SetSP3_40(idAtm);}
          }
          else {SetSP3_40(idAtm);}
        }
        else if (atoms[idAtm] == 7) {SetSP3_31(idAtm);}
        else if (atoms[idAtm] == 8) {SetSP3_22(idAtm);}
        else if (atoms[idAtm] == 15) {SetSP3_31(idAtm);}
        else if (atoms[idAtm] == 16) {SetSP3_22(idAtm);}
        continue;
      }
      //check for distorted carboxyls and missing double bonds
      if ((atoms[idAtm] == 6)&&(neighbours[idAtm] == 3)) {
        nO = 0;           //terminal oxygens
        nC = 0;           //number of sp3 carbons
        nN = 0;           //planar nitrogens
        int n22 = 0;      //number of 2,2 atoms
        for (size_t idneigh = 0; idneigh < neighbours[idAtm]; ++idneigh) {
          atomB = connectivity(idAtm + 1,idneigh + 1) - 1;
          if ((atoms[atomB] == 8)&&(neighbours[atomB] == 1)) {++nO;}
          else if ((atoms[atomB] == 6)&&(modVSEPR(atomB + 1,1) == 4.0)&&(modVSEPR(atomB + 1,2) == 0.0)) {++nC;}
          else if ((atoms[atomB] == 7)&&(modVSEPR(atomB + 1,1) == 3.0)&&(modVSEPR(atomB + 1,2) == 0.0)) {++nN;}
          else if ((modVSEPR(atomB + 1,1) == 2.0)&&(modVSEPR(atomB + 1,2) == 2.0)) {++n22;}
        }
        set3Neighbours(idAtm,atomB,atomC,atomD);
        nsp2 = CheckOOP(idAtm + 1,atomB + 1,atomC + 1,atomD + 1,25.0*(nO == 2) + 10.0*(nO == 1)*(nC == 2),0.0);
        if ((nO == 2)&&(nsp2 == 4)) {
          MissedCarbonyl(idAtm,atomB,atomC,atomD);
          formalcharges[idAtm] = -1;
        }
        else if ((nO == 1)&&(nC == 2)&&(nsp2 == 4)) {MissedCarbonyl(idAtm,atomB,atomC,atomD);}
        else if (nC == 3) {
          auxvar = 360.0 - fabs(radian2degree*(Angle(connectivity(idAtm + 1,1),idAtm + 1,connectivity(idAtm + 1,2),geom) + Angle(connectivity(idAtm + 1,1),idAtm + 1,connectivity(idAtm + 1,3),geom) + Angle(connectivity(idAtm + 1,2),idAtm + 1,connectivity(idAtm + 1,3),geom)));
          if (fabs(auxvar) < 1.0) {
            bool modified = false;
            for (size_t idneigh = 0; idneigh < neighbours[idAtm]; ++idneigh) {
              atomB = connectivity(idAtm + 1,idneigh + 1) - 1;
              if ((neighbours[atomB] == 3)&&(modVSEPR(atomB + 1,1) == 4.0)&&(modVSEPR(atomB + 1,2) == 0.0)) {
                //check the angles
                auxvar = 360.0 - fabs(radian2degree*(Angle(connectivity(atomB + 1,1),atomB + 1,connectivity(atomB + 1,2),geom) + Angle(connectivity(atomB + 1,1),atomB + 1,connectivity(atomB + 1,3),geom) + Angle(connectivity(atomB + 1,2),atomB + 1,connectivity(atomB + 1,3),geom)));
                if (fabs(auxvar) < 1.0) {
                  SetPlanar(idAtm);
                  SetPlanar(atomB);
                  modified = true;
                  break;
                }
              }
            }
            if (!modified) {
              if (neighbours[connectivity(idAtm + 1,1) - 1] == 1) {
                SetPlanar(idAtm);
                SetPlanar(connectivity(idAtm + 1,1) - 1);
              }
              else if (neighbours[connectivity(idAtm + 1,2) - 1] == 1) {
                SetPlanar(idAtm);
                SetPlanar(connectivity(idAtm + 1,2) - 1);
              }
              else if (neighbours[connectivity(idAtm + 1,3) - 1] == 1) {
                SetPlanar(idAtm);
                SetPlanar(connectivity(idAtm + 1,3) - 1);
              }
            }
          }
        }
        else if ((modVSEPR(idAtm + 1,1) == 3.0)&&(modVSEPR(idAtm + 1,2) == 0.0)&&(nN + n22 == 3)) {
          formalcharges[idAtm] = 1;
        }
      }
      //bent double bonds are not supposed to be there
      if (((modVSEPR(idAtm + 1,1) == 3.0)&&(modVSEPR(idAtm + 1,2) == 0.0)&&(atoms[idAtm] == 6))||((modVSEPR(idAtm + 1,1) == 2.0)&&(modVSEPR(idAtm + 1,2) == 1.0)&&(atoms[idAtm] == 7))) {
        if (!IsCarbonyl(idAtm,neighbours[idAtm])) {
          nsp2 = 0;                                           //number of sp2 neighbours
          atomB = -1;                                         //sp2 neighbour
          for (size_t idneigh = 0; idneigh < neighbours[idAtm]; ++idneigh) {
            atomC = connectivity(idAtm + 1,idneigh + 1) - 1;
            if (atomC < 0) {break;}
            if (((modVSEPR(atomC + 1,1) == 3.0)&&(modVSEPR(atomC + 1,2) == 0.0)&&(atoms[atomC] == 6))||((modVSEPR(atomC + 1,1) == 2.0)&&(modVSEPR(atomC + 1,2) == 1.0)&&(atoms[atomC] == 7))) {
              if (IsCarbonyl(atomC,neighbours[atomC])) {continue;}
              else if (isringatom[atomC]) {continue;}
              ++nsp2;
              atomB = atomC;
            }
          }
          if ((nsp2 == 1)&&(!isringatom[idAtm])&&(!isringatom[atomB])) {
            //check that atomB only has one sp2 neighbour
            nsp2 = 0;
            for (size_t idneigh = 0; idneigh < neighbours[atomB]; ++idneigh) {
              atomC = connectivity(atomB + 1,idneigh + 1) - 1;
              if (atomC < 0) {break;}
              if (IsCarbonyl(atomC,neighbours[atomC])) {continue;}
              else if (isringatom[atomC]) {continue;}
              if (((modVSEPR(atomC + 1,1) == 3.0)&&(modVSEPR(atomC + 1,2) == 0.0)&&(atoms[atomC] == 6))||((modVSEPR(atomC + 1,1) == 2.0)&&(modVSEPR(atomC + 1,2) == 1.0)&&(atoms[atomC] == 7))) {++nsp2;}
            }
            if (nsp2 > 1) {continue;}
            //check the angles then
            nsp2 = 0;                 //distorted angles
            theringend = 0;           //angles checked
            int nOsp3A = 0;           //number of sp3 oxygens bound to A that can be sp2
            int nOsp3B = 0;           //number of sp3 oxygens bound to B that can be sp2
            int ringneighA = 0;       //number of ring neighbours of atom A that are also sp2
            int ringneighB = 0;       //number of ring neighbours of atom B that are also sp2
            int carbonylneighA = 0;   //number of carbonyl neighbours of A
            int carbonylneighB = 0;   //number of carbonyl neighbours of B
            for (size_t idneigh = 0; idneigh < neighbours[idAtm]; ++idneigh) {
              atomC = connectivity(idAtm + 1,idneigh + 1) - 1;
              if (atomC < 0) {break;}
              else if (atomC == atomB) {continue;}
              if ((atoms[atomC] == 8)&&(modVSEPR(atomC + 1,1) == 2.0)&&(modVSEPR(atomC + 1,2) == 2.0)&&(connectivity(atomC + 1,2) == 0)) {++nOsp3A;}
              if (isringatom[atomC]) {++ringneighA;}
              if (IsCarbonyl(atomC,neighbours[atomC])) {++carbonylneighA;}
              nOsp3B = 0;
              ringneighB = 0;
              carbonylneighB = 0;
              for (size_t idneigh2 = 0; idneigh2 < neighbours[atomB]; ++idneigh2) {
                atomD = connectivity(atomB + 1,idneigh2 + 1) - 1;
                if (atomD < 0) {break;}
                else if (atomD == atomA) {continue;}
                if ((atoms[atomD] == 8)&&(modVSEPR(atomD + 1,1) == 2.0)&&(modVSEPR(atomD + 1,2) == 2.0)&&(connectivity(atomD + 1,2) == 0)) {++nOsp3B;}
                if (isringatom[atomD]) {++ringneighB;}
                if (IsCarbonyl(atomD,neighbours[atomD])) {++carbonylneighB;}
                auxvar = fabs(radian2degree*Torsion(atomC + 1,idAtm + 1,atomB + 1,atomD + 1,geom));
                if (auxvar > 180.0) {auxvar -= 180.0;}
                if ((auxvar > 45.0)&&(auxvar < 180.0 - 45.0)) {++nsp2;}
                ++theringend;
              }
            }
            if (double(nsp2)/double(theringend) >= 0.5) {
              if ((nOsp3A > 0)&&(nOsp3B > 0)) {
                //make the oxygens sp2
                LoopTerminals(idAtm,neighbours[idAtm]);
                LoopTerminals(atomB,neighbours[atomB]);
              }
              else if ((nOsp3A == 0)&&(nOsp3B > 0)) {
                if (atoms[idAtm] == 6) {SetSP3_40(idAtm);}
                else if (atoms[idAtm] == 7) {
                  if (ringneighA + carbonylneighA > 0) {SetPlanar(idAtm);}
                  else {SetSP3_31(idAtm);}
                }
                LoopTerminals(atomB,neighbours[atomB]);
              }
              else if ((nOsp3A > 0)&&(nOsp3B == 0)) {
                if (atoms[atomB] == 6) {SetSP3_40(atomB);}
                else if (atoms[atomB] == 7) {
                  if (ringneighB + carbonylneighB > 0) {SetPlanar(atomB);}
                  else {SetSP3_31(atomB);}
                }
                LoopTerminals(idAtm,neighbours[idAtm]);
              }
              else {
                //make them sp3
                if (atoms[idAtm] == 6) {SetSP3_40(idAtm);}
                else if (atoms[idAtm] == 7) {
                  if (ringneighA + carbonylneighA > 0) {SetPlanar(idAtm);}
                  else {SetSP3_31(idAtm);}
                }
                if (atoms[atomB] == 6) {SetSP3_40(atomB);}
                else if (atoms[atomB] == 7) {
                  if (ringneighB + carbonylneighB > 0) {SetPlanar(atomB);}
                  else {SetSP3_31(atomB);}
                }
              }
            }
          }
        }
      }
    }
    if (verbose > 0) {PrintElectronPairing();}
  }
  void ProtonateFromVSEPR(double threshang = 6.3, int distance_filter = 1, int verbose = 1) {
    //function that protonates a ligand based on geometric arrangements
    this->RemoveProtons();
    this->Check4Metals();
    if (heavy_metal) {IsolateMetals();}
    this->VSEPRComponents(threshang,distance_filter,verbose);
    for (size_t idatm = 0; idatm < Natoms; ++idatm) {
      std::cout << atoms[idatm] << " -> " << modVSEPR(idatm + 1,1) << "," << modVSEPR(idatm + 1,2) << std::endl;
    }
    AddProtons(TheSystem,modVSEPR,connectivity);
    if (heavy_metal) {
      AddBackMetal();
      //remove clashing protons
      ReorderMoshPit(1.25);
      //get the VSEPR code for metals
      ClassifyMetals();
    }
    else {AddCompleteStructure();}
    ClassifyProtons();
    //for (size_t idatm = 0; idatm < Natoms; ++idatm) {
    //  std::cout << atoms[idatm] << " -> " << modVSEPR(idatm + 1,1) << "," << modVSEPR(idatm + 1,2) << "   " << formalcharges[idatm] << std::endl;
    //  //std::cout << atoms[idatm] << " -> " << modVSEPR(idatm + 1,1) << "," << modVSEPR(idatm + 1,2) << std::endl;
    //}
  }
};

#endif //_PATTERNS_
