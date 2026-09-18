/*ULYSSES, a semi-empirical package
    Copyright (C) 2023- Filipe Menezes (filipe.menezes@helmholtz-munich.de)
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

#ifndef _Energy_Decomposition_Analysis_
#define _Energy_Decomposition_Analysis_
#include <vector>
#include <string>
#include <math.h>
#include "QC.hpp"

//description:
//set of functions that decompose interaction energies in complexes

class EnergyDecomposition {
  Molecule SuperMolecule;                  //the molecular adduct
  std::vector<size_t> NoccMOs;
  std::vector<size_t> MonomerNAOs;
  std::vector<size_t> MonomerNatoms;
  std::vector<size_t> SuperMoleculeAtoms;
  std::vector<size_t> SuperMoleculeAOs;
  std::vector<int> MonomerNelectrons;
  std::vector<int> iauxvec;
  std::vector<Molecule> Monomers;
  std::vector<double> auxvec;
  std::vector<double> auxvec2;
  matrixE SuperMoleculeGeometry;
  //auxiliary matrices
  matrixE auxMat1;
  matrixE auxMat2;
  matrixE auxMat3;
  matrixE auxMat4;
  matrixE auxMat5;
  matrixE auxMat6;
  matrixE auxMat7;
  matrixE auxMat8;
  matrixE auxMat9;
  //matrices for infinitely separated molecules, this is what the zero stands for
  matrixE FockAB0_alpha;
  matrixE FockAB0_beta;
  matrixE HcoreAB0;
  matrixE DensAB0_alpha;
  matrixE DensAB0_beta;
  matrixE CMOAB0_alpha;
  matrixE CMOAB0_beta;
  matrixE EnuclearAB0;
  matrixE SAO0;
  std::vector<double> occMOAB0_alpha;
  std::vector<double> occMOAB0_beta;
  //matrices for equilibrium complex, this is what the int stands for
  matrixE FockABint_alpha;
  matrixE FockABint_beta;
  matrixE FockABexcint_alpha;
  matrixE FockABexcint_beta;
  matrixE HcoreABint;
  matrixE DensABint_alpha;
  matrixE DensABint_beta;
  matrixE EnuclearABint;
  matrixE SAOint;
  //matrices for the electrostatic contribution (es)
  matrixE FockABes_alpha;
  matrixE FockABes_beta;
  //HcoreABes = HcoreABint
  //DensABes_alpha = DensAB0_alpha
  //DensABes_beta = DensAB0_beta
  //EnuclearABes = EnuclearABint
  //matrices for the polarization contribution (pol)
  matrixE FockABpol_alpha;
  matrixE FockABpol_beta;
  matrixE FockABexcpol_alpha;
  matrixE FockABexcpol_beta;
  matrixE DensABpol_alpha;
  matrixE DensABpol_beta;
  //HcoreABpol = HcoreABint
  //EnuclearABpol = EnuclearABint
  //matrices for the overlap-repulsion contribution
  matrixE FprimeAB_alpha;
  matrixE FprimeAB_beta;
  matrixE PbarAB_alpha;
  matrixE PbarAB_beta;
  //quantities of interest
  matrixE Ees;
  matrixE Epol;
  matrixE Ect;
  matrixE Eexc_AB;
  matrixE Erep;
  matrixE Edisp;
  matrixE Esolv;
  //other stuff
  size_t Nmolecules;
  size_t SuperMoleculeNAtoms;
  int SuperMoleculeNelectrons;
  double electrostatic;
  double polarization;
  double chargetransfer;
  double interexchange;
  double overlaprepulsion;
  double dispersion;
  double solvation;
  bool openshell;
  bool usesolvation;
  bool usepopcontrain;           //determine whether to use population or physical constrain of electronic density for the polarization states
  bool fastedda;                 //skip polarisation terms and speed up?
  QCbasis * ElectronicStr;
public:
  EnergyDecomposition() {
    Nmolecules = 0;
    InitializeData();
  }
  EnergyDecomposition(QCbasis & ElecStruct) {
    Nmolecules = 0;
    setElectronicModel(ElecStruct);
    InitializeData();
  }
  //getters
  bool getPopulationConstraint() {return usepopcontrain;}
  bool useFastEDDA() {return fastedda;}
  void getVanDerVaartMerz(matrixE & esterms, matrixE & polterms, matrixE & ctterms) {
    //function that gets the van der Vaart and Merz interaction terms
    esterms = electrostatic;
    polterms = Epol + Erep;
    ctterms = Ect + Eexc_AB;
  }
  double getElectrostatic(matrixE & ESmat) {
    ESmat = Ees;
    return electrostatic;
  }
  double getPolarization(matrixE & POLmat) {
    POLmat = Epol;
    return polarization;
  }
  double getChargeTransfer(matrixE & CTmat) {
    CTmat = Ect;
    return chargetransfer;
  }
  double getExchangeInter(matrixE & EXCmat) {
    EXCmat = Eexc_AB;
    return interexchange;
  }
  double getExchange(matrixE & EXCmat) {
    EXCmat = Eexc_AB;
    return interexchange;
  }
  double getOvRep(matrixE & ORmat) {
    ORmat = Erep;
    return overlaprepulsion;
  }
  double getDispersion(matrixE & Dmat) {
    Dmat = Edisp;
    return dispersion;
  }
  double getSolvation(matrixE & Smat) {
    Smat = Esolv;
    return solvation;
  }
  //setters
  void setPopulationConstraint(bool newconstr) {usepopcontrain = newconstr;}
  void setFastEDDA(bool fast) {fastedda = fast;}
  void setElectronicModel(QCbasis & ElecStruct) {
    ElectronicStr = & ElecStruct;
    usesolvation = ElectronicStr->UseSolvation();
    ElecStruct.Calculate(0,0);                       //just to ensure we have meaningful shell types
    std::string shelltype = ElectronicStr->ShellType();
    openshell = ((shelltype == "open")||(shelltype == "UHF")||(shelltype == "uhf"));
  }
  //initialize quantities
  void InitializeData() {
    usepopcontrain = true;
    fastedda = false;
  }
  //resizing functions
  void ResetInteractionMatrices() {
    Ees.resize(SuperMoleculeNAtoms,SuperMoleculeNAtoms);
    Epol.resize(SuperMoleculeNAtoms,SuperMoleculeNAtoms);
    Ect.resize(SuperMoleculeNAtoms,SuperMoleculeNAtoms);
    Eexc_AB.resize(SuperMoleculeNAtoms,SuperMoleculeNAtoms);
    Erep.resize(SuperMoleculeNAtoms,SuperMoleculeNAtoms);
    Edisp.resize(SuperMoleculeNAtoms,SuperMoleculeNAtoms);
    for (size_t idAtm = 0; idAtm < SuperMoleculeNAtoms; ++idAtm) {
      Ees(idAtm + 1,idAtm + 1) = 0.0;
      Epol(idAtm + 1,idAtm + 1) = 0.0;
      Ect(idAtm + 1,idAtm + 1) = 0.0;
      Eexc_AB(idAtm + 1,idAtm + 1) = 0.0;
      Erep(idAtm + 1,idAtm + 1) = 0.0;
      Edisp(idAtm + 1,idAtm + 1) = 0.0;
      for (size_t idBtm = 0; idBtm < idAtm; ++idBtm) {
        Ees(idBtm + 1,idAtm + 1) = 0.0;
        Epol(idBtm + 1,idAtm + 1) = 0.0;
        Ect(idBtm + 1,idAtm + 1) = 0.0;
        Eexc_AB(idBtm + 1,idAtm + 1) = 0.0;
        Erep(idBtm + 1,idAtm + 1) = 0.0;
        Edisp(idBtm + 1,idAtm + 1) = 0.0;
        Ees(idAtm + 1,idBtm + 1) = 0.0;
        Epol(idAtm + 1,idBtm + 1) = 0.0;
        Ect(idAtm + 1,idBtm + 1) = 0.0;
        Eexc_AB(idAtm + 1,idBtm + 1) = 0.0;
        Erep(idAtm + 1,idBtm + 1) = 0.0;
        Edisp(idAtm + 1,idBtm + 1) = 0.0;
      }
    }
    if (usesolvation) {
      Esolv.resize(SuperMoleculeNAtoms,SuperMoleculeNAtoms);
      for (size_t idAtm = 0; idAtm < SuperMoleculeNAtoms; ++idAtm) {
        Esolv(idAtm + 1,idAtm + 1) = 0.0;
        for (size_t idBtm = 0; idBtm < idAtm; ++idBtm) {
          Esolv(idBtm + 1,idAtm + 1) = 0.0;
          Esolv(idAtm + 1,idBtm + 1) = 0.0;
        }
      }
    }
  }
  void ZeroPolarisation() {
    //function that zeroes the matrices related to the polarisation energy
    FockABpol_alpha = FockABes_alpha;
    DensABpol_alpha = DensAB0_alpha;
    //size_t numbrows = FockABes_alpha.rows();
    //FockABexcpol_alpha.resize(numbrows,numbrows);
    //for (size_t idrow = 0; idrow < numbrows; ++idrow) {
    //  FockABexcpol_alpha(idrow + 1,idrow + 1) = 0.0;
    //  for (size_t idcol = idrow + 1; idcol < numbrows; ++idcol) {
    //    FockABexcpol_alpha(idrow + 1,idcol + 1) = 0.0;
    //    FockABexcpol_alpha(idcol + 1,idrow + 1) = 0.0;
    //  }
    //}
    if (openshell) {
      FockABpol_beta = FockABes_beta;
      DensABpol_beta = DensAB0_beta;
      //numbrows = FockABes_beta.rows();
      //FockABexcpol_beta.resize(numbrows,numbrows);
      //for (size_t idrow = 0; idrow < numbrows; ++idrow) {
      //  FockABexcpol_beta(idrow + 1,idrow + 1) = 0.0;
      //  for (size_t idcol = idrow + 1; idcol < numbrows; ++idcol) {
      //    FockABexcpol_beta(idrow + 1,idcol + 1) = 0.0;
      //    FockABexcpol_beta(idcol + 1,idrow + 1) = 0.0;
      //  }
      //}
    }
  }
  void DecomposeEnergy(bool splitmonomers = true) {
    //main driver for the calculation
    std::string method = ElectronicStr->Type();
    bool NDDOmethod = ((method == "MNDO")||(method == "AM1")||(method == "PM3")||(method == "PM3PDDG")||(method == "MNDOPDDG")||(method == "PM3BP")||(method == "RM1")||(method == "MNDOd")||(method == "PM6"));
    bool GFNmethod = (method == "GFN2");
    this->ResetInteractionMatrices();
    if (NDDOmethod) {
      //for HF and derived methods
      std::cout << "calculating monomer wavefunctions" << std::endl;
      this->getQMMonomersHF(splitmonomers);
      std::cout << "calculating supermolecule wavefunction" << std::endl;
      this->getQMSuperMoleculeHF();
      std::cout << "calculating electrostatic contributions" << std::endl;
      this->getElectrostaticContributionsHF();
      if (!fastedda) {
        std::cout << "calculating polarization terms" << std::endl;
        this->getPolarizationContributionsHF();
      }
      else {ZeroPolarisation();}
      std::cout << "calculating overlap-repulsion" << std::endl;
      this->getOverlapRepulsionHF();
      std::cout << "finalizing with interaction energy calculation" << std::endl;
      this->getInteractionEnergiesHF();
    }
    else if (GFNmethod) {
      //for DFT derived methods
      std::cout << "calculating monomer wavefunctions" << std::endl;
      this->getQMMonomersDFT(splitmonomers);
      std::cout << "calculating supermolecule wavefunction" << std::endl;
      this->getQMSuperMoleculeDFT();
      std::cout << "calculating electrostatic contributions" << std::endl;
      this->getElectrostaticContributionsDFT();
      std::cout << "calculating overlap-repulsion" << std::endl;
      this->getOverlapRepulsionDFT();
      if (!fastedda) {
        std::cout << "calculating polarization terms" << std::endl;
        this->getPolarizationContributionsDFT();
      }
      else {ZeroPolarisation();}
      std::cout << "finalizing with interaction energy calculation" << std::endl;
      this->getInteractionEnergiesDFT();
    }
    this->calcTotalDecomposedEnergies();
    this->OutputEnergies();
  }
  //functions to get the molecules
  void getMoleculesFromQC() {
    //function that generates the molecules from the QC method supplied
    SuperMolecule = ElectronicStr->Component();
    SuperMoleculeAtoms = SuperMolecule.Atoms();
    auxMat1 = SuperMolecule.Geometry();
    int iatm;
    int counter = 0;
    int SuperCharge = SuperMolecule.Charge();
    int SuperMultiplicity = SuperMolecule.Multiplicity();
    std::vector<size_t> moleculesplitvec;
    std::vector<size_t> MonomerAtomList;
    std::vector<size_t> atoms = SuperMoleculeAtoms;
    matrixE SuperMoleculeGeom = auxMat1;
    SuperMoleculeNAtoms = SuperMolecule.Natoms();
    Nmolecules = SuperMolecule.MoleculeCounter(moleculesplitvec);
    if (Nmolecules > 2) {
      std::cout << "WARNING: EnergyDecompositionAnalysis.hpp: getMoleculesFromQC(): currently supported for two molecules only" << std::endl;
      Nmolecules = 2;
    }
    //resize vectors
    MonomerNAOs.resize(Nmolecules);
    MonomerNelectrons.resize(Nmolecules);
    NoccMOs.resize(Nmolecules);
    Monomers.resize(Nmolecules);
    MonomerNatoms.resize(Nmolecules);
    //assign the atoms to each monomer
    for (size_t idAB = 0; idAB < SuperMoleculeNAtoms; ++idAB) {
      if (moleculesplitvec[idAB] < 3) {++MonomerNatoms[moleculesplitvec[idAB] - 1];}
    }
    //now assemble the monomers
    for (size_t idmol = 0; idmol < Nmolecules; ++idmol) {
      iatm = 1;
      //resize containers
      MonomerAtomList.resize(MonomerNatoms[idmol]);
      auxMat2.resize(MonomerNatoms[idmol],3);
      //now build the monomer's atom list and geometry matrix
      for (size_t idatm = 0; idatm < SuperMoleculeNAtoms; ++idatm) {
        if (moleculesplitvec[idatm] == idmol + 1) {
          auxMat2(iatm,1) = auxMat1(idatm + 1,1);
          auxMat2(iatm,2) = auxMat1(idatm + 1,2);
          auxMat2(iatm,3) = auxMat1(idatm + 1,3);
          MonomerAtomList[iatm - 1] = atoms[idatm];
          ++iatm;
        }
      }
      //now get the supermolecule back in place, so that the supermolecule is basically ordered according to monomers
      for (size_t idatm = 0; idatm < MonomerNatoms[idmol]; ++idatm) {
        SuperMoleculeAtoms[idatm + counter] = MonomerAtomList[idatm];
        SuperMoleculeGeom(idatm + counter + 1,1) = auxMat2(idatm + 1,1);
        SuperMoleculeGeom(idatm + counter + 1,2) = auxMat2(idatm + 1,2);
        SuperMoleculeGeom(idatm + counter + 1,3) = auxMat2(idatm + 1,3);
      }
      counter += MonomerNatoms[idmol];
      Monomers[idmol].setGeometry(auxMat2);
      Monomers[idmol].setAtoms(MonomerAtomList);
      //extend here to include non-zero charges and more than singlet!
      Monomers[idmol].setCharge(0);
      Monomers[idmol].setMultiplicity(1);
      Monomers[idmol].calcNatoms();
      Monomers[idmol].CountElectrons();
      Monomers[idmol].masses();
    }
    //now define the full molecule
    SuperMolecule.setGeometry(SuperMoleculeGeom);
    SuperMolecule.setAtoms(SuperMoleculeAtoms);
    SuperMolecule.setCharge(SuperCharge);
    SuperMolecule.setMultiplicity(SuperMultiplicity);
    SuperMolecule.calcNatoms();
    SuperMolecule.CountElectrons();
    SuperMolecule.masses();
  }
  void addSuperMolecule(Molecule & newSuperMolecule) {
    //function that adds a new supermolecule
    SuperMolecule = newSuperMolecule;
    SuperMoleculeNAtoms = SuperMolecule.Natoms();
  }
  void addMonomer(Molecule & newMonomer) {
    //function that adds a new supermolecule
    Monomers.push_back(newMonomer);
    Nmolecules += 1;
    MonomerNAOs.resize(Nmolecules);
    MonomerNelectrons.resize(Nmolecules);
    NoccMOs.resize(Nmolecules);
    MonomerNatoms.resize(Nmolecules);
    MonomerNatoms[Nmolecules - 1] = newMonomer.Natoms();
  }
  void addMolecules(Molecule & newSuperMolecule, Molecule & newMonomerA, Molecule & newMonomerB) {
    this->addSuperMolecule(newSuperMolecule);
    this->addMonomer(newMonomerA);
    this->addMonomer(newMonomerB);
  }
  void getQMMonomersHF(bool split = true) {
    //function that does the QM calculations on the monomers and collects all the relevant data from these
    //initialization of some relevant information
    int nocc;
    int totalNAOs = 0;
    int totalATMs = 0;
    int numberAOmonomer;
    int numberATMmonomer;
    if (split) {
      //in this mode the monomers are completely split; this is advantageous if the monomers are very large, and calculation time might be saved
      for (size_t idmol = 0; idmol < Nmolecules; ++idmol) {
        //get the monomer wavefunction
        ElectronicStr->setMolecule(Monomers[idmol]);
        ElectronicStr->Calculate(0);
        //collect the relevant quantities
        auxMat1 = ElectronicStr->getFock();
        auxMat3 = ElectronicStr->getHcore();
        ElectronicStr->fetchDens(auxMat4);
        auxMat6 = ElectronicStr->EnuclearMat();
        auxMat7 = ElectronicStr->MOs();
        auxMat9 = ElectronicStr->getSAO();
        auxvec = ElectronicStr->occ();
        //count the AOs for this subsystem
        MonomerNAOs[idmol] = ElectronicStr->NAO();
        MonomerNelectrons[idmol] = Monomers[idmol].Nelectrons();
        iauxvec = ElectronicStr->CoreCharges();
        for (size_t idAtm = 0; idAtm < MonomerNatoms[idmol]; ++idAtm) {
          MonomerNelectrons[idmol] -= iauxvec[idAtm];
        }
        numberAOmonomer = totalNAOs + 1;                  //this is the initial position when we start to write into
        numberATMmonomer = totalATMs + 1;
        totalNAOs += MonomerNAOs[idmol];
        totalATMs += MonomerNatoms[idmol];
        FockAB0_alpha.resize(totalNAOs,totalNAOs);
        DensAB0_alpha.resize(totalNAOs,totalNAOs);
        CMOAB0_alpha.resize(totalNAOs,totalNAOs);
        HcoreAB0.resize(totalNAOs,totalNAOs);
        SAO0.resize(totalNAOs,totalNAOs);
        occMOAB0_alpha.resize(totalNAOs);
        nocc = auxvec.size();
        NoccMOs[idmol] = nocc;
        //now copy these quantities into their respective containers
        for (size_t idAO = 0; idAO < MonomerNAOs[idmol]; ++idAO) {
          //technically I could save some calculation time by doing only half of the matrices and then symmetrize
          //however, the MO matrix requires full loops
          for (size_t idBO = 0; idBO < MonomerNAOs[idmol]; ++idBO) {
            FockAB0_alpha(numberAOmonomer + idAO,numberAOmonomer + idBO) = auxMat1(idAO + 1,idBO + 1);
            DensAB0_alpha(numberAOmonomer + idAO,numberAOmonomer + idBO) = auxMat4(idAO + 1,idBO + 1);
            CMOAB0_alpha(numberAOmonomer + idAO,numberAOmonomer + idBO) = auxMat7(idAO + 1,idBO + 1);
            HcoreAB0(numberAOmonomer + idAO,numberAOmonomer + idBO) = auxMat3(idAO + 1,idBO + 1);
            SAO0(numberAOmonomer + idAO,numberAOmonomer + idBO) = auxMat9(idAO + 1,idBO + 1);
          }
          if (idAO + 1 <= nocc) {occMOAB0_alpha[numberAOmonomer + idAO - 1] = auxvec[idAO];}
          else {occMOAB0_alpha[numberAOmonomer + idAO - 1] = 0.0;}
        }
        EnuclearAB0.resize(totalATMs,totalATMs);
        for (size_t idAtm = 0; idAtm < MonomerNatoms[idmol]; ++idAtm) {
          for (size_t idBtm = 0; idBtm < MonomerNatoms[idmol]; ++idBtm) {
            EnuclearAB0(numberATMmonomer + idAtm,numberATMmonomer + idBtm) = auxMat6(idAtm + 1,idBtm + 1);
          }
        }
        if (openshell) {
          auxMat2 = ElectronicStr->getFockb();
          ElectronicStr->fetchbDens(auxMat5);
          auxMat8 = ElectronicStr->bMOs();
          FockAB0_beta.resize(totalNAOs,totalNAOs);
          DensAB0_beta.resize(totalNAOs,totalNAOs);
          CMOAB0_beta.resize(totalNAOs,totalNAOs);
          occMOAB0_beta.resize(totalNAOs);
          auxvec = ElectronicStr->bocc();
          for (size_t idAO = 0; idAO < MonomerNAOs[idmol]; ++idAO) {
            for (size_t idBO = 0; idBO < MonomerNAOs[idmol]; ++idBO) {
              FockAB0_beta(numberAOmonomer + idAO,numberAOmonomer + idBO) = auxMat2(idAO + 1,idBO + 1);
              DensAB0_beta(numberAOmonomer + idAO,numberAOmonomer + idBO) = auxMat5(idAO + 1,idBO + 1);
              CMOAB0_beta(numberAOmonomer + idAO,numberAOmonomer + idBO) = auxMat8(idAO + 1,idBO + 1);
            }
            occMOAB0_beta[numberAOmonomer + idAO - 1] = auxvec[idAO];
          }
        }
      }
    }
    else {
      //in this calculation mode the monomers are placed at a very large distance
      //this is to be used in case of non-singlet states or non-trivial spin multiplicities
      double shift = 1000.0;
      int idBtm;
      Molecule InfdistComplex = SuperMolecule;                //note that charge and spin come from this species
      auxMat1 = SuperMolecule.Geometry();
      auxMat2 = Monomers[0].Geometry();
      //set the first set of atoms
      for (size_t idAtm = 0; idAtm < MonomerNatoms[0]; ++idAtm) {
        auxMat1(idAtm + 1,1) = auxMat2(idAtm + 1,1);
        auxMat1(idAtm + 1,2) = auxMat2(idAtm + 1,2);
        auxMat1(idAtm + 1,3) = auxMat2(idAtm + 1,3);
      }
      ElectronicStr->setMolecule(Monomers[0]);
      ElectronicStr->Calculate(0,0);
      MonomerNAOs[0] = ElectronicStr->NAO();
      //now assemble the monomers
      for (size_t idmol = 1; idmol < Nmolecules; ++idmol) {
        ElectronicStr->setMolecule(Monomers[idmol]);
        ElectronicStr->Calculate(0,0);
        MonomerNAOs[idmol] = ElectronicStr->NAO();
        auxMat2 = Monomers[idmol].Geometry();
        idBtm = 0;
        for (size_t idAtm = MonomerNatoms[idmol - 1]; idAtm < SuperMoleculeNAtoms; ++idAtm,++idBtm) {
          auxMat1(idAtm + 1,1) = auxMat2(idBtm + 1,1) + shift;
          auxMat1(idAtm + 1,2) = auxMat2(idBtm + 1,2) + shift;
          auxMat1(idAtm + 1,3) = auxMat2(idBtm + 1,3) + shift;
        }
      }
      InfdistComplex.setGeometry(auxMat1);
      ElectronicStr->setMolecule(InfdistComplex);
      ElectronicStr->Calculate(0);
      //collect the relevant quantities
      MonomerNelectrons[0] = Monomers[0].Nelectrons();
      MonomerNelectrons[1] = Monomers[1].Nelectrons();
      iauxvec = ElectronicStr->CoreCharges();
      for (size_t idAtm = 0; idAtm < MonomerNatoms[0]; ++idAtm) {
        MonomerNelectrons[0] -= iauxvec[idAtm];
      }
      for (int idAtm = 0; idAtm < MonomerNatoms[1]; ++idAtm) {
        MonomerNelectrons[1] -= iauxvec[MonomerNatoms[0] + idAtm];
      }
      FockAB0_alpha = ElectronicStr->getFock();
      HcoreAB0 = ElectronicStr->getHcore();
      ElectronicStr->fetchDens(DensAB0_alpha);
      EnuclearAB0 = ElectronicStr->EnuclearMat();
      SAO0 = ElectronicStr->getSAO();
      CMOAB0_alpha = ElectronicStr->MOs();
      occMOAB0_alpha = ElectronicStr->occ();
      nocc = occMOAB0_alpha.size();
      if (nocc != ElectronicStr->NAO()) {
        occMOAB0_alpha.resize(ElectronicStr->NAO());
        for (size_t idAO = nocc; idAO < ElectronicStr->NAO(); ++idAO) {
          occMOAB0_alpha[idAO] = 0.0;
        }
      }
      if (openshell) {
        FockAB0_beta = ElectronicStr->getFockb();
        ElectronicStr->fetchbDens(DensAB0_beta);
        CMOAB0_beta = ElectronicStr->bMOs();
        occMOAB0_beta = ElectronicStr->bocc();
        nocc = occMOAB0_beta.size();
        if (nocc != ElectronicStr->NAO()) {
          occMOAB0_beta.resize(ElectronicStr->NAO());
          for (size_t idAO = nocc; idAO < ElectronicStr->NAO(); ++idAO) {
            occMOAB0_beta[idAO] = 0.0;
          }
        }
      }
    }
  }
  void getQMMonomersDFT(bool split = true) {
    //function that does the QM calculations on the monomers and collects all the relevant data from these
    //initialization of some relevant information
    int nocc;
    int totalNAOs = 0;
    int totalATMs = 0;
    int numberAOmonomer;
    int numberATMmonomer;
    int counter;
    int numberShellsmonomer = 1;
    int nshells;
    if (split) {
      //in this mode the monomers are completely split
      //this is advantageous if the monomers are very large, and calculation time might be saved
      for (size_t idmol = 0; idmol < Nmolecules; ++idmol) {
        //get the monomer wavefunction
        ElectronicStr->setMolecule(Monomers[idmol]);
        ElectronicStr->Calculate(0);
        //collect the relevant quantities
        auxMat1 = Monomers[idmol].Geometry();
        ElectronicStr->EnergyDispersion(auxMat1,auxMat2,1600.0);
        ElectronicStr->Enuclear(auxMat1,auxMat3,1600.0);
        ElectronicStr->EnergyMatrix(auxMat1,auxMat4,false);
        if (usesolvation) {auxMat8 = ElectronicStr->CurrentSolvationEnergyMatrix(auxMat1);}
        ElectronicStr->fetchDens(auxMat5);
        auxMat3 = ElectronicStr->getHcore();
        auxMat6 = ElectronicStr->MOs();
        auxMat7 = ElectronicStr->getSAO();
        auxvec = ElectronicStr->occ();
        //count the AOs for this subsystem
        MonomerNAOs[idmol] = ElectronicStr->NAO();
        MonomerNelectrons[idmol] = Monomers[idmol].Nelectrons();
        iauxvec = ElectronicStr->CoreCharges();
        for (size_t idAtm = 0; idAtm < MonomerNatoms[idmol]; ++idAtm) {
          MonomerNelectrons[idmol] -= iauxvec[idAtm];
        }
        numberAOmonomer = totalNAOs + 1;                  //this is the initial position when we start to write into
        numberATMmonomer = totalATMs + 1;
        totalNAOs += MonomerNAOs[idmol];
        totalATMs += MonomerNatoms[idmol];
        DensAB0_alpha.resize(totalNAOs,totalNAOs);
        HcoreAB0.resize(totalNAOs,totalNAOs);
        CMOAB0_alpha.resize(totalNAOs,totalNAOs);
        SAO0.resize(totalNAOs,totalNAOs);
        occMOAB0_alpha.resize(totalNAOs);
        nocc = auxvec.size();
        NoccMOs[idmol] = nocc;
        //now copy these quantities into their respective containers
        for (size_t idAO = 0; idAO < MonomerNAOs[idmol]; ++idAO) {
          //technically I could save some calculation time by doing only half of the matrices and then symmetrize
          //however, the MO matrix requires full loops
          for (size_t idBO = 0; idBO < MonomerNAOs[idmol]; ++idBO) {
            DensAB0_alpha(numberAOmonomer + idAO,numberAOmonomer + idBO) = auxMat5(idAO + 1,idBO + 1);
            HcoreAB0(numberAOmonomer + idAO,numberAOmonomer + idBO) = auxMat3(idAO + 1,idBO + 1);
            CMOAB0_alpha(numberAOmonomer + idAO,numberAOmonomer + idBO) = auxMat6(idAO + 1,idBO + 1);
            SAO0(numberAOmonomer + idAO,numberAOmonomer + idBO) = auxMat7(idAO + 1,idBO + 1);
          }
          if (idAO + 1 <= nocc) {occMOAB0_alpha[numberAOmonomer + idAO - 1] = auxvec[idAO];}
          else {occMOAB0_alpha[numberAOmonomer + idAO - 1] = 0.0;}
        }
        EnuclearAB0.resize(totalATMs,totalATMs);
        FockAB0_alpha.resize(totalATMs,totalATMs);    //store the energy matrix
        if (usesolvation) {
          for (size_t idAtm = 0; idAtm < MonomerNatoms[idmol]; ++idAtm) {
            for (size_t idBtm = 0; idBtm < MonomerNatoms[idmol]; ++idBtm) {
              Esolv(numberATMmonomer + idAtm,numberATMmonomer + idBtm) -= auxMat8(idAtm + 1,idBtm + 1);
            }
          }
        }
        nshells = 0;
        iauxvec = ElectronicStr->AtomicShells();
        for (size_t idAtm = 0; idAtm < MonomerNatoms[idmol]; ++idAtm) {
          for (size_t idBtm = 0; idBtm < MonomerNatoms[idmol]; ++idBtm) {
            EnuclearAB0(numberATMmonomer + idAtm,numberATMmonomer + idBtm) = auxMat3(idAtm + 1,idBtm + 1);
            FockAB0_alpha(numberATMmonomer + idAtm,numberATMmonomer + idBtm) = auxMat4(idAtm + 1,idBtm + 1);
            Edisp(numberATMmonomer + idAtm,numberATMmonomer + idBtm) -= auxMat2(idAtm + 1,idBtm + 1);
          }
          nshells += iauxvec[idAtm];
        }
        //now the first moments of multipolar expansion; fetch data
        auxvec = ElectronicStr->getQAtoms();
        auxMat3 = ElectronicStr->getQShells();
        auxMat6 = ElectronicStr->getmuAtoms();
        auxMat7 = ElectronicStr->getthetaAtoms();
        //resize arrays
        auxvec2.resize(totalATMs);                                         //charge/atom
        FockABexcpol_alpha.resize(numberShellsmonomer + nshells - 1,1);    //charge/shell
        FockABexcpol_beta.resize(totalATMs,3);                             //dipole moment
        DensABint_beta.resize(totalATMs,6);                                //quadrupole moment
        //write
        counter = 0;
        for (size_t idAtm = 0; idAtm < MonomerNatoms[idmol]; ++idAtm) {
          auxvec2[numberATMmonomer + idAtm - 1] = auxvec[idAtm];
          for (size_t idshell = 0; idshell < iauxvec[idAtm]; ++idshell,++counter) {
            FockABexcpol_alpha(numberShellsmonomer + counter,1) = auxMat3(counter + 1,1);
          }
          for (size_t idcoord = 1; idcoord < 4; ++idcoord) {
            FockABexcpol_beta(numberATMmonomer + idAtm,idcoord) = auxMat6(idAtm + 1,idcoord);
            DensABint_beta(numberATMmonomer + idAtm,idcoord) = auxMat7(idAtm + 1,idcoord);
            DensABint_beta(numberATMmonomer + idAtm,idcoord + 3) = auxMat7(idAtm + 1,idcoord + 3);
          }
        }
        numberShellsmonomer += nshells;
      }
    }
    else {
      //in this calculation mode the monomers are placed at a very large distance
      //this is to be used in case of non-singlet states or non-trivial spin multiplicities
      double shift = 1000.0;
      int idBtm;
      Molecule InfdistComplex = SuperMolecule;                //note that charge and spin come from this species
      auxMat1 = SuperMolecule.Geometry();
      auxMat2 = Monomers[0].Geometry();
      //set the first set of atoms
      for (size_t idAtm = 0; idAtm < MonomerNatoms[0]; ++idAtm) {
        auxMat1(idAtm + 1,1) = auxMat2(idAtm + 1,1);
        auxMat1(idAtm + 1,2) = auxMat2(idAtm + 1,2);
        auxMat1(idAtm + 1,3) = auxMat2(idAtm + 1,3);
      }
      ElectronicStr->setMolecule(Monomers[0]);
      ElectronicStr->Calculate(0,0);
      MonomerNAOs[0] = ElectronicStr->NAO();
      //now assemble the monomers
      for (size_t idmol = 1; idmol < Nmolecules; ++idmol) {
        ElectronicStr->setMolecule(Monomers[idmol]);
        ElectronicStr->Calculate(0,0);
        MonomerNAOs[idmol] = ElectronicStr->NAO();
        auxMat2 = Monomers[idmol].Geometry();
        idBtm = 0;
        for (size_t idAtm = MonomerNatoms[idmol - 1]; idAtm < SuperMoleculeNAtoms; ++idAtm,++idBtm) {
          auxMat1(idAtm + 1,1) = auxMat2(idBtm + 1,1) + shift;
          auxMat1(idAtm + 1,2) = auxMat2(idBtm + 1,2) + shift;
          auxMat1(idAtm + 1,3) = auxMat2(idBtm + 1,3) + shift;
        }
      }
      InfdistComplex.setGeometry(auxMat1);
      ElectronicStr->setMolecule(InfdistComplex);
      ElectronicStr->Calculate(0);
      //collect the relevant quantities
      MonomerNelectrons[0] = Monomers[0].Nelectrons();
      MonomerNelectrons[1] = Monomers[1].Nelectrons();
      iauxvec = ElectronicStr->CoreCharges();
      for (size_t idAtm = 0; idAtm < MonomerNatoms[0]; ++idAtm) {
        MonomerNelectrons[0] -= iauxvec[idAtm];
      }
      for (size_t idAtm = 0; idAtm < MonomerNatoms[1]; ++idAtm) {
        MonomerNelectrons[1] -= iauxvec[MonomerNatoms[0] + idAtm];
      }
      ElectronicStr->EnergyDispersion(auxMat1,auxMat2,1600.0);
      Edisp = auxMat2*(-1.0);
      if (usesolvation) {
        auxMat8 = ElectronicStr->CurrentSolvationEnergyMatrix(auxMat1);
        Esolv = auxMat8*(-1.0);
      }
      ElectronicStr->Enuclear(auxMat1,EnuclearAB0,1600.0);
      ElectronicStr->EnergyMatrix(auxMat1,FockAB0_alpha,false);
      HcoreAB0 = ElectronicStr->getHcore();
      ElectronicStr->fetchDens(DensAB0_alpha);
      SAO0 = ElectronicStr->getSAO();
      CMOAB0_alpha = ElectronicStr->MOs();
      occMOAB0_alpha = ElectronicStr->occ();
      nocc = occMOAB0_alpha.size();
      if (nocc != ElectronicStr->NAO()) {
        occMOAB0_alpha.resize(ElectronicStr->NAO());
        for (size_t idAO = nocc; idAO < ElectronicStr->NAO(); ++idAO) {
          occMOAB0_alpha[idAO] = 0.0;
        }
      }
      auxvec2 = ElectronicStr->getQAtoms();
      FockABexcpol_alpha = ElectronicStr->getQShells();
      FockABexcpol_beta = ElectronicStr->getmuAtoms();
      DensABint_beta = ElectronicStr->getthetaAtoms();
    }
  }
  void getQMSuperMoleculeHF() {
    //function running the full QM calculation on the supermolecule
    ElectronicStr->setMolecule(SuperMolecule);
    ElectronicStr->Calculate(0);
    //collect the relevant quantities
    SuperMoleculeNelectrons = SuperMolecule.Nelectrons();
    iauxvec = ElectronicStr->CoreCharges();
    for (size_t idAtm = 0; idAtm < SuperMoleculeNAtoms; ++idAtm) {
      SuperMoleculeNelectrons -= iauxvec[idAtm];
    }
    SuperMoleculeAOs = ElectronicStr->getAtomicOrbitals();
    SAOint = ElectronicStr->getSAO();
    FockABint_alpha = ElectronicStr->getFock();
    HcoreABint = ElectronicStr->getHcore();
    ElectronicStr->fetchDens(DensABint_alpha);
    ElectronicStr->calcExchangeFock(FockABexcint_alpha,FockABexcint_beta,SAOint);
    if (openshell) {
      FockABint_beta = ElectronicStr->getFockb();
      ElectronicStr->fetchbDens(DensABint_beta);
    }
    EnuclearABint = ElectronicStr->EnuclearMat();
  }
  void getQMSuperMoleculeDFT() {
    //function running the full QM calculation on the supermolecule
    auxMat1 = SuperMolecule.Geometry();
    ElectronicStr->setMolecule(SuperMolecule);
    ElectronicStr->Calculate(0);
    //collect the relevant quantities
    SuperMoleculeNelectrons = SuperMolecule.Nelectrons();
    iauxvec = ElectronicStr->CoreCharges();
    for (size_t idAtm = 0; idAtm < SuperMoleculeNAtoms; ++idAtm) {
      SuperMoleculeNelectrons -= iauxvec[idAtm];
    }
    SuperMoleculeAOs = ElectronicStr->getAtomicOrbitals();
    SAOint = ElectronicStr->getSAO();
    ElectronicStr->EnergyDispersion(auxMat1,auxMat2,1600.0);
    Edisp += auxMat2;
    if (usesolvation) {
      auxMat8 = ElectronicStr->CurrentSolvationEnergyMatrix(auxMat1);
      Esolv += auxMat8;
    }
    ElectronicStr->Enuclear(auxMat1,EnuclearABint,1600.0);
    ElectronicStr->EnergyMatrix(auxMat1,FockABint_alpha,false);
    HcoreABint = ElectronicStr->getHcore();
    ElectronicStr->fetchDens(DensABint_alpha);
    auxvec = ElectronicStr->getQAtoms();
    auxMat7 = ElectronicStr->getQShells();
    auxMat8 = ElectronicStr->getmuAtoms();
    auxMat9 = ElectronicStr->getthetaAtoms();
  }
  void getElectrostaticContributionsHF() {
    //function to determine the electrostatic contributions to the interaction energy
    ElectronicStr->setDens(DensAB0_alpha);
    if (openshell) {ElectronicStr->setbDens(DensAB0_beta);}
    ElectronicStr->calcFock(SAO0);
    FockABes_alpha = ElectronicStr->getFock();
    if (openshell) {FockABes_beta = ElectronicStr->getFockb();}
  }
  void getElectrostaticContributionsDFT() {
    //function to determine the electrostatic contributions to the interaction energy
    auxMat1 = SuperMolecule.Geometry();
    ElectronicStr->setDens(DensAB0_alpha);
    ElectronicStr->setQAtoms(auxvec2);
    ElectronicStr->setQShells(FockABexcpol_alpha);
    ElectronicStr->setmuAtoms(FockABexcpol_beta);
    ElectronicStr->setthetaAtoms(DensABint_beta);
    auxMat1 = SuperMolecule.Geometry();
    ElectronicStr->EnergyMatrix(auxMat1,FockABes_alpha,false,MonomerNAOs[0]);
  }
  void getPolarizationContributionsHF() {
    //function that relaxes the monomer parts of the density in SCF to get the polarization contribution to the interaction
    ElectronicStr->ConstrainedSCF(MonomerNAOs[0],MonomerNelectrons[0],MonomerNelectrons[1],false,usepopcontrain,0,200);
    ElectronicStr->fetchDens(DensABpol_alpha);
    if (openshell) {ElectronicStr->fetchbDens(DensABpol_beta);}
    FockABpol_alpha = ElectronicStr->getFock();
    ElectronicStr->calcExchangeFock(FockABexcpol_alpha,FockABexcpol_beta,SAO0);
    if (openshell) {FockABpol_beta = ElectronicStr->getFockb();}
  }
  void getPolarizationContributionsDFT() {
    //function that relaxes the monomer parts of the density in SCF to get the polarization contribution to the interaction
    auxMat1 = SuperMolecule.Geometry();
    ElectronicStr->setRestart(1);
    ElectronicStr->setQAtoms(auxvec);
    ElectronicStr->setQShells(auxMat7);
    ElectronicStr->setmuAtoms(auxMat8);
    ElectronicStr->setthetaAtoms(auxMat9);
    ElectronicStr->setDens(DensABint_alpha);
    ElectronicStr->ConstrainedSCF(MonomerNAOs[0],MonomerNelectrons[0],MonomerNelectrons[1],false,usepopcontrain,0,200);
    ElectronicStr->fetchDens(DensABpol_alpha);
    ElectronicStr->EnergyMatrix(auxMat1,FockABpol_alpha,false,MonomerNAOs[0]);
  }
  void getOverlapRepulsionHF() {
    //function that calculates the overlap-repulsion contributions
    ElectronicStr->setOcc(occMOAB0_alpha);
    ElectronicStr->setbOcc(occMOAB0_beta);
    ElectronicStr->setMOs(CMOAB0_alpha);
    ElectronicStr->setbMOs(CMOAB0_beta);
    ElectronicStr->setHcore(HcoreABint);
    auxMat1 = SuperMolecule.Geometry();
    ElectronicStr->CalcOverlap(PbarAB_alpha,PbarAB_beta,openshell,MonomerNAOs[0],true);
    ElectronicStr->calcFock(SAO0);
    FprimeAB_alpha = ElectronicStr->getFock();
    if (openshell) {FprimeAB_beta = ElectronicStr->getFockb();}
  }
  void getOverlapRepulsionDFT() {
    //function that calculates the overlap-repulsion contributions
    ElectronicStr->setOcc(occMOAB0_alpha);
    ElectronicStr->setMOs(CMOAB0_alpha);
    ElectronicStr->CalcOverlap(PbarAB_alpha,PbarAB_beta,openshell,MonomerNAOs[0],false);
    ElectronicStr->setDens(PbarAB_alpha);
    if (openshell) {ElectronicStr->setbDens(PbarAB_beta);}
    auxMat1 = SuperMolecule.Geometry();
    ElectronicStr->EnergyMatrix(auxMat1,FprimeAB_alpha,false,MonomerNAOs[0]);
  }
  void getInteractionEnergiesHF() {
    int countA = 0;
    int countB = 0;
    //function that assembles the atom-pair interaction energies for HF-like methods
    //for the electrostatic contribution
    auxMat2 = HcoreABint - HcoreAB0 + FockABes_alpha - FockAB0_alpha;
    //for the polarization contribution
    auxMat3 = DensABpol_alpha - DensAB0_alpha;
    //for the charge transfer contribution
    auxMat4 = DensABint_alpha - DensABpol_alpha;
    //for the overlap repulsion
    auxMat5 = PbarAB_alpha - DensAB0_alpha;
    for (size_t idAtm = 0; idAtm < SuperMoleculeNAtoms; ++idAtm) {
      countB = 0;
      for (size_t idBtm = 0; idBtm < SuperMoleculeNAtoms; ++idBtm) {
        for (size_t idAO = 0; idAO < SuperMoleculeAOs[idAtm]; ++idAO) {
          for (size_t idBO = 0; idBO < SuperMoleculeAOs[idBtm]; ++idBO) {
            Ees(idAtm + 1,idBtm + 1) += 0.5*auxMat2(countA + idAO + 1,countB + idBO + 1)*DensAB0_alpha(countB + idBO + 1,countA + idAO + 1);
            Epol(idAtm + 1,idBtm + 1) += 0.5*(HcoreABint(countA + idAO + 1,countB + idBO + 1)*auxMat3(countB + idBO + 1,countA + idAO + 1) + FockABpol_alpha(countA + idAO + 1,countB + idBO + 1)*DensABpol_alpha(countB + idBO + 1,countA + idAO + 1) - FockABes_alpha(countA + idAO + 1,countB + idBO + 1)*DensAB0_alpha(countB + idBO + 1,countA + idAO + 1));
            Ect(idAtm + 1,idBtm + 1) += 0.5*(HcoreABint(countA + idAO + 1,countB + idBO + 1)*auxMat4(countB + idBO + 1,countA + idAO + 1) + FockABint_alpha(countA + idAO + 1,countB + idBO + 1)*DensABint_alpha(countB + idBO + 1,countA + idAO + 1) - FockABpol_alpha(countA + idAO + 1,countB + idBO + 1)*DensABpol_alpha(countB + idBO + 1,countA + idAO + 1));
            Eexc_AB(idAtm + 1,idBtm + 1) += 0.5*(FockABexcint_alpha(countA + idAO + 1,countB + idBO + 1)*DensABint_alpha(countB + idBO + 1,countA + idAO + 1) - FockABexcpol_alpha(countA + idAO + 1,countB + idBO + 1)*DensABpol_alpha(countB + idBO + 1,countA + idAO + 1));
            Erep(idAtm + 1,idBtm + 1) -= 0.5*(HcoreABint(countA + idAO + 1,countB + idBO + 1)*auxMat5(countB + idBO + 1,countA + idAO + 1) + FprimeAB_alpha(countA + idAO + 1,countB + idBO + 1)*PbarAB_alpha(countB + idBO + 1,countA + idAO + 1) - FockABes_alpha(countA + idAO + 1,countB + idBO + 1)*DensAB0_alpha(countB + idBO + 1,countA + idAO + 1));
          }
        }
        countB += SuperMoleculeAOs[idBtm];
      }
      countA += SuperMoleculeAOs[idAtm];
    }
    if (openshell) {
      //rest of electrostatic
      auxMat2 = HcoreABint - HcoreAB0 + FockABes_beta - FockAB0_beta;
      //rest of polarization
      auxMat3 = DensABpol_beta - DensAB0_beta;
      //rest of charge transfer
      auxMat4 = DensABint_beta - DensABpol_beta;
      //for the overlap repulsion
      auxMat5 = PbarAB_beta - DensAB0_beta;
      countA = 0;
      countB = 0;
      for (size_t idAtm = 0; idAtm < SuperMoleculeNAtoms; ++idAtm) {
        countB = 0;
        for (size_t idBtm = 0; idBtm < SuperMoleculeNAtoms; ++idBtm) {
          for (size_t idAO = 0; idAO < SuperMoleculeAOs[idAtm]; ++idAO) {
            for (size_t idBO = 0; idBO < SuperMoleculeAOs[idBtm]; ++idBO) {
              Ees(idAtm + 1,idBtm + 1) += 0.5*auxMat2(countA + idAO + 1,countB + idBO + 1)*DensAB0_beta(countB + idBO + 1,countA + idAO + 1);
              Epol(idAtm + 1,idBtm + 1) += 0.5*(HcoreABint(countA + idAO + 1,countB + idBO + 1)*auxMat3(countB + idBO + 1,countA + idAO + 1) + FockABpol_beta(countA + idAO + 1,countB + idBO + 1)*DensABpol_beta(countB + idBO + 1,countA + idAO + 1) - FockABes_beta(countA + idAO + 1,countB + idBO + 1)*DensAB0_beta(countB + idBO + 1,countA + idAO + 1));
              Ect(idAtm + 1,idBtm + 1) += 0.5*(HcoreABint(countA + idAO + 1,countB + idBO + 1)*auxMat4(countB + idBO + 1,countA + idAO + 1) + FockABint_beta(countA + idAO + 1,countB + idBO + 1)*DensABint_beta(countB + idBO + 1,countA + idAO + 1) - FockABpol_beta(countA + idAO + 1,countB + idBO + 1)*DensABpol_beta(countB + idBO + 1,countA + idAO + 1));
              Eexc_AB(idAtm + 1,idBtm + 1) += 0.5*(FockABexcint_beta(countA + idAO + 1,countB + idBO + 1)*DensABint_beta(countB + idBO + 1,countA + idAO + 1) - FockABexcpol_beta(countA + idAO + 1,countB + idBO + 1)*DensABpol_beta(countB + idBO + 1,countA + idAO + 1));
              Erep(idAtm + 1,idBtm + 1) -= 0.5*(HcoreABint(countA + idAO + 1,countB + idBO + 1)*auxMat5(countB + idBO + 1,countA + idAO + 1) + FprimeAB_beta(countA + idAO + 1,countB + idBO + 1)*PbarAB_beta(countB + idBO + 1,countA + idAO + 1) - FockABes_beta(countA + idAO + 1,countB + idBO + 1)*DensAB0_beta(countB + idBO + 1,countA + idAO + 1));
            }
          }
          countB += SuperMoleculeAOs[idBtm];
        }
        countA += SuperMoleculeAOs[idAtm];
      }
    }
    //symmetrize
    for (size_t idAtm = 0; idAtm < SuperMoleculeNAtoms; ++idAtm) {
      for (size_t idBtm = 0; idBtm < idAtm; ++idBtm) {
        Ees(idAtm + 1,idBtm + 1) += Ees(idBtm + 1,idAtm + 1);
        Ees(idAtm + 1,idBtm + 1) += EnuclearABint(idAtm + 1,idBtm + 1) - EnuclearAB0(idAtm + 1,idBtm + 1);
        Epol(idAtm + 1,idBtm + 1) += Epol(idBtm + 1,idAtm + 1);
        Ect(idAtm + 1,idBtm + 1) += Ect(idBtm + 1,idAtm + 1);
        Eexc_AB(idAtm + 1,idBtm + 1) += Eexc_AB(idBtm + 1,idAtm + 1);
        Erep(idAtm + 1,idBtm + 1) += Erep(idBtm + 1,idAtm + 1);
        Ees(idBtm + 1,idAtm + 1) = Ees(idAtm + 1,idBtm + 1);
        Epol(idBtm + 1,idAtm + 1) = Epol(idAtm + 1,idBtm + 1);
        Ect(idBtm + 1,idAtm + 1) = Ect(idAtm + 1,idBtm + 1);
        Eexc_AB(idBtm + 1,idAtm + 1) = Eexc_AB(idAtm + 1,idBtm + 1);
        Erep(idBtm + 1,idAtm + 1) = Erep(idAtm + 1,idBtm + 1);
      }
    }
    //subtract intramolecular polarization-exchange from polarization
    Ees -= Erep;
    //subtract intermolecular exchange from charge transfer
    Ect -= Eexc_AB;
    (Epol + Ect).Print();
    //Edisp is absent here
    //Esolv is absent here
  }
  void getInteractionEnergiesDFT() {
    //function that assembles the atom-pair interaction energies for HF-like methods
    //electrostatic interactions
    Ees = FockABes_alpha - FockAB0_alpha;
    Epol = FockABpol_alpha - FockABes_alpha;
    Ect = FockABint_alpha - FockABpol_alpha;
    //Eexc_AB is absent here
    Erep = FprimeAB_alpha - FockABes_alpha;
    if (fastedda) {Ect -= Erep;}
    else {Epol -= Erep;}
    //Edisp is processed in getQMSuperMoleculeDFT
    //Esolv as well
  }
  void calcTotalDecomposedEnergies() {
    //function that sums each contribution to the interaction energy
    electrostatic = 0.0;
    polarization = 0.0;
    chargetransfer = 0.0;
    interexchange = 0.0;
    overlaprepulsion = 0.0;
    dispersion = 0.0;
    for (size_t idatm = 0; idatm < SuperMoleculeNAtoms; ++idatm) {
      for (size_t idbtm = 0; idbtm < idatm + 1; ++idbtm) {
        electrostatic += Ees(idatm + 1,idbtm + 1);
        polarization += Epol(idatm + 1,idbtm + 1);
        chargetransfer += Ect(idatm + 1,idbtm + 1);
        interexchange += Eexc_AB(idatm + 1,idbtm + 1);
        overlaprepulsion += Erep(idatm + 1,idbtm + 1);
        dispersion += Edisp(idatm + 1,idbtm + 1);
      }
    }
    if (usesolvation) {
      solvation = 0.0;
      for (size_t idatm = 0; idatm < SuperMoleculeNAtoms; ++idatm) {
        for (size_t idbtm = 0; idbtm < idatm + 1; ++idbtm) {
          solvation += Esolv(idatm + 1,idbtm + 1);
        }
      }
    }
  }
  void OutputEnergies() {
    std::cout << "electrostatic energies:" << std::endl;
    Ees.Print();
    std::cout << "polarization energies:" << std::endl;
    Epol.Print();
    std::cout << "charge transfer energies:" << std::endl;
    Ect.Print();
    std::cout << "intermolecular exchange:" << std::endl;
    Eexc_AB.Print();
    std::cout << "overlap repulsion:" << std::endl;
    Erep.Print();
    std::cout << "dispersion:" << std::endl;
    Edisp.Print();
    if (usesolvation) {
      std::cout << "solvation:" << std::endl;
      Esolv.Print();
    }
    std::cout << "Ees        << " << electrostatic*627.5096080305927 << std::endl;
    std::cout << "Epol       << " << polarization*627.5096080305927 << std::endl;
    std::cout << "Ect        << " << chargetransfer*627.5096080305927 << std::endl;
    std::cout << "Einterexc  << " << interexchange*627.5096080305927 << std::endl;
    std::cout << "Erep       << " << overlaprepulsion*627.5096080305927 << std::endl;
    std::cout << "Edisp      << " << dispersion*627.5096080305927 << std::endl;
    if (usesolvation) {std::cout << "Esolv      << " << solvation*627.5096080305927 << std::endl;}
  }
};

#endif //_Energy_Decomposition_Analysis_
