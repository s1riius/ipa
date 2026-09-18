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

#ifndef _MNDO_
#define _MNDO_
#include "QC.hpp"
#include "Intermolecular.hpp"
#include "atoms/AtomPackage.hpp"
#include "math/ConvergenceAccelaration.hpp"
#include "math/SolverPackage.hpp"
#include "math/FunctionPackage.hpp"
#include "parameters/EISOLpar.hpp"
#include "parameters/ZeroOvpar.hpp"
#include "basissets/2ElectronDewar.hpp"

//descrition:
//The MNDO method with sp basis

class MNDO: public QCbasis {
  //this is the implementation of Dewar's MNDO
  //M. J. S. Dewar, W. Thiel, J. Am. Chem. Soc., 99(15), 4899, 1977
  //M. J. S. Dewar, W. Thiel, Theor. Chim. Acta (Berl.), 46, 89, 1977
  //bX quantities are "barred" tensors, which are used only for the open-shell case
protected:
  matrixE gammaSE;
  matrixE SProt;
  matrixE PProt;
  matrixE SDrot;
  matrixE PDrot;
  matrixE DDrot;
  matrixE olddens;
  matrixE bolddens;
  matrixE Bij;
  matrixE AUXmat;
  matrixE integrals;
  matrixE iBlockdR;
  matrixE iBlockdt;
  matrixE iBlockdp;
  matrixE iBlockdRVC;
  matrixE iBlockdtVC;
  matrixE iBlockdpVC;
  matrixE iBlockdRVD;
  matrixE iBlockdtVD;
  matrixE iBlockdpVD;
  std::vector<std::vector<int> > chg1;
  std::vector<std::vector<int> > chg2;
  std::vector<double> AUXvec;
  std::vector<double> AUXvec2;
  std::vector<double> D;
  std::vector<double> Dd;
  std::vector<double> KA;
  std::vector<double> KB;
  std::vector<double> LA;
  std::vector<double> LB;
  std::vector<double> MA;
  std::vector<double> MB;
  std::vector<double> intn;
  std::vector<double> enuc_dR;
  std::vector<double> enuc_dR2;
  std::vector<double> chg;
  std::vector<double> p1;
  std::vector<double> bp1;
  std::vector<int> pos;
  std::vector<int> positions_pp;
  std::vector<size_t> nintegrals;
  std::vector<size_t> atomslocal;
  std::vector<size_t> ipos;
  std::vector<size_t> AOs_local;
  std::string corecorr;                 //this is a control variable to tell which empirical correction should be used for core-core terms
  size_t DIISstatus;
  bool doDIIS;
  bool doRCA;
  bool RCAstatus;
  bool RCAdone;
  bool doPDDG;
  bool oscillationdamp;
  int sizepos_pp;
  int print;
  int ncoreelectrons;
  double oscillationdampthresh;
  double natocc;
  double factorA;
  double factorB;
  double factor;
  double sumA;
  double sumB;
  double deltaA;
  double dval;
  double odval;
  double dampfactor;
public:
  MNDO(BSet _bset, Molecule _mol, const std::string & _openclosed = "0", std::string _corecorrection = "0"): QCbasis(_bset,_mol) {
    openclosed = _openclosed;
    corecorr = _corecorrection;
    pos.resize(7);
    D.resize(4);
    intn.resize(10);
    sizepos_pp = 3;
    positions_pp.resize(sizepos_pp);
    positions_pp[0] = 5;
    positions_pp[1] = 8;
    positions_pp[2] = 10;
    nintegrals.resize(2);
    atomslocal.resize(2);
    AOs_local.resize(2);
    ipos.resize(2);
    chg.resize(2);
    doPDDG = false;
    oscillationdamp = false;
    oscillationdampthresh = 0.05;
  }
  ~MNDO() {}
  //getters
  virtual std::string Type() {return "MNDO";}
  virtual std::string TypeOverlap() {return "Orthogonal";}
  double getHeatFormation() {return EnthForm(etot);}
  virtual bool AtomWithDOrbitals(size_t atomicnr) {
    //function that determines whether atom has d orbitals
    return false;
  }
  bool OscillationDamping() {return oscillationdamp;}
  double ThresholdOscillationDamping() {return oscillationdampthresh;}
  //setters
  void setOscillationDamping(bool newoscdamp) {oscillationdamp = newoscdamp;}
  void setThresholdOscillationDamping(double newthresh) {oscillationdampthresh = newthresh;}
  //Calculation
  void Calculate(int _print = 1, size_t maxiter = 200, bool _DIIS = true, bool _RCA = true, double d2threshold = 1600.0) {
    //initialize necessary quantities
    matrixE geometry = mol.Geometry();
    print = _print;
    DIISstatus = 0;
    doDIIS = _DIIS;
    doRCA = _RCA;
    RCAdone = false;
    RCAstatus = false;
    double energy = 0.0;
    double oldenergy;
    double denconv;
    double aux;
    //required for DIIS; the dimensions are variable, so it is not worth to save these globally for the object
    Bij.resize(1,1);
    Bij(1,1) = 0.0;
    std::vector<matrixE> vError;        //if not declared here, these must be resized, which is the same
    std::vector<matrixE> vFock;         //if not declared here, these must be resized, which is the same
    std::vector<matrixE> vbFock;        //if not declared here, these must be resized, which is the same
    std::vector<matrixE> vDens;         //if not declared here, these must be resized, which is the same
    std::vector<matrixE> vbDens;        //if not declared here, these must be resized, which is the same
    //check whether calculation is possible
    checkAtoms();
    //determine whether system is closed- or open-shell
    AOs = basis.AtomNAOs(atoms);
    OpenClosed();
    if (openclosed != "0") {
      //in this case we want to override what was previously determined by OpenClosed()
      if ((openclosed == "open")||(openclosed == "UHF")||(openclosed == "uhf")) {shell = "open";}
      else if ((openclosed == "closed")||(openclosed == "RHF")||(openclosed == "rhf")) {shell = "closed";}
      else {throw("ERROR: MNDOd.hpp: MNDOd: Calculate(): calculation type not possible to determine");}
    }
    //get occupation vector
    Occ();
    if (print > 0) {
      std::cout << shell << "-shell calculation" << std::endl;
      std::cout << "occupation vector " << occupancy.size() << "/" << NAOs << std::endl;
      for (size_t id = 0; id < occupancy.size(); ++id) {
        std::cout << occupancy[id] << " ";
      }
      std::cout << std::endl;
      if (shell == "open") {
        std::cout << "occupation vector (bar) " << boccupancy.size() << "/" << NAOs << std::endl;
        for (size_t id = 0; id < boccupancy.size(); ++id) {
          std::cout << boccupancy[id] << " ";
        }
        std::cout << std::endl;
      }
    }
    se_integrals(geometry);
    sao = basis.SAO();
    //get Hcore
    calcHcore(sao,geometry);
    //initiliaze Fock Matrix with the density initialization
    if (restart < 1) {
      dens.resize(NAOs,NAOs);
      dens.zero();
      olddens = dens;
      if (shell == "open") {
        bdens = dens;
        bolddens = bdens;
      }
      initDens(geometry);
      calcFock(sao);
      calcMOs();
      olddens = dens;
      Dens(dens,occupancy,CMO,DensAUX);
      if (shell == "open") {
        bolddens = bdens;
        Dens(bdens,boccupancy,bCMO,DensAUX);
      }
      if (restart == 0) {restart = 1;}
    }
    oldenergy = Energy();
    if (print > 0) {std::cout << "starting energy: " << oldenergy << std::endl;}
    vDens.push_back(dens);
    if (shell == "open") {vbDens.push_back(bdens);}
    //SCF
    for (size_t iter = 1; iter < maxiter; ++iter) {
      if (print  > 0) {std::cout << "iteration " << iter << std::endl;}
      //get new MOs
      calcFock(sao);
      ShiftFock(shell == "open",1.0e-6);
      if ((!doDIIS)||(DIISstatus == 0)) {
        if (iter > startpseudodiag) {
          PseudoDiagonalization(CMO,Fock,EMOs,occupancy.size());
          if (shell == "open") {PseudoDiagonalization(bCMO,bFock,bEMOs,boccupancy.size());}
        }
        else {calcMOs();}
      }
      //save old density
      olddens = dens;
      if (shell == "open") {bolddens = bdens;}
      //get new energy
      energy = Energy();
      if (print > 0) {
        std::cout << "CMO " << std::endl;
        CMO.Print();
        if (shell == "open") {
          std::cout << "CMO (bar) " << std::endl;
          bCMO.Print();
        }
        if (print > 1) {
          std::cout << "density " << std::endl;
          dens.Print();
          if (shell == "open") {
            std::cout << "density (bar) " << std::endl;
            bdens.Print();
          }
        }
        std::cout << "Energy = " << energy << std::endl;
      }
      //DIIS
      if (doDIIS) {
        //preparing DIIS
        vFock.push_back(Fock);
        if (shell == "open") {vbFock.push_back(bFock);}
        //DIIS step
        DIISse(vError,vFock,vbFock,dens,bdens,Bij,occupancy,boccupancy,DIISstatus,AUXmat,AUXvec,AUXvec2,DensAUX);
        if (DIISstatus == 2) {
          if (print > 0) {std::cout << "DIIS convergence after " << iter << " iterations" << std::endl;}
          break;
        }
        if ((RCAstatus)&&(DIISstatus == 0)) {
          //If DIIS still does not start converging, then RCA
          RCA(vFock,vbFock,vDens,vbDens,occupancy,boccupancy,energy - enuc,oldenergy - enuc,AUXmat,AUXvec,DensAUX,1.e-7);
          dens = vDens[1]; 
          vDens.clear();
          vDens.push_back(dens);
          if (shell == "open") {
            bdens = vbDens[1]; 
            vbDens.clear();
            vbDens.push_back(bdens);
          }
          RCAdone = true;
        }
      }
      //assume that in open-shell cases RCA is on or off simultaneously for both alpha and beta cases
      if (((!doDIIS)||(DIISstatus == 0))&&(!(RCAstatus))) {
        //get new density
        Dens(dens,occupancy,CMO,DensAUX);
        if (doDIIS) {vDens.push_back(dens);}
        if ((vDens.size() > 1)&&(doRCA)) {RCAstatus = true;}
        if (shell == "open") {
          Dens(bdens,boccupancy,bCMO,DensAUX);
          if (doDIIS) {vbDens.push_back(bdens);}
          if ((vbDens.size() > 1)&&(doRCA)&&(!(RCAstatus))) {RCAstatus = true;}
        }
        if (print > 0) {
          std::cout << "density" << std::endl;
          dens.Print();
          std::cout << "Energy = " << energy << std::endl;
          std::cout << "orbital energies: ";
          for (size_t id = 0; id < EMOs.size(); ++id) {
            std::cout << EMOs[id] << " ";
          }
          std::cout << std::endl;
          if (shell == "open") {
            std::cout << "density (bar)" << std::endl;
            bdens.Print();
            std::cout << "orbital energies (bar): ";
            for (size_t id = 0; id < bEMOs.size(); ++id) {
              std::cout << bEMOs[id] << " ";
            }
            std::cout << std::endl;
          }
        }
        //check for convergence on density
        if (shell == "open") {
          denconv = 0.0;
          for (size_t idc = 0; idc < NAOs; ++idc) {
            for (size_t idr = 0; idr < NAOs; ++idr) {
              aux = dens(idc + 1, idr + 1) + bdens(idc + 1, idr + 1) - olddens(idc + 1, idr + 1) - bolddens(idc + 1, idr + 1);
              denconv += aux*aux;
            }
          }
        }
        else {
          denconv = 0.0;
          for (size_t idc = 0; idc < NAOs; ++idc) {
            for (size_t idr = 0; idr < NAOs; ++idr) {
              aux = dens(idc + 1, idr + 1) - olddens(idc + 1, idr + 1);
              denconv += aux*aux;
            }
          }
        }
        if (denconv < thresh_dens*thresh_dens) {
          if (print > 0) {std::cout << "convergence in density after " << iter << " iterations" << std::endl;}
          break;
        }
      }
      if ((oscillationdamp)&&((iter + 1)%3 == 0)) {
        if (p1.size() != NAOs) {p1.resize(NAOs);}
        DampOscillations(iter,dens,olddens,p1);
        if (shell == "open") {
          if (bp1.size() != NAOs) {bp1.resize(NAOs);}
          DampOscillations(iter,bdens,bolddens,bp1);
        }
      }
      //check for convergence on energy
      if (fabs(oldenergy - energy) < thresh_en) {
        if (print > 0) {std::cout << "convergence in energy after " << iter << " iterations" << std::endl;}
        break;
      }
      oldenergy = energy;
      if ((RCAdone)&&(RCAstatus)) {
        RCAdone = false;
        RCAstatus = false;
      }
    }
    //almost all the techniques to avoid directly solving the SCF require calculating the MOs at the end;
    //only the "normal" way does not require this, but that should have ampered convergence as well, 
    //so one last matrix diagonalization does not affect the calculation times by much
    calcMOs();
    dens = olddens;
    if (shell == "open") {bdens = bolddens;}
    if (corecorr != "0"){
      if ((corecorr == "D3H4X")||(corecorr == "d3h4x")) {corecorrection = D3H4X(atoms,geometry,this->Type());}
      else if ((corecorr == "D3*H4X")||(corecorr == "d3*h4x")) {corecorrection = D3sH4X(atoms,geometry,this->Type(),true,true,false,1.5,5.5,9000.0);}
      else if ((corecorr == "D3RX")||(corecorr == "d3rx")) {corecorrection = D3RX(atoms,geometry,this->Type());}
      else if ((corecorr == "D3H+")||(corecorr == "d3h+")) {corecorrection = D3Hplus(atoms,geometry,this->Type(),true);}
      else if ((corecorr == "D3*H+")||(corecorr == "d3*h+")) {corecorrection = D3sHplus(atoms,geometry,this->Type(),true,true,false,9000.0,10.5,1.4);}
    }
    if (print > 0) {
      std::cout << "Optimized MOs" << std::endl;
      CMO.Print();
      std::cout << "Optimized orbital energies: ";
      for (size_t id = 0; id < EMOs.size(); ++id) {
        std::cout << EMOs[id] << " ";
      }
      std::cout << std::endl;
      if (shell == "open") {
        std::cout << "Optimized MOs (bar)" << std::endl;
        bCMO.Print();
        std::cout << "Optimized orbital energies (bar): ";
        for (size_t id = 0; id < bEMOs.size(); ++id) {
          std::cout << bEMOs[id] << " ";
        }
        std::cout << std::endl;
      }
      std::cout << "core correction: " << corecorrection << std::endl;
      std::cout << "Energy with correction = " << energy + corecorrection << std::endl;
      std::cout << "Energy = " << energy << std::endl;
      std::cout << "Nuclear Repulsion = " << enuc << std::endl;
      std::cout << "Electronic Energy = " << energy - enuc << std::endl;
      std::cout << "Enthalpy Formation = " << EnthForm(energy) << std::endl;
    }
    etot = energy;
  }
  void ConstrainedSCF(int endmonomerA, double NelecA, double NelecB, bool fromscratch, bool popconstrain, int _print = 1, size_t maxiter = 200, double thresholdVk = 1.0e-6) {
    //the constrained SCF procedure of Wu and van Voorhis
    //Q. Wu, T. van Voorhis, J. Chem. Theory Comput., 2, 765, 2006
    matrixE geometry = mol.Geometry();
    print = _print;
    double denconv;
    double aux;
    double Vk = 0.0;
    double dWdV = 0.0;
    double d2WdV2 = 0.0;
    double shift = 0.01;
    double dWdVk2;
    double mlimit = 0.0;
    double plimit = 0.0;
    double energy = 0.0;
    double oldenergy = 0.0;
    double dqcontrol = 0.0;
    double dx;
    double dxold;
    double du;
    double dl;
    bool found;
    matrixE auxmatA(NAOs,NAOs);
    matrixE auxmatB(NAOs,NAOs);
    matrixE wFock;
    matrixE wA;
    sao = Identity(NAOs);
    if (popconstrain) {BuildElectronConstrainPotential(sao,wA,endmonomerA,0.0);}
    matrixE wB = wA;
    if (popconstrain) {BuildElectronConstrainPotential(sao,wB,-endmonomerA,0.0);}
    if ((print > 0)&&(popconstrain)) {
      std::cout << "number of electrons: " << NelecA << " " << NelecB << std::endl;
      std::cout << "SAO" << std::endl;
      sao.Print();
      std::cout << "wA" << std::endl;
      wA.Print();
      std::cout << "wB" << std::endl;
      wB.Print();
    }
    //check whether calculation is possible
    checkAtoms();
    //determine whether system is closed- or open-shell
    AOs = basis.AtomNAOs(atoms);
    OpenClosed();
    if (openclosed != "0") {
      //in this case we want to override what was previously determined by OpenClosed()
      if ((openclosed == "open")||(openclosed == "UHF")||(openclosed == "uhf")) {shell = "open";}
      else if ((openclosed == "closed")||(openclosed == "RHF")||(openclosed == "rhf")) {shell = "closed";}
      else {throw("ERROR: MNDOd.hpp: MNDOd: Calculate(): calculation type not possible to determine");}
    }
    //get occupation vector
    Occ();
    if (print > 0) {
      std::cout << shell << "-shell calculation" << std::endl;
      std::cout << "occupation vector " << occupancy.size() << "/" << NAOs << std::endl;
      for (size_t id = 0; id < occupancy.size(); ++id) {
        std::cout << occupancy[id] << " ";
      }
      std::cout << std::endl;
      if (shell == "open") {
        std::cout << "occupation vector (bar) " << boccupancy.size() << "/" << NAOs << std::endl;
        for (size_t id = 0; id < boccupancy.size(); ++id) {
          std::cout << boccupancy[id] << " ";
        }
        std::cout << std::endl;
      }
    }
    se_integrals(geometry);
    sao = basis.SAO();
    //get Hcore
    calcHcore(sao,geometry);
    //initiliaze Fock Matrix with the density initialization
    if (restart < 1) {
      dens.resize(NAOs,NAOs);
      dens.zero();
      olddens = dens;
      if (shell == "open") {
        bdens = dens;
        bolddens = bdens;
      }
      initDens(geometry);
      calcFock(sao);
      calcMOs();
      olddens = dens;
      Dens(dens,occupancy,CMO,DensAUX);
      if (shell == "open") {
        bolddens = bdens;
        Dens(bdens,boccupancy,bCMO,DensAUX);
      }
      if (restart == 0) {restart = 1;}
    }
    if (!popconstrain) {
      CleanMatrix(Hcore,endmonomerA);
      std::cout << "new Hcore" << std::endl;
      Hcore.Print();
    }
    oldenergy = 0.0;
    if (print > 0) {std::cout << "starting energy: " << oldenergy << std::endl;}
    //SCF
    for (size_t iter = 1; iter < maxiter; ++iter) {
      if (print  > 0) {std::cout << "iteration " << iter << std::endl;}
      //get new MOs
      calcFock(sao);
      if (popconstrain) {
        wFock = Fock;
        Fock += (wA - wB)*Vk;
      }
      else {
        CleanMatrix(Fock,endmonomerA);
        if (shell == "open") {CleanMatrix(bFock,endmonomerA);}
      }
      calcMOs();
      olddens = dens;
      if (shell == "open") {bolddens = bdens;}
      //get new density
      Dens(dens,occupancy,CMO,DensAUX);
      if (shell == "open") {Dens(bdens,boccupancy,bCMO,DensAUX);}
      if (!popconstrain) {
        CleanMatrix(dens,endmonomerA);
        if (shell == "open") {CleanMatrix(bdens,endmonomerA);}
      }
      else {
        shift = 0.01;
        this->gradElectronConstraintPotential(dWdV,wA,wB,NelecA,NelecB,false);
        found = false;
        //bracket the minimum
        if (fabs(dWdV) > thresholdVk) {
redobracketing:
          found = BracketMinimum(wFock,auxmatA,wA,wB,Vk,dWdV,dWdVk2,mlimit,plimit,shift,NelecA,NelecB,1.0e-7,1.0e-10,1.0e-15,1.0e-8);
          if (!found) {
            shift *= 10.0;
            goto redobracketing;
          }
        }
        dx = fabs(plimit - mlimit);
        dxold = dx;
        Fock = wFock + (wA - wB)*Vk;
        calcMOs();
        Dens(dens,occupancy,CMO,DensAUX,true,1.0e-8);
        //now derivatives of the potential
        this->gradElectronConstraintPotential(dWdV,wA,wB,NelecA,NelecB,false);
        this->HessianElectronConstraintPotential(d2WdV2,wA,wB,auxmatA,auxmatB,false);
        found = ConstraintMicroIter(wFock,wA,wB,auxmatA,auxmatA,auxmatB,Vk,dWdV,d2WdV2,mlimit,plimit,du,dl,dx,dxold,NelecA,NelecB,thresholdVk,maxiter);
        if (!found) {ConstraintBissectionSearch(wFock,wA,wB,auxmatA,Vk,dWdV,dx,mlimit,plimit,NelecA,NelecB,thresholdVk,maxiter);}
        Fock = wFock;
      }
      if (print > 0) {
        std::cout << "density" << std::endl;
        dens.Print();
        std::cout << "Energy = " << energy << std::endl;
        std::cout << "orbital energies: ";
        for (size_t id = 0; id < EMOs.size(); ++id) {
          std::cout << EMOs[id] << " ";
        }
        std::cout << std::endl;
        if (shell == "open") {
          std::cout << "density (bar)" << std::endl;
          bdens.Print();
          std::cout << "orbital energies (bar): ";
          for (size_t id = 0; id < bEMOs.size(); ++id) {
            std::cout << bEMOs[id] << " ";
          }
          std::cout << std::endl;
        }
      }
      //get new energy
      energy = Energy();
      //check for convergence on density
      if (shell == "open") {
        denconv = 0.0;
        for (size_t idc = 0; idc < NAOs; ++idc) {
          for (size_t idr = 0; idr < NAOs; ++idr) {
            aux = dens(idc + 1, idr + 1) + bdens(idc + 1, idr + 1) - olddens(idc + 1, idr + 1) - bolddens(idc + 1, idr + 1);
            denconv += aux*aux;
          }
        }
      }
      else {
        denconv = 0.0;
        for (size_t idc = 0; idc < NAOs; ++idc) {
          for (size_t idr = 0; idr < NAOs; ++idr) {
            aux = dens(idc + 1, idr + 1) - olddens(idc + 1, idr + 1);
            denconv += aux*aux;
          }
        }
      }
      if (denconv < thresh_dens*thresh_dens) {
        if (print > 0) {std::cout << "convergence in density after " << iter << " iterations" << std::endl;}
        break;
      }
      //check for convergence on energy
      if (fabs(oldenergy - energy) < thresh_en) {
        if (print > 0) {std::cout << "convergence in energy after " << iter << " iterations" << std::endl;}
        break;
      }
      oldenergy = energy;
    }
    etot = energy;
  }
  void calcExchangeFock(matrixE & xFockA, matrixE & xFockB, matrixE & SAO) {
    //function that calculates the Fock matrix
    xFockA = Hcore;
    xFockA.zero();
    if (shell == "open") {xFockB = xFockA;}
    //atom position control for Fock and density
    int cntA = 0;
    int cntB = 0;
    //atom position control for gamma matrix
    int gposA = 0;
    int gposB = 0;
    int im;
    int in;
    int il;
    int is;
    bool transitionmetalA;
    bool transitionmetalB;
    for (size_t idA = 0; idA < Natoms; ++idA) {
      transitionmetalA = (((atoms[idA] > 20)&&(atoms[idA] < 30))||((atoms[idA] > 38)&&(atoms[idA] < 48))||((atoms[idA] > 56)&&(atoms[idA] < 80)));
      for (size_t imu = 0; imu < AOs[idA]; ++imu) {
        //get position in gamma matrix
        im = ShiftAO(imu,transitionmetalA);
        pos[0] = posgamma(im,im);
        for (size_t inu = 0; inu < AOs[idA]; ++inu) {
          in = ShiftAO(inu,transitionmetalA);
          pos[1] = posgamma(im,in);
          for (size_t ilambda = 0; ilambda < AOs[idA]; ++ ilambda) {
            il = ShiftAO(ilambda,transitionmetalA);
            pos[2] = posgamma(in,il);
            pos[3] = posgamma(im,il);
            xFockA(cntA + imu + 1,cntA + imu + 1) -= 0.5*dens(cntA + inu + 1,cntA + ilambda + 1)*gammaSE(gposA + pos[1] + 1,gposA + pos[3] + 1);
            if (shell == "open") {
              xFockA(cntA + imu + 1,cntA + imu + 1) -= 0.5*dens(cntA + inu + 1,cntA + ilambda + 1)*gammaSE(gposA + pos[1] + 1,gposA + pos[3] + 1);
              xFockB(cntA + imu + 1,cntA + imu + 1) -= bdens(cntA + inu + 1,cntA + ilambda + 1)*gammaSE(gposA + pos[1] + 1,gposA + pos[3] + 1);
            }
            if (imu != inu) {
              for (size_t isigma = 0; isigma < AOs[idA]; ++ isigma) {
                is = ShiftAO(isigma,transitionmetalA);
                pos[2] = posgamma(il,is);
                pos[4] = posgamma(is,in);
                xFockA(cntA + imu + 1,cntA + inu + 1) -= 0.5*dens(cntA + ilambda + 1,cntA + isigma + 1)*gammaSE(gposA + pos[3] + 1,gposA + pos[4] + 1);
                if (shell == "open") {
                  xFockA(cntA + imu + 1,cntA + inu + 1) -= 0.5*dens(cntA + ilambda + 1,cntA + isigma + 1)*gammaSE(gposA + pos[3] + 1,gposA + pos[4] + 1);
                  xFockB(cntA + imu + 1,cntA + inu + 1) -= bdens(cntA + ilambda + 1,cntA + isigma + 1)*gammaSE(gposA + pos[3] + 1,gposA + pos[4] + 1);
                }
              }
            }
          }
          cntB = 0;
          gposB = 0;
          for (size_t idB = 0; idB < idA; ++idB) {
            transitionmetalB = (((atoms[idB] > 20)&&(atoms[idB] < 30))||((atoms[idB] > 38)&&(atoms[idB] < 48))||((atoms[idB] > 56)&&(atoms[idB] < 80)));
            for (size_t ilambda = 0; ilambda < AOs[idB]; ++ilambda) {
              il = ShiftAO(ilambda,transitionmetalB);
              for (size_t isigma = 0; isigma < AOs[idB]; ++isigma) {
                is = ShiftAO(isigma,transitionmetalB);
                pos[2] = posgamma(il,is);
                xFockA(cntA + imu + 1,cntB + ilambda + 1) -= 0.5*dens(cntA + inu + 1,cntB + isigma + 1)*gammaSE(gposA + pos[1] + 1,gposB + pos[2] + 1);
                if (shell == "open") {
                  xFockA(cntA + imu + 1,cntB + ilambda + 1) -= 0.5*dens(cntA + inu + 1,cntB + isigma + 1)*gammaSE(gposA + pos[1] + 1,gposB + pos[2] + 1);
                  xFockB(cntA + imu + 1,cntB + ilambda + 1) -= bdens(cntA + inu + 1,cntB + isigma + 1)*gammaSE(gposA + pos[1] + 1,gposB + pos[2] + 1);
                }
                xFockA(cntB + ilambda + 1,cntA + imu + 1) -= 0.5*dens(cntB + isigma + 1,cntA + inu + 1)*gammaSE(gposB + pos[2] + 1,gposA + pos[1] + 1);
                if (shell == "open") {
                  xFockA(cntB + ilambda + 1,cntA + imu + 1) -= 0.5*dens(cntB + isigma + 1,cntA + inu + 1)*gammaSE(gposB + pos[2] + 1,gposA + pos[1] + 1);
                  xFockB(cntB + ilambda + 1,cntA + imu + 1) -= bdens(cntB + isigma + 1,cntA + inu + 1)*gammaSE(gposB + pos[2] + 1,gposA + pos[1] + 1);
                }
              }
            }
            cntB += AOs[idB];
            gposB += AOs[idB]*(AOs[idB] + 1)/2;
          }
        }
      }
      cntA += AOs[idA];
      gposA += AOs[idA]*(AOs[idA] + 1)/2;
    }
  }
  //calculate MO occupancy
  void Occ() {
    //function to calculate the occupancy vector for a system
    occupancy.clear();
    boccupancy.clear();
    size_t Nelec = size_t(mol.Nelectrons());
    size_t Nelec_aux = mol.Multiplicity() - 1;             //total number of unpaired electrons
    if (print > 0) {std::cout << "multiplicity " << mol.Multiplicity() << std::endl;}
    Nelec -= Nelec_aux + ncoreelectrons;
    Nelec /= 2;
    for (size_t idx = 0; idx < Nelec; ++idx) {
      occupancy.push_back(2.0);
    }
    for (size_t idx = 0; idx < Nelec_aux; ++idx) {
      occupancy.push_back(1.0);
    }
    if (shell == "open") {
      for (size_t idx = 0; idx < occupancy.size(); ++idx) {
        if (int(occupancy[idx]) == 2) {
          occupancy[idx] -= 1.0;
          boccupancy.push_back(1.0);
        }
      }
    }
  }
  //calculate energy
  double Energy() {
    double energy = enuc;
    //other terms
    iBlockdR = Hcore + Fock;
    energy += 0.5*trace(iBlockdR,dens);
    if (shell == "open") {
      iBlockdR = Hcore + bFock;
      energy += 0.5*trace(iBlockdR,bdens);
    }
    return energy;
  }
  void calcEnergy() {etot = this->Energy();}
  double EnthForm(double energy) {
    //function calculating heats of formation
    double enth = energy;
    for (size_t iatm = 0; iatm < Natoms; ++iatm) {
      enth += ElementFormExp(atoms[iatm]) - ElementFormTheo(atoms[iatm]);
    }
    return enth*au2kcalmol;
  }
  //functions required to build the Hamiltonian
  int ShiftAO(int iorb, bool transmetal) {
    //just a shift function to simplify code
    int im = iorb;
    if (transmetal) {
      if (iorb > 4) {im = iorb - 5;}
      else if (iorb < 5) {im = iorb + 4;}
    }
    return im;
  }
  void calcHcore(matrixE & sao, matrixE & geom) {
    //function that calculates the core Hamiltonian
    Hcore.resize(NAOs,NAOs);
    Hcore.zero();
    int counter = 0;
    int counter2 = 0;
    int counter3 = 0;
    int pgamma;
    int L = 0;
    int Lp = 0;
    int im;
    int in;
    bool transitionmetalA;
    bool transitionmetalB;
    double RAAmax;
    double Raux;
    double RAB;
    for (size_t idatm = 0; idatm < Natoms; ++idatm) {
      transitionmetalA = TransitionMetal(atoms[idatm]);
      L = 2*int(transitionmetalA);
      rAB[0] = geom(idatm + 1,1);
      rAB[1] = geom(idatm + 1,2);
      rAB[2] = geom(idatm + 1,3);
      RAAmax = ZeroOverlap(atoms[idatm]);
      for (size_t imu = 0; imu < AOs[idatm]; ++imu) {
        Hcore(counter + imu + 1, counter + imu + 1) = UlX(atoms[idatm],L);
        im = ShiftAO(imu,transitionmetalA);
        Hcore(counter + imu + 1, counter + imu + 1) += VAB(counter2 + posgamma(im,im) + 1,1);
        for (size_t inu = imu + 1; inu < AOs[idatm]; ++inu) {
          //determine the position in Vmunu required
          in = ShiftAO(inu,transitionmetalA);
          pgamma = posgamma(im,in);
          Hcore(counter + imu + 1, counter + inu + 1) += VAB(counter2 + pgamma + 1,1);
          Hcore(counter + inu + 1, counter + imu + 1) += VAB(counter2 + pgamma + 1,1);
        }
        counter3 = 0;
        for (size_t idbtm = 0; idbtm < idatm; ++idbtm) {
          RAB = 0.0;
          for (size_t idcoord = 0; idcoord < 3; ++idcoord){
            Raux = rAB[idcoord] - geom(idbtm + 1,idcoord + 1);
            RAB += Raux*Raux;
          }
          Raux = fmax(RAAmax,ZeroOverlap(atoms[idbtm]));         //this here is an approximation, but works conservatively well
          if (RAB < Raux*Raux) {
            transitionmetalB = TransitionMetal(atoms[idbtm]);
            Lp = 2*int(transitionmetalB);
            for (size_t jorb = 0; jorb < AOs[idbtm]; ++jorb) {
              Hcore(counter + imu + 1, counter3 + jorb + 1) = 0.5*(betaA0(atoms[idatm],L) + betaA0(atoms[idbtm],Lp))*sao(counter + imu + 1, counter3 + jorb + 1);
              Hcore(counter3 + jorb + 1,counter + imu + 1) = Hcore(counter + imu + 1, counter3 + jorb + 1);
              if (transitionmetalB) {
                if (jorb == 4) {Lp = 0;}
                if (jorb == 5) {++Lp;}
              }
              else {
                if ((jorb == 0)||(jorb == 3)) {++Lp;}
              }
            }
          }
          counter3 += AOs[idbtm];
        }
        if (transitionmetalA) {
          if (imu == 4) {L = 0;}
          if (imu == 5) {++L;}
        }
        else {
          if ((imu == 0)||(imu == 3)) {++L;}
        }
      }
      counter += AOs[idatm];
      counter2 += AOs[idatm]*(AOs[idatm] + 1)/2;
    }
    if (print > 0) {
      std::cout << "Hcore matrix" << std::endl;
      Hcore.Print();
    }
  }
  void HcoreCD_dX(matrixE & SAO_dX, size_t atmnrC, size_t atmnrD, size_t naoC, size_t naoD) {
    //function that calculates the off-diagonal block of the core Hamiltonian or any of its derivatives, involving the pair of atoms C and D
    //input is a matrix containing already the overlap matrix SAO (or its derivatives), which is modified to get the respective Hcore block
    //note that this applies only to off-diagonal blocks in terms of atoms, i.e. C != D
    bool transitionmetalC = TransitionMetal(atmnrC);
    bool transitionmetalD = TransitionMetal(atmnrD);
    int L = 2*int(transitionmetalC);
    int Lp = 0;
    for (size_t iorb = 0; iorb < naoC; ++iorb) {
      Lp = 2*int(transitionmetalD);
      for (size_t jorb = 0; jorb < naoD; ++jorb) {
        SAO_dX(iorb + 1,jorb + 1) *= 0.5*(betaA0(atmnrC,L) + betaA0(atmnrD,Lp));
        if (transitionmetalD) {
          if (jorb == 4) {Lp = 0;}
          if (jorb == 5) {++Lp;}
        }
        else {
          if ((jorb == 0)||(jorb == 3)) {++Lp;}
        }
      }
      if (transitionmetalC) {
        if (iorb == 4) {L = 0;}
        if (iorb == 5) {++L;}
      }
      else {
        if ((iorb == 0)||(iorb == 3)) {++L;}
      }
    }
  }
  void initDens(matrixE & geometry) {
    //function that calculates starting values for the density matrix according to what is done in MOPAC
    //only diagonal elements changed
    double dens_el;
    double bdens_el;
    double factorSym = 1.0;
    double factorShell = 1.0;
    bool transitionmetal;
    int cntA = 0;
    //generate asymmetry and open-shell conditions
    if (shell == "open") {
      factorShell = double(occupancy.size())/double(occupancy.size() + boccupancy.size());
      if (occupancy.size() != boccupancy.size()) {factorSym = 1.1;}
    }
    GoedeckerCharges(QAtoms,atoms,geometry,mol.Charge());
    for (size_t idA = 0; idA < Natoms; ++idA) {
      QAtoms[idA] *= -1.0;
      QAtoms[idA] += atoms[idA] - CoreCharge[idA];
      transitionmetal = (((atoms[idA] > 20)&&(atoms[idA] < 30))||((atoms[idA] > 38)&&(atoms[idA] < 48))||((atoms[idA] > 56)&&(atoms[idA] < 80)));
      for (size_t imu = 0; imu < AOs[idA]; ++imu) {
        if (!transitionmetal) {
          if (imu > 3) {break;}
          dens_el = factorShell*factorSym*QAtoms[idA]/4.0;
          bdens_el = (1.0 - factorShell)*QAtoms[idA]/(4.0*factorSym);
        }
        else {
          dens_el = factorShell*factorSym*QAtoms[idA]/9.0;
          bdens_el = (1.0 - factorShell)*QAtoms[idA]/(9.0*factorSym);
        }
        if (CoreCharge[idA] == 0) {
          dens_el *= 4.0;
          bdens_el *= 4.0;
        }
        dens(cntA + imu + 1,cntA + imu + 1) = dens_el;
        if (shell == "open") {bdens(cntA + imu + 1,cntA + imu + 1) = bdens_el;}
      }
      cntA += AOs[idA];
    }
  }
  void calcFock(matrixE & SAO) {
    //function that calculates the Fock matrix
    Fock = Hcore;
    if (shell == "open") {bFock = Hcore;}
    //atom position control for Fock and density
    int cntA = 0;
    int cntB = 0;
    //atom position control for gamma matrix
    int gposA = 0;
    int gposB = 0;
    int im;
    int in;
    int il;
    int is;
    bool transitionmetalA;
    bool transitionmetalB;
    for (size_t idA = 0; idA < Natoms; ++idA) {
      transitionmetalA = (((atoms[idA] > 20)&&(atoms[idA] < 30))||((atoms[idA] > 38)&&(atoms[idA] < 48))||((atoms[idA] > 56)&&(atoms[idA] < 80)));
      for (size_t imu = 0; imu < AOs[idA]; ++imu) {
        //get position in gamma matrix
        im = ShiftAO(imu,transitionmetalA);
        pos[0] = posgamma(im,im);
        for (size_t inu = 0; inu < AOs[idA]; ++inu) {
          in = ShiftAO(inu,transitionmetalA);
          pos[1] = posgamma(im,in);
          for (size_t ilambda = 0; ilambda < AOs[idA]; ++ ilambda) {
            il = ShiftAO(ilambda,transitionmetalA);
            pos[2] = posgamma(in,il);
            pos[3] = posgamma(im,il);
            Fock(cntA + imu + 1,cntA + imu + 1) += dens(cntA + inu + 1,cntA + ilambda + 1)*(gammaSE(gposA + pos[0] + 1,gposA + pos[2] + 1) - 0.5*gammaSE(gposA + pos[1] + 1,gposA + pos[3] + 1));
            if (shell == "open") {
              Fock(cntA + imu + 1,cntA + imu + 1) += bdens(cntA + inu + 1,cntA + ilambda + 1)*gammaSE(gposA + pos[0] + 1,gposA + pos[2] + 1) - 0.5*dens(cntA + inu + 1,cntA + ilambda + 1)*gammaSE(gposA + pos[1] + 1,gposA + pos[3] + 1);
              bFock(cntA + imu + 1,cntA + imu + 1) += bdens(cntA + inu + 1,cntA + ilambda + 1)*(gammaSE(gposA + pos[0] + 1,gposA + pos[2] + 1) - gammaSE(gposA + pos[1] + 1,gposA + pos[3] + 1));
              bFock(cntA + imu + 1,cntA + imu + 1) += dens(cntA + inu + 1,cntA + ilambda + 1)*gammaSE(gposA + pos[0] + 1,gposA + pos[2] + 1);
            }
            if (imu != inu) {
              for (size_t isigma = 0; isigma < AOs[idA]; ++ isigma) {
                is = ShiftAO(isigma,transitionmetalA);
                pos[2] = posgamma(il,is);
                pos[4] = posgamma(is,in);
                Fock(cntA + imu + 1,cntA + inu + 1) += dens(cntA + ilambda + 1,cntA + isigma + 1)*(gammaSE(gposA + pos[1] + 1,gposA + pos[2] + 1) - 0.5*gammaSE(gposA + pos[3] + 1,gposA + pos[4] + 1));
                if (shell == "open") {
                  Fock(cntA + imu + 1,cntA + inu + 1) += bdens(cntA + ilambda + 1,cntA + isigma + 1)*gammaSE(gposA + pos[1] + 1,gposA + pos[2] + 1);
                  Fock(cntA + imu + 1,cntA + inu + 1) -= 0.5*dens(cntA + ilambda + 1,cntA + isigma + 1)*gammaSE(gposA + pos[3] + 1,gposA + pos[4] + 1);
                  bFock(cntA + imu + 1,cntA + inu + 1) += (bdens(cntA + ilambda + 1,cntA + isigma + 1) + dens(cntA + ilambda + 1,cntA + isigma + 1))*gammaSE(gposA + pos[1] + 1,gposA + pos[2] + 1);
                  bFock(cntA + imu + 1,cntA + inu + 1) -= bdens(cntA + ilambda + 1,cntA + isigma + 1)*gammaSE(gposA + pos[3] + 1,gposA + pos[4] + 1);
                }
              }
            }
          }
          cntB = 0;
          gposB = 0;
          for (size_t idB = 0; idB < idA; ++idB) {
            transitionmetalB = (((atoms[idB] > 20)&&(atoms[idB] < 30))||((atoms[idB] > 38)&&(atoms[idB] < 48))||((atoms[idB] > 56)&&(atoms[idB] < 80)));
            for (size_t ilambda = 0; ilambda < AOs[idB]; ++ilambda) {
              il = ShiftAO(ilambda,transitionmetalB);
              for (size_t isigma = 0; isigma < AOs[idB]; ++isigma) {
                is = ShiftAO(isigma,transitionmetalB);
                pos[2] = posgamma(il,is);
                Fock(cntA + imu + 1,cntB + ilambda + 1) -= 0.5*dens(cntA + inu + 1,cntB + isigma + 1)*gammaSE(gposA + pos[1] + 1,gposB + pos[2] + 1);
                if (shell == "open") {
                  Fock(cntA + imu + 1,cntB + ilambda + 1) -= 0.5*dens(cntA + inu + 1,cntB + isigma + 1)*gammaSE(gposA + pos[1] + 1,gposB + pos[2] + 1);
                  bFock(cntA + imu + 1,cntB + ilambda + 1) -= bdens(cntA + inu + 1,cntB + isigma + 1)*gammaSE(gposA + pos[1] + 1,gposB + pos[2] + 1);
                }
                Fock(cntB + ilambda + 1,cntA + imu + 1) -= 0.5*dens(cntB + isigma + 1,cntA + inu + 1)*gammaSE(gposB + pos[2] + 1,gposA + pos[1] + 1);
                if (shell == "open") {
                  Fock(cntB + ilambda + 1,cntA + imu + 1) -= 0.5*dens(cntB + isigma + 1,cntA + inu + 1)*gammaSE(gposB + pos[2] + 1,gposA + pos[1] + 1);
                  bFock(cntB + ilambda + 1,cntA + imu + 1) -= bdens(cntB + isigma + 1,cntA + inu + 1)*gammaSE(gposB + pos[2] + 1,gposA + pos[1] + 1);
                }
                Fock(cntA + imu + 1,cntA + inu + 1) += dens(cntB + ilambda + 1,cntB + isigma + 1)*gammaSE(gposA + pos[1] + 1,gposB + pos[2] + 1);
                if (shell == "open") {
                  Fock(cntA + imu + 1,cntA + inu + 1) += bdens(cntB + ilambda + 1,cntB + isigma + 1)*gammaSE(gposA + pos[1] + 1,gposB + pos[2] + 1);
                  bFock(cntA + imu + 1,cntA + inu + 1) += (dens(cntB + ilambda + 1,cntB + isigma + 1) + bdens(cntB + ilambda + 1,cntB + isigma + 1))*gammaSE(gposA + pos[1] + 1,gposB + pos[2] + 1);
                }
                Fock(cntB + ilambda + 1,cntB + isigma + 1) += dens(cntA + imu + 1,cntA + inu + 1)*gammaSE(gposB + pos[2] + 1,gposA + pos[1] + 1);
                if (shell == "open") {
                  Fock(cntB + ilambda + 1,cntB + isigma + 1) += bdens(cntA + imu + 1,cntA + inu + 1)*gammaSE(gposB + pos[2] + 1,gposA + pos[1] + 1);
                  bFock(cntB + ilambda + 1,cntB + isigma + 1) += (dens(cntA + imu + 1,cntA + inu + 1) + bdens(cntA + imu + 1,cntA + inu + 1))*gammaSE(gposB + pos[2] + 1,gposA + pos[1] + 1);
                }
              }
            }
            cntB += AOs[idB];
            gposB += AOs[idB]*(AOs[idB] + 1)/2;
          }
        }
      }
      cntA += AOs[idA];
      gposA += AOs[idA]*(AOs[idA] + 1)/2;
    }
    if (print > 0) {
      std::cout << "Fock matrix" << std::endl;
      Fock.Print();
      if (shell == "open") {
        std::cout << "Fock matrix (bar)" << std::endl;
        bFock.Print();
      }
    }
  }
  void AnalyticalGrad(matrixE & gen, double tolerance = 1.0e-6) {
    //function calculating the gradients of the energy for MNDO methods
    //all final derivatives come for nuclear coordinates in atomic units!
    this->getDens(olddens);                //direct use of density matrix can lead to problems with open-shells
    if (shell == "open") {
      this->getbDens(bolddens);
      AUXmat = olddens - bolddens;
      olddens += bolddens;
      bolddens = AUXmat;
    }
    AUXmat = mol.Geometry();
    AOs = basis.AtomNAOs(atoms);
    Enuclear_dR(enuc_dR);
    int icntC = 0;
    int icntD = 0;
    int icnt = 0;
    int jcnt = 0;
    int kcnt = 0;
    int naoC;
    int naoD;
    int nintC;
    int nintD;
    int mu;
    int nu;
    int lambda;
    int sigma;
    bool transitionmetalC;
    bool transitionmetalD;
    gen.resize(3*Natoms,1);
    gen.zero();
    double enel_dR = 0.0;                            //derivatives of energy with respect to {R}_{CD}
    double enel_dt = 0.0;                            //derivatives of energy with respect to {\theta}_{CD}
    double enel_dp = 0.0;                            //derivatives of energy with respect to {\phi}_{CD}
    double RCD;
    double cost;
    double sint;
    double cosp;
    double sinp;
    double auxD;
    double auxV;
    double dRCDdxD;
    double dRCDdyD;
    double dRCDdzD;
    double dTCDdxD;
    double dTCDdyD;
    double dTCDdzD;
    double dPCDdxD;
    double dPCDdyD;
    double dPCDdzD;
    double chgC;
    double chgD;
    for (size_t idC = 0; idC < Natoms; ++idC) {
      chgC = double(int(atoms[idC]) - int(CoreCharge[idC]));
      naoC = AOs[idC];
      nintC = naoC*(naoC + 1)/2;
      icntD = icntC + naoC;
      transitionmetalC = TransitionMetal(atoms[idC]);
      for (size_t idD = idC + 1; idD < Natoms; ++idD) {
        chgD = double(int(atoms[idD]) - int(CoreCharge[idD]));
        naoD = AOs[idD];
        nintD = naoD*(naoD + 1)/2;
        transitionmetalD = TransitionMetal(atoms[idD]);
        //getting the orientation vector
        rAB[0] = (AUXmat(idD + 1,1) - AUXmat(idC + 1,1))*dist_Angstrom2aum1;                  //Delta x
        rAB[1] = (AUXmat(idD + 1,2) - AUXmat(idC + 1,2))*dist_Angstrom2aum1;                  //Delta y
        rAB[2] = (AUXmat(idD + 1,3) - AUXmat(idC + 1,3))*dist_Angstrom2aum1;                  //Delta z
        //normalizing it
        RCD = sqrt(rAB[0]*rAB[0] + rAB[1]*rAB[1] + rAB[2]*rAB[2]);
        //getting chain rule terms
        dRCDdxD = rAB[0]/RCD;
        dRCDdyD = rAB[1]/RCD;
        dRCDdzD = rAB[2]/RCD;
        dTCDdxD = 0.0;
        dTCDdyD = 0.0;
        dTCDdzD = -sqrt(rAB[0]*rAB[0] + rAB[1]*rAB[1])/(RCD*RCD);
        dPCDdxD = 0.0;
        dPCDdyD = 0.0;
        dPCDdzD = 0.0;
        if (fabs(rAB[0] + rAB[1]) > tolerance) {
          dTCDdyD = rAB[1]*rAB[2]/(sqrt(rAB[0]*rAB[0] + rAB[1]*rAB[1])*RCD*RCD);
          dTCDdxD = rAB[0]*rAB[2]/(sqrt(rAB[0]*rAB[0] + rAB[1]*rAB[1])*RCD*RCD);
          dPCDdyD = rAB[0]/(rAB[0]*rAB[0] + rAB[1]*rAB[1]);
          dPCDdxD = -rAB[1]/(rAB[0]*rAB[0] + rAB[1]*rAB[1]);
        }
        rAB[0] /= RCD;
        rAB[1] /= RCD;
        rAB[2] /= RCD;
        //getting trigonometric functions for rotations
        cost = rAB[2];
        sint = sqrt(1.0 - cost*cost);
        cosp = 1.0;
        sinp = 0.0;
        if (fabs(sint) > tolerance) {
          cosp = rAB[0]/sint;
          sinp = rAB[1]/sint;
        }
        //integrals
        iBlockdR.resize(nintC,nintD);              //resize the blocks containing the integral first-derivatives with respect to internuclear distance
        iBlockdt.resize(nintC,nintD);              //resize the blocks containing the integral first-derivatives with respect to angle theta
        iBlockdp.resize(nintC,nintD);              //resize the blocks containing the integral first-derivatives with respect to angle phi
        iBlockdRVC.resize(nintC,1);                //resize the blocks containing the integral first-derivatives with respect to internuclear distance
        iBlockdtVC.resize(nintC,1);                //resize the blocks containing the integral first-derivatives with respect to angle theta
        iBlockdpVC.resize(nintC,1);                //resize the blocks containing the integral first-derivatives with respect to angle phi
        iBlockdRVD.resize(nintD,1);                //resize the blocks containing the integral first-derivatives with respect to internuclear distance
        iBlockdtVD.resize(nintD,1);                //resize the blocks containing the integral first-derivatives with respect to angle theta
        iBlockdpVD.resize(nintD,1);                //resize the blocks containing the integral first-derivatives with respect to angle phi
        IntegralBlock2C_dR(0,iBlockdR,atoms[idC],atoms[idD],RCD,cost,sint,cosp,sinp);
        IntegralBlock2C_dA(0,iBlockdt,2,atoms[idC],atoms[idD],RCD,cost,sint,cosp,sinp);
        IntegralBlock2C_dA(0,iBlockdp,3,atoms[idC],atoms[idD],RCD,cost,sint,cosp,sinp);
        IntegralBlock2C_dR(1,iBlockdRVC,atoms[idC],atoms[idD],RCD,cost,sint,cosp,sinp);
        IntegralBlock2C_dA(1,iBlockdtVC,2,atoms[idC],atoms[idD],RCD,cost,sint,cosp,sinp);
        IntegralBlock2C_dA(1,iBlockdpVC,3,atoms[idC],atoms[idD],RCD,cost,sint,cosp,sinp);
        IntegralBlock2C_dR(1,iBlockdRVD,atoms[idD],atoms[idC],RCD,-cost,sint,-cosp,-sinp);
        IntegralBlock2C_dA(1,iBlockdtVD,2,atoms[idD],atoms[idC],RCD,-cost,sint,-cosp,-sinp);
        IntegralBlock2C_dA(1,iBlockdpVD,3,atoms[idD],atoms[idC],RCD,-cost,sint,-cosp,-sinp);
        enel_dR = 0.0;
        enel_dt = 0.0;
        enel_dp = 0.0;
        for (size_t idmu = 0; idmu < naoC; ++idmu) {
          mu = ShiftAO(idmu,transitionmetalC);
          for (size_t idnu = idmu; idnu < naoC; ++idnu) {
            nu = ShiftAO(idnu,transitionmetalC);
            jcnt = posgamma(mu,nu);
            auxV = olddens(icntC + idmu + 1,icntC + idnu + 1);
            if (idmu != idnu) {auxV *= 2.0;}
            enel_dR -= auxV*chgD*iBlockdRVC(jcnt + 1,1);
            enel_dt -= auxV*chgD*iBlockdtVC(jcnt + 1,1);
            enel_dp -= auxV*chgD*iBlockdpVC(jcnt + 1,1);
            for (size_t idlambda = 0; idlambda < naoD; ++idlambda) {
              lambda = ShiftAO(idlambda,transitionmetalD);
              for (size_t idsigma = idlambda; idsigma < naoD; ++idsigma) {
                sigma = ShiftAO(idsigma,transitionmetalD);
                kcnt = posgamma(lambda,sigma);
                if ((mu == 0)&&(nu == 0)) {
                  auxV = olddens(icntD + idlambda + 1,icntD + idsigma + 1);
                  if (idlambda != idsigma) {auxV *= 2.0;}
                  enel_dR -= auxV*chgC*iBlockdRVD(kcnt + 1,1);
                  enel_dt += auxV*chgC*iBlockdtVD(kcnt + 1,1);
                  enel_dp -= auxV*chgC*iBlockdpVD(kcnt + 1,1);
                }
                auxD = olddens(icntC + idmu + 1,icntC + idnu + 1)*olddens(icntD + idlambda + 1,icntD + idsigma + 1) - 0.5*olddens(icntC + idmu + 1,icntD + idlambda + 1)*olddens(icntC + idnu + 1,icntD + idsigma + 1);
                if (shell == "open") {auxD -= 0.5*bolddens(icntC + idmu + 1,icntD + idlambda + 1)*bolddens(icntC + idnu + 1,icntD + idsigma + 1);}
                if (idlambda != idsigma) {
                  auxD += olddens(icntC + idmu + 1,icntC + idnu + 1)*olddens(icntD + idsigma + 1,icntD + idlambda + 1) - 0.5*olddens(icntC + idmu + 1,icntD + idsigma + 1)*olddens(icntC + idnu + 1,icntD + idlambda + 1);
                  if (shell == "open") {auxD -= 0.5*bolddens(icntC + idmu + 1,icntD + idsigma + 1)*bolddens(icntC + idnu + 1,icntD + idlambda + 1);}
                }
                if (idmu != idnu) {
                  auxD += olddens(icntC + idnu + 1,icntC + idmu + 1)*olddens(icntD + idlambda + 1,icntD + idsigma + 1) - 0.5*olddens(icntC + idnu + 1,icntD + idlambda + 1)*olddens(icntC + idmu + 1,icntD + idsigma + 1);
                  if (shell == "open") {auxD -= 0.5*bolddens(icntC + idnu + 1,icntD + idlambda + 1)*bolddens(icntC + idmu + 1,icntD + idsigma + 1);}
                }
                if ((idlambda != idsigma)&&(idmu != idnu)) {
                  auxD += olddens(icntC + idnu + 1,icntC + idmu + 1)*olddens(icntD + idsigma + 1,icntD + idlambda + 1) - 0.5*olddens(icntC + idnu + 1,icntD + idsigma + 1)*olddens(icntC + idmu + 1,icntD + idlambda + 1);
                  if (shell == "open") {auxD -= 0.5*bolddens(icntC + idnu + 1,icntD + idsigma + 1)*bolddens(icntC + idmu + 1,icntD + idlambda + 1);}
                }
                enel_dR += auxD*iBlockdR(jcnt + 1,kcnt + 1);
                enel_dt += auxD*iBlockdt(jcnt + 1,kcnt + 1);
                enel_dp += auxD*iBlockdp(jcnt + 1,kcnt + 1);
              }
            }
          }
        }
        //Overlap terms
        basis.Overlap_dXCD(iBlockdR,"R",atoms[idC],atoms[idD],RCD,rAB);
        basis.Overlap_dXCD(iBlockdt,"t",atoms[idC],atoms[idD],RCD,rAB);
        basis.Overlap_dXCD(iBlockdp,"p",atoms[idC],atoms[idD],RCD,rAB);
        //transform to Hcore
        HcoreCD_dX(iBlockdR,atoms[idC],atoms[idD],naoC,naoD);
        HcoreCD_dX(iBlockdt,atoms[idC],atoms[idD],naoC,naoD);
        HcoreCD_dX(iBlockdp,atoms[idC],atoms[idD],naoC,naoD);
        //even though I could solve this with trace, I merge into a single loop
        for (size_t idmu = 0; idmu < naoC; ++idmu) {
          for (size_t idlambda = 0; idlambda < naoD; ++idlambda) {
            enel_dR += 2.0*olddens(icntC + idmu + 1,icntD + idlambda + 1)*iBlockdR(idmu + 1,idlambda + 1);
            enel_dt += 2.0*olddens(icntC + idmu + 1,icntD + idlambda + 1)*iBlockdt(idmu + 1,idlambda + 1);
            enel_dp += 2.0*olddens(icntC + idmu + 1,icntD + idlambda + 1)*iBlockdp(idmu + 1,idlambda + 1);
          }
        }
        gen(3*idC + 1,1) -= dist_Angstrom2aum1*((enel_dR + enuc_dR[icnt])*dRCDdxD + enel_dt*dTCDdxD + enel_dp*dPCDdxD);        //xC
        gen(3*idC + 2,1) -= dist_Angstrom2aum1*((enel_dR + enuc_dR[icnt])*dRCDdyD + enel_dt*dTCDdyD + enel_dp*dPCDdyD);        //yC
        gen(3*idC + 3,1) -= dist_Angstrom2aum1*((enel_dR + enuc_dR[icnt])*dRCDdzD + enel_dt*dTCDdzD + enel_dp*dPCDdzD);        //zC
        gen(3*idD + 1,1) += dist_Angstrom2aum1*((enel_dR + enuc_dR[icnt])*dRCDdxD + enel_dt*dTCDdxD + enel_dp*dPCDdxD);        //xD
        gen(3*idD + 2,1) += dist_Angstrom2aum1*((enel_dR + enuc_dR[icnt])*dRCDdyD + enel_dt*dTCDdyD + enel_dp*dPCDdyD);        //yD
        gen(3*idD + 3,1) += dist_Angstrom2aum1*((enel_dR + enuc_dR[icnt])*dRCDdzD + enel_dt*dTCDdzD + enel_dp*dPCDdzD);        //zD
        ++icnt;
        icntD += naoD;
      }
      icntC += naoC;
    }
    //add the empirical corrections
    if (corecorr != "0"){
      std::vector<double> dispersioncorrection;
      if ((corecorr == "D3H4X")||(corecorr == "d3h4x")) {gD3H4X(dispersioncorrection,atoms,AUXmat,this->Type());}
      else if ((corecorr == "D3*H4X")||(corecorr == "d3*h4x")) {gD3sH4X(dispersioncorrection,atoms,AUXmat,this->Type(),true,true,false,1.5,5.5,9000.0,1.0e-7);}
      else if ((corecorr == "D3RX")||(corecorr == "d3rx")) {gD3RX(dispersioncorrection,atoms,AUXmat,this->Type());}
      else if ((corecorr == "D3H+")||(corecorr == "d3h+")) {gD3Hplus(dispersioncorrection,atoms,AUXmat,this->Type(),true,9000.0,10.5,1.4);}
      else if ((corecorr == "D3*H+")||(corecorr == "d3*h+")) {gD3sHplus(dispersioncorrection,atoms,AUXmat,this->Type(),true,true,false,9000.0,10.5,1.4);}
      for (size_t idx = 0; idx < 3*Natoms; ++idx) {
        gen(idx + 1,1) += dispersioncorrection[idx];
      }
    }
  }
  int AvailableGradients() {return 1;}
  size_t derindex(size_t ic1, size_t ic2) {
    //function returning the index of the derivatives for a pair of internal coordinates ic1,ic2
    size_t index = 0;
    if (ic1 == 0) {
      if (ic2 == 0) {index = 3;}
      else if (ic2 == 1) {index = 6;}
      else if (ic2 == 2) {index = 7;}
    }
    else if (ic1 == 1) {
      if (ic2 == 0) {index = 9;}
      else if (ic2 == 1) {index = 4;}
      else if (ic2 == 2) {index = 8;}
    }
    else if (ic1 == 2) {
      if (ic2 == 0) {index = 10;}
      else if (ic2 == 1) {index = 11;}
      else if (ic2 == 2) {index = 5;}
    }
    return index;
  }
  //THIS IS BUGGY, DO NOT USE
  //void AnalyticalHess(matrixE & hessian, int algorithm = 2, double tolerance = 1.0e-7, int maxiter = 50) {
  //  //function calculating the hessian of the energy
  //  //all final derivatives come for nuclear coordinates in atomic units!
  //  //algorithm defines the algorithm used: 1 -> Head-Zerner V1 approximation, 2 -> full analytical calculation
  //  matrixE geometry = mol.Geometry();
  //  matrixE Amatrix_aa(1,1);                       //the full matrix that couples Uai vector elements in the CPHF equations to get derivatives of the density matrix: alpha-alpha matrix or closed-shell
  //  matrixE Amatrix_bb(1,1);                       //the full matrix that couples Uai vector elements in the CPHF equations to get derivatives of the density matrix: beta-beta matrix
  //  matrixE Amatrix_ab(1,1);                       //the full matrix that couples Uai vector elements in the CPHF equations to get derivatives of the density matrix: alpha-beta matrix
  //  matrixE Amatrix_ba(1,1);                       //the full matrix that couples Uai vector elements in the CPHF equations to get derivatives of the density matrix: beta-alpha matrix
  //  std::vector<matrixE> block_sao;                //derivatives of the overlap matrix, calculated for off-diagonal blocks between two (different) atoms
  //  std::vector<matrixE> iBlock;                   //derivatives of two-electron-two-center integrals, calculated in blocks of two atoms
  //  std::vector<matrixE> VBlockD;                  //derivatives of nuclear attraction integrals with sDsD as ket, calculated in blocks of two atoms
  //  std::vector<matrixE> VBlockC;                  //derivatives of nuclear attraction integrals with sCsC as ket, calculated in blocks of two atoms
  //  std::vector<matrixE> dPdX;                     //derivatives of the density matrix
  //  std::vector<matrixE> dQdX;                     //derivatives of the spin-density matrix
  //  this->getDens(olddens);                        //direct use of density matrix can lead to problems with open-shells
  //  if (shell == "open") {
  //    this->getbDens(bolddens);
  //    AUXmat = olddens - bolddens;
  //    olddens += bolddens;
  //    bolddens = AUXmat;
  //  }
  //  for (size_t idx = 0; idx < 9; ++idx) {
  //    block_sao.push_back(olddens);
  //    iBlock.push_back(olddens);
  //    VBlockC.push_back(olddens);
  //    VBlockD.push_back(olddens);
  //    if (idx < 3) {
  //      dPdX.push_back(olddens);
  //      dQdX.push_back(olddens);
  //    }
  //  }
  //  Enuclear_dR(enuc_dR);
  //  Enuclear_dR2(enuc_dR2);
  //  std::vector<double> dRCD(3,0.0);                  //first derivatives of {R}_{CD} with respect to nuclear coordinates
  //  std::vector<double> dtCD(3,0.0);                  //first derivatives of {\theta}_{CD} with respect to nuclear coordinates
  //  std::vector<double> dpCD(3,0.0);                  //first derivatives of {\phi}_{CD} with respect to nuclear coordinates
  //  std::vector<double> d2RCD(6,0.0);                 //second derivatives of {R}_{CD} with respect to nuclear coordinates
  //  std::vector<double> d2tCD(6,0.0);                 //second derivatives of theta_CD with respect to nuclear coordinates
  //  std::vector<double> d2pCD(6,0.0);                 //second derivatives of phi_CD with respect to nuclear coordinates
  //  std::vector<double> enel_dX(12,0.0);              //all the required derivatives: d{R}_{CD}, d{\theta}_{CD}, d{\phi}_{CD}, d2{R}_{CD}, d2{\theta}_{CD}, d2{\phi}_{CD}, d{R}_{CD}d{\theta}_{CD}, d{R}_{CD}d{\phi}_{CD}, d{\theta}_{CD}d{\phi}_{CD}
  //  std::vector<double> aux(4,0.0);
  //  std::vector<double> Qtilde_R(1,0.0);
  //  std::vector<double> bQtilde_R(1,0.0);
  //  std::vector<double> Qtilde_t(1,0.0);
  //  std::vector<double> bQtilde_t(1,0.0);
  //  std::vector<double> Qtilde_p(1,0.0);
  //  std::vector<double> bQtilde_p(1,0.0);
  //  std::vector<double> Uaa;                          //only for open-shell systems
  //  std::vector<double> Ubb;                          //only for open-shell systems
  //  std::vector<double> CROSSab;                      //only for open-shell systems
  //  std::vector<double> CROSSba;                      //only for open-shell systems
  //  std::vector<double> HZa;                          //Vector with Head-Zerner coefficients
  //  std::vector<double> HZb;                          //Vector with Head-Zerner coefficients
  //  AOs = basis.AtomNAOs(atoms);
  //  std::vector<size_t> twoatom(2,0);
  //  std::vector<size_t> twoatoms(2,0);
  //  std::vector<size_t> twoatomAOs(2,0);
  //  int naosC;
  //  int naosD;
  //  int nintC;
  //  int nintD;
  //  int spos;
  //  int icnt = 0;
  //  int icntC = 0;
  //  int icntD;
  //  int jcnt;
  //  int kcnt;
  //  int maxsze;
  //  int aA;
  //  int iA;
  //  int bA;
  //  int jA;
  //  int aB;
  //  int iB;
  //  int bB;
  //  int jB;
  //  int mu;
  //  int nu;
  //  int lambda;
  //  int sigma;
  //  bool transitionmetalC;
  //  bool transitionmetalD;
  //  int noccA = occupancy.size();          //number of occupied orbitals of the alpha component
  //  int nvirA = NAOs - noccA;              //number of virtual orbitals of the alpha component
  //  int noccB = boccupancy.size();         //number of occupied orbitals of the beta component
  //  int nvirB = NAOs - noccB;              //number of virtual orbitals of the beta component
  //  int corr = 1;                          //just auxiliary to avoid overflow problem
  //  if (noccB == 0) {
  //    noccB = 1;
  //    corr = 0;
  //  }
  //  double RCD;
  //  double cost;
  //  double sint;
  //  double cosp;
  //  double sinp;
  //  double chgC;
  //  double chgD;
  //  double unitconversion = dist_Angstrom2aum1*dist_Angstrom2aum1;
  //  iA = noccA*nvirA;                       //matrix size for alpha electrons
  //  iB = corr*noccB*nvirB;                  //matrix size for beta electrons
  //  maxsze = std::max(iA,iB);
  //  matrixE megamat(2*maxsze,2*maxsze);
  //  std::vector<double> megavec(2*maxsze,0.0);
  //  hessian.resize(3*Natoms,3*Natoms);
  //  hessian.zero();
  //  Qtilde_R.resize(maxsze);
  //  Qtilde_t.resize(maxsze);
  //  Qtilde_p.resize(maxsze);
  //  HZa.resize(maxsze);
  //  if (algorithm == 2) {Amatrix_aa.resize(maxsze,maxsze);}
  //  if (shell == "open") {
  //    HZb.resize(maxsze);
  //    bQtilde_R.resize(maxsze);
  //    bQtilde_t.resize(maxsze);
  //    bQtilde_p.resize(maxsze);
  //    if (algorithm == 2) {
  //      Uaa.resize(maxsze);
  //      Ubb.resize(maxsze);
  //      CROSSab.resize(maxsze);
  //      CROSSba.resize(maxsze);
  //      Amatrix_bb.resize(maxsze,maxsze);
  //      Amatrix_ab.resize(maxsze,maxsze);
  //    }
  //  }
  //  Amatrix_aa.zero();
  //  Amatrix_bb.zero();
  //  Amatrix_ab.zero();
  //  icntC = 0;
  //  for (size_t idC = 0; idC < Natoms; ++idC) {
  //    naosC = AOs[idC];
  //    nintC = naosC*(naosC + 1)/2;
  //    transitionmetalC = TransitionMetal(atoms[idC]);
  //    icntD = 0;
  //    for (size_t idD = 0; idD < Natoms; ++idD) {
  //      naosD = AOs[idD];
  //      nintD = naosD*(naosD + 1)/2;
  //      iBlock[0].resize(nintC,nintD);
  //      transitionmetalD = TransitionMetal(atoms[idD]);
  //      if (idC == idD) {IntegralBlock1C(iBlock[0],atoms[idC]);}
  //      else {
  //        //getting the orientation vector
  //        rAB[0] = (geometry(idD + 1,1) - geometry(idC + 1,1))*dist_Angstrom2aum1;                  //Delta x
  //        rAB[1] = (geometry(idD + 1,2) - geometry(idC + 1,2))*dist_Angstrom2aum1;                  //Delta y
  //        rAB[2] = (geometry(idD + 1,3) - geometry(idC + 1,3))*dist_Angstrom2aum1;                  //Delta z
  //        //normalizing it
  //        RCD = sqrt(rAB[0]*rAB[0] + rAB[1]*rAB[1] + rAB[2]*rAB[2]);
  //        rAB[0] /= RCD;
  //        rAB[1] /= RCD;
  //        rAB[2] /= RCD;
  //        //getting trigonometric functions for rotations
  //        cost = rAB[2];
  //        sint = sqrt(1.0 - cost*cost);
  //        cosp = 1.0;
  //        sinp = 0.0;
  //        if (fabs(sint) > tolerance) {
  //          cosp = rAB[0]/sint;
  //          sinp = rAB[1]/sint;
  //        }
  //        IntegralBlock2C(0,iBlock[0],atoms[idC],atoms[idD],RCD,cost,sint,cosp,sinp);
  //      }
  //      for (int ai = 0; ai < maxsze; ++ai) {
  //        aA = int(ai/noccA);
  //        iA = int(ai%noccA);
  //        aB = int(corr*ai/noccB);
  //        iB = int(corr*ai%noccB);
  //        for (int bj = ai; bj < maxsze; ++bj) {
  //          if ((bj > ai)&&(algorithm == 1)) {break;}
  //          bA = int(bj/noccA);
  //          jA = int(bj%noccA);
  //          bB = int(corr*bj/noccB);
  //          jB = int(corr*bj%noccB);
  //          for (size_t idmu = 0; idmu < naosC; ++idmu) {
  //            mu = ShiftAO(idmu,transitionmetalC);
  //            for (size_t idnu = idmu; idnu < naosC; ++idnu) {
  //              nu = ShiftAO(idnu,transitionmetalC);
  //              jcnt = posgamma(mu,nu);
  //              sint = 1.0;
  //              if (idmu == idnu) {sint = 0.5;}       //this is just a correction factor to ensure no double counting takes place
  //              d2RCD[0] = 0.0;
  //              d2RCD[1] = 0.0;
  //              d2RCD[3] = 0.0;
  //              if (ai < noccA*nvirA) {
  //                d2RCD[0] = sint*(CMO(icntC + idmu + 1,noccA + aA + 1)*CMO(icntC + idnu + 1,iA + 1) + CMO(icntC + idmu + 1,iA + 1)*CMO(icntC + idnu + 1,noccA + aA + 1));
  //                if (bj < noccA*nvirA) {d2RCD[1] = 0.5*sint*(CMO(icntC + idmu + 1,noccA + aA + 1)*CMO(icntC + idnu + 1,noccA + bA + 1) + CMO(icntC + idmu + 1,noccA + bA + 1)*CMO(icntC + idnu + 1,noccA + aA + 1));}
  //                else {d2RCD[1] = 0.0;}
  //                d2RCD[3] = 0.5*sint*(CMO(icntC + idmu + 1,noccA + aA + 1)*CMO(icntC + idnu + 1,jA + 1) + CMO(icntC + idmu + 1,jA + 1)*CMO(icntC + idnu + 1,noccA + aA + 1));
  //              }
  //              if (shell == "open") {
  //                d2RCD[1] *= 2.0;
  //                d2RCD[3] *= 2.0;
  //                aux[3] = 0.0;
  //                d2tCD[0] = 0.0;
  //                dRCD[2] = 0.0;
  //                aux[2] = 0.0;
  //                if (ai < noccB*nvirB) {
  //                  aux[3] = sint*(bCMO(icntC + idmu + 1,noccB + aB + 1)*bCMO(icntC + idnu + 1,iB + 1) + bCMO(icntC + idmu + 1,iB + 1)*bCMO(icntC + idnu + 1,noccB + aB + 1));
  //                  d2tCD[0] = sint*(bCMO(icntC + idmu + 1,noccB + aB + 1)*bCMO(icntC + idnu + 1,jB + 1) + bCMO(icntC + idmu + 1,jB + 1)*bCMO(icntC + idnu + 1,noccB + aB + 1));
  //                  if (bj < noccB*nvirB) {dRCD[2] = sint*(bCMO(icntC + idmu + 1,noccB + aB + 1)*bCMO(icntC + idnu + 1,noccB + bB + 1) + bCMO(icntC + idmu + 1,noccB + bB + 1)*bCMO(icntC + idnu + 1,noccB + aB + 1));}
  //                  else {dRCD[2] = 0.0;}
  //                  aux[2] = 2.0*sint*(bCMO(icntC + idmu + 1,noccB + aB + 1)*bCMO(icntC + idnu + 1,iB + 1) + bCMO(icntC + idmu + 1,iB + 1)*bCMO(icntC + idnu + 1,noccB + aB + 1));
  //                }
  //              }
  //              for (size_t idlambda = 0; idlambda < naosD; ++idlambda) {
  //                lambda = ShiftAO(idlambda,transitionmetalD);
  //                for (size_t idsigma = idlambda; idsigma < naosD; ++idsigma) {
  //                  sigma = ShiftAO(idsigma,transitionmetalD);
  //                  kcnt = posgamma(lambda,sigma);
  //                  cost = 1.0;
  //                  if (idsigma == idlambda) {cost = 0.5;}       //this is just a correction factor to ensure no double counting takes place
  //                  d2RCD[5] = 0.0;
  //                  dRCD[0] = 0.0;
  //                  dRCD[1] = 0.0;
  //                  d2RCD[4] = 0.0;
  //                  if (bj < noccA*nvirA) {
  //                    d2RCD[5] = cost*(CMO(icntD + idlambda + 1,noccA + bA + 1)*CMO(icntD + idsigma + 1,jA + 1) + CMO(icntD + idlambda + 1,jA + 1)*CMO(icntD + idsigma + 1,noccA + bA + 1));
  //                    dRCD[0] = cost*(CMO(icntD + idlambda + 1,iA + 1)*CMO(icntD + idsigma + 1,jA + 1) + CMO(icntD + idlambda + 1,jA + 1)*CMO(icntD + idsigma + 1,iA + 1));
  //                    dRCD[1] = cost*(CMO(icntD + idlambda + 1,noccA + bA + 1)*CMO(icntD + idsigma + 1,iA + 1) + CMO(icntD + idlambda + 1,iA + 1)*CMO(icntD + idsigma + 1,noccA + bA + 1));
  //                    d2RCD[4] = cost*(CMO(icntD + idlambda + 1,noccA + bA + 1)*CMO(icntD + idsigma + 1,jA + 1) + CMO(icntD + idlambda + 1,jA + 1)*CMO(icntD + idsigma + 1,noccA + bA + 1));
  //                  }
  //                  if ((ai < noccA*nvirA)&&(ai == bj)) {HZa[ai] += iBlock[0](jcnt + 1,kcnt + 1)*(d2RCD[0]*d2RCD[5] - d2RCD[1]*dRCD[0]);}
  //                  if (algorithm == 2) {Amatrix_aa(ai + 1,bj + 1) += iBlock[0](jcnt + 1,kcnt + 1)*(2.0*d2RCD[0]*d2RCD[5] - d2RCD[3]*dRCD[1] - d2RCD[1]*dRCD[0]);}
  //                  if (shell == "open") {
  //                    d2tCD[1] = 0.0;
  //                    d2tCD[3] = 0.0;
  //                    d2tCD[2] = 0.0;
  //                    if (bj < noccB*nvirB) {
  //                      d2tCD[1] = cost*(bCMO(icntD + idlambda + 1,noccB + bB + 1)*bCMO(icntD + idsigma + 1,jB + 1) + bCMO(icntD + idlambda + 1,jB + 1)*bCMO(icntD + idsigma + 1,noccB + bB + 1));
  //                      d2tCD[3] = cost*(bCMO(icntD + idlambda + 1,noccB + bB + 1)*bCMO(icntD + idsigma + 1,iB + 1) + bCMO(icntD + idlambda + 1,iB + 1)*bCMO(icntD + idsigma + 1,noccB + bB + 1));
  //                      d2tCD[2] = cost*(bCMO(icntD + idlambda + 1,iB + 1)*bCMO(icntD + idsigma + 1,jB + 1) + bCMO(icntD + idlambda + 1,jB + 1)*bCMO(icntD + idsigma + 1,iB + 1));
  //                    }
  //                    if ((ai < noccB*nvirB)&&(ai == bj)) {HZb[ai] += iBlock[0](jcnt + 1,kcnt + 1)*(aux[3]*d2tCD[1] - dRCD[2]*d2tCD[2]);}
  //                    if (algorithm == 2) {
  //                      Amatrix_bb(ai + 1,bj + 1) += iBlock[0](jcnt + 1,kcnt + 1)*(2.0*aux[3]*d2tCD[1] - d2tCD[0]*d2tCD[3] - dRCD[2]*d2tCD[2]);
  //                      Amatrix_ab(ai + 1,bj + 1) += 2.0*iBlock[0](jcnt + 1,kcnt + 1)*d2RCD[0]*d2tCD[1];
  //                      if (ai != bj) {Amatrix_ab(bj + 1,ai + 1) += iBlock[0](jcnt + 1,kcnt + 1)*aux[2]*d2RCD[4];}
  //                    }
  //                  }
  //                }
  //              }
  //            }
  //          }
  //        }
  //      }
  //      icntD += naosD;
  //    }
  //    icntC += naosC;
  //  }
  //  if (algorithm == 2) {
  //    //symmetrization and division by (ei - ea) factor
  //    for (size_t idr = 0; idr < maxsze; ++idr) {
  //      for (size_t idc = idr + 1; idc < maxsze; ++idc) {
  //        Amatrix_aa(idc + 1, idr + 1) = Amatrix_aa(idr + 1, idc + 1);
  //        if (shell == "open") {Amatrix_bb(idc + 1, idr + 1) = Amatrix_bb(idr + 1, idc + 1);}
  //      }
  //    }
  //    cosp = 2.0;                         //number of electrons per occupied orbital
  //    if (shell == "open") {
  //      Amatrix_ba = Amatrix_ab.trans();
  //      cosp = 1.0;
  //    }
  //    for (size_t idr = 0; idr < maxsze; ++idr) {
  //      cost = 0.0;
  //      sint = 0.0;
  //      if (idr < noccA*nvirA) {cost = 1.0/(EMOs[int(idr%noccA)] - EMOs[noccA + int(idr/noccA)]);}
  //      if ((shell == "open")&&(idr < noccB*nvirB)) {sint = 1.0/(bEMOs[int(corr*idr%noccB)] - bEMOs[noccB + int(corr*idr/noccB)]);}
  //      for (size_t idc = 0; idc < maxsze; ++idc) {
  //        Amatrix_aa(idr + 1,idc + 1) *= -cosp*cost;
  //        if (shell == "open") {
  //          Amatrix_bb(idr + 1,idc + 1) *= -sint;
  //          Amatrix_ab(idr + 1,idc + 1) *= cost;
  //          Amatrix_ba(idr + 1,idc + 1) *= sint;
  //        }
  //      }
  //      Amatrix_aa(idr + 1,idr + 1) += 1.0;
  //      if (shell == "open") {Amatrix_bb(idr + 1,idr + 1) += 1.0;}
  //    }
  //  }
  //  //main loop for Hessian
  //  icntC = 0;
  //  for (size_t idC = 0; idC < Natoms; ++idC) {
  //    chgC = double(int(atoms[idC]) - int(CoreCharge[idC]));
  //    naosC = AOs[idC];
  //    nintC = naosC*(naosC + 1)/2;
  //    icntD = icntC + naosC;
  //    transitionmetalC = TransitionMetal(atoms[idC]);
  //    twoatom[0] = icntC;
  //    twoatoms[1] = naosC;
  //    twoatomAOs[0] = naosC;
  //    for (size_t idD = idC + 1; idD < Natoms; ++idD) {
  //      transitionmetalD = TransitionMetal(atoms[idD]);
  //      //getting the orientation vector
  //      rAB[0] = (geometry(idD + 1,1) - geometry(idC + 1,1))*dist_Angstrom2aum1;                  //Delta x
  //      rAB[1] = (geometry(idD + 1,2) - geometry(idC + 1,2))*dist_Angstrom2aum1;                  //Delta y
  //      rAB[2] = (geometry(idD + 1,3) - geometry(idC + 1,3))*dist_Angstrom2aum1;                  //Delta z
  //      //normalizing it
  //      RCD = sqrt(rAB[0]*rAB[0] + rAB[1]*rAB[1] + rAB[2]*rAB[2]);
  //      aux[2] = RCD*RCD;
  //      aux[0] = rAB[0]*rAB[0] + rAB[1]*rAB[1];
  //      aux[1] = sqrt(aux[0]);
  //      //getting chain rule terms
  //      for (int idx = 0; idx < 3; ++idx) {
  //        dRCD[idx] = rAB[idx]/RCD;
  //        d2RCD[idx] = (RCD - rAB[idx]*dRCD[idx])/aux[2]; 
  //      }
  //      d2RCD[3] = -rAB[0]*dRCD[1]/aux[2];                                                                                                             //\frac{{d}^{2}{R}_{CD}}{dxdy}
  //      d2RCD[4] = -rAB[0]*dRCD[2]/aux[2];                                                                                                             //\frac{{d}^{2}{R}_{CD}}{dxdz}
  //      d2RCD[5] = -rAB[1]*dRCD[2]/aux[2];                                                                                                             //\frac{{d}^{2}{R}_{CD}}{dydz}
  //      for (size_t idx = 0; idx < dtCD.size(); ++idx) {
  //        d2pCD[idx] = 0.0;
  //        d2pCD[idx + 3] = 0.0;
  //        d2tCD[idx] = 0.0;
  //        d2tCD[idx + 3] = 0.0;
  //        dtCD[idx] = 0.0;
  //        dpCD[idx] = 0.0;
  //      }
  //      dtCD[2] = -sqrt(rAB[0]*rAB[0] + rAB[1]*rAB[1])/(RCD*RCD);                                                                                      //\frac{{d}{\theta}_{CD}}{dz}
  //      d2tCD[2] = 2.0*aux[1]*dRCD[2]/(aux[2]*RCD);                                                                                                    //\frac{{d}^{2}{\theta}_{CD}}{d{z}^{2}}
  //      if (fabs(rAB[0] + rAB[1]) > tolerance) {
  //        dtCD[1] = rAB[1]*rAB[2]/(sqrt(rAB[0]*rAB[0] + rAB[1]*rAB[1])*RCD*RCD);
  //        dtCD[0] = rAB[0]*rAB[2]/(sqrt(rAB[0]*rAB[0] + rAB[1]*rAB[1])*RCD*RCD);
  //        dpCD[1] = rAB[0]/(rAB[0]*rAB[0] + rAB[1]*rAB[1]);
  //        dpCD[0] = -rAB[1]/(rAB[0]*rAB[0] + rAB[1]*rAB[1]);
  //        d2tCD[0] = rAB[2]/(aux[2]*aux[1]) - 2.0*rAB[0]*rAB[2]*dRCD[0]/(aux[2]*aux[1]*RCD) - rAB[0]*rAB[0]*rAB[2]/(aux[2]*aux[0]*aux[1]);             //\frac{{d}^{2}{\theta}_{CD}}{d{x}^{2}}
  //        d2tCD[1] = rAB[2]/(aux[2]*aux[1]) - 2.0*rAB[1]*rAB[2]*dRCD[1]/(aux[2]*aux[1]*RCD) - rAB[1]*rAB[1]*rAB[2]/(aux[2]*aux[0]*aux[1]);             //\frac{{d}^{2}{\theta}_{CD}}{d{y}^{2}}
  //        d2tCD[3] = -2.0*rAB[2]*rAB[1]*dRCD[0]/(RCD*aux[2]*aux[1]) - rAB[2]*rAB[1]*rAB[0]/(aux[2]*aux[0]*aux[1]);                                     //\frac{{d}^{2}{\theta}_{CD}}{dxdy}
  //        d2tCD[4] = -rAB[0]/(aux[2]*aux[1]) + 2.0*aux[1]*dRCD[0]/(RCD*aux[2]);                                                                        //\frac{{d}^{2}{\theta}_{CD}}{dxdz}
  //        d2tCD[5] = -rAB[1]/(aux[2]*aux[1]) + 2.0*aux[1]*dRCD[1]/(RCD*aux[2]);                                                                        //\frac{{d}^{2}{\theta}_{CD}}{dydz}
  //        d2pCD[0] = 2.0*rAB[0]*rAB[1]/(aux[0]*aux[0]);                                                                                                //\frac{{d}^{2}{\phi}_{CD}}{d{x}^{2}}
  //        d2pCD[1] = -d2pCD[0];                                                                                                                        //\frac{{d}^{2}{\phi}_{CD}}{d{y}^{2}}
  //        d2pCD[3] = (rAB[1]*rAB[1] - rAB[0]*rAB[0])/(aux[0]*aux[0]);                                                                                  //\frac{{d}^{2}{\phi}_{CD}}{dxdy}
  //      }
  //      else {
  //        d2tCD[4] = -1.0/aux[2];
  //        d2tCD[5] = -1.0/aux[2];
  //      }
  //      rAB[0] /= RCD;
  //      rAB[1] /= RCD;
  //      rAB[2] /= RCD;
  //      //getting trigonometric functions for rotations
  //      cost = rAB[2];
  //      sint = sqrt(1.0 - cost*cost);
  //      cosp = 1.0;
  //      sinp = 0.0;
  //      if (fabs(sint) > tolerance) {
  //        cosp = rAB[0]/sint;
  //        sinp = rAB[1]/sint;
  //      }
  //      chgD = double(int(atoms[idD]) - int(CoreCharge[idD]));
  //      naosD = AOs[idD];
  //      nintD = naosD*(naosD + 1)/2;
  //      twoatom[1] = icntD;
  //      twoatomAOs[1] = naosD;
  //      for (size_t idx = 0; idx < 12; ++idx) {
  //        enel_dX[idx] = 0.0;
  //      }
  //      //Overlap terms
  //      basis.Overlap_dXCD(block_sao[0],"1",atoms[idC],atoms[idD],RCD,rAB);           //dS/dR
  //      basis.Overlap_dXCD(block_sao[1],"2",atoms[idC],atoms[idD],RCD,rAB);           //dS/dt
  //      basis.Overlap_dXCD(block_sao[2],"3",atoms[idC],atoms[idD],RCD,rAB);           //dS/dp
  //      basis.Overlap_dXCD(block_sao[3],"4",atoms[idC],atoms[idD],RCD,rAB);           //d2S/dR2
  //      basis.Overlap_dXCD(block_sao[4],"5",atoms[idC],atoms[idD],RCD,rAB);           //d2S/dt2
  //      basis.Overlap_dXCD(block_sao[5],"6",atoms[idC],atoms[idD],RCD,rAB);           //d2S/dp2
  //      basis.Overlap_dXCD(block_sao[6],"7",atoms[idC],atoms[idD],RCD,rAB);           //d2S/dRdt
  //      basis.Overlap_dXCD(block_sao[7],"8",atoms[idC],atoms[idD],RCD,rAB);           //d2S/dRdp
  //      basis.Overlap_dXCD(block_sao[8],"9",atoms[idC],atoms[idD],RCD,rAB);           //d2S/dtdp
  //      //even though I could solve this with trace, I merge into a single set of loops
  //      aA = 0;
  //      for (size_t idx = 0; idx < 12; ++idx) {                                                             //static terms
  //        if (idx == 9) {aA = 6;}
  //        if (idx < 9) {HcoreCD_dX(block_sao[aA],atoms[idC],atoms[idD],naosC,naosD);}                       //transform to Hcore
  //        for (size_t idmu = 0; idmu < naosC; ++idmu) {
  //          for (size_t idlambda = 0; idlambda < naosD; ++idlambda) {
  //            enel_dX[idx] += 2.0*olddens(icntC + idmu + 1,icntD + idlambda + 1)*block_sao[aA](idmu + 1,idlambda + 1);     //Note that Hcore/SAO and the density matrix are in the same "orbital space", meaning the orbitals are in the same order
  //          }
  //        }
  //        ++aA;
  //      }
  //      //getting derivatives of the Fock matrix so that CPHF may be solved
  //      dFockdX(1,iBlock[3],iBlock[6],olddens,bolddens,block_sao[0],RCD,cost,sint,cosp,sinp,idC,idD,atoms[idC],atoms[idD],icntC,icntD);
  //      dFockdX(2,iBlock[4],iBlock[7],olddens,bolddens,block_sao[1],RCD,cost,sint,cosp,sinp,idC,idD,atoms[idC],atoms[idD],icntC,icntD);
  //      dFockdX(3,iBlock[5],iBlock[8],olddens,bolddens,block_sao[2],RCD,cost,sint,cosp,sinp,idC,idD,atoms[idC],atoms[idD],icntC,icntD);
  //      //putting derivatives of Fock matrix into MO basis, these are in atomic units
  //      for (size_t ai = 0; ai < maxsze; ++ai) {
  //        if (ai < nvirA*noccA) {
  //          aA = int(ai/noccA);
  //          iA = int(ai%noccA);
  //          Qtilde_R[ai] = 0.0;
  //          Qtilde_t[ai] = 0.0;
  //          Qtilde_p[ai] = 0.0;
  //        }
  //        if ((shell == "open")&&(ai < nvirB*noccB)) {
  //          aB = int(corr*ai/noccB);
  //          iB = int(corr*ai%noccB);
  //          bQtilde_R[ai] = 0.0;
  //          bQtilde_t[ai] = 0.0;
  //          bQtilde_p[ai] = 0.0;
  //        }
  //        for (size_t idA = 0; idA < 2; ++idA) {
  //          for (size_t idB = 0; idB < 2; ++idB) {
  //            for (size_t idmu = 0; idmu < twoatomAOs[idA]; ++idmu) {
  //              for (size_t idnu = 0; idnu < twoatomAOs[idB]; ++idnu) {              //Fock and CMO matrices are in the same "orbital space", i.e., the orbitals are in the same order in both matrices
  //                if (ai < nvirA*noccA) {
  //                  Qtilde_R[ai] += CMO(twoatom[idA] + idmu + 1,noccA + aA + 1)*iBlock[3](twoatoms[idA] + idmu + 1,twoatoms[idB] + idnu + 1)*CMO(twoatom[idB] + idnu + 1,iA + 1);
  //                  Qtilde_t[ai] += CMO(twoatom[idA] + idmu + 1,noccA + aA + 1)*iBlock[4](twoatoms[idA] + idmu + 1,twoatoms[idB] + idnu + 1)*CMO(twoatom[idB] + idnu + 1,iA + 1);
  //                  Qtilde_p[ai] += CMO(twoatom[idA] + idmu + 1,noccA + aA + 1)*iBlock[5](twoatoms[idA] + idmu + 1,twoatoms[idB] + idnu + 1)*CMO(twoatom[idB] + idnu + 1,iA + 1);
  //                }
  //                if ((shell == "open")&&(ai < nvirB*noccB)) {
  //                  bQtilde_R[ai] += bCMO(twoatom[idA] + idmu + 1,noccB + aB + 1)*iBlock[6](twoatoms[idA] + idmu + 1,twoatoms[idB] + idnu + 1)*bCMO(twoatom[idB] + idnu + 1,iB + 1);
  //                  bQtilde_t[ai] += bCMO(twoatom[idA] + idmu + 1,noccB + aB + 1)*iBlock[7](twoatoms[idA] + idmu + 1,twoatoms[idB] + idnu + 1)*bCMO(twoatom[idB] + idnu + 1,iB + 1);
  //                  bQtilde_p[ai] += bCMO(twoatom[idA] + idmu + 1,noccB + aB + 1)*iBlock[8](twoatoms[idA] + idmu + 1,twoatoms[idB] + idnu + 1)*bCMO(twoatom[idB] + idnu + 1,iB + 1);
  //                }
  //              }
  //            }
  //          }
  //        }
  //        aux[0] = 0.0;
  //        if (algorithm == 1) {aux[0] = HZa[ai];}
  //        if (shell == "open") {
  //          if ((ai < nvirA*noccA)&&(algorithm == 2)) {Uaa[ai] = Qtilde_R[ai]/(EMOs[iA] - EMOs[noccA + aA] - HZa[ai]);}
  //          if (ai < nvirB*noccB) {
  //            aux[3] = 0.0;
  //            if (algorithm == 1) {aux[3] = HZb[ai];}
  //            else {Ubb[ai] = bQtilde_R[ai]/(bEMOs[iB] - bEMOs[noccB + aB] - HZb[ai]);}
  //            bQtilde_R[ai] /= (bEMOs[iB] - bEMOs[noccB + aB] - aux[3]);
  //            bQtilde_t[ai] /= (bEMOs[iB] - bEMOs[noccB + aB] - aux[3]);
  //            bQtilde_p[ai] /= (bEMOs[iB] - bEMOs[noccB + aB] - aux[3]);
  //          }
  //        }
  //        if (ai < nvirA*noccA) {
  //          Qtilde_R[ai] /= (EMOs[iA] - EMOs[noccA + aA] - aux[0]);
  //          Qtilde_t[ai] /= (EMOs[iA] - EMOs[noccA + aA] - aux[0]);
  //          Qtilde_p[ai] /= (EMOs[iA] - EMOs[noccA + aA] - aux[0]);
  //        }
  //      }
  //      if (algorithm == 2) {                     //At this point I have the A matrices and the Q vectors; solve (I - A)U = Q
  //        if (shell == "open") {
  //          for (size_t idr = 0; idr < maxsze; ++idr) {
  //            megavec[idr] = Qtilde_R[idr];
  //            megavec[maxsze + idr] = bQtilde_R[idr];
  //            for (size_t idc = 0; idc < maxsze; ++idc) {
  //              megamat(idr + 1,idc + 1) = Amatrix_aa(idr + 1,idc + 1);
  //              megamat(maxsze + idr + 1,idc + 1) = -Amatrix_ba(idr + 1,idc + 1);
  //              megamat(idr + 1,maxsze + idc + 1) = -Amatrix_ab(idr + 1,idc + 1);
  //              megamat(maxsze + idr + 1,maxsze + idc + 1) = Amatrix_bb(idr + 1,idc + 1);
  //            }
  //          }
  //          Solve_Ax_eq_b(megamat,megavec);
  //          for (size_t idr = 0; idr < maxsze; ++idr) {
  //            Qtilde_R[idr] = megavec[idr];
  //            bQtilde_R[idr] = megavec[maxsze + idr];
  //          }
  //          for (size_t idr = 0; idr < maxsze; ++idr) {
  //            megavec[idr] = Qtilde_t[idr];
  //            megavec[maxsze + idr] = bQtilde_t[idr];
  //          }
  //          Solve_Ax_eq_b(megamat,megavec);
  //          for (size_t idr = 0; idr < maxsze; ++idr) {
  //            Qtilde_t[idr] = megavec[idr];
  //            bQtilde_t[idr] = megavec[maxsze + idr];
  //          }
  //          for (size_t idr = 0; idr < maxsze; ++idr) {
  //            megavec[idr] = Qtilde_p[idr];
  //            megavec[maxsze + idr] = bQtilde_p[idr];
  //          }
  //          Solve_Ax_eq_b(megamat,megavec);
  //          for (size_t idr = 0; idr < maxsze; ++idr) {
  //            Qtilde_p[idr] = megavec[idr];
  //            bQtilde_p[idr] = megavec[maxsze + idr];
  //          }
  //        }
  //        else {
  //          Solve_Ax_eq_b(Amatrix_aa,Qtilde_R);             //derivative with respect to {R}_{CD}; result stored to vector
  //          Solve_Ax_eq_b(Amatrix_aa,Qtilde_t);             //derivative with respect to {\theta}_{CD}; result stored to vector
  //          Solve_Ax_eq_b(Amatrix_aa,Qtilde_p);             //derivative with respect to {\phi}_{CD}; result stored to vector
  //        }
  //      }
  //      //convert back to AO basis
  //      iBlock[0].resize(naosC + naosD,naosC + naosD);
  //      iBlock[1].resize(naosC + naosD,naosC + naosD);
  //      iBlock[2].resize(naosC + naosD,naosC + naosD);
  //      aux[0] = 2.0;
  //      if (shell == "open") {
  //        aux[0] = 1.0;
  //        iBlock[3].resize(naosC + naosD,naosC + naosD);
  //        iBlock[4].resize(naosC + naosD,naosC + naosD);
  //        iBlock[5].resize(naosC + naosD,naosC + naosD);
  //      }
  //      for (size_t idA = 0; idA < 2; ++idA) {
  //        for (size_t idB = 0; idB < 2; ++idB) {
  //          for (size_t idmu = 0; idmu < twoatomAOs[idA]; ++idmu) {
  //            for (size_t idnu = 0; idnu < twoatomAOs[idB]; ++idnu) {
  //              iBlock[0](twoatoms[idA] + idmu + 1,twoatoms[idB] + idnu + 1) = 0.0;
  //              iBlock[1](twoatoms[idA] + idmu + 1,twoatoms[idB] + idnu + 1) = 0.0;
  //              iBlock[2](twoatoms[idA] + idmu + 1,twoatoms[idB] + idnu + 1) = 0.0;
  //              if (shell == "open") {
  //                iBlock[3](twoatoms[idA] + idmu + 1,twoatoms[idB] + idnu + 1) = 0.0;
  //                iBlock[4](twoatoms[idA] + idmu + 1,twoatoms[idB] + idnu + 1) = 0.0;
  //                iBlock[5](twoatoms[idA] + idmu + 1,twoatoms[idB] + idnu + 1) = 0.0;
  //              }
  //              for (size_t ai = 0; ai < maxsze; ++ai) {
  //                if (ai < nvirA*noccA) {
  //                  aA = int(ai/noccA);
  //                  iA = int(ai%noccA);
  //                  iBlock[0](twoatoms[idA] + idmu + 1,twoatoms[idB] + idnu + 1) += aux[0]*Qtilde_R[ai]*(CMO(twoatom[idA] + idmu + 1,noccA + aA + 1)*CMO(twoatom[idB] + idnu + 1,iA + 1) + CMO(twoatom[idA] + idmu + 1,iA + 1)*CMO(twoatom[idB] + idnu + 1,noccA + aA + 1));
  //                  iBlock[1](twoatoms[idA] + idmu + 1,twoatoms[idB] + idnu + 1) += aux[0]*Qtilde_t[ai]*(CMO(twoatom[idA] + idmu + 1,noccA + aA + 1)*CMO(twoatom[idB] + idnu + 1,iA + 1) + CMO(twoatom[idA] + idmu + 1,iA + 1)*CMO(twoatom[idB] + idnu + 1,noccA + aA + 1));
  //                  iBlock[2](twoatoms[idA] + idmu + 1,twoatoms[idB] + idnu + 1) += aux[0]*Qtilde_p[ai]*(CMO(twoatom[idA] + idmu + 1,noccA + aA + 1)*CMO(twoatom[idB] + idnu + 1,iA + 1) + CMO(twoatom[idA] + idmu + 1,iA + 1)*CMO(twoatom[idB] + idnu + 1,noccA + aA + 1));
  //                }
  //                if ((shell == "open")&&(ai < nvirB*noccB)) {
  //                  aB = int(corr*ai/noccB);
  //                  iB = int(corr*ai%noccB);
  //                  iBlock[3](twoatoms[idA] + idmu + 1,twoatoms[idB] + idnu + 1) += bQtilde_R[ai]*(bCMO(twoatom[idA] + idmu + 1,noccB + aB + 1)*bCMO(twoatom[idB] + idnu + 1,iB + 1) + bCMO(twoatom[idA] + idmu + 1,iB + 1)*bCMO(twoatom[idB] + idnu + 1,noccB + aB + 1));
  //                  iBlock[4](twoatoms[idA] + idmu + 1,twoatoms[idB] + idnu + 1) += bQtilde_t[ai]*(bCMO(twoatom[idA] + idmu + 1,noccB + aB + 1)*bCMO(twoatom[idB] + idnu + 1,iB + 1) + bCMO(twoatom[idA] + idmu + 1,iB + 1)*bCMO(twoatom[idB] + idnu + 1,noccB + aB + 1));
  //                  iBlock[5](twoatoms[idA] + idmu + 1,twoatoms[idB] + idnu + 1) += bQtilde_p[ai]*(bCMO(twoatom[idA] + idmu + 1,noccB + aB + 1)*bCMO(twoatom[idB] + idnu + 1,iB + 1) + bCMO(twoatom[idA] + idmu + 1,iB + 1)*bCMO(twoatom[idB] + idnu + 1,noccB + aB + 1));
  //                }
  //              }
  //            }
  //          }
  //        }
  //      }
  //      dPdX[0] = iBlock[0];
  //      dPdX[1] = iBlock[1];
  //      dPdX[2] = iBlock[2];
  //      if (shell == "open") {             //if open-shell then we must build the derivatives of the P and Q matrices
  //        for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
  //          dPdX[idcoord] += iBlock[idcoord + 3];
  //          dQdX[idcoord] = iBlock[idcoord] - iBlock[idcoord + 3];
  //        }
  //      }
  //      for (size_t idmu = 0; idmu < naosC; ++idmu) {                                                     //response terms
  //        for (size_t idlambda = 0; idlambda < naosD; ++idlambda) {
  //          for (size_t idc1 = 0; idc1 < 3; ++idc1) {
  //            for (size_t idc2 = 0; idc2 < 3; ++idc2) {
  //              enel_dX[derindex(idc1,idc2)] += 2.0*dPdX[idc1](idmu + 1,naosC + idlambda + 1)*block_sao[idc2](idmu + 1,idlambda + 1);
  //            }
  //          }
  //        }
  //      }
  //      for (size_t idcoord = 0; idcoord < 9; ++idcoord) {
  //        iBlock[idcoord].resize(nintC,nintD);
  //        VBlockC[idcoord].resize(nintD,1);
  //        VBlockD[idcoord].resize(nintC,1);
  //      }
  //      IntegralBlock2C_dR(0,iBlock[0],atoms[idC],atoms[idD],RCD,cost,sint,cosp,sinp);             //d(mn|ls)/dR
  //      IntegralBlock2C_dA(0,iBlock[1],2,atoms[idC],atoms[idD],RCD,cost,sint,cosp,sinp);           //d(mn|ls)/dt
  //      IntegralBlock2C_dA(0,iBlock[2],3,atoms[idC],atoms[idD],RCD,cost,sint,cosp,sinp);           //d(mn|ls)/dp
  //      IntegralBlock2C_dR2(0,iBlock[3],atoms[idC],atoms[idD],RCD,cost,sint,cosp,sinp);            //d2(mn|ls)/dR2
  //      IntegralBlock2C_dRA(0,iBlock[6],2,atoms[idC],atoms[idD],RCD,cost,sint,cosp,sinp);          //d2(mn|ls)/dRdt
  //      IntegralBlock2C_dRA(0,iBlock[7],3,atoms[idC],atoms[idD],RCD,cost,sint,cosp,sinp);          //d2(mn|ls)/dRdp
  //      IntegralBlock2C_dAB(0,iBlock[4],2,2,atoms[idC],atoms[idD],RCD,cost,sint,cosp,sinp);        //d2(mn|ls)/dt2
  //      IntegralBlock2C_dAB(0,iBlock[5],3,3,atoms[idC],atoms[idD],RCD,cost,sint,cosp,sinp);        //d2(mn|ls)/dp2
  //      IntegralBlock2C_dAB(0,iBlock[8],2,3,atoms[idC],atoms[idD],RCD,cost,sint,cosp,sinp);        //d2(mn|ls)/dtdp
  //      IntegralBlock2C_dR(1,VBlockD[0],atoms[idC],atoms[idD],RCD,cost,sint,cosp,sinp);            //sCsC
  //      IntegralBlock2C_dA(1,VBlockD[1],2,atoms[idC],atoms[idD],RCD,cost,sint,cosp,sinp);          //sCsC
  //      IntegralBlock2C_dA(1,VBlockD[2],3,atoms[idC],atoms[idD],RCD,cost,sint,cosp,sinp);          //sCsC
  //      IntegralBlock2C_dR2(1,VBlockD[3],atoms[idC],atoms[idD],RCD,cost,sint,cosp,sinp);           //sCsC
  //      IntegralBlock2C_dRA(1,VBlockD[6],2,atoms[idC],atoms[idD],RCD,cost,sint,cosp,sinp);         //sCsC
  //      IntegralBlock2C_dRA(1,VBlockD[7],3,atoms[idC],atoms[idD],RCD,cost,sint,cosp,sinp);         //sCsC
  //      IntegralBlock2C_dAB(1,VBlockD[4],2,2,atoms[idC],atoms[idD],RCD,cost,sint,cosp,sinp);       //sCsC
  //      IntegralBlock2C_dAB(1,VBlockD[5],3,3,atoms[idC],atoms[idD],RCD,cost,sint,cosp,sinp);       //sCsC
  //      IntegralBlock2C_dAB(1,VBlockD[8],2,3,atoms[idC],atoms[idD],RCD,cost,sint,cosp,sinp);       //sCsC
  //      IntegralBlock2C_dR(1,VBlockC[0],atoms[idD],atoms[idC],RCD,-cost,sint,-cosp,-sinp);         //sDsD
  //      IntegralBlock2C_dA(1,VBlockC[1],2,atoms[idD],atoms[idC],RCD,-cost,sint,-cosp,-sinp);       //sDsD
  //      IntegralBlock2C_dA(1,VBlockC[2],3,atoms[idD],atoms[idC],RCD,-cost,sint,-cosp,-sinp);       //sDsD
  //      IntegralBlock2C_dR2(1,VBlockC[3],atoms[idD],atoms[idC],RCD,-cost,sint,-cosp,-sinp);        //sDsD
  //      IntegralBlock2C_dRA(1,VBlockC[6],2,atoms[idD],atoms[idC],RCD,-cost,sint,-cosp,-sinp);      //sDsD
  //      IntegralBlock2C_dRA(1,VBlockC[7],3,atoms[idD],atoms[idC],RCD,-cost,sint,-cosp,-sinp);      //sDsD
  //      IntegralBlock2C_dAB(1,VBlockC[4],2,2,atoms[idD],atoms[idC],RCD,-cost,sint,-cosp,-sinp);    //sDsD
  //      IntegralBlock2C_dAB(1,VBlockC[5],3,3,atoms[idD],atoms[idC],RCD,-cost,sint,-cosp,-sinp);    //sDsD
  //      IntegralBlock2C_dAB(1,VBlockC[8],2,3,atoms[idD],atoms[idC],RCD,-cost,sint,-cosp,-sinp);    //sDsD
  //      VBlockC[1] *= -1.0;
  //      VBlockC[6] *= -1.0;
  //      VBlockC[8] *= -1.0;
  //      for (size_t idmu = 0; idmu < AOs[idC]; ++idmu) {
  //        mu = ShiftAO(idmu,transitionmetalC);
  //        for (size_t idnu = idmu; idnu < AOs[idC]; ++idnu) {
  //          nu = ShiftAO(idnu,transitionmetalC);
  //          RCD = 1.0 + double(idmu != idnu);            //factor to account for restricted loops
  //          jcnt = posgamma(mu,nu);
  //          aA = 0;
  //          for (size_t idx = 0; idx < 12; ++idx) {
  //            if (idx == 9) {aA = 6;}
  //            enel_dX[idx] -= RCD*olddens(icntC + idmu + 1,icntC + idnu + 1)*chgD*VBlockD[aA](jcnt + 1,1);
  //            ++aA;
  //          }
  //          for (size_t idc1 = 0; idc1 < 3; ++idc1) {
  //            for (size_t idc2 = 0; idc2 < 3; ++idc2) {
  //              enel_dX[derindex(idc1,idc2)] -= RCD*dPdX[idc1](idmu + 1,idnu + 1)*chgD*VBlockD[idc2](jcnt + 1,1);
  //            }
  //          }
  //          for (size_t idlambda = 0; idlambda < AOs[idD]; ++idlambda) {
  //            lambda = ShiftAO(idlambda,transitionmetalD);
  //            for (size_t idsigma = idlambda; idsigma < AOs[idD]; ++idsigma) {
  //              sigma = ShiftAO(idsigma,transitionmetalD);
  //              kcnt = posgamma(lambda,sigma);
  //              if ((idmu == 0)&&(idnu == 0)) {
  //                RCD = 1.0 + double(idlambda != idsigma);            //factor to account for restricted loops
  //                aA = 0;
  //                for (size_t idx = 0; idx < 12; ++idx) {
  //                  if (idx == 9) {aA = 6;}
  //                  enel_dX[idx] -= RCD*olddens(icntD + idlambda + 1,icntD + idsigma + 1)*chgC*VBlockC[aA](kcnt + 1,1);
  //                  ++aA;
  //                }
  //                for (size_t idc1 = 0; idc1 < 3; ++idc1) {
  //                  for (size_t idc2 = 0; idc2 < 3; ++idc2) {
  //                    enel_dX[derindex(idc1,idc2)] -= RCD*dPdX[idc1](naosC + idlambda + 1,naosC + idsigma + 1)*chgC*VBlockC[idc2](kcnt + 1,1);
  //                  }
  //                }
  //              }
  //              //terms with second derivatives
  //              aux[0] = olddens(icntC + idmu + 1,icntC + idnu + 1)*olddens(icntD + idlambda + 1,icntD + idsigma + 1) - 0.5*olddens(icntC + idmu + 1,icntD + idlambda + 1)*olddens(icntC + idnu + 1,icntD + idsigma + 1);
  //              if (shell == "open") {aux[0] -= 0.5*bolddens(icntC + idmu + 1,icntD + idlambda + 1)*bolddens(icntC + idnu + 1,icntD + idsigma + 1);}
  //              if (idlambda != idsigma) {
  //                aux[0] += olddens(icntC + idmu + 1,icntC + idnu + 1)*olddens(icntD + idsigma + 1,icntD + idlambda + 1) - 0.5*olddens(icntC + idmu + 1,icntD + idsigma + 1)*olddens(icntC + idnu + 1,icntD + idlambda + 1);
  //                if (shell == "open") {aux[0] -= 0.5*bolddens(icntC + idmu + 1,icntD + idsigma + 1)*bolddens(icntC + idnu + 1,icntD + idlambda + 1);}
  //              }
  //              if (idmu != idnu) {
  //                aux[0] += olddens(icntC + idnu + 1,icntC + idmu + 1)*olddens(icntD + idlambda + 1,icntD + idsigma + 1) - 0.5*olddens(icntC + idnu + 1,icntD + idlambda + 1)*olddens(icntC + idmu + 1,icntD + idsigma + 1);
  //                if (shell == "open") {aux[0] -= 0.5*bolddens(icntC + idnu + 1,icntD + idlambda + 1)*bolddens(icntC + idmu + 1,icntD + idsigma + 1);}
  //              }
  //              if ((idlambda != idsigma)&&(idmu != idnu)) {
  //                aux[0] += olddens(icntC + idnu + 1,icntC + idmu + 1)*olddens(icntD + idsigma + 1,icntD + idlambda + 1) - 0.5*olddens(icntC + idnu + 1,icntD + idsigma + 1)*olddens(icntC + idmu + 1,icntD + idlambda + 1);
  //                if (shell == "open") {aux[0] -= 0.5*bolddens(icntC + idnu + 1,icntD + idsigma + 1)*bolddens(icntC + idmu + 1,icntD + idlambda + 1);}
  //              }
  //              aA = 0;
  //              for (size_t idx = 0; idx < 12; ++idx) {
  //                if (idx == 9) {aA = 6;}
  //                enel_dX[idx] += aux[0]*iBlock[aA](jcnt + 1,kcnt + 1);
  //                ++aA;
  //              }
  //              //terms with derivative of density matrix
  //              for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
  //                aux[idcoord] = dPdX[idcoord](idmu + 1,idnu + 1)*olddens(icntD + idlambda + 1,icntD + idsigma + 1) - 0.5*dPdX[idcoord](idmu + 1,naosC + idlambda + 1)*olddens(icntC + idnu + 1,icntD + idsigma + 1);
  //                aux[idcoord] += olddens(icntC + idmu + 1,icntC + idnu + 1)*dPdX[idcoord](naosC + idlambda + 1,naosC + idsigma + 1) - 0.5*olddens(icntC + idmu + 1,icntD + idlambda + 1)*dPdX[idcoord](idnu + 1,naosC + idsigma + 1);
  //                if (shell == "open") {
  //                  aux[idcoord] -= 0.5*dQdX[idcoord](idmu + 1,naosC + idlambda + 1)*bolddens(icntC + idnu + 1,icntD + idsigma + 1);
  //                  aux[idcoord] -= 0.5*bolddens(icntC + idmu + 1,icntD + idlambda + 1)*dQdX[idcoord](idnu + 1,naosC + idsigma + 1);
  //                }
  //                if (idlambda != idsigma) {
  //                  aux[idcoord] += dPdX[idcoord](idmu + 1,idnu + 1)*olddens(icntD + idsigma + 1,icntD + idlambda + 1) - 0.5*dPdX[idcoord](idmu + 1,naosC + idsigma + 1)*olddens(icntC + idnu + 1,icntD + idlambda + 1);
  //                  aux[idcoord] += olddens(icntC + idmu + 1,icntC + idnu + 1)*dPdX[idcoord](naosC + idsigma + 1,naosC + idlambda + 1) - 0.5*olddens(icntC + idmu + 1,icntD + idsigma + 1)*dPdX[idcoord](idnu + 1,naosC + idlambda + 1);
  //                  if (shell == "open") {
  //                    aux[idcoord] -= 0.5*dQdX[idcoord](idmu + 1,naosC + idsigma + 1)*bolddens(icntC + idnu + 1,icntD + idlambda + 1);
  //                    aux[idcoord] -= 0.5*bolddens(icntC + idmu + 1,icntD + idsigma + 1)*dQdX[idcoord](idnu + 1,naosC + idlambda + 1);
  //                  }
  //                }
  //                if (idmu != idnu) {
  //                  aux[idcoord] += dPdX[idcoord](idnu + 1,idmu + 1)*olddens(icntD + idlambda + 1,icntD + idsigma + 1) - 0.5*dPdX[idcoord](idnu + 1,naosC + idlambda + 1)*olddens(icntC + idmu + 1,icntD + idsigma + 1);
  //                  aux[idcoord] += olddens(icntC + idnu + 1,icntC + idmu + 1)*dPdX[idcoord](naosC + idlambda + 1,naosC + idsigma + 1) - 0.5*olddens(icntC + idnu + 1,icntD + idlambda + 1)*dPdX[idcoord](idmu + 1,naosC + idsigma + 1);
  //                  if (shell == "open") {
  //                    aux[idcoord] -= 0.5*dQdX[idcoord](idnu + 1,naosC + idlambda + 1)*bolddens(icntC + idmu + 1,icntD + idsigma + 1);
  //                    aux[idcoord] -= 0.5*bolddens(icntC + idnu + 1,icntD + idlambda + 1)*dQdX[idcoord](idmu + 1,naosC + idsigma + 1);
  //                  }
  //                }
  //                if ((idlambda != idsigma)&&(idmu != idnu)) {
  //                  aux[idcoord] += dPdX[idcoord](idnu + 1,idmu + 1)*olddens(icntD + idsigma + 1,icntD + idlambda + 1) - 0.5*dPdX[idcoord](idnu + 1,naosC + idsigma + 1)*olddens(icntC + idmu + 1,icntD + idlambda + 1);
  //                  aux[idcoord] += olddens(icntC + idnu + 1,icntC + idmu + 1)*dPdX[idcoord](naosC + idsigma + 1,naosC + idlambda + 1) - 0.5*olddens(icntC + idnu + 1,icntD + idsigma + 1)*dPdX[idcoord](idmu + 1,naosC + idlambda + 1);
  //                  if (shell == "open") {
  //                    aux[idcoord] -= 0.5*dQdX[idcoord](idnu + 1,naosC + idsigma + 1)*bolddens(icntC + idmu + 1,icntD + idlambda + 1);
  //                    aux[idcoord] -= 0.5*bolddens(icntC + idnu + 1,icntD + idsigma + 1)*dQdX[idcoord](idmu + 1,naosC + idlambda + 1);
  //                  }
  //                }
  //              }
  //              for (size_t idc1 = 0; idc1 < 3; ++idc1) {
  //                for (size_t idc2 = 0; idc2 < 3; ++idc2) {
  //                  enel_dX[derindex(idc1,idc2)] += aux[idc1]*iBlock[idc2](jcnt + 1,kcnt + 1);
  //                }
  //              }
  //            }
  //          }
  //        }
  //      }
  //      //nuclear energy terms
  //      enel_dX[0] += enuc_dR[icnt];
  //      enel_dX[3] += enuc_dR2[icnt];
  //      //get Hessian contributions
  //      for (size_t idrow = 0; idrow < 3; ++idrow) {
  //        for (size_t idcol = 0; idcol < 3; ++idcol) {
  //          if (idrow == idcol) {spos = idrow;}
  //          else {spos = 2 + idrow + idcol;}
  //          //diagonal terms
  //          hessian(3*idC + idrow + 1,3*idC + idcol + 1) += unitconversion*(dRCD[idrow]*enel_dX[3]*dRCD[idcol] + dtCD[idrow]*enel_dX[4]*dtCD[idcol] + dpCD[idrow]*enel_dX[5]*dpCD[idcol]);
  //          hessian(3*idD + idrow + 1,3*idD + idcol + 1) += unitconversion*(dRCD[idrow]*enel_dX[3]*dRCD[idcol] + dtCD[idrow]*enel_dX[4]*dtCD[idcol] + dpCD[idrow]*enel_dX[5]*dpCD[idcol]);
  //          hessian(3*idC + idrow + 1,3*idD + idcol + 1) -= unitconversion*(dRCD[idrow]*enel_dX[3]*dRCD[idcol] + dtCD[idrow]*enel_dX[4]*dtCD[idcol] + dpCD[idrow]*enel_dX[5]*dpCD[idcol]);
  //          hessian(3*idD + idrow + 1,3*idC + idcol + 1) -= unitconversion*(dRCD[idrow]*enel_dX[3]*dRCD[idcol] + dtCD[idrow]*enel_dX[4]*dtCD[idcol] + dpCD[idrow]*enel_dX[5]*dpCD[idcol]);
  //          //cross terms
  //          hessian(3*idC + idrow + 1,3*idC + idcol + 1) += unitconversion*(dRCD[idrow]*enel_dX[6]*dtCD[idcol] + dtCD[idrow]*enel_dX[9]*dRCD[idcol] + dRCD[idrow]*enel_dX[7]*dpCD[idcol] + dpCD[idrow]*enel_dX[10]*dRCD[idcol] + dtCD[idrow]*enel_dX[8]*dpCD[idcol] + dpCD[idrow]*enel_dX[11]*dtCD[idcol]);
  //          hessian(3*idD + idrow + 1,3*idD + idcol + 1) += unitconversion*(dRCD[idrow]*enel_dX[6]*dtCD[idcol] + dtCD[idrow]*enel_dX[9]*dRCD[idcol] + dRCD[idrow]*enel_dX[7]*dpCD[idcol] + dpCD[idrow]*enel_dX[10]*dRCD[idcol] + dtCD[idrow]*enel_dX[8]*dpCD[idcol] + dpCD[idrow]*enel_dX[11]*dtCD[idcol]);
  //          hessian(3*idC + idrow + 1,3*idD + idcol + 1) -= unitconversion*(dRCD[idrow]*enel_dX[6]*dtCD[idcol] + dtCD[idrow]*enel_dX[9]*dRCD[idcol] + dRCD[idrow]*enel_dX[7]*dpCD[idcol] + dpCD[idrow]*enel_dX[10]*dRCD[idcol] + dtCD[idrow]*enel_dX[8]*dpCD[idcol] + dpCD[idrow]*enel_dX[11]*dtCD[idcol]);
  //          hessian(3*idD + idrow + 1,3*idC + idcol + 1) -= unitconversion*(dRCD[idrow]*enel_dX[6]*dtCD[idcol] + dtCD[idrow]*enel_dX[9]*dRCD[idcol] + dRCD[idrow]*enel_dX[7]*dpCD[idcol] + dpCD[idrow]*enel_dX[10]*dRCD[idcol] + dtCD[idrow]*enel_dX[8]*dpCD[idcol] + dpCD[idrow]*enel_dX[11]*dtCD[idcol]);
  //          //gradient terms
  //          hessian(3*idC + idrow + 1,3*idC + idcol + 1) += unitconversion*(enel_dX[0]*d2RCD[spos] + enel_dX[1]*d2tCD[spos] + enel_dX[2]*d2pCD[spos]);
  //          hessian(3*idD + idrow + 1,3*idD + idcol + 1) += unitconversion*(enel_dX[0]*d2RCD[spos] + enel_dX[1]*d2tCD[spos] + enel_dX[2]*d2pCD[spos]);
  //          hessian(3*idC + idrow + 1,3*idD + idcol + 1) -= unitconversion*(enel_dX[0]*d2RCD[spos] + enel_dX[1]*d2tCD[spos] + enel_dX[2]*d2pCD[spos]);
  //          hessian(3*idD + idrow + 1,3*idC + idcol + 1) -= unitconversion*(enel_dX[0]*d2RCD[spos] + enel_dX[1]*d2tCD[spos] + enel_dX[2]*d2pCD[spos]);
  //        }
  //      }
  //      ++icnt;
  //      icntD += naosD;
  //    }
  //    icntC += naosC;
  //  }
  //  //add the empirical corrections
  //  if (corecorr != "0"){
  //    matrixE dispersioncorrection(1,1);
  //    if ((corecorr == "D3H4X")||(corecorr == "d3h4x")) {hD3H4X(dispersioncorrection,atoms,geometry,this->Type());}
  //    else if ((corecorr == "D3RX")||(corecorr == "d3rx")) {hD3RX(dispersioncorrection,atoms,geometry,this->Type());}
  //    hessian += dispersioncorrection;
  //  }
  //}
  void IntegralBlock1C(matrixE & block, int atom) {
    //function calculating a block of one-center integrals
    block.zero();
    block(1,1) = eri1Center(atom,0,0);
    size_t nrows = block.rows();
    if (nrows > 1) {
      intn[0] = eri1Center(atom,1,1);
      for (size_t idsp = 2; idsp < 5; ++idsp) {
        block(idsp,idsp) = intn[0];                                //<spx|spx>,<spy|spy>,<spz|spz>
      }
      intn[0] = eri1Center(atom,2,2);
      intn[1] = eri1Center(atom,0,2);
      intn[2] = eri1Center(atom,2,-2);
      for (size_t idx = 0; idx < sizepos_pp; ++idx) {
        block(positions_pp[idx],positions_pp[idx]) = intn[0];      //<pxpx|pxpx>,<pypy|pypy>,<pzpz|pzpz>
        block(1,positions_pp[idx]) = intn[1];                      //<ss|pxpx>,<ss|pypy>,<ss|pzpz>
        for (size_t idy = idx + 1; idy < sizepos_pp; ++idy) {
          block(positions_pp[idx],positions_pp[idy]) = intn[2];    //<pxpx|pypy>,<pxpx|pzpz>,<pypy|pzpz>
        }
      }
      intn[0] = 0.5*(eri1Center(atom,2,2) - eri1Center(atom,2,-2));
      block(6,6) = intn[0];                                        //<pxpy|pxpy>
      block(7,7) = intn[0];                                        //<pxpz|pxpz>
      block(9,9) = intn[0];                                        //<pypz|pypz>
      if (nrows > 10) {                                    //then d orbitals
        //basic quantities
        intn[0] = SlaterCondonRadialIntegral(atom,0,1,6);        //F0sd
        intn[1] = SlaterCondonRadialIntegral(atom,0,3,6);        //F0pd
        intn[2] = SlaterCondonRadialIntegral(atom,0,6,6);        //F0dd
        intn[3] = SlaterCondonRadialIntegral(atom,1,5,5);        //G1pd
        intn[4] = SlaterCondonRadialIntegral(atom,1,2,5);        //r1sppd
        intn[5] = SlaterCondonRadialIntegral(atom,2,4,4);        //G2sd
        intn[6] = SlaterCondonRadialIntegral(atom,2,3,6);        //F2pd
        intn[7] = SlaterCondonRadialIntegral(atom,2,6,6);        //F2dd
        intn[8] = SlaterCondonRadialIntegral(atom,2,3,4);        //r2sdpp
        intn[9] = SlaterCondonRadialIntegral(atom,2,4,6);        //r2sddd
        intn[10] = SlaterCondonRadialIntegral(atom,3,5,5);       //G3pd
        intn[11] = SlaterCondonRadialIntegral(atom,4,6,6);       //F4dd
        for (size_t idp1 = 1; idp1 < 4; ++idp1) {
          for (size_t idp2 = 1; idp2 < 4; ++idp2) {
            if (idp1 == idp2) {
              if (idp1 < 3) {
                block(1 + idp1,5*(2 + idp2) + 3) = -sqrt(1.0/45.0)*intn[4];                                           //(spx|pxdz2),(spy|pydz2)
                block(1 + idp1,5*(2 + idp2) + 5) = pow(-1.0,idp1 - 1)*sqrt(3.0/45.0)*intn[4];                         //(spx|pxdx2-y2),(spy|pydx2-y2)
              }
            }
            else if ((idp1 < 3)&&(idp2 == 3)) {
              block(1 + idp2,5*(2 + idp1) + idp1) = sqrt(3.0/45.0)*intn[4];                                           //(spz|pxdxz),(spz|pydyz)
              block(1 + idp1,25 + idp1) = sqrt(3.0/45.0)*intn[4];                                                     //(spx|pzdxz),(spy|pydxz)
            }
          }
        }
        int cnt1 = 1;
        int cnt2 = 1;
        for (size_t ippb1 = 1; ippb1 < 4; ++ippb1) {
          for (size_t ippb2 = ippb1; ippb2 < 4; ++ippb2, ++cnt1) {
            cnt2 = 1;
            for (size_t iddk1 = 1; iddk1 < 6; ++iddk1) {
              if ((ippb1 == ippb2)&&(iddk1 == 3)) {
                if (ippb1 < 3) {block(4 + cnt1,10 + iddk1) = -intn[8]/sqrt(125.0);}                                   //(pxpx|sdz2),(pypy|sdz2)
              }
              else if ((ippb1 == ippb2)&&(ippb1 < 3)&&(iddk1 == 5)) {block(4 + cnt1,10 + iddk1) = pow(-1.0,ippb1 - 1)*sqrt(3.0/125.0)*intn[8];} //(pxpx|sdx2-y2),(pypy|sdx2-y2)
              else if ((ippb1 == iddk1)&&(ippb2 == 3)) {block(4 + cnt1,10 + iddk1) = sqrt(3.0/125.0)*intn[8];}        //(pxpz|sdxz),(pypz|sdyz)
              for (size_t iddk2 = iddk1; iddk2 < 6; ++iddk2, ++cnt2) {
                if ((ippb1 == ippb2)&&(iddk1 == iddk2)) {
                  if (ippb1 == 3) {
                    if (iddk1 < 3) {intn[12] = 2.0;}
                    else if (iddk1 == 3) {intn[12] = 4.0;}
                    else if (iddk1 > 3) {intn[12] = -4.0;}
                  }
                  else if (ippb1 < 3) {
                    if ((iddk1 > 3)||(ippb1 == iddk1)) {intn[12] = 2.0;}
                    else if ((iddk1 < 3)&&(ippb1 != iddk1)) {intn[12] = -4.0;}
                    else if (iddk1 == 3) {intn[12] = -2.0;}
                  }
                  block(4 + cnt1,30 + cnt2) = intn[1] + intn[12]*intn[6]/35.0;                                         //(pp|dd)
                }
                else if ((ippb1 != ippb2)&&(iddk1 != iddk2)) {
                  if (ippb2 == 3) {
                    if ((ippb1 == iddk1)&&(iddk2 == 5)) {block(4 + cnt1,30 + cnt2) = pow(-1.0,iddk1 - 1)*3.0*intn[6]/35.0;} //(pxpz|dxzdx2-y2),(pypz|dyzdx2-y2)
                    else if ((iddk2 == 4)&&(ippb1 + iddk1 == 3)) {block(4 + cnt1,30 + cnt2) = 3.0*intn[6]/35.0;}      //(pxpz|dyzdxy),(pypz|dxzdxy)
                    else if ((iddk2 == 3)&&(ippb1 == iddk1)) {block(4 + cnt1,30 + cnt2) = sqrt(3.0)*intn[6]/35.0;}    //(pxpz|dxzdz2),(pypz|dyzdz2)
                  }
                }
              }
            }
          }
        }
        cnt1 = 1;
        for (int idddb1 = 1; idddb1 < 6; ++idddb1) {
          block(10 + idddb1,10 + idddb1) = 0.2*intn[5];                                                               //(sdxz|sdxz),(sdyz|sdyz),(sdz2|sdz2),(sdxy|sdxy),(sdx2-y2|sdx2-y2)
          if (idddb1 < 3) {block(10 + idddb1,37 + idddb1) = pow(-1.0,idddb1 - 1)*sqrt(3.0/245.0)*intn[9];}            //(sdxz|dyzdxy),(sdyz|dyzdx2-y2)
          if (idddb1 != 3) {block(10 + idddb1,36 - idddb1) = sqrt(3.0/245.0)*intn[9];}                                //(sdxz|dxzdx2-y2),(sdyz|dxzdxy),(sdxy|dxzdyz),(sdx2-y2|dxzdxz)
          for (int idddb2 = idddb1; idddb2 < 6; ++idddb2, ++cnt1) {
            //bra and ket alike
            if (idddb1 == idddb2) {
              block(30 + cnt1,30 + cnt1) = intn[2] + 4.0*intn[7]/49.0 + 36.0*intn[11]/441.0;                          //(dxzdxz|dxzdxz),(dyzdyz|dyzdyz),(dz2dz2|dz2dz2),(dxydxy|dxydxy),(dx2-y2dx2-y2|dx2-y2dx2-y2)
              block(1,30 + cnt1) = intn[0];                                                                           //(ss|dxzdxz),(ss|dyzdyz),(ss|dz2dxz2),(ss|dxydxy),(ss|dx2-y2dx2-y2)
              intn[12] = 1.0;
              if (idddb1 == 3) {intn[12] = 2.0;}
              else if (idddb1 > 3) {intn[12] = -2.0;}
              block(13,30 + cnt1) = intn[12]*intn[9]/sqrt(245.0);                                                      //(sdz2|dd)
            }
            else {
              if ((idddb1 < 3)&&(idddb2 == 3)) {
                block(30 + cnt1,30 + cnt1) = intn[7]/49.0 + 30.0*intn[11]/441.0;                                      //(dxzdz2|dxzdz2),(dyzdz2|dyzdz2)
                block(10 + idddb1,30 + cnt1) = intn[9]/sqrt(245.0);                                                   //(sdxz|dz2dxz),(sdyz|dz2dyz)
              }
              else if (idddb1 == 3) {
                block(30 + cnt1,30 + cnt1) = 4.0*intn[7]/49.0 + 15.0*intn[11]/441.0;                                  //(dz2dxy|dz2dxy),(dz2dx2-y2|dz2dx2-y2)
                block(10 + idddb2,30 + cnt1) = -2.0*intn[9]/sqrt(245.0);                                              //(sdxy|dz2dxy),(sdx2-y2|dz2dx2-y2)
              }
              else if ((idddb1 < 3)&&(idddb2 != 3)) {block(30 + cnt1,30 + cnt1) = 3.0*intn[7]/49.0 + 20.0*intn[11]/441.0;}//(dxzdxy|dxzdxy),(dyzdxy|dyzdxy),(dxzdyz|dxzdyz),(dxzdx2-y2|dxzdx2-y2),(dyzdx2-y2|dyzdx2-y2)
            }
            cnt2 = 1;
            for (int idddk1 = 1; idddk1 < 6; ++idddk1) {
              for (int idddk2 = idddk1; idddk2 < 6; ++idddk2, ++cnt2) {
                if (cnt2 < cnt1) continue;
                if ((idddb1 == idddb2)&&(idddk1 == idddk2)&&(idddb1 == 3)&&((idddk1 > 3))) {block(30 + cnt1,30 + cnt2) = intn[2] - 4.0*intn[7]/49.0 + 6.0*intn[11]/441.0;}       //(dz2dz2|dxydxy),(dz2dz2|dx2-y2dx2-y2)
                else if ((idddb1 == idddb2)&&(idddk1 == idddk2)&&(idddk1 == 3)&&((idddb1 < 3))) {block(30 + cnt1,30 + cnt2) = intn[2] + 2.0*intn[7]/49.0 - 24.0*intn[11]/441.0;} //(dxzdxz|dz2dz2),(dyzdyz|dz2dz2)
                else if ((idddb1 == idddb2)&&(idddk1 == idddk2)&&(idddb1 != idddk1)&&(idddb1 != 3)&&((idddk1 != 3))) {
                  if (((idddb1 == 4)&&(idddk1 == 5))||((idddb1 == 5)&&(idddk1 == 4))) {block(30 + cnt1,30 + cnt2) = intn[2] + 4.0*intn[7]/49.0 - 34.0*intn[11]/441.0;}           //(dxydxy|dx2-y2dx2-y2)
                  else {block(30 + cnt1,30 + cnt2) = intn[2] - 2.0*intn[7]/49.0 - 4.0*intn[11]/441.0;}                                                                           //(dxzdxz|dyzdyz),(dxzdxz|dxydxy),(dxzdxz|dx2-y2dx2-y2),(dyzdyz|dxydxy),(dyzdyz|dx2-y2dx2-y2),(dyzdyz|dxzdxz),(dxydxy|dxzdxz),(dx2-y2dx2-y2|dxzdxz),(dxydxy|dyzdyz),(dx2-y2dx2-y2|dyzdyz)
                }
                else if ((idddk1 == 3)&&(idddk2 == 5)&&(idddb1 < 3)&&(idddb1 == idddb2)) {block(30 + cnt1,30 + cnt2) = pow(-1,idddb1+1)*(-sqrt(12.0)*intn[7]/49.0 + sqrt(300.0)*intn[11]/441.0);} //(dxzdxz|dz2dx2-y2),(dyzdyz|dz2dx2-y2)
                else if (((idddb1 < 3)&&(idddb2 == 3)&&(idddk1 < 3)&&(idddk2 > 3))||((idddb1 < 3)&&(idddb2 > 3)&&(idddk1 < 3)&&(idddk2 == 3))) {
                  if ((idddb1 + idddb2 + idddk1 + idddk2)%2 != 0) {continue;}
                  block(30 + cnt1,30 + cnt2) = pow(-1.0,DeltaDirac(idddk1+idddk2,7))*((sqrt(3.0)/49.0)*intn[7] - (sqrt(75.0)/441.0)*intn[11]);                                   //(dxzdz2|dxzdx2-y2),(dxzdz2|dyzdxy),(dyzdz2|dyzdx2-y2),(dxzdxy|dyzdz2)
                }
                else if (((idddb1 < 3)&&(idddb2 == 5)&&(idddk1 < 3)&&(idddk2 == 4))||((idddb1 < 3)&&(idddb2 == 4)&&(idddk1 < 3)&&(idddk2 == 5))) {
                  if (idddb1 == idddk1) {continue;}
                  block(30 + cnt1,30 + cnt2) = pow(-1.0,idddb1 + idddb2)*((3.0/49.0)*intn[7] - (15.0/441.0)*intn[11]);                                                           //(dxzdx2-y2|dyzdxy),(dxzdxy|dyzdx2-y2)
                }
              }
            }
          }
        }
        for (size_t idx = 0; idx < 12; ++idx) {
          if ((idx == 1)||(idx == 2)||(idx == 5)||(idx == 7)) {continue;}
          block(16 + idx,16 + idx) = 0.2*intn[3] + 24.0*intn[10]/245.0;                                               //(pxdxy|pxdxy),(pxdx2-y2|pxdx2-y2),(pydyz|pydyz),(pydxy|pydxy),(pydx2-y2|pydx2-y2),(pzdxz|pzdxz),(pzdyz|pzdyz)
        }
        cnt1 = 0;
        for (size_t idx = 4; idx < 13; cnt1 += idx, idx += 4) {
          cnt2 = cnt1;
          for (size_t idy = idx; idy < 13; cnt2 += idy, idy += 4) {
            block(17 + cnt1,17 + cnt2) = 3.0*intn[10]/49.0;                                                           //(pxdyz|pzdxy),(pydxz|pzdxy),(pxdyz|pydxz),(pxdyz|pxdyz),(pydxz|pydxz),(pzdxy|pzdxy)
          }
        }
        block(6,32) = 3.0*intn[6]/35.0;                                                                               //(pxpy|dxzdyz)
        block(5,42) = -sqrt(12.0)*intn[6]/35.0;                                                                       //(pxpx|dz2dx2-y2)
        block(6,41) = -sqrt(12.0)*intn[6]/35.0;                                                                       //(pxpy|dz2dxy)
        block(8,42) =  sqrt(12.0)*intn[6]/35.0;                                                                       //(pypy|dz2dx2-y2)
        block(15,36) = -sqrt(3.0/245.0)*intn[9];                                                                      //(sdx2-y2|dyzdyz)
        block(6,14) =  sqrt(3.0/125.0)*intn[8];                                                                       //(pxpy|sdxy)
        block(2,24) =  sqrt(3.0/45.0)*intn[4];                                                                        //(spx|pydxy)
        block(3,19) =  sqrt(3.0/45.0)*intn[4];                                                                        //(spy|pxdxy)
        block(4,28) =  2.0*intn[4]/sqrt(45.0);                                                                        //(spz|pzdz2)
        block(10,13) = 2.0*intn[8]/sqrt(125.0);                                                                       //(pzpz|sdz2)
        block(44,44) = 35.0*intn[11]/441.0;                                                                           //(dxydx2-y2|dxydx2-y2)
        block(32,41) = -sqrt(12.0)*intn[7]/49.0 + sqrt(300.0)*intn[11]/441.0;                                         //(dxzdyz|dz2dxy)
        block(18,18) = intn[3]/15.0 + 18.0*intn[10]/245.0;                                                            //(pxdz2|pxdz2)
        block(23,23) = intn[3]/15.0 + 18.0*intn[10]/245.0;                                                            //(pydz2|pydz2)
        block(28,28) = 4.0*intn[3]/15.0 + 27.0*intn[10]/245.0;                                                        //(pzdz2|pzdz2)
        block(18,20) = -sqrt(3.0)*intn[3]/15.0 - sqrt(27.0)*intn[10]/245.0;                                           //(pxdz2|pxdx2-y2)
        block(18,24) = -sqrt(3.0)*intn[3]/15.0 - sqrt(27.0)*intn[10]/245.0;                                           //(pxdz2|pydxy)
        block(19,23) = -sqrt(3.0)*intn[3]/15.0 - sqrt(27.0)*intn[10]/245.0;                                           //(pxdxy|pydz2)
        block(23,25) =  sqrt(3.0)*intn[3]/15.0 + sqrt(27.0)*intn[10]/245.0;                                           //(pydz2|pydx2-y2)
        block(16,30) =  3.0*intn[10]/49.0;                                                                            //(pxdxz|pzdx2-y2)
        block(22,30) = -3.0*intn[10]/49.0;                                                                            //(pydyz|pzdx2-y2)
        block(30,30) =  3.0*intn[10]/49.0;                                                                            //(pzdx2-y2|pzdx2-y2)
        block(20,24) =  0.2*intn[3] - 21.0*intn[10]/245.0;                                                            //(pxdx2-y2|pydxy)
        block(19,25) = -0.2*intn[3] + 21.0*intn[10]/245.0;                                                            //(pxdxy|pydx2-y2)
        block(16,28) = sqrt(12.0)*intn[3]/15.0 - sqrt(243.0)*intn[10]/245.0;                                          //(pxdxz|pzdz2)
        block(22,28) = sqrt(12.0)*intn[3]/15.0 - sqrt(243.0)*intn[10]/245.0;                                          //(pydyz|pzdz2)
        block(16,22) =  0.2*intn[3] - 6.0*intn[10]/245.0;                                                             //(pxdxz|pydyz)
        block(20,26) =  0.2*intn[3] - 6.0*intn[10]/245.0;                                                             //(pxdx2-y2|pzdxz)
        block(24,26) =  0.2*intn[3] - 6.0*intn[10]/245.0;                                                             //(pydxy|pzdxz)
        block(19,27) =  0.2*intn[3] - 6.0*intn[10]/245.0;                                                             //(pxdxy|pzdyz)
        block(25,27) = -0.2*intn[3] + 6.0*intn[10]/245.0;                                                             //(pydx2-y2|pzdyz)
        block(23,27) = -sqrt(3.0)*intn[3]/15.0 + sqrt(432.0)*intn[10]/245.0;                                          //(pydz2|pzdyz)
        block(18,26) = -sqrt(3.0)*intn[3]/15.0 + sqrt(432.0)*intn[10]/245.0;                                          //(pxdz2|pzdxz)
      }
      //symmetrize
      for (size_t idr = 1; idr < nrows + 1; ++idr) {
        for (size_t idc = idr + 1; idc < nrows + 1; ++idc) {
          block(idc,idr) = block(idr,idc);
        }
      }
    }
  }
  void IntegralBlock2C(int type, matrixE & block, int atmC, int atmD, double RCD, double cost, double sint, double cosp, double sinp) {
    //function calculating a block of two-center integrals
    int cnt1;
    int cnt2;
    size_t nrows = block.rows();
    size_t ncols = block.cols();
    if ((nrows > 1)||(ncols > 1)) {
      SProt = SPtransf(cost,sint,cosp,sinp);
      PProt = PPtransf(cost,sint,cosp,sinp);
      if ((nrows > 10)||(ncols > 10)) {
        SDrot = SDtransf(cost,sint,cosp,sinp);
        PDrot = PDtransf(cost,sint,cosp,sinp);
        DDrot = DDtransf(cost,sint,cosp,sinp);
      }
    }
    D[0] = Dvalue(atmC,1);
    D[1] = Dvalue(atmC,2);
    D[2] = Dvalue(atmD,1);
    D[3] = Dvalue(atmD,2);
    //integral calculation
    block(1,1) = eri2Center(0,0,0,0,0,0,0,0,RCD,D,atmC,atmD,type);       //(ss|ss)
    if (nrows > 1) {
      intn[0] = eri2Center(0,0,1,0,0,0,0,0,RCD,D,atmC,atmD,type);
      intn[1] = eri2Center(1,0,1,0,0,0,0,0,RCD,D,atmC,atmD,type);
      intn[2] = eri2Center(1,1,1,1,0,0,0,0,RCD,D,atmC,atmD,type);
      for (size_t idx = 0; idx < 6; ++idx) {
        if (idx < 3) {block(idx + 2,1) = -intn[0]*SProt(idx + 1,1);}                                   //(sp|ss) integrals
        block(5 + idx,1) = intn[1]*PProt(idx + 1,1) + intn[2]*(PProt(idx + 1,2) + PProt(idx + 1,3));   //(pp|ss) integrals
      }
      if (nrows > 10) {
        //(sd|ss) integrals
        Dd[0] = Dvalue(atmC,4);
        Dd[1] = Dvalue(atmD,1);
        intn[0] = d_eri2Center(0,0,2,0,0,0,0,0,RCD,Dd,atmC,atmD,type);
        for (size_t idsd = 1; idsd < 6; ++idsd) {
          block(10 + idsd,1) = intn[0]*SDrot(idsd,1);
        }
        //(pd|ss) integrals
        Dd[0] = Dvalue(atmC,5);
        intn[0] = d_eri2Center(1,0,2,0,0,0,0,0,RCD,Dd,atmC,atmD,type);
        intn[1] = d_eri2Center(1,1,2,1,0,0,0,0,RCD,Dd,atmC,atmD,type);
        for (int idpd1 = 1; idpd1 < 6; ++idpd1) {
          for (int idpd2 = 1; idpd2 < 4; ++idpd2) {
            block(5*(idpd2+2) + idpd1,1) = -(intn[0]*PDrot((idpd1 - 1)*3 + idpd2,1) + intn[1]*(PDrot((idpd1 - 1)*3 + idpd2,2) + PDrot((idpd1 - 1)*3 + idpd2,3)));
          }
        }
        //(dd|ss) integrals
        Dd[0] = Dvalue(atmC,6);
        intn[0] = d_eri2Center(2,0,2,0,0,0,0,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(2,1,2,1,0,0,0,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[2] = d_eri2Center(2,2,2,2,0,0,0,0,RCD,Dd,atmC,atmD,type);    //(dl,dl|sg,sg)=(-dl,-dl|sg,sg)
        for (size_t iddd = 1; iddd < 16; ++iddd) {
          block(30 + iddd,1) = intn[0]*DDrot(iddd,1) + intn[1]*(DDrot(iddd,2) + DDrot(iddd,3)) + intn[2]*(DDrot(iddd,4) + DDrot(iddd,5));
        }
      }
    }
    if ((ncols > 1)&&(type == 0)) {
      intn[0] = eri2Center(0,0,0,0,0,0,1,0,RCD,D,atmC,atmD,type);
      intn[1] = eri2Center(0,0,0,0,1,0,1,0,RCD,D,atmC,atmD,type);
      intn[2] = eri2Center(0,0,0,0,1,1,1,1,RCD,D,atmC,atmD,type);
      for (size_t idx = 0; idx < 6; ++idx) {
        if (idx < 3) {block(1,idx + 2) = -intn[0]*SProt(idx + 1,1);}                                   //(ss|sp) integrals
        block(1,5 + idx) = intn[1]*PProt(idx + 1,1) + intn[2]*(PProt(idx + 1,2) + PProt(idx + 1,3));   //(ss|pp) integrals
      }
      if (ncols > 10) {
        //(ss|sd) integrals
        Dd[0] = Dvalue(atmC,1);
        Dd[1] = Dvalue(atmD,4);
        intn[0] = d_eri2Center(0,0,0,0,0,0,2,0,RCD,Dd,atmC,atmD,type);
        for (size_t idsd = 1; idsd < 6; ++idsd) {
          block(1,10 + idsd) = intn[0]*SDrot(idsd,1);
        }
        //(ss|pd) integrals
        Dd[1] = Dvalue(atmD,5);
        intn[0] = d_eri2Center(0,0,0,0,1,0,2,0,RCD,Dd,atmC,atmD,type);
        intn[1] = d_eri2Center(0,0,0,0,1,1,2,1,RCD,Dd,atmC,atmD,type);
        for (int idpd1 = 1; idpd1 < 6; ++idpd1) {
          for (int idpd2 = 1; idpd2 < 4; ++idpd2) {
            block(1,5*(idpd2 + 2) + idpd1) = -intn[0]*PDrot((idpd1 - 1)*3 + idpd2,1) - intn[1]*(PDrot((idpd1 - 1)*3 + idpd2,2) + PDrot((idpd1 - 1)*3 + idpd2,3));
          }
        }
        //(ss|dd) integrals
        Dd[1] = Dvalue(atmD,6);
        intn[0] = d_eri2Center(0,0,0,0,2,0,2,0,RCD,Dd,atmC,atmD,type);
        intn[1] = d_eri2Center(0,0,0,0,2,1,2,1,RCD,Dd,atmC,atmD,type);
        intn[2] = d_eri2Center(0,0,0,0,2,2,2,2,RCD,Dd,atmC,atmD,type);
        for (size_t iddd = 1; iddd < 16; ++iddd) {
          block(1,30 + iddd) = intn[0]*DDrot(iddd,1) + intn[1]*(DDrot(iddd,2) + DDrot(iddd,3)) + intn[2]*(DDrot(iddd,4) + DDrot(iddd,5));
        }
      }
    }
    if ((nrows > 1)&&(ncols > 1)&&(type == 0)) {
      intn[0] = eri2Center(0,0,1,0,0,0,1,0,RCD,D,atmC,atmD,type);
      intn[1] = eri2Center(0,0,1,1,0,0,1,1,RCD,D,atmC,atmD,type);
      intn[2] = eri2Center(0,0,1,0,1,0,1,0,RCD,D,atmC,atmD,type);
      intn[3] = eri2Center(0,0,1,0,1,1,1,1,RCD,D,atmC,atmD,type);
      intn[4] = eri2Center(0,0,1,1,1,0,1,1,RCD,D,atmC,atmD,type);
      for (size_t idpb = 1; idpb < 4; ++idpb) {
        for (size_t idpk = 1; idpk < 4; ++idpk) {           //(sp|sp)
          block(idpb + 1,idpk + 1) = intn[0]*SProt(idpb,1)*SProt(idpk,1) + intn[1]*(SProt(idpb,2)*SProt(idpk,2) + SProt(idpb,3)*SProt(idpk,3));
        }
        for (size_t idc = 1; idc < 7; ++idc) {              //(sp|pp)
          block(idpb + 1,idc + 4) = -SProt(idpb,1)*(intn[2]*PProt(idc,1) + intn[3]*(PProt(idc,2) + PProt(idc,3))) - intn[4]*(SProt(idpb,2)*PProt(idc,4) + SProt(idpb,3)*PProt(idc,5));
        }
      }
      intn[0] = eri2Center(1,0,1,0,1,0,1,0,RCD,D,atmC,atmD,type);                          //(ps,ps|ps,ps)
      intn[1] = eri2Center(1,1,1,1,1,0,1,0,RCD,D,atmC,atmD,type);                          //(pp,pp|ps,ps)
      intn[2] = eri2Center(1,0,1,0,1,1,1,1,RCD,D,atmC,atmD,type);                          //(ps,ps|pp,pp)
      intn[3] = eri2Center(1,1,1,1,1,1,1,1,RCD,D,atmC,atmD,type);                          //(pp,pp|pp,pp)
      intn[4] = eri2Center(1,1,1,0,1,1,1,0,RCD,D,atmC,atmD,type);                          //(pp,ps|pp,ps)
      intn[5] = eri2Center(1,1,1,1,1,-1,1,-1,RCD,D,atmC,atmD,type);                        //(pp,pp|pp*,pp*)
      intn[6] = eri2Center(1,-1,1,1,1,-1,1,1,RCD,D,atmC,atmD,type);                        //(pp*,pp|pp*,pp)
      intn[7] = eri2Center(1,0,1,0,0,0,1,0,RCD,D,atmC,atmD,type);
      intn[8] = eri2Center(1,1,1,1,0,0,1,0,RCD,D,atmC,atmD,type);
      intn[9] = eri2Center(1,0,1,1,0,0,1,1,RCD,D,atmC,atmD,type);
      for (size_t idr = 1; idr < 7; ++idr) {
        for (size_t idc = 1; idc < 4; ++idc) {                //(pp|sp)
          block(idr + 4,idc + 1) = -SProt(idc,1)*(intn[7]*PProt(idr,1) + intn[8]*(PProt(idr,2) + PProt(idr,3))) - intn[9]*(SProt(idc,2)*PProt(idr,4) + SProt(idc,3)*PProt(idr,5));
        }
        for (size_t idc = 1; idc < 7; ++idc) {                //(pp|pp)
          block(idr + 4,idc + 4) = intn[0]*PProt(idr,1)*PProt(idc,1) + intn[1]*(PProt(idr,2) + PProt(idr,3))*PProt(idc,1) + intn[2]*PProt(idr,1)*(PProt(idc,2) + PProt(idc,3)) + intn[3]*(PProt(idr,2)*PProt(idc,2) + PProt(idr,3)*PProt(idc,3)) + intn[4]*(PProt(idr,4)*PProt(idc,4) + PProt(idr,5)*PProt(idc,5)) + intn[5]*(PProt(idr,2)*PProt(idc,3) + PProt(idr,3)*PProt(idc,2)) + intn[6]*PProt(idr,6)*PProt(idc,6);
        }
      }
      if (nrows > 10) {
        //(sd|sp) integrals
        Dd[0] = Dvalue(atmC,4);
        Dd[1] = Dvalue(atmD,1);
        intn[0] = d_eri2Center(0,0,2,0,0,0,1,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(0,0,2,1,0,0,1,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        //(sd|sp) rotate
        for (size_t idsd = 1; idsd < 6; ++idsd) {
          for (size_t idsp = 1; idsp < 4; ++idsp) {
            block(10 + idsd,1 + idsp) = -intn[0]*SDrot(idsd,1)*SProt(idsp,1) - intn[1]*(SDrot(idsd,2)*SProt(idsp,2) + SDrot(idsd,3)*SProt(idsp,3));
          }
        }
        //(sd|pp) integrals
        Dd[0] = Dvalue(atmC,4);
        Dd[1] = Dvalue(atmD,3);
        intn[0] = d_eri2Center(0,0,2,0,1,0,1,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(0,0,2,0,1,1,1,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center(0,0,2,1,1,0,1,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[3] = d_eri2Center(0,0,2,2,1,1,1,1,RCD,Dd,atmC,atmD,type);    //(sg,dl|pi,pi)=-(sg,dl|-pi,-pi)=(sg,-dl|pi,-pi)
        //(sd|pp) rotate
        for (size_t idsd = 1; idsd < 6; ++idsd) {
          for (size_t idpp = 1; idpp < 7; ++idpp) {
            block(10 + idsd,4 + idpp) = intn[0]*SDrot(idsd,1)*PProt(idpp,1) + intn[1]*SDrot(idsd,1)*(PProt(idpp,2) + PProt(idpp,3));
            block(10 + idsd,4 + idpp) += intn[2]*(SDrot(idsd,2)*PProt(idpp,4) + SDrot(idsd,3)*PProt(idpp,5));
            block(10 + idsd,4 + idpp) += intn[3]*(SDrot(idsd,4)*PProt(idpp,2) - SDrot(idsd,4)*PProt(idpp,3) + SDrot(idsd,5)*PProt(idpp,6));
          }
        }
        //(pd|sp) integrals
        Dd[0] = Dvalue(atmC,5);
        Dd[1] = Dvalue(atmD,1);
        intn[0] = d_eri2Center(1,0,2,0,0,0,1,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(1,1,2,1,0,0,1,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[2] = d_eri2Center(1,0,2,1,0,0,1,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[3] = d_eri2Center(1,1,2,0,0,0,1,1,RCD,Dd,atmC,atmD,type);    //(pi,sg|sg,pi)=(-pi,sg|sg,-pi)
        intn[4] = d_eri2Center(1,1,2,2,0,0,1,1,RCD,Dd,atmC,atmD,type);    //(pi,dl|sg,pi)=(-pi,-dl|sg,pi)=-(-pi,dl|sg,-pi)=(pi,-dl|sg,-pi)
        //(pd|sp) rotate
        for (size_t idsp = 1; idsp < 4; ++idsp) {
          for (size_t idpd1 = 1; idpd1 < 6; ++idpd1) {
            for (size_t idpd2 = 1; idpd2 < 4; ++idpd2) {
              block(5*(idpd2+2) + idpd1,1 + idsp) = intn[0]*PDrot((idpd1 - 1)*3 + idpd2,1)*SProt(idsp,1) + intn[1]*(PDrot((idpd1 - 1)*3 + idpd2,2) + PDrot((idpd1 - 1)*3 + idpd2,3))*SProt(idsp,1);
              block(5*(idpd2+2) + idpd1,1 + idsp) += intn[2]*(PDrot((idpd1 - 1)*3 + idpd2,4)*SProt(idsp,2) + PDrot((idpd1 - 1)*3 + idpd2,5)*SProt(idsp,3));
              block(5*(idpd2+2) + idpd1,1 + idsp) += intn[3]*(PDrot((idpd1 - 1)*3 + idpd2,8)*SProt(idsp,2) + PDrot((idpd1 - 1)*3 + idpd2,12)*SProt(idsp,3));
              block(5*(idpd2+2) + idpd1,1 + idsp) += intn[4]*((PDrot((idpd1 - 1)*3 + idpd2,10) + PDrot((idpd1 - 1)*3 + idpd2,15))*SProt(idsp,2) + (PDrot((idpd1 - 1)*3 + idpd2,11) - PDrot((idpd1 - 1)*3 + idpd2,14))*SProt(idsp,3));
            }
          }
        }
        //(dd|sp) integrals
        Dd[0] = Dvalue(atmC,6);
        Dd[1] = Dvalue(atmD,1);
        intn[0] = d_eri2Center(2,0,2,0,0,0,1,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(2,1,2,1,0,0,1,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[2] = d_eri2Center(2,2,2,2,0,0,1,0,RCD,Dd,atmC,atmD,type);    //(dl,dl|sg,sg)=(-dl,-dl|sg,sg)
        intn[3] = d_eri2Center(2,0,2,1,0,0,1,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[4] = d_eri2Center(2,1,2,2,0,0,1,1,RCD,Dd,atmC,atmD,type);    //(pi,dl|sg,pi)=(-pi,-dl|sg,pi)=-(-pi,dl|sg,-pi)=(pi,-dl|sg,-pi)
        //(dd|sp) rotate
        for (size_t idsp = 1; idsp < 4; ++idsp) {
          for (size_t iddd = 1; iddd < 16; ++iddd) {
            block(30 + iddd,1 + idsp) = -intn[0]*DDrot(iddd,1)*SProt(idsp,1) - intn[1]*(DDrot(iddd,2) + DDrot(iddd,3))*SProt(idsp,1) - intn[2]*(DDrot(iddd,4) + DDrot(iddd,5))*SProt(idsp,1);
            block(30 + iddd,1 + idsp) -= intn[3]*(DDrot(iddd,6)*SProt(idsp,2) + DDrot(iddd,7)*SProt(idsp,3));
            block(30 + iddd,1 + idsp) -= intn[4]*((DDrot(iddd,11) + DDrot(iddd,14))*SProt(idsp,2) + (DDrot(iddd,12) - DDrot(iddd,13))*SProt(idsp,3));
          }
        }
        //(pd|pp) integrals
        Dd[0] = Dvalue(atmC,5);
        Dd[1] = Dvalue(atmD,3);
        intn[0] = d_eri2Center(1,0,2,0,1,0,1,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(1,1,2,1,1,0,1,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[2] = d_eri2Center(1,0,2,0,1,1,1,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[3] = d_eri2Center(1,1,2,1,1,1,1,1,RCD,Dd,atmC,atmD,type);    //(pi,pi|pi,pi)=(pi,pi|-pi,-pi)=(-pi,-pi|pi,pi)=(-pi,-pi|-pi,-pi)
        intn[4] = d_eri2Center(1,0,2,1,1,0,1,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[5] = d_eri2Center(1,1,2,0,1,1,1,0,RCD,Dd,atmC,atmD,type);    //(pi,sg|pi,sg)=(-pi,sg|-pi,sg)
        intn[6] = d_eri2Center(1,1,2,2,1,0,1,1,RCD,Dd,atmC,atmD,type);    //(pi,dl|sg,pi)=(-pi,-dl|sg,pi)=-(-pi,dl|sg,-pi)=(pi,-dl|sg,-pi)
        //(pd|pp) rotate
        for (size_t idpp = 1; idpp < 7; ++idpp) {
          for (size_t idpd1 = 1; idpd1 < 6; ++idpd1) {
            for (size_t idpd2 = 1; idpd2 < 4; ++idpd2) {
              block(5*(idpd2+2) + idpd1,4 + idpp) = -intn[0]*PDrot((idpd1 - 1)*3 + idpd2,1)*PProt(idpp,1) - intn[1]*(PDrot((idpd1 - 1)*3 + idpd2,2) + PDrot((idpd1 - 1)*3 + idpd2,3))*PProt(idpp,1);
              block(5*(idpd2+2) + idpd1,4 + idpp) -= intn[2]*PDrot((idpd1 - 1)*3 + idpd2,1)*(PProt(idpp,2) + PProt(idpp,3));
              block(5*(idpd2+2) + idpd1,4 + idpp) -= intn[3]*(PDrot((idpd1 - 1)*3 + idpd2,2) + PDrot((idpd1 - 1)*3 + idpd2,3))*(PProt(idpp,2) + PProt(idpp,3));
              block(5*(idpd2+2) + idpd1,4 + idpp) -= intn[4]*(PDrot((idpd1 - 1)*3 + idpd2,4)*PProt(idpp,4) + PDrot((idpd1 - 1)*3 + idpd2,5)*PProt(idpp,5));
              block(5*(idpd2+2) + idpd1,4 + idpp) -= intn[5]*(PDrot((idpd1 - 1)*3 + idpd2,8)*PProt(idpp,4) + PDrot((idpd1 - 1)*3 + idpd2,12)*PProt(idpp,5));
              block(5*(idpd2+2) + idpd1,4 + idpp) -= intn[6]*((PDrot((idpd1 - 1)*3 + idpd2,10) + PDrot((idpd1 - 1)*3 + idpd2,15))*PProt(idpp,4) + (PDrot((idpd1 - 1)*3 + idpd2,11) - PDrot((idpd1 - 1)*3 + idpd2,14))*PProt(idpp,5));
            }
          }
        }
        //(dd|pp) integrals
        Dd[0] = Dvalue(atmC,6);
        Dd[1] = Dvalue(atmD,3);
        intn[0] = d_eri2Center(2,0,2,0,1,0,1,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(2,1,2,1,1,0,1,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[2] = d_eri2Center(2,2,2,2,1,0,1,0,RCD,Dd,atmC,atmD,type);    //(dl,dl|sg,sg)=(-dl,-dl|sg,sg)
        intn[3] = d_eri2Center(2,0,2,0,1,1,1,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[4] = d_eri2Center(2,1,2,1,1,1,1,1,RCD,Dd,atmC,atmD,type);    //(pi,pi|pi,pi)=(-pi,-pi|-pi,-pi)
        intn[5] = d_eri2Center(2,2,2,2,1,1,1,1,RCD,Dd,atmC,atmD,type);    //(dl,dl|pi,pi)=(-dl,-dl|pi,pi)=(dl,dl|-pi,-pi)=(-dl,-dl|-pi,-pi)
        intn[6] = d_eri2Center(2,0,2,1,1,0,1,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[7] = d_eri2Center(2,-1,2,-1,1,1,1,1,RCD,Dd,atmC,atmD,type);  //(-pi,-pi|pi,pi)=(pi,pi|-pi,-pi)
        intn[8] = d_eri2Center(2,1,2,-1,1,1,1,-1,RCD,Dd,atmC,atmD,type);  //(pi,-pi|pi,-pi)
        intn[9] = d_eri2Center(2,0,2,2,1,1,1,1,RCD,Dd,atmC,atmD,type);    //(sg,dl|pi,pi)=(sg,-dl|pi,-pi)=-(sg,dl|-pi,-pi)
        intn[10] = d_eri2Center(2,1,2,2,1,0,1,1,RCD,Dd,atmC,atmD,type);   //(pi,dl|sg,pi)=(-pi,-dl|sg,pi)=-(-pi,dl|sg,-pi)=(pi,-dl|sg,-pi)
        //(dd|pp) rotate
        for (size_t iddd = 1; iddd < 16; ++iddd) {
          for (size_t idpp = 1; idpp < 7; ++idpp) {
            block(30 + iddd,4 + idpp) = intn[0]*DDrot(iddd,1)*PProt(idpp,1) + intn[1]*(DDrot(iddd,2) + DDrot(iddd,3))*PProt(idpp,1) + intn[2]*(DDrot(iddd,4) + DDrot(iddd,5))*PProt(idpp,1);
            block(30 + iddd,4 + idpp) += intn[3]*DDrot(iddd,1)*(PProt(idpp,2) + PProt(idpp,3)) + intn[4]*(DDrot(iddd,2)*PProt(idpp,2) + DDrot(iddd,3)*PProt(idpp,3));
            block(30 + iddd,4 + idpp) += intn[5]*(DDrot(iddd,4) + DDrot(iddd,5))*(PProt(idpp,2) + PProt(idpp,3)) + intn[6]*(DDrot(iddd,6)*PProt(idpp,4) + DDrot(iddd,7)*PProt(idpp,5));
            block(30 + iddd,4 + idpp) += intn[7]*(DDrot(iddd,2)*PProt(idpp,3) + DDrot(iddd,3)*PProt(idpp,2)) + intn[8]*DDrot(iddd,10)*PProt(idpp,6);
            block(30 + iddd,4 + idpp) += intn[9]*(DDrot(iddd,8)*(PProt(idpp,2) - PProt(idpp,3)) + DDrot(iddd,9)*PProt(idpp,6));
            block(30 + iddd,4 + idpp) += intn[10]*((DDrot(iddd,11) + DDrot(iddd,14))*PProt(idpp,4) + (DDrot(iddd,12) - DDrot(iddd,13))*PProt(idpp,5));
          }
        }
      }
      if (ncols > 10) {
        //(sp|sd) integrals
        Dd[0] = Dvalue(atmC,1);
        Dd[1] = Dvalue(atmD,4);
        intn[0] = d_eri2Center(0,0,1,0,0,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(0,0,1,1,0,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        //(sp|sd) rotate
        for (size_t idsp = 1; idsp < 4; ++idsp) {
          for (size_t idsd = 1; idsd < 6; ++idsd) {
            block(1 + idsp,10 + idsd) = -intn[0]*SDrot(idsd,1)*SProt(idsp,1) - intn[1]*(SDrot(idsd,2)*SProt(idsp,2) + SDrot(idsd,3)*SProt(idsp,3));
          }
        }
        //(pp|sd) integrals
        Dd[0] = Dvalue(atmC,3);
        Dd[1] = Dvalue(atmD,4);
        intn[0] = d_eri2Center(1,0,1,0,0,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(1,1,1,1,0,0,2,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[2] = d_eri2Center(1,0,1,1,0,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[3] = d_eri2Center(1,1,1,1,0,0,2,2,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,dl)=-(-pi,-pi|sg,dl)=(pi,-pi|sg,-dl)
        //(pp|sd) rotate
        for (size_t idsd = 1; idsd < 6; ++idsd) {
          for (size_t idpp = 1; idpp < 7; ++idpp) {
            block(4 + idpp,10 + idsd) = intn[0]*SDrot(idsd,1)*PProt(idpp,1) + intn[1]*SDrot(idsd,1)*(PProt(idpp,2) + PProt(idpp,3));
            block(4 + idpp,10 + idsd) += intn[2]*(SDrot(idsd,2)*PProt(idpp,4) + SDrot(idsd,3)*PProt(idpp,5));
            block(4 + idpp,10 + idsd) += intn[3]*(SDrot(idsd,4)*PProt(idpp,2) - SDrot(idsd,4)*PProt(idpp,3) + SDrot(idsd,5)*PProt(idpp,6));
          }
        }
        //(sp|pd) integrals
        Dd[0] = Dvalue(atmC,1);
        Dd[1] = Dvalue(atmD,5);
        intn[0] = d_eri2Center(0,0,1,0,1,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(0,0,1,0,1,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center(0,0,1,1,1,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[3] = d_eri2Center(0,0,1,1,1,1,2,0,RCD,Dd,atmC,atmD,type);    //(sg,pi|pi,sg)=(sg,-pi|-pi,sg)
        intn[4] = d_eri2Center(0,0,1,1,1,1,2,2,RCD,Dd,atmC,atmD,type);    //(sg,pi|pi,dl)=(sg,pi|-pi,-dl)=-(sg,-pi|-pi,dl)=(sg,-pi|pi,-dl)
        //(sp|pd) rotate
        for (size_t idsp = 1; idsp < 4; ++idsp) {
          for (size_t idpd1 = 1; idpd1 < 6; ++idpd1) {
            for (size_t idpd2 = 1; idpd2 < 4; ++idpd2) {
              block(1 + idsp,5*(idpd2+2) + idpd1) = intn[0]*PDrot((idpd1 - 1)*3 + idpd2,1)*SProt(idsp,1) + intn[1]*(PDrot((idpd1 - 1)*3 + idpd2,2) + PDrot((idpd1 - 1)*3 + idpd2,3))*SProt(idsp,1);
              block(1 + idsp,5*(idpd2+2) + idpd1) += intn[2]*(PDrot((idpd1 - 1)*3 + idpd2,4)*SProt(idsp,2) + PDrot((idpd1 - 1)*3 + idpd2,5)*SProt(idsp,3));
              block(1 + idsp,5*(idpd2+2) + idpd1) += intn[3]*(PDrot((idpd1 - 1)*3 + idpd2,8)*SProt(idsp,2) + PDrot((idpd1 - 1)*3 + idpd2,12)*SProt(idsp,3));
              block(1 + idsp,5*(idpd2+2) + idpd1) += intn[4]*((PDrot((idpd1 - 1)*3 + idpd2,10) + PDrot((idpd1 - 1)*3 + idpd2,15))*SProt(idsp,2) + (PDrot((idpd1 - 1)*3 + idpd2,11) - PDrot((idpd1 - 1)*3 + idpd2,14))*SProt(idsp,3));
            }
          }
        }
        //(sp|dd) integrals
        Dd[0] = Dvalue(atmC,1);
        Dd[1] = Dvalue(atmD,6);
        intn[0] = d_eri2Center(0,0,1,0,2,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(0,0,1,0,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center(0,0,1,0,2,2,2,2,RCD,Dd,atmC,atmD,type);    //(sg,sg|dl,dl)=(sg,sg|-dl,-dl)
        intn[3] = d_eri2Center(0,0,1,1,2,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[4] = d_eri2Center(0,0,1,1,2,1,2,2,RCD,Dd,atmC,atmD,type);    //(sg,pi|pi,dl)=(sg,pi|-pi,-dl)=-(sg,-pi|-pi,dl)=(sg,-pi|pi,-dl)
        //(sp|dd) rotate
        for (size_t idsp = 1; idsp < 4; ++idsp) {
          for (size_t iddd = 1; iddd < 16; ++iddd) {
            block(1 + idsp,30 + iddd) = -intn[0]*DDrot(iddd,1)*SProt(idsp,1) - intn[1]*(DDrot(iddd,2) + DDrot(iddd,3))*SProt(idsp,1) - intn[2]*(DDrot(iddd,4) + DDrot(iddd,5))*SProt(idsp,1);
            block(1 + idsp,30 + iddd) -= intn[3]*(DDrot(iddd,6)*SProt(idsp,2) + DDrot(iddd,7)*SProt(idsp,3));
            block(1 + idsp,30 + iddd) -= intn[4]*((DDrot(iddd,11) + DDrot(iddd,14))*SProt(idsp,2) + (DDrot(iddd,12) - DDrot(iddd,13))*SProt(idsp,3));
          }
        }
        //(pp|pd) integrals
        Dd[0] = Dvalue(atmC,3);
        Dd[1] = Dvalue(atmD,5);
        intn[0] = d_eri2Center(1,0,1,0,1,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(1,0,1,0,1,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center(1,1,1,1,1,0,2,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[3] = d_eri2Center(1,1,1,1,1,1,2,1,RCD,Dd,atmC,atmD,type);    //(pi,pi|pi,pi)=(pi,pi|-pi,-pi)=(-pi,-pi|pi,pi)=(-pi,-pi|-pi,-pi)
        intn[4] = d_eri2Center(1,0,1,1,1,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[5] = d_eri2Center(1,1,1,0,1,1,2,0,RCD,Dd,atmC,atmD,type);    //(pi,sg|pi,sg)=(-pi,sg|-pi,sg)
        intn[6] = d_eri2Center(1,0,1,1,1,1,2,2,RCD,Dd,atmC,atmD,type);    //(sg,pi|pi,dl)=(sg,pi|-pi,-dl)=-(sg,-pi|-pi,dl)=(sg,-pi|pi,-dl)
        //(pp|pd) rotate
        for (size_t idpp = 1; idpp < 7; ++idpp) {
          for (size_t idpd1 = 1; idpd1 < 6; ++idpd1) {
            for (size_t idpd2 = 1; idpd2 < 4; ++idpd2) {
              block(4 + idpp,5*(idpd2+2) + idpd1) = -intn[0]*PDrot((idpd1 - 1)*3 + idpd2,1)*PProt(idpp,1) - intn[1]*(PDrot((idpd1 - 1)*3 + idpd2,2) + PDrot((idpd1 - 1)*3 + idpd2,3))*PProt(idpp,1);
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[2]*PDrot((idpd1 - 1)*3 + idpd2,1)*(PProt(idpp,2) + PProt(idpp,3));
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[3]*(PDrot((idpd1 - 1)*3 + idpd2,2) + PDrot((idpd1 - 1)*3 + idpd2,3))*(PProt(idpp,2) + PProt(idpp,3));
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[4]*(PDrot((idpd1 - 1)*3 + idpd2,4)*PProt(idpp,4) + PDrot((idpd1 - 1)*3 + idpd2,5)*PProt(idpp,5));
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[5]*(PDrot((idpd1 - 1)*3 + idpd2,8)*PProt(idpp,4) + PDrot((idpd1 - 1)*3 + idpd2,12)*PProt(idpp,5));
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[6]*((PDrot((idpd1 - 1)*3 + idpd2,10) + PDrot((idpd1 - 1)*3 + idpd2,15))*PProt(idpp,4) + (PDrot((idpd1 - 1)*3 + idpd2,11) - PDrot((idpd1 - 1)*3 + idpd2,14))*PProt(idpp,5));
            }
          }
        }
        //(pp|dd) integrals
        Dd[0] = Dvalue(atmC,3);
        Dd[1] = Dvalue(atmD,6);
        intn[0] = d_eri2Center(1,0,1,0,2,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(1,0,1,0,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center(1,0,1,0,2,2,2,2,RCD,Dd,atmC,atmD,type);    //(sg,sg|dl,dl)=(sg,sg|-dl,-dl)
        intn[3] = d_eri2Center(1,1,1,1,2,0,2,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[4] = d_eri2Center(1,1,1,1,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(pi,pi|pi,pi)=(-pi,-pi|-pi,-pi)
        intn[5] = d_eri2Center(1,1,1,1,2,2,2,2,RCD,Dd,atmC,atmD,type);    //(pi,pi|dl,dl)=(pi,pi|-dl,-dl)=(-pi,-pi|dl,dl)=(-pi,-pi|-dl,-dl)
        intn[6] = d_eri2Center(1,0,1,1,2,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[7] = d_eri2Center(1,-1,1,-1,2,1,2,1,RCD,Dd,atmC,atmD,type);  //(-pi,-pi|pi,pi)=(pi,pi|-pi,-pi)
        intn[8] = d_eri2Center(1,1,1,-1,2,1,2,-1,RCD,Dd,atmC,atmD,type);  //(pi,-pi|pi,-pi)
        intn[9] = d_eri2Center(1,1,1,1,2,0,2,2,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,dl)=(pi,-pi|sg,-dl)=-(-pi,-pi|sg,dl)
        intn[10] = d_eri2Center(1,0,1,1,2,1,2,2,RCD,Dd,atmC,atmD,type);   //(sg,pi|pi,dl)=(sg,pi|-pi,-dl)=-(sg,-pi|-pi,dl)=(sg,-pi|pi,-dl)
        //(pp|dd) rotate
        for (size_t iddd = 1; iddd < 16; ++iddd) {
          for (size_t idpp = 1; idpp < 7; ++idpp) {
            block(4 + idpp,30 + iddd) = intn[0]*DDrot(iddd,1)*PProt(idpp,1) + intn[1]*(DDrot(iddd,2) + DDrot(iddd,3))*PProt(idpp,1) + intn[2]*(DDrot(iddd,4) + DDrot(iddd,5))*PProt(idpp,1);
            block(4 + idpp,30 + iddd) += intn[3]*DDrot(iddd,1)*(PProt(idpp,2) + PProt(idpp,3)) + intn[4]*(DDrot(iddd,2)*PProt(idpp,2) + DDrot(iddd,3)*PProt(idpp,3));
            block(4 + idpp,30 + iddd) += intn[5]*(DDrot(iddd,4) + DDrot(iddd,5))*(PProt(idpp,2) + PProt(idpp,3)) + intn[6]*(DDrot(iddd,6)*PProt(idpp,4) + DDrot(iddd,7)*PProt(idpp,5));
            block(4 + idpp,30 + iddd) += intn[7]*(DDrot(iddd,2)*PProt(idpp,3) + DDrot(iddd,3)*PProt(idpp,2)) + intn[8]*DDrot(iddd,10)*PProt(idpp,6);
            block(4 + idpp,30 + iddd) += intn[9]*(DDrot(iddd,8)*(PProt(idpp,2) - PProt(idpp,3)) + DDrot(iddd,9)*PProt(idpp,6));
            block(4 + idpp,30 + iddd) += intn[10]*((DDrot(iddd,11) + DDrot(iddd,14))*PProt(idpp,4) + (DDrot(iddd,12) - DDrot(iddd,13))*PProt(idpp,5));
          }
        }
      }
      if ((nrows > 10)&&(ncols > 10)) {
        //(sd|sd) integrals
        Dd[0] = Dvalue(atmC,4);
        Dd[1] = Dvalue(atmD,4);
        intn[0] = d_eri2Center(0,0,2,0,0,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(0,0,2,1,0,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[2] = d_eri2Center(0,0,2,2,0,0,2,2,RCD,Dd,atmC,atmD,type);    //(sg,dl|sg,dl)=(sg,-dl|sg,-dl)
        //(sd|sd) rotate
        for (size_t idsdb = 1; idsdb < 6; ++idsdb) {
          for (size_t idsdk = 1; idsdk < 6; ++idsdk) {
            block(10 + idsdb,10 + idsdk) = intn[0]*SDrot(idsdb,1)*SDrot(idsdk,1) + intn[1]*(SDrot(idsdb,2)*SDrot(idsdk,2) + SDrot(idsdb,3)*SDrot(idsdk,3)) + intn[2]*(SDrot(idsdb,4)*SDrot(idsdk,4) + SDrot(idsdb,5)*SDrot(idsdk,5));
          }
        }
        //(sd|pd) integrals
        Dd[1] = Dvalue(atmD,5);
        intn[0] = d_eri2Center(0,0,2,0,1,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(0,0,2,0,1,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center(0,0,2,1,1,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[3] = d_eri2Center(0,0,2,1,1,1,2,2,RCD,Dd,atmC,atmD,type);    //(sg,pi|pi,dl)=-(sg,-pi|-pi,dl)=(sg,-pi|pi,-dl)=(sg,pi|-pi,-dl)
        intn[4] = d_eri2Center(0,0,2,1,1,1,2,0,RCD,Dd,atmC,atmD,type);    //(sg,pi|pi,sg)=(sg,-pi|-pi,sg)
        //(sd|pd) rotate
        for (size_t idsd = 1; idsd < 6; ++idsd) {
          for (size_t idpd1 = 1; idpd1 < 6; ++idpd1) {
            for (size_t idpd2 = 1; idpd2 < 4; ++idpd2) {
              block(10 + idsd,5*(idpd2+2) + idpd1) = -intn[0]*PDrot((idpd1 - 1)*3 + idpd2,1)*SDrot(idsd,1) - intn[1]*(PDrot((idpd1 - 1)*3 + idpd2,2) + PDrot((idpd1 - 1)*3 + idpd2,3))*SDrot(idsd,1);
              block(10 + idsd,5*(idpd2+2) + idpd1) -= intn[2]*(PDrot((idpd1 - 1)*3 + idpd2,4)*SDrot(idsd,2) + PDrot((idpd1 - 1)*3 + idpd2,5)*SDrot(idsd,3));
              block(10 + idsd,5*(idpd2+2) + idpd1) -= intn[3]*((PDrot((idpd1 - 1)*3 + idpd2,10) + PDrot((idpd1 - 1)*3 + idpd2,15))*SDrot(idsd,2) + (PDrot((idpd1 - 1)*3 + idpd2,11) - PDrot((idpd1 - 1)*3 + idpd2,14))*SDrot(idsd,3));
              block(10 + idsd,5*(idpd2+2) + idpd1) -= intn[4]*(PDrot((idpd1 - 1)*3 + idpd2,8)*SDrot(idsd,2) + PDrot((idpd1 - 1)*3 + idpd2,12)*SDrot(idsd,3));
            }
          }
        }
        //(pd|sd) integrals
        Dd[0] = Dvalue(atmC,5);
        Dd[1] = Dvalue(atmD,4);
        intn[0] = d_eri2Center(1,0,2,0,0,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(1,1,2,1,0,0,2,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[2] = d_eri2Center(1,0,2,1,0,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[3] = d_eri2Center(1,1,2,2,0,0,2,1,RCD,Dd,atmC,atmD,type);    //(pi,dl|sg,pi)=-(-pi,dl|sg,-pi)=(pi,-dl|sg,-pi)=(-pi,-dl|sg,pi)
        intn[4] = d_eri2Center(1,1,2,0,0,0,2,1,RCD,Dd,atmC,atmD,type);    //(pi,sg|sg,pi)=(-pi,sg|sg,-pi)
        //(pd|sd) rotate
        for (size_t idpd1 = 1; idpd1 < 6; ++idpd1) {
          for (size_t idpd2 = 1; idpd2 < 4; ++idpd2) {
            for (size_t idsd = 1; idsd < 6; ++idsd) {
              block(5*(idpd2+2) + idpd1,10 + idsd) = -intn[0]*PDrot((idpd1 - 1)*3 + idpd2,1)*SDrot(idsd,1) - intn[1]*(PDrot((idpd1 - 1)*3 + idpd2,2) + PDrot((idpd1 - 1)*3 + idpd2,3))*SDrot(idsd,1);
              block(5*(idpd2+2) + idpd1,10 + idsd) -= intn[2]*(PDrot((idpd1 - 1)*3 + idpd2,4)*SDrot(idsd,2) + PDrot((idpd1 - 1)*3 + idpd2,5)*SDrot(idsd,3));
              block(5*(idpd2+2) + idpd1,10 + idsd) -= intn[3]*((PDrot((idpd1 - 1)*3 + idpd2,10) + PDrot((idpd1 - 1)*3 + idpd2,15))*SDrot(idsd,2) + (PDrot((idpd1 - 1)*3 + idpd2,11) - PDrot((idpd1 - 1)*3 + idpd2,14))*SDrot(idsd,3));
              block(5*(idpd2+2) + idpd1,10 + idsd) -= intn[4]*(PDrot((idpd1 - 1)*3 + idpd2,8)*SDrot(idsd,2) + PDrot((idpd1 - 1)*3 + idpd2,12)*SDrot(idsd,3));
            }
          }
        }
        //(dd|sd) integrals
        Dd[0] = Dvalue(atmC,6);
        Dd[1] = Dvalue(atmD,4);
        intn[0] = d_eri2Center(2,0,2,0,0,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(2,1,2,1,0,0,2,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[2] = d_eri2Center(2,2,2,2,0,0,2,0,RCD,Dd,atmC,atmD,type);    //(dl,dl|sg,sg)=(-dl,-dl|sg,sg)
        intn[3] = d_eri2Center(2,0,2,1,0,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[4] = d_eri2Center(2,0,2,2,0,0,2,2,RCD,Dd,atmC,atmD,type);    //(sg,dl|sg,dl)=(sg,-dl|sg,-dl)
        intn[5] = d_eri2Center(2,1,2,1,0,0,2,2,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,dl)=-(-pi,-pi|sg,dl)=(pi,-pi|sg,-dl)
        intn[6] = d_eri2Center(2,1,2,2,0,0,2,1,RCD,Dd,atmC,atmD,type);    //(pi,dl|sg,pi)=(-pi,-dl|sg,pi)=(pi,-dl|sg,-pi)=-(-pi,dl|sg,-pi)
        //(dd|sd) rotate
        for (size_t iddd = 1; iddd < 16; ++iddd) {
          for (size_t idsd = 1; idsd < 6; ++idsd) {
            block(30 + iddd,10 + idsd) = intn[0]*DDrot(iddd,1)*SDrot(idsd,1) + intn[1]*(DDrot(iddd,2) + DDrot(iddd,3))*SDrot(idsd,1) + intn[2]*(DDrot(iddd,4) + DDrot(iddd,5))*SDrot(idsd,1);
            block(30 + iddd,10 + idsd) += intn[3]*(DDrot(iddd,6)*SDrot(idsd,2) + DDrot(iddd,7)*SDrot(idsd,3)) + intn[4]*(DDrot(iddd,8)*SDrot(idsd,4) + DDrot(iddd,9)*SDrot(idsd,5));
            block(30 + iddd,10 + idsd) += intn[5]*((DDrot(iddd,2) - DDrot(iddd,3))*SDrot(idsd,4) + DDrot(iddd,10)*SDrot(idsd,5));
            block(30 + iddd,10 + idsd) += intn[6]*((DDrot(iddd,11) + DDrot(iddd,14))*SDrot(idsd,2) + (DDrot(iddd,12) - DDrot(iddd,13))*SDrot(idsd,3));
          }
        }
        //(sd|dd) integrals
        Dd[0] = Dvalue(atmC,4);
        Dd[1] = Dvalue(atmD,6);
        intn[0] = d_eri2Center(0,0,2,0,2,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(0,0,2,0,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center(0,0,2,0,2,2,2,2,RCD,Dd,atmC,atmD,type);    //(sg,sg|dl,dl)=(sg,sg|-dl,-dl)
        intn[3] = d_eri2Center(0,0,2,1,2,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[4] = d_eri2Center(0,0,2,2,2,0,2,2,RCD,Dd,atmC,atmD,type);    //(sg,dl|sg,dl)=(sg,-dl|sg,-dl)
        intn[5] = d_eri2Center(0,0,2,2,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,dl|pi,pi)=-(sg,dl|-pi,-pi)=(sg,-dl|pi,-pi)
        intn[6] = d_eri2Center(0,0,2,1,2,1,2,2,RCD,Dd,atmC,atmD,type);    //(sg,pi|pi,dl)=(sg,pi|-pi,-dl)=(sg,-pi|pi,-dl)=-(sg,-pi|-pi,dl)
        //(sd|dd) rotate
        for (size_t idsd = 1; idsd < 6; ++idsd) {
          for (size_t iddd = 1; iddd < 16; ++iddd) {
            block(10 + idsd,30 + iddd) = intn[0]*DDrot(iddd,1)*SDrot(idsd,1) + intn[1]*(DDrot(iddd,2) + DDrot(iddd,3))*SDrot(idsd,1) + intn[2]*(DDrot(iddd,4) + DDrot(iddd,5))*SDrot(idsd,1);
            block(10 + idsd,30 + iddd) += intn[3]*(DDrot(iddd,6)*SDrot(idsd,2) + DDrot(iddd,7)*SDrot(idsd,3)) + intn[4]*(DDrot(iddd,8)*SDrot(idsd,4) + DDrot(iddd,9)*SDrot(idsd,5));
            block(10 + idsd,30 + iddd) += intn[5]*((DDrot(iddd,2) - DDrot(iddd,3))*SDrot(idsd,4) + DDrot(iddd,10)*SDrot(idsd,5));
            block(10 + idsd,30 + iddd) += intn[6]*((DDrot(iddd,11) + DDrot(iddd,14))*SDrot(idsd,2) + (DDrot(iddd,12) - DDrot(iddd,13))*SDrot(idsd,3));
          }
        }
        //(pd|pd) integrals
        Dd[0] = Dvalue(atmC,5);
        Dd[1] = Dvalue(atmD,5);
        intn[0] = d_eri2Center(1,0,2,0,1,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(1,1,2,1,1,0,2,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)=(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center(1,1,2,1,1,1,2,1,RCD,Dd,atmC,atmD,type);    //(pi,pi|pi,pi)=(-pi,-pi|pi,pi)=(pi,pi|-pi,-pi)=(-pi,-pi|-pi,-pi)
        intn[3] = d_eri2Center(1,0,2,1,1,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[4] = d_eri2Center(1,1,2,0,1,1,2,0,RCD,Dd,atmC,atmD,type);    //(pi,sg|pi,sg)=(-pi,sg|-pi,sg)
        intn[5] = d_eri2Center(1,1,2,2,1,1,2,2,RCD,Dd,atmC,atmD,type);    //(pi,dl|pi,dl)=(pi,dl|-pi,-dl)=(-pi,dl|-pi,dl)=-(-pi,dl|pi,-dl)=-(pi,-dl|-pi,dl)=(pi,-dl|pi,-dl)=(-pi,-dl|pi,dl)=(-pi,-dl|-pi,-dl)
        intn[6] = d_eri2Center(1,1,2,0,1,0,2,1,RCD,Dd,atmC,atmD,type);    //(pi,sg|sg,pi)=(-pi,sg|sg,-pi)=(sg,pi|pi,sg)=(sg,-pi|-pi,sg)
        intn[7] = d_eri2Center(1,1,2,2,1,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|pi,dl)=(sg,pi|-pi,-dl)=-(sg,-pi|-pi,dl)=(sg,-pi|pi,-dl)=(pi,dl|sg,pi)=-(-pi,dl|sg,-pi)=(pi,-dl|sg,-pi)=(-pi,-dl|sg,pi)
        intn[8] = d_eri2Center(1,1,2,2,1,1,2,0,RCD,Dd,atmC,atmD,type);    //(pi,sg|pi,dl)=(pi,sg|-pi,-dl)=-(-pi,sg|-pi,dl)=(-pi,sg|pi,-dl)=(pi,dl|pi,sg)=-(-pi,dl|-pi,sg)=(pi,-dl|-pi,sg)=(-pi,-dl|pi,sg)
        //(pd|pd) rotate
        for (size_t idpdb1 = 1; idpdb1 < 6; ++idpdb1) {            //d orbital on bra
          for (size_t idpdb2 = 1; idpdb2 < 4; ++idpdb2) {          //p orbital on bra
            for (size_t idpdk1 = 1; idpdk1 < 6; ++idpdk1) {        //d orbital on ket
              for (size_t idpdk2 = 1; idpdk2 < 4; ++idpdk2) {      //p orbital on ket
                cnt1 = (idpdb1 - 1)*3 + idpdb2;
                cnt2 = (idpdk1 - 1)*3 + idpdk2;
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) = intn[0]*PDrot(cnt1,1)*PDrot(cnt2,1) + intn[1]*((PDrot(cnt1,2) + PDrot(cnt1,3))*PDrot(cnt2,1) + PDrot(cnt1,1)*(PDrot(cnt2,2) + PDrot(cnt2,3)));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[2]*(PDrot(cnt1,2) + PDrot(cnt1,3))*(PDrot(cnt2,2) + PDrot(cnt2,3)) + intn[3]*(PDrot(cnt1,4)*PDrot(cnt2,4) + PDrot(cnt1,5)*PDrot(cnt2,5));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[4]*(PDrot(cnt1,8)*PDrot(cnt2,8) + PDrot(cnt1,12)*PDrot(cnt2,12));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[5]*((PDrot(cnt1,10) + PDrot(cnt1,15))*(PDrot(cnt2,10) + PDrot(cnt2,15)) + (PDrot(cnt1,14) - PDrot(cnt1,11))*(PDrot(cnt2,14) - PDrot(cnt2,11)));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[6]*(PDrot(cnt1,4)*PDrot(cnt2,8) + PDrot(cnt1,8)*PDrot(cnt2,4) + PDrot(cnt1,5)*PDrot(cnt2,12) + PDrot(cnt1,12)*PDrot(cnt2,5));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[7]*((PDrot(cnt1,10) + PDrot(cnt1,15))*PDrot(cnt2,4) + PDrot(cnt1,4)*(PDrot(cnt2,10) + PDrot(cnt2,15)) + (PDrot(cnt1,11) - PDrot(cnt1,14))*PDrot(cnt2,5) + PDrot(cnt1,5)*(PDrot(cnt2,11) - PDrot(cnt2,14)));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[8]*((PDrot(cnt1,10) + PDrot(cnt1,15))*PDrot(cnt2,8) + PDrot(cnt1,8)*(PDrot(cnt2,10) + PDrot(cnt2,15)) + (PDrot(cnt1,11) - PDrot(cnt1,14))*PDrot(cnt2,12) + PDrot(cnt1,12)*(PDrot(cnt2,11) - PDrot(cnt2,14)));
              }
            }
          }
        }
        //(pd|dd) integrals
        Dd[0] = Dvalue(atmC,5);
        Dd[1] = Dvalue(atmD,6);
        intn[0] = d_eri2Center(1,0,2,0,2,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(1,0,2,0,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center(1,1,2,1,2,0,2,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[3] = d_eri2Center(1,1,2,1,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(pi,pi|pi,pi)=(-pi,-pi|-pi,-pi)=(pi,pi|-pi,-pi)=(-pi,-pi|pi,pi)
        intn[4] = d_eri2Center(1,0,2,0,2,2,2,2,RCD,Dd,atmC,atmD,type);    //(sg,sg|dl,dl)=(sg,sg|-dl,-dl)
        intn[5] = d_eri2Center(1,1,2,1,2,2,2,2,RCD,Dd,atmC,atmD,type);    //(pi,pi|dl,dl)=(-pi,-pi|-dl,-dl)=(pi,pi|-dl,-dl)=(-pi,-pi|dl,dl)
        intn[6] = d_eri2Center(1,0,2,1,2,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[7] = d_eri2Center(1,1,2,0,2,1,2,0,RCD,Dd,atmC,atmD,type);    //(pi,sg|pi,sg)=(-pi,sg|-pi,sg)
        intn[8] = d_eri2Center(1,1,2,2,2,1,2,2,RCD,Dd,atmC,atmD,type);    //(pi,dl|pi,dl)=(pi,dl|-pi,-dl)=(-pi,-dl|-pi,-dl)=(-pi,-dl|pi,dl)=(-pi,dl|-pi,dl)=-(-pi,dl|pi,-dl)=(pi,-dl|pi,-dl)=-(pi,-dl|-pi,dl)
        intn[9] = d_eri2Center(1,1,2,0,2,1,2,2,RCD,Dd,atmC,atmD,type);    //(pi,sg|pi,dl)=(pi,sg|-pi,-dl)=-(-pi,sg|-pi,dl)=(-pi,sg|pi,-dl)
        intn[10] = d_eri2Center(1,0,2,1,2,1,2,2,RCD,Dd,atmC,atmD,type);   //(sg,pi|pi,dl)=(sg,pi|-pi,-dl)=-(sg,-pi|-pi,dl)=(sg,-pi|pi,-dl)
        intn[11] = d_eri2Center(1,1,2,2,2,0,2,1,RCD,Dd,atmC,atmD,type);   //(pi,dl|sg,pi)=(-pi,-dl|sg,pi)=-(-pi,dl|sg,-pi)=(pi,-dl|sg,-pi)
        //(pd|dd) rotate
        for (size_t iddd = 1; iddd < 16; ++iddd) {
          for (size_t idpd1 = 1; idpd1 < 6; ++idpd1) {
            for (size_t idpd2 = 1; idpd2 < 4; ++idpd2) {
              cnt1 = (idpd1 - 1)*3 + idpd2;
              block(5*(idpd2+2) + idpd1,30 + iddd) = -intn[0]*PDrot(cnt1,1)*DDrot(iddd,1) - intn[1]*PDrot(cnt1,1)*(DDrot(iddd,2) + DDrot(iddd,3));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[2]*(PDrot(cnt1,2) + PDrot(cnt1,3))*DDrot(iddd,1) + intn[3]*(PDrot(cnt1,2) + PDrot(cnt1,3))*(DDrot(iddd,2) + DDrot(iddd,3));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[4]*PDrot(cnt1,1)*(DDrot(iddd,4) + DDrot(iddd,5)) + intn[5]*(PDrot(cnt1,2) + PDrot(cnt1,3))*(DDrot(iddd,4) + DDrot(iddd,5));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[6]*(PDrot(cnt1,4)*DDrot(iddd,6) + PDrot(cnt1,5)*DDrot(iddd,7)) + intn[7]*(PDrot(cnt1,8)*DDrot(iddd,6) + PDrot(cnt1,12)*DDrot(iddd,7));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[8]*((PDrot(cnt1,10) + PDrot(cnt1,15))*(DDrot(iddd,11) + DDrot(iddd,14)) + (PDrot(cnt1,14) - PDrot(cnt1,11))*(DDrot(iddd,13) - DDrot(iddd,12)));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[9]*(PDrot(cnt1,8)*(DDrot(iddd,11) + DDrot(iddd,14)) + PDrot(cnt1,12)*(DDrot(iddd,12) - DDrot(iddd,13)));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[10]*(PDrot(cnt1,4)*(DDrot(iddd,11) + DDrot(iddd,14)) + PDrot(cnt1,5)*(DDrot(iddd,12) - DDrot(iddd,13)));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[11]*((PDrot(cnt1,10) + PDrot(cnt1,15))*DDrot(iddd,6) + (PDrot(cnt1,11) - PDrot(cnt1,14))*DDrot(iddd,7));
            }
          }
        }
        //(dd|pd) integrals
        Dd[0] = Dvalue(atmC,6);
        Dd[1] = Dvalue(atmD,5);
        intn[0] = d_eri2Center(2,0,2,0,1,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(2,1,2,1,1,0,2,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[2] = d_eri2Center(2,0,2,0,1,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[3] = d_eri2Center(2,1,2,1,1,1,2,1,RCD,Dd,atmC,atmD,type);    //(pi,pi|pi,pi)=(-pi,-pi|-pi,-pi)=(pi,pi|-pi,-pi)=(-pi,-pi|pi,pi)
        intn[4] = d_eri2Center(2,2,2,2,1,0,2,0,RCD,Dd,atmC,atmD,type);    //(dl,dl|sg,sg)=(-dl,-dl|sg,sg)
        intn[5] = d_eri2Center(2,2,2,2,1,1,2,1,RCD,Dd,atmC,atmD,type);    //(dl,dl|pi,pi)=(-dl,-dl|-pi,-pi)=(-dl,-dl|pi,pi)=(dl,dl|-pi,-pi)
        intn[6] = d_eri2Center(2,0,2,1,1,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[7] = d_eri2Center(2,1,2,0,1,1,2,0,RCD,Dd,atmC,atmD,type);    //(pi,sg|pi,sg)=(-pi,sg|-pi,sg)
        intn[8] = d_eri2Center(2,1,2,2,1,1,2,2,RCD,Dd,atmC,atmD,type);    //(pi,dl|pi,dl)=(-pi,-dl|-pi,-dl)=(-pi,-dl|pi,dl)=(pi,dl|-pi,-dl)=(-pi,dl|-pi,dl)=(pi,-dl|pi,-dl)=-(-pi,dl|pi,-dl)=-(pi,-dl|-pi,dl)
        intn[9] = d_eri2Center(2,1,2,2,1,1,2,0,RCD,Dd,atmC,atmD,type);    //(pi,dl|pi,sg)=-(-pi,dl|-pi,sg)=(pi,-dl|-pi,sg)=(-pi,-dl|pi,sg)
        intn[10] = d_eri2Center(2,1,2,2,1,0,2,1,RCD,Dd,atmC,atmD,type);   //(pi,dl|sg,pi)=-(-pi,dl|sg,-pi)=(pi,-dl|sg,-pi)=(-pi,-dl|sg,pi)
        intn[11] = d_eri2Center(2,0,2,1,1,1,2,2,RCD,Dd,atmC,atmD,type);   //(sg,pi|pi,dl)=(sg,pi|-pi,-dl)=-(sg,-pi|-pi,dl)=(sg,-pi|pi,-dl)
        //(dd|pd) rotate
        for (size_t iddd = 1; iddd < 16; ++iddd) {
          for (size_t idpd1 = 1; idpd1 < 6; ++idpd1) {
            for (size_t idpd2 = 1; idpd2 < 4; ++idpd2) {
              cnt1 = (idpd1 - 1)*3 + idpd2;
              block(30 + iddd,5*(idpd2+2) + idpd1) = -intn[0]*PDrot(cnt1,1)*DDrot(iddd,1) - intn[1]*PDrot(cnt1,1)*(DDrot(iddd,2) + DDrot(iddd,3));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[2]*(PDrot(cnt1,2) + PDrot(cnt1,3))*DDrot(iddd,1) + intn[3]*(PDrot(cnt1,2) + PDrot(cnt1,3))*(DDrot(iddd,2) + DDrot(iddd,3));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[4]*PDrot(cnt1,1)*(DDrot(iddd,4) + DDrot(iddd,5)) + intn[5]*(PDrot(cnt1,2) + PDrot(cnt1,3))*(DDrot(iddd,4) + DDrot(iddd,5));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[6]*(PDrot(cnt1,4)*DDrot(iddd,6) + PDrot(cnt1,5)*DDrot(iddd,7)) + intn[7]*(PDrot(cnt1,8)*DDrot(iddd,6) + PDrot(cnt1,12)*DDrot(iddd,7));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[8]*((PDrot(cnt1,10) + PDrot(cnt1,15))*(DDrot(iddd,11) + DDrot(iddd,14)) + (PDrot(cnt1,14) - PDrot(cnt1,11))*(DDrot(iddd,13) - DDrot(iddd,12)));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[9]*(PDrot(cnt1,8)*(DDrot(iddd,11) + DDrot(iddd,14)) + PDrot(cnt1,12)*(DDrot(iddd,12) - DDrot(iddd,13)));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[10]*(PDrot(cnt1,4)*(DDrot(iddd,11) + DDrot(iddd,14)) + PDrot(cnt1,5)*(DDrot(iddd,12) - DDrot(iddd,13)));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[11]*((PDrot(cnt1,10) + PDrot(cnt1,15))*DDrot(iddd,6) + (PDrot(cnt1,11) - PDrot(cnt1,14))*DDrot(iddd,7));
            }
          }
        }
        //(dd|dd) integrals
        Dd[0] = Dvalue(atmC,6);
        Dd[1] = Dvalue(atmD,6);
        intn[0] = d_eri2Center(2,0,2,0,2,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(2,0,2,0,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center(2,0,2,0,2,2,2,2,RCD,Dd,atmC,atmD,type);    //(sg,sg|dl,dl)=(sg,sg|-dl,-dl)
        intn[3] = d_eri2Center(2,1,2,1,2,0,2,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[4] = d_eri2Center(2,2,2,2,2,0,2,0,RCD,Dd,atmC,atmD,type);    //(dl,dl|sg,sg)=(-dl,-dl|sg,sg)
        intn[5] = d_eri2Center(2,2,2,2,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(dl,dl|pi,pi)=(-dl,-dl|pi,pi)=(dl,dl|-pi,-pi)=(-dl,-dl|-pi,-pi)
        intn[6] = d_eri2Center(2,1,2,1,2,2,2,2,RCD,Dd,atmC,atmD,type);    //(pi,pi|dl,dl)=(-pi,-pi|dl,dl)=(pi,pi|-dl,-dl)=(-pi,-pi|-dl,-dl)
        intn[7] = d_eri2Center(2,1,2,1,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(pi,pi|pi,pi)=(-pi,-pi|-pi,-pi)
        intn[8] = d_eri2Center(2,2,2,2,2,2,2,2,RCD,Dd,atmC,atmD,type);    //(dl,dl|dl,dl)=(-dl,-dl|-dl,-dl)
        intn[9] = d_eri2Center(2,1,2,1,2,-1,2,-1,RCD,Dd,atmC,atmD,type);  //(pi,pi|-pi,-pi)=(-pi,-pi|pi,pi)
        intn[10] = d_eri2Center(2,-2,2,-2,2,2,2,2,RCD,Dd,atmC,atmD,type); //(-dl,-dl|dl,dl)=(dl,dl|-dl,-dl)
        intn[11] = d_eri2Center(2,0,2,1,2,0,2,1,RCD,Dd,atmC,atmD,type);   //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[12] = d_eri2Center(2,0,2,2,2,0,2,2,RCD,Dd,atmC,atmD,type);   //(sg,dl|sg,dl)=(sg,-dl|sg,-dl)
        intn[13] = d_eri2Center(2,1,2,2,2,1,2,2,RCD,Dd,atmC,atmD,type);   //(pi,dl|pi,dl)=(pi,dl|-pi,-dl)=(-pi,-dl|pi,dl)=(-pi,-dl|-pi,-dl)=-(-pi,dl|pi,-dl)=-(pi,-dl|-pi,dl)=(-pi,dl|-pi,dl)=(pi,-dl|pi,-dl)
        intn[14] = d_eri2Center(2,1,2,-1,2,1,2,-1,RCD,Dd,atmC,atmD,type); //(pi,-pi|pi,-pi)
        intn[15] = d_eri2Center(2,1,2,-2,2,0,2,-1,RCD,Dd,atmC,atmD,type); //(pi,dl|sg,pi)=(-pi,-dl|sg,pi)=(pi,-dl|sg,-pi)=-(-pi,dl|sg,-pi)=(sg,pi|pi,dl)=(sg,pi|-pi,-dl)=(sg,-pi|pi,-dl)=-(sg,-pi|-pi,dl)
        intn[16] = d_eri2Center(2,1,2,1,2,0,2,2,RCD,Dd,atmC,atmD,type);   //(pi,pi|sg,dl)=-(-pi,-pi|sg,dl)=(pi,-pi|sg,-dl)=(sg,dl|pi,pi)=-(sg,dl|-pi,-pi)=(sg,-dl|pi,-pi)
        //(dd|dd) rotate
        for (size_t iddd1 = 1; iddd1 < 16; ++iddd1) {
          for (size_t iddd2 = 1; iddd2 < 16; ++iddd2) {
            block(30 + iddd1,30 + iddd2) = intn[0]*DDrot(iddd1,1)*DDrot(iddd2,1) + intn[1]*DDrot(iddd1,1)*(DDrot(iddd2,2) + DDrot(iddd2,3)) + intn[2]*DDrot(iddd1,1)*(DDrot(iddd2,4) + DDrot(iddd2,5));
            block(30 + iddd1,30 + iddd2) += intn[3]*(DDrot(iddd1,2) + DDrot(iddd1,3))*DDrot(iddd2,1) + intn[4]*(DDrot(iddd1,4) + DDrot(iddd1,5))*DDrot(iddd2,1);
            block(30 + iddd1,30 + iddd2) += intn[5]*(DDrot(iddd1,4) + DDrot(iddd1,5))*(DDrot(iddd2,2) + DDrot(iddd2,3)) + intn[6]*(DDrot(iddd1,2) + DDrot(iddd1,3))*(DDrot(iddd2,4) + DDrot(iddd2,5));
            block(30 + iddd1,30 + iddd2) += intn[7]*(DDrot(iddd1,2)*DDrot(iddd2,2) + DDrot(iddd1,3)*DDrot(iddd2,3)) + intn[8]*(DDrot(iddd1,4)*DDrot(iddd2,4) + DDrot(iddd1,5)*DDrot(iddd2,5));
            block(30 + iddd1,30 + iddd2) += intn[9]*(DDrot(iddd1,2)*DDrot(iddd2,3) + DDrot(iddd1,3)*DDrot(iddd2,2)) + intn[10]*(DDrot(iddd1,4)*DDrot(iddd2,5) + DDrot(iddd1,5)*DDrot(iddd2,4));
            block(30 + iddd1,30 + iddd2) += intn[11]*(DDrot(iddd1,6)*DDrot(iddd2,6) + DDrot(iddd1,7)*DDrot(iddd2,7)) + intn[12]*(DDrot(iddd1,8)*DDrot(iddd2,8) + DDrot(iddd1,9)*DDrot(iddd2,9));
            block(30 + iddd1,30 + iddd2) += intn[13]*((DDrot(iddd1,11) + DDrot(iddd1,14))*(DDrot(iddd2,11) + DDrot(iddd2,14)) + (DDrot(iddd1,12) - DDrot(iddd1,13))*(DDrot(iddd2,12) - DDrot(iddd2,13)));
            block(30 + iddd1,30 + iddd2) += intn[14]*DDrot(iddd1,10)*DDrot(iddd2,10);
            block(30 + iddd1,30 + iddd2) += intn[15]*((DDrot(iddd1,11) + DDrot(iddd1,14))*DDrot(iddd2,6) + (DDrot(iddd1,12) - DDrot(iddd1,13))*DDrot(iddd2,7) + DDrot(iddd1,6)*(DDrot(iddd2,11) + DDrot(iddd2,14)) + DDrot(iddd1,7)*(DDrot(iddd2,12) - DDrot(iddd2,13)));
            block(30 + iddd1,30 + iddd2) += intn[16]*((DDrot(iddd1,2) - DDrot(iddd1,3))*DDrot(iddd2,8) + DDrot(iddd1,10)*DDrot(iddd2,9) + DDrot(iddd1,8)*(DDrot(iddd2,2) - DDrot(iddd2,3)) + DDrot(iddd1,9)*DDrot(iddd2,10));
          }
        }
      }
    }
  }
  void IntegralBlock2C_dR(int type, matrixE & block, int atmC, int atmD, double RCD, double cost, double sint, double cosp, double sinp) {
    //function calculating a block of first-derivatives of two-center integrals with respect to the internuclear distance
    int cnt1;
    int cnt2;
    size_t nrows = block.rows();
    size_t ncols = block.cols();
    if ((nrows > 1)||(ncols > 1)) {
      SProt = SPtransf(cost,sint,cosp,sinp);
      PProt = PPtransf(cost,sint,cosp,sinp);
      if ((nrows > 10)||(ncols > 10)) {
        SDrot = SDtransf(cost,sint,cosp,sinp);
        PDrot = PDtransf(cost,sint,cosp,sinp);
        DDrot = DDtransf(cost,sint,cosp,sinp);
      }
    }
    D[0] = Dvalue(atmC,1);
    D[1] = Dvalue(atmC,2);
    D[2] = Dvalue(atmD,1);
    D[3] = Dvalue(atmD,2);
    //integral calculation
    block(1,1) = eri2Center_dR(0,0,0,0,0,0,0,0,RCD,D,atmC,atmD,type);       //(ss|ss)
    if (nrows > 1) {
      intn[0] = eri2Center_dR(0,0,1,0,0,0,0,0,RCD,D,atmC,atmD,type);
      intn[1] = eri2Center_dR(1,0,1,0,0,0,0,0,RCD,D,atmC,atmD,type);
      intn[2] = eri2Center_dR(1,1,1,1,0,0,0,0,RCD,D,atmC,atmD,type);
      for (size_t idx = 0; idx < 6; ++idx) {
        if (idx < 3) {block(idx + 2,1) = -intn[0]*SProt(idx + 1,1);}                                                   //(sp|ss) integrals
        block(5 + idx,1) = intn[1]*PProt(idx + 1,1) + intn[2]*(PProt(idx + 1,2) + PProt(idx + 1,3));   //(pp|ss) integrals
      }
      if (nrows > 10) {
        //(sd|ss) integrals
        Dd[0] = Dvalue(atmC,4);
        Dd[1] = Dvalue(atmD,1);
        intn[0] = d_eri2Center_dR(0,0,2,0,0,0,0,0,RCD,Dd,atmC,atmD,type);
        for (size_t idsd = 1; idsd < 6; ++idsd) {
          block(10 + idsd,1) = intn[0]*SDrot(idsd,1);
        }
        //(pd|ss) integrals
        Dd[0] = Dvalue(atmC,5);
        intn[0] = d_eri2Center_dR(1,0,2,0,0,0,0,0,RCD,Dd,atmC,atmD,type);
        intn[1] = d_eri2Center_dR(1,1,2,1,0,0,0,0,RCD,Dd,atmC,atmD,type);
        for (int idpd1 = 1; idpd1 < 6; ++idpd1) {
          for (int idpd2 = 1; idpd2 < 4; ++idpd2) {
            block(5*(idpd2+2) + idpd1,1) = -(intn[0]*PDrot((idpd1 - 1)*3 + idpd2,1) + intn[1]*(PDrot((idpd1 - 1)*3 + idpd2,2) + PDrot((idpd1 - 1)*3 + idpd2,3)));
          }
        }
        //(dd|ss) integrals
        Dd[0] = Dvalue(atmC,6);
        intn[0] = d_eri2Center_dR(2,0,2,0,0,0,0,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR(2,1,2,1,0,0,0,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[2] = d_eri2Center_dR(2,2,2,2,0,0,0,0,RCD,Dd,atmC,atmD,type);    //(dl,dl|sg,sg)=(-dl,-dl|sg,sg)
        for (size_t iddd = 1; iddd < 16; ++iddd) {
          block(30 + iddd,1) = intn[0]*DDrot(iddd,1) + intn[1]*(DDrot(iddd,2) + DDrot(iddd,3)) + intn[2]*(DDrot(iddd,4) + DDrot(iddd,5));
        }
      }
    }
    if ((ncols > 1)&&(type == 0)) {
      intn[0] = eri2Center_dR(0,0,0,0,0,0,1,0,RCD,D,atmC,atmD,type);
      intn[1] = eri2Center_dR(0,0,0,0,1,0,1,0,RCD,D,atmC,atmD,type);
      intn[2] = eri2Center_dR(0,0,0,0,1,1,1,1,RCD,D,atmC,atmD,type);
      for (size_t idx = 0; idx < 6; ++idx) {
        if (idx < 3) {block(1,idx + 2) = -intn[0]*SProt(idx + 1,1);}                                                   //(ss|sp) integrals
        block(1,5 + idx) = intn[1]*PProt(idx + 1,1) + intn[2]*(PProt(idx + 1,2) + PProt(idx + 1,3));   //(ss|pp) integrals
      }
      if (ncols > 10) {
        //(ss|sd) integrals
        Dd[0] = Dvalue(atmC,1);
        Dd[1] = Dvalue(atmD,4);
        intn[0] = d_eri2Center_dR(0,0,0,0,0,0,2,0,RCD,Dd,atmC,atmD,type);
        for (size_t idsd = 1; idsd < 6; ++idsd) {
          block(1,10 + idsd) = intn[0]*SDrot(idsd,1);
        }
        //(ss|pd) integrals
        Dd[1] = Dvalue(atmD,5);
        intn[0] = d_eri2Center_dR(0,0,0,0,1,0,2,0,RCD,Dd,atmC,atmD,type);
        intn[1] = d_eri2Center_dR(0,0,0,0,1,1,2,1,RCD,Dd,atmC,atmD,type);
        for (int idpd1 = 1; idpd1 < 6; ++idpd1) {
          for (int idpd2 = 1; idpd2 < 4; ++idpd2) {
            block(1,5*(idpd2 + 2) + idpd1) = -intn[0]*PDrot((idpd1 - 1)*3 + idpd2,1) - intn[1]*(PDrot((idpd1 - 1)*3 + idpd2,2) + PDrot((idpd1 - 1)*3 + idpd2,3));
          }
        }
        //(ss|dd) integrals
        Dd[1] = Dvalue(atmD,6);
        intn[0] = d_eri2Center_dR(0,0,0,0,2,0,2,0,RCD,Dd,atmC,atmD,type);
        intn[1] = d_eri2Center_dR(0,0,0,0,2,1,2,1,RCD,Dd,atmC,atmD,type);
        intn[2] = d_eri2Center_dR(0,0,0,0,2,2,2,2,RCD,Dd,atmC,atmD,type);
        for (size_t iddd = 1; iddd < 16; ++iddd) {
          block(1,30 + iddd) = intn[0]*DDrot(iddd,1) + intn[1]*(DDrot(iddd,2) + DDrot(iddd,3)) + intn[2]*(DDrot(iddd,4) + DDrot(iddd,5));
        }
      }
    }
    if ((nrows > 1)&&(ncols > 1)&&(type == 0)) {
      intn[0] = eri2Center_dR(0,0,1,0,0,0,1,0,RCD,D,atmC,atmD,type);
      intn[1] = eri2Center_dR(0,0,1,1,0,0,1,1,RCD,D,atmC,atmD,type);
      intn[2] = eri2Center_dR(0,0,1,0,1,0,1,0,RCD,D,atmC,atmD,type);
      intn[3] = eri2Center_dR(0,0,1,0,1,1,1,1,RCD,D,atmC,atmD,type);
      intn[4] = eri2Center_dR(0,0,1,1,1,0,1,1,RCD,D,atmC,atmD,type);
      for (size_t idpb = 1; idpb < 4; ++idpb) {
        for (size_t idpk = 1; idpk < 4; ++idpk) {           //(sp|sp)
          block(idpb + 1,idpk + 1) = intn[0]*SProt(idpb,1)*SProt(idpk,1) + intn[1]*(SProt(idpb,2)*SProt(idpk,2) + SProt(idpb,3)*SProt(idpk,3));
        }
        for (size_t idc = 1; idc < 7; ++idc) {              //(sp|pp)
          block(idpb + 1,idc + 4) = -SProt(idpb,1)*(intn[2]*PProt(idc,1) + intn[3]*(PProt(idc,2) + PProt(idc,3))) - intn[4]*(SProt(idpb,2)*PProt(idc,4) + SProt(idpb,3)*PProt(idc,5));
        }
      }
      intn[0] = eri2Center_dR(1,0,1,0,1,0,1,0,RCD,D,atmC,atmD,type);                          //(ps,ps|ps,ps)
      intn[1] = eri2Center_dR(1,1,1,1,1,0,1,0,RCD,D,atmC,atmD,type);                          //(pp,pp|ps,ps)
      intn[2] = eri2Center_dR(1,0,1,0,1,1,1,1,RCD,D,atmC,atmD,type);                          //(ps,ps|pp,pp)
      intn[3] = eri2Center_dR(1,1,1,1,1,1,1,1,RCD,D,atmC,atmD,type);                          //(pp,pp|pp,pp)
      intn[4] = eri2Center_dR(1,1,1,0,1,1,1,0,RCD,D,atmC,atmD,type);                          //(pp,ps|pp,ps)
      intn[5] = eri2Center_dR(1,1,1,1,1,-1,1,-1,RCD,D,atmC,atmD,type);                        //(pp,pp|pp*,pp*)
      intn[6] = eri2Center_dR(1,-1,1,1,1,-1,1,1,RCD,D,atmC,atmD,type);                        //(pp*,pp|pp*,pp)
      intn[7] = eri2Center_dR(1,0,1,0,0,0,1,0,RCD,D,atmC,atmD,type);
      intn[8] = eri2Center_dR(1,1,1,1,0,0,1,0,RCD,D,atmC,atmD,type);
      intn[9] = eri2Center_dR(1,0,1,1,0,0,1,1,RCD,D,atmC,atmD,type);
      for (size_t idr = 1; idr < 7; ++idr) {
        for (size_t idc = 1; idc < 4; ++idc) {                //(pp|sp)
          block(idr + 4,idc + 1) = -SProt(idc,1)*(intn[7]*PProt(idr,1) + intn[8]*(PProt(idr,2) + PProt(idr,3))) - intn[9]*(SProt(idc,2)*PProt(idr,4) + SProt(idc,3)*PProt(idr,5));
        }
        for (size_t idc = 1; idc < 7; ++idc) {                //(pp|pp)
          block(idr + 4,idc + 4) = intn[0]*PProt(idr,1)*PProt(idc,1) + intn[1]*(PProt(idr,2) + PProt(idr,3))*PProt(idc,1) + intn[2]*PProt(idr,1)*(PProt(idc,2) + PProt(idc,3)) + intn[3]*(PProt(idr,2)*PProt(idc,2) + PProt(idr,3)*PProt(idc,3)) + intn[4]*(PProt(idr,4)*PProt(idc,4) + PProt(idr,5)*PProt(idc,5)) + intn[5]*(PProt(idr,2)*PProt(idc,3) + PProt(idr,3)*PProt(idc,2)) + intn[6]*PProt(idr,6)*PProt(idc,6);
        }
      }
      if (nrows > 10) {
        //(sd|sp) integrals
        Dd[0] = Dvalue(atmC,4);
        Dd[1] = Dvalue(atmD,1);
        intn[0] = d_eri2Center_dR(0,0,2,0,0,0,1,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR(0,0,2,1,0,0,1,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        //(sd|sp) rotate
        for (size_t idsd = 1; idsd < 6; ++idsd) {
          for (size_t idsp = 1; idsp < 4; ++idsp) {
            block(10 + idsd,1 + idsp) = -intn[0]*SDrot(idsd,1)*SProt(idsp,1) - intn[1]*(SDrot(idsd,2)*SProt(idsp,2) + SDrot(idsd,3)*SProt(idsp,3));
          }
        }
        //(sd|pp) integrals
        Dd[0] = Dvalue(atmC,4);
        Dd[1] = Dvalue(atmD,3);
        intn[0] = d_eri2Center_dR(0,0,2,0,1,0,1,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR(0,0,2,0,1,1,1,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center_dR(0,0,2,1,1,0,1,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[3] = d_eri2Center_dR(0,0,2,2,1,1,1,1,RCD,Dd,atmC,atmD,type);    //(sg,dl|pi,pi)=-(sg,dl|-pi,-pi)=(sg,-dl|pi,-pi)
        //(sd|pp) rotate
        for (size_t idsd = 1; idsd < 6; ++idsd) {
          for (size_t idpp = 1; idpp < 7; ++idpp) {
            block(10 + idsd,4 + idpp) = intn[0]*SDrot(idsd,1)*PProt(idpp,1) + intn[1]*SDrot(idsd,1)*(PProt(idpp,2) + PProt(idpp,3));
            block(10 + idsd,4 + idpp) += intn[2]*(SDrot(idsd,2)*PProt(idpp,4) + SDrot(idsd,3)*PProt(idpp,5));
            block(10 + idsd,4 + idpp) += intn[3]*(SDrot(idsd,4)*PProt(idpp,2) - SDrot(idsd,4)*PProt(idpp,3) + SDrot(idsd,5)*PProt(idpp,6));
          }
        }
        //(pd|sp) integrals
        Dd[0] = Dvalue(atmC,5);
        Dd[1] = Dvalue(atmD,1);
        intn[0] = d_eri2Center_dR(1,0,2,0,0,0,1,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR(1,1,2,1,0,0,1,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[2] = d_eri2Center_dR(1,0,2,1,0,0,1,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[3] = d_eri2Center_dR(1,1,2,0,0,0,1,1,RCD,Dd,atmC,atmD,type);    //(pi,sg|sg,pi)=(-pi,sg|sg,-pi)
        intn[4] = d_eri2Center_dR(1,1,2,2,0,0,1,1,RCD,Dd,atmC,atmD,type);    //(pi,dl|sg,pi)=(-pi,-dl|sg,pi)=-(-pi,dl|sg,-pi)=(pi,-dl|sg,-pi)
        //(pd|sp) rotate
        for (size_t idsp = 1; idsp < 4; ++idsp) {
          for (size_t idpd1 = 1; idpd1 < 6; ++idpd1) {
            for (size_t idpd2 = 1; idpd2 < 4; ++idpd2) {
              block(5*(idpd2+2) + idpd1,1 + idsp) = intn[0]*PDrot((idpd1 - 1)*3 + idpd2,1)*SProt(idsp,1) + intn[1]*(PDrot((idpd1 - 1)*3 + idpd2,2) + PDrot((idpd1 - 1)*3 + idpd2,3))*SProt(idsp,1);
              block(5*(idpd2+2) + idpd1,1 + idsp) += intn[2]*(PDrot((idpd1 - 1)*3 + idpd2,4)*SProt(idsp,2) + PDrot((idpd1 - 1)*3 + idpd2,5)*SProt(idsp,3));
              block(5*(idpd2+2) + idpd1,1 + idsp) += intn[3]*(PDrot((idpd1 - 1)*3 + idpd2,8)*SProt(idsp,2) + PDrot((idpd1 - 1)*3 + idpd2,12)*SProt(idsp,3));
              block(5*(idpd2+2) + idpd1,1 + idsp) += intn[4]*((PDrot((idpd1 - 1)*3 + idpd2,10) + PDrot((idpd1 - 1)*3 + idpd2,15))*SProt(idsp,2) + (PDrot((idpd1 - 1)*3 + idpd2,11) - PDrot((idpd1 - 1)*3 + idpd2,14))*SProt(idsp,3));
            }
          }
        }
        //(dd|sp) integrals
        Dd[0] = Dvalue(atmC,6);
        Dd[1] = Dvalue(atmD,1);
        intn[0] = d_eri2Center_dR(2,0,2,0,0,0,1,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR(2,1,2,1,0,0,1,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[2] = d_eri2Center_dR(2,2,2,2,0,0,1,0,RCD,Dd,atmC,atmD,type);    //(dl,dl|sg,sg)=(-dl,-dl|sg,sg)
        intn[3] = d_eri2Center_dR(2,0,2,1,0,0,1,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[4] = d_eri2Center_dR(2,1,2,2,0,0,1,1,RCD,Dd,atmC,atmD,type);    //(pi,dl|sg,pi)=(-pi,-dl|sg,pi)=-(-pi,dl|sg,-pi)=(pi,-dl|sg,-pi)
        //(dd|sp) rotate
        for (size_t idsp = 1; idsp < 4; ++idsp) {
          for (size_t iddd = 1; iddd < 16; ++iddd) {
            block(30 + iddd,1 + idsp) = -intn[0]*DDrot(iddd,1)*SProt(idsp,1) - intn[1]*(DDrot(iddd,2) + DDrot(iddd,3))*SProt(idsp,1) - intn[2]*(DDrot(iddd,4) + DDrot(iddd,5))*SProt(idsp,1);
            block(30 + iddd,1 + idsp) -= intn[3]*(DDrot(iddd,6)*SProt(idsp,2) + DDrot(iddd,7)*SProt(idsp,3));
            block(30 + iddd,1 + idsp) -= intn[4]*((DDrot(iddd,11) + DDrot(iddd,14))*SProt(idsp,2) + (DDrot(iddd,12) - DDrot(iddd,13))*SProt(idsp,3));
          }
        }
        //(pd|pp) integrals
        Dd[0] = Dvalue(atmC,5);
        Dd[1] = Dvalue(atmD,3);
        intn[0] = d_eri2Center_dR(1,0,2,0,1,0,1,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR(1,1,2,1,1,0,1,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[2] = d_eri2Center_dR(1,0,2,0,1,1,1,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[3] = d_eri2Center_dR(1,1,2,1,1,1,1,1,RCD,Dd,atmC,atmD,type);    //(pi,pi|pi,pi)=(pi,pi|-pi,-pi)=(-pi,-pi|pi,pi)=(-pi,-pi|-pi,-pi)
        intn[4] = d_eri2Center_dR(1,0,2,1,1,0,1,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[5] = d_eri2Center_dR(1,1,2,0,1,1,1,0,RCD,Dd,atmC,atmD,type);    //(pi,sg|pi,sg)=(-pi,sg|-pi,sg)
        intn[6] = d_eri2Center_dR(1,1,2,2,1,0,1,1,RCD,Dd,atmC,atmD,type);    //(pi,dl|sg,pi)=(-pi,-dl|sg,pi)=-(-pi,dl|sg,-pi)=(pi,-dl|sg,-pi)
        //(pd|pp) rotate
        for (size_t idpp = 1; idpp < 7; ++idpp) {
          for (size_t idpd1 = 1; idpd1 < 6; ++idpd1) {
            for (size_t idpd2 = 1; idpd2 < 4; ++idpd2) {
              block(5*(idpd2+2) + idpd1,4 + idpp) = -intn[0]*PDrot((idpd1 - 1)*3 + idpd2,1)*PProt(idpp,1) - intn[1]*(PDrot((idpd1 - 1)*3 + idpd2,2) + PDrot((idpd1 - 1)*3 + idpd2,3))*PProt(idpp,1);
              block(5*(idpd2+2) + idpd1,4 + idpp) -= intn[2]*PDrot((idpd1 - 1)*3 + idpd2,1)*(PProt(idpp,2) + PProt(idpp,3));
              block(5*(idpd2+2) + idpd1,4 + idpp) -= intn[3]*(PDrot((idpd1 - 1)*3 + idpd2,2) + PDrot((idpd1 - 1)*3 + idpd2,3))*(PProt(idpp,2) + PProt(idpp,3));
              block(5*(idpd2+2) + idpd1,4 + idpp) -= intn[4]*(PDrot((idpd1 - 1)*3 + idpd2,4)*PProt(idpp,4) + PDrot((idpd1 - 1)*3 + idpd2,5)*PProt(idpp,5));
              block(5*(idpd2+2) + idpd1,4 + idpp) -= intn[5]*(PDrot((idpd1 - 1)*3 + idpd2,8)*PProt(idpp,4) + PDrot((idpd1 - 1)*3 + idpd2,12)*PProt(idpp,5));
              block(5*(idpd2+2) + idpd1,4 + idpp) -= intn[6]*((PDrot((idpd1 - 1)*3 + idpd2,10) + PDrot((idpd1 - 1)*3 + idpd2,15))*PProt(idpp,4) + (PDrot((idpd1 - 1)*3 + idpd2,11) - PDrot((idpd1 - 1)*3 + idpd2,14))*PProt(idpp,5));
            }
          }
        }
        //(dd|pp) integrals
        Dd[0] = Dvalue(atmC,6);
        Dd[1] = Dvalue(atmD,3);
        intn[0] = d_eri2Center_dR(2,0,2,0,1,0,1,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR(2,1,2,1,1,0,1,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[2] = d_eri2Center_dR(2,2,2,2,1,0,1,0,RCD,Dd,atmC,atmD,type);    //(dl,dl|sg,sg)=(-dl,-dl|sg,sg)
        intn[3] = d_eri2Center_dR(2,0,2,0,1,1,1,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[4] = d_eri2Center_dR(2,1,2,1,1,1,1,1,RCD,Dd,atmC,atmD,type);    //(pi,pi|pi,pi)=(-pi,-pi|-pi,-pi)
        intn[5] = d_eri2Center_dR(2,2,2,2,1,1,1,1,RCD,Dd,atmC,atmD,type);    //(dl,dl|pi,pi)=(-dl,-dl|pi,pi)=(dl,dl|-pi,-pi)=(-dl,-dl|-pi,-pi)
        intn[6] = d_eri2Center_dR(2,0,2,1,1,0,1,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[7] = d_eri2Center_dR(2,-1,2,-1,1,1,1,1,RCD,Dd,atmC,atmD,type);  //(-pi,-pi|pi,pi)=(pi,pi|-pi,-pi)
        intn[8] = d_eri2Center_dR(2,1,2,-1,1,1,1,-1,RCD,Dd,atmC,atmD,type);  //(pi,-pi|pi,-pi)
        intn[9] = d_eri2Center_dR(2,0,2,2,1,1,1,1,RCD,Dd,atmC,atmD,type);    //(sg,dl|pi,pi)=(sg,-dl|pi,-pi)=-(sg,dl|-pi,-pi)
        intn[10] = d_eri2Center_dR(2,1,2,2,1,0,1,1,RCD,Dd,atmC,atmD,type);   //(pi,dl|sg,pi)=(-pi,-dl|sg,pi)=-(-pi,dl|sg,-pi)=(pi,-dl|sg,-pi)
        //(dd|pp) rotate
        for (size_t iddd = 1; iddd < 16; ++iddd) {
          for (size_t idpp = 1; idpp < 7; ++idpp) {
            block(30 + iddd,4 + idpp) = intn[0]*DDrot(iddd,1)*PProt(idpp,1) + intn[1]*(DDrot(iddd,2) + DDrot(iddd,3))*PProt(idpp,1) + intn[2]*(DDrot(iddd,4) + DDrot(iddd,5))*PProt(idpp,1);
            block(30 + iddd,4 + idpp) += intn[3]*DDrot(iddd,1)*(PProt(idpp,2) + PProt(idpp,3)) + intn[4]*(DDrot(iddd,2)*PProt(idpp,2) + DDrot(iddd,3)*PProt(idpp,3));
            block(30 + iddd,4 + idpp) += intn[5]*(DDrot(iddd,4) + DDrot(iddd,5))*(PProt(idpp,2) + PProt(idpp,3)) + intn[6]*(DDrot(iddd,6)*PProt(idpp,4) + DDrot(iddd,7)*PProt(idpp,5));
            block(30 + iddd,4 + idpp) += intn[7]*(DDrot(iddd,2)*PProt(idpp,3) + DDrot(iddd,3)*PProt(idpp,2)) + intn[8]*DDrot(iddd,10)*PProt(idpp,6);
            block(30 + iddd,4 + idpp) += intn[9]*(DDrot(iddd,8)*(PProt(idpp,2) - PProt(idpp,3)) + DDrot(iddd,9)*PProt(idpp,6));
            block(30 + iddd,4 + idpp) += intn[10]*((DDrot(iddd,11) + DDrot(iddd,14))*PProt(idpp,4) + (DDrot(iddd,12) - DDrot(iddd,13))*PProt(idpp,5));
          }
        }
      }
      if (ncols > 10) {
        //(sp|sd) integrals
        Dd[0] = Dvalue(atmC,1);
        Dd[1] = Dvalue(atmD,4);
        intn[0] = d_eri2Center_dR(0,0,1,0,0,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR(0,0,1,1,0,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        //(sp|sd) rotate
        for (size_t idsp = 1; idsp < 4; ++idsp) {
          for (size_t idsd = 1; idsd < 6; ++idsd) {
            block(1 + idsp,10 + idsd) = -intn[0]*SDrot(idsd,1)*SProt(idsp,1) - intn[1]*(SDrot(idsd,2)*SProt(idsp,2) + SDrot(idsd,3)*SProt(idsp,3));
          }
        }
        //(pp|sd) integrals
        Dd[0] = Dvalue(atmC,3);
        Dd[1] = Dvalue(atmD,4);
        intn[0] = d_eri2Center_dR(1,0,1,0,0,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR(1,1,1,1,0,0,2,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[2] = d_eri2Center_dR(1,0,1,1,0,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[3] = d_eri2Center_dR(1,1,1,1,0,0,2,2,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,dl)=-(-pi,-pi|sg,dl)=(pi,-pi|sg,-dl)
        //(pp|sd) rotate
        for (size_t idsd = 1; idsd < 6; ++idsd) {
          for (size_t idpp = 1; idpp < 7; ++idpp) {
            block(4 + idpp,10 + idsd) = intn[0]*SDrot(idsd,1)*PProt(idpp,1) + intn[1]*SDrot(idsd,1)*(PProt(idpp,2) + PProt(idpp,3));
            block(4 + idpp,10 + idsd) += intn[2]*(SDrot(idsd,2)*PProt(idpp,4) + SDrot(idsd,3)*PProt(idpp,5));
            block(4 + idpp,10 + idsd) += intn[3]*(SDrot(idsd,4)*PProt(idpp,2) - SDrot(idsd,4)*PProt(idpp,3) + SDrot(idsd,5)*PProt(idpp,6));
          }
        }
        //(sp|pd) integrals
        Dd[0] = Dvalue(atmC,1);
        Dd[1] = Dvalue(atmD,5);
        intn[0] = d_eri2Center_dR(0,0,1,0,1,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR(0,0,1,0,1,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center_dR(0,0,1,1,1,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[3] = d_eri2Center_dR(0,0,1,1,1,1,2,0,RCD,Dd,atmC,atmD,type);    //(sg,pi|pi,sg)=(sg,-pi|-pi,sg)
        intn[4] = d_eri2Center_dR(0,0,1,1,1,1,2,2,RCD,Dd,atmC,atmD,type);    //(sg,pi|pi,dl)=(sg,pi|-pi,-dl)=-(sg,-pi|-pi,dl)=(sg,-pi|pi,-dl)
        //(sp|pd) rotate
        for (size_t idsp = 1; idsp < 4; ++idsp) {
          for (size_t idpd1 = 1; idpd1 < 6; ++idpd1) {
            for (size_t idpd2 = 1; idpd2 < 4; ++idpd2) {
              block(1 + idsp,5*(idpd2+2) + idpd1) = intn[0]*PDrot((idpd1 - 1)*3 + idpd2,1)*SProt(idsp,1) + intn[1]*(PDrot((idpd1 - 1)*3 + idpd2,2) + PDrot((idpd1 - 1)*3 + idpd2,3))*SProt(idsp,1);
              block(1 + idsp,5*(idpd2+2) + idpd1) += intn[2]*(PDrot((idpd1 - 1)*3 + idpd2,4)*SProt(idsp,2) + PDrot((idpd1 - 1)*3 + idpd2,5)*SProt(idsp,3));
              block(1 + idsp,5*(idpd2+2) + idpd1) += intn[3]*(PDrot((idpd1 - 1)*3 + idpd2,8)*SProt(idsp,2) + PDrot((idpd1 - 1)*3 + idpd2,12)*SProt(idsp,3));
              block(1 + idsp,5*(idpd2+2) + idpd1) += intn[4]*((PDrot((idpd1 - 1)*3 + idpd2,10) + PDrot((idpd1 - 1)*3 + idpd2,15))*SProt(idsp,2) + (PDrot((idpd1 - 1)*3 + idpd2,11) - PDrot((idpd1 - 1)*3 + idpd2,14))*SProt(idsp,3));
            }
          }
        }
        //(sp|dd) integrals
        Dd[0] = Dvalue(atmC,1);
        Dd[1] = Dvalue(atmD,6);
        intn[0] = d_eri2Center_dR(0,0,1,0,2,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR(0,0,1,0,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center_dR(0,0,1,0,2,2,2,2,RCD,Dd,atmC,atmD,type);    //(sg,sg|dl,dl)=(sg,sg|-dl,-dl)
        intn[3] = d_eri2Center_dR(0,0,1,1,2,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[4] = d_eri2Center_dR(0,0,1,1,2,1,2,2,RCD,Dd,atmC,atmD,type);    //(sg,pi|pi,dl)=(sg,pi|-pi,-dl)=-(sg,-pi|-pi,dl)=(sg,-pi|pi,-dl)
        //(sp|dd) rotate
        for (size_t idsp = 1; idsp < 4; ++idsp) {
          for (size_t iddd = 1; iddd < 16; ++iddd) {
            block(1 + idsp,30 + iddd) = -intn[0]*DDrot(iddd,1)*SProt(idsp,1) - intn[1]*(DDrot(iddd,2) + DDrot(iddd,3))*SProt(idsp,1) - intn[2]*(DDrot(iddd,4) + DDrot(iddd,5))*SProt(idsp,1);
            block(1 + idsp,30 + iddd) -= intn[3]*(DDrot(iddd,6)*SProt(idsp,2) + DDrot(iddd,7)*SProt(idsp,3));
            block(1 + idsp,30 + iddd) -= intn[4]*((DDrot(iddd,11) + DDrot(iddd,14))*SProt(idsp,2) + (DDrot(iddd,12) - DDrot(iddd,13))*SProt(idsp,3));
          }
        }
        //(pp|pd) integrals
        Dd[0] = Dvalue(atmC,3);
        Dd[1] = Dvalue(atmD,5);
        intn[0] = d_eri2Center_dR(1,0,1,0,1,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR(1,0,1,0,1,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center_dR(1,1,1,1,1,0,2,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[3] = d_eri2Center_dR(1,1,1,1,1,1,2,1,RCD,Dd,atmC,atmD,type);    //(pi,pi|pi,pi)=(pi,pi|-pi,-pi)=(-pi,-pi|pi,pi)=(-pi,-pi|-pi,-pi)
        intn[4] = d_eri2Center_dR(1,0,1,1,1,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[5] = d_eri2Center_dR(1,1,1,0,1,1,2,0,RCD,Dd,atmC,atmD,type);    //(pi,sg|pi,sg)=(-pi,sg|-pi,sg)
        intn[6] = d_eri2Center_dR(1,0,1,1,1,1,2,2,RCD,Dd,atmC,atmD,type);    //(sg,pi|pi,dl)=(sg,pi|-pi,-dl)=-(sg,-pi|-pi,dl)=(sg,-pi|pi,-dl)
        //(pp|pd) rotate
        for (size_t idpp = 1; idpp < 7; ++idpp) {
          for (size_t idpd1 = 1; idpd1 < 6; ++idpd1) {
            for (size_t idpd2 = 1; idpd2 < 4; ++idpd2) {
              block(4 + idpp,5*(idpd2+2) + idpd1) = -intn[0]*PDrot((idpd1 - 1)*3 + idpd2,1)*PProt(idpp,1) - intn[1]*(PDrot((idpd1 - 1)*3 + idpd2,2) + PDrot((idpd1 - 1)*3 + idpd2,3))*PProt(idpp,1);
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[2]*PDrot((idpd1 - 1)*3 + idpd2,1)*(PProt(idpp,2) + PProt(idpp,3));
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[3]*(PDrot((idpd1 - 1)*3 + idpd2,2) + PDrot((idpd1 - 1)*3 + idpd2,3))*(PProt(idpp,2) + PProt(idpp,3));
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[4]*(PDrot((idpd1 - 1)*3 + idpd2,4)*PProt(idpp,4) + PDrot((idpd1 - 1)*3 + idpd2,5)*PProt(idpp,5));
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[5]*(PDrot((idpd1 - 1)*3 + idpd2,8)*PProt(idpp,4) + PDrot((idpd1 - 1)*3 + idpd2,12)*PProt(idpp,5));
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[6]*((PDrot((idpd1 - 1)*3 + idpd2,10) + PDrot((idpd1 - 1)*3 + idpd2,15))*PProt(idpp,4) + (PDrot((idpd1 - 1)*3 + idpd2,11) - PDrot((idpd1 - 1)*3 + idpd2,14))*PProt(idpp,5));
            }
          }
        }
        //(pp|dd) integrals
        Dd[0] = Dvalue(atmC,3);
        Dd[1] = Dvalue(atmD,6);
        intn[0] = d_eri2Center_dR(1,0,1,0,2,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR(1,0,1,0,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center_dR(1,0,1,0,2,2,2,2,RCD,Dd,atmC,atmD,type);    //(sg,sg|dl,dl)=(sg,sg|-dl,-dl)
        intn[3] = d_eri2Center_dR(1,1,1,1,2,0,2,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[4] = d_eri2Center_dR(1,1,1,1,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(pi,pi|pi,pi)=(-pi,-pi|-pi,-pi)
        intn[5] = d_eri2Center_dR(1,1,1,1,2,2,2,2,RCD,Dd,atmC,atmD,type);    //(pi,pi|dl,dl)=(pi,pi|-dl,-dl)=(-pi,-pi|dl,dl)=(-pi,-pi|-dl,-dl)
        intn[6] = d_eri2Center_dR(1,0,1,1,2,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[7] = d_eri2Center_dR(1,-1,1,-1,2,1,2,1,RCD,Dd,atmC,atmD,type);  //(-pi,-pi|pi,pi)=(pi,pi|-pi,-pi)
        intn[8] = d_eri2Center_dR(1,1,1,-1,2,1,2,-1,RCD,Dd,atmC,atmD,type);  //(pi,-pi|pi,-pi)
        intn[9] = d_eri2Center_dR(1,1,1,1,2,0,2,2,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,dl)=(pi,-pi|sg,-dl)=-(-pi,-pi|sg,dl)
        intn[10] = d_eri2Center_dR(1,0,1,1,2,1,2,2,RCD,Dd,atmC,atmD,type);   //(sg,pi|pi,dl)=(sg,pi|-pi,-dl)=-(sg,-pi|-pi,dl)=(sg,-pi|pi,-dl)
        //(pp|dd) rotate
        for (size_t iddd = 1; iddd < 16; ++iddd) {
          for (size_t idpp = 1; idpp < 7; ++idpp) {
            block(4 + idpp,30 + iddd) = intn[0]*DDrot(iddd,1)*PProt(idpp,1) + intn[1]*(DDrot(iddd,2) + DDrot(iddd,3))*PProt(idpp,1) + intn[2]*(DDrot(iddd,4) + DDrot(iddd,5))*PProt(idpp,1);
            block(4 + idpp,30 + iddd) += intn[3]*DDrot(iddd,1)*(PProt(idpp,2) + PProt(idpp,3)) + intn[4]*(DDrot(iddd,2)*PProt(idpp,2) + DDrot(iddd,3)*PProt(idpp,3));
            block(4 + idpp,30 + iddd) += intn[5]*(DDrot(iddd,4) + DDrot(iddd,5))*(PProt(idpp,2) + PProt(idpp,3)) + intn[6]*(DDrot(iddd,6)*PProt(idpp,4) + DDrot(iddd,7)*PProt(idpp,5));
            block(4 + idpp,30 + iddd) += intn[7]*(DDrot(iddd,2)*PProt(idpp,3) + DDrot(iddd,3)*PProt(idpp,2)) + intn[8]*DDrot(iddd,10)*PProt(idpp,6);
            block(4 + idpp,30 + iddd) += intn[9]*(DDrot(iddd,8)*(PProt(idpp,2) - PProt(idpp,3)) + DDrot(iddd,9)*PProt(idpp,6));
            block(4 + idpp,30 + iddd) += intn[10]*((DDrot(iddd,11) + DDrot(iddd,14))*PProt(idpp,4) + (DDrot(iddd,12) - DDrot(iddd,13))*PProt(idpp,5));
          }
        }
      }
      if ((nrows > 10)&&(ncols > 10)) {
        //(sd|sd) integrals
        Dd[0] = Dvalue(atmC,4);
        Dd[1] = Dvalue(atmD,4);
        intn[0] = d_eri2Center_dR(0,0,2,0,0,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR(0,0,2,1,0,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[2] = d_eri2Center_dR(0,0,2,2,0,0,2,2,RCD,Dd,atmC,atmD,type);    //(sg,dl|sg,dl)=(sg,-dl|sg,-dl)
        //(sd|sd) rotate
        for (size_t idsdb = 1; idsdb < 6; ++idsdb) {
          for (size_t idsdk = 1; idsdk < 6; ++idsdk) {
            block(10 + idsdb,10 + idsdk) = intn[0]*SDrot(idsdb,1)*SDrot(idsdk,1) + intn[1]*(SDrot(idsdb,2)*SDrot(idsdk,2) + SDrot(idsdb,3)*SDrot(idsdk,3)) + intn[2]*(SDrot(idsdb,4)*SDrot(idsdk,4) + SDrot(idsdb,5)*SDrot(idsdk,5));
          }
        }
        //(sd|pd) integrals
        Dd[1] = Dvalue(atmD,5);
        intn[0] = d_eri2Center_dR(0,0,2,0,1,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR(0,0,2,0,1,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center_dR(0,0,2,1,1,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[3] = d_eri2Center_dR(0,0,2,1,1,1,2,2,RCD,Dd,atmC,atmD,type);    //(sg,pi|pi,dl)=-(sg,-pi|-pi,dl)=(sg,-pi|pi,-dl)=(sg,pi|-pi,-dl)
        intn[4] = d_eri2Center_dR(0,0,2,1,1,1,2,0,RCD,Dd,atmC,atmD,type);    //(sg,pi|pi,sg)=(sg,-pi|-pi,sg)
        //(sd|pd) rotate
        for (size_t idsd = 1; idsd < 6; ++idsd) {
          for (size_t idpd1 = 1; idpd1 < 6; ++idpd1) {
            for (size_t idpd2 = 1; idpd2 < 4; ++idpd2) {
              block(10 + idsd,5*(idpd2+2) + idpd1) = -intn[0]*PDrot((idpd1 - 1)*3 + idpd2,1)*SDrot(idsd,1) - intn[1]*(PDrot((idpd1 - 1)*3 + idpd2,2) + PDrot((idpd1 - 1)*3 + idpd2,3))*SDrot(idsd,1);
              block(10 + idsd,5*(idpd2+2) + idpd1) -= intn[2]*(PDrot((idpd1 - 1)*3 + idpd2,4)*SDrot(idsd,2) + PDrot((idpd1 - 1)*3 + idpd2,5)*SDrot(idsd,3));
              block(10 + idsd,5*(idpd2+2) + idpd1) -= intn[3]*((PDrot((idpd1 - 1)*3 + idpd2,10) + PDrot((idpd1 - 1)*3 + idpd2,15))*SDrot(idsd,2) + (PDrot((idpd1 - 1)*3 + idpd2,11) - PDrot((idpd1 - 1)*3 + idpd2,14))*SDrot(idsd,3));
              block(10 + idsd,5*(idpd2+2) + idpd1) -= intn[4]*(PDrot((idpd1 - 1)*3 + idpd2,8)*SDrot(idsd,2) + PDrot((idpd1 - 1)*3 + idpd2,12)*SDrot(idsd,3));
            }
          }
        }
        //(pd|sd) integrals
        Dd[0] = Dvalue(atmC,5);
        Dd[1] = Dvalue(atmD,4);
        intn[0] = d_eri2Center_dR(1,0,2,0,0,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR(1,1,2,1,0,0,2,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[2] = d_eri2Center_dR(1,0,2,1,0,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[3] = d_eri2Center_dR(1,1,2,2,0,0,2,1,RCD,Dd,atmC,atmD,type);    //(pi,dl|sg,pi)=-(-pi,dl|sg,-pi)=(pi,-dl|sg,-pi)=(-pi,-dl|sg,pi)
        intn[4] = d_eri2Center_dR(1,1,2,0,0,0,2,1,RCD,Dd,atmC,atmD,type);    //(pi,sg|sg,pi)=(-pi,sg|sg,-pi)
        //(pd|sd) rotate
        for (size_t idpd1 = 1; idpd1 < 6; ++idpd1) {
          for (size_t idpd2 = 1; idpd2 < 4; ++idpd2) {
            for (size_t idsd = 1; idsd < 6; ++idsd) {
              block(5*(idpd2+2) + idpd1,10 + idsd) = -intn[0]*PDrot((idpd1 - 1)*3 + idpd2,1)*SDrot(idsd,1) - intn[1]*(PDrot((idpd1 - 1)*3 + idpd2,2) + PDrot((idpd1 - 1)*3 + idpd2,3))*SDrot(idsd,1);
              block(5*(idpd2+2) + idpd1,10 + idsd) -= intn[2]*(PDrot((idpd1 - 1)*3 + idpd2,4)*SDrot(idsd,2) + PDrot((idpd1 - 1)*3 + idpd2,5)*SDrot(idsd,3));
              block(5*(idpd2+2) + idpd1,10 + idsd) -= intn[3]*((PDrot((idpd1 - 1)*3 + idpd2,10) + PDrot((idpd1 - 1)*3 + idpd2,15))*SDrot(idsd,2) + (PDrot((idpd1 - 1)*3 + idpd2,11) - PDrot((idpd1 - 1)*3 + idpd2,14))*SDrot(idsd,3));
              block(5*(idpd2+2) + idpd1,10 + idsd) -= intn[4]*(PDrot((idpd1 - 1)*3 + idpd2,8)*SDrot(idsd,2) + PDrot((idpd1 - 1)*3 + idpd2,12)*SDrot(idsd,3));
            }
          }
        }
        //(dd|sd) integrals
        Dd[0] = Dvalue(atmC,6);
        Dd[1] = Dvalue(atmD,4);
        intn[0] = d_eri2Center_dR(2,0,2,0,0,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR(2,1,2,1,0,0,2,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[2] = d_eri2Center_dR(2,2,2,2,0,0,2,0,RCD,Dd,atmC,atmD,type);    //(dl,dl|sg,sg)=(-dl,-dl|sg,sg)
        intn[3] = d_eri2Center_dR(2,0,2,1,0,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[4] = d_eri2Center_dR(2,0,2,2,0,0,2,2,RCD,Dd,atmC,atmD,type);    //(sg,dl|sg,dl)=(sg,-dl|sg,-dl)
        intn[5] = d_eri2Center_dR(2,1,2,1,0,0,2,2,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,dl)=-(-pi,-pi|sg,dl)=(pi,-pi|sg,-dl)
        intn[6] = d_eri2Center_dR(2,1,2,2,0,0,2,1,RCD,Dd,atmC,atmD,type);    //(pi,dl|sg,pi)=(-pi,-dl|sg,pi)=(pi,-dl|sg,-pi)=-(-pi,dl|sg,-pi)
        //(dd|sd) rotate
        for (size_t iddd = 1; iddd < 16; ++iddd) {
          for (size_t idsd = 1; idsd < 6; ++idsd) {
            block(30 + iddd,10 + idsd) = intn[0]*DDrot(iddd,1)*SDrot(idsd,1) + intn[1]*(DDrot(iddd,2) + DDrot(iddd,3))*SDrot(idsd,1) + intn[2]*(DDrot(iddd,4) + DDrot(iddd,5))*SDrot(idsd,1);
            block(30 + iddd,10 + idsd) += intn[3]*(DDrot(iddd,6)*SDrot(idsd,2) + DDrot(iddd,7)*SDrot(idsd,3)) + intn[4]*(DDrot(iddd,8)*SDrot(idsd,4) + DDrot(iddd,9)*SDrot(idsd,5));
            block(30 + iddd,10 + idsd) += intn[5]*((DDrot(iddd,2) - DDrot(iddd,3))*SDrot(idsd,4) + DDrot(iddd,10)*SDrot(idsd,5));
            block(30 + iddd,10 + idsd) += intn[6]*((DDrot(iddd,11) + DDrot(iddd,14))*SDrot(idsd,2) + (DDrot(iddd,12) - DDrot(iddd,13))*SDrot(idsd,3));
          }
        }
        //(sd|dd) integrals
        Dd[0] = Dvalue(atmC,4);
        Dd[1] = Dvalue(atmD,6);
        intn[0] = d_eri2Center_dR(0,0,2,0,2,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR(0,0,2,0,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center_dR(0,0,2,0,2,2,2,2,RCD,Dd,atmC,atmD,type);    //(sg,sg|dl,dl)=(sg,sg|-dl,-dl)
        intn[3] = d_eri2Center_dR(0,0,2,1,2,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[4] = d_eri2Center_dR(0,0,2,2,2,0,2,2,RCD,Dd,atmC,atmD,type);    //(sg,dl|sg,dl)=(sg,-dl|sg,-dl)
        intn[5] = d_eri2Center_dR(0,0,2,2,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,dl|pi,pi)=-(sg,dl|-pi,-pi)=(sg,-dl|pi,-pi)
        intn[6] = d_eri2Center_dR(0,0,2,1,2,1,2,2,RCD,Dd,atmC,atmD,type);    //(sg,pi|pi,dl)=(sg,pi|-pi,-dl)=(sg,-pi|pi,-dl)=-(sg,-pi|-pi,dl)
        //(sd|dd) rotate
        for (size_t idsd = 1; idsd < 6; ++idsd) {
          for (size_t iddd = 1; iddd < 16; ++iddd) {
            block(10 + idsd,30 + iddd) = intn[0]*DDrot(iddd,1)*SDrot(idsd,1) + intn[1]*(DDrot(iddd,2) + DDrot(iddd,3))*SDrot(idsd,1) + intn[2]*(DDrot(iddd,4) + DDrot(iddd,5))*SDrot(idsd,1);
            block(10 + idsd,30 + iddd) += intn[3]*(DDrot(iddd,6)*SDrot(idsd,2) + DDrot(iddd,7)*SDrot(idsd,3)) + intn[4]*(DDrot(iddd,8)*SDrot(idsd,4) + DDrot(iddd,9)*SDrot(idsd,5));
            block(10 + idsd,30 + iddd) += intn[5]*((DDrot(iddd,2) - DDrot(iddd,3))*SDrot(idsd,4) + DDrot(iddd,10)*SDrot(idsd,5));
            block(10 + idsd,30 + iddd) += intn[6]*((DDrot(iddd,11) + DDrot(iddd,14))*SDrot(idsd,2) + (DDrot(iddd,12) - DDrot(iddd,13))*SDrot(idsd,3));
          }
        }
        //(pd|pd) integrals
        Dd[0] = Dvalue(atmC,5);
        Dd[1] = Dvalue(atmD,5);
        intn[0] = d_eri2Center_dR(1,0,2,0,1,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR(1,1,2,1,1,0,2,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)=(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center_dR(1,1,2,1,1,1,2,1,RCD,Dd,atmC,atmD,type);    //(pi,pi|pi,pi)=(-pi,-pi|pi,pi)=(pi,pi|-pi,-pi)=(-pi,-pi|-pi,-pi)
        intn[3] = d_eri2Center_dR(1,0,2,1,1,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[4] = d_eri2Center_dR(1,1,2,0,1,1,2,0,RCD,Dd,atmC,atmD,type);    //(pi,sg|pi,sg)=(-pi,sg|-pi,sg)
        intn[5] = d_eri2Center_dR(1,1,2,2,1,1,2,2,RCD,Dd,atmC,atmD,type);    //(pi,dl|pi,dl)=(pi,dl|-pi,-dl)=(-pi,dl|-pi,dl)=-(-pi,dl|pi,-dl)=-(pi,-dl|-pi,dl)=(pi,-dl|pi,-dl)=(-pi,-dl|pi,dl)=(-pi,-dl|-pi,-dl)
        intn[6] = d_eri2Center_dR(1,1,2,0,1,0,2,1,RCD,Dd,atmC,atmD,type);    //(pi,sg|sg,pi)=(-pi,sg|sg,-pi)=(sg,pi|pi,sg)=(sg,-pi|-pi,sg)
        intn[7] = d_eri2Center_dR(1,1,2,2,1,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|pi,dl)=(sg,pi|-pi,-dl)=-(sg,-pi|-pi,dl)=(sg,-pi|pi,-dl)=(pi,dl|sg,pi)=-(-pi,dl|sg,-pi)=(pi,-dl|sg,-pi)=(-pi,-dl|sg,pi)
        intn[8] = d_eri2Center_dR(1,1,2,2,1,1,2,0,RCD,Dd,atmC,atmD,type);    //(pi,sg|pi,dl)=(pi,sg|-pi,-dl)=-(-pi,sg|-pi,dl)=(-pi,sg|pi,-dl)=(pi,dl|pi,sg)=-(-pi,dl|-pi,sg)=(pi,-dl|-pi,sg)=(-pi,-dl|pi,sg)
        //(pd|pd) rotate
        for (size_t idpdb1 = 1; idpdb1 < 6; ++idpdb1) {            //d orbital on bra
          for (size_t idpdb2 = 1; idpdb2 < 4; ++idpdb2) {          //p orbital on bra
            for (size_t idpdk1 = 1; idpdk1 < 6; ++idpdk1) {        //d orbital on ket
              for (size_t idpdk2 = 1; idpdk2 < 4; ++idpdk2) {      //p orbital on ket
                cnt1 = (idpdb1 - 1)*3 + idpdb2;
                cnt2 = (idpdk1 - 1)*3 + idpdk2;
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) = intn[0]*PDrot(cnt1,1)*PDrot(cnt2,1) + intn[1]*((PDrot(cnt1,2) + PDrot(cnt1,3))*PDrot(cnt2,1) + PDrot(cnt1,1)*(PDrot(cnt2,2) + PDrot(cnt2,3)));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[2]*(PDrot(cnt1,2) + PDrot(cnt1,3))*(PDrot(cnt2,2) + PDrot(cnt2,3)) + intn[3]*(PDrot(cnt1,4)*PDrot(cnt2,4) + PDrot(cnt1,5)*PDrot(cnt2,5));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[4]*(PDrot(cnt1,8)*PDrot(cnt2,8) + PDrot(cnt1,12)*PDrot(cnt2,12));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[5]*((PDrot(cnt1,10) + PDrot(cnt1,15))*(PDrot(cnt2,10) + PDrot(cnt2,15)) + (PDrot(cnt1,14) - PDrot(cnt1,11))*(PDrot(cnt2,14) - PDrot(cnt2,11)));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[6]*(PDrot(cnt1,4)*PDrot(cnt2,8) + PDrot(cnt1,8)*PDrot(cnt2,4) + PDrot(cnt1,5)*PDrot(cnt2,12) + PDrot(cnt1,12)*PDrot(cnt2,5));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[7]*((PDrot(cnt1,10) + PDrot(cnt1,15))*PDrot(cnt2,4) + PDrot(cnt1,4)*(PDrot(cnt2,10) + PDrot(cnt2,15)) + (PDrot(cnt1,11) - PDrot(cnt1,14))*PDrot(cnt2,5) + PDrot(cnt1,5)*(PDrot(cnt2,11) - PDrot(cnt2,14)));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[8]*((PDrot(cnt1,10) + PDrot(cnt1,15))*PDrot(cnt2,8) + PDrot(cnt1,8)*(PDrot(cnt2,10) + PDrot(cnt2,15)) + (PDrot(cnt1,11) - PDrot(cnt1,14))*PDrot(cnt2,12) + PDrot(cnt1,12)*(PDrot(cnt2,11) - PDrot(cnt2,14)));
              }
            }
          }
        }
        //(pd|dd) integrals
        Dd[0] = Dvalue(atmC,5);
        Dd[1] = Dvalue(atmD,6);
        intn[0] = d_eri2Center_dR(1,0,2,0,2,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR(1,0,2,0,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center_dR(1,1,2,1,2,0,2,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[3] = d_eri2Center_dR(1,1,2,1,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(pi,pi|pi,pi)=(-pi,-pi|-pi,-pi)=(pi,pi|-pi,-pi)=(-pi,-pi|pi,pi)
        intn[4] = d_eri2Center_dR(1,0,2,0,2,2,2,2,RCD,Dd,atmC,atmD,type);    //(sg,sg|dl,dl)=(sg,sg|-dl,-dl)
        intn[5] = d_eri2Center_dR(1,1,2,1,2,2,2,2,RCD,Dd,atmC,atmD,type);    //(pi,pi|dl,dl)=(-pi,-pi|-dl,-dl)=(pi,pi|-dl,-dl)=(-pi,-pi|dl,dl)
        intn[6] = d_eri2Center_dR(1,0,2,1,2,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[7] = d_eri2Center_dR(1,1,2,0,2,1,2,0,RCD,Dd,atmC,atmD,type);    //(pi,sg|pi,sg)=(-pi,sg|-pi,sg)
        intn[8] = d_eri2Center_dR(1,1,2,2,2,1,2,2,RCD,Dd,atmC,atmD,type);    //(pi,dl|pi,dl)=(pi,dl|-pi,-dl)=(-pi,-dl|-pi,-dl)=(-pi,-dl|pi,dl)=(-pi,dl|-pi,dl)=-(-pi,dl|pi,-dl)=(pi,-dl|pi,-dl)=-(pi,-dl|-pi,dl)
        intn[9] = d_eri2Center_dR(1,1,2,0,2,1,2,2,RCD,Dd,atmC,atmD,type);    //(pi,sg|pi,dl)=(pi,sg|-pi,-dl)=-(-pi,sg|-pi,dl)=(-pi,sg|pi,-dl)
        intn[10] = d_eri2Center_dR(1,0,2,1,2,1,2,2,RCD,Dd,atmC,atmD,type);   //(sg,pi|pi,dl)=(sg,pi|-pi,-dl)=-(sg,-pi|-pi,dl)=(sg,-pi|pi,-dl)
        intn[11] = d_eri2Center_dR(1,1,2,2,2,0,2,1,RCD,Dd,atmC,atmD,type);   //(pi,dl|sg,pi)=(-pi,-dl|sg,pi)=-(-pi,dl|sg,-pi)=(pi,-dl|sg,-pi)
        //(pd|dd) rotate
        for (size_t iddd = 1; iddd < 16; ++iddd) {
          for (size_t idpd1 = 1; idpd1 < 6; ++idpd1) {
            for (size_t idpd2 = 1; idpd2 < 4; ++idpd2) {
              cnt1 = (idpd1 - 1)*3 + idpd2;
              block(5*(idpd2+2) + idpd1,30 + iddd) = -intn[0]*PDrot(cnt1,1)*DDrot(iddd,1) - intn[1]*PDrot(cnt1,1)*(DDrot(iddd,2) + DDrot(iddd,3));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[2]*(PDrot(cnt1,2) + PDrot(cnt1,3))*DDrot(iddd,1) + intn[3]*(PDrot(cnt1,2) + PDrot(cnt1,3))*(DDrot(iddd,2) + DDrot(iddd,3));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[4]*PDrot(cnt1,1)*(DDrot(iddd,4) + DDrot(iddd,5)) + intn[5]*(PDrot(cnt1,2) + PDrot(cnt1,3))*(DDrot(iddd,4) + DDrot(iddd,5));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[6]*(PDrot(cnt1,4)*DDrot(iddd,6) + PDrot(cnt1,5)*DDrot(iddd,7)) + intn[7]*(PDrot(cnt1,8)*DDrot(iddd,6) + PDrot(cnt1,12)*DDrot(iddd,7));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[8]*((PDrot(cnt1,10) + PDrot(cnt1,15))*(DDrot(iddd,11) + DDrot(iddd,14)) + (PDrot(cnt1,14) - PDrot(cnt1,11))*(DDrot(iddd,13) - DDrot(iddd,12)));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[9]*(PDrot(cnt1,8)*(DDrot(iddd,11) + DDrot(iddd,14)) + PDrot(cnt1,12)*(DDrot(iddd,12) - DDrot(iddd,13)));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[10]*(PDrot(cnt1,4)*(DDrot(iddd,11) + DDrot(iddd,14)) + PDrot(cnt1,5)*(DDrot(iddd,12) - DDrot(iddd,13)));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[11]*((PDrot(cnt1,10) + PDrot(cnt1,15))*DDrot(iddd,6) + (PDrot(cnt1,11) - PDrot(cnt1,14))*DDrot(iddd,7));
            }
          }
        }
        //(dd|pd) integrals
        Dd[0] = Dvalue(atmC,6);
        Dd[1] = Dvalue(atmD,5);
        intn[0] = d_eri2Center_dR(2,0,2,0,1,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR(2,1,2,1,1,0,2,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[2] = d_eri2Center_dR(2,0,2,0,1,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[3] = d_eri2Center_dR(2,1,2,1,1,1,2,1,RCD,Dd,atmC,atmD,type);    //(pi,pi|pi,pi)=(-pi,-pi|-pi,-pi)=(pi,pi|-pi,-pi)=(-pi,-pi|pi,pi)
        intn[4] = d_eri2Center_dR(2,2,2,2,1,0,2,0,RCD,Dd,atmC,atmD,type);    //(dl,dl|sg,sg)=(-dl,-dl|sg,sg)
        intn[5] = d_eri2Center_dR(2,2,2,2,1,1,2,1,RCD,Dd,atmC,atmD,type);    //(dl,dl|pi,pi)=(-dl,-dl|-pi,-pi)=(-dl,-dl|pi,pi)=(dl,dl|-pi,-pi)
        intn[6] = d_eri2Center_dR(2,0,2,1,1,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[7] = d_eri2Center_dR(2,1,2,0,1,1,2,0,RCD,Dd,atmC,atmD,type);    //(pi,sg|pi,sg)=(-pi,sg|-pi,sg)
        intn[8] = d_eri2Center_dR(2,1,2,2,1,1,2,2,RCD,Dd,atmC,atmD,type);    //(pi,dl|pi,dl)=(-pi,-dl|-pi,-dl)=(-pi,-dl|pi,dl)=(pi,dl|-pi,-dl)=(-pi,dl|-pi,dl)=(pi,-dl|pi,-dl)=-(-pi,dl|pi,-dl)=-(pi,-dl|-pi,dl)
        intn[9] = d_eri2Center_dR(2,1,2,2,1,1,2,0,RCD,Dd,atmC,atmD,type);    //(pi,dl|pi,sg)=-(-pi,dl|-pi,sg)=(pi,-dl|-pi,sg)=(-pi,-dl|pi,sg)
        intn[10] = d_eri2Center_dR(2,1,2,2,1,0,2,1,RCD,Dd,atmC,atmD,type);   //(pi,dl|sg,pi)=-(-pi,dl|sg,-pi)=(pi,-dl|sg,-pi)=(-pi,-dl|sg,pi)
        intn[11] = d_eri2Center_dR(2,0,2,1,1,1,2,2,RCD,Dd,atmC,atmD,type);   //(sg,pi|pi,dl)=(sg,pi|-pi,-dl)=-(sg,-pi|-pi,dl)=(sg,-pi|pi,-dl)
        //(dd|pd) rotate
        for (size_t iddd = 1; iddd < 16; ++iddd) {
          for (size_t idpd1 = 1; idpd1 < 6; ++idpd1) {
            for (size_t idpd2 = 1; idpd2 < 4; ++idpd2) {
              cnt1 = (idpd1 - 1)*3 + idpd2;
              block(30 + iddd,5*(idpd2+2) + idpd1) = -intn[0]*PDrot(cnt1,1)*DDrot(iddd,1) - intn[1]*PDrot(cnt1,1)*(DDrot(iddd,2) + DDrot(iddd,3));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[2]*(PDrot(cnt1,2) + PDrot(cnt1,3))*DDrot(iddd,1) + intn[3]*(PDrot(cnt1,2) + PDrot(cnt1,3))*(DDrot(iddd,2) + DDrot(iddd,3));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[4]*PDrot(cnt1,1)*(DDrot(iddd,4) + DDrot(iddd,5)) + intn[5]*(PDrot(cnt1,2) + PDrot(cnt1,3))*(DDrot(iddd,4) + DDrot(iddd,5));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[6]*(PDrot(cnt1,4)*DDrot(iddd,6) + PDrot(cnt1,5)*DDrot(iddd,7)) + intn[7]*(PDrot(cnt1,8)*DDrot(iddd,6) + PDrot(cnt1,12)*DDrot(iddd,7));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[8]*((PDrot(cnt1,10) + PDrot(cnt1,15))*(DDrot(iddd,11) + DDrot(iddd,14)) + (PDrot(cnt1,14) - PDrot(cnt1,11))*(DDrot(iddd,13) - DDrot(iddd,12)));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[9]*(PDrot(cnt1,8)*(DDrot(iddd,11) + DDrot(iddd,14)) + PDrot(cnt1,12)*(DDrot(iddd,12) - DDrot(iddd,13)));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[10]*(PDrot(cnt1,4)*(DDrot(iddd,11) + DDrot(iddd,14)) + PDrot(cnt1,5)*(DDrot(iddd,12) - DDrot(iddd,13)));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[11]*((PDrot(cnt1,10) + PDrot(cnt1,15))*DDrot(iddd,6) + (PDrot(cnt1,11) - PDrot(cnt1,14))*DDrot(iddd,7));
            }
          }
        }
        //(dd|dd) integrals
        Dd[0] = Dvalue(atmC,6);
        Dd[1] = Dvalue(atmD,6);
        intn[0] = d_eri2Center_dR(2,0,2,0,2,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR(2,0,2,0,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center_dR(2,0,2,0,2,2,2,2,RCD,Dd,atmC,atmD,type);    //(sg,sg|dl,dl)=(sg,sg|-dl,-dl)
        intn[3] = d_eri2Center_dR(2,1,2,1,2,0,2,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[4] = d_eri2Center_dR(2,2,2,2,2,0,2,0,RCD,Dd,atmC,atmD,type);    //(dl,dl|sg,sg)=(-dl,-dl|sg,sg)
        intn[5] = d_eri2Center_dR(2,2,2,2,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(dl,dl|pi,pi)=(-dl,-dl|pi,pi)=(dl,dl|-pi,-pi)=(-dl,-dl|-pi,-pi)
        intn[6] = d_eri2Center_dR(2,1,2,1,2,2,2,2,RCD,Dd,atmC,atmD,type);    //(pi,pi|dl,dl)=(-pi,-pi|dl,dl)=(pi,pi|-dl,-dl)=(-pi,-pi|-dl,-dl)
        intn[7] = d_eri2Center_dR(2,1,2,1,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(pi,pi|pi,pi)=(-pi,-pi|-pi,-pi)
        intn[8] = d_eri2Center_dR(2,2,2,2,2,2,2,2,RCD,Dd,atmC,atmD,type);    //(dl,dl|dl,dl)=(-dl,-dl|-dl,-dl)
        intn[9] = d_eri2Center_dR(2,1,2,1,2,-1,2,-1,RCD,Dd,atmC,atmD,type);  //(pi,pi|-pi,-pi)=(-pi,-pi|pi,pi)
        intn[10] = d_eri2Center_dR(2,-2,2,-2,2,2,2,2,RCD,Dd,atmC,atmD,type); //(-dl,-dl|dl,dl)=(dl,dl|-dl,-dl)
        intn[11] = d_eri2Center_dR(2,0,2,1,2,0,2,1,RCD,Dd,atmC,atmD,type);   //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[12] = d_eri2Center_dR(2,0,2,2,2,0,2,2,RCD,Dd,atmC,atmD,type);   //(sg,dl|sg,dl)=(sg,-dl|sg,-dl)
        intn[13] = d_eri2Center_dR(2,1,2,2,2,1,2,2,RCD,Dd,atmC,atmD,type);   //(pi,dl|pi,dl)=(pi,dl|-pi,-dl)=(-pi,-dl|pi,dl)=(-pi,-dl|-pi,-dl)=-(-pi,dl|pi,-dl)=-(pi,-dl|-pi,dl)=(-pi,dl|-pi,dl)=(pi,-dl|pi,-dl)
        intn[14] = d_eri2Center_dR(2,1,2,-1,2,1,2,-1,RCD,Dd,atmC,atmD,type); //(pi,-pi|pi,-pi)
        intn[15] = d_eri2Center_dR(2,1,2,-2,2,0,2,-1,RCD,Dd,atmC,atmD,type); //(pi,dl|sg,pi)=(-pi,-dl|sg,pi)=(pi,-dl|sg,-pi)=-(-pi,dl|sg,-pi)=(sg,pi|pi,dl)=(sg,pi|-pi,-dl)=(sg,-pi|pi,-dl)=-(sg,-pi|-pi,dl)
        intn[16] = d_eri2Center_dR(2,1,2,1,2,0,2,2,RCD,Dd,atmC,atmD,type);   //(pi,pi|sg,dl)=-(-pi,-pi|sg,dl)=(pi,-pi|sg,-dl)=(sg,dl|pi,pi)=-(sg,dl|-pi,-pi)=(sg,-dl|pi,-pi)
        //(dd|dd) rotate
        for (size_t iddd1 = 1; iddd1 < 16; ++iddd1) {
          for (size_t iddd2 = 1; iddd2 < 16; ++iddd2) {
            block(30 + iddd1,30 + iddd2) = intn[0]*DDrot(iddd1,1)*DDrot(iddd2,1) + intn[1]*DDrot(iddd1,1)*(DDrot(iddd2,2) + DDrot(iddd2,3)) + intn[2]*DDrot(iddd1,1)*(DDrot(iddd2,4) + DDrot(iddd2,5));
            block(30 + iddd1,30 + iddd2) += intn[3]*(DDrot(iddd1,2) + DDrot(iddd1,3))*DDrot(iddd2,1) + intn[4]*(DDrot(iddd1,4) + DDrot(iddd1,5))*DDrot(iddd2,1);
            block(30 + iddd1,30 + iddd2) += intn[5]*(DDrot(iddd1,4) + DDrot(iddd1,5))*(DDrot(iddd2,2) + DDrot(iddd2,3)) + intn[6]*(DDrot(iddd1,2) + DDrot(iddd1,3))*(DDrot(iddd2,4) + DDrot(iddd2,5));
            block(30 + iddd1,30 + iddd2) += intn[7]*(DDrot(iddd1,2)*DDrot(iddd2,2) + DDrot(iddd1,3)*DDrot(iddd2,3)) + intn[8]*(DDrot(iddd1,4)*DDrot(iddd2,4) + DDrot(iddd1,5)*DDrot(iddd2,5));
            block(30 + iddd1,30 + iddd2) += intn[9]*(DDrot(iddd1,2)*DDrot(iddd2,3) + DDrot(iddd1,3)*DDrot(iddd2,2)) + intn[10]*(DDrot(iddd1,4)*DDrot(iddd2,5) + DDrot(iddd1,5)*DDrot(iddd2,4));
            block(30 + iddd1,30 + iddd2) += intn[11]*(DDrot(iddd1,6)*DDrot(iddd2,6) + DDrot(iddd1,7)*DDrot(iddd2,7)) + intn[12]*(DDrot(iddd1,8)*DDrot(iddd2,8) + DDrot(iddd1,9)*DDrot(iddd2,9));
            block(30 + iddd1,30 + iddd2) += intn[13]*((DDrot(iddd1,11) + DDrot(iddd1,14))*(DDrot(iddd2,11) + DDrot(iddd2,14)) + (DDrot(iddd1,12) - DDrot(iddd1,13))*(DDrot(iddd2,12) - DDrot(iddd2,13)));
            block(30 + iddd1,30 + iddd2) += intn[14]*DDrot(iddd1,10)*DDrot(iddd2,10);
            block(30 + iddd1,30 + iddd2) += intn[15]*((DDrot(iddd1,11) + DDrot(iddd1,14))*DDrot(iddd2,6) + (DDrot(iddd1,12) - DDrot(iddd1,13))*DDrot(iddd2,7) + DDrot(iddd1,6)*(DDrot(iddd2,11) + DDrot(iddd2,14)) + DDrot(iddd1,7)*(DDrot(iddd2,12) - DDrot(iddd2,13)));
            block(30 + iddd1,30 + iddd2) += intn[16]*((DDrot(iddd1,2) - DDrot(iddd1,3))*DDrot(iddd2,8) + DDrot(iddd1,10)*DDrot(iddd2,9) + DDrot(iddd1,8)*(DDrot(iddd2,2) - DDrot(iddd2,3)) + DDrot(iddd1,9)*DDrot(iddd2,10));
          }
        }
      }
    }
  }
  void IntegralBlock2C_dA(int type, matrixE & block, int der, int atmC, int atmD, double RCD, double cost, double sint, double cosp, double sinp) {
    //function calculating a block of first-derivatives of two-center integrals with respect to angles
    //2 -> theta; 3 -> phi
    matrixE SProt_dA(1,1);
    matrixE PProt_dA(1,1);
    matrixE SDrot_dA(1,1);
    matrixE PDrot_dA(1,1);
    matrixE DDrot_dA(1,1);
    int cnt1;
    int cnt2;
    size_t nrows = block.rows();
    size_t ncols = block.cols();
    if ((nrows > 1)||(ncols > 1)) {
      SProt = SPtransf(cost,sint,cosp,sinp);
      PProt = PPtransf(cost,sint,cosp,sinp);
      if (der == 2) {
        SProt_dA = SPtransf_dt(cost,sint,cosp,sinp);
        PProt_dA = PPtransf_dt(cost,sint,cosp,sinp);
      }
      else if (der == 3) {
        SProt_dA = SPtransf_dp(cost,sint,cosp,sinp);
        PProt_dA = PPtransf_dp(cost,sint,cosp,sinp);
      }
      if ((nrows > 10)||(ncols > 10)) {
        SDrot = SDtransf(cost,sint,cosp,sinp);
        PDrot = PDtransf(cost,sint,cosp,sinp);
        DDrot = DDtransf(cost,sint,cosp,sinp);
        if (der == 2) {
          SDrot_dA = SDtransf_dt(cost,sint,cosp,sinp);
          PDrot_dA = PDtransf_dt(cost,sint,cosp,sinp);
          DDrot_dA = DDtransf_dt(cost,sint,cosp,sinp);
        }
        else if (der == 3) {
          SDrot_dA = SDtransf_dp(cost,sint,cosp,sinp);
          PDrot_dA = PDtransf_dp(cost,sint,cosp,sinp);
          DDrot_dA = DDtransf_dp(cost,sint,cosp,sinp);
        }
      }
    }
    D[0] = Dvalue(atmC,1);
    D[1] = Dvalue(atmC,2);
    D[2] = Dvalue(atmD,1);
    D[3] = Dvalue(atmD,2);
    //integral calculation
    block(1,1) = 0.0;                                                                                                                   //(ss|ss)
    if (nrows > 1) {
      intn[0] = eri2Center(0,0,1,0,0,0,0,0,RCD,D,atmC,atmD,type);
      intn[1] = eri2Center(1,0,1,0,0,0,0,0,RCD,D,atmC,atmD,type);
      intn[2] = eri2Center(1,1,1,1,0,0,0,0,RCD,D,atmC,atmD,type);
      for (size_t idx = 0; idx < 6; ++idx) {
        if (idx < 3) {block(idx + 2,1) = -intn[0]*SProt_dA(idx + 1,1);}                                                         //(sp|ss) integrals
        block(5 + idx,1) = intn[1]*PProt_dA(idx + 1,1) + intn[2]*(PProt_dA(idx + 1,2) + PProt_dA(idx + 1,3));   //(pp|ss) integrals
      }
      if (nrows > 10) {
        //(sd|ss) integrals
        Dd[0] = Dvalue(atmC,4);
        Dd[1] = Dvalue(atmD,1);
        intn[0] = d_eri2Center(0,0,2,0,0,0,0,0,RCD,Dd,atmC,atmD,type);
        for (size_t idsd = 1; idsd < 6; ++idsd) {
          block(10 + idsd,1) = intn[0]*SDrot_dA(idsd,1);
        }
        //(pd|ss) integrals
        Dd[0] = Dvalue(atmC,5);
        intn[0] = d_eri2Center(1,0,2,0,0,0,0,0,RCD,Dd,atmC,atmD,type);
        intn[1] = d_eri2Center(1,1,2,1,0,0,0,0,RCD,Dd,atmC,atmD,type);
        for (int idpd1 = 1; idpd1 < 6; ++idpd1) {
          for (int idpd2 = 1; idpd2 < 4; ++idpd2) {
            block(5*(idpd2+2) + idpd1,1) = -(intn[0]*PDrot_dA((idpd1 - 1)*3 + idpd2,1) + intn[1]*(PDrot_dA((idpd1 - 1)*3 + idpd2,2) + PDrot_dA((idpd1 - 1)*3 + idpd2,3)));
          }
        }
        //(dd|ss) integrals
        Dd[0] = Dvalue(atmC,6);
        intn[0] = d_eri2Center(2,0,2,0,0,0,0,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(2,1,2,1,0,0,0,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[2] = d_eri2Center(2,2,2,2,0,0,0,0,RCD,Dd,atmC,atmD,type);    //(dl,dl|sg,sg)=(-dl,-dl|sg,sg)
        for (size_t iddd = 1; iddd < 16; ++iddd) {
          block(30 + iddd,1) = intn[0]*DDrot_dA(iddd,1) + intn[1]*(DDrot_dA(iddd,2) + DDrot_dA(iddd,3)) + intn[2]*(DDrot_dA(iddd,4) + DDrot_dA(iddd,5));
        }
      }
    }
    if ((ncols > 1)&&(type == 0)) {
      intn[0] = eri2Center(0,0,0,0,0,0,1,0,RCD,D,atmC,atmD,type);
      intn[1] = eri2Center(0,0,0,0,1,0,1,0,RCD,D,atmC,atmD,type);
      intn[2] = eri2Center(0,0,0,0,1,1,1,1,RCD,D,atmC,atmD,type);
      for (size_t idx = 0; idx < 6; ++idx) {
        if (idx < 3) {block(1,idx + 2) = -intn[0]*SProt_dA(idx + 1,1);}                                                         //(ss|sp) integrals
        block(1,5 + idx) = intn[1]*PProt_dA(idx + 1,1) + intn[2]*(PProt_dA(idx + 1,2) + PProt_dA(idx + 1,3));   //(ss|pp) integrals
      }
      if (ncols > 10) {
        //(ss|sd) integrals
        Dd[0] = Dvalue(atmC,1);
        Dd[1] = Dvalue(atmD,4);
        intn[0] = d_eri2Center(0,0,0,0,0,0,2,0,RCD,Dd,atmC,atmD,type);
        for (size_t idsd = 1; idsd < 6; ++idsd) {
          block(1,10 + idsd) = intn[0]*SDrot_dA(idsd,1);
        }
        //(ss|pd) integrals
        Dd[1] = Dvalue(atmD,5);
        intn[0] = d_eri2Center(0,0,0,0,1,0,2,0,RCD,Dd,atmC,atmD,type);
        intn[1] = d_eri2Center(0,0,0,0,1,1,2,1,RCD,Dd,atmC,atmD,type);
        for (int idpd1 = 1; idpd1 < 6; ++idpd1) {
          for (int idpd2 = 1; idpd2 < 4; ++idpd2) {
            block(1,5*(idpd2 + 2) + idpd1) = -intn[0]*PDrot_dA((idpd1 - 1)*3 + idpd2,1) - intn[1]*(PDrot_dA((idpd1 - 1)*3 + idpd2,2) + PDrot_dA((idpd1 - 1)*3 + idpd2,3));
          }
        }
        //(ss|dd) integrals
        Dd[1] = Dvalue(atmD,6);
        intn[0] = d_eri2Center(0,0,0,0,2,0,2,0,RCD,Dd,atmC,atmD,type);
        intn[1] = d_eri2Center(0,0,0,0,2,1,2,1,RCD,Dd,atmC,atmD,type);
        intn[2] = d_eri2Center(0,0,0,0,2,2,2,2,RCD,Dd,atmC,atmD,type);
        for (size_t iddd = 1; iddd < 16; ++iddd) {
          block(1,30 + iddd) = intn[0]*DDrot_dA(iddd,1) + intn[1]*(DDrot_dA(iddd,2) + DDrot_dA(iddd,3)) + intn[2]*(DDrot_dA(iddd,4) + DDrot_dA(iddd,5));
        }
      }
    }
    if ((nrows > 1)&&(ncols > 1)&&(type == 0)) {
      intn[0] = eri2Center(0,0,1,0,0,0,1,0,RCD,D,atmC,atmD,type);
      intn[1] = eri2Center(0,0,1,1,0,0,1,1,RCD,D,atmC,atmD,type);
      intn[2] = eri2Center(0,0,1,0,1,0,1,0,RCD,D,atmC,atmD,type);
      intn[3] = eri2Center(0,0,1,0,1,1,1,1,RCD,D,atmC,atmD,type);
      intn[4] = eri2Center(0,0,1,1,1,0,1,1,RCD,D,atmC,atmD,type);
      for (size_t idpb = 1; idpb < 4; ++idpb) {
        for (size_t idpk = 1; idpk < 4; ++idpk) {           //(sp|sp)
          block(idpb + 1,idpk + 1) = intn[0]*SProt_dA(idpb,1)*SProt(idpk,1) + intn[1]*(SProt_dA(idpb,2)*SProt(idpk,2) + SProt_dA(idpb,3)*SProt(idpk,3));
          block(idpb + 1,idpk + 1) += intn[0]*SProt(idpb,1)*SProt_dA(idpk,1) + intn[1]*(SProt(idpb,2)*SProt_dA(idpk,2) + SProt(idpb,3)*SProt_dA(idpk,3));
        }
        for (size_t idc = 1; idc < 7; ++idc) {              //(sp|pp)
          block(idpb + 1,idc + 4) = -SProt_dA(idpb,1)*(intn[2]*PProt(idc,1) + intn[3]*(PProt(idc,2) + PProt(idc,3))) - intn[4]*(SProt_dA(idpb,2)*PProt(idc,4) + SProt_dA(idpb,3)*PProt(idc,5));
          block(idpb + 1,idc + 4) += -SProt(idpb,1)*(intn[2]*PProt_dA(idc,1) + intn[3]*(PProt_dA(idc,2) + PProt_dA(idc,3))) - intn[4]*(SProt(idpb,2)*PProt_dA(idc,4) + SProt(idpb,3)*PProt_dA(idc,5));
        }
      }
      intn[0] = eri2Center(1,0,1,0,1,0,1,0,RCD,D,atmC,atmD,type);                          //(ps,ps|ps,ps)
      intn[1] = eri2Center(1,1,1,1,1,0,1,0,RCD,D,atmC,atmD,type);                          //(pp,pp|ps,ps)
      intn[2] = eri2Center(1,0,1,0,1,1,1,1,RCD,D,atmC,atmD,type);                          //(ps,ps|pp,pp)
      intn[3] = eri2Center(1,1,1,1,1,1,1,1,RCD,D,atmC,atmD,type);                          //(pp,pp|pp,pp)
      intn[4] = eri2Center(1,1,1,0,1,1,1,0,RCD,D,atmC,atmD,type);                          //(pp,ps|pp,ps)
      intn[5] = eri2Center(1,1,1,1,1,-1,1,-1,RCD,D,atmC,atmD,type);                        //(pp,pp|pp*,pp*)
      intn[6] = eri2Center(1,-1,1,1,1,-1,1,1,RCD,D,atmC,atmD,type);                        //(pp*,pp|pp*,pp)
      intn[7] = eri2Center(1,0,1,0,0,0,1,0,RCD,D,atmC,atmD,type);
      intn[8] = eri2Center(1,1,1,1,0,0,1,0,RCD,D,atmC,atmD,type);
      intn[9] = eri2Center(1,0,1,1,0,0,1,1,RCD,D,atmC,atmD,type);
      for (size_t idr = 1; idr < 7; ++idr) {
        for (size_t idc = 1; idc < 4; ++idc) {                //(pp|sp)
          block(idr + 4,idc + 1) = -SProt_dA(idc,1)*(intn[7]*PProt(idr,1) + intn[8]*(PProt(idr,2) + PProt(idr,3))) - intn[9]*(SProt_dA(idc,2)*PProt(idr,4) + SProt_dA(idc,3)*PProt(idr,5));
          block(idr + 4,idc + 1) += -SProt(idc,1)*(intn[7]*PProt_dA(idr,1) + intn[8]*(PProt_dA(idr,2) + PProt_dA(idr,3))) - intn[9]*(SProt(idc,2)*PProt_dA(idr,4) + SProt(idc,3)*PProt_dA(idr,5));
        }
        for (size_t idc = 1; idc < 7; ++idc) {                //(pp|pp)
          block(idr + 4,idc + 4) = intn[0]*PProt_dA(idr,1)*PProt(idc,1) + intn[1]*(PProt_dA(idr,2) + PProt_dA(idr,3))*PProt(idc,1) + intn[2]*PProt_dA(idr,1)*(PProt(idc,2) + PProt(idc,3)) + intn[3]*(PProt_dA(idr,2)*PProt(idc,2) + PProt_dA(idr,3)*PProt(idc,3)) + intn[4]*(PProt_dA(idr,4)*PProt(idc,4) + PProt_dA(idr,5)*PProt(idc,5)) + intn[5]*(PProt_dA(idr,2)*PProt(idc,3) + PProt_dA(idr,3)*PProt(idc,2)) + intn[6]*PProt_dA(idr,6)*PProt(idc,6);
          block(idr + 4,idc + 4) += intn[0]*PProt(idr,1)*PProt_dA(idc,1) + intn[1]*(PProt(idr,2) + PProt(idr,3))*PProt_dA(idc,1) + intn[2]*PProt(idr,1)*(PProt_dA(idc,2) + PProt_dA(idc,3)) + intn[3]*(PProt(idr,2)*PProt_dA(idc,2) + PProt(idr,3)*PProt_dA(idc,3)) + intn[4]*(PProt(idr,4)*PProt_dA(idc,4) + PProt(idr,5)*PProt_dA(idc,5)) + intn[5]*(PProt(idr,2)*PProt_dA(idc,3) + PProt(idr,3)*PProt_dA(idc,2)) + intn[6]*PProt(idr,6)*PProt_dA(idc,6);
        }
      }
      if (nrows > 10) {
        //(sd|sp) integrals
        Dd[0] = Dvalue(atmC,4);
        Dd[1] = Dvalue(atmD,1);
        intn[0] = d_eri2Center(0,0,2,0,0,0,1,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(0,0,2,1,0,0,1,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        //(sd|sp) rotate
        for (size_t idsd = 1; idsd < 6; ++idsd) {
          for (size_t idsp = 1; idsp < 4; ++idsp) {
            block(10 + idsd,1 + idsp) = -intn[0]*SDrot_dA(idsd,1)*SProt(idsp,1) - intn[1]*(SDrot_dA(idsd,2)*SProt(idsp,2) + SDrot_dA(idsd,3)*SProt(idsp,3));
            block(10 + idsd,1 + idsp) += -intn[0]*SDrot(idsd,1)*SProt_dA(idsp,1) - intn[1]*(SDrot(idsd,2)*SProt_dA(idsp,2) + SDrot(idsd,3)*SProt_dA(idsp,3));
          }
        }
        //(sd|pp) integrals
        Dd[0] = Dvalue(atmC,4);
        Dd[1] = Dvalue(atmD,3);
        intn[0] = d_eri2Center(0,0,2,0,1,0,1,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(0,0,2,0,1,1,1,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center(0,0,2,1,1,0,1,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[3] = d_eri2Center(0,0,2,2,1,1,1,1,RCD,Dd,atmC,atmD,type);    //(sg,dl|pi,pi)=-(sg,dl|-pi,-pi)=(sg,-dl|pi,-pi)
        //(sd|pp) rotate
        for (size_t idsd = 1; idsd < 6; ++idsd) {
          for (size_t idpp = 1; idpp < 7; ++idpp) {
            block(10 + idsd,4 + idpp) = intn[0]*SDrot_dA(idsd,1)*PProt(idpp,1) + intn[1]*SDrot_dA(idsd,1)*(PProt(idpp,2) + PProt(idpp,3));
            block(10 + idsd,4 + idpp) += intn[0]*SDrot(idsd,1)*PProt_dA(idpp,1) + intn[1]*SDrot(idsd,1)*(PProt_dA(idpp,2) + PProt_dA(idpp,3));
            block(10 + idsd,4 + idpp) += intn[2]*(SDrot_dA(idsd,2)*PProt(idpp,4) + SDrot_dA(idsd,3)*PProt(idpp,5));
            block(10 + idsd,4 + idpp) += intn[2]*(SDrot(idsd,2)*PProt_dA(idpp,4) + SDrot(idsd,3)*PProt_dA(idpp,5));
            block(10 + idsd,4 + idpp) += intn[3]*(SDrot_dA(idsd,4)*PProt(idpp,2) - SDrot_dA(idsd,4)*PProt(idpp,3) + SDrot_dA(idsd,5)*PProt(idpp,6));
            block(10 + idsd,4 + idpp) += intn[3]*(SDrot(idsd,4)*PProt_dA(idpp,2) - SDrot(idsd,4)*PProt_dA(idpp,3) + SDrot(idsd,5)*PProt_dA(idpp,6));
          }
        }
        //(pd|sp) integrals
        Dd[0] = Dvalue(atmC,5);
        Dd[1] = Dvalue(atmD,1);
        intn[0] = d_eri2Center(1,0,2,0,0,0,1,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(1,1,2,1,0,0,1,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[2] = d_eri2Center(1,0,2,1,0,0,1,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[3] = d_eri2Center(1,1,2,0,0,0,1,1,RCD,Dd,atmC,atmD,type);    //(pi,sg|sg,pi)=(-pi,sg|sg,-pi)
        intn[4] = d_eri2Center(1,1,2,2,0,0,1,1,RCD,Dd,atmC,atmD,type);    //(pi,dl|sg,pi)=(-pi,-dl|sg,pi)=-(-pi,dl|sg,-pi)=(pi,-dl|sg,-pi)
        //(pd|sp) rotate
        for (size_t idsp = 1; idsp < 4; ++idsp) {
          for (size_t idpd1 = 1; idpd1 < 6; ++idpd1) {
            for (size_t idpd2 = 1; idpd2 < 4; ++idpd2) {
              block(5*(idpd2+2) + idpd1,1 + idsp) = intn[0]*PDrot_dA((idpd1 - 1)*3 + idpd2,1)*SProt(idsp,1) + intn[1]*(PDrot_dA((idpd1 - 1)*3 + idpd2,2) + PDrot_dA((idpd1 - 1)*3 + idpd2,3))*SProt(idsp,1);
              block(5*(idpd2+2) + idpd1,1 + idsp) += intn[0]*PDrot((idpd1 - 1)*3 + idpd2,1)*SProt_dA(idsp,1) + intn[1]*(PDrot((idpd1 - 1)*3 + idpd2,2) + PDrot((idpd1 - 1)*3 + idpd2,3))*SProt_dA(idsp,1);
              block(5*(idpd2+2) + idpd1,1 + idsp) += intn[2]*(PDrot_dA((idpd1 - 1)*3 + idpd2,4)*SProt(idsp,2) + PDrot_dA((idpd1 - 1)*3 + idpd2,5)*SProt(idsp,3));
              block(5*(idpd2+2) + idpd1,1 + idsp) += intn[2]*(PDrot((idpd1 - 1)*3 + idpd2,4)*SProt_dA(idsp,2) + PDrot((idpd1 - 1)*3 + idpd2,5)*SProt_dA(idsp,3));
              block(5*(idpd2+2) + idpd1,1 + idsp) += intn[3]*(PDrot_dA((idpd1 - 1)*3 + idpd2,8)*SProt(idsp,2) + PDrot_dA((idpd1 - 1)*3 + idpd2,12)*SProt(idsp,3));
              block(5*(idpd2+2) + idpd1,1 + idsp) += intn[3]*(PDrot((idpd1 - 1)*3 + idpd2,8)*SProt_dA(idsp,2) + PDrot((idpd1 - 1)*3 + idpd2,12)*SProt_dA(idsp,3));
              block(5*(idpd2+2) + idpd1,1 + idsp) += intn[4]*((PDrot_dA((idpd1 - 1)*3 + idpd2,10) + PDrot_dA((idpd1 - 1)*3 + idpd2,15))*SProt(idsp,2) + (PDrot_dA((idpd1 - 1)*3 + idpd2,11) - PDrot_dA((idpd1 - 1)*3 + idpd2,14))*SProt(idsp,3));
              block(5*(idpd2+2) + idpd1,1 + idsp) += intn[4]*((PDrot((idpd1 - 1)*3 + idpd2,10) + PDrot((idpd1 - 1)*3 + idpd2,15))*SProt_dA(idsp,2) + (PDrot((idpd1 - 1)*3 + idpd2,11) - PDrot((idpd1 - 1)*3 + idpd2,14))*SProt_dA(idsp,3));
            }
          }
        }
        //(dd|sp) integrals
        Dd[0] = Dvalue(atmC,6);
        Dd[1] = Dvalue(atmD,1);
        intn[0] = d_eri2Center(2,0,2,0,0,0,1,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(2,1,2,1,0,0,1,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[2] = d_eri2Center(2,2,2,2,0,0,1,0,RCD,Dd,atmC,atmD,type);    //(dl,dl|sg,sg)=(-dl,-dl|sg,sg)
        intn[3] = d_eri2Center(2,0,2,1,0,0,1,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[4] = d_eri2Center(2,1,2,2,0,0,1,1,RCD,Dd,atmC,atmD,type);    //(pi,dl|sg,pi)=(-pi,-dl|sg,pi)=-(-pi,dl|sg,-pi)=(pi,-dl|sg,-pi)
        //(dd|sp) rotate
        for (size_t idsp = 1; idsp < 4; ++idsp) {
          for (size_t iddd = 1; iddd < 16; ++iddd) {
            block(30 + iddd,1 + idsp) = -intn[0]*DDrot_dA(iddd,1)*SProt(idsp,1) - intn[1]*(DDrot_dA(iddd,2) + DDrot_dA(iddd,3))*SProt(idsp,1) - intn[2]*(DDrot_dA(iddd,4) + DDrot_dA(iddd,5))*SProt(idsp,1);
            block(30 + iddd,1 + idsp) += -intn[0]*DDrot(iddd,1)*SProt_dA(idsp,1) - intn[1]*(DDrot(iddd,2) + DDrot(iddd,3))*SProt_dA(idsp,1) - intn[2]*(DDrot(iddd,4) + DDrot(iddd,5))*SProt_dA(idsp,1);
            block(30 + iddd,1 + idsp) -= intn[3]*(DDrot_dA(iddd,6)*SProt(idsp,2) + DDrot_dA(iddd,7)*SProt(idsp,3));
            block(30 + iddd,1 + idsp) -= intn[3]*(DDrot(iddd,6)*SProt_dA(idsp,2) + DDrot(iddd,7)*SProt_dA(idsp,3));
            block(30 + iddd,1 + idsp) -= intn[4]*((DDrot_dA(iddd,11) + DDrot_dA(iddd,14))*SProt(idsp,2) + (DDrot_dA(iddd,12) - DDrot_dA(iddd,13))*SProt(idsp,3));
            block(30 + iddd,1 + idsp) -= intn[4]*((DDrot(iddd,11) + DDrot(iddd,14))*SProt_dA(idsp,2) + (DDrot(iddd,12) - DDrot(iddd,13))*SProt_dA(idsp,3));
          }
        }
        //(pd|pp) integrals
        Dd[0] = Dvalue(atmC,5);
        Dd[1] = Dvalue(atmD,3);
        intn[0] = d_eri2Center(1,0,2,0,1,0,1,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(1,1,2,1,1,0,1,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[2] = d_eri2Center(1,0,2,0,1,1,1,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[3] = d_eri2Center(1,1,2,1,1,1,1,1,RCD,Dd,atmC,atmD,type);    //(pi,pi|pi,pi)=(pi,pi|-pi,-pi)=(-pi,-pi|pi,pi)=(-pi,-pi|-pi,-pi)
        intn[4] = d_eri2Center(1,0,2,1,1,0,1,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[5] = d_eri2Center(1,1,2,0,1,1,1,0,RCD,Dd,atmC,atmD,type);    //(pi,sg|pi,sg)=(-pi,sg|-pi,sg)
        intn[6] = d_eri2Center(1,1,2,2,1,0,1,1,RCD,Dd,atmC,atmD,type);    //(pi,dl|sg,pi)=(-pi,-dl|sg,pi)=-(-pi,dl|sg,-pi)=(pi,-dl|sg,-pi)
        //(pd|pp) rotate
        for (size_t idpp = 1; idpp < 7; ++idpp) {
          for (size_t idpd1 = 1; idpd1 < 6; ++idpd1) {
            for (size_t idpd2 = 1; idpd2 < 4; ++idpd2) {
              block(5*(idpd2+2) + idpd1,4 + idpp) = -intn[0]*PDrot_dA((idpd1 - 1)*3 + idpd2,1)*PProt(idpp,1) - intn[1]*(PDrot_dA((idpd1 - 1)*3 + idpd2,2) + PDrot_dA((idpd1 - 1)*3 + idpd2,3))*PProt(idpp,1);
              block(5*(idpd2+2) + idpd1,4 + idpp) += -intn[0]*PDrot((idpd1 - 1)*3 + idpd2,1)*PProt_dA(idpp,1) - intn[1]*(PDrot((idpd1 - 1)*3 + idpd2,2) + PDrot((idpd1 - 1)*3 + idpd2,3))*PProt_dA(idpp,1);
              block(5*(idpd2+2) + idpd1,4 + idpp) -= intn[2]*PDrot_dA((idpd1 - 1)*3 + idpd2,1)*(PProt(idpp,2) + PProt(idpp,3));
              block(5*(idpd2+2) + idpd1,4 + idpp) -= intn[2]*PDrot((idpd1 - 1)*3 + idpd2,1)*(PProt_dA(idpp,2) + PProt_dA(idpp,3));
              block(5*(idpd2+2) + idpd1,4 + idpp) -= intn[3]*(PDrot_dA((idpd1 - 1)*3 + idpd2,2) + PDrot_dA((idpd1 - 1)*3 + idpd2,3))*(PProt(idpp,2) + PProt(idpp,3));
              block(5*(idpd2+2) + idpd1,4 + idpp) -= intn[3]*(PDrot((idpd1 - 1)*3 + idpd2,2) + PDrot((idpd1 - 1)*3 + idpd2,3))*(PProt_dA(idpp,2) + PProt_dA(idpp,3));
              block(5*(idpd2+2) + idpd1,4 + idpp) -= intn[4]*(PDrot_dA((idpd1 - 1)*3 + idpd2,4)*PProt(idpp,4) + PDrot_dA((idpd1 - 1)*3 + idpd2,5)*PProt(idpp,5));
              block(5*(idpd2+2) + idpd1,4 + idpp) -= intn[4]*(PDrot((idpd1 - 1)*3 + idpd2,4)*PProt_dA(idpp,4) + PDrot((idpd1 - 1)*3 + idpd2,5)*PProt_dA(idpp,5));
              block(5*(idpd2+2) + idpd1,4 + idpp) -= intn[5]*(PDrot_dA((idpd1 - 1)*3 + idpd2,8)*PProt(idpp,4) + PDrot_dA((idpd1 - 1)*3 + idpd2,12)*PProt(idpp,5));
              block(5*(idpd2+2) + idpd1,4 + idpp) -= intn[5]*(PDrot((idpd1 - 1)*3 + idpd2,8)*PProt_dA(idpp,4) + PDrot((idpd1 - 1)*3 + idpd2,12)*PProt_dA(idpp,5));
              block(5*(idpd2+2) + idpd1,4 + idpp) -= intn[6]*((PDrot_dA((idpd1 - 1)*3 + idpd2,10) + PDrot_dA((idpd1 - 1)*3 + idpd2,15))*PProt(idpp,4) + (PDrot_dA((idpd1 - 1)*3 + idpd2,11) - PDrot_dA((idpd1 - 1)*3 + idpd2,14))*PProt(idpp,5));
              block(5*(idpd2+2) + idpd1,4 + idpp) -= intn[6]*((PDrot((idpd1 - 1)*3 + idpd2,10) + PDrot((idpd1 - 1)*3 + idpd2,15))*PProt_dA(idpp,4) + (PDrot((idpd1 - 1)*3 + idpd2,11) - PDrot((idpd1 - 1)*3 + idpd2,14))*PProt_dA(idpp,5));
            }
          }
        }
        //(dd|pp) integrals
        Dd[0] = Dvalue(atmC,6);
        Dd[1] = Dvalue(atmD,3);
        intn[0] = d_eri2Center(2,0,2,0,1,0,1,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(2,1,2,1,1,0,1,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[2] = d_eri2Center(2,2,2,2,1,0,1,0,RCD,Dd,atmC,atmD,type);    //(dl,dl|sg,sg)=(-dl,-dl|sg,sg)
        intn[3] = d_eri2Center(2,0,2,0,1,1,1,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[4] = d_eri2Center(2,1,2,1,1,1,1,1,RCD,Dd,atmC,atmD,type);    //(pi,pi|pi,pi)=(-pi,-pi|-pi,-pi)
        intn[5] = d_eri2Center(2,2,2,2,1,1,1,1,RCD,Dd,atmC,atmD,type);    //(dl,dl|pi,pi)=(-dl,-dl|pi,pi)=(dl,dl|-pi,-pi)=(-dl,-dl|-pi,-pi)
        intn[6] = d_eri2Center(2,0,2,1,1,0,1,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[7] = d_eri2Center(2,-1,2,-1,1,1,1,1,RCD,Dd,atmC,atmD,type);  //(-pi,-pi|pi,pi)=(pi,pi|-pi,-pi)
        intn[8] = d_eri2Center(2,1,2,-1,1,1,1,-1,RCD,Dd,atmC,atmD,type);  //(pi,-pi|pi,-pi)
        intn[9] = d_eri2Center(2,0,2,2,1,1,1,1,RCD,Dd,atmC,atmD,type);    //(sg,dl|pi,pi)=(sg,-dl|pi,-pi)=-(sg,dl|-pi,-pi)
        intn[10] = d_eri2Center(2,1,2,2,1,0,1,1,RCD,Dd,atmC,atmD,type);   //(pi,dl|sg,pi)=(-pi,-dl|sg,pi)=-(-pi,dl|sg,-pi)=(pi,-dl|sg,-pi)
        //(dd|pp) rotate
        for (size_t iddd = 1; iddd < 16; ++iddd) {
          for (size_t idpp = 1; idpp < 7; ++idpp) {
            block(30 + iddd,4 + idpp) = intn[0]*DDrot_dA(iddd,1)*PProt(idpp,1) + intn[1]*(DDrot_dA(iddd,2) + DDrot_dA(iddd,3))*PProt(idpp,1) + intn[2]*(DDrot_dA(iddd,4) + DDrot_dA(iddd,5))*PProt(idpp,1);
            block(30 + iddd,4 + idpp) += intn[0]*DDrot(iddd,1)*PProt_dA(idpp,1) + intn[1]*(DDrot(iddd,2) + DDrot(iddd,3))*PProt_dA(idpp,1) + intn[2]*(DDrot(iddd,4) + DDrot(iddd,5))*PProt_dA(idpp,1);
            block(30 + iddd,4 + idpp) += intn[3]*DDrot_dA(iddd,1)*(PProt(idpp,2) + PProt(idpp,3)) + intn[4]*(DDrot_dA(iddd,2)*PProt(idpp,2) + DDrot_dA(iddd,3)*PProt(idpp,3));
            block(30 + iddd,4 + idpp) += intn[3]*DDrot(iddd,1)*(PProt_dA(idpp,2) + PProt_dA(idpp,3)) + intn[4]*(DDrot(iddd,2)*PProt_dA(idpp,2) + DDrot(iddd,3)*PProt_dA(idpp,3));
            block(30 + iddd,4 + idpp) += intn[5]*(DDrot_dA(iddd,4) + DDrot_dA(iddd,5))*(PProt(idpp,2) + PProt(idpp,3)) + intn[6]*(DDrot_dA(iddd,6)*PProt(idpp,4) + DDrot_dA(iddd,7)*PProt(idpp,5));
            block(30 + iddd,4 + idpp) += intn[5]*(DDrot(iddd,4) + DDrot(iddd,5))*(PProt_dA(idpp,2) + PProt_dA(idpp,3)) + intn[6]*(DDrot(iddd,6)*PProt_dA(idpp,4) + DDrot(iddd,7)*PProt_dA(idpp,5));
            block(30 + iddd,4 + idpp) += intn[7]*(DDrot_dA(iddd,2)*PProt(idpp,3) + DDrot_dA(iddd,3)*PProt(idpp,2)) + intn[8]*DDrot_dA(iddd,10)*PProt(idpp,6);
            block(30 + iddd,4 + idpp) += intn[7]*(DDrot(iddd,2)*PProt_dA(idpp,3) + DDrot(iddd,3)*PProt_dA(idpp,2)) + intn[8]*DDrot(iddd,10)*PProt_dA(idpp,6);
            block(30 + iddd,4 + idpp) += intn[9]*(DDrot_dA(iddd,8)*(PProt(idpp,2) - PProt(idpp,3)) + DDrot_dA(iddd,9)*PProt(idpp,6));
            block(30 + iddd,4 + idpp) += intn[9]*(DDrot(iddd,8)*(PProt_dA(idpp,2) - PProt_dA(idpp,3)) + DDrot(iddd,9)*PProt_dA(idpp,6));
            block(30 + iddd,4 + idpp) += intn[10]*((DDrot_dA(iddd,11) + DDrot_dA(iddd,14))*PProt(idpp,4) + (DDrot_dA(iddd,12) - DDrot_dA(iddd,13))*PProt(idpp,5));
            block(30 + iddd,4 + idpp) += intn[10]*((DDrot(iddd,11) + DDrot(iddd,14))*PProt_dA(idpp,4) + (DDrot(iddd,12) - DDrot(iddd,13))*PProt_dA(idpp,5));
          }
        }
      }
      if (ncols > 10) {
        //(sp|sd) integrals
        Dd[0] = Dvalue(atmC,1);
        Dd[1] = Dvalue(atmD,4);
        intn[0] = d_eri2Center(0,0,1,0,0,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(0,0,1,1,0,0,2,1,RCD,Dd,atmC,atmD,type);   //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        //(sp|sd) rotate
        for (size_t idsp = 1; idsp < 4; ++idsp) {
          for (size_t idsd = 1; idsd < 6; ++idsd) {
            block(1 + idsp,10 + idsd) = -intn[0]*SDrot_dA(idsd,1)*SProt(idsp,1) - intn[1]*(SDrot_dA(idsd,2)*SProt(idsp,2) + SDrot_dA(idsd,3)*SProt(idsp,3));
            block(1 + idsp,10 + idsd) -= intn[0]*SDrot(idsd,1)*SProt_dA(idsp,1) + intn[1]*(SDrot(idsd,2)*SProt_dA(idsp,2) + SDrot(idsd,3)*SProt_dA(idsp,3));
          }
        }
        //(pp|sd) integrals
        Dd[0] = Dvalue(atmC,3);
        Dd[1] = Dvalue(atmD,4);
        intn[0] = d_eri2Center(1,0,1,0,0,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(1,1,1,1,0,0,2,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[2] = d_eri2Center(1,0,1,1,0,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[3] = d_eri2Center(1,1,1,1,0,0,2,2,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,dl)=-(-pi,-pi|sg,dl)=(pi,-pi|sg,-dl)
        //(pp|sd) rotate
        for (size_t idsd = 1; idsd < 6; ++idsd) {
          for (size_t idpp = 1; idpp < 7; ++idpp) {
            block(4 + idpp,10 + idsd) = intn[0]*SDrot_dA(idsd,1)*PProt(idpp,1) + intn[1]*SDrot_dA(idsd,1)*(PProt(idpp,2) + PProt(idpp,3));
            block(4 + idpp,10 + idsd) += intn[0]*SDrot(idsd,1)*PProt_dA(idpp,1) + intn[1]*SDrot(idsd,1)*(PProt_dA(idpp,2) + PProt_dA(idpp,3));
            block(4 + idpp,10 + idsd) += intn[2]*(SDrot_dA(idsd,2)*PProt(idpp,4) + SDrot_dA(idsd,3)*PProt(idpp,5));
            block(4 + idpp,10 + idsd) += intn[2]*(SDrot(idsd,2)*PProt_dA(idpp,4) + SDrot(idsd,3)*PProt_dA(idpp,5));
            block(4 + idpp,10 + idsd) += intn[3]*(SDrot_dA(idsd,4)*PProt(idpp,2) - SDrot_dA(idsd,4)*PProt(idpp,3) + SDrot_dA(idsd,5)*PProt(idpp,6));
            block(4 + idpp,10 + idsd) += intn[3]*(SDrot(idsd,4)*PProt_dA(idpp,2) - SDrot(idsd,4)*PProt_dA(idpp,3) + SDrot(idsd,5)*PProt_dA(idpp,6));
          }
        }
        //(sp|pd) integrals
        Dd[0] = Dvalue(atmC,1);
        Dd[1] = Dvalue(atmD,5);
        intn[0] = d_eri2Center(0,0,1,0,1,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(0,0,1,0,1,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center(0,0,1,1,1,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[3] = d_eri2Center(0,0,1,1,1,1,2,0,RCD,Dd,atmC,atmD,type);    //(sg,pi|pi,sg)=(sg,-pi|-pi,sg)
        intn[4] = d_eri2Center(0,0,1,1,1,1,2,2,RCD,Dd,atmC,atmD,type);    //(sg,pi|pi,dl)=(sg,pi|-pi,-dl)=-(sg,-pi|-pi,dl)=(sg,-pi|pi,-dl)
        //(sp|pd) rotate
        for (size_t idsp = 1; idsp < 4; ++idsp) {
          for (size_t idpd1 = 1; idpd1 < 6; ++idpd1) {
            for (size_t idpd2 = 1; idpd2 < 4; ++idpd2) {
              block(1 + idsp,5*(idpd2+2) + idpd1) = intn[0]*PDrot_dA((idpd1 - 1)*3 + idpd2,1)*SProt(idsp,1) + intn[1]*(PDrot_dA((idpd1 - 1)*3 + idpd2,2) + PDrot_dA((idpd1 - 1)*3 + idpd2,3))*SProt(idsp,1);
              block(1 + idsp,5*(idpd2+2) + idpd1) += intn[0]*PDrot((idpd1 - 1)*3 + idpd2,1)*SProt_dA(idsp,1) + intn[1]*(PDrot((idpd1 - 1)*3 + idpd2,2) + PDrot((idpd1 - 1)*3 + idpd2,3))*SProt_dA(idsp,1);
              block(1 + idsp,5*(idpd2+2) + idpd1) += intn[2]*(PDrot_dA((idpd1 - 1)*3 + idpd2,4)*SProt(idsp,2) + PDrot_dA((idpd1 - 1)*3 + idpd2,5)*SProt(idsp,3));
              block(1 + idsp,5*(idpd2+2) + idpd1) += intn[2]*(PDrot((idpd1 - 1)*3 + idpd2,4)*SProt_dA(idsp,2) + PDrot((idpd1 - 1)*3 + idpd2,5)*SProt_dA(idsp,3));
              block(1 + idsp,5*(idpd2+2) + idpd1) += intn[3]*(PDrot_dA((idpd1 - 1)*3 + idpd2,8)*SProt(idsp,2) + PDrot_dA((idpd1 - 1)*3 + idpd2,12)*SProt(idsp,3));
              block(1 + idsp,5*(idpd2+2) + idpd1) += intn[3]*(PDrot((idpd1 - 1)*3 + idpd2,8)*SProt_dA(idsp,2) + PDrot((idpd1 - 1)*3 + idpd2,12)*SProt_dA(idsp,3));
              block(1 + idsp,5*(idpd2+2) + idpd1) += intn[4]*((PDrot_dA((idpd1 - 1)*3 + idpd2,10) + PDrot_dA((idpd1 - 1)*3 + idpd2,15))*SProt(idsp,2) + (PDrot_dA((idpd1 - 1)*3 + idpd2,11) - PDrot_dA((idpd1 - 1)*3 + idpd2,14))*SProt(idsp,3));
              block(1 + idsp,5*(idpd2+2) + idpd1) += intn[4]*((PDrot((idpd1 - 1)*3 + idpd2,10) + PDrot((idpd1 - 1)*3 + idpd2,15))*SProt_dA(idsp,2) + (PDrot((idpd1 - 1)*3 + idpd2,11) - PDrot((idpd1 - 1)*3 + idpd2,14))*SProt_dA(idsp,3));
            }
          }
        }
        //(sp|dd) integrals
        Dd[0] = Dvalue(atmC,1);
        Dd[1] = Dvalue(atmD,6);
        intn[0] = d_eri2Center(0,0,1,0,2,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(0,0,1,0,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center(0,0,1,0,2,2,2,2,RCD,Dd,atmC,atmD,type);    //(sg,sg|dl,dl)=(sg,sg|-dl,-dl)
        intn[3] = d_eri2Center(0,0,1,1,2,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[4] = d_eri2Center(0,0,1,1,2,1,2,2,RCD,Dd,atmC,atmD,type);    //(sg,pi|pi,dl)=(sg,pi|-pi,-dl)=-(sg,-pi|-pi,dl)=(sg,-pi|pi,-dl)
        //(sp|dd) rotate
        for (size_t idsp = 1; idsp < 4; ++idsp) {
          for (size_t iddd = 1; iddd < 16; ++iddd) {
            block(1 + idsp,30 + iddd) = -intn[0]*DDrot_dA(iddd,1)*SProt(idsp,1) - intn[1]*(DDrot_dA(iddd,2) + DDrot_dA(iddd,3))*SProt(idsp,1) - intn[2]*(DDrot_dA(iddd,4) + DDrot_dA(iddd,5))*SProt(idsp,1);
            block(1 + idsp,30 + iddd) -= intn[0]*DDrot(iddd,1)*SProt_dA(idsp,1) + intn[1]*(DDrot(iddd,2) + DDrot(iddd,3))*SProt_dA(idsp,1) + intn[2]*(DDrot(iddd,4) + DDrot(iddd,5))*SProt_dA(idsp,1);
            block(1 + idsp,30 + iddd) -= intn[3]*(DDrot_dA(iddd,6)*SProt(idsp,2) + DDrot_dA(iddd,7)*SProt(idsp,3));
            block(1 + idsp,30 + iddd) -= intn[3]*(DDrot(iddd,6)*SProt_dA(idsp,2) + DDrot(iddd,7)*SProt_dA(idsp,3));
            block(1 + idsp,30 + iddd) -= intn[4]*((DDrot_dA(iddd,11) + DDrot_dA(iddd,14))*SProt(idsp,2) + (DDrot_dA(iddd,12) - DDrot_dA(iddd,13))*SProt(idsp,3));
            block(1 + idsp,30 + iddd) -= intn[4]*((DDrot(iddd,11) + DDrot(iddd,14))*SProt_dA(idsp,2) + (DDrot(iddd,12) - DDrot(iddd,13))*SProt_dA(idsp,3));
          }
        }
        //(pp|pd) integrals
        Dd[0] = Dvalue(atmC,3);
        Dd[1] = Dvalue(atmD,5);
        intn[0] = d_eri2Center(1,0,1,0,1,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(1,0,1,0,1,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center(1,1,1,1,1,0,2,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[3] = d_eri2Center(1,1,1,1,1,1,2,1,RCD,Dd,atmC,atmD,type);    //(pi,pi|pi,pi)=(pi,pi|-pi,-pi)=(-pi,-pi|pi,pi)=(-pi,-pi|-pi,-pi)
        intn[4] = d_eri2Center(1,0,1,1,1,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[5] = d_eri2Center(1,1,1,0,1,1,2,0,RCD,Dd,atmC,atmD,type);    //(pi,sg|pi,sg)=(-pi,sg|-pi,sg)
        intn[6] = d_eri2Center(1,0,1,1,1,1,2,2,RCD,Dd,atmC,atmD,type);    //(sg,pi|pi,dl)=(sg,pi|-pi,-dl)=-(sg,-pi|-pi,dl)=(sg,-pi|pi,-dl)
        //(pp|pd) rotate
        for (size_t idpp = 1; idpp < 7; ++idpp) {
          for (size_t idpd1 = 1; idpd1 < 6; ++idpd1) {
            for (size_t idpd2 = 1; idpd2 < 4; ++idpd2) {
              block(4 + idpp,5*(idpd2+2) + idpd1) = -intn[0]*PDrot_dA((idpd1 - 1)*3 + idpd2,1)*PProt(idpp,1) - intn[1]*(PDrot_dA((idpd1 - 1)*3 + idpd2,2) + PDrot_dA((idpd1 - 1)*3 + idpd2,3))*PProt(idpp,1);
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[0]*PDrot((idpd1 - 1)*3 + idpd2,1)*PProt_dA(idpp,1) + intn[1]*(PDrot((idpd1 - 1)*3 + idpd2,2) + PDrot((idpd1 - 1)*3 + idpd2,3))*PProt_dA(idpp,1);
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[2]*PDrot_dA((idpd1 - 1)*3 + idpd2,1)*(PProt(idpp,2) + PProt(idpp,3));
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[2]*PDrot((idpd1 - 1)*3 + idpd2,1)*(PProt_dA(idpp,2) + PProt_dA(idpp,3));
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[3]*(PDrot_dA((idpd1 - 1)*3 + idpd2,2) + PDrot_dA((idpd1 - 1)*3 + idpd2,3))*(PProt(idpp,2) + PProt(idpp,3));
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[3]*(PDrot((idpd1 - 1)*3 + idpd2,2) + PDrot((idpd1 - 1)*3 + idpd2,3))*(PProt_dA(idpp,2) + PProt_dA(idpp,3));
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[4]*(PDrot_dA((idpd1 - 1)*3 + idpd2,4)*PProt(idpp,4) + PDrot_dA((idpd1 - 1)*3 + idpd2,5)*PProt(idpp,5));
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[4]*(PDrot((idpd1 - 1)*3 + idpd2,4)*PProt_dA(idpp,4) + PDrot((idpd1 - 1)*3 + idpd2,5)*PProt_dA(idpp,5));
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[5]*(PDrot_dA((idpd1 - 1)*3 + idpd2,8)*PProt(idpp,4) + PDrot_dA((idpd1 - 1)*3 + idpd2,12)*PProt(idpp,5));
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[5]*(PDrot((idpd1 - 1)*3 + idpd2,8)*PProt_dA(idpp,4) + PDrot((idpd1 - 1)*3 + idpd2,12)*PProt_dA(idpp,5));
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[6]*((PDrot_dA((idpd1 - 1)*3 + idpd2,10) + PDrot_dA((idpd1 - 1)*3 + idpd2,15))*PProt(idpp,4) + (PDrot_dA((idpd1 - 1)*3 + idpd2,11) - PDrot_dA((idpd1 - 1)*3 + idpd2,14))*PProt(idpp,5));
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[6]*((PDrot((idpd1 - 1)*3 + idpd2,10) + PDrot((idpd1 - 1)*3 + idpd2,15))*PProt_dA(idpp,4) + (PDrot((idpd1 - 1)*3 + idpd2,11) - PDrot((idpd1 - 1)*3 + idpd2,14))*PProt_dA(idpp,5));
            }
          }
        }
        //(pp|dd) integrals
        Dd[0] = Dvalue(atmC,3);
        Dd[1] = Dvalue(atmD,6);
        intn[0] = d_eri2Center(1,0,1,0,2,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(1,0,1,0,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center(1,0,1,0,2,2,2,2,RCD,Dd,atmC,atmD,type);    //(sg,sg|dl,dl)=(sg,sg|-dl,-dl)
        intn[3] = d_eri2Center(1,1,1,1,2,0,2,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[4] = d_eri2Center(1,1,1,1,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(pi,pi|pi,pi)=(-pi,-pi|-pi,-pi)
        intn[5] = d_eri2Center(1,1,1,1,2,2,2,2,RCD,Dd,atmC,atmD,type);    //(pi,pi|dl,dl)=(pi,pi|-dl,-dl)=(-pi,-pi|dl,dl)=(-pi,-pi|-dl,-dl)
        intn[6] = d_eri2Center(1,0,1,1,2,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[7] = d_eri2Center(1,-1,1,-1,2,1,2,1,RCD,Dd,atmC,atmD,type);  //(-pi,-pi|pi,pi)=(pi,pi|-pi,-pi)
        intn[8] = d_eri2Center(1,1,1,-1,2,1,2,-1,RCD,Dd,atmC,atmD,type);  //(pi,-pi|pi,-pi)
        intn[9] = d_eri2Center(1,1,1,1,2,0,2,2,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,dl)=(pi,-pi|sg,-dl)=-(-pi,-pi|sg,dl)
        intn[10] = d_eri2Center(1,0,1,1,2,1,2,2,RCD,Dd,atmC,atmD,type);   //(sg,pi|pi,dl)=(sg,pi|-pi,-dl)=-(sg,-pi|-pi,dl)=(sg,-pi|pi,-dl)
        //(pp|dd) rotate
        for (size_t iddd = 1; iddd < 16; ++iddd) {
          for (size_t idpp = 1; idpp < 7; ++idpp) {
            block(4 + idpp,30 + iddd) = intn[0]*DDrot_dA(iddd,1)*PProt(idpp,1) + intn[1]*(DDrot_dA(iddd,2) + DDrot_dA(iddd,3))*PProt(idpp,1) + intn[2]*(DDrot_dA(iddd,4) + DDrot_dA(iddd,5))*PProt(idpp,1);
            block(4 + idpp,30 + iddd) += intn[0]*DDrot(iddd,1)*PProt_dA(idpp,1) + intn[1]*(DDrot(iddd,2) + DDrot(iddd,3))*PProt_dA(idpp,1) + intn[2]*(DDrot(iddd,4) + DDrot(iddd,5))*PProt_dA(idpp,1);
            block(4 + idpp,30 + iddd) += intn[3]*DDrot_dA(iddd,1)*(PProt(idpp,2) + PProt(idpp,3)) + intn[4]*(DDrot_dA(iddd,2)*PProt(idpp,2) + DDrot_dA(iddd,3)*PProt(idpp,3));
            block(4 + idpp,30 + iddd) += intn[3]*DDrot(iddd,1)*(PProt_dA(idpp,2) + PProt_dA(idpp,3)) + intn[4]*(DDrot(iddd,2)*PProt_dA(idpp,2) + DDrot(iddd,3)*PProt_dA(idpp,3));
            block(4 + idpp,30 + iddd) += intn[5]*(DDrot_dA(iddd,4) + DDrot_dA(iddd,5))*(PProt(idpp,2) + PProt(idpp,3)) + intn[6]*(DDrot_dA(iddd,6)*PProt(idpp,4) + DDrot_dA(iddd,7)*PProt(idpp,5));
            block(4 + idpp,30 + iddd) += intn[5]*(DDrot(iddd,4) + DDrot(iddd,5))*(PProt_dA(idpp,2) + PProt_dA(idpp,3)) + intn[6]*(DDrot(iddd,6)*PProt_dA(idpp,4) + DDrot(iddd,7)*PProt_dA(idpp,5));
            block(4 + idpp,30 + iddd) += intn[7]*(DDrot_dA(iddd,2)*PProt(idpp,3) + DDrot_dA(iddd,3)*PProt(idpp,2)) + intn[8]*DDrot_dA(iddd,10)*PProt(idpp,6);
            block(4 + idpp,30 + iddd) += intn[7]*(DDrot(iddd,2)*PProt_dA(idpp,3) + DDrot(iddd,3)*PProt_dA(idpp,2)) + intn[8]*DDrot(iddd,10)*PProt_dA(idpp,6);
            block(4 + idpp,30 + iddd) += intn[9]*(DDrot_dA(iddd,8)*(PProt(idpp,2) - PProt(idpp,3)) + DDrot_dA(iddd,9)*PProt(idpp,6));
            block(4 + idpp,30 + iddd) += intn[9]*(DDrot(iddd,8)*(PProt_dA(idpp,2) - PProt_dA(idpp,3)) + DDrot(iddd,9)*PProt_dA(idpp,6));
            block(4 + idpp,30 + iddd) += intn[10]*((DDrot_dA(iddd,11) + DDrot_dA(iddd,14))*PProt(idpp,4) + (DDrot_dA(iddd,12) - DDrot_dA(iddd,13))*PProt(idpp,5));
            block(4 + idpp,30 + iddd) += intn[10]*((DDrot(iddd,11) + DDrot(iddd,14))*PProt_dA(idpp,4) + (DDrot(iddd,12) - DDrot(iddd,13))*PProt_dA(idpp,5));
          }
        }
      }
      if ((nrows > 10)&&(ncols > 10)) {
        //(sd|sd) integrals
        Dd[0] = Dvalue(atmC,4);
        Dd[1] = Dvalue(atmD,4);
        intn[0] = d_eri2Center(0,0,2,0,0,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(0,0,2,1,0,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[2] = d_eri2Center(0,0,2,2,0,0,2,2,RCD,Dd,atmC,atmD,type);    //(sg,dl|sg,dl)=(sg,-dl|sg,-dl)
        //(sd|sd) rotate
        for (size_t idsdb = 1; idsdb < 6; ++idsdb) {
          for (size_t idsdk = 1; idsdk < 6; ++idsdk) {
            block(10 + idsdb,10 + idsdk) = intn[0]*SDrot_dA(idsdb,1)*SDrot(idsdk,1) + intn[1]*(SDrot_dA(idsdb,2)*SDrot(idsdk,2) + SDrot_dA(idsdb,3)*SDrot(idsdk,3)) + intn[2]*(SDrot_dA(idsdb,4)*SDrot(idsdk,4) + SDrot_dA(idsdb,5)*SDrot(idsdk,5));
            block(10 + idsdb,10 + idsdk) += intn[0]*SDrot(idsdb,1)*SDrot_dA(idsdk,1) + intn[1]*(SDrot(idsdb,2)*SDrot_dA(idsdk,2) + SDrot(idsdb,3)*SDrot_dA(idsdk,3)) + intn[2]*(SDrot(idsdb,4)*SDrot_dA(idsdk,4) + SDrot(idsdb,5)*SDrot_dA(idsdk,5));
          }
        }
        //(sd|pd) integrals
        Dd[1] = Dvalue(atmD,5);
        intn[0] = d_eri2Center(0,0,2,0,1,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(0,0,2,0,1,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center(0,0,2,1,1,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[3] = d_eri2Center(0,0,2,1,1,1,2,2,RCD,Dd,atmC,atmD,type);    //(sg,pi|pi,dl)=-(sg,-pi|-pi,dl)=(sg,-pi|pi,-dl)=(sg,pi|-pi,-dl)
        intn[4] = d_eri2Center(0,0,2,1,1,1,2,0,RCD,Dd,atmC,atmD,type);    //(sg,pi|pi,sg)=(sg,-pi|-pi,sg)
        //(sd|pd) rotate
        for (size_t idsd = 1; idsd < 6; ++idsd) {
          for (size_t idpd1 = 1; idpd1 < 6; ++idpd1) {
            for (size_t idpd2 = 1; idpd2 < 4; ++idpd2) {
              block(10 + idsd,5*(idpd2+2) + idpd1) = -intn[0]*PDrot_dA((idpd1 - 1)*3 + idpd2,1)*SDrot(idsd,1) - intn[1]*(PDrot_dA((idpd1 - 1)*3 + idpd2,2) + PDrot_dA((idpd1 - 1)*3 + idpd2,3))*SDrot(idsd,1);
              block(10 + idsd,5*(idpd2+2) + idpd1) -= intn[0]*PDrot((idpd1 - 1)*3 + idpd2,1)*SDrot_dA(idsd,1) + intn[1]*(PDrot((idpd1 - 1)*3 + idpd2,2) + PDrot((idpd1 - 1)*3 + idpd2,3))*SDrot_dA(idsd,1);
              block(10 + idsd,5*(idpd2+2) + idpd1) -= intn[2]*(PDrot_dA((idpd1 - 1)*3 + idpd2,4)*SDrot(idsd,2) + PDrot_dA((idpd1 - 1)*3 + idpd2,5)*SDrot(idsd,3));
              block(10 + idsd,5*(idpd2+2) + idpd1) -= intn[2]*(PDrot((idpd1 - 1)*3 + idpd2,4)*SDrot_dA(idsd,2) + PDrot((idpd1 - 1)*3 + idpd2,5)*SDrot_dA(idsd,3));
              block(10 + idsd,5*(idpd2+2) + idpd1) -= intn[3]*((PDrot_dA((idpd1 - 1)*3 + idpd2,10) + PDrot_dA((idpd1 - 1)*3 + idpd2,15))*SDrot(idsd,2) + (PDrot_dA((idpd1 - 1)*3 + idpd2,11) - PDrot_dA((idpd1 - 1)*3 + idpd2,14))*SDrot(idsd,3));
              block(10 + idsd,5*(idpd2+2) + idpd1) -= intn[3]*((PDrot((idpd1 - 1)*3 + idpd2,10) + PDrot((idpd1 - 1)*3 + idpd2,15))*SDrot_dA(idsd,2) + (PDrot((idpd1 - 1)*3 + idpd2,11) - PDrot((idpd1 - 1)*3 + idpd2,14))*SDrot_dA(idsd,3));
              block(10 + idsd,5*(idpd2+2) + idpd1) -= intn[4]*(PDrot_dA((idpd1 - 1)*3 + idpd2,8)*SDrot(idsd,2) + PDrot_dA((idpd1 - 1)*3 + idpd2,12)*SDrot(idsd,3));
              block(10 + idsd,5*(idpd2+2) + idpd1) -= intn[4]*(PDrot((idpd1 - 1)*3 + idpd2,8)*SDrot_dA(idsd,2) + PDrot((idpd1 - 1)*3 + idpd2,12)*SDrot_dA(idsd,3));
            }
          }
        }
        //(pd|sd) integrals
        Dd[0] = Dvalue(atmC,5);
        Dd[1] = Dvalue(atmD,4);
        intn[0] = d_eri2Center(1,0,2,0,0,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(1,1,2,1,0,0,2,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[2] = d_eri2Center(1,0,2,1,0,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[3] = d_eri2Center(1,1,2,2,0,0,2,1,RCD,Dd,atmC,atmD,type);    //(pi,dl|sg,pi)=-(-pi,dl|sg,-pi)=(pi,-dl|sg,-pi)=(-pi,-dl|sg,pi)
        intn[4] = d_eri2Center(1,1,2,0,0,0,2,1,RCD,Dd,atmC,atmD,type);    //(pi,sg|sg,pi)=(-pi,sg|sg,-pi)
        //(pd|sd) rotate
        for (size_t idpd1 = 1; idpd1 < 6; ++idpd1) {
          for (size_t idpd2 = 1; idpd2 < 4; ++idpd2) {
            for (size_t idsd = 1; idsd < 6; ++idsd) {
              block(5*(idpd2+2) + idpd1,10 + idsd) = -intn[0]*PDrot_dA((idpd1 - 1)*3 + idpd2,1)*SDrot(idsd,1) - intn[1]*(PDrot_dA((idpd1 - 1)*3 + idpd2,2) + PDrot_dA((idpd1 - 1)*3 + idpd2,3))*SDrot(idsd,1);
              block(5*(idpd2+2) + idpd1,10 + idsd) -= intn[0]*PDrot((idpd1 - 1)*3 + idpd2,1)*SDrot_dA(idsd,1) + intn[1]*(PDrot((idpd1 - 1)*3 + idpd2,2) + PDrot((idpd1 - 1)*3 + idpd2,3))*SDrot_dA(idsd,1);
              block(5*(idpd2+2) + idpd1,10 + idsd) -= intn[2]*(PDrot_dA((idpd1 - 1)*3 + idpd2,4)*SDrot(idsd,2) + PDrot_dA((idpd1 - 1)*3 + idpd2,5)*SDrot(idsd,3));
              block(5*(idpd2+2) + idpd1,10 + idsd) -= intn[2]*(PDrot((idpd1 - 1)*3 + idpd2,4)*SDrot_dA(idsd,2) + PDrot((idpd1 - 1)*3 + idpd2,5)*SDrot_dA(idsd,3));
              block(5*(idpd2+2) + idpd1,10 + idsd) -= intn[3]*((PDrot_dA((idpd1 - 1)*3 + idpd2,10) + PDrot_dA((idpd1 - 1)*3 + idpd2,15))*SDrot(idsd,2) + (PDrot_dA((idpd1 - 1)*3 + idpd2,11) - PDrot_dA((idpd1 - 1)*3 + idpd2,14))*SDrot(idsd,3));
              block(5*(idpd2+2) + idpd1,10 + idsd) -= intn[3]*((PDrot((idpd1 - 1)*3 + idpd2,10) + PDrot((idpd1 - 1)*3 + idpd2,15))*SDrot_dA(idsd,2) + (PDrot((idpd1 - 1)*3 + idpd2,11) - PDrot((idpd1 - 1)*3 + idpd2,14))*SDrot_dA(idsd,3));
              block(5*(idpd2+2) + idpd1,10 + idsd) -= intn[4]*(PDrot_dA((idpd1 - 1)*3 + idpd2,8)*SDrot(idsd,2) + PDrot_dA((idpd1 - 1)*3 + idpd2,12)*SDrot(idsd,3));
              block(5*(idpd2+2) + idpd1,10 + idsd) -= intn[4]*(PDrot((idpd1 - 1)*3 + idpd2,8)*SDrot_dA(idsd,2) + PDrot((idpd1 - 1)*3 + idpd2,12)*SDrot_dA(idsd,3));
            }
          }
        }
        //(dd|sd) integrals
        Dd[0] = Dvalue(atmC,6);
        Dd[1] = Dvalue(atmD,4);
        intn[0] = d_eri2Center(2,0,2,0,0,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(2,1,2,1,0,0,2,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[2] = d_eri2Center(2,2,2,2,0,0,2,0,RCD,Dd,atmC,atmD,type);    //(dl,dl|sg,sg)=(-dl,-dl|sg,sg)
        intn[3] = d_eri2Center(2,0,2,1,0,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[4] = d_eri2Center(2,0,2,2,0,0,2,2,RCD,Dd,atmC,atmD,type);    //(sg,dl|sg,dl)=(sg,-dl|sg,-dl)
        intn[5] = d_eri2Center(2,1,2,1,0,0,2,2,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,dl)=-(-pi,-pi|sg,dl)=(pi,-pi|sg,-dl)
        intn[6] = d_eri2Center(2,1,2,2,0,0,2,1,RCD,Dd,atmC,atmD,type);    //(pi,dl|sg,pi)=(-pi,-dl|sg,pi)=(pi,-dl|sg,-pi)=-(-pi,dl|sg,-pi)
        //(dd|sd) rotate
        for (size_t iddd = 1; iddd < 16; ++iddd) {
          for (size_t idsd = 1; idsd < 6; ++idsd) {
            block(30 + iddd,10 + idsd) = intn[0]*DDrot_dA(iddd,1)*SDrot(idsd,1) + intn[1]*(DDrot_dA(iddd,2) + DDrot_dA(iddd,3))*SDrot(idsd,1) + intn[2]*(DDrot(iddd,4) + DDrot(iddd,5))*SDrot_dA(idsd,1);
            block(30 + iddd,10 + idsd) += intn[0]*DDrot(iddd,1)*SDrot_dA(idsd,1) + intn[1]*(DDrot(iddd,2) + DDrot(iddd,3))*SDrot_dA(idsd,1) + intn[2]*(DDrot_dA(iddd,4) + DDrot_dA(iddd,5))*SDrot(idsd,1);
            block(30 + iddd,10 + idsd) += intn[3]*(DDrot_dA(iddd,6)*SDrot(idsd,2) + DDrot_dA(iddd,7)*SDrot(idsd,3)) + intn[4]*(DDrot_dA(iddd,8)*SDrot(idsd,4) + DDrot_dA(iddd,9)*SDrot(idsd,5));
            block(30 + iddd,10 + idsd) += intn[3]*(DDrot(iddd,6)*SDrot_dA(idsd,2) + DDrot(iddd,7)*SDrot_dA(idsd,3)) + intn[4]*(DDrot(iddd,8)*SDrot_dA(idsd,4) + DDrot(iddd,9)*SDrot_dA(idsd,5));
            block(30 + iddd,10 + idsd) += intn[5]*((DDrot_dA(iddd,2) - DDrot_dA(iddd,3))*SDrot(idsd,4) + DDrot_dA(iddd,10)*SDrot(idsd,5));
            block(30 + iddd,10 + idsd) += intn[5]*((DDrot(iddd,2) - DDrot(iddd,3))*SDrot_dA(idsd,4) + DDrot(iddd,10)*SDrot_dA(idsd,5));
            block(30 + iddd,10 + idsd) += intn[6]*((DDrot_dA(iddd,11) + DDrot_dA(iddd,14))*SDrot(idsd,2) + (DDrot(iddd,12) - DDrot(iddd,13))*SDrot_dA(idsd,3));
            block(30 + iddd,10 + idsd) += intn[6]*((DDrot(iddd,11) + DDrot(iddd,14))*SDrot_dA(idsd,2) + (DDrot_dA(iddd,12) - DDrot_dA(iddd,13))*SDrot(idsd,3));
          }
        }
        //(sd|dd) integrals
        Dd[0] = Dvalue(atmC,4);
        Dd[1] = Dvalue(atmD,6);
        intn[0] = d_eri2Center(0,0,2,0,2,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(0,0,2,0,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center(0,0,2,0,2,2,2,2,RCD,Dd,atmC,atmD,type);    //(sg,sg|dl,dl)=(sg,sg|-dl,-dl)
        intn[3] = d_eri2Center(0,0,2,1,2,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[4] = d_eri2Center(0,0,2,2,2,0,2,2,RCD,Dd,atmC,atmD,type);    //(sg,dl|sg,dl)=(sg,-dl|sg,-dl)
        intn[5] = d_eri2Center(0,0,2,2,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,dl|pi,pi)=-(sg,dl|-pi,-pi)=(sg,-dl|pi,-pi)
        intn[6] = d_eri2Center(0,0,2,1,2,1,2,2,RCD,Dd,atmC,atmD,type);    //(sg,pi|pi,dl)=(sg,pi|-pi,-dl)=(sg,-pi|pi,-dl)=-(sg,-pi|-pi,dl)
        //(sd|dd) rotate
        for (size_t idsd = 1; idsd < 6; ++idsd) {
          for (size_t iddd = 1; iddd < 16; ++iddd) {
            block(10 + idsd,30 + iddd) = intn[0]*DDrot_dA(iddd,1)*SDrot(idsd,1) + intn[1]*(DDrot_dA(iddd,2) + DDrot_dA(iddd,3))*SDrot(idsd,1) + intn[2]*(DDrot(iddd,4) + DDrot(iddd,5))*SDrot_dA(idsd,1);
            block(10 + idsd,30 + iddd) += intn[0]*DDrot(iddd,1)*SDrot_dA(idsd,1) + intn[1]*(DDrot(iddd,2) + DDrot(iddd,3))*SDrot_dA(idsd,1) + intn[2]*(DDrot_dA(iddd,4) + DDrot_dA(iddd,5))*SDrot(idsd,1);
            block(10 + idsd,30 + iddd) += intn[3]*(DDrot_dA(iddd,6)*SDrot(idsd,2) + DDrot_dA(iddd,7)*SDrot(idsd,3)) + intn[4]*(DDrot_dA(iddd,8)*SDrot(idsd,4) + DDrot_dA(iddd,9)*SDrot(idsd,5));
            block(10 + idsd,30 + iddd) += intn[3]*(DDrot(iddd,6)*SDrot_dA(idsd,2) + DDrot(iddd,7)*SDrot_dA(idsd,3)) + intn[4]*(DDrot(iddd,8)*SDrot_dA(idsd,4) + DDrot(iddd,9)*SDrot_dA(idsd,5));
            block(10 + idsd,30 + iddd) += intn[5]*((DDrot_dA(iddd,2) - DDrot_dA(iddd,3))*SDrot(idsd,4) + DDrot_dA(iddd,10)*SDrot(idsd,5));
            block(10 + idsd,30 + iddd) += intn[5]*((DDrot(iddd,2) - DDrot(iddd,3))*SDrot_dA(idsd,4) + DDrot(iddd,10)*SDrot_dA(idsd,5));
            block(10 + idsd,30 + iddd) += intn[6]*((DDrot_dA(iddd,11) + DDrot_dA(iddd,14))*SDrot(idsd,2) + (DDrot(iddd,12) - DDrot(iddd,13))*SDrot_dA(idsd,3));
            block(10 + idsd,30 + iddd) += intn[6]*((DDrot(iddd,11) + DDrot(iddd,14))*SDrot_dA(idsd,2) + (DDrot_dA(iddd,12) - DDrot_dA(iddd,13))*SDrot(idsd,3));
          }
        }
        //(pd|pd) integrals
        Dd[0] = Dvalue(atmC,5);
        Dd[1] = Dvalue(atmD,5);
        intn[0] = d_eri2Center(1,0,2,0,1,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(1,1,2,1,1,0,2,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)=(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center(1,1,2,1,1,1,2,1,RCD,Dd,atmC,atmD,type);    //(pi,pi|pi,pi)=(-pi,-pi|pi,pi)=(pi,pi|-pi,-pi)=(-pi,-pi|-pi,-pi)
        intn[3] = d_eri2Center(1,0,2,1,1,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[4] = d_eri2Center(1,1,2,0,1,1,2,0,RCD,Dd,atmC,atmD,type);    //(pi,sg|pi,sg)=(-pi,sg|-pi,sg)
        intn[5] = d_eri2Center(1,1,2,2,1,1,2,2,RCD,Dd,atmC,atmD,type);    //(pi,dl|pi,dl)=(pi,dl|-pi,-dl)=(-pi,dl|-pi,dl)=-(-pi,dl|pi,-dl)=-(pi,-dl|-pi,dl)=(pi,-dl|pi,-dl)=(-pi,-dl|pi,dl)=(-pi,-dl|-pi,-dl)
        intn[6] = d_eri2Center(1,1,2,0,1,0,2,1,RCD,Dd,atmC,atmD,type);    //(pi,sg|sg,pi)=(-pi,sg|sg,-pi)=(sg,pi|pi,sg)=(sg,-pi|-pi,sg)
        intn[7] = d_eri2Center(1,1,2,2,1,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|pi,dl)=(sg,pi|-pi,-dl)=-(sg,-pi|-pi,dl)=(sg,-pi|pi,-dl)=(pi,dl|sg,pi)=-(-pi,dl|sg,-pi)=(pi,-dl|sg,-pi)=(-pi,-dl|sg,pi)
        intn[8] = d_eri2Center(1,1,2,2,1,1,2,0,RCD,Dd,atmC,atmD,type);    //(pi,sg|pi,dl)=(pi,sg|-pi,-dl)=-(-pi,sg|-pi,dl)=(-pi,sg|pi,-dl)=(pi,dl|pi,sg)=-(-pi,dl|-pi,sg)=(pi,-dl|-pi,sg)=(-pi,-dl|pi,sg)
        //(pd|pd) rotate
        for (size_t idpdb1 = 1; idpdb1 < 6; ++idpdb1) {            //d orbital on bra
          for (size_t idpdb2 = 1; idpdb2 < 4; ++idpdb2) {          //p orbital on bra
            for (size_t idpdk1 = 1; idpdk1 < 6; ++idpdk1) {        //d orbital on ket
              for (size_t idpdk2 = 1; idpdk2 < 4; ++idpdk2) {      //p orbital on ket
                cnt1 = (idpdb1 - 1)*3 + idpdb2;
                cnt2 = (idpdk1 - 1)*3 + idpdk2;
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) = intn[0]*PDrot_dA(cnt1,1)*PDrot(cnt2,1) + intn[1]*((PDrot_dA(cnt1,2) + PDrot_dA(cnt1,3))*PDrot(cnt2,1) + PDrot_dA(cnt1,1)*(PDrot(cnt2,2) + PDrot(cnt2,3)));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[0]*PDrot(cnt1,1)*PDrot_dA(cnt2,1) + intn[1]*((PDrot(cnt1,2) + PDrot(cnt1,3))*PDrot_dA(cnt2,1) + PDrot(cnt1,1)*(PDrot_dA(cnt2,2) + PDrot_dA(cnt2,3)));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[2]*(PDrot_dA(cnt1,2) + PDrot_dA(cnt1,3))*(PDrot(cnt2,2) + PDrot(cnt2,3)) + intn[3]*(PDrot_dA(cnt1,4)*PDrot(cnt2,4) + PDrot_dA(cnt1,5)*PDrot(cnt2,5));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[2]*(PDrot(cnt1,2) + PDrot(cnt1,3))*(PDrot_dA(cnt2,2) + PDrot_dA(cnt2,3)) + intn[3]*(PDrot(cnt1,4)*PDrot_dA(cnt2,4) + PDrot(cnt1,5)*PDrot_dA(cnt2,5));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[4]*(PDrot_dA(cnt1,8)*PDrot(cnt2,8) + PDrot_dA(cnt1,12)*PDrot(cnt2,12));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[4]*(PDrot(cnt1,8)*PDrot_dA(cnt2,8) + PDrot(cnt1,12)*PDrot_dA(cnt2,12));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[5]*((PDrot_dA(cnt1,10) + PDrot_dA(cnt1,15))*(PDrot(cnt2,10) + PDrot(cnt2,15)) + (PDrot_dA(cnt1,14) - PDrot_dA(cnt1,11))*(PDrot(cnt2,14) - PDrot(cnt2,11)));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[5]*((PDrot(cnt1,10) + PDrot(cnt1,15))*(PDrot_dA(cnt2,10) + PDrot_dA(cnt2,15)) + (PDrot(cnt1,14) - PDrot(cnt1,11))*(PDrot_dA(cnt2,14) - PDrot_dA(cnt2,11)));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[6]*(PDrot_dA(cnt1,4)*PDrot(cnt2,8) + PDrot_dA(cnt1,8)*PDrot(cnt2,4) + PDrot_dA(cnt1,5)*PDrot(cnt2,12) + PDrot_dA(cnt1,12)*PDrot(cnt2,5));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[6]*(PDrot(cnt1,4)*PDrot_dA(cnt2,8) + PDrot(cnt1,8)*PDrot_dA(cnt2,4) + PDrot(cnt1,5)*PDrot_dA(cnt2,12) + PDrot(cnt1,12)*PDrot_dA(cnt2,5));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[7]*((PDrot_dA(cnt1,10) + PDrot_dA(cnt1,15))*PDrot(cnt2,4) + PDrot_dA(cnt1,4)*(PDrot(cnt2,10) + PDrot(cnt2,15)) + (PDrot_dA(cnt1,11) - PDrot_dA(cnt1,14))*PDrot(cnt2,5) + PDrot_dA(cnt1,5)*(PDrot(cnt2,11) - PDrot(cnt2,14)));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[7]*((PDrot(cnt1,10) + PDrot(cnt1,15))*PDrot_dA(cnt2,4) + PDrot(cnt1,4)*(PDrot_dA(cnt2,10) + PDrot_dA(cnt2,15)) + (PDrot(cnt1,11) - PDrot(cnt1,14))*PDrot_dA(cnt2,5) + PDrot(cnt1,5)*(PDrot_dA(cnt2,11) - PDrot_dA(cnt2,14)));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[8]*((PDrot_dA(cnt1,10) + PDrot_dA(cnt1,15))*PDrot(cnt2,8) + PDrot_dA(cnt1,8)*(PDrot(cnt2,10) + PDrot(cnt2,15)) + (PDrot_dA(cnt1,11) - PDrot_dA(cnt1,14))*PDrot(cnt2,12) + PDrot_dA(cnt1,12)*(PDrot(cnt2,11) - PDrot(cnt2,14)));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[8]*((PDrot(cnt1,10) + PDrot(cnt1,15))*PDrot_dA(cnt2,8) + PDrot(cnt1,8)*(PDrot_dA(cnt2,10) + PDrot_dA(cnt2,15)) + (PDrot(cnt1,11) - PDrot(cnt1,14))*PDrot_dA(cnt2,12) + PDrot(cnt1,12)*(PDrot_dA(cnt2,11) - PDrot_dA(cnt2,14)));
              }
            }
          }
        }
        //(pd|dd) integrals
        Dd[0] = Dvalue(atmC,5);
        Dd[1] = Dvalue(atmD,6);
        intn[0] = d_eri2Center(1,0,2,0,2,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(1,0,2,0,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center(1,1,2,1,2,0,2,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[3] = d_eri2Center(1,1,2,1,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(pi,pi|pi,pi)=(-pi,-pi|-pi,-pi)=(pi,pi|-pi,-pi)=(-pi,-pi|pi,pi)
        intn[4] = d_eri2Center(1,0,2,0,2,2,2,2,RCD,Dd,atmC,atmD,type);    //(sg,sg|dl,dl)=(sg,sg|-dl,-dl)
        intn[5] = d_eri2Center(1,1,2,1,2,2,2,2,RCD,Dd,atmC,atmD,type);    //(pi,pi|dl,dl)=(-pi,-pi|-dl,-dl)=(pi,pi|-dl,-dl)=(-pi,-pi|dl,dl)
        intn[6] = d_eri2Center(1,0,2,1,2,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[7] = d_eri2Center(1,1,2,0,2,1,2,0,RCD,Dd,atmC,atmD,type);    //(pi,sg|pi,sg)=(-pi,sg|-pi,sg)
        intn[8] = d_eri2Center(1,1,2,2,2,1,2,2,RCD,Dd,atmC,atmD,type);    //(pi,dl|pi,dl)=(pi,dl|-pi,-dl)=(-pi,-dl|-pi,-dl)=(-pi,-dl|pi,dl)=(-pi,dl|-pi,dl)=-(-pi,dl|pi,-dl)=(pi,-dl|pi,-dl)=-(pi,-dl|-pi,dl)
        intn[9] = d_eri2Center(1,1,2,0,2,1,2,2,RCD,Dd,atmC,atmD,type);    //(pi,sg|pi,dl)=(pi,sg|-pi,-dl)=-(-pi,sg|-pi,dl)=(-pi,sg|pi,-dl)
        intn[10] = d_eri2Center(1,0,2,1,2,1,2,2,RCD,Dd,atmC,atmD,type);   //(sg,pi|pi,dl)=(sg,pi|-pi,-dl)=-(sg,-pi|-pi,dl)=(sg,-pi|pi,-dl)
        intn[11] = d_eri2Center(1,1,2,2,2,0,2,1,RCD,Dd,atmC,atmD,type);   //(pi,dl|sg,pi)=(-pi,-dl|sg,pi)=-(-pi,dl|sg,-pi)=(pi,-dl|sg,-pi)
        //(pd|dd) rotate
        for (size_t iddd = 1; iddd < 16; ++iddd) {
          for (size_t idpd1 = 1; idpd1 < 6; ++idpd1) {
            for (size_t idpd2 = 1; idpd2 < 4; ++idpd2) {
              cnt1 = (idpd1 - 1)*3 + idpd2;
              block(5*(idpd2+2) + idpd1,30 + iddd) = -intn[0]*PDrot_dA(cnt1,1)*DDrot(iddd,1) - intn[1]*PDrot_dA(cnt1,1)*(DDrot(iddd,2) + DDrot(iddd,3));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[0]*PDrot(cnt1,1)*DDrot_dA(iddd,1) + intn[1]*PDrot(cnt1,1)*(DDrot_dA(iddd,2) + DDrot_dA(iddd,3));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[2]*(PDrot_dA(cnt1,2) + PDrot_dA(cnt1,3))*DDrot(iddd,1) + intn[3]*(PDrot_dA(cnt1,2) + PDrot_dA(cnt1,3))*(DDrot(iddd,2) + DDrot(iddd,3));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[2]*(PDrot(cnt1,2) + PDrot(cnt1,3))*DDrot_dA(iddd,1) + intn[3]*(PDrot(cnt1,2) + PDrot(cnt1,3))*(DDrot_dA(iddd,2) + DDrot_dA(iddd,3));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[4]*PDrot_dA(cnt1,1)*(DDrot(iddd,4) + DDrot(iddd,5)) + intn[5]*(PDrot_dA(cnt1,2) + PDrot_dA(cnt1,3))*(DDrot(iddd,4) + DDrot(iddd,5));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[4]*PDrot(cnt1,1)*(DDrot_dA(iddd,4) + DDrot_dA(iddd,5)) + intn[5]*(PDrot(cnt1,2) + PDrot(cnt1,3))*(DDrot_dA(iddd,4) + DDrot_dA(iddd,5));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[6]*(PDrot_dA(cnt1,4)*DDrot(iddd,6) + PDrot_dA(cnt1,5)*DDrot(iddd,7)) + intn[7]*(PDrot_dA(cnt1,8)*DDrot(iddd,6) + PDrot_dA(cnt1,12)*DDrot(iddd,7));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[6]*(PDrot(cnt1,4)*DDrot_dA(iddd,6) + PDrot(cnt1,5)*DDrot_dA(iddd,7)) + intn[7]*(PDrot(cnt1,8)*DDrot_dA(iddd,6) + PDrot(cnt1,12)*DDrot_dA(iddd,7));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[8]*((PDrot_dA(cnt1,10) + PDrot_dA(cnt1,15))*(DDrot(iddd,11) + DDrot(iddd,14)) + (PDrot_dA(cnt1,14) - PDrot_dA(cnt1,11))*(DDrot(iddd,13) - DDrot(iddd,12)));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[8]*((PDrot(cnt1,10) + PDrot(cnt1,15))*(DDrot_dA(iddd,11) + DDrot_dA(iddd,14)) + (PDrot(cnt1,14) - PDrot(cnt1,11))*(DDrot_dA(iddd,13) - DDrot_dA(iddd,12)));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[9]*(PDrot_dA(cnt1,8)*(DDrot(iddd,11) + DDrot(iddd,14)) + PDrot_dA(cnt1,12)*(DDrot(iddd,12) - DDrot(iddd,13)));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[9]*(PDrot(cnt1,8)*(DDrot_dA(iddd,11) + DDrot_dA(iddd,14)) + PDrot(cnt1,12)*(DDrot_dA(iddd,12) - DDrot_dA(iddd,13)));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[10]*(PDrot_dA(cnt1,4)*(DDrot(iddd,11) + DDrot(iddd,14)) + PDrot(cnt1,5)*(DDrot_dA(iddd,12) - DDrot_dA(iddd,13)));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[10]*(PDrot(cnt1,4)*(DDrot_dA(iddd,11) + DDrot_dA(iddd,14)) + PDrot_dA(cnt1,5)*(DDrot(iddd,12) - DDrot(iddd,13)));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[11]*((PDrot_dA(cnt1,10) + PDrot_dA(cnt1,15))*DDrot(iddd,6) + (PDrot(cnt1,11) - PDrot(cnt1,14))*DDrot_dA(iddd,7));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[11]*((PDrot(cnt1,10) + PDrot(cnt1,15))*DDrot_dA(iddd,6) + (PDrot_dA(cnt1,11) - PDrot_dA(cnt1,14))*DDrot(iddd,7));
            }
          }
        }
        //(dd|pd) integrals
        Dd[0] = Dvalue(atmC,6);
        Dd[1] = Dvalue(atmD,5);
        intn[0] = d_eri2Center(2,0,2,0,1,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(2,1,2,1,1,0,2,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[2] = d_eri2Center(2,0,2,0,1,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[3] = d_eri2Center(2,1,2,1,1,1,2,1,RCD,Dd,atmC,atmD,type);    //(pi,pi|pi,pi)=(-pi,-pi|-pi,-pi)=(pi,pi|-pi,-pi)=(-pi,-pi|pi,pi)
        intn[4] = d_eri2Center(2,2,2,2,1,0,2,0,RCD,Dd,atmC,atmD,type);    //(dl,dl|sg,sg)=(-dl,-dl|sg,sg)
        intn[5] = d_eri2Center(2,2,2,2,1,1,2,1,RCD,Dd,atmC,atmD,type);    //(dl,dl|pi,pi)=(-dl,-dl|-pi,-pi)=(-dl,-dl|pi,pi)=(dl,dl|-pi,-pi)
        intn[6] = d_eri2Center(2,0,2,1,1,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[7] = d_eri2Center(2,1,2,0,1,1,2,0,RCD,Dd,atmC,atmD,type);    //(pi,sg|pi,sg)=(-pi,sg|-pi,sg)
        intn[8] = d_eri2Center(2,1,2,2,1,1,2,2,RCD,Dd,atmC,atmD,type);    //(pi,dl|pi,dl)=(-pi,-dl|-pi,-dl)=(-pi,-dl|pi,dl)=(pi,dl|-pi,-dl)=(-pi,dl|-pi,dl)=(pi,-dl|pi,-dl)=-(-pi,dl|pi,-dl)=-(pi,-dl|-pi,dl)
        intn[9] = d_eri2Center(2,1,2,2,1,1,2,0,RCD,Dd,atmC,atmD,type);    //(pi,dl|pi,sg)=-(-pi,dl|-pi,sg)=(pi,-dl|-pi,sg)=(-pi,-dl|pi,sg)
        intn[10] = d_eri2Center(2,1,2,2,1,0,2,1,RCD,Dd,atmC,atmD,type);   //(pi,dl|sg,pi)=-(-pi,dl|sg,-pi)=(pi,-dl|sg,-pi)=(-pi,-dl|sg,pi)
        intn[11] = d_eri2Center(2,0,2,1,1,1,2,2,RCD,Dd,atmC,atmD,type);   //(sg,pi|pi,dl)=(sg,pi|-pi,-dl)=-(sg,-pi|-pi,dl)=(sg,-pi|pi,-dl)
        //(dd|pd) rotate
        for (size_t iddd = 1; iddd < 16; ++iddd) {
          for (size_t idpd1 = 1; idpd1 < 6; ++idpd1) {
            for (size_t idpd2 = 1; idpd2 < 4; ++idpd2) {
              cnt1 = (idpd1 - 1)*3 + idpd2;
              block(30 + iddd,5*(idpd2+2) + idpd1) = -intn[0]*PDrot_dA(cnt1,1)*DDrot(iddd,1) - intn[1]*PDrot_dA(cnt1,1)*(DDrot(iddd,2) + DDrot(iddd,3));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[0]*PDrot(cnt1,1)*DDrot_dA(iddd,1) + intn[1]*PDrot(cnt1,1)*(DDrot_dA(iddd,2) + DDrot_dA(iddd,3));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[2]*(PDrot_dA(cnt1,2) + PDrot_dA(cnt1,3))*DDrot(iddd,1) + intn[3]*(PDrot(cnt1,2) + PDrot(cnt1,3))*(DDrot_dA(iddd,2) + DDrot_dA(iddd,3));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[2]*(PDrot(cnt1,2) + PDrot(cnt1,3))*DDrot_dA(iddd,1) + intn[3]*(PDrot_dA(cnt1,2) + PDrot_dA(cnt1,3))*(DDrot(iddd,2) + DDrot(iddd,3));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[4]*PDrot_dA(cnt1,1)*(DDrot(iddd,4) + DDrot(iddd,5)) + intn[5]*(PDrot_dA(cnt1,2) + PDrot_dA(cnt1,3))*(DDrot(iddd,4) + DDrot(iddd,5));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[4]*PDrot(cnt1,1)*(DDrot_dA(iddd,4) + DDrot_dA(iddd,5)) + intn[5]*(PDrot(cnt1,2) + PDrot(cnt1,3))*(DDrot_dA(iddd,4) + DDrot_dA(iddd,5));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[6]*(PDrot_dA(cnt1,4)*DDrot(iddd,6) + PDrot_dA(cnt1,5)*DDrot(iddd,7)) + intn[7]*(PDrot_dA(cnt1,8)*DDrot(iddd,6) + PDrot_dA(cnt1,12)*DDrot(iddd,7));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[6]*(PDrot(cnt1,4)*DDrot_dA(iddd,6) + PDrot(cnt1,5)*DDrot_dA(iddd,7)) + intn[7]*(PDrot(cnt1,8)*DDrot_dA(iddd,6) + PDrot(cnt1,12)*DDrot_dA(iddd,7));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[8]*((PDrot_dA(cnt1,10) + PDrot_dA(cnt1,15))*(DDrot(iddd,11) + DDrot(iddd,14)) + (PDrot_dA(cnt1,14) - PDrot_dA(cnt1,11))*(DDrot(iddd,13) - DDrot(iddd,12)));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[8]*((PDrot(cnt1,10) + PDrot(cnt1,15))*(DDrot_dA(iddd,11) + DDrot_dA(iddd,14)) + (PDrot(cnt1,14) - PDrot(cnt1,11))*(DDrot_dA(iddd,13) - DDrot_dA(iddd,12)));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[9]*(PDrot_dA(cnt1,8)*(DDrot(iddd,11) + DDrot(iddd,14)) + PDrot(cnt1,12)*(DDrot_dA(iddd,12) - DDrot_dA(iddd,13)));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[9]*(PDrot(cnt1,8)*(DDrot_dA(iddd,11) + DDrot_dA(iddd,14)) + PDrot_dA(cnt1,12)*(DDrot(iddd,12) - DDrot(iddd,13)));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[10]*(PDrot_dA(cnt1,4)*(DDrot(iddd,11) + DDrot(iddd,14)) + PDrot(cnt1,5)*(DDrot_dA(iddd,12) - DDrot_dA(iddd,13)));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[10]*(PDrot(cnt1,4)*(DDrot_dA(iddd,11) + DDrot_dA(iddd,14)) + PDrot_dA(cnt1,5)*(DDrot(iddd,12) - DDrot(iddd,13)));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[11]*((PDrot_dA(cnt1,10) + PDrot_dA(cnt1,15))*DDrot(iddd,6) + (PDrot(cnt1,11) - PDrot(cnt1,14))*DDrot_dA(iddd,7));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[11]*((PDrot(cnt1,10) + PDrot(cnt1,15))*DDrot_dA(iddd,6) + (PDrot_dA(cnt1,11) - PDrot_dA(cnt1,14))*DDrot(iddd,7));
            }
          }
        }
        //(dd|dd) integrals
        Dd[0] = Dvalue(atmC,6);
        Dd[1] = Dvalue(atmD,6);
        intn[0] = d_eri2Center(2,0,2,0,2,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(2,0,2,0,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center(2,0,2,0,2,2,2,2,RCD,Dd,atmC,atmD,type);    //(sg,sg|dl,dl)=(sg,sg|-dl,-dl)
        intn[3] = d_eri2Center(2,1,2,1,2,0,2,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[4] = d_eri2Center(2,2,2,2,2,0,2,0,RCD,Dd,atmC,atmD,type);    //(dl,dl|sg,sg)=(-dl,-dl|sg,sg)
        intn[5] = d_eri2Center(2,2,2,2,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(dl,dl|pi,pi)=(-dl,-dl|pi,pi)=(dl,dl|-pi,-pi)=(-dl,-dl|-pi,-pi)
        intn[6] = d_eri2Center(2,1,2,1,2,2,2,2,RCD,Dd,atmC,atmD,type);    //(pi,pi|dl,dl)=(-pi,-pi|dl,dl)=(pi,pi|-dl,-dl)=(-pi,-pi|-dl,-dl)
        intn[7] = d_eri2Center(2,1,2,1,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(pi,pi|pi,pi)=(-pi,-pi|-pi,-pi)
        intn[8] = d_eri2Center(2,2,2,2,2,2,2,2,RCD,Dd,atmC,atmD,type);    //(dl,dl|dl,dl)=(-dl,-dl|-dl,-dl)
        intn[9] = d_eri2Center(2,1,2,1,2,-1,2,-1,RCD,Dd,atmC,atmD,type);  //(pi,pi|-pi,-pi)=(-pi,-pi|pi,pi)
        intn[10] = d_eri2Center(2,-2,2,-2,2,2,2,2,RCD,Dd,atmC,atmD,type); //(-dl,-dl|dl,dl)=(dl,dl|-dl,-dl)
        intn[11] = d_eri2Center(2,0,2,1,2,0,2,1,RCD,Dd,atmC,atmD,type);   //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[12] = d_eri2Center(2,0,2,2,2,0,2,2,RCD,Dd,atmC,atmD,type);   //(sg,dl|sg,dl)=(sg,-dl|sg,-dl)
        intn[13] = d_eri2Center(2,1,2,2,2,1,2,2,RCD,Dd,atmC,atmD,type);   //(pi,dl|pi,dl)=(pi,dl|-pi,-dl)=(-pi,-dl|pi,dl)=(-pi,-dl|-pi,-dl)=-(-pi,dl|pi,-dl)=-(pi,-dl|-pi,dl)=(-pi,dl|-pi,dl)=(pi,-dl|pi,-dl)
        intn[14] = d_eri2Center(2,1,2,-1,2,1,2,-1,RCD,Dd,atmC,atmD,type); //(pi,-pi|pi,-pi)
        intn[15] = d_eri2Center(2,1,2,-2,2,0,2,-1,RCD,Dd,atmC,atmD,type); //(pi,dl|sg,pi)=(-pi,-dl|sg,pi)=(pi,-dl|sg,-pi)=-(-pi,dl|sg,-pi)=(sg,pi|pi,dl)=(sg,pi|-pi,-dl)=(sg,-pi|pi,-dl)=-(sg,-pi|-pi,dl)
        intn[16] = d_eri2Center(2,1,2,1,2,0,2,2,RCD,Dd,atmC,atmD,type);   //(pi,pi|sg,dl)=-(-pi,-pi|sg,dl)=(pi,-pi|sg,-dl)=(sg,dl|pi,pi)=-(sg,dl|-pi,-pi)=(sg,-dl|pi,-pi)
        //(dd|dd) rotate
        for (size_t iddd1 = 1; iddd1 < 16; ++iddd1) {
          for (size_t iddd2 = 1; iddd2 < 16; ++iddd2) {
            block(30 + iddd1,30 + iddd2) = intn[0]*DDrot_dA(iddd1,1)*DDrot(iddd2,1) + intn[1]*DDrot_dA(iddd1,1)*(DDrot(iddd2,2) + DDrot(iddd2,3)) + intn[2]*DDrot(iddd1,1)*(DDrot_dA(iddd2,4) + DDrot_dA(iddd2,5));
            block(30 + iddd1,30 + iddd2) += intn[0]*DDrot(iddd1,1)*DDrot_dA(iddd2,1) + intn[1]*DDrot(iddd1,1)*(DDrot_dA(iddd2,2) + DDrot_dA(iddd2,3)) + intn[2]*DDrot_dA(iddd1,1)*(DDrot(iddd2,4) + DDrot(iddd2,5));
            block(30 + iddd1,30 + iddd2) += intn[3]*(DDrot_dA(iddd1,2) + DDrot_dA(iddd1,3))*DDrot(iddd2,1) + intn[4]*(DDrot(iddd1,4) + DDrot(iddd1,5))*DDrot_dA(iddd2,1);
            block(30 + iddd1,30 + iddd2) += intn[3]*(DDrot(iddd1,2) + DDrot(iddd1,3))*DDrot_dA(iddd2,1) + intn[4]*(DDrot_dA(iddd1,4) + DDrot_dA(iddd1,5))*DDrot(iddd2,1);
            block(30 + iddd1,30 + iddd2) += intn[5]*(DDrot_dA(iddd1,4) + DDrot_dA(iddd1,5))*(DDrot(iddd2,2) + DDrot(iddd2,3)) + intn[6]*(DDrot_dA(iddd1,2) + DDrot_dA(iddd1,3))*(DDrot(iddd2,4) + DDrot(iddd2,5));
            block(30 + iddd1,30 + iddd2) += intn[5]*(DDrot(iddd1,4) + DDrot(iddd1,5))*(DDrot_dA(iddd2,2) + DDrot_dA(iddd2,3)) + intn[6]*(DDrot(iddd1,2) + DDrot(iddd1,3))*(DDrot_dA(iddd2,4) + DDrot_dA(iddd2,5));
            block(30 + iddd1,30 + iddd2) += intn[7]*(DDrot_dA(iddd1,2)*DDrot(iddd2,2) + DDrot_dA(iddd1,3)*DDrot(iddd2,3)) + intn[8]*(DDrot_dA(iddd1,4)*DDrot(iddd2,4) + DDrot_dA(iddd1,5)*DDrot(iddd2,5));
            block(30 + iddd1,30 + iddd2) += intn[7]*(DDrot(iddd1,2)*DDrot_dA(iddd2,2) + DDrot(iddd1,3)*DDrot_dA(iddd2,3)) + intn[8]*(DDrot(iddd1,4)*DDrot_dA(iddd2,4) + DDrot(iddd1,5)*DDrot_dA(iddd2,5));
            block(30 + iddd1,30 + iddd2) += intn[9]*(DDrot_dA(iddd1,2)*DDrot(iddd2,3) + DDrot_dA(iddd1,3)*DDrot(iddd2,2)) + intn[10]*(DDrot_dA(iddd1,4)*DDrot(iddd2,5) + DDrot_dA(iddd1,5)*DDrot(iddd2,4));
            block(30 + iddd1,30 + iddd2) += intn[9]*(DDrot(iddd1,2)*DDrot_dA(iddd2,3) + DDrot(iddd1,3)*DDrot_dA(iddd2,2)) + intn[10]*(DDrot(iddd1,4)*DDrot_dA(iddd2,5) + DDrot(iddd1,5)*DDrot_dA(iddd2,4));
            block(30 + iddd1,30 + iddd2) += intn[11]*(DDrot_dA(iddd1,6)*DDrot(iddd2,6) + DDrot_dA(iddd1,7)*DDrot(iddd2,7)) + intn[12]*(DDrot_dA(iddd1,8)*DDrot(iddd2,8) + DDrot_dA(iddd1,9)*DDrot(iddd2,9));
            block(30 + iddd1,30 + iddd2) += intn[11]*(DDrot(iddd1,6)*DDrot_dA(iddd2,6) + DDrot(iddd1,7)*DDrot_dA(iddd2,7)) + intn[12]*(DDrot(iddd1,8)*DDrot_dA(iddd2,8) + DDrot(iddd1,9)*DDrot_dA(iddd2,9));
            block(30 + iddd1,30 + iddd2) += intn[13]*((DDrot_dA(iddd1,11) + DDrot_dA(iddd1,14))*(DDrot(iddd2,11) + DDrot(iddd2,14)) + (DDrot_dA(iddd1,12) - DDrot_dA(iddd1,13))*(DDrot(iddd2,12) - DDrot(iddd2,13)));
            block(30 + iddd1,30 + iddd2) += intn[13]*((DDrot(iddd1,11) + DDrot(iddd1,14))*(DDrot_dA(iddd2,11) + DDrot_dA(iddd2,14)) + (DDrot(iddd1,12) - DDrot(iddd1,13))*(DDrot_dA(iddd2,12) - DDrot_dA(iddd2,13)));
            block(30 + iddd1,30 + iddd2) += intn[14]*DDrot_dA(iddd1,10)*DDrot(iddd2,10) + intn[14]*DDrot(iddd1,10)*DDrot_dA(iddd2,10);
            block(30 + iddd1,30 + iddd2) += intn[15]*((DDrot_dA(iddd1,11) + DDrot_dA(iddd1,14))*DDrot(iddd2,6) + (DDrot(iddd1,12) - DDrot(iddd1,13))*DDrot_dA(iddd2,7) + DDrot_dA(iddd1,6)*(DDrot(iddd2,11) + DDrot(iddd2,14)) + DDrot(iddd1,7)*(DDrot_dA(iddd2,12) - DDrot_dA(iddd2,13)));
            block(30 + iddd1,30 + iddd2) += intn[15]*((DDrot(iddd1,11) + DDrot(iddd1,14))*DDrot_dA(iddd2,6) + (DDrot_dA(iddd1,12) - DDrot_dA(iddd1,13))*DDrot(iddd2,7) + DDrot(iddd1,6)*(DDrot_dA(iddd2,11) + DDrot_dA(iddd2,14)) + DDrot_dA(iddd1,7)*(DDrot(iddd2,12) - DDrot(iddd2,13)));
            block(30 + iddd1,30 + iddd2) += intn[16]*((DDrot_dA(iddd1,2) - DDrot_dA(iddd1,3))*DDrot(iddd2,8) + DDrot_dA(iddd1,10)*DDrot(iddd2,9) + DDrot_dA(iddd1,8)*(DDrot(iddd2,2) - DDrot(iddd2,3)) + DDrot_dA(iddd1,9)*DDrot(iddd2,10));
            block(30 + iddd1,30 + iddd2) += intn[16]*((DDrot(iddd1,2) - DDrot(iddd1,3))*DDrot_dA(iddd2,8) + DDrot(iddd1,10)*DDrot_dA(iddd2,9) + DDrot(iddd1,8)*(DDrot_dA(iddd2,2) - DDrot_dA(iddd2,3)) + DDrot(iddd1,9)*DDrot_dA(iddd2,10));
          }
        }
      }
    }
  }
  void IntegralBlock2C_dR2(int type, matrixE & block, int atmC, int atmD, double RCD, double cost, double sint, double cosp, double sinp) {
    //function calculating a block of second-derivatives of two-center integrals with respect to the internuclear distance
    int cnt1;
    int cnt2;
    size_t nrows = block.rows();
    size_t ncols = block.cols();
    if ((nrows > 1)||(ncols > 1)) {
      SProt = SPtransf(cost,sint,cosp,sinp);
      PProt = PPtransf(cost,sint,cosp,sinp);
      if ((nrows > 10)||(ncols > 10)) {
        SDrot = SDtransf(cost,sint,cosp,sinp);
        PDrot = PDtransf(cost,sint,cosp,sinp);
        DDrot = DDtransf(cost,sint,cosp,sinp);
      }
    }
    D[0] = Dvalue(atmC,1);
    D[1] = Dvalue(atmC,2);
    D[2] = Dvalue(atmD,1);
    D[3] = Dvalue(atmD,2);
    //integral calculation
    block(1,1) = eri2Center_dR2(0,0,0,0,0,0,0,0,RCD,D,atmC,atmD,type);       //(ss|ss)
    if (nrows > 1) {
      intn[0] = eri2Center_dR2(0,0,1,0,0,0,0,0,RCD,D,atmC,atmD,type);
      intn[1] = eri2Center_dR2(1,0,1,0,0,0,0,0,RCD,D,atmC,atmD,type);
      intn[2] = eri2Center_dR2(1,1,1,1,0,0,0,0,RCD,D,atmC,atmD,type);
      for (size_t idx = 0; idx < 6; ++idx) {
        if (idx < 3) {block(idx + 2,1) = -intn[0]*SProt(idx + 1,1);}                                                   //(sp|ss) integrals
        block(5 + idx,1) = intn[1]*PProt(idx + 1,1) + intn[2]*(PProt(idx + 1,2) + PProt(idx + 1,3));   //(pp|ss) integrals
      }
      if (nrows > 10) {
        //(sd|ss) integrals
        Dd[0] = Dvalue(atmC,4);
        Dd[1] = Dvalue(atmD,1);
        intn[0] = d_eri2Center_dR2(0,0,2,0,0,0,0,0,RCD,Dd,atmC,atmD,type);
        for (size_t idsd = 1; idsd < 6; ++idsd) {
          block(10 + idsd,1) = intn[0]*SDrot(idsd,1);
        }
        //(pd|ss) integrals
        Dd[0] = Dvalue(atmC,5);
        intn[0] = d_eri2Center_dR2(1,0,2,0,0,0,0,0,RCD,Dd,atmC,atmD,type);
        intn[1] = d_eri2Center_dR2(1,1,2,1,0,0,0,0,RCD,Dd,atmC,atmD,type);
        for (int idpd1 = 1; idpd1 < 6; ++idpd1) {
          for (int idpd2 = 1; idpd2 < 4; ++idpd2) {
            block(5*(idpd2+2) + idpd1,1) = -(intn[0]*PDrot((idpd1 - 1)*3 + idpd2,1) + intn[1]*(PDrot((idpd1 - 1)*3 + idpd2,2) + PDrot((idpd1 - 1)*3 + idpd2,3)));
          }
        }
        //(dd|ss) integrals
        Dd[0] = Dvalue(atmC,6);
        intn[0] = d_eri2Center_dR2(2,0,2,0,0,0,0,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR2(2,1,2,1,0,0,0,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[2] = d_eri2Center_dR2(2,2,2,2,0,0,0,0,RCD,Dd,atmC,atmD,type);    //(dl,dl|sg,sg)=(-dl,-dl|sg,sg)
        for (size_t iddd = 1; iddd < 16; ++iddd) {
          block(30 + iddd,1) = intn[0]*DDrot(iddd,1) + intn[1]*(DDrot(iddd,2) + DDrot(iddd,3)) + intn[2]*(DDrot(iddd,4) + DDrot(iddd,5));
        }
      }
    }
    if ((ncols > 1)&&(type == 0)) {
      intn[0] = eri2Center_dR2(0,0,0,0,0,0,1,0,RCD,D,atmC,atmD,type);
      intn[1] = eri2Center_dR2(0,0,0,0,1,0,1,0,RCD,D,atmC,atmD,type);
      intn[2] = eri2Center_dR2(0,0,0,0,1,1,1,1,RCD,D,atmC,atmD,type);
      for (size_t idx = 0; idx < 6; ++idx) {
        if (idx < 3) {block(1,idx + 2) = -intn[0]*SProt(idx + 1,1);}                                                   //(ss|sp) integrals
        block(1,5 + idx) = intn[1]*PProt(idx + 1,1) + intn[2]*(PProt(idx + 1,2) + PProt(idx + 1,3));   //(ss|pp) integrals
      }
      if (ncols > 10) {
        //(ss|sd) integrals
        Dd[0] = Dvalue(atmC,1);
        Dd[1] = Dvalue(atmD,4);
        intn[0] = d_eri2Center_dR2(0,0,0,0,0,0,2,0,RCD,Dd,atmC,atmD,type);
        for (size_t idsd = 1; idsd < 6; ++idsd) {
          block(1,10 + idsd) = intn[0]*SDrot(idsd,1);
        }
        //(ss|pd) integrals
        Dd[1] = Dvalue(atmD,5);
        intn[0] = d_eri2Center_dR2(0,0,0,0,1,0,2,0,RCD,Dd,atmC,atmD,type);
        intn[1] = d_eri2Center_dR2(0,0,0,0,1,1,2,1,RCD,Dd,atmC,atmD,type);
        for (int idpd1 = 1; idpd1 < 6; ++idpd1) {
          for (int idpd2 = 1; idpd2 < 4; ++idpd2) {
            block(1,5*(idpd2 + 2) + idpd1) = -intn[0]*PDrot((idpd1 - 1)*3 + idpd2,1) - intn[1]*(PDrot((idpd1 - 1)*3 + idpd2,2) + PDrot((idpd1 - 1)*3 + idpd2,3));
          }
        }
        //(ss|dd) integrals
        Dd[1] = Dvalue(atmD,6);
        intn[0] = d_eri2Center_dR2(0,0,0,0,2,0,2,0,RCD,Dd,atmC,atmD,type);
        intn[1] = d_eri2Center_dR2(0,0,0,0,2,1,2,1,RCD,Dd,atmC,atmD,type);
        intn[2] = d_eri2Center_dR2(0,0,0,0,2,2,2,2,RCD,Dd,atmC,atmD,type);
        for (size_t iddd = 1; iddd < 16; ++iddd) {
          block(1,30 + iddd) = intn[0]*DDrot(iddd,1) + intn[1]*(DDrot(iddd,2) + DDrot(iddd,3)) + intn[2]*(DDrot(iddd,4) + DDrot(iddd,5));
        }
      }
    }
    if ((nrows > 1)&&(ncols > 1)&&(type == 0)) {
      intn[0] = eri2Center_dR2(0,0,1,0,0,0,1,0,RCD,D,atmC,atmD,type);
      intn[1] = eri2Center_dR2(0,0,1,1,0,0,1,1,RCD,D,atmC,atmD,type);
      intn[2] = eri2Center_dR2(0,0,1,0,1,0,1,0,RCD,D,atmC,atmD,type);
      intn[3] = eri2Center_dR2(0,0,1,0,1,1,1,1,RCD,D,atmC,atmD,type);
      intn[4] = eri2Center_dR2(0,0,1,1,1,0,1,1,RCD,D,atmC,atmD,type);
      for (size_t idpb = 1; idpb < 4; ++idpb) {
        for (size_t idpk = 1; idpk < 4; ++idpk) {           //(sp|sp)
          block(idpb + 1,idpk + 1) = intn[0]*SProt(idpb,1)*SProt(idpk,1) + intn[1]*(SProt(idpb,2)*SProt(idpk,2) + SProt(idpb,3)*SProt(idpk,3));
        }
        for (size_t idc = 1; idc < 7; ++idc) {              //(sp|pp)
          block(idpb + 1,idc + 4) = -SProt(idpb,1)*(intn[2]*PProt(idc,1) + intn[3]*(PProt(idc,2) + PProt(idc,3))) - intn[4]*(SProt(idpb,2)*PProt(idc,4) + SProt(idpb,3)*PProt(idc,5));
        }
      }
      intn[0] = eri2Center_dR2(1,0,1,0,1,0,1,0,RCD,D,atmC,atmD,type);                          //(ps,ps|ps,ps)
      intn[1] = eri2Center_dR2(1,1,1,1,1,0,1,0,RCD,D,atmC,atmD,type);                          //(pp,pp|ps,ps)
      intn[2] = eri2Center_dR2(1,0,1,0,1,1,1,1,RCD,D,atmC,atmD,type);                          //(ps,ps|pp,pp)
      intn[3] = eri2Center_dR2(1,1,1,1,1,1,1,1,RCD,D,atmC,atmD,type);                          //(pp,pp|pp,pp)
      intn[4] = eri2Center_dR2(1,1,1,0,1,1,1,0,RCD,D,atmC,atmD,type);                          //(pp,ps|pp,ps)
      intn[5] = eri2Center_dR2(1,1,1,1,1,-1,1,-1,RCD,D,atmC,atmD,type);                        //(pp,pp|pp*,pp*)
      intn[6] = eri2Center_dR2(1,-1,1,1,1,-1,1,1,RCD,D,atmC,atmD,type);                        //(pp*,pp|pp*,pp)
      intn[7] = eri2Center_dR2(1,0,1,0,0,0,1,0,RCD,D,atmC,atmD,type);
      intn[8] = eri2Center_dR2(1,1,1,1,0,0,1,0,RCD,D,atmC,atmD,type);
      intn[9] = eri2Center_dR2(1,0,1,1,0,0,1,1,RCD,D,atmC,atmD,type);
      for (size_t idr = 1; idr < 7; ++idr) {
        for (size_t idc = 1; idc < 4; ++idc) {                //(pp|sp)
          block(idr + 4,idc + 1) = -SProt(idc,1)*(intn[7]*PProt(idr,1) + intn[8]*(PProt(idr,2) + PProt(idr,3))) - intn[9]*(SProt(idc,2)*PProt(idr,4) + SProt(idc,3)*PProt(idr,5));
        }
        for (size_t idc = 1; idc < 7; ++idc) {                //(pp|pp)
          block(idr + 4,idc + 4) = intn[0]*PProt(idr,1)*PProt(idc,1) + intn[1]*(PProt(idr,2) + PProt(idr,3))*PProt(idc,1) + intn[2]*PProt(idr,1)*(PProt(idc,2) + PProt(idc,3)) + intn[3]*(PProt(idr,2)*PProt(idc,2) + PProt(idr,3)*PProt(idc,3)) + intn[4]*(PProt(idr,4)*PProt(idc,4) + PProt(idr,5)*PProt(idc,5)) + intn[5]*(PProt(idr,2)*PProt(idc,3) + PProt(idr,3)*PProt(idc,2)) + intn[6]*PProt(idr,6)*PProt(idc,6);
        }
      }
      if (nrows > 10) {
        //(sd|sp) integrals
        Dd[0] = Dvalue(atmC,4);
        Dd[1] = Dvalue(atmD,1);
        intn[0] = d_eri2Center_dR2(0,0,2,0,0,0,1,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR2(0,0,2,1,0,0,1,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        //(sd|sp) rotate
        for (size_t idsd = 1; idsd < 6; ++idsd) {
          for (size_t idsp = 1; idsp < 4; ++idsp) {
            block(10 + idsd,1 + idsp) = -intn[0]*SDrot(idsd,1)*SProt(idsp,1) - intn[1]*(SDrot(idsd,2)*SProt(idsp,2) + SDrot(idsd,3)*SProt(idsp,3));
          }
        }
        //(sd|pp) integrals
        Dd[0] = Dvalue(atmC,4);
        Dd[1] = Dvalue(atmD,3);
        intn[0] = d_eri2Center_dR2(0,0,2,0,1,0,1,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR2(0,0,2,0,1,1,1,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center_dR2(0,0,2,1,1,0,1,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[3] = d_eri2Center_dR2(0,0,2,2,1,1,1,1,RCD,Dd,atmC,atmD,type);    //(sg,dl|pi,pi)=-(sg,dl|-pi,-pi)=(sg,-dl|pi,-pi)
        //(sd|pp) rotate
        for (size_t idsd = 1; idsd < 6; ++idsd) {
          for (size_t idpp = 1; idpp < 7; ++idpp) {
            block(10 + idsd,4 + idpp) = intn[0]*SDrot(idsd,1)*PProt(idpp,1) + intn[1]*SDrot(idsd,1)*(PProt(idpp,2) + PProt(idpp,3));
            block(10 + idsd,4 + idpp) += intn[2]*(SDrot(idsd,2)*PProt(idpp,4) + SDrot(idsd,3)*PProt(idpp,5));
            block(10 + idsd,4 + idpp) += intn[3]*(SDrot(idsd,4)*PProt(idpp,2) - SDrot(idsd,4)*PProt(idpp,3) + SDrot(idsd,5)*PProt(idpp,6));
          }
        }
        //(pd|sp) integrals
        Dd[0] = Dvalue(atmC,5);
        Dd[1] = Dvalue(atmD,1);
        intn[0] = d_eri2Center_dR2(1,0,2,0,0,0,1,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR2(1,1,2,1,0,0,1,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[2] = d_eri2Center_dR2(1,0,2,1,0,0,1,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[3] = d_eri2Center_dR2(1,1,2,0,0,0,1,1,RCD,Dd,atmC,atmD,type);    //(pi,sg|sg,pi)=(-pi,sg|sg,-pi)
        intn[4] = d_eri2Center_dR2(1,1,2,2,0,0,1,1,RCD,Dd,atmC,atmD,type);    //(pi,dl|sg,pi)=(-pi,-dl|sg,pi)=-(-pi,dl|sg,-pi)=(pi,-dl|sg,-pi)
        //(pd|sp) rotate
        for (size_t idsp = 1; idsp < 4; ++idsp) {
          for (size_t idpd1 = 1; idpd1 < 6; ++idpd1) {
            for (size_t idpd2 = 1; idpd2 < 4; ++idpd2) {
              block(5*(idpd2+2) + idpd1,1 + idsp) = intn[0]*PDrot((idpd1 - 1)*3 + idpd2,1)*SProt(idsp,1) + intn[1]*(PDrot((idpd1 - 1)*3 + idpd2,2) + PDrot((idpd1 - 1)*3 + idpd2,3))*SProt(idsp,1);
              block(5*(idpd2+2) + idpd1,1 + idsp) += intn[2]*(PDrot((idpd1 - 1)*3 + idpd2,4)*SProt(idsp,2) + PDrot((idpd1 - 1)*3 + idpd2,5)*SProt(idsp,3));
              block(5*(idpd2+2) + idpd1,1 + idsp) += intn[3]*(PDrot((idpd1 - 1)*3 + idpd2,8)*SProt(idsp,2) + PDrot((idpd1 - 1)*3 + idpd2,12)*SProt(idsp,3));
              block(5*(idpd2+2) + idpd1,1 + idsp) += intn[4]*((PDrot((idpd1 - 1)*3 + idpd2,10) + PDrot((idpd1 - 1)*3 + idpd2,15))*SProt(idsp,2) + (PDrot((idpd1 - 1)*3 + idpd2,11) - PDrot((idpd1 - 1)*3 + idpd2,14))*SProt(idsp,3));
            }
          }
        }
        //(dd|sp) integrals
        Dd[0] = Dvalue(atmC,6);
        Dd[1] = Dvalue(atmD,1);
        intn[0] = d_eri2Center_dR2(2,0,2,0,0,0,1,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR2(2,1,2,1,0,0,1,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[2] = d_eri2Center_dR2(2,2,2,2,0,0,1,0,RCD,Dd,atmC,atmD,type);    //(dl,dl|sg,sg)=(-dl,-dl|sg,sg)
        intn[3] = d_eri2Center_dR2(2,0,2,1,0,0,1,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[4] = d_eri2Center_dR2(2,1,2,2,0,0,1,1,RCD,Dd,atmC,atmD,type);    //(pi,dl|sg,pi)=(-pi,-dl|sg,pi)=-(-pi,dl|sg,-pi)=(pi,-dl|sg,-pi)
        //(dd|sp) rotate
        for (size_t idsp = 1; idsp < 4; ++idsp) {
          for (size_t iddd = 1; iddd < 16; ++iddd) {
            block(30 + iddd,1 + idsp) = -intn[0]*DDrot(iddd,1)*SProt(idsp,1) - intn[1]*(DDrot(iddd,2) + DDrot(iddd,3))*SProt(idsp,1) - intn[2]*(DDrot(iddd,4) + DDrot(iddd,5))*SProt(idsp,1);
            block(30 + iddd,1 + idsp) -= intn[3]*(DDrot(iddd,6)*SProt(idsp,2) + DDrot(iddd,7)*SProt(idsp,3));
            block(30 + iddd,1 + idsp) -= intn[4]*((DDrot(iddd,11) + DDrot(iddd,14))*SProt(idsp,2) + (DDrot(iddd,12) - DDrot(iddd,13))*SProt(idsp,3));
          }
        }
        //(pd|pp) integrals
        Dd[0] = Dvalue(atmC,5);
        Dd[1] = Dvalue(atmD,3);
        intn[0] = d_eri2Center_dR2(1,0,2,0,1,0,1,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR2(1,1,2,1,1,0,1,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[2] = d_eri2Center_dR2(1,0,2,0,1,1,1,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[3] = d_eri2Center_dR2(1,1,2,1,1,1,1,1,RCD,Dd,atmC,atmD,type);    //(pi,pi|pi,pi)=(pi,pi|-pi,-pi)=(-pi,-pi|pi,pi)=(-pi,-pi|-pi,-pi)
        intn[4] = d_eri2Center_dR2(1,0,2,1,1,0,1,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[5] = d_eri2Center_dR2(1,1,2,0,1,1,1,0,RCD,Dd,atmC,atmD,type);    //(pi,sg|pi,sg)=(-pi,sg|-pi,sg)
        intn[6] = d_eri2Center_dR2(1,1,2,2,1,0,1,1,RCD,Dd,atmC,atmD,type);    //(pi,dl|sg,pi)=(-pi,-dl|sg,pi)=-(-pi,dl|sg,-pi)=(pi,-dl|sg,-pi)
        //(pd|pp) rotate
        for (size_t idpp = 1; idpp < 7; ++idpp) {
          for (size_t idpd1 = 1; idpd1 < 6; ++idpd1) {
            for (size_t idpd2 = 1; idpd2 < 4; ++idpd2) {
              block(5*(idpd2+2) + idpd1,4 + idpp) = -intn[0]*PDrot((idpd1 - 1)*3 + idpd2,1)*PProt(idpp,1) - intn[1]*(PDrot((idpd1 - 1)*3 + idpd2,2) + PDrot((idpd1 - 1)*3 + idpd2,3))*PProt(idpp,1);
              block(5*(idpd2+2) + idpd1,4 + idpp) -= intn[2]*PDrot((idpd1 - 1)*3 + idpd2,1)*(PProt(idpp,2) + PProt(idpp,3));
              block(5*(idpd2+2) + idpd1,4 + idpp) -= intn[3]*(PDrot((idpd1 - 1)*3 + idpd2,2) + PDrot((idpd1 - 1)*3 + idpd2,3))*(PProt(idpp,2) + PProt(idpp,3));
              block(5*(idpd2+2) + idpd1,4 + idpp) -= intn[4]*(PDrot((idpd1 - 1)*3 + idpd2,4)*PProt(idpp,4) + PDrot((idpd1 - 1)*3 + idpd2,5)*PProt(idpp,5));
              block(5*(idpd2+2) + idpd1,4 + idpp) -= intn[5]*(PDrot((idpd1 - 1)*3 + idpd2,8)*PProt(idpp,4) + PDrot((idpd1 - 1)*3 + idpd2,12)*PProt(idpp,5));
              block(5*(idpd2+2) + idpd1,4 + idpp) -= intn[6]*((PDrot((idpd1 - 1)*3 + idpd2,10) + PDrot((idpd1 - 1)*3 + idpd2,15))*PProt(idpp,4) + (PDrot((idpd1 - 1)*3 + idpd2,11) - PDrot((idpd1 - 1)*3 + idpd2,14))*PProt(idpp,5));
            }
          }
        }
        //(dd|pp) integrals
        Dd[0] = Dvalue(atmC,6);
        Dd[1] = Dvalue(atmD,3);
        intn[0] = d_eri2Center_dR2(2,0,2,0,1,0,1,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR2(2,1,2,1,1,0,1,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[2] = d_eri2Center_dR2(2,2,2,2,1,0,1,0,RCD,Dd,atmC,atmD,type);    //(dl,dl|sg,sg)=(-dl,-dl|sg,sg)
        intn[3] = d_eri2Center_dR2(2,0,2,0,1,1,1,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[4] = d_eri2Center_dR2(2,1,2,1,1,1,1,1,RCD,Dd,atmC,atmD,type);    //(pi,pi|pi,pi)=(-pi,-pi|-pi,-pi)
        intn[5] = d_eri2Center_dR2(2,2,2,2,1,1,1,1,RCD,Dd,atmC,atmD,type);    //(dl,dl|pi,pi)=(-dl,-dl|pi,pi)=(dl,dl|-pi,-pi)=(-dl,-dl|-pi,-pi)
        intn[6] = d_eri2Center_dR2(2,0,2,1,1,0,1,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[7] = d_eri2Center_dR2(2,-1,2,-1,1,1,1,1,RCD,Dd,atmC,atmD,type);  //(-pi,-pi|pi,pi)=(pi,pi|-pi,-pi)
        intn[8] = d_eri2Center_dR2(2,1,2,-1,1,1,1,-1,RCD,Dd,atmC,atmD,type);  //(pi,-pi|pi,-pi)
        intn[9] = d_eri2Center_dR2(2,0,2,2,1,1,1,1,RCD,Dd,atmC,atmD,type);    //(sg,dl|pi,pi)=(sg,-dl|pi,-pi)=-(sg,dl|-pi,-pi)
        intn[10] = d_eri2Center_dR2(2,1,2,2,1,0,1,1,RCD,Dd,atmC,atmD,type);   //(pi,dl|sg,pi)=(-pi,-dl|sg,pi)=-(-pi,dl|sg,-pi)=(pi,-dl|sg,-pi)
        //(dd|pp) rotate
        for (size_t iddd = 1; iddd < 16; ++iddd) {
          for (size_t idpp = 1; idpp < 7; ++idpp) {
            block(30 + iddd,4 + idpp) = intn[0]*DDrot(iddd,1)*PProt(idpp,1) + intn[1]*(DDrot(iddd,2) + DDrot(iddd,3))*PProt(idpp,1) + intn[2]*(DDrot(iddd,4) + DDrot(iddd,5))*PProt(idpp,1);
            block(30 + iddd,4 + idpp) += intn[3]*DDrot(iddd,1)*(PProt(idpp,2) + PProt(idpp,3)) + intn[4]*(DDrot(iddd,2)*PProt(idpp,2) + DDrot(iddd,3)*PProt(idpp,3));
            block(30 + iddd,4 + idpp) += intn[5]*(DDrot(iddd,4) + DDrot(iddd,5))*(PProt(idpp,2) + PProt(idpp,3)) + intn[6]*(DDrot(iddd,6)*PProt(idpp,4) + DDrot(iddd,7)*PProt(idpp,5));
            block(30 + iddd,4 + idpp) += intn[7]*(DDrot(iddd,2)*PProt(idpp,3) + DDrot(iddd,3)*PProt(idpp,2)) + intn[8]*DDrot(iddd,10)*PProt(idpp,6);
            block(30 + iddd,4 + idpp) += intn[9]*(DDrot(iddd,8)*(PProt(idpp,2) - PProt(idpp,3)) + DDrot(iddd,9)*PProt(idpp,6));
            block(30 + iddd,4 + idpp) += intn[10]*((DDrot(iddd,11) + DDrot(iddd,14))*PProt(idpp,4) + (DDrot(iddd,12) - DDrot(iddd,13))*PProt(idpp,5));
          }
        }
      }
      if (ncols > 10) {
        //(sp|sd) integrals
        Dd[0] = Dvalue(atmC,1);
        Dd[1] = Dvalue(atmD,4);
        intn[0] = d_eri2Center_dR2(0,0,1,0,0,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR2(0,0,1,1,0,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        //(sp|sd) rotate
        for (size_t idsp = 1; idsp < 4; ++idsp) {
          for (size_t idsd = 1; idsd < 6; ++idsd) {
            block(1 + idsp,10 + idsd) = -intn[0]*SDrot(idsd,1)*SProt(idsp,1) - intn[1]*(SDrot(idsd,2)*SProt(idsp,2) + SDrot(idsd,3)*SProt(idsp,3));
          }
        }
        //(pp|sd) integrals
        Dd[0] = Dvalue(atmC,3);
        Dd[1] = Dvalue(atmD,4);
        intn[0] = d_eri2Center_dR2(1,0,1,0,0,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR2(1,1,1,1,0,0,2,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[2] = d_eri2Center_dR2(1,0,1,1,0,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[3] = d_eri2Center_dR2(1,1,1,1,0,0,2,2,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,dl)=-(-pi,-pi|sg,dl)=(pi,-pi|sg,-dl)
        //(pp|sd) rotate
        for (size_t idsd = 1; idsd < 6; ++idsd) {
          for (size_t idpp = 1; idpp < 7; ++idpp) {
            block(4 + idpp,10 + idsd) = intn[0]*SDrot(idsd,1)*PProt(idpp,1) + intn[1]*SDrot(idsd,1)*(PProt(idpp,2) + PProt(idpp,3));
            block(4 + idpp,10 + idsd) += intn[2]*(SDrot(idsd,2)*PProt(idpp,4) + SDrot(idsd,3)*PProt(idpp,5));
            block(4 + idpp,10 + idsd) += intn[3]*(SDrot(idsd,4)*PProt(idpp,2) - SDrot(idsd,4)*PProt(idpp,3) + SDrot(idsd,5)*PProt(idpp,6));
          }
        }
        //(sp|pd) integrals
        Dd[0] = Dvalue(atmC,1);
        Dd[1] = Dvalue(atmD,5);
        intn[0] = d_eri2Center_dR2(0,0,1,0,1,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR2(0,0,1,0,1,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center_dR2(0,0,1,1,1,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[3] = d_eri2Center_dR2(0,0,1,1,1,1,2,0,RCD,Dd,atmC,atmD,type);    //(sg,pi|pi,sg)=(sg,-pi|-pi,sg)
        intn[4] = d_eri2Center_dR2(0,0,1,1,1,1,2,2,RCD,Dd,atmC,atmD,type);    //(sg,pi|pi,dl)=(sg,pi|-pi,-dl)=-(sg,-pi|-pi,dl)=(sg,-pi|pi,-dl)
        //(sp|pd) rotate
        for (size_t idsp = 1; idsp < 4; ++idsp) {
          for (size_t idpd1 = 1; idpd1 < 6; ++idpd1) {
            for (size_t idpd2 = 1; idpd2 < 4; ++idpd2) {
              block(1 + idsp,5*(idpd2+2) + idpd1) = intn[0]*PDrot((idpd1 - 1)*3 + idpd2,1)*SProt(idsp,1) + intn[1]*(PDrot((idpd1 - 1)*3 + idpd2,2) + PDrot((idpd1 - 1)*3 + idpd2,3))*SProt(idsp,1);
              block(1 + idsp,5*(idpd2+2) + idpd1) += intn[2]*(PDrot((idpd1 - 1)*3 + idpd2,4)*SProt(idsp,2) + PDrot((idpd1 - 1)*3 + idpd2,5)*SProt(idsp,3));
              block(1 + idsp,5*(idpd2+2) + idpd1) += intn[3]*(PDrot((idpd1 - 1)*3 + idpd2,8)*SProt(idsp,2) + PDrot((idpd1 - 1)*3 + idpd2,12)*SProt(idsp,3));
              block(1 + idsp,5*(idpd2+2) + idpd1) += intn[4]*((PDrot((idpd1 - 1)*3 + idpd2,10) + PDrot((idpd1 - 1)*3 + idpd2,15))*SProt(idsp,2) + (PDrot((idpd1 - 1)*3 + idpd2,11) - PDrot((idpd1 - 1)*3 + idpd2,14))*SProt(idsp,3));
            }
          }
        }
        //(sp|dd) integrals
        Dd[0] = Dvalue(atmC,1);
        Dd[1] = Dvalue(atmD,6);
        intn[0] = d_eri2Center_dR2(0,0,1,0,2,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR2(0,0,1,0,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center_dR2(0,0,1,0,2,2,2,2,RCD,Dd,atmC,atmD,type);    //(sg,sg|dl,dl)=(sg,sg|-dl,-dl)
        intn[3] = d_eri2Center_dR2(0,0,1,1,2,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[4] = d_eri2Center_dR2(0,0,1,1,2,1,2,2,RCD,Dd,atmC,atmD,type);    //(sg,pi|pi,dl)=(sg,pi|-pi,-dl)=-(sg,-pi|-pi,dl)=(sg,-pi|pi,-dl)
        //(sp|dd) rotate
        for (size_t idsp = 1; idsp < 4; ++idsp) {
          for (size_t iddd = 1; iddd < 16; ++iddd) {
            block(1 + idsp,30 + iddd) = -intn[0]*DDrot(iddd,1)*SProt(idsp,1) - intn[1]*(DDrot(iddd,2) + DDrot(iddd,3))*SProt(idsp,1) - intn[2]*(DDrot(iddd,4) + DDrot(iddd,5))*SProt(idsp,1);
            block(1 + idsp,30 + iddd) -= intn[3]*(DDrot(iddd,6)*SProt(idsp,2) + DDrot(iddd,7)*SProt(idsp,3));
            block(1 + idsp,30 + iddd) -= intn[4]*((DDrot(iddd,11) + DDrot(iddd,14))*SProt(idsp,2) + (DDrot(iddd,12) - DDrot(iddd,13))*SProt(idsp,3));
          }
        }
        //(pp|pd) integrals
        Dd[0] = Dvalue(atmC,3);
        Dd[1] = Dvalue(atmD,5);
        intn[0] = d_eri2Center_dR2(1,0,1,0,1,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR2(1,0,1,0,1,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center_dR2(1,1,1,1,1,0,2,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[3] = d_eri2Center_dR2(1,1,1,1,1,1,2,1,RCD,Dd,atmC,atmD,type);    //(pi,pi|pi,pi)=(pi,pi|-pi,-pi)=(-pi,-pi|pi,pi)=(-pi,-pi|-pi,-pi)
        intn[4] = d_eri2Center_dR2(1,0,1,1,1,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[5] = d_eri2Center_dR2(1,1,1,0,1,1,2,0,RCD,Dd,atmC,atmD,type);    //(pi,sg|pi,sg)=(-pi,sg|-pi,sg)
        intn[6] = d_eri2Center_dR2(1,0,1,1,1,1,2,2,RCD,Dd,atmC,atmD,type);    //(sg,pi|pi,dl)=(sg,pi|-pi,-dl)=-(sg,-pi|-pi,dl)=(sg,-pi|pi,-dl)
        //(pp|pd) rotate
        for (size_t idpp = 1; idpp < 7; ++idpp) {
          for (size_t idpd1 = 1; idpd1 < 6; ++idpd1) {
            for (size_t idpd2 = 1; idpd2 < 4; ++idpd2) {
              block(4 + idpp,5*(idpd2+2) + idpd1) = -intn[0]*PDrot((idpd1 - 1)*3 + idpd2,1)*PProt(idpp,1) - intn[1]*(PDrot((idpd1 - 1)*3 + idpd2,2) + PDrot((idpd1 - 1)*3 + idpd2,3))*PProt(idpp,1);
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[2]*PDrot((idpd1 - 1)*3 + idpd2,1)*(PProt(idpp,2) + PProt(idpp,3));
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[3]*(PDrot((idpd1 - 1)*3 + idpd2,2) + PDrot((idpd1 - 1)*3 + idpd2,3))*(PProt(idpp,2) + PProt(idpp,3));
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[4]*(PDrot((idpd1 - 1)*3 + idpd2,4)*PProt(idpp,4) + PDrot((idpd1 - 1)*3 + idpd2,5)*PProt(idpp,5));
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[5]*(PDrot((idpd1 - 1)*3 + idpd2,8)*PProt(idpp,4) + PDrot((idpd1 - 1)*3 + idpd2,12)*PProt(idpp,5));
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[6]*((PDrot((idpd1 - 1)*3 + idpd2,10) + PDrot((idpd1 - 1)*3 + idpd2,15))*PProt(idpp,4) + (PDrot((idpd1 - 1)*3 + idpd2,11) - PDrot((idpd1 - 1)*3 + idpd2,14))*PProt(idpp,5));
            }
          }
        }
        //(pp|dd) integrals
        Dd[0] = Dvalue(atmC,3);
        Dd[1] = Dvalue(atmD,6);
        intn[0] = d_eri2Center_dR2(1,0,1,0,2,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR2(1,0,1,0,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center_dR2(1,0,1,0,2,2,2,2,RCD,Dd,atmC,atmD,type);    //(sg,sg|dl,dl)=(sg,sg|-dl,-dl)
        intn[3] = d_eri2Center_dR2(1,1,1,1,2,0,2,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[4] = d_eri2Center_dR2(1,1,1,1,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(pi,pi|pi,pi)=(-pi,-pi|-pi,-pi)
        intn[5] = d_eri2Center_dR2(1,1,1,1,2,2,2,2,RCD,Dd,atmC,atmD,type);    //(pi,pi|dl,dl)=(pi,pi|-dl,-dl)=(-pi,-pi|dl,dl)=(-pi,-pi|-dl,-dl)
        intn[6] = d_eri2Center_dR2(1,0,1,1,2,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[7] = d_eri2Center_dR2(1,-1,1,-1,2,1,2,1,RCD,Dd,atmC,atmD,type);  //(-pi,-pi|pi,pi)=(pi,pi|-pi,-pi)
        intn[8] = d_eri2Center_dR2(1,1,1,-1,2,1,2,-1,RCD,Dd,atmC,atmD,type);  //(pi,-pi|pi,-pi)
        intn[9] = d_eri2Center_dR2(1,1,1,1,2,0,2,2,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,dl)=(pi,-pi|sg,-dl)=-(-pi,-pi|sg,dl)
        intn[10] = d_eri2Center_dR2(1,0,1,1,2,1,2,2,RCD,Dd,atmC,atmD,type);   //(sg,pi|pi,dl)=(sg,pi|-pi,-dl)=-(sg,-pi|-pi,dl)=(sg,-pi|pi,-dl)
        //(pp|dd) rotate
        for (size_t iddd = 1; iddd < 16; ++iddd) {
          for (size_t idpp = 1; idpp < 7; ++idpp) {
            block(4 + idpp,30 + iddd) = intn[0]*DDrot(iddd,1)*PProt(idpp,1) + intn[1]*(DDrot(iddd,2) + DDrot(iddd,3))*PProt(idpp,1) + intn[2]*(DDrot(iddd,4) + DDrot(iddd,5))*PProt(idpp,1);
            block(4 + idpp,30 + iddd) += intn[3]*DDrot(iddd,1)*(PProt(idpp,2) + PProt(idpp,3)) + intn[4]*(DDrot(iddd,2)*PProt(idpp,2) + DDrot(iddd,3)*PProt(idpp,3));
            block(4 + idpp,30 + iddd) += intn[5]*(DDrot(iddd,4) + DDrot(iddd,5))*(PProt(idpp,2) + PProt(idpp,3)) + intn[6]*(DDrot(iddd,6)*PProt(idpp,4) + DDrot(iddd,7)*PProt(idpp,5));
            block(4 + idpp,30 + iddd) += intn[7]*(DDrot(iddd,2)*PProt(idpp,3) + DDrot(iddd,3)*PProt(idpp,2)) + intn[8]*DDrot(iddd,10)*PProt(idpp,6);
            block(4 + idpp,30 + iddd) += intn[9]*(DDrot(iddd,8)*(PProt(idpp,2) - PProt(idpp,3)) + DDrot(iddd,9)*PProt(idpp,6));
            block(4 + idpp,30 + iddd) += intn[10]*((DDrot(iddd,11) + DDrot(iddd,14))*PProt(idpp,4) + (DDrot(iddd,12) - DDrot(iddd,13))*PProt(idpp,5));
          }
        }
      }
      if ((nrows > 10)&&(ncols > 10)) {
        //(sd|sd) integrals
        Dd[0] = Dvalue(atmC,4);
        Dd[1] = Dvalue(atmD,4);
        intn[0] = d_eri2Center_dR2(0,0,2,0,0,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR2(0,0,2,1,0,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[2] = d_eri2Center_dR2(0,0,2,2,0,0,2,2,RCD,Dd,atmC,atmD,type);    //(sg,dl|sg,dl)=(sg,-dl|sg,-dl)
        //(sd|sd) rotate
        for (size_t idsdb = 1; idsdb < 6; ++idsdb) {
          for (size_t idsdk = 1; idsdk < 6; ++idsdk) {
            block(10 + idsdb,10 + idsdk) = intn[0]*SDrot(idsdb,1)*SDrot(idsdk,1) + intn[1]*(SDrot(idsdb,2)*SDrot(idsdk,2) + SDrot(idsdb,3)*SDrot(idsdk,3)) + intn[2]*(SDrot(idsdb,4)*SDrot(idsdk,4) + SDrot(idsdb,5)*SDrot(idsdk,5));
          }
        }
        //(sd|pd) integrals
        Dd[1] = Dvalue(atmD,5);
        intn[0] = d_eri2Center_dR2(0,0,2,0,1,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR2(0,0,2,0,1,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center_dR2(0,0,2,1,1,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[3] = d_eri2Center_dR2(0,0,2,1,1,1,2,2,RCD,Dd,atmC,atmD,type);    //(sg,pi|pi,dl)=-(sg,-pi|-pi,dl)=(sg,-pi|pi,-dl)=(sg,pi|-pi,-dl)
        intn[4] = d_eri2Center_dR2(0,0,2,1,1,1,2,0,RCD,Dd,atmC,atmD,type);    //(sg,pi|pi,sg)=(sg,-pi|-pi,sg)
        //(sd|pd) rotate
        for (size_t idsd = 1; idsd < 6; ++idsd) {
          for (size_t idpd1 = 1; idpd1 < 6; ++idpd1) {
            for (size_t idpd2 = 1; idpd2 < 4; ++idpd2) {
              block(10 + idsd,5*(idpd2+2) + idpd1) = -intn[0]*PDrot((idpd1 - 1)*3 + idpd2,1)*SDrot(idsd,1) - intn[1]*(PDrot((idpd1 - 1)*3 + idpd2,2) + PDrot((idpd1 - 1)*3 + idpd2,3))*SDrot(idsd,1);
              block(10 + idsd,5*(idpd2+2) + idpd1) -= intn[2]*(PDrot((idpd1 - 1)*3 + idpd2,4)*SDrot(idsd,2) + PDrot((idpd1 - 1)*3 + idpd2,5)*SDrot(idsd,3));
              block(10 + idsd,5*(idpd2+2) + idpd1) -= intn[3]*((PDrot((idpd1 - 1)*3 + idpd2,10) + PDrot((idpd1 - 1)*3 + idpd2,15))*SDrot(idsd,2) + (PDrot((idpd1 - 1)*3 + idpd2,11) - PDrot((idpd1 - 1)*3 + idpd2,14))*SDrot(idsd,3));
              block(10 + idsd,5*(idpd2+2) + idpd1) -= intn[4]*(PDrot((idpd1 - 1)*3 + idpd2,8)*SDrot(idsd,2) + PDrot((idpd1 - 1)*3 + idpd2,12)*SDrot(idsd,3));
            }
          }
        }
        //(pd|sd) integrals
        Dd[0] = Dvalue(atmC,5);
        Dd[1] = Dvalue(atmD,4);
        intn[0] = d_eri2Center_dR2(1,0,2,0,0,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR2(1,1,2,1,0,0,2,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[2] = d_eri2Center_dR2(1,0,2,1,0,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[3] = d_eri2Center_dR2(1,1,2,2,0,0,2,1,RCD,Dd,atmC,atmD,type);    //(pi,dl|sg,pi)=-(-pi,dl|sg,-pi)=(pi,-dl|sg,-pi)=(-pi,-dl|sg,pi)
        intn[4] = d_eri2Center_dR2(1,1,2,0,0,0,2,1,RCD,Dd,atmC,atmD,type);    //(pi,sg|sg,pi)=(-pi,sg|sg,-pi)
        //(pd|sd) rotate
        for (size_t idpd1 = 1; idpd1 < 6; ++idpd1) {
          for (size_t idpd2 = 1; idpd2 < 4; ++idpd2) {
            for (size_t idsd = 1; idsd < 6; ++idsd) {
              block(5*(idpd2+2) + idpd1,10 + idsd) = -intn[0]*PDrot((idpd1 - 1)*3 + idpd2,1)*SDrot(idsd,1) - intn[1]*(PDrot((idpd1 - 1)*3 + idpd2,2) + PDrot((idpd1 - 1)*3 + idpd2,3))*SDrot(idsd,1);
              block(5*(idpd2+2) + idpd1,10 + idsd) -= intn[2]*(PDrot((idpd1 - 1)*3 + idpd2,4)*SDrot(idsd,2) + PDrot((idpd1 - 1)*3 + idpd2,5)*SDrot(idsd,3));
              block(5*(idpd2+2) + idpd1,10 + idsd) -= intn[3]*((PDrot((idpd1 - 1)*3 + idpd2,10) + PDrot((idpd1 - 1)*3 + idpd2,15))*SDrot(idsd,2) + (PDrot((idpd1 - 1)*3 + idpd2,11) - PDrot((idpd1 - 1)*3 + idpd2,14))*SDrot(idsd,3));
              block(5*(idpd2+2) + idpd1,10 + idsd) -= intn[4]*(PDrot((idpd1 - 1)*3 + idpd2,8)*SDrot(idsd,2) + PDrot((idpd1 - 1)*3 + idpd2,12)*SDrot(idsd,3));
            }
          }
        }
        //(dd|sd) integrals
        Dd[0] = Dvalue(atmC,6);
        Dd[1] = Dvalue(atmD,4);
        intn[0] = d_eri2Center_dR2(2,0,2,0,0,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR2(2,1,2,1,0,0,2,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[2] = d_eri2Center_dR2(2,2,2,2,0,0,2,0,RCD,Dd,atmC,atmD,type);    //(dl,dl|sg,sg)=(-dl,-dl|sg,sg)
        intn[3] = d_eri2Center_dR2(2,0,2,1,0,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[4] = d_eri2Center_dR2(2,0,2,2,0,0,2,2,RCD,Dd,atmC,atmD,type);    //(sg,dl|sg,dl)=(sg,-dl|sg,-dl)
        intn[5] = d_eri2Center_dR2(2,1,2,1,0,0,2,2,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,dl)=-(-pi,-pi|sg,dl)=(pi,-pi|sg,-dl)
        intn[6] = d_eri2Center_dR2(2,1,2,2,0,0,2,1,RCD,Dd,atmC,atmD,type);    //(pi,dl|sg,pi)=(-pi,-dl|sg,pi)=(pi,-dl|sg,-pi)=-(-pi,dl|sg,-pi)
        //(dd|sd) rotate
        for (size_t iddd = 1; iddd < 16; ++iddd) {
          for (size_t idsd = 1; idsd < 6; ++idsd) {
            block(30 + iddd,10 + idsd) = intn[0]*DDrot(iddd,1)*SDrot(idsd,1) + intn[1]*(DDrot(iddd,2) + DDrot(iddd,3))*SDrot(idsd,1) + intn[2]*(DDrot(iddd,4) + DDrot(iddd,5))*SDrot(idsd,1);
            block(30 + iddd,10 + idsd) += intn[3]*(DDrot(iddd,6)*SDrot(idsd,2) + DDrot(iddd,7)*SDrot(idsd,3)) + intn[4]*(DDrot(iddd,8)*SDrot(idsd,4) + DDrot(iddd,9)*SDrot(idsd,5));
            block(30 + iddd,10 + idsd) += intn[5]*((DDrot(iddd,2) - DDrot(iddd,3))*SDrot(idsd,4) + DDrot(iddd,10)*SDrot(idsd,5));
            block(30 + iddd,10 + idsd) += intn[6]*((DDrot(iddd,11) + DDrot(iddd,14))*SDrot(idsd,2) + (DDrot(iddd,12) - DDrot(iddd,13))*SDrot(idsd,3));
          }
        }
        //(sd|dd) integrals
        Dd[0] = Dvalue(atmC,4);
        Dd[1] = Dvalue(atmD,6);
        intn[0] = d_eri2Center_dR2(0,0,2,0,2,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR2(0,0,2,0,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center_dR2(0,0,2,0,2,2,2,2,RCD,Dd,atmC,atmD,type);    //(sg,sg|dl,dl)=(sg,sg|-dl,-dl)
        intn[3] = d_eri2Center_dR2(0,0,2,1,2,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[4] = d_eri2Center_dR2(0,0,2,2,2,0,2,2,RCD,Dd,atmC,atmD,type);    //(sg,dl|sg,dl)=(sg,-dl|sg,-dl)
        intn[5] = d_eri2Center_dR2(0,0,2,2,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,dl|pi,pi)=-(sg,dl|-pi,-pi)=(sg,-dl|pi,-pi)
        intn[6] = d_eri2Center_dR2(0,0,2,1,2,1,2,2,RCD,Dd,atmC,atmD,type);    //(sg,pi|pi,dl)=(sg,pi|-pi,-dl)=(sg,-pi|pi,-dl)=-(sg,-pi|-pi,dl)
        //(sd|dd) rotate
        for (size_t idsd = 1; idsd < 6; ++idsd) {
          for (size_t iddd = 1; iddd < 16; ++iddd) {
            block(10 + idsd,30 + iddd) = intn[0]*DDrot(iddd,1)*SDrot(idsd,1) + intn[1]*(DDrot(iddd,2) + DDrot(iddd,3))*SDrot(idsd,1) + intn[2]*(DDrot(iddd,4) + DDrot(iddd,5))*SDrot(idsd,1);
            block(10 + idsd,30 + iddd) += intn[3]*(DDrot(iddd,6)*SDrot(idsd,2) + DDrot(iddd,7)*SDrot(idsd,3)) + intn[4]*(DDrot(iddd,8)*SDrot(idsd,4) + DDrot(iddd,9)*SDrot(idsd,5));
            block(10 + idsd,30 + iddd) += intn[5]*((DDrot(iddd,2) - DDrot(iddd,3))*SDrot(idsd,4) + DDrot(iddd,10)*SDrot(idsd,5));
            block(10 + idsd,30 + iddd) += intn[6]*((DDrot(iddd,11) + DDrot(iddd,14))*SDrot(idsd,2) + (DDrot(iddd,12) - DDrot(iddd,13))*SDrot(idsd,3));
          }
        }
        //(pd|pd) integrals
        Dd[0] = Dvalue(atmC,5);
        Dd[1] = Dvalue(atmD,5);
        intn[0] = d_eri2Center_dR2(1,0,2,0,1,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR2(1,1,2,1,1,0,2,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)=(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center_dR2(1,1,2,1,1,1,2,1,RCD,Dd,atmC,atmD,type);    //(pi,pi|pi,pi)=(-pi,-pi|pi,pi)=(pi,pi|-pi,-pi)=(-pi,-pi|-pi,-pi)
        intn[3] = d_eri2Center_dR2(1,0,2,1,1,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[4] = d_eri2Center_dR2(1,1,2,0,1,1,2,0,RCD,Dd,atmC,atmD,type);    //(pi,sg|pi,sg)=(-pi,sg|-pi,sg)
        intn[5] = d_eri2Center_dR2(1,1,2,2,1,1,2,2,RCD,Dd,atmC,atmD,type);    //(pi,dl|pi,dl)=(pi,dl|-pi,-dl)=(-pi,dl|-pi,dl)=-(-pi,dl|pi,-dl)=-(pi,-dl|-pi,dl)=(pi,-dl|pi,-dl)=(-pi,-dl|pi,dl)=(-pi,-dl|-pi,-dl)
        intn[6] = d_eri2Center_dR2(1,1,2,0,1,0,2,1,RCD,Dd,atmC,atmD,type);    //(pi,sg|sg,pi)=(-pi,sg|sg,-pi)=(sg,pi|pi,sg)=(sg,-pi|-pi,sg)
        intn[7] = d_eri2Center_dR2(1,1,2,2,1,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|pi,dl)=(sg,pi|-pi,-dl)=-(sg,-pi|-pi,dl)=(sg,-pi|pi,-dl)=(pi,dl|sg,pi)=-(-pi,dl|sg,-pi)=(pi,-dl|sg,-pi)=(-pi,-dl|sg,pi)
        intn[8] = d_eri2Center_dR2(1,1,2,2,1,1,2,0,RCD,Dd,atmC,atmD,type);    //(pi,sg|pi,dl)=(pi,sg|-pi,-dl)=-(-pi,sg|-pi,dl)=(-pi,sg|pi,-dl)=(pi,dl|pi,sg)=-(-pi,dl|-pi,sg)=(pi,-dl|-pi,sg)=(-pi,-dl|pi,sg)
        //(pd|pd) rotate
        for (size_t idpdb1 = 1; idpdb1 < 6; ++idpdb1) {            //d orbital on bra
          for (size_t idpdb2 = 1; idpdb2 < 4; ++idpdb2) {          //p orbital on bra
            for (size_t idpdk1 = 1; idpdk1 < 6; ++idpdk1) {        //d orbital on ket
              for (size_t idpdk2 = 1; idpdk2 < 4; ++idpdk2) {      //p orbital on ket
                cnt1 = (idpdb1 - 1)*3 + idpdb2;
                cnt2 = (idpdk1 - 1)*3 + idpdk2;
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) = intn[0]*PDrot(cnt1,1)*PDrot(cnt2,1) + intn[1]*((PDrot(cnt1,2) + PDrot(cnt1,3))*PDrot(cnt2,1) + PDrot(cnt1,1)*(PDrot(cnt2,2) + PDrot(cnt2,3)));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[2]*(PDrot(cnt1,2) + PDrot(cnt1,3))*(PDrot(cnt2,2) + PDrot(cnt2,3)) + intn[3]*(PDrot(cnt1,4)*PDrot(cnt2,4) + PDrot(cnt1,5)*PDrot(cnt2,5));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[4]*(PDrot(cnt1,8)*PDrot(cnt2,8) + PDrot(cnt1,12)*PDrot(cnt2,12));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[5]*((PDrot(cnt1,10) + PDrot(cnt1,15))*(PDrot(cnt2,10) + PDrot(cnt2,15)) + (PDrot(cnt1,14) - PDrot(cnt1,11))*(PDrot(cnt2,14) - PDrot(cnt2,11)));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[6]*(PDrot(cnt1,4)*PDrot(cnt2,8) + PDrot(cnt1,8)*PDrot(cnt2,4) + PDrot(cnt1,5)*PDrot(cnt2,12) + PDrot(cnt1,12)*PDrot(cnt2,5));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[7]*((PDrot(cnt1,10) + PDrot(cnt1,15))*PDrot(cnt2,4) + PDrot(cnt1,4)*(PDrot(cnt2,10) + PDrot(cnt2,15)) + (PDrot(cnt1,11) - PDrot(cnt1,14))*PDrot(cnt2,5) + PDrot(cnt1,5)*(PDrot(cnt2,11) - PDrot(cnt2,14)));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[8]*((PDrot(cnt1,10) + PDrot(cnt1,15))*PDrot(cnt2,8) + PDrot(cnt1,8)*(PDrot(cnt2,10) + PDrot(cnt2,15)) + (PDrot(cnt1,11) - PDrot(cnt1,14))*PDrot(cnt2,12) + PDrot(cnt1,12)*(PDrot(cnt2,11) - PDrot(cnt2,14)));
              }
            }
          }
        }
        //(pd|dd) integrals
        Dd[0] = Dvalue(atmC,5);
        Dd[1] = Dvalue(atmD,6);
        intn[0] = d_eri2Center_dR2(1,0,2,0,2,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR2(1,0,2,0,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center_dR2(1,1,2,1,2,0,2,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[3] = d_eri2Center_dR2(1,1,2,1,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(pi,pi|pi,pi)=(-pi,-pi|-pi,-pi)=(pi,pi|-pi,-pi)=(-pi,-pi|pi,pi)
        intn[4] = d_eri2Center_dR2(1,0,2,0,2,2,2,2,RCD,Dd,atmC,atmD,type);    //(sg,sg|dl,dl)=(sg,sg|-dl,-dl)
        intn[5] = d_eri2Center_dR2(1,1,2,1,2,2,2,2,RCD,Dd,atmC,atmD,type);    //(pi,pi|dl,dl)=(-pi,-pi|-dl,-dl)=(pi,pi|-dl,-dl)=(-pi,-pi|dl,dl)
        intn[6] = d_eri2Center_dR2(1,0,2,1,2,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[7] = d_eri2Center_dR2(1,1,2,0,2,1,2,0,RCD,Dd,atmC,atmD,type);    //(pi,sg|pi,sg)=(-pi,sg|-pi,sg)
        intn[8] = d_eri2Center_dR2(1,1,2,2,2,1,2,2,RCD,Dd,atmC,atmD,type);    //(pi,dl|pi,dl)=(pi,dl|-pi,-dl)=(-pi,-dl|-pi,-dl)=(-pi,-dl|pi,dl)=(-pi,dl|-pi,dl)=-(-pi,dl|pi,-dl)=(pi,-dl|pi,-dl)=-(pi,-dl|-pi,dl)
        intn[9] = d_eri2Center_dR2(1,1,2,0,2,1,2,2,RCD,Dd,atmC,atmD,type);    //(pi,sg|pi,dl)=(pi,sg|-pi,-dl)=-(-pi,sg|-pi,dl)=(-pi,sg|pi,-dl)
        intn[10] = d_eri2Center_dR2(1,0,2,1,2,1,2,2,RCD,Dd,atmC,atmD,type);   //(sg,pi|pi,dl)=(sg,pi|-pi,-dl)=-(sg,-pi|-pi,dl)=(sg,-pi|pi,-dl)
        intn[11] = d_eri2Center_dR2(1,1,2,2,2,0,2,1,RCD,Dd,atmC,atmD,type);   //(pi,dl|sg,pi)=(-pi,-dl|sg,pi)=-(-pi,dl|sg,-pi)=(pi,-dl|sg,-pi)
        //(pd|dd) rotate
        for (size_t iddd = 1; iddd < 16; ++iddd) {
          for (size_t idpd1 = 1; idpd1 < 6; ++idpd1) {
            for (size_t idpd2 = 1; idpd2 < 4; ++idpd2) {
              cnt1 = (idpd1 - 1)*3 + idpd2;
              block(5*(idpd2+2) + idpd1,30 + iddd) = -intn[0]*PDrot(cnt1,1)*DDrot(iddd,1) - intn[1]*PDrot(cnt1,1)*(DDrot(iddd,2) + DDrot(iddd,3));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[2]*(PDrot(cnt1,2) + PDrot(cnt1,3))*DDrot(iddd,1) + intn[3]*(PDrot(cnt1,2) + PDrot(cnt1,3))*(DDrot(iddd,2) + DDrot(iddd,3));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[4]*PDrot(cnt1,1)*(DDrot(iddd,4) + DDrot(iddd,5)) + intn[5]*(PDrot(cnt1,2) + PDrot(cnt1,3))*(DDrot(iddd,4) + DDrot(iddd,5));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[6]*(PDrot(cnt1,4)*DDrot(iddd,6) + PDrot(cnt1,5)*DDrot(iddd,7)) + intn[7]*(PDrot(cnt1,8)*DDrot(iddd,6) + PDrot(cnt1,12)*DDrot(iddd,7));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[8]*((PDrot(cnt1,10) + PDrot(cnt1,15))*(DDrot(iddd,11) + DDrot(iddd,14)) + (PDrot(cnt1,14) - PDrot(cnt1,11))*(DDrot(iddd,13) - DDrot(iddd,12)));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[9]*(PDrot(cnt1,8)*(DDrot(iddd,11) + DDrot(iddd,14)) + PDrot(cnt1,12)*(DDrot(iddd,12) - DDrot(iddd,13)));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[10]*(PDrot(cnt1,4)*(DDrot(iddd,11) + DDrot(iddd,14)) + PDrot(cnt1,5)*(DDrot(iddd,12) - DDrot(iddd,13)));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[11]*((PDrot(cnt1,10) + PDrot(cnt1,15))*DDrot(iddd,6) + (PDrot(cnt1,11) - PDrot(cnt1,14))*DDrot(iddd,7));
            }
          }
        }
        //(dd|pd) integrals
        Dd[0] = Dvalue(atmC,6);
        Dd[1] = Dvalue(atmD,5);
        intn[0] = d_eri2Center_dR2(2,0,2,0,1,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR2(2,1,2,1,1,0,2,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[2] = d_eri2Center_dR2(2,0,2,0,1,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[3] = d_eri2Center_dR2(2,1,2,1,1,1,2,1,RCD,Dd,atmC,atmD,type);    //(pi,pi|pi,pi)=(-pi,-pi|-pi,-pi)=(pi,pi|-pi,-pi)=(-pi,-pi|pi,pi)
        intn[4] = d_eri2Center_dR2(2,2,2,2,1,0,2,0,RCD,Dd,atmC,atmD,type);    //(dl,dl|sg,sg)=(-dl,-dl|sg,sg)
        intn[5] = d_eri2Center_dR2(2,2,2,2,1,1,2,1,RCD,Dd,atmC,atmD,type);    //(dl,dl|pi,pi)=(-dl,-dl|-pi,-pi)=(-dl,-dl|pi,pi)=(dl,dl|-pi,-pi)
        intn[6] = d_eri2Center_dR2(2,0,2,1,1,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[7] = d_eri2Center_dR2(2,1,2,0,1,1,2,0,RCD,Dd,atmC,atmD,type);    //(pi,sg|pi,sg)=(-pi,sg|-pi,sg)
        intn[8] = d_eri2Center_dR2(2,1,2,2,1,1,2,2,RCD,Dd,atmC,atmD,type);    //(pi,dl|pi,dl)=(-pi,-dl|-pi,-dl)=(-pi,-dl|pi,dl)=(pi,dl|-pi,-dl)=(-pi,dl|-pi,dl)=(pi,-dl|pi,-dl)=-(-pi,dl|pi,-dl)=-(pi,-dl|-pi,dl)
        intn[9] = d_eri2Center_dR2(2,1,2,2,1,1,2,0,RCD,Dd,atmC,atmD,type);    //(pi,dl|pi,sg)=-(-pi,dl|-pi,sg)=(pi,-dl|-pi,sg)=(-pi,-dl|pi,sg)
        intn[10] = d_eri2Center_dR2(2,1,2,2,1,0,2,1,RCD,Dd,atmC,atmD,type);   //(pi,dl|sg,pi)=-(-pi,dl|sg,-pi)=(pi,-dl|sg,-pi)=(-pi,-dl|sg,pi)
        intn[11] = d_eri2Center_dR2(2,0,2,1,1,1,2,2,RCD,Dd,atmC,atmD,type);   //(sg,pi|pi,dl)=(sg,pi|-pi,-dl)=-(sg,-pi|-pi,dl)=(sg,-pi|pi,-dl)
        //(dd|pd) rotate
        for (size_t iddd = 1; iddd < 16; ++iddd) {
          for (size_t idpd1 = 1; idpd1 < 6; ++idpd1) {
            for (size_t idpd2 = 1; idpd2 < 4; ++idpd2) {
              cnt1 = (idpd1 - 1)*3 + idpd2;
              block(30 + iddd,5*(idpd2+2) + idpd1) = -intn[0]*PDrot(cnt1,1)*DDrot(iddd,1) - intn[1]*PDrot(cnt1,1)*(DDrot(iddd,2) + DDrot(iddd,3));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[2]*(PDrot(cnt1,2) + PDrot(cnt1,3))*DDrot(iddd,1) + intn[3]*(PDrot(cnt1,2) + PDrot(cnt1,3))*(DDrot(iddd,2) + DDrot(iddd,3));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[4]*PDrot(cnt1,1)*(DDrot(iddd,4) + DDrot(iddd,5)) + intn[5]*(PDrot(cnt1,2) + PDrot(cnt1,3))*(DDrot(iddd,4) + DDrot(iddd,5));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[6]*(PDrot(cnt1,4)*DDrot(iddd,6) + PDrot(cnt1,5)*DDrot(iddd,7)) + intn[7]*(PDrot(cnt1,8)*DDrot(iddd,6) + PDrot(cnt1,12)*DDrot(iddd,7));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[8]*((PDrot(cnt1,10) + PDrot(cnt1,15))*(DDrot(iddd,11) + DDrot(iddd,14)) + (PDrot(cnt1,14) - PDrot(cnt1,11))*(DDrot(iddd,13) - DDrot(iddd,12)));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[9]*(PDrot(cnt1,8)*(DDrot(iddd,11) + DDrot(iddd,14)) + PDrot(cnt1,12)*(DDrot(iddd,12) - DDrot(iddd,13)));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[10]*(PDrot(cnt1,4)*(DDrot(iddd,11) + DDrot(iddd,14)) + PDrot(cnt1,5)*(DDrot(iddd,12) - DDrot(iddd,13)));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[11]*((PDrot(cnt1,10) + PDrot(cnt1,15))*DDrot(iddd,6) + (PDrot(cnt1,11) - PDrot(cnt1,14))*DDrot(iddd,7));
            }
          }
        }
        //(dd|dd) integrals
        Dd[0] = Dvalue(atmC,6);
        Dd[1] = Dvalue(atmD,6);
        intn[0] = d_eri2Center_dR2(2,0,2,0,2,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR2(2,0,2,0,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center_dR2(2,0,2,0,2,2,2,2,RCD,Dd,atmC,atmD,type);    //(sg,sg|dl,dl)=(sg,sg|-dl,-dl)
        intn[3] = d_eri2Center_dR2(2,1,2,1,2,0,2,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[4] = d_eri2Center_dR2(2,2,2,2,2,0,2,0,RCD,Dd,atmC,atmD,type);    //(dl,dl|sg,sg)=(-dl,-dl|sg,sg)
        intn[5] = d_eri2Center_dR2(2,2,2,2,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(dl,dl|pi,pi)=(-dl,-dl|pi,pi)=(dl,dl|-pi,-pi)=(-dl,-dl|-pi,-pi)
        intn[6] = d_eri2Center_dR2(2,1,2,1,2,2,2,2,RCD,Dd,atmC,atmD,type);    //(pi,pi|dl,dl)=(-pi,-pi|dl,dl)=(pi,pi|-dl,-dl)=(-pi,-pi|-dl,-dl)
        intn[7] = d_eri2Center_dR2(2,1,2,1,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(pi,pi|pi,pi)=(-pi,-pi|-pi,-pi)
        intn[8] = d_eri2Center_dR2(2,2,2,2,2,2,2,2,RCD,Dd,atmC,atmD,type);    //(dl,dl|dl,dl)=(-dl,-dl|-dl,-dl)
        intn[9] = d_eri2Center_dR2(2,1,2,1,2,-1,2,-1,RCD,Dd,atmC,atmD,type);  //(pi,pi|-pi,-pi)=(-pi,-pi|pi,pi)
        intn[10] = d_eri2Center_dR2(2,-2,2,-2,2,2,2,2,RCD,Dd,atmC,atmD,type); //(-dl,-dl|dl,dl)=(dl,dl|-dl,-dl)
        intn[11] = d_eri2Center_dR2(2,0,2,1,2,0,2,1,RCD,Dd,atmC,atmD,type);   //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[12] = d_eri2Center_dR2(2,0,2,2,2,0,2,2,RCD,Dd,atmC,atmD,type);   //(sg,dl|sg,dl)=(sg,-dl|sg,-dl)
        intn[13] = d_eri2Center_dR2(2,1,2,2,2,1,2,2,RCD,Dd,atmC,atmD,type);   //(pi,dl|pi,dl)=(pi,dl|-pi,-dl)=(-pi,-dl|pi,dl)=(-pi,-dl|-pi,-dl)=-(-pi,dl|pi,-dl)=-(pi,-dl|-pi,dl)=(-pi,dl|-pi,dl)=(pi,-dl|pi,-dl)
        intn[14] = d_eri2Center_dR2(2,1,2,-1,2,1,2,-1,RCD,Dd,atmC,atmD,type); //(pi,-pi|pi,-pi)
        intn[15] = d_eri2Center_dR2(2,1,2,-2,2,0,2,-1,RCD,Dd,atmC,atmD,type); //(pi,dl|sg,pi)=(-pi,-dl|sg,pi)=(pi,-dl|sg,-pi)=-(-pi,dl|sg,-pi)=(sg,pi|pi,dl)=(sg,pi|-pi,-dl)=(sg,-pi|pi,-dl)=-(sg,-pi|-pi,dl)
        intn[16] = d_eri2Center_dR2(2,1,2,1,2,0,2,2,RCD,Dd,atmC,atmD,type);   //(pi,pi|sg,dl)=-(-pi,-pi|sg,dl)=(pi,-pi|sg,-dl)=(sg,dl|pi,pi)=-(sg,dl|-pi,-pi)=(sg,-dl|pi,-pi)
        //(dd|dd) rotate
        for (size_t iddd1 = 1; iddd1 < 16; ++iddd1) {
          for (size_t iddd2 = 1; iddd2 < 16; ++iddd2) {
            block(30 + iddd1,30 + iddd2) = intn[0]*DDrot(iddd1,1)*DDrot(iddd2,1) + intn[1]*DDrot(iddd1,1)*(DDrot(iddd2,2) + DDrot(iddd2,3)) + intn[2]*DDrot(iddd1,1)*(DDrot(iddd2,4) + DDrot(iddd2,5));
            block(30 + iddd1,30 + iddd2) += intn[3]*(DDrot(iddd1,2) + DDrot(iddd1,3))*DDrot(iddd2,1) + intn[4]*(DDrot(iddd1,4) + DDrot(iddd1,5))*DDrot(iddd2,1);
            block(30 + iddd1,30 + iddd2) += intn[5]*(DDrot(iddd1,4) + DDrot(iddd1,5))*(DDrot(iddd2,2) + DDrot(iddd2,3)) + intn[6]*(DDrot(iddd1,2) + DDrot(iddd1,3))*(DDrot(iddd2,4) + DDrot(iddd2,5));
            block(30 + iddd1,30 + iddd2) += intn[7]*(DDrot(iddd1,2)*DDrot(iddd2,2) + DDrot(iddd1,3)*DDrot(iddd2,3)) + intn[8]*(DDrot(iddd1,4)*DDrot(iddd2,4) + DDrot(iddd1,5)*DDrot(iddd2,5));
            block(30 + iddd1,30 + iddd2) += intn[9]*(DDrot(iddd1,2)*DDrot(iddd2,3) + DDrot(iddd1,3)*DDrot(iddd2,2)) + intn[10]*(DDrot(iddd1,4)*DDrot(iddd2,5) + DDrot(iddd1,5)*DDrot(iddd2,4));
            block(30 + iddd1,30 + iddd2) += intn[11]*(DDrot(iddd1,6)*DDrot(iddd2,6) + DDrot(iddd1,7)*DDrot(iddd2,7)) + intn[12]*(DDrot(iddd1,8)*DDrot(iddd2,8) + DDrot(iddd1,9)*DDrot(iddd2,9));
            block(30 + iddd1,30 + iddd2) += intn[13]*((DDrot(iddd1,11) + DDrot(iddd1,14))*(DDrot(iddd2,11) + DDrot(iddd2,14)) + (DDrot(iddd1,12) - DDrot(iddd1,13))*(DDrot(iddd2,12) - DDrot(iddd2,13)));
            block(30 + iddd1,30 + iddd2) += intn[14]*DDrot(iddd1,10)*DDrot(iddd2,10);
            block(30 + iddd1,30 + iddd2) += intn[15]*((DDrot(iddd1,11) + DDrot(iddd1,14))*DDrot(iddd2,6) + (DDrot(iddd1,12) - DDrot(iddd1,13))*DDrot(iddd2,7) + DDrot(iddd1,6)*(DDrot(iddd2,11) + DDrot(iddd2,14)) + DDrot(iddd1,7)*(DDrot(iddd2,12) - DDrot(iddd2,13)));
            block(30 + iddd1,30 + iddd2) += intn[16]*((DDrot(iddd1,2) - DDrot(iddd1,3))*DDrot(iddd2,8) + DDrot(iddd1,10)*DDrot(iddd2,9) + DDrot(iddd1,8)*(DDrot(iddd2,2) - DDrot(iddd2,3)) + DDrot(iddd1,9)*DDrot(iddd2,10));
          }
        }
      }
    }
  }
  void IntegralBlock2C_dAB(int type, matrixE & block, int derA, int derB, int atmC, int atmD, double RCD, double cost, double sint, double cosp, double sinp) {
    //function calculating a block of second-derivatives with respect to pair of angular variables of two-center integrals
    //2 -> theta; 3 -> phi
    matrixE SProt_dA(1,1);
    matrixE PProt_dA(1,1);
    matrixE SDrot_dA(1,1);
    matrixE PDrot_dA(1,1);
    matrixE DDrot_dA(1,1);
    matrixE SProt_dB(1,1);
    matrixE PProt_dB(1,1);
    matrixE SDrot_dB(1,1);
    matrixE PDrot_dB(1,1);
    matrixE DDrot_dB(1,1);
    matrixE SProt_dAB(1,1);
    matrixE PProt_dAB(1,1);
    matrixE SDrot_dAB(1,1);
    matrixE PDrot_dAB(1,1);
    matrixE DDrot_dAB(1,1);
    int cnt1;
    int cnt2;
    size_t nrows = block.rows();
    size_t ncols = block.cols();
    if ((nrows > 1)||(ncols > 1)) {
      SProt = SPtransf(cost,sint,cosp,sinp);
      PProt = PPtransf(cost,sint,cosp,sinp);
      if ((nrows > 10)||(ncols > 10)) {
        SDrot = SDtransf(cost,sint,cosp,sinp);
        PDrot = PDtransf(cost,sint,cosp,sinp);
        DDrot = DDtransf(cost,sint,cosp,sinp);
      }
      if (derA == 2) {
        SProt_dA = SPtransf_dt(cost,sint,cosp,sinp);
        PProt_dA = PPtransf_dt(cost,sint,cosp,sinp);
        if ((nrows > 10)||(ncols > 10)) {
          SDrot_dA = SDtransf_dt(cost,sint,cosp,sinp);
          PDrot_dA = PDtransf_dt(cost,sint,cosp,sinp);
          DDrot_dA = DDtransf_dt(cost,sint,cosp,sinp);
        }
      }
      else if (derA == 3) {
        SProt_dA = SPtransf_dp(cost,sint,cosp,sinp);
        PProt_dA = PPtransf_dp(cost,sint,cosp,sinp);
        if ((nrows > 10)||(ncols > 10)) {
          SDrot_dA = SDtransf_dp(cost,sint,cosp,sinp);
          PDrot_dA = PDtransf_dp(cost,sint,cosp,sinp);
          DDrot_dA = DDtransf_dp(cost,sint,cosp,sinp);
        }
      }
      if (derB == 2) {
        SProt_dB = SPtransf_dt(cost,sint,cosp,sinp);
        PProt_dB = PPtransf_dt(cost,sint,cosp,sinp);
        if ((nrows > 10)||(ncols > 10)) {
          SDrot_dB = SDtransf_dt(cost,sint,cosp,sinp);
          PDrot_dB = PDtransf_dt(cost,sint,cosp,sinp);
          DDrot_dB = DDtransf_dt(cost,sint,cosp,sinp);
        }
      }
      else if (derB == 3) {
        SProt_dB = SPtransf_dp(cost,sint,cosp,sinp);
        PProt_dB = PPtransf_dp(cost,sint,cosp,sinp);
        if ((nrows > 10)||(ncols > 10)) {
          SDrot_dB = SDtransf_dp(cost,sint,cosp,sinp);
          PDrot_dB = PDtransf_dp(cost,sint,cosp,sinp);
          DDrot_dB = DDtransf_dp(cost,sint,cosp,sinp);
        }
      }
      if (derA == derB) {
        if (derA == 2) {
          SProt_dAB = SPtransf_dt2(cost,sint,cosp,sinp);
          PProt_dAB = PPtransf_dt2(cost,sint,cosp,sinp);
          if ((nrows > 10)||(ncols > 10)) {
            SDrot_dAB = SDtransf_dt2(cost,sint,cosp,sinp);
            PDrot_dAB = PDtransf_dt2(cost,sint,cosp,sinp);
            DDrot_dAB = DDtransf_dt2(cost,sint,cosp,sinp);
          }
        }
        else if (derA == 3) {
          SProt_dAB = SPtransf_dp2(cost,sint,cosp,sinp);
          PProt_dAB = PPtransf_dp2(cost,sint,cosp,sinp);
          if ((nrows > 10)||(ncols > 10)) {
            SDrot_dAB = SDtransf_dp2(cost,sint,cosp,sinp);
            PDrot_dAB = PDtransf_dp2(cost,sint,cosp,sinp);
            DDrot_dAB = DDtransf_dp2(cost,sint,cosp,sinp);
          }
        }
      }
      else {
        SProt_dAB = SPtransf_dtdp(cost,sint,cosp,sinp);
        PProt_dAB = PPtransf_dtdp(cost,sint,cosp,sinp);
        if ((nrows > 10)||(ncols > 10)) {
          SDrot_dAB = SDtransf_dtdp(cost,sint,cosp,sinp);
          PDrot_dAB = PDtransf_dtdp(cost,sint,cosp,sinp);
          DDrot_dAB = DDtransf_dtdp(cost,sint,cosp,sinp);
        }
      }
    }
    D[0] = Dvalue(atmC,1);
    D[1] = Dvalue(atmC,2);
    D[2] = Dvalue(atmD,1);
    D[3] = Dvalue(atmD,2);
    //integral calculation
    block(1,1) = 0.0;                                                                                                                      //(ss|ss)
    if (nrows > 1) {
      intn[0] = eri2Center(0,0,1,0,0,0,0,0,RCD,D,atmC,atmD,type);
      intn[1] = eri2Center(1,0,1,0,0,0,0,0,RCD,D,atmC,atmD,type);
      intn[2] = eri2Center(1,1,1,1,0,0,0,0,RCD,D,atmC,atmD,type);
      for (size_t idx = 0; idx < 6; ++idx) {
        if (idx < 3) {block(idx + 2,1) = -intn[0]*SProt_dAB(idx + 1,1);}                                                           //(sp|ss) integrals
        block(5 + idx,1) = intn[1]*PProt_dAB(idx + 1,1) + intn[2]*(PProt_dAB(idx + 1,2) + PProt_dAB(idx + 1,3));   //(pp|ss) integrals
      }
      if (nrows > 10) {
        //(sd|ss) integrals
        Dd[0] = Dvalue(atmC,4);
        Dd[1] = Dvalue(atmD,1);
        intn[0] = d_eri2Center(0,0,2,0,0,0,0,0,RCD,Dd,atmC,atmD,type);
        for (size_t idsd = 1; idsd < 6; ++idsd) {
          block(10 + idsd,1) = intn[0]*SDrot_dAB(idsd,1);
        }
        //(pd|ss) integrals
        Dd[0] = Dvalue(atmC,5);
        intn[0] = d_eri2Center(1,0,2,0,0,0,0,0,RCD,Dd,atmC,atmD,type);
        intn[1] = d_eri2Center(1,1,2,1,0,0,0,0,RCD,Dd,atmC,atmD,type);
        for (int idpd1 = 1; idpd1 < 6; ++idpd1) {
          for (int idpd2 = 1; idpd2 < 4; ++idpd2) {
            block(5*(idpd2+2) + idpd1,1) = -(intn[0]*PDrot_dAB((idpd1 - 1)*3 + idpd2,1) + intn[1]*(PDrot_dAB((idpd1 - 1)*3 + idpd2,2) + PDrot_dAB((idpd1 - 1)*3 + idpd2,3)));
          }
        }
        //(dd|ss) integrals
        Dd[0] = Dvalue(atmC,6);
        intn[0] = d_eri2Center(2,0,2,0,0,0,0,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(2,1,2,1,0,0,0,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[2] = d_eri2Center(2,2,2,2,0,0,0,0,RCD,Dd,atmC,atmD,type);    //(dl,dl|sg,sg)=(-dl,-dl|sg,sg)
        for (size_t iddd = 1; iddd < 16; ++iddd) {
          block(30 + iddd,1) = intn[0]*DDrot_dAB(iddd,1) + intn[1]*(DDrot_dAB(iddd,2) + DDrot_dAB(iddd,3)) + intn[2]*(DDrot_dAB(iddd,4) + DDrot_dAB(iddd,5));
        }
      }
    }
    if ((ncols > 1)&&(type == 0)) {
      intn[0] = eri2Center(0,0,0,0,0,0,1,0,RCD,D,atmC,atmD,type);
      intn[1] = eri2Center(0,0,0,0,1,0,1,0,RCD,D,atmC,atmD,type);
      intn[2] = eri2Center(0,0,0,0,1,1,1,1,RCD,D,atmC,atmD,type);
      for (size_t idx = 0; idx < 6; ++idx) {
        if (idx < 3) {block(1,idx + 2) = -intn[0]*SProt_dAB(idx + 1,1);}                                                           //(ss|sp) integrals
        block(1,5 + idx) = intn[1]*PProt_dAB(idx + 1,1) + intn[2]*(PProt_dAB(idx + 1,2) + PProt_dAB(idx + 1,3));   //(ss|pp) integrals
      }
      if (ncols > 10) {
        //(ss|sd) integrals
        Dd[0] = Dvalue(atmC,1);
        Dd[1] = Dvalue(atmD,4);
        intn[0] = d_eri2Center(0,0,0,0,0,0,2,0,RCD,Dd,atmC,atmD,type);
        for (size_t idsd = 1; idsd < 6; ++idsd) {
          block(1,10 + idsd) = intn[0]*SDrot_dAB(idsd,1);
        }
        //(ss|pd) integrals
        Dd[1] = Dvalue(atmD,5);
        intn[0] = d_eri2Center(0,0,0,0,1,0,2,0,RCD,Dd,atmC,atmD,type);
        intn[1] = d_eri2Center(0,0,0,0,1,1,2,1,RCD,Dd,atmC,atmD,type);
        for (int idpd1 = 1; idpd1 < 6; ++idpd1) {
          for (int idpd2 = 1; idpd2 < 4; ++idpd2) {
            block(1,5*(idpd2 + 2) + idpd1) = -intn[0]*PDrot_dAB((idpd1 - 1)*3 + idpd2,1) - intn[1]*(PDrot_dAB((idpd1 - 1)*3 + idpd2,2) + PDrot_dAB((idpd1 - 1)*3 + idpd2,3));
          }
        }
        //(ss|dd) integrals
        Dd[1] = Dvalue(atmD,6);
        intn[0] = d_eri2Center(0,0,0,0,2,0,2,0,RCD,Dd,atmC,atmD,type);
        intn[1] = d_eri2Center(0,0,0,0,2,1,2,1,RCD,Dd,atmC,atmD,type);
        intn[2] = d_eri2Center(0,0,0,0,2,2,2,2,RCD,Dd,atmC,atmD,type);
        for (size_t iddd = 1; iddd < 16; ++iddd) {
          block(1,30 + iddd) = intn[0]*DDrot_dAB(iddd,1) + intn[1]*(DDrot_dAB(iddd,2) + DDrot_dAB(iddd,3)) + intn[2]*(DDrot_dAB(iddd,4) + DDrot_dAB(iddd,5));
        }
      }
    }
    if ((nrows > 1)&&(ncols > 1)&&(type == 0)) {
      intn[0] = eri2Center(0,0,1,0,0,0,1,0,RCD,D,atmC,atmD,type);
      intn[1] = eri2Center(0,0,1,1,0,0,1,1,RCD,D,atmC,atmD,type);
      intn[2] = eri2Center(0,0,1,0,1,0,1,0,RCD,D,atmC,atmD,type);
      intn[3] = eri2Center(0,0,1,0,1,1,1,1,RCD,D,atmC,atmD,type);
      intn[4] = eri2Center(0,0,1,1,1,0,1,1,RCD,D,atmC,atmD,type);
      for (size_t idpb = 1; idpb < 4; ++idpb) {
        for (size_t idpk = 1; idpk < 4; ++idpk) {           //(sp|sp)
          block(idpb + 1,idpk + 1) = intn[0]*SProt_dAB(idpb,1)*SProt(idpk,1) + intn[1]*(SProt_dAB(idpb,2)*SProt(idpk,2) + SProt_dAB(idpb,3)*SProt(idpk,3));
          block(idpb + 1,idpk + 1) += intn[0]*SProt_dA(idpb,1)*SProt_dB(idpk,1) + intn[1]*(SProt_dA(idpb,2)*SProt_dB(idpk,2) + SProt_dA(idpb,3)*SProt_dB(idpk,3));
          block(idpb + 1,idpk + 1) += intn[0]*SProt_dB(idpb,1)*SProt_dA(idpk,1) + intn[1]*(SProt_dB(idpb,2)*SProt_dA(idpk,2) + SProt_dB(idpb,3)*SProt_dA(idpk,3));
          block(idpb + 1,idpk + 1) += intn[0]*SProt(idpb,1)*SProt_dAB(idpk,1) + intn[1]*(SProt(idpb,2)*SProt_dAB(idpk,2) + SProt(idpb,3)*SProt_dAB(idpk,3));
        }
        for (size_t idc = 1; idc < 7; ++idc) {              //(sp|pp)
          block(idpb + 1,idc + 4) = -SProt_dAB(idpb,1)*(intn[2]*PProt(idc,1) + intn[3]*(PProt(idc,2) + PProt(idc,3))) - intn[4]*(SProt_dAB(idpb,2)*PProt(idc,4) + SProt_dAB(idpb,3)*PProt(idc,5));
          block(idpb + 1,idc + 4) -= SProt_dA(idpb,1)*(intn[2]*PProt_dB(idc,1) + intn[3]*(PProt_dB(idc,2) + PProt_dB(idc,3))) + intn[4]*(SProt_dA(idpb,2)*PProt_dB(idc,4) + SProt_dA(idpb,3)*PProt_dB(idc,5));
          block(idpb + 1,idc + 4) -= SProt_dB(idpb,1)*(intn[2]*PProt_dA(idc,1) + intn[3]*(PProt_dA(idc,2) + PProt_dA(idc,3))) + intn[4]*(SProt_dB(idpb,2)*PProt_dA(idc,4) + SProt_dB(idpb,3)*PProt_dA(idc,5));
          block(idpb + 1,idc + 4) -= SProt(idpb,1)*(intn[2]*PProt_dAB(idc,1) + intn[3]*(PProt_dAB(idc,2) + PProt_dAB(idc,3))) + intn[4]*(SProt(idpb,2)*PProt_dAB(idc,4) + SProt(idpb,3)*PProt_dAB(idc,5));
        }
      }
      intn[0] = eri2Center(1,0,1,0,1,0,1,0,RCD,D,atmC,atmD,type);                          //(ps,ps|ps,ps)
      intn[1] = eri2Center(1,1,1,1,1,0,1,0,RCD,D,atmC,atmD,type);                          //(pp,pp|ps,ps)
      intn[2] = eri2Center(1,0,1,0,1,1,1,1,RCD,D,atmC,atmD,type);                          //(ps,ps|pp,pp)
      intn[3] = eri2Center(1,1,1,1,1,1,1,1,RCD,D,atmC,atmD,type);                          //(pp,pp|pp,pp)
      intn[4] = eri2Center(1,1,1,0,1,1,1,0,RCD,D,atmC,atmD,type);                          //(pp,ps|pp,ps)
      intn[5] = eri2Center(1,1,1,1,1,-1,1,-1,RCD,D,atmC,atmD,type);                        //(pp,pp|pp*,pp*)
      intn[6] = eri2Center(1,-1,1,1,1,-1,1,1,RCD,D,atmC,atmD,type);                        //(pp*,pp|pp*,pp)
      intn[7] = eri2Center(1,0,1,0,0,0,1,0,RCD,D,atmC,atmD,type);
      intn[8] = eri2Center(1,1,1,1,0,0,1,0,RCD,D,atmC,atmD,type);
      intn[9] = eri2Center(1,0,1,1,0,0,1,1,RCD,D,atmC,atmD,type);
      for (size_t idr = 1; idr < 7; ++idr) {
        for (size_t idc = 1; idc < 4; ++idc) {                //(pp|sp)
          block(idr + 4,idc + 1) = -SProt_dAB(idc,1)*(intn[7]*PProt(idr,1) + intn[8]*(PProt(idr,2) + PProt(idr,3))) - intn[9]*(SProt_dAB(idc,2)*PProt(idr,4) + SProt_dAB(idc,3)*PProt(idr,5));
          block(idr + 4,idc + 1) -= SProt_dA(idc,1)*(intn[7]*PProt_dB(idr,1) + intn[8]*(PProt_dB(idr,2) + PProt_dB(idr,3))) + intn[9]*(SProt_dA(idc,2)*PProt_dB(idr,4) + SProt_dA(idc,3)*PProt_dB(idr,5));
          block(idr + 4,idc + 1) -= SProt_dB(idc,1)*(intn[7]*PProt_dA(idr,1) + intn[8]*(PProt_dA(idr,2) + PProt_dA(idr,3))) + intn[9]*(SProt_dB(idc,2)*PProt_dA(idr,4) + SProt_dB(idc,3)*PProt_dA(idr,5));
          block(idr + 4,idc + 1) -= SProt(idc,1)*(intn[7]*PProt_dAB(idr,1) + intn[8]*(PProt_dAB(idr,2) + PProt_dAB(idr,3))) + intn[9]*(SProt(idc,2)*PProt_dAB(idr,4) + SProt(idc,3)*PProt_dAB(idr,5));
        }
        for (size_t idc = 1; idc < 7; ++idc) {                //(pp|pp)
          block(idr + 4,idc + 4) = intn[0]*PProt_dAB(idr,1)*PProt(idc,1) + intn[1]*(PProt_dAB(idr,2) + PProt_dAB(idr,3))*PProt(idc,1) + intn[2]*PProt_dAB(idr,1)*(PProt(idc,2) + PProt(idc,3)) + intn[3]*(PProt_dAB(idr,2)*PProt(idc,2) + PProt_dAB(idr,3)*PProt(idc,3)) + intn[4]*(PProt_dAB(idr,4)*PProt(idc,4) + PProt_dAB(idr,5)*PProt(idc,5)) + intn[5]*(PProt_dAB(idr,2)*PProt(idc,3) + PProt_dAB(idr,3)*PProt(idc,2)) + intn[6]*PProt_dAB(idr,6)*PProt(idc,6);
          block(idr + 4,idc + 4) += intn[0]*PProt_dA(idr,1)*PProt_dB(idc,1) + intn[1]*(PProt_dA(idr,2) + PProt_dA(idr,3))*PProt_dB(idc,1) + intn[2]*PProt_dA(idr,1)*(PProt_dB(idc,2) + PProt_dB(idc,3)) + intn[3]*(PProt_dA(idr,2)*PProt_dB(idc,2) + PProt_dA(idr,3)*PProt_dB(idc,3)) + intn[4]*(PProt_dA(idr,4)*PProt_dB(idc,4) + PProt_dA(idr,5)*PProt_dB(idc,5)) + intn[5]*(PProt_dA(idr,2)*PProt_dB(idc,3) + PProt_dA(idr,3)*PProt_dB(idc,2)) + intn[6]*PProt_dA(idr,6)*PProt_dB(idc,6);
          block(idr + 4,idc + 4) += intn[0]*PProt_dB(idr,1)*PProt_dA(idc,1) + intn[1]*(PProt_dB(idr,2) + PProt_dB(idr,3))*PProt_dA(idc,1) + intn[2]*PProt_dB(idr,1)*(PProt_dA(idc,2) + PProt_dA(idc,3)) + intn[3]*(PProt_dB(idr,2)*PProt_dA(idc,2) + PProt_dB(idr,3)*PProt_dA(idc,3)) + intn[4]*(PProt_dB(idr,4)*PProt_dA(idc,4) + PProt_dB(idr,5)*PProt_dA(idc,5)) + intn[5]*(PProt_dB(idr,2)*PProt_dA(idc,3) + PProt_dB(idr,3)*PProt_dA(idc,2)) + intn[6]*PProt_dB(idr,6)*PProt_dA(idc,6);
          block(idr + 4,idc + 4) += intn[0]*PProt(idr,1)*PProt_dAB(idc,1) + intn[1]*(PProt(idr,2) + PProt(idr,3))*PProt_dAB(idc,1) + intn[2]*PProt(idr,1)*(PProt_dAB(idc,2) + PProt_dAB(idc,3)) + intn[3]*(PProt(idr,2)*PProt_dAB(idc,2) + PProt(idr,3)*PProt_dAB(idc,3)) + intn[4]*(PProt(idr,4)*PProt_dAB(idc,4) + PProt(idr,5)*PProt_dAB(idc,5)) + intn[5]*(PProt(idr,2)*PProt_dAB(idc,3) + PProt(idr,3)*PProt_dAB(idc,2)) + intn[6]*PProt(idr,6)*PProt_dAB(idc,6);
        }
      }
      if (nrows > 10) {
        //(sd|sp) integrals
        Dd[0] = Dvalue(atmC,4);
        Dd[1] = Dvalue(atmD,1);
        intn[0] = d_eri2Center(0,0,2,0,0,0,1,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(0,0,2,1,0,0,1,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        //(sd|sp) rotate
        for (size_t idsd = 1; idsd < 6; ++idsd) {
          for (size_t idsp = 1; idsp < 4; ++idsp) {
            block(10 + idsd,1 + idsp) = -intn[0]*SDrot_dAB(idsd,1)*SProt(idsp,1) - intn[1]*(SDrot_dAB(idsd,2)*SProt(idsp,2) + SDrot_dAB(idsd,3)*SProt(idsp,3));
            block(10 + idsd,1 + idsp) -= intn[0]*SDrot_dA(idsd,1)*SProt_dB(idsp,1) + intn[1]*(SDrot_dA(idsd,2)*SProt_dB(idsp,2) + SDrot_dA(idsd,3)*SProt_dB(idsp,3));
            block(10 + idsd,1 + idsp) -= intn[0]*SDrot_dB(idsd,1)*SProt_dA(idsp,1) + intn[1]*(SDrot_dB(idsd,2)*SProt_dA(idsp,2) + SDrot_dB(idsd,3)*SProt_dA(idsp,3));
            block(10 + idsd,1 + idsp) -= intn[0]*SDrot(idsd,1)*SProt_dAB(idsp,1) + intn[1]*(SDrot(idsd,2)*SProt_dAB(idsp,2) + SDrot(idsd,3)*SProt_dAB(idsp,3));
          }
        }
        //(sd|pp) integrals
        Dd[0] = Dvalue(atmC,4);
        Dd[1] = Dvalue(atmD,3);
        intn[0] = d_eri2Center(0,0,2,0,1,0,1,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(0,0,2,0,1,1,1,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center(0,0,2,1,1,0,1,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[3] = d_eri2Center(0,0,2,2,1,1,1,1,RCD,Dd,atmC,atmD,type);    //(sg,dl|pi,pi)=-(sg,dl|-pi,-pi)=(sg,-dl|pi,-pi)
        //(sd|pp) rotate
        for (size_t idsd = 1; idsd < 6; ++idsd) {
          for (size_t idpp = 1; idpp < 7; ++idpp) {
            block(10 + idsd,4 + idpp) = intn[0]*SDrot_dAB(idsd,1)*PProt(idpp,1) + intn[1]*SDrot_dAB(idsd,1)*(PProt(idpp,2) + PProt(idpp,3));
            block(10 + idsd,4 + idpp) += intn[0]*SDrot_dA(idsd,1)*PProt_dB(idpp,1) + intn[1]*SDrot_dA(idsd,1)*(PProt_dB(idpp,2) + PProt_dB(idpp,3));
            block(10 + idsd,4 + idpp) += intn[0]*SDrot_dB(idsd,1)*PProt_dA(idpp,1) + intn[1]*SDrot_dB(idsd,1)*(PProt_dA(idpp,2) + PProt_dA(idpp,3));
            block(10 + idsd,4 + idpp) += intn[0]*SDrot(idsd,1)*PProt_dAB(idpp,1) + intn[1]*SDrot(idsd,1)*(PProt_dAB(idpp,2) + PProt_dAB(idpp,3));
            block(10 + idsd,4 + idpp) += intn[2]*(SDrot_dAB(idsd,2)*PProt(idpp,4) + SDrot_dAB(idsd,3)*PProt(idpp,5));
            block(10 + idsd,4 + idpp) += intn[2]*(SDrot_dA(idsd,2)*PProt_dB(idpp,4) + SDrot_dA(idsd,3)*PProt_dB(idpp,5));
            block(10 + idsd,4 + idpp) += intn[2]*(SDrot_dB(idsd,2)*PProt_dA(idpp,4) + SDrot_dB(idsd,3)*PProt_dA(idpp,5));
            block(10 + idsd,4 + idpp) += intn[2]*(SDrot(idsd,2)*PProt_dAB(idpp,4) + SDrot(idsd,3)*PProt_dAB(idpp,5));
            block(10 + idsd,4 + idpp) += intn[3]*(SDrot_dAB(idsd,4)*PProt(idpp,2) - SDrot_dAB(idsd,4)*PProt(idpp,3) + SDrot_dAB(idsd,5)*PProt(idpp,6));
            block(10 + idsd,4 + idpp) += intn[3]*(SDrot_dA(idsd,4)*PProt_dB(idpp,2) - SDrot_dA(idsd,4)*PProt_dB(idpp,3) + SDrot_dA(idsd,5)*PProt_dB(idpp,6));
            block(10 + idsd,4 + idpp) += intn[3]*(SDrot_dB(idsd,4)*PProt_dA(idpp,2) - SDrot_dB(idsd,4)*PProt_dA(idpp,3) + SDrot_dB(idsd,5)*PProt_dA(idpp,6));
            block(10 + idsd,4 + idpp) += intn[3]*(SDrot(idsd,4)*PProt_dAB(idpp,2) - SDrot(idsd,4)*PProt_dAB(idpp,3) + SDrot(idsd,5)*PProt_dAB(idpp,6));
          }
        }
        //(pd|sp) integrals
        Dd[0] = Dvalue(atmC,5);
        Dd[1] = Dvalue(atmD,1);
        intn[0] = d_eri2Center(1,0,2,0,0,0,1,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(1,1,2,1,0,0,1,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[2] = d_eri2Center(1,0,2,1,0,0,1,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[3] = d_eri2Center(1,1,2,0,0,0,1,1,RCD,Dd,atmC,atmD,type);    //(pi,sg|sg,pi)=(-pi,sg|sg,-pi)
        intn[4] = d_eri2Center(1,1,2,2,0,0,1,1,RCD,Dd,atmC,atmD,type);    //(pi,dl|sg,pi)=(-pi,-dl|sg,pi)=-(-pi,dl|sg,-pi)=(pi,-dl|sg,-pi)
        //(pd|sp) rotate
        for (size_t idsp = 1; idsp < 4; ++idsp) {
          for (size_t idpd1 = 1; idpd1 < 6; ++idpd1) {
            for (size_t idpd2 = 1; idpd2 < 4; ++idpd2) {
              block(5*(idpd2+2) + idpd1,1 + idsp) = intn[0]*PDrot_dAB((idpd1 - 1)*3 + idpd2,1)*SProt(idsp,1) + intn[1]*(PDrot_dAB((idpd1 - 1)*3 + idpd2,2) + PDrot_dAB((idpd1 - 1)*3 + idpd2,3))*SProt(idsp,1);
              block(5*(idpd2+2) + idpd1,1 + idsp) += intn[0]*PDrot_dA((idpd1 - 1)*3 + idpd2,1)*SProt_dB(idsp,1) + intn[1]*(PDrot_dA((idpd1 - 1)*3 + idpd2,2) + PDrot_dA((idpd1 - 1)*3 + idpd2,3))*SProt_dB(idsp,1);
              block(5*(idpd2+2) + idpd1,1 + idsp) += intn[0]*PDrot_dB((idpd1 - 1)*3 + idpd2,1)*SProt_dA(idsp,1) + intn[1]*(PDrot_dB((idpd1 - 1)*3 + idpd2,2) + PDrot_dB((idpd1 - 1)*3 + idpd2,3))*SProt_dA(idsp,1);
              block(5*(idpd2+2) + idpd1,1 + idsp) += intn[0]*PDrot((idpd1 - 1)*3 + idpd2,1)*SProt_dAB(idsp,1) + intn[1]*(PDrot((idpd1 - 1)*3 + idpd2,2) + PDrot((idpd1 - 1)*3 + idpd2,3))*SProt_dAB(idsp,1);
              block(5*(idpd2+2) + idpd1,1 + idsp) += intn[2]*(PDrot_dAB((idpd1 - 1)*3 + idpd2,4)*SProt(idsp,2) + PDrot_dAB((idpd1 - 1)*3 + idpd2,5)*SProt(idsp,3));
              block(5*(idpd2+2) + idpd1,1 + idsp) += intn[2]*(PDrot_dA((idpd1 - 1)*3 + idpd2,4)*SProt_dB(idsp,2) + PDrot_dA((idpd1 - 1)*3 + idpd2,5)*SProt_dB(idsp,3));
              block(5*(idpd2+2) + idpd1,1 + idsp) += intn[2]*(PDrot_dB((idpd1 - 1)*3 + idpd2,4)*SProt_dA(idsp,2) + PDrot_dB((idpd1 - 1)*3 + idpd2,5)*SProt_dA(idsp,3));
              block(5*(idpd2+2) + idpd1,1 + idsp) += intn[2]*(PDrot((idpd1 - 1)*3 + idpd2,4)*SProt_dAB(idsp,2) + PDrot((idpd1 - 1)*3 + idpd2,5)*SProt_dAB(idsp,3));
              block(5*(idpd2+2) + idpd1,1 + idsp) += intn[3]*(PDrot_dAB((idpd1 - 1)*3 + idpd2,8)*SProt(idsp,2) + PDrot_dAB((idpd1 - 1)*3 + idpd2,12)*SProt(idsp,3));
              block(5*(idpd2+2) + idpd1,1 + idsp) += intn[3]*(PDrot_dA((idpd1 - 1)*3 + idpd2,8)*SProt_dB(idsp,2) + PDrot_dA((idpd1 - 1)*3 + idpd2,12)*SProt_dB(idsp,3));
              block(5*(idpd2+2) + idpd1,1 + idsp) += intn[3]*(PDrot_dB((idpd1 - 1)*3 + idpd2,8)*SProt_dA(idsp,2) + PDrot_dB((idpd1 - 1)*3 + idpd2,12)*SProt_dA(idsp,3));
              block(5*(idpd2+2) + idpd1,1 + idsp) += intn[3]*(PDrot((idpd1 - 1)*3 + idpd2,8)*SProt_dAB(idsp,2) + PDrot((idpd1 - 1)*3 + idpd2,12)*SProt_dAB(idsp,3));
              block(5*(idpd2+2) + idpd1,1 + idsp) += intn[4]*((PDrot_dAB((idpd1 - 1)*3 + idpd2,10) + PDrot_dAB((idpd1 - 1)*3 + idpd2,15))*SProt(idsp,2) + (PDrot_dAB((idpd1 - 1)*3 + idpd2,11) - PDrot_dAB((idpd1 - 1)*3 + idpd2,14))*SProt(idsp,3));
              block(5*(idpd2+2) + idpd1,1 + idsp) += intn[4]*((PDrot_dA((idpd1 - 1)*3 + idpd2,10) + PDrot_dA((idpd1 - 1)*3 + idpd2,15))*SProt_dB(idsp,2) + (PDrot_dA((idpd1 - 1)*3 + idpd2,11) - PDrot_dA((idpd1 - 1)*3 + idpd2,14))*SProt_dB(idsp,3));
              block(5*(idpd2+2) + idpd1,1 + idsp) += intn[4]*((PDrot_dB((idpd1 - 1)*3 + idpd2,10) + PDrot_dB((idpd1 - 1)*3 + idpd2,15))*SProt_dA(idsp,2) + (PDrot_dB((idpd1 - 1)*3 + idpd2,11) - PDrot_dB((idpd1 - 1)*3 + idpd2,14))*SProt_dA(idsp,3));
              block(5*(idpd2+2) + idpd1,1 + idsp) += intn[4]*((PDrot((idpd1 - 1)*3 + idpd2,10) + PDrot((idpd1 - 1)*3 + idpd2,15))*SProt_dAB(idsp,2) + (PDrot((idpd1 - 1)*3 + idpd2,11) - PDrot((idpd1 - 1)*3 + idpd2,14))*SProt_dAB(idsp,3));
            }
          }
        }
        //(dd|sp) integrals
        Dd[0] = Dvalue(atmC,6);
        Dd[1] = Dvalue(atmD,1);
        intn[0] = d_eri2Center(2,0,2,0,0,0,1,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(2,1,2,1,0,0,1,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[2] = d_eri2Center(2,2,2,2,0,0,1,0,RCD,Dd,atmC,atmD,type);    //(dl,dl|sg,sg)=(-dl,-dl|sg,sg)
        intn[3] = d_eri2Center(2,0,2,1,0,0,1,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[4] = d_eri2Center(2,1,2,2,0,0,1,1,RCD,Dd,atmC,atmD,type);    //(pi,dl|sg,pi)=(-pi,-dl|sg,pi)=-(-pi,dl|sg,-pi)=(pi,-dl|sg,-pi)
        //(dd|sp) rotate
        for (size_t idsp = 1; idsp < 4; ++idsp) {
          for (size_t iddd = 1; iddd < 16; ++iddd) {
            block(30 + iddd,1 + idsp) = -intn[0]*DDrot_dAB(iddd,1)*SProt(idsp,1) - intn[1]*(DDrot_dAB(iddd,2) + DDrot_dAB(iddd,3))*SProt(idsp,1) - intn[2]*(DDrot_dAB(iddd,4) + DDrot_dAB(iddd,5))*SProt(idsp,1);
            block(30 + iddd,1 + idsp) -= intn[0]*DDrot_dA(iddd,1)*SProt_dB(idsp,1) + intn[1]*(DDrot_dA(iddd,2) + DDrot_dA(iddd,3))*SProt_dB(idsp,1) + intn[2]*(DDrot_dA(iddd,4) + DDrot_dA(iddd,5))*SProt_dB(idsp,1);
            block(30 + iddd,1 + idsp) -= intn[0]*DDrot_dB(iddd,1)*SProt_dA(idsp,1) + intn[1]*(DDrot_dB(iddd,2) + DDrot_dB(iddd,3))*SProt_dA(idsp,1) + intn[2]*(DDrot_dB(iddd,4) + DDrot_dB(iddd,5))*SProt_dA(idsp,1);
            block(30 + iddd,1 + idsp) -= intn[0]*DDrot(iddd,1)*SProt_dAB(idsp,1) + intn[1]*(DDrot(iddd,2) + DDrot(iddd,3))*SProt_dAB(idsp,1) + intn[2]*(DDrot(iddd,4) + DDrot(iddd,5))*SProt_dAB(idsp,1);
            block(30 + iddd,1 + idsp) -= intn[3]*(DDrot_dAB(iddd,6)*SProt(idsp,2) + DDrot_dAB(iddd,7)*SProt(idsp,3));
            block(30 + iddd,1 + idsp) -= intn[3]*(DDrot_dA(iddd,6)*SProt_dB(idsp,2) + DDrot_dA(iddd,7)*SProt_dB(idsp,3));
            block(30 + iddd,1 + idsp) -= intn[3]*(DDrot_dB(iddd,6)*SProt_dA(idsp,2) + DDrot_dB(iddd,7)*SProt_dA(idsp,3));
            block(30 + iddd,1 + idsp) -= intn[3]*(DDrot(iddd,6)*SProt_dAB(idsp,2) + DDrot(iddd,7)*SProt_dAB(idsp,3));
            block(30 + iddd,1 + idsp) -= intn[4]*((DDrot_dAB(iddd,11) + DDrot_dAB(iddd,14))*SProt(idsp,2) + (DDrot_dAB(iddd,12) - DDrot_dAB(iddd,13))*SProt(idsp,3));
            block(30 + iddd,1 + idsp) -= intn[4]*((DDrot_dA(iddd,11) + DDrot_dA(iddd,14))*SProt_dB(idsp,2) + (DDrot_dA(iddd,12) - DDrot_dA(iddd,13))*SProt_dB(idsp,3));
            block(30 + iddd,1 + idsp) -= intn[4]*((DDrot_dB(iddd,11) + DDrot_dB(iddd,14))*SProt_dA(idsp,2) + (DDrot_dB(iddd,12) - DDrot_dB(iddd,13))*SProt_dA(idsp,3));
            block(30 + iddd,1 + idsp) -= intn[4]*((DDrot(iddd,11) + DDrot(iddd,14))*SProt_dAB(idsp,2) + (DDrot(iddd,12) - DDrot(iddd,13))*SProt_dAB(idsp,3));
          }
        }
        //(pd|pp) integrals
        Dd[0] = Dvalue(atmC,5);
        Dd[1] = Dvalue(atmD,3);
        intn[0] = d_eri2Center(1,0,2,0,1,0,1,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(1,1,2,1,1,0,1,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[2] = d_eri2Center(1,0,2,0,1,1,1,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[3] = d_eri2Center(1,1,2,1,1,1,1,1,RCD,Dd,atmC,atmD,type);    //(pi,pi|pi,pi)=(pi,pi|-pi,-pi)=(-pi,-pi|pi,pi)=(-pi,-pi|-pi,-pi)
        intn[4] = d_eri2Center(1,0,2,1,1,0,1,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[5] = d_eri2Center(1,1,2,0,1,1,1,0,RCD,Dd,atmC,atmD,type);    //(pi,sg|pi,sg)=(-pi,sg|-pi,sg)
        intn[6] = d_eri2Center(1,1,2,2,1,0,1,1,RCD,Dd,atmC,atmD,type);    //(pi,dl|sg,pi)=(-pi,-dl|sg,pi)=-(-pi,dl|sg,-pi)=(pi,-dl|sg,-pi)
        //(pd|pp) rotate
        for (size_t idpp = 1; idpp < 7; ++idpp) {
          for (size_t idpd1 = 1; idpd1 < 6; ++idpd1) {
            for (size_t idpd2 = 1; idpd2 < 4; ++idpd2) {
              block(5*(idpd2+2) + idpd1,4 + idpp) = -intn[0]*PDrot_dAB((idpd1 - 1)*3 + idpd2,1)*PProt(idpp,1) - intn[1]*(PDrot_dAB((idpd1 - 1)*3 + idpd2,2) + PDrot_dAB((idpd1 - 1)*3 + idpd2,3))*PProt(idpp,1);
              block(5*(idpd2+2) + idpd1,4 + idpp) -= intn[0]*PDrot_dA((idpd1 - 1)*3 + idpd2,1)*PProt_dB(idpp,1) + intn[1]*(PDrot_dA((idpd1 - 1)*3 + idpd2,2) + PDrot_dA((idpd1 - 1)*3 + idpd2,3))*PProt_dB(idpp,1);
              block(5*(idpd2+2) + idpd1,4 + idpp) -= intn[0]*PDrot_dB((idpd1 - 1)*3 + idpd2,1)*PProt_dA(idpp,1) + intn[1]*(PDrot_dB((idpd1 - 1)*3 + idpd2,2) + PDrot_dB((idpd1 - 1)*3 + idpd2,3))*PProt_dA(idpp,1);
              block(5*(idpd2+2) + idpd1,4 + idpp) -= intn[0]*PDrot((idpd1 - 1)*3 + idpd2,1)*PProt_dAB(idpp,1) + intn[1]*(PDrot((idpd1 - 1)*3 + idpd2,2) + PDrot((idpd1 - 1)*3 + idpd2,3))*PProt_dAB(idpp,1);
              block(5*(idpd2+2) + idpd1,4 + idpp) -= intn[2]*PDrot_dAB((idpd1 - 1)*3 + idpd2,1)*(PProt(idpp,2) + PProt(idpp,3));
              block(5*(idpd2+2) + idpd1,4 + idpp) -= intn[2]*PDrot_dA((idpd1 - 1)*3 + idpd2,1)*(PProt_dB(idpp,2) + PProt_dB(idpp,3));
              block(5*(idpd2+2) + idpd1,4 + idpp) -= intn[2]*PDrot_dB((idpd1 - 1)*3 + idpd2,1)*(PProt_dA(idpp,2) + PProt_dA(idpp,3));
              block(5*(idpd2+2) + idpd1,4 + idpp) -= intn[2]*PDrot((idpd1 - 1)*3 + idpd2,1)*(PProt_dAB(idpp,2) + PProt_dAB(idpp,3));
              block(5*(idpd2+2) + idpd1,4 + idpp) -= intn[3]*(PDrot_dAB((idpd1 - 1)*3 + idpd2,2) + PDrot_dAB((idpd1 - 1)*3 + idpd2,3))*(PProt(idpp,2) + PProt(idpp,3));
              block(5*(idpd2+2) + idpd1,4 + idpp) -= intn[3]*(PDrot_dA((idpd1 - 1)*3 + idpd2,2) + PDrot_dA((idpd1 - 1)*3 + idpd2,3))*(PProt_dB(idpp,2) + PProt_dB(idpp,3));
              block(5*(idpd2+2) + idpd1,4 + idpp) -= intn[3]*(PDrot_dB((idpd1 - 1)*3 + idpd2,2) + PDrot_dB((idpd1 - 1)*3 + idpd2,3))*(PProt_dA(idpp,2) + PProt_dA(idpp,3));
              block(5*(idpd2+2) + idpd1,4 + idpp) -= intn[3]*(PDrot((idpd1 - 1)*3 + idpd2,2) + PDrot((idpd1 - 1)*3 + idpd2,3))*(PProt_dAB(idpp,2) + PProt_dAB(idpp,3));
              block(5*(idpd2+2) + idpd1,4 + idpp) -= intn[4]*(PDrot_dAB((idpd1 - 1)*3 + idpd2,4)*PProt(idpp,4) + PDrot_dAB((idpd1 - 1)*3 + idpd2,5)*PProt(idpp,5));
              block(5*(idpd2+2) + idpd1,4 + idpp) -= intn[4]*(PDrot_dA((idpd1 - 1)*3 + idpd2,4)*PProt_dB(idpp,4) + PDrot_dA((idpd1 - 1)*3 + idpd2,5)*PProt_dB(idpp,5));
              block(5*(idpd2+2) + idpd1,4 + idpp) -= intn[4]*(PDrot_dB((idpd1 - 1)*3 + idpd2,4)*PProt_dA(idpp,4) + PDrot_dB((idpd1 - 1)*3 + idpd2,5)*PProt_dA(idpp,5));
              block(5*(idpd2+2) + idpd1,4 + idpp) -= intn[4]*(PDrot((idpd1 - 1)*3 + idpd2,4)*PProt_dAB(idpp,4) + PDrot((idpd1 - 1)*3 + idpd2,5)*PProt_dAB(idpp,5));
              block(5*(idpd2+2) + idpd1,4 + idpp) -= intn[5]*(PDrot_dAB((idpd1 - 1)*3 + idpd2,8)*PProt(idpp,4) + PDrot_dAB((idpd1 - 1)*3 + idpd2,12)*PProt(idpp,5));
              block(5*(idpd2+2) + idpd1,4 + idpp) -= intn[5]*(PDrot_dA((idpd1 - 1)*3 + idpd2,8)*PProt_dB(idpp,4) + PDrot_dA((idpd1 - 1)*3 + idpd2,12)*PProt_dB(idpp,5));
              block(5*(idpd2+2) + idpd1,4 + idpp) -= intn[5]*(PDrot_dB((idpd1 - 1)*3 + idpd2,8)*PProt_dA(idpp,4) + PDrot_dB((idpd1 - 1)*3 + idpd2,12)*PProt_dA(idpp,5));
              block(5*(idpd2+2) + idpd1,4 + idpp) -= intn[5]*(PDrot((idpd1 - 1)*3 + idpd2,8)*PProt_dAB(idpp,4) + PDrot((idpd1 - 1)*3 + idpd2,12)*PProt_dAB(idpp,5));
              block(5*(idpd2+2) + idpd1,4 + idpp) -= intn[6]*((PDrot_dAB((idpd1 - 1)*3 + idpd2,10) + PDrot_dAB((idpd1 - 1)*3 + idpd2,15))*PProt(idpp,4) + (PDrot_dAB((idpd1 - 1)*3 + idpd2,11) - PDrot_dAB((idpd1 - 1)*3 + idpd2,14))*PProt(idpp,5));
              block(5*(idpd2+2) + idpd1,4 + idpp) -= intn[6]*((PDrot_dA((idpd1 - 1)*3 + idpd2,10) + PDrot_dA((idpd1 - 1)*3 + idpd2,15))*PProt_dB(idpp,4) + (PDrot_dA((idpd1 - 1)*3 + idpd2,11) - PDrot_dA((idpd1 - 1)*3 + idpd2,14))*PProt_dB(idpp,5));
              block(5*(idpd2+2) + idpd1,4 + idpp) -= intn[6]*((PDrot_dB((idpd1 - 1)*3 + idpd2,10) + PDrot_dB((idpd1 - 1)*3 + idpd2,15))*PProt_dA(idpp,4) + (PDrot_dB((idpd1 - 1)*3 + idpd2,11) - PDrot_dB((idpd1 - 1)*3 + idpd2,14))*PProt_dA(idpp,5));
              block(5*(idpd2+2) + idpd1,4 + idpp) -= intn[6]*((PDrot((idpd1 - 1)*3 + idpd2,10) + PDrot((idpd1 - 1)*3 + idpd2,15))*PProt_dAB(idpp,4) + (PDrot((idpd1 - 1)*3 + idpd2,11) - PDrot((idpd1 - 1)*3 + idpd2,14))*PProt_dAB(idpp,5));
            }
          }
        }
        //(dd|pp) integrals
        Dd[0] = Dvalue(atmC,6);
        Dd[1] = Dvalue(atmD,3);
        intn[0] = d_eri2Center(2,0,2,0,1,0,1,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(2,1,2,1,1,0,1,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[2] = d_eri2Center(2,2,2,2,1,0,1,0,RCD,Dd,atmC,atmD,type);    //(dl,dl|sg,sg)=(-dl,-dl|sg,sg)
        intn[3] = d_eri2Center(2,0,2,0,1,1,1,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[4] = d_eri2Center(2,1,2,1,1,1,1,1,RCD,Dd,atmC,atmD,type);    //(pi,pi|pi,pi)=(-pi,-pi|-pi,-pi)
        intn[5] = d_eri2Center(2,2,2,2,1,1,1,1,RCD,Dd,atmC,atmD,type);    //(dl,dl|pi,pi)=(-dl,-dl|pi,pi)=(dl,dl|-pi,-pi)=(-dl,-dl|-pi,-pi)
        intn[6] = d_eri2Center(2,0,2,1,1,0,1,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[7] = d_eri2Center(2,-1,2,-1,1,1,1,1,RCD,Dd,atmC,atmD,type);  //(-pi,-pi|pi,pi)=(pi,pi|-pi,-pi)
        intn[8] = d_eri2Center(2,1,2,-1,1,1,1,-1,RCD,Dd,atmC,atmD,type);  //(pi,-pi|pi,-pi)
        intn[9] = d_eri2Center(2,0,2,2,1,1,1,1,RCD,Dd,atmC,atmD,type);    //(sg,dl|pi,pi)=(sg,-dl|pi,-pi)=-(sg,dl|-pi,-pi)
        intn[10] = d_eri2Center(2,1,2,2,1,0,1,1,RCD,Dd,atmC,atmD,type);   //(pi,dl|sg,pi)=(-pi,-dl|sg,pi)=-(-pi,dl|sg,-pi)=(pi,-dl|sg,-pi)
        //(dd|pp) rotate
        for (size_t iddd = 1; iddd < 16; ++iddd) {
          for (size_t idpp = 1; idpp < 7; ++idpp) {
            block(30 + iddd,4 + idpp) = intn[0]*DDrot_dAB(iddd,1)*PProt(idpp,1) + intn[1]*(DDrot_dAB(iddd,2) + DDrot_dAB(iddd,3))*PProt(idpp,1) + intn[2]*(DDrot_dAB(iddd,4) + DDrot_dAB(iddd,5))*PProt(idpp,1);
            block(30 + iddd,4 + idpp) += intn[0]*DDrot_dA(iddd,1)*PProt_dB(idpp,1) + intn[1]*(DDrot_dA(iddd,2) + DDrot_dA(iddd,3))*PProt_dB(idpp,1) + intn[2]*(DDrot_dA(iddd,4) + DDrot_dA(iddd,5))*PProt_dB(idpp,1);
            block(30 + iddd,4 + idpp) += intn[0]*DDrot_dB(iddd,1)*PProt_dA(idpp,1) + intn[1]*(DDrot_dB(iddd,2) + DDrot_dB(iddd,3))*PProt_dA(idpp,1) + intn[2]*(DDrot_dB(iddd,4) + DDrot_dB(iddd,5))*PProt_dA(idpp,1);
            block(30 + iddd,4 + idpp) += intn[0]*DDrot(iddd,1)*PProt_dAB(idpp,1) + intn[1]*(DDrot(iddd,2) + DDrot(iddd,3))*PProt_dAB(idpp,1) + intn[2]*(DDrot(iddd,4) + DDrot(iddd,5))*PProt_dAB(idpp,1);
            block(30 + iddd,4 + idpp) += intn[3]*DDrot_dAB(iddd,1)*(PProt(idpp,2) + PProt(idpp,3)) + intn[4]*(DDrot_dAB(iddd,2)*PProt(idpp,2) + DDrot_dAB(iddd,3)*PProt(idpp,3));
            block(30 + iddd,4 + idpp) += intn[3]*DDrot_dA(iddd,1)*(PProt_dB(idpp,2) + PProt_dB(idpp,3)) + intn[4]*(DDrot_dA(iddd,2)*PProt_dB(idpp,2) + DDrot_dA(iddd,3)*PProt_dB(idpp,3));
            block(30 + iddd,4 + idpp) += intn[3]*DDrot_dB(iddd,1)*(PProt_dA(idpp,2) + PProt_dA(idpp,3)) + intn[4]*(DDrot_dB(iddd,2)*PProt_dA(idpp,2) + DDrot_dB(iddd,3)*PProt_dA(idpp,3));
            block(30 + iddd,4 + idpp) += intn[3]*DDrot(iddd,1)*(PProt_dAB(idpp,2) + PProt_dAB(idpp,3)) + intn[4]*(DDrot(iddd,2)*PProt_dAB(idpp,2) + DDrot(iddd,3)*PProt_dAB(idpp,3));
            block(30 + iddd,4 + idpp) += intn[5]*(DDrot_dAB(iddd,4) + DDrot_dAB(iddd,5))*(PProt(idpp,2) + PProt(idpp,3)) + intn[6]*(DDrot_dAB(iddd,6)*PProt(idpp,4) + DDrot_dAB(iddd,7)*PProt(idpp,5));
            block(30 + iddd,4 + idpp) += intn[5]*(DDrot_dA(iddd,4) + DDrot_dA(iddd,5))*(PProt_dB(idpp,2) + PProt_dB(idpp,3)) + intn[6]*(DDrot_dA(iddd,6)*PProt_dB(idpp,4) + DDrot_dA(iddd,7)*PProt_dB(idpp,5));
            block(30 + iddd,4 + idpp) += intn[5]*(DDrot_dB(iddd,4) + DDrot_dB(iddd,5))*(PProt_dA(idpp,2) + PProt_dA(idpp,3)) + intn[6]*(DDrot_dB(iddd,6)*PProt_dA(idpp,4) + DDrot_dB(iddd,7)*PProt_dA(idpp,5));
            block(30 + iddd,4 + idpp) += intn[5]*(DDrot(iddd,4) + DDrot(iddd,5))*(PProt_dAB(idpp,2) + PProt_dAB(idpp,3)) + intn[6]*(DDrot(iddd,6)*PProt_dAB(idpp,4) + DDrot(iddd,7)*PProt_dAB(idpp,5));
            block(30 + iddd,4 + idpp) += intn[7]*(DDrot_dAB(iddd,2)*PProt(idpp,3) + DDrot_dAB(iddd,3)*PProt(idpp,2)) + intn[8]*DDrot_dAB(iddd,10)*PProt(idpp,6);
            block(30 + iddd,4 + idpp) += intn[7]*(DDrot_dA(iddd,2)*PProt_dB(idpp,3) + DDrot_dA(iddd,3)*PProt_dB(idpp,2)) + intn[8]*DDrot_dA(iddd,10)*PProt_dB(idpp,6);
            block(30 + iddd,4 + idpp) += intn[7]*(DDrot_dB(iddd,2)*PProt_dA(idpp,3) + DDrot_dB(iddd,3)*PProt_dA(idpp,2)) + intn[8]*DDrot_dB(iddd,10)*PProt_dA(idpp,6);
            block(30 + iddd,4 + idpp) += intn[7]*(DDrot(iddd,2)*PProt_dAB(idpp,3) + DDrot(iddd,3)*PProt_dAB(idpp,2)) + intn[8]*DDrot(iddd,10)*PProt_dAB(idpp,6);
            block(30 + iddd,4 + idpp) += intn[9]*(DDrot_dAB(iddd,8)*(PProt(idpp,2) - PProt(idpp,3)) + DDrot_dAB(iddd,9)*PProt(idpp,6));
            block(30 + iddd,4 + idpp) += intn[9]*(DDrot_dA(iddd,8)*(PProt_dB(idpp,2) - PProt_dB(idpp,3)) + DDrot_dA(iddd,9)*PProt_dB(idpp,6));
            block(30 + iddd,4 + idpp) += intn[9]*(DDrot_dB(iddd,8)*(PProt_dA(idpp,2) - PProt_dA(idpp,3)) + DDrot_dB(iddd,9)*PProt_dA(idpp,6));
            block(30 + iddd,4 + idpp) += intn[9]*(DDrot(iddd,8)*(PProt_dAB(idpp,2) - PProt_dAB(idpp,3)) + DDrot(iddd,9)*PProt_dAB(idpp,6));
            block(30 + iddd,4 + idpp) += intn[10]*((DDrot_dAB(iddd,11) + DDrot_dAB(iddd,14))*PProt(idpp,4) + (DDrot_dAB(iddd,12) - DDrot_dAB(iddd,13))*PProt(idpp,5));
            block(30 + iddd,4 + idpp) += intn[10]*((DDrot_dA(iddd,11) + DDrot_dA(iddd,14))*PProt_dB(idpp,4) + (DDrot_dA(iddd,12) - DDrot_dA(iddd,13))*PProt_dB(idpp,5));
            block(30 + iddd,4 + idpp) += intn[10]*((DDrot_dB(iddd,11) + DDrot_dB(iddd,14))*PProt_dA(idpp,4) + (DDrot_dB(iddd,12) - DDrot_dB(iddd,13))*PProt_dA(idpp,5));
            block(30 + iddd,4 + idpp) += intn[10]*((DDrot(iddd,11) + DDrot(iddd,14))*PProt_dAB(idpp,4) + (DDrot(iddd,12) - DDrot(iddd,13))*PProt_dAB(idpp,5));
          }
        }
      }
      if (ncols > 10) {
        //(sp|sd) integrals
        Dd[0] = Dvalue(atmC,1);
        Dd[1] = Dvalue(atmD,4);
        intn[0] = d_eri2Center(0,0,1,0,0,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(0,0,1,1,0,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        //(sp|sd) rotate
        for (size_t idsp = 1; idsp < 4; ++idsp) {
          for (size_t idsd = 1; idsd < 6; ++idsd) {
            block(1 + idsp,10 + idsd) = -intn[0]*SDrot_dAB(idsd,1)*SProt(idsp,1) - intn[1]*(SDrot_dAB(idsd,2)*SProt(idsp,2) + SDrot_dAB(idsd,3)*SProt(idsp,3));
            block(1 + idsp,10 + idsd) -= intn[0]*SDrot_dA(idsd,1)*SProt_dB(idsp,1) + intn[1]*(SDrot_dA(idsd,2)*SProt_dB(idsp,2) + SDrot_dA(idsd,3)*SProt_dB(idsp,3));
            block(1 + idsp,10 + idsd) -= intn[0]*SDrot_dB(idsd,1)*SProt_dA(idsp,1) + intn[1]*(SDrot_dB(idsd,2)*SProt_dA(idsp,2) + SDrot_dB(idsd,3)*SProt_dA(idsp,3));
            block(1 + idsp,10 + idsd) -= intn[0]*SDrot(idsd,1)*SProt_dAB(idsp,1) + intn[1]*(SDrot(idsd,2)*SProt_dAB(idsp,2) + SDrot(idsd,3)*SProt_dAB(idsp,3));
          }
        }
        //(pp|sd) integrals
        Dd[0] = Dvalue(atmC,3);
        Dd[1] = Dvalue(atmD,4);
        intn[0] = d_eri2Center(1,0,1,0,0,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(1,1,1,1,0,0,2,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[2] = d_eri2Center(1,0,1,1,0,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[3] = d_eri2Center(1,1,1,1,0,0,2,2,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,dl)=-(-pi,-pi|sg,dl)=(pi,-pi|sg,-dl)
        //(pp|sd) rotate
        for (size_t idsd = 1; idsd < 6; ++idsd) {
          for (size_t idpp = 1; idpp < 7; ++idpp) {
            block(4 + idpp,10 + idsd) = intn[0]*SDrot_dAB(idsd,1)*PProt(idpp,1) + intn[1]*SDrot_dAB(idsd,1)*(PProt(idpp,2) + PProt(idpp,3));
            block(4 + idpp,10 + idsd) += intn[0]*SDrot_dA(idsd,1)*PProt_dB(idpp,1) + intn[1]*SDrot_dA(idsd,1)*(PProt_dB(idpp,2) + PProt_dB(idpp,3));
            block(4 + idpp,10 + idsd) += intn[0]*SDrot_dB(idsd,1)*PProt_dA(idpp,1) + intn[1]*SDrot_dB(idsd,1)*(PProt_dA(idpp,2) + PProt_dA(idpp,3));
            block(4 + idpp,10 + idsd) += intn[0]*SDrot(idsd,1)*PProt_dAB(idpp,1) + intn[1]*SDrot(idsd,1)*(PProt_dAB(idpp,2) + PProt_dAB(idpp,3));
            block(4 + idpp,10 + idsd) += intn[2]*(SDrot_dAB(idsd,2)*PProt(idpp,4) + SDrot_dAB(idsd,3)*PProt(idpp,5));
            block(4 + idpp,10 + idsd) += intn[2]*(SDrot_dA(idsd,2)*PProt_dB(idpp,4) + SDrot_dA(idsd,3)*PProt_dB(idpp,5));
            block(4 + idpp,10 + idsd) += intn[2]*(SDrot_dB(idsd,2)*PProt_dA(idpp,4) + SDrot_dB(idsd,3)*PProt_dA(idpp,5));
            block(4 + idpp,10 + idsd) += intn[2]*(SDrot(idsd,2)*PProt_dAB(idpp,4) + SDrot(idsd,3)*PProt_dAB(idpp,5));
            block(4 + idpp,10 + idsd) += intn[3]*(SDrot_dAB(idsd,4)*PProt(idpp,2) - SDrot_dAB(idsd,4)*PProt(idpp,3) + SDrot_dAB(idsd,5)*PProt(idpp,6));
            block(4 + idpp,10 + idsd) += intn[3]*(SDrot_dA(idsd,4)*PProt_dB(idpp,2) - SDrot_dA(idsd,4)*PProt_dB(idpp,3) + SDrot_dA(idsd,5)*PProt_dB(idpp,6));
            block(4 + idpp,10 + idsd) += intn[3]*(SDrot_dB(idsd,4)*PProt_dA(idpp,2) - SDrot_dB(idsd,4)*PProt_dA(idpp,3) + SDrot_dB(idsd,5)*PProt_dA(idpp,6));
            block(4 + idpp,10 + idsd) += intn[3]*(SDrot(idsd,4)*PProt_dAB(idpp,2) - SDrot(idsd,4)*PProt_dAB(idpp,3) + SDrot(idsd,5)*PProt_dAB(idpp,6));
          }
        }
        //(sp|pd) integrals
        Dd[0] = Dvalue(atmC,1);
        Dd[1] = Dvalue(atmD,5);
        intn[0] = d_eri2Center(0,0,1,0,1,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(0,0,1,0,1,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center(0,0,1,1,1,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[3] = d_eri2Center(0,0,1,1,1,1,2,0,RCD,Dd,atmC,atmD,type);    //(sg,pi|pi,sg)=(sg,-pi|-pi,sg)
        intn[4] = d_eri2Center(0,0,1,1,1,1,2,2,RCD,Dd,atmC,atmD,type);    //(sg,pi|pi,dl)=(sg,pi|-pi,-dl)=-(sg,-pi|-pi,dl)=(sg,-pi|pi,-dl)
        //(sp|pd) rotate
        for (size_t idsp = 1; idsp < 4; ++idsp) {
          for (size_t idpd1 = 1; idpd1 < 6; ++idpd1) {
            for (size_t idpd2 = 1; idpd2 < 4; ++idpd2) {
              block(1 + idsp,5*(idpd2+2) + idpd1) = intn[0]*PDrot_dAB((idpd1 - 1)*3 + idpd2,1)*SProt(idsp,1) + intn[1]*(PDrot_dAB((idpd1 - 1)*3 + idpd2,2) + PDrot_dAB((idpd1 - 1)*3 + idpd2,3))*SProt(idsp,1);
              block(1 + idsp,5*(idpd2+2) + idpd1) += intn[0]*PDrot_dA((idpd1 - 1)*3 + idpd2,1)*SProt_dB(idsp,1) + intn[1]*(PDrot_dA((idpd1 - 1)*3 + idpd2,2) + PDrot_dA((idpd1 - 1)*3 + idpd2,3))*SProt_dB(idsp,1);
              block(1 + idsp,5*(idpd2+2) + idpd1) += intn[0]*PDrot_dB((idpd1 - 1)*3 + idpd2,1)*SProt_dA(idsp,1) + intn[1]*(PDrot_dB((idpd1 - 1)*3 + idpd2,2) + PDrot_dB((idpd1 - 1)*3 + idpd2,3))*SProt_dA(idsp,1);
              block(1 + idsp,5*(idpd2+2) + idpd1) += intn[0]*PDrot((idpd1 - 1)*3 + idpd2,1)*SProt_dAB(idsp,1) + intn[1]*(PDrot((idpd1 - 1)*3 + idpd2,2) + PDrot((idpd1 - 1)*3 + idpd2,3))*SProt_dAB(idsp,1);
              block(1 + idsp,5*(idpd2+2) + idpd1) += intn[2]*(PDrot_dAB((idpd1 - 1)*3 + idpd2,4)*SProt(idsp,2) + PDrot_dAB((idpd1 - 1)*3 + idpd2,5)*SProt(idsp,3));
              block(1 + idsp,5*(idpd2+2) + idpd1) += intn[2]*(PDrot_dA((idpd1 - 1)*3 + idpd2,4)*SProt_dB(idsp,2) + PDrot_dA((idpd1 - 1)*3 + idpd2,5)*SProt_dB(idsp,3));
              block(1 + idsp,5*(idpd2+2) + idpd1) += intn[2]*(PDrot_dB((idpd1 - 1)*3 + idpd2,4)*SProt_dA(idsp,2) + PDrot_dB((idpd1 - 1)*3 + idpd2,5)*SProt_dA(idsp,3));
              block(1 + idsp,5*(idpd2+2) + idpd1) += intn[2]*(PDrot((idpd1 - 1)*3 + idpd2,4)*SProt_dAB(idsp,2) + PDrot((idpd1 - 1)*3 + idpd2,5)*SProt_dAB(idsp,3));
              block(1 + idsp,5*(idpd2+2) + idpd1) += intn[3]*(PDrot_dAB((idpd1 - 1)*3 + idpd2,8)*SProt(idsp,2) + PDrot_dAB((idpd1 - 1)*3 + idpd2,12)*SProt(idsp,3));
              block(1 + idsp,5*(idpd2+2) + idpd1) += intn[3]*(PDrot_dA((idpd1 - 1)*3 + idpd2,8)*SProt_dB(idsp,2) + PDrot_dA((idpd1 - 1)*3 + idpd2,12)*SProt_dB(idsp,3));
              block(1 + idsp,5*(idpd2+2) + idpd1) += intn[3]*(PDrot_dB((idpd1 - 1)*3 + idpd2,8)*SProt_dA(idsp,2) + PDrot_dB((idpd1 - 1)*3 + idpd2,12)*SProt_dA(idsp,3));
              block(1 + idsp,5*(idpd2+2) + idpd1) += intn[3]*(PDrot((idpd1 - 1)*3 + idpd2,8)*SProt_dAB(idsp,2) + PDrot((idpd1 - 1)*3 + idpd2,12)*SProt_dAB(idsp,3));
              block(1 + idsp,5*(idpd2+2) + idpd1) += intn[4]*((PDrot_dAB((idpd1 - 1)*3 + idpd2,10) + PDrot_dAB((idpd1 - 1)*3 + idpd2,15))*SProt(idsp,2) + (PDrot_dAB((idpd1 - 1)*3 + idpd2,11) - PDrot_dAB((idpd1 - 1)*3 + idpd2,14))*SProt(idsp,3));
              block(1 + idsp,5*(idpd2+2) + idpd1) += intn[4]*((PDrot_dA((idpd1 - 1)*3 + idpd2,10) + PDrot_dA((idpd1 - 1)*3 + idpd2,15))*SProt_dB(idsp,2) + (PDrot_dA((idpd1 - 1)*3 + idpd2,11) - PDrot_dA((idpd1 - 1)*3 + idpd2,14))*SProt_dB(idsp,3));
              block(1 + idsp,5*(idpd2+2) + idpd1) += intn[4]*((PDrot_dB((idpd1 - 1)*3 + idpd2,10) + PDrot_dB((idpd1 - 1)*3 + idpd2,15))*SProt_dA(idsp,2) + (PDrot_dB((idpd1 - 1)*3 + idpd2,11) - PDrot_dB((idpd1 - 1)*3 + idpd2,14))*SProt_dA(idsp,3));
              block(1 + idsp,5*(idpd2+2) + idpd1) += intn[4]*((PDrot((idpd1 - 1)*3 + idpd2,10) + PDrot((idpd1 - 1)*3 + idpd2,15))*SProt_dAB(idsp,2) + (PDrot((idpd1 - 1)*3 + idpd2,11) - PDrot((idpd1 - 1)*3 + idpd2,14))*SProt_dAB(idsp,3));
            }
          }
        }
        //(sp|dd) integrals
        Dd[0] = Dvalue(atmC,1);
        Dd[1] = Dvalue(atmD,6);
        intn[0] = d_eri2Center(0,0,1,0,2,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(0,0,1,0,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center(0,0,1,0,2,2,2,2,RCD,Dd,atmC,atmD,type);    //(sg,sg|dl,dl)=(sg,sg|-dl,-dl)
        intn[3] = d_eri2Center(0,0,1,1,2,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[4] = d_eri2Center(0,0,1,1,2,1,2,2,RCD,Dd,atmC,atmD,type);    //(sg,pi|pi,dl)=(sg,pi|-pi,-dl)=-(sg,-pi|-pi,dl)=(sg,-pi|pi,-dl)
        //(sp|dd) rotate
        for (size_t idsp = 1; idsp < 4; ++idsp) {
          for (size_t iddd = 1; iddd < 16; ++iddd) {
            block(1 + idsp,30 + iddd) = -intn[0]*DDrot_dAB(iddd,1)*SProt(idsp,1) - intn[1]*(DDrot_dAB(iddd,2) + DDrot_dAB(iddd,3))*SProt(idsp,1) - intn[2]*(DDrot_dAB(iddd,4) + DDrot_dAB(iddd,5))*SProt(idsp,1);
            block(1 + idsp,30 + iddd) -= intn[0]*DDrot_dA(iddd,1)*SProt_dB(idsp,1) + intn[1]*(DDrot_dA(iddd,2) + DDrot_dA(iddd,3))*SProt_dB(idsp,1) + intn[2]*(DDrot_dA(iddd,4) + DDrot_dA(iddd,5))*SProt_dB(idsp,1);
            block(1 + idsp,30 + iddd) -= intn[0]*DDrot_dB(iddd,1)*SProt_dA(idsp,1) + intn[1]*(DDrot_dB(iddd,2) + DDrot_dB(iddd,3))*SProt_dA(idsp,1) + intn[2]*(DDrot_dB(iddd,4) + DDrot_dB(iddd,5))*SProt_dA(idsp,1);
            block(1 + idsp,30 + iddd) -= intn[0]*DDrot(iddd,1)*SProt_dAB(idsp,1) + intn[1]*(DDrot(iddd,2) + DDrot(iddd,3))*SProt_dAB(idsp,1) + intn[2]*(DDrot(iddd,4) + DDrot(iddd,5))*SProt_dAB(idsp,1);
            block(1 + idsp,30 + iddd) -= intn[3]*(DDrot_dAB(iddd,6)*SProt(idsp,2) + DDrot_dAB(iddd,7)*SProt(idsp,3));
            block(1 + idsp,30 + iddd) -= intn[3]*(DDrot_dA(iddd,6)*SProt_dB(idsp,2) + DDrot_dA(iddd,7)*SProt_dB(idsp,3));
            block(1 + idsp,30 + iddd) -= intn[3]*(DDrot_dB(iddd,6)*SProt_dA(idsp,2) + DDrot_dB(iddd,7)*SProt_dA(idsp,3));
            block(1 + idsp,30 + iddd) -= intn[3]*(DDrot(iddd,6)*SProt_dAB(idsp,2) + DDrot(iddd,7)*SProt_dAB(idsp,3));
            block(1 + idsp,30 + iddd) -= intn[4]*((DDrot_dAB(iddd,11) + DDrot_dAB(iddd,14))*SProt(idsp,2) + (DDrot_dAB(iddd,12) - DDrot_dAB(iddd,13))*SProt(idsp,3));
            block(1 + idsp,30 + iddd) -= intn[4]*((DDrot_dA(iddd,11) + DDrot_dA(iddd,14))*SProt_dB(idsp,2) + (DDrot_dA(iddd,12) - DDrot_dA(iddd,13))*SProt_dB(idsp,3));
            block(1 + idsp,30 + iddd) -= intn[4]*((DDrot_dB(iddd,11) + DDrot_dB(iddd,14))*SProt_dA(idsp,2) + (DDrot_dB(iddd,12) - DDrot_dB(iddd,13))*SProt_dA(idsp,3));
            block(1 + idsp,30 + iddd) -= intn[4]*((DDrot(iddd,11) + DDrot(iddd,14))*SProt_dAB(idsp,2) + (DDrot(iddd,12) - DDrot(iddd,13))*SProt_dAB(idsp,3));
          }
        }
        //(pp|pd) integrals
        Dd[0] = Dvalue(atmC,3);
        Dd[1] = Dvalue(atmD,5);
        intn[0] = d_eri2Center(1,0,1,0,1,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(1,0,1,0,1,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center(1,1,1,1,1,0,2,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[3] = d_eri2Center(1,1,1,1,1,1,2,1,RCD,Dd,atmC,atmD,type);    //(pi,pi|pi,pi)=(pi,pi|-pi,-pi)=(-pi,-pi|pi,pi)=(-pi,-pi|-pi,-pi)
        intn[4] = d_eri2Center(1,0,1,1,1,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[5] = d_eri2Center(1,1,1,0,1,1,2,0,RCD,Dd,atmC,atmD,type);    //(pi,sg|pi,sg)=(-pi,sg|-pi,sg)
        intn[6] = d_eri2Center(1,0,1,1,1,1,2,2,RCD,Dd,atmC,atmD,type);    //(sg,pi|pi,dl)=(sg,pi|-pi,-dl)=-(sg,-pi|-pi,dl)=(sg,-pi|pi,-dl)
        //(pp|pd) rotate
        for (size_t idpp = 1; idpp < 7; ++idpp) {
          for (size_t idpd1 = 1; idpd1 < 6; ++idpd1) {
            for (size_t idpd2 = 1; idpd2 < 4; ++idpd2) {
              block(4 + idpp,5*(idpd2+2) + idpd1) = -intn[0]*PDrot_dAB((idpd1 - 1)*3 + idpd2,1)*PProt(idpp,1) - intn[1]*(PDrot_dAB((idpd1 - 1)*3 + idpd2,2) + PDrot_dAB((idpd1 - 1)*3 + idpd2,3))*PProt(idpp,1);
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[0]*PDrot_dA((idpd1 - 1)*3 + idpd2,1)*PProt_dB(idpp,1) + intn[1]*(PDrot_dA((idpd1 - 1)*3 + idpd2,2) + PDrot_dA((idpd1 - 1)*3 + idpd2,3))*PProt_dB(idpp,1);
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[0]*PDrot_dB((idpd1 - 1)*3 + idpd2,1)*PProt_dA(idpp,1) + intn[1]*(PDrot_dB((idpd1 - 1)*3 + idpd2,2) + PDrot_dB((idpd1 - 1)*3 + idpd2,3))*PProt_dA(idpp,1);
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[0]*PDrot((idpd1 - 1)*3 + idpd2,1)*PProt_dAB(idpp,1) + intn[1]*(PDrot((idpd1 - 1)*3 + idpd2,2) + PDrot((idpd1 - 1)*3 + idpd2,3))*PProt_dAB(idpp,1);
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[2]*PDrot_dAB((idpd1 - 1)*3 + idpd2,1)*(PProt(idpp,2) + PProt(idpp,3));
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[2]*PDrot_dA((idpd1 - 1)*3 + idpd2,1)*(PProt_dB(idpp,2) + PProt_dB(idpp,3));
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[2]*PDrot_dB((idpd1 - 1)*3 + idpd2,1)*(PProt_dA(idpp,2) + PProt_dA(idpp,3));
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[2]*PDrot((idpd1 - 1)*3 + idpd2,1)*(PProt_dAB(idpp,2) + PProt_dAB(idpp,3));
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[3]*(PDrot_dAB((idpd1 - 1)*3 + idpd2,2) + PDrot_dAB((idpd1 - 1)*3 + idpd2,3))*(PProt(idpp,2) + PProt(idpp,3));
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[3]*(PDrot_dA((idpd1 - 1)*3 + idpd2,2) + PDrot_dA((idpd1 - 1)*3 + idpd2,3))*(PProt_dB(idpp,2) + PProt_dB(idpp,3));
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[3]*(PDrot_dB((idpd1 - 1)*3 + idpd2,2) + PDrot_dB((idpd1 - 1)*3 + idpd2,3))*(PProt_dA(idpp,2) + PProt_dA(idpp,3));
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[3]*(PDrot((idpd1 - 1)*3 + idpd2,2) + PDrot((idpd1 - 1)*3 + idpd2,3))*(PProt_dAB(idpp,2) + PProt_dAB(idpp,3));
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[4]*(PDrot_dAB((idpd1 - 1)*3 + idpd2,4)*PProt(idpp,4) + PDrot_dAB((idpd1 - 1)*3 + idpd2,5)*PProt(idpp,5));
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[4]*(PDrot_dA((idpd1 - 1)*3 + idpd2,4)*PProt_dB(idpp,4) + PDrot_dA((idpd1 - 1)*3 + idpd2,5)*PProt_dB(idpp,5));
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[4]*(PDrot_dB((idpd1 - 1)*3 + idpd2,4)*PProt_dA(idpp,4) + PDrot_dB((idpd1 - 1)*3 + idpd2,5)*PProt_dA(idpp,5));
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[4]*(PDrot((idpd1 - 1)*3 + idpd2,4)*PProt_dAB(idpp,4) + PDrot((idpd1 - 1)*3 + idpd2,5)*PProt_dAB(idpp,5));
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[5]*(PDrot_dAB((idpd1 - 1)*3 + idpd2,8)*PProt(idpp,4) + PDrot_dAB((idpd1 - 1)*3 + idpd2,12)*PProt(idpp,5));
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[5]*(PDrot_dA((idpd1 - 1)*3 + idpd2,8)*PProt_dB(idpp,4) + PDrot_dA((idpd1 - 1)*3 + idpd2,12)*PProt_dB(idpp,5));
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[5]*(PDrot_dB((idpd1 - 1)*3 + idpd2,8)*PProt_dA(idpp,4) + PDrot_dB((idpd1 - 1)*3 + idpd2,12)*PProt_dA(idpp,5));
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[5]*(PDrot((idpd1 - 1)*3 + idpd2,8)*PProt_dAB(idpp,4) + PDrot((idpd1 - 1)*3 + idpd2,12)*PProt_dAB(idpp,5));
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[6]*((PDrot_dAB((idpd1 - 1)*3 + idpd2,10) + PDrot_dAB((idpd1 - 1)*3 + idpd2,15))*PProt(idpp,4) + (PDrot_dAB((idpd1 - 1)*3 + idpd2,11) - PDrot_dAB((idpd1 - 1)*3 + idpd2,14))*PProt(idpp,5));
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[6]*((PDrot_dA((idpd1 - 1)*3 + idpd2,10) + PDrot_dA((idpd1 - 1)*3 + idpd2,15))*PProt_dB(idpp,4) + (PDrot_dA((idpd1 - 1)*3 + idpd2,11) - PDrot_dA((idpd1 - 1)*3 + idpd2,14))*PProt_dB(idpp,5));
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[6]*((PDrot_dB((idpd1 - 1)*3 + idpd2,10) + PDrot_dB((idpd1 - 1)*3 + idpd2,15))*PProt_dA(idpp,4) + (PDrot_dB((idpd1 - 1)*3 + idpd2,11) - PDrot_dB((idpd1 - 1)*3 + idpd2,14))*PProt_dA(idpp,5));
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[6]*((PDrot((idpd1 - 1)*3 + idpd2,10) + PDrot((idpd1 - 1)*3 + idpd2,15))*PProt_dAB(idpp,4) + (PDrot((idpd1 - 1)*3 + idpd2,11) - PDrot((idpd1 - 1)*3 + idpd2,14))*PProt_dAB(idpp,5));
            }
          }
        }
        //(pp|dd) integrals
        Dd[0] = Dvalue(atmC,3);
        Dd[1] = Dvalue(atmD,6);
        intn[0] = d_eri2Center(1,0,1,0,2,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(1,0,1,0,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center(1,0,1,0,2,2,2,2,RCD,Dd,atmC,atmD,type);    //(sg,sg|dl,dl)=(sg,sg|-dl,-dl)
        intn[3] = d_eri2Center(1,1,1,1,2,0,2,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[4] = d_eri2Center(1,1,1,1,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(pi,pi|pi,pi)=(-pi,-pi|-pi,-pi)
        intn[5] = d_eri2Center(1,1,1,1,2,2,2,2,RCD,Dd,atmC,atmD,type);    //(pi,pi|dl,dl)=(pi,pi|-dl,-dl)=(-pi,-pi|dl,dl)=(-pi,-pi|-dl,-dl)
        intn[6] = d_eri2Center(1,0,1,1,2,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[7] = d_eri2Center(1,-1,1,-1,2,1,2,1,RCD,Dd,atmC,atmD,type);  //(-pi,-pi|pi,pi)=(pi,pi|-pi,-pi)
        intn[8] = d_eri2Center(1,1,1,-1,2,1,2,-1,RCD,Dd,atmC,atmD,type);  //(pi,-pi|pi,-pi)
        intn[9] = d_eri2Center(1,1,1,1,2,0,2,2,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,dl)=(pi,-pi|sg,-dl)=-(-pi,-pi|sg,dl)
        intn[10] = d_eri2Center(1,0,1,1,2,1,2,2,RCD,Dd,atmC,atmD,type);   //(sg,pi|pi,dl)=(sg,pi|-pi,-dl)=-(sg,-pi|-pi,dl)=(sg,-pi|pi,-dl)
        //(pp|dd) rotate
        for (size_t iddd = 1; iddd < 16; ++iddd) {
          for (size_t idpp = 1; idpp < 7; ++idpp) {
            block(4 + idpp,30 + iddd) = intn[0]*DDrot_dAB(iddd,1)*PProt(idpp,1) + intn[1]*(DDrot_dAB(iddd,2) + DDrot_dAB(iddd,3))*PProt(idpp,1) + intn[2]*(DDrot_dAB(iddd,4) + DDrot_dAB(iddd,5))*PProt(idpp,1);
            block(4 + idpp,30 + iddd) += intn[0]*DDrot_dA(iddd,1)*PProt_dB(idpp,1) + intn[1]*(DDrot_dA(iddd,2) + DDrot_dA(iddd,3))*PProt_dB(idpp,1) + intn[2]*(DDrot_dA(iddd,4) + DDrot_dA(iddd,5))*PProt_dB(idpp,1);
            block(4 + idpp,30 + iddd) += intn[0]*DDrot_dB(iddd,1)*PProt_dA(idpp,1) + intn[1]*(DDrot_dB(iddd,2) + DDrot_dB(iddd,3))*PProt_dA(idpp,1) + intn[2]*(DDrot_dB(iddd,4) + DDrot_dB(iddd,5))*PProt_dA(idpp,1);
            block(4 + idpp,30 + iddd) += intn[0]*DDrot(iddd,1)*PProt_dAB(idpp,1) + intn[1]*(DDrot(iddd,2) + DDrot(iddd,3))*PProt_dAB(idpp,1) + intn[2]*(DDrot(iddd,4) + DDrot(iddd,5))*PProt_dAB(idpp,1);
            block(4 + idpp,30 + iddd) += intn[3]*DDrot_dAB(iddd,1)*(PProt(idpp,2) + PProt(idpp,3)) + intn[4]*(DDrot_dAB(iddd,2)*PProt(idpp,2) + DDrot_dAB(iddd,3)*PProt(idpp,3));
            block(4 + idpp,30 + iddd) += intn[3]*DDrot_dA(iddd,1)*(PProt_dB(idpp,2) + PProt_dB(idpp,3)) + intn[4]*(DDrot_dA(iddd,2)*PProt_dB(idpp,2) + DDrot_dA(iddd,3)*PProt_dB(idpp,3));
            block(4 + idpp,30 + iddd) += intn[3]*DDrot_dB(iddd,1)*(PProt_dA(idpp,2) + PProt_dA(idpp,3)) + intn[4]*(DDrot_dB(iddd,2)*PProt_dA(idpp,2) + DDrot_dB(iddd,3)*PProt_dA(idpp,3));
            block(4 + idpp,30 + iddd) += intn[3]*DDrot(iddd,1)*(PProt_dAB(idpp,2) + PProt_dAB(idpp,3)) + intn[4]*(DDrot(iddd,2)*PProt_dAB(idpp,2) + DDrot(iddd,3)*PProt_dAB(idpp,3));
            block(4 + idpp,30 + iddd) += intn[5]*(DDrot_dAB(iddd,4) + DDrot_dAB(iddd,5))*(PProt(idpp,2) + PProt(idpp,3)) + intn[6]*(DDrot_dAB(iddd,6)*PProt(idpp,4) + DDrot_dAB(iddd,7)*PProt(idpp,5));
            block(4 + idpp,30 + iddd) += intn[5]*(DDrot_dA(iddd,4) + DDrot_dA(iddd,5))*(PProt_dB(idpp,2) + PProt_dB(idpp,3)) + intn[6]*(DDrot_dA(iddd,6)*PProt_dB(idpp,4) + DDrot_dA(iddd,7)*PProt_dB(idpp,5));
            block(4 + idpp,30 + iddd) += intn[5]*(DDrot_dB(iddd,4) + DDrot_dB(iddd,5))*(PProt_dA(idpp,2) + PProt_dA(idpp,3)) + intn[6]*(DDrot_dB(iddd,6)*PProt_dA(idpp,4) + DDrot_dB(iddd,7)*PProt_dA(idpp,5));
            block(4 + idpp,30 + iddd) += intn[5]*(DDrot(iddd,4) + DDrot(iddd,5))*(PProt_dAB(idpp,2) + PProt_dAB(idpp,3)) + intn[6]*(DDrot(iddd,6)*PProt_dAB(idpp,4) + DDrot(iddd,7)*PProt_dAB(idpp,5));
            block(4 + idpp,30 + iddd) += intn[7]*(DDrot_dAB(iddd,2)*PProt(idpp,3) + DDrot_dAB(iddd,3)*PProt(idpp,2)) + intn[8]*DDrot_dAB(iddd,10)*PProt(idpp,6);
            block(4 + idpp,30 + iddd) += intn[7]*(DDrot_dA(iddd,2)*PProt_dB(idpp,3) + DDrot_dA(iddd,3)*PProt_dB(idpp,2)) + intn[8]*DDrot_dA(iddd,10)*PProt_dB(idpp,6);
            block(4 + idpp,30 + iddd) += intn[7]*(DDrot_dB(iddd,2)*PProt_dA(idpp,3) + DDrot_dB(iddd,3)*PProt_dA(idpp,2)) + intn[8]*DDrot_dB(iddd,10)*PProt_dA(idpp,6);
            block(4 + idpp,30 + iddd) += intn[7]*(DDrot(iddd,2)*PProt_dAB(idpp,3) + DDrot(iddd,3)*PProt_dAB(idpp,2)) + intn[8]*DDrot(iddd,10)*PProt_dAB(idpp,6);
            block(4 + idpp,30 + iddd) += intn[9]*(DDrot_dAB(iddd,8)*(PProt(idpp,2) - PProt(idpp,3)) + DDrot_dAB(iddd,9)*PProt(idpp,6));
            block(4 + idpp,30 + iddd) += intn[9]*(DDrot_dA(iddd,8)*(PProt_dB(idpp,2) - PProt_dB(idpp,3)) + DDrot_dA(iddd,9)*PProt_dB(idpp,6));
            block(4 + idpp,30 + iddd) += intn[9]*(DDrot_dB(iddd,8)*(PProt_dA(idpp,2) - PProt_dA(idpp,3)) + DDrot_dB(iddd,9)*PProt_dA(idpp,6));
            block(4 + idpp,30 + iddd) += intn[9]*(DDrot(iddd,8)*(PProt_dAB(idpp,2) - PProt_dAB(idpp,3)) + DDrot(iddd,9)*PProt_dAB(idpp,6));
            block(4 + idpp,30 + iddd) += intn[10]*((DDrot_dAB(iddd,11) + DDrot_dAB(iddd,14))*PProt(idpp,4) + (DDrot_dAB(iddd,12) - DDrot_dAB(iddd,13))*PProt(idpp,5));
            block(4 + idpp,30 + iddd) += intn[10]*((DDrot_dA(iddd,11) + DDrot_dA(iddd,14))*PProt_dB(idpp,4) + (DDrot_dA(iddd,12) - DDrot_dA(iddd,13))*PProt_dB(idpp,5));
            block(4 + idpp,30 + iddd) += intn[10]*((DDrot_dB(iddd,11) + DDrot_dB(iddd,14))*PProt_dA(idpp,4) + (DDrot_dB(iddd,12) - DDrot_dB(iddd,13))*PProt_dA(idpp,5));
            block(4 + idpp,30 + iddd) += intn[10]*((DDrot(iddd,11) + DDrot(iddd,14))*PProt_dAB(idpp,4) + (DDrot(iddd,12) - DDrot(iddd,13))*PProt_dAB(idpp,5));
          }
        }
      }
      if ((nrows > 10)&&(ncols > 10)) {
        //(sd|sd) integrals
        Dd[0] = Dvalue(atmC,4);
        Dd[1] = Dvalue(atmD,4);
        intn[0] = d_eri2Center(0,0,2,0,0,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(0,0,2,1,0,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[2] = d_eri2Center(0,0,2,2,0,0,2,2,RCD,Dd,atmC,atmD,type);    //(sg,dl|sg,dl)=(sg,-dl|sg,-dl)
        //(sd|sd) rotate
        for (size_t idsdb = 1; idsdb < 6; ++idsdb) {
          for (size_t idsdk = 1; idsdk < 6; ++idsdk) {
            block(10 + idsdb,10 + idsdk) = intn[0]*SDrot_dAB(idsdb,1)*SDrot(idsdk,1) + intn[1]*(SDrot_dAB(idsdb,2)*SDrot(idsdk,2) + SDrot_dAB(idsdb,3)*SDrot(idsdk,3)) + intn[2]*(SDrot_dAB(idsdb,4)*SDrot(idsdk,4) + SDrot_dAB(idsdb,5)*SDrot(idsdk,5));
            block(10 + idsdb,10 + idsdk) += intn[0]*SDrot_dA(idsdb,1)*SDrot_dB(idsdk,1) + intn[1]*(SDrot_dA(idsdb,2)*SDrot_dB(idsdk,2) + SDrot_dA(idsdb,3)*SDrot_dB(idsdk,3)) + intn[2]*(SDrot_dA(idsdb,4)*SDrot_dB(idsdk,4) + SDrot_dA(idsdb,5)*SDrot_dB(idsdk,5));
            block(10 + idsdb,10 + idsdk) += intn[0]*SDrot_dB(idsdb,1)*SDrot_dA(idsdk,1) + intn[1]*(SDrot_dB(idsdb,2)*SDrot_dA(idsdk,2) + SDrot_dB(idsdb,3)*SDrot_dA(idsdk,3)) + intn[2]*(SDrot_dB(idsdb,4)*SDrot_dA(idsdk,4) + SDrot_dB(idsdb,5)*SDrot_dA(idsdk,5));
            block(10 + idsdb,10 + idsdk) += intn[0]*SDrot(idsdb,1)*SDrot_dAB(idsdk,1) + intn[1]*(SDrot(idsdb,2)*SDrot_dAB(idsdk,2) + SDrot(idsdb,3)*SDrot_dAB(idsdk,3)) + intn[2]*(SDrot(idsdb,4)*SDrot_dAB(idsdk,4) + SDrot(idsdb,5)*SDrot_dAB(idsdk,5));
          }
        }
        //(sd|pd) integrals
        Dd[1] = Dvalue(atmD,5);
        intn[0] = d_eri2Center(0,0,2,0,1,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(0,0,2,0,1,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center(0,0,2,1,1,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[3] = d_eri2Center(0,0,2,1,1,1,2,2,RCD,Dd,atmC,atmD,type);    //(sg,pi|pi,dl)=-(sg,-pi|-pi,dl)=(sg,-pi|pi,-dl)=(sg,pi|-pi,-dl)
        intn[4] = d_eri2Center(0,0,2,1,1,1,2,0,RCD,Dd,atmC,atmD,type);    //(sg,pi|pi,sg)=(sg,-pi|-pi,sg)
        //(sd|pd) rotate
        for (size_t idsd = 1; idsd < 6; ++idsd) {
          for (size_t idpd1 = 1; idpd1 < 6; ++idpd1) {
            for (size_t idpd2 = 1; idpd2 < 4; ++idpd2) {
              block(10 + idsd,5*(idpd2+2) + idpd1) = -intn[0]*PDrot_dAB((idpd1 - 1)*3 + idpd2,1)*SDrot(idsd,1) - intn[1]*(PDrot_dAB((idpd1 - 1)*3 + idpd2,2) + PDrot_dAB((idpd1 - 1)*3 + idpd2,3))*SDrot(idsd,1);
              block(10 + idsd,5*(idpd2+2) + idpd1) -= intn[0]*PDrot_dA((idpd1 - 1)*3 + idpd2,1)*SDrot_dB(idsd,1) + intn[1]*(PDrot_dA((idpd1 - 1)*3 + idpd2,2) + PDrot_dA((idpd1 - 1)*3 + idpd2,3))*SDrot_dB(idsd,1);
              block(10 + idsd,5*(idpd2+2) + idpd1) -= intn[0]*PDrot_dB((idpd1 - 1)*3 + idpd2,1)*SDrot_dA(idsd,1) + intn[1]*(PDrot_dB((idpd1 - 1)*3 + idpd2,2) + PDrot_dB((idpd1 - 1)*3 + idpd2,3))*SDrot_dA(idsd,1);
              block(10 + idsd,5*(idpd2+2) + idpd1) -= intn[0]*PDrot((idpd1 - 1)*3 + idpd2,1)*SDrot_dAB(idsd,1) + intn[1]*(PDrot((idpd1 - 1)*3 + idpd2,2) + PDrot((idpd1 - 1)*3 + idpd2,3))*SDrot_dAB(idsd,1);
              block(10 + idsd,5*(idpd2+2) + idpd1) -= intn[2]*(PDrot_dAB((idpd1 - 1)*3 + idpd2,4)*SDrot(idsd,2) + PDrot_dAB((idpd1 - 1)*3 + idpd2,5)*SDrot(idsd,3));
              block(10 + idsd,5*(idpd2+2) + idpd1) -= intn[2]*(PDrot_dA((idpd1 - 1)*3 + idpd2,4)*SDrot_dB(idsd,2) + PDrot_dA((idpd1 - 1)*3 + idpd2,5)*SDrot_dB(idsd,3));
              block(10 + idsd,5*(idpd2+2) + idpd1) -= intn[2]*(PDrot_dB((idpd1 - 1)*3 + idpd2,4)*SDrot_dA(idsd,2) + PDrot_dB((idpd1 - 1)*3 + idpd2,5)*SDrot_dA(idsd,3));
              block(10 + idsd,5*(idpd2+2) + idpd1) -= intn[2]*(PDrot((idpd1 - 1)*3 + idpd2,4)*SDrot_dAB(idsd,2) + PDrot((idpd1 - 1)*3 + idpd2,5)*SDrot_dAB(idsd,3));
              block(10 + idsd,5*(idpd2+2) + idpd1) -= intn[3]*((PDrot_dAB((idpd1 - 1)*3 + idpd2,10) + PDrot_dAB((idpd1 - 1)*3 + idpd2,15))*SDrot(idsd,2) + (PDrot_dAB((idpd1 - 1)*3 + idpd2,11) - PDrot_dAB((idpd1 - 1)*3 + idpd2,14))*SDrot(idsd,3));
              block(10 + idsd,5*(idpd2+2) + idpd1) -= intn[3]*((PDrot_dA((idpd1 - 1)*3 + idpd2,10) + PDrot_dA((idpd1 - 1)*3 + idpd2,15))*SDrot_dB(idsd,2) + (PDrot_dA((idpd1 - 1)*3 + idpd2,11) - PDrot_dA((idpd1 - 1)*3 + idpd2,14))*SDrot_dB(idsd,3));
              block(10 + idsd,5*(idpd2+2) + idpd1) -= intn[3]*((PDrot_dB((idpd1 - 1)*3 + idpd2,10) + PDrot_dB((idpd1 - 1)*3 + idpd2,15))*SDrot_dA(idsd,2) + (PDrot_dB((idpd1 - 1)*3 + idpd2,11) - PDrot_dB((idpd1 - 1)*3 + idpd2,14))*SDrot_dA(idsd,3));
              block(10 + idsd,5*(idpd2+2) + idpd1) -= intn[3]*((PDrot((idpd1 - 1)*3 + idpd2,10) + PDrot((idpd1 - 1)*3 + idpd2,15))*SDrot_dAB(idsd,2) + (PDrot((idpd1 - 1)*3 + idpd2,11) - PDrot((idpd1 - 1)*3 + idpd2,14))*SDrot_dAB(idsd,3));
              block(10 + idsd,5*(idpd2+2) + idpd1) -= intn[4]*(PDrot_dAB((idpd1 - 1)*3 + idpd2,8)*SDrot(idsd,2) + PDrot_dAB((idpd1 - 1)*3 + idpd2,12)*SDrot(idsd,3));
              block(10 + idsd,5*(idpd2+2) + idpd1) -= intn[4]*(PDrot_dA((idpd1 - 1)*3 + idpd2,8)*SDrot_dB(idsd,2) + PDrot_dA((idpd1 - 1)*3 + idpd2,12)*SDrot_dB(idsd,3));
              block(10 + idsd,5*(idpd2+2) + idpd1) -= intn[4]*(PDrot_dB((idpd1 - 1)*3 + idpd2,8)*SDrot_dA(idsd,2) + PDrot_dB((idpd1 - 1)*3 + idpd2,12)*SDrot_dA(idsd,3));
              block(10 + idsd,5*(idpd2+2) + idpd1) -= intn[4]*(PDrot((idpd1 - 1)*3 + idpd2,8)*SDrot_dAB(idsd,2) + PDrot((idpd1 - 1)*3 + idpd2,12)*SDrot_dAB(idsd,3));
            }
          }
        }
        //(pd|sd) integrals
        Dd[0] = Dvalue(atmC,5);
        Dd[1] = Dvalue(atmD,4);
        intn[0] = d_eri2Center(1,0,2,0,0,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(1,1,2,1,0,0,2,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[2] = d_eri2Center(1,0,2,1,0,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[3] = d_eri2Center(1,1,2,2,0,0,2,1,RCD,Dd,atmC,atmD,type);    //(pi,dl|sg,pi)=-(-pi,dl|sg,-pi)=(pi,-dl|sg,-pi)=(-pi,-dl|sg,pi)
        intn[4] = d_eri2Center(1,1,2,0,0,0,2,1,RCD,Dd,atmC,atmD,type);    //(pi,sg|sg,pi)=(-pi,sg|sg,-pi)
        //(pd|sd) rotate
        for (size_t idpd1 = 1; idpd1 < 6; ++idpd1) {
          for (size_t idpd2 = 1; idpd2 < 4; ++idpd2) {
            for (size_t idsd = 1; idsd < 6; ++idsd) {
              block(5*(idpd2+2) + idpd1,10 + idsd) = -intn[0]*PDrot_dAB((idpd1 - 1)*3 + idpd2,1)*SDrot(idsd,1) - intn[1]*(PDrot_dAB((idpd1 - 1)*3 + idpd2,2) + PDrot_dAB((idpd1 - 1)*3 + idpd2,3))*SDrot(idsd,1);
              block(5*(idpd2+2) + idpd1,10 + idsd) -= intn[0]*PDrot_dA((idpd1 - 1)*3 + idpd2,1)*SDrot_dB(idsd,1) + intn[1]*(PDrot_dA((idpd1 - 1)*3 + idpd2,2) + PDrot_dA((idpd1 - 1)*3 + idpd2,3))*SDrot_dB(idsd,1);
              block(5*(idpd2+2) + idpd1,10 + idsd) -= intn[0]*PDrot_dB((idpd1 - 1)*3 + idpd2,1)*SDrot_dA(idsd,1) + intn[1]*(PDrot_dB((idpd1 - 1)*3 + idpd2,2) + PDrot_dB((idpd1 - 1)*3 + idpd2,3))*SDrot_dA(idsd,1);
              block(5*(idpd2+2) + idpd1,10 + idsd) -= intn[0]*PDrot((idpd1 - 1)*3 + idpd2,1)*SDrot_dAB(idsd,1) + intn[1]*(PDrot((idpd1 - 1)*3 + idpd2,2) + PDrot((idpd1 - 1)*3 + idpd2,3))*SDrot_dAB(idsd,1);
              block(5*(idpd2+2) + idpd1,10 + idsd) -= intn[2]*(PDrot_dAB((idpd1 - 1)*3 + idpd2,4)*SDrot(idsd,2) + PDrot_dAB((idpd1 - 1)*3 + idpd2,5)*SDrot(idsd,3));
              block(5*(idpd2+2) + idpd1,10 + idsd) -= intn[2]*(PDrot_dA((idpd1 - 1)*3 + idpd2,4)*SDrot_dB(idsd,2) + PDrot_dA((idpd1 - 1)*3 + idpd2,5)*SDrot_dB(idsd,3));
              block(5*(idpd2+2) + idpd1,10 + idsd) -= intn[2]*(PDrot_dB((idpd1 - 1)*3 + idpd2,4)*SDrot_dA(idsd,2) + PDrot_dB((idpd1 - 1)*3 + idpd2,5)*SDrot_dA(idsd,3));
              block(5*(idpd2+2) + idpd1,10 + idsd) -= intn[2]*(PDrot((idpd1 - 1)*3 + idpd2,4)*SDrot_dAB(idsd,2) + PDrot((idpd1 - 1)*3 + idpd2,5)*SDrot_dAB(idsd,3));
              block(5*(idpd2+2) + idpd1,10 + idsd) -= intn[3]*((PDrot_dAB((idpd1 - 1)*3 + idpd2,10) + PDrot_dAB((idpd1 - 1)*3 + idpd2,15))*SDrot(idsd,2) + (PDrot_dAB((idpd1 - 1)*3 + idpd2,11) - PDrot_dAB((idpd1 - 1)*3 + idpd2,14))*SDrot(idsd,3));
              block(5*(idpd2+2) + idpd1,10 + idsd) -= intn[3]*((PDrot_dA((idpd1 - 1)*3 + idpd2,10) + PDrot_dA((idpd1 - 1)*3 + idpd2,15))*SDrot_dB(idsd,2) + (PDrot_dA((idpd1 - 1)*3 + idpd2,11) - PDrot_dA((idpd1 - 1)*3 + idpd2,14))*SDrot_dB(idsd,3));
              block(5*(idpd2+2) + idpd1,10 + idsd) -= intn[3]*((PDrot_dB((idpd1 - 1)*3 + idpd2,10) + PDrot_dB((idpd1 - 1)*3 + idpd2,15))*SDrot_dA(idsd,2) + (PDrot_dB((idpd1 - 1)*3 + idpd2,11) - PDrot_dB((idpd1 - 1)*3 + idpd2,14))*SDrot_dA(idsd,3));
              block(5*(idpd2+2) + idpd1,10 + idsd) -= intn[3]*((PDrot((idpd1 - 1)*3 + idpd2,10) + PDrot((idpd1 - 1)*3 + idpd2,15))*SDrot_dAB(idsd,2) + (PDrot((idpd1 - 1)*3 + idpd2,11) - PDrot((idpd1 - 1)*3 + idpd2,14))*SDrot_dAB(idsd,3));
              block(5*(idpd2+2) + idpd1,10 + idsd) -= intn[4]*(PDrot_dAB((idpd1 - 1)*3 + idpd2,8)*SDrot(idsd,2) + PDrot_dAB((idpd1 - 1)*3 + idpd2,12)*SDrot(idsd,3));
              block(5*(idpd2+2) + idpd1,10 + idsd) -= intn[4]*(PDrot_dA((idpd1 - 1)*3 + idpd2,8)*SDrot_dB(idsd,2) + PDrot_dA((idpd1 - 1)*3 + idpd2,12)*SDrot_dB(idsd,3));
              block(5*(idpd2+2) + idpd1,10 + idsd) -= intn[4]*(PDrot_dB((idpd1 - 1)*3 + idpd2,8)*SDrot_dA(idsd,2) + PDrot_dB((idpd1 - 1)*3 + idpd2,12)*SDrot_dA(idsd,3));
              block(5*(idpd2+2) + idpd1,10 + idsd) -= intn[4]*(PDrot((idpd1 - 1)*3 + idpd2,8)*SDrot_dAB(idsd,2) + PDrot((idpd1 - 1)*3 + idpd2,12)*SDrot_dAB(idsd,3));
            }
          }
        }
        //(dd|sd) integrals
        Dd[0] = Dvalue(atmC,6);
        Dd[1] = Dvalue(atmD,4);
        intn[0] = d_eri2Center(2,0,2,0,0,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(2,1,2,1,0,0,2,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[2] = d_eri2Center(2,2,2,2,0,0,2,0,RCD,Dd,atmC,atmD,type);    //(dl,dl|sg,sg)=(-dl,-dl|sg,sg)
        intn[3] = d_eri2Center(2,0,2,1,0,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[4] = d_eri2Center(2,0,2,2,0,0,2,2,RCD,Dd,atmC,atmD,type);    //(sg,dl|sg,dl)=(sg,-dl|sg,-dl)
        intn[5] = d_eri2Center(2,1,2,1,0,0,2,2,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,dl)=-(-pi,-pi|sg,dl)=(pi,-pi|sg,-dl)
        intn[6] = d_eri2Center(2,1,2,2,0,0,2,1,RCD,Dd,atmC,atmD,type);    //(pi,dl|sg,pi)=(-pi,-dl|sg,pi)=(pi,-dl|sg,-pi)=-(-pi,dl|sg,-pi)
        //(dd|sd) rotate
        for (size_t iddd = 1; iddd < 16; ++iddd) {
          for (size_t idsd = 1; idsd < 6; ++idsd) {
            block(30 + iddd,10 + idsd) = intn[0]*DDrot_dAB(iddd,1)*SDrot(idsd,1) + intn[1]*(DDrot_dAB(iddd,2) + DDrot_dAB(iddd,3))*SDrot(idsd,1) + intn[2]*(DDrot_dB(iddd,4) + DDrot_dB(iddd,5))*SDrot_dA(idsd,1);
            block(30 + iddd,10 + idsd) += intn[0]*DDrot_dA(iddd,1)*SDrot_dB(idsd,1) + intn[1]*(DDrot_dA(iddd,2) + DDrot_dA(iddd,3))*SDrot_dB(idsd,1) + intn[2]*(DDrot(iddd,4) + DDrot(iddd,5))*SDrot_dAB(idsd,1);
            block(30 + iddd,10 + idsd) += intn[0]*DDrot_dB(iddd,1)*SDrot_dA(idsd,1) + intn[1]*(DDrot_dB(iddd,2) + DDrot_dB(iddd,3))*SDrot_dA(idsd,1) + intn[2]*(DDrot_dAB(iddd,4) + DDrot_dAB(iddd,5))*SDrot(idsd,1);
            block(30 + iddd,10 + idsd) += intn[0]*DDrot(iddd,1)*SDrot_dAB(idsd,1) + intn[1]*(DDrot(iddd,2) + DDrot(iddd,3))*SDrot_dAB(idsd,1) + intn[2]*(DDrot_dA(iddd,4) + DDrot_dA(iddd,5))*SDrot_dB(idsd,1);
            block(30 + iddd,10 + idsd) += intn[3]*(DDrot_dAB(iddd,6)*SDrot(idsd,2) + DDrot_dAB(iddd,7)*SDrot(idsd,3)) + intn[4]*(DDrot_dAB(iddd,8)*SDrot(idsd,4) + DDrot_dAB(iddd,9)*SDrot(idsd,5));
            block(30 + iddd,10 + idsd) += intn[3]*(DDrot_dA(iddd,6)*SDrot_dB(idsd,2) + DDrot_dA(iddd,7)*SDrot_dB(idsd,3)) + intn[4]*(DDrot_dA(iddd,8)*SDrot_dB(idsd,4) + DDrot_dA(iddd,9)*SDrot_dB(idsd,5));
            block(30 + iddd,10 + idsd) += intn[3]*(DDrot_dB(iddd,6)*SDrot_dA(idsd,2) + DDrot_dB(iddd,7)*SDrot_dA(idsd,3)) + intn[4]*(DDrot_dB(iddd,8)*SDrot_dA(idsd,4) + DDrot_dB(iddd,9)*SDrot_dA(idsd,5));
            block(30 + iddd,10 + idsd) += intn[3]*(DDrot(iddd,6)*SDrot_dAB(idsd,2) + DDrot(iddd,7)*SDrot_dAB(idsd,3)) + intn[4]*(DDrot(iddd,8)*SDrot_dAB(idsd,4) + DDrot(iddd,9)*SDrot_dAB(idsd,5));
            block(30 + iddd,10 + idsd) += intn[5]*((DDrot_dAB(iddd,2) - DDrot_dAB(iddd,3))*SDrot(idsd,4) + DDrot_dAB(iddd,10)*SDrot(idsd,5));
            block(30 + iddd,10 + idsd) += intn[5]*((DDrot_dA(iddd,2) - DDrot_dA(iddd,3))*SDrot_dB(idsd,4) + DDrot_dA(iddd,10)*SDrot_dB(idsd,5));
            block(30 + iddd,10 + idsd) += intn[5]*((DDrot_dB(iddd,2) - DDrot_dB(iddd,3))*SDrot_dA(idsd,4) + DDrot_dB(iddd,10)*SDrot_dA(idsd,5));
            block(30 + iddd,10 + idsd) += intn[5]*((DDrot(iddd,2) - DDrot(iddd,3))*SDrot_dAB(idsd,4) + DDrot(iddd,10)*SDrot_dAB(idsd,5));
            block(30 + iddd,10 + idsd) += intn[6]*((DDrot_dAB(iddd,11) + DDrot_dAB(iddd,14))*SDrot(idsd,2) + (DDrot_dB(iddd,12) - DDrot_dB(iddd,13))*SDrot_dA(idsd,3));
            block(30 + iddd,10 + idsd) += intn[6]*((DDrot_dA(iddd,11) + DDrot_dA(iddd,14))*SDrot_dB(idsd,2) + (DDrot(iddd,12) - DDrot(iddd,13))*SDrot_dAB(idsd,3));
            block(30 + iddd,10 + idsd) += intn[6]*((DDrot_dB(iddd,11) + DDrot_dB(iddd,14))*SDrot_dA(idsd,2) + (DDrot_dAB(iddd,12) - DDrot_dAB(iddd,13))*SDrot(idsd,3));
            block(30 + iddd,10 + idsd) += intn[6]*((DDrot(iddd,11) + DDrot(iddd,14))*SDrot_dAB(idsd,2) + (DDrot_dA(iddd,12) - DDrot_dA(iddd,13))*SDrot_dB(idsd,3));
          }
        }
        //(sd|dd) integrals
        Dd[0] = Dvalue(atmC,4);
        Dd[1] = Dvalue(atmD,6);
        intn[0] = d_eri2Center(0,0,2,0,2,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(0,0,2,0,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center(0,0,2,0,2,2,2,2,RCD,Dd,atmC,atmD,type);    //(sg,sg|dl,dl)=(sg,sg|-dl,-dl)
        intn[3] = d_eri2Center(0,0,2,1,2,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[4] = d_eri2Center(0,0,2,2,2,0,2,2,RCD,Dd,atmC,atmD,type);    //(sg,dl|sg,dl)=(sg,-dl|sg,-dl)
        intn[5] = d_eri2Center(0,0,2,2,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,dl|pi,pi)=-(sg,dl|-pi,-pi)=(sg,-dl|pi,-pi)
        intn[6] = d_eri2Center(0,0,2,1,2,1,2,2,RCD,Dd,atmC,atmD,type);    //(sg,pi|pi,dl)=(sg,pi|-pi,-dl)=(sg,-pi|pi,-dl)=-(sg,-pi|-pi,dl)
        //(sd|dd) rotate
        for (size_t idsd = 1; idsd < 6; ++idsd) {
          for (size_t iddd = 1; iddd < 16; ++iddd) {
            block(10 + idsd,30 + iddd) = intn[0]*DDrot_dAB(iddd,1)*SDrot(idsd,1) + intn[1]*(DDrot_dAB(iddd,2) + DDrot_dAB(iddd,3))*SDrot(idsd,1) + intn[2]*(DDrot_dB(iddd,4) + DDrot_dB(iddd,5))*SDrot_dA(idsd,1);
            block(10 + idsd,30 + iddd) += intn[0]*DDrot_dA(iddd,1)*SDrot_dB(idsd,1) + intn[1]*(DDrot_dA(iddd,2) + DDrot_dA(iddd,3))*SDrot_dB(idsd,1) + intn[2]*(DDrot(iddd,4) + DDrot(iddd,5))*SDrot_dAB(idsd,1);
            block(10 + idsd,30 + iddd) += intn[0]*DDrot_dB(iddd,1)*SDrot_dA(idsd,1) + intn[1]*(DDrot_dB(iddd,2) + DDrot_dB(iddd,3))*SDrot_dA(idsd,1) + intn[2]*(DDrot_dAB(iddd,4) + DDrot_dAB(iddd,5))*SDrot(idsd,1);
            block(10 + idsd,30 + iddd) += intn[0]*DDrot(iddd,1)*SDrot_dAB(idsd,1) + intn[1]*(DDrot(iddd,2) + DDrot(iddd,3))*SDrot_dAB(idsd,1) + intn[2]*(DDrot_dA(iddd,4) + DDrot_dA(iddd,5))*SDrot_dB(idsd,1);
            block(10 + idsd,30 + iddd) += intn[3]*(DDrot_dAB(iddd,6)*SDrot(idsd,2) + DDrot_dAB(iddd,7)*SDrot(idsd,3)) + intn[4]*(DDrot_dAB(iddd,8)*SDrot(idsd,4) + DDrot_dAB(iddd,9)*SDrot(idsd,5));
            block(10 + idsd,30 + iddd) += intn[3]*(DDrot_dA(iddd,6)*SDrot_dB(idsd,2) + DDrot_dA(iddd,7)*SDrot_dB(idsd,3)) + intn[4]*(DDrot_dA(iddd,8)*SDrot_dB(idsd,4) + DDrot_dA(iddd,9)*SDrot_dB(idsd,5));
            block(10 + idsd,30 + iddd) += intn[3]*(DDrot_dB(iddd,6)*SDrot_dA(idsd,2) + DDrot_dB(iddd,7)*SDrot_dA(idsd,3)) + intn[4]*(DDrot_dB(iddd,8)*SDrot_dA(idsd,4) + DDrot_dB(iddd,9)*SDrot_dA(idsd,5));
            block(10 + idsd,30 + iddd) += intn[3]*(DDrot(iddd,6)*SDrot_dAB(idsd,2) + DDrot(iddd,7)*SDrot_dAB(idsd,3)) + intn[4]*(DDrot(iddd,8)*SDrot_dAB(idsd,4) + DDrot(iddd,9)*SDrot_dAB(idsd,5));
            block(10 + idsd,30 + iddd) += intn[5]*((DDrot_dAB(iddd,2) - DDrot_dAB(iddd,3))*SDrot(idsd,4) + DDrot_dAB(iddd,10)*SDrot(idsd,5));
            block(10 + idsd,30 + iddd) += intn[5]*((DDrot_dA(iddd,2) - DDrot_dA(iddd,3))*SDrot_dB(idsd,4) + DDrot_dA(iddd,10)*SDrot_dB(idsd,5));
            block(10 + idsd,30 + iddd) += intn[5]*((DDrot_dB(iddd,2) - DDrot_dB(iddd,3))*SDrot_dA(idsd,4) + DDrot_dB(iddd,10)*SDrot_dA(idsd,5));
            block(10 + idsd,30 + iddd) += intn[5]*((DDrot(iddd,2) - DDrot(iddd,3))*SDrot_dAB(idsd,4) + DDrot(iddd,10)*SDrot_dAB(idsd,5));
            block(10 + idsd,30 + iddd) += intn[6]*((DDrot_dAB(iddd,11) + DDrot_dAB(iddd,14))*SDrot(idsd,2) + (DDrot_dB(iddd,12) - DDrot_dB(iddd,13))*SDrot_dA(idsd,3));
            block(10 + idsd,30 + iddd) += intn[6]*((DDrot_dA(iddd,11) + DDrot_dA(iddd,14))*SDrot_dB(idsd,2) + (DDrot(iddd,12) - DDrot(iddd,13))*SDrot_dAB(idsd,3));
            block(10 + idsd,30 + iddd) += intn[6]*((DDrot_dB(iddd,11) + DDrot_dB(iddd,14))*SDrot_dA(idsd,2) + (DDrot_dAB(iddd,12) - DDrot_dAB(iddd,13))*SDrot(idsd,3));
            block(10 + idsd,30 + iddd) += intn[6]*((DDrot(iddd,11) + DDrot(iddd,14))*SDrot_dAB(idsd,2) + (DDrot_dA(iddd,12) - DDrot_dA(iddd,13))*SDrot_dB(idsd,3));
          }
        }
        //(pd|pd) integrals
        Dd[0] = Dvalue(atmC,5);
        Dd[1] = Dvalue(atmD,5);
        intn[0] = d_eri2Center(1,0,2,0,1,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(1,1,2,1,1,0,2,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)=(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center(1,1,2,1,1,1,2,1,RCD,Dd,atmC,atmD,type);    //(pi,pi|pi,pi)=(-pi,-pi|pi,pi)=(pi,pi|-pi,-pi)=(-pi,-pi|-pi,-pi)
        intn[3] = d_eri2Center(1,0,2,1,1,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[4] = d_eri2Center(1,1,2,0,1,1,2,0,RCD,Dd,atmC,atmD,type);    //(pi,sg|pi,sg)=(-pi,sg|-pi,sg)
        intn[5] = d_eri2Center(1,1,2,2,1,1,2,2,RCD,Dd,atmC,atmD,type);    //(pi,dl|pi,dl)=(pi,dl|-pi,-dl)=(-pi,dl|-pi,dl)=-(-pi,dl|pi,-dl)=-(pi,-dl|-pi,dl)=(pi,-dl|pi,-dl)=(-pi,-dl|pi,dl)=(-pi,-dl|-pi,-dl)
        intn[6] = d_eri2Center(1,1,2,0,1,0,2,1,RCD,Dd,atmC,atmD,type);    //(pi,sg|sg,pi)=(-pi,sg|sg,-pi)=(sg,pi|pi,sg)=(sg,-pi|-pi,sg)
        intn[7] = d_eri2Center(1,1,2,2,1,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|pi,dl)=(sg,pi|-pi,-dl)=-(sg,-pi|-pi,dl)=(sg,-pi|pi,-dl)=(pi,dl|sg,pi)=-(-pi,dl|sg,-pi)=(pi,-dl|sg,-pi)=(-pi,-dl|sg,pi)
        intn[8] = d_eri2Center(1,1,2,2,1,1,2,0,RCD,Dd,atmC,atmD,type);    //(pi,sg|pi,dl)=(pi,sg|-pi,-dl)=-(-pi,sg|-pi,dl)=(-pi,sg|pi,-dl)=(pi,dl|pi,sg)=-(-pi,dl|-pi,sg)=(pi,-dl|-pi,sg)=(-pi,-dl|pi,sg)
        //(pd|pd) rotate
        for (size_t idpdb1 = 1; idpdb1 < 6; ++idpdb1) {            //d orbital on bra
          for (size_t idpdb2 = 1; idpdb2 < 4; ++idpdb2) {          //p orbital on bra
            for (size_t idpdk1 = 1; idpdk1 < 6; ++idpdk1) {        //d orbital on ket
              for (size_t idpdk2 = 1; idpdk2 < 4; ++idpdk2) {      //p orbital on ket
                cnt1 = (idpdb1 - 1)*3 + idpdb2;
                cnt2 = (idpdk1 - 1)*3 + idpdk2;
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) = intn[0]*PDrot_dAB(cnt1,1)*PDrot(cnt2,1) + intn[1]*((PDrot_dAB(cnt1,2) + PDrot_dAB(cnt1,3))*PDrot(cnt2,1) + PDrot_dAB(cnt1,1)*(PDrot(cnt2,2) + PDrot(cnt2,3)));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[0]*PDrot_dA(cnt1,1)*PDrot_dB(cnt2,1) + intn[1]*((PDrot_dA(cnt1,2) + PDrot_dA(cnt1,3))*PDrot_dB(cnt2,1) + PDrot_dA(cnt1,1)*(PDrot_dB(cnt2,2) + PDrot_dB(cnt2,3)));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[0]*PDrot_dB(cnt1,1)*PDrot_dA(cnt2,1) + intn[1]*((PDrot_dB(cnt1,2) + PDrot_dB(cnt1,3))*PDrot_dA(cnt2,1) + PDrot_dB(cnt1,1)*(PDrot_dA(cnt2,2) + PDrot_dA(cnt2,3)));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[0]*PDrot(cnt1,1)*PDrot_dAB(cnt2,1) + intn[1]*((PDrot(cnt1,2) + PDrot(cnt1,3))*PDrot_dAB(cnt2,1) + PDrot(cnt1,1)*(PDrot_dAB(cnt2,2) + PDrot_dAB(cnt2,3)));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[2]*(PDrot_dAB(cnt1,2) + PDrot_dAB(cnt1,3))*(PDrot(cnt2,2) + PDrot(cnt2,3)) + intn[3]*(PDrot_dAB(cnt1,4)*PDrot(cnt2,4) + PDrot_dAB(cnt1,5)*PDrot(cnt2,5));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[2]*(PDrot_dA(cnt1,2) + PDrot_dA(cnt1,3))*(PDrot_dB(cnt2,2) + PDrot_dB(cnt2,3)) + intn[3]*(PDrot_dA(cnt1,4)*PDrot_dB(cnt2,4) + PDrot_dA(cnt1,5)*PDrot_dB(cnt2,5));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[2]*(PDrot_dB(cnt1,2) + PDrot_dB(cnt1,3))*(PDrot_dA(cnt2,2) + PDrot_dA(cnt2,3)) + intn[3]*(PDrot_dB(cnt1,4)*PDrot_dA(cnt2,4) + PDrot_dB(cnt1,5)*PDrot_dA(cnt2,5));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[2]*(PDrot(cnt1,2) + PDrot(cnt1,3))*(PDrot_dAB(cnt2,2) + PDrot_dAB(cnt2,3)) + intn[3]*(PDrot(cnt1,4)*PDrot_dAB(cnt2,4) + PDrot(cnt1,5)*PDrot_dAB(cnt2,5));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[4]*(PDrot_dAB(cnt1,8)*PDrot(cnt2,8) + PDrot_dAB(cnt1,12)*PDrot(cnt2,12));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[4]*(PDrot_dA(cnt1,8)*PDrot_dB(cnt2,8) + PDrot_dA(cnt1,12)*PDrot_dB(cnt2,12));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[4]*(PDrot_dB(cnt1,8)*PDrot_dA(cnt2,8) + PDrot_dB(cnt1,12)*PDrot_dA(cnt2,12));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[4]*(PDrot(cnt1,8)*PDrot_dAB(cnt2,8) + PDrot(cnt1,12)*PDrot_dAB(cnt2,12));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[5]*((PDrot_dAB(cnt1,10) + PDrot_dAB(cnt1,15))*(PDrot(cnt2,10) + PDrot(cnt2,15)) + (PDrot_dAB(cnt1,14) - PDrot_dAB(cnt1,11))*(PDrot(cnt2,14) - PDrot(cnt2,11)));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[5]*((PDrot_dA(cnt1,10) + PDrot_dA(cnt1,15))*(PDrot_dB(cnt2,10) + PDrot_dB(cnt2,15)) + (PDrot_dA(cnt1,14) - PDrot_dA(cnt1,11))*(PDrot_dB(cnt2,14) - PDrot_dB(cnt2,11)));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[5]*((PDrot_dB(cnt1,10) + PDrot_dB(cnt1,15))*(PDrot_dA(cnt2,10) + PDrot_dA(cnt2,15)) + (PDrot_dB(cnt1,14) - PDrot_dB(cnt1,11))*(PDrot_dA(cnt2,14) - PDrot_dA(cnt2,11)));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[5]*((PDrot(cnt1,10) + PDrot(cnt1,15))*(PDrot_dAB(cnt2,10) + PDrot_dAB(cnt2,15)) + (PDrot(cnt1,14) - PDrot(cnt1,11))*(PDrot_dAB(cnt2,14) - PDrot_dAB(cnt2,11)));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[6]*(PDrot_dAB(cnt1,4)*PDrot(cnt2,8) + PDrot_dAB(cnt1,8)*PDrot(cnt2,4) + PDrot_dAB(cnt1,5)*PDrot(cnt2,12) + PDrot_dAB(cnt1,12)*PDrot(cnt2,5));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[6]*(PDrot_dA(cnt1,4)*PDrot_dB(cnt2,8) + PDrot_dA(cnt1,8)*PDrot_dB(cnt2,4) + PDrot_dA(cnt1,5)*PDrot_dB(cnt2,12) + PDrot_dA(cnt1,12)*PDrot_dB(cnt2,5));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[6]*(PDrot_dB(cnt1,4)*PDrot_dA(cnt2,8) + PDrot_dB(cnt1,8)*PDrot_dA(cnt2,4) + PDrot_dB(cnt1,5)*PDrot_dA(cnt2,12) + PDrot_dB(cnt1,12)*PDrot_dA(cnt2,5));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[6]*(PDrot(cnt1,4)*PDrot_dAB(cnt2,8) + PDrot(cnt1,8)*PDrot_dAB(cnt2,4) + PDrot(cnt1,5)*PDrot_dAB(cnt2,12) + PDrot(cnt1,12)*PDrot_dAB(cnt2,5));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[7]*((PDrot_dAB(cnt1,10) + PDrot_dAB(cnt1,15))*PDrot(cnt2,4) + PDrot_dAB(cnt1,4)*(PDrot(cnt2,10) + PDrot(cnt2,15)) + (PDrot_dAB(cnt1,11) - PDrot_dAB(cnt1,14))*PDrot(cnt2,5) + PDrot_dAB(cnt1,5)*(PDrot(cnt2,11) - PDrot(cnt2,14)));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[7]*((PDrot_dA(cnt1,10) + PDrot_dA(cnt1,15))*PDrot_dB(cnt2,4) + PDrot_dA(cnt1,4)*(PDrot_dB(cnt2,10) + PDrot_dB(cnt2,15)) + (PDrot_dA(cnt1,11) - PDrot_dA(cnt1,14))*PDrot_dB(cnt2,5) + PDrot_dA(cnt1,5)*(PDrot_dB(cnt2,11) - PDrot_dB(cnt2,14)));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[7]*((PDrot_dB(cnt1,10) + PDrot_dB(cnt1,15))*PDrot_dA(cnt2,4) + PDrot_dB(cnt1,4)*(PDrot_dA(cnt2,10) + PDrot_dA(cnt2,15)) + (PDrot_dB(cnt1,11) - PDrot_dB(cnt1,14))*PDrot_dA(cnt2,5) + PDrot_dB(cnt1,5)*(PDrot_dA(cnt2,11) - PDrot_dA(cnt2,14)));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[7]*((PDrot(cnt1,10) + PDrot(cnt1,15))*PDrot_dAB(cnt2,4) + PDrot(cnt1,4)*(PDrot_dAB(cnt2,10) + PDrot_dAB(cnt2,15)) + (PDrot(cnt1,11) - PDrot(cnt1,14))*PDrot_dAB(cnt2,5) + PDrot(cnt1,5)*(PDrot_dAB(cnt2,11) - PDrot_dAB(cnt2,14)));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[8]*((PDrot_dAB(cnt1,10) + PDrot_dAB(cnt1,15))*PDrot(cnt2,8) + PDrot_dAB(cnt1,8)*(PDrot(cnt2,10) + PDrot(cnt2,15)) + (PDrot_dAB(cnt1,11) - PDrot_dAB(cnt1,14))*PDrot(cnt2,12) + PDrot_dAB(cnt1,12)*(PDrot(cnt2,11) - PDrot(cnt2,14)));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[8]*((PDrot_dA(cnt1,10) + PDrot_dA(cnt1,15))*PDrot_dB(cnt2,8) + PDrot_dA(cnt1,8)*(PDrot_dB(cnt2,10) + PDrot_dB(cnt2,15)) + (PDrot_dA(cnt1,11) - PDrot_dA(cnt1,14))*PDrot_dB(cnt2,12) + PDrot_dA(cnt1,12)*(PDrot_dB(cnt2,11) - PDrot_dB(cnt2,14)));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[8]*((PDrot_dB(cnt1,10) + PDrot_dB(cnt1,15))*PDrot_dA(cnt2,8) + PDrot_dB(cnt1,8)*(PDrot_dA(cnt2,10) + PDrot_dA(cnt2,15)) + (PDrot_dB(cnt1,11) - PDrot_dB(cnt1,14))*PDrot_dA(cnt2,12) + PDrot_dB(cnt1,12)*(PDrot_dA(cnt2,11) - PDrot_dA(cnt2,14)));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[8]*((PDrot(cnt1,10) + PDrot(cnt1,15))*PDrot_dAB(cnt2,8) + PDrot(cnt1,8)*(PDrot_dAB(cnt2,10) + PDrot_dAB(cnt2,15)) + (PDrot(cnt1,11) - PDrot(cnt1,14))*PDrot_dAB(cnt2,12) + PDrot(cnt1,12)*(PDrot_dAB(cnt2,11) - PDrot_dAB(cnt2,14)));
              }
            }
          }
        }
        //(pd|dd) integrals
        Dd[0] = Dvalue(atmC,5);
        Dd[1] = Dvalue(atmD,6);
        intn[0] = d_eri2Center(1,0,2,0,2,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(1,0,2,0,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center(1,1,2,1,2,0,2,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[3] = d_eri2Center(1,1,2,1,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(pi,pi|pi,pi)=(-pi,-pi|-pi,-pi)=(pi,pi|-pi,-pi)=(-pi,-pi|pi,pi)
        intn[4] = d_eri2Center(1,0,2,0,2,2,2,2,RCD,Dd,atmC,atmD,type);    //(sg,sg|dl,dl)=(sg,sg|-dl,-dl)
        intn[5] = d_eri2Center(1,1,2,1,2,2,2,2,RCD,Dd,atmC,atmD,type);    //(pi,pi|dl,dl)=(-pi,-pi|-dl,-dl)=(pi,pi|-dl,-dl)=(-pi,-pi|dl,dl)
        intn[6] = d_eri2Center(1,0,2,1,2,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[7] = d_eri2Center(1,1,2,0,2,1,2,0,RCD,Dd,atmC,atmD,type);    //(pi,sg|pi,sg)=(-pi,sg|-pi,sg)
        intn[8] = d_eri2Center(1,1,2,2,2,1,2,2,RCD,Dd,atmC,atmD,type);    //(pi,dl|pi,dl)=(pi,dl|-pi,-dl)=(-pi,-dl|-pi,-dl)=(-pi,-dl|pi,dl)=(-pi,dl|-pi,dl)=-(-pi,dl|pi,-dl)=(pi,-dl|pi,-dl)=-(pi,-dl|-pi,dl)
        intn[9] = d_eri2Center(1,1,2,0,2,1,2,2,RCD,Dd,atmC,atmD,type);    //(pi,sg|pi,dl)=(pi,sg|-pi,-dl)=-(-pi,sg|-pi,dl)=(-pi,sg|pi,-dl)
        intn[10] = d_eri2Center(1,0,2,1,2,1,2,2,RCD,Dd,atmC,atmD,type);   //(sg,pi|pi,dl)=(sg,pi|-pi,-dl)=-(sg,-pi|-pi,dl)=(sg,-pi|pi,-dl)
        intn[11] = d_eri2Center(1,1,2,2,2,0,2,1,RCD,Dd,atmC,atmD,type);   //(pi,dl|sg,pi)=(-pi,-dl|sg,pi)=-(-pi,dl|sg,-pi)=(pi,-dl|sg,-pi)
        //(pd|dd) rotate
        for (size_t iddd = 1; iddd < 16; ++iddd) {
          for (size_t idpd1 = 1; idpd1 < 6; ++idpd1) {
            for (size_t idpd2 = 1; idpd2 < 4; ++idpd2) {
              cnt1 = (idpd1 - 1)*3 + idpd2;
              block(5*(idpd2+2) + idpd1,30 + iddd) = -intn[0]*PDrot_dAB(cnt1,1)*DDrot(iddd,1) - intn[1]*PDrot_dAB(cnt1,1)*(DDrot(iddd,2) + DDrot(iddd,3));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[0]*PDrot_dA(cnt1,1)*DDrot_dB(iddd,1) + intn[1]*PDrot_dA(cnt1,1)*(DDrot_dB(iddd,2) + DDrot_dB(iddd,3));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[0]*PDrot_dB(cnt1,1)*DDrot_dA(iddd,1) + intn[1]*PDrot_dB(cnt1,1)*(DDrot_dA(iddd,2) + DDrot_dA(iddd,3));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[0]*PDrot(cnt1,1)*DDrot_dAB(iddd,1) + intn[1]*PDrot(cnt1,1)*(DDrot_dAB(iddd,2) + DDrot_dAB(iddd,3));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[2]*(PDrot_dAB(cnt1,2) + PDrot_dAB(cnt1,3))*DDrot(iddd,1) + intn[3]*(PDrot_dAB(cnt1,2) + PDrot_dAB(cnt1,3))*(DDrot(iddd,2) + DDrot(iddd,3));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[2]*(PDrot_dA(cnt1,2) + PDrot_dA(cnt1,3))*DDrot_dB(iddd,1) + intn[3]*(PDrot_dA(cnt1,2) + PDrot_dA(cnt1,3))*(DDrot_dB(iddd,2) + DDrot_dB(iddd,3));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[2]*(PDrot_dB(cnt1,2) + PDrot_dB(cnt1,3))*DDrot_dA(iddd,1) + intn[3]*(PDrot_dB(cnt1,2) + PDrot_dB(cnt1,3))*(DDrot_dA(iddd,2) + DDrot_dA(iddd,3));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[2]*(PDrot(cnt1,2) + PDrot(cnt1,3))*DDrot_dAB(iddd,1) + intn[3]*(PDrot(cnt1,2) + PDrot(cnt1,3))*(DDrot_dAB(iddd,2) + DDrot_dAB(iddd,3));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[4]*PDrot_dAB(cnt1,1)*(DDrot(iddd,4) + DDrot(iddd,5)) + intn[5]*(PDrot_dAB(cnt1,2) + PDrot_dAB(cnt1,3))*(DDrot(iddd,4) + DDrot(iddd,5));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[4]*PDrot_dA(cnt1,1)*(DDrot_dB(iddd,4) + DDrot_dB(iddd,5)) + intn[5]*(PDrot_dA(cnt1,2) + PDrot_dA(cnt1,3))*(DDrot_dB(iddd,4) + DDrot_dB(iddd,5));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[4]*PDrot_dB(cnt1,1)*(DDrot_dA(iddd,4) + DDrot_dA(iddd,5)) + intn[5]*(PDrot_dB(cnt1,2) + PDrot_dB(cnt1,3))*(DDrot_dA(iddd,4) + DDrot_dA(iddd,5));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[4]*PDrot(cnt1,1)*(DDrot_dAB(iddd,4) + DDrot_dAB(iddd,5)) + intn[5]*(PDrot(cnt1,2) + PDrot(cnt1,3))*(DDrot_dAB(iddd,4) + DDrot_dAB(iddd,5));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[6]*(PDrot_dAB(cnt1,4)*DDrot(iddd,6) + PDrot_dAB(cnt1,5)*DDrot(iddd,7)) + intn[7]*(PDrot_dAB(cnt1,8)*DDrot(iddd,6) + PDrot_dAB(cnt1,12)*DDrot(iddd,7));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[6]*(PDrot_dA(cnt1,4)*DDrot_dB(iddd,6) + PDrot_dA(cnt1,5)*DDrot_dB(iddd,7)) + intn[7]*(PDrot_dA(cnt1,8)*DDrot_dB(iddd,6) + PDrot_dA(cnt1,12)*DDrot_dB(iddd,7));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[6]*(PDrot_dB(cnt1,4)*DDrot_dA(iddd,6) + PDrot_dB(cnt1,5)*DDrot_dA(iddd,7)) + intn[7]*(PDrot_dB(cnt1,8)*DDrot_dA(iddd,6) + PDrot_dB(cnt1,12)*DDrot_dA(iddd,7));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[6]*(PDrot(cnt1,4)*DDrot_dAB(iddd,6) + PDrot(cnt1,5)*DDrot_dAB(iddd,7)) + intn[7]*(PDrot(cnt1,8)*DDrot_dAB(iddd,6) + PDrot(cnt1,12)*DDrot_dAB(iddd,7));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[8]*((PDrot_dAB(cnt1,10) + PDrot_dAB(cnt1,15))*(DDrot(iddd,11) + DDrot(iddd,14)) + (PDrot_dAB(cnt1,14) - PDrot_dAB(cnt1,11))*(DDrot(iddd,13) - DDrot(iddd,12)));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[8]*((PDrot_dA(cnt1,10) + PDrot_dA(cnt1,15))*(DDrot_dB(iddd,11) + DDrot_dB(iddd,14)) + (PDrot_dA(cnt1,14) - PDrot_dA(cnt1,11))*(DDrot_dB(iddd,13) - DDrot_dB(iddd,12)));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[8]*((PDrot_dB(cnt1,10) + PDrot_dB(cnt1,15))*(DDrot_dA(iddd,11) + DDrot_dA(iddd,14)) + (PDrot_dB(cnt1,14) - PDrot_dB(cnt1,11))*(DDrot_dA(iddd,13) - DDrot_dA(iddd,12)));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[8]*((PDrot(cnt1,10) + PDrot(cnt1,15))*(DDrot_dAB(iddd,11) + DDrot_dAB(iddd,14)) + (PDrot(cnt1,14) - PDrot(cnt1,11))*(DDrot_dAB(iddd,13) - DDrot_dAB(iddd,12)));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[9]*(PDrot_dAB(cnt1,8)*(DDrot(iddd,11) + DDrot(iddd,14)) + PDrot_dAB(cnt1,12)*(DDrot(iddd,12) - DDrot(iddd,13)));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[9]*(PDrot_dA(cnt1,8)*(DDrot_dB(iddd,11) + DDrot_dB(iddd,14)) + PDrot_dA(cnt1,12)*(DDrot_dB(iddd,12) - DDrot_dB(iddd,13)));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[9]*(PDrot_dB(cnt1,8)*(DDrot_dA(iddd,11) + DDrot_dA(iddd,14)) + PDrot_dB(cnt1,12)*(DDrot_dA(iddd,12) - DDrot_dA(iddd,13)));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[9]*(PDrot(cnt1,8)*(DDrot_dAB(iddd,11) + DDrot_dAB(iddd,14)) + PDrot(cnt1,12)*(DDrot_dAB(iddd,12) - DDrot_dAB(iddd,13)));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[10]*(PDrot_dAB(cnt1,4)*(DDrot(iddd,11) + DDrot(iddd,14)) + PDrot_dB(cnt1,5)*(DDrot_dA(iddd,12) - DDrot_dA(iddd,13)));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[10]*(PDrot_dA(cnt1,4)*(DDrot_dB(iddd,11) + DDrot_dB(iddd,14)) + PDrot(cnt1,5)*(DDrot_dAB(iddd,12) - DDrot_dAB(iddd,13)));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[10]*(PDrot_dB(cnt1,4)*(DDrot_dA(iddd,11) + DDrot_dA(iddd,14)) + PDrot_dAB(cnt1,5)*(DDrot(iddd,12) - DDrot(iddd,13)));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[10]*(PDrot(cnt1,4)*(DDrot_dAB(iddd,11) + DDrot_dAB(iddd,14)) + PDrot_dA(cnt1,5)*(DDrot_dB(iddd,12) - DDrot_dB(iddd,13)));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[11]*((PDrot_dAB(cnt1,10) + PDrot_dAB(cnt1,15))*DDrot(iddd,6) + (PDrot_dB(cnt1,11) - PDrot_dB(cnt1,14))*DDrot_dA(iddd,7));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[11]*((PDrot_dA(cnt1,10) + PDrot_dA(cnt1,15))*DDrot_dB(iddd,6) + (PDrot(cnt1,11) - PDrot(cnt1,14))*DDrot_dAB(iddd,7));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[11]*((PDrot_dB(cnt1,10) + PDrot_dB(cnt1,15))*DDrot_dA(iddd,6) + (PDrot_dAB(cnt1,11) - PDrot_dAB(cnt1,14))*DDrot(iddd,7));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[11]*((PDrot(cnt1,10) + PDrot(cnt1,15))*DDrot_dAB(iddd,6) + (PDrot_dA(cnt1,11) - PDrot_dA(cnt1,14))*DDrot_dB(iddd,7));
            }
          }
        }
        //(dd|pd) integrals
        Dd[0] = Dvalue(atmC,6);
        Dd[1] = Dvalue(atmD,5);
        intn[0] = d_eri2Center(2,0,2,0,1,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(2,1,2,1,1,0,2,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[2] = d_eri2Center(2,0,2,0,1,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[3] = d_eri2Center(2,1,2,1,1,1,2,1,RCD,Dd,atmC,atmD,type);    //(pi,pi|pi,pi)=(-pi,-pi|-pi,-pi)=(pi,pi|-pi,-pi)=(-pi,-pi|pi,pi)
        intn[4] = d_eri2Center(2,2,2,2,1,0,2,0,RCD,Dd,atmC,atmD,type);    //(dl,dl|sg,sg)=(-dl,-dl|sg,sg)
        intn[5] = d_eri2Center(2,2,2,2,1,1,2,1,RCD,Dd,atmC,atmD,type);    //(dl,dl|pi,pi)=(-dl,-dl|-pi,-pi)=(-dl,-dl|pi,pi)=(dl,dl|-pi,-pi)
        intn[6] = d_eri2Center(2,0,2,1,1,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[7] = d_eri2Center(2,1,2,0,1,1,2,0,RCD,Dd,atmC,atmD,type);    //(pi,sg|pi,sg)=(-pi,sg|-pi,sg)
        intn[8] = d_eri2Center(2,1,2,2,1,1,2,2,RCD,Dd,atmC,atmD,type);    //(pi,dl|pi,dl)=(-pi,-dl|-pi,-dl)=(-pi,-dl|pi,dl)=(pi,dl|-pi,-dl)=(-pi,dl|-pi,dl)=(pi,-dl|pi,-dl)=-(-pi,dl|pi,-dl)=-(pi,-dl|-pi,dl)
        intn[9] = d_eri2Center(2,1,2,2,1,1,2,0,RCD,Dd,atmC,atmD,type);    //(pi,dl|pi,sg)=-(-pi,dl|-pi,sg)=(pi,-dl|-pi,sg)=(-pi,-dl|pi,sg)
        intn[10] = d_eri2Center(2,1,2,2,1,0,2,1,RCD,Dd,atmC,atmD,type);   //(pi,dl|sg,pi)=-(-pi,dl|sg,-pi)=(pi,-dl|sg,-pi)=(-pi,-dl|sg,pi)
        intn[11] = d_eri2Center(2,0,2,1,1,1,2,2,RCD,Dd,atmC,atmD,type);   //(sg,pi|pi,dl)=(sg,pi|-pi,-dl)=-(sg,-pi|-pi,dl)=(sg,-pi|pi,-dl)
        //(dd|pd) rotate
        for (size_t iddd = 1; iddd < 16; ++iddd) {
          for (size_t idpd1 = 1; idpd1 < 6; ++idpd1) {
            for (size_t idpd2 = 1; idpd2 < 4; ++idpd2) {
              cnt1 = (idpd1 - 1)*3 + idpd2;
              block(30 + iddd,5*(idpd2+2) + idpd1) = -intn[0]*PDrot_dAB(cnt1,1)*DDrot(iddd,1) - intn[1]*PDrot_dAB(cnt1,1)*(DDrot(iddd,2) + DDrot(iddd,3));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[0]*PDrot_dA(cnt1,1)*DDrot_dB(iddd,1) + intn[1]*PDrot_dA(cnt1,1)*(DDrot_dB(iddd,2) + DDrot_dB(iddd,3));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[0]*PDrot_dB(cnt1,1)*DDrot_dA(iddd,1) + intn[1]*PDrot_dB(cnt1,1)*(DDrot_dA(iddd,2) + DDrot_dA(iddd,3));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[0]*PDrot(cnt1,1)*DDrot_dAB(iddd,1) + intn[1]*PDrot(cnt1,1)*(DDrot_dAB(iddd,2) + DDrot_dAB(iddd,3));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[2]*(PDrot_dAB(cnt1,2) + PDrot_dAB(cnt1,3))*DDrot(iddd,1) + intn[3]*(PDrot_dB(cnt1,2) + PDrot_dB(cnt1,3))*(DDrot_dA(iddd,2) + DDrot_dA(iddd,3));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[2]*(PDrot_dA(cnt1,2) + PDrot_dA(cnt1,3))*DDrot_dB(iddd,1) + intn[3]*(PDrot(cnt1,2) + PDrot(cnt1,3))*(DDrot_dAB(iddd,2) + DDrot_dAB(iddd,3));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[2]*(PDrot_dB(cnt1,2) + PDrot_dB(cnt1,3))*DDrot_dA(iddd,1) + intn[3]*(PDrot_dAB(cnt1,2) + PDrot_dAB(cnt1,3))*(DDrot(iddd,2) + DDrot(iddd,3));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[2]*(PDrot(cnt1,2) + PDrot(cnt1,3))*DDrot_dAB(iddd,1) + intn[3]*(PDrot_dA(cnt1,2) + PDrot_dA(cnt1,3))*(DDrot_dB(iddd,2) + DDrot_dB(iddd,3));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[4]*PDrot_dAB(cnt1,1)*(DDrot(iddd,4) + DDrot(iddd,5)) + intn[5]*(PDrot_dAB(cnt1,2) + PDrot_dAB(cnt1,3))*(DDrot(iddd,4) + DDrot(iddd,5));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[4]*PDrot_dA(cnt1,1)*(DDrot_dB(iddd,4) + DDrot_dB(iddd,5)) + intn[5]*(PDrot_dA(cnt1,2) + PDrot_dA(cnt1,3))*(DDrot_dB(iddd,4) + DDrot_dB(iddd,5));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[4]*PDrot_dB(cnt1,1)*(DDrot_dA(iddd,4) + DDrot_dA(iddd,5)) + intn[5]*(PDrot_dB(cnt1,2) + PDrot_dB(cnt1,3))*(DDrot_dA(iddd,4) + DDrot_dA(iddd,5));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[4]*PDrot(cnt1,1)*(DDrot_dAB(iddd,4) + DDrot_dAB(iddd,5)) + intn[5]*(PDrot(cnt1,2) + PDrot(cnt1,3))*(DDrot_dAB(iddd,4) + DDrot_dAB(iddd,5));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[6]*(PDrot_dAB(cnt1,4)*DDrot(iddd,6) + PDrot_dAB(cnt1,5)*DDrot(iddd,7)) + intn[7]*(PDrot_dAB(cnt1,8)*DDrot(iddd,6) + PDrot_dAB(cnt1,12)*DDrot(iddd,7));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[6]*(PDrot_dA(cnt1,4)*DDrot_dB(iddd,6) + PDrot_dA(cnt1,5)*DDrot_dB(iddd,7)) + intn[7]*(PDrot_dA(cnt1,8)*DDrot_dB(iddd,6) + PDrot_dA(cnt1,12)*DDrot_dB(iddd,7));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[6]*(PDrot_dB(cnt1,4)*DDrot_dA(iddd,6) + PDrot_dB(cnt1,5)*DDrot_dA(iddd,7)) + intn[7]*(PDrot_dB(cnt1,8)*DDrot_dA(iddd,6) + PDrot_dB(cnt1,12)*DDrot_dA(iddd,7));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[6]*(PDrot(cnt1,4)*DDrot_dAB(iddd,6) + PDrot(cnt1,5)*DDrot_dAB(iddd,7)) + intn[7]*(PDrot(cnt1,8)*DDrot_dAB(iddd,6) + PDrot(cnt1,12)*DDrot_dAB(iddd,7));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[8]*((PDrot_dAB(cnt1,10) + PDrot_dAB(cnt1,15))*(DDrot(iddd,11) + DDrot(iddd,14)) + (PDrot_dAB(cnt1,14) - PDrot_dAB(cnt1,11))*(DDrot(iddd,13) - DDrot(iddd,12)));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[8]*((PDrot_dA(cnt1,10) + PDrot_dA(cnt1,15))*(DDrot_dB(iddd,11) + DDrot_dB(iddd,14)) + (PDrot_dA(cnt1,14) - PDrot_dA(cnt1,11))*(DDrot_dB(iddd,13) - DDrot_dB(iddd,12)));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[8]*((PDrot_dB(cnt1,10) + PDrot_dB(cnt1,15))*(DDrot_dA(iddd,11) + DDrot_dA(iddd,14)) + (PDrot_dB(cnt1,14) - PDrot_dB(cnt1,11))*(DDrot_dA(iddd,13) - DDrot_dA(iddd,12)));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[8]*((PDrot(cnt1,10) + PDrot(cnt1,15))*(DDrot_dAB(iddd,11) + DDrot_dAB(iddd,14)) + (PDrot(cnt1,14) - PDrot(cnt1,11))*(DDrot_dAB(iddd,13) - DDrot_dAB(iddd,12)));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[9]*(PDrot_dAB(cnt1,8)*(DDrot(iddd,11) + DDrot(iddd,14)) + PDrot_dB(cnt1,12)*(DDrot_dA(iddd,12) - DDrot_dA(iddd,13)));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[9]*(PDrot_dA(cnt1,8)*(DDrot_dB(iddd,11) + DDrot_dB(iddd,14)) + PDrot(cnt1,12)*(DDrot_dAB(iddd,12) - DDrot_dAB(iddd,13)));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[9]*(PDrot_dB(cnt1,8)*(DDrot_dA(iddd,11) + DDrot_dA(iddd,14)) + PDrot_dAB(cnt1,12)*(DDrot(iddd,12) - DDrot(iddd,13)));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[9]*(PDrot(cnt1,8)*(DDrot_dAB(iddd,11) + DDrot_dAB(iddd,14)) + PDrot_dA(cnt1,12)*(DDrot_dB(iddd,12) - DDrot_dB(iddd,13)));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[10]*(PDrot_dAB(cnt1,4)*(DDrot(iddd,11) + DDrot(iddd,14)) + PDrot_dB(cnt1,5)*(DDrot_dA(iddd,12) - DDrot_dA(iddd,13)));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[10]*(PDrot_dA(cnt1,4)*(DDrot_dB(iddd,11) + DDrot_dB(iddd,14)) + PDrot(cnt1,5)*(DDrot_dAB(iddd,12) - DDrot_dAB(iddd,13)));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[10]*(PDrot_dB(cnt1,4)*(DDrot_dA(iddd,11) + DDrot_dA(iddd,14)) + PDrot_dAB(cnt1,5)*(DDrot(iddd,12) - DDrot(iddd,13)));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[10]*(PDrot(cnt1,4)*(DDrot_dAB(iddd,11) + DDrot_dAB(iddd,14)) + PDrot_dA(cnt1,5)*(DDrot_dB(iddd,12) - DDrot_dB(iddd,13)));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[11]*((PDrot_dAB(cnt1,10) + PDrot_dAB(cnt1,15))*DDrot(iddd,6) + (PDrot_dB(cnt1,11) - PDrot_dB(cnt1,14))*DDrot_dA(iddd,7));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[11]*((PDrot_dA(cnt1,10) + PDrot_dA(cnt1,15))*DDrot_dB(iddd,6) + (PDrot(cnt1,11) - PDrot(cnt1,14))*DDrot_dAB(iddd,7));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[11]*((PDrot_dB(cnt1,10) + PDrot_dB(cnt1,15))*DDrot_dA(iddd,6) + (PDrot_dAB(cnt1,11) - PDrot_dAB(cnt1,14))*DDrot(iddd,7));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[11]*((PDrot(cnt1,10) + PDrot(cnt1,15))*DDrot_dAB(iddd,6) + (PDrot_dA(cnt1,11) - PDrot_dA(cnt1,14))*DDrot_dB(iddd,7));
            }
          }
        }
        //(dd|dd) integrals
        Dd[0] = Dvalue(atmC,6);
        Dd[1] = Dvalue(atmD,6);
        intn[0] = d_eri2Center(2,0,2,0,2,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(2,0,2,0,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center(2,0,2,0,2,2,2,2,RCD,Dd,atmC,atmD,type);    //(sg,sg|dl,dl)=(sg,sg|-dl,-dl)
        intn[3] = d_eri2Center(2,1,2,1,2,0,2,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[4] = d_eri2Center(2,2,2,2,2,0,2,0,RCD,Dd,atmC,atmD,type);    //(dl,dl|sg,sg)=(-dl,-dl|sg,sg)
        intn[5] = d_eri2Center(2,2,2,2,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(dl,dl|pi,pi)=(-dl,-dl|pi,pi)=(dl,dl|-pi,-pi)=(-dl,-dl|-pi,-pi)
        intn[6] = d_eri2Center(2,1,2,1,2,2,2,2,RCD,Dd,atmC,atmD,type);    //(pi,pi|dl,dl)=(-pi,-pi|dl,dl)=(pi,pi|-dl,-dl)=(-pi,-pi|-dl,-dl)
        intn[7] = d_eri2Center(2,1,2,1,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(pi,pi|pi,pi)=(-pi,-pi|-pi,-pi)
        intn[8] = d_eri2Center(2,2,2,2,2,2,2,2,RCD,Dd,atmC,atmD,type);    //(dl,dl|dl,dl)=(-dl,-dl|-dl,-dl)
        intn[9] = d_eri2Center(2,1,2,1,2,-1,2,-1,RCD,Dd,atmC,atmD,type);  //(pi,pi|-pi,-pi)=(-pi,-pi|pi,pi)
        intn[10] = d_eri2Center(2,-2,2,-2,2,2,2,2,RCD,Dd,atmC,atmD,type); //(-dl,-dl|dl,dl)=(dl,dl|-dl,-dl)
        intn[11] = d_eri2Center(2,0,2,1,2,0,2,1,RCD,Dd,atmC,atmD,type);   //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[12] = d_eri2Center(2,0,2,2,2,0,2,2,RCD,Dd,atmC,atmD,type);   //(sg,dl|sg,dl)=(sg,-dl|sg,-dl)
        intn[13] = d_eri2Center(2,1,2,2,2,1,2,2,RCD,Dd,atmC,atmD,type);   //(pi,dl|pi,dl)=(pi,dl|-pi,-dl)=(-pi,-dl|pi,dl)=(-pi,-dl|-pi,-dl)=-(-pi,dl|pi,-dl)=-(pi,-dl|-pi,dl)=(-pi,dl|-pi,dl)=(pi,-dl|pi,-dl)
        intn[14] = d_eri2Center(2,1,2,-1,2,1,2,-1,RCD,Dd,atmC,atmD,type); //(pi,-pi|pi,-pi)
        intn[15] = d_eri2Center(2,1,2,-2,2,0,2,-1,RCD,Dd,atmC,atmD,type); //(pi,dl|sg,pi)=(-pi,-dl|sg,pi)=(pi,-dl|sg,-pi)=-(-pi,dl|sg,-pi)=(sg,pi|pi,dl)=(sg,pi|-pi,-dl)=(sg,-pi|pi,-dl)=-(sg,-pi|-pi,dl)
        intn[16] = d_eri2Center(2,1,2,1,2,0,2,2,RCD,Dd,atmC,atmD,type);   //(pi,pi|sg,dl)=-(-pi,-pi|sg,dl)=(pi,-pi|sg,-dl)=(sg,dl|pi,pi)=-(sg,dl|-pi,-pi)=(sg,-dl|pi,-pi)
        //(dd|dd) rotate
        for (size_t iddd1 = 1; iddd1 < 16; ++iddd1) {
          for (size_t iddd2 = 1; iddd2 < 16; ++iddd2) {
            block(30 + iddd1,30 + iddd2) = intn[0]*DDrot_dAB(iddd1,1)*DDrot(iddd2,1) + intn[1]*DDrot_dAB(iddd1,1)*(DDrot(iddd2,2) + DDrot(iddd2,3)) + intn[2]*DDrot_dB(iddd1,1)*(DDrot_dA(iddd2,4) + DDrot_dA(iddd2,5));
            block(30 + iddd1,30 + iddd2) += intn[0]*DDrot_dA(iddd1,1)*DDrot_dB(iddd2,1) + intn[1]*DDrot_dA(iddd1,1)*(DDrot_dB(iddd2,2) + DDrot_dB(iddd2,3)) + intn[2]*DDrot(iddd1,1)*(DDrot_dAB(iddd2,4) + DDrot_dAB(iddd2,5));
            block(30 + iddd1,30 + iddd2) += intn[0]*DDrot_dB(iddd1,1)*DDrot_dA(iddd2,1) + intn[1]*DDrot_dB(iddd1,1)*(DDrot_dA(iddd2,2) + DDrot_dA(iddd2,3)) + intn[2]*DDrot_dAB(iddd1,1)*(DDrot(iddd2,4) + DDrot(iddd2,5));
            block(30 + iddd1,30 + iddd2) += intn[0]*DDrot(iddd1,1)*DDrot_dAB(iddd2,1) + intn[1]*DDrot(iddd1,1)*(DDrot_dAB(iddd2,2) + DDrot_dAB(iddd2,3)) + intn[2]*DDrot_dA(iddd1,1)*(DDrot_dB(iddd2,4) + DDrot_dB(iddd2,5));
            block(30 + iddd1,30 + iddd2) += intn[3]*(DDrot_dAB(iddd1,2) + DDrot_dAB(iddd1,3))*DDrot(iddd2,1) + intn[4]*(DDrot_dB(iddd1,4) + DDrot_dB(iddd1,5))*DDrot_dA(iddd2,1);
            block(30 + iddd1,30 + iddd2) += intn[3]*(DDrot_dA(iddd1,2) + DDrot_dA(iddd1,3))*DDrot_dB(iddd2,1) + intn[4]*(DDrot(iddd1,4) + DDrot(iddd1,5))*DDrot_dAB(iddd2,1);
            block(30 + iddd1,30 + iddd2) += intn[3]*(DDrot_dB(iddd1,2) + DDrot_dB(iddd1,3))*DDrot_dA(iddd2,1) + intn[4]*(DDrot_dAB(iddd1,4) + DDrot_dAB(iddd1,5))*DDrot(iddd2,1);
            block(30 + iddd1,30 + iddd2) += intn[3]*(DDrot(iddd1,2) + DDrot(iddd1,3))*DDrot_dAB(iddd2,1) + intn[4]*(DDrot_dA(iddd1,4) + DDrot_dA(iddd1,5))*DDrot_dB(iddd2,1);
            block(30 + iddd1,30 + iddd2) += intn[5]*(DDrot_dAB(iddd1,4) + DDrot_dAB(iddd1,5))*(DDrot(iddd2,2) + DDrot(iddd2,3)) + intn[6]*(DDrot_dAB(iddd1,2) + DDrot_dAB(iddd1,3))*(DDrot(iddd2,4) + DDrot(iddd2,5));
            block(30 + iddd1,30 + iddd2) += intn[5]*(DDrot_dA(iddd1,4) + DDrot_dA(iddd1,5))*(DDrot_dB(iddd2,2) + DDrot_dB(iddd2,3)) + intn[6]*(DDrot_dA(iddd1,2) + DDrot_dA(iddd1,3))*(DDrot_dB(iddd2,4) + DDrot_dB(iddd2,5));
            block(30 + iddd1,30 + iddd2) += intn[5]*(DDrot_dB(iddd1,4) + DDrot_dB(iddd1,5))*(DDrot_dA(iddd2,2) + DDrot_dA(iddd2,3)) + intn[6]*(DDrot_dB(iddd1,2) + DDrot_dB(iddd1,3))*(DDrot_dA(iddd2,4) + DDrot_dA(iddd2,5));
            block(30 + iddd1,30 + iddd2) += intn[5]*(DDrot(iddd1,4) + DDrot(iddd1,5))*(DDrot_dAB(iddd2,2) + DDrot_dAB(iddd2,3)) + intn[6]*(DDrot(iddd1,2) + DDrot(iddd1,3))*(DDrot_dAB(iddd2,4) + DDrot_dAB(iddd2,5));
            block(30 + iddd1,30 + iddd2) += intn[7]*(DDrot_dAB(iddd1,2)*DDrot(iddd2,2) + DDrot_dAB(iddd1,3)*DDrot(iddd2,3)) + intn[8]*(DDrot_dAB(iddd1,4)*DDrot(iddd2,4) + DDrot_dAB(iddd1,5)*DDrot(iddd2,5));
            block(30 + iddd1,30 + iddd2) += intn[7]*(DDrot_dA(iddd1,2)*DDrot_dB(iddd2,2) + DDrot_dA(iddd1,3)*DDrot_dB(iddd2,3)) + intn[8]*(DDrot_dA(iddd1,4)*DDrot_dB(iddd2,4) + DDrot_dA(iddd1,5)*DDrot_dB(iddd2,5));
            block(30 + iddd1,30 + iddd2) += intn[7]*(DDrot_dB(iddd1,2)*DDrot_dA(iddd2,2) + DDrot_dB(iddd1,3)*DDrot_dA(iddd2,3)) + intn[8]*(DDrot_dB(iddd1,4)*DDrot_dA(iddd2,4) + DDrot_dB(iddd1,5)*DDrot_dA(iddd2,5));
            block(30 + iddd1,30 + iddd2) += intn[7]*(DDrot(iddd1,2)*DDrot_dAB(iddd2,2) + DDrot(iddd1,3)*DDrot_dAB(iddd2,3)) + intn[8]*(DDrot(iddd1,4)*DDrot_dAB(iddd2,4) + DDrot(iddd1,5)*DDrot_dAB(iddd2,5));
            block(30 + iddd1,30 + iddd2) += intn[9]*(DDrot_dAB(iddd1,2)*DDrot(iddd2,3) + DDrot_dAB(iddd1,3)*DDrot(iddd2,2)) + intn[10]*(DDrot_dAB(iddd1,4)*DDrot(iddd2,5) + DDrot_dAB(iddd1,5)*DDrot(iddd2,4));
            block(30 + iddd1,30 + iddd2) += intn[9]*(DDrot_dA(iddd1,2)*DDrot_dB(iddd2,3) + DDrot_dA(iddd1,3)*DDrot_dB(iddd2,2)) + intn[10]*(DDrot_dA(iddd1,4)*DDrot_dB(iddd2,5) + DDrot_dA(iddd1,5)*DDrot_dB(iddd2,4));
            block(30 + iddd1,30 + iddd2) += intn[9]*(DDrot_dB(iddd1,2)*DDrot_dA(iddd2,3) + DDrot_dB(iddd1,3)*DDrot_dA(iddd2,2)) + intn[10]*(DDrot_dB(iddd1,4)*DDrot_dA(iddd2,5) + DDrot_dB(iddd1,5)*DDrot_dA(iddd2,4));
            block(30 + iddd1,30 + iddd2) += intn[9]*(DDrot(iddd1,2)*DDrot_dAB(iddd2,3) + DDrot(iddd1,3)*DDrot_dAB(iddd2,2)) + intn[10]*(DDrot(iddd1,4)*DDrot_dAB(iddd2,5) + DDrot(iddd1,5)*DDrot_dAB(iddd2,4));
            block(30 + iddd1,30 + iddd2) += intn[11]*(DDrot_dAB(iddd1,6)*DDrot(iddd2,6) + DDrot_dAB(iddd1,7)*DDrot(iddd2,7)) + intn[12]*(DDrot_dAB(iddd1,8)*DDrot(iddd2,8) + DDrot_dAB(iddd1,9)*DDrot(iddd2,9));
            block(30 + iddd1,30 + iddd2) += intn[11]*(DDrot_dA(iddd1,6)*DDrot_dB(iddd2,6) + DDrot_dA(iddd1,7)*DDrot_dB(iddd2,7)) + intn[12]*(DDrot_dA(iddd1,8)*DDrot_dB(iddd2,8) + DDrot_dA(iddd1,9)*DDrot_dB(iddd2,9));
            block(30 + iddd1,30 + iddd2) += intn[11]*(DDrot_dB(iddd1,6)*DDrot_dA(iddd2,6) + DDrot_dB(iddd1,7)*DDrot_dA(iddd2,7)) + intn[12]*(DDrot_dB(iddd1,8)*DDrot_dA(iddd2,8) + DDrot_dB(iddd1,9)*DDrot_dA(iddd2,9));
            block(30 + iddd1,30 + iddd2) += intn[11]*(DDrot(iddd1,6)*DDrot_dAB(iddd2,6) + DDrot(iddd1,7)*DDrot_dAB(iddd2,7)) + intn[12]*(DDrot(iddd1,8)*DDrot_dAB(iddd2,8) + DDrot(iddd1,9)*DDrot_dAB(iddd2,9));
            block(30 + iddd1,30 + iddd2) += intn[13]*((DDrot_dAB(iddd1,11) + DDrot_dAB(iddd1,14))*(DDrot(iddd2,11) + DDrot(iddd2,14)) + (DDrot_dAB(iddd1,12) - DDrot_dAB(iddd1,13))*(DDrot(iddd2,12) - DDrot(iddd2,13)));
            block(30 + iddd1,30 + iddd2) += intn[13]*((DDrot_dA(iddd1,11) + DDrot_dA(iddd1,14))*(DDrot_dB(iddd2,11) + DDrot_dB(iddd2,14)) + (DDrot_dA(iddd1,12) - DDrot_dA(iddd1,13))*(DDrot_dB(iddd2,12) - DDrot_dB(iddd2,13)));
            block(30 + iddd1,30 + iddd2) += intn[13]*((DDrot_dB(iddd1,11) + DDrot_dB(iddd1,14))*(DDrot_dA(iddd2,11) + DDrot_dA(iddd2,14)) + (DDrot_dB(iddd1,12) - DDrot_dB(iddd1,13))*(DDrot_dA(iddd2,12) - DDrot_dA(iddd2,13)));
            block(30 + iddd1,30 + iddd2) += intn[13]*((DDrot(iddd1,11) + DDrot(iddd1,14))*(DDrot_dAB(iddd2,11) + DDrot_dAB(iddd2,14)) + (DDrot(iddd1,12) - DDrot(iddd1,13))*(DDrot_dAB(iddd2,12) - DDrot_dAB(iddd2,13)));
            block(30 + iddd1,30 + iddd2) += intn[14]*DDrot_dAB(iddd1,10)*DDrot(iddd2,10) + intn[14]*DDrot_dB(iddd1,10)*DDrot_dA(iddd2,10);
            block(30 + iddd1,30 + iddd2) += intn[14]*DDrot_dA(iddd1,10)*DDrot_dB(iddd2,10) + intn[14]*DDrot(iddd1,10)*DDrot_dAB(iddd2,10);
            block(30 + iddd1,30 + iddd2) += intn[15]*((DDrot_dAB(iddd1,11) + DDrot_dAB(iddd1,14))*DDrot(iddd2,6) + (DDrot_dB(iddd1,12) - DDrot_dB(iddd1,13))*DDrot_dA(iddd2,7) + DDrot_dAB(iddd1,6)*(DDrot(iddd2,11) + DDrot(iddd2,14)) + DDrot_dB(iddd1,7)*(DDrot_dA(iddd2,12) - DDrot_dA(iddd2,13)));
            block(30 + iddd1,30 + iddd2) += intn[15]*((DDrot_dA(iddd1,11) + DDrot_dA(iddd1,14))*DDrot_dB(iddd2,6) + (DDrot(iddd1,12) - DDrot(iddd1,13))*DDrot_dAB(iddd2,7) + DDrot_dA(iddd1,6)*(DDrot_dB(iddd2,11) + DDrot_dB(iddd2,14)) + DDrot(iddd1,7)*(DDrot_dAB(iddd2,12) - DDrot_dAB(iddd2,13)));
            block(30 + iddd1,30 + iddd2) += intn[15]*((DDrot_dB(iddd1,11) + DDrot_dB(iddd1,14))*DDrot_dA(iddd2,6) + (DDrot_dAB(iddd1,12) - DDrot_dAB(iddd1,13))*DDrot(iddd2,7) + DDrot_dB(iddd1,6)*(DDrot_dA(iddd2,11) + DDrot_dA(iddd2,14)) + DDrot_dAB(iddd1,7)*(DDrot(iddd2,12) - DDrot(iddd2,13)));
            block(30 + iddd1,30 + iddd2) += intn[15]*((DDrot(iddd1,11) + DDrot(iddd1,14))*DDrot_dAB(iddd2,6) + (DDrot_dA(iddd1,12) - DDrot_dA(iddd1,13))*DDrot_dB(iddd2,7) + DDrot(iddd1,6)*(DDrot_dAB(iddd2,11) + DDrot_dAB(iddd2,14)) + DDrot_dA(iddd1,7)*(DDrot_dB(iddd2,12) - DDrot_dB(iddd2,13)));
            block(30 + iddd1,30 + iddd2) += intn[16]*((DDrot_dAB(iddd1,2) - DDrot_dAB(iddd1,3))*DDrot(iddd2,8) + DDrot_dAB(iddd1,10)*DDrot(iddd2,9) + DDrot_dAB(iddd1,8)*(DDrot(iddd2,2) - DDrot(iddd2,3)) + DDrot_dAB(iddd1,9)*DDrot(iddd2,10));
            block(30 + iddd1,30 + iddd2) += intn[16]*((DDrot_dA(iddd1,2) - DDrot_dA(iddd1,3))*DDrot_dB(iddd2,8) + DDrot_dA(iddd1,10)*DDrot_dB(iddd2,9) + DDrot_dA(iddd1,8)*(DDrot_dB(iddd2,2) - DDrot_dB(iddd2,3)) + DDrot_dA(iddd1,9)*DDrot_dB(iddd2,10));
            block(30 + iddd1,30 + iddd2) += intn[16]*((DDrot_dB(iddd1,2) - DDrot_dB(iddd1,3))*DDrot_dA(iddd2,8) + DDrot_dB(iddd1,10)*DDrot_dA(iddd2,9) + DDrot_dB(iddd1,8)*(DDrot_dA(iddd2,2) - DDrot_dA(iddd2,3)) + DDrot_dB(iddd1,9)*DDrot_dA(iddd2,10));
            block(30 + iddd1,30 + iddd2) += intn[16]*((DDrot(iddd1,2) - DDrot(iddd1,3))*DDrot_dAB(iddd2,8) + DDrot(iddd1,10)*DDrot_dAB(iddd2,9) + DDrot(iddd1,8)*(DDrot_dAB(iddd2,2) - DDrot_dAB(iddd2,3)) + DDrot(iddd1,9)*DDrot_dAB(iddd2,10));
          }
        }
      }
    }
  }
  void IntegralBlock2C_dRA(int type, matrixE & block, int der, int atmC, int atmD, double RCD, double cost, double sint, double cosp, double sinp) {
    //function calculating a block of second-derivatives of two-center integrals with respect to an angle and the internuclear distance
    //2 -> theta; 3 -> phi
    matrixE SProt_dA(1,1);
    matrixE PProt_dA(1,1);
    matrixE SDrot_dA(1,1);
    matrixE PDrot_dA(1,1);
    matrixE DDrot_dA(1,1);
    int cnt1;
    int cnt2;
    size_t nrows = block.rows();
    size_t ncols = block.cols();
    if ((nrows > 1)||(ncols > 1)) {
      SProt = SPtransf(cost,sint,cosp,sinp);
      PProt = PPtransf(cost,sint,cosp,sinp);
      if (der == 2) {
        SProt_dA = SPtransf_dt(cost,sint,cosp,sinp);
        PProt_dA = PPtransf_dt(cost,sint,cosp,sinp);
      }
      else if (der == 3) {
        SProt_dA = SPtransf_dp(cost,sint,cosp,sinp);
        PProt_dA = PPtransf_dp(cost,sint,cosp,sinp);
      }
      if ((nrows > 10)||(ncols > 10)) {
        SDrot = SDtransf(cost,sint,cosp,sinp);
        PDrot = PDtransf(cost,sint,cosp,sinp);
        DDrot = DDtransf(cost,sint,cosp,sinp);
        if (der == 2) {
          SDrot_dA = SDtransf_dt(cost,sint,cosp,sinp);
          PDrot_dA = PDtransf_dt(cost,sint,cosp,sinp);
          DDrot_dA = DDtransf_dt(cost,sint,cosp,sinp);
        }
        else if (der == 3) {
          SDrot_dA = SDtransf_dp(cost,sint,cosp,sinp);
          PDrot_dA = PDtransf_dp(cost,sint,cosp,sinp);
          DDrot_dA = DDtransf_dp(cost,sint,cosp,sinp);
        }
      }
    }
    D[0] = Dvalue(atmC,1);
    D[1] = Dvalue(atmC,2);
    D[2] = Dvalue(atmD,1);
    D[3] = Dvalue(atmD,2);
    //integral calculation
    block(1,1) = 0.0;                                                                                                                   //(ss|ss)
    if (nrows > 1) {
      intn[0] = eri2Center_dR(0,0,1,0,0,0,0,0,RCD,D,atmC,atmD,type);
      intn[1] = eri2Center_dR(1,0,1,0,0,0,0,0,RCD,D,atmC,atmD,type);
      intn[2] = eri2Center_dR(1,1,1,1,0,0,0,0,RCD,D,atmC,atmD,type);
      for (size_t idx = 0; idx < 6; ++idx) {
        if (idx < 3) {block(idx + 2,1) = -intn[0]*SProt_dA(idx + 1,1);}                                                         //(sp|ss) integrals
        block(5 + idx,1) = intn[1]*PProt_dA(idx + 1,1) + intn[2]*(PProt_dA(idx + 1,2) + PProt_dA(idx + 1,3));   //(pp|ss) integrals
      }
      if (nrows > 10) {
        //(sd|ss) integrals
        Dd[0] = Dvalue(atmC,4);
        Dd[1] = Dvalue(atmD,1);
        intn[0] = d_eri2Center_dR(0,0,2,0,0,0,0,0,RCD,Dd,atmC,atmD,type);
        for (size_t idsd = 1; idsd < 6; ++idsd) {
          block(10 + idsd,1) = intn[0]*SDrot_dA(idsd,1);
        }
        //(pd|ss) integrals
        Dd[0] = Dvalue(atmC,5);
        intn[0] = d_eri2Center_dR(1,0,2,0,0,0,0,0,RCD,Dd,atmC,atmD,type);
        intn[1] = d_eri2Center_dR(1,1,2,1,0,0,0,0,RCD,Dd,atmC,atmD,type);
        for (int idpd1 = 1; idpd1 < 6; ++idpd1) {
          for (int idpd2 = 1; idpd2 < 4; ++idpd2) {
            block(5*(idpd2+2) + idpd1,1) = -(intn[0]*PDrot_dA((idpd1 - 1)*3 + idpd2,1) + intn[1]*(PDrot_dA((idpd1 - 1)*3 + idpd2,2) + PDrot_dA((idpd1 - 1)*3 + idpd2,3)));
          }
        }
        //(dd|ss) integrals
        Dd[0] = Dvalue(atmC,6);
        intn[0] = d_eri2Center_dR(2,0,2,0,0,0,0,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR(2,1,2,1,0,0,0,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[2] = d_eri2Center_dR(2,2,2,2,0,0,0,0,RCD,Dd,atmC,atmD,type);    //(dl,dl|sg,sg)=(-dl,-dl|sg,sg)
        for (size_t iddd = 1; iddd < 16; ++iddd) {
          block(30 + iddd,1) = intn[0]*DDrot_dA(iddd,1) + intn[1]*(DDrot_dA(iddd,2) + DDrot_dA(iddd,3)) + intn[2]*(DDrot_dA(iddd,4) + DDrot_dA(iddd,5));
        }
      }
    }
    if ((ncols > 1)&&(type == 0)) {
      intn[0] = eri2Center_dR(0,0,0,0,0,0,1,0,RCD,D,atmC,atmD,type);
      intn[1] = eri2Center_dR(0,0,0,0,1,0,1,0,RCD,D,atmC,atmD,type);
      intn[2] = eri2Center_dR(0,0,0,0,1,1,1,1,RCD,D,atmC,atmD,type);
      for (size_t idx = 0; idx < 6; ++idx) {
        if (idx < 3) {block(1,idx + 2) = -intn[0]*SProt_dA(idx + 1,1);}                                                         //(ss|sp) integrals
        block(1,5 + idx) = intn[1]*PProt_dA(idx + 1,1) + intn[2]*(PProt_dA(idx + 1,2) + PProt_dA(idx + 1,3));   //(ss|pp) integrals
      }
      if (ncols > 10) {
        //(ss|sd) integrals
        Dd[0] = Dvalue(atmC,1);
        Dd[1] = Dvalue(atmD,4);
        intn[0] = d_eri2Center_dR(0,0,0,0,0,0,2,0,RCD,Dd,atmC,atmD,type);
        for (size_t idsd = 1; idsd < 6; ++idsd) {
          block(1,10 + idsd) = intn[0]*SDrot_dA(idsd,1);
        }
        //(ss|pd) integrals
        Dd[1] = Dvalue(atmD,5);
        intn[0] = d_eri2Center_dR(0,0,0,0,1,0,2,0,RCD,Dd,atmC,atmD,type);
        intn[1] = d_eri2Center_dR(0,0,0,0,1,1,2,1,RCD,Dd,atmC,atmD,type);
        for (int idpd1 = 1; idpd1 < 6; ++idpd1) {
          for (int idpd2 = 1; idpd2 < 4; ++idpd2) {
            block(1,5*(idpd2 + 2) + idpd1) = -intn[0]*PDrot_dA((idpd1 - 1)*3 + idpd2,1) - intn[1]*(PDrot_dA((idpd1 - 1)*3 + idpd2,2) + PDrot_dA((idpd1 - 1)*3 + idpd2,3));
          }
        }
        //(ss|dd) integrals
        Dd[1] = Dvalue(atmD,6);
        intn[0] = d_eri2Center_dR(0,0,0,0,2,0,2,0,RCD,Dd,atmC,atmD,type);
        intn[1] = d_eri2Center_dR(0,0,0,0,2,1,2,1,RCD,Dd,atmC,atmD,type);
        intn[2] = d_eri2Center_dR(0,0,0,0,2,2,2,2,RCD,Dd,atmC,atmD,type);
        for (size_t iddd = 1; iddd < 16; ++iddd) {
          block(1,30 + iddd) = intn[0]*DDrot_dA(iddd,1) + intn[1]*(DDrot_dA(iddd,2) + DDrot_dA(iddd,3)) + intn[2]*(DDrot_dA(iddd,4) + DDrot_dA(iddd,5));
        }
      }
    }
    if ((nrows > 1)&&(ncols > 1)&&(type == 0)) {
      intn[0] = eri2Center_dR(0,0,1,0,0,0,1,0,RCD,D,atmC,atmD,type);
      intn[1] = eri2Center_dR(0,0,1,1,0,0,1,1,RCD,D,atmC,atmD,type);
      intn[2] = eri2Center_dR(0,0,1,0,1,0,1,0,RCD,D,atmC,atmD,type);
      intn[3] = eri2Center_dR(0,0,1,0,1,1,1,1,RCD,D,atmC,atmD,type);
      intn[4] = eri2Center_dR(0,0,1,1,1,0,1,1,RCD,D,atmC,atmD,type);
      for (size_t idpb = 1; idpb < 4; ++idpb) {
        for (size_t idpk = 1; idpk < 4; ++idpk) {           //(sp|sp)
          block(idpb + 1,idpk + 1) = intn[0]*SProt_dA(idpb,1)*SProt(idpk,1) + intn[1]*(SProt_dA(idpb,2)*SProt(idpk,2) + SProt_dA(idpb,3)*SProt(idpk,3));
          block(idpb + 1,idpk + 1) += intn[0]*SProt(idpb,1)*SProt_dA(idpk,1) + intn[1]*(SProt(idpb,2)*SProt_dA(idpk,2) + SProt(idpb,3)*SProt_dA(idpk,3));
        }
        for (size_t idc = 1; idc < 7; ++idc) {              //(sp|pp)
          block(idpb + 1,idc + 4) = -SProt_dA(idpb,1)*(intn[2]*PProt(idc,1) + intn[3]*(PProt(idc,2) + PProt(idc,3))) - intn[4]*(SProt_dA(idpb,2)*PProt(idc,4) + SProt_dA(idpb,3)*PProt(idc,5));
          block(idpb + 1,idc + 4) += -SProt(idpb,1)*(intn[2]*PProt_dA(idc,1) + intn[3]*(PProt_dA(idc,2) + PProt_dA(idc,3))) - intn[4]*(SProt(idpb,2)*PProt_dA(idc,4) + SProt(idpb,3)*PProt_dA(idc,5));
        }
      }
      intn[0] = eri2Center_dR(1,0,1,0,1,0,1,0,RCD,D,atmC,atmD,type);                          //(ps,ps|ps,ps)
      intn[1] = eri2Center_dR(1,1,1,1,1,0,1,0,RCD,D,atmC,atmD,type);                          //(pp,pp|ps,ps)
      intn[2] = eri2Center_dR(1,0,1,0,1,1,1,1,RCD,D,atmC,atmD,type);                          //(ps,ps|pp,pp)
      intn[3] = eri2Center_dR(1,1,1,1,1,1,1,1,RCD,D,atmC,atmD,type);                          //(pp,pp|pp,pp)
      intn[4] = eri2Center_dR(1,1,1,0,1,1,1,0,RCD,D,atmC,atmD,type);                          //(pp,ps|pp,ps)
      intn[5] = eri2Center_dR(1,1,1,1,1,-1,1,-1,RCD,D,atmC,atmD,type);                        //(pp,pp|pp*,pp*)
      intn[6] = eri2Center_dR(1,-1,1,1,1,-1,1,1,RCD,D,atmC,atmD,type);                        //(pp*,pp|pp*,pp)
      intn[7] = eri2Center_dR(1,0,1,0,0,0,1,0,RCD,D,atmC,atmD,type);
      intn[8] = eri2Center_dR(1,1,1,1,0,0,1,0,RCD,D,atmC,atmD,type);
      intn[9] = eri2Center_dR(1,0,1,1,0,0,1,1,RCD,D,atmC,atmD,type);
      for (size_t idr = 1; idr < 7; ++idr) {
        for (size_t idc = 1; idc < 4; ++idc) {                //(pp|sp)
          block(idr + 4,idc + 1) = -SProt_dA(idc,1)*(intn[7]*PProt(idr,1) + intn[8]*(PProt(idr,2) + PProt(idr,3))) - intn[9]*(SProt_dA(idc,2)*PProt(idr,4) + SProt_dA(idc,3)*PProt(idr,5));
          block(idr + 4,idc + 1) += -SProt(idc,1)*(intn[7]*PProt_dA(idr,1) + intn[8]*(PProt_dA(idr,2) + PProt_dA(idr,3))) - intn[9]*(SProt(idc,2)*PProt_dA(idr,4) + SProt(idc,3)*PProt_dA(idr,5));
        }
        for (size_t idc = 1; idc < 7; ++idc) {                //(pp|pp)
          block(idr + 4,idc + 4) = intn[0]*PProt_dA(idr,1)*PProt(idc,1) + intn[1]*(PProt_dA(idr,2) + PProt_dA(idr,3))*PProt(idc,1) + intn[2]*PProt_dA(idr,1)*(PProt(idc,2) + PProt(idc,3)) + intn[3]*(PProt_dA(idr,2)*PProt(idc,2) + PProt_dA(idr,3)*PProt(idc,3)) + intn[4]*(PProt_dA(idr,4)*PProt(idc,4) + PProt_dA(idr,5)*PProt(idc,5)) + intn[5]*(PProt_dA(idr,2)*PProt(idc,3) + PProt_dA(idr,3)*PProt(idc,2)) + intn[6]*PProt_dA(idr,6)*PProt(idc,6);
          block(idr + 4,idc + 4) += intn[0]*PProt(idr,1)*PProt_dA(idc,1) + intn[1]*(PProt(idr,2) + PProt(idr,3))*PProt_dA(idc,1) + intn[2]*PProt(idr,1)*(PProt_dA(idc,2) + PProt_dA(idc,3)) + intn[3]*(PProt(idr,2)*PProt_dA(idc,2) + PProt(idr,3)*PProt_dA(idc,3)) + intn[4]*(PProt(idr,4)*PProt_dA(idc,4) + PProt(idr,5)*PProt_dA(idc,5)) + intn[5]*(PProt(idr,2)*PProt_dA(idc,3) + PProt(idr,3)*PProt_dA(idc,2)) + intn[6]*PProt(idr,6)*PProt_dA(idc,6);
        }
      }
      if (nrows > 10) {
        //(sd|sp) integrals
        Dd[0] = Dvalue(atmC,4);
        Dd[1] = Dvalue(atmD,1);
        intn[0] = d_eri2Center_dR(0,0,2,0,0,0,1,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR(0,0,2,1,0,0,1,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        //(sd|sp) rotate
        for (size_t idsd = 1; idsd < 6; ++idsd) {
          for (size_t idsp = 1; idsp < 4; ++idsp) {
            block(10 + idsd,1 + idsp) = -intn[0]*SDrot_dA(idsd,1)*SProt(idsp,1) - intn[1]*(SDrot_dA(idsd,2)*SProt(idsp,2) + SDrot_dA(idsd,3)*SProt(idsp,3));
            block(10 + idsd,1 + idsp) += -intn[0]*SDrot(idsd,1)*SProt_dA(idsp,1) - intn[1]*(SDrot(idsd,2)*SProt_dA(idsp,2) + SDrot(idsd,3)*SProt_dA(idsp,3));
          }
        }
        //(sd|pp) integrals
        Dd[0] = Dvalue(atmC,4);
        Dd[1] = Dvalue(atmD,3);
        intn[0] = d_eri2Center_dR(0,0,2,0,1,0,1,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR(0,0,2,0,1,1,1,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center_dR(0,0,2,1,1,0,1,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[3] = d_eri2Center_dR(0,0,2,2,1,1,1,1,RCD,Dd,atmC,atmD,type);    //(sg,dl|pi,pi)=-(sg,dl|-pi,-pi)=(sg,-dl|pi,-pi)
        //(sd|pp) rotate
        for (size_t idsd = 1; idsd < 6; ++idsd) {
          for (size_t idpp = 1; idpp < 7; ++idpp) {
            block(10 + idsd,4 + idpp) = intn[0]*SDrot_dA(idsd,1)*PProt(idpp,1) + intn[1]*SDrot_dA(idsd,1)*(PProt(idpp,2) + PProt(idpp,3));
            block(10 + idsd,4 + idpp) += intn[0]*SDrot(idsd,1)*PProt_dA(idpp,1) + intn[1]*SDrot(idsd,1)*(PProt_dA(idpp,2) + PProt_dA(idpp,3));
            block(10 + idsd,4 + idpp) += intn[2]*(SDrot_dA(idsd,2)*PProt(idpp,4) + SDrot_dA(idsd,3)*PProt(idpp,5));
            block(10 + idsd,4 + idpp) += intn[2]*(SDrot(idsd,2)*PProt_dA(idpp,4) + SDrot(idsd,3)*PProt_dA(idpp,5));
            block(10 + idsd,4 + idpp) += intn[3]*(SDrot_dA(idsd,4)*PProt(idpp,2) - SDrot_dA(idsd,4)*PProt(idpp,3) + SDrot_dA(idsd,5)*PProt(idpp,6));
            block(10 + idsd,4 + idpp) += intn[3]*(SDrot(idsd,4)*PProt_dA(idpp,2) - SDrot(idsd,4)*PProt_dA(idpp,3) + SDrot(idsd,5)*PProt_dA(idpp,6));
          }
        }
        //(pd|sp) integrals
        Dd[0] = Dvalue(atmC,5);
        Dd[1] = Dvalue(atmD,1);
        intn[0] = d_eri2Center_dR(1,0,2,0,0,0,1,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR(1,1,2,1,0,0,1,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[2] = d_eri2Center_dR(1,0,2,1,0,0,1,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[3] = d_eri2Center_dR(1,1,2,0,0,0,1,1,RCD,Dd,atmC,atmD,type);    //(pi,sg|sg,pi)=(-pi,sg|sg,-pi)
        intn[4] = d_eri2Center_dR(1,1,2,2,0,0,1,1,RCD,Dd,atmC,atmD,type);    //(pi,dl|sg,pi)=(-pi,-dl|sg,pi)=-(-pi,dl|sg,-pi)=(pi,-dl|sg,-pi)
        //(pd|sp) rotate
        for (size_t idsp = 1; idsp < 4; ++idsp) {
          for (size_t idpd1 = 1; idpd1 < 6; ++idpd1) {
            for (size_t idpd2 = 1; idpd2 < 4; ++idpd2) {
              block(5*(idpd2+2) + idpd1,1 + idsp) = intn[0]*PDrot_dA((idpd1 - 1)*3 + idpd2,1)*SProt(idsp,1) + intn[1]*(PDrot_dA((idpd1 - 1)*3 + idpd2,2) + PDrot_dA((idpd1 - 1)*3 + idpd2,3))*SProt(idsp,1);
              block(5*(idpd2+2) + idpd1,1 + idsp) += intn[0]*PDrot((idpd1 - 1)*3 + idpd2,1)*SProt_dA(idsp,1) + intn[1]*(PDrot((idpd1 - 1)*3 + idpd2,2) + PDrot((idpd1 - 1)*3 + idpd2,3))*SProt_dA(idsp,1);
              block(5*(idpd2+2) + idpd1,1 + idsp) += intn[2]*(PDrot_dA((idpd1 - 1)*3 + idpd2,4)*SProt(idsp,2) + PDrot_dA((idpd1 - 1)*3 + idpd2,5)*SProt(idsp,3));
              block(5*(idpd2+2) + idpd1,1 + idsp) += intn[2]*(PDrot((idpd1 - 1)*3 + idpd2,4)*SProt_dA(idsp,2) + PDrot((idpd1 - 1)*3 + idpd2,5)*SProt_dA(idsp,3));
              block(5*(idpd2+2) + idpd1,1 + idsp) += intn[3]*(PDrot_dA((idpd1 - 1)*3 + idpd2,8)*SProt(idsp,2) + PDrot_dA((idpd1 - 1)*3 + idpd2,12)*SProt(idsp,3));
              block(5*(idpd2+2) + idpd1,1 + idsp) += intn[3]*(PDrot((idpd1 - 1)*3 + idpd2,8)*SProt_dA(idsp,2) + PDrot((idpd1 - 1)*3 + idpd2,12)*SProt_dA(idsp,3));
              block(5*(idpd2+2) + idpd1,1 + idsp) += intn[4]*((PDrot_dA((idpd1 - 1)*3 + idpd2,10) + PDrot_dA((idpd1 - 1)*3 + idpd2,15))*SProt(idsp,2) + (PDrot_dA((idpd1 - 1)*3 + idpd2,11) - PDrot_dA((idpd1 - 1)*3 + idpd2,14))*SProt(idsp,3));
              block(5*(idpd2+2) + idpd1,1 + idsp) += intn[4]*((PDrot((idpd1 - 1)*3 + idpd2,10) + PDrot((idpd1 - 1)*3 + idpd2,15))*SProt_dA(idsp,2) + (PDrot((idpd1 - 1)*3 + idpd2,11) - PDrot((idpd1 - 1)*3 + idpd2,14))*SProt_dA(idsp,3));
            }
          }
        }
        //(dd|sp) integrals
        Dd[0] = Dvalue(atmC,6);
        Dd[1] = Dvalue(atmD,1);
        intn[0] = d_eri2Center_dR(2,0,2,0,0,0,1,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR(2,1,2,1,0,0,1,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[2] = d_eri2Center_dR(2,2,2,2,0,0,1,0,RCD,Dd,atmC,atmD,type);    //(dl,dl|sg,sg)=(-dl,-dl|sg,sg)
        intn[3] = d_eri2Center_dR(2,0,2,1,0,0,1,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[4] = d_eri2Center_dR(2,1,2,2,0,0,1,1,RCD,Dd,atmC,atmD,type);    //(pi,dl|sg,pi)=(-pi,-dl|sg,pi)=-(-pi,dl|sg,-pi)=(pi,-dl|sg,-pi)
        //(dd|sp) rotate
        for (size_t idsp = 1; idsp < 4; ++idsp) {
          for (size_t iddd = 1; iddd < 16; ++iddd) {
            block(30 + iddd,1 + idsp) = -intn[0]*DDrot_dA(iddd,1)*SProt(idsp,1) - intn[1]*(DDrot_dA(iddd,2) + DDrot_dA(iddd,3))*SProt(idsp,1) - intn[2]*(DDrot_dA(iddd,4) + DDrot_dA(iddd,5))*SProt(idsp,1);
            block(30 + iddd,1 + idsp) -= intn[0]*DDrot(iddd,1)*SProt_dA(idsp,1) + intn[1]*(DDrot(iddd,2) + DDrot(iddd,3))*SProt_dA(idsp,1) + intn[2]*(DDrot(iddd,4) + DDrot(iddd,5))*SProt_dA(idsp,1);
            block(30 + iddd,1 + idsp) -= intn[3]*(DDrot_dA(iddd,6)*SProt(idsp,2) + DDrot_dA(iddd,7)*SProt(idsp,3));
            block(30 + iddd,1 + idsp) -= intn[3]*(DDrot(iddd,6)*SProt_dA(idsp,2) + DDrot(iddd,7)*SProt_dA(idsp,3));
            block(30 + iddd,1 + idsp) -= intn[4]*((DDrot_dA(iddd,11) + DDrot_dA(iddd,14))*SProt(idsp,2) + (DDrot_dA(iddd,12) - DDrot_dA(iddd,13))*SProt(idsp,3));
            block(30 + iddd,1 + idsp) -= intn[4]*((DDrot(iddd,11) + DDrot(iddd,14))*SProt_dA(idsp,2) + (DDrot(iddd,12) - DDrot(iddd,13))*SProt_dA(idsp,3));
          }
        }
        //(pd|pp) integrals
        Dd[0] = Dvalue(atmC,5);
        Dd[1] = Dvalue(atmD,3);
        intn[0] = d_eri2Center_dR(1,0,2,0,1,0,1,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR(1,1,2,1,1,0,1,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[2] = d_eri2Center_dR(1,0,2,0,1,1,1,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[3] = d_eri2Center_dR(1,1,2,1,1,1,1,1,RCD,Dd,atmC,atmD,type);    //(pi,pi|pi,pi)=(pi,pi|-pi,-pi)=(-pi,-pi|pi,pi)=(-pi,-pi|-pi,-pi)
        intn[4] = d_eri2Center_dR(1,0,2,1,1,0,1,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[5] = d_eri2Center_dR(1,1,2,0,1,1,1,0,RCD,Dd,atmC,atmD,type);    //(pi,sg|pi,sg)=(-pi,sg|-pi,sg)
        intn[6] = d_eri2Center_dR(1,1,2,2,1,0,1,1,RCD,Dd,atmC,atmD,type);    //(pi,dl|sg,pi)=(-pi,-dl|sg,pi)=-(-pi,dl|sg,-pi)=(pi,-dl|sg,-pi)
        //(pd|pp) rotate
        for (size_t idpp = 1; idpp < 7; ++idpp) {
          for (size_t idpd1 = 1; idpd1 < 6; ++idpd1) {
            for (size_t idpd2 = 1; idpd2 < 4; ++idpd2) {
              block(5*(idpd2+2) + idpd1,4 + idpp) = -intn[0]*PDrot_dA((idpd1 - 1)*3 + idpd2,1)*PProt(idpp,1) - intn[1]*(PDrot_dA((idpd1 - 1)*3 + idpd2,2) + PDrot_dA((idpd1 - 1)*3 + idpd2,3))*PProt(idpp,1);
              block(5*(idpd2+2) + idpd1,4 + idpp) += -intn[0]*PDrot((idpd1 - 1)*3 + idpd2,1)*PProt_dA(idpp,1) - intn[1]*(PDrot((idpd1 - 1)*3 + idpd2,2) + PDrot((idpd1 - 1)*3 + idpd2,3))*PProt_dA(idpp,1);
              block(5*(idpd2+2) + idpd1,4 + idpp) -= intn[2]*PDrot_dA((idpd1 - 1)*3 + idpd2,1)*(PProt(idpp,2) + PProt(idpp,3));
              block(5*(idpd2+2) + idpd1,4 + idpp) -= intn[2]*PDrot((idpd1 - 1)*3 + idpd2,1)*(PProt_dA(idpp,2) + PProt_dA(idpp,3));
              block(5*(idpd2+2) + idpd1,4 + idpp) -= intn[3]*(PDrot_dA((idpd1 - 1)*3 + idpd2,2) + PDrot_dA((idpd1 - 1)*3 + idpd2,3))*(PProt(idpp,2) + PProt(idpp,3));
              block(5*(idpd2+2) + idpd1,4 + idpp) -= intn[3]*(PDrot((idpd1 - 1)*3 + idpd2,2) + PDrot((idpd1 - 1)*3 + idpd2,3))*(PProt_dA(idpp,2) + PProt_dA(idpp,3));
              block(5*(idpd2+2) + idpd1,4 + idpp) -= intn[4]*(PDrot_dA((idpd1 - 1)*3 + idpd2,4)*PProt(idpp,4) + PDrot_dA((idpd1 - 1)*3 + idpd2,5)*PProt(idpp,5));
              block(5*(idpd2+2) + idpd1,4 + idpp) -= intn[4]*(PDrot((idpd1 - 1)*3 + idpd2,4)*PProt_dA(idpp,4) + PDrot((idpd1 - 1)*3 + idpd2,5)*PProt_dA(idpp,5));
              block(5*(idpd2+2) + idpd1,4 + idpp) -= intn[5]*(PDrot_dA((idpd1 - 1)*3 + idpd2,8)*PProt(idpp,4) + PDrot_dA((idpd1 - 1)*3 + idpd2,12)*PProt(idpp,5));
              block(5*(idpd2+2) + idpd1,4 + idpp) -= intn[5]*(PDrot((idpd1 - 1)*3 + idpd2,8)*PProt_dA(idpp,4) + PDrot((idpd1 - 1)*3 + idpd2,12)*PProt_dA(idpp,5));
              block(5*(idpd2+2) + idpd1,4 + idpp) -= intn[6]*((PDrot_dA((idpd1 - 1)*3 + idpd2,10) + PDrot_dA((idpd1 - 1)*3 + idpd2,15))*PProt(idpp,4) + (PDrot_dA((idpd1 - 1)*3 + idpd2,11) - PDrot_dA((idpd1 - 1)*3 + idpd2,14))*PProt(idpp,5));
              block(5*(idpd2+2) + idpd1,4 + idpp) -= intn[6]*((PDrot((idpd1 - 1)*3 + idpd2,10) + PDrot((idpd1 - 1)*3 + idpd2,15))*PProt_dA(idpp,4) + (PDrot((idpd1 - 1)*3 + idpd2,11) - PDrot((idpd1 - 1)*3 + idpd2,14))*PProt_dA(idpp,5));
            }
          }
        }
        //(dd|pp) integrals
        Dd[0] = Dvalue(atmC,6);
        Dd[1] = Dvalue(atmD,3);
        intn[0] = d_eri2Center_dR(2,0,2,0,1,0,1,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR(2,1,2,1,1,0,1,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[2] = d_eri2Center_dR(2,2,2,2,1,0,1,0,RCD,Dd,atmC,atmD,type);    //(dl,dl|sg,sg)=(-dl,-dl|sg,sg)
        intn[3] = d_eri2Center_dR(2,0,2,0,1,1,1,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[4] = d_eri2Center_dR(2,1,2,1,1,1,1,1,RCD,Dd,atmC,atmD,type);    //(pi,pi|pi,pi)=(-pi,-pi|-pi,-pi)
        intn[5] = d_eri2Center_dR(2,2,2,2,1,1,1,1,RCD,Dd,atmC,atmD,type);    //(dl,dl|pi,pi)=(-dl,-dl|pi,pi)=(dl,dl|-pi,-pi)=(-dl,-dl|-pi,-pi)
        intn[6] = d_eri2Center_dR(2,0,2,1,1,0,1,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[7] = d_eri2Center_dR(2,-1,2,-1,1,1,1,1,RCD,Dd,atmC,atmD,type);  //(-pi,-pi|pi,pi)=(pi,pi|-pi,-pi)
        intn[8] = d_eri2Center_dR(2,1,2,-1,1,1,1,-1,RCD,Dd,atmC,atmD,type);  //(pi,-pi|pi,-pi)
        intn[9] = d_eri2Center_dR(2,0,2,2,1,1,1,1,RCD,Dd,atmC,atmD,type);    //(sg,dl|pi,pi)=(sg,-dl|pi,-pi)=-(sg,dl|-pi,-pi)
        intn[10] = d_eri2Center_dR(2,1,2,2,1,0,1,1,RCD,Dd,atmC,atmD,type);   //(pi,dl|sg,pi)=(-pi,-dl|sg,pi)=-(-pi,dl|sg,-pi)=(pi,-dl|sg,-pi)
        //(dd|pp) rotate
        for (size_t iddd = 1; iddd < 16; ++iddd) {
          for (size_t idpp = 1; idpp < 7; ++idpp) {
            block(30 + iddd,4 + idpp) = intn[0]*DDrot_dA(iddd,1)*PProt(idpp,1) + intn[1]*(DDrot_dA(iddd,2) + DDrot_dA(iddd,3))*PProt(idpp,1) + intn[2]*(DDrot_dA(iddd,4) + DDrot_dA(iddd,5))*PProt(idpp,1);
            block(30 + iddd,4 + idpp) += intn[0]*DDrot(iddd,1)*PProt_dA(idpp,1) + intn[1]*(DDrot(iddd,2) + DDrot(iddd,3))*PProt_dA(idpp,1) + intn[2]*(DDrot(iddd,4) + DDrot(iddd,5))*PProt_dA(idpp,1);
            block(30 + iddd,4 + idpp) += intn[3]*DDrot_dA(iddd,1)*(PProt(idpp,2) + PProt(idpp,3)) + intn[4]*(DDrot_dA(iddd,2)*PProt(idpp,2) + DDrot_dA(iddd,3)*PProt(idpp,3));
            block(30 + iddd,4 + idpp) += intn[3]*DDrot(iddd,1)*(PProt_dA(idpp,2) + PProt_dA(idpp,3)) + intn[4]*(DDrot(iddd,2)*PProt_dA(idpp,2) + DDrot(iddd,3)*PProt_dA(idpp,3));
            block(30 + iddd,4 + idpp) += intn[5]*(DDrot_dA(iddd,4) + DDrot_dA(iddd,5))*(PProt(idpp,2) + PProt(idpp,3)) + intn[6]*(DDrot_dA(iddd,6)*PProt(idpp,4) + DDrot_dA(iddd,7)*PProt(idpp,5));
            block(30 + iddd,4 + idpp) += intn[5]*(DDrot(iddd,4) + DDrot(iddd,5))*(PProt_dA(idpp,2) + PProt_dA(idpp,3)) + intn[6]*(DDrot(iddd,6)*PProt_dA(idpp,4) + DDrot(iddd,7)*PProt_dA(idpp,5));
            block(30 + iddd,4 + idpp) += intn[7]*(DDrot_dA(iddd,2)*PProt(idpp,3) + DDrot_dA(iddd,3)*PProt(idpp,2)) + intn[8]*DDrot_dA(iddd,10)*PProt(idpp,6);
            block(30 + iddd,4 + idpp) += intn[7]*(DDrot(iddd,2)*PProt_dA(idpp,3) + DDrot(iddd,3)*PProt_dA(idpp,2)) + intn[8]*DDrot(iddd,10)*PProt_dA(idpp,6);
            block(30 + iddd,4 + idpp) += intn[9]*(DDrot_dA(iddd,8)*(PProt(idpp,2) - PProt(idpp,3)) + DDrot_dA(iddd,9)*PProt(idpp,6));
            block(30 + iddd,4 + idpp) += intn[9]*(DDrot(iddd,8)*(PProt_dA(idpp,2) - PProt_dA(idpp,3)) + DDrot(iddd,9)*PProt_dA(idpp,6));
            block(30 + iddd,4 + idpp) += intn[10]*((DDrot_dA(iddd,11) + DDrot_dA(iddd,14))*PProt(idpp,4) + (DDrot_dA(iddd,12) - DDrot_dA(iddd,13))*PProt(idpp,5));
            block(30 + iddd,4 + idpp) += intn[10]*((DDrot(iddd,11) + DDrot(iddd,14))*PProt_dA(idpp,4) + (DDrot(iddd,12) - DDrot(iddd,13))*PProt_dA(idpp,5));
          }
        }
      }
      if (ncols > 10) {
        //(sp|sd) integrals
        Dd[0] = Dvalue(atmC,1);
        Dd[1] = Dvalue(atmD,4);
        intn[0] = d_eri2Center_dR(0,0,1,0,0,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR(0,0,1,1,0,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        //(sp|sd) rotate
        for (size_t idsp = 1; idsp < 4; ++idsp) {
          for (size_t idsd = 1; idsd < 6; ++idsd) {
            block(1 + idsp,10 + idsd) = -intn[0]*SDrot_dA(idsd,1)*SProt(idsp,1) - intn[1]*(SDrot_dA(idsd,2)*SProt(idsp,2) + SDrot_dA(idsd,3)*SProt(idsp,3));
            block(1 + idsp,10 + idsd) -= intn[0]*SDrot(idsd,1)*SProt_dA(idsp,1) + intn[1]*(SDrot(idsd,2)*SProt_dA(idsp,2) + SDrot(idsd,3)*SProt_dA(idsp,3));
          }
        }
        //(pp|sd) integrals
        Dd[0] = Dvalue(atmC,3);
        Dd[1] = Dvalue(atmD,4);
        intn[0] = d_eri2Center_dR(1,0,1,0,0,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR(1,1,1,1,0,0,2,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[2] = d_eri2Center_dR(1,0,1,1,0,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[3] = d_eri2Center_dR(1,1,1,1,0,0,2,2,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,dl)=-(-pi,-pi|sg,dl)=(pi,-pi|sg,-dl)
        //(pp|sd) rotate
        for (size_t idsd = 1; idsd < 6; ++idsd) {
          for (size_t idpp = 1; idpp < 7; ++idpp) {
            block(4 + idpp,10 + idsd) = intn[0]*SDrot_dA(idsd,1)*PProt(idpp,1) + intn[1]*SDrot_dA(idsd,1)*(PProt(idpp,2) + PProt(idpp,3));
            block(4 + idpp,10 + idsd) += intn[0]*SDrot(idsd,1)*PProt_dA(idpp,1) + intn[1]*SDrot(idsd,1)*(PProt_dA(idpp,2) + PProt_dA(idpp,3));
            block(4 + idpp,10 + idsd) += intn[2]*(SDrot_dA(idsd,2)*PProt(idpp,4) + SDrot_dA(idsd,3)*PProt(idpp,5));
            block(4 + idpp,10 + idsd) += intn[2]*(SDrot(idsd,2)*PProt_dA(idpp,4) + SDrot(idsd,3)*PProt_dA(idpp,5));
            block(4 + idpp,10 + idsd) += intn[3]*(SDrot_dA(idsd,4)*PProt(idpp,2) - SDrot_dA(idsd,4)*PProt(idpp,3) + SDrot_dA(idsd,5)*PProt(idpp,6));
            block(4 + idpp,10 + idsd) += intn[3]*(SDrot(idsd,4)*PProt_dA(idpp,2) - SDrot(idsd,4)*PProt_dA(idpp,3) + SDrot(idsd,5)*PProt_dA(idpp,6));
          }
        }
        //(sp|pd) integrals
        Dd[0] = Dvalue(atmC,1);
        Dd[1] = Dvalue(atmD,5);
        intn[0] = d_eri2Center_dR(0,0,1,0,1,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR(0,0,1,0,1,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center_dR(0,0,1,1,1,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[3] = d_eri2Center_dR(0,0,1,1,1,1,2,0,RCD,Dd,atmC,atmD,type);    //(sg,pi|pi,sg)=(sg,-pi|-pi,sg)
        intn[4] = d_eri2Center_dR(0,0,1,1,1,1,2,2,RCD,Dd,atmC,atmD,type);    //(sg,pi|pi,dl)=(sg,pi|-pi,-dl)=-(sg,-pi|-pi,dl)=(sg,-pi|pi,-dl)
        //(sp|pd) rotate
        for (size_t idsp = 1; idsp < 4; ++idsp) {
          for (size_t idpd1 = 1; idpd1 < 6; ++idpd1) {
            for (size_t idpd2 = 1; idpd2 < 4; ++idpd2) {
              block(1 + idsp,5*(idpd2+2) + idpd1) = intn[0]*PDrot_dA((idpd1 - 1)*3 + idpd2,1)*SProt(idsp,1) + intn[1]*(PDrot_dA((idpd1 - 1)*3 + idpd2,2) + PDrot_dA((idpd1 - 1)*3 + idpd2,3))*SProt(idsp,1);
              block(1 + idsp,5*(idpd2+2) + idpd1) += intn[0]*PDrot((idpd1 - 1)*3 + idpd2,1)*SProt_dA(idsp,1) + intn[1]*(PDrot((idpd1 - 1)*3 + idpd2,2) + PDrot((idpd1 - 1)*3 + idpd2,3))*SProt_dA(idsp,1);
              block(1 + idsp,5*(idpd2+2) + idpd1) += intn[2]*(PDrot_dA((idpd1 - 1)*3 + idpd2,4)*SProt(idsp,2) + PDrot_dA((idpd1 - 1)*3 + idpd2,5)*SProt(idsp,3));
              block(1 + idsp,5*(idpd2+2) + idpd1) += intn[2]*(PDrot((idpd1 - 1)*3 + idpd2,4)*SProt_dA(idsp,2) + PDrot((idpd1 - 1)*3 + idpd2,5)*SProt_dA(idsp,3));
              block(1 + idsp,5*(idpd2+2) + idpd1) += intn[3]*(PDrot_dA((idpd1 - 1)*3 + idpd2,8)*SProt(idsp,2) + PDrot_dA((idpd1 - 1)*3 + idpd2,12)*SProt(idsp,3));
              block(1 + idsp,5*(idpd2+2) + idpd1) += intn[3]*(PDrot((idpd1 - 1)*3 + idpd2,8)*SProt_dA(idsp,2) + PDrot((idpd1 - 1)*3 + idpd2,12)*SProt_dA(idsp,3));
              block(1 + idsp,5*(idpd2+2) + idpd1) += intn[4]*((PDrot_dA((idpd1 - 1)*3 + idpd2,10) + PDrot_dA((idpd1 - 1)*3 + idpd2,15))*SProt(idsp,2) + (PDrot_dA((idpd1 - 1)*3 + idpd2,11) - PDrot_dA((idpd1 - 1)*3 + idpd2,14))*SProt(idsp,3));
              block(1 + idsp,5*(idpd2+2) + idpd1) += intn[4]*((PDrot((idpd1 - 1)*3 + idpd2,10) + PDrot((idpd1 - 1)*3 + idpd2,15))*SProt_dA(idsp,2) + (PDrot((idpd1 - 1)*3 + idpd2,11) - PDrot((idpd1 - 1)*3 + idpd2,14))*SProt_dA(idsp,3));
            }
          }
        }
        //(sp|dd) integrals
        Dd[0] = Dvalue(atmC,1);
        Dd[1] = Dvalue(atmD,6);
        intn[0] = d_eri2Center_dR(0,0,1,0,2,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR(0,0,1,0,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center_dR(0,0,1,0,2,2,2,2,RCD,Dd,atmC,atmD,type);    //(sg,sg|dl,dl)=(sg,sg|-dl,-dl)
        intn[3] = d_eri2Center_dR(0,0,1,1,2,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[4] = d_eri2Center_dR(0,0,1,1,2,1,2,2,RCD,Dd,atmC,atmD,type);    //(sg,pi|pi,dl)=(sg,pi|-pi,-dl)=-(sg,-pi|-pi,dl)=(sg,-pi|pi,-dl)
        //(sp|dd) rotate
        for (size_t idsp = 1; idsp < 4; ++idsp) {
          for (size_t iddd = 1; iddd < 16; ++iddd) {
            block(1 + idsp,30 + iddd) = -intn[0]*DDrot_dA(iddd,1)*SProt(idsp,1) - intn[1]*(DDrot_dA(iddd,2) + DDrot_dA(iddd,3))*SProt(idsp,1) - intn[2]*(DDrot_dA(iddd,4) + DDrot_dA(iddd,5))*SProt(idsp,1);
            block(1 + idsp,30 + iddd) -= intn[0]*DDrot(iddd,1)*SProt_dA(idsp,1) + intn[1]*(DDrot(iddd,2) + DDrot(iddd,3))*SProt_dA(idsp,1) + intn[2]*(DDrot(iddd,4) + DDrot(iddd,5))*SProt_dA(idsp,1);
            block(1 + idsp,30 + iddd) -= intn[3]*(DDrot_dA(iddd,6)*SProt(idsp,2) + DDrot_dA(iddd,7)*SProt(idsp,3));
            block(1 + idsp,30 + iddd) -= intn[3]*(DDrot(iddd,6)*SProt_dA(idsp,2) + DDrot(iddd,7)*SProt_dA(idsp,3));
            block(1 + idsp,30 + iddd) -= intn[4]*((DDrot_dA(iddd,11) + DDrot_dA(iddd,14))*SProt(idsp,2) + (DDrot_dA(iddd,12) - DDrot_dA(iddd,13))*SProt(idsp,3));
            block(1 + idsp,30 + iddd) -= intn[4]*((DDrot(iddd,11) + DDrot(iddd,14))*SProt_dA(idsp,2) + (DDrot(iddd,12) - DDrot(iddd,13))*SProt_dA(idsp,3));
          }
        }
        //(pp|pd) integrals
        Dd[0] = Dvalue(atmC,3);
        Dd[1] = Dvalue(atmD,5);
        intn[0] = d_eri2Center_dR(1,0,1,0,1,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR(1,0,1,0,1,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center_dR(1,1,1,1,1,0,2,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[3] = d_eri2Center_dR(1,1,1,1,1,1,2,1,RCD,Dd,atmC,atmD,type);    //(pi,pi|pi,pi)=(pi,pi|-pi,-pi)=(-pi,-pi|pi,pi)=(-pi,-pi|-pi,-pi)
        intn[4] = d_eri2Center_dR(1,0,1,1,1,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[5] = d_eri2Center_dR(1,1,1,0,1,1,2,0,RCD,Dd,atmC,atmD,type);    //(pi,sg|pi,sg)=(-pi,sg|-pi,sg)
        intn[6] = d_eri2Center_dR(1,0,1,1,1,1,2,2,RCD,Dd,atmC,atmD,type);    //(sg,pi|pi,dl)=(sg,pi|-pi,-dl)=-(sg,-pi|-pi,dl)=(sg,-pi|pi,-dl)
        //(pp|pd) rotate
        for (size_t idpp = 1; idpp < 7; ++idpp) {
          for (size_t idpd1 = 1; idpd1 < 6; ++idpd1) {
            for (size_t idpd2 = 1; idpd2 < 4; ++idpd2) {
              block(4 + idpp,5*(idpd2+2) + idpd1) = -intn[0]*PDrot_dA((idpd1 - 1)*3 + idpd2,1)*PProt(idpp,1) - intn[1]*(PDrot_dA((idpd1 - 1)*3 + idpd2,2) + PDrot_dA((idpd1 - 1)*3 + idpd2,3))*PProt(idpp,1);
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[0]*PDrot((idpd1 - 1)*3 + idpd2,1)*PProt_dA(idpp,1) + intn[1]*(PDrot((idpd1 - 1)*3 + idpd2,2) + PDrot((idpd1 - 1)*3 + idpd2,3))*PProt_dA(idpp,1);
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[2]*PDrot_dA((idpd1 - 1)*3 + idpd2,1)*(PProt(idpp,2) + PProt(idpp,3));
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[2]*PDrot((idpd1 - 1)*3 + idpd2,1)*(PProt_dA(idpp,2) + PProt_dA(idpp,3));
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[3]*(PDrot_dA((idpd1 - 1)*3 + idpd2,2) + PDrot_dA((idpd1 - 1)*3 + idpd2,3))*(PProt(idpp,2) + PProt(idpp,3));
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[3]*(PDrot((idpd1 - 1)*3 + idpd2,2) + PDrot((idpd1 - 1)*3 + idpd2,3))*(PProt_dA(idpp,2) + PProt_dA(idpp,3));
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[4]*(PDrot_dA((idpd1 - 1)*3 + idpd2,4)*PProt(idpp,4) + PDrot_dA((idpd1 - 1)*3 + idpd2,5)*PProt(idpp,5));
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[4]*(PDrot((idpd1 - 1)*3 + idpd2,4)*PProt_dA(idpp,4) + PDrot((idpd1 - 1)*3 + idpd2,5)*PProt_dA(idpp,5));
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[5]*(PDrot_dA((idpd1 - 1)*3 + idpd2,8)*PProt(idpp,4) + PDrot_dA((idpd1 - 1)*3 + idpd2,12)*PProt(idpp,5));
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[5]*(PDrot((idpd1 - 1)*3 + idpd2,8)*PProt_dA(idpp,4) + PDrot((idpd1 - 1)*3 + idpd2,12)*PProt_dA(idpp,5));
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[6]*((PDrot_dA((idpd1 - 1)*3 + idpd2,10) + PDrot_dA((idpd1 - 1)*3 + idpd2,15))*PProt(idpp,4) + (PDrot_dA((idpd1 - 1)*3 + idpd2,11) - PDrot_dA((idpd1 - 1)*3 + idpd2,14))*PProt(idpp,5));
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[6]*((PDrot((idpd1 - 1)*3 + idpd2,10) + PDrot((idpd1 - 1)*3 + idpd2,15))*PProt_dA(idpp,4) + (PDrot((idpd1 - 1)*3 + idpd2,11) - PDrot((idpd1 - 1)*3 + idpd2,14))*PProt_dA(idpp,5));
            }
          }
        }
        //(pp|dd) integrals
        Dd[0] = Dvalue(atmC,3);
        Dd[1] = Dvalue(atmD,6);
        intn[0] = d_eri2Center_dR(1,0,1,0,2,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR(1,0,1,0,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center_dR(1,0,1,0,2,2,2,2,RCD,Dd,atmC,atmD,type);    //(sg,sg|dl,dl)=(sg,sg|-dl,-dl)
        intn[3] = d_eri2Center_dR(1,1,1,1,2,0,2,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[4] = d_eri2Center_dR(1,1,1,1,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(pi,pi|pi,pi)=(-pi,-pi|-pi,-pi)
        intn[5] = d_eri2Center_dR(1,1,1,1,2,2,2,2,RCD,Dd,atmC,atmD,type);    //(pi,pi|dl,dl)=(pi,pi|-dl,-dl)=(-pi,-pi|dl,dl)=(-pi,-pi|-dl,-dl)
        intn[6] = d_eri2Center_dR(1,0,1,1,2,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[7] = d_eri2Center_dR(1,-1,1,-1,2,1,2,1,RCD,Dd,atmC,atmD,type);  //(-pi,-pi|pi,pi)=(pi,pi|-pi,-pi)
        intn[8] = d_eri2Center_dR(1,1,1,-1,2,1,2,-1,RCD,Dd,atmC,atmD,type);  //(pi,-pi|pi,-pi)
        intn[9] = d_eri2Center_dR(1,1,1,1,2,0,2,2,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,dl)=(pi,-pi|sg,-dl)=-(-pi,-pi|sg,dl)
        intn[10] = d_eri2Center_dR(1,0,1,1,2,1,2,2,RCD,Dd,atmC,atmD,type);   //(sg,pi|pi,dl)=(sg,pi|-pi,-dl)=-(sg,-pi|-pi,dl)=(sg,-pi|pi,-dl)
        //(pp|dd) rotate
        for (size_t iddd = 1; iddd < 16; ++iddd) {
          for (size_t idpp = 1; idpp < 7; ++idpp) {
            block(4 + idpp,30 + iddd) = intn[0]*DDrot_dA(iddd,1)*PProt(idpp,1) + intn[1]*(DDrot_dA(iddd,2) + DDrot_dA(iddd,3))*PProt(idpp,1) + intn[2]*(DDrot_dA(iddd,4) + DDrot_dA(iddd,5))*PProt(idpp,1);
            block(4 + idpp,30 + iddd) += intn[0]*DDrot(iddd,1)*PProt_dA(idpp,1) + intn[1]*(DDrot(iddd,2) + DDrot(iddd,3))*PProt_dA(idpp,1) + intn[2]*(DDrot(iddd,4) + DDrot(iddd,5))*PProt_dA(idpp,1);
            block(4 + idpp,30 + iddd) += intn[3]*DDrot_dA(iddd,1)*(PProt(idpp,2) + PProt(idpp,3)) + intn[4]*(DDrot_dA(iddd,2)*PProt(idpp,2) + DDrot_dA(iddd,3)*PProt(idpp,3));
            block(4 + idpp,30 + iddd) += intn[3]*DDrot(iddd,1)*(PProt_dA(idpp,2) + PProt_dA(idpp,3)) + intn[4]*(DDrot(iddd,2)*PProt_dA(idpp,2) + DDrot(iddd,3)*PProt_dA(idpp,3));
            block(4 + idpp,30 + iddd) += intn[5]*(DDrot_dA(iddd,4) + DDrot_dA(iddd,5))*(PProt(idpp,2) + PProt(idpp,3)) + intn[6]*(DDrot_dA(iddd,6)*PProt(idpp,4) + DDrot_dA(iddd,7)*PProt(idpp,5));
            block(4 + idpp,30 + iddd) += intn[5]*(DDrot(iddd,4) + DDrot(iddd,5))*(PProt_dA(idpp,2) + PProt_dA(idpp,3)) + intn[6]*(DDrot(iddd,6)*PProt_dA(idpp,4) + DDrot(iddd,7)*PProt_dA(idpp,5));
            block(4 + idpp,30 + iddd) += intn[7]*(DDrot_dA(iddd,2)*PProt(idpp,3) + DDrot_dA(iddd,3)*PProt(idpp,2)) + intn[8]*DDrot_dA(iddd,10)*PProt(idpp,6);
            block(4 + idpp,30 + iddd) += intn[7]*(DDrot(iddd,2)*PProt_dA(idpp,3) + DDrot(iddd,3)*PProt_dA(idpp,2)) + intn[8]*DDrot(iddd,10)*PProt_dA(idpp,6);
            block(4 + idpp,30 + iddd) += intn[9]*(DDrot_dA(iddd,8)*(PProt(idpp,2) - PProt(idpp,3)) + DDrot_dA(iddd,9)*PProt(idpp,6));
            block(4 + idpp,30 + iddd) += intn[9]*(DDrot(iddd,8)*(PProt_dA(idpp,2) - PProt_dA(idpp,3)) + DDrot(iddd,9)*PProt_dA(idpp,6));
            block(4 + idpp,30 + iddd) += intn[10]*((DDrot_dA(iddd,11) + DDrot_dA(iddd,14))*PProt(idpp,4) + (DDrot_dA(iddd,12) - DDrot_dA(iddd,13))*PProt(idpp,5));
            block(4 + idpp,30 + iddd) += intn[10]*((DDrot(iddd,11) + DDrot(iddd,14))*PProt_dA(idpp,4) + (DDrot(iddd,12) - DDrot(iddd,13))*PProt_dA(idpp,5));
          }
        }
      }
      if ((nrows > 10)&&(ncols > 10)) {
        //(sd|sd) integrals
        Dd[0] = Dvalue(atmC,4);
        Dd[1] = Dvalue(atmD,4);
        intn[0] = d_eri2Center_dR(0,0,2,0,0,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR(0,0,2,1,0,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[2] = d_eri2Center_dR(0,0,2,2,0,0,2,2,RCD,Dd,atmC,atmD,type);    //(sg,dl|sg,dl)=(sg,-dl|sg,-dl)
        //(sd|sd) rotate
        for (size_t idsdb = 1; idsdb < 6; ++idsdb) {
          for (size_t idsdk = 1; idsdk < 6; ++idsdk) {
            block(10 + idsdb,10 + idsdk) = intn[0]*SDrot_dA(idsdb,1)*SDrot(idsdk,1) + intn[1]*(SDrot_dA(idsdb,2)*SDrot(idsdk,2) + SDrot_dA(idsdb,3)*SDrot(idsdk,3)) + intn[2]*(SDrot_dA(idsdb,4)*SDrot(idsdk,4) + SDrot_dA(idsdb,5)*SDrot(idsdk,5));
            block(10 + idsdb,10 + idsdk) += intn[0]*SDrot(idsdb,1)*SDrot_dA(idsdk,1) + intn[1]*(SDrot(idsdb,2)*SDrot_dA(idsdk,2) + SDrot(idsdb,3)*SDrot_dA(idsdk,3)) + intn[2]*(SDrot(idsdb,4)*SDrot_dA(idsdk,4) + SDrot(idsdb,5)*SDrot_dA(idsdk,5));
          }
        }
        //(sd|pd) integrals
        Dd[1] = Dvalue(atmD,5);
        intn[0] = d_eri2Center_dR(0,0,2,0,1,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR(0,0,2,0,1,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center_dR(0,0,2,1,1,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[3] = d_eri2Center_dR(0,0,2,1,1,1,2,2,RCD,Dd,atmC,atmD,type);    //(sg,pi|pi,dl)=-(sg,-pi|-pi,dl)=(sg,-pi|pi,-dl)=(sg,pi|-pi,-dl)
        intn[4] = d_eri2Center_dR(0,0,2,1,1,1,2,0,RCD,Dd,atmC,atmD,type);    //(sg,pi|pi,sg)=(sg,-pi|-pi,sg)
        //(sd|pd) rotate
        for (size_t idsd = 1; idsd < 6; ++idsd) {
          for (size_t idpd1 = 1; idpd1 < 6; ++idpd1) {
            for (size_t idpd2 = 1; idpd2 < 4; ++idpd2) {
              block(10 + idsd,5*(idpd2+2) + idpd1) = -intn[0]*PDrot_dA((idpd1 - 1)*3 + idpd2,1)*SDrot(idsd,1) - intn[1]*(PDrot_dA((idpd1 - 1)*3 + idpd2,2) + PDrot_dA((idpd1 - 1)*3 + idpd2,3))*SDrot(idsd,1);
              block(10 + idsd,5*(idpd2+2) + idpd1) -= intn[0]*PDrot((idpd1 - 1)*3 + idpd2,1)*SDrot_dA(idsd,1) + intn[1]*(PDrot((idpd1 - 1)*3 + idpd2,2) + PDrot((idpd1 - 1)*3 + idpd2,3))*SDrot_dA(idsd,1);
              block(10 + idsd,5*(idpd2+2) + idpd1) -= intn[2]*(PDrot_dA((idpd1 - 1)*3 + idpd2,4)*SDrot(idsd,2) + PDrot_dA((idpd1 - 1)*3 + idpd2,5)*SDrot(idsd,3));
              block(10 + idsd,5*(idpd2+2) + idpd1) -= intn[2]*(PDrot((idpd1 - 1)*3 + idpd2,4)*SDrot_dA(idsd,2) + PDrot((idpd1 - 1)*3 + idpd2,5)*SDrot_dA(idsd,3));
              block(10 + idsd,5*(idpd2+2) + idpd1) -= intn[3]*((PDrot_dA((idpd1 - 1)*3 + idpd2,10) + PDrot_dA((idpd1 - 1)*3 + idpd2,15))*SDrot(idsd,2) + (PDrot_dA((idpd1 - 1)*3 + idpd2,11) - PDrot_dA((idpd1 - 1)*3 + idpd2,14))*SDrot(idsd,3));
              block(10 + idsd,5*(idpd2+2) + idpd1) -= intn[3]*((PDrot((idpd1 - 1)*3 + idpd2,10) + PDrot((idpd1 - 1)*3 + idpd2,15))*SDrot_dA(idsd,2) + (PDrot((idpd1 - 1)*3 + idpd2,11) - PDrot((idpd1 - 1)*3 + idpd2,14))*SDrot_dA(idsd,3));
              block(10 + idsd,5*(idpd2+2) + idpd1) -= intn[4]*(PDrot_dA((idpd1 - 1)*3 + idpd2,8)*SDrot(idsd,2) + PDrot_dA((idpd1 - 1)*3 + idpd2,12)*SDrot(idsd,3));
              block(10 + idsd,5*(idpd2+2) + idpd1) -= intn[4]*(PDrot((idpd1 - 1)*3 + idpd2,8)*SDrot_dA(idsd,2) + PDrot((idpd1 - 1)*3 + idpd2,12)*SDrot_dA(idsd,3));
            }
          }
        }
        //(pd|sd) integrals
        Dd[0] = Dvalue(atmC,5);
        Dd[1] = Dvalue(atmD,4);
        intn[0] = d_eri2Center_dR(1,0,2,0,0,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR(1,1,2,1,0,0,2,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[2] = d_eri2Center_dR(1,0,2,1,0,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[3] = d_eri2Center_dR(1,1,2,2,0,0,2,1,RCD,Dd,atmC,atmD,type);    //(pi,dl|sg,pi)=-(-pi,dl|sg,-pi)=(pi,-dl|sg,-pi)=(-pi,-dl|sg,pi)
        intn[4] = d_eri2Center_dR(1,1,2,0,0,0,2,1,RCD,Dd,atmC,atmD,type);    //(pi,sg|sg,pi)=(-pi,sg|sg,-pi)
        //(pd|sd) rotate
        for (size_t idpd1 = 1; idpd1 < 6; ++idpd1) {
          for (size_t idpd2 = 1; idpd2 < 4; ++idpd2) {
            for (size_t idsd = 1; idsd < 6; ++idsd) {
              block(5*(idpd2+2) + idpd1,10 + idsd) = -intn[0]*PDrot_dA((idpd1 - 1)*3 + idpd2,1)*SDrot(idsd,1) - intn[1]*(PDrot_dA((idpd1 - 1)*3 + idpd2,2) + PDrot_dA((idpd1 - 1)*3 + idpd2,3))*SDrot(idsd,1);
              block(5*(idpd2+2) + idpd1,10 + idsd) -= intn[0]*PDrot((idpd1 - 1)*3 + idpd2,1)*SDrot_dA(idsd,1) + intn[1]*(PDrot((idpd1 - 1)*3 + idpd2,2) + PDrot((idpd1 - 1)*3 + idpd2,3))*SDrot_dA(idsd,1);
              block(5*(idpd2+2) + idpd1,10 + idsd) -= intn[2]*(PDrot_dA((idpd1 - 1)*3 + idpd2,4)*SDrot(idsd,2) + PDrot_dA((idpd1 - 1)*3 + idpd2,5)*SDrot(idsd,3));
              block(5*(idpd2+2) + idpd1,10 + idsd) -= intn[2]*(PDrot((idpd1 - 1)*3 + idpd2,4)*SDrot_dA(idsd,2) + PDrot((idpd1 - 1)*3 + idpd2,5)*SDrot_dA(idsd,3));
              block(5*(idpd2+2) + idpd1,10 + idsd) -= intn[3]*((PDrot_dA((idpd1 - 1)*3 + idpd2,10) + PDrot_dA((idpd1 - 1)*3 + idpd2,15))*SDrot(idsd,2) + (PDrot_dA((idpd1 - 1)*3 + idpd2,11) - PDrot_dA((idpd1 - 1)*3 + idpd2,14))*SDrot(idsd,3));
              block(5*(idpd2+2) + idpd1,10 + idsd) -= intn[3]*((PDrot((idpd1 - 1)*3 + idpd2,10) + PDrot((idpd1 - 1)*3 + idpd2,15))*SDrot_dA(idsd,2) + (PDrot((idpd1 - 1)*3 + idpd2,11) - PDrot((idpd1 - 1)*3 + idpd2,14))*SDrot_dA(idsd,3));
              block(5*(idpd2+2) + idpd1,10 + idsd) -= intn[4]*(PDrot_dA((idpd1 - 1)*3 + idpd2,8)*SDrot(idsd,2) + PDrot_dA((idpd1 - 1)*3 + idpd2,12)*SDrot(idsd,3));
              block(5*(idpd2+2) + idpd1,10 + idsd) -= intn[4]*(PDrot((idpd1 - 1)*3 + idpd2,8)*SDrot_dA(idsd,2) + PDrot((idpd1 - 1)*3 + idpd2,12)*SDrot_dA(idsd,3));
            }
          }
        }
        //(dd|sd) integrals
        Dd[0] = Dvalue(atmC,6);
        Dd[1] = Dvalue(atmD,4);
        intn[0] = d_eri2Center_dR(2,0,2,0,0,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR(2,1,2,1,0,0,2,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[2] = d_eri2Center_dR(2,2,2,2,0,0,2,0,RCD,Dd,atmC,atmD,type);    //(dl,dl|sg,sg)=(-dl,-dl|sg,sg)
        intn[3] = d_eri2Center_dR(2,0,2,1,0,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[4] = d_eri2Center_dR(2,0,2,2,0,0,2,2,RCD,Dd,atmC,atmD,type);    //(sg,dl|sg,dl)=(sg,-dl|sg,-dl)
        intn[5] = d_eri2Center_dR(2,1,2,1,0,0,2,2,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,dl)=-(-pi,-pi|sg,dl)=(pi,-pi|sg,-dl)
        intn[6] = d_eri2Center_dR(2,1,2,2,0,0,2,1,RCD,Dd,atmC,atmD,type);    //(pi,dl|sg,pi)=(-pi,-dl|sg,pi)=(pi,-dl|sg,-pi)=-(-pi,dl|sg,-pi)
        //(dd|sd) rotate
        for (size_t iddd = 1; iddd < 16; ++iddd) {
          for (size_t idsd = 1; idsd < 6; ++idsd) {
            block(30 + iddd,10 + idsd) = intn[0]*DDrot_dA(iddd,1)*SDrot(idsd,1) + intn[1]*(DDrot_dA(iddd,2) + DDrot_dA(iddd,3))*SDrot(idsd,1) + intn[2]*(DDrot(iddd,4) + DDrot(iddd,5))*SDrot_dA(idsd,1);
            block(30 + iddd,10 + idsd) += intn[0]*DDrot(iddd,1)*SDrot_dA(idsd,1) + intn[1]*(DDrot(iddd,2) + DDrot(iddd,3))*SDrot_dA(idsd,1) + intn[2]*(DDrot_dA(iddd,4) + DDrot_dA(iddd,5))*SDrot(idsd,1);
            block(30 + iddd,10 + idsd) += intn[3]*(DDrot_dA(iddd,6)*SDrot(idsd,2) + DDrot_dA(iddd,7)*SDrot(idsd,3)) + intn[4]*(DDrot_dA(iddd,8)*SDrot(idsd,4) + DDrot_dA(iddd,9)*SDrot(idsd,5));
            block(30 + iddd,10 + idsd) += intn[3]*(DDrot(iddd,6)*SDrot_dA(idsd,2) + DDrot(iddd,7)*SDrot_dA(idsd,3)) + intn[4]*(DDrot(iddd,8)*SDrot_dA(idsd,4) + DDrot(iddd,9)*SDrot_dA(idsd,5));
            block(30 + iddd,10 + idsd) += intn[5]*((DDrot_dA(iddd,2) - DDrot_dA(iddd,3))*SDrot(idsd,4) + DDrot_dA(iddd,10)*SDrot(idsd,5));
            block(30 + iddd,10 + idsd) += intn[5]*((DDrot(iddd,2) - DDrot(iddd,3))*SDrot_dA(idsd,4) + DDrot(iddd,10)*SDrot_dA(idsd,5));
            block(30 + iddd,10 + idsd) += intn[6]*((DDrot_dA(iddd,11) + DDrot_dA(iddd,14))*SDrot(idsd,2) + (DDrot(iddd,12) - DDrot(iddd,13))*SDrot_dA(idsd,3));
            block(30 + iddd,10 + idsd) += intn[6]*((DDrot(iddd,11) + DDrot(iddd,14))*SDrot_dA(idsd,2) + (DDrot_dA(iddd,12) - DDrot_dA(iddd,13))*SDrot(idsd,3));
          }
        }
        //(sd|dd) integrals
        Dd[0] = Dvalue(atmC,4);
        Dd[1] = Dvalue(atmD,6);
        intn[0] = d_eri2Center_dR(0,0,2,0,2,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR(0,0,2,0,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center_dR(0,0,2,0,2,2,2,2,RCD,Dd,atmC,atmD,type);    //(sg,sg|dl,dl)=(sg,sg|-dl,-dl)
        intn[3] = d_eri2Center_dR(0,0,2,1,2,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[4] = d_eri2Center_dR(0,0,2,2,2,0,2,2,RCD,Dd,atmC,atmD,type);    //(sg,dl|sg,dl)=(sg,-dl|sg,-dl)
        intn[5] = d_eri2Center_dR(0,0,2,2,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,dl|pi,pi)=-(sg,dl|-pi,-pi)=(sg,-dl|pi,-pi)
        intn[6] = d_eri2Center_dR(0,0,2,1,2,1,2,2,RCD,Dd,atmC,atmD,type);    //(sg,pi|pi,dl)=(sg,pi|-pi,-dl)=(sg,-pi|pi,-dl)=-(sg,-pi|-pi,dl)
        //(sd|dd) rotate
        for (size_t idsd = 1; idsd < 6; ++idsd) {
          for (size_t iddd = 1; iddd < 16; ++iddd) {
            block(10 + idsd,30 + iddd) = intn[0]*DDrot_dA(iddd,1)*SDrot(idsd,1) + intn[1]*(DDrot_dA(iddd,2) + DDrot_dA(iddd,3))*SDrot(idsd,1) + intn[2]*(DDrot(iddd,4) + DDrot(iddd,5))*SDrot_dA(idsd,1);
            block(10 + idsd,30 + iddd) += intn[0]*DDrot(iddd,1)*SDrot_dA(idsd,1) + intn[1]*(DDrot(iddd,2) + DDrot(iddd,3))*SDrot_dA(idsd,1) + intn[2]*(DDrot_dA(iddd,4) + DDrot_dA(iddd,5))*SDrot(idsd,1);
            block(10 + idsd,30 + iddd) += intn[3]*(DDrot_dA(iddd,6)*SDrot(idsd,2) + DDrot_dA(iddd,7)*SDrot(idsd,3)) + intn[4]*(DDrot_dA(iddd,8)*SDrot(idsd,4) + DDrot_dA(iddd,9)*SDrot(idsd,5));
            block(10 + idsd,30 + iddd) += intn[3]*(DDrot(iddd,6)*SDrot_dA(idsd,2) + DDrot(iddd,7)*SDrot_dA(idsd,3)) + intn[4]*(DDrot(iddd,8)*SDrot_dA(idsd,4) + DDrot(iddd,9)*SDrot_dA(idsd,5));
            block(10 + idsd,30 + iddd) += intn[5]*((DDrot_dA(iddd,2) - DDrot_dA(iddd,3))*SDrot(idsd,4) + DDrot_dA(iddd,10)*SDrot(idsd,5));
            block(10 + idsd,30 + iddd) += intn[5]*((DDrot(iddd,2) - DDrot(iddd,3))*SDrot_dA(idsd,4) + DDrot(iddd,10)*SDrot_dA(idsd,5));
            block(10 + idsd,30 + iddd) += intn[6]*((DDrot_dA(iddd,11) + DDrot_dA(iddd,14))*SDrot(idsd,2) + (DDrot(iddd,12) - DDrot(iddd,13))*SDrot_dA(idsd,3));
            block(10 + idsd,30 + iddd) += intn[6]*((DDrot(iddd,11) + DDrot(iddd,14))*SDrot_dA(idsd,2) + (DDrot_dA(iddd,12) - DDrot_dA(iddd,13))*SDrot(idsd,3));
          }
        }
        //(pd|pd) integrals
        Dd[0] = Dvalue(atmC,5);
        Dd[1] = Dvalue(atmD,5);
        intn[0] = d_eri2Center_dR(1,0,2,0,1,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR(1,1,2,1,1,0,2,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)=(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center_dR(1,1,2,1,1,1,2,1,RCD,Dd,atmC,atmD,type);    //(pi,pi|pi,pi)=(-pi,-pi|pi,pi)=(pi,pi|-pi,-pi)=(-pi,-pi|-pi,-pi)
        intn[3] = d_eri2Center_dR(1,0,2,1,1,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[4] = d_eri2Center_dR(1,1,2,0,1,1,2,0,RCD,Dd,atmC,atmD,type);    //(pi,sg|pi,sg)=(-pi,sg|-pi,sg)
        intn[5] = d_eri2Center_dR(1,1,2,2,1,1,2,2,RCD,Dd,atmC,atmD,type);    //(pi,dl|pi,dl)=(pi,dl|-pi,-dl)=(-pi,dl|-pi,dl)=-(-pi,dl|pi,-dl)=-(pi,-dl|-pi,dl)=(pi,-dl|pi,-dl)=(-pi,-dl|pi,dl)=(-pi,-dl|-pi,-dl)
        intn[6] = d_eri2Center_dR(1,1,2,0,1,0,2,1,RCD,Dd,atmC,atmD,type);    //(pi,sg|sg,pi)=(-pi,sg|sg,-pi)=(sg,pi|pi,sg)=(sg,-pi|-pi,sg)
        intn[7] = d_eri2Center_dR(1,1,2,2,1,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|pi,dl)=(sg,pi|-pi,-dl)=-(sg,-pi|-pi,dl)=(sg,-pi|pi,-dl)=(pi,dl|sg,pi)=-(-pi,dl|sg,-pi)=(pi,-dl|sg,-pi)=(-pi,-dl|sg,pi)
        intn[8] = d_eri2Center_dR(1,1,2,2,1,1,2,0,RCD,Dd,atmC,atmD,type);    //(pi,sg|pi,dl)=(pi,sg|-pi,-dl)=-(-pi,sg|-pi,dl)=(-pi,sg|pi,-dl)=(pi,dl|pi,sg)=-(-pi,dl|-pi,sg)=(pi,-dl|-pi,sg)=(-pi,-dl|pi,sg)
        //(pd|pd) rotate
        for (size_t idpdb1 = 1; idpdb1 < 6; ++idpdb1) {            //d orbital on bra
          for (size_t idpdb2 = 1; idpdb2 < 4; ++idpdb2) {          //p orbital on bra
            for (size_t idpdk1 = 1; idpdk1 < 6; ++idpdk1) {        //d orbital on ket
              for (size_t idpdk2 = 1; idpdk2 < 4; ++idpdk2) {      //p orbital on ket
                cnt1 = (idpdb1 - 1)*3 + idpdb2;
                cnt2 = (idpdk1 - 1)*3 + idpdk2;
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) = intn[0]*PDrot_dA(cnt1,1)*PDrot(cnt2,1) + intn[1]*((PDrot_dA(cnt1,2) + PDrot_dA(cnt1,3))*PDrot(cnt2,1) + PDrot_dA(cnt1,1)*(PDrot(cnt2,2) + PDrot(cnt2,3)));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[0]*PDrot(cnt1,1)*PDrot_dA(cnt2,1) + intn[1]*((PDrot(cnt1,2) + PDrot(cnt1,3))*PDrot_dA(cnt2,1) + PDrot(cnt1,1)*(PDrot_dA(cnt2,2) + PDrot_dA(cnt2,3)));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[2]*(PDrot_dA(cnt1,2) + PDrot_dA(cnt1,3))*(PDrot(cnt2,2) + PDrot(cnt2,3)) + intn[3]*(PDrot_dA(cnt1,4)*PDrot(cnt2,4) + PDrot_dA(cnt1,5)*PDrot(cnt2,5));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[2]*(PDrot(cnt1,2) + PDrot(cnt1,3))*(PDrot_dA(cnt2,2) + PDrot_dA(cnt2,3)) + intn[3]*(PDrot(cnt1,4)*PDrot_dA(cnt2,4) + PDrot(cnt1,5)*PDrot_dA(cnt2,5));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[4]*(PDrot_dA(cnt1,8)*PDrot(cnt2,8) + PDrot_dA(cnt1,12)*PDrot(cnt2,12));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[4]*(PDrot(cnt1,8)*PDrot_dA(cnt2,8) + PDrot(cnt1,12)*PDrot_dA(cnt2,12));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[5]*((PDrot_dA(cnt1,10) + PDrot_dA(cnt1,15))*(PDrot(cnt2,10) + PDrot(cnt2,15)) + (PDrot_dA(cnt1,14) - PDrot_dA(cnt1,11))*(PDrot(cnt2,14) - PDrot(cnt2,11)));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[5]*((PDrot(cnt1,10) + PDrot(cnt1,15))*(PDrot_dA(cnt2,10) + PDrot_dA(cnt2,15)) + (PDrot(cnt1,14) - PDrot(cnt1,11))*(PDrot_dA(cnt2,14) - PDrot_dA(cnt2,11)));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[6]*(PDrot_dA(cnt1,4)*PDrot(cnt2,8) + PDrot_dA(cnt1,8)*PDrot(cnt2,4) + PDrot_dA(cnt1,5)*PDrot(cnt2,12) + PDrot_dA(cnt1,12)*PDrot(cnt2,5));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[6]*(PDrot(cnt1,4)*PDrot_dA(cnt2,8) + PDrot(cnt1,8)*PDrot_dA(cnt2,4) + PDrot(cnt1,5)*PDrot_dA(cnt2,12) + PDrot(cnt1,12)*PDrot_dA(cnt2,5));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[7]*((PDrot_dA(cnt1,10) + PDrot_dA(cnt1,15))*PDrot(cnt2,4) + PDrot_dA(cnt1,4)*(PDrot(cnt2,10) + PDrot(cnt2,15)) + (PDrot_dA(cnt1,11) - PDrot_dA(cnt1,14))*PDrot(cnt2,5) + PDrot_dA(cnt1,5)*(PDrot(cnt2,11) - PDrot(cnt2,14)));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[7]*((PDrot(cnt1,10) + PDrot(cnt1,15))*PDrot_dA(cnt2,4) + PDrot(cnt1,4)*(PDrot_dA(cnt2,10) + PDrot_dA(cnt2,15)) + (PDrot(cnt1,11) - PDrot(cnt1,14))*PDrot_dA(cnt2,5) + PDrot(cnt1,5)*(PDrot_dA(cnt2,11) - PDrot_dA(cnt2,14)));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[8]*((PDrot_dA(cnt1,10) + PDrot_dA(cnt1,15))*PDrot(cnt2,8) + PDrot_dA(cnt1,8)*(PDrot(cnt2,10) + PDrot(cnt2,15)) + (PDrot_dA(cnt1,11) - PDrot_dA(cnt1,14))*PDrot(cnt2,12) + PDrot_dA(cnt1,12)*(PDrot(cnt2,11) - PDrot(cnt2,14)));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[8]*((PDrot(cnt1,10) + PDrot(cnt1,15))*PDrot_dA(cnt2,8) + PDrot(cnt1,8)*(PDrot_dA(cnt2,10) + PDrot_dA(cnt2,15)) + (PDrot(cnt1,11) - PDrot(cnt1,14))*PDrot_dA(cnt2,12) + PDrot(cnt1,12)*(PDrot_dA(cnt2,11) - PDrot_dA(cnt2,14)));
              }
            }
          }
        }
        //(pd|dd) integrals
        Dd[0] = Dvalue(atmC,5);
        Dd[1] = Dvalue(atmD,6);
        intn[0] = d_eri2Center_dR(1,0,2,0,2,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR(1,0,2,0,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center_dR(1,1,2,1,2,0,2,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[3] = d_eri2Center_dR(1,1,2,1,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(pi,pi|pi,pi)=(-pi,-pi|-pi,-pi)=(pi,pi|-pi,-pi)=(-pi,-pi|pi,pi)
        intn[4] = d_eri2Center_dR(1,0,2,0,2,2,2,2,RCD,Dd,atmC,atmD,type);    //(sg,sg|dl,dl)=(sg,sg|-dl,-dl)
        intn[5] = d_eri2Center_dR(1,1,2,1,2,2,2,2,RCD,Dd,atmC,atmD,type);    //(pi,pi|dl,dl)=(-pi,-pi|-dl,-dl)=(pi,pi|-dl,-dl)=(-pi,-pi|dl,dl)
        intn[6] = d_eri2Center_dR(1,0,2,1,2,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[7] = d_eri2Center_dR(1,1,2,0,2,1,2,0,RCD,Dd,atmC,atmD,type);    //(pi,sg|pi,sg)=(-pi,sg|-pi,sg)
        intn[8] = d_eri2Center_dR(1,1,2,2,2,1,2,2,RCD,Dd,atmC,atmD,type);    //(pi,dl|pi,dl)=(pi,dl|-pi,-dl)=(-pi,-dl|-pi,-dl)=(-pi,-dl|pi,dl)=(-pi,dl|-pi,dl)=-(-pi,dl|pi,-dl)=(pi,-dl|pi,-dl)=-(pi,-dl|-pi,dl)
        intn[9] = d_eri2Center_dR(1,1,2,0,2,1,2,2,RCD,Dd,atmC,atmD,type);    //(pi,sg|pi,dl)=(pi,sg|-pi,-dl)=-(-pi,sg|-pi,dl)=(-pi,sg|pi,-dl)
        intn[10] = d_eri2Center_dR(1,0,2,1,2,1,2,2,RCD,Dd,atmC,atmD,type);   //(sg,pi|pi,dl)=(sg,pi|-pi,-dl)=-(sg,-pi|-pi,dl)=(sg,-pi|pi,-dl)
        intn[11] = d_eri2Center_dR(1,1,2,2,2,0,2,1,RCD,Dd,atmC,atmD,type);   //(pi,dl|sg,pi)=(-pi,-dl|sg,pi)=-(-pi,dl|sg,-pi)=(pi,-dl|sg,-pi)
        //(pd|dd) rotate
        for (size_t iddd = 1; iddd < 16; ++iddd) {
          for (size_t idpd1 = 1; idpd1 < 6; ++idpd1) {
            for (size_t idpd2 = 1; idpd2 < 4; ++idpd2) {
              cnt1 = (idpd1 - 1)*3 + idpd2;
              block(5*(idpd2+2) + idpd1,30 + iddd) = -intn[0]*PDrot_dA(cnt1,1)*DDrot(iddd,1) - intn[1]*PDrot_dA(cnt1,1)*(DDrot(iddd,2) + DDrot(iddd,3));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[0]*PDrot(cnt1,1)*DDrot_dA(iddd,1) + intn[1]*PDrot(cnt1,1)*(DDrot_dA(iddd,2) + DDrot_dA(iddd,3));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[2]*(PDrot_dA(cnt1,2) + PDrot_dA(cnt1,3))*DDrot(iddd,1) + intn[3]*(PDrot_dA(cnt1,2) + PDrot_dA(cnt1,3))*(DDrot(iddd,2) + DDrot(iddd,3));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[2]*(PDrot(cnt1,2) + PDrot(cnt1,3))*DDrot_dA(iddd,1) + intn[3]*(PDrot(cnt1,2) + PDrot(cnt1,3))*(DDrot_dA(iddd,2) + DDrot_dA(iddd,3));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[4]*PDrot_dA(cnt1,1)*(DDrot(iddd,4) + DDrot(iddd,5)) + intn[5]*(PDrot_dA(cnt1,2) + PDrot_dA(cnt1,3))*(DDrot(iddd,4) + DDrot(iddd,5));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[4]*PDrot(cnt1,1)*(DDrot_dA(iddd,4) + DDrot_dA(iddd,5)) + intn[5]*(PDrot(cnt1,2) + PDrot(cnt1,3))*(DDrot_dA(iddd,4) + DDrot_dA(iddd,5));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[6]*(PDrot_dA(cnt1,4)*DDrot(iddd,6) + PDrot_dA(cnt1,5)*DDrot(iddd,7)) + intn[7]*(PDrot_dA(cnt1,8)*DDrot(iddd,6) + PDrot_dA(cnt1,12)*DDrot(iddd,7));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[6]*(PDrot(cnt1,4)*DDrot_dA(iddd,6) + PDrot(cnt1,5)*DDrot_dA(iddd,7)) + intn[7]*(PDrot(cnt1,8)*DDrot_dA(iddd,6) + PDrot(cnt1,12)*DDrot_dA(iddd,7));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[8]*((PDrot_dA(cnt1,10) + PDrot_dA(cnt1,15))*(DDrot(iddd,11) + DDrot(iddd,14)) + (PDrot_dA(cnt1,14) - PDrot_dA(cnt1,11))*(DDrot(iddd,13) - DDrot(iddd,12)));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[8]*((PDrot(cnt1,10) + PDrot(cnt1,15))*(DDrot_dA(iddd,11) + DDrot_dA(iddd,14)) + (PDrot(cnt1,14) - PDrot(cnt1,11))*(DDrot_dA(iddd,13) - DDrot_dA(iddd,12)));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[9]*(PDrot_dA(cnt1,8)*(DDrot(iddd,11) + DDrot(iddd,14)) + PDrot_dA(cnt1,12)*(DDrot(iddd,12) - DDrot(iddd,13)));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[9]*(PDrot(cnt1,8)*(DDrot_dA(iddd,11) + DDrot_dA(iddd,14)) + PDrot(cnt1,12)*(DDrot_dA(iddd,12) - DDrot_dA(iddd,13)));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[10]*(PDrot_dA(cnt1,4)*(DDrot(iddd,11) + DDrot(iddd,14)) + PDrot(cnt1,5)*(DDrot_dA(iddd,12) - DDrot_dA(iddd,13)));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[10]*(PDrot(cnt1,4)*(DDrot_dA(iddd,11) + DDrot_dA(iddd,14)) + PDrot_dA(cnt1,5)*(DDrot(iddd,12) - DDrot(iddd,13)));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[11]*((PDrot_dA(cnt1,10) + PDrot_dA(cnt1,15))*DDrot(iddd,6) + (PDrot(cnt1,11) - PDrot(cnt1,14))*DDrot_dA(iddd,7));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[11]*((PDrot(cnt1,10) + PDrot(cnt1,15))*DDrot_dA(iddd,6) + (PDrot_dA(cnt1,11) - PDrot_dA(cnt1,14))*DDrot(iddd,7));
            }
          }
        }
        //(dd|pd) integrals
        Dd[0] = Dvalue(atmC,6);
        Dd[1] = Dvalue(atmD,5);
        intn[0] = d_eri2Center_dR(2,0,2,0,1,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR(2,1,2,1,1,0,2,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[2] = d_eri2Center_dR(2,0,2,0,1,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[3] = d_eri2Center_dR(2,1,2,1,1,1,2,1,RCD,Dd,atmC,atmD,type);    //(pi,pi|pi,pi)=(-pi,-pi|-pi,-pi)=(pi,pi|-pi,-pi)=(-pi,-pi|pi,pi)
        intn[4] = d_eri2Center_dR(2,2,2,2,1,0,2,0,RCD,Dd,atmC,atmD,type);    //(dl,dl|sg,sg)=(-dl,-dl|sg,sg)
        intn[5] = d_eri2Center_dR(2,2,2,2,1,1,2,1,RCD,Dd,atmC,atmD,type);    //(dl,dl|pi,pi)=(-dl,-dl|-pi,-pi)=(-dl,-dl|pi,pi)=(dl,dl|-pi,-pi)
        intn[6] = d_eri2Center_dR(2,0,2,1,1,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[7] = d_eri2Center_dR(2,1,2,0,1,1,2,0,RCD,Dd,atmC,atmD,type);    //(pi,sg|pi,sg)=(-pi,sg|-pi,sg)
        intn[8] = d_eri2Center_dR(2,1,2,2,1,1,2,2,RCD,Dd,atmC,atmD,type);    //(pi,dl|pi,dl)=(-pi,-dl|-pi,-dl)=(-pi,-dl|pi,dl)=(pi,dl|-pi,-dl)=(-pi,dl|-pi,dl)=(pi,-dl|pi,-dl)=-(-pi,dl|pi,-dl)=-(pi,-dl|-pi,dl)
        intn[9] = d_eri2Center_dR(2,1,2,2,1,1,2,0,RCD,Dd,atmC,atmD,type);    //(pi,dl|pi,sg)=-(-pi,dl|-pi,sg)=(pi,-dl|-pi,sg)=(-pi,-dl|pi,sg)
        intn[10] = d_eri2Center_dR(2,1,2,2,1,0,2,1,RCD,Dd,atmC,atmD,type);   //(pi,dl|sg,pi)=-(-pi,dl|sg,-pi)=(pi,-dl|sg,-pi)=(-pi,-dl|sg,pi)
        intn[11] = d_eri2Center_dR(2,0,2,1,1,1,2,2,RCD,Dd,atmC,atmD,type);   //(sg,pi|pi,dl)=(sg,pi|-pi,-dl)=-(sg,-pi|-pi,dl)=(sg,-pi|pi,-dl)
        //(dd|pd) rotate
        for (size_t iddd = 1; iddd < 16; ++iddd) {
          for (size_t idpd1 = 1; idpd1 < 6; ++idpd1) {
            for (size_t idpd2 = 1; idpd2 < 4; ++idpd2) {
              cnt1 = (idpd1 - 1)*3 + idpd2;
              block(30 + iddd,5*(idpd2+2) + idpd1) = -intn[0]*PDrot_dA(cnt1,1)*DDrot(iddd,1) - intn[1]*PDrot_dA(cnt1,1)*(DDrot(iddd,2) + DDrot(iddd,3));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[0]*PDrot(cnt1,1)*DDrot_dA(iddd,1) + intn[1]*PDrot(cnt1,1)*(DDrot_dA(iddd,2) + DDrot_dA(iddd,3));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[2]*(PDrot_dA(cnt1,2) + PDrot_dA(cnt1,3))*DDrot(iddd,1) + intn[3]*(PDrot(cnt1,2) + PDrot(cnt1,3))*(DDrot_dA(iddd,2) + DDrot_dA(iddd,3));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[2]*(PDrot(cnt1,2) + PDrot(cnt1,3))*DDrot_dA(iddd,1) + intn[3]*(PDrot_dA(cnt1,2) + PDrot_dA(cnt1,3))*(DDrot(iddd,2) + DDrot(iddd,3));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[4]*PDrot_dA(cnt1,1)*(DDrot(iddd,4) + DDrot(iddd,5)) + intn[5]*(PDrot_dA(cnt1,2) + PDrot_dA(cnt1,3))*(DDrot(iddd,4) + DDrot(iddd,5));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[4]*PDrot(cnt1,1)*(DDrot_dA(iddd,4) + DDrot_dA(iddd,5)) + intn[5]*(PDrot(cnt1,2) + PDrot(cnt1,3))*(DDrot_dA(iddd,4) + DDrot_dA(iddd,5));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[6]*(PDrot_dA(cnt1,4)*DDrot(iddd,6) + PDrot_dA(cnt1,5)*DDrot(iddd,7)) + intn[7]*(PDrot_dA(cnt1,8)*DDrot(iddd,6) + PDrot_dA(cnt1,12)*DDrot(iddd,7));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[6]*(PDrot(cnt1,4)*DDrot_dA(iddd,6) + PDrot(cnt1,5)*DDrot_dA(iddd,7)) + intn[7]*(PDrot(cnt1,8)*DDrot_dA(iddd,6) + PDrot(cnt1,12)*DDrot_dA(iddd,7));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[8]*((PDrot_dA(cnt1,10) + PDrot_dA(cnt1,15))*(DDrot(iddd,11) + DDrot(iddd,14)) + (PDrot_dA(cnt1,14) - PDrot_dA(cnt1,11))*(DDrot(iddd,13) - DDrot(iddd,12)));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[8]*((PDrot(cnt1,10) + PDrot(cnt1,15))*(DDrot_dA(iddd,11) + DDrot_dA(iddd,14)) + (PDrot(cnt1,14) - PDrot(cnt1,11))*(DDrot_dA(iddd,13) - DDrot_dA(iddd,12)));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[9]*(PDrot_dA(cnt1,8)*(DDrot(iddd,11) + DDrot(iddd,14)) + PDrot(cnt1,12)*(DDrot_dA(iddd,12) - DDrot_dA(iddd,13)));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[9]*(PDrot(cnt1,8)*(DDrot_dA(iddd,11) + DDrot_dA(iddd,14)) + PDrot_dA(cnt1,12)*(DDrot(iddd,12) - DDrot(iddd,13)));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[10]*(PDrot_dA(cnt1,4)*(DDrot(iddd,11) + DDrot(iddd,14)) + PDrot(cnt1,5)*(DDrot_dA(iddd,12) - DDrot_dA(iddd,13)));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[10]*(PDrot(cnt1,4)*(DDrot_dA(iddd,11) + DDrot_dA(iddd,14)) + PDrot_dA(cnt1,5)*(DDrot(iddd,12) - DDrot(iddd,13)));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[11]*((PDrot_dA(cnt1,10) + PDrot_dA(cnt1,15))*DDrot(iddd,6) + (PDrot(cnt1,11) - PDrot(cnt1,14))*DDrot_dA(iddd,7));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[11]*((PDrot(cnt1,10) + PDrot(cnt1,15))*DDrot_dA(iddd,6) + (PDrot_dA(cnt1,11) - PDrot_dA(cnt1,14))*DDrot(iddd,7));
            }
          }
        }
        //(dd|dd) integrals
        Dd[0] = Dvalue(atmC,6);
        Dd[1] = Dvalue(atmD,6);
        intn[0] = d_eri2Center_dR(2,0,2,0,2,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR(2,0,2,0,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center_dR(2,0,2,0,2,2,2,2,RCD,Dd,atmC,atmD,type);    //(sg,sg|dl,dl)=(sg,sg|-dl,-dl)
        intn[3] = d_eri2Center_dR(2,1,2,1,2,0,2,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[4] = d_eri2Center_dR(2,2,2,2,2,0,2,0,RCD,Dd,atmC,atmD,type);    //(dl,dl|sg,sg)=(-dl,-dl|sg,sg)
        intn[5] = d_eri2Center_dR(2,2,2,2,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(dl,dl|pi,pi)=(-dl,-dl|pi,pi)=(dl,dl|-pi,-pi)=(-dl,-dl|-pi,-pi)
        intn[6] = d_eri2Center_dR(2,1,2,1,2,2,2,2,RCD,Dd,atmC,atmD,type);    //(pi,pi|dl,dl)=(-pi,-pi|dl,dl)=(pi,pi|-dl,-dl)=(-pi,-pi|-dl,-dl)
        intn[7] = d_eri2Center_dR(2,1,2,1,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(pi,pi|pi,pi)=(-pi,-pi|-pi,-pi)
        intn[8] = d_eri2Center_dR(2,2,2,2,2,2,2,2,RCD,Dd,atmC,atmD,type);    //(dl,dl|dl,dl)=(-dl,-dl|-dl,-dl)
        intn[9] = d_eri2Center_dR(2,1,2,1,2,-1,2,-1,RCD,Dd,atmC,atmD,type);  //(pi,pi|-pi,-pi)=(-pi,-pi|pi,pi)
        intn[10] = d_eri2Center_dR(2,-2,2,-2,2,2,2,2,RCD,Dd,atmC,atmD,type); //(-dl,-dl|dl,dl)=(dl,dl|-dl,-dl)
        intn[11] = d_eri2Center_dR(2,0,2,1,2,0,2,1,RCD,Dd,atmC,atmD,type);   //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[12] = d_eri2Center_dR(2,0,2,2,2,0,2,2,RCD,Dd,atmC,atmD,type);   //(sg,dl|sg,dl)=(sg,-dl|sg,-dl)
        intn[13] = d_eri2Center_dR(2,1,2,2,2,1,2,2,RCD,Dd,atmC,atmD,type);   //(pi,dl|pi,dl)=(pi,dl|-pi,-dl)=(-pi,-dl|pi,dl)=(-pi,-dl|-pi,-dl)=-(-pi,dl|pi,-dl)=-(pi,-dl|-pi,dl)=(-pi,dl|-pi,dl)=(pi,-dl|pi,-dl)
        intn[14] = d_eri2Center_dR(2,1,2,-1,2,1,2,-1,RCD,Dd,atmC,atmD,type); //(pi,-pi|pi,-pi)
        intn[15] = d_eri2Center_dR(2,1,2,-2,2,0,2,-1,RCD,Dd,atmC,atmD,type); //(pi,dl|sg,pi)=(-pi,-dl|sg,pi)=(pi,-dl|sg,-pi)=-(-pi,dl|sg,-pi)=(sg,pi|pi,dl)=(sg,pi|-pi,-dl)=(sg,-pi|pi,-dl)=-(sg,-pi|-pi,dl)
        intn[16] = d_eri2Center_dR(2,1,2,1,2,0,2,2,RCD,Dd,atmC,atmD,type);   //(pi,pi|sg,dl)=-(-pi,-pi|sg,dl)=(pi,-pi|sg,-dl)=(sg,dl|pi,pi)=-(sg,dl|-pi,-pi)=(sg,-dl|pi,-pi)
        //(dd|dd) rotate
        for (size_t iddd1 = 1; iddd1 < 16; ++iddd1) {
          for (size_t iddd2 = 1; iddd2 < 16; ++iddd2) {
            block(30 + iddd1,30 + iddd2) = intn[0]*DDrot_dA(iddd1,1)*DDrot(iddd2,1) + intn[1]*DDrot_dA(iddd1,1)*(DDrot(iddd2,2) + DDrot(iddd2,3)) + intn[2]*DDrot(iddd1,1)*(DDrot_dA(iddd2,4) + DDrot_dA(iddd2,5));
            block(30 + iddd1,30 + iddd2) += intn[0]*DDrot(iddd1,1)*DDrot_dA(iddd2,1) + intn[1]*DDrot(iddd1,1)*(DDrot_dA(iddd2,2) + DDrot_dA(iddd2,3)) + intn[2]*DDrot_dA(iddd1,1)*(DDrot(iddd2,4) + DDrot(iddd2,5));
            block(30 + iddd1,30 + iddd2) += intn[3]*(DDrot_dA(iddd1,2) + DDrot_dA(iddd1,3))*DDrot(iddd2,1) + intn[4]*(DDrot(iddd1,4) + DDrot(iddd1,5))*DDrot_dA(iddd2,1);
            block(30 + iddd1,30 + iddd2) += intn[3]*(DDrot(iddd1,2) + DDrot(iddd1,3))*DDrot_dA(iddd2,1) + intn[4]*(DDrot_dA(iddd1,4) + DDrot_dA(iddd1,5))*DDrot(iddd2,1);
            block(30 + iddd1,30 + iddd2) += intn[5]*(DDrot_dA(iddd1,4) + DDrot_dA(iddd1,5))*(DDrot(iddd2,2) + DDrot(iddd2,3)) + intn[6]*(DDrot_dA(iddd1,2) + DDrot_dA(iddd1,3))*(DDrot(iddd2,4) + DDrot(iddd2,5));
            block(30 + iddd1,30 + iddd2) += intn[5]*(DDrot(iddd1,4) + DDrot(iddd1,5))*(DDrot_dA(iddd2,2) + DDrot_dA(iddd2,3)) + intn[6]*(DDrot(iddd1,2) + DDrot(iddd1,3))*(DDrot_dA(iddd2,4) + DDrot_dA(iddd2,5));
            block(30 + iddd1,30 + iddd2) += intn[7]*(DDrot_dA(iddd1,2)*DDrot(iddd2,2) + DDrot_dA(iddd1,3)*DDrot(iddd2,3)) + intn[8]*(DDrot_dA(iddd1,4)*DDrot(iddd2,4) + DDrot_dA(iddd1,5)*DDrot(iddd2,5));
            block(30 + iddd1,30 + iddd2) += intn[7]*(DDrot(iddd1,2)*DDrot_dA(iddd2,2) + DDrot(iddd1,3)*DDrot_dA(iddd2,3)) + intn[8]*(DDrot(iddd1,4)*DDrot_dA(iddd2,4) + DDrot(iddd1,5)*DDrot_dA(iddd2,5));
            block(30 + iddd1,30 + iddd2) += intn[9]*(DDrot_dA(iddd1,2)*DDrot(iddd2,3) + DDrot_dA(iddd1,3)*DDrot(iddd2,2)) + intn[10]*(DDrot_dA(iddd1,4)*DDrot(iddd2,5) + DDrot_dA(iddd1,5)*DDrot(iddd2,4));
            block(30 + iddd1,30 + iddd2) += intn[9]*(DDrot(iddd1,2)*DDrot_dA(iddd2,3) + DDrot(iddd1,3)*DDrot_dA(iddd2,2)) + intn[10]*(DDrot(iddd1,4)*DDrot_dA(iddd2,5) + DDrot(iddd1,5)*DDrot_dA(iddd2,4));
            block(30 + iddd1,30 + iddd2) += intn[11]*(DDrot_dA(iddd1,6)*DDrot(iddd2,6) + DDrot_dA(iddd1,7)*DDrot(iddd2,7)) + intn[12]*(DDrot_dA(iddd1,8)*DDrot(iddd2,8) + DDrot_dA(iddd1,9)*DDrot(iddd2,9));
            block(30 + iddd1,30 + iddd2) += intn[11]*(DDrot(iddd1,6)*DDrot_dA(iddd2,6) + DDrot(iddd1,7)*DDrot_dA(iddd2,7)) + intn[12]*(DDrot(iddd1,8)*DDrot_dA(iddd2,8) + DDrot(iddd1,9)*DDrot_dA(iddd2,9));
            block(30 + iddd1,30 + iddd2) += intn[13]*((DDrot_dA(iddd1,11) + DDrot_dA(iddd1,14))*(DDrot(iddd2,11) + DDrot(iddd2,14)) + (DDrot_dA(iddd1,12) - DDrot_dA(iddd1,13))*(DDrot(iddd2,12) - DDrot(iddd2,13)));
            block(30 + iddd1,30 + iddd2) += intn[13]*((DDrot(iddd1,11) + DDrot(iddd1,14))*(DDrot_dA(iddd2,11) + DDrot_dA(iddd2,14)) + (DDrot(iddd1,12) - DDrot(iddd1,13))*(DDrot_dA(iddd2,12) - DDrot_dA(iddd2,13)));
            block(30 + iddd1,30 + iddd2) += intn[14]*DDrot_dA(iddd1,10)*DDrot(iddd2,10) + intn[14]*DDrot(iddd1,10)*DDrot_dA(iddd2,10);
            block(30 + iddd1,30 + iddd2) += intn[15]*((DDrot_dA(iddd1,11) + DDrot_dA(iddd1,14))*DDrot(iddd2,6) + (DDrot(iddd1,12) - DDrot(iddd1,13))*DDrot_dA(iddd2,7) + DDrot_dA(iddd1,6)*(DDrot(iddd2,11) + DDrot(iddd2,14)) + DDrot(iddd1,7)*(DDrot_dA(iddd2,12) - DDrot_dA(iddd2,13)));
            block(30 + iddd1,30 + iddd2) += intn[15]*((DDrot(iddd1,11) + DDrot(iddd1,14))*DDrot_dA(iddd2,6) + (DDrot_dA(iddd1,12) - DDrot_dA(iddd1,13))*DDrot(iddd2,7) + DDrot(iddd1,6)*(DDrot_dA(iddd2,11) + DDrot_dA(iddd2,14)) + DDrot_dA(iddd1,7)*(DDrot(iddd2,12) - DDrot(iddd2,13)));
            block(30 + iddd1,30 + iddd2) += intn[16]*((DDrot_dA(iddd1,2) - DDrot_dA(iddd1,3))*DDrot(iddd2,8) + DDrot_dA(iddd1,10)*DDrot(iddd2,9) + DDrot_dA(iddd1,8)*(DDrot(iddd2,2) - DDrot(iddd2,3)) + DDrot_dA(iddd1,9)*DDrot(iddd2,10));
            block(30 + iddd1,30 + iddd2) += intn[16]*((DDrot(iddd1,2) - DDrot(iddd1,3))*DDrot_dA(iddd2,8) + DDrot(iddd1,10)*DDrot_dA(iddd2,9) + DDrot(iddd1,8)*(DDrot_dA(iddd2,2) - DDrot_dA(iddd2,3)) + DDrot(iddd1,9)*DDrot_dA(iddd2,10));
          }
        }
      }
    }
  }
  void dFockdX(int der, matrixE & FockA, matrixE & FockB, matrixE & dens, matrixE & qdens, matrixE & hcoreCD, double RCD, double cost, double sint, double cosp, double sinp, int iC, int iD, int atmC, int atmD, int posC, int posD) {
    //function that calculates the first-static-derivatives of the Fock matrix; note that this function cannot be used for the full calculation of the Fock matrix because it only loops over two atoms
    //the derivative variable controls which derivatives are calculated, 1 -> with respect to RCD, 2 -> with respect to thetaCD, 3 -> with respect to phiCD
    //results in atomic units
    size_t natoms = 2;
    int cntA = 0;
    int cntB = 0;
    int im;
    int in;
    int il;
    int is;
    int nAOsC = 1;
    int nAOsD = 1;
    if (atmC > 1) {
      nAOsC = 4;
      if (AtomWithDOrbitals(atmC)) {nAOsC = 9;}
    }
    if (atmD > 1) {
      nAOsD = 4;
      if (AtomWithDOrbitals(atmD)) {nAOsD = 9;}
    }
    chg[0] = double(atmD) - double(CoreCharge[iD]);
    chg[1] = double(atmC) - double(CoreCharge[iC]);
    nintegrals[0] = nAOsC*(nAOsC + 1)/2;
    nintegrals[1] = nAOsD*(nAOsD + 1)/2;
    atomslocal[0] = atmC;
    atomslocal[1] = atmD;
    ipos[0] = posC;
    ipos[1] = posD;
    AOs_local[0] = nAOsC;
    AOs_local[1] = nAOsD;
    std::vector<matrixE> integrals;
    matrixE block(nintegrals[0],1);
    integrals.push_back(block);
    block.resize(nintegrals[1],1);
    integrals.push_back(block);
    block.resize(nintegrals[0],nintegrals[1]);
    bool transitionmetalA;
    bool transitionmetalB;
    double factor = 1.0;
    if (der == 1) {
      IntegralBlock2C_dR(0,block,atmC,atmD,RCD,cost,sint,cosp,sinp);
      IntegralBlock2C_dR(1,integrals[0],atmC,atmD,RCD,cost,sint,cosp,sinp);
      IntegralBlock2C_dR(1,integrals[1],atmD,atmC,RCD,-cost,sint,-cosp,-sinp);
    }
    else if (der > 1) {
      IntegralBlock2C_dA(0,block,der,atmC,atmD,RCD,cost,sint,cosp,sinp);
      IntegralBlock2C_dA(1,integrals[0],der,atmC,atmD,RCD,cost,sint,cosp,sinp);
      IntegralBlock2C_dA(1,integrals[1],der,atmD,atmC,RCD,-cost,sint,-cosp,-sinp);
    }
    //zero Fock matrix and resize
    FockA.resize(nAOsC + nAOsD,nAOsC + nAOsD);
    FockA.zero();
    if (shell == "open") {
      FockB.resize(nAOsC + nAOsD,nAOsC + nAOsD);
      FockB.zero();
    }
    for (size_t idA = 0; idA < natoms; ++idA) {
      transitionmetalA = TransitionMetal(atomslocal[idA]);
      factor = 1.0;
      if ((idA == 1)&&(der == 2)) {factor = -1.0;}
      for (size_t imu = 0; imu < AOs_local[idA]; ++imu) {
        //get position in gamma matrix
        im = imu;
        if (transitionmetalA) {
          if (imu > 4) {im = imu - 5;}
          else if (imu < 5) {im = imu + 4;}
        }
        pos[0] = posgamma(im,im);
        for (size_t inu = 0; inu < AOs_local[idA]; ++inu) {
          in = inu;
          if (transitionmetalA) {
            if (inu > 4) {in = inu - 5;}
            else if (inu < 5) {in = inu + 4;}
          }
          pos[1] = posgamma(im,in);
          FockA(cntA + imu + 1,cntA + inu + 1) -= factor*chg[idA]*integrals[idA](pos[1] + 1,1);
          if (shell == "open") {FockB(cntA + imu + 1,cntA + inu + 1) -= factor*chg[idA]*integrals[idA](pos[1] + 1,1);}
          cntB = 0;
          for (size_t idB = 0; idB < natoms; ++idB) {
            transitionmetalB = TransitionMetal(atomslocal[idB]);
            if (idB != idA) {
              for (size_t ilambda = 0; ilambda < AOs_local[idB]; ++ilambda) {
                il = ilambda;
                if (transitionmetalB) {
                  if (ilambda > 4) {il = ilambda - 5;}
                  else if (ilambda < 5) {il = ilambda + 4;}
                }
                if (imu == inu) {
                  FockA(cntA + imu + 1,cntB + ilambda + 1) += hcoreCD(imu*(idA == 0) + ilambda*(idA == 1) + 1,imu*(idA == 1) + ilambda*(idA == 0) + 1);
                  if (shell == "open") {FockB(cntA + imu + 1,cntB + ilambda + 1) += hcoreCD(imu*(idA == 0) + ilambda*(idA == 1) + 1,imu*(idA == 1) + ilambda*(idA == 0) + 1);}
                }
                for (size_t isigma = 0; isigma < AOs_local[idB]; ++isigma) {
                  is = isigma;
                  if (transitionmetalB) {
                    if (isigma > 4) {is = isigma - 5;}
                    else if (isigma < 5) {is = isigma + 4;}
                  }
                  pos[2] = posgamma(il,is);
                  FockA(cntA + imu + 1,cntB + ilambda + 1) -= 0.5*dens(ipos[idA] + inu + 1,ipos[idB] + isigma + 1)*block(pos[1]*(idA == 0) + pos[2]*(idA == 1) + 1,pos[2]*(idA == 0) + pos[1]*(idA == 1) + 1);
                  if (shell == "open") {
                    FockA(cntA + imu + 1,cntB + ilambda + 1) -= 0.5*qdens(ipos[idA] + inu + 1,ipos[idB] + isigma + 1)*block(pos[1]*(idA == 0) + pos[2]*(idA == 1) + 1,pos[2]*(idA == 0) + pos[1]*(idA == 1) + 1);
                    FockB(cntA + imu + 1,cntB + ilambda + 1) -= 0.5*(dens(ipos[idA] + inu + 1,ipos[idB] + isigma + 1) - qdens(ipos[idA] + inu + 1,ipos[idB] + isigma + 1))*block(pos[1]*(idA == 0) + pos[2]*(idA == 1) + 1,pos[2]*(idA == 0) + pos[1]*(idA == 1) + 1);
                  }
                  FockA(cntA + imu + 1,cntA + inu + 1) += dens(ipos[idB] + ilambda + 1,ipos[idB] + isigma + 1)*block(pos[1]*(idA == 0) + pos[2]*(idA == 1) + 1,pos[2]*(idA == 0) + pos[1]*(idA == 1) + 1);
                  if (shell == "open") {
                    FockB(cntA + imu + 1,cntA + inu + 1) += dens(ipos[idB] + ilambda + 1,ipos[idB] + isigma + 1)*block(pos[1]*(idA == 0) + pos[2]*(idA == 1) + 1,pos[2]*(idA == 0) + pos[1]*(idA == 1) + 1);
                  }
                }
              }
            }
            cntB += AOs_local[idB];
          }
        }
      }
      cntA += AOs_local[idA];
    }
  }
  void dFockdXTEST(int der, matrixE & FockA, matrixE & FockB, matrixE & dens, matrixE & qdens, matrixE & hcoreCD, double RCD, double cost, double sint, double cosp, double sinp, int iC, int iD, int atmC, int atmD, int posC, int posD) {
    size_t natoms = 2;
    int cntA = 0;
    int cntB = 0;
    int im;
    int in;
    int il;
    int is;
    int nAOsC = 1;
    int nAOsD = 1;
    if (atmC > 1) {
      nAOsC = 4;
      if (AtomWithDOrbitals(atmC)) {nAOsC = 9;}
    }
    if (atmD > 1) {
      nAOsD = 4;
      if (AtomWithDOrbitals(atmD)) {nAOsD = 9;}
    }
    chg[0] = double(atmD) - double(CoreCharge[iD]);
    chg[1] = double(atmC) - double(CoreCharge[iC]);
    nintegrals[0] = nAOsC*(nAOsC + 1)/2;
    nintegrals[1] = nAOsD*(nAOsD + 1)/2;
    atomslocal[0] = atmC;
    atomslocal[1] = atmD;
    ipos[0] = posC;
    ipos[1] = posD;
    AOs_local[0] = nAOsC;
    AOs_local[1] = nAOsD;
    std::vector<matrixE> integrals;
    matrixE block(nintegrals[0],1);
    integrals.push_back(block);
    block.resize(nintegrals[1],1);
    integrals.push_back(block);
    block.resize(nintegrals[0],nintegrals[1]);
    bool transitionmetalA;
    bool transitionmetalB;
    double factor = 1.0;
    IntegralBlock2C(0,block,atmC,atmD,RCD,cost,sint,cosp,sinp);
    IntegralBlock2C(1,integrals[0],atmC,atmD,RCD,cost,sint,cosp,sinp);
    IntegralBlock2C(1,integrals[1],atmD,atmC,RCD,-cost,sint,-cosp,-sinp);
    //zero Fock matrix and resize
    FockA.resize(nAOsC + nAOsD,nAOsC + nAOsD);
    FockA.zero();
    if (shell == "open") {
      FockB.resize(nAOsC + nAOsD,nAOsC + nAOsD);
      FockB.zero();
    }
    for (size_t idA = 0; idA < natoms; ++idA) {
      transitionmetalA = TransitionMetal(atomslocal[idA]);
      factor = 1.0;
      if ((idA == 1)&&(der == 2)) {factor = -1.0;}
      for (size_t imu = 0; imu < AOs_local[idA]; ++imu) {
        //get position in gamma matrix
        im = imu;
        if (transitionmetalA) {
          if (imu > 4) {im = imu - 5;}
          else if (imu < 5) {im = imu + 4;}
        }
        pos[0] = posgamma(im,im);
        for (size_t inu = 0; inu < AOs_local[idA]; ++inu) {
          in = inu;
          if (transitionmetalA) {
            if (inu > 4) {in = inu - 5;}
            else if (inu < 5) {in = inu + 4;}
          }
          pos[1] = posgamma(im,in);
          FockA(cntA + imu + 1,cntA + inu + 1) -= factor*chg[idA]*integrals[idA](pos[1] + 1,1);
          if (shell == "open") {FockB(cntA + imu + 1,cntA + inu + 1) -= factor*chg[idA]*integrals[idA](pos[1] + 1,1);}
          cntB = 0;
          for (size_t idB = 0; idB < natoms; ++idB) {
            transitionmetalB = TransitionMetal(atomslocal[idB]);
            if (idB != idA) {
              for (size_t ilambda = 0; ilambda < AOs_local[idB]; ++ilambda) {
                il = ilambda;
                if (transitionmetalB) {
                  if (ilambda > 4) {il = ilambda - 5;}
                  else if (ilambda < 5) {il = ilambda + 4;}
                }
                if (imu == inu) {
                  FockA(cntA + imu + 1,cntB + ilambda + 1) += hcoreCD(imu*(idA == 0) + ilambda*(idA == 1) + 1,imu*(idA == 1) + ilambda*(idA == 0) + 1);
                  if (shell == "open") {FockB(cntA + imu + 1,cntB + ilambda + 1) += hcoreCD(imu*(idA == 0) + ilambda*(idA == 1) + 1,imu*(idA == 1) + ilambda*(idA == 0) + 1);}
                }
                for (size_t isigma = 0; isigma < AOs_local[idB]; ++isigma) {
                  is = isigma;
                  if (transitionmetalB) {
                    if (isigma > 4) {is = isigma - 5;}
                    else if (isigma < 5) {is = isigma + 4;}
                  }
                  pos[2] = posgamma(il,is);
                  FockA(cntA + imu + 1,cntB + ilambda + 1) -= 0.5*dens(ipos[idA] + inu + 1,ipos[idB] + isigma + 1)*block(pos[1]*(idA == 0) + pos[2]*(idA == 1) + 1,pos[2]*(idA == 0) + pos[1]*(idA == 1) + 1);
                  if (shell == "open") {
                    FockA(cntA + imu + 1,cntB + ilambda + 1) -= 0.5*qdens(ipos[idA] + inu + 1,ipos[idB] + isigma + 1)*block(pos[1]*(idA == 0) + pos[2]*(idA == 1) + 1,pos[2]*(idA == 0) + pos[1]*(idA == 1) + 1);
                    FockB(cntA + imu + 1,cntB + ilambda + 1) -= 0.5*(dens(ipos[idA] + inu + 1,ipos[idB] + isigma + 1) - qdens(ipos[idA] + inu + 1,ipos[idB] + isigma + 1))*block(pos[1]*(idA == 0) + pos[2]*(idA == 1) + 1,pos[2]*(idA == 0) + pos[1]*(idA == 1) + 1);
                  }
                  FockA(cntA + imu + 1,cntA + inu + 1) += dens(ipos[idB] + ilambda + 1,ipos[idB] + isigma + 1)*block(pos[1]*(idA == 0) + pos[2]*(idA == 1) + 1,pos[2]*(idA == 0) + pos[1]*(idA == 1) + 1);
                  if (shell == "open") {
                    FockB(cntA + imu + 1,cntA + inu + 1) += dens(ipos[idB] + ilambda + 1,ipos[idB] + isigma + 1)*block(pos[1]*(idA == 0) + pos[2]*(idA == 1) + 1,pos[2]*(idA == 0) + pos[1]*(idA == 1) + 1);
                  }
                }
              }
            }
            cntB += AOs_local[idB];
          }
        }
      }
      cntA += AOs_local[idA];
    }
  }
  int posgamma(int mu1, int nu1) {
    //function returning the gamma matrix element for the orbital pair (mu,nu)
    int mu = std::min(mu1,nu1);
    int nu = std::max(mu1,nu1);
    int pos = 0;
    if (mu == 0) {
      if (nu == 1) {pos = 1;}       //s,px
      else if (nu == 2) {pos = 2;}  //s,py
      else if (nu == 3) {pos = 3;}  //s,pz
      else if (nu == 4) {pos = 10;} //s,dxz
      else if (nu == 5) {pos = 11;} //s,dyz
      else if (nu == 6) {pos = 12;} //s,dz2
      else if (nu == 7) {pos = 13;} //s,dxy
      else if (nu == 8) {pos = 14;} //s,dx2-y2
    }
    else if (mu == 1) {
      if (nu == 1) {pos = 4;}       //px,px
      else if (nu == 2) {pos = 5;}  //px,py
      else if (nu == 3) {pos = 6;}  //px,pz
      else if (nu == 4) {pos = 15;} //px,dxz
      else if (nu == 5) {pos = 16;} //px,dyz
      else if (nu == 6) {pos = 17;} //px,dz2
      else if (nu == 7) {pos = 18;} //px,dxy
      else if (nu == 8) {pos = 19;} //px,dx2-y2
    }
    else if (mu == 2) {
      if (nu == 2) {pos = 7;}       //py,py
      else if (nu == 3) {pos = 8;}  //py,pz
      else if (nu == 4) {pos = 20;} //py,dxz
      else if (nu == 5) {pos = 21;} //py,dyz
      else if (nu == 6) {pos = 22;} //py,dz2
      else if (nu == 7) {pos = 23;} //py,dxy
      else if (nu == 8) {pos = 24;} //py,dx2-y2
    }
    else if (mu == 3) {
      if (nu == 3) {pos = 9;}       //pz,pz
      else if (nu == 4) {pos = 25;} //pz,dxz
      else if (nu == 5) {pos = 26;} //pz,dyz
      else if (nu == 6) {pos = 27;} //pz,dz2
      else if (nu == 7) {pos = 28;} //pz,dxy
      else if (nu == 8) {pos = 29;} //pz,dx2-y2
    }
    else if (mu == 4) {
      if (nu == 4) {pos = 30;}      //dxz,dxz
      else if (nu == 5) {pos = 31;} //dxz,dyz
      else if (nu == 6) {pos = 32;} //dxz,dz2
      else if (nu == 7) {pos = 33;} //dxz,dxy
      else if (nu == 8) {pos = 34;} //dxz,dx2-y2
    }
    else if (mu == 5) {
      if (nu == 5) {pos = 35;}      //dyz,dyz
      else if (nu == 6) {pos = 36;} //dyz,dz2
      else if (nu == 7) {pos = 37;} //dyz,dxy
      else if (nu == 8) {pos = 38;} //dyz,dx2-y2
    }
    else if (mu == 6) {
      if (nu == 6) {pos = 39;}      //dz2,dz2
      else if (nu == 7) {pos = 40;} //dz2,dxy
      else if (nu == 8) {pos = 41;} //dz2,dx2-y2
    }
    else if (mu == 7) {
      if (nu == 7) {pos = 42;}      //dxy,dxy
      else if (nu == 8) {pos = 43;} //dxy,dx2-y2
    }
    else if (mu == 8) {
      if (nu == 8) {pos = 44;}      //dx2-y2,dx2-y2
    }
    return pos;
  }
  virtual double gfactor(size_t iatm1, size_t iatm2, double RAB) {
    //function determining special factors in the nuclear repulsion terms of MNDO
    double gfac = 1.0;
    if (iatm2 == 1) {
      if ((iatm1 == 7)||(iatm1 == 8)) {gfac = RAB*dist_Angstrom2au;}
    }
    return gfac;
  }
  virtual double gfactor_dR(size_t iatm1, size_t iatm2) {
    //function determining the derivative of the special factors in the nuclear repulsion terms of MNDO
    //with respect to internuclear distance
    double gfac = 0.0;
    if (iatm2 == 1) {
      if ((iatm1 == 7)||(iatm1 == 8)) {gfac = dist_Angstrom2au;}
    }
    return gfac;
  }
  virtual double AM1factor(double RAB, std::vector<double> & K1, std::vector<double> & L1, std::vector<double> & M1) {return 0.0;}
  virtual double AM1factor_dR(double RAB, std::vector<double> & K1, std::vector<double> & L1, std::vector<double> & M1) {return 0.0;}
  virtual double AM1factor_dR2(double RAB, std::vector<double> & K1, std::vector<double> & L1, std::vector<double> & M1) {return 0.0;}
  double EPDDG(int atmA, int atmB, double RAB) {
    //function calculating the PDDG contribution to the nuclear energy
    double nA = ValenceElectrons(atmA);
    double nB = ValenceElectrons(atmB);
    double PAi;
    double DAi;
    double epddg = 0.0;
    double aux;
    for (size_t idi = 1; idi < 3; ++idi) {
      PAi = PAPDDG(atmA,idi);
      DAi = DAPDDG(atmA,idi);
      for (size_t idj = 1; idj < 3; ++idj) {
        aux = RAB - DAi - DAPDDG(atmB,idj);
        epddg += (PAi*nA + PAPDDG(atmB,idj)*nB)*exp(-10.0*aux*aux);
      }
    }
    return epddg/(nA + nB);
  }
  double dEPDDG(int atmA, int atmB, double RAB) {
    //function calculating the PDDG contribution to the nuclear energy
    double nA = ValenceElectrons(atmA);
    double nB = ValenceElectrons(atmB);
    double PAi;
    double DAi;
    double depddg = 0.0;
    double aux;
    for (size_t idi = 1; idi < 3; ++idi) {
      PAi = PAPDDG(atmA,idi);
      DAi = DAPDDG(atmA,idi);
      for (size_t idj = 1; idj < 3; ++idj) {
        aux = RAB - DAi - DAPDDG(atmB,idj);
        depddg += (PAi*nA + PAPDDG(atmB,idj)*nB)*exp(-10.0*aux*aux)*aux;
      }
    }
    return -20.0*depddg/(nA + nB);
  }
  double d2EPDDG(int atmA, int atmB, double RAB) {
    //function calculating the PDDG contribution to the nuclear energy
    double nA = ValenceElectrons(atmA);
    double nB = ValenceElectrons(atmB);
    double PAi;
    double DAi;
    double d2epddg = 0.0;
    double aux;
    for (size_t idi = 1; idi < 3; ++idi) {
      PAi = PAPDDG(atmA,idi);
      DAi = DAPDDG(atmA,idi);
      for (size_t idj = 1; idj < 3; ++idj) {
        aux = RAB - DAi - DAPDDG(atmB,idj);
        d2epddg += (PAi*nA + PAPDDG(atmB,idj)*nB)*exp(-10.0*aux*aux)*(1.0 - 20.0*aux*aux);
      }
    }
    return -20.0*d2epddg/(nA + nB);
  }
  matrixE EnuclearMat() {
    //function that calculates the nuclear energy
    double XAB;
    double alphaA;
    double chgA;
    double alphaB;
    double chgB;
    double RAB;
    double integral;
    double aux;
    double enucaux = 0.0;
    double gA;
    double gB;
    double factorA;
    double factorB;
    double factorC;
    double unitconv = 1.0/(au2eV*dist_Angstrom2au);
    matrixE Enuclear(Natoms,Natoms);
    for (size_t iatm = 0; iatm < Natoms; ++iatm) {
      //charge and nuclear energy variables
      KA = AM1K(atoms[iatm]);
      LA = AM1L(atoms[iatm]);
      MA = AM1M(atoms[iatm]);
      chgA = double(int(atoms[iatm]) - int(CoreCharge[iatm]));
      Enuclear(iatm + 1,iatm + 1) = 0.0;
      for (size_t ibtm = iatm + 1; ibtm < Natoms; ++ibtm) {
        alphaA = alpha(atoms[iatm],atoms[ibtm]);
        XAB = xAB(atoms[iatm],atoms[ibtm]);
        //charge and nuclear energy variables
        KB = AM1K(atoms[ibtm],atoms[iatm]);
        LB = AM1L(atoms[ibtm],atoms[iatm]);
        MB = AM1M(atoms[ibtm],atoms[iatm]);
        if (atoms[iatm] == 5) {
          KA = AM1K(atoms[iatm],atoms[ibtm]);
          LA = AM1L(atoms[iatm],atoms[ibtm]);
          MA = AM1M(atoms[iatm],atoms[ibtm]);
        }
        chgB = double(int(atoms[ibtm]) - int(CoreCharge[ibtm]));
        //interatomic distance
        RAB = mol.AUdistance(iatm + 1,ibtm + 1);
        //integral
        integral = eri2Center(0,0,0,0,0,0,0,0,RAB,D,atoms[iatm],atoms[ibtm],2);
        factorA = 1.0;
        factorB = 0.0;
        factorC = 0.0;
        if (XAB > 0.0) {           //PM6
          //assume here that these parameters take values considerably different from zero
          if (((atoms[iatm] == 1)&&(atoms[ibtm] == 6))||((atoms[iatm] == 6)&&(atoms[ibtm] == 1))||((atoms[iatm] == 1)&&(atoms[ibtm] == 7))||((atoms[iatm] == 7)&&(atoms[ibtm] == 1))||((atoms[iatm] == 1)&&(atoms[ibtm] == 8))||((atoms[iatm] == 8)&&(atoms[ibtm] == 1))) {
            factorA += 2.0*XAB*exp(-alphaA*RAB*RAB*dist_Angstrom2au);
          }
          else {
            factorA += 2.0*XAB*exp(-alphaA*RAB*(1.0 + 0.0003*RAB*RAB*RAB*RAB*RAB*dist_Angstrom2au*dist_Angstrom2au*dist_Angstrom2au*dist_Angstrom2au*dist_Angstrom2au));
          }
          if ((atoms[iatm] == 6)&&(atoms[ibtm] == 6)) {factorA += 9.28*exp(-5.98*RAB*dist_Angstrom2au);}
          else if (((atoms[iatm] == 8)&&(atoms[ibtm] == 14))||((atoms[iatm] == 14)&&(atoms[ibtm] == 8))) {
            factorA -= 0.0007*exp(-(RAB*dist_Angstrom2au - 2.9)*(RAB*dist_Angstrom2au - 2.9));
          }
          factorB = 1.0e-8*pow((pow(chgA,1.0/3.0) + pow(chgB,1.0/3.0))/(RAB*dist_Angstrom2au),12.0)/au2eV;
        }
        gA = gfactor(atoms[iatm],atoms[ibtm],RAB);
        gB = gfactor(atoms[ibtm],atoms[iatm],RAB);
        alphaB = alpha(atoms[ibtm],atoms[iatm]);
        factorA += gA*exp(-alphaA*RAB*dist_Angstrom2au) + gB*exp(-alphaB*RAB*dist_Angstrom2au);
        if (fabs(KA[0]) > 1.0e-5) {factorC += AM1factor(RAB,KA,LA,MA);}
        if (fabs(KB[0]) > 1.0e-5) {factorC += AM1factor(RAB,KB,LB,MB);}
        aux = chgA*chgB*integral*factorA + factorB + chgA*chgB*factorC;
        if (doPDDG) {aux += EPDDG(atoms[iatm],atoms[ibtm],RAB*dist_Angstrom2au);}
        enucaux += aux;
        Enuclear(iatm + 1,ibtm + 1) = aux;
        Enuclear(ibtm + 1,iatm + 1) = aux;
      }
    }
    return Enuclear;
  }
  void Enuclear_dR(std::vector<double> & enuclear_dr) {
    //function that calculates the derivative of the nuclear energy with respect to internuclear distance
    size_t index = 0;
    enuclear_dr.resize(Natoms*(Natoms - 1)/2);
    for (size_t idx = 0; idx < Natoms*(Natoms - 1)/2; ++idx) {
      enuclear_dr[idx] = 0.0;
    }
    double alphaA;
    double chgA;
    double alphaB;
    double chgB;
    double XAB;
    double RAB;
    double RABangstroem;
    double integral;
    double integral_dR;
    double aux;
    double gA;
    double gB;
    double gAdR;
    double gBdR;
    double expA;
    double expB;
    double fac1;
    double fac2;
    double fac3;
    double fac4;
    double A2au5 = dist_Angstrom2au*dist_Angstrom2au*dist_Angstrom2au*dist_Angstrom2au*dist_Angstrom2au;
    double RAB5;
    double unitconv = 1.0/(au2eV*dist_Angstrom2au);
    for (size_t iatm = 0; iatm < Natoms; ++iatm) {
      //charge and nuclear energy variables
      KA = AM1K(atoms[iatm]);
      LA = AM1L(atoms[iatm]);
      MA = AM1M(atoms[iatm]);
      chgA = double(int(atoms[iatm]) - int(CoreCharge[iatm]));
      for (size_t ibtm = iatm + 1; ibtm < Natoms; ++ibtm) {
        alphaA = alpha(atoms[iatm],atoms[ibtm]);
        alphaB = alpha(atoms[ibtm],atoms[iatm]);
        //charge and nuclear energy variables
        XAB = xAB(atoms[iatm],atoms[ibtm]);
        KB = AM1K(atoms[ibtm],atoms[iatm]);
        LB = AM1L(atoms[ibtm],atoms[iatm]);
        MB = AM1M(atoms[ibtm],atoms[iatm]);
        if (atoms[iatm] == 5) {
          KA = AM1K(atoms[iatm],atoms[ibtm]);
          LA = AM1L(atoms[iatm],atoms[ibtm]);
          MA = AM1M(atoms[iatm],atoms[ibtm]);
        }
        chgB = double(int(atoms[ibtm]) - int(CoreCharge[ibtm]));
        //interatomic distance
        RAB = mol.AUdistance(iatm + 1,ibtm + 1);
        RABangstroem = RAB*dist_Angstrom2au;
        //integrals required
        integral = eri2Center(0,0,0,0,0,0,0,0,RAB,D,atoms[iatm],atoms[ibtm],2);
        integral_dR = eri2Center_dR(0,0,0,0,0,0,0,0,RAB,D,atoms[iatm],atoms[ibtm],2);
        //auxiliaries
        fac1 = 1.0;
        fac2 = 0.0;
        fac3 = 0.0;
        fac4 = 0.0;
        if (XAB > 0.0) {           //PM6
          if (((atoms[iatm] == 1)&&(atoms[ibtm] == 6))||((atoms[iatm] == 6)&&(atoms[ibtm] == 1))||((atoms[iatm] == 1)&&(atoms[ibtm] == 7))||((atoms[iatm] == 7)&&(atoms[ibtm] == 1))||((atoms[iatm] == 1)&&(atoms[ibtm] == 8))||((atoms[iatm] == 8)&&(atoms[ibtm] == 1))) {
            expA = 2.0*exp(-alphaA*RAB*RABangstroem);
            fac1 += XAB*expA;
            fac2 += 2.0*XAB*alphaA*expA*RABangstroem;
          }
          else {
            RAB5 = RAB*RAB*RAB*RAB*RAB*A2au5;
            expA = 2.0*exp(-alphaA*RAB*(1.0 + 0.0003*RAB5));
            fac1 += XAB*expA;
            fac2 += XAB*alphaA*expA*(1.0 + 0.0018*RAB5);
          }
          if ((atoms[iatm] == 6)&&(atoms[ibtm] == 6)) {
            expA = exp(-5.98*RABangstroem);
            fac1 += 9.28*expA;
            fac2 += 55.4944*expA*dist_Angstrom2au;
          }
          else if (((atoms[iatm] == 8)&&(atoms[ibtm] == 14))||((atoms[iatm] == 14)&&(atoms[ibtm] == 8))) {
            aux = RABangstroem - 2.9;
            expA = exp(-aux*aux);
            fac1 -= 0.0007*expA;
            fac2 -= 0.0014*aux*expA*dist_Angstrom2au;
          }
          fac3 = 1.2e-7*pow((pow(chgA,1.0/3.0) + pow(chgB,1.0/3.0))/(RABangstroem),12.0)/(au2eV*RAB);
        }
        gA = gfactor(atoms[iatm],atoms[ibtm],RAB);
        gB = gfactor(atoms[ibtm],atoms[iatm],RAB);
        if ((gA > 0.0)||(gB > 0.0)) {
          expA = exp(-alphaA*RABangstroem);
          expB = exp(-alphaB*RABangstroem);
          gAdR = gfactor_dR(atoms[iatm],atoms[ibtm]);
          gBdR = gfactor_dR(atoms[ibtm],atoms[iatm]);
          fac1 += gA*expA + gB*expB;
          fac2 -= gAdR*expA - gA*alphaA*dist_Angstrom2au*expA + gBdR*expB - gB*alphaB*dist_Angstrom2au*expB;
        }
        if (fabs(KA[0]) > 1.0e-5) {fac4 += AM1factor_dR(RAB,KA,LA,MA);}
        if (fabs(KB[0]) > 1.0e-5) {fac4 += AM1factor_dR(RAB,KB,LB,MB);}
        aux = chgA*chgB*(integral_dR*fac1 - integral*fac2) - fac3 - chgA*chgB*fac4;
        if (doPDDG) {aux += dEPDDG(atoms[iatm],atoms[ibtm],RABangstroem);}
        enuclear_dr[index] = aux;
        ++index;
      }
    }
  }
  virtual void Enuclear_dR2(std::vector<double> & enuclear_dr2) {
    //function that calculates the second-derivative of the nuclear energy with respect to internuclear distance
    size_t index = 0;
    enuclear_dr2.resize(Natoms*(Natoms - 1)/2);
    for (size_t idx = 0; idx < Natoms*(Natoms - 1)/2; ++idx) {
      enuclear_dr2[idx] = 0.0;
    }
    double alphaA;
    double chgA;
    double alphaB;
    double chgB;
    double RAB;
    double XAB;
    double integral;
    double integral_dR;
    double integral_dR2;
    double aux;
    double gA;
    double gB;
    double gAdR;
    double gBdR;
    double expA;
    double expB;
    double fac1;
    double fac2;
    double fac3;
    double fac4;
    double fac5;
    double RAB5;
    double RR;
    for (size_t iatm = 0; iatm < Natoms; ++iatm) {
      //charge and nuclear energy variables
      KA = AM1K(atoms[iatm]);
      LA = AM1L(atoms[iatm]);
      MA = AM1M(atoms[iatm]);
      chgA = double(int(atoms[iatm]) - int(CoreCharge[iatm]));
      for (size_t ibtm = iatm + 1; ibtm < Natoms; ++ibtm) {
        alphaA = alpha(atoms[iatm],atoms[ibtm]);
        //charge and nuclear energy variables
        XAB = xAB(atoms[iatm],atoms[ibtm]);
        alphaB = alpha(atoms[ibtm],atoms[iatm]);
        KB = AM1K(atoms[ibtm],atoms[iatm]);
        LB = AM1L(atoms[ibtm],atoms[iatm]);
        MB = AM1M(atoms[ibtm],atoms[iatm]);
        if (atoms[iatm] == 5) {
          KA = AM1K(atoms[iatm],atoms[ibtm]);
          LA = AM1L(atoms[iatm],atoms[ibtm]);
          MA = AM1M(atoms[iatm],atoms[ibtm]);
        }
        chgB = double(int(atoms[ibtm]) - int(CoreCharge[ibtm]));
        //interatomic distance
        RAB = mol.AUdistance(iatm + 1,ibtm + 1);
        RR = RAB*dist_Angstrom2au;
        //integrals required
        integral = eri2Center(0,0,0,0,0,0,0,0,RAB,D,atoms[iatm],atoms[ibtm],2);
        integral_dR = eri2Center_dR(0,0,0,0,0,0,0,0,RAB,D,atoms[iatm],atoms[ibtm],2);
        integral_dR2 = eri2Center_dR2(0,0,0,0,0,0,0,0,RAB,D,atoms[iatm],atoms[ibtm],2);
        //auxiliaries
        fac1 = 1.0;                                     //factor assiciated to (mn|sl) dR2
        fac2 = 0.0;                                     //factor assiciated to (mn|sl) dR
        fac3 = 0.0;                                     //factor assiciated to (mn|sl)
        fac4 = 0.0;
        fac5 = 0.0;
        if (XAB > 0.0) {           //PM6
          if (((atoms[iatm] == 1)&&(atoms[ibtm] == 6))||((atoms[iatm] == 6)&&(atoms[ibtm] == 1))||((atoms[iatm] == 1)&&(atoms[ibtm] == 7))||((atoms[iatm] == 7)&&(atoms[ibtm] == 1))||((atoms[iatm] == 1)&&(atoms[ibtm] == 8))||((atoms[iatm] == 8)&&(atoms[ibtm] == 1))) {
            expA = 2.0*exp(-alphaA*RAB*RR);
            fac1 += XAB*expA;
            fac2 += 4.0*XAB*alphaA*expA*RR;
            fac3 += 2.0*XAB*alphaA*expA*(dist_Angstrom2au - 2.0*alphaA*RR*RR);
          }
          else {
            RAB5 = RR*RR*RR*RR*RR;
            expA = 2.0*exp(-alphaA*RAB*(1.0 + 0.0003*RAB5));
            fac1 += XAB*expA;
            fac2 += 2.0*XAB*alphaA*expA*(1.0 + 0.0018*RAB5);
            fac3 += XAB*alphaA*expA*(0.009*RAB5/RAB - alphaA*(1.0 + 0.0018*RAB5)*(1.0 + 0.0018*RAB5));
          }
          if ((atoms[iatm] == 6)&&(atoms[ibtm] == 6)) {
            expA = exp(-5.98*RR);
            fac1 += 9.28*expA;
            fac2 += 110.9888*expA*dist_Angstrom2au;
            fac3 -= 331.856512*expA*dist_Angstrom2au*dist_Angstrom2au;
          }
          else if (((atoms[iatm] == 8)&&(atoms[ibtm] == 14))||((atoms[iatm] == 14)&&(atoms[ibtm] == 8))) {
            aux = RR - 2.9;
            expA = exp(-aux*aux);
            fac1 -= 0.0007*expA;
            fac2 -= 0.0028*aux*expA*dist_Angstrom2au;
            fac3 -= 0.0014*expA*(1.0 - 2.0*aux*aux)*dist_Angstrom2au*dist_Angstrom2au;
          }
          fac4 = 1.56e-6*pow((pow(chgA,1.0/3.0) + pow(chgB,1.0/3.0))/RR,12.0)/(au2eV*RAB*RAB);
        }
        gA = gfactor(atoms[iatm],atoms[ibtm],RAB);
        gB = gfactor(atoms[ibtm],atoms[iatm],RAB);
        if ((gA > 0.0)||(gB > 0.0)) {
          gAdR = gfactor_dR(atoms[iatm],atoms[ibtm]);
          gBdR = gfactor_dR(atoms[ibtm],atoms[iatm]);
          expA = exp(-alphaA*RR);
          expB = exp(-alphaB*RR);
          fac1 += gA*expA + gB*expB;
          fac2 -= 2.0*((gAdR - alphaA*gA*dist_Angstrom2au)*expA + (gBdR - alphaB*gB*dist_Angstrom2au)*expB);
          fac3 -= (alphaA*(alphaA*gA*dist_Angstrom2au - 2.0*gAdR)*expA + alphaB*(alphaB*gB*dist_Angstrom2au - 2.0*gBdR)*expB)*dist_Angstrom2au;
        }
        if (fabs(KA[0]) > 1.0e-5) {fac5 += AM1factor_dR2(RAB,KA,LA,MA);}
        if (fabs(KB[0]) > 1.0e-5) {fac5 += AM1factor_dR2(RAB,KB,LB,MB);}
        aux = chgA*chgB*(integral_dR2*fac1 - integral_dR*fac2 - integral*fac3) + fac4 - chgA*chgB*fac5;
        if (doPDDG) {aux += d2EPDDG(atoms[iatm],atoms[ibtm],RR);}
        enuclear_dr2[index] = aux;
        ++index;
      }
    }
  }
  void se_integrals(matrixE & geometry, double tolerance = 1e-6) {
    //function calculating ALL the required integrals according to Dewar's semiempirical formulas
    //get total number of orbital pairs in system
    int irow = 0;
    int icol = 0;
    int naosA;
    int naosB;
    int nintA;
    int nintB;
    int AOpair = 0;
    for (size_t iatm = 0; iatm < Natoms; ++ iatm) {
      AOpair += AOs[iatm]*(AOs[iatm] + 1)/2;                    //saving solely unique pairs
    }
    gammaSE.resize(AOpair,AOpair);
    VAB.resize(AOpair,1);
    enuc = 0.0;
    //auxiliaries
    double intn;
    double RAB;
    double RR;
    double cost;
    double sint;
    double cosp;
    double sinp;
    double chgA;
    double chgB;
    double alphaA;
    double alphaB;
    double XAB;
    double gA;
    double gB;
    double factorA;
    double factorB;
    double factorC;
    for (size_t iatm = 0; iatm < Natoms; ++iatm) {
      //charge and nuclear energy variables
      KA = AM1K(atoms[iatm]);
      LA = AM1L(atoms[iatm]);
      MA = AM1M(atoms[iatm]);
      chgA = double(int(atoms[iatm]) - int(CoreCharge[iatm]));
      naosA = AOs[iatm];
      nintA = naosA*(naosA + 1)/2;
      icol = 0;
      for (size_t idr = 1; idr < nintA + 1; ++idr) {
        VAB(irow + idr,1) = 0.0;
      }
      for (size_t ibtm = 0; ibtm < Natoms; ++ibtm) {
        naosB = AOs[ibtm];
        nintB = naosB*(naosB + 1)/2;
        iBlockdR.resize(nintA,nintB);
        if (iatm == ibtm) {
          IntegralBlock1C(iBlockdR,atoms[iatm]);
          for (size_t idr = 1; idr < nintA + 1; ++idr) {
            gammaSE(irow + idr,icol + idr) = iBlockdR(idr,idr);
            for (size_t idc = idr + 1; idc < nintA + 1; ++idc) {
              gammaSE(irow + idr,icol + idc) = iBlockdR(idr,idc);
              gammaSE(icol + idc,irow + idr) = iBlockdR(idr,idc);
            }
          }
        }
        else {
          chgB = double(int(atoms[ibtm]) - int(CoreCharge[ibtm]));
          //getting the orientation vector
          rAB[0] = (geometry(ibtm + 1,1) - geometry(iatm + 1,1));                  //Delta x
          rAB[1] = (geometry(ibtm + 1,2) - geometry(iatm + 1,2));                  //Delta y
          rAB[2] = (geometry(ibtm + 1,3) - geometry(iatm + 1,3));                  //Delta z
          //normalizing it
          RAB = sqrt(rAB[0]*rAB[0] + rAB[1]*rAB[1] + rAB[2]*rAB[2]);
          rAB[0] /= RAB;
          rAB[1] /= RAB;
          rAB[2] /= RAB;
          //getting trigonometric functions for rotations
          cost = rAB[2];
          sint = sqrt(1.0 - cost*cost);
          cosp = 1.0;
          sinp = 0.0;
          if (fabs(sint) > tolerance) {
            cosp = rAB[0]/sint;
            sinp = rAB[1]/sint;
          }
          RR = RAB;
          RAB *= dist_Angstrom2aum1;
          iBlockdRVC.resize(nintA,1);
          IntegralBlock2C(1,iBlockdRVC,atoms[iatm],atoms[ibtm],RAB,cost,sint,cosp,sinp);
          for (size_t idr = 1; idr < nintA + 1; ++idr) {
            VAB(irow + idr,1) -= chgB*iBlockdRVC(idr,1);
          }
          if (ibtm > iatm) {
            //nuclear energy terms
            alphaA = alpha(atoms[iatm],atoms[ibtm]);
            XAB = xAB(atoms[iatm],atoms[ibtm]);
            KB = AM1K(atoms[ibtm],atoms[iatm]);
            LB = AM1L(atoms[ibtm],atoms[iatm]);
            MB = AM1M(atoms[ibtm],atoms[iatm]);
            if (atoms[iatm] == 5) {
              KA = AM1K(atoms[iatm],atoms[ibtm]);
              LA = AM1L(atoms[iatm],atoms[ibtm]);
              MA = AM1M(atoms[iatm],atoms[ibtm]);
            }
            factorA = 1.0;
            factorB = 0.0;
            factorC = 0.0;
            if (XAB > 0.0) {
              //assume here that these parameters take values considerably different from zero
              if (((atoms[iatm] == 1)&&(atoms[ibtm] == 6))||((atoms[iatm] == 6)&&(atoms[ibtm] == 1))||((atoms[iatm] == 1)&&(atoms[ibtm] == 7))||((atoms[iatm] == 7)&&(atoms[ibtm] == 1))||((atoms[iatm] == 1)&&(atoms[ibtm] == 8))||((atoms[iatm] == 8)&&(atoms[ibtm] == 1))) {
                factorA += 2.0*XAB*exp(-alphaA*RAB*RR);
              }
              else {
                factorA += 2.0*XAB*exp(-alphaA*RAB*(1.0 + 0.0003*RR*RR*RR*RR*RR));
              }
              if ((atoms[iatm] == 6)&&(atoms[ibtm] == 6)) {factorA += 9.28*exp(-5.98*RR);}
              else if (((atoms[iatm] == 8)&&(atoms[ibtm] == 14))||((atoms[iatm] == 14)&&(atoms[ibtm] == 8))) {
                factorA -= 0.0007*exp(-(RR - 2.9)*(RR - 2.9));
              }
              factorB = 1.0e-8*pow((pow(chgA,1.0/3.0) + pow(chgB,1.0/3.0))/(RR),12.0)/au2eV;
            }
            gA = gfactor(atoms[iatm],atoms[ibtm],RAB);
            gB = gfactor(atoms[ibtm],atoms[iatm],RAB);
            alphaB = alpha(atoms[ibtm],atoms[iatm]);
            factorA += gA*exp(-alphaA*RR) + gB*exp(-alphaB*RR);
            if (fabs(KA[0]) > 1.0e-5) {factorC += AM1factor(RAB,KA,LA,MA);}
            if (fabs(KB[0]) > 1.0e-5) {factorC += AM1factor(RAB,KB,LB,MB);}
            intn = eri2Center(0,0,0,0,0,0,0,0,RAB,D,atoms[iatm],atoms[ibtm],2);
            enuc += chgA*chgB*intn*factorA + factorB + chgA*chgB*factorC;
            if (doPDDG) {enuc += EPDDG(atoms[iatm],atoms[ibtm],RR);}
            //2-electron integrals
            IntegralBlock2C(0,iBlockdR,atoms[iatm],atoms[ibtm],RAB,cost,sint,cosp,sinp);
            for (size_t idr = 1; idr < nintA + 1; ++idr) {
              for (size_t idc = 1; idc < nintB + 1; ++idc) {
                gammaSE(irow + idr,icol + idc) = iBlockdR(idr,idc);
                gammaSE(icol + idc,irow + idr) = iBlockdR(idr,idc);
              }
            }
          }
        }
        icol += nintB;
      }
      irow += nintA;
    }
    if (print > 0) {
      std::cout << "Matrix Semiempirical Gammas:" << std::endl;
      gammaSE.Print();
      std::cout << "VAB:" << std::endl;
      VAB.Print();
    }
  }
  virtual double SlaterCondonRadialIntegral(int atmnr, int L, int II, int JJ, double threshold = 1e-7) {return 0.0;}
  virtual std::vector<double> AM1K(int atomicnr, int atm2 = 0) {
    std::vector<double> am1k(1,0.0);
    return am1k;
  }
  virtual std::vector<double> AM1L(int atomicnr, int atm2 = 0) {
    std::vector<double> am1l;
    return am1l;
  }
  virtual std::vector<double> AM1M(int atomicnr, int atm2 = 0) {
    std::vector<double> am1m;
    return am1m;
  }
  double eri2Center(int L1, int M1, int L2, int M2, int L3, int M3, int L4, int M4, double RAB, std::vector<double> D,size_t atmnrA,size_t atmnrB, int core = 0) {
    //function that calculates 2 electron repulsion integrals as in the NDDO approximation as applied in MNDO
    //AOs contains AO information on the 4 orbitals involved in the integral, meaning that this is only for s,p basis sets/atoms
    //RAB is distance between centers
    //D is a vector with all necessary D parameters: {D1A, D2A, D1B, D2B}
    //M. J. S. Dewar, W. Thiel, Theoret. Chim. Acta (Berl.), 46, 89, 1977
    double eri = 0.0;
    double rhoA;
    double rhoB;
    int lA;
    int lB;
    int mA;
    int mB;
    chg1 = NonZeroMultipoleMom_sp(L1,M1,L2,M2);
    chg2 = NonZeroMultipoleMom_sp(L3,M3,L4,M4);
    for (size_t mp1 = 0; mp1 < chg1.size(); ++mp1) {
      lA = chg1[mp1][0];
      mA = chg1[mp1][1];
      for (size_t mp2 = 0; mp2 < chg2.size(); ++mp2) {
        lB = chg2[mp2][0];
        mB = chg2[mp2][1];
        if (((abs(mA) % 2 == 0)&&(abs(mB) % 2 == 0))||((abs(mA) % 2 != 0)&&(abs(mB) % 2 != 0))) {
          //even-even or odd-odd only
          if (abs(mA) != abs(mB)) {
            mA = abs(mA);
            mB = abs(mB);
          }
          rhoA = rho(atmnrA,lA);
          rhoB = rho(atmnrB,lB);
          if (core > 0) {rhoB = rhocore(atmnrB);}                        //in this case use the correct core rho for the ss atom
          if ((core > 1)||(core < 0)) {rhoA = rhocore(atmnrA);}          //in this case use the correct core rho for the ss atom
          if (mA != 3) {eri += se_multipole(lA,mA,lB,mB,RAB,rhoA + rhoB,D);}
          else {eri += 0.5*(se_multipole(2,2,2,2,RAB,rhoA + rhoB,D) - se_multipole(2,2,2,-2,RAB,rhoA + rhoB,D));}
        }
      }
    }
    return eri;
  }
  double eri2Center_dR(int L1, int M1, int L2, int M2, int L3, int M3, int L4, int M4, double RAB, std::vector<double> D,size_t atmnrA,size_t atmnrB, int core = 0) {
    //function that calculates the first-derivatives of 2 electron repulsion integrals as in the NDDO approximation as applied in MNDO
    //AOs contains AO information on the 4 orbitals involved in the integral, meaning that this is only for s,p basis sets/atoms
    //RAB is distance between centers
    //D is a vector with all necessary D parameters: {D1A, D2A, D1B, D2B}
    double eri = 0.0;
    double rhoA;
    double rhoB;
    int lA;
    int lB;
    int mA;
    int mB;
    chg1 = NonZeroMultipoleMom_sp(L1,M1,L2,M2);
    chg2 = NonZeroMultipoleMom_sp(L3,M3,L4,M4);
    for (size_t mp1 = 0; mp1 < chg1.size(); ++mp1) {
      lA = chg1[mp1][0];
      mA = chg1[mp1][1];
      for (size_t mp2 = 0; mp2 < chg2.size(); ++mp2) {
        lB = chg2[mp2][0];
        mB = chg2[mp2][1];
        if (((abs(mA) % 2 == 0)&&(abs(mB) % 2 == 0))||((abs(mA) % 2 != 0)&&(abs(mB) % 2 != 0))) {
          //even-even or odd-odd only
          if (abs(mA) != abs(mB)) {
            mA = abs(mA);
            mB = abs(mB);
          }
          rhoA = rho(atmnrA,lA);
          rhoB = rho(atmnrB,lB);
          if (core > 0) {rhoB = rhocore(atmnrB);}                        //in this case use the correct core rho for the ss atom
          if ((core > 1)||(core < 0)) {rhoA = rhocore(atmnrA);}          //in this case use the correct core rho for the ss atom
          if (mA != 3) {eri += se_multipole_dR(lA,mA,lB,mB,RAB,rhoA + rhoB,D);}
          else {eri += 0.5*(se_multipole_dR(2,2,2,2,RAB,rhoA + rhoB,D) - se_multipole_dR(2,2,2,-2,RAB,rhoA + rhoB,D));}
        }
      }
    }
    return eri;
  }
  double eri2Center_dR2(int L1, int M1, int L2, int M2, int L3, int M3, int L4, int M4, double RAB, std::vector<double> D,size_t atmnrA,size_t atmnrB, int core = 0) {
    //function that calculates the second-derivatives of 2 electron repulsion integrals as in the NDDO approximation as applied in MNDO
    //AOs contains AO information on the 4 orbitals involved in the integral, meaning that this is only for s,p basis sets/atoms
    //RAB is distance between centers
    //D is a vector with all necessary D parameters: {D1A, D2A, D1B, D2B}
    double eri = 0.0;
    double rhoA;
    double rhoB;
    int lA;
    int lB;
    int mA;
    int mB;
    chg1 = NonZeroMultipoleMom_sp(L1,M1,L2,M2);
    chg2 = NonZeroMultipoleMom_sp(L3,M3,L4,M4);
    for (size_t mp1 = 0; mp1 < chg1.size(); ++mp1) {
      lA = chg1[mp1][0];
      mA = chg1[mp1][1];
      for (size_t mp2 = 0; mp2 < chg2.size(); ++mp2) {
        lB = chg2[mp2][0];
        mB = chg2[mp2][1];
        if (((abs(mA) % 2 == 0)&&(abs(mB) % 2 == 0))||((abs(mA) % 2 != 0)&&(abs(mB) % 2 != 0))) {
          //even-even or odd-odd only
          if (abs(mA) != abs(mB)) {
            mA = abs(mA);
            mB = abs(mB);
          }
          rhoA = rho(atmnrA,lA);
          rhoB = rho(atmnrB,lB);
          if (core > 0) {rhoB = rhocore(atmnrB);}                        //in this case use the correct core rho for the ss atom
          if ((core > 1)||(core < 0)) {rhoA = rhocore(atmnrA);}          //in this case use the correct core rho for the ss atom
          if (mA != 3) {eri += se_multipole_dR2(lA,mA,lB,mB,RAB,rhoA + rhoB,D);}
          else {eri += 0.5*(se_multipole_dR2(2,2,2,2,RAB,rhoA + rhoB,D) - se_multipole_dR2(2,2,2,-2,RAB,rhoA + rhoB,D));}
        }
      }
    }
    return eri;
  }
  virtual double d_eri2Center(int L1, int M1, int L2, int M2, int L3, int M3, int L4, int M4, double RAB, std::vector<double> D, size_t atmnrA, size_t atmnrB, int core = 0) {
    return 0.0;
  }
  virtual double d_eri2Center_dR(int L1, int M1, int L2, int M2, int L3, int M3, int L4, int M4, double RAB, std::vector<double> D, size_t atmnrA, size_t atmnrB, int core = 0) {
    return 0.0;
  }
  virtual double d_eri2Center_dR2(int L1, int M1, int L2, int M2, int L3, int M3, int L4, int M4, double RAB, std::vector<double> D, size_t atmnrA, size_t atmnrB, int core = 0) {
    return 0.0;
  }
  //other auxiliary functions
  virtual void checkAtoms() {
    //function that checks whether the atoms in the molecule are fully defined in the respective MNDO theories
    bool def;
    for (size_t iatm = 0; iatm < Natoms; ++iatm) {
      def = false;
      if (atoms[iatm] == 1) {def = true;}                                                       //H
      else if ((atoms[iatm] > 2)&&(atoms[iatm] < 10)) {def = true;}                             //Li,Be,B,C,N,O,F
      else if ((atoms[iatm] > 12)&&(atoms[iatm] < 18)) {def = true;}                            //Al,Si,P,S,Cl
      else if ((atoms[iatm] == 30)||(atoms[iatm] == 32)) {def = true;}                          //Zn,Ge
      else if ((atoms[iatm] == 35)||(atoms[iatm] == 53)) {def = true;}                          //Br,I
      else if ((atoms[iatm] == 50)||(atoms[iatm] == 80)||(atoms[iatm] == 82)) {def = true;}     //Sn,Hg,Pb
      if (!def) {throw("ERROR: MNDO.hpp: MNDO: checkAtoms(): atom not fully specified for MNDO-theory");}
    }
  }
  void OpenClosed() {
    //function to determine whether system is open- or closed-shell; also gets the number of core electrons per atom
    shell = "closed";
    ncoreelectrons = 0;
    if (mol.Multiplicity() > 1) {shell = "open";}
    else if (mol.Nelectrons()%2 != 0) {shell = "open";}
    if (CoreCharge.size() != Natoms) {CoreCharge.resize(Natoms);}
    for (size_t idx = 0; idx < Natoms; ++idx) {
      if (atoms[idx] < 3) {CoreCharge[idx] = 0;}                                //H,He
      else if ((atoms[idx] > 2)&&(atoms[idx] < 11))  {CoreCharge[idx] = 2;}     //Li-Ne
      else if ((atoms[idx] > 10)&&(atoms[idx] < 19)) {CoreCharge[idx] = 10;}    //Na-Ar
      else if ((atoms[idx] > 18)&&(atoms[idx] < 30)) {CoreCharge[idx] = 18;}    //K-Cu
      else if ((atoms[idx] > 29)&&(atoms[idx] < 37)) {CoreCharge[idx] = 28;}    //Zn-Kr
      else if ((atoms[idx] > 36)&&(atoms[idx] < 48)) {CoreCharge[idx] = 36;}    //Rb-Ag
      else if ((atoms[idx] > 47)&&(atoms[idx] < 55)) {CoreCharge[idx] = 46;}    //Cd-Xe
      else if ((atoms[idx] > 54)&&(atoms[idx] < 58)) {CoreCharge[idx] = 54;}    //Cs-La
      else if ((atoms[idx] > 70)&&(atoms[idx] < 80)) {CoreCharge[idx] = 68;}    //Lu-Au
      else if ((atoms[idx] > 79)&&(atoms[idx] < 84)) {CoreCharge[idx] = 78;}    //Hg-Bi
      if ((atoms[idx] == 10)||(atoms[idx] == 18)||(atoms[idx] == 36)||(atoms[idx] == 54)) {CoreCharge[idx] += 2;}
      ncoreelectrons += CoreCharge[idx];
    }
    if (print > 0) {
      std::cout << "core charge" << std::endl;
      for (size_t idx = 0; idx < CoreCharge.size(); ++idx) {
        std::cout << CoreCharge[idx] << " ";
      }
      std::cout << std::endl;
    }
  }
  void DampOscillations(int iteration, matrixE & DM, matrixE & oDM, std::vector<double> & pvec) {
    //apply mopac's oscillation damping, i.e., check whether the density matrix is changing by more than a certain value between iterations
    natocc = 1.0 + double(shell != "open");
    factorA = 0.0;
    factorB = 0.0;
    factor = 0.0;
    dampfactor = 10.0;
    if (iteration > 3) {dampfactor = oscillationdampthresh;}
    oscillationdamp = ((iteration + 1)%3 != 0);
    sumA = 0.0;
    sumB = 0.0;
    for (size_t iorb = 0; iorb < NAOs; ++iorb) {
      dval = DM(iorb + 1,iorb + 1);
      sumA += dval;
      deltaA = fabs(dval - oDM(iorb + 1,iorb + 1));
      if (!oscillationdamp) {
        factorA += deltaA*deltaA;
        deltaA = (dval - 2.0*oDM(iorb + 1,iorb + 1) + pvec[iorb]);
        factorB += deltaA*deltaA;
      }
      pvec[iorb] = oDM(iorb + 1,iorb + 1);
      oDM(iorb + 1,iorb + 1) = dval;
    }
    if (factorB >= 0.0) {
      if (factorA < 100.0*factorB) {factor = sqrt(factorA/factorB);}
    }
    for (size_t iorb = 0; iorb < NAOs; ++iorb) {
      for (size_t jorb = 0; jorb < iorb; ++jorb) {
        dval = DM(iorb + 1,jorb + 1);
        odval = oDM(iorb + 1,jorb + 1);
        oDM(iorb + 1,jorb + 1) = dval + factor*(dval - odval);
        DM(iorb + 1,jorb + 1) = oDM(iorb + 1,jorb + 1);
        DM(jorb + 1,iorb + 1) = oDM(iorb + 1,jorb + 1);
      }
      if (fabs(oDM(iorb + 1,iorb + 1) - pvec[iorb]) > dampfactor) {
        odval = oDM(iorb + 1,iorb + 1) - pvec[iorb];
        if (odval != 0.0) {
          dval = fabs(odval);
          odval /= dval;
        }
        oDM(iorb + 1,iorb + 1) = pvec[iorb] + odval*dampfactor;
      }
      else {oDM(iorb + 1,iorb + 1) += factor*(oDM(iorb + 1,iorb + 1) - pvec[iorb]);}
      odval = oDM(iorb + 1,iorb + 1);
      oDM(iorb + 1,iorb + 1) = fmin(natocc,fmax(0.0,odval));
      sumB += oDM(iorb + 1,iorb + 1);
      DM(iorb + 1,iorb + 1) = oDM(iorb + 1,iorb + 1);
    }
    deltaA = sumA;
label:
    dval = 0.0;
    if (sumB > 1.0e-3) {dval = sumA/sumB;}
    sumA = deltaA;
    if ((sumB > 1.0e-3)&&(fabs(dval - 1.0) > 1.0e-5)) {
      sumB = 0.0;
      for (size_t iorb = 0; iorb < NAOs; ++iorb) {
        oDM(iorb + 1,iorb + 1) *= dval;
        odval = oDM(iorb + 1,iorb + 1) + 1.0e-20;
        oDM(iorb + 1,iorb + 1) = fmax(odval,0.0);
        if (oDM(iorb + 1,iorb + 1) > natocc) {
          oDM(iorb + 1,iorb + 1) = natocc;
          sumA -= natocc;
        }
        else {sumB += oDM(iorb + 1,iorb + 1);}
        DM(iorb + 1,iorb + 1) = oDM(iorb + 1,iorb + 1);
      }
      goto label;
    }
    oscillationdamp = true;
  }
  //charge model stuff
  std::vector<double> CM1charge() {
    //function returning the CM1 partial charges
    //J. W. Storer, D. J. Giesen, C. J. Cramer, D. G. Truhlar, J. Comp.-Aid. Mol. Des., 9, 87, 1995
    //get Mulliken charges; note that actual population analysis not required here
    AOs = basis.AtomNAOs(atoms);
    this->getDens(dens);
    if (this->Shell() == "open") {this->getbDens(bdens);}
    std::vector<double> q0k;
    std::vector<double> Dqk(Natoms,0.0);
    std::vector<double> qk(Natoms,0.0);
    std::vector<double> Bkkp(Natoms,0.0);
    for (size_t idatm = 0; idatm < Natoms; ++idatm) {
      qk[idatm] = double(CoreCharge[idatm]);
    }
    LMcharges(q0k,dens,AOs,qk);
    matrixE BKKp(Natoms,Natoms);
    ArmstrongBondOrder(BKKp);
    for (size_t irow = 0; irow < Natoms; ++irow) {
      qk[irow] = q0k[irow];
      Dqk[irow] += CM1ck(atoms[irow])*q0k[irow] + CM1dk(atoms[irow]);
      for (size_t icol = 0; icol < Natoms; ++icol) {
        if (irow == icol) {continue;}
        Bkkp[irow] += BKKp(irow + 1,icol + 1);
        if (atoms[irow] == 1) {Dqk[irow] += BKKp(irow + 1,icol + 1)*CM1dkkp(atoms[irow],atoms[icol]);}
        else if (atoms[irow] == 7) {
          if (atoms[icol] == 6) {Dqk[irow] += 0.5*(CM1ckkp(atoms[irow],atoms[icol])*q0k[irow] + CM1dkkp(atoms[irow],atoms[icol]))*(tanh(10.0*(BKKp(irow + 1,icol + 1) - 2.3)) + 1.0);}
          else if (atoms[icol] == 8) {Dqk[irow] += BKKp(irow + 1,icol + 1)*CM1dkkp(atoms[irow],atoms[icol]);}
        }
        else if (atoms[irow] == 8) {
          if (atoms[icol] == 16) {Dqk[irow] += BKKp(irow + 1,icol + 1)*CM1dkkp(atoms[irow],atoms[icol]);}
        }
      }
    }
    for (size_t iatm = 0; iatm < Natoms; ++iatm) {
      qk[iatm] += Bkkp[iatm]*Dqk[iatm];
      for (size_t ibtm = 0; ibtm < Natoms; ++ibtm) {
        if (iatm == ibtm) {continue;}
        qk[iatm] -= BKKp(iatm + 1,ibtm + 1)*Dqk[ibtm];
      }
    }
    return qk;
  }
  std::vector<double> CM2charge() {
    //function returning the CM2 partial charges
    //J. Li, T. Zhu, C. J. Cramer, D. G. Truhlar, J. Phys. Chem. A, 102, 1820, 1998
    AOs = basis.AtomNAOs(atoms);
    this->getDens(dens);
    if (this->Shell() == "open") {this->getbDens(bdens);}
    sao = Identity(NAOs);                                   //for NDDO methods the overlap is supposed to be the identity matrix
    std::vector<double> q0k(Natoms,0.0);
    std::vector<double> qk(Natoms,0.0);
    for (size_t idatm = 0; idatm < Natoms; ++idatm) {
      q0k[idatm] = double(CoreCharge[idatm]);
    }
    LMcharges(qk,dens,AOs,q0k);
    matrixE BKKp(Natoms,Natoms);
    MayerBondOrder(BKKp);
    for (size_t iatm = 0; iatm < Natoms; ++iatm) {
      for (size_t ibtm = 0; ibtm < Natoms; ++ibtm) {
        if (iatm == ibtm) {continue;}
        qk[iatm] += BKKp(iatm + 1,ibtm + 1)*(CM2dkkp(atoms[iatm],atoms[ibtm]) + CM2ckkp(atoms[iatm],atoms[ibtm])*BKKp(iatm + 1,ibtm + 1));
      }
    }
    return qk;
  }
  std::vector<double> CM3charge() {
    //function returning the CM3 partial charges
    //P. Winget, J. D. Thompson, J. D. Xidos, C. J. Cramer, D. G. Truhlar, J. Phys. Chem. A, 106(44), 10707 2002
    //J. D. Thompson, C. J. Cramer, D. G. Truhlar, J. Comp. Chem., 24(11), 1291, 2003
    int other;
    AOs = basis.AtomNAOs(atoms);
    this->getDens(dens);
    if (this->Shell() == "open") {this->getbDens(bdens);}
    sao = Identity(NAOs);                                   //for NDDO methods the overlap is supposed to be the identity matrix
    std::vector<double> q0k(Natoms,0.0);
    std::vector<double> qk(Natoms,0.0);
    for (size_t idatm = 0; idatm < Natoms; ++idatm) {
      q0k[idatm] = double(CoreCharge[idatm]);
    }
    LMcharges(qk,dens,AOs,q0k);
    matrixE BKKp(Natoms,Natoms);
    MayerBondOrder(BKKp);
    for (size_t iatm = 0; iatm < Natoms; ++iatm) {
      if ((atoms[iatm] != 7)&&(atoms[iatm] != 8)) {
        for (size_t ibtm = 0; ibtm < Natoms; ++ibtm) {
          if (iatm == ibtm) {continue;}
          qk[iatm] += BKKp(iatm + 1,ibtm + 1)*(CM3dkkp(atoms[iatm],atoms[ibtm]) + CM3ckkp(atoms[iatm],atoms[ibtm])*BKKp(iatm + 1,ibtm + 1));
        }
      }
      else {
        if (atoms[iatm] == 8) {other = 7;}
        else {other = 8;}
        for (size_t ibtm = 0; ibtm < Natoms; ++ibtm) {
          if (iatm == ibtm) {continue;}
          if ((atoms[iatm] != 7)||(atoms[iatm] != 8)) {
            qk[iatm] += BKKp(iatm + 1,ibtm + 1)*(CM3dkkp(atoms[iatm],atoms[ibtm]) + CM3ckkp(atoms[iatm],atoms[ibtm])*BKKp(iatm + 1,ibtm + 1));}
          if (atoms[iatm] == other) {
            qk[iatm] += BKKp(iatm + 1,ibtm + 1)*(CM3dkkp(atoms[iatm],atoms[ibtm]) + CM3ckkp(atoms[iatm],atoms[ibtm])*exp(-(BKKp(iatm + 1,ibtm + 1)*BKKp(iatm + 1,ibtm + 1))/(B0()*B0())));}
        }
      }
    }
    return qk;
  }
  //estimate parameters for heat of formation
  virtual double calcEISOL(int atmnr) {
    //function that calculates the theoretical heats of formation for atoms
    //this function is not used directly by the methods, its purpose is the generation of the parameters that are then tabulated
    //integrals
    double Uss = UlX(atmnr,0);
    double Upp = UlX(atmnr,1);
    double Udd = UlX(atmnr,2);
    double Gss = eri1Center(atmnr,0,0);
    double Gpp = eri1Center(atmnr,2,2);
    double Gsp = eri1Center(atmnr,0,2);
    double Gp2 = eri1Center(atmnr,2,-2);
    double Hsp = eri1Center(atmnr,1,1);
    if (Hsp < 1.0e-7) {Hsp = 1.0e-7;}
    //Take into account constraints on the values of the integrals
    double Hpp = 0.5*(Gpp - Gp2);
    if (Hpp < 0.1) {Hpp = 0.1;}
    //counting how often integrals show up
    int nrs = NrS(atmnr);
    int nrp = NrP(atmnr);
    int nrd = NrD(atmnr);
    int ngss = (nrs - 1)*int(nrs > 1);                     //number of integrals of type <ss|ss>
    int ngsp = nrs*nrp;                                    //number of integrals of type <ss|pp>
    int iaux = std::min(nrp,6 - nrp);
    int npp = -iaux*(iaux - 1)/4;
    int ngp2 = nrp*(nrp - 1)/2 - npp;                      //number of integrals of type <pp|pp>; note that hpp is replaced by 0.5*(gpp - gp2) to insure rotational invariance
    int nhsp = -nrp;                                       //number of integrals of type <sp|sp>; note that if nrp != 0 then nrs = 2 due to Aufbau principle
    double eisol = Uss*nrs + Upp*nrp + Udd*nrd + Gss*ngss + Gpp*npp + Gsp*ngsp + Gp2*ngp2 + Hsp*nhsp;
    return eisol;
  }
  //parameters
  virtual double B0() {return 1.0;}
  virtual double ZeroOverlap(size_t atm) {return ZeroOverlapMNDO(atm);}
  virtual double ElementFormTheo(size_t atomicnr) {
    //function returning theoretical heats of formation for elements; values directly in a.u.
    double enth = 0.0;
    switch (atomicnr) {
      case 1:            //H
        enth = -0.437547617217;
        break;
      case 3:            //Li
        enth = -0.188450543318;
        break;
      case 4:            //Be
        enth = -0.889508467058;
        break;
      case 5:            //B
        enth = -2.3635678084;
        break;
      case 6:            //C
        enth = -4.40553146138;
        break;
      case 7:            //N
        enth = -7.41863055941;
        break;
      case 8:            //O
        enth = -11.6531563457;
        break;
      case 9:            //F
        enth = -17.5178076287;
        break;
      case 13:           //Al
        enth = -1.63475966018;
        break;
      case 14:           //Si
        enth = -3.03014442762;
        break;
      case 15:           //P
        enth = -5.60353574631;
        break;
      case 16:           //S
        enth = -8.2861420927;
        break;
      case 17:           //Cl
        enth = -12.9768362327;
        break;
      case 30:           //Zn
        enth = -1.09804898487;
        break;
      case 32:           //Ge
        enth = -2.78739756374;
        break;
      case 35:           //Br
        enth = -12.7403019342;
        break;
      case 50:           //Sn
        enth = -3.3781487843;
        break;
      case 53:           //I
        enth = -12.5167626458;
        break;
      case 80:           //Hg
        enth = -1.05908426259;
        break;
      case 82:           //Pb
        enth = -3.87464695071;
        break;
    }
    return enth;
  }
  virtual double betaA0(size_t atomicnr, int L) {
    //function returning the bonding parameters; values stored in eV, returned however in a.u.
    double betaa0 = 0.0;
    switch (atomicnr) {
      case 1:                   //H
        betaa0 = -6.989064;  
        break;
      case 3:                   //Li
        betaa0 = -1.3500400;
        break;
      case 4:                   //Be
        betaa0 = -4.0170960;
        break;
      case 5:                   //B
        betaa0 = -8.252054;
        break;
      case 6:                   //C
        if (L == 0) {betaa0 = -18.985044;}
        else if (L == 1) {betaa0 = -7.934122;}
        break;
      case 7:                   //N
        betaa0 = -20.495758;
        break;
      case 8:                   //O
        betaa0 = -32.688082;
        break;
      case 9:                   //F
        if (L == 0) {betaa0 = -48.290460;}
        else if (L == 1) {betaa0 = -36.508540;}
        break;
      case 13:                  //Al
        betaa0 = -2.6702840;
        break;
      case 14:                  //Si
        if (L == 0) {betaa0 = -9.0868040;}
        else if (L == 1) {betaa0 = -1.0758270;}
        break;
      case 15:                  //P
        betaa0 = -6.7916000;
        break;
      case 16:                  //S
        if (L == 0) {betaa0 = -10.7616700;}
        else if (L == 1) {betaa0 = -10.1084330;}
        break;
      case 17:                  //Cl
        betaa0 = -14.2623200;
        break;
      case 30:                  //Zn
        if (L == 0) {betaa0 = -1.000;}
        else if (L == 1) {betaa0 = -2.000;}
        break;
      case 32:                  //Ge
        if (L == 0) {betaa0 = -4.5164790;}
        else if (L == 1) {betaa0 = -1.7555170;}
        break;
      case 35:                  //Br
        if (L == 0) {betaa0 = -8.91710680;}
        else if (L == 1) {betaa0 = -9.9437398;}
        break;
      case 50:                  //Sn
        if (L == 0) {betaa0 = -3.2351470;}
        else if (L == 1) {betaa0 = -4.2904160;}
        break;
      case 53:                  //I
        if (L == 0) {betaa0 = -7.4144514;}
        else if (L == 1) {betaa0 = -6.1967812;}
        break;
      case 80:                  //Hg
        if (L == 0) {betaa0 = -0.4045250;}
        else if (L == 1) {betaa0 = -6.2066830;}
        break;
      case 82:                  //Pb
        if (L == 0) {betaa0 = -8.0423870;}
        else if (L == 1) {betaa0 = -3.00000;}
        break;
    }
    return betaa0/au2eV;
  }
  virtual double UlX(size_t atomicnr, int L) {
    //function that returns the atomic potential U
    double ulx = 0.0;
    switch (atomicnr) {
      case 1:                 //H
        ulx = -11.906276;
        break;
      case 3:                 //Li
        if (L == 0) {ulx = -5.1280000;}
        else if (L == 1) {ulx = -2.7212000;}
        break;
      case 4:                 //Be
        if (L == 0) {ulx = -16.6023780;}
        else if (L == 1) {ulx = -10.7037710;}
        break;
      case 5:                 //B
        if (L == 0) {ulx = -34.547130;}
        else if (L == 1) {ulx = -23.121690;}
        break;
      case 6:                 //C
        if (L == 0) {ulx = -52.279745;}
        else if (L == 1) {ulx = -39.205558;}
        break;
      case 7:                 //N
        if (L == 0) {ulx = -71.932122;}
        else if (L == 1) {ulx = -57.172319;}
        break;
      case 8:                 //O
        if (L == 0) {ulx = -99.644309;}
        else if (L == 1) {ulx = -77.797472;}
        break;
      case 9:                 //F
        if (L == 0) {ulx = -131.071548;}
        else if (L == 1) {ulx = -105.782137;}
        break;
      case 13:                //Al
        if (L == 0) {ulx = -23.8070970;}
        else if (L == 1) {ulx = -17.5198780;}
        break;
      case 14:                //Si
        if (L == 0) {ulx = -37.0375330;}
        else if (L == 1) {ulx = -27.76967800;}
        break;
      case 15:                //P
        if (L == 0) {ulx = -56.1433600;}
        else if (L == 1) {ulx = -42.8510800;}
        break;
      case 16:                //S
        if (L == 0) {ulx = -72.2422810;}
        else if (L == 1) {ulx = -56.9732070;}
        break;
      case 17:                //Cl
        if (L == 0) {ulx = -100.2271660;}
        else if (L == 1) {ulx = -77.3786670;}
        break;
      case 30:                //Zn
        if (L == 0) {ulx = -20.8397160;}
        else if (L == 1) {ulx = -19.6252240;}
        break;
      case 32:                //Ge
        if (L == 0) {ulx = -33.9493670;}
        else if (L == 1) {ulx = -27.4251050;}
        break;
      case 35:                //Br
        if (L == 0) {ulx = -99.98644054;}
        else if (L == 1) {ulx = -75.67130754;}
        break;
      case 50:                //Sn
        if (L == 0) {ulx = -40.8518020;}
        else if (L == 1) {ulx = -28.5602490;}
        break;
      case 53:                //I
        if (L == 0) {ulx = -100.0030538;}
        else if (L == 1) {ulx = -74.61146919;}
        break;
      case 80:                //Hg
        if (L == 0) {ulx = -19.8095740;}
        else if (L == 1) {ulx = -13.1025300;}
        break;
      case 82:                //Pb
        if (L == 0) {ulx = -47.3196920;}
        else if (L == 1) {ulx = -28.8475600;}
        break;
    }
    return ulx/au2eV;
  }
  virtual double alpha(int atomicnr, int atm2) {
    //function that returns alpha values for MNDO
    //alphas stored and returned in 1/Angstrom
    double alpha = 0.0;
    switch (atomicnr) {
      case 1:     //H
        alpha = 2.544134;
        break;
      case 3:     //Li
        alpha = 1.250140;
        break;
      case 4:     //Be
        alpha = 1.669434;
        break;
      case 5:     //B
        alpha = 2.134993;
        break;
      case 6:     //C
        alpha = 2.546380;
        break;
      case 7:     //N
        alpha = 2.861342;
        break;
      case 8:     //O
        alpha = 3.160604;
        break;
      case 9:     //F
        alpha = 3.419661;
        break;
      case 13:    //Al
        alpha = 1.868839;
        break;
      case 14:    //Si
        alpha = 2.205316;
        break;
      case 15:    //P
        alpha = 2.415280;
        break;
      case 16:    //S
        alpha = 2.478026;
        break;
      case 17:    //Cl
        alpha = 2.542201;
        break;
      case 30:    //Zn
        alpha = 1.506457;
        break;
      case 32:    //Ge
        alpha = 1.978498;
        break;
      case 35:    //Br
        alpha = 2.44570512;
        break;
      case 50:    //Sn
        alpha = 1.800814;
        break;
      case 53:    //I
        alpha = 2.20732001;
        break;
      case 80:    //Hg
        alpha = 1.335641;
        break;
      case 82:    //Pb
        alpha = 1.728333;
        break;
    }
    return alpha;
  }
  virtual double Dvalue(size_t atmnr, size_t idx) {
    //function returning the D values needed to calculate eris; values stored in Angstrom, but returned in a.u.
    //note that this function returns both D1 and D2; idx is then either 1 or 2
    double dval = 0.0;
    switch (atmnr) {
      case 1:                    //H
        dval = 0.0;
        break;
      case 3:                    //Li
        if (idx == 1) {dval = 1.0874267;}
        else if (idx == 2) {dval = 0.9227121;}
        break;
      case 4:                    //Be
        if (idx == 1) {dval = 0.7605847;}
        else if (idx == 2) {dval = 0.6453775;}
        break;
      case 5:                    //B
        if (idx == 1) {dval = 0.5068929;}
        else if (idx == 2) {dval = 0.4301129;}
        break;
      case 6:                    //C
        if (idx == 1) {dval = 0.4272845;}
        else if (idx == 2) {dval = 0.3625629;}
        break;
      case 7:                    //N
        if (idx == 1) {dval = 0.3386159;}
        else if (idx == 2) {dval = 0.2873251;}
        break;
      case 8:                    //O
        if (idx == 1) {dval = 0.2828939;}
        else if (idx == 2) {dval = 0.2400435;}
        break;
      case 9:                    //F
        if (idx == 1) {dval = 0.2681377;}
        else if (idx == 2) {dval = 0.2275224;}
        break;
      case 13:                   //Al
        if (idx == 1) {dval = 0.7404309;}
        else if (idx == 2) {dval = 0.6131351;}
        break;
      case 14:                   //Si
        if (idx == 1) {dval = 0.6657106;}
        else if (idx == 2) {dval = 0.5178335;}
        break;
      case 15:                   //P
        if (idx == 1) {dval = 0.5360302;}
        else if (idx == 2) {dval = 0.4958342;}
        break;
      case 16:                   //S
        if (idx == 1) {dval = 0.4863010;}
        else if (idx == 2) {dval = 0.4407175;}
        break;
      case 17:                   //Cl
        if (idx == 1) {dval = 0.2638887;}
        else if (idx == 2) {dval = 0.4348484;}
        break;
      case 30:                   //Zn
        if (idx == 1) {dval = 0.6899187;}
        else if (idx == 2) {dval = 0.7683602;}
        break;
      case 32:                   //Ge
        if (idx == 1) {dval = 0.6644269;}
        else if (idx == 2) {dval = 0.5555542;}
        break;
      case 35:                   //Br
        if (idx == 1) {dval = 0.3202029;}
        else if (idx == 2) {dval = 0.5104278;}
        break;
      case 50:                   //Sn
        if (idx == 1) {dval = 0.8306740;}
        else if (idx == 2) {dval = 0.7017967;}
        break;
      case 53:                   //I
        if (idx == 1) {dval = 0.7542341;}
        else if (idx == 2) {dval = 0.6266241;}
        break;
      case 80:                   //Hg
        if (idx == 1) {dval = 0.9195890;}
        else if (idx == 2) {dval = 0.7730105;}
        break;
      case 82:                   //Pb
        if (idx == 1) {dval = 0.8216177;}
        else if (idx == 2) {dval = 0.7666867;}
        break;
    }
    return dval*dist_Angstrom2aum1;
  }
  virtual double rhocore(size_t atmnr) {
    //function returning rho values used to calculate VAB integrals; values already in atomic units
    double rho_ = rho(atmnr,0);
    return rho_;
  }
  virtual double rho(size_t atmnr, size_t l) {
    //function returning the rho values needed to calculate eris; values stored in Angstrom but returned in atomic units
    double rho = 0.0;
    switch (atmnr) {
      case 1:                    //H
        if (l == 0) {rho = 0.560345403;}
        else if (l == 1) {rho = 0.0;}
        else if (l == 2) {rho = 0.0;}
        break;
      case 3:                    //Li
        if (l == 0) {rho = 0.986207884;}
        else if (l == 1) {rho = 1.165672376;}
        else if (l == 2) {rho = 1.011953732;}
        break;
      case 4:                    //Be
        if (l == 0) {rho = 0.799924235;}
        else if (l == 1) {rho = 0.788356095;}
        else if (l == 2) {rho = 0.687877905;}
        break;
      case 5:                    //B
        if (l == 0) {rho = 0.679822279;}
        else if (l == 1) {rho = 0.539445596;}
        else if (l == 2) {rho = 0.476128306;}
        break;
      case 6:                    //C
        if (l == 0) {rho = 0.588660438;}
        else if (l == 1) {rho = 0.430253872;}
        else if (l == 2) {rho = 0.395733736;}
        break;
      case 7:                    //N
        if (l == 0) {rho = 0.529751104;}
        else if (l == 1) {rho = 0.337322211;}
        else if (l == 2) {rho = 0.325583327;}
        break;
      case 8:                    //O
        if (l == 0) {rho = 0.466881794;}
        else if (l == 1) {rho = 0.275821517;}
        else if (l == 2) {rho = 0.2786282;}
        break;
      case 9:                    //F
        if (l == 0) {rho = 0.425491557;}
        else if (l == 1) {rho = 0.243848996;}
        else if (l == 2) {rho = 0.255793341;}
        break;
      case 13:                   //Al
        if (l == 0) {rho = 0.889903107;}
        else if (l == 1) {rho = 1.00389327;}
        else if (l == 2) {rho = 0.720237317;}
        break;
      case 14:                   //Si
        if (l == 0) {rho = 0.733128067;}
        else if (l == 1) {rho = 0.722068454;}
        else if (l == 2) {rho = 0.587084012;}
        break;
      case 15:                   //P
        if (l == 0) {rho = 0.622778301;}
        else if (l == 1) {rho = 0.541910569;}
        else if (l == 2) {rho = 0.531355547;}
        break;
      case 16:                   //S
        if (l == 0) {rho = 0.558953167;}
        else if (l == 1) {rho = 0.477199756;}
        else if (l == 2) {rho = 0.473718785;}
        break;
      case 17:                   //Cl
        if (l == 0) {rho = 0.478996434;}
        else if (l == 1) {rho = 0.328217689;}
        else if (l == 2) {rho = 0.435988;}
        break;
      case 30:                   //Zn
        if (l == 0) {rho = 0.610111605;}
        else if (l == 1) {rho = 1.113608164;}
        else if (l == 2) {rho = 0.966035844;}
        break;
      case 32:                   //Ge
        if (l == 0) {rho = 0.734624198;}
        else if (l == 1) {rho = 0.726135254;}
        else if (l == 2) {rho = 0.608610513;}
        break;
      case 35:                   //Br
        if (l == 0) {rho = 0.478791611;}
        else if (l == 1) {rho = 0.364523933;}
        else if (l == 2) {rho = 0.474624228;}
        break;
      case 50:                   //Sn
        if (l == 0) {rho = 0.734624198;}
        else if (l == 1) {rho = 0.821688399;}
        else if (l == 2) {rho = 0.712428177;}
        break;
      case 53:                   //I
        if (l == 0) {rho = 0.478664021;}
        else if (l == 1) {rho = 0.576001573;}
        else if (l == 2) {rho = 0.577015931;}
        break;
      case 80:                   //Hg
        if (l == 0) {rho = 0.666603429;}
        else if (l == 1) {rho = 0.868143258;}
        else if (l == 2) {rho = 0.759620304;}
        break;
      case 82:                   //Pb
        if (l == 0) {rho = 0.734624198;}
        else if (l == 1) {rho = 0.816789939;}
        else if (l == 2) {rho = 0.755508834;}
        break;
    }
    return rho*dist_Angstrom2aum1;
  }
  virtual double eri1Center(int atmnr, int Lbra, int Lket) {
    //function that gives back the semi-empirical 1-center eris
    //values stored in eV but returned in a.u.
    //Lbra is the sum of azimuthal quantum numbers for bra (ss = 0; pp = 2; sp = 1)
    //Lket is the sum of azimuthal quantum numbers for ket (ss = 0; pp = 2; sp = 1; p*p* = -2)
    double eri = 0.0;
    switch (atmnr) {
      case 1:            //H
        eri = 12.848;                                                                            //(ss|ss)
        break;
      case 3:           //Li
        if ((Lbra == 0)&&(Lket == 0)) {eri = 7.3;}                                               //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 5.0;}                                          //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 5.42;}           //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 4.52;}         //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 0.83;}                                         //(sp|sp)||(ps|ps)
        break;
      case 4:          //Be
        if ((Lbra == 0)&&(Lket == 0)) {eri = 9.0;}                                               //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 6.97;}                                         //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 7.43;}           //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 6.22;}         //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 1.28;}                                         //(sp|sp)||(ps|ps)
        break;
      case 5:          //B
        if ((Lbra == 0)&&(Lket == 0)) {eri = 10.59;}                                             //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 8.86;}                                         //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 9.56;}           //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 7.86;}         //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 1.81;}                                         //(sp|sp)||(ps|ps)
        break;
      case 6:          //C
        if ((Lbra == 0)&&(Lket == 0)) {eri = 12.23;}                                             //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 11.08;}                                        //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 11.47;}          //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 9.84;}         //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 2.43;}                                         //(sp|sp)||(ps|ps)
        break;
      case 7:          //N
        if ((Lbra == 0)&&(Lket == 0)) {eri = 13.59;}                                             //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 12.98;}                                        //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 12.66;}          //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 11.59;}        //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 3.14;}                                         //(sp|sp)||(ps|ps)
        break;
      case 8:          //O
        if ((Lbra == 0)&&(Lket == 0)) {eri = 15.42;}                                             //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 14.52;}                                        //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 14.48;}          //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 12.98;}        //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 3.94;}                                         //(sp|sp)||(ps|ps)
        break;
      case 9:          //F
        if ((Lbra == 0)&&(Lket == 0)) {eri = 16.92;}                                             //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 16.71;}                                        //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 17.25;}          //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 14.91;}        //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 4.83;}                                         //(sp|sp)||(ps|ps)
        break;
      case 13:         //Al
        if ((Lbra == 0)&&(Lket == 0)) {eri = 8.09;}                                              //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 5.98;}                                         //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 6.63;}           //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 5.40;}         //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 0.70;}                                         //(sp|sp)||(ps|ps)
        break;
      case 14:         //Si
        if ((Lbra == 0)&&(Lket == 0)) {eri = 9.82;}                                              //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 7.31;}                                         //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 8.36;}           //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 6.54;}         //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 1.32;}                                         //(sp|sp)||(ps|ps)
        break;
      case 15:         //P
        if ((Lbra == 0)&&(Lket == 0)) {eri = 11.56;}                                             //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 8.64;}                                         //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 10.08;}          //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 7.68;}         //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 1.92;}                                         //(sp|sp)||(ps|ps)
        break;
      case 16:         //S
        if ((Lbra == 0)&&(Lket == 0)) {eri = 12.88;}                                             //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 9.90;}                                         //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 11.26;}          //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 8.83;}         //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 2.26;}                                         //(sp|sp)||(ps|ps)
        break;
      case 17:         //Cl
        if ((Lbra == 0)&&(Lket == 0)) {eri = 15.03;}                                             //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 11.30;}                                        //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 13.16;}          //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 9.97;}         //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 2.42;}                                         //(sp|sp)||(ps|ps)
        break;
      case 30:         //Zn
        if ((Lbra == 0)&&(Lket == 0)) {eri = 11.8;}                                              //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 13.3;}                                         //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 11.182018;}      //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 12.93052;}     //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 0.484606;}                                     //(sp|sp)||(ps|ps)
        break;
      case 32:         //Ge
        if ((Lbra == 0)&&(Lket == 0)) {eri = 9.8;}                                               //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 7.3;}                                          //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 8.3;}            //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 6.5;}          //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 1.3;}                                          //(sp|sp)||(ps|ps)
        break;
      case 35:         //Br
        if ((Lbra == 0)&&(Lket == 0)) {eri = 15.03643948;}                                       //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 11.27632539;}                                  //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 13.03468242;}    //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 9.85442552;}   //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 2.45586832;}                                   //(sp|sp)||(ps|ps)
        break;
      case 50:         //Sn
        if ((Lbra == 0)&&(Lket == 0)) {eri = 9.8;}                                               //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 7.3;}                                          //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 8.3;}            //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 6.5;}          //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 1.3;}                                          //(sp|sp)||(ps|ps)
        break;
      case 53:         //I
        if ((Lbra == 0)&&(Lket == 0)) {eri = 15.04044855;}                                       //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 11.14778369;}                                  //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 13.05655798;}    //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 9.91409071;}   //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 2.45638202;}                                   //(sp|sp)||(ps|ps)
        break;
      case 80:         //Hg
        if ((Lbra == 0)&&(Lket == 0)) {eri = 10.8;}                                              //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 14.3;}                                         //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 9.3;}            //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 13.5;}         //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 1.3;}                                          //(sp|sp)||(ps|ps)
        break;
      case 82:         //Pb
        if ((Lbra == 0)&&(Lket == 0)) {eri = 9.8;}                                               //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 7.3;}                                          //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 8.3;}            //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 6.5;}          //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 1.3;}                                          //(sp|sp)||(ps|ps)
        break;
    }
    return eri/au2eV;
  }
  virtual double xAB(int atm1, int atm2) {
    //function returning xAB values for MNDOd, all zero
    return 0.0;
  }
  virtual double CM1ck(size_t atomicnr) {
    //function returning the atomic parameters ck for CM1
    return 0.0;
  }
  virtual double CM1dk(size_t atomicnr) {
    //function returning the atomic parameters dk for CM1
    return 0.0;
  }
  virtual double CM1ckkp(size_t atom1, size_t atom2) {
    //function returning the atom-pair parameters ckk' for CM1
    return 0.0;
  }
  virtual double CM1dkkp(size_t atom1, size_t atom2) {
    //function returning the atom-pair parameters dkk' for CM1
    return 0.0;
  }
  virtual double CM2ckkp(size_t atom1, size_t atom2) {
    //function returning the atom-pair parameters ckk' for CM2
    return 0.0;
  }
  virtual double CM2dkkp(size_t atom1, size_t atom2) {
    //function returning the atom-pair parameters dkk' for CM2
    return 0.0;
  }
  virtual double CM3ckkp(size_t atom1, size_t atom2) {
    //function returning the atom-pair parameters ckk' for CM3
    return 0.0;
  }
  virtual double CM3dkkp(size_t atom1, size_t atom2) {
    //function returning the atom-pair parameters dkk' for CM3
    return 0.0;
  }
  virtual double PAPDDG(int atmnr, int index) {return 0.0;}
  virtual double DAPDDG(int atmnr, int index) {return 0.0;}
  double ValenceElectrons(int atmA) {
    //function returning the number of valence electrons of a certain atom
    double valelectron = 0.0;
    switch (atmA) {
      case 1:      //H
        valelectron = 1.0;
        break;
      case 6:      //C
        valelectron = 4.0;
        break;
      case 7:      //N
        valelectron = 5.0;
        break;
      case 8:      //O
        valelectron = 6.0;
        break;
      case 9:      //F
        valelectron = 7.0;
        break;
      case 14:     //Si
        valelectron = 4.0;
        break;
      case 15:     //P
        valelectron = 5.0;
        break;
      case 16:     //S
        valelectron = 6.0;
        break;
      case 17:     //Cl
        valelectron = 7.0;
        break;
      case 35:     //Br
        valelectron = 7.0;
        break;
      case 53:     //I
        valelectron = 7.0;
        break;
    }
    return valelectron;
  }
};

#endif //_MNDO_
