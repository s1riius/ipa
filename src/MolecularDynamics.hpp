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

#ifndef _Molecular_Dynamics_
#define _Molecular_Dynamics_
#include <vector>
#include <string>
#include <math.h>
#include <time.h>
#include "QC.hpp"
#include "BSet.hpp"
#include "Molecule.hpp"
#include "Potentials.hpp"
#include "ConstantsPackage.hpp"
#include "atoms/AtomPackage.hpp"
#include "math/MatrixPackage.hpp"
#include "math/RandomPackage.hpp"
#include "math/SolverPackage.hpp"

#include "other/auxMCMD.hpp"
#include "other/output.hpp"
#include "other/python.hpp"

//description:
//molecular dynamics package

std::vector<int> getNumbersFromFile(std::string filename) {
  //function reading a series of comma separated numbers from a file
  std::vector<int> numbers;
  std::ifstream file(filename);
  if (!file.is_open()) {throw std::string("ERROR: MolecularDynamics.hpp: getNumbersFromFile(): file cannot be open");}
  std::string line;
  while (std::getline(file,line)) {
      std::stringstream ss(line);
      std::string numberStr;
      while (std::getline(ss, numberStr,',')) {
          try {
            int number = std::stoi(numberStr);
            numbers.push_back(number);
          } 
          catch (std::invalid_argument & e) {std::cout << "WARNING: MolecularDynamics.hpp: getNumbersFromFile(): skipping invalid number" << std::endl;}
      }
  }
  file.close();
  return numbers;
}

class MetaDynamics {
  //class for metadynamics objects, based on
  //S. Grimme, J. Chem. Theory Comput., 15, 5, 2847, 2019
protected:
  matrixE gauxREF;
  matrixE gaux0;
  matrixE Rmat;
  matrixE Smat;
  matrixE gRMSD;                                          //RMSD gradients
  std::vector<matrixE> metaset;                           //meta geometries
  double alpha;                                           //general alpha to use everywhere
  double kpush;                                           //general k for pushing
  double kpull;                                           //general k for pulling
  double ppull;                                           //pull force
  double qi[4];                                           //auxiliary
  size_t nmtdstructs;                                     //number of meta-structures (size of metaset)
  std::vector<size_t> metaatoms;                          //list of atoms for which metadynamics applies
  std::vector<double> mtdfactors;                         //kfactors for all declared structures
  std::vector<double> mtdalphas;                          //alphas for all declared structures
  std::vector<double> auxV;
public:
  MetaDynamics(double alp = 1.2, double push = 0.003, double pull = -0.015) {
    setOptions(push,pull,0.05,alp);
    Rmat.resize(3,3);
    Smat.resize(4,4);
  }
  ~MetaDynamics() {}
  //getters
  std::vector<size_t> MetaAtoms() {return metaatoms;}
  std::vector<double> MetaFactors() {return mtdfactors;}
  std::vector<double> MetaAlphas() {return mtdalphas;}
  double Alpha() {return alpha;}
  double kPush() {return kpush;}
  double kPull() {return kpull;}
  double pPull() {return ppull;}
  std::vector<matrixE> MetaSet() {return metaset;}
  matrixE MetaGeometry(size_t istruc) {return metaset[istruc - 1];}
  size_t TotalMetaAtoms() {return metaatoms.size();}
  //setters
  void setOptions(double push = 0.003, double pull = -0.015, double p_pull = 0.05, double alp = 1.2) {
    //general object configurator
    kpush = push;
    kpull = pull;
    alpha = alp;
    ppull = p_pull;
    srand(time(NULL));
    std::cout << "temporarily using fixed seed\n";
    srand(100);
  }
  void setAlpha(double alp) {alpha = alp;}
  void setkPush(double kappa) {kpush = kappa;}
  void setkPull(double kappa) {kpull = kappa;}
  void setpPull(double pparam) {ppull = pparam;}
  void setMetaFactors(const std::vector<double> & newfactors) {mtdfactors = newfactors;}
  void setMetaAlphas(const std::vector<double> & newfactors) {mtdalphas = newfactors;}
  void setMetaAlphas2Alpha(double newalpha = 0.0) {
    if (fabs(newalpha) > 1.0e-12) {alpha = newalpha;}
    for (size_t istruct = 0; istruct < nmtdstructs; ++istruct) {
      mtdalphas[istruct] = alpha;
    }
  }
  void setMetaFactors2kPush(double newk = 0.0) {
    if (fabs(newk) > 1.0e-12) {kpush = newk;}
    for (size_t istruct = 0; istruct < nmtdstructs; ++istruct) {
      mtdfactors[istruct] = kpush;
    }
  }
  void setMetaFactors2kPull(double newk = 0.0) {
    if (fabs(newk) > 1.0e-12) {kpull = newk;}
    for (size_t istruct = 0; istruct < nmtdstructs; ++istruct) {
      mtdfactors[istruct] = kpull;
    }
  }
  //functions for setting the metaset
  void setMetaSet(const std::vector<matrixE> & newset) {
    metaset = newset;
    nmtdstructs = metaset.size();
    mtdfactors.resize(nmtdstructs);
    mtdalphas.resize(nmtdstructs);
    for (size_t idx = 0; idx < nmtdstructs; ++idx) {
      mtdfactors[idx] = kpush;
      mtdalphas[idx] = alpha;
    }
  }
  void setMetaSet2Structure(const matrixE & newgeom) {
    metaset.resize(1);
    metaset[0] = newgeom;
    nmtdstructs = 1;
    mtdfactors.resize(1);
    mtdfactors[0] = kpush;
    mtdalphas.resize(1);
    mtdalphas[0] = alpha;
  }
  void addMetaStructure(const matrixE & newgeom) {
    metaset.push_back(newgeom);
    nmtdstructs = metaset.size();
    mtdfactors.resize(nmtdstructs + 1);
    mtdalphas.resize(nmtdstructs + 1);
    mtdfactors[nmtdstructs] = kpush;
    mtdalphas[nmtdstructs] = alpha;
  }
  void clearMetaSet() {
    metaset.clear();
    mtdfactors.clear();
    mtdalphas.clear();
    nmtdstructs = 0;
  }
  //functions to define list of meta-atoms; these functions are based on index 1 (input), so internally shifted by one
  void setMetaAtoms(const std::vector<size_t> & atmmeta, bool shift = true) {
    //function that takes a vector of atomic indices to add to the meta-atom list
    metaatoms = atmmeta;
    size_t nconstr = metaatoms.size();
    gRMSD.resize(nconstr,3);
    gauxREF.resize(nconstr,3);
    gaux0.resize(nconstr,3);
    if (shift) {
      for (size_t idx = 0; idx < nconstr; ++idx) {
        --metaatoms[idx];
      }
    }
  }
  void addMetaAtom(size_t atmidx) {
    //function that adds atom atmidx - 1 to the meta-atoms list
    metaatoms.push_back(atmidx - 1);
    gRMSD.resize(metaatoms.size(),3);
    gauxREF.resize(metaatoms.size(),3);
    gaux0.resize(metaatoms.size(),3);
  }
  void addMetaElements(const std::string & atmsymbol, const std::vector<size_t> & atoms) {
    //function that adds specific elements (e.g., all protons) to meta-atoms
    int elnumber = Symbol2AtomNr(atmsymbol);
    size_t Natoms = atoms.size();
    for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
      if (atoms[idAtm] == elnumber) {metaatoms.push_back(idAtm);}
    }
    gRMSD.resize(metaatoms.size(),3);
    gauxREF.resize(metaatoms.size(),3);
    gaux0.resize(metaatoms.size(),3);
  }
  void addMetaRange(size_t init, size_t end) {
    //adds atoms between positions init and end, inclusively, to the meta-atom list
    --init;
    int delta = int(end) - int(init);
    int currentsize = metaatoms.size();
    metaatoms.resize(currentsize + delta);
    for (size_t idpos = 0; idpos < delta; ++idpos) {
      metaatoms[currentsize + idpos] = init + idpos;
    }
    gRMSD.resize(currentsize + delta,3);
    gauxREF.resize(currentsize + delta,3);
    gaux0.resize(currentsize + delta,3);
  }
  void clearMetaAtoms() {
    //function that clears the meta-atom list
    metaatoms.clear();
    gRMSD.clear();
    gauxREF.clear();
    gaux0.clear();
  }
  //RMSD functions
  void Quaternion2Rotation(matrixE & Quat, matrixE & Rot, int pos) {
    //function that converts a quaternion in position pos of matrix Quat into a rotation matrix in Rot
    qi[0] = Quat(1,pos);
    qi[1] = Quat(2,pos);
    qi[2] = Quat(3,pos);
    qi[3] = Quat(4,pos);
    for (size_t idr = 1; idr < 5; ++idr) {
      for (size_t idc = idr; idc < 5; ++idc) {
        Quat(idr,idc) = 2.0*qi[idr - 1]*qi[idc - 1];
      }
    }
    Quat(1,1) -= 1.0;
    Rot(1,1) = Quat(1,1) + Quat(2,2);
    Rot(1,2) = Quat(2,3) - Quat(1,4);
    Rot(1,3) = Quat(2,4) + Quat(1,3);
    Rot(2,1) = Quat(2,3) + Quat(1,4);
    Rot(2,2) = Quat(1,1) + Quat(3,3);
    Rot(2,3) = Quat(3,4) - Quat(1,2);
    Rot(3,1) = Quat(2,4) - Quat(1,3);
    Rot(3,2) = Quat(3,4) + Quat(1,2);
    Rot(3,3) = Quat(1,1) + Quat(4,4);
  }
  double RMSD(matrixE & geomA, matrixE & geomB, matrixE & grad) {
    //calculation of least square root-mean-square-deviation between two sets of points, using a quaternion based-method
    double centerA[3];
    double centerB[3];
    double anorm = 0.0;
    double bnorm = 0.0;
    double ermsd;
    size_t nconstraint = geomA.rows();
    //get barycenters, centroidal coordinates
    for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
      centerA[idcoord] = 0.0;
      centerB[idcoord] = 0.0;
      for (size_t idAtm = 0; idAtm < nconstraint; ++idAtm) {
        centerA[idcoord] += geomA(idAtm + 1,idcoord + 1);
        centerB[idcoord] += geomB(idAtm + 1,idcoord + 1);
      }
      centerA[idcoord] /= double(nconstraint);
      centerB[idcoord] /= double(nconstraint);
    }
    //norms
    for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
      for (size_t idAtm = 0; idAtm < nconstraint; ++idAtm) {
        geomA(idAtm + 1,idcoord + 1) -= centerA[idcoord];
        anorm += geomA(idAtm + 1,idcoord + 1)*geomA(idAtm + 1,idcoord + 1);
        geomB(idAtm + 1,idcoord + 1) -= centerB[idcoord];
        bnorm += geomB(idAtm + 1,idcoord + 1)*geomB(idAtm + 1,idcoord + 1);
      }
    }
    //Rmatrix
    for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
      for (size_t jdcoord = 0; jdcoord < 3; ++jdcoord) {
        Rmat(idcoord + 1,jdcoord + 1) = 0.0;
        for (size_t idAtm = 0; idAtm < nconstraint; ++idAtm) {
          Rmat(idcoord + 1,jdcoord + 1) += geomA(idAtm + 1,idcoord + 1)*geomB(idAtm + 1,jdcoord + 1);
        }
      }
    }
    //s matrix
    Smat(1,1) = Rmat(1,1) + Rmat(2,2) + Rmat(3,3);
    Smat(2,1) = Rmat(2,3) - Rmat(3,2);
    Smat(3,1) = Rmat(3,1) - Rmat(1,3);
    Smat(4,1) = Rmat(1,2) - Rmat(2,1);
    Smat(1,2) = Smat(2,1);
    Smat(2,2) = Rmat(1,1) - Rmat(2,2) - Rmat(3,3);
    Smat(3,2) = Rmat(1,2) + Rmat(2,1);
    Smat(4,2) = Rmat(1,3) + Rmat(3,1);
    Smat(1,3) = Smat(3,1);
    Smat(2,3) = Smat(3,2);
    Smat(3,3) = -Rmat(1,1) + Rmat(2,2) - Rmat(3,3);
    Smat(4,3) = Rmat(2,3) + Rmat(3,2);
    Smat(1,4) = Smat(4,1);
    Smat(2,4) = Smat(4,2);
    Smat(3,4) = Smat(4,3);
    Smat(4,4) = -Rmat(1,1) - Rmat(2,2) + Rmat(3,3);
    //get eigenvalues
    auxV = MatDiag(Smat);
    //convert quaternion q to rotation matrix U
    Quaternion2Rotation(Smat,Rmat,4);
    //RMSD + small number to avoid division by 0
    ermsd = sqrt(fmax(0.0,(anorm + bnorm - 2.0*auxV[3]))/double(nconstraint)) + 1.e-9;
    //now gradients
    anorm = dist_Angstrom2aum1/(ermsd*double(nconstraint));
    for (size_t idAtm = 0; idAtm < nconstraint; ++idAtm) {
      for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
        centerA[idcoord] = 0.0;
        for (size_t jdcoord = 0; jdcoord < 3; ++jdcoord) {
          centerA[idcoord] += Rmat(jdcoord + 1,idcoord + 1)*geomB(idAtm + 1,jdcoord + 1);
        }
        grad(idAtm + 1,idcoord + 1) = anorm*(geomA(idAtm + 1,idcoord + 1) - centerA[idcoord]);
      }
    }
    return ermsd;
  }
  void MTDPrepare(const matrixE & geometry) {
    //function to prepare atom constraint if user did not apply any constraint specifically
    size_t natm = geometry.rows();
    size_t nsets = metaset.size();
    if (metaatoms.size() == 0) {
      metaatoms.resize(natm);
      for (size_t idx = 0; idx < natm; ++idx) {
        metaatoms[idx] = idx;
      }
    }
    //else -> user defined it
    gauxREF.resize(metaatoms.size(),3);
    gaux0.resize(metaatoms.size(),3);
    metaset[0] = geometry;
    for (size_t idset = 1; idset < nsets; ++idset) {
      metaset[idset].resize(1,1);
    }
    //randomize structure
    Randomize(metaset[0],1.0e-6);
    gRMSD.resize(natm,3);
  }
  //other functions
  void PreparePathFinder(const matrixE & Reag, const matrixE & Prod, double k_push = 0.0, double k_pull = 0.0, double alph = 0.0, double pp = 0.0) {
    //function that sets everything for metadynamics path search
    metaset.resize(2);
    nmtdstructs = 2;
    metaset[0] = Reag;
    metaset[1] = Prod;
    mtdfactors.resize(2);
    if (fabs(k_push) > 1.0e-12) {kpush = k_push;}
    if (fabs(k_pull) > 1.0e-12) {kpull = k_pull;}
    if (fabs(pp) > 1.0e-12) {ppull = pp;}
    mtdfactors[0] = kpush;
    mtdfactors[1] = kpull;
    mtdalphas.resize(2);
    if (fabs(alph) > 1.0e-12) {alpha = alph;}
    mtdalphas[0] = alpha;
    mtdalphas[1] = alpha;
  }
  void RescalePathFinder(double delta_k, double delta_alpha, bool limitalpha = true, double alpha_lim = 0.2) {
    //function rescaling parameters for path-finder
    kpush *= delta_k;
    kpull *= delta_k;
    mtdfactors[0] = kpush;
    mtdfactors[1] = kpull;
    alpha -= delta_alpha;
    if (limitalpha) {
      if (alpha < alpha_lim) {alpha = alpha_lim;}
    }
    mtdalphas[0] = alpha;
    mtdalphas[1] = alpha;
  }
  void PathFinderResetStructures(const matrixE & Reag, const matrixE & Prod) {
    //to reintroduce the structures
    metaset[0] = Reag;
    metaset[1] = Prod;
  }
  //the main function
  double Metadynamic(const matrixE & geometry, matrixE & gradients) {
    double Ebias = 0.0;
    nmtdstructs = metaset.size();
    if (nmtdstructs > 0) {
      double rmsd;
      double Ei;
      size_t masz = metaatoms.size();
      size_t idAtm;
      size_t icnt;
      for (size_t istruct = 0; istruct < nmtdstructs; ++istruct) {
        //build constrained geometry and the one at stake
        if ((istruct > 0)&&(metaset[istruct].rows() != metaset[istruct - 1].rows())) {break;}
        for (size_t idconstr = 0; idconstr < masz; ++idconstr) {
          idAtm = metaatoms[idconstr];
          for (size_t idcoord = 1; idcoord < 4; ++idcoord) {
            gauxREF(idconstr + 1,idcoord) = metaset[istruct](idAtm + 1,idcoord)*dist_Angstrom2aum1;
            gaux0(idconstr + 1,idcoord) = geometry(idAtm + 1,idcoord)*dist_Angstrom2aum1;
          }
        }
        rmsd = RMSD(gaux0,gauxREF,gRMSD);
        Ei = mtdfactors[istruct]*exp(-mtdalphas[istruct]*rmsd*rmsd);
        Ebias += Ei;
        Ei *= -2.0*mtdalphas[istruct]*rmsd;
        for (size_t idconstr = 0; idconstr < masz; ++idconstr) {
          idAtm = metaatoms[idconstr];
          icnt = 3*idAtm;
          for (size_t idcoord = 1; idcoord < 4; ++idcoord) {
            gradients(icnt + idcoord,1) += Ei*gRMSD(idconstr + 1,idcoord);
          }
        }
      }
    }
    return Ebias;
  }
  void NumericalHessian(matrixE & hessian, const matrixE & geometry, double step = 1.0e-3, bool resetHessian = true) {
    //numerical Hessian with central differences
    int Natoms = geometry.rows();
    matrixE gplus(3*Natoms,1);
    matrixE gminus(3*Natoms,1);
    matrixE stepGeom = geometry;               //matrix containing the actual step
    double invstep = 0.5/step;
    double emtd;
    size_t icounter = 1;
    size_t jcounter = 1;
    if (resetHessian) {
      hessian.resize(3*Natoms,3*Natoms);
      hessian.zero();
    }
    for (size_t idatm = 0; idatm < Natoms; ++idatm) {
      for (size_t idcoord = 0; idcoord < 3; ++idcoord,++icounter) {
        //zero gradients
        gplus.zero();
        gminus.zero();
        //do geometry + step
        stepGeom(idatm + 1,idcoord + 1) = geometry(idatm + 1,idcoord + 1) + step;
        emtd = this->Metadynamic(stepGeom,gplus);
        //do geometry - step
        stepGeom(idatm + 1,idcoord + 1) = geometry(idatm + 1,idcoord + 1) - step;
        emtd = this->Metadynamic(stepGeom,gminus);
        //reset geometry
        stepGeom(idatm + 1,idcoord + 1) = geometry(idatm + 1,idcoord + 1);
        //now Hessian
        jcounter = 1;
        for (size_t idbtm = 0; idbtm < Natoms; ++idbtm) {
          for (size_t jdcoord = 0; jdcoord < 3; ++jdcoord,++jcounter) {
            hessian(icounter,jcounter) += (gplus(jcounter,1) - gminus(jcounter,1))*invstep;
          }
        }
      }
    }
    for (size_t idatm = 0; idatm < 3*Natoms; ++idatm) {
      for (size_t idbtm = 0; idbtm < idatm; ++idbtm) {
        invstep = 0.5*(hessian(idatm + 1,idbtm + 1) + hessian(idbtm + 1,idatm + 1));
        hessian(idatm + 1,idbtm + 1) = invstep;
        hessian(idbtm + 1,idatm + 1) = invstep;
      }
    }
  }
  void Randomize(matrixE & geom, double scale = 1.0) {
    //function that randomizes the geometry prior to metadynamics calculations in order to get the first structure for bias potential
    size_t Natoms = geom.rows();
    for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
      for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
        geom(idAtm + 1,idcoord + 1) += fRandom(-1.0,1.0)*scale;
      }
    }
  }
};
class Dynamics: public MetaDynamics {
  //class for molecular dynamics on a certain macro-structure
  //note that though there is only one geometry, this might include several molecules
  Molecule system;
  matrixE geometry;
  matrixE geometryNew;
  matrixE inertia;
  matrixE velocity;
  matrixE acceleration;
  matrixE velocityAux;
  matrixE velocityN;
  matrixE gradients;
  matrixE block;
  matrixE property;                                       //property to be measured along time
  matrixE propauxmat;                                     //auxiliary matrix for property calculation
  matrixE constraintList;                                 //SHAKE constraints
  matrixE DRorig;
  matrixE DR;
  std::vector<double> masses;
  std::vector<double> Lmoment;
  std::vector<double> propaux;
  std::vector<double> propaux2;
  std::vector<double> auxV;
  double mtdramp;                                         //ramp parameter
  double cmass[3];
  double velcmass[3];
  double velang[3];
  double maxsimulationtime;                               //maximum simulation time for each run (enter in ps, working in a.u.)
  double minequilibrationtime;                            //minimum simulation time for equilibration (enter in ps, working in a.u.)
  double timestep;                                        //timestep (enter in ps, working in a.u.)
  double totalmass;                                       //total mass of the system
  double scaleprotons;                                    //scaling factor for protons so that they don't just dance everywhere
  double timedumpgeom;                                    //frequency at which geometry is dumped (enter in ps, working in a.u.)
  double timedumptraj;                                    //frequency at which trajectories are dumped (enter in ps, working in a.u.)
  //FMCM Add this also to MSMC
  double SCFEaccuracy;                                    //accuracy with which SCF energy is converged
  double SCFDaccuracy;                                    //accuracy with which SCF density is converged
  double temperature;                                     //temperature in Kelvin
  double Telectron;                                       //electronic temperature for fractional occupation numbers in DFT(B) methods
  double Tinit;
  double KBEh;                                            //Boltzmann constant in Hartree
  double Ekin;                                            //kinetic energy
  double Epot;                                            //potential energy
  double Etotal;                                          //total energy
  double Eatm;                                            //energy per atom to check whether something is broken
  double aux;
  std::vector<size_t> atoms;
  matrix<int> neighbours;
  int Natoms;
  int nsteps;                                             //number of steps in integration
  int ndumptraj;                                          //number of steps between trajectory dumps to file
  int ndumpgeom;                                          //number of steps between geometry dumping
  int degreesoffreedom;                                   //how many degrees of freedom allowed to play with
  int SHAKE;                                              //apply SHAKE?
                                                          // 0 -> no
                                                          // 1 -> only on protons
                                                          // 2 -> all atoms
                                                          // 3 -> user defined
  int nmeasurements;                                      //how many measurements were done?
  int measurementfreq;                                    //how frequently should property be measured?
  int blocksz;                                            //size of the block matrix
  int dumpcounter;                                        //counter for dumped geometries
  int nrowprops;
  int szprpaux;
  int mtdmaxstruct;                                       //number of structures considered in metadynamics
  int mtddumpfreq;                                        //frequency with which metadynamics structures are dumped
  bool equilibration;                                     //bool to determine whether to equilibrate
  bool optgeom;                                           //determine whether to optimize geometries at all
  bool dumpgeom;                                          //determines whether to write to file
  bool dumpvelocity;                                      //determines whether to dump velocities as well
  bool restart;                                           //control whether to restart
  bool restartF;                                          //control whether to restart from file
  bool transrotConstr;                                    //translational and rotational constraints
  bool usethermostat;                                     //using thermostat?
  bool measure;                                           //measure properties?
  bool doMetaDyn;                                         //control variable that determines whether to do metadynamics
  bool useContrainPot;                                    //use contrain potential?
  bool forceMD;                                           //this will force the MD run no matter what
  EllipsoidalConstraint ConstrPot;                        //contraining potential
  std::string ensemble;                                   //type of ensemble to be used
                                                          //NVT (canonical) -> default
                                                          //NVE (microcanonical)
  std::string restartFile;                                //file to read restart data from
  std::string restartOutput;                              //file to write at the end of the run to allow for restart from data
  std::string equilibrationFile;                          //file to write trajectories during equilibration into
  std::string trajectoryFile;                             //file to write trajectories into
  std::string geometryFile;                               //file to write geometries into
  std::string thermostat;                                 //type of thermostat to use
  std::string integration;                                //type of integration to use
  std::string propmeasure;                                //property to measure
  std::string auxstrng;
public:
  Dynamics(const Molecule & _system, double _timemax = 50, bool _equilib = true, bool optimizegeom = true, double _timestep = 1.0e-3, double tdumpg = 1.0) {
    setOptions(_timemax,_equilib,optimizegeom,_timestep,tdumpg);
    setNewSystem(_system);
  }
  Dynamics(double _timemax = 50, bool _equilib = true, bool optimizegeom = true, double _timestep = 1.0e-3, double tdumpg = 1.0) {
    setOptions(_timemax,_equilib,optimizegeom,_timestep,tdumpg);
  }
  ~Dynamics() {}
  //getters
  bool Equilibration() {return equilibration;}
  bool GeometryOptBool() {return optgeom;}
  bool GeometryDumpBool() {return dumpgeom;}
  bool VelocityDumpBool() {return dumpvelocity;}
  bool Restart() {return restart;}
  bool RestartFromFile() {return restartF;}
  bool UseThermostat() {return usethermostat;}
  bool Measure() {return measure;}
  bool DoMTD() {return doMetaDyn;}
  bool ConstraintPotential() {return useContrainPot;}
  bool ForceMD() {return forceMD;}
  int DegreesOfFreedom() {return degreesoffreedom;}
  int SHAKEType() {return SHAKE;}
  int NSteps() {return nsteps;}
  int NDumpGeometries() {return ndumpgeom;}
  int NDumpTrajectories() {return ndumptraj;}
  int DumpCounter() {return dumpcounter;}
  int MeasurementFrequency() {return measurementfreq;}
  int MaxNumberMTDStructures() {return mtdmaxstruct;}
  int MTDDumpingFrequency() {return mtddumpfreq;}
  double MaxSimulationTime() {return 0.001*maxsimulationtime*au2fs;}          //returned in ps
  double MinEquilibrationTime() {return 0.001*minequilibrationtime*au2fs;}    //returned in ps
  double TimeStep() {return 0.001*timestep*au2fs;}                            //returned in ps
  double ScaleProtons() {return scaleprotons;}
  double TimeDumpGeometry() {return 0.001*timedumpgeom*au2fs;}                //returned in ps
  double TimeDumpTrajectories() {return 0.001*timedumptraj*au2fs;}            //returned in ps
  double SCFEnergyAccuracy() {return SCFEaccuracy;}
  double SCFDensityAccuracy() {return SCFDaccuracy;}
  double Temperature() {return temperature;}
  double TElectron() {return Telectron;}
  double MTDRamp() {return mtdramp;}
  Molecule System() {return system;}
  std::string EnsembleType() {return ensemble;}
  std::string RestartFile() {return restartFile;}
  std::string RestartOutput() {return restartOutput;}
  std::string TrajectoryFile() {return trajectoryFile;}
  std::string EquilibrationFile() {return equilibrationFile;}
  std::string GeometryFile() {return geometryFile;}
  std::string Thermostat() {return thermostat;}
  std::string Integration() {return integration;}
  std::string PropertyToMeasure() {return propmeasure;}
  matrixE Property() {return property;}
  matrixE ConstraintList() {return constraintList;}
  //setters
  void setOptions(double _tmax, bool _equilib, bool optimizegeom, double _timestep, double tdumpg) {
    maxsimulationtime = 1000.0*_tmax*fs2au;          // <-> tmax
    minequilibrationtime = 5000.0*fs2au;             // <-> mintime
    equilibration = _equilib;
    optgeom = optimizegeom;
    dumpgeom = true;
    dumpvelocity = true;
    timestep = 1000.0*_timestep*fs2au;                      // <-> tstep_md,tstep0
    nsteps = this->calcNstep(maxsimulationtime,timestep);
    timedumpgeom = 1000.0*tdumpg*fs2au;                     // <-> dump_md
    timedumptraj = 1000.0*tdumpg*fs2au;                     // <-> dump_md
    ndumpgeom = this->calcNstep(timedumpgeom,timestep);     // <-> cdump0
    ndumptraj = this->calcNstep(timedumptraj,timestep);     // <-> dumpstep
    scaleprotons = 3.96847889203;
    ensemble = "NVT";                         // <-> thermostat,nvt_md
    usethermostat = true;
    thermostat = "Berendsen";
    SCFEaccuracy = 1.0e-7;                    //well converged
    SCFDaccuracy = 1.0e-4;                    //well-converged
    temperature = 298.15;
    restart = false;
    restartF = false;
    restartFile = "";
    restartOutput = "";
    SHAKE = 0;
    transrotConstr = false;
    useContrainPot = false;
    KBEh = KB/au2J;
    trajectoryFile = "trajectories.trj";
    equilibrationFile = "equilibration.trj";
    geometryFile = "geom";
    dumpcounter = 0;
    measurementfreq = 50;
    Telectron = 300.0;
    inertia.resize(3,3);
    Lmoment.resize(3);
    integration = "Verlet";
    measure = false;
    doMetaDyn = false;
    mtddumpfreq = 2*nsteps;
    //setting seed for random numbers
    srand(time(NULL));
    std::cout << "temporarily using fixed seed\n";
    srand(100);
    forceMD = false;
  }
  void setDumpCounter(int _dumpcounter) {dumpcounter = _dumpcounter;}
  void resetDumpCounter() {dumpcounter = 0;}
  void setMeasurementFrequency(int freq) {measurementfreq = freq;}
  void setEquilibration(bool _equilib) {equilibration = _equilib;}
  void setGeometryOptBool(bool optimizegeom) {optgeom = optimizegeom;}
  void setGeometryDumpBool(bool dumpgeometry) {dumpgeom = dumpgeometry;}
  void setVelocityDumpBool(bool dumpvel) {dumpvelocity = dumpvel;}
  void setRestart(bool rstrt) {restart = rstrt;}
  void setRestartFromFile(bool rstrt) {
    restartF = rstrt;
    restart = rstrt;
  }
  void ApplyContraintPotential(bool constr, double r1 = 1.0, double r2 = -1.0, double r3 = -1.0, int pottype = 0, int beta = 6.0, int alpha = 30.0) {
    //setting up a constraining potential; note that this has nothing to do with the constraints defined in the functions below
    //by default center potential with molecule
    useContrainPot = constr;
    ConstrPot.setRadii(r1,r2,r3);
    std::vector<double> cmass = system.CM();
    ConstrPot.setCenter(cmass[0],cmass[1],cmass[2]);
    ConstrPot.setTemperature(Telectron);
    ConstrPot.setTypeOfPotential(pottype);
    ConstrPot.setExponent(alpha);
    ConstrPot.setBeta(beta);
  }
  void setTransRotConstraint(bool transrot) {ApplyConstraints(transrot,SHAKE);}
  void ApplyConstraints(bool transrot, int SHAKE_ = 0) {
    transrotConstr = transrot;
    SHAKE = SHAKE_;                                  //these are applied later
    degreesoffreedom = 3*Natoms - 6*transrotConstr;
  }
  void setForceMD(bool newfmd) {forceMD = newfmd;}
  void setSHAKE(int _SHAKE) {ApplyConstraints(transrotConstr,_SHAKE);}
  void setMaxSimulationTime(double _timemax) {maxsimulationtime = _timemax*fs2au;}
  void setMinEquilibrationTime(double _timemin) {minequilibrationtime = _timemin*fs2au;}
  void setTimeStep(double _timestep) {timestep = _timestep*fs2au;}
  void setTimeDumpGeometry(double topt) {timedumpgeom = topt*fs2au;}
  void setTimeDumpTrajectory(double tdump) {timedumptraj = tdump*fs2au;}
  void setScaleProtons(double scalefactor) {scaleprotons = scalefactor;}
  void setSCFEnergyAccuracy(double acc) {SCFEaccuracy = acc;}
  void setSCFDensityAccuracy(double acc) {SCFDaccuracy = acc;}
  void setTemperature(double newT) {temperature = newT;}
  void setTElectron(double newT) {Telectron = newT;}
  void setEnsembleType(std::string newensemble) {
    ensemble = newensemble;
    usethermostat = (ensemble == "NVT");
  }
  void setRestartFile(std::string filename) {
    restartFile = filename + ".restart";
    if (filename + "blahblah" != "blahblah") {
      restartF = true;
      restart = true;
    }
    else {
      restartF = false;
      restart = false;
    }
  }
  void setRestartOutput(std::string filename) {restartOutput = filename + ".restart";}
  void setTrajectoryFile(std::string filename) {trajectoryFile = filename + ".trj";}
  void setEquilibrationFile(std::string filename) {equilibrationFile = filename + ".trj";}
  void setGeometryFile(std::string filename) {geometryFile = filename;}
  void setThermostat(std::string newthermo) {
    if ((newthermo == "Berendsen")||(newthermo == "BERENDSEN")||(newthermo == "berendsen")) {thermostat = "Berendsen";}
    else if ((newthermo == "Nose-Hoover")||(newthermo == "NOSE-HOOVER")||(newthermo == "nose-hoover")) {thermostat = "Nose-Hoover";}
    else if ((newthermo == "Maxwell-Boltzmann")||(newthermo == "MAXWELL-BOLTZMANN")||(newthermo == "maxwell-boltzmann")) {thermostat = "Maxwell-Boltzmann";}
    usethermostat = true;
    ensemble = "NVT";
  }
  void RemoveThermostat() {usethermostat = false;}
  void setIntegration(std::string method) {
    if ((method == "Verlet")||(method == "VERLET")||(method == "verlet")) {integration = "Verlet";}
    else if ((method == "LeapFrog")||(method == "LEAPFROG")||(method == "leapfrog")) {integration = "LeapFrog";}
  }
  void setPropertyToMeasure(std::string followprop) {
    if ((followprop == "harmonic-frequencies")||(followprop == "Harmonic-Frequencies")) {
      propmeasure = "harmonic-frequencies";
      nrowprops = 3*Natoms;
    }
    else if ((followprop == "atomic-polarizabilities")||(followprop == "Atomic-Polarizabilities")) {
      propmeasure = "atomic-polarizabilities";
      nrowprops = Natoms;
    }
    else if ((followprop == "partial-charges")||(followprop == "Partial-Charges")) {
      propmeasure = "partial-charges";
      nrowprops = Natoms;
    }
    else if ((followprop == "molecular-polarizabilities")||(followprop == "Molecular-Polarizabilities")) {
      propmeasure = "molecular-polarizabilities";
      propaux.resize(1);
      nrowprops = 1;
    }
    else if ((followprop == "molecular-dispersion")||(followprop == "Molecular-Dispersion")) {
      propmeasure = "molecular-dispersion";
      propaux.resize(2);
      nrowprops = 2;
    }
    else if ((followprop == "heat-of-formation")||(followprop == "Heat-Of-Formation")) {
      propmeasure = "heat-of-formation";
      propaux.resize(1);
      nrowprops = 1;
    }
    else if ((followprop == "ionization-potential")||(followprop == "IP")||(followprop == "ip")) {
      propmeasure = "ip";
      propaux.resize(1);
      nrowprops = 1;
    }
    else if ((followprop == "electron-affinity")||(followprop == "EA")||(followprop == "ea")) {
      propmeasure = "ea";
      propaux.resize(1);
      nrowprops = 1;
    }
    else if ((followprop == "HSAB")||(followprop == "hsab")) {
      propmeasure = "hsab";
      propaux.resize(2);
      nrowprops = 2;
    }
    else if ((followprop == "Fukui-Electrophiles")||(followprop == "Fukui-Indices+")||(followprop == "fukui-indices+")) {
      propmeasure = "fukui+";
      nrowprops = Natoms;
      propaux.resize(Natoms);
    }
    else if ((followprop == "Fukui-Nucleophiles")||(followprop == "Fukui-Indices-")||(followprop == "fukui-indices-")) {
      propmeasure = "fukui-";
      nrowprops = Natoms;
      propaux.resize(Natoms);
    }
    else if ((followprop == "Fukui-Radicals")||(followprop == "Fukui-Indices0")||(followprop == "fukui-indices0")) {
      propmeasure = "fukui0";
      nrowprops = Natoms;
      propaux.resize(Natoms);
    }
    else if ((followprop == "Softness-Electrophiles")||(followprop == "Softness+")||(followprop == "softness+")) {
      propmeasure = "softness+";
      nrowprops = Natoms;
      propaux.resize(Natoms);
    }
    else if ((followprop == "Softness-Nucleophiles")||(followprop == "Softness-")||(followprop == "softness-")) {
      propmeasure = "softness-";
      nrowprops = Natoms;
      propaux.resize(Natoms);
    }
    else if ((followprop == "Softness-Radicals")||(followprop == "Softness0")||(followprop == "softness0")) {
      propmeasure = "softness0";
      nrowprops = Natoms;
      propaux.resize(Natoms);
    }
    measure = true;
  }
  void setPropertyToMeasureOff() {measure = false;}
  void setMTDRamp(double _ramp) {mtdramp = _ramp;}
  void setMaxMTDStructures(int numbstruct) {mtdmaxstruct = numbstruct;}
  void setMTDDumpFrequency(int mtddpfrq) {mtddumpfreq = mtddpfrq;}
  void setMetaDynamics(bool doit, int numbstruct = 10, double k_factor = 0.02, double alpha_exponent = 1.2, int mtddpfrq = 1000, double _ramp = 0.03) {
    //function that sets up the metadynamics calculation
    doMetaDyn = doit;
    mtdmaxstruct = numbstruct;
    kpush = k_factor;
    alpha = alpha_exponent;
    mtdfactors.resize(mtdmaxstruct);
    mtdalphas.resize(mtdmaxstruct);
    for (size_t idx = 0; idx < mtdmaxstruct; ++idx) {
      mtdfactors[idx] = k_factor;
      mtdalphas[idx] = alpha;
    }
    mtdramp = _ramp;
    mtddumpfreq = mtddpfrq;
    metaset.resize(mtdmaxstruct);
  }
  void setNewSystem(const Molecule & newsystem) {
    //function that sets the system's coordinates from an already existing structure
    //it is assumed that the coordinates were already processed by the molecule class
    system = newsystem;
    Natoms = system.Natoms();
    atoms = system.Atoms();
    geometry = system.Geometry();
    system.ConnectivityMatrix(neighbours);    //getting neighbours
    masses.resize(Natoms);
    totalmass = system.AtomicMasses(masses,true,scaleprotons);
    velocity.resize(Natoms,3);
    acceleration.resize(Natoms,3);
    velocityAux.resize(Natoms,3);
    velocityN.resize(Natoms,3);
    degreesoffreedom = 3*Natoms;
  }
  //thermostat functions
  void MaxwellBoltzmann(matrixE & vel, double Tgoal) {
    //Maxwell-Boltzmann thermostat
    //T. A. Andrea, W. C. Swope, H. C. Andersen, J. Chem. Phys., 79, 4576, 1983
    double sigma;
    double velcomp;
    double totalmass = 0.0;
    double vcmass[3];
    vcmass[0] = 0.0;
    vcmass[1] = 0.0;
    vcmass[2] = 0.0;
    for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
      sigma = sqrt(KBEh*Tgoal/masses[idAtm]);
      for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
        velcomp = GaussianRandom(0.0,sigma,1.0e-14);
        vel(idAtm + 1,idcoord + 1) = velcomp;
        vcmass[idcoord] += velcomp*masses[idAtm];
        totalmass += masses[idAtm];
      }
    }
    vcmass[0] /= totalmass;
    vcmass[1] /= totalmass;
    vcmass[2] /= totalmass;
    velcomp = 0.0;
    for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
      for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
        vel(idAtm + 1,idcoord + 1) -= vcmass[idcoord];
        sigma = vel(idAtm + 1,idcoord + 1);
        velcomp += sigma*sigma*masses[idAtm];
      }
    }
    sigma = sqrt(Tgoal*KBEh*double(degreesoffreedom)/velcomp);
    vel *= sigma;
  }
  void Berendsen(matrixE & vel, double step, double tau, double Ti, double Ta) {
    //Berendsen thermostat as declared in 
    //H. J. C. Berendsen, J. P. M. Postma, W. F. van Gunsteren, A. DiNola, J. R. Haak, J. Chem. Phys., 81(8), 3684, 1984
    double Tscal = sqrt(1.0 + (step/tau)*(Ti/Ta - 1.0));
    vel *= Tscal;
  }
  void NoseHoover() {
    //Nose'-Hoover thermostat
  }
  //SHAKE/RATTLE functions; note that RATTLE is a variant of SHAKE where velocities are also corrected
  void AddSHAKEPair(size_t atmA, size_t atmB) {
    //function that adds a particular atom pair to the SHAKE constraint list
    int ccol = constraintList.cols();
    if (ccol < 2) {
      if (constraintList.rows() == 1) {ccol = 0;}
    }
    constraintList.resize(2,ccol + 1);
    constraintList(1,ccol + 1) = atmA;
    constraintList(2,ccol + 1) = atmB;
    SHAKE = 3;
  }
  void RemoveSHAKEPair(size_t atmA, size_t atmB) {
    //function that adds a particular atom pair to the SHAKE constraint list
    size_t ccol = constraintList.cols();
    bool foundpair = false;
    for (size_t idcons = 0; idcons < ccol; ++idcons) {
      if ((atmA == constraintList(1,idcons + 1))&&(atmB == constraintList(2,idcons + 1))) {foundpair = true;}
      else if ((atmA == constraintList(2,idcons + 1))&&(atmB == constraintList(1,idcons + 1))) {foundpair = true;}
      if (foundpair) {
        constraintList(1,idcons + 1) = constraintList(1,ccol);
        constraintList(2,idcons + 1) = constraintList(2,ccol);
        --ccol;
        if (ccol > 0) {constraintList.resize(2,ccol);}
        else {constraintList.clear();}
        break;
      }
    }
    SHAKE = 3;
  }
  void SetSHAKE2Lists(std::vector<size_t> atmsA, std::vector<size_t> atmsB) {
    //function that sets the SHAKE constraint list to predefined lists supplied by the user
    size_t szA = atmsA.size();
    if (szA == atmsB.size()) {
      constraintList.resize(2,szA);
      for (size_t idcons = 0; idcons < szA; ++idcons) {
        constraintList(1,idcons + 1) = atmsA[idcons];
        constraintList(2,idcons + 1) = atmsB[idcons];
      }
      SHAKE = 3;
    }
    else {std::cout << "WARNING: MolecularDynamics.hpp: Dynamics: SetSHAKE2Lists(): user supplied lists are inconsistent\n";}
  }
  void SetSHAKE2List(matrixE & list) {
    //function that sets the SHAKE constraint list to predefined list supplied by the user
    constraintList = list.trans();
    for (size_t idcons = 0; idcons < constraintList.cols(); ++idcons) {
      constraintList(1,idcons + 1) += 1.0;
      constraintList(2,idcons + 1) += 1.0;
    }
    SHAKE = 3;
  }
  template<class T>
  void SHAKEinit(T & ElecStruct, bool simplify = true) {
    //initialization of SHAKE; this structure avoids doing anything for SHAKE = 0
    //simplify applies for SHAKE = 2, where we don't look again over all atoms and simply use the neighbour list already available
    double r2min;
    double rAB2;
    double aux;
    int Bclosest;
    int consrows = 0;
    int idAtm;
    int idBtm;
    size_t nneighbours = neighbours.cols();
    //count the number of protons
    for (idAtm = 0; idAtm < Natoms; ++idAtm) {
      if (atoms[idAtm] == 1) {++consrows;}
    }
    if (SHAKE == 2) {                     //in this case must consider other atoms as well
      Bclosest = consrows;
      //protons -> 1 neighbour; other atoms up to nneighbours
      consrows = Bclosest + (Natoms - Bclosest)*nneighbours;
    }
    else if (SHAKE == 3) {consrows = constraintList.cols();}
    constraintList.resize(2,consrows);
    consrows = 1;
    if (SHAKE == 1) {
      //protons only
      for (idAtm = 0; idAtm < Natoms; ++idAtm) {
        if (atoms[idAtm] != 1) {continue;}
        r2min = 1000.0;
        for (size_t idB = 0; idB < nneighbours; ++idB) {
          idBtm = neighbours(idAtm + 1,idB + 1) - 1;
          if (idBtm < 0) {break;}
          Distance2(idAtm + 1,idBtm + 1,geometry,aux,rAB2);
          if (rAB2 < r2min) {
            r2min = rAB2;
            Bclosest = idBtm;
          }
        }
        if (atoms[Bclosest] == 1) {
          if (Bclosest > idAtm) {
            constraintList(1,consrows) = idAtm + 1;
            constraintList(2,consrows) = Bclosest + 1;
            ++consrows;
          }
        }
        else {
          constraintList(1,consrows) = idAtm + 1;
          constraintList(2,consrows) = Bclosest + 1;
          ++consrows;
        }
      }
    }
    else if (SHAKE == 2) {
      //everyone on board
      bool metalA;
      bool metalB;
      if (!simplify) {                  //this is the way consistent with xTB
        double Rcov;
        double wbothresh;
        //get Mayer/Wiberg bond indices
        matrixE WBO;
        ElecStruct.MayerBondOrder(WBO);
        for (idAtm = 0; idAtm < Natoms; ++idAtm) {
          metalA = MetalHead(atoms[idAtm]);
          if (metalA) {
            if ((atoms[idAtm] != 3)&&(atoms[idAtm] != 4)) {continue;}
          }
          for (idBtm = idAtm + 1; idBtm < Natoms; ++idBtm) {
            metalB = MetalHead(atoms[idBtm]);
            if (metalB) {
              if ((atoms[idBtm] != 3)&&(atoms[idBtm] != 4)) {continue;}
            }
            Distance2(idAtm + 1,idBtm + 1,geometry,aux,rAB2);
            Rcov = AtmRadii(atoms[idAtm]) + AtmRadii(atoms[idBtm]);
            if (rAB2 > 1.44*Rcov*Rcov) {continue;}
            wbothresh = 0.5;
            if ((metalA)||(metalB)) {wbothresh = 0.1;}
            if (WBO(idAtm + 1,idBtm + 1) < wbothresh) {continue;}
            constraintList(1,consrows) = idAtm + 1;
            constraintList(2,consrows) = idBtm + 1;
            ++consrows;
          }
        }
      }
      else {                            //this is a faster way of doing the same
        for (idAtm = 0; idAtm < Natoms; ++idAtm) {
          for (size_t idB = 0; idB < nneighbours; ++idB) {
            idBtm = neighbours(idAtm + 1,idB + 1) - 1;
            if (idBtm < 0) {break;}
            if (idBtm < idAtm) {continue;}
            constraintList(1,consrows) = idAtm + 1;
            constraintList(2,consrows) = idBtm + 1;
            ++consrows;
          }
        }
      }
      constraintList.resize(2,consrows - 1);
    }
    consrows = constraintList.cols();
    DRorig.resize(consrows,4);
    DR.resize(consrows,4);
    for (size_t idpair = 0; idpair < consrows; ++idpair) {
      idAtm = constraintList(1,idpair + 1);
      idBtm = constraintList(2,idpair + 1);
      for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
        DRorig(idpair + 1,idcoord + 1) = geometry(idAtm,idcoord + 1) - geometry(idBtm,idcoord + 1);
      }
      DRorig(idpair + 1,4) = DRorig(idpair + 1,1)*DRorig(idpair + 1,1) + DRorig(idpair + 1,2)*DRorig(idpair + 1,2) + DRorig(idpair + 1,3)*DRorig(idpair + 1,3);
    }
  }
  bool SHAKEit(matrixE & velocity, matrixE & geomA, matrixE & geomB, double SHAKEtol = 1.0e-7, size_t maxiter = 250) {
    //apply SHAKE geometric constraints
    bool converged = false;
    double tau1 = 1.0/timestep;
    double tau2 = tau1*tau1;
    double maxdeviation;
    double dev;
    double delta;
    double mA;
    double mB;
    double lambda;
    size_t nconstraints = constraintList.cols();
    size_t idAtm;
    size_t idBtm;
    for (size_t idpair = 0; idpair < nconstraints; ++idpair) {
      idAtm = constraintList(1,idpair + 1);
      idBtm = constraintList(2,idpair + 1);
      for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
        DRorig(idpair + 1,idcoord + 1) = geomA(idAtm,idcoord + 1) - geomA(idBtm,idcoord + 1);
      }
    }
    geomA = geomB;
    for (size_t ishake = 0; ishake < maxiter; ++ishake) {
      maxdeviation = 0.0;
      for (size_t idpair = 0; idpair < nconstraints; ++idpair) {
        idAtm = constraintList(1,idpair + 1);
        idBtm = constraintList(2,idpair + 1);
        for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
          DR(idpair + 1,idcoord + 1) = geomA(idAtm,idcoord + 1) - geomA(idBtm,idcoord + 1);
        }
        DR(idpair + 1,4) = DR(idpair + 1,1)*DR(idpair + 1,1) + DR(idpair + 1,2)*DR(idpair + 1,2) + DR(idpair + 1,3)*DR(idpair + 1,3);
        delta = DRorig(idpair + 1,4);
        dev = fabs(DR(idpair + 1,4) - delta)/delta;
        if (dev > maxdeviation) {maxdeviation = dev;}
      }
      //is largest deviation below tolerance?
      if (maxdeviation < SHAKEtol) {
        converged = true;
        break;
      }
      for (size_t idpair = 0; idpair < nconstraints; ++idpair) {
        idAtm = constraintList(1,idpair + 1);
        idBtm = constraintList(2,idpair + 1);
        delta = DRorig(idpair + 1,4);
        mA = 1.0/masses[idAtm - 1];
        mB = 1.0/masses[idBtm - 1];
        //get the Lagrange multiplier
        dev = 2.0*(mA + mB)*(DR(idpair + 1,1)*DRorig(idpair + 1,1) + DR(idpair + 1,2)*DRorig(idpair + 1,2) + DR(idpair + 1,3)*DRorig(idpair + 1,3));
        lambda = (delta - DR(idpair + 1,4))/dev;
        for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
          geomA(idAtm,idcoord + 1) += mA*lambda*DRorig(idpair + 1,idcoord + 1);
          geomA(idBtm,idcoord + 1) -= mB*lambda*DRorig(idpair + 1,idcoord + 1);
        }
      }
    }
    if (converged) {
      velocity += (geomA - geomB)*tau1*dist_Angstrom2aum1;
      acceleration += (geomA - geomB)*tau2*dist_Angstrom2aum1;
      geomB = geomA;
    }
    else {std::cout << "WARNING: MolecularDynamics.hpp: Dynamics: SHAKEit(): SHAKE did not converge\n";}
    return converged;
  }
  bool RATTLEit(double tstep, matrixE & velocity, double RATTLEtol = 1.0e-7, size_t maxiter = 250) {
    //apply RATTLE velocity constraints
    //note that positions don't have to be updated, it was already done in SHAKEit (geometric part of RATTLE)
    bool converged = false;
    size_t nconstraints = constraintList.cols();
    size_t idAtm;
    size_t idBtm;
    double maxdeviation;
    double dev;
    double mA;
    double mB;
    double dRdV;
    double mu;
    double delta;
    for (size_t irattle = 0; irattle < maxiter; ++irattle) {
      maxdeviation = 0.0;
      for (size_t idpair = 0; idpair < nconstraints; ++idpair) {
        idAtm = constraintList(1,idpair + 1);
        idBtm = constraintList(2,idpair + 1);
        dRdV = 0.0;
        for (size_t idcoord = 1; idcoord < 4; ++idcoord) {
          dRdV += DR(idpair + 1,idcoord)*(velocity(idAtm,idcoord) - velocity(idBtm,idcoord));
        }
        dev = dRdV*tstep/DR(idpair + 1,4);
        if (dev > maxdeviation) {maxdeviation = dev;}
      }
      //is largest deviation below tolerance?
      if (maxdeviation < RATTLEtol) {
        converged = true;
        break;
      }
      for (size_t idpair = 0; idpair < nconstraints; ++idpair) {
        idAtm = constraintList(1,idpair + 1);
        idBtm = constraintList(2,idpair + 1);
        delta = DRorig(idpair + 1,4);
        mA = 1.0/masses[idAtm - 1];
        mB = 1.0/masses[idBtm - 1];
        dRdV = 0.0;
        for (size_t idcoord = 1; idcoord < 4; ++idcoord) {
          dRdV += DR(idpair + 1,idcoord)*(velocity(idAtm,idcoord) - velocity(idBtm,idcoord));
        }
        //get the Lagrange multiplier
        mu = dRdV/((mA + mB)*DR(idpair + 1,4));
        mA *= mu;
        mB *= mu;
        //update velocities
        for (size_t idcoord = 1; idcoord < 4; ++idcoord) {
          velocity(idAtm,idcoord) -= mA*DR(idpair + 1,idcoord);
          velocity(idBtm,idcoord) += mB*DR(idpair + 1,idcoord);
        }
      }
    }
    return converged;
  }
  //other functions
  int calcNstep(double ttime, double tstep) {
    //function that calculates the number of steps
    return int(ttime/tstep);
  }
  void EstimateTimeStep() {
    //function to estimate a reasonable time step
    double minmass = 100000.0;
    for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
      if (minmass > Weight(atoms[idAtm])) {minmass = Weight(atoms[idAtm]);}
    }
    timestep = fs2au*pow(minmass/Weight(1),1.0/3.0);
  }
  void DumpTrajectory(std::ofstream & gfile, std::string comment = "", bool dumpvel = true, int prc = 7) {
    //function that dumps trajectory into file; the file must be open previously
    if (!gfile.is_open()) {throw std::string("ERROR: MolecularDynamics.hpp: Dynamics: DumpTrajectory(): file for writing not open");}
    double precision = pow(10.0,-prc);
    double number;
    std::string AtmID;
    gfile << std::fixed;
    gfile << std::setprecision(prc);
    //dump geometry
    gfile << Natoms << "\n";
    gfile << comment << "\n";
    for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
      AtmID = AtomNr2Symbol(atoms[idAtm]);
      gfile << AtmID << "    ";
      if (AtmID.length() == 1) {gfile << " ";}
      for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
        number = geometry(idAtm + 1,idcoord + 1);
        if (fabs(number) < precision) {number = 0.0;}
        if (number >= -precision) {gfile << " ";}
        gfile << number << "    ";
      }
      gfile << "\n";
    }
    //dump velocity
    if (dumpvel) {
      for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
        gfile << "    ";
        for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
          number = velocity(idAtm + 1,idcoord + 1);
          if (fabs(number) < precision) {number = 0.0;}
          if (number >= -precision) {gfile << " ";}
          gfile << number << "    ";
        }
        gfile << "\n";
      }
    }
    gfile.flush();
  }
  void ReadRestart(matrixE & vel) {
    //function to read geometries from a given xyz file
    std::ifstream gfile(restartFile, std::ios::in);
    if (!gfile.is_open()) {throw std::string("ERROR: MolecularDynamics.hpp: Dynamics: ReadRestart(): restart file could not be open");}
    std::string read;
    gfile >> Natoms;
    geometry.resize(Natoms,3);
    vel.resize(Natoms,3);
    atoms.resize(Natoms);
    //taking care of the comment
    int atmnr;
    for (size_t idx = 0; idx < 500; ++idx) {
      gfile >> read;
      atmnr = Symbol2AtomNr(read);
      if (atmnr != 0) {
        atoms[0] = atmnr;
        gfile >> geometry(1,1);
        gfile >> geometry(1,2);
        gfile >> geometry(1,3);
        break;
      }
    }
    for (size_t idAtm = 1; idAtm < Natoms; ++idAtm) {             //read geometry and atoms
      gfile >> read;
      atoms[idAtm] = Symbol2AtomNr(read);
      for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
        gfile >> geometry(idAtm + 1,idcoord + 1);
      }
    }
    for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {             //read velocities
      for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
        gfile >> vel(idAtm + 1,idcoord + 1);
      }
    }
    gfile.close();
    system.setGeometry(geometry);
    system.setAtoms(atoms);
  }
  double Ekinetic(const matrixE & velocity, const std::vector<double> & mass) {
    //function that calculates the total kinetic energy of a system
    int Natoms = velocity.rows();
    double auxV[3];
    double energy = 0.0;
    for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
      auxV[0] = velocity(idAtm + 1,1);
      auxV[1] = velocity(idAtm + 1,2);
      auxV[2] = velocity(idAtm + 1,3);
      energy += mass[idAtm]*(auxV[0]*auxV[0] + auxV[1]*auxV[1] + auxV[2]*auxV[2]);
    }
    return 0.5*energy;
  }
  void InitVel(double Ekin) {
    //function that initializes velocities
    double Eatom = 2.0*Ekin/double(3*Natoms);         //this is no temperature, but a kinetic energy per atomic coordinate
    double velA;
    double randomNum;
    double factor;
    double factor2;
    for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
      factor2 = 1.0 + double(atoms[idAtm] == 1);
      velA = sqrt(Eatom/masses[idAtm]);
      randomNum = iRandom01();
      factor = 1.0;
      if (randomNum > 0.5) {factor = -1.0;}
      velocity(idAtm + 1,1) = velA*factor*factor2;
      randomNum = iRandom01();
      factor = 1.0;
      if (randomNum > 0.5) {factor = -1.0;}
      velocity(idAtm + 1,2) = velA*factor*factor2;
      randomNum = iRandom01();
      factor = 1.0;
      if (randomNum > 0.5) {factor = -1.0;}
      velocity(idAtm + 1,3) = velA*factor*factor2;
    }
  }
  void AverageET(double & Eavg, double & Tavg, int sizebl) {
    //function that averages temperatures and potential energies from the block matrix
    Eavg = 0.0;
    Tavg = 0.0;
    for (size_t idx = 0; idx < sizebl; ++idx) {
      Eavg += block(idx + 1,1);
      Tavg += block(idx + 1,2);
    }
    Eavg /= double(sizebl);
    Tavg /= double(sizebl);
  }
  double Regression(int istart, int iend, const std::vector<double> & regressions) {
    //function that calculates the slope of energy to check whether it is drifting
    double numberelements = double(iend - istart) + 1.0;
    double sx = 0.0;
    double sy = 0.0;
    double sxx = 0.0;
    double sxy = 0.0;
    double xx = 0.0;
    for (size_t idx = istart; idx < iend; ++idx) {
      xx += 1.0;
      sx += xx;
      sy += regressions[idx];
      sxx += xx*xx;
      sxy += xx*regressions[idx];
    }
    return (numberelements*sxy - sx*sy)/(numberelements*sxx - sx*sx);
  }
  void RemoveTransRot(const std::vector<double> & masses, matrixE & velocities) {
    //remove translations and rotations from velocities
    for (size_t idx = 0; idx < 3; ++idx) {
      cmass[idx] = 0.0;
      velcmass[idx] = 0.0;
      Lmoment[idx] = 0.0;
      inertia(idx + 1,idx + 1) = 0.0;
      for (size_t idy = 0; idy < idx; ++idy) {
        inertia(idx + 1,idy + 1) = 0.0;
        inertia(idy + 1,idx + 1) = 0.0;
      }
    }
    //get center of mass and the geometry
    aux = 0.0;
    for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
      cmass[0] += geometryNew(idAtm + 1,1)*masses[idAtm];
      cmass[1] += geometryNew(idAtm + 1,2)*masses[idAtm];
      cmass[2] += geometryNew(idAtm + 1,3)*masses[idAtm];
      velcmass[0] += velocities(idAtm + 1,1)*masses[idAtm];
      velcmass[1] += velocities(idAtm + 1,2)*masses[idAtm];
      velcmass[2] += velocities(idAtm + 1,3)*masses[idAtm];
      aux += masses[idAtm];
    }
    cmass[0] /= aux;
    cmass[1] /= aux;
    cmass[2] /= aux;
    //get the angular moment and the inertia matrix
    for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
      //recenter geometry
      geometry(idAtm + 1,1) = dist_Angstrom2aum1*(geometryNew(idAtm + 1,1) - cmass[0]);
      geometry(idAtm + 1,2) = dist_Angstrom2aum1*(geometryNew(idAtm + 1,2) - cmass[1]);
      geometry(idAtm + 1,3) = dist_Angstrom2aum1*(geometryNew(idAtm + 1,3) - cmass[2]);
      Lmoment[0] += masses[idAtm]*(geometry(idAtm + 1,2)*velocities(idAtm + 1,3) - geometry(idAtm + 1,3)*velocities(idAtm + 1,2));
      Lmoment[1] += masses[idAtm]*(geometry(idAtm + 1,3)*velocities(idAtm + 1,1) - geometry(idAtm + 1,1)*velocities(idAtm + 1,3));
      Lmoment[2] += masses[idAtm]*(geometry(idAtm + 1,1)*velocities(idAtm + 1,2) - geometry(idAtm + 1,2)*velocities(idAtm + 1,1));
      inertia(1,1) += masses[idAtm]*(geometry(idAtm + 1,2)*geometry(idAtm + 1,2) + geometry(idAtm + 1,3)*geometry(idAtm + 1,3));
      inertia(2,2) += masses[idAtm]*(geometry(idAtm + 1,3)*geometry(idAtm + 1,3) + geometry(idAtm + 1,1)*geometry(idAtm + 1,1));
      inertia(3,3) += masses[idAtm]*(geometry(idAtm + 1,1)*geometry(idAtm + 1,1) + geometry(idAtm + 1,2)*geometry(idAtm + 1,2));
      inertia(1,2) -= masses[idAtm]*geometry(idAtm + 1,1)*geometry(idAtm + 1,2);
      inertia(1,3) -= masses[idAtm]*geometry(idAtm + 1,1)*geometry(idAtm + 1,3);
      inertia(2,3) -= masses[idAtm]*geometry(idAtm + 1,2)*geometry(idAtm + 1,3);
    }
    inertia(2,1) = inertia(1,2);
    inertia(3,1) = inertia(1,3);
    inertia(3,2) = inertia(2,3);
    //get angular velocity omega
    Solve_Ax_eq_b(inertia,Lmoment);
    for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
      velang[0] = (Lmoment[1]*geometry(idAtm + 1,3) - Lmoment[2]*geometry(idAtm + 1,2));
      velang[1] = (Lmoment[2]*geometry(idAtm + 1,1) - Lmoment[0]*geometry(idAtm + 1,3));
      velang[2] = (Lmoment[0]*geometry(idAtm + 1,2) - Lmoment[1]*geometry(idAtm + 1,1));
      velocities(idAtm + 1,1) -= velcmass[0]/aux + velang[0];
      velocities(idAtm + 1,2) -= velcmass[1]/aux + velang[1];
      velocities(idAtm + 1,3) -= velcmass[2]/aux + velang[2];
    }
  }
  //measurement related functions
  template<class T>
  void Measure(T & ElecStruct, int position) {
    //function measuring a property; always write it to propaux and then process it
    if (property.cols() <= position) {property.resize(nrowprops,position + 20);}
    if (propmeasure == "harmonic-frequencies") {propaux = ElecStruct.CalcVibrFrequencies();}
    else if (propmeasure == "atomic-polarizabilities") {
      if (ElecStruct.Type() != "GFN2") {ElecStruct.getMullikenCharges();}
      propaux2 = ElecStruct.getQAtoms();
      ElecStruct.AtomicPolarizabilities(propaux,propaux2);
    }
    else if (propmeasure == "partial-charges") {
      if (ElecStruct.Type() != "GFN2") {ElecStruct.getMullikenCharges();}
      propaux = ElecStruct.getCharges();
    }
    else if (propmeasure == "molecular-polarizabilities") {
      if (ElecStruct.Type() != "GFN2") {ElecStruct.getMullikenCharges();}
      propaux2 = ElecStruct.getQAtoms();
      ElecStruct.TotalPolarizability(propaux[0],propaux2);
    }
    else if (propmeasure == "molecular-dispersion") {
      if (ElecStruct.Type() != "GFN2") {ElecStruct.getMullikenCharges();}
      propaux2 = ElecStruct.getQAtoms();
      ElecStruct.TotalDispersion(propaux[0],propaux[1],propaux2);
    }
    else if (propmeasure == "heat-of-formation") {
      propaux[0] = ElecStruct.getHeatFormation();
    }
    else if (propmeasure == "ip") {
      propaux[0] = ElecStruct.IonizationPotential(false);
    }
    else if (propmeasure == "ea") {
      propaux[0] = ElecStruct.ElectronAffinity();
    }
    else if (propmeasure == "hsab") {
      ElecStruct.HSABdata(propaux[0],propaux[1]);
    }
    else if (propmeasure == "fukui+") {
      ElecStruct.ReactivityIndices(propauxmat,false);
      for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
        propaux[idAtm] = propauxmat(idAtm + 1,1);
      }
    }
    else if (propmeasure == "fukui-") {
      ElecStruct.ReactivityIndices(propauxmat,false);
      for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
        propaux[idAtm] = propauxmat(idAtm + 1,2);
      }
    }
    else if (propmeasure == "fukui0") {
      ElecStruct.ReactivityIndices(propauxmat,false);
      for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
        propaux[idAtm] = propauxmat(idAtm + 1,3);
      }
    }
    else if (propmeasure == "softness+") {
      ElecStruct.ReactivityIndices(propauxmat,false);
      for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
        propaux[idAtm] = propauxmat(idAtm + 1,4);
      }
    }
    else if (propmeasure == "softness-") {
      ElecStruct.ReactivityIndices(propauxmat,false);
      for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
        propaux[idAtm] = propauxmat(idAtm + 1,5);
      }
    }
    else if (propmeasure == "softness0") {
      ElecStruct.ReactivityIndices(propauxmat,false);
      for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
        propaux[idAtm] = propauxmat(idAtm + 1,6);
      }
    }
    szprpaux = propaux.size();
    for (size_t idcnt = 0; idcnt < szprpaux; ++idcnt) {
      property(idcnt + 1,position) = propaux[idcnt];
    }
  }
  //actual MD functions
  template<class T>
  void runMD(T & ElecStruct, double newtemp = -1.0, double driftthreshold = 2.0e-3, int printfrequency = 200) {
    //function to run the MD simulation
    if (int(newtemp) > 0) {temperature = newtemp;}
    ElecStruct.setElectronTemp(Telectron);
    //read from restart file?
    if ((restart)&&(restartF)) {                    //restart from file
      ReadRestart(velocity);
    }
    //making sure the setup is correct
    ElecStruct.setThresholdDensity(SCFDaccuracy);
    ElecStruct.setThresholdEnergy(SCFEaccuracy);
    ElecStruct.setMolecule(system);
    ElecStruct.setRestart(0);                       //use restart, otherwise this will be a killer
    ElecStruct.Calculate(0,200,true,true,1600.0);
    if (optgeom) {
      //run geometry optimization
      geometry = system.Geometry();
    }
    if (SHAKE > 0) {                                //SHAKE initialization?
      SHAKEinit(ElecStruct);
      degreesoffreedom -= constraintList.cols();
    }
    Tinit = temperature;
    if (equilibration) {
      if (!restart) {                               //initialize velocities
        Ekin = 0.5*(1.0 + double(!usethermostat))*Tinit*KBEh*double(degreesoffreedom);
        InitVel(Ekin);
      }
      //else -> already considered above
      std::cout << "Running equilibration \n";
      if (integration == "Verlet") {Verlet(ElecStruct,100.0*fs2au,true,driftthreshold,printfrequency);}
      else if (integration == "LeapFrog") {LeapFrog(ElecStruct,100.0*fs2au,true,driftthreshold,printfrequency);}
    }
    else {
      if (!restart) {            //initialize velocities
        Ekin = 0.5*(1.0 + double(!usethermostat))*Tinit*KBEh*double(degreesoffreedom);
        InitVel(Ekin);
      }
      //else -> already considered above
    }
    restart = true;
    restartF = false;
    std::cout << "Running productive calculation \n";
    if (integration == "Verlet") {Verlet(ElecStruct,500.0*fs2au,false,driftthreshold,printfrequency);}
    else if (integration == "LeapFrog") {LeapFrog(ElecStruct,500.0*fs2au,false,driftthreshold,printfrequency);}
    if (SHAKE > 0) {this->ApplyConstraints(transrotConstr,SHAKE);}        //reset SHAKE constraints
  }
  template<class T>
  void Verlet(T & ElecStruct, double Taut, bool equilibrate, double driftthreshold = 2.0e-3, int printfrequency = 200) {
    //the actual MD loop function
    //printfrequency is the frequency with which output is printed
    //equilibrate determines whether step is of equilibration
    int gradtype = ElecStruct.AvailableGradients();
    int iblock = 0;
    int nblocks = 0;
    int nregressions = -1;
    int gdump = 1;
    int imeasure = 1;
    int tdump = 1;
    int iprintoutput = 1;
    int icount;
    int mtddump = 0;
    int mtdnstruct = 1;
    auxstrng = trajectoryFile;
    if (equilibrate) {auxstrng = equilibrationFile;}
    std::ofstream trjfile(auxstrng, std::ios::out);
    double Taux = 0.0;
    double aux;
    double blockavgE;
    double blockavgT;
    double slope;
    double Epavg = 0.0;
    double Ekavg = 0.0;
    double Tavg = 0.0;
    double gnorm;
    double Edev;
    double mtdtime = 0.0;
    double Emtd;
    double Ecp;
    size_t szregr = 50;
    std::string strngaux;
    std::vector<double> regressions(szregr,0.0);
    geometry = system.Geometry();
    if (doMetaDyn) {MTDPrepare(geometry);}
    Ekin = Ekinetic(velocity,masses);
    std::cout << "initial MD velocity:\n";
    velocity.Print(12);
    std::cout << "initial MD geometry:\n";
    geometry.Print();
    ElecStruct.Calculate(0,200,true,true,1600.0);
    Epot = ElecStruct.getEnergy(1);
    ElecStruct.gEnergy(gradients,gradtype,0,1.0e-8);
    if (useContrainPot) {
      //constraint potential for the molecule
      Ecp = ConstrPot.Calculate(geometry,gradients);
      Epot += Ecp;
    }
    //get acceleration (F = m.a)
    icount = 1;
    for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
      for (size_t idcoord = 0; idcoord < 3; ++idcoord, ++icount) {
        acceleration(idAtm + 1,idcoord + 1) = -gradients(icount,1)*dist_Angstrom2au/masses[idAtm];
      }
    }
    Etotal = Ekin + Epot;
    nmeasurements = 0;
    std::cout << "nsteps           = " << nsteps << "\n";
    std::cout << "ndumptraj        = " << ndumptraj << "\n";
    std::cout << "ndumpgeom        = " << ndumpgeom << "\n";
    std::cout << "degreesoffreedom = " << degreesoffreedom << "\n";
    std::cout << "SHAKE            = " << SHAKE << "\n";
    std::cout << "nmeasurements    = " << nmeasurements << "\n";
    std::cout << "measurementfreq  = " << measurementfreq << "\n";
    std::cout << "dumpcounter      = " << dumpcounter << "\n";
    blocksz = std::min(5000,int(5000.0*fs2au/timestep));
    block.resize(blocksz,2);
    std::cout << "blocksz          = " << blocksz << "\n";
    //main loop
    for (int idstep = 1; idstep < nsteps + 1; ++idstep) {
      if ((imeasure == measurementfreq)&&(measure)&&(!equilibrate)) {
        //measure whatever
        imeasure = 1;
        ++nmeasurements;
        this->Measure(ElecStruct,nmeasurements);
      }
      if (iblock == blocksz) {                         //average
        if (equilibrate) {                             //increase temperature?
          aux = 1.5*Tinit;
          Tinit = fmin(temperature,aux);
        }
        ++nblocks;
        iblock = 0;
        AverageET(blockavgE,blockavgT,blocksz);
        ++nregressions;
        if (nregressions > szregr) {                   //make sure there is always enough space
          szregr += 50;
          regressions.resize(szregr);
        }
        regressions[nregressions] = blockavgE;
        if (nregressions > 2) {slope = Regression(nregressions - 3,nregressions,regressions);}
        else {slope = 99.0;}
      }
      else {                                           //accumulate
        ++iblock;
        block(iblock,1) = Epot;
        block(iblock,2) = Taux;
      }
      geometry = system.Geometry();
      if ((gdump == ndumpgeom)&&(dumpgeom)) {          //dump geometry?
        system.WriteXYZ(geometryFile,dumpcounter,7);
        gdump = 1;
        ++dumpcounter;
      }
      if (iprintoutput == printfrequency) {            //write details?
        std::cout << idstep << ", " << 0.001*idstep*timestep*au2fs << ", " << Epavg/double(idstep);
        std::cout << ", " << Ekin << ", " << Tavg/double(idstep) << ", " << Taux << ", " << Ekin + Epot;
        std::cout << ", " << Etotal/double(idstep) - Ekin - Epot << "\n";
        iprintoutput = 1;
      }
      if (tdump == ndumptraj) {                        //dump trajectory?
        tdump = 0;
        gnorm = 0.0;
        for (size_t idx = 0; idx < 3*Natoms; ++idx) {
          gnorm += gradients(idx + 1,1)*gradients(idx + 1,1);
        }
        strngaux = "Energy = " + to_string(Epot) + " Eh      |g| = " + to_string(sqrt(gnorm)) + " Eh/A";
        this->DumpTrajectory(trjfile,strngaux,dumpvelocity,7);
      }
      //new positions:         drift                                kick 1/2 step
      geometryNew = geometry + velocity*timestep*dist_Angstrom2au + acceleration*0.5*timestep*timestep*dist_Angstrom2au;
      //Apply SHAKE constraints on positions
      if (SHAKE > 0) {SHAKEit(velocity,geometry,geometryNew,1.0e-7,250);}
      //upd velocities
      velocityN = velocity + acceleration*0.5*timestep;
      //update system and electronic structure
      system.setGeometry(geometryNew);
      ElecStruct.setMolecule(system);
      ElecStruct.Calculate(0,200,true,true,1600.0);
      Epot = ElecStruct.getEnergy(1);
      ElecStruct.gEnergy(gradients,gradtype,0,1.0e-8);
      if (useContrainPot) {
        //constraint potential for the molecule
        Ecp = ConstrPot.Calculate(geometryNew,gradients);
        Epot += Ecp;
      }
      if (doMetaDyn) {                  //do metadynamics
        mtdtime += 1.0;
        for (size_t idx = 0; idx < mtdnstruct; ++idx) {
          mtdfactors[idx] = kpush;
        }
        mtdfactors[mtdnstruct - 1] *= (2.0/(1.0 + exp(-mtdramp*mtdtime)) - 1.0);
        if (mtddump > mtddumpfreq) {
          mtddump = 0;
          if (mtdnstruct < mtdmaxstruct) {          //store next structure
            mtdtime = 0.0;
            ++mtdnstruct;
            metaset[mtdnstruct - 1] = geometry;
          }
          else {                                    //shift all previous structures by one
            for (size_t istruct = 1; istruct < mtdmaxstruct; ++istruct) {
              metaset[istruct - 1] = metaset[istruct];
            }
            metaset[mtdmaxstruct - 1] = geometry;
          }
        }
        Emtd = Metadynamic(geometry,gradients);
        Epot += Emtd;
      }
      Eatm = fabs(Epot/double(Natoms));
      if ((Eatm > 1.0e5)||(Taux > 10000.0)) {
        if (!forceMD) {
          std::cout << "WARNING: MolecularDynamics.hpp: Dynamics: MDcycle(): system's energy is absurd; aborting MD run\n";
          break;
        }
      }
      //get acceleration (F = m.a)
      icount = 1;
      for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
        for (size_t idcoord = 0; idcoord < 3; ++idcoord, ++icount) {
          acceleration(idAtm + 1,idcoord + 1) = -gradients(icount,1)*dist_Angstrom2au/masses[idAtm];
        }
      }
      //kinetic energy
      Ekin = Ekinetic(velocityN,masses);
      //temperature
      Taux = 2.0*Ekin/(KBEh*double(degreesoffreedom));
      if (usethermostat) {                                      //thermostat?
        if (thermostat == "Berendsen") {this->Berendsen(velocityN,timestep,Taut,Tinit,Taux);}
        else if (thermostat == "Nose-Hoover") {}
        else if (thermostat == "Maxwell-Boltzmann") {this->MaxwellBoltzmann(velocityN,Tinit);}
      }
      //upd velocity
      velocity = velocityN + 0.5*acceleration*timestep;           //vmiddle
      //apply RATTLE constraints on velocities
      if (SHAKE > 0) {RATTLEit(timestep,velocity,1.0e-7,250);}
      //leave only internal motion
      RemoveTransRot(masses,velocity);
      //kinetic energy
      Ekin = Ekinetic(velocity,masses);
      //temperature
      Taux = 2.0*Ekin/(KBEh*double(degreesoffreedom));
      ++tdump;
      ++gdump;
      ++imeasure;
      ++iprintoutput;
      Etotal += Epot + Ekin;
      Edev = Etotal/double(idstep) - Epot - Ekin;
      Tavg += Taux;
      Epavg += Epot;
      Ekavg += Ekin;
      if ((equilibrate)&&(0.001*double(idstep)*timestep > minequilibrationtime)) {
        //can we prematurely exit?
        if ((nblocks > 1)&&(fabs(blockavgT - Tinit)/Tinit < 0.02)&&(fabs(slope) < driftthreshold)) {
          std::cout << "early equilibration of the system; prematurely exiting the equilibration\n";
          break;
        }
      }
    }
    geometry = system.Geometry();
    std::cout << "final MD geometry:\n";
    geometry.Print();
    std::cout << "final MD velocities:\n";
    velocity.Print(10);
    std::cout << "average properties " << std::endl;
    std::cout << "Epot               :" << Epavg/double(nsteps) << std::endl;
    std::cout << "Ekin               :" << Ekavg/double(nsteps) << std::endl;
    std::cout << "Etot               :" << (Ekavg + Epavg)/double(nsteps) << std::endl;
    std::cout << "T                  :" << Tavg/double(nsteps) << std::endl;
    if ((fabs(Tavg/double(nsteps) - Tinit) > 0.02*Tinit)&&(nsteps > 500)&&(usethermostat)&&(!equilibrate)) {std::cout << "WARNING: MolecularDynamics.hpp: Dynamics: MDcycle(): problem with thermostat" << std::endl;}
    trjfile.close();
    if (restartOutput + "blahblah" != "blahblah") {           //write last state to file for restart from file
      trjfile.open(restartOutput, std::ios::out);
      DumpTrajectory(trjfile,"",true,10);
      trjfile.close();
    }
  }
  template<class T>
  void LeapFrog(T & ElecStruct, double Taut, bool equilibrate, double driftthreshold = 2.0e-3, int printfrequency = 200) {
    //the actual MD loop function
    //printfrequency is the frequency with which output is printed
    //equilibrate determines whether step is of equilibration
    int gradtype = ElecStruct.AvailableGradients();
    int iblock = 0;
    int nblocks = 0;
    int nregressions = -1;
    int gdump = 1;
    int imeasure = 1;
    int tdump = 1;
    int iprintoutput = 1;
    int mtddump = 0;
    int icount;
    int mtdnstruct = 1;
    auxstrng = trajectoryFile;
    if (equilibrate) {auxstrng = equilibrationFile;}
    std::ofstream trjfile(auxstrng, std::ios::out);
    double Taux = 0.0;
    double aux;
    double blockavgE;
    double blockavgT;
    double slope;
    double Epavg = 0.0;
    double Ekavg = 0.0;
    double Tavg = 0.0;
    double gnorm;
    double Edev;
    double mtdtime = 0.0;
    double Emtd;
    double Ecp;
    size_t szregr = 50;
    std::string strngaux;
    std::vector<double> regressions(szregr,0.0);
    geometry = system.Geometry();
    Ekin = Ekinetic(velocity,masses);
    std::cout << "initial MD velocity:\n";
    velocity.Print(12);
    std::cout << "initial MD geometry:\n";
    geometry.Print();
    ElecStruct.Calculate(0,200,true,true,1600.0);
    Epot = ElecStruct.getEnergy(1);
    Etotal = Ekin + Epot;
    nmeasurements = 0;
    std::cout << "nsteps           = " << nsteps << "\n";
    std::cout << "ndumptraj        = " << ndumptraj << "\n";
    std::cout << "ndumpgeom        = " << ndumpgeom << "\n";
    std::cout << "degreesoffreedom = " << degreesoffreedom << "\n";
    std::cout << "SHAKE            = " << SHAKE << "\n";
    std::cout << "nmeasurements    = " << nmeasurements << "\n";
    std::cout << "measurementfreq  = " << measurementfreq << "\n";
    std::cout << "dumpcounter      = " << dumpcounter << "\n";
    std::cout << "mtddpfrq         = " << mtddumpfreq << "\n";
    blocksz = std::min(5000,int(5000.0*fs2au/timestep));
    block.resize(blocksz,2);
    std::cout << "blocksz          = " << blocksz << "\n";
    if (doMetaDyn) {MTDPrepare(geometry);}
    //main loop
    for (int idstep = 1; idstep < nsteps + 1; ++idstep) {
      ElecStruct.Calculate(0,200,true,true,1600.0);
      Epot = ElecStruct.getEnergy(1);
      geometry = system.Geometry();
      ElecStruct.gEnergy(gradients,gradtype,0,1.0e-8);
      if (useContrainPot) {
        //constraint potential for the molecule
        Ecp = ConstrPot.Calculate(geometry,gradients);
        Epot += Ecp;
      }
      if (doMetaDyn) {                  //do metadynamics
        mtdtime += 1.0;
        for (size_t idx = 0; idx < mtdnstruct; ++idx) {
          mtdfactors[idx] = kpush;
        }
        mtdfactors[mtdnstruct - 1] *= (2.0/(1.0 + exp(-mtdramp*mtdtime)) - 1.0);
        if (mtddump > mtddumpfreq) {
          mtddump = 0;
          if (mtdnstruct < mtdmaxstruct) {          //store next structure
            mtdtime = 0.0;
            ++mtdnstruct;
            metaset[mtdnstruct - 1] = geometry;
          }
          else {                                    //shift all previous structures by one
            for (size_t istruct = 1; istruct < mtdmaxstruct; ++istruct) {
              metaset[istruct - 1] = metaset[istruct];
            }
            metaset[mtdmaxstruct - 1] = geometry;
          }
        }
        
        Emtd = Metadynamic(geometry,gradients);
        Epot += Emtd;
      }
      if ((imeasure == measurementfreq)&&(measure)&&(!equilibrate)) {
        //measure whatever
        imeasure = 1;
        ++nmeasurements;
        this->Measure(ElecStruct,nmeasurements);
      }
      Eatm = fabs(Epot/double(Natoms));
      if ((Eatm > 1.0e5)||(Taux > 10000.0)) {
        if (!forceMD) {
          std::cout << "WARNING: MolecularDynamics.hpp: Dynamics: MDcycle(): system's energy is absurd; aborting MD run\n";
          break;
        }
      }
      if (iblock == blocksz) {                         //average
        if (equilibrate) {                             //increase temperature?
          aux = 1.5*Tinit;
          Tinit = fmin(temperature,aux);
        }
        ++nblocks;
        iblock = 0;
        AverageET(blockavgE,blockavgT,blocksz);
        ++nregressions;
        if (nregressions > szregr) {                   //make sure there is always enough space
          szregr += 50;
          regressions.resize(szregr);
        }
        regressions[nregressions] = blockavgE;
        if (nregressions > 2) {slope = Regression(nregressions - 3,nregressions,regressions);}
        else {slope = 99.0;}
      }
      else {                                           //accumulate
        ++iblock;
        block(iblock,1) = Epot;
        block(iblock,2) = Taux;
      }
      if ((gdump == ndumpgeom)&&(dumpgeom)) {          //dump geometry?
        system.WriteXYZ(geometryFile,dumpcounter,7);
        gdump = 1;
        ++dumpcounter;
      }
      if (iprintoutput == printfrequency) {            //write details?
        std::cout << idstep << ", " << 0.001*idstep*timestep*au2fs << ", " << Epavg/double(idstep);
        std::cout << ", " << Ekin << ", " << Tavg/double(idstep) << ", " << Taux << ", " << Ekin + Epot;
        std::cout << ", " << Etotal/double(idstep) - Ekin - Epot << "\n";
        iprintoutput = 1;
      }
      if (tdump == ndumptraj) {                        //dump trajectory?
        tdump = 0;
        gnorm = 0.0;
        for (size_t idx = 0; idx < 3*Natoms; ++idx) {
          gnorm += gradients(idx + 1,1)*gradients(idx + 1,1);
        }
        strngaux = "Energy = " + to_string(Epot) + " Eh      |g| = " + to_string(sqrt(gnorm)) + " Eh/A";
        this->DumpTrajectory(trjfile,strngaux,dumpvelocity,7);
      }
      //get acceleration (F = m.a)
      icount = 1;
      for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
        for (size_t idcoord = 0; idcoord < 3; ++idcoord, ++icount) {
          acceleration(idAtm + 1,idcoord + 1) = -gradients(icount,1)*dist_Angstrom2au/masses[idAtm];
        }
      }
      //upd velocities
      velocityN = velocity + acceleration*0.5*timestep;
      //kinetic energy
      Ekin = Ekinetic(velocityN,masses);
      //temperature
      Taux = 2.0*Ekin/(KBEh*double(degreesoffreedom));
      velocityAux = velocity + acceleration*timestep;
      if (usethermostat) {                                      //thermostat?
        if (thermostat == "Berendsen") {this->Berendsen(velocityAux,timestep,Taut,Tinit,Taux);}
        else if (thermostat == "Nose-Hoover") {}
        else if (thermostat == "Maxwell-Boltzmann") {this->MaxwellBoltzmann(velocityAux,Tinit);}
      }
      //positions at t + dt
      geometryNew = geometry + velocityAux*timestep*dist_Angstrom2au;
      //apply SHAKE constraints
      if (SHAKE > 0) {SHAKEit(velocityAux,geometry,geometryNew,1.0e-7,250);}
      velocity = velocityAux;
      //leave only internal motion
      RemoveTransRot(masses,velocity);
      //update system and electronic structure
      system.setGeometry(geometryNew);
      ElecStruct.setMolecule(system);
      ++tdump;
      ++gdump;
      ++imeasure;
      ++iprintoutput;
      ++mtddump;
      Etotal += Epot + Ekin;
      Edev = Etotal/double(idstep) - Epot - Ekin;
      Tavg += Taux;
      Epavg += Epot;
      Ekavg += Ekin;
      if ((equilibrate)&&(0.001*double(idstep)*timestep > minequilibrationtime)) {
        //can we exit prematurely?
        if ((nblocks > 1)&&(fabs(blockavgT - Tinit)/Tinit < 0.02)&&(fabs(slope) < driftthreshold)) {
          std::cout << "early equilibration of the system; prematurely exiting the equilibration\n";
          break;
        }
      }
    }
    geometry = system.Geometry();
    std::cout << "final MD geometry:\n";
    geometry.Print();
    std::cout << "final MD velocities:\n";
    velocity.Print(10);
    std::cout << "average properties " << std::endl;
    std::cout << "Epot               :" << Epavg/double(nsteps) << std::endl;
    std::cout << "Ekin               :" << Ekavg/double(nsteps) << std::endl;
    std::cout << "Etot               :" << (Ekavg + Epavg)/double(nsteps) << std::endl;
    std::cout << "T                  :" << Tavg/double(nsteps) << std::endl;
    if ((fabs(Tavg/double(nsteps) - Tinit) > 0.02*Tinit)&&(nsteps > 500)&&(usethermostat)&&(!equilibrate)) {std::cout << "WARNING: MolecularDynamics.hpp: Dynamics: MDcycle(): problem with thermostat" << std::endl;}
    trjfile.close();
    if (restartOutput + "blahblah" != "blahblah") {           //write last state to file for restart from file
      trjfile.open(restartOutput, std::ios::out);
      DumpTrajectory(trjfile,"",true,10);
      trjfile.close();
    }
  }
};
class CompositeMTDQM {
  //class that merges a metadynamics object with a quantum chemical method so that the main solver may be used with composite methods
  //this is not supposed to configure any object, it simply couples two objects for calculation purposes
  MetaDynamics * MTDobj;                //pointer to metadynamics class
  Method * ElectronicStr;              //pointer to quantum chemical method
  double energy;
  matrixE geometry;
  matrixE aux;
  matrixE auxHess;
public:
  CompositeMTDQM(MetaDynamics & MTD, Method & EStrct, const matrixE & geom) {
    MTDobj = & MTD;
    ElectronicStr = & EStrct;
    geometry = geom;
    aux.resize(3*geom.rows(),1);
  }
  ~CompositeMTDQM() {}
  std::string Type() {return "CompositeMTDQM";}
  std::vector<matrixE> MetaSet() {return MTDobj->MetaSet();}
  size_t NAtoms() {return ElectronicStr->NAtoms();}
  Molecule & Component() {return ElectronicStr->Component();}
  double getEnergy(bool complete = 0) {
    energy = ElectronicStr->getEnergy(complete);
    energy += MTDobj->Metadynamic(geometry,aux);
    return energy;
  }
  void gEnergy(matrixE & gx, int type = 0, bool project = 0, double threshold = 1.0e-7) {
    //first plain gradients without projection
    ElectronicStr->gEnergy(aux,type,false,threshold);
    energy = MTDobj->Metadynamic(geometry,aux);
    gx = aux;
    if (project) {ElectronicStr->ProjectGradients(gx);}
  }
  void hEnergy(matrixE & hx, int type = 0, bool project = 0, bool setrestart1 = false, double damp = 0.03, double threshold = 1.0e-7, int maxiter = 100) {
    ElectronicStr->hEnergy(hx,type,0,setrestart1,damp,threshold,maxiter);
    MTDobj->NumericalHessian(auxHess,geometry,1.0e-3,true);
    hx += auxHess;
    if (project) {
      auxHess = geometry;
      ProjectForceConstants(hx,auxHess);
      auxHess = geometry;
      ElectronicStr->WilsonBmatrix(hx,auxHess,damp);              //this modifies the Hessian as well -> diagonal matrix
    }
  }
  std::vector<double> CalcVibrFrequencies(bool restartcalcs = true) {
    //wrapper to get vibrational frequencies without headaches
    this->hEnergy(auxHess,0,0,restartcalcs,0.0,1.0e-7,50);
    std::vector<double> Vfreq = ElectronicStr->Vibrations(auxHess,1,1,true);
    return Vfreq;
  }
  void setGeometry(matrixE & _x) {
    geometry = _x;
    ElectronicStr->setGeometry(_x);
  }
  void updGeometry(matrixE & step, bool project = 0) {
    ElectronicStr->updGeometry(step,project);
    geometry = ElectronicStr->Geometry();
  }
  matrixE Geometry() {return geometry;}
  void Calculate(int _prnt = 0, size_t maxiter = 200, bool _DIIS = true, bool _RCA = true, double d2threshold = 1600.0) {
    ElectronicStr->Calculate(_prnt,maxiter,_DIIS,_RCA,d2threshold);
  }
  std::vector<size_t> Atoms() {return ElectronicStr->Atoms();}
};
template <class slvr> class PathFinder {
  //class for determining reaction path based on metadynamics, based on
  //S. Grimme, J. Chem. Theory Comput., 15, 5, 2847, 2019
  slvr solver;                                            //solver method chosen to optimize geometry
  MetaDynamics MTD;
  Molecule reagent;
  Molecule product;
  size_t Natoms;
  size_t nrefinements;
  size_t MaxNpoints;                                      //maximum number of points per path
  size_t tsindex;                                         //index number of estimated transition state structure
  size_t startindex;                                      //if given geometry is not the minimum on the reagent side, update this
  int hessiantype;
  int gradtype;
  double aline;                                           //alpha for the line search/trust region
  double Telectron;                                       //electronic temperature for fractional occupation numbers in DFT(B) methods
  double alpha;                                           //general alpha to use everywhere
  double kpush;                                           //k for pushing
  double kpull;                                           //k for pulling
  double ppull;                                           //pull force
  double deltak;                                          //increase factor for kpush and kpull between different runs (multiplicative)
  double deltaalpha;                                      //increase factor for alpha between different runs (additive or shift)
  double SCFEaccuracy;                                    //accuracy with which SCF energy is converged
  double SCFDaccuracy;                                    //accuracy with which SCF density is converged
  double rmsdthreshold_product;                           //RMSD-distance threshold from product to terminate scan
  double econvthreshold;                                  //energy convergence threshold
  double gconvthreshold;                                  //gradient convergence threshold
  std::vector<size_t> atoms;
  std::vector<size_t> molvector;
  std::vector<double> inertiaR;
  std::vector<double> inertiaP;
  std::vector<double> bestpathRMSD;
  std::vector<double> bestpathRMSDtemp;
  std::vector<double> bestpathenergies;                   //stored energies from the most promising path found in run
  std::vector<double> bestpathgradient;                   //stored estimated gradients from the most promising path found in run
  std::vector<matrixE> bestpathRP;                        //the geometries associated to each point in the previous vector
public:
  PathFinder(const Molecule & reag, const Molecule & prod, double alp = 1.2, double push = 0.003, double pull = -0.015) {
    setNewSystem(reag,prod);
    setOptions(alp,push,pull,0.05);
  }
  PathFinder(double alp = 1.2, double push = 0.003, double pull = -0.015) {
    setOptions(alp,push,pull,0.05);
  }
  ~PathFinder() {}
  //getters
  Molecule Reagent() {return reagent;}
  Molecule Product() {return product;}
  double LineSearch() {return aline;}
  double Alpha() {return alpha;}
  double kPush() {return kpush;}
  double kPull() {return kpull;}
  double pPull() {return ppull;}
  double SCFEnergyAccuracy() {return SCFEaccuracy;}
  double SCFDensityAccuracy() {return SCFDaccuracy;}
  double EnergyConvergenceThreshold() {return econvthreshold;}
  double GradientConvergenceThreshold() {return gconvthreshold;}
  double DeltaK() {return deltak;}
  double DeltaAlpha() {return deltaalpha;}
  double ThresholRMSDProduct() {return rmsdthreshold_product;}
  double ElectronTemperature() {return Telectron;}
  int HessianType() {return hessiantype;}
  size_t NumberRefinementSteps() {return nrefinements;}
  size_t MaxNumberPointsPerPath() {return MaxNpoints;}
  size_t TSIndex() {return tsindex;}
  size_t StartIndex() {return startindex;}
  std::vector<size_t> MoleculeVector() {return molvector;}
  std::vector<size_t> MetaAtoms() {return MTD.MetaAtoms();}
  std::vector<double> OptPathEnergies() {return bestpathenergies;}
  std::vector<double> OptPathRMSD() {return bestpathRMSD;}
  std::vector<double> OptPathGradients() {return bestpathgradient;}
  std::vector<matrixE> OptPathStructures() {return bestpathRP;}
  void WritePath2XYZ() {
    Molecule auxie;
    auxie.setAtoms(atoms);
    size_t szpath = bestpathRP.size();
    for (size_t idpath = 0; idpath < szpath; ++idpath) {
      //assume that the atom's list is unchanged from the calculation of the best path
      auxie.setGeometry(bestpathRP[idpath]);
      auxie.WriteXYZ("bestpathg",idpath + 1);
    }
  }
  //setters
  void setOptions(double alp = 1.2, double push = 0.003, double pull = -0.015, double p_pull = 0.05) {
    //general object configurator
    kpush = push;
    kpull = pull;
    alpha = alp;
    ppull = p_pull;
    nrefinements = 1;
    hessiantype = 0;
    MaxNpoints = 25;
    rmsdthreshold_product = 0.2;
    deltak = 1.25;
    deltaalpha = 0.2;
    SCFEaccuracy = 1.0e-7;                    //well converged
    SCFDaccuracy = 1.0e-4;                    //well-converged
    econvthreshold = 5.0e-6;                  //Hartree
    gconvthreshold = 5.0e-3;                  //Eh/A
    aline = 0.05;
  }
  void setLineSearch(double alin) {aline = alin;}
  void setAlpha(double alp) {alpha = alp;}
  void setkPush(double kappa) {kpush = kappa;}
  void setkPull(double kappa) {kpull = kappa;}
  void setpPull(double pparam) {ppull = pparam;}
  void setSCFEnergyAccuracy(double acc) {SCFEaccuracy = acc;}
  void setSCFDensityAccuracy(double acc) {SCFDaccuracy = acc;}
  void setEnergyConvergenceThreshold(double thresh) {econvthreshold = thresh;}
  void setGradientConvergenceThreshold(double thresh) {gconvthreshold = thresh;}
  void setRMSDThresholdProduct(double trmsd) {rmsdthreshold_product = trmsd;}
  void setDeltaK(double delta) {deltak = delta;}
  void setDeltaAlpha(double delta) {deltaalpha = delta;}
  void setElectronTemperature(double Te) {Telectron = Te;}
  void setHessianType(int typeH) {hessiantype = typeH;}
  void setNumberRefinementSteps(size_t npathopt) {nrefinements = npathopt;}
  void setMaxNumberPointsPerPath(size_t maxn) {MaxNpoints = maxn;}
  void setNewSystem(const Molecule & reag, const Molecule & prod) {
    //function that prepares reagent and product information
    reagent = reag;
    product = prod;
    Natoms = reag.Natoms();
    molvector.resize(Natoms);
    if (Natoms != product.Natoms()) {std::cout << "WARNING: MolecularDynamics.hpp: PathFinder: setNewSystem(): reagent and product do not have the same number of atoms\n";}
    atoms = reag.Atoms();
    bool match = true;
    //check whether atom list matches
    for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
      if (atoms[idAtm] != product.Atom(idAtm + 1)) {match = false;}
      if (!match) {
        std::cout << "WARNING: MolecularDynamics.hpp: PathFinder: setNewSystem(): reagent and product do not have the same atom list\n";
        break;
      }
    }
  }
  void setReagent(const Molecule & reag) {
    //function that sets the reagent
    reagent = reag;
    Natoms = reagent.Natoms();
    molvector.resize(Natoms);
    if (Natoms != product.Natoms()) {std::cout << "WARNING: MolecularDynamics.hpp: PathFinder: setReagent(): reagent and product do not have the same number of atoms\n";}
    atoms = reagent.Atoms();
    bool match = true;
    //check whether atom list matches
    for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
      if (atoms[idAtm] != product.Atom(idAtm + 1)) {match = false;}
      if (!match) {
        std::cout << "WARNING: MolecularDynamics.hpp: PathFinder: setReagent(): reagent and product do not have the same atom list\n";
        break;
      }
    }
  }
  void setProduct(const Molecule & prod) {
    //function that prepares the product
    product = prod;
    if (Natoms != product.Natoms()) {std::cout << "WARNING: MolecularDynamics.hpp: PathFinder: setProduct(): reagent and product do not have the same number of atoms\n";}
    bool match = true;
    //check whether atom list matches
    for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
      if (atoms[idAtm] != product.Atom(idAtm + 1)) {match = false;}
      if (!match) {
        std::cout << "WARNING: MolecularDynamics.hpp: PathFinder: setNewSystem(): reagent and product do not have the same atom list\n";
        break;
      }
    }
  }
  void setMetaAtoms(const std::vector<size_t> & atmmeta, bool shift = true) {
    //function that takes a vector of atomic indices to add to the meta-atom list
    MTD.setMetaAtoms(atmmeta,shift);
  }
  void addMetaAtom(size_t atmidx) {
    //function that adds atom atmidx - 1 to the meta-atoms list
    MTD.addMetaAtom(atmidx);
  }
  void addMetaElements(const std::string & atmsymbol) {
    //function that adds specific elements (e.g., all protons) to meta-atoms
    MTD.addMetaElements(atmsymbol,atoms);
  }
  void addMetaRange(size_t init, size_t end) {
    //adds atoms between positions init and end, inclusively, to the meta-atom list
    MTD.addMetaRange(init,end);
  }
  void clearMetaAtoms() {
    //function that clears the meta-atom list
    MTD.clearMetaAtoms();
  }
  //other functions
  int MoleculeCounter(Molecule & structure) {
    int counter = structure.MoleculeCounter(molvector);
    return counter;
  }
  void PrintReactionProfile(double re, double te, double pe, bool processed, std::string header, int prc = 7) {
    //printing information on the reaction
    std::cout << std::setprecision(prc);
    if (tsindex == 0) {std::cout << "WARNING: MolecularDynamics.hpp: PathFinder: PrintReactionProfile(): no transition state found\n";}
    double Erx;
    double ETSforward;
    double ETSbackward;
    if (processed) {                  //if the data has already been processed, look here to know the order of arguments
      Erx = re*au2kcalmol;
      ETSforward = te*au2kcalmol;
      ETSbackward = pe*au2kcalmol;
    }
    else {
      Erx = (pe - re)*au2kcalmol;
      ETSforward = (te - re)*au2kcalmol;
      ETSbackward = (te - pe)*au2kcalmol;
    }
    std::cout << header << "\n";
    std::cout << "E(reaction) = " << Erx << " kcal/mol \n";
    std::cout << "Ea(direct)  = " << ETSforward << " kcal/mol \n";
    std::cout << "Ea(reverse) = " << ETSbackward << " kcal/mol \n";
  }
  //path finder functions
  void PathCleaner(std::vector<matrixE> & thepath, std::vector<double> & Epath, std::vector<double> & Eaux, size_t maxnpoints) {
    //function that removes points along the path that are not really needed
    size_t lengthpath = thepath.size();
    double dgeom;
    double de;
    double thresholdA = 0.0001;
    double thresholdB = 0.00001;
    double tsenergy;
    size_t NA = thepath[0].rows();
    size_t tsposition;
Remove_More:
    tsenergy = Epath[0];
    for (size_t ipath = 1; ipath < lengthpath - 1; ++ipath) {
      dgeom = 1.0;
      for (size_t idAtm = 0; idAtm < NA; ++idAtm) {
        for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
          de = (thepath[ipath](idAtm + 1,idcoord + 1) - thepath[ipath - 1](idAtm + 1,idcoord + 1))*dist_Angstrom2aum1;
          dgeom += de*de;
        }
      }
      if (Epath[ipath] > tsenergy) {       //store the position of the transition state not to delete this
        tsposition = ipath;
        tsenergy = Epath[ipath];
      }
      de = Epath[ipath] - Epath[ipath - 1];
      Eaux[ipath] = de/dgeom;
    }
    for (size_t ipath = 1; ipath < lengthpath - 1; ++ipath) {
      if ((fabs(Eaux[ipath]) < thresholdA)&&(ipath < tsposition)) {
        Eaux[ipath + 1] += Eaux[ipath];
        thepath.erase(thepath.begin() + ipath);
        Epath.erase(Epath.begin() + ipath);
        Eaux.erase(Eaux.begin() + ipath);
        --lengthpath;
        --ipath;
      }
      else if ((fabs(Eaux[ipath]) < thresholdB)&&(ipath > tsposition)) {
        Eaux[ipath + 1] += Eaux[ipath];
        thepath.erase(thepath.begin() + ipath);
        Epath.erase(Epath.begin() + ipath);
        Eaux.erase(Eaux.begin() + ipath);
        --lengthpath;
        --ipath;
      }
    }
    lengthpath = thepath.size();
    if (lengthpath >= maxnpoints) {
      thresholdA *= 1.3;
      thresholdB *= 1.2;
      goto Remove_More;
    }
  }
  void Equidistance(const std::vector<double> & vin, std::vector<double> & vout) {
    //function that places all points equidistant
    size_t vsz = vin.size();
    if (vout.size() != vsz) {vout.resize(vsz);}
    double start = vin[0];
    double end = vin[vsz - 1];
    double step = (end - start)/double(vsz - 1);
    vout[0] = start;
    for (size_t idpos = 1; idpos < vsz; ++idpos) {
      start += step;
      vout[idpos] = start;
    }
  }
  void CubicSpline(const std::vector<double> & xvec, const std::vector<double> & xvec2, std::vector<double> & yvec) {
    //function that applies a cubic spline, based on xTB's implementation, which is itself based on C. J. C. Kruger, Constrained Cubic Spline Interpolation for Chemical Engineering Application
    size_t numberpoints = xvec.size();
    matrixE coeff(numberpoints,4);
    getCubicCoefficients(xvec,yvec,coeff);
    for (size_t idpoint = 1; idpoint < numberpoints - 1; ++idpoint) {
      EvalCubicSpline(coeff,xvec,xvec2[idpoint],yvec[idpoint]);
    }
  }
  void getCubicCoefficients(const std::vector<double> & xvec, const std::vector<double> & yvec, matrixE & coeff) {
    //function to calculate the coefficients of a cubic polynomial for applying a spline
    size_t szx = xvec.size();
    std::vector<double> deltax(szx - 1);
    std::vector<double> deltay(szx - 1);
    std::vector<double> function(szx);
    double slope;
    double dx;
    double dx2;
    double dy;
    double xkm2;
    double xkm3;
    double xk2;
    double fxkm;
    double fxk;
    deltax[0] = xvec[1] - xvec[0];
    deltay[0] = yvec[1] - yvec[0];
    for (size_t idpos = 1; idpos < szx - 1; ++idpos) {
      deltax[idpos] = xvec[idpos + 1] - xvec[idpos];
      deltay[idpos] = yvec[idpos + 1] - yvec[idpos];
      slope = deltay[idpos]*deltay[idpos - 1];
      function[idpos] = double(slope > 0.0)*2.0/((deltax[idpos]/deltay[idpos]) + (deltax[idpos - 1]/deltay[idpos - 1]));
    }
    function[0] = 1.5*deltay[0]/(deltax[0]) - 0.5*function[1];
    function[szx - 1] = 1.5*deltay[szx - 2]/(deltax[szx - 2]) - 0.5*function[szx - 2];
    for (size_t idpos = 1; idpos < szx; ++idpos) {
      dx = deltax[idpos - 1];
      dy = deltay[idpos - 1];
      xkm2 = xvec[idpos - 1]*xvec[idpos - 1];
      xkm3 = xkm2*xvec[idpos - 1];
      dx2 = dx*dx;
      xk2 = xvec[idpos]*xvec[idpos];
      fxkm = -2.0*(function[idpos] + 2.0*function[idpos - 1])/dx + 6.0*dy/dx2;
      fxk = 2.0*(2.0*function[idpos] + function[idpos - 1])/dx - 6.0*dy/dx2;
      coeff(idpos + 1,4) = (fxk - fxkm)/(6.0*dx);
      coeff(idpos + 1,3) = 0.5*(xvec[idpos]*fxkm - xvec[idpos - 1]*fxk)/dx;
      coeff(idpos + 1,2) = (dy - coeff(idpos + 1,3)*(xk2 - xkm2) - coeff(idpos + 1,4)*(xk2*xvec[idpos] - xkm3))/dx;
      coeff(idpos + 1,1) = yvec[idpos - 1] - coeff(idpos + 1,2)*xvec[idpos - 1] - coeff(idpos + 1,3)*xkm2 - coeff(idpos + 1,4)*xkm3;
    }
  }
  void EvalCubicSpline(const matrixE & coeff,const std::vector<double> xvec, double xval, double & yval) {
    //function that evaluates a function based on constrained cubic spline
    double xv2;
    size_t totalpoints = xvec.size();
    for (size_t idpoint = 1; idpoint < totalpoints; ++idpoint) {
      if (xval < xvec[idpoint]) {
        if (xval >= xvec[idpoint - 1]) {
          xv2 = xval*xval;
          yval = coeff(idpoint + 1,1) + coeff(idpoint + 1,2)*xval + coeff(idpoint + 1,3)*xv2 + coeff(idpoint + 1,4)*xv2*xval;
        }
      }
    }
  }
  void PathFit(std::vector<matrixE> & thepath) {
    //function fitting a certain path to a cubic function
    size_t Natom = thepath[0].rows();
    size_t npointspath = thepath.size();
    matrixE atomPath(npointspath,3);
    std::vector<double> atomPathSpline(npointspath,0.0);
    std::vector<double> atomPathSplineFix(npointspath,0.0);
    std::vector<double> yk(npointspath);
    double geok[3];
    double geokp1[3];
    double geodiff2;
    //for each atom, build its path
    for (size_t idAtm = 0; idAtm < Natom; ++idAtm) {
      for (size_t ipath = 0; ipath < npointspath; ++ipath) {
        atomPath(ipath + 1,1) = thepath[ipath](idAtm + 1,1);
        atomPath(ipath + 1,2) = thepath[ipath](idAtm + 1,2);
        atomPath(ipath + 1,3) = thepath[ipath](idAtm + 1,3);
      }
      for (size_t ipath = 1; ipath < npointspath; ++ipath) {
        geodiff2 = 0.0;
        for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
          geok[idcoord] = atomPath(ipath,idcoord + 1);
          geokp1[idcoord] = atomPath(ipath + 1,idcoord + 1);
          geodiff2 += (geokp1[idcoord] - geok[idcoord])*(geokp1[idcoord] - geok[idcoord]);
        }
        atomPathSpline[ipath] = atomPathSpline[ipath - 1] + sqrt(geodiff2);
      }
      for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
        for (size_t idpoint = 0; idpoint < npointspath; ++idpoint) {
          yk[idpoint] = atomPath(idpoint + 1,idcoord + 1);
        }
        //place all points equidistantly along the path
        Equidistance(atomPathSpline,atomPathSplineFix);
        CubicSpline(atomPathSpline,atomPathSplineFix,yk);
        for (size_t idpoint = 0; idpoint < npointspath; ++idpoint) {
          atomPath(idpoint + 1,idcoord + 1) = yk[idpoint];
        }
      }
      for (size_t ipath = 0; ipath < npointspath; ++ipath) {
        thepath[ipath](idAtm + 1,1) = atomPath(ipath + 1,1);
        thepath[ipath](idAtm + 1,2) = atomPath(ipath + 1,2);
        thepath[ipath](idAtm + 1,3) = atomPath(ipath + 1,3);
      }
    }
  }
  void Bias2Point(const matrixE & geomTS) {
    MTD.setMetaSet2Structure(geomTS);
    MTD.setMetaAlphas2Alpha(1.3);
    MTD.setMetaFactors2kPull(-5.0*fabs(kpull)*double(Natoms));       //make sure it is negative
  }
  double GuessTSMode(size_t totalpoints, size_t ipoint, const std::vector<matrixE> & path, std::vector<double> & pv1, std::vector<double> & pv2) {
    //function that estimates the harmonic mode of the best guess of the transition state
    size_t ipos = 0;
    double rms1 = 0.0;
    double rms2 = 0.0;
    double f1 = 1.0/3.0;
    double f2 = 2.0/3.0;
    double auxA;
    double auxB;
    double p1p2;
    for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
      for (size_t idcoord = 0; idcoord < 3; ++idcoord, ++ipos) {
        auxA = path[ipoint + 1](idAtm + 1,idcoord + 1) - path[ipoint](idAtm + 1,idcoord + 1);
        if (ipoint + 2 < totalpoints) {
          auxB = path[ipoint + 2](idAtm + 1,idcoord + 1) - path[ipoint + 1](idAtm + 1,idcoord + 1);
          auxA = f2*auxA + f1*auxB;
        }
        pv1[ipos] = auxA;
        rms1 += auxA*auxA;
        auxA = path[ipoint](idAtm + 1,idcoord + 1) - path[ipoint - 1](idAtm + 1,idcoord + 1);
        if (ipoint > 2) {
          auxB = path[ipoint - 1](idAtm + 1,idcoord + 1) - path[ipoint - 2](idAtm + 1,idcoord + 1);
          auxA = f2*auxA + f1*auxB;
        }
        pv2[ipos] = auxA;
        rms2 += auxA*auxA;
      }
    }
    f1 = 1.0/sqrt(rms1);
    f2 = 1.0/sqrt(rms2);
    p1p2 = 0.0;
    for (size_t idpos = 0; idpos < 3*Natoms; ++idpos) {
      pv1[idpos] *= f1;
      pv2[idpos] *= f2;
      p1p2 += pv1[idpos]*pv2[idpos];
    }
    return p1p2;
  }
  bool DetermineRXCoordinates(const std::vector<double> & energies, size_t pathlength, double & maximum, double & rmin, double & pmin) {
    //function that finds the minima and maxima along the reaction path
    //determine the positions to start counting
    bool successful = true;
    int count_start = 0;
    int count_end = pathlength;
    int counter = 0;
    std::vector<int> maxima(1,-1);
    for (size_t idpath = 0; idpath < pathlength - 1; ++idpath) {   //check the first position at which the energy "stabilized" along the path
      if (fabs(energies[idpath] - energies[idpath + 1]) < 0.5) {   //if there are big jumps, ignore these
        count_start = idpath;
        break;
      }
    }
    for (size_t idpath = pathlength; idpath > 1; --idpath) {       //check the last position at which the energy "stabilized" along the path
      if (fabs(energies[idpath] - energies[idpath - 1]) < 0.5) {
        count_end = idpath;
        break;
      }
    }
    //find the maxima
    maximum = -10000000000000000000.0;
    for (size_t idpath = count_start + 1; idpath < count_end; ++idpath) {
      if ((energies[idpath] > energies[idpath - 1])&&(energies[idpath] > energies[idpath + 1])) {
        maxima[counter] = idpath;
        ++counter;
        maxima.push_back(-1);
        if (energies[idpath] > maximum) {
          maximum = energies[idpath];
          tsindex = idpath;
        }
      }
    }
    if (maxima.size() > 1) {
      //find the first minimum (reagent)
      rmin = 10000000000000000000.0;
      for (size_t idpath = count_start + 1; idpath < maxima[0]; ++idpath) {
        if ((energies[idpath] < energies[idpath - 1])&&(energies[idpath] < energies[idpath + 1])&&(energies[idpath] < rmin)) {rmin = energies[idpath];}
      }
      if (energies[count_start] < rmin) {rmin = energies[count_start];}
      //find the last minimum (product)
      pmin = 10000000000000000000.0;
      for (size_t idpath = maxima[counter - 1]; idpath < count_end; ++idpath) {
        if ((energies[idpath] < energies[idpath - 1])&&(energies[idpath] < energies[idpath + 1])&&(energies[idpath] < pmin)) {pmin = energies[idpath];}
      }
      if (energies[count_end] < pmin) {pmin = energies[count_end];}
    }
    else {successful = false;}
    return successful;
  }
  template <class func>
  bool SolverOpt(func & ObjFunc, std::vector<matrixE> & Xiter, matrixE & aux, matrixE & goalgeom, double thresh_energy = 1.0e-7, double thresh_grad = 1.0e-7, int maxiter = 500, bool checkrmsd = false) {
    //geometry optimization algorithm based on the main solvers available, but adapted to this specific case
    std::string type = solver.Type();
    size_t szgrad;
    size_t szrmsdvec = 20;
    int szxiter;
    int degree = solver.Degree();
    int linesearch = solver.LineSearch();
    int doNewton = solver.Newton();
    int failedsteps = 0;
    int totalnrsteps = 0;
    int hessupd = 1;                                                                 //algorithm with normal (1) or inverted (-1) Hessian
    //assume that there will be at most 50 iterations
    Xiter.clear();
    Xiter.resize(50);
    Xiter[0] = ObjFunc.Geometry();
    szxiter = Xiter.size() - 1;
    Xiter[szxiter] = goalgeom;
    ObjFunc.setGeometry(Xiter[0]);
    matrixE sk(1,1);                                                                 //step; in vector form
    matrixE grad_new(1,1);                                                           //new gradient; in vector form
    matrixE grad_old(1,1);                                                           //old gradient; in vector form
    matrixE delta_g(1,1);                                                            //new gradient; in vector form
    matrixE Hessian(1,1);                                                            //hessian matrix; in vector form
    double normA;
    double normB;
    double normG;                                                                    //control variable for checking whether gradient norm is below threshold
    double maxG;                                                                     //control variable for maximum absolute value of gradient
    double normE;                                                                    //control variable for change in energy
    double fold = 0.0;
    double fpredicted;
    double rmsd = 10000.0;
    double rmsd_initial = 100001.0;;
    double thren = thresh_energy;
    if (thren > 1.0) {thren = econvthreshold;}
    double thrgr = thresh_grad;
    if (thrgr > 1.0) {thrgr = gconvthreshold;}
    ObjFunc.Calculate(0,200,true,true,1600.0);
    double fnew = ObjFunc.getEnergy(1);
    double alp = aline;
    //control variables to determine when convergence obtained
    bool converged = false;
    bool convergedG = false;                                      //convergence of gradient
    bool convergedE = false;                                      //convergence of energy
    bool TROK;
    bool skip;
    bool recalculateHess;
    bool RMSDonlyincreases;
    if (degree == 2) {                           //get Hessian or the approximation to it
      ObjFunc.hEnergy(delta_g,hessiantype,1,1,0.03,thrgr,maxiter);
      if ((type == "BFGSi")||(type == "CGBFGS")) {                                                         //this is not really the hessian but its inverse
        Hessian = delta_g.inverse();
        hessupd = -1;
      }
      else if ((type == "Newton-Raphson")||(type == "BFGSd")||(type == "BakerRFO")||(type == "Wales-NR")) {Hessian = delta_g;}
    }
    ObjFunc.gEnergy(grad_new,gradtype,degree == 2,thrgr);
    szgrad = grad_new.rows();
    grad_old.resize(szgrad,1);
    sk.resize(szgrad,1);
    bestpathRMSDtemp.resize(szrmsdvec);                     //control how RMSD is evolving; if it is just increasing, then discard the trajectory
    if ((linesearch == 6)||(linesearch == 7)) {solver.pstore = Hessian;}
    for (int iter = 0; iter < maxiter; ++iter, ++totalnrsteps) {
      std::cout << "iter = " << iter + 1 << "\n";
      //check how RMSD has been evolving?
      if ((rmsd > 10.0*rmsd_initial)&&(checkrmsd)) {               //not worth continuing
        Xiter.resize(iter + 2);
        szxiter = Xiter.size() - 1;
        Xiter[szxiter] = goalgeom;
        break;
      }
      if ((iter == szrmsdvec)&&(checkrmsd)) {
        RMSDonlyincreases = true;
        for (size_t idrmsd = 1; idrmsd < szrmsdvec; ++idrmsd) {
          if (bestpathRMSDtemp[idrmsd] < bestpathRMSDtemp[idrmsd - 1]) {
            RMSDonlyincreases = false;
            break;
          }
        }
        if (RMSDonlyincreases) {                    //if the RMSD has only been increasing, then this is anyway a bad trajectory, so just move on
          Xiter.resize(iter + 2);
          szxiter = Xiter.size() - 1;
          Xiter[szxiter] = goalgeom;
          break;
        }
      }
      if ((iter + 1)%50 == 0) {          //increase storage
        Xiter.resize(iter + 51);
        szxiter = Xiter.size() - 1;
        Xiter[szxiter] = goalgeom;
      }
      TROK = true;
      skip = false;
      recalculateHess = false;
      //get search direction
      fpredicted = alp;
      if (degree == 1) {solver.Step(sk,grad_new,grad_old,fpredicted);}
      else if (degree == 2) {solver.Step(sk,grad_new,Hessian,fpredicted);}     //fpredicted contains the predicted change in energy for trust region methods implemented as such
      normG = 0.0;
      maxG = 0.0;
      normB = 0.0;                                                             //this is to check whether there is trouble with the calculated step
      for (size_t idcoord = 0; idcoord < szgrad; ++idcoord) {
        normA = fabs(grad_new(idcoord + 1,1));
        if (maxG < normA) {maxG = normA;}
        normG += normA*normA;
        normA = fabs(sk(idcoord + 1,1));
        normB += normA*normA;
      }
      normA = sqrt(normB/double(Natoms));
      if ((normA > 50.0)||(normA < 1.0e-6)) {         //if the step is absurd, do something else; this something else should avoid Hessian, since we are here because it is ill-conditioned
        recalculateHess = true;
        sk = grad_new;
      }
      normA = sqrt(normG/double(Natoms));
      normG = normA;
      if ((normG/thresh_grad < 1.0 + thresh_grad)&&(0.5*maxG/thresh_grad < 1.0 + thresh_grad)) {convergedG = true;}                  //converged on gradients?
      else {convergedG = false;}
      if (normG > 10000) {break;}            //something too off
      std::cout << "        |g(k)| = " << normG << "         max|g(k)| = " << maxG << "        congG = " << convergedG << "\n";
      if (linesearch == 1) {alp = solver.Davidon(ObjFunc,sk,grad_new,Hessian,fnew,degree == 2,skip);}                                //Davidon line search
      else if (linesearch == 2) {alp = solver.BarzilaiBorwein(sk,grad_new,grad_old,iter%2,alp);}                                     //Barzilai-Borwein
      else if (linesearch == 3) {alp = solver.Armijo(ObjFunc,sk,grad_new,fold,degree == 2,skip,0.5,0.75,maxiter);}                   //Armijo rule
      else if (linesearch == 4) {alp = solver.MoreThuente(ObjFunc,sk,grad_new,fnew,alp,skip,degree == 2,maxiter);}                   //More'-Thuente
      else if (linesearch == 5) {alp = solver.Fletcher(ObjFunc,sk,grad_new,fnew,alp,degree == 2,skip,maxiter);}                      //Fletcher
      else if (linesearch == 6) {alp = solver.DogLeg(ObjFunc,sk,Hessian,grad_new,fnew,alp,hessupd,TROK,skip);}                       //Dogleg
      else if (linesearch == 7) {alp = solver.Trustregion(ObjFunc,sk,fnew,10000.0,10000.0,fpredicted,alp,thren,TROK,skip);}          //Trust-region things for such methods
      //corrective action needed?
      solver.UpdateAlphaControl(alp,iter);
      if (solver.SteepestDescentStep()) {                                                                                            //if there are too many bad steps, then do steepest descent
        std::cout << "WARNING: SolverPackage.hpp: SolverOpt(): doing steepest descent step\n";
        solver.ResetAlphaControl();                                                                                                  //reset counter
        sk = grad_new;
        alp = aline;
        skip = false;
      }
      else if (!TROK) {
        std::cout << "     -> iteration rejected\n";
        normE *= 0.5;
        ++failedsteps;
        if (failedsteps > 10) {                                              //do steepest descent
          sk = grad_new;
          if (linesearch != 6) {alp = aline;}
          skip = false;
          recalculateHess = true;
          failedsteps = 0;
          TROK = true;
        }
        else {--iter;}
      }
      else {                                                                                                                  //restrict alpha?
        if (alp < solver.AlphaMin()) {alp = solver.AlphaMin();}
        if (alp > solver.AlphaMax()) {alp = solver.AlphaMax();}
        if (alp > aline) {alp = aline;}
      }
      //update x and get norm of step
      if (TROK) {
        failedsteps = 0;
        //update function, gradient
        fold = fnew;
        grad_old = grad_new;
        if (!skip) {
          delta_g = sk*alp;
          ObjFunc.updGeometry(delta_g,degree == 2);
          ObjFunc.Calculate(0,200,true,true,1600.0);
        }
        Xiter[iter + 1] = ObjFunc.Geometry();
        fnew = ObjFunc.getEnergy(1);
        normE = fabs(fnew - fold);
        std::cout << "        |Delta E| = " << normE << "\n";
        ObjFunc.gEnergy(grad_new,gradtype,degree == 2,thrgr);
        if (degree == 2) {
          if (!recalculateHess) {
            recalculateHess = (int(fpredicted) == -1) + (type == "Newton-Raphson");
            if (((type == "BakerRFO")||(type == "Wales-NR"))&&(((iter + 1)%5 == 0)&&(doNewton))) {
              recalculateHess = true;
              alp = aline;                            //reset trust region
            }
          }
          //update Hessian or recalculate it?
          if (recalculateHess) {ObjFunc.hEnergy(Hessian,hessiantype,1,1,0.03,thrgr,maxiter);}
          else {
            delta_g = grad_new - grad_old;
            solver.HessianUpdate(Hessian,sk,delta_g,alp);
          }
        }
        alp = aline;
        rmsd = MTD.RMSD(Xiter[iter + 1],goalgeom,aux);
      }
      if (normE/thresh_energy < 1.0 + thresh_energy) {convergedE = true;}                                             //converged?
      else {convergedE = false;}
      if ((iter < szrmsdvec)&&(checkrmsd)) {
        if (iter == 0) {rmsd_initial = rmsd;}
        bestpathRMSDtemp[iter] = rmsd;
      }
      std::cout << "current rmsd = " << rmsd << "     " << rmsdthreshold_product << "\n";
      if ((convergedG)&&(convergedE)) {
        if (checkrmsd) {
          if (rmsd < rmsdthreshold_product) {converged = true;}
        }
        else {converged = true;}
      }
      if ((converged)&&(iter > 1)) {
        Xiter.resize(iter + 2);
        szxiter = Xiter.size() - 1;
        Xiter[szxiter] = goalgeom;
        std::cout << "Convergence in " << iter + 1 + (!TROK) << " iterations.\n";
        break;
      }
      if (totalnrsteps > maxiter) {
        Xiter.resize(maxiter + 1);
        szxiter = Xiter.size() - 1;
        Xiter[szxiter] = goalgeom;
      }
    }
    if (!converged) {
      Xiter.resize(maxiter + 2);
      szxiter = Xiter.size() - 1;
      Xiter[szxiter] = goalgeom;
    }
    return converged;
  }
  template<class T>
  void FindPath(T & ElecStruct, int nalpha = 4, int nkappa = 6, bool simplify = false) {
    //main function to find a path
    //nalpha is the number of runs with different values for alpha exponent in metadynamics
    //nkappa is the number of runs with different values for the push/pull forces
    //variables
    gradtype = ElecStruct.AvailableGradients();
    int totalnmol;
    size_t maxnpaths = nalpha*nkappa;
    size_t pathlength;
    size_t irun = 0;
    double rmsd;
    double Afactor;
    double Kfactor;
    double auxC;
    double auxD;
    double NMatms;
    double barrierdirect;
    double barrierreverse;
    double deltaE;
    double bestpathdirectbarrier = 1.0e100;
    double bestpathreversebarrier = 1.0e100;
    double bestpathErx;
    double bestpathrmsd;
    double p1dotp2;
    double gradpath;
    bool directbarrierOK;
    bool reversebarrierOK;
    bool RPdegenerate = false;
    bool pathfound = false;
    bool converged;
    bool energypathgood;
    matrixE geomR = reagent.Geometry();
    matrixE geomP = product.Geometry();
    matrixE gradient = geomR;
    std::vector<double> energies(50);
    std::vector<double> pvec1(3*Natoms);
    std::vector<double> pvec2(3*Natoms);
    bestpathenergies.resize(50);
    std::vector<matrixE> pathRP;
    CompositeMTDQM mixmethod(MTD,ElecStruct,gradient);          //this contains pointers to the two methods, so the object should be always up-to-date
    ElecStruct.setElectronTemp(Telectron);
    ElecStruct.setMolecule(reagent);
    ElecStruct.setRestart(0);                       //use restart, otherwise this will be a killer
    ElecStruct.Calculate(0,200,true,true,1600.0);
    MTD.PreparePathFinder(geomR,geomP,kpush,kpull,alpha,ppull);
    rmsd = MTD.RMSD(geomR,geomP,gradient);                //gradients not really needed here, so geometry is dummy
    std::cout << std::setprecision(10);
    std::cout << "RMSD(reagent,product) = " << rmsd << "\n";
    //if user did not supply anything, then consider all atoms meta
    if (MTD.TotalMetaAtoms() == 0) {MTD.addMetaRange(1,Natoms);}
    NMatms = double(MTD.TotalMetaAtoms());
    totalnmol = MoleculeCounter(product);          //count the number of sub-molecules in the product state
    rmsdthreshold_product *= double(totalnmol);    //in these cases the RMSD from product might be large and the algorithm does not terminate
    //check for degenerate structures
    Kfactor = h*1.0e18/(8.0*pi*pi*c);
    auxC = 0.0;
    auxD = 0.0;
    inertiaR = reagent.InertiaEigenvalues(false);
    inertiaP = product.InertiaEigenvalues(false);
    for (size_t idx = 0; idx < 3; ++idx) {
      Afactor = inertiaR[idx];
      inertiaR[idx] = Kfactor/Afactor;
      Afactor = inertiaP[idx];
      inertiaP[idx] = Kfactor/Afactor;
      auxC += inertiaR[idx] - inertiaP[idx];
      auxD += inertiaR[idx] + inertiaP[idx];
    }
    Kfactor = 2.0*fabs(auxC)/auxD;
    if (Kfactor < 0.01) {
      RPdegenerate = true;         // <-> degen
      std::cout << "WARNING: MolecularDynamics.hpp: PathFinder: FindPath(): reagent and product are too similar\n";
    }
    Afactor = 0.0;
    tsindex = 0;
    //find paths
    for (size_t ialpha = 0; ialpha < nalpha; ++ialpha) {
      std::cout << "alpha = " << ialpha << "\n";
      Kfactor = 1.0;
      for (size_t ikappa = 0; ikappa < nkappa; ++ikappa, ++irun) {
        std::cout << "  kappa = " << ikappa <<  "\n";
        MTD.RescalePathFinder(Kfactor*NMatms,Afactor,true,0.2);
        //now get the path
        geomR = reagent.Geometry();
        geomP = product.Geometry();
        if (RPdegenerate) {MTD.Randomize(geomR,0.05);}                        //randomize structure! avoid symmetry trapping
        //preparing structures
        mixmethod.setGeometry(geomR);
        MTD.PathFinderResetStructures(geomR,geomP);
        //run optimization path from current reagent
        converged = SolverOpt(mixmethod,pathRP,gradient,geomP,10000.0,10000.0,200,true);  //here we do not care about energy nor gradients, only about RMSD
        if (!converged) {
          if (pathfound) {break;}
          else {continue;}                                                                //if the product was not found, then just move on with life
        }
        pathlength = pathRP.size();
        //check whether what happened was in the right direction
        geomR = pathRP[pathlength - 1];
        rmsd = MTD.RMSD(geomR,geomP,gradient);
        if (energies.size() < pathlength) {energies.resize(pathlength);}
        barrierdirect = 0.0;
        barrierreverse = 0.0;
        std::cout << "    energies along current path\n";
        for (size_t idpath = 0; idpath < pathlength; ++idpath) {
          ElecStruct.setGeometry(pathRP[idpath]);
          ElecStruct.Calculate(0,200,true,true,1600.0);
          energies[idpath] = ElecStruct.getEnergy(1);
          std::cout << "      " << idpath << "  " << energies[idpath] << "\n";
        }
        energypathgood = DetermineRXCoordinates(energies,pathlength,auxD,auxC,deltaE);
        if (!energypathgood) {
          if (pathfound) {break;}
          else {continue;}
        }
        barrierdirect = auxD - auxC;
        barrierreverse = auxD - deltaE;
        deltaE -= auxC;
        PrintReactionProfile(deltaE,barrierdirect,barrierreverse,true,"current reaction energies",7);
        //check the barriers
        directbarrierOK = ((barrierdirect > 0.0)&&(barrierdirect < 100.0)&&(barrierdirect < bestpathdirectbarrier));          //make sure things are within reasonable bounds
        reversebarrierOK = ((barrierreverse > 0.0)&&(barrierreverse < 100.0)&&(barrierreverse < bestpathreversebarrier));     //make sure things are within reasonable bounds
        //is this the best path so far?
        if ((directbarrierOK)&&(reversebarrierOK)&&(rmsd < rmsdthreshold_product)) {
          bestpathdirectbarrier = barrierdirect;
          bestpathreversebarrier = barrierreverse;
          bestpathenergies = energies;
          bestpathErx = energies[pathlength - 1] - energies[0];
          bestpathrmsd = rmsd;
          bestpathRP = pathRP;
          std::cout << "update of current best path with barrier: " << bestpathdirectbarrier << " Hartree\n";
        }
        //rescale
        Kfactor *= deltak;
        if (pathfound) {break;}
        //again, check the current barriers
        directbarrierOK = ((barrierdirect > 0.0)&&(barrierdirect < 100.0));
        reversebarrierOK = ((barrierreverse > 0.0)&&(barrierreverse < 100.0));
        if ((rmsd < rmsdthreshold_product)&&(directbarrierOK)&&(reversebarrierOK)) {              //if product was found try a softer push/pull to see if something is improved
          Kfactor /= deltak*1.2;
          pathfound = true;
        }
      }
      if (pathfound) {break;}
      Afactor += deltaalpha;
    }
    if (!pathfound) {std::cout << "WARNING: MolecularDynamics.hpp: PathFinder: FindPath(): no product found; try smaller alphas\n";}
    else {
      auxC = MTD.Alpha()/1.2;
      ppull *= auxC*auxC;                 //adjust ppull to the value of alpha used
      //refine the reaction path 
      pathlength = bestpathRP.size();
      energies.resize(pathlength);
      bestpathenergies.resize(pathlength);
      if (pathlength > MaxNpoints) {
        PathCleaner(bestpathRP,bestpathenergies,energies,MaxNpoints);
        pathlength = bestpathRP.size();
      }
      bestpathgradient.resize(pathlength);
      bestpathRMSD.resize(pathlength);
      bestpathgradient[0] = 0.0;
      bestpathgradient[pathlength - 1] = 0.0;
      geomR = reagent.Geometry();
      barrierdirect = 0.0;
      pathRP = bestpathRP;
      Afactor = bestpathenergies[0];
      startindex = 0;
      if (simplify) {
        //getting the RMSD from reagent to define reaction axis
        for (size_t idpstep = 0; idpstep < pathlength; ++idpstep) {
          mixmethod.setGeometry(bestpathRP[idpstep]);
          //get RMSD from reagent
          bestpathRMSD[idpstep] = MTD.RMSD(bestpathRP[idpstep],geomR,gradient);
          deltaE = bestpathenergies[idpstep] - bestpathenergies[0];
          if (deltaE > barrierdirect) {
            barrierdirect = deltaE;
            tsindex = idpstep;
          }
          if (Afactor > bestpathenergies[idpstep]) {
            startindex = idpstep;
            Afactor = bestpathenergies[idpstep];
          }
        }
      }
      else {
        //refining the path found
        //std::cout << "refining the best path\n";
        std::cout << "refining the best path" << std::endl;
        for (size_t ipopt = 0; ipopt < nrefinements; ++ipopt) {
          //PathFit(pathRP);                   //apply cubic spline to make the path nicer
          barrierdirect = 0.0;
          //std::cout << "    energies along best path\n";
          std::cout << "    energies along best path" << std::endl;
          for (size_t idpstep = 0; idpstep < pathlength; ++idpstep) {
            mixmethod.setGeometry(bestpathRP[idpstep]);
            //mixmethod.setGeometry(pathRP[idpstep]);
            Bias2Point(bestpathRP[idpstep]);
            converged = SolverOpt(mixmethod,pathRP,gradient,bestpathRP[idpstep],econvthreshold,gconvthreshold,200,false);
            bestpathRP[idpstep] = pathRP[pathRP.size() - 2];
            ElecStruct.setGeometry(bestpathRP[idpstep]);
            ElecStruct.Calculate(0);
            bestpathenergies[idpstep] = ElecStruct.getEnergy(1);
            std::cout << idpstep << "  " << bestpathenergies[idpstep] << std::endl;
            //get RMSD from reagent
            bestpathRMSD[idpstep] = MTD.RMSD(bestpathRP[idpstep],geomR,gradient);
            deltaE = bestpathenergies[idpstep] - bestpathenergies[0];
            if (deltaE > barrierdirect) {
              barrierdirect = deltaE;
              tsindex = idpstep;
            }
            if (Afactor > bestpathenergies[idpstep]) {
              startindex = idpstep;
              Afactor = bestpathenergies[idpstep];
            }
          }
        }
      }
      energypathgood = DetermineRXCoordinates(bestpathenergies,pathlength,auxD,auxC,deltaE);
      PrintReactionProfile(auxC,auxD,deltaE,false,"optimized reaction energies",7);
      //trying to guess the TS mode, which estimates the gradient along the PES
      for (size_t ipath = 1; ipath < pathlength - 1; ++ipath) {
        p1dotp2 = GuessTSMode(pathlength,ipath,bestpathRP,pvec1,pvec2);
        geomR = bestpathRP[ipath];
        geomP = bestpathRP[ipath - 1];
        rmsd = MTD.RMSD(geomR,geomP,gradient);
        gradpath = (bestpathenergies[ipath] - bestpathenergies[ipath - 1])/(rmsd*double(Natoms));   //gradient on point
        bestpathgradient[ipath] = gradpath;
      }
    }
  }
};
class MolDynamicsHS {
  //class for hard-sphere molecular dynamics
  int Nparticle;                                   //total number of particles in system
  int ndim;                                        //dimension of box
  int nevents;                                     //total of events to record
  double mass;                                     //particle mass
  double sigma;                                    //particle radius
  double time_sim;                                 //time ran in simulation
  double time_max;                                 //max time to run simulation
  matrixE CollisionW;                              //matrix with collision times with walls; index is particle; 
                                                   //particle i [time; wall coordinate]
  matrixE CollisionP;                              //matrix with pair collision times            
  std::vector<std::vector<double> > particles;     //particle positions
  std::vector<std::vector<double> > velocities;    //their velocities
  std::vector<std::vector<double> > accellerations;//and their accelerations
  std::vector<double> boxdim;                      //maximum delimiter of box
  double maxtime;                                  //definition of infinity time for minimization purposes
  //folder information to save states
  std::string folderevents;
  std::string foldertimes;
  std::string fileevents;
  std::string filetimes;
  bool PBC;                     //use of periodic boundary conditions
public:
  MolDynamicsHS(int _Nparticle, int _ndim, int _nevents, double _time_max, double _mass, double _sigma, std::vector<double> _boxdim, bool _pbc = true) {
    if (_boxdim.size() != _ndim) {throw("MolDynamicsHS - MolDynamicsHS: defining box with different number of dimensions than space.");}
    Nparticle = _Nparticle;
    ndim = _ndim;
    nevents = _nevents;
    mass = _mass;
    sigma = _sigma;
    boxdim = _boxdim;
    time_max = _time_max;
    maxtime = 10000000.0;
    PBC = _pbc;
    //plant random seed
    srand(time(NULL));
    //folders and files to save data
    folderevents = "events";
    foldertimes = "times";
    fileevents = "events.out";
    filetimes = "times.out";
  }
  ~MolDynamicsHS() {}
  //getters ------------------------------------------------------------------------------------------------------------------------------
  int Nparticles() {return Nparticle;}
  double RadiusParticle() {return sigma;}
  double MaxTime() {return time_max;}
  int Ndimensions() {return ndim;}
  int Nevents() {return nevents;}
  bool PeriodicConditions() {return PBC;}
  std::vector<double> BoxDimensions() {return boxdim;}
  std::vector<std::vector<double> > PosParticles() {return particles;}
  std::vector<std::vector<double> > VelParticles() {return velocities;}
  //setters ------------------------------------------------------------------------------------------------------------------------------
  void setNparticles(int _npart) {Nparticle = _npart;}
  void setParticlePositions(const std::vector<std::vector<double> > & newpos) {particles = newpos;}
  void setRadius(double _sigma) {sigma = _sigma;}
  void setMaxTime(double _time_max) {time_max = _time_max;}
  void setDimensions(int _dim) {ndim = _dim;}
  void setNevents(int _nevents) {nevents = _nevents;}
  void setBoxDimensions(std::vector<double> dimensions) {
    boxdim = dimensions;
    if (boxdim.size() != ndim) {throw("MonteCarlo - MonteCarlo: defining box with more dimensions than space.");}
  }
  void setFileEvents(std::string _file) {fileevents = _file;}
  void setFileTimes(std::string _file) {filetimes = _file;}
  void setFolderEvents(std::string _folder) {folderevents = _folder;}
  void setFolderTimes(std::string _folder) {foldertimes = _folder;}
  void setPBC(bool _pbc) {PBC = _pbc;}
  //other functions ----------------------------------------------------------------------------------------------------------------------
  void PlaceParticles(int maxiter = 200) {
    //interface to particle distribution function
    std::vector<std::vector<double> > * pParticles = & particles;
    std::vector<double> * pboxdim = & boxdim;
    DistributeParticles(pParticles,pboxdim,Nparticle,ndim,sigma,maxiter);
  }
  void setVelocities(std::vector<std::vector<double> > * _velocities) {
    if (_velocities->size() != Nparticle) {throw("MolDynamicsHS - setVelocities: vector with velocities does not have correct number of particles.");}
    for (size_t idpart = 0; idpart < Nparticle; ++idpart) {
      for (size_t idcoord = 0; idcoord < ndim; ++idcoord) {
        velocities[idpart][idcoord] = _velocities->at(idpart)[idcoord];
      }
    }
  }
  void CanonicalVelocities(double Temp, double deltaT = 0.01) {
    //funtion that generates particle velocities according to system's temperature, using Gaussian distribution
    //temperatures in Kelvin, naturally
    //generate average velocity per coordinate
    double velocity2_avg = 2.0*KB*Temp/(mass*double(ndim));
    double velocity2_max = 2.0*KB*(Temp + deltaT)/(mass*double(ndim));
    double deltaV;
    std::vector<double> velocity(ndim,0.0);
    for (size_t idpart = 0; idpart < Nparticle; ++idpart) {
      for (size_t idcoord = 0; idcoord < ndim; ++idcoord) {
        deltaV = sqrt(velocity2_max) - sqrt(velocity2_avg);
        velocity[idcoord] = pow(-1.0,double(iRandom01()))*fRandom(sqrt(velocity2_avg) - deltaV,sqrt(velocity2_avg) + deltaV);
      }
      velocities.push_back(velocity);
    }
  }
  //functions to determine next event
  double CollisionWall_time(size_t iparticle, size_t icoord) {
    //function to calculate time necessary to collide against wall
    double tcollision = maxtime;                                           //case velocity = 0
    if (!PBC) {
      if (velocities[iparticle][icoord] > 0.0) {
        //particle goes to upper limit of wall
        tcollision = (boxdim[icoord] - sigma - particles[iparticle][icoord])/velocities[iparticle][icoord];
      }
      else if (velocities[iparticle][icoord] < 0.0) {
        //particle goes to lower limit of wall
        tcollision = (sigma - particles[iparticle][icoord])/velocities[iparticle][icoord];
      }
    }
    else {
      if (velocities[iparticle][icoord] > 0.0) {
        //particle goes to upper limit of wall
        tcollision = (boxdim[icoord] - particles[iparticle][icoord])/velocities[iparticle][icoord];
      }
      else if (velocities[iparticle][icoord] < 0.0) {
        //particle goes to lower limit of wall
        tcollision = -particles[iparticle][icoord]/velocities[iparticle][icoord];
      }
    }
    return tcollision;
  }
  double PairCollision_time(int ipart1, int ipart2, double threshold = 1.0e-10) {
    //function calculating collision time between two particles
    double t_paircollision = maxtime; 
    if (accellerations.size() == 0) {
      //solve quadratic polynomial
      std::vector<double> delta_r(ndim,0.0);
      std::vector<double> delta_s(ndim,0.0);
      std::vector<double> delta_v(ndim,0.0);
      std::vector<double> delta_r_norm(pow(2.0,ndim),0.0);
      std::vector<double> crossprod(pow(2.0,ndim),0.0);
      std::vector<double> sqrtterm(pow(2.0,ndim),0.0);
      std::vector<double> collision_times(pow(2.0,ndim),maxtime);
      double delta_v_norm = 0.0;
      std::vector<bool> collisiontype(ndim,0.0);      //determine whether collision takes place by going inside or outside of box
                                                      //0 -> internal; 1 -> external
      int upper;                                      //which is "upper" particle, i.e., the particle with larger position vector component in the box?
      int lower;                                      //which is "upper" particle, i.e., the particle with larger position vector component in the box?
      for (size_t idx = 0; idx < ndim; ++idx) {
        std::cout << particles[ipart1][idx] << " " << particles[ipart2][idx] << " " << particles[ipart1][idx] - particles[ipart2][idx] << " " << boxdim[idx] - particles[ipart1][idx] + particles[ipart2][idx] << std::endl;
        //determine which is the "upper" particle
        if (particles[ipart1][idx] > particles[ipart2][idx]) {       //box defined between 0 and L>0, therefore coordinates are supposed to be always positive
          upper = ipart1;
          lower = ipart2;
        }
        else {
          upper = ipart2;
          lower = ipart1;
        }
        //determine whether collision is internal or external
        collisiontype[idx] = 0;         //default is internal collision
        //zero or threshold?
        if (PBC) {                      //external collisions only allowed with periodic boundary conditions
          if ((velocities[upper][idx] > 0.0)&&(velocities[lower][idx] < 0.0)) {collisiontype[idx] = 1;}
          else if ((velocities[upper][idx] > 0.0)&&(velocities[lower][idx] > 0.0)&&(velocities[upper][idx] > velocities[lower][idx])) {collisiontype[idx] = 1;}
          else if ((velocities[upper][idx] < 0.0)&&(velocities[lower][idx] < 0.0)&&(velocities[upper][idx] > velocities[lower][idx])) {collisiontype[idx] = 1;}
        }
        //get distance contribution
        delta_r[idx] = particles[upper][idx] - particles[lower][idx];
        delta_s[idx] = delta_r[idx] - boxdim[idx];
        std::cout << "trials  " << particles[upper][idx] - particles[lower][idx] << "   " << particles[upper][idx] - particles[lower][idx] - boxdim[idx] << std::endl;
        delta_v[idx] = velocities[upper][idx] - velocities[lower][idx];
        delta_v_norm += delta_v[idx]*delta_v[idx];
        delta_r_norm[0] += delta_r[idx]*delta_r[idx];
        crossprod[0] += delta_r[idx]*delta_v[idx];
      }
      std::cout << "Dr     Dv" << std::endl;
      for (size_t idx = 0; idx < delta_r.size(); ++idx) {
        std::cout << delta_r[idx] << "  " << delta_v[idx] << std::endl;
      }
      sqrtterm[0] = crossprod[0]*crossprod[0] - delta_v_norm*(delta_r_norm[0] - 4.0*sigma*sigma);
      std::cout << "<dr|dv>      |dv|2      |dr|2      sqrt" << std::endl;
      std::cout << crossprod[0] << " " << delta_v_norm << " " << delta_r_norm[0] << " " << sqrtterm[0] << std::endl;
      if ((sqrtterm[0] > 0.0)&&(crossprod[0] < 0.0)) {collision_times[0] = -(crossprod[0] + sqrt(sqrtterm[0]))/delta_v_norm;}
      t_paircollision = collision_times[0];
      if (PBC) {
        //calculate others
        delta_r_norm[1] = delta_r[0]*delta_r[0] + delta_s[1]*delta_s[1];
        delta_r_norm[2] = delta_s[0]*delta_s[0] + delta_r[1]*delta_r[1];
        delta_r_norm[3] = delta_s[0]*delta_s[0] + delta_s[1]*delta_s[1];
        crossprod[1] = delta_r[0]*delta_v[0] + delta_s[1]*delta_v[1];
        crossprod[2] = delta_s[0]*delta_v[0] + delta_r[1]*delta_v[1];
        crossprod[3] = delta_s[0]*delta_v[0] + delta_s[1]*delta_v[1];
        if (delta_r_norm.size() > 4) {               //3 case
          delta_r_norm[4] = delta_r_norm[0] - delta_r[2]*delta_r[2] + delta_s[2]*delta_s[2];
          delta_r_norm[5] = delta_r_norm[1] + delta_s[2]*delta_s[2];
          delta_r_norm[6] = delta_r_norm[2] + delta_s[2]*delta_s[2];
          delta_r_norm[7] = delta_r_norm[3] + delta_s[2]*delta_s[2];
          delta_r_norm[1] += delta_r[2]*delta_r[2];
          delta_r_norm[2] += delta_r[2]*delta_r[2];
          delta_r_norm[3] += delta_r[2]*delta_r[2];
          crossprod[4] = crossprod[0] - delta_r[2]*delta_v[2] + delta_s[2]*delta_v[2];
          crossprod[5] = crossprod[1] + delta_s[2]*delta_v[2];
          crossprod[6] = crossprod[2] + delta_s[2]*delta_v[2];
          crossprod[7] = crossprod[3] + delta_s[2]*delta_v[2];
          crossprod[1] += delta_r[2]*delta_v[2];
          crossprod[2] += delta_r[2]*delta_v[2];
          crossprod[3] += delta_r[2]*delta_v[2];
        }
        for (size_t idx = 1; idx < sqrtterm.size(); ++idx) {
          sqrtterm[idx] = crossprod[idx]*crossprod[idx] - delta_v_norm*(delta_r_norm[idx] - 4.0*sigma*sigma);
          std::cout << crossprod[idx] << " " << delta_v_norm << " " << delta_r_norm[idx] << " " << sqrtterm[idx] << std::endl;
          if ((sqrtterm[idx] > 0.0)&&(crossprod[idx] < 0.0)) {collision_times[idx] = -(crossprod[idx] + sqrt(sqrtterm[idx]))/delta_v_norm;}
          if ((collision_times[idx] < t_paircollision)&&(collision_times[idx] > threshold)&&(fabs(delta_r_norm[idx] - 4.0*sigma*sigma) > threshold)) {t_paircollision = collision_times[idx];}
        }
        std::cout << "collision times" << std::endl;
        for (size_t idx = 0; idx < collision_times.size(); ++idx) {
          std::cout << collision_times[idx] << std::endl;
        }
      }
    }
    else {}
    return t_paircollision;
  }
  //functions to solve next event
  void CollisionWall(size_t iparticle, size_t icoord, double threshold = 1.0e-7) {
    //function that performs collision with wall for particle iparticle along coordinate icoord
    if (!PBC) {velocities[iparticle][icoord] *= -1.0;}
    else {
      double factor = 1.0;
      if (fabs(particles[iparticle][icoord] - boxdim[icoord]) < threshold) {factor = -1.0;}
      std::cout << "old " << particles[iparticle][icoord] << "   new " << factor*boxdim[icoord] << std::endl;
      particles[iparticle][icoord] += factor*boxdim[icoord];}
  }
  void PairCollision(size_t ipart1, size_t ipart2, double threshold = 1.0e-7) {
    //function that performs collision between particles ipart1 and ipart2
    std::vector<double> delta_r(ndim,0.0);
    std::vector<double> delta_v(ndim,0.0);
    double delta_r_norm = 0.0;
    double crossterm = 0.0;
    double aux1;
    double aux2;
    //bool collisiontype;                           //determine whether collision takes place by going inside or outside of box
                                                    //0 -> internal; 1 -> external
    std::vector<int> upper(ndim,0);                 //which is "upper" particle, i.e., the particle with larger position vector component in the box?
    std::vector<int> lower(ndim,0);                 //which is "upper" particle, i.e., the particle with larger position vector component in the box?
    //get deltas and norms
    for (size_t idx = 0; idx < ndim; ++idx) {
      std::cout << particles[ipart1][idx] << " " << particles[ipart2][idx] << std::endl;
      //determine which is the "upper" particle
      if (particles[ipart1][idx] > particles[ipart2][idx]) {       //box defined between 0 and L>0, therefore coordinates are supposed to be always positive
        upper[idx] = ipart1;
        lower[idx] = ipart2;
      }
      else {
        upper[idx] = ipart2;
        lower[idx] = ipart1;
      }
      //get distance contribution
      aux1 = particles[upper[idx]][idx] - particles[lower[idx]][idx];
      aux2 = aux1 - boxdim[idx];
      if (fabs(aux1) < fabs(aux2)) {delta_r[idx] = aux1;}
      else {delta_r[idx] = aux2;}
      std::cout << "trials  " << particles[upper[idx]][idx] - particles[lower[idx]][idx] << "   " << particles[upper[idx]][idx] - particles[lower[idx]][idx] - boxdim[idx] << std::endl;
      delta_v[idx] = velocities[upper[idx]][idx] - velocities[lower[idx]][idx];
      delta_r_norm += delta_r[idx]*delta_r[idx];
      crossterm += delta_r[idx]*delta_v[idx];
    }
    std::cout << "distance = " << delta_r_norm << std::endl;
    for (size_t idx = 0; idx < ndim; ++idx) {
      std::cout << " v(" << upper[idx]+1 << "," << idx << ") = " << velocities[upper[idx]][idx] << "-0.25*" << delta_r[idx] << "*" << crossterm << "/" << sigma*sigma << std::endl;
      std::cout << " v(" << lower[idx]+1 << "," << idx << ") = " << velocities[lower[idx]][idx] << "+0.25*" << delta_r[idx] << "*" << crossterm << "/" << sigma*sigma << std::endl;
      velocities[upper[idx]][idx] -= 0.25*delta_r[idx]*crossterm/(sigma*sigma);
      velocities[lower[idx]][idx] += 0.25*delta_r[idx]*crossterm/(sigma*sigma);
    }
  }
  void PairCollision_new(size_t ipart1, size_t ipart2, double threshold = 1.0e-10) {
    //function that performs collision between particles ipart1 and ipart2
    std::vector<double> delta_r(ndim,0.0);
    std::vector<double> delta_v(ndim,0.0);
    double delta_r_norm = 0.0;
    bool collisiontype;                             //determine whether collision takes place by going inside or outside of box
                                                    //0 -> internal; 1 -> external
    int upper;                                      //which is "upper" particle, i.e., the particle with larger position vector component in the box?
    int lower;                                      //which is "upper" particle, i.e., the particle with larger position vector component in the box?
    //get deltas and norms
    for (size_t idx = 0; idx < ndim; ++idx) {
      //determine which is the "upper" and the "lower" particles
      if (particles[ipart1][idx] > particles[ipart2][idx]) {       //box defined between 0 and L>0, therefore coordinates are supposed to be always positive
        upper = ipart1;
        lower = ipart2;
      }
      else {
        upper = ipart2;
        lower = ipart1;
      }
      //determine whether collision is internal or external
      collisiontype = 0;              //default is internal collision
      //zero or threshold?
      if (PBC) {                      //external collisions only allowed with periodic boundary conditions
        if ((velocities[upper][idx] > 0.0)&&(velocities[lower][idx] < 0.0)) {collisiontype = 1;}
        else if ((velocities[upper][idx] > 0.0)&&(velocities[lower][idx] > 0.0)&&(velocities[upper][idx] > velocities[lower][idx])) {collisiontype = 1;}
        else if ((velocities[upper][idx] < 0.0)&&(velocities[lower][idx] < 0.0)&&(velocities[upper][idx] > velocities[lower][idx])) {collisiontype = 1;}
      }
      delta_r[idx] = particles[upper][idx] - particles[lower][idx];
      if (collisiontype == 1) {delta_r[idx] -= boxdim[idx];}
      delta_v[idx] = velocities[upper][idx] - velocities[lower][idx];
      delta_r_norm += delta_r[idx]*delta_r[idx];
    }
    std::cout << "distance = " << delta_r_norm << std::endl;
    //get normalized vector
    double sqrt_delta_r_norm = sqrt(delta_r_norm);
    if (sqrt_delta_r_norm < threshold) {sqrt_delta_r_norm = 1.0;}
    delta_r_norm = 0.0;
    for (size_t idx = 0; idx < ndim; ++idx) {
      delta_r[idx] /= sqrt_delta_r_norm;
      delta_r_norm += delta_v[idx]*delta_r[idx];
    }
    for (size_t idx = 0; idx < ndim; ++idx) {
      velocities[upper][idx] -= delta_r[idx]*delta_r_norm;
      velocities[lower][idx] += delta_r[idx]*delta_r_norm;
    }
  }
  void UpdatePositions(double time) {
    //function that updates the positions of all particles except particles ipart1 and ipart2
    for (size_t idpart = 0; idpart < Nparticle; ++idpart) {
      for (size_t idcoord = 0; idcoord < ndim; ++idcoord) {
        particles[idpart][idcoord] += time*velocities[idpart][idcoord];
      }
    }
  }
  //run simulations
  void RunMDEvents(bool continue_sim = false, bool snap = true, bool write2file = true) {
    //run event driven molecular dynamics simulation
    if (particles.size() != Nparticle) {
      std::cout << "Particle array inconsistent with simulation constructor. Recalculating particle array." << std::endl;
      PlaceParticles();
    }
    if (velocities.size() != Nparticle) {
      std::cout << "Velocity array inconsistent with simulation constructor. Recalculating velocity array with canonical distribution at 298.15 K." << std::endl;
      CanonicalVelocities(298.15);
    }
    std::vector<std::vector<double> > * pparticles = & particles;
    std::vector<std::vector<double> > * pvelocities = & velocities;
    std::vector<int> mintpaircoll(2,0);
    double timecollision = 0.0;
    double mintimecoll = maxtime;
    double timeevent;
    double coordinate;
    double dmintpaircoll;
    int pmintwall = 0;
    if (!continue_sim) {
      time_sim = 0.0;
      //determine collision times
      CollisionW.resize(Nparticle,2);
      for (size_t idpart = 0; idpart < Nparticle; ++idpart) {
        mintimecoll = maxtime;
        for (size_t idcoord = 0; idcoord < ndim; ++idcoord) {
          timecollision = CollisionWall_time(idpart,idcoord);
          if (mintimecoll > timecollision) {
            mintimecoll = timecollision;
            coordinate = double(idcoord);
          }
        }
        CollisionW(idpart + 1, 1) = mintimecoll;
        CollisionW(idpart + 1, 2) = coordinate;
        if (mintimecoll < CollisionW(pmintwall + 1,1)) {pmintwall = idpart;}
      }
      CollisionP.resize(Nparticle,Nparticle);
      dmintpaircoll = maxtime;
      for (size_t idpart1 = 0; idpart1 < Nparticle; ++idpart1){
        for (size_t idpart2 = idpart1 + 1; idpart2 < Nparticle; ++idpart2){
          CollisionP(idpart1 + 1,idpart2 + 1) = PairCollision_time(idpart1,idpart2);
          CollisionP(idpart2 + 1,idpart1 + 1) = CollisionP(idpart1 + 1,idpart2 + 1);
          if (CollisionP(idpart1 + 1,idpart2 + 1) < dmintpaircoll) {
            dmintpaircoll = CollisionP(idpart1 + 1,idpart2 + 1);
            mintpaircoll[0] = idpart1;
            mintpaircoll[1] = idpart2;
          }
        }
      }
    }
    //snap initial state?
    if ((snap)||(write2file)) {
      if (!continue_sim) {MD2file(true,Nparticle,ndim,sigma,true,time_sim,pparticles,pvelocities,boxdim,fileevents,folderevents);}
      MD2file(false,Nparticle,ndim,sigma,true,time_sim,pparticles,pvelocities,boxdim,fileevents,folderevents);
    }
    std::cout << "wall collisions" << std::endl;
    CollisionW.Print();
    std::cout << "minimum wall collision with particle " << pmintwall + 1 << std::endl;
    std::cout << "pair collisions" << std::endl;
    CollisionP.Print();
    std::cout << "pair with lowest collision time: " << mintpaircoll[0] + 1 << "," << mintpaircoll[1] + 1 << std::endl << std::endl;
    for (size_t ievent = 0; ievent < nevents; ++ievent) {
      std::cout << "event " << ievent + 1 << std::endl;
      if ((CollisionW(pmintwall + 1,1) < CollisionP(mintpaircoll[0] + 1,mintpaircoll[1] + 1))||(mintpaircoll[0] == mintpaircoll[1])) {
        //update time
        timeevent = CollisionW(pmintwall + 1,1);
        time_sim += timeevent;
        std::cout << "Wall collision     " << time_sim << std::endl;
        UpdatePositions(timeevent);
        CollisionWall(pmintwall,size_t(CollisionW(pmintwall + 1,2)));
        //update collisions against wall
        mintimecoll = maxtime;
        for (size_t idpart = 0; idpart < Nparticle; ++idpart) {
          if (idpart == pmintwall) {
            for (size_t idcoord = 0; idcoord < ndim; ++idcoord) {
              timecollision = CollisionWall_time(pmintwall,idcoord);
              if (mintimecoll > timecollision) {
                mintimecoll = timecollision;
                coordinate = double(idcoord);
              }
            }
            CollisionW(pmintwall + 1, 1) = mintimecoll;
            CollisionW(pmintwall + 1, 2) = coordinate;
          }
          else {CollisionW(idpart + 1, 1) -= timeevent;}
        }
        //update pair collisions
        for (size_t idpart1 = 0; idpart1 < Nparticle; ++idpart1){
          if (idpart1 == pmintwall) {
            for (size_t idpart2 = 0; idpart2 < Nparticle; ++idpart2){
              if (idpart2 == pmintwall) {continue;}
              CollisionP(pmintwall + 1,idpart2 + 1) = PairCollision_time(pmintwall,idpart2);
              CollisionP(idpart2 + 1,pmintwall + 1) = CollisionP(pmintwall + 1,idpart2 + 1);
            }
          }
          else {
            for (size_t idpart2 = idpart1 + 1; idpart2 < Nparticle; ++idpart2){
              if (idpart2 == pmintwall) {continue;}
              if (CollisionP(idpart1 + 1,idpart2 + 1) == maxtime) {
                CollisionP(idpart1 + 1,idpart2 + 1) = PairCollision_time(idpart1,idpart2);
                CollisionP(idpart2 + 1,idpart1 + 1) = CollisionP(idpart1 + 1,idpart2 + 1);
              }
              else {
                CollisionP(idpart1 + 1,idpart2 + 1) -= timeevent;
                CollisionP(idpart2 + 1,idpart1 + 1) -= timeevent;
              }
            }
          }
        }
        //get new minima: Wall collisions
        pmintwall = 0;
        for (size_t idpart = 1; idpart < Nparticle; ++idpart) {
          if (CollisionW(idpart + 1, 1) < CollisionW(pmintwall + 1,1)) {pmintwall = idpart;}
        }
        std::cout << "updated wall collisions" << std::endl;
        CollisionW.Print();
        std::cout << "new minimum wall collision with particle " << pmintwall + 1 << std::endl;
        //get new minima: Pair collisions
        dmintpaircoll = maxtime;
        for (size_t idpart1 = 0; idpart1 < Nparticle; ++idpart1){
          for (size_t idpart2 = idpart1 + 1; idpart2 < Nparticle; ++idpart2){
            if (CollisionP(idpart1 + 1,idpart2 + 1) < dmintpaircoll) {
              dmintpaircoll = CollisionP(idpart1 + 1,idpart2 + 1);
              mintpaircoll[0] = idpart1;
              mintpaircoll[1] = idpart2;
            }
          }
        }
        std::cout << "updated pair collisions" << std::endl;
        CollisionP.Print();
        std::cout << "new pair with lowest collision time: " << mintpaircoll[0] + 1 << "," << mintpaircoll[1] + 1 << std::endl;
        std::cout << std::endl;
      }
      else {
        //update time
        timeevent = CollisionP(mintpaircoll[0] + 1,mintpaircoll[1] + 1);
        time_sim += timeevent;
        std::cout << "Particle collision     " << time_sim << std::endl;
        UpdatePositions(timeevent);
        PairCollision(mintpaircoll[0],mintpaircoll[1]);
        //update collisions against wall
        for (size_t idpart = 0; idpart < Nparticle; ++idpart) {
          if ((idpart == mintpaircoll[0])||(idpart == mintpaircoll[1])) {
            mintimecoll = maxtime;
            for (size_t idcoord = 0; idcoord < ndim; ++idcoord) {
              timecollision = CollisionWall_time(idpart,idcoord);
              if (mintimecoll > timecollision) {
                mintimecoll = timecollision;
                coordinate = double(idcoord);
              }
            }
            CollisionW(idpart + 1, 1) = mintimecoll;
            CollisionW(idpart + 1, 2) = coordinate;
          }
          else {CollisionW(idpart + 1, 1) -= timeevent;}
        }
        //update pair collisions
        for (size_t idpart1 = 0; idpart1 < Nparticle; ++idpart1){
          if ((idpart1 == mintpaircoll[0])||(idpart1 == mintpaircoll[1])) {
            for (size_t idpart2 = 0; idpart2 < Nparticle; ++idpart2){
              if (idpart1 == idpart2) {continue;}
              CollisionP(idpart1 + 1,idpart2 + 1) = PairCollision_time(idpart1,idpart2);
              CollisionP(idpart2 + 1,idpart1 + 1) = CollisionP(idpart1 + 1,idpart2 + 1);
            }
          }
          else {
            for (size_t idpart2 = idpart1 + 1; idpart2 < Nparticle; ++idpart2){
              if ((idpart2 == mintpaircoll[0])||(idpart2 == mintpaircoll[1])) {continue;}
              if (CollisionP(idpart1 + 1,idpart2 + 1) == maxtime) {
                CollisionP(idpart1 + 1,idpart2 + 1) = PairCollision_time(idpart1,idpart2);
                CollisionP(idpart2 + 1,idpart1 + 1) = CollisionP(idpart1 + 1,idpart2 + 1);
              }
              else {
                CollisionP(idpart1 + 1,idpart2 + 1) -= timeevent;
                CollisionP(idpart2 + 1,idpart1 + 1) -= timeevent;
              }
            }
          }
        }
        //get new minima: Wall collisions
        pmintwall = 0;
        for (size_t idpart = 1; idpart < Nparticle; ++idpart) {
          if (CollisionW(idpart + 1, 1) < CollisionW(pmintwall + 1,1)) {pmintwall = idpart;}
        }
        std::cout << "updated wall collisions" << std::endl;
        CollisionW.Print();
        std::cout << "new minimum wall collision with particle " << pmintwall + 1 << std::endl;
        //get new minima: Pair collisions
        dmintpaircoll = maxtime;
        for (size_t idpart1 = 0; idpart1 < Nparticle; ++idpart1){
          for (size_t idpart2 = idpart1 + 1; idpart2 < Nparticle; ++idpart2){
            if (CollisionP(idpart1 + 1,idpart2 + 1) < dmintpaircoll) {
              dmintpaircoll = CollisionP(idpart1 + 1,idpart2 + 1);
              mintpaircoll[0] = idpart1;
              mintpaircoll[1] = idpart2;
            }
          }
        }
        std::cout << "updated pair collisions" << std::endl;
        CollisionP.Print();
        std::cout << "new pair with lowest collision time: " << mintpaircoll[0] + 1 << "," << mintpaircoll[1] + 1 << std::endl;
        std::cout << "===============" << std::endl;
      }
      if ((snap)||(write2file)) {MD2file(false,Nparticle,ndim,sigma,true,time_sim,pparticles,pvelocities,boxdim,fileevents,folderevents);}
    }
    if (snap) {SnapEvents(fileevents,folderevents);}
  }
  void RunMDTime(double step = 1.0, bool continue_sim = false, bool snap = true, bool write2file = true) {
    //run time driven molecular dynamics simulation
    if (particles.size() != Nparticle) {
      std::cout << "Particle array inconsistent with simulation constructor. Recalculating particle array." << std::endl;
      PlaceParticles();
    }
    if (velocities.size() != Nparticle) {
      std::cout << "Velocity array inconsistent with simulation constructor. Recalculating velocity array with canonical distribution at 298.15 K." << std::endl;
      CanonicalVelocities(298.15);
    }
    std::vector<std::vector<double> > * pparticles = & particles;
    std::vector<std::vector<double> > * pvelocities = & velocities;
    std::vector<int> mintpaircoll(2,0);
    double timecollision = 0.0;
    double mintimecoll = maxtime;
    double timeevent;
    double coordinate;
    double dmintpaircoll;
    int pmintwall = 0;
    int nsteps = int(time_max/step) + 1;
    double initime;
    if (!continue_sim) {
      time_sim = 0.0;
      //determine collision times
      CollisionW.resize(Nparticle,2);
      for (size_t idpart = 0; idpart < Nparticle; ++idpart) {
        mintimecoll = maxtime;
        for (size_t idcoord = 0; idcoord < ndim; ++idcoord) {
          timecollision = CollisionWall_time(idpart,idcoord);
          if (mintimecoll > timecollision) {
            mintimecoll = timecollision;
            coordinate = double(idcoord);
          }
        }
        CollisionW(idpart + 1, 1) = mintimecoll;
        CollisionW(idpart + 1, 2) = coordinate;
        if (mintimecoll < CollisionW(pmintwall + 1,1)) {pmintwall = idpart;}
      }
      CollisionP.resize(Nparticle,Nparticle);
      dmintpaircoll = maxtime;
      for (size_t idpart1 = 0; idpart1 < Nparticle; ++idpart1){
        for (size_t idpart2 = idpart1 + 1; idpart2 < Nparticle; ++idpart2){
          CollisionP(idpart1 + 1,idpart2 + 1) = PairCollision_time(idpart1,idpart2);
          CollisionP(idpart2 + 1,idpart1 + 1) = CollisionP(idpart1 + 1,idpart2 + 1);
          if (CollisionP(idpart1 + 1,idpart2 + 1) < dmintpaircoll) {
            dmintpaircoll = CollisionP(idpart1 + 1,idpart2 + 1);
            mintpaircoll[0] = idpart1;
            mintpaircoll[1] = idpart2;
          }
        }
      }
    }
    //snap initial state?
    if ((snap)||(write2file)) {
      if (!continue_sim) {MD2file(true,Nparticle,ndim,sigma,true,time_sim,pparticles,pvelocities,boxdim,filetimes,foldertimes);}
      MD2file(false,Nparticle,ndim,sigma,true,time_sim,pparticles,pvelocities,boxdim,filetimes,foldertimes);
    }
    std::cout << "wall collisions" << std::endl;
    CollisionW.Print();
    std::cout << "minimum wall collision with particle " << pmintwall + 1 << std::endl;
    std::cout << "pair collisions" << std::endl;
    CollisionP.Print();
    std::cout << "pair with lowest collision time: " << mintpaircoll[0] + 1 << "," << mintpaircoll[1] + 1 << std::endl << std::endl;
    initime = time_sim;
    for (size_t istep = 0; istep < nsteps; ++istep) {
      std::cout << "clock " << time_sim << std::endl;
      std::cout << step << " " << CollisionW(pmintwall + 1,1) << " " << CollisionP(mintpaircoll[0] + 1,mintpaircoll[1] + 1) << std::endl;
      if ((time_sim - initime) > time_max) {break;}
      if ((step < CollisionW(pmintwall + 1,1))&&(step < CollisionP(mintpaircoll[0] + 1,mintpaircoll[1] + 1))) {
        std::cout << "step" << std::endl;
        //no collision
        UpdatePositions(step);
        time_sim += step;
        //update collisions against wall
        for (size_t idpart = 0; idpart < Nparticle; ++idpart) {
          CollisionW(idpart + 1, 1) -= step;
        }
        //update pair collisions
        for (size_t idpart1 = 0; idpart1 < Nparticle; ++idpart1) {
          for (size_t idpart2 = idpart1 + 1; idpart2 < Nparticle; ++idpart2){
            if (CollisionP(idpart1 + 1,idpart2 + 1) == maxtime) {
              CollisionP(idpart1 + 1,idpart2 + 1) = PairCollision_time(idpart1,idpart2);
              CollisionP(idpart2 + 1,idpart1 + 1) = CollisionP(idpart1 + 1,idpart2 + 1);
            }
            else {
              CollisionP(idpart1 + 1,idpart2 + 1) -= step;
              CollisionP(idpart2 + 1,idpart1 + 1) -= step;
            }
          }
        }
        //get new minima: Wall collisions
        pmintwall = 0;
        for (size_t idpart = 1; idpart < Nparticle; ++idpart) {
          if (CollisionW(idpart + 1, 1) < CollisionW(pmintwall + 1,1)) {pmintwall = idpart;}
        }
        std::cout << "updated wall collisions" << std::endl;
        CollisionW.Print();
        std::cout << "new minimum wall collision with particle " << pmintwall + 1 << std::endl;
        //get new minima: Pair collisions
        dmintpaircoll = maxtime;
        for (size_t idpart1 = 0; idpart1 < Nparticle; ++idpart1){
          for (size_t idpart2 = idpart1 + 1; idpart2 < Nparticle; ++idpart2){
            if (CollisionP(idpart1 + 1,idpart2 + 1) < dmintpaircoll) {
              dmintpaircoll = CollisionP(idpart1 + 1,idpart2 + 1);
              mintpaircoll[0] = idpart1;
              mintpaircoll[1] = idpart2;
            }
          }
        }
        std::cout << "updated pair collisions" << std::endl;
        CollisionP.Print();
        std::cout << "new pair with lowest collision time: " << mintpaircoll[0] + 1 << "," << mintpaircoll[1] + 1 << std::endl;
        std::cout << std::endl;
      }
      else {
        std::cout << "collision/crossing" << std::endl;
        --istep;
        if ((CollisionW(pmintwall + 1,1) < CollisionP(mintpaircoll[0] + 1,mintpaircoll[1] + 1))||(mintpaircoll[0] == mintpaircoll[1])) {
          //update time
          timeevent = CollisionW(pmintwall + 1,1);
          time_sim += timeevent;
          std::cout << "Wall collision     " << time_sim << std::endl;
          UpdatePositions(timeevent);
          CollisionWall(pmintwall,size_t(CollisionW(pmintwall + 1,2)));
          //update collisions against wall
          mintimecoll = maxtime;
          for (size_t idpart = 0; idpart < Nparticle; ++idpart) {
            if (idpart == pmintwall) {
              for (size_t idcoord = 0; idcoord < ndim; ++idcoord) {
                timecollision = CollisionWall_time(pmintwall,idcoord);
                if (mintimecoll > timecollision) {
                  mintimecoll = timecollision;
                  coordinate = double(idcoord);
                }
              }
              CollisionW(pmintwall + 1, 1) = mintimecoll;
              CollisionW(pmintwall + 1, 2) = coordinate;
            }
            else {CollisionW(idpart + 1, 1) -= timeevent;}
          }
          //update pair collisions
          for (size_t idpart1 = 0; idpart1 < Nparticle; ++idpart1){
            if (idpart1 == pmintwall) {
              for (size_t idpart2 = 0; idpart2 < Nparticle; ++idpart2){
                if (idpart2 == pmintwall) {continue;}
                CollisionP(pmintwall + 1,idpart2 + 1) = PairCollision_time(pmintwall,idpart2);
                CollisionP(idpart2 + 1,pmintwall + 1) = CollisionP(pmintwall + 1,idpart2 + 1);
              }
            }
            else {
              for (size_t idpart2 = idpart1 + 1; idpart2 < Nparticle; ++idpart2){
                if (idpart2 == pmintwall) {continue;}
                if (CollisionP(idpart1 + 1,idpart2 + 1) == maxtime) {
                  CollisionP(idpart1 + 1,idpart2 + 1) = PairCollision_time(idpart1,idpart2);
                  CollisionP(idpart2 + 1,idpart1 + 1) = CollisionP(idpart1 + 1,idpart2 + 1);
                }
                else {
                  CollisionP(idpart1 + 1,idpart2 + 1) -= timeevent;
                  CollisionP(idpart2 + 1,idpart1 + 1) -= timeevent;
                }
              }
            }
          }
          //get new minima: Wall collisions
          pmintwall = 0;
          for (size_t idpart = 1; idpart < Nparticle; ++idpart) {
            if (CollisionW(idpart + 1, 1) < CollisionW(pmintwall + 1,1)) {pmintwall = idpart;}
          }
          std::cout << "updated wall collisions" << std::endl;
          CollisionW.Print();
          std::cout << "new minimum wall collision with particle " << pmintwall + 1 << std::endl;
          //get new minima: Pair collisions
          dmintpaircoll = maxtime;
          for (size_t idpart1 = 0; idpart1 < Nparticle; ++idpart1){
            for (size_t idpart2 = idpart1 + 1; idpart2 < Nparticle; ++idpart2){
              if (CollisionP(idpart1 + 1,idpart2 + 1) < dmintpaircoll) {
                dmintpaircoll = CollisionP(idpart1 + 1,idpart2 + 1);
                mintpaircoll[0] = idpart1;
                mintpaircoll[1] = idpart2;
              }
            }
          }
          std::cout << "updated pair collisions" << std::endl;
          CollisionP.Print();
          std::cout << "new pair with lowest collision time: " << mintpaircoll[0] + 1 << "," << mintpaircoll[1] + 1 << std::endl;
          std::cout << std::endl;
        }
        else {
          //update time
          timeevent = CollisionP(mintpaircoll[0] + 1,mintpaircoll[1] + 1);
          time_sim += timeevent;
          std::cout << "Particle collision     " << time_sim << std::endl;
          UpdatePositions(timeevent);
          PairCollision(mintpaircoll[0],mintpaircoll[1]);
          //update collisions against wall
          for (size_t idpart = 0; idpart < Nparticle; ++idpart) {
            if ((idpart == mintpaircoll[0])||(idpart == mintpaircoll[1])) {
              mintimecoll = maxtime;
              for (size_t idcoord = 0; idcoord < ndim; ++idcoord) {
                timecollision = CollisionWall_time(idpart,idcoord);
                if (mintimecoll > timecollision) {
                  mintimecoll = timecollision;
                  coordinate = double(idcoord);
                }
              }
              CollisionW(idpart + 1, 1) = mintimecoll;
              CollisionW(idpart + 1, 2) = coordinate;
            }
            else {CollisionW(idpart + 1, 1) -= timeevent;}
          }
          //update pair collisions
          for (size_t idpart1 = 0; idpart1 < Nparticle; ++idpart1){
            if ((idpart1 == mintpaircoll[0])||(idpart1 == mintpaircoll[1])) {
              for (size_t idpart2 = 0; idpart2 < Nparticle; ++idpart2){
                if (idpart1 == idpart2) {continue;}
                CollisionP(idpart1 + 1,idpart2 + 1) = PairCollision_time(idpart1,idpart2);
                CollisionP(idpart2 + 1,idpart1 + 1) = CollisionP(idpart1 + 1,idpart2 + 1);
              }
            }
            else {
              for (size_t idpart2 = idpart1 + 1; idpart2 < Nparticle; ++idpart2){
                if ((idpart2 == mintpaircoll[0])||(idpart2 == mintpaircoll[1])) {continue;}
                if (CollisionP(idpart1 + 1,idpart2 + 1) == maxtime) {
                  CollisionP(idpart1 + 1,idpart2 + 1) = PairCollision_time(idpart1,idpart2);
                  CollisionP(idpart2 + 1,idpart1 + 1) = CollisionP(idpart1 + 1,idpart2 + 1);
                }
                else {
                  CollisionP(idpart1 + 1,idpart2 + 1) -= timeevent;
                  CollisionP(idpart2 + 1,idpart1 + 1) -= timeevent;
                }
              }
            }
          }
          //get new minima: Wall collisions
          pmintwall = 0;
          for (size_t idpart = 1; idpart < Nparticle; ++idpart) {
            if (CollisionW(idpart + 1, 1) < CollisionW(pmintwall + 1,1)) {pmintwall = idpart;}
          }
          std::cout << "updated wall collisions" << std::endl;
          CollisionW.Print();
          std::cout << "new minimum wall collision with particle " << pmintwall + 1 << std::endl;
          //get new minima: Pair collisions
          dmintpaircoll = maxtime;
          for (size_t idpart1 = 0; idpart1 < Nparticle; ++idpart1){
            for (size_t idpart2 = idpart1 + 1; idpart2 < Nparticle; ++idpart2){
              if (CollisionP(idpart1 + 1,idpart2 + 1) < dmintpaircoll) {
                dmintpaircoll = CollisionP(idpart1 + 1,idpart2 + 1);
                mintpaircoll[0] = idpart1;
                mintpaircoll[1] = idpart2;
              }
            }
          }
          std::cout << "updated pair collisions" << std::endl;
          CollisionP.Print();
          std::cout << "new pair with lowest collision time: " << mintpaircoll[0] + 1 << "," << mintpaircoll[1] + 1 << std::endl;
          std::cout << "===============" << std::endl;
        }
      }
      if ((snap)||(write2file)) {MD2file(false,Nparticle,ndim,sigma,true,time_sim,pparticles,pvelocities,boxdim,filetimes,foldertimes);}
    }
    if (snap) {SnapEvents(filetimes,foldertimes);}
  }
};

#endif //_Molecular_Dynamics_
