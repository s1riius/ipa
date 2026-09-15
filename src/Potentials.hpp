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

#ifndef _Model_Potentials_
#define _Model_Potentials_
#include <vector>
#include <string>
#include <math.h>
#include <iostream>
#include "BSet.hpp"
#include "Molecule.hpp"
#include "ConstantsPackage.hpp"
#include "UnitConversion.hpp"
#include "math/FunctionPackage.hpp"

//description:
//classes to calculate molecular energies according to a given model potential

class HardSphere {
  //class for the hard-sphere model potential between point particles
protected:
  double epsilon;                           //for the hard-sphere this is the energy for distances below sigma; in Joule
  double etot;                              //in Joule
  double sigma;                             //in Angstrom
  Molecule mol;
  std::vector<size_t> atoms;                //container with all atoms
  std::vector<double> pairenergies;         //container with pair energies
public:
  HardSphere(double _sigma = 1.0, double _epsilon = 1.0e100) {
    epsilon = _epsilon;
    sigma = _sigma;
    etot = 0.0;
  }
  HardSphere(Molecule _mol, double _sigma = 1.0, double _epsilon = 1.0e100): mol(_mol) {
    epsilon = _epsilon;
    sigma = _sigma;
    atoms = mol.Atoms();
    int natoms = atoms.size();
    pairenergies.resize(natoms*(natoms - 1)/2);
    etot = 0.0;
  }
  ~HardSphere() {}
  //getters
  virtual std::string Type() {return "HS";}
  std::string MechanicsType() {return "classical";}
  Molecule & Component() {return mol;}
  double Epsilon() {return epsilon;}
  double Sigma() {return sigma;}
  double getEnergy(bool complete = 0) {return etot;}
  std::vector<double> getPairEnergies() {return pairenergies;}
  double getMayer(double beta) {
    //function to get products of Mayer functions
    double mayer = 1.0;
    for (size_t idx = 0; idx < pairenergies.size(); ++idx) {
      mayer *= exp(-beta*pairenergies[idx]) - 1.0;
    }
    return mayer;
  }
  void setAOBasis(const BSet & _bset) {}
  //setters
  void setMolecule(const Molecule & _component) {
    mol = _component;
    atoms = mol.Atoms();
    int natoms = atoms.size();
    pairenergies.resize(natoms*(natoms - 1)/2);
  }
  void setSigma(double _sigma) {sigma = _sigma;}
  void setEpsilon(double _epsilon) {epsilon = _epsilon;}
  virtual double Potential(int atm1, int atm2, std::vector<double> & pos1, std::vector<double> & pos2) {
    //calculation of energy
    double dist2 = 0.0;
    for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
      dist2 += (pos1[idcoord] - pos2[idcoord])*(pos1[idcoord] - pos2[idcoord]);
    }
    double dist = sqrt(dist2);
    double potentialenergy = 0.0;
    if (dist < sigma) {potentialenergy = epsilon;}
    return potentialenergy;
  }
  void Calculate(int print = 1, double thresh_dens = 1e-7, double thresh_en = 1e-7, size_t maxiter = 200) {
    //calculation of the energy of a hard-sphere system
    matrixE geometry = mol.Geometry();
    std::vector<double> posA(3,0.0);
    std::vector<double> posB(3,0.0);
    size_t Natoms = geometry.rows();
    etot = 0.0;
    int counter = 0;
    for (size_t idatm = 0; idatm < Natoms; ++idatm) {
      //set posA
      for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
        posA[idcoord] = geometry(idatm + 1,idcoord + 1);
      }
      for (size_t idbtm = idatm + 1; idbtm < Natoms; ++idbtm, ++counter) {
        //set posB
        for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
          posB[idcoord] = geometry(idbtm + 1,idcoord + 1);
        }
        pairenergies[counter] = Potential(atoms[idatm],atoms[idbtm],posA,posB);
        etot += pairenergies[counter];
      }
    }
    if (print == 1) {std::cout << "New total energy: " << etot << std::endl;}
  }
  void Collapse2Point() {
    //function that collapses the molecule to a point particle with coordinates of center of mass and pseudoparticle is used for atomic number
    matrixE geometry(1,3);
    std::vector<double> com = mol.CM();
    geometry(1,1) = com[0];
    geometry(1,2) = com[1];
    geometry(1,3) = com[2];
    mol.setGeometry(geometry);
    atoms.resize(1);
    atoms[0] = 0;
    mol.setAtoms(atoms);
  }
  void CalcSigma(int level = 0) {
    //function calculating the molecular radius, only make sense if the molecule is a single molecule, not a composite one
    //level determines how accurate we want to calculate the radius: 
    //         level = 0 looks for the largest distance to center of mass and double it
    //         level = 1 looks for the largest distance between two atoms
    matrixE geometry = mol.Geometry();
    size_t Natoms = geometry.rows();
    double aux;
    double dist;
    sigma = 0.0;
    if (level == 0) {
      std::vector<double> com = mol.CM();
      for (size_t idatm = 0; idatm < Natoms; ++idatm) {
        aux = 0.0;
        for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
          aux += (geometry(idatm + 1,idcoord + 1) - com[idcoord])*(geometry(idatm + 1,idcoord + 1) - com[idcoord]);
        }
        dist = sqrt(aux);
        if (dist > sigma) {sigma = dist;}
      }
      sigma *= 2.0;
    }
    else if (level == 1) {
      for (size_t idatm = 0; idatm < Natoms; ++idatm) {
        for (size_t idbtm = idatm; idbtm < Natoms; ++idbtm) {
          aux = 0.0;
          for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
            aux += (geometry(idatm + 1,idcoord + 1) - geometry(idbtm + 1,idcoord + 1))*(geometry(idatm + 1,idcoord + 1) - geometry(idbtm + 1,idcoord + 1));
          }
          dist = sqrt(aux);
          if (dist > sigma) {sigma = dist;}
        }
      }
    }
  }
  virtual void CalcParam(matrixE & Potential, matrixE & param, bool setparam = true, int mode = 0, double temperature = 298.15) {
    //function that calculates model potential parameters and eventually sets parameters
    //assume that the matrix Potential contains coordinates in the first column and energies (or interaction energies) in the other coordinates
    //that is why param is also a matrix, organized in a way that each column contains a set of parameters
    //the potential energy should come in Eh, distances are assumed to be in Angstrom
    //mode = 0 -> determine sigma by inversion of signal of model potential
    //mode = 1 -> determine sigma by checking largest distance towards center of mass
    //mode = 2 -> determine sigma by determining the largest interatomic distance in a molecule
    int npotentials = Potential.cols() - 1;
    if (npotentials > 0) {
      //there is something to adjust
      param.resize(1,npotentials);
      size_t npts = Potential.rows();
      matrixE Emin(1,npotentials);
      if (mode == 0) {
        for (size_t ipot = 0; ipot < npotentials; ++ipot) {
          param(1,ipot + 1) = 0.0;
          Emin(1,ipot + 1) = 10000000000000000.0;   //just a very large number
          for (size_t ipt = 0; ipt < npts - 1; ++ipt) {
            if (Potential(ipt + 1,ipot + 2)*Potential(ipt + 2,ipot + 2) < 0.0) {
              param(1,ipot + 1) = 0.5*(Potential(ipt + 1,1) + Potential(ipt + 2,1));
            }
            if (Potential(ipt + 1,ipot + 2) < Potential(ipt + 2,ipot + 2)) {
              //if the energy starts going up, then we are at the minimum
              Emin(1,ipot + 1) = Potential(ipt + 1,ipot + 2);
              break;
            }
          }
        }
        if (setparam) {
          //there is only one sigma, therefore we make a Boltzmann average over all the sigma
          double emin = Emin(1,1);
          //get the absolute minimum of all potentials
          for (size_t ipot = 1; ipot < npotentials; ++ipot) {
            if (Emin(1,ipot + 1) < emin) {emin = Emin(1,ipot + 1);}
          }
          //now Boltzmann average
          double conversionfactor = Energy(1.0,"J/mol")/(R*temperature);
          sigma = 0.0;
          double norm = 0.0;
          double auxval;
          for (size_t ipot = 0; ipot < npotentials; ++ipot) {
            auxval = exp(-(Emin(1,ipot + 1) - emin)*conversionfactor);
            sigma += param(1,ipot + 1)*auxval;
            norm += auxval;
          }
          sigma /= norm;
        }
      }
    }
    else if ((mode == 1)||(mode == 2)) {CalcSigma(mode - 1);}
  }
  virtual double B2(double temp = 298.15, double threshold = 1.0e-7, int maxiter = 500) {
    //function that calculates the second-virial coefficient in m^3/mol
    double sigma3 = sigma*sigma*sigma*1.0e-30*NA;
    double b2 = 2.0*pi*sigma3/3.0;
    return b2;
  }
  virtual double B3(double temp = 298.15, double threshold = 1.0e-7, int maxiter = 500) {
    //function that calculates the third-virial coefficient in m^6/mol^2
    double sigma3 = sigma*sigma*sigma*1.0e-30*NA;
    double b3 = 5.0*pi*pi*sigma3*sigma3/18.0;
    return b3;
  }
  virtual double B4(double temp = 298.15, double threshold = 1.0e-7, int maxiter = 500) {
    //function that calculates the fourth-virial coefficient in m^9/mol^3
    double sigma3 = sigma*sigma*sigma*1.0e-30*NA;
    double b4 = 0.2869495*8.0*pi*pi*pi*sigma3*sigma3*sigma3/27.0;
    return b4;
  }
  virtual double B5(double temp = 298.15, double threshold = 1.0e-7, int maxiter = 500) {
    //function that calculates the fifth-virial coefficient in m^12/mol^4
    double sigma3 = sigma*sigma*sigma*1.0e-30*NA;
    double sigma6 = sigma3*sigma3;
    double b5 = 0.1102521*16.0*pi*pi*pi*pi*sigma6*sigma6/81.0;
    return b5;
  }
  virtual double B6(double temp = 298.15, double threshold = 1.0e-7, int maxiter = 500) {
    //function that calculates the sixth-virial coefficient in m^15/mol^5
    double sigma3 = sigma*sigma*sigma*1.0e-30*NA;
    double sigma6 = sigma3*sigma3;
    double b6 = 0.03880855*32.0*pi*pi*pi*pi*pi*sigma6*sigma6*sigma3/243.0;
    return b6;
  }
  virtual double B7(double temp = 298.15, double threshold = 1.0e-7, int maxiter = 500) {
    //function that calculates the seventh-virial coefficient in m^18/mol^6
    double sigma3 = sigma*sigma*sigma*1.0e-30*NA;
    double sigma9 = sigma3*sigma3*sigma3;
    double b7 = 0.01307170*64.0*pi*pi*pi*pi*pi*pi*sigma9*sigma9/729.0;
    return b7;
  }
  virtual double B8(double temp = 298.15, double threshold = 1.0e-7, int maxiter = 500) {
    //function that calculates the eighth-virial coefficient in m^21/mol^7
    double sigma3 = sigma*sigma*sigma*1.0e-30*NA;
    double sigma9 = sigma3*sigma3*sigma3;
    double b8 = 0.0043210*128.0*pi*pi*pi*pi*pi*pi*pi*sigma9*sigma9*sigma3/2187.0;
    return b8;
  }
  virtual double B9(double temp = 298.15, double threshold = 1.0e-7, int maxiter = 500) {
    //function that calculates the nineth-virial coefficient in m^24/mol^8
    double sigma3 = sigma*sigma*sigma*1.0e-30*NA;
    double sigma12 = sigma3*sigma3*sigma3*sigma3;
    double b9 = 0.0013092*256.0*pi*pi*pi*pi*pi*pi*pi*pi*sigma12*sigma12/6561.0;
    return b9;
  }
  virtual double B10(double temp = 298.15, double threshold = 1.0e-7, int maxiter = 500) {
    //function that calculates the eighth-virial coefficient in m^27/mol^9
    double sigma3 = sigma*sigma*sigma*1.0e-30*NA;
    double sigma12 = sigma3*sigma3*sigma3*sigma3;
    double b10 = 0.000403515*512.0*pi*pi*pi*pi*pi*pi*pi*pi*pi*sigma12*sigma12*sigma3/19683.0;
    return b10;
  }
};
class LennardJonesSimple: public HardSphere {
  //class containing the Lennard-Jones potential between point particles
protected:
  int aattrct;           //exponent for attractive part of the potential
  int arepuls;           //exponent for repulsive part of the potential
public:
  LennardJonesSimple(Molecule _mol, double _sigma = 1.0, double _epsilon = 1.0e100, int _arepuls = 12, int _aattrct = 6): HardSphere(_mol,_sigma,_epsilon) {
    aattrct = _aattrct;
    arepuls = _arepuls;
  }
  ~LennardJonesSimple() {}
  std::string Type() {return "LJ";}
  int AttractionExp() {return aattrct;}
  int RepulsionExp() {return arepuls;}
  void setAttractionExp(int _aa) {aattrct = _aa;}
  void setRepulsionExp(int _aa) {arepuls = _aa;}
  double Potential(int atm1, int atm2, std::vector<double> & pos1, std::vector<double> & pos2) {
    //calculation of energy
    double dist = 0.0;
    for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
      dist += (pos1[idcoord] - pos2[idcoord])*(pos1[idcoord] - pos2[idcoord]);
    }
    double aux = sigma*sigma/dist;                //\frac{{\sigma}^{2}}{{r}^{2}}
    double aux2 = aux*aux*aux;                    //\frac{{\sigma}^{6}}{{r}^{6}}
    double potentialenergy = 4.0*epsilon*aux2*(aux2 - 1.0);
    return potentialenergy;
  }
  double LJsimpleEval(double RAB, double sig, double eps) {
    //simple evaluation of the potential
    double x = sig/RAB;
    double xx = x*x;
    x = xx*xx*xx;
    return 4.0*eps*(x*x - x);
  }
  void CalcParam(matrixE & Potential, matrixE & param, bool setparam = true, int mode = 0, double temperature = 298.15) {
    //further documentation comments on mother class
    //mode = 0 -> determine sigma by inversion of signal of model potential and epsilon by mathing area of attractive part
    //mode = 1 -> determine sigma by checking largest distance towards center of mass
    //mode = 2 -> determine sigma by determining the largest interatomic distance in a molecule
    int npotentials = Potential.cols() - 1;
    if (npotentials > 0) {
      //there is something to adjust
      param.resize(2,npotentials);        //row 1 -> sigma; row 2 -> epsilon
      size_t npts = Potential.rows();
      int firstpos;               //first position for integration
      int lastpos;                //last position for integration
      double LJint;
      double LJfp;
      double LJlp;
      matrixE Emin(1,npotentials);
      if (mode == 0) {
        //step is supposed to be constant
        double step = Potential(2,1) - Potential(1,1);
        for (size_t ipot = 0; ipot < npotentials; ++ipot) {
          param(1,ipot + 1) = 0.0;
          param(2,ipot + 1) = 0.0;
          Emin(1,ipot + 1) = 10000000000000000.0;   //just a very large number
          firstpos = 0;
          LJint = 0.0;
          LJfp = 0.0;
          LJlp = 0.0;
          for (size_t ipt = 0; ipt < npts - 1; ++ipt) {
            //sigma
            if (Potential(ipt + 1,ipot + 2)*Potential(ipt + 2,ipot + 2) < 0.0) {
              param(1,ipot + 1) = 0.5*(Potential(ipt + 1,1) + Potential(ipt + 2,1));
            }
            //epsilon
            if (Potential(ipt + 1,ipot + 2) < 0.0) {
              //if the energy is negative, then use for fitting; note that Lennard-Jones has sigma in agreement with
              //the model potential, meaning that the zero of the potentials match
              param(2,ipot + 1) += Potential(ipt + 1,ipot + 2);
              LJlp = LJsimpleEval(Potential(ipt + 1,1),param(1,ipot + 1),1.0);
              LJint += LJlp;
              lastpos = ipt;
              if (ipt > 0) {
                if (Potential(ipt,ipot + 2) >= 0.0) {
                  firstpos = ipt;
                  LJfp = LJlp;
                }
              }
              if ((Potential(ipt + 1,ipot + 2) < Potential(ipt + 2,ipot + 2))&&(Potential(ipt + 1,ipot + 2) < Emin(1,ipot + 1))) {Emin(1,ipot + 1) = Potential(ipt + 1,ipot + 2);}
            }
          }
          //fix the first and last positions
          if ((lastpos == npts - 2)&&(Potential(npts,ipot + 2) < 0.0)) {
            lastpos = npts - 1;
            LJlp = LJsimpleEval(Potential(lastpos + 1,1),param(1,ipot + 1),1.0);
            LJint += LJlp;
            param(2,ipot + 1) += Potential(lastpos + 1,ipot + 2);
          }
          //finalize integration
          param(2,ipot + 1) *= 2.0;
          LJint *= 2.0;
          param(2,ipot + 1) -= (Potential(firstpos + 1,ipot + 2) + Potential(lastpos + 1,ipot + 2));
          LJint -= (LJfp + LJlp); 
          param(2,ipot + 1) *= 0.5*step;
          LJint *= 0.5*step;
          //now get epsilon
          param(2,ipot + 1) /= LJint;
        }
        std::cout << "----" << std::endl;
        Emin.Print(12);
        std::cout << "----" << std::endl;
        if (setparam) {
          //there is only one sigma, therefore we make a Boltzmann average over all the sigma
          double emin = Emin(1,1);
          //get the absolute minimum of all potentials
          for (size_t ipot = 1; ipot < npotentials; ++ipot) {
            if (Emin(1,ipot + 1) < emin) {emin = Emin(1,ipot + 1);}
          }
          //now Boltzmann average
          std::cout << std::setprecision(26);
          double conversionfactor = Energy(1.0,"J/mol")/(R*temperature);
          sigma = 0.0;
          epsilon = 0.0;
          double norm = 0.0;
          double auxval;
          for (size_t ipot = 0; ipot < npotentials; ++ipot) {
            auxval = exp(-(Emin(1,ipot + 1) - emin)*conversionfactor);
            std::cout << (Emin(1,ipot + 1) - emin) << " " << conversionfactor << std::endl;
            sigma += param(1,ipot + 1)*auxval;
            epsilon += param(2,ipot + 1)*auxval;
            std::cout << param(1,ipot + 1) << "*" << auxval << "       " << param(2,ipot + 1) << "*" << auxval << std::endl;
            norm += auxval;
          }
          sigma /= norm;
          epsilon *= Energy(1.0,"J")/norm;
        }
      }
    }
    else if ((mode == 1)||(mode == 2)) {CalcSigma(mode - 1);}
  }
  double B2(double temp = 298.15, double threshold = 1.0e-7, int maxiter = 500) {
    //function that calculates the second-virial coefficient in m^3/mol
    double sigma3 = sigma*sigma*sigma*1.0e-30*NA;
    double hs = 2.0*pi*sigma3/3.0;
    double ebeta = epsilon/(KB*temp);
    double ebeta14 = pow(ebeta,0.25);
    double ebeta34 = ebeta14*ebeta14*ebeta14;
    double gamma14 = 3.6256099082219083119;
    double gamma34 = sqrt(2.0)*pi/gamma14;
    double b2 = hs*sqrt(2.0)*(ebeta14*gamma34*((1.0 + 4.0*ebeta)*hyper1F1(0.75,0.5,ebeta) - 6.0*ebeta*hyper1F1(7.0/4.0,1.5,ebeta)) + 0.5*ebeta34*gamma14*((4.0*ebeta - 1.0)*hyper1F1(5.0/4.0,1.5,ebeta) - 10.0*ebeta*hyper1F1(9.0/4.0,2.5,ebeta)/3.0));
    return b2;
  }
  double B3(double temp = 298.15, double threshold = 1.0e-7, int maxiter = 500) {
    //function that calculates the third-virial coefficient in m^6/mol^2
    return 0.0;
  }
  double B4(double temp = 298.15, double threshold = 1.0e-7, int maxiter = 500) {
    //function that calculates the fourth-virial coefficient in m^9/mol^3
    return 0.0;
  }
  double B5(double temp = 298.15, double threshold = 1.0e-7, int maxiter = 500) {
    //function that calculates the fifth-virial coefficient in m^12/mol^4
    return 0.0;
  }
  double B6(double temp = 298.15, double threshold = 1.0e-7, int maxiter = 500) {
    //function that calculates the sixth-virial coefficient in m^15/mol^5
    return 0.0;
  }
  double B7(double temp = 298.15, double threshold = 1.0e-7, int maxiter = 500) {
    //function that calculates the seventh-virial coefficient in m^18/mol^6
    return 0.0;
  }
  double B8(double temp = 298.15, double threshold = 1.0e-7, int maxiter = 500) {
    //function that calculates the eighth-virial coefficient in m^21/mol^7
    return 0.0;
  }
  double B9(double temp = 298.15, double threshold = 1.0e-7, int maxiter = 500) {
    //function that calculates the nineth-virial coefficient in m^24/mol^8
    return 0.0;
  }
  double B10(double temp = 298.15, double threshold = 1.0e-7, int maxiter = 500) {
    //function that calculates the eighth-virial coefficient in m^27/mol^9
    return 0.0;
  }
};
class Sutherland: public HardSphere {
  //class containing the Sutherland potential between point particles
public:
  Sutherland(Molecule _mol, double _sigma = 1.0, double _epsilon = 1.0e100): HardSphere(_mol,_sigma,_epsilon) {}
  ~Sutherland() {}
  std::string Type() {return "Sutherland";}
  double Potential(int atm1, int atm2, std::vector<double> & pos1, std::vector<double> & pos2, double infinity = 1.0e100) {
    //calculation of energy
    double dist = 0.0;
    for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
      dist += (pos1[idcoord] - pos2[idcoord])*(pos1[idcoord] - pos2[idcoord]);
    }
    double potentialenergy = infinity;
    if (dist >= sigma*sigma) {
      double aux = sigma*sigma/dist;                //\frac{{\sigma}^{2}}{{r}^{2}}
      double aux2 = aux*aux*aux;                    //\frac{{\sigma}^{6}}{{r}^{6}}
      potentialenergy = -epsilon*aux2;
    }
    return potentialenergy;
  }
  double B2(double temp = 298.15, double threshold = 1.0e-7, int maxiter = 500) {
    //function that calculates the second-virial coefficient in m^3/mol
    double sigma3 = sigma*sigma*sigma*1.0e-30*NA;
    double hs = 2.0*pi*sigma3/3.0;
    double ebeta = epsilon/(KB*temp);
    double expebeta = exp(ebeta);
    double b2 = hs*(expebeta - 2.0*ebeta*hyper1F1(0.5,1.5,ebeta));
    return b2;
  }
  double B3(double temp = 298.15, double threshold = 1.0e-7, int maxiter = 500) {
    //function that calculates the third-virial coefficient in m^6/mol^2
    return 0.0;
  }
  double B4(double temp = 298.15, double threshold = 1.0e-7, int maxiter = 500) {
    //function that calculates the fourth-virial coefficient in m^9/mol^3
    return 0.0;
  }
  double B5(double temp = 298.15, double threshold = 1.0e-7, int maxiter = 500) {
    //function that calculates the fifth-virial coefficient in m^12/mol^4
    return 0.0;
  }
  double B6(double temp = 298.15, double threshold = 1.0e-7, int maxiter = 500) {
    //function that calculates the sixth-virial coefficient in m^15/mol^5
    return 0.0;
  }
  double B7(double temp = 298.15, double threshold = 1.0e-7, int maxiter = 500) {
    //function that calculates the seventh-virial coefficient in m^18/mol^6
    return 0.0;
  }
  double B8(double temp = 298.15, double threshold = 1.0e-7, int maxiter = 500) {
    //function that calculates the eighth-virial coefficient in m^21/mol^7
    return 0.0;
  }
  double B9(double temp = 298.15, double threshold = 1.0e-7, int maxiter = 500) {
    //function that calculates the nineth-virial coefficient in m^24/mol^8
    return 0.0;
  }
  double B10(double temp = 298.15, double threshold = 1.0e-7, int maxiter = 500) {
    //function that calculates the eighth-virial coefficient in m^27/mol^9
    return 0.0;
  }
};
class EllipsoidalConstraint {
  //class that evaluates constraint potential of ellipsoidal form
  //only parameters evaluated
  double radii[3];                        //ellipsoid's radii in Angstroem
  double maxrad;
  double anisotropy[3];
  double exponent;                        //exponent in the case of power potential
  double beta;                            //beta for log-Fermi
  double temperature;                     //temperature scaling parameter
  double center[3];                       //origin of the potential
  double KBEh;                            //Boltzmann constant in Hartree
  int potentialtype;                      //type of potential to be used: 0 -> log-Fermi (default)
                                          //                              1 -> power
public:
  //radii must be positive
  EllipsoidalConstraint(double r1 = 1.0, double r2 = -1.0, double r3 = -1.0) {
    radii[0] = r1;
    if (r2 < 0.0) {radii[1] = r1;}
    else {radii[1] = r2;}
    if (r3 < 0.0) {radii[2] = r1;}
    else {radii[2] = r3;}
    maxrad = r1;
    if (r2 > maxrad) {maxrad = r2;}
    if (r3 > maxrad) {maxrad = r3;}
    for (size_t idx = 0; idx < 3; ++idx) {
      anisotropy[idx] = maxrad/radii[idx];
    }
    maxrad *= dist_Angstrom2aum1;
    exponent = 30.0;
    beta = 6.0;
    temperature = 300.0;
    center[0] = 0.0;
    center[1] = 0.0;
    center[2] = 0.0;
    potentialtype = 0;
    KBEh = KB/au2J;
  }
  ~EllipsoidalConstraint() {}
  //getters
  std::vector<double> Radii() {
    std::vector<double> rad(3,radii[0]);
    rad[1] = radii[1];
    rad[2] = radii[2];
    return rad;
  }
  double Exponent() {return exponent;}
  double Beta() {return beta;}
  double Temperature() {return temperature;}
  std::vector<double> Center() {
    std::vector<double> centre(3,center[0]);
    centre[1] = center[1];
    centre[2] = center[2];
    return centre;
  }
  int TypeOfPotential() {return potentialtype;}
  //setters
  void setRadii(double r1, double r2 = -1.0, double r3 = -1.0) {
    radii[0] = r1;
    if (r2 < 0.0) {radii[1] = r1;}
    else {radii[1] = r2;}
    if (r3 < 0.0) {radii[2] = r1;}
    else {radii[2] = r3;}
    maxrad = r1;
    if (r2 > maxrad) {maxrad = r2;}
    if (r3 > maxrad) {maxrad = r3;}
    for (size_t idx = 0; idx < 3; ++idx) {
      anisotropy[idx] = maxrad/radii[idx];
    }
    maxrad *= dist_Angstrom2aum1;
  }
  void setExponent(double newexp) {exponent = newexp;}
  void setBeta(double newb) {beta = newb;}
  void setTemperature(double Temp) {temperature = Temp;}
  void setCenter(double c1, double c2, double c3) {
    center[0] = c1;
    center[1] = c2;
    center[2] = c3;
  }
  void setTypeOfPotential(int ptype) {potentialtype = ptype;}
  //other functions
  void LogFermiPotential(double & result, const matrixE & geometry, matrixE & grad) {
    //the log-Fermi potential
    double rAO[3];
    double RAO;
    double expterm;
    double RAOANIS[3];
    double Fermi;
    size_t Natoms = geometry.rows();
    result = 0.0;
    for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
      expterm = 0.0;
      for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
        rAO[idcoord] = anisotropy[idcoord]*(geometry(idAtm + 1,idcoord + 1) - center[idcoord])*dist_Angstrom2aum1;
        RAOANIS[idcoord] = rAO[idcoord]*anisotropy[idcoord];
        expterm += rAO[idcoord]*rAO[idcoord];
      }
      RAO = sqrt(expterm);
      expterm = exp(beta*(RAO - maxrad));
      Fermi = dist_Angstrom2aum1*KBEh*temperature/(1.0 + expterm);
      result += log(1.0 + expterm);
      grad(3*idAtm + 1,1) += beta*expterm*Fermi*RAOANIS[0]/(RAO + 1.0e-14);
      grad(3*idAtm + 2,1) += beta*expterm*Fermi*RAOANIS[1]/(RAO + 1.0e-14);
      grad(3*idAtm + 3,1) += beta*expterm*Fermi*RAOANIS[2]/(RAO + 1.0e-14);
    }
    result *= KBEh*temperature;
  }
  void PowerPotential(double & result, const matrixE & geometry, matrixE & grad) {
    //the power potential function
    double RAO2;
    double polyterm;
    double RAOANIS[3];
    size_t Natoms = geometry.rows();
    result = 0.0;
    for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
      RAO2 = 0.0;
      for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
        polyterm = anisotropy[idcoord]*(geometry(idAtm + 1,idcoord + 1) - center[idcoord])*dist_Angstrom2aum1;
        RAOANIS[idcoord] = polyterm*anisotropy[idcoord];
        RAO2 += polyterm*polyterm;
      }
      polyterm = pow(RAO2/(maxrad*maxrad),0.5*exponent);
      result += polyterm;
      grad(3*idAtm + 1,1) += exponent*polyterm*RAOANIS[0]/(RAO2 + 1.0e-14);
      grad(3*idAtm + 2,1) += exponent*polyterm*RAOANIS[1]/(RAO2 + 1.0e-14);
      grad(3*idAtm + 3,1) += exponent*polyterm*RAOANIS[2]/(RAO2 + 1.0e-14);
    }
  }
  double Calculate(const matrixE & geometry, matrixE & grad) {
    //function that evaluates potential according to setup
    double result;
    if (potentialtype == 0) {LogFermiPotential(result,geometry,grad);}
    else if (potentialtype == 1) {PowerPotential(result,geometry,grad);}
    return result;
  }
};

#endif //_Model_Potentials_
