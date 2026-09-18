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

#ifndef _PM3BP_
#define _PM3BP_
#include "MNDO.hpp"

//descrition:
//The PM3-BasePair method

class PM3BP: public MNDO {
  //this is the implementation of PM3BP for DNA/RNA base pairs
  //T. J. Giese, E. C. Sherer, C. J. Cramer, D. M. York, J. Chem. Theory Comput., 1, 1275, 2005
  //bX quantities are "barred" tensors, which are used only for the open-shell case
public:
  PM3BP(BSet _bset, Molecule _mol, std::string _openclosed = "0", std::string _corecorrection = "0"): MNDO(_bset,_mol,_openclosed,_corecorrection) {}
  ~PM3BP() {}
  std::string Type() {return "PM3BP";}
  double gfactor(size_t iatm1, size_t iatm2, double RAB) {return 1.0;}
  virtual double gfactor_dR(size_t iatm1, size_t iatm2) {return 0.0;}
  double AM1factor(double RAB, std::vector<double> & K1, std::vector<double> & L1, std::vector<double> & M1) {
    //RAB in a.u.
    double AM1factor = 0.0;
    RAB *= dist_Angstrom2au;
    for (size_t idx = 0; idx < K1.size(); ++idx) {
      AM1factor += K1[idx]*exp(-L1[idx]*(RAB - M1[idx])*(RAB - M1[idx]))/RAB;
    }
    return AM1factor;
  }
  virtual double AM1factor_dR(double RAB, std::vector<double> & K1, std::vector<double> & L1, std::vector<double> & M1) {
    //RAB in a.u.
    double am1factor_dr = 0.0;
    double RR = RAB*dist_Angstrom2au;
    for (size_t idx = 0; idx < K1.size(); ++idx) {
      am1factor_dr += K1[idx]*(2.0*RR*RR*L1[idx] - 2.0*L1[idx]*M1[idx]*RR + 1.0)*exp(-L1[idx]*(RR - M1[idx])*(RR - M1[idx]))/(RAB*RR);
    }
    return am1factor_dr;  
  }
  virtual double AM1factor_dR2(double RAB, std::vector<double> & K1, std::vector<double> & L1, std::vector<double> & M1) {
    //RAB in a.u.
    double am1factor_dr2 = 0.0;
    double RR = RAB*dist_Angstrom2au;
    for (size_t idx = 0; idx < K1.size(); ++idx) {
      am1factor_dr2 += 2.0*K1[idx]*exp(-L1[idx]*(RR - M1[idx])*(RR - M1[idx]))*(L1[idx]*M1[idx]*RR - 1.0 - RR*L1[idx]*(RR - M1[idx])*(2.0*RR*RR*L1[idx] - 2.0*L1[idx]*M1[idx]*RR + 1.0))/(RAB*RAB*RR);
    }
    return am1factor_dr2;  
  }
  //other auxiliary functions
  void checkAtoms() {
    //function that checks whether the atoms in the molecule are fully defined in PM3BP
    bool def;
    for (size_t iatm = 0; iatm < Natoms; ++iatm) {
      def = false;
      if (atoms[iatm] == 1) {def = true;}                            //H
      else if ((atoms[iatm] > 5)&&(atoms[iatm] < 9)) {def = true;}   //C,N,O
      if (!def) {throw("ERROR: PM3BP.hpp: PM3BP: checkAtoms(): atom not fully specified for PM3BP-theory");}
    }
  }
  double ZeroOverlap(size_t atm) {return ZeroOverlapPM3(atm);}
  double ElementFormTheo(size_t atomicnr) {
    //function returning theoretical heats of formation for elements; values directly in a.u.
    double enth = 0.0;
    switch (atomicnr) {
      case 1:      //H
        enth = -0.46874201793;
        break;
      case 6:      //C
        enth = -4.05540034771;
        break;
      case 7:      //N
        enth = -5.75406193855;
        break;
      case 8:      //O
        enth = -10.6892450762;
        break;
    }
    return enth;
  }
  double betaA0(size_t atomicnr, int L) {
    //function returning the bonding parameters; values stored in eV, returned however in a.u.
    double betaa0 = 0.0;
    switch (atomicnr) {
      case 1:         //H
        betaa0 = -4.87823460;
        break;
      case 6:         //C
        if (L == 0) {betaa0 = -11.91001500;}
        else if (L == 1) {betaa0 = -9.80275500;}
        break;
      case 7:         //N
        if (L == 0) {betaa0 = -14.33884665;}
        else if (L == 1) {betaa0 = -19.30862853;}
        break;
      case 8:         //O
        if (L == 0) {betaa0 = -46.87741023;}
        else if (L == 1) {betaa0 = -24.74232518;}
        break;
    }
    return betaa0/au2eV;
  }
  double UlX(size_t atomicnr, int L) {
    //function that returns the atomic potential U
    double ulx = 0.0;
    switch (atomicnr) {
      case 1:                 //H
        ulx = -12.75511880;
        break;
      case 6:                 //C
        if (L == 0) {ulx = -47.27032000;}
        else if (L == 1) {ulx = -36.2669180;}
        break;
      case 7:                 //N
        if (L == 0) {ulx = -48.79493385;}
        else if (L == 1) {ulx = -46.57945903;}
        break;
      case 8:                 //O
        if (L == 0) {ulx = -87.38709324;}
        else if (L == 1) {ulx = -71.70268570;}
        break;
    }
    return ulx/au2eV;
  }
  double alpha(int atomicnr, int atm2) {
    //function that returns alpha values for PM3BP
    //alphas stored and returned in 1/Angstrom
    double alpha = 0.0;
    switch (atomicnr) {
      case 1:     //H
        alpha = 3.35638600;
        break;
      case 6:     //C
        alpha = 2.70780700;
        break;
      case 7:     //N
        alpha = 2.83054500;
        break;
      case 8:     //O
        alpha = 3.21710200;
        break;
    }
    return alpha;
  }
  std::vector<double> AM1K(int atomicnr, int atm2 = 0) {
    std::vector<double> am1k;
    switch (atomicnr) {
      case 1:     //H
        am1k.push_back(1.12172594/au2eV);
        am1k.push_back(-1.06492525/au2eV);
        break;
      case 6:     //C
        am1k.push_back(0.05010700/au2eV);
        am1k.push_back(0.05073300/au2eV);
        break;
      case 7:     //N
        am1k.push_back(1.50167153/au2eV);
        am1k.push_back(-1.51571618/au2eV);
        break;
      case 8:     //O
        am1k.push_back(-1.13117677/au2eV);
        am1k.push_back(1.13098909/au2eV);
        break;
    }
    return am1k;
  }
  std::vector<double> AM1L(int atomicnr, int atm2 = 0) {
    std::vector<double> am1l;
    switch (atomicnr) {
      case 1:     //H
        am1l.push_back(5.09516707);
        am1l.push_back(6.02315366);
        break;
      case 6:     //C
        am1l.push_back(6.00316500);
        am1l.push_back(6.00297900);
        break;
      case 7:     //N
        am1l.push_back(5.90399175);
        am1l.push_back(5.97579498);
        break;
      case 8:     //O
        am1l.push_back(6.00999815);
        am1l.push_back(5.87216545);
        break;
    }
    return am1l;
  }
  std::vector<double> AM1M(int atomicnr, int atm2 = 0) {
    std::vector<double> am1m;
    switch (atomicnr) {
      case 1:     //H
        am1m.push_back(1.53693700);
        am1m.push_back(1.57130732);
        break;
      case 6:     //C
        am1m.push_back(1.64221400);
        am1m.push_back(0.89248800);
        break;
      case 7:     //N
        am1m.push_back(1.71042669);
        am1m.push_back(1.71093513);
        break;
      case 8:     //O
        am1m.push_back(1.60731100);
        am1m.push_back(1.60347421);
        break;
    }
    return am1m;
  }
  double Dvalue(size_t atmnr, size_t idx) {
    //function returning the D values needed to calculate eris; values stored in Angstrom, but returned in a.u.
    //note that this function returns both D1 and D2; idx is then either 1 or 2
    double dval = 0.0;
    switch (atmnr) {
      case 1:                  //H
        dval = 0.0;
      case 6:                  //C
        if (idx == 1) {dval = 0.4409314053780;}
        else if (idx == 2) {dval = 0.3517837786516;}
        break;
      case 7:                  //N
        if (idx == 1) {dval = 0.3480401674091;}
        else if (idx == 2) {dval = 0.2801137638434;}
        break;
      case 8:                  //O
        if (idx == 1) {dval = 0.2162309620795;}
        else if (idx == 2) {dval = 0.2712423725348;}
        break;
    }
    return dval*dist_Angstrom2aum1;
  }
  double rho(size_t atmnr, size_t l) {
    //function returning the rho values needed to calculate eris; values stored in Angstrom but returned in atomic units
    double rho = 0.0;
    switch (atmnr) {
      case 1:                  //H
        if (l == 0) {rho = 0.4792425114721178;}
        else if (l == 1) {rho = 0.0;}
        else if (l == 2) {rho = 0.0;}
        break;
      case 6:                  //C
        if (l == 0) {rho = 0.6428006131604468;}
        else if (l == 1) {rho = 0.4495426970001134;}
        else if (l == 2) {rho = 0.34597896187068383;}
        break;
      case 7:                  //N
        if (l == 0) {rho = 0.5799189343842014;}
        else if (l == 1) {rho = 0.5259276440914536;}
        else if (l == 2) {rho = 0.2321791302258449;}
        break;
      case 8:                  //O
        if (l == 0) {rho = 0.4717592829251376;}
        else if (l == 1) {rho = 0.4992677822133302;}
        else if (l == 2) {rho = 0.32420165686630137;}
        break;
    }
    return rho*dist_Angstrom2aum1;
  }
  double eri1Center(int atmnr, int Lbra, int Lket) {
    //function that gives back the semi-empirical 1-center eris
    //values stored in eV but returned in a.u.
    //Lbra is the sum of azimuthal quantum numbers for bra (ss = 0; pp = 2; sp = 1)
    //Lket is the sum of azimuthal quantum numbers for ket (ss = 0; pp = 2; sp = 1; p*p* = -2)
    double eri = 0.0;
    switch (atmnr) {
      case 1:       //H
        eri = 15.02333745;                                                                        //(ss|ss)
        break;
      case 6:       //C
        if ((Lbra == 0)&&(Lket == 0)) {eri = 11.20070800;}                                        //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 10.79629200;}                                   //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 10.26502700;}     //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 9.04256600;}    //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 2.29098000;}                                    //(sp|sp)||(ps|ps)
        break;
      case 7:       //N
        if ((Lbra == 0)&&(Lket == 0)) {eri = 12.41522141;}                                        //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 13.96611412;}                                   //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 7.34540226;}      //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 10.41021925;}   //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 1.13671300;}                                    //(sp|sp)||(ps|ps)
        break;
      case 8:       //O
        if ((Lbra == 0)&&(Lket == 0)) {eri = 15.26164345;}                                        //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 13.65930075;}                                   //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 10.41332625;}     //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 12.42051017;}   //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 0.59388300;}                                    //(sp|sp)||(ps|ps)
        break;
    }
    return eri/au2eV;
  }
};

#endif //_PM3BP_
