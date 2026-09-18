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

#ifndef _RM1_
#define _RM1_
#include "MNDO.hpp"

//descrition:
//The RM1 method

class RM1: public MNDO {
  //this is the implementation of RM1
  //G. B. Rocha, R. O. Freire, A. M. Simas, J. J. P. Stewart, J. Comp. Chem., 27(10), 1101, 2005
  //bX quantities are "barred" tensors, which are used only for the open-shell case
public:
  RM1(BSet _bset, Molecule _mol, std::string _openclosed = "0", std::string _corecorrection = "0"): MNDO(_bset,_mol,_openclosed,_corecorrection) {}
  ~RM1() {}
  std::string Type() {return "RM1";}
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
    //function that checks whether the atoms in the molecule are fully defined in the respective MNDO theories
    bool def;
    for (size_t iatm = 0; iatm < Natoms; ++iatm) {
      def = false;
      if (atoms[iatm] == 1) {def = true;}                                   //H
      else if ((atoms[iatm] > 5)&&(atoms[iatm] < 10)) {def = true;}         //C,N,O,F
      else if ((atoms[iatm] > 14)&&(atoms[iatm] < 18)) {def = true;}        //P,S,Cl
      else if ((atoms[iatm] == 35)||(atoms[iatm] == 53)) {def = true;}      //Br,I
      if (!def) {throw("ERROR: RM1.hpp: RM1: checkAtoms(): atom not fully specified for RM1-theory");}
    }
  }
  double ZeroOverlap(size_t atm) {return ZeroOverlapRM1(atm);}
  double ElementFormTheo(size_t atomicnr) {
    //function returning theoretical heats of formation for elements; values directly in a.u.
    double enth = 0.0;
    switch (atomicnr) {
      case 1:      //H
        enth = -0.439546816194;
        break;
      case 6:      //C
        enth = -4.30899631856;
        break;
      case 7:      //N
        enth = -7.50452619284;
        break;
      case 8:      //O
        enth = -11.440778309;
        break;
      case 9:      //F
        enth = -17.8085654062;
        break;
      case 15:     //P
        enth = -4.52307090781;
        break;
      case 16:     //S
        enth = -6.79711669201;
        break;
      case 17:     //Cl
        enth = -14.0555181302;
        break;
      case 35:     //Br
        enth = -13.1237879735;
        break;
      case 53:     //I
        enth = -9.13196215556;
        break;
    }
    return enth;
  }
  double betaA0(size_t atomicnr, int L) {
    //function returning the bonding parameters; values stored in eV, returned however in a.u.
    double betaa0 = 0.0;
    switch (atomicnr) {
      case 1:                         //H
        betaa0 = -5.76544469;
        break;
      case 6:                         //C
        if (L == 0) {betaa0 = -15.45932428;}
        else if (L == 1) {betaa0 = -8.23608638;}
        break;
      case 7:                         //N
        if (L == 0) {betaa0 = -20.87124548;}
        else if (L == 1) {betaa0 = -16.67171853;}
        break;
      case 8:                         //O
        if (L == 0) {betaa0 = -29.85101212;}
        else if (L == 1) {betaa0 = -29.15101314;}
        break;
      case 9:                         //F
        if (L == 0) {betaa0 = -70.00000512;}
        else if (L == 1) {betaa0 = -32.67982711;}
        break;
      case 15:                        //P
        if (L == 0) {betaa0 = -6.13514969;}
        else if (L == 1) {betaa0 = -5.94442127;}
        break;
      case 16:                        //S
        if (L == 0) {betaa0 = -1.95910719;}
        else if (L == 1) {betaa0 = -8.77430652;}
        break;
      case 17:                        //Cl
        if (L == 0) {betaa0 = -19.92430432;}
        else if (L == 1) {betaa0 = -11.52935197;}
        break;
      case 35:                        //Br
        if (L == 0) {betaa0 = -1.34139841;}
        else if (L == 1) {betaa0 = -8.20225991;}
        break;
      case 53:                        //I
        if (L == 0) {betaa0 = -4.19316149;}
        else if (L == 1) {betaa0 = -4.40038412;}
        break;
    }
    return betaa0/au2eV;
  }
  double UlX(size_t atomicnr, int L) {
    //function that returns the atomic potential U
    double ulx = 0.0;
    switch (atomicnr) {
      case 1:                 //H
        ulx = -11.96067697;
        break;
      case 6:                 //C
        if (L == 0) {ulx = -51.72556032;}
        else if (L == 1) {ulx = -39.40728943;}
        break;
      case 7:                 //N
        if (L == 0) {ulx = -70.85123715;}
        else if (L == 1) {ulx = -57.97730920;}
        break;
      case 8:                 //O
        if (L == 0) {ulx = -96.94948069;}
        else if (L == 1) {ulx = -77.89092978;}
        break;
      case 9:                 //F
        if (L == 0) {ulx = -134.18369591;}
        else if (L == 1) {ulx = -107.84660920;}
        break;
      case 15:                //P
        if (L == 0) {ulx = -41.81533184;}
        else if (L == 1) {ulx = -34.38342529;}
        break;
      case 16:                //S
        if (L == 0) {ulx = -55.16775121;}
        else if (L == 1) {ulx = -46.52930422;}
        break;
      case 17:                //Cl
        if (L == 0) {ulx = -118.47306918;}
        else if (L == 1) {ulx = -76.35330340;}
        break;
      case 35:                //Br
        if (L == 0) {ulx = -113.48398183;}
        else if (L == 1) {ulx = -76.18720023;}
        break;
      case 53:                //I
        if (L == 0) {ulx = -74.89997837;}
        else if (L == 1) {ulx = -51.41023805;}
        break;
    }
    return ulx/au2eV;
  }
  double alpha(int atomicnr, int atm2) {
    //function that returns alpha values for AM1
    //alphas stored and returned in 1/Angstrom
    double alpha = 0.0;
    switch (atomicnr) {
      case 1:     //H
        alpha = 3.06835947;
        break;
      case 6:     //C
        alpha = 2.79282078;
        break;
      case 7:     //N
        alpha = 2.96422542;
        break;
      case 8:     //O
        alpha = 4.17196717;
        break;
      case 9:     //F
        alpha = 6.00000062;
        break;
      case 15:    //P
        alpha = 1.90993294;
        break;
      case 16:    //S
        alpha = 2.44015636;
        break;
      case 17:    //Cl
        alpha = 3.69358828;
        break;
      case 35:    //Br
        alpha = 2.86710532;
        break;
      case 53:    //I
        alpha = 2.14157092;
        break;
    }
    return alpha;
  }
  std::vector<double> AM1K(int atomicnr, int atm2 = 0) {
    std::vector<double> am1k;
    switch (atomicnr) {
      case 1:     //H
        am1k.push_back(0.10288875/au2eV);
        am1k.push_back(0.06457449/au2eV);
        am1k.push_back(-0.03567387/au2eV);
        break;
      case 6:     //C
        am1k.push_back(0.07462271/au2eV);
        am1k.push_back(0.01177053/au2eV);
        am1k.push_back(0.03720662/au2eV);
        am1k.push_back(-0.00270657/au2eV);
        break;
      case 7:     //N
        am1k.push_back(0.06073380/au2eV);
        am1k.push_back(0.02438558/au2eV);
        am1k.push_back(-0.02283430/au2eV);
        break;
      case 8:     //O
        am1k.push_back(0.23093552/au2eV);
        am1k.push_back(0.05859873/au2eV);
        break;
      case 9:     //F
        am1k.push_back(0.40302025/au2eV);
        am1k.push_back(0.07085831/au2eV);
        break;
      case 15:    //P
        am1k.push_back(-0.41063467/au2eV);
        am1k.push_back(-0.16299288/au2eV);
        am1k.push_back(-0.04887125/au2eV);
        break;
      case 16:    //S
        am1k.push_back(-0.74601055/au2eV);
        am1k.push_back(-0.06519286/au2eV);
        am1k.push_back(-0.00655977/au2eV);
        break;
      case 17:    //Cl
        am1k.push_back(0.12947108/au2eV);
        am1k.push_back(0.00288899/au2eV);
        break;
      case 35:    //Br
        am1k.push_back(0.98689937/au2eV);
        am1k.push_back(-0.92731247/au2eV);
        break;
      case 53:    //I
        am1k.push_back(-0.08147724/au2eV);
        am1k.push_back(0.05914991/au2eV);
        break;
    }
    return am1k;
  }
  std::vector<double> AM1L(int atomicnr, int atm2 = 0) {
    std::vector<double> am1l;
    switch (atomicnr) {
      case 1:     //H
        am1l.push_back(5.90172268);
        am1l.push_back(6.41785671);
        am1l.push_back(2.80473127);
        break;
      case 6:     //C
        am1l.push_back(5.73921605);
        am1l.push_back(6.92401726);
        am1l.push_back(6.26158944);
        am1l.push_back(9.00003735);
        break;
      case 7:     //N
        am1l.push_back(4.58892946);
        am1l.push_back(4.62730519);
        am1l.push_back(2.05274659);
        break;
      case 8:     //O
        am1l.push_back(5.21828736);
        am1l.push_back(7.42932932);
        break;
      case 9:     //F
        am1l.push_back(7.20441959);
        am1l.push_back(9.00001562);
        break;
      case 15:    //P
        am1l.push_back(6.08752832);
        am1l.push_back(7.09472602);
        am1l.push_back(8.99979308);
        break;
      case 16:    //S
        am1l.push_back(4.81038002);
        am1l.push_back(7.20760864);
        am1l.push_back(9.00000180);
        break;
      case 17:    //Cl
        am1l.push_back(2.97724424);
        am1l.push_back(7.09827589);
        break;
      case 35:    //Br
        am1l.push_back(4.28484191);
        am1l.push_back(4.54005910);
        break;
      case 53:    //I
        am1l.push_back(1.56065072);
        am1l.push_back(5.76111270);
        break;
    }
    return am1l;
  }
  std::vector<double> AM1M(int atomicnr, int atm2 = 0) {
    std::vector<double> am1m;
    switch (atomicnr) {
      case 1:     //H
        am1m.push_back(1.17501185);
        am1m.push_back(1.93844484);
        am1m.push_back(1.63655241);
        break;
      case 6:     //C
        am1m.push_back(1.04396983);
        am1m.push_back(1.66159571);
        am1m.push_back(1.63158721);
        am1m.push_back(2.79557901);
        break;
      case 7:     //N
        am1m.push_back(1.37873881);
        am1m.push_back(2.08370698);
        am1m.push_back(1.86763816);
        break;
      case 8:     //O
        am1m.push_back(0.90363555);
        am1m.push_back(1.51754610);
        break;
      case 9:     //F
        am1m.push_back(0.81653013);
        am1m.push_back(1.43802381);
        break;
      case 15:    //P
        am1m.push_back(1.31650261);
        am1m.push_back(1.90721319);
        am1m.push_back(2.65857780);
        break;
      case 16:    //S
        am1m.push_back(0.59380129);
        am1m.push_back(1.29492008);
        am1m.push_back(1.80060151);
        break;
      case 17:    //Cl
        am1m.push_back(1.46749784);
        am1m.push_back(2.50002723);
        break;
      case 35:    //Br
        am1m.push_back(2.00019696);
        am1m.push_back(2.01617695);
        break;
      case 53:    //I
        am1m.push_back(2.00002063);
        am1m.push_back(2.20488800);
        break;
    }
    return am1m;
  }
  double Dvalue(size_t atmnr, size_t idx) {
    //function returning the D values needed to calculate eris; values stored in Angstrom, but returned in a.u.
    //note that this function returns both D1 and D2; idx is then either 1 or 2
    double dval = 0.0;
    switch (atmnr) {
      case 1:                   //H
        dval = 0.0;
        break;
      case 6:                   //C
        if (idx == 1) {dval = 0.421625705;}
        else if (idx == 2) {dval = 0.366514014;}
        break;
      case 7:                  //N
        if (idx == 1) {dval = 0.343733389;}
        else if (idx == 2) {dval = 0.327636952;}
        break;
      case 8:                  //O
        if (idx == 1) {dval = 0.258593059;}
        else if (idx == 2) {dval = 0.253799433;}
        break;
      case 9:                  //F
        if (idx == 1) {dval = 0.184626083;}
        else if (idx == 2) {dval = 0.244715022;}
        break;
      case 15:                 //P
        if (idx == 1) {dval = 0.534837011;}
        else if (idx == 2) {dval = 0.507940817;}
        break;
      case 16:                 //S
        if (idx == 1) {dval = 0.525839200;}
        else if (idx == 2) {dval = 0.472356644;}
        break;
      case 17:                 //Cl
        if (idx == 1) {dval = 0.240341095;}
        else if (idx == 2) {dval = 0.467043700;}
        break;
      case 35:                 //Br
        if (idx == 1) {dval = 0.111074578;}
        else if (idx == 2) {dval = 0.552580723;}
        break;
      case 53:                 //I
        if (idx == 1) {dval = 0.685994904;}
        else if (idx == 2) {dval = 0.586643919;}
        break;
    }
    return dval*dist_Angstrom2aum1;
  }
  double rho(size_t atmnr, size_t l) {
    //function returning the rho values needed to calculate eris; values stored in Angstrom but returned in atomic units
    double rho = 0.0;
    switch (atmnr) {
      case 1:                  //H
        if (l == 0) {rho = 0.51488046;}
        else if (l == 1) {rho = 0.0;}
        else if (l == 2) {rho = 0.0;}
        break;
      case 6:                  //C
        if (l == 0) {rho = 0.55156780;}
        else if (l == 1) {rho = 0.51905507;}
        else if (l == 2) {rho = 0.39998984;}
        break;
      case 7:                  //N
        if (l == 0) {rho = 0.55012484;}
        else if (l == 1) {rho = 0.27259716;}
        else if (l == 2) {rho = 0.32975880;}
        break;
      case 8:                  //O
        if (l == 0) {rho = 0.51417391;}
        else if (l == 1) {rho = 0.26285897;}
        else if (l == 2) {rho = 0.29551626;}
        break;
      case 9:                  //F
        if (l == 0) {rho = 0.43057954;}
        else if (l == 1) {rho = 0.28680176;}
        else if (l == 2) {rho = 0.42425093;}
        break;
      case 15:                 //P
        if (l == 0) {rho = 0.64975614;}
        else if (l == 1) {rho = 0.67287406;}
        else if (l == 2) {rho = 0.84709086;}
        break;
      case 16:                 //S
        if (l == 0) {rho = 0.57651500;}
        else if (l == 1) {rho = 0.38175294;}
        else if (l == 2) {rho = 0.53350588;}
        break;
      case 17:                 //Cl
        if (l == 0) {rho = 0.46872232;}
        else if (l == 1) {rho = 0.35326854;}
        else if (l == 2) {rho = 0.34329568;}
        break;
      case 35:                 //Br
        if (l == 0) {rho = 0.42064959;}
        else if (l == 1) {rho = 0.20096728;}
        else if (l == 2) {rho = 0.45124871;}
        break;
      case 53:                 //I
        if (l == 0) {rho = 0.35998881;}
        else if (l == 1) {rho = 0.71175720;}
        else if (l == 2) {rho = 0.75323162;}
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
        eri = 13.98321296;                                                                       //(ss|ss)
        break;
      case 6:       //C
        if ((Lbra == 0)&&(Lket == 0)) {eri = 13.05312440;}                                       //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 10.95113739;}                                  //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 11.33479389;}    //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 9.72395099;}   //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 1.55215133;}                                   //(sp|sp)||(ps|ps)
        break;
      case 7:       //N
        if ((Lbra == 0)&&(Lket == 0)) {eri = 13.08736234;}                                       //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 13.69924324;}                                  //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 13.21226834;}    //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 11.94103953;}  //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 5.00000846;}                                   //(sp|sp)||(ps|ps)
        break;
      case 8:       //O
        if ((Lbra == 0)&&(Lket == 0)) {eri = 14.00242788;}                                       //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 14.14515138;}                                  //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 14.95625043;}    //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 12.70325497;}  //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 3.93217161;}                                   //(sp|sp)||(ps|ps)
        break;
      case 9:       //F
        if ((Lbra == 0)&&(Lket == 0)) {eri = 16.72091319;}                                       //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 15.22581028;}                                  //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 16.76142629;}    //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 14.86578679;}  //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 1.99766171;}                                   //(sp|sp)||(ps|ps)
        break;
      case 15:      //P
        if ((Lbra == 0)&&(Lket == 0)) {eri = 11.08059265;}                                       //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 7.60417563;}                                   //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 5.68339201;}     //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 7.40265182;}   //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 1.16181792;}                                   //(sp|sp)||(ps|ps)
        break;
      case 16:      //S
        if ((Lbra == 0)&&(Lket == 0)) {eri = 12.48828408;}                                       //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 8.52301167;}                                   //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 8.56910574;}     //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 7.66863296;}   //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 3.88978932;}                                   //(sp|sp)||(ps|ps)
        break;
      case 17:      //Cl
        if ((Lbra == 0)&&(Lket == 0)) {eri = 15.36023105;}                                       //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 12.56502640;}                                  //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 13.30671171;}    //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 9.66397083;}   //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 1.76489897;}                                   //(sp|sp)||(ps|ps)
        break;
      case 35:      //Br
        if ((Lbra == 0)&&(Lket == 0)) {eri = 17.11563074;}                                       //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 10.73546293;}                                  //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 15.62419251;}    //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 8.86056199;}   //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 2.23512762;}                                   //(sp|sp)||(ps|ps)
        break;
      case 53:      //I
        if ((Lbra == 0)&&(Lket == 0)) {eri = 19.99974131;}                                       //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 7.30488343;}                                   //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 7.68957672;}     //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 6.85424614;}   //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 1.41602940;}                                   //(sp|sp)||(ps|ps)
        break;
    }
    return eri/au2eV;
  }
};

#endif //_RM1_
