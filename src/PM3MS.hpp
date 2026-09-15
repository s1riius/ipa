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

#ifndef _PM3_
#define _PM3_
#include "MNDO.hpp"

//description:
//The PM3-MS method

class PM3ms: public MNDO {
  //this is the implementation of PM3-MS for sugars
  //J. Y. Mane, M. Klobukowski, Chem. Phys. Lett., 500, 140, 2010
  //bX quantities are "barred" tensors, which are used only for the open-shell case
public:
  PM3ms(BSet _bset, Molecule _mol, std::string _openclosed = "0", std::string _corecorrection = "0"): MNDO(_bset,_mol,_openclosed,_corecorrection) {}
  ~PM3ms() {}
  std::string Type() {return "PM3ms";}
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
      if (atoms[iatm] == 1) {def = true;}                            //H
      else if ((atoms[iatm] > 5)&&(atoms[iatm] < 9)) {def = true;}   //C,N,O
      if (!def) {throw("ERROR: PM3ms.hpp: PM3ms: checkAtoms(): atom not fully specified for PM3ms-theory");}
    }
  }
  double ZeroOverlap(size_t atm) {return ZeroOverlapPM3(atm);}
  double ElementFormTheo(size_t atomicnr) {
    //function returning theoretical heats of formation for elements; values directly in a.u.
    double enth = 0.0;
    switch (atomicnr) {
      case 1:      //H
        enth = -0.480435734286;
        break;
      case 6:      //C
        enth = -4.05540034771;
        break;
      case 7:      //N
        enth = -5.77479194152;
        break;
      case 8:      //O
        enth = -10.6102009213;
        break;
    }
    return enth;
  }
  double betaA0(size_t atomicnr, int L) {
    //function returning the bonding parameters; values stored in eV, returned however in a.u.
    double betaa0 = 0.0;
    switch (atomicnr) {
      case 1:         //H
        betaa0 = -4.289542;
        break;
      case 6:         //C
        if (L == 0) {betaa0 = -12.110500;}
        else if (L == 1) {betaa0 = -9.549442;}
        break;
      case 7:         //N
        if (L == 0) {betaa0 = -14.0625210;}
        else if (L == 1) {betaa0 = -20.0438480;}
        break;
      case 8:         //O
        if (L == 0) {betaa0 = -41.899362;}
        else if (L == 1) {betaa0 = -28.450886;}
        break;
    }
    return betaa0/au2eV;
  }
  double UlX(size_t atomicnr, int L) {
    //function that returns the atomic potential U
    double ulx = 0.0;
    switch (atomicnr) {
      case 1:                 //H
        ulx = -12.436734;
        break;
      case 6:                 //C
        if (L == 0) {ulx = -44.415413;}
        else if (L == 1) {ulx = -37.396380;}
        break;
      case 7:                 //N
        if (L == 0) {ulx = -49.3356720;}
        else if (L == 1) {ulx = -47.5097360;}
        break;
      case 8:                 //O
        if (L == 0) {ulx = -82.154826;}
        else if (L == 1) {ulx = -74.660258;}
        break;
    }
    return ulx/au2eV;
  }
  double alpha(int atomicnr, int atm2) {
    //function that returns alpha values for PM3ms
    //alphas stored and returned in 1/Angstrom
    double alpha = 0.0;
    switch (atomicnr) {
      case 1:     //H
        alpha = 3.118553;
        break;
      case 6:     //C
        alpha = 2.814684;
        break;
      case 7:     //N
        alpha = 2.830545;
        break;
      case 8:     //O
        alpha = 3.071855;
        break;
    }
    return alpha;
  }
  std::vector<double> AM1K(int atomicnr, int atm2 = 0) {
    std::vector<double> am1k;
    switch (atomicnr) {
      case 1:     //H
        am1k.push_back(1.127411/au2eV);
        am1k.push_back(-1.0603290/au2eV);
        break;
      case 6:     //C
        am1k.push_back(0.048141/au2eV);
        am1k.push_back(0.050727/au2eV);
        break;
      case 7:     //N
        am1k.push_back(1.5016740/au2eV);
        am1k.push_back(-1.5057720/au2eV);
        break;
      case 8:     //O
        am1k.push_back(-1.1311280/au2eV);
        am1k.push_back(1.1378910/au2eV);
        break;
    }
    return am1k;
  }
  std::vector<double> AM1L(int atomicnr, int atm2 = 0) {
    std::vector<double> am1l;
    switch (atomicnr) {
      case 1:     //H
        am1l.push_back(5.0941340);
        am1l.push_back(6.0068720);
        break;
      case 6:     //C
        am1l.push_back(4.1429530);
        am1l.push_back(5.2812140);
        break;
      case 7:     //N
        am1l.push_back(5.9011480);
        am1l.push_back(6.0046580);
        break;
      case 8:     //O
        am1l.push_back(6.0024770);
        am1l.push_back(5.9505120);
        break;
    }
    return am1l;
  }
  std::vector<double> AM1M(int atomicnr, int atm2 = 0) {
    std::vector<double> am1m;
    switch (atomicnr) {
      case 1:     //H
        am1m.push_back(1.5374650);
        am1m.push_back(1.5701890);
        break;
      case 6:     //C
        am1m.push_back(1.5810810);
        am1m.push_back(0.8905620);
        break;
      case 7:     //N
        am1m.push_back(1.7107400);
        am1m.push_back(1.7161490);
        break;
      case 8:     //O
        am1m.push_back(1.6073110);
        am1m.push_back(1.5983950);
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
        if (l == 0) {rho = 0.486640184440;}
        else if (l == 1) {rho = 0.0;}
        else if (l == 2) {rho = 0.0;}
        break;
      case 6:                  //C
        if (l == 0) {rho = 0.642767927834;}
        else if (l == 1) {rho = 0.449532462965;}
        else if (l == 2) {rho = 0.345972966858;}
        break;
      case 7:                  //N
        if (l == 0) {rho = 0.604753079729;}
        else if (l == 1) {rho = 0.525917048522;}
        else if (l == 2) {rho = 0.359254597636;}
        break;
      case 8:                  //O
        if (l == 0) {rho = 0.456941201522;}
        else if (l == 1) {rho = 0.499258635325;}
        else if (l == 2) {rho = 0.323472582688;}
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
        eri = 14.7942080;                                                                        //(ss|ss)
        break;
      case 6:       //C
        if ((Lbra == 0)&&(Lket == 0)) {eri = 11.2007080;}                                        //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 10.7962920;}                                   //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 10.2650270;}     //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 9.0425660;}    //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 2.2909800;}                                    //(sp|sp)||(ps|ps)
        break;
      case 7:       //N
        if ((Lbra == 0)&&(Lket == 0)) {eri = 11.9047870;}                                        //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 11.7546720;}                                   //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 7.3485650;}      //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 10.8072770;}   //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 1.1367130;}                                    //(sp|sp)||(ps|ps)
        break;
      case 8:       //O
        if ((Lbra == 0)&&(Lket == 0)) {eri = 15.7557600;}                                        //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 13.6540160;}                                   //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 10.6211600;}     //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 12.4060950;}   //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 0.5938830;}                                    //(sp|sp)||(ps|ps)
        break;
    }
    return eri/au2eV;
  }
};

#endif //_PM3_
