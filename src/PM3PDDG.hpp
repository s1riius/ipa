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

#ifndef _PM3PDDG_
#define _PM3PDDG_
#include "MNDO.hpp"

//descrition:
//The PM3-PDDG method

class PM3PDDG: public MNDO {
  //this is the implementation of Jorgensen's PM3PDDG
  //M. P. Repasky, J. Chandrasekhar, W. L. Jorgensen, J. Comput. Chem., 23, 1601, 2002
  //I. Tubert-Brohman, C. R. W. Guimaraes, M. P. Repasky, W. L. Jorgensen, J. Comput. Chem., 25, 138, 2004
  //I. Tubert-Brohman, C. R. W. Guimaraes, W. L. Jorgensen, J. Chem. Theory Comput., 1, 817, 2005
  //bX quantities are "barred" tensors, which are used only for the open-shell case
public:
  PM3PDDG(BSet _bset, Molecule _mol, std::string _openclosed = "0", std::string _corecorrection = "0"): MNDO(_bset,_mol,_openclosed,_corecorrection) {
    doPDDG = true;
  }
  ~PM3PDDG() {}
  std::string Type() {return "PM3PDDG";}
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
      else if ((atoms[iatm] > 13)&&(atoms[iatm] < 18)) {def = true;}        //Si,P,S,Cl
      else if ((atoms[iatm] == 35)||(atoms[iatm] == 53)) {def = true;}      //Br,I
      if (!def) {throw("ERROR: PM3PDDG.hpp: PM3PDDG: checkAtoms(): atom not fully specified for PM3PDDG-theory");}
    }
  }
  double ZeroOverlap(size_t atm) {return ZeroOverlapPM3PDDG(atm);}
  double ElementFormTheo(size_t atomicnr) {
    //function returning theoretical heats of formation for elements; values directly in a.u.
    double enth = 0.0;
    switch (atomicnr) {
      case 1:      //H
        enth = -13.120566;
        break;
      case 6:      //C
        enth = -113.428242;
        break;
      case 7:      //N
        enth = -158.416205;
        break;
      case 8:      //O
        enth = -292.188766;
        break;
      case 9:      //F
        enth = -442.457133;
        break;
      case 14:     //Si
        enth = -66.839000;
        break;
      case 15:     //P
        enth = -117.212854;
        break;
      case 16:     //S
        enth = -166.336554;
        break;
      case 17:     //Cl
        enth = -305.715201;
        break;
      case 35:     //Br
        enth = -351.013887;
        break;
      case 53:     //I
        enth = -291.537869;
        break;
    }
    return enth/au2eV;
  }
  double betaA0(size_t atomicnr, int L) {
    //function returning the bonding parameters; values stored in eV, returned however in a.u.
    double betaa0 = 0.0;
    switch (atomicnr) {
      case 1:         //H
        betaa0 = -6.152654;
        break;
      case 6:         //C
        if (L == 0) {betaa0 = -11.952818;}
        else if (L == 1) {betaa0 = -9.922411;}
        break;
      case 7:         //N
        if (L == 0) {betaa0 = -14.117230;}
        else if (L == 1) {betaa0 = -19.938509;}
        break;
      case 8:         //O
        if (L == 0) {betaa0 = -44.874553;}
        else if (L == 1) {betaa0 = -24.601939;}
        break;
      case 9:         //F
        if (L == 0) {betaa0 = -50.937301;}
        else if (L == 1) {betaa0 = -31.636976;}
        break;
      case 14:        //Si
        if (L == 0) {betaa0 = -3.376445;}
        else if (L == 1) {betaa0 = -3.620969;}
        break;
      case 15:        //P
        if (L == 0) {betaa0 = -12.676297;}
        else if (L == 1) {betaa0 = -7.093318;}
        break;
      case 16:        //S
        if (L == 0) {betaa0 = -2.953912;}
        else if (L == 1) {betaa0 = -8.507779;}
        break;
      case 17:        //Cl
        if (L == 0) {betaa0 = -26.913129;}
        else if (L == 1) {betaa0 = -14.991178;}
        break;
      case 35:        //Br
        if (L == 0) {betaa0 = -21.538044;}
        else if (L == 1) {betaa0 = -8.524764;}
        break;
      case 53:        //I
        if (L == 0) {betaa0 = -16.592621;}
        else if (L == 1) {betaa0 = -6.599816;}
        break;
    }
    return betaa0/au2eV;
  }
  double UlX(size_t atomicnr, int L) {
    //function that returns the atomic potential U
    double ulx = 0.0;
    switch (atomicnr) {
      case 1:                 //H
        ulx = -12.893272;
        break;
      case 6:                 //C
        if (L == 0) {ulx = -48.241241;}
        else if (L == 1) {ulx = -36.461256;}
        break;
      case 7:                 //N
        if (L == 0) {ulx = -49.454546;}
        else if (L == 1) {ulx = -47.757406;}
        break;
      case 8:                 //O
        if (L == 0) {ulx = -87.412505;}
        else if (L == 1) {ulx = -72.183070;}
        break;
      case 9:                 //F
        if (L == 0) {ulx = -111.400432;}
        else if (L == 1) {ulx = -106.395264;}
        break;
      case 14:                //Si
        if (L == 0) {ulx = -26.332522;}
        else if (L == 1) {ulx = -22.602540;}
        break;
      case 15:                //P
        if (L == 0) {ulx = -37.882113;}
        else if (L == 1) {ulx = -30.312979;}
        break;
      case 16:                //S
        if (L == 0) {ulx = -43.906366;}
        else if (L == 1) {ulx = -43.461348;}
        break;
      case 17:                //Cl
        if (L == 0) {ulx = -95.094434;}
        else if (L == 1) {ulx = -53.921651;}
        break;
      case 35:                //Br
        if (L == 0) {ulx = -115.841963;}
        else if (L == 1) {ulx = -74.205146;}
        break;
      case 53:                //I
        if (L == 0) {ulx = -97.664174;}
        else if (L == 1) {ulx = -61.167137;}
        break;
    }
    return ulx/au2eV;
  }
  double alpha(int atomicnr, int atm2) {
    //function that returns alpha values for PM3PDDG
    //alphas stored and returned in 1/Angstrom
    double alpha = 0.0;
    switch (atomicnr) {
      case 1:     //H
        alpha = 3.381686;
        break;
      case 6:     //C
        alpha = 2.725772;
        break;
      case 7:     //N
        alpha = 2.849124;
        break;
      case 8:     //O
        alpha = 3.225309;
        break;
      case 9:     //F
        alpha = 3.200571;
        break;
      case 14:    //Si
        alpha = 2.215157;
        break;
      case 15:    //P
        alpha = 2.005294;
        break;
      case 16:    //S
        alpha = 2.539751;
        break;
      case 17:    //Cl
        alpha = 2.497617;
        break;
      case 35:    //Br
        alpha = 2.424673;
        break;
      case 53:    //I
        alpha = 1.978170;
        break;
    }
    return alpha;
  }
  std::vector<double> AM1K(int atomicnr, int atm2 = 0) {
    std::vector<double> am1k;
    switch (atomicnr) {
      case 1:     //H
        am1k.push_back(1.122244/au2eV);
        am1k.push_back(-1.069737/au2eV);
        break;
      case 6:     //C
        am1k.push_back(0.048906/au2eV);
        am1k.push_back(0.047697/au2eV);
        break;
      case 7:     //N
        am1k.push_back(1.513320/au2eV);
        am1k.push_back(-1.511892/au2eV);
        break;
      case 8:     //O
        am1k.push_back(-1.138455/au2eV);
        am1k.push_back(1.146007/au2eV);
        break;
      case 9:     //F
        am1k.push_back(-0.008079/au2eV);
        am1k.push_back(-0.002659/au2eV);
        break;
      case 14:    //Si
        am1k.push_back(-0.071314/au2eV);
        am1k.push_back(0.089451/au2eV);
        break;
      case 15:    //P
        am1k.push_back(-0.398055/au2eV);
        am1k.push_back(-0.079653/au2eV);
        break;
      case 16:    //S
        am1k.push_back(-0.330692/au2eV);
        am1k.push_back(0.024171/au2eV);
        break;
      case 17:    //Cl
        am1k.push_back(-0.112222/au2eV);
        am1k.push_back(-0.013061/au2eV);
        break;
      case 35:    //Br
        am1k.push_back(0.961362/au2eV);
        am1k.push_back(-0.948834/au2eV);
        break;
      case 53:    //I
        am1k.push_back(-0.136003/au2eV);
        am1k.push_back(-0.037287/au2eV);
        break;
    }
    return am1k;
  }
  std::vector<double> AM1L(int atomicnr, int atm2 = 0) {
    std::vector<double> am1l;
    switch (atomicnr) {
      case 1:     //H
        am1l.push_back(4.707790);
        am1l.push_back(5.857995);
        break;
      case 6:     //C
        am1l.push_back(5.765340);
        am1l.push_back(5.973721);
        break;
      case 7:     //N
        am1l.push_back(5.904394);
        am1l.push_back(6.030014);
        break;
      case 8:     //O
        am1l.push_back(6.000043);
        am1l.push_back(5.963494);
        break;
      case 9:     //F
        am1l.push_back(5.938969);
        am1l.push_back(5.925105);
        break;
      case 14:    //Si
        am1l.push_back(6.000000);
        am1l.push_back(6.000000);
        break;
      case 15:    //P
        am1l.push_back(1.997272);
        am1l.push_back(1.998360);
        break;
      case 16:    //S
        am1l.push_back(6.000000);
        am1l.push_back(6.000000);
        break;
      case 17:    //Cl
        am1l.push_back(5.963719);
        am1l.push_back(1.999556);
        break;
      case 35:    //Br
        am1l.push_back(6.013600);
        am1l.push_back(5.976329);
        break;
      case 53:    //I
        am1l.push_back(3.852912);
        am1l.push_back(5.229264);
        break;
    }
    return am1l;
  }
  std::vector<double> AM1M(int atomicnr, int atm2 = 0) {
    std::vector<double> am1m;
    switch (atomicnr) {
      case 1:     //H
        am1m.push_back(1.547099);
        am1m.push_back(1.567893);
        break;
      case 6:     //C
        am1m.push_back(1.682232);
        am1m.push_back(0.894406);
        break;
      case 7:     //N
        am1m.push_back(1.728376);
        am1m.push_back(1.734108);
        break;
      case 8:     //O
        am1m.push_back(1.622362);
        am1m.push_back(1.614788);
        break;
      case 9:     //F
        am1m.push_back(1.863949);
        am1m.push_back(2.388864);
        break;
      case 14:    //Si
        am1m.push_back(0.237995);
        am1m.push_back(1.897728);
        break;
      case 15:    //P
        am1m.push_back(0.950073);
        am1m.push_back(2.336959);
        break;
      case 16:    //S
        am1m.push_back(0.823837);
        am1m.push_back(2.017756);
        break;
      case 17:    //Cl
        am1m.push_back(1.027719);
        am1m.push_back(2.286377);
        break;
      case 35:    //Br
        am1m.push_back(2.340445);
        am1m.push_back(2.348745);
        break;
      case 53:    //I
        am1m.push_back(1.697455);
        am1m.push_back(2.768669);
        break;
    }
    return am1m;
  }
  double Dvalue(size_t atmnr, size_t idx) {
    //function returning the D values needed to calculate eris; values stored and returned in a.u.
    //note that this function returns both D1 and D2; idx is then either 1 or 2
    double dval = 0.0;
    switch (atmnr) {
      case 1:                  //H
        dval = 0.0;
      case 6:                  //C
        if (idx == 1) {dval = 0.831413;}
        else if (idx == 2) {dval = 0.663222;}
        break;
      case 7:                  //N
        if (idx == 1) {dval = 0.654855;}
        else if (idx == 2) {dval = 0.526924;}
        break;
      case 8:                  //O
        if (idx == 1) {dval = 0.403741;}
        else if (idx == 2) {dval = 0.528360;}
        break;
      case 9:                  //F
        if (idx == 1) {dval = 0.246601;}
        else if (idx == 2) {dval = 0.482551;}
        break;
      case 14:                 //Si
        if (idx == 1) {dval = 1.310515;}
        else if (idx == 2) {dval = 1.126089;}
        break;
      case 15:                 //P
        if (idx == 1) {dval = 0.893978;}
        else if (idx == 2) {dval = 0.960457;}
        break;
      case 16:                 //S
        if (idx == 1) {dval = 1.006989;}
        else if (idx == 2) {dval = 0.891487;}
        break;
      case 17:                 //Cl
        if (idx == 1) {dval = 0.827561;}
        else if (idx == 2) {dval = 0.732427;}
        break;
      case 35:                 //Br
        if (idx == 1) {dval = 0.473860;}
        else if (idx == 2) {dval = 0.968214;}
        break;
      case 53:                 //I
        if (idx == 1) {dval = 0.407261;}
        else if (idx == 2) {dval = 1.062574;}
        break;
    }
    return dval;
  }
  double rho(size_t atmnr, size_t l) {
    //function returning the rho values needed to calculate eris; values stored and returned in atomic units
    double rho = 0.0;
    switch (atmnr) {
      case 1:                  //H
        if (l == 0) {rho = 0.919616;}
        else if (l == 1) {rho = 0.0;}
        else if (l == 2) {rho = 0.0;}
        break;
      case 6:                  //C
        if (l == 0) {rho = 1.214657;}
        else if (l == 1) {rho = 0.848467;}
        else if (l == 2) {rho = 0.652785;}
        break;
      case 7:                  //N
        if (l == 0) {rho = 1.142818;}
        else if (l == 1) {rho = 0.991235;}
        else if (l == 2) {rho = 0.676704;}
        break;
      case 8:                  //O
        if (l == 0) {rho = 0.863494;}
        else if (l == 1) {rho = 0.936266;}
        else if (l == 2) {rho = 0.624291;}
        break;
      case 9:                  //F
        if (l == 0) {rho = 1.296126;}
        else if (l == 1) {rho = 0.634633;}
        else if (l == 2) {rho = 0.805802;}
        break;
      case 14:                 //Si
        if (l == 0) {rho = 2.695556;}
        else if (l == 1) {rho = 1.630757;}
        else if (l == 2) {rho = 0.949200;}
        break;
      case 15:                 //P
        if (l == 0) {rho = 1.743870;}
        else if (l == 1) {rho = 1.050851;}
        else if (l == 2) {rho = 1.208907;}
        break;
      case 16:                 //S
        if (l == 0) {rho = 1.517625;}
        else if (l == 1) {rho = 0.711672;}
        else if (l == 2) {rho = 0.754336;}
        break;
      case 17:                 //Cl
        if (l == 0) {rho = 0.849590;}
        else if (l == 1) {rho = 0.696164;}
        else if (l == 2) {rho = 2.299104;}
        break;
      case 35:                 //Br
        if (l == 0) {rho = 0.853330;}
        else if (l == 1) {rho = 1.046430;}
        else if (l == 2) {rho = 1.280643;}
        break;
      case 53:                 //I
        if (l == 0) {rho = 0.998024;}
        else if (l == 1) {rho = 0.532290;}
        else if (l == 2) {rho = 0.979783;}
        break;
    }
    return rho;
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
      case 9:       //F
        if ((Lbra == 0)&&(Lket == 0)) {eri = 10.4966670;}                                        //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 14.8172560;}                                   //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 16.0736890;}     //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 14.4183930;}   //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 0.7277630;}                                    //(sp|sp)||(ps|ps)
        break;
      case 14:      //Si
        if ((Lbra == 0)&&(Lket == 0)) {eri = 5.0471960;}                                         //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 6.7593670;}                                    //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 5.9490570;}      //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 5.1612970;}    //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 0.9198320;}                                    //(sp|sp)||(ps|ps)
        break;
      case 15:      //P
        if ((Lbra == 0)&&(Lket == 0)) {eri = 7.8016150;}                                         //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 6.6184780;}                                    //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 5.1869490;}      //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 6.0620020;}    //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 1.5428090;}                                    //(sp|sp)||(ps|ps)
        break;
      case 16:      //S
        if ((Lbra == 0)&&(Lket == 0)) {eri = 8.9646670;}                                         //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 9.9681640;}                                    //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 6.7859360;}      //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 7.9702470;}    //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 4.0418360;}                                    //(sp|sp)||(ps|ps)
        break;
      case 17:      //Cl
        if ((Lbra == 0)&&(Lket == 0)) {eri = 16.0136010;}                                        //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 7.5222150;}                                    //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 8.0481150;}      //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 7.5041540;}    //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 3.4811530;}                                    //(sp|sp)||(ps|ps)
        break;
      case 35:      //Br
        if ((Lbra == 0)&&(Lket == 0)) {eri = 15.9434250;}                                        //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 8.2827630;}                                    //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 16.0616800;}     //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 7.8168490;}    //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 0.5788690;}                                    //(sp|sp)||(ps|ps)
        break;
      case 53:      //I
        if ((Lbra == 0)&&(Lket == 0)) {eri = 13.6319430;}                                        //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 7.2883300;}                                    //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 14.9904060;}     //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 5.9664070;}    //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 2.6300350;}                                    //(sp|sp)||(ps|ps)
        break;
    }
    return eri/au2eV;
  }
  double PAPDDG(int atmnr, int index) {
    //function returning the PDDG PAi parameters in Hartree
    double pa = 0.0;
    switch (atmnr) {
      case 1:      //H
        if (index == 1) {pa = 0.057193;}
        else if (index == 2) {pa = -0.034823;}
        break;
      case 6:      //C
        if (index == 1) {pa = -0.000743;}
        else if (index == 2) {pa = 0.000985;}
        break;
      case 7:      //N
        if (index == 1) {pa = -0.003160;}
        else if (index == 2) {pa = 0.012501;}
        break;
      case 8:      //O
        if (index == 1) {pa = -0.001000;}
        else if (index == 2) {pa = -0.001522;}
        break;
      case 9:      //F
        if (index == 1) {pa = -0.012866;}
        else if (index == 2) {pa = 0.007315;}
        break;
      case 14:     //Si
        if (index == 1) {pa = -0.091928;}
        else if (index == 2) {pa = -0.040753;}
        break;
      case 15:     //P
        if (index == 1) {pa = 0.462741;}
        else if (index == 2) {pa = -0.020444;}
        break;
      case 16:     //S
        if (index == 1) {pa = 0.120434;}
        else if (index == 2) {pa = -0.002663;}
        break;
      case 17:     //Cl
        if (index == 1) {pa = -0.016552;}
        else if (index == 2) {pa = -0.016646;}
        break;
      case 35:     //Br
        if (index == 1) {pa = -0.013772;}
        else if (index == 2) {pa = 0.008849;}
        break;
      case 53:     //I
        if (index == 1) {pa = 0.012901;}
        else if (index == 2) {pa = -0.012825;}
        break;
    }
    return pa/au2eV;
  }
  double DAPDDG(int atmnr, int index) {
    //function returning the PDDG DAi parameters in inverse Angstroem
    double da = 0.0;
    switch (atmnr) {
      case 1:      //H
        if (index == 1) {da = 0.663395;}
        else if (index == 2) {da = 1.081901;}
        break;
      case 6:      //C
        if (index == 1) {da = 0.836915;}
        else if (index == 2) {da = 1.585236;}
        break;
      case 7:      //N
        if (index == 1) {da = 1.004172;}
        else if (index == 2) {da = 1.516336;}
        break;
      case 8:      //O
        if (index == 1) {da = 1.360685;}
        else if (index == 2) {da = 1.366407;}
        break;
      case 9:      //F
        if (index == 1) {da = 1.305681;}
        else if (index == 2) {da = 1.842572;}
        break;
      case 14:     //Si
        if (index == 1) {da = 1.163190;}
        else if (index == 2) {da = 2.190526;}
        break;
      case 15:     //P
        if (index == 1) {da = 0.714296;}
        else if (index == 2) {da = 2.041209;}
        break;
      case 16:     //S
        if (index == 1) {da = 0.672870;}
        else if (index == 2) {da = 2.032340;}
        break;
      case 17:     //Cl
        if (index == 1) {da = 1.727690;}
        else if (index == 2) {da = 1.784655;}
        break;
      case 35:     //Br
        if (index == 1) {da = 1.852030;}
        else if (index == 2) {da = 2.338958;}
        break;
      case 53:     //I
        if (index == 1) {da = 1.994299;}
        else if (index == 2) {da = 2.263417;}
        break;
    }
    return da;
  }
  //charge model information; this is separate calculation
  double CM1ck(size_t atomicnr) {
    //function returning the atomic parameters ck for CM1
    double cm1ck = 0.0;
    if ((atomicnr == 1)||((atomicnr > 5)&&(atomicnr < 10))||(atomicnr == 14)||(atomicnr == 16)||(atomicnr == 17)||(atomicnr == 35)||(atomicnr == 53)) {
      if (atomicnr == 9) {cm1ck = 0.3381;}                 //F
      else if (atomicnr == 16) {cm1ck = -0.0834;}          //S
      else if (atomicnr == 17) {cm1ck = -0.1080;}          //Cl
      else if (atomicnr == 35) {cm1ck = -0.0116;}          //Br
      else if (atomicnr == 53) {cm1ck = -0.3213;}          //I
    }
    else {std::cout << "ERROR: PM3PDDG.hpp: PM3PDDG: CM1ck(): partial charge calculation involving non-parametrized atom";}
    return cm1ck;
  }
  double CM1dk(size_t atomicnr) {
    //function returning the atomic parameters dk for CM1
    double cm1dk = 0.0;
    if ((atomicnr == 1)||((atomicnr > 5)&&(atomicnr < 10))||(atomicnr == 14)||(atomicnr == 16)||(atomicnr == 17)||(atomicnr == 35)||(atomicnr == 53)) {
      if (atomicnr == 7) {cm1dk = -0.0909;}                //N
      else if (atomicnr == 8) {cm1dk = -0.0449;}           //O
      else if (atomicnr == 9) {cm1dk = 0.0148;}            //F
      else if (atomicnr == 16) {cm1dk = -0.0848;}          //S
      else if (atomicnr == 17) {cm1dk = -0.1168;}          //Cl
      else if (atomicnr == 35) {cm1dk = -0.0338;}          //Br
      else if (atomicnr == 53) {cm1dk = -0.0636;}          //I
    }
    else {std::cout << "ERROR: PM3PDDG.hpp: PM3PDDG: CM1dk(): partial charge calculation involving non-parametrized atom";}
    return cm1dk;
  }
  double CM1ckkp(size_t atom1, size_t atom2) {
    //function returning the atom-pair parameters ckk' for CM1; implemented the function just in case a refit is ever going to take place
    return 0.0;
  }
  double CM1dkkp(size_t atom1, size_t atom2) {
    //function returning the atom-pair parameters dkk' for CM1
    double cm1dkkp = 0.0;
    bool valid1 = false;
    bool valid2 = false;
    if ((atom1 == 1)||((atom1 > 5)&&(atom1 < 10))||(atom1 == 14)||(atom1 == 16)||(atom1 == 17)||(atom1 == 35)||(atom1 == 53)) {valid1 = true;}
    if ((atom2 == 1)||((atom2 > 5)&&(atom2 < 10))||(atom2 == 14)||(atom2 == 16)||(atom2 == 17)||(atom2 == 35)||(atom2 == 53)) {valid2 = true;}
    if ((valid1)&&(valid2)) {
      size_t atomA = std::min(atom1,atom2);
      size_t atomB = std::max(atom1,atom2);
      if ((atomA == 1)&&(atomB == 7)) {cm1dkkp = 0.1854;}           //H-N
      else if ((atomA == 1)&&(atomB == 8)) {cm1dkkp = 0.1434;}      //H-O
      else if ((atomA == 1)&&(atomB == 14)) {cm1dkkp = -0.1004;}    //H-Si
    }
    else {std::cout << "ERROR: PM3PDDG.hpp: PM3PDDG: CM1dkkp(): partial charge calculation involving non-parametrized atom";}
    return cm1dkkp;
  }
  double CM2ckkp(size_t atom1, size_t atom2) {
    //function returning the atom-pair parameters ckk' for CM2
    double cm2ckkp = 0.0;
    size_t atomA = std::min(atom1,atom2);
    size_t atomB = std::max(atom1,atom2);
    if ((atomA == 1)&&(atomB == 6)) {cm2ckkp = 0.003;}            //H-C
    else if ((atomA == 1)&&(atomB == 7)) {cm2ckkp = 0.274;}       //H-N
    else if ((atomA == 1)&&(atomB == 8)) {cm2ckkp = 0.185;}       //H-O
    else if ((atomA == 1)&&(atomB == 14)) {cm2ckkp = -0.021;}     //H-Si
    else if ((atomA == 1)&&(atomB == 16)) {cm2ckkp = 0.089;}      //H-S
    else if ((atomA == 6)&&(atomB == 7)) {cm2ckkp = -0.022;}      //C-N
    else if ((atomA == 6)&&(atomB == 8)) {cm2ckkp = 0.025;}       //C-O
    else if ((atomA == 6)&&(atomB == 14)) {cm2ckkp = -0.107;}     //C-Si
    else if ((atomA == 6)&&(atomB == 16)) {cm2ckkp = -0.033;}     //C-S
    else if ((atomA == 7)&&(atomB == 8)) {cm2ckkp = -0.030;}      //N-O
    if (atom1 != atomA) {cm2ckkp *= -1.0;}
    return cm2ckkp;
  }
  double CM2dkkp(size_t atom1, size_t atom2) {
    //function returning the atom-pair parameters dkk' for CM2
    double cm2dkkp = 0.0;
    bool valid1 = false;
    bool valid2 = false;
    if ((atom1 == 1)||((atom1 > 5)&&(atom1 < 10))||((atom1 > 13)&&(atom1 < 18))||(atom1 == 35)||(atom1 == 53)) {valid1 = true;}
    if ((atom2 == 1)||((atom2 > 5)&&(atom2 < 10))||((atom2 > 13)&&(atom2 < 18))||(atom2 == 35)||(atom2 == 53)) {valid2 = true;}
    if ((valid1)&&(valid2)) {
      size_t atomA = std::min(atom1,atom2);
      size_t atomB = std::max(atom1,atom2);
      if ((atomA == 1)&&(atomB == 15)) {cm2dkkp = 0.253;}           //H-P
      else if ((atomA == 6)&&(atomB == 7)) {cm2dkkp = 0.156;}       //C-N
      else if ((atomA == 6)&&(atomB == 8)) {cm2dkkp = 0.016;}       //C-O
      else if ((atomA == 6)&&(atomB == 9)) {cm2dkkp = 0.025;}       //C-F
      else if ((atomA == 6)&&(atomB == 15)) {cm2dkkp = 0.082;}      //C-P
      else if ((atomA == 6)&&(atomB == 16)) {cm2dkkp = 0.112;}      //C-S
      else if ((atomA == 6)&&(atomB == 17)) {cm2dkkp = 0.117;}      //C-Cl
      else if ((atomA == 6)&&(atomB == 35)) {cm2dkkp = 0.040;}      //C-Br
      else if ((atomA == 6)&&(atomB == 53)) {cm2dkkp = -0.032;}     //C-I
      else if ((atomA == 7)&&(atomB == 8)) {cm2dkkp = -0.043;}      //N-O
      else if ((atomA == 8)&&(atomB == 15)) {cm2dkkp = 0.181;}      //O-P
      else if ((atomA == 8)&&(atomB == 16)) {cm2dkkp = 0.056;}      //O-S
      else if ((atomA == 9)&&(atomB == 15)) {cm2dkkp = 0.244;}      //F-P
      else if ((atomA == 15)&&(atomB == 16)) {cm2dkkp = -0.087;}    //P-S
      if (atom1 != atomA) {cm2dkkp *= -1.0;}
    }
    else {std::cout << "ERROR: PM3PDDG.hpp: PM3PDDG: CM2dkkp(): partial charge calculation involving non-parametrized atom";}
    return cm2dkkp;
  }
  double CM3ckkp(size_t atom1, size_t atom2) {
    //function returning the atom-pair parameters ckk' for CM3
    double cm3ckkp = 0.0;
    size_t atomA = std::min(atom1,atom2);
    size_t atomB = std::max(atom1,atom2);
    if ((atomA == 6)&&(atomB == 7)) {cm3ckkp = -0.017;}           //C-N
    else if ((atomA == 6)&&(atomB == 8)) {cm3ckkp = 0.010;}       //C-O
    else if ((atomA == 7)&&(atomB == 8)) {cm3ckkp = 0.274;}       //N-O
    else if ((atomA == 8)&&(atomB == 14)) {cm3ckkp = -0.044;}     //O-Si
    else if ((atomA == 8)&&(atomB == 15)) {cm3ckkp = -0.076;}     //O-P
    else if ((atomA == 15)&&(atomB == 16)) {cm3ckkp = 0.245;}     //P-S
    if (atom1 != atomA) {cm3ckkp *= -1.0;}
    return cm3ckkp;
  }
  double CM3dkkp(size_t atom1, size_t atom2) {
    //function returning the atom-pair parameters dkk' for CM3
    double cm3dkkp = 0.0;
    bool valid1 = false;
    bool valid2 = false;
    if ((atom1 == 1)||(atom1 == 3)||((atom1 > 5)&&(atom1 < 10))||((atom1 > 13)&&(atom1 < 18))||(atom1 == 35)) {valid1 = true;}
    if ((atom2 == 1)||(atom2 == 3)||((atom2 > 5)&&(atom2 < 10))||((atom2 > 13)&&(atom2 < 18))||(atom2 == 35)) {valid2 = true;}
    if ((valid1)&&(valid2)) {
      size_t atomA = std::min(atom1,atom2);
      size_t atomB = std::max(atom1,atom2);
      if ((atomA == 1)&&(atomB == 6)) {cm3dkkp = 0.021;}            //H-C
      else if ((atomA == 1)&&(atomB == 7)) {cm3dkkp = 0.243;}       //H-N
      else if ((atomA == 1)&&(atomB == 8)) {cm3dkkp = 0.153;}       //H-O
      else if ((atomA == 1)&&(atomB == 14)) {cm3dkkp = 0.191;}      //H-Si
      else if ((atomA == 1)&&(atomB == 15)) {cm3dkkp = 0.233;}      //H-P
      else if ((atomA == 1)&&(atomB == 16)) {cm3dkkp = 0.118;}      //H-S
      else if ((atomA == 3)&&(atomB == 6)) {cm3dkkp = 0.192;}       //Li-C
      else if ((atomA == 3)&&(atomB == 7)) {cm3dkkp = 0.403;}       //Li-N
      else if ((atomA == 3)&&(atomB == 8)) {cm3dkkp = 0.390;}       //Li-O
      else if ((atomA == 3)&&(atomB == 9)) {cm3dkkp = 0.430;}       //Li-F
      else if ((atomA == 3)&&(atomB == 16)) {cm3dkkp = 0.268;}      //Li-S
      else if ((atomA == 3)&&(atomB == 17)) {cm3dkkp = 0.117;}      //Li-Cl
      else if ((atomA == 6)&&(atomB == 7)) {cm3dkkp = 0.126;}       //C-N
      else if ((atomA == 6)&&(atomB == 8)) {cm3dkkp = 0.018;}       //C-O
      else if ((atomA == 6)&&(atomB == 9)) {cm3dkkp = 0.019;}       //C-F
      else if ((atomA == 6)&&(atomB == 14)) {cm3dkkp = 0.069;}      //C-Si
      else if ((atomA == 6)&&(atomB == 15)) {cm3dkkp = 0.197;}      //C-P
      else if ((atomA == 6)&&(atomB == 16)) {cm3dkkp = 0.076;}      //C-S
      else if ((atomA == 6)&&(atomB == 17)) {cm3dkkp = 0.097;}      //C-Cl
      else if ((atomA == 6)&&(atomB == 35)) {cm3dkkp = 0.028;}      //C-Br
      else if ((atomA == 7)&&(atomB == 8)) {cm3dkkp = -0.085;}      //N-O
      else if ((atomA == 7)&&(atomB == 15)) {cm3dkkp = 0.018;}      //N-P
      else if ((atomA == 8)&&(atomB == 14)) {cm3dkkp = 0.141;}      //O-Si
      else if ((atomA == 8)&&(atomB == 15)) {cm3dkkp = 0.231;}      //O-P
      else if ((atomA == 8)&&(atomB == 16)) {cm3dkkp = 0.125;}      //O-S
      else if ((atomA == 9)&&(atomB == 14)) {cm3dkkp = 0.099;}      //F-Si
      else if ((atomA == 9)&&(atomB == 15)) {cm3dkkp = 0.175;}      //F-P
      else if ((atomA == 14)&&(atomB == 17)) {cm3dkkp = -0.180;}    //Si-Cl
      else if ((atomA == 15)&&(atomB == 16)) {cm3dkkp = -0.534;}    //P-S
      else if ((atomA == 15)&&(atomB == 17)) {cm3dkkp = -0.394;}    //P-Cl
      if (atom1 != atomA) {cm3dkkp *= -1.0;}
    }
    else {std::cout << "WARNING: PM3PDDG.hpp: PM3PDDG: CM3dkkp(): partial charge calculation involving non-parametrized atom";}
    return cm3dkkp;
  }
};

#endif //_PM3PDDG_
