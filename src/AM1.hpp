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

#ifndef _AM1_
#define _AM1_
#include "MNDO.hpp"

//descrition:
//AM1 method

class AM1: public MNDO {
  //this is the implementation of Dewar's AM1
  //M. J. S. Dewar, E. G. Zoebisch, E. F. Healy, J. J. P. Stewart, J. Am. Chem. Soc., 107(13), 3902, 1985
  //for charge and dipole moment model: J. W. Storer, D. J. Giesen, C. J. Cramer, D. G. Truhlar, J. Comp.-Aid. Mol. Des., 9, 87, 1995   -> CM1
  //                                    J. Li, T. Zhu, C. J. Cramer, D. G. Truhlar, J. Phys. Chem. A, 102, 1820, 1998                   -> CM2
  //bX quantities are "barred" tensors, which are used only for the open-shell case
public:
  AM1(BSet _bset, Molecule _mol, std::string _openclosed = "0", std::string _corecorrection = "0"): MNDO(_bset,_mol,_openclosed,_corecorrection) {}
  ~AM1() {}
  std::string Type() {return "AM1";}
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
      else if ((atoms[iatm] > 2)&&(atoms[iatm] < 10)) {def = true;}         //Li,Be,B,C,N,O,F
      else if ((atoms[iatm] > 12)&&(atoms[iatm] < 18)) {def = true;}        //Al,Si,P,S,Cl
      else if ((atoms[iatm] == 30)||(atoms[iatm] == 32)) {def = true;}      //Zn,Ge
      else if ((atoms[iatm] == 33)||(atoms[iatm] == 34)) {def = true;}      //As,Se
      else if ((atoms[iatm] == 35)||(atoms[iatm] == 53)) {def = true;}      //Br,I
      else if ((atoms[iatm] == 51)||(atoms[iatm] == 52)) {def = true;}      //Sb,Te
      else if ((atoms[iatm] == 11)||(atoms[iatm] == 80)) {def = true;}      //Na,Hg
      if (!def) {throw("ERROR: AM1.hpp: AM1: checkAtoms(): atom not fully specified for AM1-theory");}
    }
  }
  double ZeroOverlap(size_t atm) {return ZeroOverlapAM1(atm);}
  double ElementFormTheo(size_t atomicnr) {
    //function returning theoretical heats of formation for elements; values directly in a.u.
    double enth = 0.0;
    switch (atomicnr) {
      case 1:       //H
        enth = -0.418811010146;
        break;
      case 3:       //Li
        enth = -0.188450543318;
        break;
      case 4:       //Be
        enth = -0.889508467058;
        break;
      case 5:       //B
        enth = -2.34156653822;
        break;
      case 6:       //C
        enth = -4.41711440792;
        break;
      case 7:       //N
        enth = -7.41280733472;
        break;
      case 8:       //O
        enth = -11.5881473027;
        break;
      case 9:       //F
        enth = -17.7238538228;
        break;
      case 11:      //Na
        enth = -0.19313741269;
        break;
      case 13:      //Al
        enth = -1.70593365991;
        break;
      case 14:      //Si
        enth = -2.88911227453;
        break;
      case 15:      //P
        enth = -4.56249735046;
        break;
      case 16:      //S
        enth = -7.00454584332;
        break;
      case 17:      //Cl
        enth = -13.6780414477;
        break;
      case 30:      //Zn
        enth = -1.11277017685;
        break;
      case 32:      //Ge
        enth = -2.88509340499;
        break;
      case 33:      //As
        enth = -4.49926063551;
        break;
      case 34:      //Se
        enth = -5.25247215157;
        break;
      case 35:      //Br
        enth = -12.9473096825;
        break;
      case 51:      //Sb
        enth = -4.43488928513;
        break;
      case 52:      //Te
        enth = -5.12393240263;
        break;
      case 53:      //I
        enth = -12.7470286834;
        break;
      case 80:      //Hg
        enth = -1.06878637863;
        break;
    }
    return enth;
  }
  double betaA0(size_t atomicnr, int L) {
    //function returning the bonding parameters; values stored in eV, returned however in a.u.
    double betaa0 = 0.0;
    switch (atomicnr) {
      case 1:      //H
        betaa0 = -6.173787;
        break;
      case 3:      //Li
        betaa0 = -1.35004;
        break;
      case 4:      //Be
        betaa0 = -4.017096;
        break;
      case 5:      //B
        if (L == 0) {betaa0 = -9.599114;}
        else if (L == 1) {betaa0 = -6.273757;}
        break;
      case 6:      //C
        if (L == 0) {betaa0 = -15.715783;}
        else if (L == 1) {betaa0 = -7.719283;}
        break;
      case 7:      //N
        if (L == 0) {betaa0 = -20.299110;}
        else if (L == 1) {betaa0 = -18.238666;}
        break;
      case 8:      //O
        betaa0 = -29.272773;
        break;
      case 9:      //F
        if (L == 0) {betaa0 = -69.590277;}
        else if (L == 1) {betaa0 = -27.92236;}
        break;
      case 11:     //Na
        if (L == 0) {betaa0 = -1.4536944;}
        else if (L == 1) {betaa0 = -0.2298064;}
        break;
      case 13:     //Al
        if (L == 0) {betaa0 = -3.866822;}
        else if (L == 1) {betaa0 = -2.317146;}
        break;
      case 14:     //Si
        if (L == 0) {betaa0 = -3.784852;}
        else if (L == 1) {betaa0 = -1.968123;}
        break;
      case 15:     //P
        if (L == 0) {betaa0 = -6.353764;}
        else if (L == 1) {betaa0 = -6.590709;}
        break;
      case 16:     //S
        if (L == 0) {betaa0 = -3.920566;}
        else if (L == 1) {betaa0 = -7.905278;}
        break;
      case 17:     //Cl
        if (L == 0) {betaa0 = -24.59467;}
        else if (L == 1) {betaa0 = -14.637216;}
        break;
      case 30:     //Zn
        if (L == 0) {betaa0 = -1.997429;}
        else if (L == 1) {betaa0 = -4.758119;}
        break;
      case 32:     //Ge
        if (L == 0) {betaa0 = -4.356607;}
        else if (L == 1) {betaa0 = -0.991091;}
        break;
      case 33:     //As
        if (L == 0) {betaa0 = -5.6481504;}
        else if (L == 1) {betaa0 = -4.9979109;}
        break;
      case 34:     //Se
        if (L == 0) {betaa0 = -3.1470826;}
        else if (L == 1) {betaa0 = -6.1468406;}
        break;
      case 35:     //Br
        if (L == 0) {betaa0 = -19.39988;}
        else if (L == 1) {betaa0 = -8.957195;}
        break;
      case 51:     //Sb
        if (L == 0) {betaa0 = -7.38233;}
        else if (L == 1) {betaa0 = -3.633119;}
        break;
      case 52:     //Te
        if (L == 0) {betaa0 = -8.3897294;}
        else if (L == 1) {betaa0 = -5.1065429;}
        break;
      case 53:     //I
        if (L == 0) {betaa0 = -8.443327;}
        else if (L == 1) {betaa0 = -6.323405;}
        break;
      case 80:     //Hg
        if (L == 0) {betaa0 = -0.908657;}
        else if (L == 1) {betaa0 = -4.909384;}
        break;
    }
    return betaa0/au2eV;
  }
  double UlX(size_t atomicnr, int L) {
    //function that returns the atomic potential U
    double ulx = 0.0;
    switch (atomicnr) {
      case 1:                 //H
        ulx = -11.396427;
        break;
      case 3:                 //Li
        if (L == 0) {ulx = -5.128;}
        else if (L == 1) {ulx = -2.7212;}
        break;
      case 4:                 //Be
        if (L == 0) {ulx = -16.602378;}
        else if (L == 1) {ulx = -10.703771;}
        break;
      case 5:                 //B
        if (L == 0) {ulx = -34.492870;}
        else if (L == 1) {ulx = -22.631525;}
        break;
      case 6:                 //C
        if (L == 0) {ulx = -52.028658;}
        else if (L == 1) {ulx = -39.614239;}
        break;
      case 7:                 //N
        if (L == 0) {ulx = -71.860000;}
        else if (L == 1) {ulx = -57.167581;}
        break;
      case 8:                 //O
        if (L == 0) {ulx = -97.830000;}
        else if (L == 1) {ulx = -78.262380;}
        break;
      case 9:                 //F
        if (L == 0) {ulx = -136.105579;}
        else if (L == 1) {ulx = -104.889885;}
        break;
      case 11:                //Na
        if (L == 0) {ulx = -5.2555362;}
        else if (L == 1) {ulx = -2.0812781;}
        break;
      case 13:                //Al
        if (L == 0) {ulx = -24.353585;}
        else if (L == 1) {ulx = -18.363645;}
        break;
      case 14:                //Si
        if (L == 0) {ulx = -33.953622;}
        else if (L == 1) {ulx = -28.934749;}
        break;
      case 15:                //P
        if (L == 0) {ulx = -42.029863;}
        else if (L == 1) {ulx = -34.030709;}
        break;
      case 16:                //S
        if (L == 0) {ulx = -56.694056;}
        else if (L == 1) {ulx = -48.717049;}
        break;
      case 17:                //Cl
        if (L == 0) {ulx = -111.613948;}
        else if (L == 1) {ulx = -76.640107;}
        break;
      case 30:                //Zn
        if (L == 0) {ulx = -21.040008;}
        else if (L == 1) {ulx = -17.655574;}
        break;
      case 32:                //Ge
        if (L == 0) {ulx = -34.183889;}
        else if (L == 1) {ulx = -28.640811;}
        break;
      case 33:                //As
        if (L == 0) {ulx = -41.681751;}
        else if (L == 1) {ulx = -33.4506152;}
        break;
      case 34:                //Se
        if (L == 0) {ulx = -41.9984056;}
        else if (L == 1) {ulx = -32.8575485;}
        break;
      case 35:                //Br
        if (L == 0) {ulx = -104.656063;}
        else if (L == 1) {ulx = -74.930052;}
        break;
      case 51:                //Sb
        if (L == 0) {ulx = -44.438162;}
        else if (L == 1) {ulx = -32.389514;}
        break;
      case 52:                //Te
        if (L == 0) {ulx = -39.245423;}
        else if (L == 1) {ulx = -30.8515845;}
        break;
      case 53:                //I
        if (L == 0) {ulx = -103.589663;}
        else if (L == 1) {ulx = -74.429997;}
        break;
      case 80:                //Hg
        if (L == 0) {ulx = -19.941578;}
        else if (L == 1) {ulx = -11.11087;}
        break;
    }
    return ulx/au2eV;
  }
  double alpha(int atomicnr, int atm2) {
    //function that returns alpha values for AM1
    //alphas stored and returned in 1/Angstrom
    double alpha = 0.0;
    switch (atomicnr) {
      case 1:      //H
        alpha = 2.882324;
        break;
      case 3:      //Li
        alpha = 1.25014;
        break;
      case 4:      //Be
        alpha = 1.669434;
        break;
      case 5:      //B
        alpha = 2.446909;
        break;
      case 6:      //C
        alpha = 2.648274;
        break;
      case 7:      //N
        alpha = 2.947286;
        break;
      case 8:      //O
        alpha = 4.455371;
        break;
      case 9:      //F
        alpha = 5.5178;
        break;
      case 11:     //Na
        alpha = 2.2487164;
        break;
      case 13:     //Al
        alpha = 1.976586;
        break;
      case 14:     //Si
        alpha = 2.257816;
        break;
      case 15:     //P
        alpha = 2.455322;
        break;
      case 16:     //S
        alpha = 2.461648;
        break;
      case 17:     //Cl
        alpha = 2.919368;
        break;
      case 30:     //Zn
        alpha = 1.484563;
        break;
      case 32:     //Ge
        alpha = 2.136405;
        break;
      case 33:     //As
        alpha = 2.240538;
        break;
      case 34:     //Se
        alpha = 2.6375694;
        break;
      case 35:     //Br
        alpha = 2.576546;
        break;
      case 51:     //Sb
        alpha = 2.276331;
        break;
      case 52:     //Te          this value is totally unphysical!!!
        alpha = 6.0171167;
        break;
      case 53:     //I
        alpha = 2.299424;
        break;
      case 80:     //Hg
        alpha = 1.484734;
        break;
    }
    return alpha;
  }
  std::vector<double> AM1K(int atomicnr, int atm2 = 0) {
    std::vector<double> am1k;
    switch (atomicnr) {
      case 1:     //H
        am1k.push_back(0.122796/au2eV);
        am1k.push_back(0.005090/au2eV);
        am1k.push_back(-0.018336/au2eV);
        break;
      case 3:     //Li
        am1k.push_back(0.0);
        break;
      case 4:     //Be
        am1k.push_back(0.0);
        break;
      case 5:     //B, boron has different interactions with different atoms, as a unique set of parameters could not be fitted
        if (atm2 == 1) {                                                       //H
          am1k.push_back(0.412253/au2eV);
          am1k.push_back(-0.149917/au2eV);
        }
        else if (atm2 == 6) {                                                  //C
          am1k.push_back(0.261751/au2eV);
          am1k.push_back(0.050275/au2eV);
        }
        else if ((atm2 == 9)||(atm2 == 17)||(atm2 == 35)||(atm2 == 53)) {      //Halogen
          am1k.push_back(0.359244/au2eV);
          am1k.push_back(0.074729/au2eV);
        }
        else {                                                                 //anything else
          am1k.push_back(0.182613/au2eV);
          am1k.push_back(0.118587/au2eV);
          am1k.push_back(-0.073280/au2eV);
        }
        break;
      case 6:     //C
        am1k.push_back(0.011355/au2eV);
        am1k.push_back(0.045924/au2eV);
        am1k.push_back(-0.020061/au2eV);
        am1k.push_back(-0.001260/au2eV);
        break;
      case 7:     //N
        am1k.push_back(0.025251/au2eV);
        am1k.push_back(0.028953/au2eV);
        am1k.push_back(-0.005806/au2eV);
        break;
      case 8:     //O
        am1k.push_back(0.280962/au2eV);
        am1k.push_back(0.081430/au2eV);
        break;
      case 9:     //F
        am1k.push_back(0.242079/au2eV);
        am1k.push_back(0.003607/au2eV);
        break;
      case 11:    //Na
        am1k.push_back(0.5322668/au2eV);
        am1k.push_back(0.9223598/au2eV);
        break;
      case 13:    //Al
        am1k.push_back(0.09/au2eV);
        break;
      case 14:    //Si
        am1k.push_back(0.25/au2eV);
        am1k.push_back(0.061513/au2eV);
        am1k.push_back(0.020789/au2eV);
        break;
      case 15:    //P
        am1k.push_back(-0.031827/au2eV);
        am1k.push_back(0.01847/au2eV);
        am1k.push_back(0.03329/au2eV);
        break;
      case 16:    //S
        am1k.push_back(-0.509195/au2eV);
        am1k.push_back(-0.011863/au2eV);
        am1k.push_back(0.012334/au2eV);
        break;
      case 17:    //Cl
        am1k.push_back(0.094243/au2eV);
        am1k.push_back(0.027168/au2eV);
        break;
      case 30:    //Zn
        am1k.push_back(0.0);
        break;
      case 32:    //Ge
        am1k.push_back(0.0);
        break;
      case 33:    //As
        am1k.push_back(-0.0073614/au2eV);
        am1k.push_back(0.0437629/au2eV);
        break;
      case 34:    //Se
        am1k.push_back(0.1116681/au2eV);
        am1k.push_back(0.0396143/au2eV);
        break;
      case 35:    //Br
        am1k.push_back(0.066685/au2eV);
        am1k.push_back(0.025568/au2eV);
        break;
      case 51:    //Sb
        am1k.push_back(-0.596447/au2eV);
        am1k.push_back(0.895513/au2eV);
        break;
      case 52:    //Te
        am1k.push_back(0.4873378/au2eV);
        am1k.push_back(0.1520464/au2eV);
        break;
      case 53:    //I
        am1k.push_back(0.004361/au2eV);
        am1k.push_back(0.015706/au2eV);
        break;
      case 80:    //Ge
        am1k.push_back(0.0);
        break;
    }
    return am1k;
  }
  std::vector<double> AM1L(int atomicnr, int atm2 = 0) {
    std::vector<double> am1l;
    switch (atomicnr) {
      case 1:     //H
        am1l.push_back(5.0);
        am1l.push_back(5.0);
        am1l.push_back(2.0);
        break;
      case 3:     //Li
        am1l.push_back(0.0);
        break;
      case 4:     //Be
        am1l.push_back(0.0);
        break;
      case 5:     //B
        if (atm2 == 1) {
          am1l.push_back(10.0);
          am1l.push_back(6.0);
        }
        else if (atm2 == 6) {
          am1l.push_back(8.0);
          am1l.push_back(5.0);
        }
        else if ((atm2 == 9)||(atm2 == 17)||(atm2 == 35)||(atm2 == 53)) {
          am1l.push_back(9.0);
          am1l.push_back(9.0);
        }
        else {
          am1l.push_back(6.0);
          am1l.push_back(6.0);
          am1l.push_back(5.0);
        }
        break;
      case 6:     //C
        am1l.push_back(5.0);
        am1l.push_back(5.0);
        am1l.push_back(5.0);
        am1l.push_back(5.0);
        break;
      case 7:     //N
        am1l.push_back(5.0);
        am1l.push_back(5.0);
        am1l.push_back(2.0);
        break;
      case 8:     //O
        am1l.push_back(5.0);
        am1l.push_back(7.0);
        break;
      case 9:     //F
        am1l.push_back(4.8);
        am1l.push_back(4.6);
        break;
      case 11:    //Na
        am1l.push_back(0.4800304);
        am1l.push_back(1.9076776);
        break;
      case 13:    //Al
        am1l.push_back(12.392443);
        break;
      case 14:    //Si
        am1l.push_back(9.0);
        am1l.push_back(5.0);
        am1l.push_back(5.0);
        break;
      case 15:    //P
        am1l.push_back(6.0);
        am1l.push_back(7.0);
        am1l.push_back(9.0);
        break;
      case 16:    //S
        am1l.push_back(4.593691);
        am1l.push_back(5.865731);
        am1l.push_back(13.557336);
        break;
      case 17:    //Cl
        am1l.push_back(4.0);
        am1l.push_back(4.0);
        break;
      case 30:    //Zn
        am1l.push_back(0.0);
        break;
      case 32:    //Ge
        am1l.push_back(0.0);
        break;
      case 33:    //As
        am1l.push_back(4.9433993);
        am1l.push_back(3.1944613);
        break;
      case 34:    //Se
        am1l.push_back(6.5086644);
        am1l.push_back(6.5241228);
        break;
      case 35:    //Br
        am1l.push_back(4.0);
        am1l.push_back(4.0);
        break;
      case 51:    //Sb
        am1l.push_back(6.02795);
        am1l.push_back(3.028109);
        break;
      case 52:    //Te
        am1l.push_back(6.0519413);
        am1l.push_back(3.8304067);
        break;
      case 53:    //I
        am1l.push_back(2.3);
        am1l.push_back(3.0);
        break;
    }
    return am1l;
  }
  std::vector<double> AM1M(int atomicnr, int atm2 = 0) {
    std::vector<double> am1m;
    switch (atomicnr) {
      case 1:     //H
        am1m.push_back(1.2);
        am1m.push_back(1.8);
        am1m.push_back(2.1);
        break;
      case 3:     //Li
        am1m.push_back(0.0);
        break;
      case 4:     //Be
        am1m.push_back(0.0);
        break;
      case 5:     //B
        if (atm2 == 1) {
          am1m.push_back(0.832586);
          am1m.push_back(1.186220);
        }
        else if (atm2 == 6) {
          am1m.push_back(1.063995);
          am1m.push_back(1.936492);
        }
        else if ((atm2 == 9)||(atm2 == 17)||(atm2 == 35)||(atm2 == 53)) {
          am1m.push_back(0.819351);
          am1m.push_back(1.574414);
        }
        else {
          am1m.push_back(0.727592);
          am1m.push_back(1.466639);
          am1m.push_back(1.570975);
        }
        break;
      case 6:     //C
        am1m.push_back(1.60);
        am1m.push_back(1.85);
        am1m.push_back(2.05);
        am1m.push_back(2.65);
        break;
      case 7:     //N
        am1m.push_back(1.5);
        am1m.push_back(2.1);
        am1m.push_back(2.4);
        break;
      case 8:     //O
        am1m.push_back(0.847918);
        am1m.push_back(1.445071);
        break;
      case 9:     //F
        am1m.push_back(0.93);
        am1m.push_back(1.66);
        break;
      case 11:    //Na
        am1m.push_back(1.1681055);
        am1m.push_back(1.1537670);
        break;
      case 13:    //Al
        am1m.push_back(2.050394);
        break;
      case 14:    //Si
        am1m.push_back(0.911453);
        am1m.push_back(1.995569);
        am1m.push_back(2.990610);
        break;
      case 15:    //P
        am1m.push_back(1.474323);
        am1m.push_back(1.779354);
        am1m.push_back(3.006576);
        break;
      case 16:    //S
        am1m.push_back(0.770665);
        am1m.push_back(1.503313);
        am1m.push_back(2.009173);
        break;
      case 17:    //Cl
        am1m.push_back(1.3);
        am1m.push_back(2.1);
        break;
      case 30:    //Zn
        am1m.push_back(0.0);
        break;
      case 32:    //Ge
        am1m.push_back(0.0);
        break;
      case 33:    //As
        am1m.push_back(1.4544264);
        am1m.push_back(2.0144939);
        break;
      case 34:    //Se
        am1m.push_back(1.4981077);
        am1m.push_back(2.0751916);
        break;
      case 35:    //Br
        am1m.push_back(1.5);
        am1m.push_back(2.3);
        break;
      case 51:    //Sb
        am1m.push_back(1.710367);
        am1m.push_back(1.538318);
        break;
      case 52:    //Te
        am1m.push_back(1.3079857);
        am1m.push_back(2.0899707);
        break;
      case 53:    //I
        am1m.push_back(1.8);
        am1m.push_back(2.24);
        break;
    }
    return am1m;
  }
  double Dvalue(size_t atmnr, size_t idx) {
    //function returning the D values needed to calculate eris; values stored in Angstrom, but returned in a.u.
    //note that this function returns both D1 and D2; idx is then either 1 or 2
    double dval = 0.0;
    switch (atmnr) {
      case 1:                    //H
        dval = 0.0;
        break;
      case 3:                  //Li
        if (idx == 1) {dval = 1.0874267;}
        else if (idx == 2) {dval = 0.9227121;}
        break;
      case 4:                  //Be
        if (idx == 1) {dval = 0.7605847;}
        else if (idx == 2) {dval = 0.6453775;}
        break;
      case 5:                  //B
        if (idx == 1) {dval = 0.4819453;}
        else if (idx == 2) {dval = 0.4166779;}
        break;
      case 6:                  //C
        if (idx == 1) {dval = 0.4358609;}
        else if (idx == 2) {dval = 0.3845994;}
        break;
      case 7:                  //N
        if (idx == 1) {dval = 0.3404262;}
        else if (idx == 2) {dval = 0.3003302;}
        break;
      case 8:                  //O
        if (idx == 1) {dval = 0.2639959;}
        else if (idx == 2) {dval = 0.2567689;}
        break;
      case 9:                  //F
        if (idx == 1) {dval = 0.2193505;}
        else if (idx == 2) {dval = 0.2597917;}
        break;
      case 11:                 //Na
        if (idx == 1) {dval = 0.84137879635;}
        else if (idx == 2) {dval = 0.72756674912;}
        break;
      case 13:                 //Al
        if (idx == 1) {dval = 0.7429739;}
        else if (idx == 2) {dval = 0.6778182;}
        break;
      case 14:                 //Si
        if (idx == 1) {dval = 0.6154798;}
        else if (idx == 2) {dval = 0.6891036;}
        break;
      case 15:                 //P
        if (idx == 1) {dval = 0.5530865;}
        else if (idx == 2) {dval = 0.4722106;}
        break;
      case 16:                 //S
        if (idx == 1) {dval = 0.4764760;}
        else if (idx == 2) {dval = 0.5310894;}
        break;
      case 17:                 //Cl
        if (idx == 1) {dval = 0.2860828;}
        else if (idx == 2) {dval = 0.4263609;}
        break;
      case 30:                 //Zn
        if (idx == 1) {dval = 0.7186677;}
        else if (idx == 2) {dval = 0.8179549;}
        break;
      case 32:                 //Ge
        if (idx == 1) {dval = 0.6599821;}
        else if (idx == 2) {dval = 0.5661368;}
        break;
      case 33:                 //As
        if (idx == 1) {dval = 0.6365298;}
        else if (idx == 2) {dval = 0.6507545;}
        break;
      case 34:                 //Se
        if (idx == 1) {dval = 0.5353874;}
        else if (idx == 2) {dval = 0.5474124;}
        break;
      case 35:                 //Br
        if (idx == 1) {dval = 0.4475750;}
        else if (idx == 2) {dval = 0.5507111;}
        break;
      case 51:                 //Sb
        if (idx == 1) {dval = 0.7509809;}
        else if (idx == 2) {dval = 0.6127556;}
        break;
      case 52:                 //Te
        if (idx == 1) {dval = 0.8121000;}
        else if (idx == 2) {dval = 0.6896346;}
        break;
      case 53:                 //I
        if (idx == 1) {dval = 0.7873358;}
        else if (idx == 2) {dval = 0.6290413;}
        break;
      case 80:                 //Hg
        if (idx == 1) {dval = 0.9922320;}
        else if (idx == 2) {dval = 0.8161999;}
        break;
    }
    return dval*dist_Angstrom2aum1;
  }
  double rho(size_t atmnr, size_t l) {
    //function returning the rho values needed to calculate eris; values stored in Angstrom but returned in atomic units
    double rho = 0.0;
    switch (atmnr) {
      case 1:                  //H
        if (l == 0) {rho = 0.560345403;}
        else if (l == 1) {rho = 0.0;}
        else if (l == 2) {rho = 0.0;}
        break;
      case 3:                  //Li
        if (l == 0) {rho = 0.986207884;}
        else if (l == 1) {rho = 1.165672376;}
        else if (l == 2) {rho = 1.011953732;}
        break;
      case 4:                  //Be
        if (l == 0) {rho = 0.799924235;}
        else if (l == 1) {rho = 0.788356095;}
        else if (l == 2) {rho = 0.687877905;}
        break;
      case 5:                  //B
        if (l == 0) {rho = 0.679822279;}
        else if (l == 1) {rho = 0.524431177;}
        else if (l == 2) {rho = 0.465909859;}
        break;
      case 6:                  //C
        if (l == 0) {rho = 0.588660438;}
        else if (l == 1) {rho = 0.434959475;}
        else if (l == 2) {rho = 0.411899789;}
        break;
      case 7:                  //N
        if (l == 0) {rho = 0.529751104;}
        else if (l == 1) {rho = 0.338305732;}
        else if (l == 2) {rho = 0.335616880;}
        break;
      case 8:                  //O
        if (l == 0) {rho = 0.466881794;}
        else if (l == 1) {rho = 0.265617656;}
        else if (l == 2) {rho = 0.291866510;}
        break;
      case 9:                  //F
        if (l == 0) {rho = 0.425491557;}
        else if (l == 1) {rho = 0.218866782;}
        else if (l == 2) {rho = 0.280001651;}
        break;
      case 11:                 //Na
        if (l == 0) {rho = 0.9800619701;}
        else if (l == 1) {rho = 1.63523435444;}
        else if (l == 2) {rho = 1.41814447768;}
        break;
      case 13:                 //Al
        if (l == 0) {rho = 0.889903107;}
        else if (l == 1) {rho = 1.005933324;}
        else if (l == 2) {rho = 0.771868341;}
        break;
      case 14:                 //Si
        if (l == 0) {rho = 0.733128067;}
        else if (l == 1) {rho = 0.690852269;}
        else if (l == 2) {rho = 0.712758472;}
        break;
      case 15:                 //P
        if (l == 0) {rho = 0.622778008;}
        else if (l == 1) {rho = 0.807809865;}
        else if (l == 2) {rho = 0.603128119;}
        break;
      case 16:                 //S
        if (l == 0) {rho = 0.610819239;}
        else if (l == 1) {rho = 0.447906465;}
        else if (l == 2) {rho = 0.409892853;}
        break;
      case 17:                 //Cl
        if (l == 0) {rho = 0.478996434;}
        else if (l == 1) {rho = 0.343918655;}
        else if (l == 2) {rho = 0.431383633;}
        break;
      case 30:                 //Zn
        if (l == 0) {rho = 0.610111605;}
        else if (l == 1) {rho = 1.141714309;}
        else if (l == 2) {rho = 1.009411845;}
        break;
      case 32:                 //Ge
        if (l == 0) {rho = 0.707994522;}
        else if (l == 1) {rho = 0.831942745;}
        else if (l == 2) {rho = 0.759073299;}
        break;
      case 33:                 //As
        if (l == 0) {rho = 0.64880765;}
        else if (l == 1) {rho = 0.95995755;}
        else if (l == 2) {rho = 0.927371034;}
        break;
      case 34:                 //Se
        if (l == 0) {rho = 1.060143453;}
        else if (l == 1) {rho = 0.358473620;}
        else if (l == 2) {rho = 0.527366991;}
        break;
      case 35:                 //Br
        if (l == 0) {rho = 0.478791351;}
        else if (l == 1) {rho = 0.439172041;}
        else if (l == 2) {rho = 0.498503548;}
        break;
      case 51:                 //Sb
        if (l == 0) {rho = 0.629847642;}
        else if (l == 1) {rho = 1.086004527;}
        else if (l == 2) {rho = 0.974630761;}
        break;
      case 52:                 //Te
        if (l == 0) {rho = 1.442019673;}
        else if (l == 1) {rho = 0.457437268;}
        else if (l == 2) {rho = 0.551264090;}
        break;
      case 53:                 //I
        if (l == 0) {rho = 0.478663761;}
        else if (l == 1) {rho = 0.588287153;}
        else if (l == 2) {rho = 0.571235649;}
        break;
      case 80:                 //Hg
        if (l == 0) {rho = 0.666603429;}
        else if (l == 1) {rho = 0.904062899;}
        else if (l == 2) {rho = 0.787310536;}
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
      case 1:         //H
        eri = 12.848;                                                                            //(ss|ss)
        break;
      case 3:         //Li
        if ((Lbra == 0)&&(Lket == 0)) {eri = 7.3;}                                               //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 5.0;}                                          //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 5.42;}           //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 4.52;}         //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 0.83;}                                         //(sp|sp)||(ps|ps)
        break;
      case 4:         //Be
        if ((Lbra == 0)&&(Lket == 0)) {eri = 9.0;}                                               //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 6.97;}                                         //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 7.43;}           //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 6.22;}         //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 1.28;}                                         //(sp|sp)||(ps|ps)
        break;
      case 5:         //B
        if ((Lbra == 0)&&(Lket == 0)) {eri = 10.59;}                                             //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 8.86;}                                         //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 9.56;}           //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 7.86;}         //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 1.81;}                                         //(sp|sp)||(ps|ps)
        break;
      case 6:         //C
        if ((Lbra == 0)&&(Lket == 0)) {eri = 12.23;}                                             //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 11.08;}                                        //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 11.47;}          //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 9.84;}         //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 2.43;}                                         //(sp|sp)||(ps|ps)
        break;
      case 7:         //N
        if ((Lbra == 0)&&(Lket == 0)) {eri = 13.59;}                                             //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 12.98;}                                        //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 12.66;}          //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 11.59;}        //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 3.14;}                                         //(sp|sp)||(ps|ps)
        break;
      case 8:         //O
        if ((Lbra == 0)&&(Lket == 0)) {eri = 15.42;}                                             //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 14.52;}                                        //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 14.48;}          //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 12.98;}        //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 3.94;}                                         //(sp|sp)||(ps|ps)
        break;
      case 9:         //F
        if ((Lbra == 0)&&(Lket == 0)) {eri = 16.92;}                                             //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 16.71;}                                        //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 17.25;}          //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 14.91;}        //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 4.83;}                                         //(sp|sp)||(ps|ps)
        break;
      case 11:        //Na
        if ((Lbra == 0)&&(Lket == 0)) {eri = 7.3459178;}                                         //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 4.1130516;}                                    //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 8.4042550;}      //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 4.0370957;}    //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 0.2487699;}                                    //(sp|sp)||(ps|ps)
        break;
      case 13:        //Al
        if ((Lbra == 0)&&(Lket == 0)) {eri = 8.09;}                                              //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 5.98;}                                         //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 6.63;}           //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 5.40;}         //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 0.70;}                                         //(sp|sp)||(ps|ps)
        break;
      case 14:        //Si
        if ((Lbra == 0)&&(Lket == 0)) {eri = 9.82;}                                              //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 7.31;}                                         //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 8.36;}           //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 6.54;}         //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 1.32;}                                         //(sp|sp)||(ps|ps)
        break;
      case 15:        //P
        if ((Lbra == 0)&&(Lket == 0)) {eri = 11.560005;}                                         //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 7.877589;}                                     //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 5.237449;}       //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 7.307648;}     //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 0.779238;}                                     //(sp|sp)||(ps|ps)
        break;
      case 16:        //S
        if ((Lbra == 0)&&(Lket == 0)) {eri = 11.786329;}                                         //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 10.039308;}                                    //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 8.663127;}       //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 7.781688;}     //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 2.532137;}                                     //(sp|sp)||(ps|ps)
        break;
      case 17:        //Cl
        if ((Lbra == 0)&&(Lket == 0)) {eri = 15.03;}                                             //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 11.30;}                                        //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 13.16;}          //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 9.97;}         //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 2.42;}                                         //(sp|sp)||(ps|ps)
        break;
      case 30:        //Zn
        if ((Lbra == 0)&&(Lket == 0)) {eri = 11.8;}                                              //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 13.3;}                                         //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 11.182018;}      //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 12.93052;}     //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 0.484606;}                                     //(sp|sp)||(ps|ps)
        break;
      case 32:        //Ge
        if ((Lbra == 0)&&(Lket == 0)) {eri = 10.168605;}                                       //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 6.671902;}                                   //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 8.144473;}     //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 6.269706;}   //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 0.937093;}                                   //(sp|sp)||(ps|ps)
        break;
      case 33:        //As
        if ((Lbra == 0)&&(Lket == 0)) {eri = 11.0962258;}                                        //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 7.8781648;}                                    //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 4.9259328;}      //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 7.5961088;}    //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 0.6246173;}                                    //(sp|sp)||(ps|ps)
        break;
      case 34:        //Se
        if ((Lbra == 0)&&(Lket == 0)) {eri = 6.7908891;}                                         //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 6.4769273;}                                    //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 6.4812786;}      //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 5.2796993;}    //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 4.4548356;}                                    //(sp|sp)||(ps|ps)
        break;
      case 35:        //Br
        if ((Lbra == 0)&&(Lket == 0)) {eri = 15.03643948;}                                       //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 11.27632539;}                                  //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 13.03468242;}    //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 9.85442552;}   //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 2.45586832;}                                   //(sp|sp)||(ps|ps)
        break;
      case 51:        //Sb
        if ((Lbra == 0)&&(Lket == 0)) {eri = 11.430251;}                                         //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 6.424094;}                                     //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 5.787922;}       //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 6.849181;}     //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 0.58834;}                                      //(sp|sp)||(ps|ps)
        break;
      case 52:        //Te
        if ((Lbra == 0)&&(Lket == 0)) {eri = 4.9925231;}                                         //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 7.2097852;}                                    //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 4.9721484;}      //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 5.6211521;}    //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 4.0071821;}                                    //(sp|sp)||(ps|ps)
        break;
      case 53:        //I
        if ((Lbra == 0)&&(Lket == 0)) {eri = 15.04044855;}                                       //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 11.14778369;}                                  //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 13.05655798;}    //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 9.91409071;}   //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 2.45638202;}                                   //(sp|sp)||(ps|ps)
        break;
      case 80:        //Hg
        if ((Lbra == 0)&&(Lket == 0)) {eri = 10.8;}                                              //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 14.3;}                                         //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 9.3;}            //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 13.5;}         //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 1.3;}                                          //(sp|sp)||(ps|ps)
        break;
    }
    return eri/au2eV;
  }
  //charge model information; this is separate calculation
  double CM1ck(size_t atomicnr) {
    //function returning the atomic parameters ck for CM1
    double cm1ck = 0.0;
    if ((atomicnr == 1)||((atomicnr > 5)&&(atomicnr < 10))||(atomicnr == 14)||(atomicnr == 16)||(atomicnr == 17)||(atomicnr == 35)||(atomicnr == 53)) {
      if (atomicnr == 7) {cm1ck = 0.3846;}                 //N
      else if (atomicnr == 9) {cm1ck = 0.1468;}            //F
      else if (atomicnr == 16) {cm1ck = -0.1311;}          //S
      else if (atomicnr == 17) {cm1ck = 0.0405;}           //Cl
      else if (atomicnr == 35) {cm1ck = 0.1761;}           //Br
      else if (atomicnr == 53) {cm1ck = 0.2380;}           //I
    }
    else {std::cout << "WARNING: AM1.hpp: AM1: CM1ck(): partial charge calculation involving non-parametrized atom";}
    return cm1ck;
  }
  double CM1dk(size_t atomicnr) {
    //function returning the atomic parameters dk for CM1
    double cm1dk = 0.0;
    if ((atomicnr == 1)||((atomicnr > 5)&&(atomicnr < 10))||(atomicnr == 14)||(atomicnr == 16)||(atomicnr == 17)||(atomicnr == 35)||(atomicnr == 53)) {
      if (atomicnr == 8) {cm1dk = -0.0283;}                //O
      else if (atomicnr == 9) {cm1dk = 0.0399;}            //F
      else if (atomicnr == 16) {cm1dk = -0.0956;}          //S
      else if (atomicnr == 17) {cm1dk = -0.0276;}          //Cl
      else if (atomicnr == 35) {cm1dk = -0.0802;}          //Br
      else if (atomicnr == 53) {cm1dk = -0.1819;}          //I
    }
    else {std::cout << "WARNING: AM1.hpp: AM1: CM1dk(): partial charge calculation involving non-parametrized atom";}
    return cm1dk;
  }
  double CM1ckkp(size_t atom1, size_t atom2) {
    //function returning the atom-pair parameters ckk' for CM1
    size_t atomA = std::min(atom1,atom2);
    size_t atomB = std::max(atom1,atom2);
    double cm1ckkp = 0.0;
    if ((atomA == 6)&&(atomB == 7)) {cm1ckkp = 0.3846;}           //C-N
    return cm1ckkp;
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
      if ((atomA == 1)&&(atomB == 7)) {cm1dkkp = 0.0850;}           //H-N
      else if ((atomA == 1)&&(atomB == 8)) {cm1dkkp = 0.1447;}      //H-O
      else if ((atomA == 1)&&(atomB == 14)) {cm1dkkp = 0.0640;}     //H-Si
      else if ((atomA == 6)&&(atomB == 7)) {cm1dkkp = -0.0880;}     //C-N
      else if ((atomA == 7)&&(atomB == 8)) {cm1dkkp = -0.0630;}     //N-O
      else if ((atomA == 8)&&(atomB == 16)) {cm1dkkp = -0.0600;}    //O-S
    }
    else {std::cout << "WARNING: AM1.hpp: AM1: CM1dkkp(): partial charge calculation involving non-parametrized atom";}
    return cm1dkkp;
  }
  double CM2ckkp(size_t atom1, size_t atom2) {
    //function returning the atom-pair parameters ckk' for CM2
    double cm2ckkp = 0.0;
    size_t atomA = std::min(atom1,atom2);
    size_t atomB = std::max(atom1,atom2);
    if ((atomA == 1)&&(atomB == 6)) {cm2ckkp = -0.020;}           //H-C
    else if ((atomA == 1)&&(atomB == 7)) {cm2ckkp = 0.207;}       //H-N
    else if ((atomA == 1)&&(atomB == 8)) {cm2ckkp = 0.177;}       //H-O
    else if ((atomA == 1)&&(atomB == 14)) {cm2ckkp = -0.083;}     //H-Si
    else if ((atomA == 1)&&(atomB == 16)) {cm2ckkp = 0.038;}      //H-S
    else if ((atomA == 6)&&(atomB == 7)) {cm2ckkp = 0.008;}       //C-N
    else if ((atomA == 6)&&(atomB == 8)) {cm2ckkp = 0.026;}       //C-O
    else if ((atomA == 6)&&(atomB == 14)) {cm2ckkp = 0.062;}      //C-Si
    else if ((atomA == 6)&&(atomB == 16)) {cm2ckkp = -0.059;}     //C-S
    else if ((atomA == 7)&&(atomB == 8)) {cm2ckkp = -0.197;}      //N-O
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
      if ((atomA == 1)&&(atomB == 15)) {cm2dkkp = 0.103;}           //H-P
      else if ((atomA == 6)&&(atomB == 7)) {cm2dkkp = 0.086;}       //C-N
      else if ((atomA == 6)&&(atomB == 8)) {cm2dkkp = 0.016;}       //C-O
      else if ((atomA == 6)&&(atomB == 9)) {cm2dkkp = 0.019;}       //C-F
      else if ((atomA == 6)&&(atomB == 15)) {cm2dkkp = -0.019;}     //C-P
      else if ((atomA == 6)&&(atomB == 16)) {cm2dkkp = 0.171;}      //C-S
      else if ((atomA == 6)&&(atomB == 17)) {cm2dkkp = 0.027;}      //C-Cl
      else if ((atomA == 6)&&(atomB == 35)) {cm2dkkp = 0.081;}      //C-Br
      else if ((atomA == 6)&&(atomB == 53)) {cm2dkkp = 0.147;}      //C-I
      else if ((atomA == 7)&&(atomB == 8)) {cm2dkkp = 0.134;}       //N-O
      else if ((atomA == 8)&&(atomB == 15)) {cm2dkkp = 0.088;}      //O-P
      else if ((atomA == 9)&&(atomB == 15)) {cm2dkkp = 0.252;}      //F-P
      else if ((atomA == 15)&&(atomB == 16)) {cm2dkkp = -0.080;}    //P-S
      if (atom1 != atomA) {cm2dkkp *= -1.0;}
    }
    else {std::cout << "WARNING: AM1.hpp: AM1: CM2dkkp(): partial charge calculation involving non-parametrized atom";}
    return cm2dkkp;
  }
  double CM3ckkp(size_t atom1, size_t atom2) {
    //function returning the atom-pair parameters ckk' for CM3
    double cm3ckkp = 0.0;
    size_t atomA = std::min(atom1,atom2);
    size_t atomB = std::max(atom1,atom2);
    if ((atomA == 6)&&(atomB == 8)) {cm3ckkp = 0.029;}            //C-O
    else if ((atomA == 7)&&(atomB == 8)) {cm3ckkp = 0.439;}       //N-O
    else if ((atomA == 8)&&(atomB == 14)) {cm3ckkp = 0.201;}      //O-Si
    else if ((atomA == 8)&&(atomB == 15)) {cm3ckkp = 0.004;}      //O-P
    else if ((atomA == 15)&&(atomB == 16)) {cm3ckkp = 0.341;}     //P-S
    if (atom1 != atomA) {cm3ckkp *= -1.0;}
    return cm3ckkp;
  }
  double CM3dkkp(size_t atom1, size_t atom2) {
    //function returning the atom-pair parameters dkk' for CM3
    double cm3dkkp = 0.0;
    bool valid1 = false;
    bool valid2 = false;
    if ((atom1 == 1)||((atom1 > 5)&&(atom1 < 10))||((atom1 > 13)&&(atom1 < 18))||(atom1 == 35)) {valid1 = true;}
    if ((atom2 == 1)||((atom2 > 5)&&(atom2 < 10))||((atom2 > 13)&&(atom2 < 18))||(atom2 == 35)) {valid2 = true;}
    if ((valid1)&&(valid2)) {
      size_t atomA = std::min(atom1,atom2);
      size_t atomB = std::max(atom1,atom2);
      if ((atomA == 1)&&(atomB == 6)) {cm3dkkp = -0.009;}           //H-C
      else if ((atomA == 1)&&(atomB == 7)) {cm3dkkp = 0.200;}       //H-N
      else if ((atomA == 1)&&(atomB == 8)) {cm3dkkp = 0.154;}       //H-O
      else if ((atomA == 1)&&(atomB == 14)) {cm3dkkp = 0.231;}      //H-Si
      else if ((atomA == 1)&&(atomB == 15)) {cm3dkkp = 0.078;}      //H-P
      else if ((atomA == 1)&&(atomB == 16)) {cm3dkkp = 0.072;}      //H-S
      else if ((atomA == 6)&&(atomB == 7)) {cm3dkkp = 0.091;}       //C-N
      else if ((atomA == 6)&&(atomB == 8)) {cm3dkkp = -0.012;}      //C-O
      else if ((atomA == 6)&&(atomB == 9)) {cm3dkkp = 0.006;}       //C-F
      else if ((atomA == 6)&&(atomB == 14)) {cm3dkkp = 0.163;}      //C-Si
      else if ((atomA == 6)&&(atomB == 15)) {cm3dkkp = 0.106;}      //C-P
      else if ((atomA == 6)&&(atomB == 16)) {cm3dkkp = 0.103;}      //C-S
      else if ((atomA == 6)&&(atomB == 17)) {cm3dkkp = 0.013;}      //C-Cl
      else if ((atomA == 6)&&(atomB == 35)) {cm3dkkp = 0.073;}      //C-Br
      else if ((atomA == 7)&&(atomB == 8)) {cm3dkkp = -0.132;}      //N-O
      else if ((atomA == 7)&&(atomB == 15)) {cm3dkkp = -0.067;}     //N-P
      else if ((atomA == 8)&&(atomB == 14)) {cm3dkkp = 0.169;}      //O-Si
      else if ((atomA == 8)&&(atomB == 15)) {cm3dkkp = 0.033;}      //O-P
      else if ((atomA == 8)&&(atomB == 16)) {cm3dkkp = 0.103;}      //O-S
      else if ((atomA == 9)&&(atomB == 14)) {cm3dkkp = 0.228;}      //F-Si
      else if ((atomA == 9)&&(atomB == 15)) {cm3dkkp = 0.110;}      //F-P
      else if ((atomA == 14)&&(atomB == 17)) {cm3dkkp = -0.171;}    //Si-Cl
      else if ((atomA == 15)&&(atomB == 16)) {cm3dkkp = -0.685;}    //P-S
      else if ((atomA == 15)&&(atomB == 17)) {cm3dkkp = -0.252;}    //P-Cl
      if (atom1 != atomA) {cm3dkkp *= -1.0;}
    }
    else {std::cout << "WARNING: AM1.hpp: AM1: CM3dkkp(): partial charge calculation involving non-parametrized atom";}
    return cm3dkkp;
  }
  double B0() {return 0.482;}
};

#endif //_AM1_
