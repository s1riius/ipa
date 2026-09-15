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

//descrition:
//The PM3 method

class PM3: public MNDO {
  //this is the implementation of Stewart's PM3
  //J. J. P. Stewart, J. Comput. Chem., 10(2), 209, 1989
  //J. J. P. Stewart, J. Comput. Chem., 12(3), 320, 1991
  //E. Anders, R. Koch, P. Freunscht, J. Comput. Chem., 14(11), 1301, 1993
  //for charge and dipole moment model: J. W. Storer, D. J. Giesen, C. J. Cramer, D. G. Truhlar, J. Comp.-Aid. Mol. Des., 9, 87, 1995
  //                                    J. Li, T. Zhu, C. J. Cramer, D. G. Truhlar, J. Phys. Chem. A, 102, 1820, 1998                   -> CM2
  //bX quantities are "barred" tensors, which are used only for the open-shell case
public:
  PM3(BSet _bset, Molecule _mol, std::string _openclosed = "0", std::string _corecorrection = "0"): MNDO(_bset,_mol,_openclosed,_corecorrection) {}
  ~PM3() {}
  std::string Type() {return "PM3";}
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
      else if ((atoms[iatm] > 10)&&(atoms[iatm] < 18)) {def = true;}        //Na,Mg,Al,Si,P,S,Cl
      else if ((atoms[iatm] > 29)||(atoms[iatm] < 36)) {def = true;}        //Zn,Ga,Ge,As,Se,Br
      else if ((atoms[iatm] > 47)||(atoms[iatm] < 54)) {def = true;}        //Cd-I
      else if ((atoms[iatm] > 79)||(atoms[iatm] < 84)) {def = true;}        //Hg,Tl,Pb,Bi
      if (atoms[iatm] == 5) {def = false;}                                  //remove B
      if (!def) {throw("ERROR: PM3.hpp: PM3: checkAtoms(): atom not fully specified for PM3-theory");}
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
      case 3:      //Li
        enth = -0.194771427376;
        break;
      case 4:      //Be
        enth = -0.937719797485;
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
      case 9:      //F
        enth = -16.0784610392;
        break;
      case 11:     //Na
        enth = -0.19457272009;
        break;
      case 12:     //Mg
        enth = -0.828810340423;
        break;
      case 13:     //Al
        enth = -1.72224844965;
        break;
      case 14:     //Si
        enth = -2.46180717218;
        break;
      case 15:     //P
        enth = -4.32469507163;
        break;
      case 16:     //S
        enth = -6.70509021263;
        break;
      case 17:     //Cl
        enth = -11.5832018724;
        break;
      case 30:     //Zn
        enth = -1.00646113884;
        break;
      case 31:     //Ga
        enth = -2.10676627507;
        break;
      case 32:     //Ge
        enth = -3.07377982018;
        break;
      case 33:     //As
        enth = -4.50525281523;
        break;
      case 34:     //Se
        enth = -7.05046081862;
        break;
      case 35:     //Br
        enth = -12.9556035747;
        break;
      case 48:     //Cd
        enth = -0.825030010764;
        break;
      case 49:     //In
        enth = -1.91004794191;
        break;
      case 50:     //Sn
        enth = -2.89003979076;
        break;
      case 51:     //Sb
        enth = -5.47154425813;
        break;
      case 52:     //Te
        enth = -6.17695006812;
        break;
      case 53:     //I
        enth = -10.5954218289;
        break;
      case 80:     //Hg
        enth = -1.06204589076;
        break;
      case 81:     //Tl
        enth = -2.08182009765;
        break;
      case 82:     //Pb
        enth = -2.69729075229;
        break;
      case 83:     //Bi
        enth = -4.009248964;
        break;
    }
    return enth;
  }
  double betaA0(size_t atomicnr, int L) {
    //function returning the bonding parameters; values stored in eV, returned however in a.u.
    double betaa0 = 0.0;
    switch (atomicnr) {
      case 1:         //H
        betaa0 = -5.6265120;
        break;
      case 3:         //Li
        if (L == 0) {betaa0 = -0.5500000;}
        else if (L == 1) {betaa0 = -1.5000000;}
        break;
      case 4:         //Be
        if (L == 0) {betaa0 = -3.9620530;}
        else if (L == 1) {betaa0 = -2.7806840;}
        break;
      case 6:         //C
        if (L == 0) {betaa0 = -11.9100150;}
        else if (L == 1) {betaa0 = -9.8027550;}
        break;
      case 7:         //N
        if (L == 0) {betaa0 = -14.0625210;}
        else if (L == 1) {betaa0 = -20.0438480;}
        break;
      case 8:         //O
        if (L == 0) {betaa0 = -45.2026510;}
        else if (L == 1) {betaa0 = -24.7525150;}
        break;
      case 9:         //F
        if (L == 0) {betaa0 = -48.4059390;}
        else if (L == 1) {betaa0 = -27.7446600;}
        break;
      case 11:        //Na
        if (L == 0) {betaa0 = -0.1510870;}
        else if (L == 1) {betaa0 = -0.2184096;}
        break;
      case 12:        //Mg
        if (L == 0) {betaa0 = -2.0716910;}
        else if (L == 1) {betaa0 = -0.5695810;}
        break;
      case 13:        //Al
        if (L == 0) {betaa0 = -0.5943010;}
        else if (L == 1) {betaa0 = -0.9565500;}
        break;
      case 14:        //Si
        if (L == 0) {betaa0 = -2.8621450;}
        else if (L == 1) {betaa0 = -3.9331480;}
        break;
      case 15:        //P
        if (L == 0) {betaa0 = -12.6158790;}
        else if (L == 1) {betaa0 = -4.1600400;}
        break;
      case 16:        //S
        if (L == 0) {betaa0 = -8.8274650;}
        else if (L == 1) {betaa0 = -8.0914150;}
        break;
      case 17:        //Cl
        if (L == 0) {betaa0 = -27.5285600;}
        else if (L == 1) {betaa0 = -11.5939220;}
        break;
      case 30:        //Zn
        if (L == 0) {betaa0 = -0.7155780;}
        else if (L == 1) {betaa0 = -6.3518640;}
        break;
      case 31:        //Ga
        if (L == 0) {betaa0 = -4.9456180;}
        else if (L == 1) {betaa0 = -0.4070530;}
        break;
      case 32:        //Ge
        if (L == 0) {betaa0 = -5.3250024;}
        else if (L == 1) {betaa0 = -2.2501567;}
        break;
      case 33:        //As
        if (L == 0) {betaa0 = -8.2321650;}
        else if (L == 1) {betaa0 = -5.0173860;}
        break;
      case 34:        //Se
        if (L == 0) {betaa0 = -6.1578220;}
        else if (L == 1) {betaa0 = -5.4930390;}
        break;
      case 35:        //Br
        if (L == 0) {betaa0 = -31.1713420;}
        else if (L == 1) {betaa0 = -6.8140130;}
        break;
      case 48:        //Cd
        if (L == 0) {betaa0 = -8.5819440;}
        else if (L == 1) {betaa0 = -0.6010340;}
        break;
      case 49:        //In
        if (L == 0) {betaa0 = -2.9933190;}
        else if (L == 1) {betaa0 = -1.8289080;}
        break;
      case 50:        //Sn
        if (L == 0) {betaa0 = -2.7858020;}
        else if (L == 1) {betaa0 = -2.0059990;}
        break;
      case 51:        //Sb
        if (L == 0) {betaa0 = -14.7942170;}
        else if (L == 1) {betaa0 = -2.8179480;}
        break;
      case 52:        //Te
        if (L == 0) {betaa0 = -2.6651460;}
        else if (L == 1) {betaa0 = -3.8954300;}
        break;
      case 53:        //I
        if (L == 0) {betaa0 = -14.4942340;}
        else if (L == 1) {betaa0 = -5.8947030;}
        break;
      case 80:        //Hg
        if (L == 0) {betaa0 = -3.1013650;}
        else if (L == 1) {betaa0 = -3.4640310;}
        break;
      case 81:        //Tl
        if (L == 0) {betaa0 = -1.0844950;}
        else if (L == 1) {betaa0 = -7.9467990;}
        break;
      case 82:        //Pb
        if (L == 0) {betaa0 = -6.1260240;}
        else if (L == 1) {betaa0 = -1.3954300;}
        break;
      case 83:        //Bi
        if (L == 0) {betaa0 = -5.6072830;}
        else if (L == 1) {betaa0 = -5.8001520;}
        break;
    }
    return betaa0/au2eV;
  }
  double UlX(size_t atomicnr, int L) {
    //function that returns the atomic potential U
    double ulx = 0.0;
    switch (atomicnr) {
      case 1:                 //H
        ulx = -13.0733210;
        break;
      case 3:                 //Li
        if (L == 0) {ulx = -5.3000000;}
        else if (L == 1) {ulx = -3.4000000;}
        break;
      case 4:                 //Be
        if (L == 0) {ulx = -17.2647520;}
        else if (L == 1) {ulx = -11.3042430;}
        break;
      case 6:                 //C
        if (L == 0) {ulx = -47.2703200;}
        else if (L == 1) {ulx = -36.2669180;}
        break;
      case 7:                 //N
        if (L == 0) {ulx = -49.3356720;}
        else if (L == 1) {ulx = -47.5097360;}
        break;
      case 8:                 //O
        if (L == 0) {ulx = -86.9930020;}
        else if (L == 1) {ulx = -71.8795800;}
        break;
      case 9:                 //F
        if (L == 0) {ulx = -110.4353030;}
        else if (L == 1) {ulx = -105.6850470;}
        break;
      case 11:                //Na
        if (L == 0) {ulx = -5.2945929;}
        else if (L == 1) {ulx = -2.4596564;}
        break;
      case 12:                //Mg
        if (L == 0) {ulx = -14.6236880;}
        else if (L == 1) {ulx = -14.1734600;}
        break;
      case 13:                //Al
        if (L == 0) {ulx = -24.8454040;}
        else if (L == 1) {ulx = -22.2641590;}
        break;
      case 14:                //Si
        if (L == 0) {ulx = -26.7634830;}
        else if (L == 1) {ulx = -22.8136350;}
        break;
      case 15:                //P
        if (L == 0) {ulx = -40.4130960;}
        else if (L == 1) {ulx = -29.5930520;}
        break;
      case 16:                //S
        if (L == 0) {ulx = -49.8953710;}
        else if (L == 1) {ulx = -44.3925830;}
        break;
      case 17:                //Cl
        if (L == 0) {ulx = -100.6267470;}
        else if (L == 1) {ulx = -53.6143960;}
        break;
      case 30:                //Zn
        if (L == 0) {ulx = -18.5321980;}
        else if (L == 1) {ulx = -11.0474090;}
        break;
      case 31:                //Ga
        if (L == 0) {ulx = -29.8555930;}
        else if (L == 1) {ulx = -21.8753710;}
        break;
      case 32:                //Ge
        if (L == 0) {ulx = -35.4671955;}
        else if (L == 1) {ulx = -31.5863583;}
        break;
      case 33:                //As
        if (L == 0) {ulx = -38.5074240;}
        else if (L == 1) {ulx = -35.1524150;}
        break;
      case 34:                //Se
        if (L == 0) {ulx = -55.3781350;}
        else if (L == 1) {ulx = -49.8230760;}
        break;
      case 35:                //Br
        if (L == 0) {ulx = -116.6193110;}
        else if (L == 1) {ulx = -74.2271290;}
        break;
      case 48:                //Cd
        if (L == 0) {ulx = -15.8285840;}
        else if (L == 1) {ulx = 8.7497950;}
        break;
      case 49:                //In
        if (L == 0) {ulx = -26.1762050;}
        else if (L == 1) {ulx = -20.0058220;}
        break;
      case 50:                //Sn
        if (L == 0) {ulx = -34.5501920;}
        else if (L == 1) {ulx = -25.8944190;}
        break;
      case 51:                //Sb
        if (L == 0) {ulx = -56.4321960;}
        else if (L == 1) {ulx = -29.4349540;}
        break;
      case 52:                //Te
        if (L == 0) {ulx = -44.9380360;}
        else if (L == 1) {ulx = -46.3140990;}
        break;
      case 53:                //I
        if (L == 0) {ulx = -96.4540370;}
        else if (L == 1) {ulx = -61.09158209;}
        break;
      case 80:                //Hg
        if (L == 0) {ulx = -17.7622290;}
        else if (L == 1) {ulx = -18.3307510;}
        break;
      case 81:                //Tl
        if (L == 0) {ulx = -30.0531700;}
        else if (L == 1) {ulx = -26.9206370;}
        break;
      case 82:                //Pb
        if (L == 0) {ulx = -30.3227560;}
        else if (L == 1) {ulx = -24.4258340;}
        break;
      case 83:                //Bi
        if (L == 0) {ulx = -33.4959380;}
        else if (L == 1) {ulx = -35.5210260;}
        break;
    }
    return ulx/au2eV;
  }
  double alpha(int atomicnr, int atm2) {
    //function that returns alpha values for PM3
    //alphas stored and returned in 1/Angstrom
    double alpha = 0.0;
    switch (atomicnr) {
      case 1:     //H
        alpha = 3.356386;
        break;
      case 3:     //Li
        alpha = 1.255000;
        break;
      case 4:     //Be
        alpha = 1.593536;
        break;
      case 6:     //C
        alpha = 2.707807;
        break;
      case 7:     //N
        alpha = 2.830545;
        break;
      case 8:     //O
        alpha = 3.217102;
        break;
      case 9:     //F
        alpha = 3.358921;
        break;
      case 11:    //Na
        alpha = 2.3677169;
        break;
      case 12:    //Mg
        alpha = 1.329147;
        break;
      case 13:    //Al
        alpha = 1.521703;
        break;
      case 14:    //Si
        alpha = 2.135809;
        break;
      case 15:    //P
        alpha = 1.940534;
        break;
      case 16:    //S
        alpha = 2.269706;
        break;
      case 17:    //Cl
        alpha = 2.517296;
        break;
      case 30:    //Zn
        alpha = 1.350126;
        break;
      case 31:    //Ga
        alpha = 1.605115;
        break;
      case 32:    //Ge
        alpha = 1.972337;
        break;
      case 33:    //As
        alpha = 1.794477;
        break;
      case 34:    //Se
        alpha = 3.043957;
        break;
      case 35:    //Br
        alpha = 2.511842;
        break;
      case 48:    //Cd
        alpha = 1.525382;
        break;
      case 49:    //In
        alpha = 1.418385;
        break;
      case 50:    //Sn
        alpha = 1.699650;
        break;
      case 51:    //Sb
        alpha = 2.034301;
        break;
      case 52:    //Te
        alpha = 2.485019;
        break;
      case 53:    //I
        alpha = 1.990185;
        break;
      case 80:    //Hg
        alpha = 1.529377;
        break;
      case 81:    //Tl
        alpha = 1.340951;
        break;
      case 82:    //Pb
        alpha = 1.620045;
        break;
      case 83:    //Bi
        alpha = 1.857431;
        break;
    }
    return alpha;
  }
  std::vector<double> AM1K(int atomicnr, int atm2 = 0) {
    std::vector<double> am1k;
    switch (atomicnr) {
      case 1:     //H
        am1k.push_back(1.1287500/au2eV);
        am1k.push_back(-1.0603290/au2eV);
        break;
      case 3:     //Li
        am1k.push_back(-0.4500000/au2eV);
        am1k.push_back(0.8000000/au2eV);
        break;
      case 4:     //Be
        am1k.push_back(1.6315720/au2eV);
        am1k.push_back(-2.1109590/au2eV);
        break;
      case 6:     //C
        am1k.push_back(0.0501070/au2eV);
        am1k.push_back(0.0507330/au2eV);
        break;
      case 7:     //N
        am1k.push_back(1.5016740/au2eV);
        am1k.push_back(-1.5057720/au2eV);
        break;
      case 8:     //O
        am1k.push_back(-1.1311280/au2eV);
        am1k.push_back(1.1378910/au2eV);
        break;
      case 9:     //F
        am1k.push_back(-0.0121660/au2eV);
        am1k.push_back(-0.0028520/au2eV);
        break;
      case 11:    //Na
        am1k.push_back(0.6433655/au2eV);
        am1k.push_back(1.0871788/au2eV);
        break;
      case 12:    //Mg
        am1k.push_back(2.1170500/au2eV);
        am1k.push_back(-2.5477670/au2eV);
        break;
      case 13:    //Al
        am1k.push_back(-0.4730900/au2eV);
        am1k.push_back(-0.1540510/au2eV);
        break;
      case 14:    //Si
        am1k.push_back(-0.3906000/au2eV);
        am1k.push_back(0.0572590/au2eV);
        break;
      case 15:    //P
        am1k.push_back(-0.6114210/au2eV);
        am1k.push_back(-0.0939350/au2eV);
        break;
      case 16:    //S
        am1k.push_back(-0.3991910/au2eV);
        am1k.push_back(-0.0548990/au2eV);
        break;
      case 17:    //Cl
        am1k.push_back(-0.1715910/au2eV);
        am1k.push_back(-0.0134580/au2eV);
        break;
      case 30:    //Zn
        am1k.push_back(-0.1112340/au2eV);
        am1k.push_back(-0.1323700/au2eV);
        break;
      case 31:    //Ga
        am1k.push_back(-0.5601790/au2eV);
        am1k.push_back(-0.2727310/au2eV);
        break;
      case 32:    //Ge
        am1k.push_back(0.9631726/au2eV);
        am1k.push_back(-0.9593891/au2eV);
        break;
      case 33:    //As
        am1k.push_back(-0.4600950/au2eV);
        am1k.push_back(-0.0889960/au2eV);
        break;
      case 34:    //Se
        am1k.push_back(0.0478730/au2eV);
        am1k.push_back(0.1147200/au2eV);
        break;
      case 35:    //Br
        am1k.push_back(0.9604580/au2eV);
        am1k.push_back(-0.9549160/au2eV);
        break;
      case 48:    //Cd
        am1k.push_back(0.0);
        break;
      case 49:    //In
        am1k.push_back(-0.3431380/au2eV);
        am1k.push_back(-0.1095320/au2eV);
        break;
      case 50:    //Sn
        am1k.push_back(-0.1503530/au2eV);
        am1k.push_back(-0.0444170/au2eV);
        break;
      case 51:    //Sb
        am1k.push_back(3.0020280/au2eV);
        am1k.push_back(-0.0188920/au2eV);
        break;
      case 52:    //Te
        am1k.push_back(0.0333910/au2eV);
        am1k.push_back(-1.9218670/au2eV);
        break;
      case 53:    //I
        am1k.push_back(-0.1314810/au2eV);
        am1k.push_back(-0.0368970/au2eV);
        break;
      case 80:    //Hg
        am1k.push_back(1.0827200/au2eV);
        am1k.push_back(-0.0965530/au2eV);
        break;
      case 81:    //Tl
        am1k.push_back(-1.3613990/au2eV);
        am1k.push_back(-0.0454010/au2eV);
        break;
      case 82:    //Pb
        am1k.push_back(-0.1225760/au2eV);
        am1k.push_back(-0.0566480/au2eV);
        break;
      case 83:    //Bi
        am1k.push_back(2.5816930/au2eV);
        am1k.push_back(0.0603200/au2eV);
        break;
    }
    return am1k;
  }
  std::vector<double> AM1L(int atomicnr, int atm2 = 0) {
    std::vector<double> am1l;
    switch (atomicnr) {
      case 1:     //H
        am1l.push_back(5.0962820);
        am1l.push_back(6.0037880);
        break;
      case 3:     //Li
        am1l.push_back(5.0000000);
        am1l.push_back(6.5000000);
        break;
      case 4:     //Be
        am1l.push_back(2.6729620);
        am1l.push_back(1.9685940);
        break;
      case 6:     //C
        am1l.push_back(6.0031650);
        am1l.push_back(6.0029790);
        break;
      case 7:     //N
        am1l.push_back(5.9011480);
        am1l.push_back(6.0046580);
        break;
      case 8:     //O
        am1l.push_back(6.0024770);
        am1l.push_back(5.9505120);
        break;
      case 9:     //F
        am1l.push_back(6.0235740);
        am1l.push_back(6.0037170);
        break;
      case 11:    //Na
        am1l.push_back(1.5465054);
        am1l.push_back(1.4529000);
        break;
      case 12:    //Mg
        am1l.push_back(6.0094770);
        am1l.push_back(4.3953700);
        break;
      case 13:    //Al
        am1l.push_back(1.9158250);
        am1l.push_back(6.0050860);
        break;
      case 14:    //Si
        am1l.push_back(6.0000540);
        am1l.push_back(6.0071830);
        break;
      case 15:    //P
        am1l.push_back(1.9972720);
        am1l.push_back(1.9983600);
        break;
      case 16:    //S
        am1l.push_back(6.0006690);
        am1l.push_back(6.0018450);
        break;
      case 17:    //Cl
        am1l.push_back(6.0008020);
        am1l.push_back(1.9666180);
        break;
      case 30:    //Zn
        am1l.push_back(6.0014780);
        am1l.push_back(1.9958390);
        break;
      case 31:    //Ga
        am1l.push_back(5.6232730);
        am1l.push_back(1.9918430);
        break;
      case 32:    //Ge
        am1l.push_back(6.0120134);
        am1l.push_back(5.7491802);
        break;
      case 33:    //As
        am1l.push_back(1.9831150);
        am1l.push_back(1.9929440);
        break;
      case 34:    //Se
        am1l.push_back(6.0074000);
        am1l.push_back(6.0086720);
        break;
      case 35:    //Br
        am1l.push_back(5.9765080);
        am1l.push_back(5.9447030);
        break;
      case 48:    //Cd
        am1l.push_back(0.0);
        break;
      case 49:    //In
        am1l.push_back(1.9940340);
        am1l.push_back(5.6832170);
        break;
      case 50:    //Sn
        am1l.push_back(6.0056940);
        am1l.push_back(2.2573810);
        break;
      case 51:    //Sb
        am1l.push_back(6.0053420);
        am1l.push_back(6.0114780);
        break;
      case 52:    //Te
        am1l.push_back(5.9563790);
        am1l.push_back(4.9732190);
        break;
      case 53:    //I
        am1l.push_back(5.2064170);
        am1l.push_back(6.0101170);
        break;
      case 80:    //Hg
        am1l.push_back(6.4965980);
        am1l.push_back(3.9262810);
        break;
      case 81:    //Tl
        am1l.push_back(3.5572260);
        am1l.push_back(2.3069950);
        break;
      case 82:    //Pb
        am1l.push_back(6.0030620);
        am1l.push_back(4.7437050);
        break;
      case 83:    //Bi
        am1l.push_back(5.0940220);
        am1l.push_back(6.0015380);
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
      case 3:     //Li
        am1m.push_back(1.0000000);
        am1m.push_back(1.0000000);
        break;
      case 4:     //Be
        am1m.push_back(1.7916860);
        am1m.push_back(1.7558710);
        break;
      case 6:     //C
        am1m.push_back(1.6422140);
        am1m.push_back(0.8924880);
        break;
      case 7:     //N
        am1m.push_back(1.7107400);
        am1m.push_back(1.7161490);
        break;
      case 8:     //O
        am1m.push_back(1.6073110);
        am1m.push_back(1.5983950);
        break;
      case 9:     //F
        am1m.push_back(1.8568590);
        am1m.push_back(2.6361580);
        break;
      case 11:    //Na
        am1m.push_back(0.9976699);
        am1m.push_back(1.4506099);
        break;
      case 12:    //Mg
        am1m.push_back(2.0844060);
        am1m.push_back(2.0636740);
        break;
      case 13:    //Al
        am1m.push_back(1.4517280);
        am1m.push_back(2.5199970);
        break;
      case 14:    //Si
        am1m.push_back(0.6322620);
        am1m.push_back(2.0199870);
        break;
      case 15:    //P
        am1m.push_back(0.7946240);
        am1m.push_back(1.9106770);
        break;
      case 16:    //S
        am1m.push_back(0.9621230);
        am1m.push_back(1.5799440);
        break;
      case 17:    //Cl
        am1m.push_back(1.0875020);
        am1m.push_back(2.2928910);
        break;
      case 30:    //Zn
        am1m.push_back(1.5160320);
        am1m.push_back(2.5196420);
        break;
      case 31:    //Ga
        am1m.push_back(1.5317800);
        am1m.push_back(2.1838640);
        break;
      case 32:    //Ge
        am1m.push_back(2.1633655);
        am1m.push_back(2.1693724);
        break;
      case 33:    //As
        am1m.push_back(1.0867930);
        am1m.push_back(2.1400580);
        break;
      case 34:    //Se
        am1m.push_back(2.0817170);
        am1m.push_back(1.5164230);
        break;
      case 35:    //Br
        am1m.push_back(2.3216540);
        am1m.push_back(2.3281420);
        break;
      case 48:    //Cd
        am1m.push_back(0.0);
        break;
      case 49:    //In
        am1m.push_back(1.6255160);
        am1m.push_back(2.8670090);
        break;
      case 50:    //Sn
        am1m.push_back(1.7046420);
        am1m.push_back(2.4698690);
        break;
      case 51:    //Sb
        am1m.push_back(0.8530600);
        am1m.push_back(2.7933110);
        break;
      case 52:    //Te
        am1m.push_back(2.2775750);
        am1m.push_back(0.5242430);
        break;
      case 53:    //I
        am1m.push_back(1.7488240);
        am1m.push_back(2.7103730);
        break;
      case 80:    //Hg
        am1m.push_back(1.1951460);
        am1m.push_back(2.6271600);
        break;
      case 81:    //Tl
        am1m.push_back(1.0928020);
        am1m.push_back(2.9650290);
        break;
      case 82:    //Pb
        am1m.push_back(1.9015970);
        am1m.push_back(2.8618790);
        break;
      case 83:    //Bi
        am1m.push_back(0.4997870);
        am1m.push_back(2.4279700);
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
      case 3:                  //Li
        if (idx == 1) {dval = 1.0772805453024;}
        else if (idx == 2) {dval = 0.8641427827586;}
        break;
      case 4:                  //Be
        if (idx == 1) {dval = 0.5339679024905;}
        else if (idx == 2) {dval = 0.4295641497663;}
        break;
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
      case 9:                  //F
        if (idx == 1) {dval = 0.1653838587473;}
        else if (idx == 2) {dval = 0.2601608726156;}
        break;
      case 11:                 //Na
        if (idx == 1) {dval = 0.91824824182;}
        else if (idx == 2) {dval = 0.74551702212;}
        break;
      case 12:                 //Mg
        if (idx == 1) {dval = 0.6034710424661;}
        else if (idx == 2) {dval = 0.5969065462793;}
        break;
      case 13:                 //Al
        if (idx == 1) {dval = 0.6404525387508;}
        else if (idx == 2) {dval = 0.8247568110748;}
        break;
      case 14:                 //Si
        if (idx == 1) {dval = 0.6955796273438;}
        else if (idx == 2) {dval = 0.6743514719617;}
        break;
      case 15:                 //P
        if (idx == 1) {dval = 0.5633063336025;}
        else if (idx == 2) {dval = 0.5884654818765;}
        break;
      case 16:                 //S
        if (idx == 1) {dval = 0.5934358846409;}
        else if (idx == 2) {dval = 0.5337539561451;}
        break;
      case 17:                 //Cl
        if (idx == 1) {dval = 0.4855653861897;}
        else if (idx == 2) {dval = 0.4116591214170;}
        break;
      case 30:                 //Zn
        if (idx == 1) {dval = 0.7940705126955;}
        else if (idx == 2) {dval = 0.7449319638157;}
        break;
      case 31:                 //Ga
        if (idx == 1) {dval = 0.5173602579027;}
        else if (idx == 2) {dval = 1.3373119811729;}
        break;
      case 32:                 //Ge
        if (idx == 1) {dval = 0.6307953192878;}
        else if (idx == 2) {dval = 0.7049308765449;}
        break;
      case 33:                 //As
        if (idx == 1) {dval = 0.5122252810265;}
        else if (idx == 2) {dval = 0.6588189567081;}
        break;
      case 34:                 //Se
        if (idx == 1) {dval = 0.4614326300290;}
        else if (idx == 2) {dval = 0.6479255792865;}
        break;
      case 35:                 //Br
        if (idx == 1) {dval = 0.1460013147146;}
        else if (idx == 2) {dval = 0.5276178289085;}
        break;
      case 48:                 //Cd
        if (idx == 1) {dval = 0.8457670512825;}
        else if (idx == 2) {dval = 0.6578943782897;}
        break;
      case 49:                 //In
        if (idx == 1) {dval = 0.8343135397859;}
        else if (idx == 2) {dval = 0.9405893627198;}
        break;
      case 50:                 //Sn
        if (idx == 1) {dval = 0.6942825081708;}
        else if (idx == 2) {dval = 0.8298458553693;}
        break;
      case 51:                 //Sb
        if (idx == 1) {dval = 0.7457113889258;}
        else if (idx == 2) {dval = 0.7155192362237;}
        break;
      case 52:                 //Te
        if (idx == 1) {dval = 0.1843747058104;}
        else if (idx == 2) {dval = 0.8251505189177;}
        break;
      case 53:                 //I
        if (idx == 1) {dval = 0.0836877350442;}
        else if (idx == 2) {dval = 0.5539057650855;}
        break;
      case 80:                 //Hg
        if (idx == 1) {dval = 0.6518304837420;}
        else if (idx == 2) {dval = 0.6436929024681;}
        break;
      case 81:                 //Tl
        if (idx == 1) {dval = 0.0413478961836;}
        else if (idx == 2) {dval = 0.8105465509115;}
        break;
      case 82:                 //Pb
        if (idx == 1) {dval = 0.5221015798537;}
        else if (idx == 2) {dval = 0.8435382097988;}
        break;
      case 83:                 //Bi
        if (idx == 1) {dval = 0.1480960097760;}
        else if (idx == 2) {dval = 0.8250028255589;}
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
      case 3:                  //Li
        if (l == 0) {rho = 1.599878849930;}
        else if (l == 1) {rho = 2.287368480615;}
        else if (l == 2) {rho = 0.835697784494;}
        break;
      case 4:                  //Be
        if (l == 0) {rho = 0.798799045243;}
        else if (l == 1) {rho = 0.909553000943;}
        else if (l == 2) {rho = 0.749541089292;}
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
      case 9:                  //F
        if (l == 0) {rho = 0.685880274657;}
        else if (l == 1) {rho = 0.390911556297;}
        else if (l == 2) {rho = 0.432331000322;}
        break;
      case 11:                 //Na
        if (l == 0) {rho = 1.81994308925;}
        else if (l == 1) {rho = 1.24743397132;}
        else if (l == 2) {rho = 1.02719786193;}
        break;
      case 12:                 //Mg
        if (l == 0) {rho = 1.075460694417;}
        else if (l == 1) {rho = 0.981502387101;}
        else if (l == 2) {rho = 0.956048783543;}
        break;
      case 13:                 //Al
        if (l == 0) {rho = 1.246284086598;}
        else if (l == 1) {rho = 0.412222702317;}
        else if (l == 2) {rho = 1.169277713000;}
        break;
      case 14:                 //Si
        if (l == 0) {rho = 1.426426712705;}
        else if (l == 1) {rho = 0.864466649632;}
        else if (l == 2) {rho = 0.542475229082;}
        break;
      case 15:                 //P
        if (l == 0) {rho = 0.922816001025;}
        else if (l == 1) {rho = 0.613973592320;}
        else if (l == 2) {rho = 0.708874478410;}
        break;
      case 16:                 //S
        if (l == 0) {rho = 0.803092446274;}
        else if (l == 1) {rho = 0.396143029893;}
        else if (l == 2) {rho = 0.431103562107;}
        break;
      case 17:                 //Cl
        if (l == 0) {rho = 0.449583792797;}
        else if (l == 1) {rho = 0.388271699986;}
        else if (l == 2) {rho = 0.726154842182;}
        break;
      case 30:                 //Zn
        if (l == 0) {rho = 0.743960973132;}
        else if (l == 1) {rho = 1.113734180497;}
        else if (l == 2) {rho = 0.994294413835;}
        break;
      case 31:                 //Ga
        if (l == 0) {rho = 0.851144894368;}
        else if (l == 1) {rho = 0.515831612819;}
        else if (l == 2) {rho = 1.711209644194;}
        break;
      case 32:                 //Ge
        if (l == 0) {rho = 1.338944749466;}
        else if (l == 1) {rho = 0.696619077635;}
        else if (l == 2) {rho = 0.730772694640;}
        break;
      case 33:                 //As
        if (l == 0) {rho = 0.819143788072;}
        else if (l == 1) {rho = 0.524744273631;}
        else if (l == 2) {rho = 1.027840304780;}
        break;
      case 34:                 //Se
        if (l == 0) {rho = 0.968633392541;}
        else if (l == 1) {rho = 0.352329267283;}
        else if (l == 2) {rho = 0.500760359861;}
        break;
      case 35:                 //Br
        if (l == 0) {rho = 0.451562701487;}
        else if (l == 1) {rho = 0.391670767080;}
        else if (l == 2) {rho = 0.691966653807;}
        break;
      case 48:                 //Cd
        if (l == 0) {rho = 0.781957934859;}
        else if (l == 1) {rho = 0.741084349319;}
        else if (l == 2) {rho = 0.938063861121;}
        break;
      case 49:                 //In
        if (l == 0) {rho = 1.098331941968;}
        else if (l == 1) {rho = 0.583738678882;}
        else if (l == 2) {rho = 0.717078821775;}
        break;
      case 50:                 //Sn
        if (l == 0) {rho = 0.706519361502;}
        else if (l == 1) {rho = 0.822172786627;}
        else if (l == 2) {rho = 0.934107308885;}
        break;
      case 51:                 //Sb
        if (l == 0) {rho = 0.779307247171;}
        else if (l == 1) {rho = 0.576570118943;}
        else if (l == 2) {rho = 1.091774958213;}
        break;
      case 52:                 //Te
        if (l == 0) {rho = 0.702038451270;}
        else if (l == 1) {rho = 0.221214187240;}
        else if (l == 2) {rho = 1.211050437677;}
        break;
      case 53:                 //I
        if (l == 0) {rho = 0.528131297085;}
        else if (l == 1) {rho = 0.158444790902;}
        else if (l == 2) {rho = 0.513457003310;}
        break;
      case 80:                 //Hg
        if (l == 0) {rho = 1.086756136177;}
        else if (l == 1) {rho = 0.585960015150;}
        else if (l == 2) {rho = 1.010499581625;}
        break;
      case 81:                 //Tl
        if (l == 0) {rho = 0.688257458272;}
        else if (l == 1) {rho = 0.102785527809;}
        else if (l == 2) {rho = 1.195467888842;}
        break;
      case 82:                 //Pb
        if (l == 0) {rho = 1.026734684572;}
        else if (l == 1) {rho = 0.584380347173;}
        else if (l == 2) {rho = 1.230544509879;}
        break;
      case 83:                 //Bi
        if (l == 0) {rho = 1.442927491968;}
        else if (l == 1) {rho = 0.390478300666;}
        else if (l == 2) {rho = 1.022952090658;}
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
      case 3:       //Li
        if ((Lbra == 0)&&(Lket == 0)) {eri = 4.5000000;}                                         //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 5.2500000;}                                    //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 3.0000000;}      //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 4.5000000;}    //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 0.15000003;}                                   //(sp|sp)||(ps|ps)
        break;
      case 4:       //Be
        if ((Lbra == 0)&&(Lket == 0)) {eri = 9.0128510;}                                         //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 6.0571820;}                                    //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 6.5761990;}      //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 9.0052190;}    //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 0.5446790;}                                    //(sp|sp)||(ps|ps)
        break;
      case 5:       //B
        if ((Lbra == 0)&&(Lket == 0)) {eri = 10.59;}                                               //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 8.86;}                                           //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 9.56;}             //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 7.86;}           //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 1.81;}                                           //(sp|sp)||(ps|ps)
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
      case 11:      //Na
        if ((Lbra == 0)&&(Lket == 0)) {eri = 3.9558692;}                                         //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 5.3363963;}                                    //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 7.1929109;}      //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 5.0588074;}    //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 0.5687889;}                                    //(sp|sp)||(ps|ps)
        break;
      case 12:      //Mg
        if ((Lbra == 0)&&(Lket == 0)) {eri = 6.6943000;}                                         //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 6.9104460;}                                    //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 6.7939950;}      //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 7.0908230;}    //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 0.5433000;}                                    //(sp|sp)||(ps|ps)
        break;
      case 13:      //Al
        if ((Lbra == 0)&&(Lket == 0)) {eri = 5.7767370;}                                         //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 6.3477900;}                                    //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 11.6598560;}     //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 6.1210770;}    //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 4.0062450;}                                    //(sp|sp)||(ps|ps)
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
      case 30:      //Zn
        if ((Lbra == 0)&&(Lket == 0)) {eri = 9.6771960;}                                         //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 4.9801740;}                                    //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 7.7362040;}      //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 4.6696560;}    //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 0.6004130;}                                    //(sp|sp)||(ps|ps)
        break;
      case 31:      //Ga
        if ((Lbra == 0)&&(Lket == 0)) {eri = 8.4585540;}                                         //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 5.0868550;}                                    //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 8.9256190;}      //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 4.9830450;}    //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 2.0512600;}                                    //(sp|sp)||(ps|ps)
        break;
      case 32:      //Ge
        if ((Lbra == 0)&&(Lket == 0)) {eri = 5.3769635;}                                       //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 7.6718647;}                                  //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 10.2095293;}   //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 6.9242663;}  //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 1.3370204;}                                  //(sp|sp)||(ps|ps)
        break;
      case 33:      //As
        if ((Lbra == 0)&&(Lket == 0)) {eri = 8.7890010;}                                         //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 8.2872500;}                                    //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 5.3979830;}      //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 8.2103460;}    //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 1.9510340;}                                    //(sp|sp)||(ps|ps)
        break;
      case 34:      //Se
        if ((Lbra == 0)&&(Lket == 0)) {eri = 7.4325910;}                                         //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 9.5683260;}                                    //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 10.0604610;}     //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 7.7242890;}    //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 4.0165580;}                                    //(sp|sp)||(ps|ps)
        break;
      case 35:      //Br
        if ((Lbra == 0)&&(Lket == 0)) {eri = 15.9434250;}                                        //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 8.2827630;}                                    //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 16.0616800;}     //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 7.8168490;}    //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 0.5788690;}                                    //(sp|sp)||(ps|ps)
        break;
      case 48:      //Cd
        if ((Lbra == 0)&&(Lket == 0)) {eri = 9.2069600;}                                         //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 4.9481040;}                                    //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 8.2315390;}      //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 4.6696560;}    //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 1.6562340;}                                    //(sp|sp)||(ps|ps)
        break;
      case 49:      //In
        if ((Lbra == 0)&&(Lket == 0)) {eri = 6.5549000;}                                         //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 6.2992690;}                                    //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 8.2298730;}      //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 4.9842110;}    //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 2.6314610;}                                    //(sp|sp)||(ps|ps)
        break;
      case 50:      //Sn
        if ((Lbra == 0)&&(Lket == 0)) {eri = 10.1900330;}                                        //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 5.6738100;}                                    //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 7.2353270;}      //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 5.1822140;}    //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 1.0331570;}                                    //(sp|sp)||(ps|ps)
        break;
      case 51:      //Sb
        if ((Lbra == 0)&&(Lket == 0)) {eri = 9.2382770;}                                         //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 6.3500000;}                                    //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 5.2776800;}      //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 6.2500000;}    //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 2.4244640;}                                    //(sp|sp)||(ps|ps)
        break;
      case 52:      //Te
        if ((Lbra == 0)&&(Lket == 0)) {eri = 10.2550730;}                                        //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 7.7775920;}                                    //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 8.1691450;}      //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 7.7551210;}    //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 3.7724620;}                                    //(sp|sp)||(ps|ps)
        break;
      case 53:      //I
        if ((Lbra == 0)&&(Lket == 0)) {eri = 13.6319430;}                                        //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 7.2883300;}                                    //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 14.9904060;}     //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 5.9664070;}    //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 2.6300350;}                                    //(sp|sp)||(ps|ps)
        break;
      case 80:      //Hg
        if ((Lbra == 0)&&(Lket == 0)) {eri = 6.6247200;}                                         //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 14.7092830;}                                   //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 10.6392970;}     //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 16.0007400;}   //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 2.0363110;}                                    //(sp|sp)||(ps|ps)v
        break;
      case 81:      //Tl
        if ((Lbra == 0)&&(Lket == 0)) {eri = 10.4604120;}                                        //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 4.9927850;}                                    //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 11.2238830;}     //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 8.9627270;}    //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 2.5304060;}                                    //(sp|sp)||(ps|ps)
        break;
      case 82:      //Pb
        if ((Lbra == 0)&&(Lket == 0)) {eri = 7.0119920;}                                         //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 5.1837800;}                                    //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 6.7937820;}      //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 5.0456510;}    //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 1.5663020;}                                    //(sp|sp)||(ps|ps)
        break;
      case 83:      //Bi
        if ((Lbra == 0)&&(Lket == 0)) {eri = 4.9894800;}                                         //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 8.6960070;}                                    //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 6.1033080;}      //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 8.3354470;}    //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 0.5991220;}                                    //(sp|sp)||(ps|ps)
        break;
    }
    return eri/au2eV;
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
    else {std::cout << "ERROR: PM3.hpp: PM3: CM1ck(): partial charge calculation involving non-parametrized atom";}
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
    else {std::cout << "ERROR: PM3.hpp: PM3: CM1dk(): partial charge calculation involving non-parametrized atom";}
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
    else {std::cout << "ERROR: PM3.hpp: PM3: CM1dkkp(): partial charge calculation involving non-parametrized atom";}
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
    else {std::cout << "ERROR: PM3.hpp: PM3: CM2dkkp(): partial charge calculation involving non-parametrized atom";}
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
    else {std::cout << "WARNING: PM3.hpp: PM3: CM3dkkp(): partial charge calculation involving non-parametrized atom";}
    return cm3dkkp;
  }
  double B0() {return 0.460;}
};

#endif //_PM3_
