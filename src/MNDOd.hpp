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

#ifndef _MNDOD_
#define _MNDOD_
#include <string>
#include <sstream>
#include <vector>
#include "MNDO.hpp"
//FMCM??
#include "basissets/HcoreFunctions.hpp"

//descrition:
//NDDO methods with spd bases

class MNDOd: public MNDO {
  //implementation of Thiel's and Voityuk's MNDOd
  //W. Thiel, A. A. Voityuk, J. Phys. Chem., 100, 616, 1996
  //bX quantities are "barred" tensors, which are used only for the open-shell case
public:
  MNDOd(BSet _bset, Molecule _mol, std::string _openclosed = "0", std::string _corecorrection = "0"): MNDO(_bset,_mol,_openclosed,_corecorrection) {
    Dd.resize(2);
    intn.resize(31);
  }
  ~MNDOd() {}
  //getters
  virtual std::string Type() {return "MNDOd";}
  bool AtomWithDOrbitals(size_t atomicnr) {
    //function that determines whether atom has d orbitals
    return (((atomicnr > 12)&&(atomicnr < 18))||((atomicnr > 20)&&(atomicnr < 30))||(atomicnr == 33)||(atomicnr == 35)||((atomicnr > 38)&&(atomicnr < 48))||(atomicnr == 51)||(atomicnr == 53)||(atomicnr == 57)||((atomicnr > 70)&&(atomicnr < 80)));
  }
  double SlaterCondonRadialIntegral(int atmnr, int L, int II, int JJ, double threshold = 1e-7) {
    //function returning Slater Condon radial integrals used to calculate 1-Center eris
    //II and JJ refer to the orbitals in the bra and ket respectively
    //they take the values
    //        1 -> ss; 2 -> sp; 3 -> pp; 4 -> sd; 5 -> pd; 6 -> dd
    //L is the radial number of the Slater Condon integral
    //these integrals translate into
    //    L,II,JJ = 0,6,6 -> F0dd;  2,6,6 -> F2dd;  4,6,6 -> F4dd;  0,1,6 -> F0sd;  2,4,4 -> G2sd;
    //            = 0,3,6 -> F0pd;  2,3,6 -> F2pd;  1,5,5 -> G1pd;  3,5,5 -> G3pd;
    //more combinations possible though
    bool calculate = true;
    double scint = 0.0;
    //check whether there are tabulated values so that calculation is skipped
    if ((L == 0)&&(II == 1)&&(JJ == 6)) {
      scint = F0sd(atmnr);
      if (fabs(scint) > threshold) {calculate = false;}
    }
    else if ((L == 2)&&(II == 4)&&(JJ == 4)) {
      scint = G2sd(atmnr);
      if (fabs(scint) > threshold) {calculate = false;}
    }
    if (calculate) {
      //get principal quantum numbers of orbitals in the integral (mu,nu|lambda,sigma)
      size_t l_1 = 0;                //s
      size_t l_2 = 0;                //s
      if (II == 2) {l_2 = 1;}        //sp
      else if (II == 3) {            //pp
        l_1 = 1;
        l_2 = 1;
      }
      else if (II == 4) {l_2 = 2;}   //sd
      else if (II == 5) {            //pd
        l_1 = 1;
        l_2 = 2;
      }
      else if (II == 6) {            //dd
        l_1 = 2;
        l_2 = 2;
      }
      double n_mu = PrincipalQuantNumb(atmnr,l_1);
      double n_nu = PrincipalQuantNumb(atmnr,l_2);
      double z_mu = zn(atmnr,l_1);
      double z_nu = zn(atmnr,l_2);
      l_1 = 0;                       //s
      l_2 = 0;                       //s
      if (JJ == 2) {l_2 = 1;}        //sp
      else if (JJ == 3) {            //pp
        l_1 = 1;
        l_2 = 1;
      }
      else if (JJ == 4) {l_2 = 2;}   //sd
      else if (JJ == 5) {            //pd
        l_1 = 1;
        l_2 = 2;
      }
      else if (JJ == 6) {            //dd
        l_1 = 2;
        l_2 = 2;
      }
      double n_lambda = PrincipalQuantNumb(atmnr,l_1);
      double n_sigma = PrincipalQuantNumb(atmnr,l_2);
      double z_lambda = zn(atmnr,l_1);
      double z_sigma = zn(atmnr,l_2);
      int nmn = int(n_mu + n_nu);
      int nls = int(n_lambda + n_sigma);
      int nmnls = nmn + nls;
      double lzm = log(z_mu);
      double lzn = log(z_nu);
      double lzl = log(z_lambda);
      double lzs = log(z_sigma);
      double zls = z_lambda + z_sigma;
      double zmn = z_mu + z_nu;
      double zmnls = zls + zmn;
      double lzmnls = log(zmnls);
      double lzls = log(zls);
      double lzmn = log(zmn);
      double fac_a = Factorial(nmnls - 1)/sqrt(Factorial(2*int(n_mu))*Factorial(2*int(n_nu))*Factorial(2*int(n_lambda))*Factorial(2*int(n_sigma)));
      double fac0 = fac_a*exp(n_mu*lzm + n_nu*lzn + n_lambda*lzl + n_sigma*lzs + 0.5*(lzm + lzn + lzl + lzs) + log(2.0)*double(nmnls + 2) - lzmnls*double(nmnls));
      double aux = 1.0/zmnls;
      double fac1 = 0.0;
      double fac2 = 0.0;
      double fac3;
      for (size_t i = 0; i < nls - L; i++) {
        aux = aux*zmnls/zls;
        fac1 += aux*(dCombinations(nls - L - 1,i) - dCombinations(nls + L,i))/dCombinations(nmnls - 1,i);
      }
      for (size_t i = nls - L; i < nls + L + 1; i++) {
        aux = aux*zmnls/zls;
        fac2 += aux*dCombinations(nls + L,i)/dCombinations(nmnls - 1,i);
      }
      fac3 = exp(lzmnls*nmnls - lzls*(nls + L + 1) - lzmn*(nmn - L))/dCombinations(nmnls - 1,nls + L);
      scint = fac0*(fac1 - fac2 + fac3);
    }
    return scint;
  }
  double d_eri2Center(int L1, int M1, int L2, int M2, int L3, int M3, int L4, int M4, double RAB, std::vector<double> D,size_t atmnrA,size_t atmnrB, int core = 0) {
    //function that calculates 2 electron repulsion integrals as in the NDDO approximation MNDO/d
    //AOs contains AO information on the orbitals involved in the integral, which is valid for s,p,d basis/atoms
    //RAB is distance between centers
    //Thiel, Voityuk, Theor. Chim. Acta, 81, 391, 1992
    double eri = 0.0;
    double rhoA;
    double rhoB;
    int lA;
    int lB;
    int mA;
    int mB;
    int l_pairA;
    int l_pairB;
    chg1 = NonZeroMultipoleMom_d(L1,M1,L2,M2);
    chg2 = NonZeroMultipoleMom_d(L3,M3,L4,M4);
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
          if ((L1 == 0)&&(L2 == 0)) {l_pairA = 0;}
          else if (((L1 == 0)&&(L2 == 1))||((L1 == 1)&&(L2 == 0))) {l_pairA = 1;}
          else if ((L1 == 1)&&(L2 == 1)) {
            l_pairA = 2;
            if (lA == 0) {l_pairA = 0;}
          }
          else if (((L1 == 0)&&(L2 == 2))||((L1 == 2)&&(L2 == 0))) {l_pairA = 3;}
          else if (((L1 == 1)&&(L2 == 2))||((L1 == 2)&&(L2 == 1))) {l_pairA = 4;}
          else if ((L1 == 2)&&(L2 == 2)) {
            l_pairA = 5;
            if (lA == 2) {l_pairA = 6;}
            //if ((lA == 2)&&(atmnrB < 11)) {l_pairA = 4;}
          }
          if ((L3 == 0)&&(L4 == 0)) {l_pairB = 0;}
          else if (((L3 == 0)&&(L4 == 1))||((L3 == 1)&&(L4 == 0))) {l_pairB = 1;}
          else if ((L3 == 1)&&(L4 == 1)) {
            l_pairB = 2;
            if (lB == 0) {l_pairB = 0;}
          }
          else if (((L3 == 0)&&(L4 == 2))||((L3 == 2)&&(L4 == 0))) {l_pairB = 3;}
          else if (((L3 == 1)&&(L4 == 2))||((L3 == 2)&&(L4 == 1))) {l_pairB = 4;}
          else if ((L3 == 2)&&(L4 == 2)) {
            l_pairB = 5;
            if (lB == 2) {l_pairB = 6;}
            //if ((lB == 2)&&(atmnrA < 11)) {l_pairB = 4;}
          }
          rhoA = rho(atmnrA,l_pairA);
          rhoB = rho(atmnrB,l_pairB);
          if (core > 0) {rhoB = rhocore(atmnrB);}                        //in this case use the correct core rho for the ss atom
          if ((core > 1)||(core < 0)) {rhoA = rhocore(atmnrA);}          //in this case use the correct core rho for the ss atom
          eri += d_se_multipole(lA,mA,lB,mB,RAB,rhoA + rhoB,D)*clm_munu(L1,M1,L2,M2,lA,mA)*clm_munu(L3,M3,L4,M4,lB,mB);
        }
      }
    }
    return eri;
  }
  double d_eri2Center_dR(int L1, int M1, int L2, int M2, int L3, int M3, int L4, int M4, double RAB, std::vector<double> D,size_t atmnrA,size_t atmnrB, int core = 0) {
    //function that calculates the first-derivatives of 2 electron repulsion integrals as in the NDDO approximation MNDO/d with respect to internuclear distance
    //AOs contains AO information on the orbitals involved in the integral, which is valid for s,p,d basis/atoms
    //RAB is distance between centers
    double eri = 0.0;
    double rhoA;
    double rhoB;
    int lA;
    int lB;
    int mA;
    int mB;
    int l_pairA;
    int l_pairB;
    chg1 = NonZeroMultipoleMom_d(L1,M1,L2,M2);
    chg2 = NonZeroMultipoleMom_d(L3,M3,L4,M4);
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
          if ((L1 == 0)&&(L2 == 0)) {l_pairA = 0;}
          else if (((L1 == 0)&&(L2 == 1))||((L1 == 1)&&(L2 == 0))) {l_pairA = 1;}
          else if ((L1 == 1)&&(L2 == 1)) {
            l_pairA = 2;
            if (lA == 0) {l_pairA = 0;}
          }
          else if (((L1 == 0)&&(L2 == 2))||((L1 == 2)&&(L2 == 0))) {l_pairA = 3;}
          else if (((L1 == 1)&&(L2 == 2))||((L1 == 2)&&(L2 == 1))) {l_pairA = 4;}
          else if ((L1 == 2)&&(L2 == 2)) {
            l_pairA = 5;
            if (lA == 2) {l_pairA = 6;}
          }
          if ((L3 == 0)&&(L4 == 0)) {l_pairB = 0;}
          else if (((L3 == 0)&&(L4 == 1))||((L3 == 1)&&(L4 == 0))) {l_pairB = 1;}
          else if ((L3 == 1)&&(L4 == 1)) {
            l_pairB = 2;
            if (lB == 0) {l_pairB = 0;}
          }
          else if (((L3 == 0)&&(L4 == 2))||((L3 == 2)&&(L4 == 0))) {l_pairB = 3;}
          else if (((L3 == 1)&&(L4 == 2))||((L3 == 2)&&(L4 == 1))) {l_pairB = 4;}
          else if ((L3 == 2)&&(L4 == 2)) {
            l_pairB = 5;
            if (lB == 2) {l_pairB = 6;}
          }
          rhoA = rho(atmnrA,l_pairA);
          rhoB = rho(atmnrB,l_pairB);
          if (core > 0) {rhoB = rhocore(atmnrB);}                        //in this case use the correct core rho for the ss atom
          if ((core > 1)||(core < 0)) {rhoA = rhocore(atmnrA);}          //in this case use the correct core rho for the ss atom
          eri += d_se_multipole_dR(lA,mA,lB,mB,RAB,rhoA + rhoB,D)*clm_munu(L1,M1,L2,M2,lA,mA)*clm_munu(L3,M3,L4,M4,lB,mB);
        }
      }
    }
    return eri;
  }
  double d_eri2Center_dR2(int L1, int M1, int L2, int M2, int L3, int M3, int L4, int M4, double RAB, std::vector<double> D,size_t atmnrA,size_t atmnrB, int core = 0) {
    //function that calculates the second-derivatives of 2 electron repulsion integrals as in the NDDO approximation MNDO/d with respect to internuclear distance
    //AOs contains AO information on the orbitals involved in the integral, which is valid for s,p,d basis/atoms
    //RAB is distance between centers
    double eri = 0.0;
    double rhoA;
    double rhoB;
    int lA;
    int lB;
    int mA;
    int mB;
    int l_pairA;
    int l_pairB;
    chg1 = NonZeroMultipoleMom_d(L1,M1,L2,M2);
    chg2 = NonZeroMultipoleMom_d(L3,M3,L4,M4);
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
          if ((L1 == 0)&&(L2 == 0)) {l_pairA = 0;}
          else if (((L1 == 0)&&(L2 == 1))||((L1 == 1)&&(L2 == 0))) {l_pairA = 1;}
          else if ((L1 == 1)&&(L2 == 1)) {
            l_pairA = 2;
            if (lA == 0) {l_pairA = 0;}
          }
          else if (((L1 == 0)&&(L2 == 2))||((L1 == 2)&&(L2 == 0))) {l_pairA = 3;}
          else if (((L1 == 1)&&(L2 == 2))||((L1 == 2)&&(L2 == 1))) {l_pairA = 4;}
          else if ((L1 == 2)&&(L2 == 2)) {
            l_pairA = 5;
            if (lA == 2) {l_pairA = 6;}
          }
          if ((L3 == 0)&&(L4 == 0)) {l_pairB = 0;}
          else if (((L3 == 0)&&(L4 == 1))||((L3 == 1)&&(L4 == 0))) {l_pairB = 1;}
          else if ((L3 == 1)&&(L4 == 1)) {
            l_pairB = 2;
            if (lB == 0) {l_pairB = 0;}
          }
          else if (((L3 == 0)&&(L4 == 2))||((L3 == 2)&&(L4 == 0))) {l_pairB = 3;}
          else if (((L3 == 1)&&(L4 == 2))||((L3 == 2)&&(L4 == 1))) {l_pairB = 4;}
          else if ((L3 == 2)&&(L4 == 2)) {
            l_pairB = 5;
            if (lB == 2) {l_pairB = 6;}
          }
          rhoA = rho(atmnrA,l_pairA);
          rhoB = rho(atmnrB,l_pairB);
          if (core > 0) {rhoB = rhocore(atmnrB);}                        //in this case use the correct core rho for the ss atom
          if ((core > 1)||(core < 0)) {rhoA = rhocore(atmnrA);}          //in this case use the correct core rho for the ss atom
          eri += d_se_multipole_dR2(lA,mA,lB,mB,RAB,rhoA + rhoB,D)*clm_munu(L1,M1,L2,M2,lA,mA)*clm_munu(L3,M3,L4,M4,lB,mB);
        }
      }
    }
    return eri;
  }
  //other auxiliary functions
  virtual void checkAtoms() {
    //function that checks whether the atoms in the molecule are defined in MNDOd
    bool def;
    for (size_t iatm = 0; iatm < Natoms; ++iatm) {
      def = false;
      if (atoms[iatm] == 1) {def = true;}                                                       //H
      else if ((atoms[iatm] > 2)&&(atoms[iatm] < 10)) {def = true;}                             //Li,Be,B,C,N,O,F
      else if ((atoms[iatm] > 10)&&(atoms[iatm] < 18)) {def = true;}                            //Na,Mg,Al,Si,P,S,Cl
      else if ((atoms[iatm] == 30)||(atoms[iatm] == 32)||(atoms[iatm] == 48)) {def = true;}     //Zn,Ge,Cd
      else if ((atoms[iatm] == 35)||(atoms[iatm] == 53)) {def = true;}                          //Br,I
      else if ((atoms[iatm] == 50)||(atoms[iatm] == 80)||(atoms[iatm] == 82)) {def = true;}     //Sn,Hg,Pb
      if (!def) {throw("ERROR: MNDO.hpp: MNDO: checkAtoms(): atom not fully specified for MNDO-theory");}
    }
  }
  //estimate parameters for heat of formation
  double calcEISOL(int atmnr) {
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
    bool dodorbs = (((atmnr > 20)&&(atmnr < 30))||((atmnr > 38)&&(atmnr < 48))||((atmnr > 71)&&(atmnr < 80)));
    if (dodorbs) {
      //correction for atoms with partially filled d orbitals
      double r016 = SlaterCondonRadialIntegral(atmnr,0,1,6);        //F0sd
      double r066 = SlaterCondonRadialIntegral(atmnr,0,6,6);        //F0dd
      double r244 = SlaterCondonRadialIntegral(atmnr,2,4,4);        //G2sd
      double r266 = SlaterCondonRadialIntegral(atmnr,2,6,6);        //F2dd
      double r466 = SlaterCondonRadialIntegral(atmnr,4,6,6);        //F4dd
      int n016 = Nr016(atmnr);
      int n066 = Nr066(atmnr);
      int n244 = Nr244(atmnr);
      int n266 = Nr266(atmnr);
      int n466 = Nr466(atmnr);
      eisol += n016*r016 + n066*r066 - n244*r244/5.0 - n266*r266/49.0 - n466*r466/49.0;
    }
    return eisol;
  }
  //parameters
  virtual double ZeroOverlap(size_t atm) {return ZeroOverlapMNDOd(atm);}
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
        enth = -2.363567808404;
        break;
      case 6:            //C
        enth = -4.405531461381;
        break;
      case 7:            //N
        enth = -7.418630559411;
        break;
      case 8:            //O
        enth = -11.653156345739;
        break;
      case 9:            //F
        enth = -17.517807628747;
        break;
      case 11:           //Na
        enth = -0.191133244110;
        break;
      case 12:           //Mg
        enth = -0.838577980508;
        break;
      case 13:           //Al
        enth = -2.001448256136;
        break;
      case 14:           //Si
        enth = -2.990566579586;
        break;
      case 15:           //P
        enth = -4.811178862983;
        break;
      case 16:           //S
        enth = -6.805235065709;
        break;
      case 17:           //Cl
        enth = -9.676042537744;
        break;
      case 30:           //Zn
        enth = -1.010065292227;
        break;
      case 32:           //Ge
        enth = -2.787397563744;
        break;
      case 35:           //Br
        enth = -9.100628685109;
        break;
      case 48:           //Cd
        enth = -0.956767450850;
        break;
      case 50:           //Sn
        enth = -3.378148784301;
        break;
      case 53:           //I
        enth = -8.438500792289;
        break;
      case 80:           //Hg
        enth = -1.077330339657;
        break;
      case 82:           //Pb
        enth = -3.874646950710;
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
      case 11:                  //Na
        if (L == 0) {betaa0 = -1.087382;}
        else if (L == 1) {betaa0 = -0.486239;}
        break;
      case 12:                  //Mg
        if (L == 0) {betaa0 = -1.895884;}
        else if (L == 1) {betaa0 = -2.141089;}
        break;
      case 13:                  //Al
        if (L == 0) {betaa0 = -5.402719;}
        else if (L == 1) {betaa0 = -3.406734;}
        else if (L == 2) {betaa0 = -1.778777;}
        break;
      case 14:                  //Si
        if (L == 0) {betaa0 = -8.210734;}
        else if (L == 1) {betaa0 = -4.884620;}
        else if (L == 2) {betaa0 = -2.608011;}
        break;
      case 15:                  //P
        if (L == 0) {betaa0 = -8.902104;}
        else if (L == 1) {betaa0 = -9.386110;}
        else if (L == 2) {betaa0 = -2.091701;}
        break;
      case 16:                  //S
        if (L == 0) {betaa0 = -10.999545;}
        else if (L == 1) {betaa0 = -12.215437;}
        else if (L == 2) {betaa0 = -1.880669;}
        break;
      case 17:                  //Cl
        if (L == 0) {betaa0 = -6.037292;}
        else if (L == 1) {betaa0 = -19.183386;}
        else if (L == 2) {betaa0 = -1.877782;}
        break;
      case 30:                  //Zn
        if (L == 0) {betaa0 = -5.017261;}
        else if (L == 1) {betaa0 = -0.712060;}
        break;
      case 32:                  //Ge
        if (L == 0) {betaa0 = -4.5164790;}
        else if (L == 1) {betaa0 = -1.7555170;}
        break;
      case 35:                  //Br
        if (L == 0) {betaa0 = -8.314976;}
        else if (L == 1) {betaa0 = -10.507041;}
        else if (L == 2) {betaa0 = -0.962599;}
        break;
      case 48:                  //Cd
        if (L == 0) {betaa0 = -2.771544;}
        else if (L == 1) {betaa0 = -1.805650;}
        break;
      case 50:                  //Sn
        if (L == 0) {betaa0 = -3.2351470;}
        else if (L == 1) {betaa0 = -4.2904160;}
        break;
      case 53:                  //I
        if (L == 0) {betaa0 = -10.699487;}
        else if (L == 1) {betaa0 = -4.941178;}
        else if (L == 2) {betaa0 = -2.350461;}
        break;
      case 80:                  //Hg
        if (L == 0) {betaa0 = -2.218722;}
        else if (L == 1) {betaa0 = -2.909786;}
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
      case 11:                //Na
        if (L == 0) {ulx = -5.201000;}
        else if (L == 1) {ulx = -2.712573;}
        break;
      case 12:                //Mg
        if (L == 0) {ulx = -15.097000;}
        else if (L == 1) {ulx = -10.650000;}
        break;
      case 13:                //Al
        if (L == 0) {ulx = -28.961830;}
        else if (L == 1) {ulx = -22.814474;}
        else if (L == 2) {ulx = -6.690879;}
        break;
      case 14:                //Si
        if (L == 0) {ulx = -36.051529;}
        else if (L == 1) {ulx = -27.535692;}
        else if (L == 2) {ulx = -14.677439;}
        break;
      case 15:                //P
        if (L == 0) {ulx = -47.055531;}
        else if (L == 1) {ulx = -38.067059;}
        else if (L == 2) {ulx = -23.691597;}
        break;
      case 16:                //S
        if (L == 0) {ulx = -56.889130;}
        else if (L == 1) {ulx = -47.274746;}
        else if (L == 2) {ulx = -25.095118;}
        break;
      case 17:                //Cl
        if (L == 0) {ulx = -69.622971;}
        else if (L == 1) {ulx = -59.100731;}
        else if (L == 2) {ulx = -36.674572;}
        break;
      case 30:                //Zn
        if (L == 0) {ulx = -18.023001;}
        else if (L == 1) {ulx = -12.242166;}
        break;
      case 32:                //Ge
        if (L == 0) {ulx = -33.9493670;}
        else if (L == 1) {ulx = -27.4251050;}
        break;
      case 35:                //Br
        if (L == 0) {ulx = -65.402779;}
        else if (L == 1) {ulx = -54.553753;}
        else if (L == 2) {ulx = -13.728099;}
        break;
      case 48:                //Cd
        if (L == 0) {ulx = -16.969700;}
        else if (L == 1) {ulx = -12.400965;}
        break;
      case 50:                //Sn
        if (L == 0) {ulx = -40.8518020;}
        else if (L == 1) {ulx = -28.5602490;}
        break;
      case 53:                //I
        if (L == 0) {ulx = -62.765354;}
        else if (L == 1) {ulx = -50.292114;}
        else if (L == 2) {ulx = -12.248305;}
        break;
      case 80:                //Hg
        if (L == 0) {ulx = -18.815649;}
        else if (L == 1) {ulx = -13.397114;}
        break;
      case 82:                //Pb
        if (L == 0) {ulx = -47.3196920;}
        else if (L == 1) {ulx = -28.8475600;}
        break;
    }
    return ulx/au2eV;
  }
  virtual double alpha(int atm1, int atm2) {
    //function that returns alpha values for MNDOd
    //alphas stored and returned in 1/Angstrom
    double alpha = 0.0;
    switch (atm1) {
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
      case 11:    //Na
        if ((atm2 == 1)||(atm2 == 6)) {alpha = 1.052252;}
        else {alpha = 1.170102;}
        break;
      case 12:    //Mg
        if (atm2 == 1) {alpha = 1.350530;}
        else if ((atm2 == 6)||(atm2 == 16)) {alpha = 1.481721;}
        else {alpha = 1.621470;}
        break;
      case 13:    //Al
        if ((atm2 == 1)||(atm2 == 6)||(atm2 == 13)) {alpha = 1.387880;}
        else {alpha = 1.439879;}
        break;
      case 14:    //Si
        alpha = 1.660069;
        break;
      case 15:    //P
        alpha = 1.852551;
        break;
      case 16:    //S
        alpha = 2.023060;
        break;
      case 17:    //Cl
        alpha = 2.180300;
        break;
      case 30:    //Zn
        alpha = 1.517637;
        break;
      case 32:    //Ge
        alpha = 1.978498;
        break;
      case 35:    //Br
        alpha = 2.091050;
        break;
      case 48:    //Cd
        alpha = 1.424613;
        break;
      case 50:    //Sn
        alpha = 1.800814;
        break;
      case 53:    //I
        alpha = 1.906174;
        break;
      case 80:    //Hg
        alpha = 1.382242;
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
        else if (idx == 3) {dval = 1.3049119659857598;}
        break;
      case 4:                    //Be
        if (idx == 1) {dval = 0.7605847;}
        else if (idx == 2) {dval = 0.6453775;}
        else if (idx == 3) {dval = 0.9127016133504423;}
        break;
      case 5:                    //B
        if (idx == 1) {dval = 0.5068929;}
        else if (idx == 2) {dval = 0.4301129;}
        else if (idx == 3) {dval = 0.6082714965316228;}
        break;
      case 6:                    //C
        if (idx == 1) {dval = 0.4272845;}
        else if (idx == 2) {dval = 0.3625629;}
        else if (idx == 3) {dval = 0.5127413703933202;}
        break;
      case 7:                    //N
        if (idx == 1) {dval = 0.3386159;}
        else if (idx == 2) {dval = 0.2873251;}
        else if (idx == 3) {dval = 0.4063390532302058;}
        break;
      case 8:                    //O
        if (idx == 1) {dval = 0.2828939;}
        else if (idx == 2) {dval = 0.2400435;}
        else if (idx == 3) {dval = 0.33947277325950603;}
        break;
      case 9:                    //F
        if (idx == 1) {dval = 0.2681377;}
        else if (idx == 2) {dval = 0.2275224;}
        else if (idx == 3) {dval = 0.3217652638236763;}
        break;
      case 11:                   //Na
        if (idx == 1) {dval = 1.722047013848492;}
        else if (idx == 2) {dval = 1.568867318336547;}
        else if (idx == 3) {dval = 2.218713439155453;}
        break;
      case 12:                    //Mg
        if (idx == 1) {dval = 1.067674252498970;}
        else if (idx == 2) {dval = 0.956997330811926;}
        else if (idx == 3) {dval = 1.353398604389077;}
        break;
      case 13:                   //Al
        if (idx == 1) {dval = 0.832308566407991;}               //ss,sp
        else if (idx == 2) {dval = 0.952665777575407;}          //pp
        else if (idx == 3) {dval = 1.347272863055851;}          //D2pp
        else if (idx == 4) {dval = 1.177532394535999;}          //D2sd
        else if (idx == 5) {dval = 0.897301564955694;}          //D1pd
        else if (idx == 6) {dval = 1.154925851527145;}          //D2dd
        break;
      case 14:                   //Si
        if (idx == 1) {dval = 0.716671969127800;}               //ss,sp
        else if (idx == 2) {dval = 0.906995869298913;}          //pp
        else if (idx == 3) {dval = 1.282685859378898;}          //D2pp
        else if (idx == 4) {dval = 1.039958482223044;}          //D2sd
        else if (idx == 5) {dval = 0.864109122473992;}          //D1pd
        else if (idx == 6) {dval = 1.128112373316393;}          //D2dd
        break;
      case 15:                   //P
        if (idx == 1) {dval = 0.679314907108863;}               //ss,sp
        else if (idx == 2) {dval = 0.817686787361311;}          //pp
        else if (idx == 3) {dval = 1.156383744459651;}          //D2pp
        else if (idx == 4) {dval = 0.973232417552050;}          //D2sd
        else if (idx == 5) {dval = 0.789138952424778;}          //D1pd
        else if (idx == 6) {dval = 1.051503868894154;}          //D2dd
        break;
      case 16:                   //S
        if (idx == 1) {dval = 0.633873496222589;}               //ss,sp
        else if (idx == 2) {dval = 0.789776177529618;}          //pp
        else if (idx == 3) {dval = 1.116912181501567;}          //D2pp
        else if (idx == 4) {dval = 0.915348224929503;}          //D2sd
        else if (idx == 5) {dval = 0.759895375983441;}          //D1pd
        else if (idx == 6) {dval = 1.007144157617617;}          //D2dd
        break;
      case 17:                   //Cl
        if (idx == 1) {dval = 0.578881496107948;}               //ss,sp
        else if (idx == 2) {dval = 0.749747116762158;}          //pp
        else if (idx == 3) {dval = 1.060302540875168;}          //D2pp
        else if (idx == 4) {dval = 0.845123889350194;}          //D2sd
        else if (idx == 5) {dval = 0.712779146193005;}          //D1pd
        else if (idx == 6) {dval = 0.927879686590404;}          //D2dd
        break;
      case 30:                   //Zn
        if (idx == 1) {dval = 0.763507018483803;}
        else if (idx == 2) {dval = 1.300999463764879;}
        else if (idx == 3) {dval = 1.839891086296417;}
        break;
      case 32:                   //Ge
        if (idx == 1) {dval = 0.6644269;}
        else if (idx == 2) {dval = 0.5555542;}
        else if (idx == 3) {dval = 0.785672284273335;}
        break;
      case 35:                   //Br
        if (idx == 1) {dval = 0.600967020923898;}               //ss,sp
        else if (idx == 2) {dval = 0.783399557210443;}          //pp
        else if (idx == 3) {dval = 1.107894278564086;}          //D2pp
        else if (idx == 4) {dval = 0.856496644208800;}          //D2sd
        else if (idx == 5) {dval = 0.778100505006181;}          //D1pd
        else if (idx == 6) {dval = 1.079778579131701;}          //D2dd
        break;
      case 48:                   //Cd
        if (idx == 1) {dval = 0.992909075368672;}
        else if (idx == 2) {dval = 0.891159574668322;}
        else if (idx == 3) {dval = 1.260289956734580;}
        break;
      case 50:                   //Sn
        if (idx == 1) {dval = 0.8306740;}
        else if (idx == 2) {dval = 0.7017967;}
        else if (idx == 3) {dval = 0.9924904111686824;}
        break;
      case 53:                   //I
        if (idx == 1) {dval = 0.540129610588824;}               //ss,sp
        else if (idx == 2) {dval = 0.864645250800544;}          //pp
        else if (idx == 3) {dval = 1.222793040323616;}          //D2pp
        else if (idx == 4) {dval = 0.759451393810603;}          //D2sd
        else if (idx == 5) {dval = 0.859183859369850;}          //D1pd
        else if (idx == 6) {dval = 1.290775710431637;}          //D2dd
        break;
      case 80:                   //Hg
        if (idx == 1) {dval = 0.923531152348794;}
        else if (idx == 2) {dval = 0.936259354817651;}
        else if (idx == 3) {dval = 1.324070677481806;}
        break;
      case 82:                   //Pb
        if (idx == 1) {dval = 0.8216177;}
        else if (idx == 2) {dval = 0.7666867;}
        else if (idx == 3) {dval = 1.0842587292310726;}
        break;
    }
    return dval*dist_Angstrom2aum1;
  }
  virtual double rhocore(size_t atmnr) {
    //function returning rho values used to calculate VAB integrals; values already in atomic units
    double rho_ = 0.0;
    if (atmnr == 11) {rho_ = 1.530553;}                //Na
    else if (atmnr == 12) {rho_ = 1.350776;}           //Mg
    else if (atmnr == 29) {rho_ = 1.589234;}           //Zn
    else if (atmnr == 48) {rho_ = 1.721186;}           //Cd
    else if (atmnr == 80) {rho_ = 1.636072;}           //Hg
    else {rho_ = rho(atmnr,0);}                        //all other atoms
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
      case 11:                   //Na
        if (l == 0) {rho = 1.567071093635136;}
        else if (l == 1) {rho = 1.830222533372990;}
        else if (l == 2) {rho = 1.417020172941393;}
        break;
      case 12:                   //Mg
        if (l == 0) {rho = 0.976229438953369;}
        else if (l == 1) {rho = 1.222534313856169;}
        else if (l == 2) {rho = 0.864372080385308;}
        break;
      case 13:                   //Al
        if (l == 0) {rho = 0.760157668073853;}
        else if (l == 1) {rho = 1.029192272831372;}
        else if (l == 2) {rho = 0.860459625876539;}
        else if (l == 3) {rho = 0.306374315792088;}
        else if (l == 4) {rho = 0.885024468142918;}
        else if (l == 5) {rho = 1.118515830987198;}
        else if (l == 6) {rho = 0.972795025314052;}
        break;
      case 14:                   //Si
        if (l == 0) {rho = 0.670271695839240;}
        else if (l == 1) {rho = 0.897134091517244;}
        else if (l == 2) {rho = 0.819209949824937;}
        else if (l == 3) {rho = 2.645886041500000;}
        else if (l == 4) {rho = 0.853697951546307;}
        else if (l == 5) {rho = 1.092547670500620;}
        else if (l == 6) {rho = 0.950210003881803;}
        break;
      case 15:                   //P
        if (l == 0) {rho = 0.627175684946922;}
        else if (l == 1) {rho = 0.589029003839085;}
        else if (l == 2) {rho = 0.738545109729496;}
        else if (l == 3) {rho = 0.230274672062604;}
        else if (l == 4) {rho = 0.783143597305115;}
        else if (l == 5) {rho = 1.018354314389296;}
        else if (l == 6) {rho = 0.885682597109276;}
        break;
      case 16:                   //S
        if (l == 0) {rho = 0.590328325842112;}
        else if (l == 1) {rho = 0.508801745643787;}
        else if (l == 2) {rho = 0.713335789594268;}
        else if (l == 3) {rho = 2.645886041500000;}
        else if (l == 4) {rho = 0.753084041637511;}
        else if (l == 5) {rho = 0.975393080769647;}
        else if (l == 6) {rho = 0.848318378715980;}
        break;
      case 17:                   //Cl
        if (l == 0) {rho = 0.544980795824593;}
        else if (l == 1) {rho = 0.450625897868114;}
        else if (l == 2) {rho = 0.677180931836914;}
        else if (l == 3) {rho = 1.647969914212899;}
        else if (l == 4) {rho = 0.703877431305918;}
        else if (l == 5) {rho = 0.898627489661325;}
        else if (l == 6) {rho = 0.781553842697954;}
        break;
      case 30:                   //Zn
        if (l == 0) {rho = 0.841029168406762;}
        else if (l == 1) {rho = 1.140966140980049;}
        else if (l == 2) {rho = 1.155784659246568;}
        break;
      case 32:                   //Ge
        if (l == 0) {rho = 0.734624198;}
        else if (l == 1) {rho = 0.726135254;}
        else if (l == 2) {rho = 0.608610513;}
        break;
      case 35:                   //Br
        if (l == 0) {rho = 0.589069893791056;}
        else if (l == 1) {rho = 0.485997257226682;}
        else if (l == 2) {rho = 0.695958153863445;}
        else if (l == 3) {rho = 0.364713943660608;}
        else if (l == 4) {rho = 0.778171834446568;}
        else if (l == 5) {rho = 1.059880891099237;}
        else if (l == 6) {rho = 0.896749329209849;}
        break;
      case 50:                   //Sn
        if (l == 0) {rho = 0.734624198;}
        else if (l == 1) {rho = 0.821688399;}
        else if (l == 2) {rho = 0.712428177;}
        break;
      case 53:                   //I
        if (l == 0) {rho = 0.600947580115929;}
        else if (l == 1) {rho = 0.525723746001474;}
        else if (l == 2) {rho = 0.757767568929242;}
        else if (l == 3) {rho = 0.580249398341799;}
        else if (l == 4) {rho = 0.901076613080184;}
        else if (l == 5) {rho = 1.275688097573650;}
        else if (l == 6) {rho = 1.059659488896923;}
        break;
      case 80:                   //Hg
        if (l == 0) {rho = 0.865815999544531;}
        else if (l == 1) {rho = 1.059780622389568;}
        else if (l == 2) {rho = 0.811246212763297;}
        break;
      case 82:                   //Pb
        if (l == 0) {rho = 0.734624198;}
        else if (l == 1) {rho = 0.816789939;}
        else if (l == 2) {rho = 0.755508834;}
        break;
    }
    return rho*dist_Angstrom2aum1;
  }
  virtual double PrincipalQuantNumb(size_t atmnr, size_t l) {
    //function returning principal quantum numbers for atoms for the different orbitals
    double pqn = 0.0;
    if ((atmnr < 3)&&(l == 0)) {pqn = 1.0;}
    else if ((atmnr == 2)&&(l == 1)) {pqn = 2.0;}
    else if ((atmnr > 2)&&(atmnr < 10)) {pqn = 2.0;}
    else if ((atmnr == 10)&&(l == 0)) {pqn = 3.0;}
    else if ((atmnr == 10)&&(l == 1)) {pqn = 2.0;}
    else if ((atmnr > 10)&&(atmnr < 18)) {pqn = 3.0;}
    else if ((atmnr == 18)&&(l == 0)) {pqn = 4.0;}
    else if ((atmnr == 18)&&(l == 1)) {pqn = 3.0;}
    else if ((atmnr > 18)&&(atmnr < 36)&&(l < 2)) {pqn = 4.0;}
    else if ((atmnr > 18)&&(atmnr < 30)&&(l == 2)) {pqn = 3.0;}
    else if ((atmnr > 30)&&(atmnr < 36)&&(l == 2)) {pqn = 4.0;}
    else if ((atmnr == 36)&&(l == 0)) {pqn = 5.0;}
    else if ((atmnr == 36)&&(l == 1)) {pqn = 4.0;}
    else if ((atmnr > 36)&&(atmnr < 54)&&(l < 2)) {pqn = 5.0;}
    else if ((atmnr > 36)&&(atmnr < 48)&&(l == 2)) {pqn = 4.0;}
    else if ((atmnr > 48)&&(atmnr < 54)&&(l == 2)) {pqn = 5.0;}
    else if ((atmnr == 54)&&(l == 0)) {pqn = 6.0;}
    else if ((atmnr == 54)&&(l == 1)) {pqn = 5.0;}
    else if ((atmnr > 54)&&(atmnr < 84)&&(l < 2)) {pqn = 6.0;}
    else if ((atmnr > 54)&&(atmnr < 84)&&(l == 2)) {pqn = 5.0;}
    return pqn;
  }
  virtual double zn(size_t atmnr, size_t l) {
    //function returning the zn values needed to calculate 1-center eris; values stored and returned directly in atomic units
    double zn_ = 0.0;
    if (atmnr == 11) {                                 //Na
      if (l == 0) {zn_ = 0.654113;}                        //zs
      else if (l == 1) {zn_ = 0.564409;}                   //zp
    }
    else if (atmnr == 12) {                            //Mg
      if (l == 0) {zn_ = 1.050000;}                        //zs
      else if (l == 1) {zn_ = 0.925272;}                   //zp
    }
    else if (atmnr == 13) {                            //Al
      if (l == 0) {zn_ = 1.348458;}                        //zs
      else if (l == 1) {zn_ = 0.929479;}                   //zp
      else if (l == 2) {zn_ = 0.916383;}                   //zd
    }
    else if (atmnr == 14) {                            //Si
      if (l == 0) {zn_ = 1.529292;}                        //zs
      else if (l == 1) {zn_ = 0.976281;}                   //zp
      else if (l == 2) {zn_ = 0.938164;}                   //zd
    }
    else if (atmnr == 15) {                            //P
      if (l == 0) {zn_ = 1.634376;}                        //zs
      else if (l == 1) {zn_ = 1.082912;}                   //zp
      else if (l == 2) {zn_ = 1.006515;}                   //zd
    }
    else if (atmnr == 16) {                            //S
      if (l == 0) {zn_ = 1.736391;}                        //zs
      else if (l == 1) {zn_ = 1.121182;}                   //zp
      else if (l == 2) {zn_ = 1.050847;}                   //zd
    }
    else if (atmnr == 17) {                            //Cl
      if (l == 0) {zn_ = 1.880875;}                        //zs
      else if (l == 1) {zn_ = 1.181042;}                   //zp
      else if (l == 2) {zn_ = 1.140616;}                   //zd
    }
    else if (atmnr == 30) {                            //Zn
      if (l == 0) {zn_ = 1.566000;}                        //zs
      else if (l == 1) {zn_ = 0.862840;}                   //zp
    }
    else if (atmnr == 35) {                            //Br
      if (l == 0) {zn_ = 2.235816;}                        //zs
      else if (l == 1) {zn_ = 1.432927;}                   //zp
      else if (l == 2) {zn_ = 1.242578;}                   //zd
    }
    else if (atmnr == 48) {                            //Cd
      if (l == 0) {zn_ = 1.763148;}                        //zs
      else if (l == 1) {zn_ = 1.525519;}                   //zp
    }
    else if (atmnr == 53) {                            //I
      if (l == 0) {zn_ = 2.672411;}                        //zs
      else if (l == 1) {zn_ = 1.572299;}                   //zp
      else if (l == 2) {zn_ = 1.258848;}                   //zd
    }
    else if (atmnr == 80) {                            //Hg
      if (l == 0) {zn_ = 2.186000;}                        //zs
      else if (l == 1) {zn_ = 1.705005;}                   //zp
    }
    return zn_;
  }
  virtual double F0sd(size_t atmnr) {
    //function returning the F0sd 1-center eris; values stored in eV, returned in a.u.
    return 0.0;
  }
  virtual double G2sd(size_t atmnr) {
    //function returning the G2sd 1-center eris; values stored in eV, returned in a.u.
    return 0.0;
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
      case 11:         //Na
        if ((Lbra == 0)&&(Lket == 0)) {eri = 4.594445;}                                          //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 4.299198;}                                     //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 4.147574;}       //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 3.796957;}     //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 0.534409;}                                     //(sp|sp)||(ps|ps)
        break;
      case 12:         //Mg
        if ((Lbra == 0)&&(Lket == 0)) {eri = 7.375133;}                                          //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 7.047954;}                                     //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 6.888907;}       //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 6.224599;}     //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 0.726734;}                                     //(sp|sp)||(ps|ps)
        break;
      case 13:         //Al
        if ((Lbra == 0)&&(Lket == 0)) {eri = 9.471485;}                                          //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 7.080002;}                                     //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 8.791266;}       //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 6.252903;}     //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 0.778059;}                                     //(sp|sp)||(ps|ps)
        break;
      case 14:         //Si
        if ((Lbra == 0)&&(Lket == 0)) {eri = 10.741647;}                                         //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 7.436497;}                                     //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 7.560667;}       //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 6.567751;}     //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 0.877539;}                                     //(sp|sp)||(ps|ps)
        break;
      case 15:         //P
        if ((Lbra == 0)&&(Lket == 0)) {eri = 11.4797530;}                                        //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 8.248723;}                                     //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 8.557570;}       //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 7.285092;}     //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 2.107804;}                                     //(sp|sp)||(ps|ps)
        break;
      case 16:         //S
        if ((Lbra == 0)&&(Lket == 0)) {eri = 12.196301;}                                         //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 8.540233;}                                     //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 8.853901;}       //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 7.542547;}     //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 2.646352;}                                     //(sp|sp)||(ps|ps)
        break;
      case 17:         //Cl
        if ((Lbra == 0)&&(Lket == 0)) {eri = 13.211148;}                                         //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 8.996201;}                                     //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 9.419496;}       //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 7.945248;}     //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 3.081499;}                                     //(sp|sp)||(ps|ps)
        break;
      case 30:         //Zn
        if ((Lbra == 0)&&(Lket == 0)) {eri = 8.560728;}                                          //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 5.139648;}                                     //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 7.490036;}       //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 4.505403;}     //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 0.532946;}                                     //(sp|sp)||(ps|ps)
        break;
      case 32:         //Ge
        if ((Lbra == 0)&&(Lket == 0)) {eri = 9.8;}                                               //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 7.3;}                                          //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 8.3;}            //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 6.5;}          //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 1.3;}                                          //(sp|sp)||(ps|ps)
        break;
      case 35:         //Br
        if ((Lbra == 0)&&(Lket == 0)) {eri = 12.222356;}                                         //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 8.535464;}                                     //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 8.263721;}       //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 7.482167;}     //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 2.749522;}                                     //(sp|sp)||(ps|ps)
        break;
      case 48:         //Cd
        if ((Lbra == 0)&&(Lket == 0)) {eri = 7.904434;}                                          //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 7.480000;}                                     //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 7.515707;}       //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 6.518664;}     //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 0.636744;}                                     //(sp|sp)||(ps|ps)
        break;
      case 50:         //Sn
        if ((Lbra == 0)&&(Lket == 0)) {eri = 9.8;}                                               //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 7.3;}                                          //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 8.3;}            //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 6.5;}          //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 1.3;}                                          //(sp|sp)||(ps|ps)
        break;
      case 53:         //I
        if ((Lbra == 0)&&(Lket == 0)) {eri = 11.980782;}                                         //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 7.709372;}                                     //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 7.855902;}       //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 6.718557;}     //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 2.071475;}                                     //(sp|sp)||(ps|ps)
        break;
      case 80:         //Hg
        if ((Lbra == 0)&&(Lket == 0)) {eri = 8.315649;}                                          //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 7.115259;}                                     //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 8.212173;}       //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 6.171250;}     //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 0.835941;}                                     //(sp|sp)||(ps|ps)
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
};

#endif //_MNDOD_
