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

#ifndef _MNDOPDDG_
#define _MNDOPDDG_
#include "MNDO.hpp"

//descrition:
//The MNDO-PDDG method

class MNDOPDDG: public MNDO {
  //this is the implementation of Jorgensen's MNDOPDDG
  //M. P. Repasky, J. Chandrasekhar, W. L. Jorgensen, J. Comput. Chem., 23, 1601, 2002
  //I. Tubert-Brohman, C. R. W. Guimaraes, M. P. Repasky, W. L. Jorgensen, J. Comput. Chem., 25, 138, 2004
  //I. Tubert-Brohman, C. R. W. Guimaraes, W. L. Jorgensen, J. Chem. Theory Comput., 1, 817, 2005
  //bX quantities are "barred" tensors, which are used only for the open-shell case
public:
  MNDOPDDG(BSet _bset, Molecule _mol, std::string _openclosed = "0", std::string _corecorrection = "0"): MNDO(_bset,_mol,_openclosed,_corecorrection) {doPDDG = true;}
  ~MNDOPDDG() {}
  std::string Type() {return "MNDOPDDG";}
  double gfactor(size_t iatm1, size_t iatm2, double RAB) {return 1.0;}
  virtual double gfactor_dR(size_t iatm1, size_t iatm2) {return 0.0;}
  //other auxiliary functions
  void checkAtoms() {
    //function that checks whether the atoms in the molecule are fully defined in the respective MNDO theories
    bool def;
    for (size_t iatm = 0; iatm < Natoms; ++iatm) {
      def = false;
      if (atoms[iatm] == 1) {def = true;}                                                        //H
      else if ((atoms[iatm] > 5)&&(atoms[iatm] < 10)) {def = true;}                              //C,N,O,F
      else if ((atoms[iatm] == 17)||(atoms[iatm] == 35)||(atoms[iatm] == 53)) {def = true;}      //Cl,Br,I
      if (!def) {throw("ERROR: MNDOPDDG.hpp: MNDOPDDG: checkAtoms(): atom not fully specified for MNDOPDDG-theory");}
    }
  }
  double ZeroOverlap(size_t atm) {return ZeroOverlapMNDOPDDG(atm);}
  double ElementFormTheo(size_t atomicnr) {
    //function returning theoretical heats of formation for elements; values directly in a.u.
    double enth = 0.0;
    switch (atomicnr) {
      case 1:      //H
        enth = -12.015956;
        break;
      case 6:      //C
        enth = -123.864412;
        break;
      case 7:      //N
        enth = -206.466626;
        break;
      case 8:      //O
        enth = -310.879745;
        break;
      case 9:      //F
        enth = -488.703243;
        break;
      case 17:     //Cl
        enth = -378.909727;
        break;
      case 35:     //Br
        enth = -349.564096;
        break;
      case 53:     //I
        enth = -356.076398;
        break;
    }
    return enth/au2eV;
  }
  double betaA0(size_t atomicnr, int L) {
    //function returning the bonding parameters; values stored in eV, returned however in a.u.
    double betaa0 = 0.0;
    switch (atomicnr) {
      case 1:         //H
        betaa0 = -7.493504;
        break;
      case 6:         //C
        if (L == 0) {betaa0 = -18.841334;}
        else if (L == 1) {betaa0 = -7.922234;}
        break;
      case 7:         //N
        if (L == 0) {betaa0 = -20.375774;}
        else if (L == 1) {betaa0 = -21.085373;}
        break;
      case 8:         //O
        if (L == 0) {betaa0 = -33.606336;}
        else if (L == 1) {betaa0 = -27.984442;}
        break;
      case 9:         //F
        if (L == 0) {betaa0 = -67.827612;}
        else if (L == 1) {betaa0 = -40.924818;}
        break;
      case 17:        //Cl
        if (L == 0) {betaa0 = -15.663317;}
        else if (L == 1) {betaa0 = -15.399331;}
        break;
      case 35:        //Br
        if (L == 0) {betaa0 = -7.054170;}
        else if (L == 1) {betaa0 = -10.221030;}
        break;
      case 53:        //I
        if (L == 0) {betaa0 = -6.698375;}
        else if (L == 1) {betaa0 = -5.693814;}
        break;
    }
    return betaa0/au2eV;
  }
  double UlX(size_t atomicnr, int L) {
    //function that returns the atomic potential U
    double ulx = 0.0;
    switch (atomicnr) {
      case 1:                 //H
        ulx = -11.724114;
        break;
      case 6:                 //C
        if (L == 0) {ulx = -53.837582;}
        else if (L == 1) {ulx = -39.936409;}
        break;
      case 7:                 //N
        if (L == 0) {ulx = -71.871894;}
        else if (L == 1) {ulx = -58.216617;}
        break;
      case 8:                 //O
        if (L == 0) {ulx = -97.884970;}
        else if (L == 1) {ulx = -77.342674;}
        break;
      case 9:                 //F
        if (L == 0) {ulx = -134.220379;}
        else if (L == 1) {ulx = -107.155961;}
        break;
      case 17:                //Cl
        if (L == 0) {ulx = -111.133653;}
        else if (L == 1) {ulx = -78.062493;}
        break;
      case 35:                //Br
        if (L == 0) {ulx = -100.637007;}
        else if (L == 1) {ulx = -76.015735;}
        break;
      case 53:                //I
        if (L == 0) {ulx = -106.588422;}
        else if (L == 1) {ulx = -75.282605;}
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
        alpha = 2.491813;
        break;
      case 6:     //C
        alpha = 2.555522;
        break;
      case 7:     //N
        alpha = 2.843678;
        break;
      case 8:     //O
        alpha = 3.238842;
        break;
      case 9:     //F
        alpha = 3.322382;
        break;
      case 17:    //Cl
        alpha = 2.602846;
        break;
      case 35:    //Br
        alpha = 2.414265;
        break;
      case 53:    //I
        alpha = 2.242446;
        break;
    }
    return alpha;
  }
  double Dvalue(size_t atmnr, size_t idx) {
    //function returning the D values needed to calculate eris; values stored and returned in a.u.
    //note that this function returns both D1 and D2; idx is then either 1 or 2
    double dval = 0.0;
    switch (atmnr) {
      case 1:                  //H
        dval = 0.0;
      case 6:                  //C
        if (idx == 1) {dval = 0.794158;}
        else if (idx == 2) {dval = 0.671090;}
        break;
      case 7:                  //N
        if (idx == 1) {dval = 0.643624;}
        else if (idx == 2) {dval = 0.543495;}
        break;
      case 8:                  //O
        if (idx == 1) {dval = 0.547344;}
        else if (idx == 2) {dval = 0.454088;}
        break;
      case 9:                  //F
        if (idx == 1) {dval = 0.361556;}
        else if (idx == 2) {dval = 0.421593;}
        break;
      case 17:                 //Cl
        if (idx == 1) {dval = 0.411609;}
        else if (idx == 2) {dval = 0.821202;}
        break;
      case 35:                 //Br
        if (idx == 1) {dval = 0.574623;}
        else if (idx == 2) {dval = 0.944892;}
        break;
      case 53:                 //I
        if (idx == 1) {dval = 1.209529;}
        else if (idx == 2) {dval = 1.043559;}
        break;
    }
    return dval;
  }
  double rho(size_t atmnr, size_t l) {
    //function returning the rho values needed to calculate eris; values stored and returned in atomic units
    double rho = 0.0;
    switch (atmnr) {
      case 1:                  //H
        if (l == 0) {rho = 1.058920;}
        else if (l == 1) {rho = 0.0;}
        else if (l == 2) {rho = 0.0;}
        break;
      case 6:                  //C
        if (l == 0) {rho = 1.112429;}
        else if (l == 1) {rho = 0.805697;}
        else if (l == 2) {rho = 0.737353;}
        break;
      case 7:                  //N
        if (l == 0) {rho = 1.001103;}
        else if (l == 1) {rho = 0.639479;}
        else if (l == 2) {rho = 0.615679;}
        break;
      case 8:                  //O
        if (l == 0) {rho = 0.882296;}
        else if (l == 1) {rho = 0.527927;}
        else if (l == 2) {rho = 0.526913;}
        break;
      case 9:                  //F
        if (l == 0) {rho = 0.804078;}
        else if (l == 1) {rho = 0.383553;}
        else if (l == 2) {rho = 0.476913;}
        break;
      case 17:                 //Cl
        if (l == 0) {rho = 0.905190;}
        else if (l == 1) {rho = 0.554248;}
        else if (l == 2) {rho = 0.825604;}
        break;
      case 35:                 //Br
        if (l == 0) {rho = 0.904802;}
        else if (l == 1) {rho = 0.668867;}
        else if (l == 2) {rho = 0.884978;}
        break;
      case 53:                 //I
        if (l == 0) {rho = 0.904562;}
        else if (l == 1) {rho = 1.002309;}
        else if (l == 2) {rho = 0.992001;}
        break;
    }
    return rho;
  }
  double PAPDDG(int atmnr, int index) {
    //function returning the PDDG PAi parameters in Hartree
    double pa = 0.0;
    switch (atmnr) {
      case 1:      //H
        if (index == 1) {pa = -0.108861;}
        else if (index == 2) {pa = -0.024706;}
        break;
      case 6:      //C
        if (index == 1) {pa = -0.006889;}
        else if (index == 2) {pa = -0.027751;}
        break;
      case 7:      //N
        if (index == 1) {pa = 0.035027;}
        else if (index == 2) {pa = -0.001721;}
        break;
      case 8:      //O
        if (index == 1) {pa = 0.086344;}
        else if (index == 2) {pa = 0.030403;}
        break;
      case 9:      //F
        if (index == 1) {pa = -0.011579;}
        else if (index == 2) {pa = -0.012943;}
        break;
      case 17:     //Cl
        if (index == 1) {pa = -0.017119;}
        else if (index == 2) {pa = 0.005497;}
        break;
      case 35:     //Br
        if (index == 1) {pa = -0.017133;}
        else if (index == 2) {pa = -0.016964;}
        break;
      case 53:     //I
        if (index == 1) {pa = 0.009616;}
        else if (index == 2) {pa = -0.007505;}
        break;
    }
    return pa/au2eV;
  }
  double DAPDDG(int atmnr, int index) {
    //function returning the PDDG DAi parameters in inverse Angstroem
    double da = 0.0;
    switch (atmnr) {
      case 1:      //H
        if (index == 1) {da = 0.460721;}
        else if (index == 2) {da = 1.298731;}
        break;
      case 6:      //C
        if (index == 1) {da = 1.192456;}
        else if (index == 2) {da = 1.329522;}
        break;
      case 7:      //N
        if (index == 1) {da = 1.011630;}
        else if (index == 2) {da = 2.278423;}
        break;
      case 8:      //O
        if (index == 1) {da = 0.725408;}
        else if (index == 2) {da = 0.709728;}
        break;
      case 9:      //F
        if (index == 1) {da = 0.834606;}
        else if (index == 2) {da = 1.875603;}
        break;
      case 17:     //Cl
        if (index == 1) {da = 1.466335;}
        else if (index == 2) {da = 2.236842;}
        break;
      case 35:     //Br
        if (index == 1) {da = 2.201539;}
        else if (index == 2) {da = 2.255764;}
        break;
      case 53:     //I
        if (index == 1) {da = 2.572332;}
        else if (index == 2) {da = 2.936456;}
        break;
    }
    return da;
  }
};

#endif //_MNDOPDDG_
