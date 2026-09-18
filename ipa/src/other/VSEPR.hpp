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

#ifndef _VSEPR_Package_
#define _VSEPR_Package_
#include <vector>
#include <fstream>
#include <iostream>
#include "../nonbonded/DFTDfunc.hpp"
#include "../Geometry.hpp"

//description:
//functions with VSEPR-related information

//int VSEPRclassifier(int numberneighbors, int numberelectronpairs) {
//  //VSEPR classification in the form of xy
//  //where coordination number = x and electron pairs = y
//  //values accepted:
//  //10 -> protons
//  //20 -> linear (CO2)                                -> angle = 180
//  //30 -> trigonal planar (BCl3)                      -> angle = 120
//  //31 -> angled (SO2,O3)                             -> angle = 115
//  //32 -> linear (O2)                                 -> angle = 120
//  //40 -> tetrahedral (CH4, phosphate, sulfate)       -> angle = 109.5
//  //41 -> trigonal pyramidal (NH3)                    -> angle = 107
//  //42 -> angled (H2O)                                -> angle = 104.5
//  //43 -> linear (HCl)                                -> angle = 180
//  //50 -> trigonal bipyramidal (PCl5)                 -> angle = 90,120,180
//  //51 -> seesaw (SF4)                                -> angle = 90,120,180
//  //52 -> t-shape (ClF3)                              -> angle = 90,180
//  //53 -> linear (I3-)                                -> angle = 180
//  //60 -> octahedral (SF6)                            -> angle = 90
//  //61 -> square pyramidal (ClF5)                     -> angle = 90,180
//  //62 -> square planar (ICl4-)                       -> angle = 90,180
//  //70 -> pentagonal bipyramidal (IF7)                -> angle = 72,90,144,180
//  //this is documentation for this and the next function
//  return 10*(numberneighbors + numberelectronpairs) + numberelectronpairs;
//}
//void VSEPRangle(std::vector<double> & angles, int vseprcode) {
//  //function returning the angles accepted
//  angles.resize(4);
//  if ((vseprcode < 50)||(vseprcode == 53)||(vseprcode == 60)) {       //only one angle
//    angles[0] = 177.5*double((vseprcode == 20)||(vseprcode == 43)||(vseprcode == 53)) + 120.0*double((vseprcode == 30)||(vseprcode == 32)) + 115.0*double(vseprcode == 31);
//    angles[0] += 109.5*double(vseprcode == 40) + 107.0*double(vseprcode == 41) + 104.5*double(vseprcode == 42);
//    angles[1] = 0.0;
//    angles[2] = 0.0;
//    angles[3] = 0.0;
//  }
//  else if ((vseprcode > 51)&&(vseprcode < 63)) {
//    angles[0] = 90.0;
//    angles[1] = 180.0;
//    angles[2] = 0.0;
//    angles[3] = 0.0;
//  }
//  else if ((vseprcode > 49)&&(vseprcode < 52)) {
//    angles[0] = 90.0;
//    angles[1] = 120.0;
//    angles[2] = 180.0;
//    angles[3] = 0.0;
//  }
//  else if (vseprcode > 69) {
//    angles[0] = 72.0;
//    angles[1] = 90.0;
//    angles[2] = 144.0;
//    angles[3] = 180.0;
//  }
//}
//void HeteroAtomAngleFix(int atomnumber, double & angle, int szring) {
//  //function that fixes the angles in rings with heteroatoms
//  if (atomnumber == 16) {
//    if (szring == 5) {angle = 91.5;}
//    else if (szring == 6) {angle = 110.0;}
//  }
//}
//void OneAtomDihedral(std::vector<double> & angles, int vsepr) {
//  //function that returns dihedrals for one atom cases, like methane, ammonia, sulfate, ...
//  angles.resize(4);
//  angles[0] = 0.0;
//  angles[1] = 0.0;
//  if ((vsepr == 30)||(vsepr == 31)||(vsepr == 32)) {angles[0] = -180.0;}
//  else if ((vsepr == 40)||(vsepr == 41)||(vsepr == 42)) {
//    angles[0] = 120.0;
//    angles[1] = -120.0;
//  }
//  angles[2] = 0.0;
//  angles[3] = 0.0;
//}
//void SelectDihedrals(std::vector<double> & angles, int maxlevel, int choice) {
//  //function that selects a set of dihedral angles based on a certain choice outside this function
//  angles.resize(4);
//  angles[3] = 0.0;
//  if (maxlevel == 2) {
//    angles[2] = 0.0;
//    if (choice == 1) {
//      angles[0] = 0.0;
//      angles[1] = 180.0;
//    }
//    else if (choice == 2) {
//      angles[0] = 180.0;
//      angles[1] = 0.0;
//    }
//  }
//  else if (maxlevel == 3) {
//    if (choice == 1) {
//      angles[0] = 60.0;
//      angles[1] = 180.0;
//      angles[2] = 300.0;
//    }
//    else if (choice == 2) {
//      angles[0] = 60.0;
//      angles[1] = 300.0;
//      angles[2] = 180.0;
//    }
//    else if (choice == 3) {
//      angles[0] = 180.0;
//      angles[1] = 60.0;
//      angles[2] = 300.0;
//    }
//    else if (choice == 4) {
//      angles[0] = 300.0;
//      angles[1] = 60.0;
//      angles[2] = 180.0;
//    }
//    else if (choice == 5) {
//      angles[0] = 300.0;
//      angles[1] = 180.0;
//      angles[2] = 60.0;
//    }
//    else if (choice == 6) {
//      angles[0] = 180.0;
//      angles[1] = 300.0;
//      angles[2] = 60.0;
//    }
//  }
//}
int CountNeighbours(matrixE & Ncoord, std::vector<int> & friendlyneighbours, int idAtm, int Natoms, double threshbond = 0.3) {
  //function that counts the number of neighbours
  int nneigh = 0;
  int iaux;
  for (size_t idBtm = 0; idBtm < Natoms; ++idBtm) {
    iaux = (Ncoord(idAtm + 1,idBtm + 1) >= threshbond);
    friendlyneighbours[nneigh] = idBtm*(iaux > 0);
    nneigh += iaux;
  }
  return nneigh;
}
void GetAngles(std::vector<int> & angles, std::vector<int> & friendlyneighbours, matrixE & geom, matrix<int> & islinear, matrix<int> & istrigonal, matrix<int> & istetrahedral, matrix<int> & isorthogonal, matrix<int> & ispentagonal, matrixE & theta, int idAtm, int nneigh, double threshang = 6.0) {
  //function that gets and classifies angles
  double onedpi = 180.0/pi;
  double auxie1;
  //limits for geometry determination
  double linearmin = 180.0 - threshang;
  double linearmax = 180.0 + threshang;
  double trigonalmin = 120.0 - threshang;
  double trigonalmax = 120.0 + threshang;
  double tetrahedralmin = 109.5 - threshang;
  double tetrahedralmax = 109.5 + threshang;
  double orthogonalmin = 90.0 - threshang;
  double orthogonalmax = 90.0 + threshang;
  double pentagonalmin = 72.0 - threshang;
  double pentagonalmax = 72.0 + threshang;
  angles[0] = 0;      //nlinear
  angles[1] = 0;      //ntrigonal
  angles[2] = 0;      //ntetrahedral
  angles[3] = 0;      //northogonal
  angles[4] = 0;      //npentagonal
  for (size_t idn1 = 0; idn1 < nneigh; ++idn1) {
    for (size_t idn2 = idn1 + 1; idn2 < nneigh; ++idn2) {
      auxie1 = fabs(onedpi*Angle(friendlyneighbours[idn1] + 1,idAtm + 1,friendlyneighbours[idn2] + 1,geom,1.0e-20));
      if (auxie1 < threshang) {auxie1 += 180.0;}
      theta(idn1 + 1,idn2 + 1) = auxie1;
      theta(idn2 + 1,idn1 + 1) = auxie1;
      //linear?
      if ((linearmin < auxie1)&&(auxie1 < linearmax)) {
        islinear(idn1 + 1,idn2 + 1) = 1;
        islinear(idn2 + 1,idn1 + 1) = 1;
        ++angles[0];
      }
      else {
        islinear(idn1 + 1,idn2 + 1) = 0;
        islinear(idn2 + 1,idn1 + 1) = 0;
      }
      //trigonal?
      if ((trigonalmin < auxie1)&&(auxie1 < trigonalmax)) {
        istrigonal(idn1 + 1,idn2 + 1) = 1;
        istrigonal(idn2 + 1,idn1 + 1) = 1;
        ++angles[1];
      }
      else {
        istrigonal(idn1 + 1,idn2 + 1) = 0;
        istrigonal(idn2 + 1,idn1 + 1) = 0;
      }
      //tetrahedral?
      if ((tetrahedralmin < auxie1)&&(auxie1 < tetrahedralmax)) {
        istetrahedral(idn1 + 1,idn2 + 1) = 1;
        istetrahedral(idn2 + 1,idn1 + 1) = 1;
        ++angles[2];
      }
      else {
        istetrahedral(idn1 + 1,idn2 + 1) = 0;
        istetrahedral(idn2 + 1,idn1 + 1) = 0;
      }
      //straight angle?
      if ((orthogonalmin < auxie1)&&(auxie1 < orthogonalmax)) {
        isorthogonal(idn1 + 1,idn2 + 1) = 1;
        isorthogonal(idn2 + 1,idn1 + 1) = 1;
        ++angles[3];
      }
      else {
        isorthogonal(idn1 + 1,idn2 + 1) = 0;
        isorthogonal(idn2 + 1,idn1 + 1) = 0;
      }
      //pentagonal?
      if ((orthogonalmin < auxie1)&&(auxie1 < orthogonalmax)) {
        ispentagonal(idn1 + 1,idn2 + 1) = 1;
        ispentagonal(idn2 + 1,idn1 + 1) = 1;
        ++angles[4];
      }
      else {
        ispentagonal(idn1 + 1,idn2 + 1) = 0;
        ispentagonal(idn2 + 1,idn1 + 1) = 0;
      }
    }
  }
}
void ClassifyVSEPR(matrixE & modVSEPR, std::vector<int> & angles, int idAtm, int atomA, int nneigh, double threshang = 6.0) {
  //function that gives the VSEPR classification of atoms
  if (nneigh < 2) {                //terminal atom
    modVSEPR(idAtm + 1,2) = -1.0;
  }
  else if (nneigh == 2) {
    //check angle
    if (angles[0] == 1) {                                         //AX2E0 or AX2E3
      if ((atomA == 17)||(atomA == 18)||(atomA == 35)||(atomA == 36)||(atomA == 53)||(atomA == 54)||(atomA == 85)||(atomA == 86)||(atomA == 117)||(atomA == 118)) {
        //assume that only "late" halogens and noble gases can take this geometry
        modVSEPR(idAtm + 1,2) = 3.0;
      }
      else {modVSEPR(idAtm + 1,2) = 0.0;}
    }
    else if (angles[1] == 1) {modVSEPR(idAtm + 1,2) = 1.0;}                                         //AX2E1
    else if (angles[2] == 1) {modVSEPR(idAtm + 1,2) = 2.0;}                                         //AX2E2
    else if (angles[3] == 1) {modVSEPR(idAtm + 1,2) = 3.5;}                                         //AX2E3 -> this is a new class I defined in order to account for angles in sulfides; it is derived from T-shaped and replacing one X by E
    else {modVSEPR(idAtm + 1,2) = -2.0;}                                                            //unidentified
  }
  else if (nneigh == 3) {
    if (angles[1] == 3) {modVSEPR(idAtm + 1,2) = 0.0;}                                              //AX3E0
    else if ((angles[1] == 2)&&(atomA == 7)&&(angles[2] == 1)) {modVSEPR(idAtm + 1,2) = 0.0;}       //AX3E0 -> those weird aromatic nitrogens
    else if ((angles[1] == 2)&&(atomA == 6)&&(angles[2] == 1)) {modVSEPR(idAtm + 1,2) = 0.0;}       //AX3E0 -> some aromatic carbons are weird too
    else if (angles[2] == 3) {modVSEPR(idAtm + 1,2) = 1.0;}                                         //AX3E1
    else if (angles[3] == 3) {modVSEPR(idAtm + 1,2) = 1.0;}                                         //AX3E1 -> wrong assignment on purpose, to ensure I get sp3
    else if ((angles[0] == 1)&&(angles[3] == 2)) {modVSEPR(idAtm + 1,2) = 2.0;}                     //AX3E2
    else {modVSEPR(idAtm + 1,2) = -2.0;}                                                            //unidentified
  }
  else if (nneigh == 4) {
    if (angles[2] == 6) {modVSEPR(idAtm + 1,2) = 0.0;}                                              //AX4E0
    else if ((angles[0] == 1)&&(angles[1] == 1)&&(angles[3] == 4)) {modVSEPR(idAtm + 1,2) = 1.0;}   //AX4E1
    else if (angles[3] == 4) {modVSEPR(idAtm + 1,2) = 2.0;}                                         //AX4E2
    else {modVSEPR(idAtm + 1,2) = -2.0;}                                                            //unidentified
  }
  //from this point on things get really tricky
  else if (nneigh == 5) {
    if ((angles[3] == 6)&&(angles[1] == 3)&&(angles[0] == 1)) {modVSEPR(idAtm + 1,2) = 0.0;}        //AX5E0
    else if ((angles[3] == 8)&&(angles[0] == 2)) {modVSEPR(idAtm + 1,2) = 1.0;}                     //AX5E1
    else if (angles[4] == 5) {modVSEPR(idAtm + 1,2) = 2.0;}                                         //AX5E2
    else {modVSEPR(idAtm + 1,2) = -2.0;}                                                            //unidentified
  }
  else if (nneigh == 6) {
    if ((angles[3] == 12)&&(angles[0] == 3)) {modVSEPR(idAtm + 1,2) = 0.0;}                         //AX6E0
    else if ((angles[4] == 5)&&(angles[3] == 5)) {modVSEPR(idAtm + 1,2) = 1.0;}                     //AX6E1
    else {modVSEPR(idAtm + 1,2) = -2.0;}                                                            //unidentified
  }
  //not going beyond this for now
}
matrixE AtomicGeometry(std::vector<size_t> & atoms, matrixE & geom, matrixE & Ncoord, double threshang = 6.0, size_t maxneigh = 8, double threshbond = 0.3) {
  //function that applies some modified form of VSEPR to determine the type of geometry around each atom
  //the return matrix has dimensions of the total number of atoms times two; for each atom we give the subscript for X (number of neighbours) and for E (which is obtained from the geometry)
  //possible outcomes for each atom:
  //    2,0  -> linear arrangement (AX2E0)
  //    3,0  -> trigonal planar arrangement (AX3E0)
  //    2,1  -> bent arrangement (AX2E1)
  //    4,0  -> tetrahedral arrangement (AX4E0)
  //    3,1  -> trigonal pyramidal arrangement (AX3E1)
  //    2,2  -> bent arrangement (AX2E2)
  //    5,0  -> trigonal bipyramidal arrangement (AX5E0)
  //    4,1  -> seesaw arrangement (AX4E1)
  //    3,2  -> t-shaped arrangement (AX3E2)
  //    2,3  -> linear arrangement (AX2E3)
  //    6,0  -> octahedral arrangement (AX6E0)
  //    5,1  -> square pyramidal arrangement (AX5E1)
  //    4,2  -> square planar arrangement (AX4E2)
  //    7,0  -> pentagonal bipyramidal arrangement (AX7E0)
  //    6,1  -> pentagonal pyramidal pyramidal arrangement (AX6E1)
  //    5,2  -> pentagonal planar arrangement (AX5E2)
  //    8,0  -> square antiprimatic arrangement (AX8E0)
  //if second position outputs -1, then terminal atom; else if second position outputs -2, then atom not identified
  //threshang is the error allowed in angles for assigning geometry
  //maxneigh is the maximum number of neighbours allowed; this is limited by the maximum allowed by VSEPR
  //threshbond is the minimum value of coordination number to consider a bond
  size_t Natoms = atoms.size();
  std::vector<int> friendlyneighbours(maxneigh);
  int nneigh;
  int atomA;
  std::vector<int> angles(5);
  matrix<int> islinear(maxneigh,maxneigh);
  matrix<int> istrigonal(maxneigh,maxneigh);
  matrix<int> istetrahedral(maxneigh,maxneigh);
  matrix<int> isorthogonal(maxneigh,maxneigh);
  matrix<int> ispentagonal(maxneigh,maxneigh);
  matrixE theta(maxneigh,maxneigh);
  matrixE modVSEPR(Natoms,2);
  for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
    atomA = atoms[idAtm];
    if (atomA == 1) {
      modVSEPR(idAtm + 1,1) = 1.0;
      modVSEPR(idAtm + 1,2) = -1.0;
      continue;
    }
    //count neighbours
    nneigh = CountNeighbours(Ncoord,friendlyneighbours,idAtm,Natoms,threshbond);
    modVSEPR(idAtm + 1,1) = double(nneigh);
    //now get all the angles and classify them
    GetAngles(angles,friendlyneighbours,geom,islinear,istrigonal,istetrahedral,isorthogonal,ispentagonal,theta,idAtm,nneigh,threshang);
    ClassifyVSEPR(modVSEPR,angles,idAtm,atomA,nneigh,threshang);
  }
  return modVSEPR;
}

#endif //_VSEPR_Package_
