/*ULYSSES, a semi-empirical package
    Copyright (C) 2023- Filipe Menezes (filipe.menezes@helmholtz-munich.de)
                        Federico Ballabio (federico.ballabio@helmholtz-munich.de)
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

#ifndef _Add_H_Package_
#define _Add_H_Package_
#include <vector>
#include <fstream>
#include <iostream>

//description:
//functions to add protons

double HXdistances(int atomX, int hybdridisation) {
  //function with tabulated bond distances
  double rHX = 0.0;
  switch(atomX) {
    case 5:
      rHX = 1.15;
      break;
    case 6:
      rHX = 1.10;                                    //sp3
      if (hybdridisation == 2) {rHX = 1.09;}         //sp2
      else if (hybdridisation == 1) {rHX = 1.08;}    //sp
      break;
    case 7:
      rHX = 1.01;                                    //sp3+
      if (hybdridisation == 3) {rHX = 1.01;}         //sp3
      else if (hybdridisation == 2) {rHX = 1.01;}    //sp2
      break;
    case 8:
      rHX = 0.98;                                    //sp3
      break;
    case 16:
      rHX = 1.38;                                    //sp3
      break;
    case 34:
      rHX = 1.46;                                    //sp3
      break;
  }
  return rHX;
}
double NormAndDistance(int idAtm, int idBtm, matrixE & geometry, std::vector<double> & dist) {
  //function calculating the distance between two points and the respective norm
  double norm = 0.0;
  for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
    dist[idcoord] = geometry(idAtm + 1,idcoord + 1) - geometry(idBtm + 1,idcoord + 1);
    norm += dist[idcoord]*dist[idcoord];
  }
  double aux = sqrt(norm);
  norm = 1.0/aux;
  return norm;
}
void NormAndDistanceTwice(int idAtm, int idBtm, int idCtm, matrixE & geometry, std::vector<double> & dist1, double & norm1, std::vector<double> & dist2, double & norm2) {
  //function calculating the distance between two points and the respective norm
  norm1 = 0.0;
  norm2 = 0.0;
  for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
    dist1[idcoord] = geometry(idAtm + 1,idcoord + 1) - geometry(idBtm + 1,idcoord + 1);
    dist2[idcoord] = geometry(idAtm + 1,idcoord + 1) - geometry(idCtm + 1,idcoord + 1);
    norm1 += dist1[idcoord]*dist1[idcoord];
    norm2 += dist2[idcoord]*dist2[idcoord];
  }
  double aux = sqrt(norm1);
  norm1 = 1.0/aux;
  aux = sqrt(norm2);
  norm2 = 1.0/aux;
}
bool CalculateRotationMatrix(matrixE & Rotation, std::vector<double> & vA, std::vector<double> & vB, std::vector<double> & vBalt, std::vector<double> & vVec, matrixE & vMat, matrixE & vMat2, double norm) {
  //function that calculates a rotation matrix to make two vectors colinear
  //the matrix is good to rotate vB into vA, i.e., vA = R.vB
  double cost = 0.0;
  for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
    vA[idcoord] *= norm;
    cost += vA[idcoord]*vB[idcoord];
  }
  bool swap = false;
  if (cost < -1.0 + 0.00001) {
    swap = true;
    cost = 0.0;
    for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
      cost += vA[idcoord]*vBalt[idcoord];
    }
    vVec = crossProd(vBalt,vA,1.0);
  }
  else {vVec = crossProd(vB,vA,1.0);}
  vMat(1,1) = 0.0;
  vMat(1,2) = -vVec[2];
  vMat(1,3) = vVec[1];
  vMat(2,1) = vVec[2];
  vMat(2,2) = 0.0;
  vMat(2,3) = -vVec[0];
  vMat(3,1) = -vVec[1];
  vMat(3,2) = vVec[0];
  vMat(3,3) = 0.0;
  vMat2 = vMat*vMat;
  vMat2 *= 1.0/(1.0 + cost);
  Rotation = Identity(3);
  Rotation += vMat;
  Rotation += vMat2;
  return swap;
}
void RotateMatrix(matrixE & ML, matrixE & Rotation, matrixE & MR, matrixE & aux) {
  //function applying a certain rotation operation to a matrix
  aux = Rotation;
  aux *= MR;
  ML = aux;
}
void AddNewNeighbour(matrix<int> & connectivity, int idAtm, int newneigh, int maxneigh) {
  //function adding a new neighbour to the connectivity matrix
  for (size_t idneigh = 0; idneigh < maxneigh; ++idneigh) {
    if (connectivity(idAtm + 1,idneigh + 1) == 0) {
      connectivity(idAtm + 1,idneigh + 1) = newneigh;
      break;
    }
  }
}
void AddProtons(Molecule & TheSystem, matrixE & VSEPRcodes, matrix<int> & connectivity) {
  //function to add protons
  double norm;
  double norm2;
  double norm3;
  double norm4;
  double bdist;
  double osqrt3 = 1.0/sqrt(3.0);
  std::vector<size_t> atoms = TheSystem.Atoms();
  std::vector<double> rAB(3);
  std::vector<double> rAC(3);
  std::vector<double> rAD(3);
  std::vector<double> rsum(3);
  std::vector<double> Tetrahedron1(3,osqrt3);           //vertices of tetrahedron      ( 1, 1, 1)
  std::vector<double> Tetrahedron2(3,osqrt3);           //vertices of tetrahedron      (-1,-1, 1)
  std::vector<double> Tetrahedron3(3,osqrt3);           //vertices of tetrahedron      (-1, 1,-1)
  std::vector<double> Tetrahedron4(3,osqrt3);           //vertices of tetrahedron      ( 1,-1,-1)
  Tetrahedron2[0] *= -1.0;
  Tetrahedron2[1] *= -1.0;
  Tetrahedron3[0] *= -1.0;
  Tetrahedron3[2] *= -1.0;
  Tetrahedron4[1] *= -1.0;
  Tetrahedron4[2] *= -1.0;
  std::vector<double> vVec(3);
  matrixE Rotation(3,3);
  matrixE vMat(3,3);
  matrixE vMat2(3,3);
  matrixE mTetrahedron(3,3);
  matrixE rTetrahedron1(3,3);
  matrixE rTetrahedron2(3,3);
  matrixE rTetrahedron3(3,3);
  matrixE rTetrahedron4(3,3);
  matrixE FakeZmat(1,6);
  matrixE geom = TheSystem.Geometry();
  size_t Natoms = atoms.size();
  size_t Natoms_new = Natoms;
  int idBtm;
  int idCtm;
  int idDtm;
  int maxneigh = connectivity.cols();
  bool swap;
  bool swap2;
  for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
    idBtm = connectivity(idAtm + 1,1) - 1;
    //sp atoms
    if ((atoms[idAtm] == 6)&&(VSEPRcodes(idAtm + 1,1) == 2.0)&&(VSEPRcodes(idAtm + 1,2) == 0.0)&&(connectivity(idAtm + 1,2) == 0)) {
      bdist = HXdistances(6,1);
      ++Natoms_new;
      geom.resize(Natoms_new,3);
      atoms.push_back(1);
      //get the distance
      norm = NormAndDistance(idAtm,idBtm,geom,rAB)*bdist;
      for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
        geom(Natoms_new,idcoord + 1) = geom(idAtm + 1,idcoord + 1) + rAB[idcoord]*norm;
      }
      AddNewNeighbour(connectivity,idAtm,Natoms_new,maxneigh);
    }
    //sp2 atoms
    else if (((VSEPRcodes(idAtm + 1,1) == 3.0)&&(VSEPRcodes(idAtm + 1,2) == 0.0))||((VSEPRcodes(idAtm + 1,1) == 2.0)&&(VSEPRcodes(idAtm + 1,2) == 1.0))) {
      bdist = HXdistances(6,2)*(atoms[idAtm] == 6) + HXdistances(7,2)*(atoms[idAtm] == 7) + HXdistances(8,2)*(atoms[idAtm] == 8);
      if (connectivity(idAtm + 1,2) == 0) {
        //one neighbour, then fetch the coordinates from the neighbour
        if (connectivity(idBtm + 1,3) != 0) {
          //neighbour has 3 neighbours, so we just mirror it
          if (connectivity(idBtm + 1,1) == idAtm + 1) {
            idCtm = connectivity(idBtm + 1,2) - 1;
            idDtm = connectivity(idBtm + 1,3) - 1;
          }
          else if (connectivity(idBtm + 1,2) == idAtm + 1) {
            idCtm = connectivity(idBtm + 1,1) - 1;
            idDtm = connectivity(idBtm + 1,3) - 1;
          }
          else if (connectivity(idBtm + 1,3) == idAtm + 1) {
            idCtm = connectivity(idBtm + 1,1) - 1;
            idDtm = connectivity(idBtm + 1,2) - 1;
          }
          NormAndDistanceTwice(idBtm,idCtm,idDtm,geom,rAC,norm,rAD,norm2);
          norm *= bdist;
          norm2 *= bdist;
          ++Natoms_new;
          geom.resize(Natoms_new,3);
          atoms.push_back(1);
          for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
            geom(Natoms_new,idcoord + 1) = geom(idAtm + 1,idcoord + 1) + rAC[idcoord]*norm;
          }
          AddNewNeighbour(connectivity,idAtm,Natoms_new,maxneigh);
          if (VSEPRcodes(idAtm + 1,1) == 3.0) {
            ++Natoms_new;
            geom.resize(Natoms_new,3);
            atoms.push_back(1);
            for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
              geom(Natoms_new,idcoord + 1) = geom(idAtm + 1,idcoord + 1) + rAD[idcoord]*norm2;
            }
            AddNewNeighbour(connectivity,idAtm,Natoms_new,maxneigh);
          }
        }
        else if (connectivity(idBtm + 1,2) != 0) {
          //neighbour only has 2 neighbours, so take the only undefined neighbour available and place it trans
          idCtm = connectivity(idBtm + 1,2) - 1;
          if (idCtm == idAtm) {idCtm = connectivity(idBtm + 1,1) - 1;}
          //prepare a pseudo-fake Zmat
          FakeZmat(1,2) = bdist;
          FakeZmat(1,4) = Angle(idAtm + 1,idBtm + 1,idCtm + 1,geom,1.0e-20);
          FakeZmat(1,6) = pi;
          //get the coordinates
          GetUnitVectors(rAB,rsum,geom,idBtm + 1,idAtm + 1,idCtm + 1);
          GetLocalAxes(rAB,rsum,rAC);
          GetBondVector(rAD,FakeZmat,1);
          //now add the new guy
          ++Natoms_new;
          geom.resize(Natoms_new,3);
          atoms.push_back(1);
          geom(Natoms_new,1) = geom(idAtm + 1,1) + dotProd(rAD,rAC);
          geom(Natoms_new,2) = geom(idAtm + 1,2) + dotProd(rAD,rsum);
          geom(Natoms_new,3) = geom(idAtm + 1,3) + dotProd(rAD,rAB);
          //now add the last neighbour
          if (VSEPRcodes(idAtm + 1,1) == 3.0) {
            FakeZmat(1,6) = 0.0;
            //get the coordinates
            GetBondVector(rAD,FakeZmat,1);
            //now add the new guy
            ++Natoms_new;
            geom.resize(Natoms_new,3);
            atoms.push_back(1);
            geom(Natoms_new,1) = geom(idAtm + 1,1) + dotProd(rAD,rAC);
            geom(Natoms_new,2) = geom(idAtm + 1,2) + dotProd(rAD,rsum);
            geom(Natoms_new,3) = geom(idAtm + 1,3) + dotProd(rAD,rAB);
          }
        }
      }
      else if ((connectivity(idAtm + 1,3) == 0)&&(VSEPRcodes(idAtm + 1,1) == 3.0)) {
        //two neighbours
        idCtm = connectivity(idAtm + 1,2) - 1;
        ++Natoms_new;
        geom.resize(Natoms_new,3);
        atoms.push_back(1);
        NormAndDistanceTwice(idAtm,idBtm,idCtm,geom,rAB,norm,rAC,norm2);
        norm3 = 0.0;
        for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
          rsum[idcoord] = rAB[idcoord]*norm + rAC[idcoord]*norm2;
          norm3 += rsum[idcoord]*rsum[idcoord];
        }
        norm = bdist/sqrt(norm3);
        for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
          geom(Natoms_new,idcoord + 1) = geom(idAtm + 1,idcoord + 1) + rsum[idcoord]*norm;
        }
        AddNewNeighbour(connectivity,idAtm,Natoms_new,maxneigh);
      }
    }
    //sp3 atoms
    else if ((VSEPRcodes(idAtm + 1,1) == 2.0)&&(VSEPRcodes(idAtm + 1,2) == 2.0)&&(connectivity(idAtm + 1,2) == 0)) {
      //these are alcohols and thiols
      bdist = HXdistances(8,3)*(atoms[idAtm] == 8) + HXdistances(16,3)*(atoms[idAtm] == 16) + HXdistances(34,3)*(atoms[idAtm] == 34);
      //look for non-hydrogens
      idCtm = -1;
      if (idBtm >= 0) {
        for (size_t idneigh = 0; idneigh < maxneigh; ++idneigh) {
          idCtm = connectivity(idBtm + 1,idneigh + 1) - 1;
          if (idCtm < 0) {break;}
          if ((idCtm != idAtm)&&(atoms[idCtm] != 1)) {break;}
        }
        //if this cannot be found, then lets take a hydrogen
        if (idCtm < 0) {
          for (size_t idneigh = 0; idneigh < maxneigh; ++idneigh) {
            idCtm = connectivity(idBtm + 1,idneigh + 1) - 1;
            if (idCtm < 0) {break;}
            if (idCtm != idAtm) {break;}
          }
        }
      }
      if (idCtm >= 0) {
        FakeZmat(1,2) = bdist;
        FakeZmat(1,4) = Angle(idAtm + 1,idBtm + 1,idCtm + 1,geom,1.0e-20);
        FakeZmat(1,6) = pi;
        //get the coordinates
        GetUnitVectors(rAB,rsum,geom,idBtm + 1,idAtm + 1,idCtm + 1);
        GetLocalAxes(rAB,rsum,rAC);
        GetBondVector(rAD,FakeZmat,1);
        //now add the new guy
        ++Natoms_new;
        geom.resize(Natoms_new,3);
        atoms.push_back(1);
        geom(Natoms_new,1) = geom(idAtm + 1,1) + dotProd(rAD,rAC);
        geom(Natoms_new,2) = geom(idAtm + 1,2) + dotProd(rAD,rsum);
        geom(Natoms_new,3) = geom(idAtm + 1,3) + dotProd(rAD,rAB);
      }
      else {
        ++Natoms_new;
        geom.resize(Natoms_new,3);
        atoms.push_back(1);
        geom(Natoms_new,1) = geom(idAtm + 1,1) - bdist;
        geom(Natoms_new,2) = geom(idAtm + 1,2);
        geom(Natoms_new,3) = geom(idAtm + 1,3);
        ++Natoms_new;
        geom.resize(Natoms_new,3);
        atoms.push_back(1);
        geom(Natoms_new,1) = geom(idAtm + 1,1);
        geom(Natoms_new,2) = geom(idAtm + 1,2) - bdist;
        geom(Natoms_new,3) = geom(idAtm + 1,3);
      }
    }
    else if (((VSEPRcodes(idAtm + 1,1) == 4.0)&&(VSEPRcodes(idAtm + 1,2) == 0.0))||((VSEPRcodes(idAtm + 1,1) == 3.0)&&(VSEPRcodes(idAtm + 1,2) == 1.0))) {
      bdist = HXdistances(5,3)*(atoms[idAtm] == 5) + HXdistances(6,3)*(atoms[idAtm] == 6) + HXdistances(7,3)*(atoms[idAtm] == 7) + HXdistances(16,3)*(atoms[idAtm] == 16);
      //                          1 neighbour                       2 neighbours
      if ((connectivity(idAtm + 1,2) == 0)||(connectivity(idAtm + 1,3) == 0)) {
        //get the bond vector towards first neighbour
        norm = NormAndDistance(idBtm,idAtm,geom,rAB);
        //get rotation matrix to align the tetrahedron to the first bond vector
        swap = CalculateRotationMatrix(Rotation,rAB,Tetrahedron1,Tetrahedron2,vVec,vMat,vMat2,norm);
        //align the tetrahedron according to the bond vector
        rTetrahedron1 = Tetrahedron1;
        rTetrahedron2 = Tetrahedron2;
        rTetrahedron3 = Tetrahedron3;
        rTetrahedron4 = Tetrahedron4;
        if (swap) {
          RotateMatrix(rTetrahedron1,Rotation,rTetrahedron2,mTetrahedron);
          RotateMatrix(rTetrahedron2,Rotation,rTetrahedron1,mTetrahedron);
        }
        else {
          RotateMatrix(rTetrahedron1,Rotation,rTetrahedron1,mTetrahedron);
          RotateMatrix(rTetrahedron2,Rotation,rTetrahedron2,mTetrahedron);
        }
        RotateMatrix(rTetrahedron3,Rotation,rTetrahedron3,mTetrahedron);
        RotateMatrix(rTetrahedron4,Rotation,rTetrahedron4,mTetrahedron);
        if (connectivity(idAtm + 1,2) != 0) {
          //two neighbours, then re-align the tetrahedron including the second bond vector
          idCtm = connectivity(idAtm + 1,2) - 1;
          norm = NormAndDistance(idCtm,idAtm,geom,rAC);
          //get rotation matrix
          if (swap) {swap2 = CalculateRotationMatrix(Rotation,rAC,Tetrahedron1,Tetrahedron2,vVec,vMat,vMat2,norm);}
          else {swap2 = CalculateRotationMatrix(Rotation,rAC,Tetrahedron2,Tetrahedron1,vVec,vMat,vMat2,norm);}
          //and rotate everything
          if (swap) {
            RotateMatrix(rTetrahedron1,Rotation,rTetrahedron2,mTetrahedron);
            RotateMatrix(rTetrahedron2,Rotation,rTetrahedron1,mTetrahedron);
          }
          else {
            RotateMatrix(rTetrahedron1,Rotation,rTetrahedron1,mTetrahedron);
            RotateMatrix(rTetrahedron2,Rotation,rTetrahedron2,mTetrahedron);
          }
          RotateMatrix(rTetrahedron3,Rotation,rTetrahedron3,mTetrahedron);
          RotateMatrix(rTetrahedron4,Rotation,rTetrahedron4,mTetrahedron);
          //the vectors might still not be well aligned due to deviations of angles from ideality
          //this requires two additional rotations
          norm2 = 0.0;
          norm4 = 0.0;
          for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
            rAD[idcoord] = rAB[idcoord] + rAC[idcoord];
            norm2 += rAD[idcoord]*rAD[idcoord];
            rsum[idcoord] = rTetrahedron3(idcoord + 1,1) + rTetrahedron4(idcoord + 1,1);
            norm4 += rsum[idcoord]*rsum[idcoord];
          }
          norm3 = 1.0/sqrt(norm2);
          norm = 1.0/sqrt(norm4);
          for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
            rAD[idcoord] *= norm3;
            rsum[idcoord] *= -norm;
          }
          //get and apply the first rotation
          swap2 = CalculateRotationMatrix(Rotation,rAD,rsum,rsum,vVec,vMat,vMat2,1.0);
          RotateMatrix(rTetrahedron3,Rotation,rTetrahedron3,mTetrahedron);
          RotateMatrix(rTetrahedron4,Rotation,rTetrahedron4,mTetrahedron);
          //now the second rotation
          rAD = crossProd(rAB,rAC,1.0);
          norm2 = 0.0;
          norm = 0.0;
          for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
            rsum[idcoord] = rTetrahedron4(idcoord + 1,1) - rTetrahedron3(idcoord + 1,1);
            norm2 += rsum[idcoord]*rsum[idcoord];
            norm += rAD[idcoord]*rAD[idcoord];
          }
          norm3 = 1.0/sqrt(norm2);
          norm2 = 1.0/sqrt(norm);
          norm4 = 0.0;
          for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
            rsum[idcoord] *= norm3;
            rAD[idcoord]*= norm2;
            norm4 += rAD[idcoord]*rsum[idcoord];
          }
          if (norm4 < 0.0000001) {
            for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
              rsum[idcoord] *= -1.0;
            }
          }
          swap2 = CalculateRotationMatrix(Rotation,rAD,rsum,rsum,vVec,vMat,vMat2,1.0);
          RotateMatrix(rTetrahedron3,Rotation,rTetrahedron3,mTetrahedron);
          RotateMatrix(rTetrahedron4,Rotation,rTetrahedron4,mTetrahedron);
        }
        else {
          //one neighbour
          ++Natoms_new;
          geom.resize(Natoms_new,3);
          atoms.push_back(1);
          for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
            geom(Natoms_new,idcoord + 1) = geom(idAtm + 1,idcoord + 1) + rTetrahedron2(idcoord + 1,1)*bdist;
          }
          AddNewNeighbour(connectivity,idAtm,Natoms_new,maxneigh);
        }
        //now we add the other neighbours
        ++Natoms_new;
        geom.resize(Natoms_new,3);
        atoms.push_back(1);
        for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
          geom(Natoms_new,idcoord + 1) = geom(idAtm + 1,idcoord + 1) + rTetrahedron3(idcoord + 1,1)*bdist;
        }
        AddNewNeighbour(connectivity,idAtm,Natoms_new,maxneigh);
        if ((VSEPRcodes(idAtm + 1,1) == 4.0)&&(VSEPRcodes(idAtm + 1,2) == 0.0)) {
          ++Natoms_new;
          geom.resize(Natoms_new,3);
          atoms.push_back(1);
          for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
            geom(Natoms_new,idcoord + 1) = geom(idAtm + 1,idcoord + 1) + rTetrahedron4(idcoord + 1,1)*bdist;
          }
          AddNewNeighbour(connectivity,idAtm,Natoms_new,maxneigh);
        }
      }
      //                              3 neighbours
      else if ((connectivity(idAtm + 1,4) == 0)&&(VSEPRcodes(idAtm + 1,1) == 4.0)) {
        //three neighbours
        idCtm = connectivity(idAtm + 1,2) - 1;
        idDtm = connectivity(idAtm + 1,3) - 1;
        ++Natoms_new;
        geom.resize(Natoms_new,3);
        atoms.push_back(1);
        //get the distances
        norm = NormAndDistance(idAtm,idBtm,geom,rAB);
        NormAndDistanceTwice(idAtm,idCtm,idDtm,geom,rAC,norm2,rAD,norm3);
        norm4 = 0.0;
        //normalise
        for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
          rsum[idcoord] = rAB[idcoord]*norm + rAC[idcoord]*norm2 + rAD[idcoord]*norm3;
          norm4 += rsum[idcoord]*rsum[idcoord];
        }
        norm = bdist/sqrt(norm4);
        for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
          geom(Natoms_new,idcoord + 1) = geom(idAtm + 1,idcoord + 1) + rsum[idcoord]*norm;
        }
        AddNewNeighbour(connectivity,idAtm,Natoms_new,maxneigh);
      }
    }
  }
  TheSystem.setGeometry(geom);
  TheSystem.setAtoms(atoms);
}

#endif //_Add_H_Package_
