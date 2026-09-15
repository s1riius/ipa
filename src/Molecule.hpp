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

#ifndef _Molecule_
#define _Molecule_
#include <vector>
#include <string>
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <ctype.h>
#include "Global.hpp"
#include "Geometry.hpp"
#include "atoms/AtomPackage.hpp"
#include "parameters/AtomicRadiipar.hpp"
#include "parameters/BondDistancepar.hpp"
#include "ConstantsPackage.hpp"
#include "other/listPackage.hpp"
#include "other/auxPackage.hpp"
#include "other/FileFormats.hpp"
#include "other/MoleculeAnalyser.hpp"
#include "math/MatrixPackage.hpp"
#include "math/RandomPackage.hpp"
#include "math/VectorPackage.hpp"

//description:
//molecule class which is the basis for quantum chemical calculations

matrixE sigma(size_t pos) {
  matrixE _sigma = Identity(3);
  _sigma(pos,pos) = -1;
  return _sigma;
}
matrixE cn(size_t pos, int nm) {
  double cs = cos(2.0*pi/nm);
  double sn = sin(2.0*pi/nm);
  size_t row;
  size_t col;
  if (pos == 1) {
    row = 2;
    col = 3;
  }
  else if (pos == 2) {
    row = 1;
    col = 3;
  }
  else if (pos == 3) {
    row = 1;
    col = 2;
  }
  return Givens(row,col,3,cs,sn);
}
class Molecule {
  matrixE geometry;
  std::vector<size_t> atoms;
  std::string geomfile;
  size_t natoms;
  size_t mult;
  int nelectrons;
  int charge;
  int maxcoordinationnumber;               //the maximum coordination number allowed in getting connectivity matrix
  double MWeight;
  double mass;
  std::string symm;
public:
  Molecule(): geometry(1,1) {
    maxcoordinationnumber = 8;
  }
  Molecule(std::string _file, int _chg = 0, size_t _mult = 1, std::string _symm = "0", bool centerstructure = true) {
    charge = _chg;
    geomfile = _file;
    std::string extension = GetFileExtension(_file);
    if (extension == "xyz") {ReadXYZ();}
    else if (extension == "mol2") {ReadMOL2();}
    else if (extension == "sdf") {ReadSDF();}
    /////////else if (extension == "gjf") {ReadGJF();}
    else if (extension == "pdb") {ReadPDB();}
    //if (Linear()) {linearGeom();}
    //else if (Planar()) {planarGeom();}
    masses();
    CountElectrons();
    mult = _mult;
    if ((nelectrons%2 != 0)&&(mult == 1)) {mult = 2;}
    symm = _symm;
    if (centerstructure) {CorrectedGeom();}
    if (symm == "0") {PointGroup();}
    maxcoordinationnumber = 8;
  }
  Molecule(const Molecule & rhs): geometry(rhs.Geometry()) {
    geometry = rhs.Geometry();
    atoms = rhs.Atoms();
    geomfile = rhs.GeomFile();
    natoms = rhs.Natoms();
    mult = rhs.Multiplicity();
    nelectrons = rhs.Nelectrons();
    charge = rhs.Charge();
    MWeight = rhs.MW();
    mass = rhs.m();
    symm = rhs.PGroup();
    maxcoordinationnumber = 8;
  }
  ~Molecule() {}
//--------------------------------- Counting functions
  int MoleculeCounter(std::vector<size_t> & molvector) {
    if (molvector.size() != natoms) {molvector.resize(natoms);}
    std::vector<bool> assigned(natoms,false);
    matrix<int> neighbours;
    ConnectivityMatrix(neighbours);                 //getting neighbours
    int counter = 1;
    bool neighbourassigned;
    size_t maxnrneighbours = neighbours.cols();
    for (size_t idAtm = 0; idAtm < natoms; ++idAtm) {
      molvector[idAtm] = 0;
    }
    for (size_t idAtm = 0; idAtm < natoms; ++idAtm) {
      if (!assigned[idAtm]) {
        neighbourassigned = false;
        for (size_t idneigh = 0; idneigh < maxnrneighbours; ++idneigh) {
          if (neighbours(idAtm + 1,idneigh + 1) - 1 < 0) {break;}
          if (assigned[neighbours(idAtm + 1,idneigh + 1) - 1]) {
            molvector[idAtm] = molvector[neighbours(idAtm + 1,idneigh + 1) - 1];
            neighbourassigned = true;
            break;
          }
        }
        if (!neighbourassigned) {
          molvector[idAtm] = counter;
          ++counter;
        }
        assigned[idAtm] = true;
      }
    }
    --counter;
    return counter;
  }
//--------------------------------- functions to add molecules
  void set2System(matrixE & geom, std::vector<size_t> & atm, int chrge, int multpl, std::string pg) {
    //function that sets the molecule to external variables
    this->setGeometry(geom);
    this->setAtoms(atm);
    this->calcNatoms();
    this->setCharge(chrge);
    this->setMultiplicity(multpl);
    this->masses();
    this->setPGroup(pg);
  }
  void AddMolecule(Molecule & extramol) {
    matrixE geomextra = extramol.Geometry();
    std::vector<size_t> extraatoms = extramol.Atoms();
    AddMolecule(geomextra,extraatoms,extramol.Charge(),extramol.Multiplicity());
  }
  void AddMolecule(matrixE & extraG, std::vector<size_t> & extraA, int extrachg, size_t extramult) {
    ConcatenateMR(geometry,extraG);
    ConcatenateV(atoms,extraA);
    charge += extrachg;
    mult *= extramult;
    this->calcNatoms();
    this->masses();
  }
//--------------------------------- Geometry Functions
  void Make3D(double maxshift = 0.5, bool complexanalysis = false, double thresholdzero = 1.0e-7) {
    //function that takes a potentially flat molecule (one dimension is zero) and generates distortions on the flat dimension
    bool isflat[3];
    double randomnumber;
    int randomsign;
    isflat[0] = true;
    isflat[1] = true;
    isflat[2] = true;
    //determine whether molecule is flat
    if (complexanalysis) {
      //check for planarity
      maxshift /= 3.0;
      isflat[0] = Planar();
      isflat[1] = isflat[0];
      isflat[2] = isflat[0];
    }
    else {
      //check for zeroed dimensions
      for (size_t idAtm = 0; idAtm < natoms; ++idAtm) {
        if (fabs(geometry(idAtm + 1,1)) > thresholdzero) {isflat[0] = false;}
        if (fabs(geometry(idAtm + 1,2)) > thresholdzero) {isflat[1] = false;}
        if (fabs(geometry(idAtm + 1,3)) > thresholdzero) {isflat[2] = false;}
        if (isflat[0] + isflat[1] + isflat[2] == 0) {break;}
      }
    }
    //remove flatness
    for (size_t idAtm = 0; idAtm < natoms; ++idAtm) {
      randomnumber = fRandom(0.0,1.0);
      randomsign = iRandom01();
      geometry(idAtm + 1,1) += (2*randomsign - 1)*randomnumber*maxshift;
      randomnumber = fRandom(0.0,1.0);
      randomsign = iRandom01();
      geometry(idAtm + 1,2) += (2*randomsign - 1)*randomnumber*maxshift;
      randomnumber = fRandom(0.0,1.0);
      randomsign = iRandom01();
      geometry(idAtm + 1,3) += (2*randomsign - 1)*randomnumber*maxshift;
    }
  }
  void AlignTo(Molecule & reference, matrixE & align, int firstN = -1) {
    //function that aligns this molecule to another one
    matrixE rgeom = reference.Geometry();
    AlignBonA(rgeom,geometry,align,firstN);
  }
  void AlignTo(Molecule & reference) {
    //function that aligns this molecule to another one
    matrixE rgeom = reference.Geometry();
    matrixE align;
    AlignBonA(rgeom,geometry,align,natoms);
  }
  void setAtoms(std::vector<size_t> * _atoms) {
    atoms = *_atoms;
    natoms = atoms.size();
  }
  void setAtoms(std::vector<size_t> & _atoms) {
    atoms = _atoms;
    natoms = atoms.size();
  }
  void setGeometry(matrixE * _geometry) {geometry = *_geometry;}
  void setGeometry(matrixE & _geometry) {geometry = _geometry;}
  std::vector<size_t> Atoms() {return atoms;}
  std::vector<size_t> Atoms() const {return atoms;}
  size_t Atom(size_t pos) {return atoms[pos - 1];}
  size_t Atom(size_t pos) const {return atoms[pos - 1];}
  matrixE Geometry() {return geometry;}
  matrixE Geometry() const {return geometry;}
  void transGeometry(std::vector<double> translation) {
    //function that translates geometry according to vector
    if (translation.size() != 3) {throw std::string("ERROR: Molecule.hpp: Molecule: transGeometry(): translation with vector not in R3");}
    for (size_t iatom = 0; iatom < natoms; ++iatom) {
      for (size_t icoord = 0; icoord < 3; ++icoord) {
        geometry(iatom + 1,icoord + 1) += translation[icoord];
      }
    }
  }
  void rotGeometry(size_t axis, double theta) {
    //function that rotates geometry according to certain Given's transformation matrix
    std::vector<double> cmass = this->CM();
    //recenter molecule
    for (size_t idatm = 0; idatm < geometry.rows(); ++idatm) {
      for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
        geometry(idatm + 1,idcoord + 1) -= cmass[idcoord];
      }
    }
    matrixE rotation;
    RotateGeomMatrix(axis,theta,geometry,rotation);
    //shift molecule back to its position
    for (size_t idatm = 0; idatm < geometry.rows(); ++idatm) {
      for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
        geometry(idatm + 1,idcoord + 1) += cmass[idcoord];
      }
    }
  }
  void rotGeometry2(size_t axis, double theta, matrixE & rotation) {
    //like function before, but quick and "dirty" as it assumes structure is centered
    RotateGeomMatrix(axis,theta,geometry,rotation);
  }
  std::string GeomFile() {return geomfile;}
  std::string GeomFile() const {return geomfile;}
  void ReadXYZ() {ReadXYZFormat(geomfile,natoms,geometry,atoms);}
  void ReadMOL2() {ReadMOL2Format(geomfile,natoms,geometry,atoms,100000000);}
  void ReadSDF() {ReadSDFFormat(geomfile,natoms,geometry,atoms,charge,100000000);}
  void ReadPDB() {ReadPDBFormat(geomfile,natoms,geometry,atoms);}
  void WriteXYZ(std::string name, int counter = -1, int prc = 7) {
    //function to write geometries into a given xyz file
    WriteXYZFormat(name,atoms,geometry,counter,prc);
  }
  std::vector<double> CM() {
    //function to calculate the center of mass
    std::vector<double> cm = CenterOfMass(geometry,atoms);
    return cm;
  }
  std::vector<double> CMg(matrixE geom) {
    //just like the previous function but it reads the geometry from somewhere else
    //the vector of atoms is the same
    std::vector<double> cm = CenterOfMass(geom,atoms);
    return cm;
  }
  void CorrectedGeom() {
    //function to correct a geometry using the center of mass, i.e., that puts geometry in the center of the referencial
    std::vector<double> cm = CM();
    for (size_t idx = 0; idx < natoms; ++idx) {
      geometry(idx+1,1) -= cm[0];
      geometry(idx+1,2) -= cm[1];
      geometry(idx+1,3) -= cm[2];
    }
  }
  bool CompGeom(matrixE MA, matrixE MB, double tolerance) {
    //function comparing 2 geometries according to a given tolerance
    bool same = true;
    std::vector<bool> symmetric(natoms,false);
    std::vector<bool> comp(3,false);
    double diff;
    double MAX1;
    double MAX2;
    double aux;
    for (size_t irw = 0; irw < MA.rows(); ++irw) {
      for (size_t irw2 = 0; irw2 < MA.rows(); ++irw2) {
        for (size_t icl = 0; icl < MA.cols(); ++icl) {
          comp[icl] = false;
          diff = fabs(MB(irw2+1,icl+1) - MA(irw+1,icl+1));
          if (diff < tolerance) {
            if (atoms[irw2] == atoms[irw]) {
              comp[icl] = true;
            }
          }
        }
        if ((comp[0] == true)&&(comp[1] == true)&&(comp[2] == true)) {
          symmetric[irw] = true;
          break;
        }
      }
      if (symmetric[irw] == false) {
        same = false;
        break;
      }
    }
    return same;
  }
  bool Linear(double tolerance = pi/180) {
    //function to check whether a molecule is Linear
    bool linear = true;
    if (natoms > 2) {
      //can just choose 2 atoms and check whether the angle is close to 0 or 180. Search by positive then
      double angle;
      linear = false;
      for (size_t atms = 3; atms < (natoms+1); ++atms) {
        angle = calcAngles(1,2,atms);
        if ((fabs(angle) < tolerance)||(fabs(angle - pi) < tolerance)) {linear = true;}
        else {
          linear = false;
          break;
        }
      }
    }
    return linear;
  }
  void linearGeom() {
    //function that rebuilds the geometry to be linear
    matrixE xgeom(natoms,3);
    if (natoms == 2) {
      xgeom(1,1) = -1e10*distance(1,2)/2;
      xgeom(1,2) = 0.0;
      xgeom(1,3) = 0.0;
      xgeom(2,1) = 1e10*distance(1,2)/2;
      xgeom(2,2) = 0.0;
      xgeom(2,3) = 0.0;
    }
    else {
      std::vector<std::vector<size_t> > connect = Connections();
      size_t begin;
      for (size_t start = 0; start < connect.size(); ++start) {
        if (connect[start].size() == 1) {
          begin = start+1;
          break;
        }
        else {
          size_t zeroes = 0;
          for (size_t idx = 0; idx < connect[start].size(); ++idx) {
            if (connect[start][idx] == 0) {++zeroes;}
          }
          if (zeroes == (connect[start].size() - 1)) {
            begin = start+1;
            break;
          }
        }
      }
      std::cout << "begin = " << begin << std::endl;
      xgeom(begin,1) = 0.0;
      xgeom(begin,2) = 0.0;
      xgeom(begin,3) = 0.0;
      size_t connected = connect[begin-1][0];
      size_t aux;
      std::cout << "B" << std::endl;
      for (size_t idx = 0; idx < natoms; ++idx) {
        xgeom(connected,1) = 1e10*distance(begin,connected)+xgeom(begin,1);
        xgeom(connected,2) = 0.0;
        xgeom(connected,3) = 0.0;
        aux = begin;
        begin = connected;
        for (size_t idx2 = 0; idx2 < connect[begin-1].size(); ++idx2) {
          if ((connect[begin-1][idx2] != aux)&&(connect[begin-1][idx2] != 0)) {
            connected = connect[begin-1][idx2];
            break;
          }
        }
      }
    }
    geometry = xgeom;
  }
  bool Planar(double tolerance = pi/180) {
    //function to check whether a molecule is planar
    bool planar = true;
    if (natoms == 4) {
      double angle = calcDihedrals(1,2,3,4);
      if (fabs(angle) < 0.5*pi) {
        if (fabs(angle) < tolerance) {planar = true;}
        else {planar = false;}
      }
      else {
        if (fabs(fabs(angle) - pi) < tolerance) {planar = true;}
        else {planar = false;}
      }
    }
    else if (natoms == 5) {
      double angle = calcDihedrals(1,2,3,4);
      if (fabs(angle) < 0.5*pi) {
        if (fabs(angle) < tolerance) {
          //planar = true, lets check the other angles
          angle = calcDihedrals(2,3,4,5);
          if (fabs(angle) < 0.5*pi) {
            if (fabs(angle) < tolerance) {planar = true;}
            else {planar = false;}
          }
          else {
            if (fabs(fabs(angle) - pi) < tolerance) {planar = true;}
            else {planar = false;}
          }
        }
        else {planar = false;}
      }
      else {
        if (fabs(fabs(angle) - pi) < tolerance) {
          //planar = true, lets check the other angles
          angle = calcDihedrals(2,3,4,5);
          if (fabs(angle) < 0.5*pi) {
            if (fabs(angle) < tolerance) {planar = true;}
            else {planar = false;}
          }
          else {
            if (fabs(fabs(angle) - pi) < tolerance) {planar = true;}
            else {planar = false;}
          }
        }
        else {planar = false;}
      }
    }
    else if (natoms > 5) {
      //loop over atoms and create lines of atoms connected
      double angle;
      std::vector<std::vector<size_t> > connect = Connections();
      size_t atm2;
      size_t atm3;
      size_t atm4;
      size_t connsz1;
      size_t connsz2;
      size_t connsz3;
      planar = true;
      for (size_t atm1 = 0; atm1 < natoms; ++atm1) {
        connsz1 = connect[atm1].size();
        if (connsz1 > 4) {
          planar = false;
          break;
        }
        for (size_t idx2 = 0; idx2 < connsz1; ++idx2) {
          atm2 = connect[atm1][idx2];
          connsz2 = connect[atm2 - 1].size();
          if (atm2 == 0) {continue;}
          for (size_t idx3 = 0; idx3 < connsz2; ++idx3) {
            atm3 = connect[atm2 - 1][idx3];
            if (((atm3 - 1) == atm1)||(atm3 == 0)) {continue;}
            connsz3 = connect[atm3 - 1].size();
            for (size_t idx4 = 0; idx4 < connsz3; ++idx4) {
              atm4 = connect[atm3 - 1][idx4];
              if (((atm4 - 1) == atm1)||(atm4 == atm2)||(atm4 == 0)) {continue;}
              angle = calcDihedrals(atm1 + 1,atm2,atm3,atm4);
              if (fabs(angle) < 0.5*pi) {
                if (fabs(angle) < tolerance) {planar = true;}
                else {
                  planar = false;
                  angle = calcDihedrals(atm1 + 1,atm3,atm4,atm2);
                  if (fabs(angle) < tolerance) {
                    planar = true;
                    continue;
                  }
                  break;
                }
              }
              else {
                if (fabs(fabs(angle) - pi) < tolerance) {planar = true;}
                else {
                  planar = false;
                  angle = calcDihedrals(atm1 + 1,atm3,atm4,atm2);
                  if ((fabs(angle) > 0.5*pi)&&(fabs(fabs(angle) - pi) < tolerance)) {
                    planar = true;
                    continue;
                  }
                  else if ((fabs(angle) < 0.5*pi)&&(fabs(angle) < tolerance)) {
                    planar = true;
                    continue;
                  }
                  break;
                }
              }
            }
            if (planar == false) {break;}
          }
          if (planar == false) {break;}
        }
        if (planar == false) {break;}
      }
    }
    return planar;
  }
  matrixE makePlanar() {
    //molecule is planar, rebuild geometry so that z component is zero for all atoms
    //this is not easily done with a couple of rotations because the molecule may have funny orientation
    //furthermore, I need to be careful with the angles, since they do not give me an absolute direction!
    matrixE ngeom = geometry;
    //defining coordinates of first 2 atoms
    for (size_t idn = 1; idn < 4; ++idn) {
      for (size_t idn2 = 1; idn2 < natoms + 1; ++idn2) {
        ngeom(idn2,idn) = 0.0;
      }
    }
    ngeom(2,1) = 1e10*distance(1,2);
    std::vector<std::vector<size_t> > connmat = Connections();
    //this is a non-elegant workaround to get correct coordinates for all atoms of the planar geometry
    size_t con1 = 0;
    size_t con2 = 0;
    for (size_t idn1 = 0; idn1 < connmat[0].size(); ++idn1) {
      if (connmat[0][idn1] > con1) {
        //check whether the atom connected to the first atom has only one connection
        size_t pos = connmat[0][idn1];
        size_t count0 = 0;
        if (connmat[pos-1].size() == 1) {con1 = connmat[0][idn1];}
        for (size_t idn2 = 0; idn2 < connmat[pos-1].size(); ++idn2) {
          if (connmat[pos-1][idn2] == 0) {++count0;}
        }
        if (connmat[pos-1].size() == (count0+1)) {con1 = connmat[0][idn1];}
      }
    }
    for (size_t idn1 = 0; idn1 < connmat[1].size(); ++idn1) {
      if (connmat[1][idn1] > con2) {
        //check whether the atom connected to the second atom has only one connection
        size_t pos = connmat[1][idn1];
        size_t count0 = 0;
        if (connmat[pos-1].size() == 1) {con2 = connmat[1][idn1];}
        for (size_t idn2 = 0; idn2 < connmat[pos-1].size(); ++idn2) {
          if (connmat[pos-1][idn2] == 0) {++count0;}
        }
        if (connmat[pos-1].size() == (count0+1)) {con2 = connmat[1][idn1];}
      }
    }
    //defining rest of coordinates
    int sign;
    double theta;
    double aux;
    for (size_t idn = 3; idn < (natoms + 1); ++idn) {
      theta = calcAngles(1,2,idn);
      //must correct theta for atoms connected to atoms 1 and 2
      aux = 1e10*distance(1,idn);
      ngeom(idn,1) = aux*cos(theta);
      if ((idn == con1)||(idn == con2)) {sign = -1;}
      else {sign = 1;}
      ngeom(idn,2) = sign*aux*sin(theta);
      //making sure there are no colisions
      for (size_t idn2 = 1; idn2 < (idn); ++idn2) {
        aux = (ngeom(idn,2) - ngeom(idn2,2))*(ngeom(idn,2) - ngeom(idn2,2));
        aux += (ngeom(idn,1) - ngeom(idn2,1))*(ngeom(idn,1) - ngeom(idn2,1));
        if (aux < 1.0) {
          ngeom(idn,2) *= -1;
          break;
        }
      }
      ngeom(idn,3) = 0.0;
    }
    //ensuring the new planar geometry is correct
    std::vector<std::vector<size_t> > connmat2 = ConnectionsG(ngeom);
    bool comparison = true;
    int position;
    size_t missing;
    if (connmat.size() != connmat2.size()) {std::cout << "something very wrong... connmat and connmat2 do not have the same size" << std::endl;}
    else {
      for (size_t idb = 0; idb < connmat.size(); ++idb) {
        comparison = true;
        //look for discrepancies
        for (size_t idc = 0; idc < connmat[idb].size(); ++idc) {
          if (connmat[idb][idc] != connmat2[idb][idc]) {
            comparison = false;
            position = idc;
            break;
          }
        }
        size_t cnt1 = 0;
        size_t cnt2 = 0;
        //where is an atom missing?
        for (size_t idc = 0; idc < connmat[idb].size(); ++idc) {
          if (connmat[idb][idc] != 0) {++cnt1;}
          if (connmat2[idb][idc] != 0) {++cnt2;}
        }
        if (cnt1 > cnt2) {missing = 2;}
        else{missing = 1;}
        //do something about it
        if (!comparison) {
          //I assume that the misplaced atoms all have single connection
          if (missing == 2) {
            //in this situation there is an atom to add to list
            double distr = 1e10*distance(1,connmat[idb][position]);
            size_t pos;
            if (position == 0) {pos = 1;}
            else {pos = 0;}
            //only do this if there is actually an atom to make an angle with
            if (connmat[idb][pos] == 0) {continue;}
            if (connmat[idb].size() == 1) {continue;}
            double angtheta = calcAngles(1,2,connmat[idb][position]);
            ngeom(connmat[idb][position],1) = distr*cos(angtheta);
            ngeom(connmat[idb][position],2) = distr*sin(angtheta);
          }
        }
      }
    }
    std::vector<double> cm = CMg(ngeom);
    for (size_t idx = 0; idx < natoms; ++idx) {
      ngeom(idx+1,1) -= cm[0];
      ngeom(idx+1,2) -= cm[1];
    }
    return ngeom;
  }
  void Cart2Zmat(matrixE & zmat, std::vector<size_t> & atomorder) {
    //function converting a geometry in cartesian coordinates into Z-matrix
    size_t atmordersz = atomorder.size();
    if (atmordersz != natoms) {               //by allowing this, I can focus the Z-matrix on some specific atoms only, and the rest is whatever
      bool present;
      //loop over atoms
      for (size_t idatm = 1; idatm < natoms + 1; ++idatm) {
        present = false;
        //check whether atom position already in atomorder vector
        for (size_t idpos = 0; idpos < atmordersz; ++idpos) {
          if (atomorder.at(idpos) == idatm) {
            present = true;
            break;
          }
        }
        if (!present) {atomorder.push_back(idatm);}
      }
    }
    zmat.resize(natoms,6);                 //reshape it; note that atomic numbers are given separately in an array
    zmat.zero();
    zmat(2,1) = 1.0;
    std::vector<double> dij = this->calcUDistances();
    zmat(2,2) = dij[LTriangMat2Array_nodiag(atomorder.at(0),atomorder.at(1),natoms)];
    if (natoms > 2) {
      zmat(3,1) = 2.0;
      zmat(3,2) = dij[LTriangMat2Array_nodiag(atomorder.at(1),atomorder.at(2),natoms)];
      zmat(3,3) = 1.0;
      zmat(3,4) = calcAngles(atomorder.at(1),atomorder.at(2),atomorder.at(0));
      if (natoms > 3) {
        for (size_t idatm = 3; idatm < natoms; ++idatm) {
          zmat(idatm + 1,1) = double(idatm);
          zmat(idatm + 1,2) = dij[LTriangMat2Array_nodiag(atomorder.at(idatm),atomorder.at(idatm - 1),natoms)];
          zmat(idatm + 1,3) = double(idatm - 1);
          zmat(idatm + 1,4) = calcAngles(atomorder.at(idatm - 1),atomorder.at(idatm),atomorder.at(idatm - 2));
          zmat(idatm + 1,5) = double(idatm - 2);
          zmat(idatm + 1,6) = calcDihedrals(atomorder.at(idatm),atomorder.at(idatm - 1),atomorder.at(idatm - 2),atomorder.at(idatm - 3));
        }
      }
    }
  }
  void Zmat2Cart(matrixE & zmat) {
    //function converting a geometry in Z-matrix into cartesian coordinates; based on TMPChem's python code
    natoms = zmat.rows();
    double rij;
    double theta;
    double factor;
    size_t atmR;
    size_t atmA;
    size_t atmD;
    std::vector<double> u21(3,0.0);
    std::vector<double> u23(3,0.0);
    std::vector<double> cp(3,0.0);
    std::vector<double> bv(3,0.0);
    matrixE cart(natoms,3);
    cart(2,3) = zmat(2,2);
    if (natoms > 1) {
      if (natoms > 2) {
        factor = 1.0;
        if (int(zmat(3,1)) == 2) {factor = -1.0;}      //if atom 3 is bonded to atom 2
        rij = zmat(3,2);
        theta = zmat(3,4);
        cart(3,2) = rij*sin(theta);
        cart(3,3) = cart(2,3) + factor*rij*cos(theta);
        for (size_t idatm = 4; idatm < natoms + 1; ++idatm) {
          //get coordinates of atoms defining the position; in my z-matrix construction this step is actually not needed; I keep this so that I am as general as possible
          atmR = size_t(zmat(idatm,1));
          atmA = size_t(zmat(idatm,3));
          atmD = size_t(zmat(idatm,5));
          GetUnitVectors(u21,u23,cart,atmA,atmR,atmD);
          GetLocalAxes(u21,u23,cp);
          GetBondVector(bv,zmat,idatm);
          //get atomic coordinates
          cart(idatm,1) = cart(atmR,1) + dotProd(bv,cp);
          cart(idatm,2) = cart(atmR,2) + dotProd(bv,u23);
          cart(idatm,3) = cart(atmR,3) + dotProd(bv,u21);
        }
      }
    }
    geometry = cart;
  }
//--------------------------------- Connectivity Functions
  int MaxCoordinationNumber() {return maxcoordinationnumber;}
  void setMaxCoordinationNumber(int _maxcdn) {maxcoordinationnumber = _maxcdn;}
  void ConnectivityMatrix(matrix<int> & ConnMat, double factor = 1.3) {
    //function that calculates connectivity between atoms as a matrix
    //factor is a scaling factor to apply on vdW radii
    double rAB2;
    double aux;
    double rABcov;
    std::vector<int> count(natoms,1);
    if ((ConnMat.rows() != natoms)&&(ConnMat.cols() != maxcoordinationnumber)) {ConnMat.resize(natoms,maxcoordinationnumber);}
    for (size_t idAtm = 0; idAtm < natoms; ++idAtm) {
      for (size_t idpos = 1; idpos < maxcoordinationnumber + 1; ++idpos) {
        ConnMat(idAtm + 1,idpos) = 0;
      }
      for (size_t idBtm = 0; idBtm < idAtm; ++idBtm) {
        rAB2 = 0.0;
        for (size_t idcoord = 1; idcoord < 4; ++idcoord) {
          aux = geometry(idAtm + 1,idcoord) - geometry(idBtm + 1,idcoord);
          rAB2 += aux*aux;
        }
        rABcov = factor*(AtmRadii(atoms[idAtm]) + AtmRadii(atoms[idBtm]));
        if (rAB2 < rABcov*rABcov) {                     //neighbour found
          ConnMat(idAtm + 1,count[idAtm]) = idBtm + 1;
          ++count[idAtm];
          ConnMat(idBtm + 1,count[idBtm]) = idAtm + 1;
          ++count[idBtm];
        }
      }
    }
  }
  bool CheckBondDistances(double minbondlen = 1.0, bool quick = false, double scalefactor = 1.3) {
    //function to check bond distances; basically whether there is a large overlap between atoms
    //boolean quick is to determine which algorithm to take
    //quick -> check only immediate neighbours
    //!quick -> loop over all atoms
    double RAB;
    bool problem = false;
    if (quick) {
      matrix<int> neighbours;
      ConnectivityMatrix(neighbours,scalefactor);
      size_t totalnghb = neighbours.cols();
      for (size_t idatm = 0; idatm < natoms; ++idatm) {
        for (size_t idneigh = 0; idneigh < totalnghb; ++idneigh) {
          if (neighbours(idatm + 1,idneigh + 1) - 1 < 0) {break;}
          RAB = Distance(idatm + 1,neighbours(idatm + 1,idneigh + 1),geometry);
          if (RAB < minbondlen) {
            problem = true;
            std::cout << "potential abnormal bond: " << idatm + 1 << "," << neighbours(idatm + 1,idneigh + 1) << "\n";
          }
        }
      }
    }
    else {
      for (size_t idatm = 0; idatm < natoms; ++idatm) {
        for (size_t idbtm = idatm + 1; idbtm < natoms; ++idbtm) {
          RAB = Distance(idatm + 1,idbtm + 1,geometry);
          if (RAB < minbondlen) {
            problem = true;
            std::cout << "potential abnormal bond: " << idatm + 1 << "," << idbtm + 1 << "\n";
          }
        }
      }
    }
    return problem;
  }
  void CheckSystem(double scalefactor = 1.25, double mindist = 1.0, bool quick = false) {
    //function that checks for particular chemical functionalities and potential inconsistencies
    matrix<int> neighbours;
    ConnectivityMatrix(neighbours,scalefactor);
    std::vector<std::string> functionalgroups = CheckGroupsByNeighborhood(atoms,neighbours,geometry);
    std::cout << "functional groups in system" << std::endl;
    for (size_t idx = 0; idx < natoms; ++idx) {
      std::cout << functionalgroups[idx] << std::endl;
    }
    std::cout << "----------------" << std::endl;
    std::cout << "running check on bond distances " << std::endl;
    bool problem = CheckBondDistances(mindist,quick,scalefactor);
    std::cout << "problems found: " << problem << std::endl;
    int totalcharge = TotalChargeEstimator(functionalgroups,atoms,neighbours);
    std::cout << "estimated total charge: " << totalcharge << std::endl;
  }
  void AnalyseNeighbours(double scalefactor = 1.3) {
    //function that analyses a system by checking it for correct number of neighbours
    matrix<int> neighbours;
    ConnectivityMatrix(neighbours,scalefactor);
    size_t totalnghb = neighbours.cols();
    size_t atomA;
    size_t nneigh;
    size_t auxie;
    //now that I have the connectivity matrix, lets check for consistency in neighbours
    for (size_t idatm = 0; idatm < natoms; ++idatm) {
      atomA = atoms[idatm];
      nneigh = 0;
      for (size_t idneigh = 0; idneigh < totalnghb; ++idneigh) {
        if (neighbours(idatm + 1,idneigh + 1) > 0) {++nneigh;}
        else {break;}
      }
      switch (atomA) {
        case 1:
          if (nneigh > 1) {std::cout << "abnormal atom H (too many bonds): " << idatm + 1 << "\n";}
          break;
        case 2:
          break;
        case 3:
          if (nneigh > 1) {std::cout << "abnormal atom Li (too many bonds): " << idatm + 1 << "\n";}
          break;
        case 4:
          break;
        case 5:
          if (nneigh < 2) {std::cout << "abnormal atom B (too few bonds): " << idatm + 1 << "\n";}
          else if (nneigh > 4) {std::cout << "abnormal atom B (too many bonds): " << idatm + 1 << "\n";}
          break;
        case 6:
          if (nneigh < 2) {std::cout << "abnormal atom C (too few bonds): " << idatm + 1 << "\n";}
          else if (nneigh > 4) {std::cout << "abnormal atom C (too many bonds): " << idatm + 1 << "\n";}
          break;
        case 7:
          if (nneigh > 4) {std::cout << "abnormal atom N (too many bonds): " << idatm + 1 << "\n";}
          break;
        case 8:
          if (nneigh > 3) {std::cout << "abnormal atom O (too many bonds): " << idatm + 1 << "\n";}
          break;
        case 9:
          if (nneigh > 1) {std::cout << "abnormal atom F (too many bonds): " << idatm + 1 << "\n";}
          break;
        case 10:
          break;
        case 11:
          if (nneigh > 1) {std::cout << "abnormal atom Na (too many bonds): " << idatm + 1 << "\n";}
          break;
        case 12:
          if ((nneigh != 2)&&(nneigh != 4)) {std::cout << "abnormal atom Mg: " << idatm + 1 << "\n";}
          break;
        case 13:
          if (nneigh != 3) {std::cout << "abnormal atom Al: " << idatm + 1 << "\n";}
          break;
        case 14:
          break;
        case 15:
          if ((nneigh != 3)&&(nneigh != 4)&&(nneigh != 5)) {std::cout << "abnormal atom P: " << idatm + 1 << "\n";}
          if (nneigh == 3) {
            //count oxygens
            auxie = 0;
            for (size_t idneigh = 0; idneigh < totalnghb; ++idneigh) {
              if (atoms[neighbours(idatm + 1,idneigh + 1) - 1] == 8) {++auxie;}
              else {break;}
            }
            if (auxie == 3) {std::cout << "potential problem on atom P: " << idatm + 1 << ". Coordination number 3 to 3 oxygens -> uncommon.\n";}
          }
          break;
        case 16:
          if (nneigh > 4) {std::cout << "abnormal atom S: " << idatm + 1 << "\n";}
          break;
        case 17:
          break;
        case 18:
          break;
        case 19:
          if (nneigh > 1) {std::cout << "abnormal atom K (too many bonds): " << idatm + 1 << "\n";}
          break;
        case 20:
          break;
        case 21:
          break;
        case 22:
          break;
        case 23:
          break;
        case 24:
          break;
        case 25:
          break;
        case 26:
          break;
        case 27:
          break;
        case 28:
          break;
        case 29:
          break;
        case 30:
          break;
        case 31:
          break;
        case 32:
          break;
        case 33:
          break;
        case 34:
          break;
        case 35:
          break;
        case 36:
          break;
        case 37:
          if (nneigh > 1) {std::cout << "abnormal atom Rb (too many bonds): " << idatm + 1 << "\n";}
          break;
        case 38:
          break;
        case 39:
          break;
        case 40:
          break;
        case 41:
          break;
        case 42:
          break;
        case 43:
          break;
        case 44:
          break;
        case 45:
          break;
        case 46:
          if ((nneigh != 4)&&(nneigh != 6)) {std::cout << "abnormal atom Pd: " << idatm + 1 << "\n";}
          break;
        case 47:
          break;
        case 48:
          break;
        case 49:
          break;
        case 50:
          break;
        case 51:
          break;
        case 52:
          break;
        case 53:
          break;
        case 54:
          break;
        case 55:
          if (nneigh > 1) {std::cout << "abnormal atom Cs (too many bonds): " << idatm + 1 << "\n";}
          break;
        case 56:
          break;
        case 57:
          break;
        case 58:
          break;
        case 59:
          break;
        case 60:
          break;
        case 61:
          break;
        case 62:
          break;
        case 63:
          break;
        case 64:
          break;
        case 65:
          break;
        case 66:
          break;
        case 67:
          break;
        case 68:
          break;
        case 69:
          break;
        case 70:
          break;
        case 71:
          break;
        case 72:
          break;
        case 73:
          break;
        case 74:
          break;
        case 75:
          break;
        case 76:
          break;
        case 77:
          break;
        case 78:
          if ((nneigh != 4)&&(nneigh != 6)) {std::cout << "abnormal atom Pt: " << idatm + 1 << "\n";}
          break;
        case 79:
          break;
        case 80:
          break;
        case 81:
          break;
        case 82:
          break;
        case 83:
          break;
        case 84:
          break;
        case 85:
          break;
        case 86:
          break;
        case 87:
          break;
        case 88:
          break;
        case 89:
          break;
        case 90:
          break;
        case 91:
          break;
        case 92:
          break;
        case 93:
          break;
        case 94:
          break;
      }
    }
  }
  std::vector<std::vector<size_t> > Connections() {
    //function that returns a vector of vectors
    //the outer vector contains the atoms as they appear in the geometry
    //the inner vector contains the numbers of the atoms a given atom is connected to
    matrixE connectivity = calcDistances();
    std::vector<std::vector<size_t> > connections;
    std::vector<double> dist;
    std::vector<size_t> conns;
    size_t _octet;
    double distmol;                     //to write the distance between atoms in molecule
    double maxbonddist;                 //to write the maximum expected bond distance between atoms
    //loop over all atoms
    for (size_t idx = 0; idx < natoms; ++idx) {
      //establishing the maximum number of bonds per atom
      _octet = size_t(octet(atoms[idx]));
      //std::cout << "octet atom " << atoms[idx] << " " << octet(atoms[idx]) << std::endl;
      //std::cout << "Atom " << idx+1 << " octet " << _octet << std::endl;
      dist.resize(_octet);
      conns.resize(_octet);
      //setting an upper bound for first search. Distances are in Angstrom
      for (size_t idx2 = 0; idx2 < _octet; ++idx2) {
        dist[idx2] = 100.0;          //do not really know anymore why I set this upper limit, since 
        conns[idx2] = 0;
      }
      //loop over other atoms to check which are connected to atom idx+1
      for (size_t idx2 = 0; idx2 < natoms; ++idx2) {
        if (idx2 != idx) {
          distmol = connectivity(idx+1,idx2+1)/dist_si2Angstrom;
          maxbonddist = MaxBondDistance(atoms[idx],atoms[idx2],1.30);
          //std::cout << "distance between " << AtomNr2Symbol(atoms[idx]) << " and " << AtomNr2Symbol(atoms[idx2]) << " is ";
          //std::cout << distmol << "; max dist is " << maxbonddist << std::endl;
          //std::cout << _octet << std::endl;
          if ((distmol < dist[_octet-1])&&(distmol < maxbonddist)) {
            for (size_t idx3 = 0; idx3 < _octet; ++idx3 ) {
              if (distmol <= dist[idx3]) {
                size_t aux;
                for (size_t idx4 = 0; idx4 < _octet; ++idx4) {
                  aux = _octet - idx4 - 1;
                  if (aux > idx3) {
                    dist[aux] = dist[aux-1];
                    conns[aux] = conns[aux-1];
                  }
                }
                dist[idx3] = distmol;
                conns[idx3] = idx2+1;
                break;
              }
            }
          }
        }
      }
      connections.push_back(conns);
    }
    return connections;
  }
  std::vector<std::vector<size_t> > ConnectionsG(matrixE geom) {
    //the general Connections function. Just like the previous function but also accepts another geometry.
    matrixE connectivity = calcDistancesG(geom);
    std::vector<std::vector<size_t> > connections;
    std::vector<double> dist;
    std::vector<size_t> conns;
    size_t _octet;
    double distmol;                     //to write the distance between atoms in molecule
    double maxbonddist;                 //to write the maximum expected bond distance between atoms.
    //loop over all atoms
    for (size_t idx = 0; idx < natoms; ++idx) {
      //establishing the maximum number of bonds per atom
      _octet = size_t(octet(atoms[idx]));
      //std::cout << "Atom " << idx+1 << " octet " << _octet << std::endl;
      dist.resize(_octet);
      conns.resize(_octet);
      //setting an upper bound for first search. Distances are in meters
      for (size_t idx2 = 0; idx2 < _octet; ++idx2) {
        dist[idx2] = 1.0;
        conns[idx2] = 0.0;
      }
      //loop over other atoms to check which are connected to atom idx+1
      for (size_t idx2 = 0; idx2 < natoms; ++idx2) {
        if (idx2 != idx) {
          distmol = connectivity(idx+1,idx2+1);
          maxbonddist = MaxBondDistance(atoms[idx],atoms[idx2],1.25);
          //std::cout << "distance between " << AtomNr2Symbol(atoms[idx]) << " and " << AtomNr2Symbol(atoms[idx2]) << " is ";
          //std::cout << distmol << "; max dist is " << maxbonddist << std::endl;
          //std::cout << _octet << std::endl;
          if ((distmol < dist[_octet-1])&&(distmol < maxbonddist)) {
            for (size_t idx3 = 0; idx3 < _octet; ++idx3 ) {
              if (distmol <= dist[idx3]) {
                size_t aux;
                for (size_t idx4 = 0; idx4 < _octet; ++idx4) {
                  aux = _octet - idx4 - 1;
                  if (aux > idx3) {
                    dist[aux] = dist[aux-1];
                    conns[aux] = conns[aux-1];
                  }
                }
                dist[idx3] = distmol;
                conns[idx3] = idx2+1;
                break;
              }
            }
          }
        }
      }
      connections.push_back(conns);
    }
    return connections;
  }
//--------------------------------- distance functions
  double MaxBondDistance(int atom1, int atom2, double fac = 1.2) {
    //function to determine the maximum reasonable bond distance between 2 atoms
    double rad1 = RadiusPsi4(size_t(atom1));
    double rad2 = RadiusPsi4(size_t(atom2));
    return fac*(rad1 + rad2);
  }
  std::vector<double> calcUDistances() {
    //function that determines all bond distances and stores them in an array
    std::vector<double> distances(natoms*(natoms - 1)/2,0.0);
    double distance;
    double distance2;
    for (size_t idatm1 = 1; idatm1 < natoms; ++idatm1) {
      for (size_t idatm2 = 0; idatm2 < idatm1; ++idatm2) {
        distance2 = (geometry(idatm1 + 1,1)-geometry(idatm2 + 1,1))*(geometry(idatm1 + 1,1)-geometry(idatm2 + 1,1));
        distance2 += (geometry(idatm1 + 1,2)-geometry(idatm2 + 1,2))*(geometry(idatm1 + 1,2)-geometry(idatm2 + 1,2));
        distance2 += (geometry(idatm1 + 1,3)-geometry(idatm2 + 1,3))*(geometry(idatm1 + 1,3)-geometry(idatm2 + 1,3));
        distance = sqrt(distance2);
        distances[LTriangMat2Array_nodiag(idatm1 + 1,idatm2 + 1,natoms)] = distance;
      }
    }
    return distances;
  }
  matrixE calcDistances() {
    //function that determines a matrix of all bond distances
    matrixE distances(natoms,natoms);
    double distance;
    double distance2;
    for (size_t idx = 1; idx < natoms; ++idx) {
      for (size_t idx2 = 0; idx2 < idx; ++idx2) {
        distance2 = (geometry(idx+1,1)-geometry(idx2+1,1))*(geometry(idx+1,1)-geometry(idx2+1,1));
        distance2 += (geometry(idx+1,2)-geometry(idx2+1,2))*(geometry(idx+1,2)-geometry(idx2+1,2));
        distance2 += (geometry(idx+1,3)-geometry(idx2+1,3))*(geometry(idx+1,3)-geometry(idx2+1,3));
        distance = sqrt(distance2);
        distances(idx+1,idx2+1) = distance*1e-10;
        distances(idx2+1,idx+1) = distance*1e-10;
      }
    }
    return distances;
  }
  matrixE calcDistancesG(matrixE geom) {
    //generalization of the previous function to other geometries
    matrixE distances(natoms,natoms);
    double distance;
    double distance2;
    for (size_t idx = 1; idx < natoms; ++idx) {
      for (size_t idx2 = 0; idx2 < idx; ++idx2) {
        distance2 = (geom(idx+1,1)-geom(idx2+1,1))*(geom(idx+1,1)-geom(idx2+1,1));
        distance2 += (geom(idx+1,2)-geom(idx2+1,2))*(geom(idx+1,2)-geom(idx2+1,2));
        distance2 += (geom(idx+1,3)-geom(idx2+1,3))*(geom(idx+1,3)-geom(idx2+1,3));
        distance = sqrt(distance2);
        distances(idx+1,idx2+1) = distance*1e-10;
        distances(idx2+1,idx+1) = distance*1e-10;
      }
    }
    return distances;
  }
  double Adistance(size_t el1, size_t el2) {
    //function to calculate the distance between 2 atoms (in Angstrom)
    double _dist = (geometry(el1,1)-geometry(el2,1))*(geometry(el1,1)-geometry(el2,1));
    _dist += (geometry(el1,2)-geometry(el2,2))*(geometry(el1,2)-geometry(el2,2));
    _dist += (geometry(el1,3)-geometry(el2,3))*(geometry(el1,3)-geometry(el2,3));
    return sqrt(_dist);
  }
  double AUdistance(size_t el1, size_t el2) {
    //function to calculate the distance between 2 atoms (in Bohr)
    return Adistance(el1,el2)/dist_Angstrom2au;
  }
  double distance(size_t el1, size_t el2) {
    //function to calculate the distance between 2 atoms (in SI)
    return dist_si2Angstrom*Adistance(el1,el2);
  }
  double distanceG(matrixE & geom, size_t el1, size_t el2) {
    //like the previous function but generallized to any geometry
    double _dist = (geom(el1,1)-geom(el2,1))*(geom(el1,1)-geom(el2,1));
    _dist += (geom(el1,2)-geom(el2,2))*(geom(el1,2)-geom(el2,2));
    _dist += (geom(el1,3)-geom(el2,3))*(geom(el1,3)-geom(el2,3));
    return dist_si2Angstrom*sqrt(_dist);
  }
//--------------------------------- Angle Functions
  std::vector<matrixE> calcAngles() {
    //function that calculates all molecular angles
    std::vector<matrixE> angles;
    matrixE anglesm(natoms,natoms);
    std::vector<double> aref(3,0.0);
    std::vector<double> a1(3,0.0);
    std::vector<double> a2(3,0.0);
    double intprod;
    double norm1;
    double norm2;
    //loop over reference atom. This is unique and has the simmetry of permutation a1 with a2
    for (size_t idxr = 0; idxr < natoms; ++idxr) {
      aref[0] = geometry(idxr+1,1);
      aref[1] = geometry(idxr+1,2);
      aref[2] = geometry(idxr+1,3);
      //loop over a1
      for (size_t idx1 = 0; idx1 < natoms; ++idx1) {
        if (idxr != idx1) {
          //getting a1 with respect to reference
          a1[0] = geometry(idx1+1,1) - aref[0];
          a1[1] = geometry(idx1+1,2) - aref[1];
          a1[2] = geometry(idx1+1,3) - aref[2];
          norm1 = sqrt(a1[0]*a1[0] + a1[1]*a1[1] + a1[2]*a1[2]);
          //loop over a2
          for (size_t idx2 = idx1; idx2 < natoms; ++idx2) {
            if ((idxr != idx2)&&(idx1 != idx2)) {
              //getting a2 with respect to reference
              a2[0] = geometry(idx2+1,1) - aref[0];
              a2[1] = geometry(idx2+1,2) - aref[1];
              a2[2] = geometry(idx2+1,3) - aref[2];
              intprod = a1[0]*a2[0] + a1[1]*a2[1] + a1[2]*a2[2];
              norm2 = sqrt(a2[0]*a2[0] + a2[1]*a2[1] + a2[2]*a2[2]);
              anglesm(idx1+1,idx2+1) = acos((intprod)/(norm1*norm2));
              //anglesm(idx1+1,idx2+1) = acos((intprod)/(norm1*norm2))*180/pi;
              anglesm(idx2+1,idx1+1) = anglesm(idx1+1,idx2+1);
            }
            else {
              anglesm(idx1+1,idx2+1) = 0;
              anglesm(idx2+1,idx1+1) = 0;
            }
          }
        }
        else {
          for (size_t idx2 = 0; idx2 < natoms; ++idx2) {
            anglesm(idx2+1,idx1+1) = 0;
            anglesm(idx1+1,idx2+1) = 0;
          }
        }
      }
      angles.push_back(anglesm);
      //std::cout << "The matrix of angles for the reference atom " << idxr+1 << std::endl;
      //anglesm.Print();
    }
    return angles;
  }
  matrixE calcAngles(size_t refa) {
    //function calculating all angles for a fixed reference atom. This atom is used for the origin of referencial
    matrixE anglesm(natoms,natoms);
    std::vector<double> aref(3,0.0);
    std::vector<double> a1(3,0.0);
    std::vector<double> a2(3,0.0);
    double intprod;
    double norm1;
    double norm2;
    aref[0] = geometry(refa,1);
    aref[1] = geometry(refa,2);
    aref[2] = geometry(refa,3);
    //loop over a1
    for (size_t idx1 = 0; idx1 < natoms; ++idx1) {
      if (idx1 != (refa-1)) {
        //getting a1 with respect to reference
        a1[0] = geometry(idx1+1,1) - aref[0];
        a1[1] = geometry(idx1+1,2) - aref[1];
        a1[2] = geometry(idx1+1,3) - aref[2];
        norm1 = sqrt(a1[0]*a1[0] + a1[1]*a1[1] + a1[2]*a1[2]);
        //loop over a2
        for (size_t idx2 = 0; idx2 < natoms; ++idx2) {
          if ((idx2 != (refa-1))&&(idx1 != idx2)) {
            //getting a2 with respect to reference
            a2[0] = geometry(idx2+1,1) - aref[0];
            a2[1] = geometry(idx2+1,2) - aref[1];
            a2[2] = geometry(idx2+1,3) - aref[2];
            intprod = a1[0]*a2[0] + a1[1]*a2[1] + a1[2]*a2[2];
            norm2 = sqrt(a2[0]*a2[0] + a2[1]*a2[1] + a2[2]*a2[2]);
            anglesm(idx1+1,idx2+1) = acos((intprod)/(norm1*norm2));
            //anglesm(idx1+1,idx2+1) = acos((intprod)/(norm1*norm2))*180/pi;
            anglesm(idx2+1,idx1+1) = anglesm(idx1+1,idx2+1);
          }
          else {
            anglesm(idx1+1,idx2+1) = 0;
            anglesm(idx2+1,idx1+1) = 0;
          }
        }
      }
      else {
        for (size_t idx2 = 0; idx2 < natoms; ++idx2) {
          anglesm(idx2+1,idx1+1) = 0;
        }
      }
    }
    return anglesm;
  }
  double calcAngles(size_t refa, size_t at1, size_t at2) {
    //function that calculates the angle made by a set of 3 atoms.
    double anglesm = 0.0;
    if ((refa != at1)&&(refa != at2)&&(at2 != at1)) {
      std::vector<double> aref(3,0.0);
      std::vector<double> a1(3,0.0);
      std::vector<double> a2(3,0.0);
      double intprod;
      double norm1;
      double norm2;
      aref[0] = geometry(refa,1);
      aref[1] = geometry(refa,2);
      aref[2] = geometry(refa,3);
      a1[0]   = geometry(at1,1) - aref[0];
      a1[1]   = geometry(at1,2) - aref[1];
      a1[2]   = geometry(at1,3) - aref[2];
      norm1   = sqrt(a1[0]*a1[0] + a1[1]*a1[1] + a1[2]*a1[2]);
      a2[0]   = geometry(at2,1) - aref[0];
      a2[1]   = geometry(at2,2) - aref[1];
      a2[2]   = geometry(at2,3) - aref[2];
      intprod = a1[0]*a2[0] + a1[1]*a2[1] + a1[2]*a2[2];
      norm2   = sqrt(a2[0]*a2[0] + a2[1]*a2[1] + a2[2]*a2[2]);
      anglesm = acos((intprod)/(norm1*norm2));
    }
    return anglesm;
  }
  double calcDihedrals(size_t at1, size_t at2, size_t at3, size_t at4, double threshold = 1.0e-7) {
    //function calculating the dihedral angle between 4 atoms; adapted from TMPChem
    double dihedral = 0.0;
    if ((at1 != at2)&&(at1 != at3)&&(at1 != at4)&&(at2 != at3)&&(at2 != at4)&&(at3 != at4)) {
      //calculate vectors for planes
      std::vector<double> atomA(3,0.0);
      std::vector<double> atomB(3,0.0);
      std::vector<double> atomC(3,0.0);
      std::vector<double> atomD(3,0.0);
      for (size_t idx = 1; idx < 4; ++idx) {
        atomA[idx - 1] = geometry(at1,idx);
        atomB[idx - 1] = geometry(at2,idx);
        atomC[idx - 1] = geometry(at3,idx);
        atomD[idx - 1] = geometry(at4,idx);
      }
      std::vector<double> a1(3,0.0);
      std::vector<double> a2(3,0.0);
      a1[0] = (atomB[1] - atomA[1])*(atomC[2] - atomA[2]) - (atomB[2] - atomA[2])*(atomC[1] - atomA[1]);
      a1[1] = (atomB[2] - atomA[2])*(atomC[0] - atomA[0]) - (atomB[0] - atomA[0])*(atomC[2] - atomA[2]);
      a1[2] = (atomB[0] - atomA[0])*(atomC[1] - atomA[1]) - (atomB[1] - atomA[1])*(atomC[0] - atomA[0]);
      a2[0] = (atomC[1] - atomB[1])*(atomD[2] - atomB[2]) - (atomC[2] - atomB[2])*(atomD[1] - atomB[1]);
      a2[1] = (atomC[2] - atomB[2])*(atomD[0] - atomB[0]) - (atomC[0] - atomB[0])*(atomD[2] - atomB[2]);
      a2[2] = (atomC[0] - atomB[0])*(atomD[1] - atomB[1]) - (atomC[1] - atomB[1])*(atomD[0] - atomB[0]);
      //calculate angle
      double n1 = sqrt(a1[0]*a1[0] + a1[1]*a1[1] + a1[2]*a1[2]);
      double n2 = sqrt(a2[0]*a2[0] + a2[1]*a2[1] + a2[2]*a2[2]);
      double intprod = a1[0]*a2[0] + a1[1]*a2[1] + a1[2]*a2[2];
      double theta = (intprod)/(n1*n2);
      if (theta > 1.0) {theta = 1.0;}
      else if (theta < -1.0) {theta = -1.0;}
      dihedral = acos(theta);
      //determine sign
      std::vector<double> u21(3,0.0);
      std::vector<double> * pu21 = & u21;
      std::vector<double> u23(3,0.0);
      std::vector<double> * pu23 = & u23;
      std::vector<double> u34(3,0.0);
      std::vector<double> * pu34 = & u34;
      double nu21 = 0.0;
      double nu23 = 0.0;
      double nu34 = 0.0;
      double aux;
      //get unit vectors between atoms
      for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
        u21[idcoord] = geometry(at1,idcoord + 1) - geometry(at2,idcoord + 1);
        nu21 += u21[idcoord]*u21[idcoord];
        u23[idcoord] = geometry(at3,idcoord + 1) - geometry(at2,idcoord + 1);
        nu23 += u23[idcoord]*u23[idcoord];
        u34[idcoord] = geometry(at4,idcoord + 1) - geometry(at3,idcoord + 1);
        nu34 += u34[idcoord]*u34[idcoord];
      }
      for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
        u21[idcoord] /= sqrt(nu21);
        u23[idcoord] /= sqrt(nu23);
        u34[idcoord] /= sqrt(nu34);
      }
      aux = dotProd(pu23,pu21);
      std::vector<double> u21xu23 = crossProd(u21,u23,sqrt(1.0 - aux*aux));
      std::vector<double> * pu21xu23 = & u21xu23;
      double sign = 1.0;
      if (dotProd(pu21xu23,pu34) > threshold) {sign = -1.0;}
      dihedral *= sign;
    }
    return dihedral;
  }
//--------------------------------- Mass Functions
  void masses() {
    //function to determine the molecular masses (mol and single molecule)
    MWeight = 0.0;
    for (size_t idx = 0; idx < natoms; ++idx) {
      MWeight += Weight(atoms[idx]);
    }
    mass = MWeight/NA;
  }
  double AtomicMasses(std::vector<double> & amass, bool convert2au = false, double scaleproton = 1.0) {
    //function that populates a vector with atomic masses
    if (amass.size() != natoms) {amass.resize(natoms);}
    double totalmass = 0.0;
    double aux;
    double unitconv = 1.0;
    if (convert2au) {unitconv = 1.0/(au2kg*NA);}         //convert to atomic units
    for (size_t idAtm = 0; idAtm < natoms; ++idAtm) {
      aux = Weight(atoms[idAtm])*unitconv;
      totalmass += aux;
      if (atoms[idAtm] == 1) {aux *= scaleproton;}
      amass[idAtm] = aux;
    }
    return totalmass;
  }
  double MW() {return MWeight;}
  double MW() const {return MWeight;}
  double m() {return mass;}
  double m() const {return mass;}
//--------------------------------- Inertia Functions
  std::vector<double> InertiaEigenvalues(bool inmeter = true) {
    //function calculating the eigenvalues of the molecular inertia matrix
    matrixE _inert = InertiaMatrix(geometry,atoms,inmeter);
    std::vector<double> eigenvalues = MatDiag(_inert);
    return eigenvalues;
  }
  matrixE Inertia(bool inmeter = true, bool scale = false) {
    //function calculating the molecular inertia matrix
    matrixE _inert = InertiaMatrix(geometry,atoms,inmeter);
    if (scale) {_inert *= 1.0e44;}
    return _inert;
  }
//--------------------------------- "particle" counting functions
  void CountElectrons() {
    //function to calculate the number of electrons
    nelectrons = 0;
    for (size_t idx = 0; idx < natoms; ++idx) {
      nelectrons += atoms[idx];
    }
    nelectrons -= charge;
  }
  size_t Natoms() {return natoms;}
  size_t Natoms() const {return natoms;}
  void calcNatoms() {natoms = atoms.size();}
  int Nelectrons() {return nelectrons;}
  int Nelectrons() const {return nelectrons;}
  size_t Multiplicity() {return mult;}
  size_t Multiplicity() const {return mult;}
  void setMultiplicity(size_t _mult) {mult = _mult;}
  int Charge() {return charge;}
  void setCharge(int _charge) {
    charge = _charge;
    CountElectrons();
  }
  int Charge() const {return charge;}
//--------------------------------- Vector Related Functions
  std::vector<double> oVector2Atom(matrixE gm, size_t atm1, size_t atm2) {
    //function that returns the vector orthogonal to the plane made by 2 atoms and the origin
    //atmi are atom numbers (start on 1)
    std::vector<double> atomB(3,0.0);
    std::vector<double> atomC(3,0.0);
    for (size_t idx = 1; idx < 4; ++idx) {
      atomB[idx-1] = gm(atm1,idx);
      atomC[idx-1] = gm(atm2,idx);
    }
    std::vector<double> a1(3,0.0);
    a1[0] = atomB[1]*atomC[2] - atomB[2]*atomC[1];
    a1[1] = atomB[2]*atomC[0] - atomB[0]*atomC[2];
    a1[2] = atomB[0]*atomC[1] - atomB[1]*atomC[0];
    double norm = sqrt(a1[0]*a1[0] + a1[1]*a1[1] + a1[2]*a1[2]);
    a1[0] /= norm;
    a1[1] /= norm;
    a1[2] /= norm;
    return a1;
  }
  std::vector<double> oVector3Atom(matrixE gm, size_t atm1, size_t atm2, size_t atm3) {
    //function that returns the vector orthogonal to the plane made by 3 atoms
    //atmi are atom numbers (start on 1)
    std::vector<double> atomA(3,0.0);
    std::vector<double> atomB(3,0.0);
    std::vector<double> atomC(3,0.0);
    for (size_t idx = 1; idx < 4; ++idx) {
      atomA[idx-1] = gm(atm1,idx);
      atomB[idx-1] = gm(atm2,idx);
      atomC[idx-1] = gm(atm3,idx);
    }
    std::vector<double> a1(3,0.0);
    a1[0] = (atomB[1]-atomA[1])*(atomC[2]-atomA[2]) - (atomB[2]-atomA[2])*(atomC[1]-atomA[1]);
    a1[1] = (atomB[2]-atomA[2])*(atomC[0]-atomA[0]) - (atomB[0]-atomA[0])*(atomC[2]-atomA[2]);
    a1[2] = (atomB[0]-atomA[0])*(atomC[1]-atomA[1]) - (atomB[1]-atomA[1])*(atomC[0]-atomA[0]);
    double norm = sqrt(a1[0]*a1[0] + a1[1]*a1[1] + a1[2]*a1[2]);
    a1[0] /= norm;
    a1[1] /= norm;
    a1[2] /= norm;
    return a1;
  }
  std::vector<double> oVectorPlane(double tolerance = pi/180) {
    //function that returns the vector orthogonal to the molecular plane
    double angle;
    bool planar = false;
    bool swtch = false;
    size_t atm1;
    size_t atm2;
    size_t atm3;
    size_t atm4;
    if (natoms == 3) {
      atm1 = 1;
      atm2 = 2;
      atm3 = 3;
    }
    else {
      //loop over atoms and create lines of atoms connected
      for (size_t idx1 = 0; idx1 < natoms; ++idx1) {
        atm1 = idx1+1;
        for (size_t idx2 = 0; idx2 < natoms; ++idx2) {
          atm2 = idx2+1;
          if (atm2 == atm1) {continue;}
          for (size_t idx3 = 0; idx3 < natoms; ++idx3) {
            atm3 = idx3+1;
            if ((atm3 == atm1)||(atm3 == atm2)) {continue;}
            for (size_t idx4 = 0; idx4 < natoms; ++idx4) {
              atm4 = idx4+1;
              if ((atm4 == atm1)||(atm4 == atm2)||(atm4 == atm3)) {continue;}
              angle = calcDihedrals(atm1,atm2,atm3,atm4);
              if (fabs(angle) < pi/2) {
                if (fabs(angle) < tolerance) {
                  planar = true;
                  swtch = false;
                  break;
                }
                else {
                  angle = calcDihedrals(idx1+1,atm3,atm4,atm2);
                  if (fabs(angle) < tolerance) {
                    planar = true;
                    swtch = true;
                    break;
                  }
                }
              }
              else {
                if (fabs(fabs(angle) - pi) < tolerance) {
                  planar = true;
                  swtch = false;
                  break;
                }
                else {
                  angle = calcDihedrals(idx1+1,atm3,atm4,atm2);
                  if ((fabs(angle) > pi/2)&&(fabs(fabs(angle) - pi) < tolerance)) {
                    planar = true;
                    swtch = true;
                    break;
                  }
                  else if ((fabs(angle) < pi/2)&&(fabs(angle) < tolerance)) {
                    planar = true;
                    swtch = true;
                    break;
                  }
                }
              }
            }
            if (planar == true) {break;}
          }
          if (planar == true) {break;}
        }
        if (planar == true) {break;}
      }
    }
    //getting one of the planes
    if (swtch) {
      atm2 = atm3;
      atm3 = atm4;
    }
    return oVector3Atom(geometry,atm1,atm2,atm3);
  }
//--------------------------------- Symmetry Functions
  bool c2lin(double tol = 0.001) {
    //function applying a C2 rotation to a linear molecule
    double minpos = geometry(1,1);
    for (size_t idx = 1; idx < natoms; ++idx) {
      if (geometry(idx+1,1) < minpos) {minpos = geometry(idx+1,1);}
    }
    double newmin = minpos+100;
    std::vector<double> nwgeom(natoms,0.0);
    for (size_t idx = 0; idx < natoms; ++idx) {
      nwgeom[idx] = minpos;
      for (size_t idx2 = 0; idx2 < natoms; ++idx2) {
        if ((geometry(idx2+1,1) < newmin)&&(geometry(idx2+1,1) > minpos)) {
          newmin = geometry(idx2+1,1);
        }
      }
      minpos = newmin;
      newmin += 100;
    }
    bool c2perp = true;
    for (size_t idx = 0; idx < natoms/2; ++idx) {
      if (fabs(nwgeom[idx]+nwgeom[natoms-idx-1]) > tol) {
        c2perp = false;
        break;
      }
    }
    return c2perp;
  }
  std::string PGroup() {return symm;}
  std::string PGroup() const {return symm;}
  void setPGroup(const std::string & _pgroup) {symm = _pgroup;}
  size_t sigmaHorizontal(size_t prcn, std::vector<std::vector<double> > CNvec, std::vector<double> SIGMAvec, double tolzero = 1e-6) {
    //function to count the sigma_h planes in the molecule
    //sigma_h have orthogonal vector with norm 1 with the principal axes
    size_t res = 0;
    //if (prcn == 0) {throw std::string("There are no Cn axes");}
    //else if (prcn > CNvec.size()) {throw std::string("Error, the principal axis is not contemplated");}
    if (prcn == 0) {return 0;}
    else if (prcn > CNvec.size()) {return 0;}
    //declaring auxiliaries
    bool foundb;
    double ax; 
    double ay;
    double az;
    double norm;
    std::vector<double> aux(3,0.0);
    std::vector<std::vector<double> > found;
    //loop over principal symmetry axes
    for (size_t idn = 1; idn < (CNvec[prcn-1].size()/3 + 1); ++idn) {
      //loop over symmetry planes
      for (size_t idx = 1; idx < (SIGMAvec.size()/3 + 1); ++idx) {
        norm = SIGMAvec[3*idx-3]*SIGMAvec[3*idx-3] + SIGMAvec[3*idx-2]*SIGMAvec[3*idx-2] + SIGMAvec[3*idx-1]*SIGMAvec[3*idx-1];
        if (norm < 100) {
          if ((SIGMAvec[3*idx-1] == 1)||(SIGMAvec[3*idx-1] == 2)||(SIGMAvec[3*idx-3] == 1)) {continue;}
          norm = CNvec[prcn-1][3*idn-3]*SIGMAvec[3*idx-3];
          norm += CNvec[prcn-1][3*idn-2]*SIGMAvec[3*idx-2];
          norm += CNvec[prcn-1][3*idn-1]*SIGMAvec[3*idx-1];
          foundb = false;
          for (size_t lp = 0; lp < found.size(); ++lp) {
            if ((SIGMAvec[3*idx-3] == found[lp][0])&&(SIGMAvec[3*idx-2] == found[lp][1])&&(SIGMAvec[3*idx-1] == found[lp][2])) {
              foundb = true;
              break;
            }
          }
          if ((fabs(fabs(norm) - 1) < tolzero)&&(!foundb)) {
            ++res;
            aux[0] = SIGMAvec[3*idx-3];
            aux[1] = SIGMAvec[3*idx-2];
            aux[2] = SIGMAvec[3*idx-1];
            found.push_back(aux);
          }
          if (SIGMAvec[3*idx-1] == 9) {
            ++res;
            aux[0] = SIGMAvec[3*idx-3];
            aux[1] = SIGMAvec[3*idx-2];
            aux[2] = SIGMAvec[3*idx-1];
            found.push_back(aux);
          }
        }
        else if ((norm > 100)&&(norm < 1000)) {
          ax = SIGMAvec[3*idx-3] - 10;
          ay = SIGMAvec[3*idx-2] - 10;
          az = SIGMAvec[3*idx-1] - 10;
          norm = CNvec[prcn-1][3*idn-3]*ax;
          norm += CNvec[prcn-1][3*idn-2]*ay;
          norm += CNvec[prcn-1][3*idn-1]*az;
          foundb = false;
          for (size_t lp = 0; lp < found.size(); ++lp) {
            if ((ax == found[lp][0])&&(ay == found[lp][1])&&(az == found[lp][2])) {
              foundb = true;
              break;
            }
          }
          if ((fabs(norm) < tolzero)&&(!foundb)) {
            ++res;
            aux[0] = ax;
            aux[1] = ay;
            aux[2] = az;
            found.push_back(aux);
          }
        }
        else if (norm > 1000) {
          ax = SIGMAvec[3*idx-3] - 20;
          ay = SIGMAvec[3*idx-2] - 20;
          az = SIGMAvec[3*idx-1] - 20;
          norm = CNvec[prcn-1][3*idn-3]*ax;
          norm += CNvec[prcn-1][3*idn-2]*ay;
          norm += CNvec[prcn-1][3*idn-1]*az;
          foundb = false;
          for (size_t lp = 0; lp < found.size(); ++lp) {
            if ((ax == found[lp][0])&&(ay == found[lp][1])&&(az == found[lp][2])) {
              foundb = true;
              break;
            }
          }
          if ((fabs(norm) < tolzero)&&(!foundb)) {
            ++res;
            aux[0] = ax;
            aux[1] = ay;
            aux[2] = az;
            found.push_back(aux);
          }
        }
      }
    }
    return res;
  }
  size_t sigmaVertical(size_t prcn, std::vector<std::vector<double> > CNvec, std::vector<double> SIGMAvec, double tolzero = 5e-2) {
    //function to count the sigma_v planes in the molecule
    //sigma_v have orthogonal vector with norm 0 with the principal axes
    size_t res = 0;
    //if (prcn == 0) {throw std::string("There are no Cn axes");}
    //else if (prcn > CNvec.size()) {throw std::string("Error, the principal axis is not contemplated");}
    if (prcn == 0) {return 0;}
    else if (prcn > CNvec.size()) {return 0;}
    //declaring auxiliaries
    bool foundb;
    double ax; 
    double ay;
    double az;
    double norm;
    std::vector<double> aux(3,0.0);
    std::vector<std::vector<double> > found;
    for (size_t idn = 1; idn < (CNvec[prcn-1].size()/3 + 1); ++idn) {
      //loop over symmetry planes
      for (size_t idx = 1; idx < (SIGMAvec.size()/3 + 1); ++idx) {
        norm = SIGMAvec[3*idx-3]*SIGMAvec[3*idx-3] + SIGMAvec[3*idx-2]*SIGMAvec[3*idx-2] + SIGMAvec[3*idx-1]*SIGMAvec[3*idx-1];
        if (norm < 100) {
          norm = CNvec[prcn-1][3*idn-3]*SIGMAvec[3*idx-3];
          norm += CNvec[prcn-1][3*idn-2]*SIGMAvec[3*idx-2];
          norm += CNvec[prcn-1][3*idn-1]*SIGMAvec[3*idx-1];
          foundb = false;
          for (size_t lp = 0; lp < found.size(); ++lp) {
            if ((SIGMAvec[3*idx-3] == found[lp][0])&&(SIGMAvec[3*idx-2] == found[lp][1])&&(SIGMAvec[3*idx-1] == found[lp][2])) {
              foundb = true;
              break;
            }
          }
          if (fabs(fabs(norm) < tolzero)&&(!foundb)) {
            ++res;
            aux[0] = SIGMAvec[3*idx-3];
            aux[1] = SIGMAvec[3*idx-2];
            aux[2] = SIGMAvec[3*idx-1];
            found.push_back(aux);
          }
          if ((SIGMAvec[3*idx-3] == 10)||(SIGMAvec[3*idx-2] == 9)) {
            ++res;
            aux[0] = SIGMAvec[3*idx-3];
            aux[1] = SIGMAvec[3*idx-2];
            aux[2] = SIGMAvec[3*idx-1];
            found.push_back(aux);
          }
          if ((SIGMAvec[3*idx-1] == 1)||(SIGMAvec[3*idx-2] == 1)||(SIGMAvec[3*idx-3] == 1)) {
            ++res;
            aux[0] = SIGMAvec[3*idx-3];
            aux[1] = SIGMAvec[3*idx-2];
            aux[2] = SIGMAvec[3*idx-1];
            found.push_back(aux);
          }
        }
        else if ((norm > 100)&&(norm < 1000)) {
          ax = SIGMAvec[3*idx-3] - 10;
          ay = SIGMAvec[3*idx-2] - 10;
          az = SIGMAvec[3*idx-1] - 10;
          norm = CNvec[prcn-1][3*idn-3]*ax;
          norm += CNvec[prcn-1][3*idn-2]*ay;
          norm += CNvec[prcn-1][3*idn-1]*az;
          foundb = false;
          for (size_t lp = 0; lp < found.size(); ++lp) {
            if ((ax == found[lp][0])&&(ay == found[lp][1])&&(az == found[lp][2])) {
              foundb = true;
              break;
            }
          }
          if ((fabs(fabs(norm) - 1) < tolzero)&&(!foundb)) {
            ++res;
            aux[0] = ax;
            aux[1] = ay;
            aux[2] = az;
            found.push_back(aux);
          }
        }
        else if (norm > 1000) {
          ax = SIGMAvec[3*idx-3] - 20;
          ay = SIGMAvec[3*idx-2] - 20;
          az = SIGMAvec[3*idx-1] - 20;
          norm = CNvec[prcn-1][3*idn-3]*ax;
          norm += CNvec[prcn-1][3*idn-2]*ay;
          norm += CNvec[prcn-1][3*idn-1]*az;
          foundb = false;
          for (size_t lp = 0; lp < found.size(); ++lp) {
            if ((ax == found[lp][0])&&(ay == found[lp][1])&&(az == found[lp][2])) {
              foundb = true;
              break;
            }
          }
          if ((fabs(fabs(norm) - 1) < tolzero)&&(!foundb)) {
            ++res;
            aux[0] = ax;
            aux[1] = ay;
            aux[2] = az;
            found.push_back(aux);
          }
        }
      }
    }
    return res;
  } 
  void PointGroup(double tol = 0.025, double tol2 = 0.01) {
    //function to determine a molecule's point group symmetry
    //tol is to control the tolerance in the inertia eigenvalues
    //tol2 is used to control tolerance in comparing geometries
    if (Linear()) {
      if (c2lin()) {symm = "DINFh";}
      else {symm = "CINFv";}
    }
    else {
      //if not linear, then change to the inertia's referencial
      matrixE InMom = Inertia(true,true);
      std::vector<double> InVal(3,0.0);
      InVal = MatDiag(InMom);
      //renormalize inertia's eigenvalues so that all systems are comparable and reorder the eigenvalues of inertia
      double aux = std::min(InVal[0],InVal[1]);
      double MIN = std::min(aux,InVal[2]);
      aux = std::max(InVal[0],InVal[1]);
      double MAX = std::max(aux,InVal[2]);
      aux = MIN*MIN + MAX*MAX;
      for (size_t idx = 0; idx < 3; ++idx) {
        if ((InVal[idx] != MIN)&&(InVal[idx] != MAX)) {
          aux += InVal[idx]*InVal[idx];
          InVal[1] = InVal[idx]/sqrt(aux);
          InVal[0] = MAX/sqrt(aux);
          InVal[2] = MIN/sqrt(aux);
        }
      }
      matrixE newgeom = geometry;
      //change to the inertia basis
      if (!InMom.isdiagonal(0)) {newgeom *= InMom;}
      //recenter with center of mass; this is important so that center of mass is the origin
      CorrectedGeom();
      //allocation of variables to detect symmetries
      //center of inversion
      bool invpoint = false;
      matrixE invert = Identity(3)*(-1.0);
      //just to have it clean and without "-0" in the print due to numerical approximations
      invert *= Identity(3);
      //control of cn axes
      int maxcn = 6;
      std::vector<int> cnaxes(maxcn,0);
      //vectors of cn axes; these vectors contain all axes, one after the other
      std::vector<std::vector<double> > CNaxes;
      std::vector<double> CNadd(3,0.0);
      CNaxes.push_back(CNadd);
      CNadd.clear();
      //rotating the molecule by 2*pi yields always the same; there are infinite possible C1 rotations
      cnaxes[0] = 1000;
      matrixE cnx(1,1);        //Cnx
      matrixE cny(1,1);        //Cny
      matrixE cnz(1,1);        //Cnz
      //control of sn axes
      int maxsn = 6;
      std::vector<int> snaxes(maxsn,0);
      //vectors of sn axes; these vectors contain all axes, one after the other
      std::vector<std::vector<double> > SNaxes;
      std::vector<double> SNadd(3,0.0);
      SNaxes.push_back(SNadd);
      SNadd.clear();
      //S1 projections are sigma planes
      snaxes[0] = 0;
      matrixE snx(1,1);        //Snx
      matrixE sny(1,1);        //Sny
      matrixE snz(1,1);        //Snz
      //control of symmetry planes
      int sigmaplanes = 0;
      //vectors of sigma planes
      std::vector<double> SIGMAplanes;
      //rotations for planes of symmetry
      matrixE sigmax = Identity(3);        //plane reflection on x (yz)
      sigmax(1,1) = -1.0;
      matrixE sigmay = Identity(3);        //plane reflection on y (xz)
      sigmay(2,2) = -1.0;
      matrixE sigmaz = Identity(3);        //plane reflection on z (xy)
      sigmaz(3,3) = -1.0;
      //tolerances
      double tolzero = tol2/1000;
      double tolz;
      double auxtol;
      double tolcn;
      double tolsn;
      double tolsigma;
      double tolinv;
      //auxiliaries
      matrixE rotgeom(1,1);
      matrixE testgeom(1,1);
      bool cnplane = false;
      bool snplane = false;
      bool sigmaplane = false;
      bool calc = false;
      bool px = false;
      bool py = false;
      bool pz = false;
      double axis_x = 0.0;
      double axis_y = 0.0;
      double axis_z = 0.0;
      double theta = 0.0;
      double phi = 0.0;
      double norm = 0.0;
      bool lookextra = false;
      //vector orthogonal to the molecular plane if molecule is planar
      std::vector<double> orthvect(3,-10.0);
      if (Planar()) {orthvect = oVectorPlane();}
      //Checking whether the vector orthogonal to the molecular plane is one of the referencial's axes.
      if (fabs(orthvect[0]*orthvect[0] - 1) < tolzero) {px = true;}
      else if (fabs(orthvect[1]*orthvect[1] - 1) < tolzero) {py = true;}
      else if (fabs(orthvect[2]*orthvect[2] - 1) < tolzero) {pz = true;}
      //Look for Planes of Symmetry and Centers of Inversion:
      tolsigma = 25*tol2;
      tolinv = 5*tol2;
      sigmaplane = false;
      //check whether default geometry has an inversion point
      if (!invpoint) {
        testgeom = newgeom*invert;
        if (CompGeom(newgeom,testgeom,tolinv)) {invpoint = true;}
      }
      //check whether default geometry has any plane of symmetry
      //test sigmax: sigmax is the plane keeping y and z unchanged
      testgeom = newgeom*sigmax;
      if (Prints[0] >= 4) {std::cout << "------------*** Sigma X" << std::endl;}
      if (CompGeom(newgeom,testgeom,tolsigma)) {
        sigmaplanes += 1;
        SIGMAplanes.push_back(1.0);
        SIGMAplanes.push_back(0.0);
        SIGMAplanes.push_back(0.0);
        if (Prints[0] >= 4) {std::cout << "Sigma X found" << std::endl;}
        if (px) {sigmaplane = true;}
      }
      else {if (Prints[0] >= 4) {std::cout << "Sigma X NOT found" << std::endl;}}
      //test sigmay: sigmay is the plane keeping x and z unchanged
      testgeom = newgeom*sigmay;
      if (Prints[0] >= 4) {std::cout << "------------*** Sigma Y" << std::endl;}
      if (CompGeom(newgeom,testgeom,tolsigma)) {
        sigmaplanes += 1;
        SIGMAplanes.push_back(0.0);
        SIGMAplanes.push_back(1.0);
        SIGMAplanes.push_back(0.0);
        if (Prints[0] >= 4) {std::cout << "Sigma Y found" << std::endl;}
        if (py) {sigmaplane = true;}
      }
      else {if (Prints[0] >= 4) {std::cout << "Sigma Y NOT found" << std::endl;}}
      //test sigmaz: sigmaz is the plane keeping x and y unchanged
      testgeom = newgeom*sigmaz;
      if (Prints[0] >= 4) {std::cout << "------------*** Sigma Z" << std::endl;}
      if (CompGeom(newgeom,testgeom,tolsigma)) {
        sigmaplanes += 1;
        SIGMAplanes.push_back(0.0);
        SIGMAplanes.push_back(0.0);
        SIGMAplanes.push_back(1.0);
        if (Prints[0] >= 4) {std::cout << "Sigma Z found" << std::endl;}
        if (pz) {sigmaplane = true;}
      }
      else {if (Prints[0] >= 4) {std::cout << "Sigma Z NOT found" << std::endl;}}
      if (Prints[0] >= 5) {
        std::cout << "plane orthogonal to molecule inspected:" << std::endl;
        std::cout << px << " " << py << " " << pz << std::endl;
        std::cout << "plane symmetries found so far: " << sigmaplane << std::endl;
      }
      //if molecule is planar, test also vector orthogonal to plane; lets make it the z axis
      if (Planar()) {
        //molecule is planar, rebuild geometry so that z component is zero for all atoms
        testgeom = makePlanar();
        //if by any chance testgeom == geometry, then skip the comparison, it was already made
        if (!(CompGeom(newgeom,testgeom,tol2))) {
          if (Prints[0] >= 5) {
            std::cout << "Geometry with zero z components " << std::endl;
            testgeom.Print();
            std::cout << std::endl;
          }
          //test new sigmaz
          if (Prints[0] >= 4) {std::cout << "------------*** Sigma Z'" << std::endl;}
          //if molecule is planar then sigma z' is guaranteed to exist.
          if (!sigmaplane) {
            sigmaplanes += 1;
            SIGMAplanes.push_back(0.0);
            SIGMAplanes.push_back(0.0);
            SIGMAplanes.push_back(9.0);
            sigmaplane = true;
            if (Prints[0] >= 4) {std::cout << "Sigma Z' found" << std::endl;}
          }
          else {if (Prints[0] >= 4) {std::cout << "Sigma Z' NOT found" << std::endl;}}
          //look for new sigma y
          if (Prints[0] >= 4) {std::cout << "------------*** Sigma Y'" << std::endl;}
          if (py) {
            if (!sigmaplane) {calc = true;}
            else {calc = false;}
          }
          else {calc = true;}
          if (calc) {
            rotgeom = testgeom*sigmay;
            if (CompGeom(testgeom,rotgeom,tolsigma)) {
              sigmaplanes += 1;
              SIGMAplanes.push_back(0.0);
              SIGMAplanes.push_back(9.0);
              SIGMAplanes.push_back(0.0);
              if (Prints[0] >= 4) {std::cout << "Sigma Y' found" << std::endl;}
            }
            else {if (Prints[0] >= 4) {std::cout << "Sigma Y' NOT found" << std::endl;}}
          }
          else {if (Prints[0] >= 4) {std::cout << "Plane previously found (Y')" << std::endl;}}
          //look for new sigma x
          if (Prints[0] >= 4) {std::cout << "------------*** Sigma X'" << std::endl;}
          if (py) {
            if (!sigmaplane) {calc = true;}
            else {calc = false;}
          }
          else {calc = true;}
          if (calc) {
            rotgeom = testgeom*sigmax;
            if (CompGeom(testgeom,rotgeom,tolsigma)) {
              sigmaplanes += 1;
              SIGMAplanes.push_back(9.0);
              SIGMAplanes.push_back(0.0);
              SIGMAplanes.push_back(0.0);
              if (Prints[0] >= 4) {std::cout << "Sigma X' found" << std::endl;}
            }
            else {if (Prints[0] >= 4) {std::cout << "Sigma X' NOT found" << std::endl;}}
          }
          else {if (Prints[0] >= 4) {std::cout << "Plane previously found (X')" << std::endl;}}
          //check again whether there is inversion point, if not previously found
          if (!invpoint) {
            rotgeom = testgeom*invert;
            if (CompGeom(testgeom,rotgeom,tolinv)) {invpoint = true;}
          }
        }
      }
      if (Prints[0] >= 5) {
        std::cout << "plane orthogonal to molecule inspected (2):" << std::endl;
        std::cout << px << " " << py << " " << pz << std::endl;
        std::cout << "plane symmetries found so far: " << sigmaplane << std::endl;
      }
      bool existSIGMA = false;
      for (size_t atm = 0; atm < natoms; ++atm) {
        //look for axes along the vector CM-atom and CM-(middle 2 atoms)
        //these are called here mix vectors
        for (size_t atm2 = atm; atm2 < natoms; ++atm2) {
          existSIGMA = false;
          //calculate vector of reorientation
          aux = newgeom(atm+1,1) + newgeom(atm2+1,1);
          norm = aux*aux;
          aux = newgeom(atm+1,2) + newgeom(atm2+1,2);
          norm += aux*aux;
          aux = newgeom(atm+1,3) + newgeom(atm2+1,3);
          norm += aux*aux;
          axis_x = newgeom(atm+1,1) + newgeom(atm2+1,1);
          axis_y = newgeom(atm+1,2) + newgeom(atm2+1,2);
          axis_z = newgeom(atm+1,3) + newgeom(atm2+1,3);
          //if the atom or their mixed vector is the center of mass then just skip this iteration
          if (norm < tolzero*tolzero) {continue;}
          axis_x /= sqrt(norm);
          axis_y /= sqrt(norm);
          axis_z /= sqrt(norm);
          //verify whether this rotation vector coincides with a natural axis of the molecule
          norm = orthvect[0]*axis_x + orthvect[1]*axis_y + orthvect[2]*axis_z;
          if (Prints[0] >= 5) {
            std::cout << "Internal product with molecule's normal vector: " << norm << std::endl;
            std::cout << "!!!Axes: " << axis_x << " " << axis_y << " " << axis_z << std::endl;
          }
          //checking whether this sigma plane was previously detected for planar molecules
          tolz = 5000*tolzero;
          if (!Planar()) {tolz = 100*tolzero;}
          //comparing new mix vector with old ones
          for (size_t idn = 1; idn < (SIGMAplanes.size()/3 + 1); ++idn) {
            if (Prints[0] >= 5) {
              std::cout << "For sigma planar";
              std::cout << axis_x << "(" << SIGMAplanes[3*idn-3] << ") " << axis_y << "(" << SIGMAplanes[3*idn-2] << ") ";
              std::cout << axis_z << "(" << SIGMAplanes[3*idn-1] << ") " << std::endl;
            }
            //check whether some similar vector was already found
            aux = fabs(axis_x - SIGMAplanes[3*idn-3]);
            MIN = fabs(axis_y - SIGMAplanes[3*idn-2]);
            MAX = fabs(axis_z - SIGMAplanes[3*idn-1]);
            if ((aux < tolz)&&(MIN < tolz)&&(MAX < tolz)) {existSIGMA = true;}
            aux = fabs(axis_x + SIGMAplanes[3*idn-3]);
            MIN = fabs(axis_y + SIGMAplanes[3*idn-2]);
            MAX = fabs(axis_z + SIGMAplanes[3*idn-1]);
            if ((aux < tolz)&&(MIN < tolz)&&(MAX < tolz)&&(!existSIGMA)) {existSIGMA = true;}
            aux = fabs(axis_x + 10 - SIGMAplanes[3*idn-3]);
            MIN = fabs(axis_y + 10 - SIGMAplanes[3*idn-2]);
            MAX = fabs(axis_z + 10 - SIGMAplanes[3*idn-1]);
            if ((aux < tolz)&&(MIN < tolz)&&(MAX < tolz)&&(!existSIGMA)) {existSIGMA = true;}
            aux = fabs(axis_x - 10 + SIGMAplanes[3*idn-3]);
            MIN = fabs(axis_y - 10 + SIGMAplanes[3*idn-2]);
            MAX = fabs(axis_z - 10 + SIGMAplanes[3*idn-1]);
            if ((aux < tolz)&&(MIN < tolz)&&(MAX < tolz)&&(!existSIGMA)) {existSIGMA = true;}
            aux = fabs(axis_x + 20 - SIGMAplanes[3*idn-3]);
            MIN = fabs(axis_y + 20 - SIGMAplanes[3*idn-2]);
            MAX = fabs(axis_z + 20 - SIGMAplanes[3*idn-1]);
            if ((aux < tolz)&&(MIN < tolz)&&(MAX < tolz)&&(!existSIGMA)) {existSIGMA = true;}
            aux = fabs(axis_x - 20 + SIGMAplanes[3*idn-3]);
            MIN = fabs(axis_y - 20 + SIGMAplanes[3*idn-2]);
            MAX = fabs(axis_z - 20 + SIGMAplanes[3*idn-1]);
            if ((aux < tolz)&&(MIN < tolz)&&(MAX < tolz)&&(!existSIGMA)) {existSIGMA = true;}
            //comparing norms now
            norm = axis_x*SIGMAplanes[3*idn-3] + axis_y*SIGMAplanes[3*idn-2] + axis_z*SIGMAplanes[3*idn-1];
            if (fabs(fabs(norm) - 1) < 2*tolz) {existSIGMA = true;}
            if (existSIGMA) {break;}
          }
          //comparing new mix vector with the plane defined by atm, atm2 and the reference (center of mass)
          //this skipping only valid for planar molecules, of course!
          if (!existSIGMA) {
            std::vector<double> atomplane = oVector2Atom(newgeom,atm+1,atm2+1);
            if (Prints[0] >= 5) {
              std::cout << "2-atom plane ";
              std::cout << "(" << atomplane[0] << " " << atomplane[1] << " " << atomplane[2] << ")" << std::endl;
            }
            for (size_t idn = 1; idn < (SIGMAplanes.size()/3 + 1); ++idn) {
              norm = atomplane[0]*SIGMAplanes[3*idn-3] + atomplane[1]*SIGMAplanes[3*idn-2] + atomplane[2]*SIGMAplanes[3*idn-1];
              if ((fabs(norm) < tolzero)||(fabs(fabs(norm) - 1) < tolzero)) {
                existSIGMA = true;
                break;
              }
            }
          }
          //there is nothing else here because non-planar molecules are calculated differently below
          //actual determination of the mix vector and symmetry planes
          //for non-planar molecules existSIGMA == false since we only checked and eventually changed for planar molecules
          if (!existSIGMA) {
            //making the mix vector rotation so that it becomes (0,0,1)
            //getting the rotation angles theta (axis_x -> 0) and phi (axis_y -> 0)
            theta = atan(axis_x/axis_y);
            aux = axis_x*sin(theta) + axis_y*cos(theta);
            phi = atan(aux/axis_z);
            //rotation theta
            testgeom = newgeom*Givens(1,2,3,cos(theta),sin(theta));
            //rotation phi
            rotgeom = testgeom*Givens(2,3,3,cos(phi),sin(phi));
            //recenter ??????????????????????
            //CMg(rotgeom);
            if (Prints[0] >= 4) {std::cout << "------------ mix: " << axis_x << " " << axis_y << " " << axis_z << std::endl;}
            size_t counter = 0;
            //check whether vector (axis_x, axis_y, axis_z) is on the molecular plane
            norm = orthvect[0]*axis_x + orthvect[1]*axis_y + orthvect[2]*axis_z;
            //test mix Sigma X
            if (Prints[0] >= 4) {std::cout << "------------*** Sigma X mix" << std::endl;}
            calc = true;
            if ((fabs(norm) < 50*tolzero)||(fabs(fabs(norm) - 1) < 50*tolzero)) {
              //calculate only if sigmaplane not yet defined, which here should be always false
              if (!sigmaplane) {calc = true;}
              else {
                calc = false;
                ++counter;
              }
            }
            if (calc) {
              //apply symmetry plane
              testgeom = rotgeom*sigmax;
              if (CompGeom(rotgeom,testgeom,tolsigma)) {
                sigmaplanes += 1;
                SIGMAplanes.push_back(axis_x);
                SIGMAplanes.push_back(axis_y);
                SIGMAplanes.push_back(axis_z);
                if ((!sigmaplane)&&((fabs(norm) < 50*tolzero)||(fabs(fabs(norm) - 1) < 50*tolzero))) {sigmaplane = true;}
                if (!Planar()) {sigmaplane = false;}
                if (Prints[0] >= 4) {std::cout << "Sigma X mix found" << std::endl;}
              }
              else {if (Prints[0] >= 4) {std::cout << "Sigma X mix NOT found" << std::endl;}}
            }
            else {if (Prints[0] >= 4) {std::cout << "Plane previously found (X mix)" << std::endl;}}
            //test mix Sigma Y
            if (Prints[0] >= 4) {std::cout << "------------*** Sigma Y mix" << std::endl;}
            calc = true;
            if ((fabs(norm) < 50*tolzero)||(fabs(fabs(norm) - 1) < 50*tolzero)) {
              //calculate only if sigmaplane not yet defined, which here should be always false
              if (!sigmaplane) {calc = true;}
              else if ((sigmaplane)&&(counter <= 1)) {calc = true;}
              else {calc = false;}
            }
            if (calc) {
              testgeom = rotgeom*sigmay;
              if (CompGeom(rotgeom,testgeom,tolsigma)) {
                sigmaplanes += 1;
                SIGMAplanes.push_back(axis_x+10);
                SIGMAplanes.push_back(axis_y+10);
                SIGMAplanes.push_back(axis_z+10);
                if (Prints[0] >= 4) {std::cout << "Sigma Y mix found" << std::endl;}
              }
              else {if (Prints[0] >= 4) {std::cout << "Sigma Y mix NOT found" << std::endl;}}
            }
            else {if (Prints[0] >= 4) {std::cout << "Plane previously found (Y mix)" << std::endl;}}
            //test mix Sigma Z
            if (Prints[0] >= 4) {std::cout << "------------*** Sigma Z mix" << std::endl;}
            //check whether vector (axis_x, axis_y, axis_z) is orthogonal to molecular plane
            calc = true;
            if ((fabs(norm) < 50*tolzero)||(fabs(fabs(norm) - 1) < 50*tolzero)) {
              //calculate only if sigmaplane not yet defined, which here should be always false
              if (!sigmaplane) {calc = true;}
              else {calc = false;}
            }
            if (calc) {
              testgeom = rotgeom*sigmaz;
              if (CompGeom(rotgeom,testgeom,tolsigma)) {
                sigmaplanes += 1;
                SIGMAplanes.push_back(axis_x+20);
                SIGMAplanes.push_back(axis_y+20);
                SIGMAplanes.push_back(axis_z+20);
                if (Prints[0] >= 4) {std::cout << "Sigma Z mix found" << std::endl;}
              }
              if (Prints[0] >= 4) {std::cout << "Sigma Z mix NOT found" << std::endl;}
            }
            else {if (Prints[0] >= 4) {std::cout << "Plane previously found (Z mix)" << std::endl;}}
            if (!Planar()) {
              for (size_t atm3 = 0; atm3 < natoms; ++atm3) {
                if ((atm3 != atm)&&(atm3 != atm2)&&(atm != atm2)) {
                  rotgeom = newgeom;
                  aux = newgeom(atm3+1,1);
                  MIN = newgeom(atm3+1,2);
                  MAX = newgeom(atm3+1,3);
                  for (size_t idp1 = 0; idp1 < natoms; ++idp1) {
                    rotgeom(idp1+1,1) -= aux;
                    rotgeom(idp1+1,2) -= MIN;
                    rotgeom(idp1+1,3) -= MAX;
                  }
                  //calculate the new middle point
                  axis_x = rotgeom(atm+1,1) + rotgeom(atm2+1,1);
                  axis_y = rotgeom(atm+1,2) + rotgeom(atm2+1,2);
                  axis_z = rotgeom(atm+1,3) + rotgeom(atm2+1,3);
                  //renormalize new middle point
                  norm = axis_x*axis_x + axis_y*axis_y + axis_z*axis_z;
                  axis_x /= sqrt(norm);
                  axis_y /= sqrt(norm);
                  axis_z /= sqrt(norm);
                  if (Prints[0] >= 4) {std::cout << "------------ mix (non-planar): " << axis_x << " " << axis_y << " " << axis_z << std::endl;}
                  //check whether the new mix is an axis of the referencial
                  if ((fabs(axis_x*axis_x - 1) < tolzero)&&(axis_y*axis_y < tolzero)&&(axis_z*axis_z < tolzero)) {continue;}
                  else if ((fabs(axis_y*axis_y - 1) < tolzero)&&(axis_x*axis_x < tolzero)&&(axis_z*axis_z < tolzero)) {continue;}
                  else if ((fabs(axis_z*axis_z - 1) < tolzero)&&(axis_y*axis_y < tolzero)&&(axis_x*axis_x < tolzero)) {continue;}
                  existSIGMA = false;
                  //comparing new mix vector with old ones
                  for (size_t idn = 1; idn < (SIGMAplanes.size()/3 + 1); ++idn) {
                    if (Prints[0] >= 5) {
                      std::cout << "For sigma non-planar";
                      std::cout << axis_x << "(" << SIGMAplanes[3*idn-3] << ") " << axis_y << "(" << SIGMAplanes[3*idn-2] << ") ";
                      std::cout << axis_z << "(" << SIGMAplanes[3*idn-1] << ") " << std::endl;
                    }
                    tolz = 1000*tolzero;
                    //std::cout << "tolerance of " << tolz << std::endl;
                    aux = fabs(axis_x - SIGMAplanes[3*idn-3]);
                    MIN = fabs(axis_y - SIGMAplanes[3*idn-2]);
                    MAX = fabs(axis_z - SIGMAplanes[3*idn-1]);
                    if ((aux < tolz)&&(MIN < tolz)&&(MAX < tolz)) {existSIGMA = true;}
                    aux = fabs(axis_x + SIGMAplanes[3*idn-3]);
                    MIN = fabs(axis_y + SIGMAplanes[3*idn-2]);
                    MAX = fabs(axis_z + SIGMAplanes[3*idn-1]);
                    if ((aux < tolz)&&(MIN < tolz)&&(MAX < tolz)&&(!existSIGMA)) {existSIGMA = true;}
                    aux = fabs(axis_x + 10 - SIGMAplanes[3*idn-3]);
                    MIN = fabs(axis_y + 10 - SIGMAplanes[3*idn-2]);
                    MAX = fabs(axis_z + 10 - SIGMAplanes[3*idn-1]);
                    if ((aux < tolz)&&(MIN < tolz)&&(MAX < tolz)&&(!existSIGMA)) {existSIGMA = true;}
                    aux = fabs(axis_x - 10 + SIGMAplanes[3*idn-3]);
                    MIN = fabs(axis_y - 10 + SIGMAplanes[3*idn-2]);
                    MAX = fabs(axis_z - 10 + SIGMAplanes[3*idn-1]);
                    if ((aux < tolz)&&(MIN < tolz)&&(MAX < tolz)&&(!existSIGMA)) {existSIGMA = true;}
                    aux = fabs(axis_x + 20 - SIGMAplanes[3*idn-3]);
                    MIN = fabs(axis_y + 20 - SIGMAplanes[3*idn-2]);
                    MAX = fabs(axis_z + 20 - SIGMAplanes[3*idn-1]);
                    if ((aux < tolz)&&(MIN < tolz)&&(MAX < tolz)&&(!existSIGMA)) {existSIGMA = true;}
                    aux = fabs(axis_x - 20 + SIGMAplanes[3*idn-3]);
                    MIN = fabs(axis_y - 20 + SIGMAplanes[3*idn-2]);
                    MAX = fabs(axis_z - 20 + SIGMAplanes[3*idn-1]);
                    if ((aux < tolz)&&(MIN < tolz)&&(MAX < tolz)&&(!existSIGMA)) {existSIGMA = true;}
                    norm = axis_x*SIGMAplanes[3*idn-3] + axis_y*SIGMAplanes[3*idn-2] + axis_z*SIGMAplanes[3*idn-1];
                    if (fabs(fabs(norm) - 1) < 2*tolz) {existSIGMA = true;}
                    if (existSIGMA) {break;}
                  }
                  if (!existSIGMA) {
                    //making the mix vector rotation so that it becomes (0,0,1)
                    //getting the rotation angles theta (axis_x -> 0) and phi (axis_y -> 0)
                    theta = atan(axis_x/axis_y);
                    aux = axis_x*sin(theta) + axis_y*cos(theta);
                    phi = atan(aux/axis_z);
                    //rotation theta
                    testgeom = rotgeom*Givens(1,2,3,cos(theta),sin(theta));
                    //rotation phi
                    rotgeom = testgeom*Givens(2,3,3,cos(phi),sin(phi));
                    //recenter ??????????????????????
                    //CMg(rotgeom);
                    //test mix Sigma X
                    if (Prints[0] >= 4) {std::cout << "------------*** Sigma X mix (non-planar)" << std::endl;}
                    testgeom = rotgeom*sigmax;
                    if (CompGeom(rotgeom,testgeom,tolsigma)) {
                      sigmaplanes += 1;
                      SIGMAplanes.push_back(axis_x);
                      SIGMAplanes.push_back(axis_y);
                      SIGMAplanes.push_back(axis_z);
                      if (Prints[0] >= 4) {std::cout << "Sigma X mix found" << std::endl;}
                    }
                    else {if (Prints[0] >= 4) {std::cout << "Sigma X mix NOT found" << std::endl;}}
                    //test mix Sigma Y
                    if (Prints[0] >= 4) {std::cout << "------------*** Sigma Y mix (non-planar)" << std::endl;}
                    testgeom = rotgeom*sigmay;
                    if (CompGeom(rotgeom,testgeom,tolsigma)) {
                      sigmaplanes += 1;
                      SIGMAplanes.push_back(axis_x+10);
                      SIGMAplanes.push_back(axis_y+10);
                      SIGMAplanes.push_back(axis_z+10);
                      if (Prints[0] >= 4) {std::cout << "Sigma Y mix found" << std::endl;}
                    }
                    else {if (Prints[0] >= 4) {std::cout << "Sigma Y mix NOT found" << std::endl;}}
                    //test mix Sigma Z
                    if (Prints[0] >= 4) {std::cout << "------------*** Sigma Z mix (non-planar)" << std::endl;}
                    testgeom = rotgeom*sigmaz;
                    if (CompGeom(rotgeom,testgeom,tolsigma)) {
                      sigmaplanes += 1;
                      SIGMAplanes.push_back(axis_x+20);
                      SIGMAplanes.push_back(axis_y+20);
                      SIGMAplanes.push_back(axis_z+20);
                      if (Prints[0] >= 4) {std::cout << "Sigma Z mix found" << std::endl;} 
                    }
                    else {if (Prints[0] >= 4) {std::cout << "Sigma Z mix NOT found" << std::endl;}}
                  }
                }
              }
            }
          }
        }
      }
      //look for Cn and Sn
      for (int idx = 1; idx < maxcn; ++idx) {
        lookextra = false;
        //generate the cn transformation matrices for main axes
        if (idx == 1) {
          tolcn = 2*tol2;
          tolsn = 1.5*tol2;
        }
        else if (idx == 2) {
          tolcn = 8*tol2;
          tolsn = 8*tol2;
        }
        else if (idx == 3) {
          tolcn = 8*tol2;
          tolsn = 8*tol2;
        }
        else if (idx == 4) {
          tolcn = 8*tol2;
          tolsn = 8*tol2;
        }
        else if (idx == 5) {
          tolcn = 8*tol2;
          tolsn = 8*tol2;
        }
        else {
          tolcn = 5*tol2;
          tolsn = 5*tol2;
        }
        cnplane = false;
        snplane = false;
        cnx = cn(1,idx+1);
        cny = cn(2,idx+1);
        cnz = cn(3,idx+1);
        snx = cnx*sigmax;
        sny = cny*sigmay;
        snz = cnz*sigmaz;
        size_t basicSN = 0;
        //test Cn X
        if (Prints[0] >= 4) {std::cout << "------------*** C" << idx+1 << " X " << std::endl;}
        testgeom = newgeom*cnx;
        if (CompGeom(newgeom,testgeom,tolcn)) {
          cnaxes[idx] += 1;
          CNadd.push_back(1.0);
          CNadd.push_back(0.0);
          CNadd.push_back(0.0);
          if (idx > 1) {lookextra = true;}
          if (Prints[0] >= 4) {std::cout << "C" << idx+1 << " X found" << std::endl;}
          if (px) {cnplane = true;}
        }
        else {if (Prints[0] >= 4) {std::cout << "C" << idx+1 << " X NOT found" << std::endl;}}
        //test Cn Y
        testgeom = newgeom*cny;
        if (Prints[0] >= 4) {std::cout << "------------*** C" << idx+1 << " Y " << std::endl;}
        if (CompGeom(newgeom,testgeom,tolcn)) {
          cnaxes[idx] += 1;
          CNadd.push_back(0.0);
          CNadd.push_back(1.0);
          CNadd.push_back(0.0);
          if (idx > 1) {lookextra = true;}
          if (Prints[0] >= 4) {std::cout << "C" << idx+1 << " Y found" << std::endl;}
          if (py) {cnplane = true;}
        }
        else {if (Prints[0] >= 4) {std::cout << "C" << idx+1 << " Y NOT found" << std::endl;}}
        //test Cn Z
        testgeom = newgeom*cnz;
        if (Prints[0] >= 4) {std::cout << "------------*** C" << idx+1 << " Z " << std::endl;}
        if (CompGeom(newgeom,testgeom,tolcn)) {
          cnaxes[idx] += 1;
          CNadd.push_back(0.0);
          CNadd.push_back(0.0);
          CNadd.push_back(1.0);
          if (idx > 1) {lookextra = true;}
          if (Prints[0] >= 4) {std::cout << "C" << idx+1 << " Z found" << std::endl;}
          if (pz) {cnplane = true;}
        }
        else {if (Prints[0] >= 4) {std::cout << "C" << idx+1 << " Z NOT found" << std::endl;}}
        //test snx
        testgeom = newgeom*snx;
        if (Prints[0] >= 4) {std::cout << "------------*** S" << idx+1 << " X " << std::endl;}
        if (CompGeom(newgeom,testgeom,tolsn)) {
          snaxes[idx] += 1;
          SNadd.push_back(1.0);
          SNadd.push_back(0.0);
          SNadd.push_back(0.0);
          ++basicSN;
          if (Prints[0] >= 4) {std::cout << "S" << idx+1 << " X found" << std::endl;}
          if (px) {snplane = true;}
        }
        else {if (Prints[0] >= 4) {std::cout << "S" << idx+1 << " X NOT found" << std::endl;}}
        //test Sn Y
        testgeom = newgeom*sny;
        if (Prints[0] >= 4) {std::cout << "------------*** S" << idx+1 << " Y " << std::endl;}
        if (CompGeom(newgeom,testgeom,tolsn)) {
          snaxes[idx] += 1;
          SNadd.push_back(0.0);
          SNadd.push_back(1.0);
          SNadd.push_back(0.0);
          ++basicSN;
          if (Prints[0] >= 4) {std::cout << "S" << idx+1 << " Y found" << std::endl;}
          if (py) {snplane = true;}
        }
        else {if (Prints[0] >= 4) {std::cout << "S" << idx+1 << " Y NOT found" << std::endl;}}
        //test Sn Z
        testgeom = newgeom*snz;
        if (Prints[0] >= 4) {std::cout << "------------*** S" << idx+1 << " Z " << std::endl;}
        if (CompGeom(newgeom,testgeom,tolsn)) {
          snaxes[idx] += 1;
          SNadd.push_back(0.0);
          SNadd.push_back(0.0);
          SNadd.push_back(1.0);
          ++basicSN;
          if (Prints[0] >= 4) {std::cout << "S" << idx+1 << " Z found" << std::endl;}
          if (pz) {snplane = true;}
        }
        else {if (Prints[0] >= 4) {std::cout << "S" << idx+1 << " Z NOT found" << std::endl;}}
        if (Prints[0] >= 5) {
          std::cout << "plane orthogonal to molecule inspected (3):" << std::endl;
          std::cout << px << " " << py << " " << pz << std::endl;
          std::cout << "plane symmetries found so far: " << cnplane << "; " << snplane  << std::endl;
        }
        if (Planar()) {
          //molecule is planar, rebuild geometry so that z component is zero for all atoms
          testgeom = makePlanar();
          //if by any chance testgeom == geometry, then skip the comparison, it was already made
          if (!(CompGeom(newgeom,testgeom,tol2))) {
            //test Cn X'
            if (Prints[0] >= 4) {std::cout << "------------*** C" << idx+1 << " X' " << std::endl;}
            if ((px)&&(!cnplane)) {calc = true;}
            else if (!px) {calc = true;}
            else {calc = false;}
            if (calc) {
              rotgeom = testgeom*cnx;
              if (CompGeom(testgeom,rotgeom,tolcn)) {
                cnaxes[idx] += 1;
                CNadd.push_back(1.0);
                CNadd.push_back(0.0);
                CNadd.push_back(0.0);
                if (idx > 1) {lookextra = true;}
                if (Prints[0] >= 4) {std::cout << "C" << idx+1 << " X' found" << std::endl;}
                if (px) {cnplane = true;}
              }
              else {if (Prints[0] >= 4) {std::cout << "C" << idx+1 << " X' NOT found" << std::endl;}}
            }
            else {if (Prints[0] >= 4) {std::cout << "Axis previously found (CX')" << std::endl;}}
            //test Cn Y'
            if (Prints[0] >= 4) {std::cout << "------------*** C" << idx+1 << " Y' " << std::endl;}
            if ((py)&&(!cnplane)) {calc = true;}
            else if (!py) {calc = true;}
            else {calc = false;}
            if (calc) {
              rotgeom = testgeom*cny;
              if (CompGeom(testgeom,rotgeom,tolcn)) {
                cnaxes[idx] += 1;
                CNadd.push_back(0.0);
                CNadd.push_back(1.0);
                CNadd.push_back(0.0);
                if (idx > 1) {lookextra = true;}
                if (Prints[0] >= 4) {std::cout << "C" << idx+1 << " Y' found" << std::endl;}
                if (py) {cnplane = true;}
              }
              else {if (Prints[0] >= 4) {std::cout << "C" << idx+1 << " Y' NOT found" << std::endl;}}
            }
            else {if (Prints[0] >= 4) {std::cout << "Axis previously found (CY')" << std::endl;}}
            //test Cn Z'
            if (Prints[0] >= 4) {std::cout << "------------*** C" << idx+1 << " Z' " << std::endl;}
            if ((pz)&&(!cnplane)) {calc = true;}
            else if (!pz) {calc = true;}
            else {calc = false;}
            if (calc) {
              rotgeom = testgeom*cnz;
              if (CompGeom(testgeom,rotgeom,tolcn)) {
                cnaxes[idx] += 1;
                CNadd.push_back(0.0);
                CNadd.push_back(0.0);
                CNadd.push_back(1.0);
                if (idx > 1) {lookextra = true;}
                if (Prints[0] >= 4) {std::cout << "C" << idx+1 << " Z' found" << std::endl;}
                if (pz) {cnplane = true;}
              }
              else {if (Prints[0] >= 4) {std::cout << "C" << idx+1 << " Z' NOT found" << std::endl;}}
            }
            else {if (Prints[0] >= 4) {std::cout << "Axis previously found (CZ')" << std::endl;}}
            //test Sn X'
            if (Prints[0] >= 4) {std::cout << "------------*** S" << idx+1 << " X' " << std::endl;}
            if ((px)&&(!snplane)) {calc = true;}
            else if (!px) {calc = true;}
            else {calc = false;}
            if (basicSN >= 3) {calc = false;}
            if (calc) {
              rotgeom = testgeom*snx;
              if (CompGeom(testgeom,rotgeom,tolsn)) {
                snaxes[idx] += 1;
                SNadd.push_back(1.0);
                SNadd.push_back(0.0);
                SNadd.push_back(0.0);
                if (Prints[0] >= 4) {std::cout << "S" << idx+1 << " X' found" << std::endl;}
                if (px) {snplane = true;}
              }
              else {if (Prints[0] >= 4) {std::cout << "S" << idx+1 << " X' NOT found" << std::endl;}}
            }
            else {if (Prints[0] >= 4) {std::cout << "Axis previously found (SX')" << std::endl;}}
            //test Sn Y'
            if (Prints[0] >= 4) {std::cout << "------------*** S" << idx+1 << " Y' " << std::endl;}
            if ((py)&&(!snplane)) {calc = true;}
            else if (!py) {calc = true;}
            else {calc = false;}
            if (basicSN >= 3) {calc = false;}
            if (calc) {
              rotgeom = testgeom*sny;
              if (CompGeom(testgeom,rotgeom,tolsn)) {
                snaxes[idx] += 1;
                SNadd.push_back(0.0);
                SNadd.push_back(1.0);
                SNadd.push_back(0.0);
                if (Prints[0] >= 4) {std::cout << "S" << idx+1 << " Y' found" << std::endl;}
                if (py) {snplane = true;}
              }
              else {if (Prints[0] >= 4) {std::cout << "S" << idx+1 << " Y' NOT found" << std::endl;}}
            }
            else {if (Prints[0] >= 4) {std::cout << "Axis previously found (SY')" << std::endl;}}
            //test Sn Z'
            if (Prints[0] >= 4) {std::cout << "------------*** S" << idx+1 << " Z' " << std::endl;}
            if ((pz)&&(!snplane)) {calc = true;}
            else if (!pz) {calc = true;}
            else {calc = false;}
            if (basicSN >= 3) {calc = false;}
            if (calc) {
              rotgeom = testgeom*snz;
              if (CompGeom(testgeom,rotgeom,tolsn)) {
                snaxes[idx] += 1;
                SNadd.push_back(0.0);
                SNadd.push_back(0.0);
                SNadd.push_back(1.0);
                if (Prints[0] >= 4) {std::cout << "S" << idx+1 << " Z' found" << std::endl;}
                if (pz) {snplane = true;}
              }
              else {if (Prints[0] >= 4) {std::cout << "S" << idx+1 << " Z' NOT found" << std::endl;}}
            }
            else {if (Prints[0] >= 4) {std::cout << "Axis previously found (SZ')" << std::endl;}}
          }
        }
        if (Prints[0] >= 5) {
          std::cout << "plane orthogonal to molecule inspected (4):" << std::endl;
          std::cout << px << " " << py << " " << pz << std::endl;
          std::cout << "plane symmetries found so far: " << cnplane << "; " << snplane  << std::endl;
        }
        bool existCN = false;
        bool existSN = false;
        for (size_t atm = 0; atm < natoms; ++atm) {
          //look for axes along the vector CM-atoms and the mixed ones
          for (size_t atm2 = atm; atm2 < natoms; ++atm2) {
            if (atoms[atm] != atoms[atm2]) {continue;}
            existCN = false;
            existSN = false;
            aux = newgeom(atm + 1,1) + newgeom(atm2 + 1,1);
            norm = aux*aux;
            aux = newgeom(atm + 1,2) + newgeom(atm2 + 1,2);
            norm += aux*aux;
            aux = newgeom(atm + 1,3) + newgeom(atm2 + 1,3);
            norm += aux*aux;
            axis_x = newgeom(atm + 1,1) + newgeom(atm2+1,1);
            axis_y = newgeom(atm + 1,2) + newgeom(atm2+1,2);
            axis_z = newgeom(atm + 1,3) + newgeom(atm2+1,3);
            //if the atom or their mixed vector is the center of mass then just skip this iteration
            if (norm < tolzero*tolzero) {continue;}
            axis_x /= sqrt(norm);
            axis_y /= sqrt(norm);
            axis_z /= sqrt(norm);
            if (Prints[0] >= 4) {std::cout << "------------ mix: " << axis_x << " " << axis_y << " " << axis_z << std::endl;}
            norm = orthvect[0]*axis_x + orthvect[1]*axis_y + orthvect[2]*axis_z;
            if (Prints[0] >= 5) {
              std::cout << "Internal product with molecule's normal vector: " << norm << std::endl;
              std::cout << "!!!Axes: " << axis_x << " " << axis_y << " " << axis_z << std::endl;
            }
            //checking whether this cn axis was previously identified
            double tolz = 250*tolzero;
            //std::cout << "tolerance of " << tolz << std::endl;
            for (size_t idn = 1; idn < (CNadd.size()/3 + 1); ++idn) {
              if (Prints[0] >= 5) {
                std::cout << "For Cn";
                std::cout << axis_x << "(" << CNadd[3*idn-3] << ") " << axis_y << "(" << CNadd[3*idn-2] << ") ";
                std::cout << axis_z << "(" << CNadd[3*idn-1] << ") " << std::endl;
              }
              norm = axis_x*CNadd[3*idn-3] + axis_y*CNadd[3*idn-2] + axis_z*CNadd[3*idn-1];
              aux = fabs(axis_x - CNadd[3*idn-3]);
              MIN = fabs(axis_y - CNadd[3*idn-2]);
              MAX = fabs(axis_z - CNadd[3*idn-1]);
              if ((aux < tolz)&&(MIN < tolz)&&(MAX < tolz)) {existCN = true;}
              aux = fabs(axis_x + CNadd[3*idn-3]);
              MIN = fabs(axis_y + CNadd[3*idn-2]);
              MAX = fabs(axis_z + CNadd[3*idn-1]);
              if ((aux < tolz)&&(MIN < tolz)&&(MAX < tolz)&&(!existCN)) {existCN = true;}
              if ((fabs(fabs(norm) - 1) < tolz)&&(!existCN)) {existCN = true;}
              if (existCN) {break;}
            }
            //checking whether this sn axis was previously identified
            for (size_t idn = 1; idn < (SNadd.size()/3 + 1); ++idn) {
              if (Prints[0] >= 5) {
                std::cout << "For Sn";
                std::cout << axis_x << "(" << SNadd[3*idn-3] << ") " << axis_y << "(" << SNadd[3*idn-2] << ") ";
                std::cout << axis_z << "(" << SNadd[3*idn-1] << ") " << std::endl;
              }
              norm = axis_x*SNadd[3*idn-3] + axis_y*SNadd[3*idn-2] + axis_z*SNadd[3*idn-1];
              aux = fabs(axis_x - SNadd[3*idn-3]);
              MIN = fabs(axis_y - SNadd[3*idn-2]);
              MAX = fabs(axis_z - SNadd[3*idn-1]);
              if ((aux < 4*tolz)&&(MIN < 4*tolz)&&(MAX < 4*tolz)) {existSN = true;}
              aux = fabs(axis_x + SNadd[3*idn-3]);
              MIN = fabs(axis_y + SNadd[3*idn-2]);
              MAX = fabs(axis_z + SNadd[3*idn-1]);
              if ((aux < 4*tolz)&&(MIN < 4*tolz)&&(MAX < 4*tolz)&&(idx == 1)) {existSN = true;}
              if ((fabs(fabs(norm) - 1) < 20*tolz)&&(!existSN)&&(idx == 1)) {existSN = true;}
              if ((fabs(fabs(norm) - 1) < 4*tolz)&&(!existSN)&&(idx != 1)) {existSN = true;}
              //one last attempt to reduce the S2 axes found by using atom planes of the new mix and the previously found axes
              if (idx == 1) {
                double auxtol = 100*tolzero;
                if (!Planar()) {auxtol *= 25;}
                if ((!existSN)&&(atm != atm2)) {
                  std::vector<double> atomplane = oVector2Atom(newgeom,atm+1,atm2+1);
                  norm = atomplane[0]*SNadd[3*idn-3] + atomplane[1]*SNadd[3*idn-2] + atomplane[2]*SNadd[3*idn-1];
                  if ((fabs(norm) < auxtol)||(fabs(fabs(norm) - 1) < auxtol)) {
                    existSN = true;
                    break;
                  }
                }
              }
              if (existSN) {break;}
            }
            if ((fabs(axis_x*axis_x - 1) < tolzero)&&(axis_y*axis_y < tolzero)&&(axis_z*axis_z < tolzero)) {continue;}
            else if ((fabs(axis_y*axis_y - 1) < tolzero)&&(axis_x*axis_x < tolzero)&&(axis_z*axis_z < tolzero)) {continue;}
            else if ((fabs(axis_z*axis_z - 1) < tolzero)&&(axis_y*axis_y < tolzero)&&(axis_x*axis_x < tolzero)) {continue;}
            if ((!existCN)||(!existSN)) {
              //mixing the axes
              //getting the rotation angles theta (axis_x -> 0) and phi (axis_y -> 0)
              //so that the vector (axis_x, axis_y, axis_z) is the unitary vector (0,0,1)
              theta = atan(axis_x/axis_y);
              aux = axis_x*sin(theta) + axis_y*cos(theta);
              phi = atan(aux/axis_z);
              //rotation theta
              testgeom = newgeom*Givens(1,2,3,cos(theta),sin(theta));
              //rotation phi
              rotgeom = testgeom*Givens(2,3,3,cos(phi),sin(phi));
            }
            if (!existCN) {
              if (Prints[0] >= 4) {std::cout << "------------*** C" << idx+1 << " Z mix" << std::endl;}
              //test Cn Z mix axis
              testgeom = rotgeom*cnz;
              //look for symmetry
              if (CompGeom(rotgeom,testgeom,tolcn)) {
                cnaxes[idx] += 1;
                CNadd.push_back(axis_x);
                CNadd.push_back(axis_y);
                CNadd.push_back(axis_z);
                if (idx > 1) {lookextra = true;}
                if (Prints[0] >= 4) {std::cout << "C" << idx+1 << " Z mix found" << std::endl;}
              }
              else {if (Prints[0] >= 4) {std::cout << "C" << idx+1 << " Z mix NOT found" << std::endl;}}
            }
            if (!existSN) {
              if (Prints[0] >= 4) {std::cout << "------------*** S" << idx+1 << " Z mix" << std::endl;}
              //test Sn Z mix axis
              testgeom = rotgeom*snz;
              //look for symmetry
              if (CompGeom(rotgeom,testgeom,tolsn)) {
                snaxes[idx] += 1;
                SNadd.push_back(axis_x);
                SNadd.push_back(axis_y);
                SNadd.push_back(axis_z);
                if (Prints[0] >= 4) {std::cout << "S" << idx+1 << " Z mix found" << std::endl;}
              }
              else {if (Prints[0] >= 4) {std::cout << "S" << idx+1 << " Z mix NOT found" << std::endl;}}
            }
            //for non-planar molecules look a bit deeper using 3 atoms
            //this part of the function should not be overdone, otherwise we detect way too many Cn>2 axes and the point group symmetry comes out wrong
            //that is why this part of the code is entered if no other Cn>2 was already found for each n
            if ((!Planar())&&(idx > 1)&&(!lookextra)) {
              for (size_t atm3 = 0; atm3 < natoms; ++atm3) {
                if ((atm3 != atm)&&(atm3 != atm2)&&(atm != atm2)) {
                  if ((atoms[atm] != atoms[atm2])||(atoms[atm3] != atoms[atm2])||(atoms[atm] != atoms[atm3])) {continue;}
                  rotgeom = newgeom;
                  //calculate the new middle point
                  axis_x = rotgeom(atm+1,1) + rotgeom(atm2+1,1) + rotgeom(atm3+1,1);
                  axis_y = rotgeom(atm+1,2) + rotgeom(atm2+1,2) + rotgeom(atm3+1,2);
                  axis_z = rotgeom(atm+1,3) + rotgeom(atm2+1,3) + rotgeom(atm3+1,3);
                  //renormalize new middle point
                  norm = axis_x*axis_x + axis_y*axis_y + axis_z*axis_z;
                  axis_x /= sqrt(norm);
                  axis_y /= sqrt(norm);
                  axis_z /= sqrt(norm);
                  if (Prints[0] >= 4) {std::cout << "------------ Cn mix (non-planar): " << axis_x << " " << axis_y << " " << axis_z << std::endl;}
                  //check whether the new mix is an axis of the referencial
                  if ((fabs(axis_x*axis_x - 1) < tolzero)&&(axis_y*axis_y < tolzero)&&(axis_z*axis_z < tolzero)) {continue;}
                  else if ((fabs(axis_y*axis_y - 1) < tolzero)&&(axis_x*axis_x < tolzero)&&(axis_z*axis_z < tolzero)) {continue;}
                  else if ((fabs(axis_z*axis_z - 1) < tolzero)&&(axis_y*axis_y < tolzero)&&(axis_x*axis_x < tolzero)) {continue;}
                  existCN = false;
                  //comparing new mix vector with old ones
                  for (size_t idn = 1; idn < (CNadd.size()/3 + 1); ++idn) {
                    if (Prints[0] >= 5) {
                      std::cout << "For Cn non-planar";
                      std::cout << axis_x << "(" << CNadd[3*idn-3] << ") " << axis_y << "(" << CNadd[3*idn-2] << ") ";
                      std::cout << axis_z << "(" << CNadd[3*idn-1] << ") " << std::endl;
                    }
                    //tolz = 1000*tolzero;
                    //std::cout << "tolerance of " << tolz << std::endl;
                    aux = fabs(axis_x - CNadd[3*idn-3]);
                    MIN = fabs(axis_y - CNadd[3*idn-2]);
                    MAX = fabs(axis_z - CNadd[3*idn-1]);
                    if ((aux < tolz)&&(MIN < tolz)&&(MAX < tolz)) {existCN = true;}
                    aux = fabs(axis_x + CNadd[3*idn-3]);
                    MIN = fabs(axis_y + CNadd[3*idn-2]);
                    MAX = fabs(axis_z + CNadd[3*idn-1]);
                    if ((aux < tolz)&&(MIN < tolz)&&(MAX < tolz)&&(!existCN)) {existCN = true;}
                    norm = axis_x*CNadd[3*idn-3] + axis_y*CNadd[3*idn-2] + axis_z*CNadd[3*idn-1];
                    if (fabs(fabs(norm) - 1) < 10*tolz) {existCN = true;}
                    if (existCN) {break;}
                  }
                  if (!existCN) {
                    //making the mix vector rotation so that it becomes (0,0,1)
                    //getting the rotation angles theta (axis_x -> 0) and phi (axis_y -> 0)
                    theta = atan(axis_x/axis_y);
                    aux = axis_x*sin(theta) + axis_y*cos(theta);
                    phi = atan(aux/axis_z);
                    //rotation theta
                    testgeom = rotgeom*Givens(1,2,3,cos(theta),sin(theta));
                    //rotation phi
                    rotgeom = testgeom*Givens(2,3,3,cos(phi),sin(phi));
                    //recenter ??????????????????????
                    //CMg(rotgeom);
                    //test mix Cn Z
                    if (Prints[0] >= 4) {std::cout << "------------*** C" << idx+1 << " Z mix (non-planar)" << std::endl;}
                    testgeom = rotgeom*cnz;
                    if (CompGeom(rotgeom,testgeom,tolcn)) {
                      cnaxes[idx] += 1;
                      CNadd.push_back(axis_x);
                      CNadd.push_back(axis_y);
                      CNadd.push_back(axis_z);
                      if (idx > 1) {lookextra = true;}
                      if (Prints[0] >= 4) {std::cout << "C" << idx+1 << " Z mix (non-planar) found" << std::endl;}
                    }
                    else {if (Prints[0] >= 4) {std::cout << "C" << idx+1 << " Z mix (non-planar) NOT found" << std::endl;}}
                    //test mix Cn Y
                    if (Prints[0] >= 4) {std::cout << "------------*** C" << idx+1 << " Y mix (non-planar)" << std::endl;}
                    testgeom = rotgeom*cny;
                    if (CompGeom(rotgeom,testgeom,tolcn)) {
                      cnaxes[idx] += 1;
                      CNadd.push_back(axis_x);
                      CNadd.push_back(axis_y);
                      CNadd.push_back(axis_z);
                      if (idx > 1) {lookextra = true;}
                      if (Prints[0] >= 4) {std::cout << "C" << idx+1 << " Y mix (non-planar) found" << std::endl;}
                    }
                    else {if (Prints[0] >= 4) {std::cout << "C" << idx+1 << " Y mix (non-planar) NOT found" << std::endl;}}
                    //test mix Cn X
                    if (Prints[0] >= 4) {std::cout << "------------*** C" << idx+1 << " X mix (non-planar)" << std::endl;}
                    testgeom = rotgeom*cnx;
                    if (CompGeom(rotgeom,testgeom,tolcn)) {
                      cnaxes[idx] += 1;
                      CNadd.push_back(axis_x);
                      CNadd.push_back(axis_y);
                      CNadd.push_back(axis_z);
                      if (idx > 1) {lookextra = true;}
                      if (Prints[0] >= 4) {std::cout << "C" << idx+1 << " X mix (non-planar) found" << std::endl;}
                    }
                    else {if (Prints[0] >= 4) {std::cout << "C" << idx+1 << " X mix (non-planar) NOT found" << std::endl;}}
                  }
                }
              }
            }
          }
        }
        CNaxes.push_back(CNadd);
        CNadd.clear();
        SNaxes.push_back(SNadd);
        SNadd.clear();
      }
      //showing axes symmetry (Cn)
      if (Prints[0] >= 1) {
        for (size_t idx = 0; idx < maxcn; ++idx) {
          std::cout << "number of C" << idx+1 << " axes is " << cnaxes[idx] << std::endl;
        }
        if (Prints[0] >= 3) {
          size_t counter = 0;
          for (size_t idx = 0; idx < CNaxes.size(); ++idx) {
            counter = 1;
            std::cout << "C" << idx+1 << " vectors:" << std::endl;
            for (size_t idx2 = 0; idx2 < CNaxes[idx].size(); ++idx2) {
              std::cout << CNaxes[idx][idx2] << " ";
              ++counter;
              if (counter == 4) {
                std::cout << std::endl;
                counter = 1;
              }
            }
          }
        }
        std::cout << "------------------------------" << std::endl << std::endl;
      }
      //showing improper axes symmetry (Sn)
      if (Prints[0] >= 1) {
        for (size_t idx = 0; idx < maxsn; ++idx) {
          std::cout << "number of S" << idx+1 << " axes is " << snaxes[idx] << std::endl;
        }
        if (Prints[0] >= 3) {
          size_t counter = 0;
          for (size_t idx = 0; idx < SNaxes.size(); ++idx) {
            counter = 1;
            std::cout << "S" << idx+1 << " vectors:" << std::endl;
            for (size_t idx2 = 0; idx2 < SNaxes[idx].size(); ++idx2) {
              std::cout << SNaxes[idx][idx2] << " ";
              ++counter;
              if (counter == 4) {
                std::cout << std::endl;
                counter = 1;
              }
            }
          }
        }
        std::cout << "------------------------------" << std::endl << std::endl;
      }
      //showing planes symmetry 
      size_t highercn = 0; //highercn is principal axis
      for (size_t idx = 1; idx < maxcn; ++idx) {
        if (cnaxes[idx] > 0) {highercn = idx+1;}
      }
      size_t sigmahs = sigmaHorizontal(highercn,CNaxes,SIGMAplanes);
      size_t sigmavs = sigmaVertical(highercn,CNaxes,SIGMAplanes);
      size_t sigmads = 0;
      if ((sigmahs + sigmavs) < sigmaplanes) {sigmads = sigmaplanes - sigmahs - sigmavs;}
      if (Prints[0] >= 1) {
        std::cout << "number of sigma planes is " << sigmaplanes << std::endl;
        std::cout << "sigma planes" << std::endl;
        size_t counter = 1;
        std::cout << "Using C" << highercn << " to look for types of sigma" << std::endl;
        std::cout << "number sigma_h is " << sigmahs << std::endl;
        std::cout << "number sigma_v is " << sigmavs << std::endl;
        std::cout << "number sigma_d is " << sigmads << std::endl;
        if (Prints[0] >= 3) {
          for (size_t idx = 0; idx < SIGMAplanes.size(); ++idx) {
            std::cout << SIGMAplanes[idx] << " ";
            ++counter;
            if (counter == 4) {
              std::cout << std::endl;
              counter = 1;
            }
          }
        }
        std::cout << "------------------------------" << std::endl << std::endl;
      }
      if (Prints[0] >= 1) {
        std::cout << "Inversion Point: " << invpoint << std::endl;
        std::cout << "------------------------------" << std::endl << std::endl;
      }
      //defining symmetry point group based on the allowed operations and their number
      //I use a hack to get rid of the excessive S2 axes; if n_S2 > n_C2 then I can ignore the S2 axes
      bool cngc2 = false;                         //are there 2 or more Cn, n > 2?
      for (size_t idx = 2; idx < maxcn; ++idx) {
        if (cnaxes[idx] > 1) {
          cngc2 = true;
          break;
        }
      }
      if (cngc2) {                                            //groups with m_Cn > m_C2, n > 2 (Td, Oh, Ih)
        if (Prints[0] >= 2) {std::cout << "There are 2 or more Cn with n > 2 (Td, Oh, Ih)" << std::endl;}
        if (invpoint) {                                           //has inversion center (Oh, Ih)
          if (Prints[0] >= 2) {std::cout << "There is an inversion point (Oh, Ih)" << std::endl;}
          if (cnaxes[4] != 0) {                                       //has C5 (Ih)
            if (Prints[0] >= 2) {std::cout << "There is a C5" << std::endl;}
            symm = "Ih";
          }
          else {                                                      //has not C5 (Oh)
            if (Prints[0] >= 2) {std::cout << "There is no C5" << std::endl;}
            symm = "Oh";
          }
        }
        else {                                                    //no inversion center
          if (Prints[0] >= 2) {std::cout << "There is no inversion point" << std::endl;}
          symm = "Td";
        }
      }
      else {                                                 //groups without m_Cn > m_C2, n > 2 (C1, Ci, Cs, Cn, Cnh, Cnv, Dn, Dnd, Dnh, S2n)
        if (Prints[0] >= 2) {std::cout << "There are not 2 or more Cn with n > 2 (C1, Ci, Cs, Cn, Cnh, Cnv, Dn, Dnd, Dnh, S2n)" << std::endl;}
        size_t hascn = 0; //hascn is principal axis
        for (size_t idx = 1; idx < maxcn; ++idx) {
          if (cnaxes[idx] > 0) {hascn = idx+1;}
        }
        if (hascn > 1) {                                          //has any Cn, n > 1 (Cn, Cnh, Cnv, Dn, Dnd, Dnh, S2n)
          if (Prints[0] >= 2) {std::cout << "There is a Cn, n >= 2" << std::endl;}
          size_t nc2eqcn = 0;
          for (size_t idx = 2; idx < maxcn; ++idx) {
            if (cnaxes[1] >= (idx+1)) {nc2eqcn = idx+1;}
          }
          if (nc2eqcn != 0) {                                         //n_C2 is at least Cn (Dn, Dnd, Dnh)
            if (Prints[0] >= 2) {std::cout << "There are nC2 and Cn" << std::endl;}
            if (sigmads >= hascn) {                                       //has "hascn" sigmas diagonal (Dnd)
              if (Prints[0] >= 2) {std::cout << "There are at least " << hascn << " sigma d" << std::endl;}
              symm = "D"+to_string(hascn)+"d";
            }
            else if (sigmahs != 0) {                                      //has sigma horizontal to principal axis (Dnh)
              if (Prints[0] >= 2) {std::cout << "There are sigma h" << std::endl;}
              symm = "D"+to_string(hascn)+"h";
            }
            else {                                                        //(Dn)
              if (Prints[0] >= 2) {std::cout << "There are no sigma h and there not at least " << hascn << " sigma d" << std::endl;}
              symm = "D"+to_string(hascn);
            }
          }
          else {                                                      //n_C2 is not at least Cn (Cn, Cnh, Cnv, S2n)
            if (Prints[0] >= 2) {std::cout << "There are mC2 and Cn, m < n" << std::endl;}
            bool hass2n = false;
            for (size_t idx = 0; idx < maxsn; ++idx) {
              if (((idx+1)%2 == 0)&&(snaxes[idx] != 0)) {
                hass2n = true;
                break;
              }
            }
            if (sigmavs >= hascn) {                                       //has "hascn" sigma verticals (Cnv)
              if (Prints[0] >= 2) {std::cout << "There are at least " << hascn << " sigma v" << std::endl;}
              symm = "C"+to_string(hascn)+"v";
            }
            else if (sigmahs != 0) {                                      //has sigma horizontal to principal axis (Cnh)
              if (Prints[0] >= 2) {std::cout << "There are sigma h" << std::endl;}
              symm = "C"+to_string(hascn)+"h";
            }
            else if (hass2n) {                                            //(S2n)
              if (Prints[0] >= 2) {std::cout << "There are S2n (even)" << std::endl;}
              symm = "S"+to_string(hascn);
            }
            else {                                                        //(Cn)
              if (Prints[0] >= 2) {std::cout << "There are no sigma h, there are no S2n and there not at least " << hascn << " sigma v" << std::endl;}
              symm = "C"+to_string(hascn);
            }
          }
        }
        else {                                                    //has no Cn, n > 1 (C1, Ci, Cs)
          if (Prints[0] >= 2) {std::cout << "There is no Cn" << std::endl;}
          if (sigmaplanes != 0) {                                     //has sigma planes (Cs)
            if (Prints[0] >= 2) {std::cout << "There is a sigma" << std::endl;}
            symm = "Cs";
          }
          else {                                                      //has no sigma planes (C1, Ci)
            if (Prints[0] >= 2) {std::cout << "There is no sigma" << std::endl;}
            if (invpoint) {                                               //has inversion center (Ci)
              if (Prints[0] >= 2) {std::cout << "There is an inversion point" << std::endl;}
              symm = "Ci";
            }
            else {                                                        //has no inversion center (C1)
              if (Prints[0] >= 2) {std::cout << "There are no elements of symmetry" << std::endl;}
              symm = "C1";
            }
          }
        }
      }
      
    }
    //std::cout << "Point group symmetry found was " << symm << std::endl;
  }
};
unsigned int PointGroup2Sigma(std::string pg) {
  //http://symmetry.jacobs-university.de
  //pg is point group as a string
  int _sigma = 0;
  if (pg == "DINFh") {_sigma = 2;}
  else if (pg == "CINFv") {_sigma = 1;}
  else if (pg == "C2h") {_sigma = 2;}
  else if (pg == "C3h") {_sigma = 3;}
  else if (pg == "C4h") {_sigma = 4;}
  else if (pg == "C5h") {_sigma = 5;}
  else if (pg == "C6h") {_sigma = 6;}
  else if (pg == "C2v") {_sigma = 2;}
  else if (pg == "C3v") {_sigma = 3;}
  else if (pg == "C4v") {_sigma = 4;}
  else if (pg == "C5v") {_sigma = 5;}
  else if (pg == "C6v") {_sigma = 6;}
  else if (pg == "C7v") {_sigma = 7;}
  else if (pg == "C8v") {_sigma = 8;}
  else if (pg == "C1") {_sigma = 1;}
  else if (pg == "C2") {_sigma = 2;}
  else if (pg == "C3") {_sigma = 3;}
  else if (pg == "C4") {_sigma = 4;}
  else if (pg == "C5") {_sigma = 5;}
  else if (pg == "C6") {_sigma = 6;}
  else if (pg == "C7") {_sigma = 7;}
  else if (pg == "C8") {_sigma = 8;}
  else if (pg == "Cs") {_sigma = 1;}
  else if (pg == "Ci") {_sigma = 1;}
  else if (pg == "D2h") {_sigma = 4;}
  else if (pg == "D3h") {_sigma = 6;}
  else if (pg == "D4h") {_sigma = 8;}
  else if (pg == "D5h") {_sigma = 10;}
  else if (pg == "D6h") {_sigma = 12;}
  else if (pg == "D7h") {_sigma = 14;}
  else if (pg == "D8h") {_sigma = 16;}
  else if (pg == "D2d") {_sigma = 4;}
  else if (pg == "D3d") {_sigma = 6;}
  else if (pg == "D4d") {_sigma = 8;}
  else if (pg == "D5d") {_sigma = 10;}
  else if (pg == "D6d") {_sigma = 12;}
  else if (pg == "D7d") {_sigma = 14;}
  else if (pg == "D8d") {_sigma = 16;}
  else if (pg == "D2") {_sigma = 4;}
  else if (pg == "D3") {_sigma = 6;}
  else if (pg == "D4") {_sigma = 8;}
  else if (pg == "D5") {_sigma = 10;}
  else if (pg == "D6") {_sigma = 12;}
  else if (pg == "D7") {_sigma = 14;}
  else if (pg == "D8") {_sigma = 16;}
  else if (pg == "T") {_sigma = 12;}
  else if (pg == "Th") {_sigma = 12;}    //use the same as T
  else if (pg == "Td") {_sigma = 12;}
  else if (pg == "O") {_sigma = 24;}     //use the same as Oh
  else if (pg == "Oh") {_sigma = 24;}
  else if (pg == "I") {_sigma = 60;}     //use the same as Ih
  else if (pg == "Ih") {_sigma = 60;}
  else if (pg == "S2") {_sigma = 1;}
  else if (pg == "S4") {_sigma = 2;}
  else if (pg == "S6") {_sigma = 3;}
  else if (pg == "S8") {_sigma = 4;}
  else if (pg == "S10") {_sigma = 5;}
  else if (pg == "S12") {_sigma = 6;}
  return _sigma;
}

#endif //_Molecule_
