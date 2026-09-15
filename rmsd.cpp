#include "src/Molecule.hpp"
#include "src/MolecularDynamics.hpp"
#include <stdlib.h>

//program to align two structures

int main(int argc, char** argv) {
  
  //parameters passed as argument
  //argv[1] reference molecule
  //argv[2] molecule to reorient
  //argv[3] where to save reoriented molecule

  Molecule Mol1(argv[1],0,1,"C1",false);
  Molecule Mol2(argv[2],0,1,"C1",false);
  
  matrixE geom1 = Mol1.Geometry();
  matrixE geom2 = Mol2.Geometry();
  matrixE gradients = geom1;
  
  MetaDynamics MTD;
  
  std::cout << "RMSD: " << MTD.RMSD(geom1,geom2,gradients);
  
  return 0;
}