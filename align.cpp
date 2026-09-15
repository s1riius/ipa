#include "src/Molecule.hpp"
#include <stdlib.h>

//program to align two structures

int main(int argc, char** argv) {
  
  //parameters passed as argument
  //argv[1] reference molecule
  //argv[2] molecule to reorient
  //argv[3] where to save reoriented molecule

  Molecule Reference(argv[1],0,1,"C1",false);
  Molecule Modify(argv[2],0,1,"C1",false);
  
  Modify.AlignTo(Reference);
  Modify.WriteXYZ(argv[3]);
  
  return 0;
}