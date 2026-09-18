#include <iostream>
#include "src/Molecule.hpp"
#include "src/Patterns.hpp"
#include "src/other/AddProtons.hpp"

int main(int argc, char** argv) {
  
  Molecule Mol(argv[1],0,1,"C1");
  std::vector<size_t> atoms = Mol.Atoms();
  
  Patterns VSEPRclass(Mol);
  VSEPRclass.ProtonateFromCDF(argv[2],argv[3]);
  VSEPRclass.WriteSDF(argv[4],argv[4],-1);
  
  Mol = VSEPRclass.GetMolecule();
  Mol.WriteXYZ(argv[4]);
  
  std::cout << VSEPRclass.MolecularCharge() << std::endl;
  
  return 0;
}
