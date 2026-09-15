#include <iostream>
#include "src/Molecule.hpp"
#include "src/Patterns.hpp"
#include "src/other/AddProtons.hpp"

int main(int argc, char** argv) {
  
  Molecule Mol(argv[1],0,1,"C1",false);
  std::vector<size_t> atoms = Mol.Atoms();
  
  Patterns VSEPRclass(Mol);
  VSEPRclass.ProtonateFromSDFMOL2(argv[1]);
  
  //write SDF
  VSEPRclass.WriteProtonatedSDF(argv[2],"UNL");
  
  Mol = VSEPRclass.GetMolecule();
  Mol.WriteXYZ(argv[2]);
  
  std::cout << VSEPRclass.MolecularCharge() << std::endl;
  
  return 0;
}
