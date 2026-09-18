#include <iostream>
#include "src/GFN.hpp"
#include "src/EnergyDecompositionAnalysis.hpp"

int main(int argc, char** argv) {

  char *p;
  int chargeAB = strtol(argv[2],&p,10);
  char *q;
  int chargeA = strtol(argv[4],&q,10);
  char *r;
  int chargeB = strtol(argv[6],&r,10);
  
  std::cout << "running " << argv[1] << "  " << chargeAB << std::endl;
  std::cout << "        " << argv[3] << "  " << chargeA << std::endl;
  std::cout << "        " << argv[5] << "  " << chargeB << std::endl;

  Molecule molecule(argv[1],chargeAB,1,"C1",false);
  Molecule monomerA(argv[3],chargeA,1,"C1",false);
  Molecule monomerB(argv[5],chargeB,1,"C1",false);

  BSet basis(molecule,"gfn2");
  GFN2 electron(basis,molecule);
  electron.setSolvent("water");
  //electron.Calculate(0);
  
  EnergyDecomposition EDA;
  EDA.setPopulationConstraint(true);
  EDA.setFastEDDA(true);
  EDA.setElectronicModel(electron);
  EDA.addMolecules(molecule,monomerA,monomerB);
  EDA.DecomposeEnergy(true);
  
  return 0;
}
