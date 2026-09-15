#include "src/GFN.hpp"
#include "src/MolecularDynamics.hpp"
#include "src/math/SolverPackage.hpp"
#include <stdlib.h>

//program to run GFN2-xTB calculations

int main(int argc, char** argv) {
  
  //arguments
  // 0 exe
  // 1 geometry
  // 2 charge
  // 3 Telec
  // 4 solvation?
  // 5 solvent name
  // 6 name of density file
  
  //parameters passed as argument
  char *p;
  int charge = strtol(argv[2],&p,10);
  char *q;
  double Telec = strtod(argv[3],&q);
  char *v;
  int solvation = strtol(argv[4],&v,10);

  std::cout << "running " << argv[1] << "\n";
  std::cout << "charge     = " << charge << std::endl;

  Molecule Mol1(argv[1],charge,1,"C1",false);

  std::vector<size_t> atoms = Mol1.Atoms();
  size_t Natoms = atoms.size();
  matrixE Geometry = Mol1.Geometry();

  BSet basis(Mol1,"gfn2");
  std::vector<size_t> AOS = basis.AtomNAOs(atoms);

  GFN2 electron(basis,Mol1);
  electron.setElectronTemp(Telec);
  
  //use ALPB solvation?
  if (solvation > 0) {
    electron.setSolvent(argv[5]);
  }

  electron.setRestart(0);                                      //use restart, otherwise this will be a killer
  electron.Calculate(0);
  
  //get the density
  electron.ElectronicDensity(argv[6]);
  
  return 0;
}