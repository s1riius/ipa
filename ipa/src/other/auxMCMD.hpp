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

#ifndef _Auxiliary_Package_Monte_Carlo_Molecular_Dynamics_
#define _Auxiliary_Package_Monte_Carlo_Molecular_Dynamics_
#include <vector>
#include "../math/RandomPackage.hpp"

//description:
//auxiliary functions for monte carlo and molecular dynamics

//selection of methods to be used in monte carlo and molecular dynamics
void DistributeParticles(std::vector<std::vector<double> > * particlepos, std::vector<double> * boxdim, int Nparticles, int ndim, double sigma, int maxiter = 200) {
  //function distributing particles in a box and returning a vector containing the positions of all particles
  bool distribute = false;          //variable to control particle distribution inside box
  bool overlap;                     //to check whether particles overlap
  int iter = 0;
  std::vector<double> position(ndim,0.0);
  double distance;
  int counter = 0;
  while (!(distribute)) {
    ++counter;
    //get position of first particle; this must always work
    for (size_t idx = 0; idx < ndim; ++idx) {
      position[idx] = fRandom(sigma,boxdim->at(idx) - sigma);
    }
    particlepos->push_back(position);
    //place other particles
    for (size_t idp = 1; idp < Nparticles; ++idp) {
      //placing particle idp
      for (size_t idx = 0; idx < ndim; ++idx) {
        position[idx] = fRandom(sigma,boxdim->at(idx) - sigma);
      }
      //check for overlap between particles
      for (size_t idx = 0; idx < particlepos->size(); ++idx) {
        distance = 0.0;
        overlap = false;
        //loop over coordinates
        for (size_t idy = 0; idy < ndim; ++idy) {
          distance += (position[idy] - particlepos->at(idx)[idy])*(position[idy] - particlepos->at(idx)[idy]);
        }
        if (distance < 4.0*sigma*sigma) {
          overlap = true;
          break;
        }
      }
      if (overlap) {
        //tabulae rasa
        distribute = false;
        particlepos->clear();
        break;
      }
      else {
        distribute = true;
        particlepos->push_back(position);
      }
    }
    ++iter;
    if (iter == maxiter) break;
  }
  //std::cout << "counter = " << counter << std::endl;
  if (iter == maxiter) {std::cout << "DistributeParticles: calculation failed, as could not place all particles in box; restart with larger box or less particles" << std::endl;}
}
double PotEnergy(int potential, std::vector<double> * param, double threshold = 1.0e-8) {
  //function that calculates potential energies, according to severa types of potentials, using for that a parameter vector
  //potential defines the type of potential:       0 -> hard-sphere
  //potential defines the type of potential:       1 -> Ising
  double Epot = 0.0;
  if (potential == 0) {
    //param[0] -> sigma, the particle diameter; param[1] -> distance
    if (param->at(1) < param->at(0)) {Epot = 1.0e30;}
  }
  else if (potential == 1) {
    //param[0] -> sigma, the particle diameter; param[1] -> distance; param[2] -> J
    if (param->at(1) < param->at(0)) {Epot = 1.0e30;}
    else if (fabs(param->at(1) - 1.0) < threshold) {Epot = -param->at(2);}
  }
  return Epot;
}
#endif //_Auxiliary_Package_Monte_Carlo_Molecular_Dynamics_
