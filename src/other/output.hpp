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

#ifndef _Output_Package_
#define _Output_Package_
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <stdio.h>

//description:
//functions for output

//interface to output functions
void MD2file(bool init, size_t Nparticles, size_t ndim, double sigma, bool monochromatic, double timestamp, std::vector<std::vector<double> > * positions, std::vector<std::vector<double> > * velocities, std::vector<double> boxlim, std::string datafile, std::string folder) {
  //function preparing a string to be written to file, passed to python so that I remain OS independent and do not use boost
  //Interface from molecular dynamics to file
  std::string write2file;
  std::string initialize;
  if (init) {
    initialize = "true";
    write2file = std::to_string(Nparticles)+"\n"+std::to_string(sigma)+"\n";
    write2file += "boxlimits: ";
    for (size_t idx = 0; idx < boxlim.size();++idx) {
      write2file += std::to_string(boxlim[idx])+" ";
    }
    write2file += "\n";
  }
  else {
    initialize = "false";
    write2file = "t = "+std::to_string(timestamp)+"\n";
    for (size_t idpart = 0; idpart < Nparticles; ++idpart) {
      write2file += std::to_string(idpart + 1)+" ; t = "+std::to_string(timestamp)+";  r = ";
	    for (size_t idcoord = 0; idcoord < ndim; ++idcoord) {
	      write2file += std::to_string(positions->at(idpart)[idcoord]);
	      if (idcoord != ndim - 1) {write2file += ",";}
	    }
	    write2file += "\n";
	    write2file += std::to_string(idpart + 1)+" ; t = "+std::to_string(timestamp)+";  v = ";
	    for (size_t idcoord = 0; idcoord < ndim; ++idcoord) {
	      write2file += std::to_string(velocities->at(idpart)[idcoord]);
	      if (idcoord != ndim - 1) {write2file += ",";}
	    }
	    write2file += "\n-------\n";
    }
    write2file += "===========\n";
  }
  //send to python
  std::string command = "python other/string2file.py "+datafile+" "+folder+" '"+write2file+"' "+initialize;
  FILE * stream = popen(command.c_str(),"r");
  pclose(stream);
}
template<class T>
void LatticeFile(bool init, size_t Nparticles, size_t ndim, double sigma, double temp, double energy, std::vector<std::vector<T> > * position, std::string datafile, std::string folder) {
  //init determines what to do with file: 1 -> initialize; 0 -> adds update
  std::string write2file = "";
  FILE * pFile;
  std::string strng = folder+"/"+datafile;
  if (init) {
    write2file += "Nparticles: "+std::to_string(Nparticles)+"\nsigma: "+std::to_string(sigma)+"\nndim: "+std::to_string(ndim)+"\ntemperature: "+std::to_string(temp)+"\n";
    pFile = fopen(strng.c_str(),"w");
  }
  else {
    pFile = fopen (strng.c_str(),"a");
  }
  write2file += "event: "+std::to_string(energy)+";";
  for (size_t ipart = 0; ipart < Nparticles; ++ipart) {
    write2file += std::to_string(ipart)+",";
    for (size_t idim = 0; idim < ndim; ++idim) {
      write2file += std::to_string(position->at(ipart)[idim]);
      if (idim < ndim - 1) {write2file += ",";}
    }
    write2file += ";";
  }
  write2file += "\n";
  if (pFile!=NULL)
  {
    fputs (write2file.c_str(),pFile);
    fclose (pFile);
  }
}
template<class T>
void WriteConfig(size_t Nparticles, size_t ndim, double sigma, std::vector<std::vector<T> > * position, std::string datafile, std::string folder) {
  std::string write2file;
  FILE * pFile;
  std::string strng = folder+"/"+datafile;
  write2file = "Nparticles: "+std::to_string(Nparticles)+"\nsigma: "+std::to_string(sigma)+"\nndim: "+std::to_string(ndim)+"\n";
  pFile = fopen(strng.c_str(),"w");
  for (size_t idpart = 0; idpart < position->size(); ++idpart) {
    for (size_t idim = 0; idim < ndim; ++idim) {
      write2file += std::to_string(position->at(idpart)[idim]);
      if (idim < ndim - 1) {write2file += ",";}
    }
    write2file += "\n";
  }
  if (pFile!=NULL)
  {
    fputs (write2file.c_str(),pFile);
    fclose (pFile);
  }
}

#endif //_Output_Package_
