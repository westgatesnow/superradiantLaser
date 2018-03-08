#ifndef __SUPERRADIANTLASER__HPP__
#define __SUPERRADIANTLASER__HPP__
//This program is used to simulate the superradiant laser 
//  --using the cumulant method 
//  --without the cavity variables
//  --using the individual variables.

//Include packages
#include <vector>
#include <complex>
#include <iostream>
#include <fstream>
#include <cmath>
#include <getopt.h>
#include <time.h>
#include <stdlib.h> 
//Include Eigen package;
#include <Eigen/Core>
#include <Eigen/StdVector>
#include <Eigen/Eigenvalues>
using namespace Eigen;

//Include and define RNG;
#include "RNG.hpp"
RNG rng(time(NULL));


//Simulation parameters
typedef struct Param {
  //simulation specification
  double dt; 
  double tmax;
  int nstore; // number of times to store observables
  //beam parameters
  int nAtom; // number of intracavity atoms
  double gammac;    //collective decay rate
  double repumping; //repumping rate w
  //Other parameters
  std::string name; //name of the directory to store results

  //Set up initial values of the parameters
  Param() : dt(1.0e-3), tmax(1), nstore(100), nAtom(10), gammac(0.1e0), repumping(10), name("abracadabra")  {}
} Param;

std::ostream& operator<< (std::ostream& o,const Param& s)
{
  o << s.dt << std::endl;
  o << s.tmax << std::endl;
  o << s.nstore << std::endl;
  o << s.nAtom << std::endl;
  o << s.gammac << std::endl;
  o << s.repumping << std::endl;

  return o;
}

typedef struct Observables {
  Observables(const int n/*,const int m*/) : intensity(n), intensityUnCor(n),
                                          inversion(n), spinSpinCor(n)
                                          //,g1(m)
  {}
  VectorXd intensity;
  VectorXd intensityUnCor;
  VectorXd inversion;
  VectorXd spinSpinCor;
  //VectorXd g1;
} Observables;

typedef struct ObservableFiles {
  ObservableFiles() : intensity("intensity.dat"), intensityUnCor("intensityUnCor.dat"), 
                      inversion("inversion.dat"), spinSpinCor("spinSpinCor.dat")
                      //,g1("g1.dat")
  {}
  ~ObservableFiles() {
    intensity.close();
    intensityUnCor.close();
    inversion.close();
    spinSpinCor.close();
    //g1.close();
  }
  std::ofstream intensity, intensityUnCor, inversion, spinSpinCor;//, g1;
} ObservableFiles;

#endif