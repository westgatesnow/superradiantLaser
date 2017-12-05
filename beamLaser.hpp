#ifndef __BEAMLASER__HPP__
#define __BEAMLASER__HPP__
//This program is used to simulate the beam laser using the cumulant expansion method.

//Avoid warnings on Eigen package
//#pragma GCC diagnostic ignored "-Wignored-attributes"
//#pragma GCC diagnostic ignored "-Wmisleading-indentation"
//#pragma GCC diagnostic ignored "-Wdeprecated-declarations"

//Include Eigen package
//Work in Eigen namespace
#include <Eigen/Core>
#include <Eigen/StdVector>
#include <Eigen/Eigenvalues>
using namespace Eigen;

//Include standard packages
#include <vector>
#include <complex>
#include <iostream>
#include <fstream>
#include <cmath>
#include <getopt.h>
#include <time.h>

#define NVAR 3 // 3 variables; sigmax, sigmay, sigmaz

//Include and define RNG
#include "RNG.hpp"
RNG rng(time(NULL));

//Define the complex I and ONE
static const std::complex<double> I = std::complex<double>(0.0,1.0);
static const std::complex<double> ONE = std::complex<double>(1.0,0.0);

//Data type for atoms
typedef struct {
  Vector3d X;     //position
  Vector3d P;     //momentum. We suppose mass is one, so momentum is velocity.
} Atom;

typedef struct {
  MatrixXd cov;       //The covariance matrix between \sigma^+_j and \sigma^-_k. It is real and symmetric.
  MatrixXd covZ;      //The covariance matrix between \sigma^z_j and \sigma^z_k. Its initial value is identity. We only care about its off-diagonal elements.
  VectorXd inv;       //The population inversion \sigma^z_j.
} Internal;

typedef struct {
  std::vector<Atom> atoms;  //The position and velocity of each atom.
  Internal internal;        //Two covariance matrices and a vector of population inversion of all the atoms in the cavity.
} Ensemble;

typedef struct {
  const char* configFile;
} CmdLineArgs;

const char* usageHeader = "\nBeam Laser Simulation.\n";
const char* usageMessage =
  "\n"
  "Usage:         "
  "beamLaser "
  "--file"
  "\n"
  "--file, -f     : Configuration file to setup the simulation\n"
  "--help, -h     : Print this usage message\n"
  "\n\n";

//Simulation parameters
typedef struct Param {
  //simulation specification
  double meanAtomNumber; //the desired average number of intracavity atoms N0
  int nstore; // number of times to store observables
  int steadyTime;  // tmax/tau
  //beam parameters
  double yWall; //position of the wall where atoms are destroyed
                //The coordinated are chosen s.t. atoms are created at -yWall.
                //The walls are assumed to be in xz plane.
  Vector2d sigmaX;  //standard deviation of position in xz. y deviation is 
                    //taken care of by the Poisson distribution.
  double transitTime;     //the transit time tau1>0, unit 1/gammaC.
  Vector3d sigmaP;  //standard deviation of momentum
  int meanAtomGeneratingNumber;   //the mean number of atoms dN generated in time dt, dN << N0;
  double gammac;    //collective decay rate

  //Other parameters
  std::string name; //name of the directory to store results

  //Set up initial values of the parameters
  Param() : meanAtomNumber(100), nstore(10), steadyTime(10), yWall(5.0e0), 
    sigmaX(0.0e0,0.0e0), transitTime(1.0e0),
    sigmaP(0.0e0,0.0e0,0.0e0), meanAtomGeneratingNumber(1), gammac(0.1e0), name("abracadabra")  {}
} Param;

std::ostream& operator<< (std::ostream& o,const Param& s)
{
  o << s.meanAtomNumber << std::endl;
  o << s.steadyTime << std::endl;
  o << s.nstore << std::endl;
  o << s.yWall << std::endl;
  o << s.sigmaX << std::endl;
  o << s.transitTime << std::endl;
  o << s.sigmaP << std::endl;
  o << s.meanAtomGeneratingNumber << std::endl;
  o << s.gammac << std::endl;

  return o;
}

typedef struct Observables {
  Observables(const int n) : nAtom(n), intensity(n), inversion(n), spinSpinCor(n)
  {}
  Matrix <unsigned long int, 1, Dynamic> nAtom; 
  VectorXd intensity;
  VectorXd inversion;
  VectorXd spinSpinCor;
} Observables;

typedef struct RawData {
  RawData() : spinSpinCovSS(), invSS()
  {}
  MatrixXd spinSpinCovSS;
  VectorXd invSS;
} RawData;

typedef struct ObservableFiles {
  ObservableFiles() : nAtom("nAtom.dat"), intensity("intensity.dat"), 
                  inversion("inversion.dat"), spinSpinCor("spinSpinCor.dat"),
                  spinSpinCovSS("spinSpinCovSS.dat"), invSS("invSS.dat")
  {}
  ~ObservableFiles() {
    nAtom.close();
    intensity.close();
    inversion.close();
    spinSpinCor.close();
    spinSpinCovSS.close();
    invSS.close();
  }
  std::ofstream nAtom, intensity, inversion, spinSpinCor, spinSpinCovSS, invSS;
} ObservableFiles;

#endif
