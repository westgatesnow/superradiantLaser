//This program is used to simulate the beam laser using the cumulant expansion method.
#include "beamLaser.hpp"

void getOptions(int argc, char** argv, CmdLineArgs* cmdLineArgs)
{
  cmdLineArgs->configFile="sampleSimulation.txt";
  while (1) {
    int c;
    static struct option long_options[] = {
      {"help", no_argument, 0, 'h'},
      {"file", required_argument, 0, 'f'},
      {0, 0, 0, 0}
    };
    int option_index = 0;
    c = getopt_long(argc, argv, "hf:", long_options, &option_index);
    if (c == -1) break;
    switch (c) {
      case 'h': std::cout << usageHeader << usageMessage;
        exit(0);
      case 'f': cmdLineArgs->configFile = optarg;
        break;
      default: exit(-1);
    }
  }
  if (optind < argc) {
    std::cout << "Error: non-option arguments: ";
    while (optind < argc) std::cout << argv[optind++] << " ";
    std::cout << std::endl;
    exit(-1);
  }
  std::cout << "Using parameters file " << cmdLineArgs->configFile << std::endl;
}

void getParam(const char* filename, Param *param) 
{
  std::ifstream configInput(filename);
  std::string dummy;

  while (!configInput.eof()) {
    configInput >> dummy;
    if (configInput.eof()) break;
    if (!configInput.good()) {
      std::cout << "Bad read in input file" << std::endl;
      exit(-1);
    }
    if (dummy.compare("meanAtomNumber") == 0)
      configInput >> param->meanAtomNumber;
    else if (dummy.compare("nstore") == 0)
      configInput >> param->nstore;
    else if (dummy.compare("steadyTime") == 0)
      configInput >> param->steadyTime;
    else if (dummy.compare("yWall") == 0)
      configInput >> param->yWall;
    else if (dummy.compare("sigmaXX") == 0)
      configInput >> param->sigmaX[0];
    else if (dummy.compare("sigmaXZ") == 0)
      configInput >> param->sigmaX[1];
    else if (dummy.compare("transitTime") == 0)
      configInput >> param->transitTime;
    else if (dummy.compare("sigmaPX") == 0)
      configInput >> param->sigmaP[0];
    else if (dummy.compare("sigmaPY") == 0)
      configInput >> param->sigmaP[1];
    else if (dummy.compare("sigmaPZ") == 0)
      configInput >> param->sigmaP[2];
    else if (dummy.compare("meanAtomGeneratingNumber") == 0)
      configInput >> param->meanAtomGeneratingNumber;
    else if (dummy.compare("gammac") == 0)
      configInput >> param->gammac;
    else if (dummy.compare("name") == 0)
      configInput >> param->name;
    else {
      std::cout << "Error: invalid label " << dummy << " in "
          << filename << std::endl;
      exit(-1);
    }
  }
}


void generateMotionalState(Ensemble& ensemble, const Param& param, const double meanP)
{
  Vector3d X (0, -param.yWall, 0);
  Vector3d P (0, meanP, 0);
/*  Vector3d X (
   rng.get_gaussian_rn (param.sigmaX [0]),
   param.yWallCreate,
    rng.get_gaussian_rn (param.sigmaX [1]));
  //Only keep atoms moving to the positive y direction.
*/ 

/*  Vector3d P (
    rng.get_gaussian_rn (param.sigmaP [0]),
    rng.get_gaussian_rn (param.sigmaP [1]) + meanP,
    rng.get_gaussian_rn (param.sigmaP [2]));
  while (P[1] <= 0)
    P[1] = rng.get_gaussian_rn (param.sigmaP [1]) + meanP;
*/
  Atom newAtom = {X,P};
  ensemble.atoms.push_back(newAtom);
}


void generateInternalState(const unsigned long int nAtom, Ensemble& ensemble)
{ 
  //generate cov
  unsigned long int size = ensemble.internal.cov.rows(); 
  unsigned long int newSize = size+nAtom;
  MatrixXd newCov = MatrixXd::Identity(newSize,newSize); 
  newCov.bottomRightCorner(size,size) = ensemble.internal.cov;
  
  //generate covZ
  MatrixXd newCovZ = MatrixXd::Identity(newSize,newSize); 
  newCovZ.bottomRightCorner(size,size) = ensemble.internal.covZ;

  //generate inv
  VectorXd newInv(newSize);
  newInv.fill(1);
  newInv.tail(size) = ensemble.internal.inv;

  //put back cov, covZ, and inv
  ensemble.internal.cov = newCov;
  ensemble.internal.covZ = newCovZ;
  ensemble.internal.inv = newInv;
}

void addAtomsFromSource(Ensemble& ensemble, const Param& param, const double meanP)
{
  unsigned long int nAtom;

/////////////////////////////////////
// No Beam Noise. N0 changes.
//   double mean = param.density*param.dt;
//   nAtom = mean;

/////////////////////////////////////
// No Beam Noise. N0 is const.
  const int N0 = param.meanAtomNumber;
  nAtom = param.meanAtomGeneratingNumber;

/////////////////////////////////////
// Beam Noise. N0 is const.
//  const int N0 = param.meanAtomNumber;
//  nAtom = rng.get_poissonian_int(param.meanAtomGeneratingNumber);      

  for (unsigned long int n = 0; n < nAtom; n++)
    generateMotionalState(ensemble, param, meanP);                //For each atom, generate its own x and p;
  
  generateInternalState(nAtom, ensemble);        //Generate the covariance matrix cov and the vector inv with expanded
                                                 // due to the new atoms.
}
 
void removeAtomsAtWalls(Ensemble& ensemble, const Param& param) 
{
  std::vector<Atom> newAtoms;
  for (std::vector<Atom>::iterator a = ensemble.atoms.begin(); a != ensemble.atoms.end(); a++)
    if (a->X[1] < param.yWall)
      newAtoms.push_back(*a);
  ensemble.atoms = newAtoms;

  unsigned long int newSize = newAtoms.size();
  ensemble.internal.cov.conservativeResize(newSize,newSize); 
  ensemble.internal.covZ.conservativeResize(newSize, newSize);
  ensemble.internal.inv.conservativeResize(newSize);
}

void advanceMotionalStateOneTimeStep(Ensemble& ensemble, const Param& param, const double dt) 
{
  for (std::vector<Atom>::iterator a = ensemble.atoms.begin(); a != ensemble.atoms.end(); a++)
    a->X += dt * a->P;
}

Internal RHS(const Internal& internal, const Param& param)
{
  Internal RHS(internal);
  double gc = param.gammac;
  
  //for convenience
  MatrixXd cov = internal.cov;
  MatrixXd covZ = internal.covZ;
  VectorXd inv = internal.inv;
  int size = inv.size();
  //right hand side of the DE set
  for (int i=0; i<size; i++) {
    //inv
    RHS.inv(i) = -2*gc*cov.colwise().sum()(i);
    //cov diagonal
    RHS.cov(i,i) = -gc*cov.colwise().sum()(i);

    for (int j=i+1; j<size; j++) {
      //cov off-diagonal
      RHS.cov(i,j) = gc/2*(inv(i)*cov.colwise().sum()(j)+inv(j)*cov.rowwise().sum()(i))
                  -gc/2*(inv(i)*inv(j)-covZ(i,j)+cov(i,j)*(inv(i)+inv(j)+2));
      RHS.cov(j,i) = RHS.cov(i,j);
      //covZ off-diagonal
      RHS.covZ(i,j) = 2*gc*(inv(i)*inv(j)-covZ(i,j))
                   +2*gc*cov(i,j)*(inv(i)+inv(j)+2)
                   -2*gc*(inv(i)*cov.colwise().sum()(j)+inv(j)*cov.rowwise().sum()(i));
      RHS.covZ(j,i) = RHS.cov(i,j);
    }
  }
  return RHS;
}

void advanceInternalStateOneTimeStep(Ensemble& ensemble, const Param& param, const double dt)
{
  double gc = param.gammac;
  
  //Define a new Internal type with initial value ensemble.internal
  Internal newInternal(ensemble.internal);
  //Using RK2 method.
  newInternal.inv += dt/2*RHS(ensemble.internal, param).inv;
  newInternal.cov += dt/2*RHS(ensemble.internal, param).cov;
  newInternal.covZ += dt/2*RHS(ensemble.internal, param).covZ;
  //Second round
  ensemble.internal.inv += dt*RHS(newInternal, param).inv;
  ensemble.internal.cov += dt*RHS(newInternal, param).cov;
  ensemble.internal.covZ += dt*RHS(newInternal, param).covZ;
}

void advanceAtomsOneTimeStep(Ensemble& ensemble, const Param& param, const double dt)
{
  advanceInternalStateOneTimeStep(ensemble, param, dt);
  advanceMotionalStateOneTimeStep(ensemble, param, dt);
}

void advanceInterval(Ensemble& ensemble, const Param& param, const double meanP, const double dt)
{
  addAtomsFromSource(ensemble, param, meanP);
  removeAtomsAtWalls(ensemble, param);
  advanceAtomsOneTimeStep(ensemble, param, dt);
}

void storeObservables(Observables& observables, int s, Ensemble& ensemble, 
    const Param& param)
{
  observables.nAtom(s) = ensemble.atoms.size();
  observables.intensity(s) = param.gammac*ensemble.internal.cov.sum();
  observables.inversion(s) = ensemble.internal.inv.sum()/ensemble.internal.inv.size();
  observables.spinSpinCor(s) = (ensemble.internal.cov.sum()-ensemble.internal.cov.trace())
                          /(ensemble.atoms.size()*(ensemble.atoms.size()-1));
}

void storeRawData(RawData& rawData, Ensemble& ensemble, const Param& param)
{
  rawData.spinSpinCovSS = ensemble.internal.cov;
  rawData.invSS = ensemble.internal.inv;
}

void evolve(Ensemble& ensemble, const Param& param, Observables& observables, RawData& rawData)
{
  //meanP
  double meanP = param.yWall*2/param.transitTime; //vy = deltay/tau

  //Integration conditions
  double dt = param.meanAtomGeneratingNumber/param.meanAtomNumber*param.transitTime; //dt = dN/N0*tau
  double tmax = param.steadyTime*param.transitTime;
  
  //evolve
  int nstep = tmax/dt+0.5;
  double tstep = dt, t=0;

  for (int n=0, s=0; n<=nstep; n++, t += tstep) {
    if ((long)(n+1)*param.nstore/(nstep+1) > s) {
      storeObservables(observables, s++, ensemble, param);
    }
    if (n == nstep-1) {
      storeRawData(rawData, ensemble, param);
    }
    if (n != nstep)
      advanceInterval(ensemble, param, meanP, dt);
  }
}

void writeObservables(ObservableFiles& observableFiles, 
    Observables& observables, RawData& rawData)
{
  observableFiles.nAtom << observables.nAtom << std::endl;
  observableFiles.intensity << observables.intensity << std::endl;
  observableFiles.inversion << observables.inversion << std::endl;
  observableFiles.spinSpinCor << observables.spinSpinCor << std::endl;
  observableFiles.spinSpinCovSS << rawData.spinSpinCovSS << std::endl;
  observableFiles.invSS << rawData.invSS << std::endl;
}

void mkdir(Param& param) {
  std::string mkdir = "mkdir "+param.name; //make a new directory to store data
  system(mkdir.c_str());
  std::string cpInput = "cp input.txt "+param.name;
  system(cpInput.c_str());  
  std::string moveparam = "mv *.dat "+param.name;
  system(moveparam.c_str());
};

int main(int argc, char *argv[])
{
  CmdLineArgs config;
  getOptions(argc, argv, &config);

  //Set up parameters
  Param param;
  getParam (config.configFile, &param);
	
  //Set up initial conditions
  Ensemble ensemble;
  Observables observables(param.nstore);
  RawData rawData;

  //Start simulation
  evolve(ensemble, param, observables, rawData);

  //Write Observables
  ObservableFiles observableFiles;
  writeObservables(observableFiles, observables, rawData);
  
  //Move .dat files into the directory named "name"
  mkdir(param);

  return 0;
}


//debug
//std::cout << ensemble.cov << std::endl << std::endl;
//debug