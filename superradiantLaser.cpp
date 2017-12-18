//This program is used to simulate the superradiant laser using the cumulant expansion method.
#include "superradiantLaser.hpp"

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
    if (dummy.compare("dt") == 0)
      configInput >> param->dt;
    else if (dummy.compare("tmax") == 0)
      configInput >> param->tmax;
    else if (dummy.compare("nstore") == 0)
      configInput >> param->nstore;
    else if (dummy.compare("nAtom") == 0)
      configInput >> param->nAtom;
    else if (dummy.compare("gammac") == 0)
      configInput >> param->gammac;
    else if (dummy.compare("repumping") == 0)
      configInput >> param->repumping;
    else if (dummy.compare("name") == 0)
      configInput >> param->name;
    else {
      std::cout << "Error: invalid label " << dummy << " in "
          << filename << std::endl;
      exit(-1);
    }
  }
}

MatrixXd RHS(const MatrixXd& cov, const Param& param)
{
  MatrixXd RHS(cov);
  double gc = param.gammac;
  double w = param.repumping;
  int nAtom = param.nAtom;

  //right hand side of the DE set
   for (int i = 0; i < nAtom; i++) {
    //diagonal
    RHS(i,i) = 2*w*(1-cov(i,i))-0.5*gc*cov.colwise().sum()(i)-0.5*gc*cov.rowwise().sum()(i);
    //off-diagonal; we used the fact that RHS is symmetric
    for (int j = i+1; j < nAtom; j++) {
      RHS(i,j) = -(2*w+gc*cov(i,i)+gc*cov(j,j))*cov(i,j)
              +gc/2.0*(2*cov(i,i)-1)*cov.colwise().sum()(j)
              +gc/2.0*(2*cov(j,j)-1)*cov.rowwise().sum()(i);
      RHS(j,i) = RHS(i,j);
    }
  }
  return RHS;
}

void advanceInterval(MatrixXd& cov, const Param& param)
{
  double gc = param.gammac;
  double dt = param.dt;
  int nAtom = param.nAtom;
  
  //Define a new covariance matrix with initial value cov;
  MatrixXd newCov(cov);
  //Using RK2 method.
  newCov += dt/2.0*RHS(cov, param);
  //Second round
  cov += dt*RHS(newCov, param);
}

void storeObservables(Observables& observables, int s, const MatrixXd& cov, 
    const Param& param)
{
  observables.intensity(s) = param.gammac*cov.sum();//
  observables.intensityUnCor(s) = param.gammac*cov.diagonal().sum();
  observables.inversion(s) = (2*cov.diagonal().sum()-param.nAtom)/param.nAtom;
  observables.spinSpinCor(s) = (cov.sum()-cov.diagonal().sum())/(param.nAtom*(param.nAtom-1));
}

void evolve(MatrixXd& cov, const Param& param, Observables& observables)
{
  //Integration conditions
  double dt = param.dt;
  double tmax = param.tmax;
  
  //evolve
  int nTimeStep = tmax/dt+0.5;
  double tstep = dt, t=0;

  for (int n = 0, s = 0; n <= nTimeStep; n++, t += tstep) {
    if ((long)(n+1)*param.nstore/(nTimeStep+1) > s) {
      storeObservables(observables, s++, cov, param);
      //debug
      std::cout << "This is timestep " << n << "/" << nTimeStep << std::endl << std::endl;
      //debug
    }
    if (n != nTimeStep)
      advanceInterval(cov, param);
  }
}

void writeObservables(ObservableFiles& observableFiles, 
    Observables& observables)
{
  observableFiles.intensity << observables.intensity << std::endl;
  observableFiles.intensityUnCor << observables.intensityUnCor << std::endl;
  observableFiles.inversion << observables.inversion << std::endl;
  observableFiles.spinSpinCor << observables.spinSpinCor << std::endl;
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
  MatrixXd cov = MatrixXd::Identity(param.nAtom, param.nAtom);
  Observables observables(param.nstore);

  //Start simulation
  evolve(cov, param, observables);

  //Write Observables
  ObservableFiles observableFiles;
  writeObservables(observableFiles, observables);
  
  //Move .dat files into the directory named "name"
  mkdir(param);

  return 0;
}