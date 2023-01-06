#include "includes/problem.h"
#include <map>
#include <string>

void Problem::initialize(map<string,double> &input) {
  // tstepping
  dt = input["dt"];
  // discrete mesh
  float a = input["Left"];
  float b = input["Right"];
  int nels = int( input["nels"] );
  mesh.initialize1DMesh( a, b, nels );

  // heat source
  mhs.radius = input["radius"];
  mhs.power = input["power"];
  if (input.count("efficiency")==1) mhs.efficiency = input["efficiency"];

  mhs.speed[0] = input["speedX"];
  mhs.speed[1] = input["speedY"];
  mhs.speed[2] = input["speedZ"];
  mhs.initialPosition[0] = input["initialPositionX"];
  mhs.initialPosition[1] = input["initialPositionY"];
  mhs.initialPosition[2] = input["initialPositionZ"];
  // set type of source term
  switch (int(input["sourceTerm"])) {
    case 91:
      { mhs.powerDensity = &forcedSolutionSource91;
        break; }
    default:
      { mhs.powerDensity = &gaussianPowerDensity;
        break; }
  }

  // initialize solution and increment
  double environmentTemperature = input["environmentTemperature"];
  solution = Eigen::VectorXd::Constant( mesh.nnodes, environmentTemperature );
  // gaussian IC
  //mhs.computePulse( solution, time, mesh );

  // dirichlet BC
  vector<int> freeNodes(mesh.nnodes);


  // material. dictionnary is not efficient + involved in assembly
  material["rho"] = input["rho"];
  material["k"] = input["conductivity"];
  material["cp"] = input["specific_heat"];

  // check for advection term
  if (input.count("isAdvection")==1) {
    isAdvection = (input["isAdvection"]==1);
    if (isAdvection) {
      advectionSpeed[0] = input["advectionSpeedX"];
      advectionSpeed[1] = input["advectionSpeedY"];
      advectionSpeed[2] = input["advectionSpeedZ"];
      cout << "advectionSpeed= " << advectionSpeed << endl;
    }
  }

  // check for time dependency
  if (input.count("steadyState")==1) {
    isSteady = (input["steadyState"]==1);
  }

  // timeIntegrator
  timeIntegrator.setRequiredSteps( input["timeIntegration"] );

  // allocate storage for previous solutions
  prevSolutions = Eigen::MatrixXd::Zero( mesh.nnodes, timeIntegrator.nstepsRequired );
  prevSolutions.col( 0 ) << solution;
  ++timeIntegrator.nstepsStored;
}
