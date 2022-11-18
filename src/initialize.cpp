#include "problem.h"
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
  mhs.speed[0] = input["speedX"];
  mhs.speed[1] = input["speedY"];
  mhs.speed[2] = input["speedZ"];
  mhs.initialPosition[0] = input["initialPositionX"];
  mhs.initialPosition[1] = input["initialPositionY"];
  mhs.initialPosition[2] = input["initialPositionZ"];

  // initialize solution and increment
  double environmentTemperature = input["environmentTemperature"];
  solution = Eigen::VectorXd::Constant( mesh.nnodes, environmentTemperature );
  deltaSolution = Eigen::VectorXd::Zero( mesh.nnodes );
  // gaussian IC
  //mhs.computePulse( solution, time, mesh );

  // dirichlet BC
  vector<int> freeNodes(mesh.nnodes);


  // material. dictionnary is not efficient + involved in assembly
  material["rho"] = input["rho"];
  material["k"] = input["conductivity"];
  material["cp"] = input["specific_heat"];

  if (input["timeIntegration"] == 0) {
    timeIntegration = "ForwardEuler";
  } else if (input["timeIntegration"] == 1) {
    timeIntegration = "BackwardEuler";
  }
  // initial condition
}
