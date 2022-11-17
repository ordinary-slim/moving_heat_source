#include "problem.h"
#include <map>
#include <string>

void Problem::initialize(map<string,float> &input) {
  // tstepping
  dt = input["dt"];
  // discrete mesh
  float a = input["Left"];
  float b = input["Right"];
  int nels = int( input["nels"] );
  mesh.initialize1DMesh( a, b, nels );
  fineMesh.initialize1DMesh( a, b, 100);

  // material. dictionnary is not efficient + involved in assembly
  material["rho"] = input["rho"];
  material["k"] = input["conductivity"];
  material["cp"] = input["specific_heat"];

  // initialize solution and increment
  double initialTemperature = input["initialTemperature"];
  solution = Eigen::VectorXd::Constant( mesh.nnodes, initialTemperature );
  deltaSolution = Eigen::VectorXd::Zero( mesh.nnodes );
}
