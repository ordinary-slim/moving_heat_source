#ifndef PROBLEM
#include <map>
#include <string>
#include "mesh.h"
#include "heatSource.h"
#include <Eigen/Core>

using namespace std;

class Problem {
  public:
    Mesh mesh;
    HeatSource mhs;
    Eigen::VectorXd solution;
    Eigen::VectorXd deltaSolution;
    map<string,float> material;
    string timeIntegration;
    double time = 0.0;
    double dt;

    void initialize(map<string,double> &input);
    void iterate();
    };
#define PROBLEM
#endif
/*
maxIter		100
L		10.0
radius		2.0
P		1000000.0
x0		5.0
speed		10
t		0.0
Tfinal		10.0
rho		4000
k		200
c_p		10
initialTemperature		25
dt		0.1
nels		100
iter		0
maxIter		100
plot    True
*/
