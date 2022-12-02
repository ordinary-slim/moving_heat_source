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
    map<string,float> material;
    Eigen::VectorXd solution;
    Eigen::MatrixXd prevSolutions;
    double time = 0.0;
    double dt;
    int iter;

    int currentIntegrator = 1;
    int desiredIntegrator = 1;
    int nstepsRequired = 1;
    int nstepsStored   = 0;

    void initialize(map<string,double> &input);
    void iterate();
    };
#define PROBLEM
#endif
