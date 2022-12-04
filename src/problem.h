#ifndef PROBLEM
#include <map>
#include <string>
#include "mesh.h"
#include "heatSource.h"
#include "timeIntegrator.h"
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

    // integrator
    TimeIntegratorHandler timeIntegrator;

    void initialize(map<string,double> &input);
    void initializeIntegrator(Eigen::MatrixXd pSols);
    void iterate();
    };
#define PROBLEM
#endif
