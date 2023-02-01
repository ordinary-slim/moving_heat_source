#ifndef PROBLEM
#include <map>
#include <string>
#include "mesh.h"
#include "heatSource.h"
#include "timeIntegrator.h"
#include <Eigen/Core>
#include <Eigen/Sparse>
#include "../../external/pybind11/include/pybind11/pybind11.h"
#include "../../external/pybind11/include/pybind11/stl.h"

namespace py = pybind11;

using namespace std;

typedef Eigen::SparseMatrix<double> SpMat; // declares a column-major sparse matrix type of double
                                           //
class Problem {
  public:
    Mesh mesh;
    HeatSource mhs;
    map<string,double> material;
    Eigen::VectorXd solution;
    Eigen::MatrixXd prevSolutions;
    double time = 0.0;
    double dt;
    int iter;

    bool isAssembled = false;
    bool isAdvection = false;
    bool isSteady    = false;
    Eigen::Vector3d advectionSpeed;
    Eigen::VectorXd rhs;
    SpMat lhs;
    SpMat M; // mass mat
    SpMat K; // stiffness mat
    SpMat A; // advection mat
    SpMat I; // inactive nodes
    Eigen::VectorXd pulse; // source term


    // integrator
    TimeIntegratorHandler timeIntegrator;

    void setTime(double newTime) {
      time = newTime;
      mhs.time = newTime;
    }
    void initialize(py::dict &input);
    void initializeIntegrator(Eigen::MatrixXd pSols);
    void iterate();
    void postIterate();
    void activateDomain(vector<int> inputActiveElements ) {
      mesh.setActiveElements( inputActiveElements );
      isAssembled = false;
    }
    };
#define PROBLEM
#endif
