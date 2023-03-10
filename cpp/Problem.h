#ifndef PROBLEM
#include <map>
#include <string>
#include "mesh/Mesh.h"
#include "FEMFunction.h"
#include "heatSource.h"
#include "timeIntegrator.h"
#include <Eigen/Core>
#include <Eigen/Sparse>
#include "../external/pybind11/include/pybind11/pybind11.h"
#include "../external/pybind11/include/pybind11/stl.h"

namespace py = pybind11;

typedef Eigen::SparseMatrix<double> SpMat; // declares a column-major sparse matrix type of double
                                           //
class Problem {
  public:
    mesh::Mesh mesh;
    FEMFunction unknown;
    HeatSource mhs;
    map<string,double> material;
    double time = 0.0;
    double dt;
    int iter;

    bool isAdvection = false;
    bool isSteady    = false;
    bool isStabilized = false;
    double SCA = 512;//stabilization cte advection
    Eigen::Vector3d advectionSpeed;
    Eigen::VectorXd rhs;
    SpMat lhs;
    SpMat M; // mass mat

    // integrator
    TimeIntegratorHandler timeIntegrator;

    void setTime(double newTime) {
      time = newTime;
      mhs.time = newTime;
    }
    void setAdvectionSpeed(Eigen::Vector3d inputAdvectionSpeed){
      advectionSpeed = inputAdvectionSpeed;
    }
    void initialize(py::dict &input);
    void initializeIntegrator(Eigen::MatrixXd pSols);
    void iterate();
    void updateFRFpos();
    void assembleSpatialLHS();
    void assembleSpatialRHS();
    void assembleStabilization();//Only P1/Q1 for the moment!
    void assembleTime();
    void forceInactiveNodes();
    void preIterate();
    void postIterate();
    void activateDomain(vector<int> inputActiveElements ) {
      mesh.setActiveElements( inputActiveElements );
      //isAssembled = false;
    }
    };
#define PROBLEM
#endif
