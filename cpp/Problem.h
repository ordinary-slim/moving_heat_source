#ifndef PROBLEM
#include <map>
#include <string>
#include "mesh/Mesh.h"
#include "FEMFunction.h"
#include "HeatSource.h"
#include "timeIntegrator.h"
#include <Eigen/Core>
#include <Eigen/Sparse>
#include "../external/pybind11/include/pybind11/pybind11.h"
#include "../external/pybind11/include/pybind11/stl.h"

namespace py = pybind11;

class Problem {
  public:
    mesh::Mesh mesh;
    FEMFunction unknown;
    HeatSource mhs;
    map<string,double> material;
    double time = 0.0;
    double dt = 0.0;
    int iter;

    bool isAdvection = false;
    bool isSteady    = false;
    bool isStabilized = false;
    Eigen::Vector3d advectionSpeed;
    Eigen::VectorXd rhs;
    Eigen::SparseMatrix<double> lhs;
    vector<Eigen::Triplet<double>> lhsCoeffs;
    Eigen::SparseMatrix<double> M; // mass mat
    vector<Eigen::Triplet<double>> massCoeffs;

    // Neumann BC
    std::vector<int> neumannFacets;
    std::vector<double> neumannFluxes;

    // integrator
    TimeIntegratorHandler timeIntegrator;

    void setTime(double newTime) {
      time = newTime;
      mhs.time = newTime;
    }
    void setPointers(){
      unknown.mesh = &mesh;
    }

    void setAdvectionSpeed(Eigen::Vector3d inputAdvectionSpeed){
      advectionSpeed = inputAdvectionSpeed;
      if (advectionSpeed.norm() > 1e-10) isAdvection = true;
    }
    void initialize(py::dict &input);
    void initializeIntegrator(Eigen::MatrixXd pSols);
    void iterate();
    void updateFRFpos();
    void assembleSpatialLHS();
    void assembleSpatialRHS();
    void assembleStabilization();//Only P1/Q1 for the moment!
    void assembleTime();
    void assembleNeumann();
    void forceDirichletNodes();
    void forceInactiveNodes();
    void preIterate();
    void postIterate();
    void activateDomain(vector<int> inputActiveElements ) {
      mesh.setActiveElements( inputActiveElements );
      //isAssembled = false;
    }
    void setStabilization(bool stabilize) {
      isStabilized = stabilize;
    }
    void setNeumann( vector<vector<int>> otherNeumannFacets, double neumannFlux );
    };
#define PROBLEM
#endif
