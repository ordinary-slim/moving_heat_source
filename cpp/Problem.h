#ifndef PROBLEM
#include <map>
#include <string>
#include <list>
#include "mesh/Mesh.h"
#include "Function.h"
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
    fem::Function unknown;
    list<fem::Function> previousValues;
    HeatSource mhs;
    map<string,double> material;
    double time = 0.0;
    double dt = 0.0;
    int iter;

    bool isAdvection = false;
    bool isSteady    = false;
    bool isStabilized = false;
    bool isConvection = false;
    Eigen::Vector3d advectionSpeed;
    Eigen::VectorXd rhs;
    Eigen::SparseMatrix<double> lhs;
    vector<Eigen::Triplet<double>> lhsCoeffs;
    Eigen::SparseMatrix<double> M; // mass mat
    vector<Eigen::Triplet<double>> massCoeffs;

    // Neumann BC
    std::vector<int> neumannFacets;
    std::vector<double> neumannFluxes;

    // Convection BC
    double Tenv;
    std::vector<int> convectionFacets;

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
    void assembleSpatialLHS();//mass, diffusion, advection
    void assembleSpatialRHS();//source term
    void assembleConvectionLHS();
    void assembleConvectionRHS();
    void assembleStabilization();//Only P1/Q1 for the moment!
    void assembleTime();
    void assembleNeumann();
    void forceDirichletNodes();
    void forceInactiveNodes();
    void preIterate();
    void postIterate();
    void setActiveElements(vector<int> otherActiveElements ) {
      mesh.setActiveElements( otherActiveElements );
      //POS ACTIVATION OPS
    }
    void setActiveNodes(vector<int> otherActiveNodes ) {
      mesh.setActiveNodes( otherActiveNodes );
      //POS ACTIVATION OPS
    }
    void setStabilization(bool stabilize) {
      isStabilized = stabilize;
    }
    void setNeumann( vector<vector<int>> otherNeumannFacets, double neumannFlux );
    void setNeumann( Eigen::Vector3d pointInPlane, Eigen::Vector3d normal, double neumannFlux );
    };
#define PROBLEM
#endif
