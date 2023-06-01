#ifndef PROBLEM
#define PROBLEM
#include <map>
#include <string>
#include <list>
#include "mesh/ActiveMesh.h"
#include "Function.h"
#include "HeatSource.h"
#include "timeIntegrator.h"
#include "linearAlgebra/LinearSystem.h"
#include "../external/pybind11/include/pybind11/pybind11.h"
#include "../external/pybind11/include/pybind11/stl.h"

namespace py = pybind11;

class Problem {
  public:
    Problem(mesh::Mesh &mesh, py::dict &input);
    mesh::ActiveMesh domain;
    fem::Function unknown;
    list<fem::Function> previousValues;
    HeatSource mhs;
    double density = 1.0, conductivity = 1.0, specificHeat = 1.0, convectionCoeff = 0.0;
    double time = 0.0;
    double dt = 0.0;
    int iter;

    bool isAdvection = false;
    bool isSteady    = false;
    bool isStabilized = false;
    bool isConvection = false;
    Eigen::Vector3d advectionSpeed;

    LinearSystem myls;
    LinearSystem* ls = NULL;
    bool assembling2external = false;
    vector<int> dofNumbering;
    mesh::MeshTag<int>    forcedDofs;

    // Dirichlet BC
    mesh::MeshTag<int>    dirichletNodes;
    mesh::MeshTag<double> dirichletValues;

    // Neumann BC
    mesh::MeshTag<int>    neumannFacets;
    mesh::MeshTag<std::vector<double>> neumannFluxes;//[ifacet][igpoint]

    // Convection BC
    mesh::MeshTag<int> convectionFacets;
    double Tenv;

    // Coupling BC
    mesh::MeshTag<int> gammaNodes;
    mesh::MeshTag<int> gammaFacets;

    // integrator
    TimeIntegratorHandler timeIntegrator;

    void setTime(double newTime) {
      time = newTime;
      mhs.time = newTime;
    }
    void setDeltaT( double newDeltaT ) {
      /*
       * Should only be used before first time-step
       * Not ready yet for change of DeltaT between time iterations
       */
      dt = newDeltaT;
    }
    void setPointers(){
      unknown.domain = &domain;
    }

    void setAdvectionSpeed(Eigen::Vector3d inputAdvectionSpeed){
      advectionSpeed = inputAdvectionSpeed;
      if (advectionSpeed.norm() > 1e-10) isAdvection = true;
    }
    void initializeIntegrator(Eigen::MatrixXd pSols);
    void iterate();
    void findGamma( const Problem &pExt );
    void findGamma( mesh::MeshTag<int> &activeInExternal );
    void assemble();
    void gather();
    void updateFRFpos();
    void assembleSpatialPDE();//mass, diffusion, advection
    void assembleWeakBcs();
    void assembleTime();
    void updateForcedDofs();
    void preAssemble();
    void preIterate();
    void postIterate();
    void setStabilization(bool stabilize) {
      isStabilized = stabilize;
    }
    void setNeumann( vector<vector<int>> otherNeumannFacets, double neumannFlux );
    void setNeumann( Eigen::Vector3d pointInPlane, Eigen::Vector3d normal, double neumannFlux );
    void setNeumann( vector<int> otherNeumannFacets, std::function<Eigen::Vector3d(Eigen::Vector3d)> fluxFunc );
    void setDirichlet( vector<int> otherDirichletFacets, std::function<double(Eigen::Vector3d)> dirichletFunc );
    void setDirichlet( const vector<int> &otherDirichletNodes, const vector<double> &otherDirichletValues );
    mesh::MeshTag<int> getActiveInExternal( const Problem &pExt, double tol=1e-7 );
    void substractExternal( const Problem &pExt, bool updateGamma = true);
    void intersectExternal( const Problem &pExt, bool updateGamma = true );
    void interpolate2dirichlet( fem::Function &extFEMFunc);

    void clearBCs() {
      dirichletNodes.setCteValue( 0 );
      dirichletValues.setCteValue( 0 );
      neumannFacets.setCteValue( 0 );
      neumannFluxes.setCteValue( vector<double>() );
      convectionFacets.setCteValue( 0 );
    }
    fem::Function project( std::function<double(Eigen::Vector3d)> func );//L2 projection onto domain attribute
};
#endif
