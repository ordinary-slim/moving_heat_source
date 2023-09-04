#ifndef PROBLEM
#define PROBLEM
#include <map>
#include <string>
#include <list>
#include <memory>
#include "Material.h"
#include "mesh/Domain.h"
#include "Function.h"
#include "heatSource/HeatSource.h"
#include "heatSource/LumpedHeatSource.h"
#include "timeIntegrator.h"
#include "linearAlgebra/LinearSystem.h"
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>

namespace py = pybind11;

class Problem {
  public:
    Problem(mesh::Mesh &mesh, py::dict &input);

    mesh::Domain domain;
    fem::Function unknown;
    list<fem::Function> previousValues;
    std::unique_ptr<heat::HeatSource> mhs;//run time polymporphism
    bool hasPreIterated = false;
    ThermalMaterial material;
    double time = 0.0;
    double dt = 0.0;

    double Tdeposition, Tenv;

    bool isAdvection = false;
    bool isSteady    = false;
    bool isStabilized = false;
    Eigen::Vector3d advectionSpeed;

    std::shared_ptr<LinearSystem> ls = NULL;
    bool assembling2external = false;
    vector<int> dofNumbering;
    vector<int> freeDofsNumbering;// Necessary to assemble Gamma Dirichlet nodes
    mesh::MeshTag<int>    forcedDofs;

    // Dirichlet BC
    mesh::MeshTag<int>    dirichletNodes;// 1 is Dirichlet
                                         // 2 is Dirichlet @ interface
    mesh::MeshTag<double> dirichletValues;

    // Neumann BC
    mesh::MeshTag<int>    weakBcFacets; // 1 is Neumann
                                        // 2 is convection
                                        // 3 is Neumann @ interface
    mesh::MeshTag<std::vector<double>> neumannFluxes;//[ifacet][igpoint]

    // Coupling BC
    bool isCoupled = false;
    mesh::MeshTag<int> gammaNodes;
    mesh::MeshTag<int> gammaFacets;

    // integrator
    TimeIntegratorHandler timeIntegrator;

    void setDt( double dt ) {
      /* Should only be used before first time-step
       * Not ready yet for change of DeltaT between time iterations */
      this->dt = dt;
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
    void updateInterface( const Problem &pExt );
    void updateInterface( mesh::MeshTag<int> &activeInExternal );
    void assemble(const Problem* externalProblem = nullptr);
    void gather();
    void assembleSpatialPDE();//mass, diffusion, advection
    void assembleWeakBcs();
    void assembleTime();
    void assembleDirichletGamma( const Problem *pExt ); 
    void assembleNeumannGamma( const Problem *pExt ); 
    void updateForcedDofs();
    void preAssemble(bool allocateLs=true);
    void preIterate(bool canPreassemble=true);
    void postIterate();
    void setStabilization(bool stabilize) {
      isStabilized = stabilize;
    }
    void setNeumann( vector<vector<unsigned int>> otherNeumannNodes, double neumannFlux );
    void setNeumann( Eigen::Vector3d pointInPlane, Eigen::Vector3d normal, double neumannFlux );
    void setNeumann( vector<int> otherNeumannFacets, std::function<Eigen::Vector3d(Eigen::Vector3d)> fluxFunc );
    void setConvection(bool resetBcs = false);
    void setDirichlet( vector<int> otherDirichletFacets, std::function<double(Eigen::Vector3d)> dirichletFunc );
    void setDirichlet( const vector<int> &otherDirichletNodes, const vector<double> &otherDirichletValues );
    void setGamma2Neumann();
    void setGamma2Dirichlet();
    mesh::MeshTag<int> getActiveInExternal( const Problem &pExt, double tol=1e-5 );
    void uniteExternal( const Problem &pExt, bool updateGamma = true);
    void substractExternal( const Problem &pExt, bool updateGamma = true);
    void intersectExternal( const Problem &pExt, bool updateGamma = true, bool interpolateIntersected = false );
    void interpolate2dirichlet( fem::Function &extFEMFunc);

    void clearBCs() {
      dirichletNodes.setCteValue( 0 );
      dirichletValues.setCteValue( 0 );
      weakBcFacets.setCteValue( 0 );
      neumannFluxes.setCteValue( vector<double>() );
    }
    void clearGamma() {
      /* Clear boundary conditions at interfacet
       * and interface itself
       */
      // Clear Dirichlet BCs at Gamma
      weakBcFacets.tag( [](int tag){ return tag==3; }, 0 );
      dirichletNodes.tag( [](int tag){ return tag==2; }, 0 );
      // Clear Gamma
      gammaNodes.setCteValue( 0 );
      gammaFacets.setCteValue( 0 );

      isCoupled = false;
    }
    fem::Function project( std::function<double(Eigen::Vector3d)> func );//L2 projection onto domain attribute
};

Eigen::VectorXd CreateEigenVector(py::array_t<double> n); // Helper function initialization. Where should I put this?
#endif
