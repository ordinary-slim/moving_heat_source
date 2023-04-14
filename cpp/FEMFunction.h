#ifndef FEMFUNC
#include "mesh/Mesh.h"
#include "mesh/Element.h"
#include <Eigen/Core>
#include <vector>

class FEMFunction{
  public:
    mesh::Mesh* mesh;
    Eigen::VectorXd values;
    Eigen::MatrixXd prevValues;//col 0 is previous it, col 1 is prev prev it, etc
    std::vector<int> dirichletNodes;
    std::vector<double> dirichletValues;
    int nStepsRequired = 0;

    FEMFunction(){
    }

    FEMFunction(mesh::Mesh &otherMesh, int otherNStepsRequired){
      mesh = &otherMesh;
      nStepsRequired = otherNStepsRequired;
      values = Eigen::VectorXd::Zero( mesh->nnodes );
      prevValues = Eigen::MatrixXd::Zero( mesh->nnodes, nStepsRequired );
    }

    double evalVal( Eigen::Vector3d point );
    Eigen::Vector3d evalGrad( Eigen::Vector3d point );
    vector<double> evalValNPrevVals( Eigen::Vector3d point );
    void getFromExternal(FEMFunction &extFEMFunc );
    void forceFromExternal(FEMFunction &extFEMFunc);
    void releaseDirichlet(){
      dirichletNodes.clear();
      dirichletValues.clear();
    }
};
#define FEMFUNC
#endif
