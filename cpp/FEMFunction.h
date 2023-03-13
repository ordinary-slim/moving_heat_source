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

    double evaluateVal( Eigen::Vector3d point );
    vector<double> evaluateValNPrevVals( Eigen::Vector3d point );
    void getFromExternal( FEMFunction &extFEMFunc );
};
#define FEMFUNC
#endif
