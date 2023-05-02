#ifndef FEMFUNC
#include "mesh/Mesh.h"
#include "mesh/Element.h"
#include <Eigen/Core>
#include <list>
#include <vector>

namespace fem
{
class Function{
  public:
    mesh::Mesh* mesh;
    Eigen::VectorXd values;
    std::vector<int> dirichletNodes;
    std::vector<double> dirichletValues;

    Function(){
    }

    Function(mesh::Mesh &otherMesh, const Eigen::VectorXd &otherValues = Eigen::VectorXd()){
      mesh = &otherMesh;
      if (otherValues.size() == 0) {
        values = Eigen::VectorXd::Zero( mesh->nnodes );
      } else {
        if (otherValues.size() != mesh->nnodes) {
          cout << "Wrong vals mesh pair in Function construction!" << endl;
          exit(-1);
        }
        values = otherValues;
      }
    }

    double evaluate( Eigen::Vector3d point ) const;
    Eigen::Vector3d evalGrad( Eigen::Vector3d point );
    void interpolate(Function &extFEMFunc );
    void interpolate2dirichlet(Function &extFEMFunc);
    void releaseDirichlet(){
      dirichletNodes.clear();
      dirichletValues.clear();
    }
    void setValues( const Eigen::VectorXd &otherValues ) {
      values = otherValues;
    }
};
void interpolate( list<Function> &targetFunctions, const list<Function> &sourceFunctions );
}
#define FEMFUNC
#endif
