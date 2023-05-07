#ifndef FEMFUNC
#define FEMFUNC
#include "mesh/ActiveMesh.h"
#include "mesh/Mesh.h"
#include "mesh/Element.h"
#include <Eigen/Core>
#include <list>
#include <vector>

namespace fem
{
class Function{
  public:
    Eigen::VectorXd values;
    const mesh::ActiveMesh* domain;

    Function(const mesh::ActiveMesh* dom, const Eigen::VectorXd &otherValues = Eigen::VectorXd())
    {
      domain = dom;
      if (otherValues.size() == 0) {
        values = Eigen::VectorXd::Zero( domain->mesh->nnodes );
      } else {
        if (otherValues.size() != domain->mesh->nnodes) {
          cout << "Wrong vals mesh pair in Function construction!" << endl;
          exit(-1);
        }
        values = otherValues;
      }
    }

    double evaluate( Eigen::Vector3d point ) const;
    Eigen::Vector3d evaluateGrad( Eigen::Vector3d point );
    void interpolate(Function &extFEMFunc );
    void setValues( const Eigen::VectorXd &otherValues ) {
      values = otherValues;
    }
};
void interpolate( list<Function> &targetFunctions, const list<Function> &sourceFunctions );
}
#endif
