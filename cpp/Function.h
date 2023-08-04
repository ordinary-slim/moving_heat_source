#ifndef FEMFUNC
#define FEMFUNC
#include "mesh/Domain.h"
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
    const mesh::Domain* domain;

    Function(const mesh::Domain* dom) {
      domain = dom;
      values = Eigen::VectorXd::Zero( domain->mesh->nnodes );
    }
    Function(const mesh::Domain* dom, const Eigen::VectorXd &values) {
      domain = dom;
      if (values.size() != domain->mesh->nnodes) {
        throw std::invalid_argument("Provided values size is not nnodes.");
      }
      this->values = values;
    }
    Function(const mesh::Domain* dom, Eigen::VectorXd &values) {
      domain = dom;
      if (values.size() != domain->mesh->nnodes) {
        throw std::invalid_argument("Provided values size is not nnodes.");
      }
      this->values = move(values);
    }
    template<typename T>
    Function(const mesh::Domain* dom, const mesh::MeshTag<T> &tag) {
      if (tag.dim()!=0) {
        throw std::invalid_argument("Expected nodal MeshTag in Function constructor.");
      }
      domain = dom;
      Eigen::VectorXd convertedVals(domain->mesh->nnodes);
      // Convert to double
      for (int inode = 0; inode < domain->mesh->nnodes; ++inode) {
        convertedVals[inode] = double(tag[inode]);
      }
      values = move(convertedVals);
    }
    Function(const Function&) = default;

    double evaluate( Eigen::Vector3d point ) const;
    Eigen::Vector3d evaluateGrad( Eigen::Vector3d point );
    void interpolate(const Function &extFEMFunc );
    void interpolateInactive( const Function &extFEMFunc, bool ignoreOutside );
    void setValues( const Eigen::VectorXd &values ) { this->values = values; }
    double getL2Norm() const;
    friend Function operator-(const Function& f1, const Function& f2);
};
void interpolate( list<Function> &targetFunctions, const list<Function> &sourceFunctions );

fem::Function interpolate( const fem::Function &extFEMFunc,
                           const mesh::Domain *domain,
                           bool ignoreOutside = false );
}
#endif
