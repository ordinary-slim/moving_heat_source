#ifndef FEMFUNC
#define FEMFUNC
#include "mesh/Domain.h"
#include "mesh/Mesh.h"
#include "mesh/Element.h"
#include <Eigen/Core>
#include <vector>

class AbstractFunction {
  public:
    const mesh::Domain* domain;
    virtual double evaluate( Eigen::Vector3d &point ) const = 0;
};

namespace fem
{
class Function : public AbstractFunction {
  public:
    Eigen::VectorXd values;// value at each node
                           // mesh is accessed through the domain

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

    double evaluate( Eigen::Vector3d &point ) const;
    Eigen::Vector3d evaluateGrad( Eigen::Vector3d &point );
    void interpolate(const AbstractFunction &extFEMFunc,
        const mesh::MeshTag<int> &nodalTag,
        std::function<bool(int)> filter = nullptr,
        bool ignoreOutside=false );
    void interpolate(const AbstractFunction &extFEMFunc,
        const mesh::MeshTag<int> &nodalTag,
        bool ignoreOutside=false );
    void interpolate(const AbstractFunction &extFEMFunc, bool ignoreOutside=false );
    void interpolateInactive( const AbstractFunction &extFEMFunc, bool ignoreOutside );
    double getL2Norm() const;
    friend Function operator-(const Function& f1, const Function& f2);
};

Function interpolate( const AbstractFunction &extFEMFunc, const mesh::Domain *domain,
                           bool ignoreOutside = false );
}

class ConstantFunction : public AbstractFunction {
  public:
    double constant = 0.0;
    ConstantFunction( const mesh::Domain* dom, double c ) {
      domain = dom;
      constant = c;
    }
    double evaluate( Eigen::Vector3d &point ) const {
      //TODO: Check if point is in domain (?)
      return constant;
    }
};

#endif
