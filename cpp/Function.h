#ifndef FEMFUNC
#define FEMFUNC
#include "mesh/Domain.h"
#include "mesh/Mesh.h"
#include "mesh/Element.h"
#include <Eigen/Core>
#include <vector>

class AbstractFunction {
  public:
    const mesh::Domain* domain;//This line should stay at the top
    Eigen::VectorXd values;
    virtual double evaluate( Eigen::Vector3d &point ) const = 0;
    AbstractFunction(size_t ndofs, const mesh::Domain* dom) {
      _ndofs = ndofs;
      domain = dom;
      values = Eigen::VectorXd::Zero( _ndofs );
    }
    AbstractFunction(size_t ndofs, const mesh::Domain* dom, const Eigen::VectorXd &values) {
      _ndofs = ndofs;
      domain = dom;
      if (values.size() != _ndofs ) {
        throw std::invalid_argument("Provided values size is not compatible.");
      }
      this->values = values;
    }
    AbstractFunction(size_t ndofs, const mesh::Domain* dom, Eigen::VectorXd &values) {
      domain = dom;
      _ndofs = ndofs;
      if (values.size() != _ndofs ) {
        throw std::invalid_argument("Provided values size is not compatible.");
      }
      this->values = move(values);
    }
    template<typename T>
    AbstractFunction(size_t ndofs, const mesh::Domain* dom, const mesh::MeshTag<T> &tag) {
      domain = dom;
      _ndofs = ndofs;
      if (tag.size()!=_ndofs) {
        throw std::invalid_argument("Tag is not compatible with Function.");
      }
      Eigen::VectorXd convertedVals(_ndofs);
      // Convert to double
      for (int idof = 0; idof < _ndofs; ++idof) {
        convertedVals[idof] = double(tag[idof]);
      }
      values = move(convertedVals);
    }

  protected:
    size_t _ndofs;
};

namespace fem
{
class Function : public AbstractFunction {
  /*
   * FEM function
   */
  public:

    Function(const mesh::Domain* dom) : AbstractFunction( dom->mesh->nnodes, dom ) {}
    Function(const mesh::Domain* dom, const Eigen::VectorXd &values) : AbstractFunction( dom->mesh->nnodes, dom, values) {}
    Function(const mesh::Domain* dom, Eigen::VectorXd &values) : AbstractFunction( dom->mesh->nnodes, dom, values) {}
    template<typename T>
    Function(const mesh::Domain* dom, const mesh::MeshTag<T> &tag) : AbstractFunction( dom->mesh->nnodes, dom, tag) {}

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
    friend Function operator/(const Function& f, const double c);
};

Function interpolate( const AbstractFunction &extFEMFunc, const mesh::Domain *domain,
                           bool ignoreOutside = false );

class DG0Function : public AbstractFunction {
  public:
    DG0Function(const mesh::Domain* dom) : AbstractFunction( dom->mesh->nels, dom ) {}
    DG0Function(const mesh::Domain* dom, const Eigen::VectorXd &values) : AbstractFunction( dom->mesh->nels, dom, values) {}
    DG0Function(const mesh::Domain* dom, Eigen::VectorXd &values) : AbstractFunction( dom->mesh->nels, dom, values) {}
    template<typename T>
    DG0Function(const mesh::Domain* dom, const mesh::MeshTag<T> &tag) : AbstractFunction( dom->mesh->nels, dom, tag) {}

    double evaluate( Eigen::Vector3d &point ) const;
    void interpolate(const AbstractFunction &extFEMFunc,
        const mesh::MeshTag<int> &cellTag,
        std::function<bool(int)> filter = nullptr,
        bool ignoreOutside=false );
};

}

class ConstantFunction : public AbstractFunction {
  public:
    ConstantFunction( const mesh::Domain* dom, double c ) :
        AbstractFunction(1, dom) {
          values(0) = c;
    }
    double evaluate( Eigen::Vector3d &point ) const {
      //TODO: Check if point is in domain (?)
      return values(0);
    }
};

#endif
