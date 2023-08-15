#include <list>
#include "Function.h"

namespace fem
{
double Function::evaluate( Eigen::Vector3d point ) const {
  /*
  Output val of Function at input point
  */
  double val = 0;

  // GET VALS OF SHAPE FUNCS AT POINT
  int idxOwnerEl = domain->findOwnerElements( point );
  mesh::Element e = domain->mesh->getElement( idxOwnerEl );//Load element containing point
  Eigen::VectorXd shaFunVals = e.evaluateShaFuns( point );

  val = values( *e.con ).dot( shaFunVals );

  return val;
}

Eigen::Vector3d Function::evaluateGrad( Eigen::Vector3d point ) {
  /*
  Output gradient of Function @ input point
  */
  Eigen::Vector3d grad = Eigen::Vector3d::Zero();

  // GET VALS OF GRAD o SHAPE FUNCS AT POINT
  int idxOwnerEl = domain->findOwnerElements( point );
  mesh::Element e = domain->mesh->getElement( idxOwnerEl );//Load element containing point
  Eigen::MatrixXd gradShaFunVals = e.evaluateGradShaFuns( point );

  for (int inode = 0; inode < e.nnodes; ++inode) {
    grad += gradShaFunVals.row(inode) * values[ (*e.con)[ inode ] ] ;
  }

  return grad;
}


void Function::interpolate(const Function &extFEMFunc, const mesh::MeshTag<int> &nodalTag, bool ignoreOutside ) {

  if (nodalTag.dim() != 0) {
    throw std::invalid_argument("Interpolate requires nodal mesh tag.");
  }
  vector<int> nodesOfInterest = nodalTag.getIndices();

  for (int inode : nodesOfInterest) {
    // Move to reference frame of external
    Eigen::Vector3d posExt = domain->posLab.row(inode) - extFEMFunc.domain->translationLab.transpose();
    try {
      values[inode] = extFEMFunc.evaluate( posExt );
    } catch ( const std::invalid_argument &e ) {
      // Point outside of domain
      if (domain->activeNodes[inode] && not(ignoreOutside)) {
        throw;
      } else {
        values[inode] = -1;// sentinel value
      }
    }
  }
}

void Function::interpolate(const Function &extFEMFunc, bool ignoreOutside ) {
  mesh::MeshTag<int> allNodes = mesh::MeshTag<int>( domain->mesh, 0, 1 );
  interpolate( extFEMFunc, allNodes, ignoreOutside );
}

void Function::interpolateInactive( const Function &extFEMFunc, bool ignoreOutside ) {
  // Shorthand
  interpolate( extFEMFunc, not( domain->activeNodes ), ignoreOutside );
}

double Function::getL2Norm() const {
  // TODO: cleanup when mass mat refactored
  // Workaround mass matrix not being accurate
  // Works for closed integration...
  double l2norm = 0.0;
  Eigen::VectorXd tmp =  domain->massMat * values;
  for (int inode : domain->activeNodes.getIndices() ) {
    l2norm += values[inode] * tmp[inode];
  }
  return sqrt( l2norm );
}

fem::Function interpolate( const fem::Function &extFEMFunc,
                           const mesh::Domain *domain,
                           bool ignoreOutside ) {
  fem::Function f = fem::Function( domain );
  f.interpolate( extFEMFunc, ignoreOutside );
  return f;
}

Function operator-(const Function& f1, const Function& f2) {
  if (f1.domain != f2.domain) {
    throw std::invalid_argument("Functions don't share domain.");
  }
  Eigen::VectorXd values = f1.values - f2.values;
  return Function( f1.domain, values );
}

}
