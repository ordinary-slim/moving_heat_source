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


void Function::interpolate( const Function &extFEMFunc, bool ignoreOutside ) {

  for (int inode = 0; inode < domain->mesh->nnodes; inode++) {
    // Move to reference frame of external
    Eigen::Vector3d posExt = domain->mesh->pos.row(inode) + (domain->translationLab - extFEMFunc.domain->translationLab).transpose();
    try {
      values[inode] = extFEMFunc.evaluate( posExt );
    } catch ( const std::invalid_argument &e ) {
      // Point outside of domain
      if (domain->activeNodes[inode] && not(ignoreOutside)) {
        throw;
      } else {
        values[inode] = -1;
      }
    }
  }
}

void Function::interpolateInactive( const Function &extFEMFunc, bool ignoreOutside ) {
  // code duplicated from interpolate
  Eigen::Vector3d posExt;

  vector<int> indicesInactive = domain->activeNodes.filterIndices( [](int v){return (v==0);});
  for (int inactiveNode : indicesInactive) {
    // Move to reference frame of external
    posExt = domain->mesh->pos.row(inactiveNode) + (domain->translationLab - extFEMFunc.domain->translationLab).transpose();
    try {
      values[inactiveNode] = extFEMFunc.evaluate( posExt );
    } catch (const std::invalid_argument &e) {
      if (not(ignoreOutside)) {
        throw;
      }
    }
  }
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
  Eigen::VectorXd vals = Eigen::VectorXd::Zero( domain->mesh->nnodes );

  for (int inode = 0; inode < domain->mesh->nnodes; inode++) {
    // Move to reference frame of external
    Eigen::Vector3d posExt = domain->posLab.row(inode) - extFEMFunc.domain->translationLab.transpose();
    try {
      vals[inode] = extFEMFunc.evaluate( posExt );
    } catch ( const std::invalid_argument &e ) {
      // Point outside of domain
      if (not( ignoreOutside )) {
        throw;
      }
    }
  }

  return fem::Function( domain, vals );
}

Function operator-(const Function& f1, const Function& f2) {
  if (f1.domain != f2.domain) {
    throw std::invalid_argument("Functions don't share domain.");
  }
  Eigen::VectorXd values = f1.values - f2.values;
  return Function( f1.domain, values );
}

}
