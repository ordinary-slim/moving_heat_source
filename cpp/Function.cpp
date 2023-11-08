#include <list>
#include "Function.h"

namespace fem
{
double Function::evaluate( Eigen::Vector3d &point ) const {
  /*
  Output val of Function at input point
  */
  // GET VALS OF SHAPE FUNCS AT POINT
  int idxOwnerEl = domain->findOwnerElements( point );
  mesh::Element e = domain->mesh->getElement( idxOwnerEl );//Load element containing point
  Eigen::VectorXd shaFunVals = e.evaluateShaFuns( point );

  return  values( *e.con ).dot( shaFunVals );
}

Eigen::Vector3d Function::evaluateGrad( Eigen::Vector3d &point ) {
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


void Function::interpolate(const AbstractFunction &extFEMFunc, const mesh::MeshTag<int> &nodalTag,
    std::function<bool(int)> filter, bool ignoreOutside ) {

  if (nodalTag.dim() != 0) { throw std::invalid_argument("Interpolate requires nodal mesh tag."); }

  vector<int> nodesOfInterest;
  if (filter == nullptr) {
     nodesOfInterest = nodalTag.getIndices();
  } else {
     nodesOfInterest = nodalTag.filterIndices( filter );
  }

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
        values[inode] = -1;// sentinel value, if it shows up in active
                           // nodes of simulation something went wrong
      }
    }
  }
}

void Function::interpolate(const AbstractFunction &extFEMFunc, const mesh::MeshTag<int> &nodalTag, bool ignoreOutside ) {
  interpolate( extFEMFunc, nodalTag, nullptr, ignoreOutside );
}

void Function::interpolate(const AbstractFunction &extFEMFunc, bool ignoreOutside ) {
  mesh::MeshTag<int> allNodes = mesh::MeshTag<int>( domain->mesh, 0, 1 );
  interpolate( extFEMFunc, allNodes, nullptr, ignoreOutside );
}

void Function::interpolateInactive( const AbstractFunction &extFEMFunc, bool ignoreOutside ) {
  // Shorthand
  interpolate( extFEMFunc,
      domain->activeNodes,
      [](int val){ return (val == 0); },
      ignoreOutside );
}

double Function::getL2Norm() const {
  Eigen::VectorXd valuesAtActiveNodes(domain->ls->getNdofs());
  for (int inode = 0; inode < domain->mesh->nnodes; ++inode) {
    int idof = domain->dofNumbering[inode];
    if (idof < 0) {
      continue;
    }
    valuesAtActiveNodes[idof] = values[inode];
  }
  return sqrt( valuesAtActiveNodes.dot( (*domain->massMat) * valuesAtActiveNodes ) );
}

fem::Function interpolate( const AbstractFunction &extFEMFunc,
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

Function operator/(const Function& f, const double c) {
  Eigen::VectorXd values = f.values / c;
  return Function( f.domain, values );
}

double DG0Function::evaluate(Eigen::Vector3d &point) const {
  int idxOwnerEl = domain->findOwnerElements( point );
  return values[idxOwnerEl];
}

void DG0Function::interpolate(const AbstractFunction &extFEMFunc,
    const mesh::MeshTag<int> &cellTag, std::function<bool(int)> filter,
    bool ignoreOutside) {

  if (cellTag.dim() != domain->mesh->dim ) { throw std::invalid_argument("Interpolate requires cell mesh tag."); }

  vector<int> cellsOfInterest;
  if (filter == nullptr) {
     cellsOfInterest = cellTag.getIndices();
  } else {
     cellsOfInterest = cellTag.filterIndices( filter );
  }

  for (int ielem : cellsOfInterest) {
    // Move to reference frame of external
    mesh::Element e = domain->getElement( ielem );
    Eigen::Vector3d posExt = e.centroid + domain->translationLab - extFEMFunc.domain->translationLab;
    try {
      values[ielem] = extFEMFunc.evaluate( posExt );
    } catch ( const std::invalid_argument &e ) {
      // Point outside of domain
      if (domain->activeElements[ielem] && not(ignoreOutside)) {
        throw;
      } else {
        values[ielem] = -1;// sentinel value, if it shows up in active
                           // elements of simulation something went wrong
      }
    }
  }
}

}
