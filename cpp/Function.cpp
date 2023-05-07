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
  int idxOwnerEl = domain->findOwnerElement( point );
  if (idxOwnerEl < 0) {// Point outside of domain->mesh
    return -1;
  }
  mesh::Element e = domain->mesh->getElement( idxOwnerEl );//Load element containing point
  Eigen::VectorXd shaFunVals = e.evalShaFuns( point );

  val = values( *e.con ).dot( shaFunVals );

  return val;
}

Eigen::Vector3d Function::evaluateGrad( Eigen::Vector3d point ) {
  /*
  Output gradient of Function @ input point
  */
  Eigen::Vector3d grad = Eigen::Vector3d::Zero();

  // GET VALS OF GRAD o SHAPE FUNCS AT POINT
  int idxOwnerEl = domain->findOwnerElement( point );
  if (idxOwnerEl < 0) {// Point outside of domain->mesh
    cout << "Point outside domain in grad FEM fun eval!" << endl;
    exit(-1);
  }
  mesh::Element e = domain->mesh->getElement( idxOwnerEl );//Load element containing point
  Eigen::MatrixXd gradShaFunVals = e.evaluateGradShaFuns( point );

  for (int inode = 0; inode < e.nnodes; ++inode) {
    grad += gradShaFunVals.row(inode) * values[ (*e.con)[ inode ] ] ;
  }

  return grad;
}


void Function::interpolate( Function &extFEMFunc ){

  Eigen::Vector3d posExt;
  values.setZero();//TODO: Rethink!

  for (int inode = 0; inode < domain->mesh->nnodes; inode++) {
    // Move to reference frame of external
    posExt = domain->mesh->pos.row(inode) + (domain->mesh->shiftFRF - extFEMFunc.domain->mesh->shiftFRF).transpose();
    values[inode] = extFEMFunc.evaluate( posExt );
  }
}

void interpolate( list<Function> &targetFunctions, const list<Function> &sourceFunctions ) {
  // Check same size
  if (targetFunctions.size() != sourceFunctions.size() ){
    cout << "Mismatch in number of funs to interpolate" << endl;
    exit(-1);
  }
  // Extract adress to target domain->mesh and source domain->mesh
  std::list<Function>::iterator targetFunsIterator = targetFunctions.begin();
  std::list<Function>::const_iterator sourceFunsIterator = sourceFunctions.begin();
  mesh::Mesh *targetMesh = targetFunsIterator->domain->mesh;
  mesh::Mesh *sourceMesh = sourceFunsIterator->domain->mesh;

  Eigen::Vector3d posExt;
  for (int inode = 0; inode < targetMesh->nnodes; inode++) {
    // Move to reference frame of external
    posExt = targetMesh->pos.row(inode) + (targetMesh->shiftFRF - sourceMesh->shiftFRF).transpose();
    while (targetFunsIterator != targetFunctions.end() ) {
      targetFunsIterator->values[inode] = sourceFunsIterator->evaluate( posExt );
      ++sourceFunsIterator;
      ++targetFunsIterator;
    }
    // Reset iterators
    targetFunsIterator = targetFunctions.begin();
    sourceFunsIterator = sourceFunctions.begin();
  }
}

}
