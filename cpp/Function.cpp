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
  int idxOwnerEl = mesh->findOwnerElement( point );
  if (idxOwnerEl < 0) {// Point outside of mesh
    return -1;
  }
  mesh::Element e = mesh->getElement( idxOwnerEl );//Load element containing point
  Eigen::VectorXd shaFunVals = e.evalShaFuns( point );

  val = values( e.con ).dot( shaFunVals );

  return val;
}

Eigen::Vector3d Function::evalGrad( Eigen::Vector3d point ) {
  /*
  Output gradient of Function @ input point
  */
  Eigen::Vector3d grad;
  grad.setZero();

  // GET VALS OF GRAD o SHAPE FUNCS AT POINT
  int idxOwnerEl = mesh->findOwnerElement( point );
  if (idxOwnerEl < 0) {// Point outside of mesh
    cout << "Point outside domain in grad FEM fun eval!" << endl;
    exit(-1);
  }
  mesh::Element e = mesh->getElement( idxOwnerEl );//Load element containing point
  Eigen::MatrixXd gradShaFunVals = e.evalGradShaFuns( point );

  for (int inode = 0; inode < e.nnodes; ++inode) {
    grad += gradShaFunVals.row(inode) * values[ e.con[ inode ] ] ;
  }

  return grad;
}


void Function::interpolate( Function &extFEMFunc ){

  Eigen::Vector3d posExt;
  values.setZero();//TODO: Rethink!

  for (int inode = 0; inode < mesh->nnodes; inode++) {
    // Move to reference frame of external
    posExt = mesh->pos.row(inode) + (mesh->shiftFRF - extFEMFunc.mesh->shiftFRF).transpose();
    values[inode] = extFEMFunc.evaluate( posExt );
  }
}

void interpolate( list<Function> &targetFunctions, const list<Function> &sourceFunctions ) {
  // Check same size
  if (targetFunctions.size() != sourceFunctions.size() ){
    cout << "Mismatch in number of funs to interpolate" << endl;
    exit(-1);
  }
  // Extract adress to target mesh and source mesh
  std::list<Function>::iterator targetFunsIterator = targetFunctions.begin();
  std::list<Function>::const_iterator sourceFunsIterator = sourceFunctions.begin();
  mesh::Mesh *targetMesh = targetFunsIterator->mesh;
  mesh::Mesh *sourceMesh = sourceFunsIterator->mesh;

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

void Function::interpolate2dirichlet( Function &extFEMFunc) {
  Function fh = Function( *mesh );
  fh.interpolate( extFEMFunc );
  for (int inode = 0; inode < mesh->nnodes; inode++) {
    if (fh.values(inode) >= 0) {//If interpolated
      values(inode) = fh.values(inode);
      dirichletNodes.push_back( inode );
      dirichletValues.push_back( fh.values(inode) );
    }
  }
}
}
