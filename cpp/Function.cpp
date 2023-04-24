#include "Function.h"

double Function::evalVal( Eigen::Vector3d point ) {
  /*
  Output val of Function at input point
  TODO: Nest this in evalValNPrevVals
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

vector<double> Function::evalValNPrevVals( Eigen::Vector3d point ) {
  /*
  Output val and previous vals of Function at input point
  */
  vector<double> vals(1+prevValues.cols());

  // Get values of shape funcs at point
  int idxOwnerEl = mesh->findOwnerElement( point );
  if (idxOwnerEl < 0) {// Point outside of mesh
    fill(vals.begin(), vals.end(), -1);
    return vals;
  }
  mesh::Element e = mesh->getElement( idxOwnerEl );//Load element containing point
  Eigen::VectorXd shaFunVals = e.evalShaFuns( point );

  vals[0] = values( e.con ).dot( shaFunVals );
  for (int icol=0; icol < prevValues.cols(); ++icol) {
    vals[icol+1] = prevValues( e.con, icol ).dot( shaFunVals );
  }

  return vals;
}

void Function::getFromExternal( Function &extFEMFunc ){
  Eigen::Vector3d posExt;
  values.setZero();
  prevValues.setZero();
  vector<double> valsAtPoint( 1+prevValues.cols() );
  for (int inode = 0; inode < mesh->nnodes; inode++) {
    // MOVE TO REFERENCE FRAME OF EXTERNAL
    posExt = mesh->pos.row(inode) + (mesh->shiftFRF - extFEMFunc.mesh->shiftFRF).transpose();
    valsAtPoint = extFEMFunc.evalValNPrevVals( posExt );
    values[inode] = valsAtPoint[0];
    for (int icol = 0; icol < prevValues.cols(); ++icol) {
      prevValues(inode, icol) = valsAtPoint[icol+1];
    }
  }
}

void Function::forceFromExternal( Function &extFEMFunc) {
  Function fh = Function( *mesh, extFEMFunc.nStepsRequired );
  fh.getFromExternal( extFEMFunc );
  for (int inode = 0; inode < mesh->nnodes; inode++) {
    if (fh.values(inode) >= 0) {//If interpolated
      values(inode) = fh.values(inode);
      dirichletNodes.push_back( inode );
      dirichletValues.push_back( fh.values(inode) );
    }
  }
}
