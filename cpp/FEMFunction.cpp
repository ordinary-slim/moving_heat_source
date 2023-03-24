#include "FEMFunction.h"

double FEMFunction::evaluateVal( Eigen::Vector3d point ) {
  // Output val of FEMFunction at input point
  double val = 0;

  // GET VALS OF SHAPE FUNCS AT POINT
  int idxOwnerEl = mesh->findOwnerElement( point );
  if (idxOwnerEl < 0) {// Point outside of mesh
    return -1;
  }
  mesh::Element e = mesh->getElement( idxOwnerEl );//Load element containing point
  Eigen::VectorXd shaFunVals = e.evaluateShaFuns( point );

  val = values( e.con ).dot( shaFunVals );

  return val;
}

vector<double> FEMFunction::evaluateValNPrevVals( Eigen::Vector3d point ) {
  // Output val and previous vals of FEMFunction at input point
  vector<double> vals(1+prevValues.cols());

  // GET VALS OF SHAPE FUNCS AT POINT
  int idxOwnerEl = mesh->findOwnerElement( point );
  if (idxOwnerEl < 0) {// Point outside of mesh
    fill(vals.begin(), vals.end(), -1);
    return vals;
  }
  mesh::Element e = mesh->getElement( idxOwnerEl );//Load element containing point
  Eigen::VectorXd shaFunVals = e.evaluateShaFuns( point );

  vals[0] = values( e.con ).dot( shaFunVals );
  for (int icol=0; icol < prevValues.cols(); ++icol) {
    vals[icol+1] = prevValues( e.con, icol ).dot( shaFunVals );
  }

  return vals;
}

void FEMFunction::getFromExternal( FEMFunction &extFEMFunc ){
  Eigen::Vector3d posExt;
  values.setZero();
  prevValues.setZero();
  vector<double> valsAtPoint( 1+prevValues.cols() );
  cout << "Hello?" << endl;
  for (int inode = 0; inode < mesh->nnodes; inode++) {
    // MOVE TO REFERENCE FRAME OF EXTERNAL
    posExt = mesh->pos.row(inode) + (mesh->shiftFRF - extFEMFunc.mesh->shiftFRF).transpose();
    valsAtPoint = extFEMFunc.evaluateValNPrevVals( posExt );
    values[inode] = valsAtPoint[0];
    for (int icol = 0; icol < prevValues.cols(); ++icol) {
      prevValues(inode, icol) = valsAtPoint[icol+1];
    }
  }
}
