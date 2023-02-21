#ifndef FEMFUNC
#include "mesh/Mesh.h"
#include "mesh/Element.h"
#include <Eigen/Core>

class FEMFunction{
  public:
    mesh::Mesh* mesh;
    Eigen::VectorXd values;
    Eigen::MatrixXd prevValues;//col 0 is previous it, col 1 is prev prev it, etc

    double evaluateVal( Eigen::Vector3d point ) {
      // Output val of FEMFunction at input point
      double val = 0;

      // GET VALS OF SHAPE FUNCS AT POINT
      int idxOwnerEl = mesh->findOwnerElement( point );
      if (idxOwnerEl < 0) {// Point outside of mesh
        return 0.0;
      }
      Element e = mesh->getElement( idxOwnerEl );//Load element containing point
      Eigen::VectorXd shaFunVals = e.evaluateShaFuns( point );

      val = values( e.con ).dot( shaFunVals );

      return val;
    }
    vector<double> evaluateValNPrevVals( Eigen::Vector3d point ) {
      // Output val and previous vals of FEMFunction at input point
      vector<double> vals(1+prevValues.cols());

      // GET VALS OF SHAPE FUNCS AT POINT
      int idxOwnerEl = mesh->findOwnerElement( point );
      if (idxOwnerEl < 0) {// Point outside of mesh
        fill(vals.begin(), vals.end(), 0.0);
        return vals;
      }
      Element e = mesh->getElement( idxOwnerEl );//Load element containing point
      Eigen::VectorXd shaFunVals = e.evaluateShaFuns( point );

      vals[0] = values( e.con ).dot( shaFunVals );
      for (int icol=0; icol < prevValues.cols(); ++icol) {
        vals[icol+1] = prevValues( e.con, icol ).dot( shaFunVals );
      }

      return vals;
    }
    void getFromExternal( FEMFunction &extFEMFunc ){
      Eigen::Vector3d posExt;
      values.setZero();
      prevValues.setZero();
      vector<double> valsAtPoint( 1+prevValues.cols() );
      for (int inode = 0; inode < mesh->nnodes; inode++) {
        posExt = mesh->pos.row(inode) + (mesh->shiftFRF - extFEMFunc.mesh->shiftFRF).transpose();
        valsAtPoint = extFEMFunc.evaluateValNPrevVals( posExt );
        values[inode] = valsAtPoint[0];
        for (int icol = 0; icol < prevValues.cols(); ++icol) {
          prevValues(inode, icol) = valsAtPoint[icol+1];
        }
      }
    }

};
#define FEMFUNC
#endif
