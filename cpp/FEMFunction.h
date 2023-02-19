#include "mesh/Mesh.h"
#include "mesh/Element.h"
#include <Eigen/Core>

class FEMFunction{
  public:
    mesh::Mesh* mesh;
    Eigen::VectorXd values;
    Eigen::MatrixXd prevValues;

    double interpolate( Eigen::Vector3d point ) {
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
};
