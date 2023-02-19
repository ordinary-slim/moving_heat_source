#include "element.h"
#include <Eigen/Core>

struct AABB {
  double bounds[3][2];
  int dim;

  AABB() {
  }

  AABB( Element e ) {
    dim = e.dim;
    for (int idim = 0; idim < dim; ++idim) {
      bounds[idim][0] = e.pos.col(idim).minCoeff();
      bounds[idim][1] = e.pos.col(idim).maxCoeff();
    }
  }

  bool isPointInside( Eigen::Vector3d point ) {
    bool isInside = true;
    for (int idim = 0; idim < dim; ++idim) {
      if ( (point[idim] < bounds[idim][0]) || (point[idim] > bounds[idim][1]) ) {
        isInside = false;
        break;
      }
    }
    return isInside;
  }
};
