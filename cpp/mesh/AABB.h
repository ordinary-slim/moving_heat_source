#ifndef AABB_H
#define AABB_H
#include "Element.h"
#include <Eigen/Core>

struct AABB {
  double bounds[3][2];
  int dim;
  double pad = 1e-5;

  AABB() {
  }

  AABB( mesh::Element e ) {
    dim = e.dim;
    for (int idim = 0; idim < dim; ++idim) {
      bounds[idim][0] = e.pos.col(idim).minCoeff() - pad;
      bounds[idim][1] = e.pos.col(idim).maxCoeff() + pad;
    }
  }

  bool isPointInside( const Eigen::Vector3d &point ) const {
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
#endif
