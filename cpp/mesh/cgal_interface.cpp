#include "cgal_interface.h"
#include <CGAL/Convex_hull_3/dual/halfspace_intersection_interior_point_3.h>

myAABB::myAABB( mesh::Element e ) {
  for (int idim = 0; idim < 3; ++idim) {
    double min = e.pos.col(idim).minCoeff();
    double max = e.pos.col(idim).maxCoeff();
    double L   = max - min;
    bounds[idim][0] = min - (stretch/2)*L - pad;
    bounds[idim][1] = max + (stretch/2)*L + pad;
  }
  this->ielem = e.ient;
}


bool myOBB::hasCollided(const mesh::Element &otherConvex) const {
  bool hasCollided = false;
  std::vector<Plane3_CGAL> halfSpaces;
  appendPlanes( halfSpaces );
  otherConvex.appendPlanes( halfSpaces );
  auto p = CGAL::halfspace_intersection_interior_point_3(halfSpaces.begin(), halfSpaces.end());
  return not(p==boost::none);
};
