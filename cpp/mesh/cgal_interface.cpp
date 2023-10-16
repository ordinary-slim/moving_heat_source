#include "cgal_interface.h"
#include <CGAL/Convex_hull_3/dual/halfspace_intersection_interior_point_3.h>

MyAABB::MyAABB( mesh::Element e, double pad ) {
  for (int idim = 0; idim < 3; ++idim) {
    double min = e.pos.col(idim).minCoeff();
    double max = e.pos.col(idim).maxCoeff();
    double L   = max - min;
    bounds[idim][0] = min - (stretch/2)*L - pad;
    bounds[idim][1] = max + (stretch/2)*L + pad;
  }
  this->ielem = e.ient;
}

MyOBB::MyOBB(Eigen::Vector3d p1, Eigen::Vector3d p2, double width, 
    double height, int dim, bool shrink){
  /*
   * 3D-printing constructor for OBB
   */
  Eigen::Vector3d step = (p2 - p1);
  xAxis = step.normalized();
  setTransverseAxes(dim);
  pos = (p2 + p1)/2.0;
  halfWidths(0) = step.norm() / 2.0;
  halfWidths(1) = width / 2.0;
  halfWidths(2) = height / 2.0;
  // Avoid touching collisions
  if (shrink) {
    halfWidths *= 0.999;
  }
}
MyOBB::MyOBB(Eigen::Vector3d p1, Eigen::Vector3d p2, double width, 
    double aboveLen, double belowLen, int dim, bool shrink) {
  /*
   * Another 3D-printing constructor for OBB
   * Height and depth for Z
   */
  Eigen::Vector3d step = (p2 - p1);
  xAxis = step.normalized();
  setTransverseAxes(dim);
  pos = (p2 + p1)/2.0;
  pos += ((aboveLen - belowLen)/2)*zAxis;
  halfWidths(0) = step.norm() / 2.0;
  halfWidths(1) = width / 2.0;
  halfWidths(2) = (aboveLen + belowLen) / 2;
  // Avoid touching collisions
  if (shrink) {
    halfWidths *= 0.999;
  }
  //std::cout << "pos = " << pos << "\nMy axes are \nxAxis=" << xAxis << "\nyAxis = "  << yAxis << "\nzAxis = "   << zAxis  << std::endl << "My halfWidths are " << halfWidths << std::endl;
}


bool MyOBB::hasCollided(const mesh::Element &otherConvex) const {
  bool hasCollided = false;
  std::vector<Plane3_CGAL> halfSpaces;
  appendPlanes( halfSpaces );
  otherConvex.appendPlanes( halfSpaces );
  auto p = CGAL::halfspace_intersection_interior_point_3(halfSpaces.begin(), halfSpaces.end());
  return not(p==boost::none);
};

void MyOBB::setTransverseAxes(int dim) {
  /*
   * Set transverse axes of OBB
   * Attempt to set 0.0, 0.0, 1.0 as zAxis
   */
  if (dim==2) {
    zAxis << 0.0, 1.0, 0.0;
  } else {
    zAxis << 0.0, 0.0, 1.0;
  }
  zAxis = (zAxis - zAxis.dot( xAxis )*xAxis).normalized();//Gram schmidt
  if ( zAxis.norm() < (1 - 1e-7) ) {
    throw std::invalid_argument("Steps in Z-axis not allowed.");
  }
  yAxis = zAxis.cross( xAxis );
}

void MyOBB::appendPlanes(std::vector<Plane3_CGAL> &v) const {
  v.reserve( std::min<int>(v.size(), 12) );
  // Mins
  v.push_back( Plane3_CGAL( -xAxis[0],
                       -xAxis[1],
                       -xAxis[2],
                       (pos.dot( xAxis ) - halfWidths[0])
                       ) );
  v.push_back( Plane3_CGAL( -yAxis[0],
                       -yAxis[1],
                       -yAxis[2],
                       (pos.dot( yAxis ) - halfWidths[1])
                       ) );
  v.push_back( Plane3_CGAL( -zAxis[0],
                       -zAxis[1],
                       -zAxis[2],
                       (pos.dot( zAxis ) - halfWidths[2])
                       ) );
  //Maxes
  v.push_back( Plane3_CGAL( +xAxis[0],
                       +xAxis[1],
                       +xAxis[2],
                       -(pos.dot( xAxis ) + halfWidths[0])
                       ) );
  v.push_back( Plane3_CGAL( +yAxis[0],
                       +yAxis[1],
                       +yAxis[2],
                       -(pos.dot( yAxis ) + halfWidths[1])
                       ) );
  v.push_back( Plane3_CGAL( +zAxis[0],
                       +zAxis[1],
                       +zAxis[2],
                       -(pos.dot( zAxis ) + halfWidths[2])
                       ) );
}
