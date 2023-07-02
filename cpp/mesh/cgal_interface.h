#ifndef cgal_interface_h
#define cgal_interface_h
#include "Element.h"
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <algorithm>

using inex_K      = CGAL::Simple_cartesian<double>;
using ex_inex_K   = CGAL::Exact_predicates_inexact_constructions_kernel;
using Plane3_CGAL = ex_inex_K::Plane_3;

struct myAABB {
  int ielem;
  double bounds[3][2];

  myAABB( mesh::Element e );

};

struct myBboxPrimitive {
public:
    typedef const myAABB* Pointer;
    // this is the type of data that the queries returns
    using Id = int;
    // CGAL types returned
    using inex_K = CGAL::Simple_cartesian<double>;
    typedef inex_K::Point_3    Point; // CGAL 3D point type
    typedef inex_K::Iso_cuboid_3 Datum;
private:
    Id ielem;// this is what the AABB tree stores internally
    Pointer pt;
public:
    myBboxPrimitive() {} // default constructor needed
    // the following constructor is the one that receives the iterators from the
    // iterator range given as input to the AABB_tree
    myBboxPrimitive(std::vector<myAABB>::const_iterator it)
        : pt(&(*it)) {
          ielem = it->ielem;
        }
    const Id& id() const { return ielem; }
    // on the fly conversion from the internal data to the CGAL types
    Datum datum() const
    {
        return Datum(pt->bounds[0][0],
                     pt->bounds[1][0],
                     pt->bounds[2][0],
                     pt->bounds[0][1],
                     pt->bounds[1][1],
                     pt->bounds[2][1]
                     );
    }
    // returns a reference point which must be on the primitive
    Point reference_point() const
    { return Point(pt->bounds[0][0], pt->bounds[1][0], pt->bounds[2][0]); }
};

struct myOBB{
  public:
    Eigen::Vector3d pos, halfWidths, xAxis, yAxis, zAxis;

    myOBB(Eigen::Vector3d p1, Eigen::Vector3d p2, double width, 
        double height){
      // Printing dir as xAxis
      Eigen::Vector3d step = (p2 - p1);
      pos = (p2 + p1)/2.0;
      xAxis = step.normalized();
      zAxis << 0.0, 0.0, 1.0;
      zAxis = (zAxis - zAxis.dot( xAxis )*xAxis).normalized();//Gram schmidt
      if (zAxis.norm() < 0.9) {
        throw std::invalid_argument("Steps in Z-axis not allowed.");
      }
      yAxis = zAxis.cross( xAxis );
      halfWidths(0) = step.norm() / 2.0;
      halfWidths(1) = width / 2.0;
      halfWidths(2) = height / 2.0;
    }
    void appendPlanes(std::vector<Plane3_CGAL> &v) const {
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
    explicit operator inex_K::Iso_cuboid_3() const
    {
        double minX = pos[0] - halfWidths[0]*xAxis[0] - halfWidths[1]*yAxis[0] - halfWidths[2]*zAxis[0];
        double maxX = pos[0] + halfWidths[0]*xAxis[0] + halfWidths[1]*yAxis[0] + halfWidths[2]*zAxis[0];
        if (minX > maxX) { std::swap( minX, maxX ); }
        double minY = pos[1] - halfWidths[0]*xAxis[1] - halfWidths[1]*yAxis[1] - halfWidths[2]*zAxis[1];
        double maxY = pos[1] + halfWidths[0]*xAxis[1] + halfWidths[1]*yAxis[1] + halfWidths[2]*zAxis[1];
        if (minY > maxY) { std::swap( minY, maxY ); }
        double minZ = pos[2] - halfWidths[0]*xAxis[2] - halfWidths[1]*yAxis[2] - halfWidths[2]*zAxis[2];
        double maxZ = pos[2] + halfWidths[0]*xAxis[2] + halfWidths[1]*yAxis[2] + halfWidths[2]*zAxis[2];
        if (minZ > maxZ) { std::swap( minZ, maxZ ); }
        return inex_K::Iso_cuboid_3( minX, minY, minZ, maxX, maxY, maxZ );
    }
    //TODO: template this
    bool hasCollided(const mesh::Element &otherConvex) const;
};

#endif
