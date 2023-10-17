#ifndef cgal_interface_h
#define cgal_interface_h
#include "Element.h"
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <algorithm>

using inex_K      = CGAL::Simple_cartesian<double>;
using Plane3_CGAL = inex_K::Plane_3;

struct MyAABB {
  int ielem = -1;
  double bounds[3][2] = {
    {1, -1},
    {1, -1},
    {1, -1}
  };
  constexpr static double stretch = 0.05;

  MyAABB() = default;
  MyAABB( mesh::Element e, double pad = 1e-7 );

};

struct MyBboxPrimitive {
public:
    typedef const MyAABB* Pointer;
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
    MyBboxPrimitive() {} // default constructor needed
    // the following constructor is the one that receives the iterators from the
    // iterator range given as input to the AABB_tree
    MyBboxPrimitive(std::vector<MyAABB>::const_iterator it)
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

struct MyOBB{
  public:
    Eigen::Vector3d pos, halfWidths, xAxis, yAxis, zAxis;

    MyOBB() {
      pos << 0.0, 0.0, 0.0;
      halfWidths << 0.0, 0.0, 0.0;
      xAxis << 1.0, 0.0, 0.0;
      yAxis << 0.0, 1.0, 0.0;
      zAxis << 0.0, 0.0, 1.0;
    }

    // 3D-printing constructors
    MyOBB(Eigen::Vector3d p1, Eigen::Vector3d p2, double width, 
        double height, int dim = 3, bool shrink = true);

    MyOBB(Eigen::Vector3d p1, Eigen::Vector3d p2, double width, 
        double aboveLen, double belowLen, int dim = 3, bool shrink = true);

    void setTransverseAxes(int dim = 3);

    void appendPlanes(std::vector<Plane3_CGAL> &v) const;

    explicit operator inex_K::Iso_cuboid_3() const
    {
        double minX = pos[0] - halfWidths[0]*abs(xAxis[0]) - halfWidths[1]*abs(yAxis[0]) - halfWidths[2]*abs(zAxis[0]);
        double maxX = pos[0] + halfWidths[0]*abs(xAxis[0]) + halfWidths[1]*abs(yAxis[0]) + halfWidths[2]*abs(zAxis[0]);
        if (minX > maxX) { std::swap( minX, maxX ); }
        double minY = pos[1] - halfWidths[0]*abs(xAxis[1]) - halfWidths[1]*abs(yAxis[1]) - halfWidths[2]*abs(zAxis[1]);
        double maxY = pos[1] + halfWidths[0]*abs(xAxis[1]) + halfWidths[1]*abs(yAxis[1]) + halfWidths[2]*abs(zAxis[1]);
        if (minY > maxY) { std::swap( minY, maxY ); }
        double minZ = pos[2] - halfWidths[0]*abs(xAxis[2]) - halfWidths[1]*abs(yAxis[2]) - halfWidths[2]*abs(zAxis[2]);
        double maxZ = pos[2] + halfWidths[0]*abs(xAxis[2]) + halfWidths[1]*abs(yAxis[2]) + halfWidths[2]*abs(zAxis[2]);
        if (minZ > maxZ) { std::swap( minZ, maxZ ); }
        return inex_K::Iso_cuboid_3( minX, minY, minZ, maxX, maxY, maxZ );
    }
    //TODO: template this
    bool hasCollided(const mesh::Element &otherConvex) const;
};

#endif
