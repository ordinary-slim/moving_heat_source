#ifndef cgal_interface_h
#define cgal_interface_h
#include "Element.h"
#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>

struct myAABB {
  int ielem;
  double bounds[3][2];

  myAABB( mesh::Element e ) {
    for (int idim = 0; idim < 3; ++idim) {
      bounds[idim][0] = e.pos.col(idim).minCoeff();
      bounds[idim][1] = e.pos.col(idim).maxCoeff();
    }
    this->ielem = e.ient;
  }

};

struct myBboxPrimitive {
public:
    typedef const myAABB* Pointer;
    // this is the type of data that the queries returns
    using Id = int;
    // CGAL types returned
    using K = CGAL::Simple_cartesian<double>;
    typedef K::Point_3    Point; // CGAL 3D point type
    typedef K::Iso_cuboid_3 Datum;
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
#endif
