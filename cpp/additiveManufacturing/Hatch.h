#ifndef HATCH
#define HATCH
#include "../mesh/cgal_interface.h"
#include "../mesh/Mesh.h"

class HatchCollider {
  private:
    mesh::Mesh* m;
  public:
    double width;
    double height, depth;
    myOBB obb;//collision detection with mesh elements

    HatchCollider( mesh::Mesh* m, double width, double height, double depth = 0.0 ) {
      this->m      = m;
      this->width  = width;
      this->height = height;
      this->depth  = depth;
    }

    std::vector<int> collide( const Eigen::Vector3d &p1, const Eigen::Vector3d &p2 ) {
      /*
       * Compute mesh elements that collide with hatch p1 *------*p2
       */
      obb = myOBB( p1, p2, width, height, depth, true );
      return this->m->findCollidingElements( obb );
    }

};
#endif
