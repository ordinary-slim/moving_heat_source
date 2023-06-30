#ifndef HATCH
#define HATCH
#include "../mesh/cgal_interface.h"
#include "../mesh/Mesh.h"

class HatchCollider {
  private:
    mesh::Mesh* m;
  public:
    double width, height;
    myOBB obb;//collision detection with mesh elements

    HatchCollider( mesh::Mesh* m, double width, double height ) :
      obb( Eigen::Vector3d(-1, 0, 0), Eigen::Vector3d(+1, 0, 0), width, height )
    {
      this->m      = m;
      this->width  = width;
      this->height = height;
    }

    std::vector<int> collide( const Eigen::Vector3d &p1, const Eigen::Vector3d &p2 ) {
      /*
       * Compute mesh elements that collide with hatch p1 *------*p2
       */
      obb = myOBB( p1, p2, width, height );
      return this->m->findCollidingElement( obb );
    }

};
#endif
