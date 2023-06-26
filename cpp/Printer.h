#include "mesh/Mesh.h"
#include "mesh/cgal_interface.h"

class Printer {
  public:
    double mdwidth, mdheight;
    mesh::Mesh *mesh;
    myOBB obb;

    Printer( mesh::Mesh *mesh, double mdwidth, double mdheight ) :
      obb( Eigen::Vector3d(-1, 0, 0), Eigen::Vector3d(+1, 0, 0), mdwidth, mdheight )
    {
      this->mesh    = mesh;
      this->mdwidth = mdwidth;
      this->mdheight = mdheight;
    }

    mesh::MeshTag<int> mark( const Eigen::Vector3d &p1, const Eigen::Vector3d &p2 ) {
      obb = myOBB( p1, p2, mdwidth, mdheight );
      vector<int> collidedEls = this->mesh->findCollidingElement( obb );
      return mesh::mark( *mesh, mesh->dim, collidedEls );
    }

};
