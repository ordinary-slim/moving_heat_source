#include "Problem.h"
#include "mesh/cgal_interface.h"

class Printer {
  private:
    Problem *p;
  public:
    double mdwidth, mdheight;
    myOBB obb;

    Printer( Problem *p, double mdwidth, double mdheight ) :
      obb( Eigen::Vector3d(-1, 0, 0), Eigen::Vector3d(+1, 0, 0), mdwidth, mdheight )
    {
      this->p    = p;
      this->mdwidth = mdwidth;
      this->mdheight = mdheight;
    }

    vector<int> collide( const Eigen::Vector3d &p1, const Eigen::Vector3d &p2 ) {
      obb = myOBB( p1, p2, mdwidth, mdheight );
      return this->p->domain.mesh->findCollidingElement( obb );
    }

    void deposit( const Eigen::Vector3d &p1, const Eigen::Vector3d &p2, mesh::MeshTag<int> *activeEls = NULL ) {
      vector<int> collidedEls = collide( p1-p->domain.mesh->shiftFRF, p2-p->domain.mesh->shiftFRF);

      if (not(activeEls)) { activeEls = &p->domain.activeElements; };
      for (int ielem : collidedEls) {
        (*activeEls)[ielem] = 1;
      }
      p->domain.setActivation((*activeEls));
      // Set deposition temperature
      vector<int> indicesJustActivated =
        p->domain.activeNodes.filterIndices( [](int inode){return (inode==2);});
      for (int inode : indicesJustActivated) {
        p->unknown.values[inode] = p->Tdeposition;
      }
      for (fem::Function& prevVal : p->previousValues) {
        for (int inode : indicesJustActivated) {
          prevVal.values[inode] = p->Tdeposition;
        }
      }
    }
};
