#ifndef PRINTER
#define PRINTER
#include "Hatch.h"
#include "../Problem.h"
class Printer : public HatchCollider {
  private:
    Problem *p;
  public:
    Printer( Problem *p, double width, double height ) :
      HatchCollider( p->domain.mesh, width, height )
    {
      this->p    = p;
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
#endif
