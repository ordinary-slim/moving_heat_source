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
      vector<int> collidedEls = collide( p1-p->domain.translationLab, p2-p->domain.translationLab);

      if (not(activeEls)) { activeEls = &p->domain.activeElements; };
      for (int ielem : collidedEls) {
        (*activeEls)[ielem] = 1;
      }
      p->domain.setActivation((*activeEls));
      // Set deposition temperature at just activated nodes
      // TODO: Maybe compute MeshTag only once to avoid filtering it
      ConstantFunction depositionTemperature = ConstantFunction( &p->domain, p->Tdeposition );
      p->unknown.interpolate( depositionTemperature, p->domain.activeNodes,  [](int inode){return (inode==2);}, false );
      for (fem::Function& prevVal : p->previousValues) {
        prevVal.interpolate( depositionTemperature, p->domain.activeNodes,  [](int inode){return (inode==2);}, false );
      }
    }
};
#endif
