#ifndef PRINTER
#define PRINTER
#include "Hatch.h"
#include "../Problem.h"
class Printer : public HatchCollider {
  private:
    Problem *p;
  public:
    Printer( Problem *p, double width, double height, double depth=0.0 ) :
      HatchCollider( p->domain.mesh, width, height, depth )
    {
      this->p    = p;
    }
    void deposit( const Eigen::Vector3d &p1, const Eigen::Vector3d &p2, mesh::MeshTag<int> *activeEls = NULL, bool modifyValues = true ) {
      /*
       * Activate elements at hatch defined by points p1 and p2 in lab reference frame
       */
      vector<int> collidedEls = collide( p1-p->domain.translationLab, p2-p->domain.translationLab);

      if (not(activeEls)) { activeEls = &p->domain.activeElements; };
      for (int ielem : collidedEls) {
        (*activeEls)[ielem] = 1;
      }
      p->domain.setActivation((*activeEls));
      if (modifyValues) {
        // Set deposition temperature at just activated nodes
        // TODO: Maybe compute MeshTag only once to avoid filtering it
        ConstantFunction depositionTemperature = ConstantFunction( &p->domain, p->Tdeposition );
        p->unknown.interpolate( depositionTemperature, p->domain.activeNodes,  [](int inode){return (inode==2);}, false );
        for (fem::Function& prevVal : p->previousValues) {
          prevVal.interpolate( depositionTemperature, p->domain.activeNodes,  [](int inode){return (inode==2);}, false );
        }
      }
    }

    void melt( const Eigen::Vector3d &p1, const Eigen::Vector3d &p2, mesh::MeshTag<int> *materialTag = NULL,
               size_t idxTargetSet = 0) {
      /*
       * Change material tag of elements at hatch in lab reference frame
       */
      vector<int> collidedEls = collide( p1-p->domain.translationLab, p2-p->domain.translationLab);

      if (not(materialTag)) { materialTag = &p->domain.materialTag; };
      for (int ielem : collidedEls) {
        (*materialTag)[ielem] = idxTargetSet;
      }
    }
};
#endif
