#include "Domain.h"
#include "../Problem.h"

namespace mesh {
Domain::Domain(Mesh *m, Problem *p) :
  activeNodes(mesh::MeshTag<int>(m, 0, 1)),
  activeElements(mesh::MeshTag<int>(m, m->dim, 1)),
  justDeactivatedElements(mesh::MeshTag<int>(m, m->dim, 0)),
  justActivatedBoundary(mesh::MeshTag<int>(m, m->dim-1, 0)),
  boundaryFacets(mesh::MeshTag<int>(m, m->dim-1)),
  boundaryFacetsParentEls(mesh::MeshTag<int>(m, m->dim-1))
{
  mesh = m;
  _dim = mesh->dim;
  posLab = mesh->pos;
  problem = p;
  computeBoundary();
}

void Domain::computeBoundary() {
  /*
   * Build array of indices of boundary facets
  */
  fill(boundaryFacets.x.begin(), boundaryFacets.x.end(), 0);
  fill(boundaryFacetsParentEls.x.begin(), boundaryFacetsParentEls.x.end(), 0);
  int activeElsPerFacet;
  bool incident2JustDeactivated;
  int lastVisitedActiveEl;
  for (int ifacet = 0; ifacet < mesh->con_FacetCell.nels_oDim; ++ifacet) {
    activeElsPerFacet = 0;
    incident2JustDeactivated = false;
    const vector<unsigned int>* incidentElements = mesh->con_FacetCell.getLocalCon( ifacet );
    for (int ielem : *incidentElements ) {
      if (activeElements.x[ielem]) {
        ++activeElsPerFacet;
        lastVisitedActiveEl = ielem;
        continue;
      }
      if (not(incident2JustDeactivated) && 
          (justDeactivatedElements.x[ielem]==1)) {
        incident2JustDeactivated = true;
        continue;
      }
    }
    if (activeElsPerFacet==1) {
      boundaryFacets[ ifacet ] = 1;
      boundaryFacetsParentEls[ ifacet ] = lastVisitedActiveEl;
      if (incident2JustDeactivated) {
        justActivatedBoundary.x[ifacet] = 1;
      }
    }
  }
}

void Domain::updateActiveNodes(const MeshTag<int> *newActiveNodes) {
  if (newActiveNodes) {
    for (int inode = 0; inode < mesh->nnodes; ++inode ) {
      if (not(activeNodes.x[inode]) && (newActiveNodes->x[inode]) ) {
        activeNodes.x[inode] = 2;
      } else {
        this->activeNodes[inode] = newActiveNodes->x[inode];
      }
    }
  } else {
    // Update activeNodes after a change in activeElements
    // If a node belongs to an active element, set it to active.
    const vector<unsigned int>* incidentElements;
    for (int inode = 0; inode < activeNodes.size(); ++inode) {
      bool isInactive = true;
      incidentElements = mesh->con_PointCell.getLocalCon( inode );
      for ( auto p_ielem = incidentElements->begin(); p_ielem != incidentElements->end(); ++p_ielem ) {
        if (activeElements.x[*p_ielem]) {
          if (not(activeNodes.x[inode])) {
            activeNodes.x[inode] = 2;
          } else {
            activeNodes.x[inode] = 1;
          }
          isInactive = false;
          break;
        }
      }
      if (isInactive) { activeNodes[inode] = 0; };
    }
  }
}

void Domain::updateActiveElements(const MeshTag<int> *newActiveEls) {
  if (newActiveEls) {
    for (int iel = 0; iel < mesh->nels; ++iel ) {
      if ((activeElements.x[iel])&& not(newActiveEls->x[iel]) ) {
        justDeactivatedElements.x[iel] = 1;
      }
    }
    activeElements = *newActiveEls;
  } else {
    //Update activeElements after a change in activeNodes
    //If all the nodes of an element are active, activate it.
    const vector<unsigned int>* incidentNodes;
    bool allNodesActive;
    for (int ielem = 0; ielem < activeElements.size(); ++ielem) {
      incidentNodes = mesh->con_CellPoint.getLocalCon( ielem );
      allNodesActive = true;
      for ( auto p_inode = incidentNodes->begin(); (p_inode != incidentNodes->end())&&(*p_inode != -1); ++p_inode ) {

        if (!activeNodes.x[*p_inode]) {
          allNodesActive = false;
          break;
        }
      }
      if (allNodesActive) {
        activeElements.x[ielem] = 1;
      } else {
        if (activeElements.x[ielem]) {
          justDeactivatedElements.x[ielem] = 1;
        }
        activeElements.x[ielem] = 0;
      }
    }
  }
}

bool checkHasInactive( const vector<int> &activeElements ) {
  return (std::find( activeElements.begin(), activeElements.end(), 0) != activeElements.end() );
}

void Domain::updateBeforeActivation() {
  fill(justActivatedBoundary.x.begin(), justActivatedBoundary.x.end(), 0);
  fill(justDeactivatedElements.x.begin(), justDeactivatedElements.x.end(), 0);
}

void Domain::updateAfterActivation() {
  hasInactive = (std::find( activeElements.x.begin(), activeElements.x.end(), false) != activeElements.x.end() );
  computeBoundary();
}

void Domain::setActivation(const MeshTag<int> &activationCriterion) {
  updateBeforeActivation();
  if (activationCriterion.dim()==0) {
    // Activate by nodes. Is this useful?
    updateActiveNodes(&activationCriterion);
    updateActiveElements();
  } else if (activationCriterion.dim()==_dim) {
    updateActiveElements(&activationCriterion);
    updateActiveNodes();
  } else {
    throw std::invalid_argument("Activation criterion must be node or element tag.");
  }
  updateAfterActivation();
}

void Domain::resetActivation() {
  MeshTag<int> allActive = MeshTag<int>( mesh, mesh->dim, 1);
  setActivation( allActive );
}

void Domain::deactivate() {
  MeshTag<int> allInactive = MeshTag<int>( mesh, mesh->dim, 0);
  setActivation( allInactive );
}

int Domain::findOwnerElements( const Eigen::Vector3d &point ) const {
  /*
   * Wrapper around mesh's findOwnerElements that returns
   * an active element owning the point.
   */
  vector<int> owners = mesh->findOwnerElements( point );
  for (int ielem : owners) {
    if ( activeElements[ielem] ) {
      return ielem;
    }
  }
  throw std::invalid_argument("Point is not owned by active element.");
}


void Domain::intersect( const MeshTag<int> activeElements) {
  this->activeElements &= activeElements;
  setActivation( activeElements );
}

void Domain::preIterate() {
  translationLab += problem->dt * speedDomain;
  for (int inode=0; inode < mesh->nnodes; inode++){
    posLab.row( inode ) += problem->dt * speedDomain;
  }
}

void Domain::setSpeed(Eigen::Vector3d speedDomain){ this->speedDomain = speedDomain; }

}
