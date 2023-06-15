#include "ActiveMesh.h"

namespace mesh {
void ActiveMesh::computeBoundary() {
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

void ActiveMesh::updateActiveNodes() {
  /*
   * Update activeNodes after a change in activeElements
   * If a node belongs to an active element, set it to active.
   */
  const vector<unsigned int>* incidentElements;
  for (int inode = 0; inode < activeNodes.size(); ++inode) {
    activeNodes.x[inode] = false;
    incidentElements = mesh->con_PointCell.getLocalCon( inode );
    for ( auto p_ielem = incidentElements->begin(); p_ielem != incidentElements->end(); ++p_ielem ) {
      if (activeElements.x[*p_ielem]) {
        activeNodes.x[inode] = true;
        break;
      }
    }
  }
}

void ActiveMesh::updateActiveElements() {
  /*
   * Update activeElements after a change in activeNodes
   * If all the nodes of an element are active, activate it.
   */
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

bool checkHasInactive( const vector<int> &activeElements ) {
  return (std::find( activeElements.begin(), activeElements.end(), 0) != activeElements.end() );
}

void ActiveMesh::updateBeforeActivation() {
  fill(justActivatedBoundary.x.begin(), justActivatedBoundary.x.end(), 0);
  fill(justDeactivatedElements.x.begin(), justDeactivatedElements.x.end(), 0);
}

void ActiveMesh::updateAfterActivation() {
  hasInactive = (std::find( activeElements.x.begin(), activeElements.x.end(), false) != activeElements.x.end() );
  if (hasInactive) {
    computeBoundary();
  }
}

void ActiveMesh::setActivation(const MeshTag<int> &activationCriterion) {
  updateBeforeActivation();
  if (activationCriterion.dim()==0) {
    // Activate by nodes. Is this useful?
    activeNodes = activationCriterion;
    updateActiveElements();
  } else if (activationCriterion.dim()==_dim) {
    for (int iel = 0; iel < mesh->nels; ++iel ) {
      if ((activeElements.x[iel])&&(activationCriterion.x[iel]==0)) {
        justDeactivatedElements.x[iel] = 1;
      }
    }
    activeElements = activationCriterion;
    updateActiveNodes();
  } else {
    cout << "Bad call to setActivation!" << endl;
    exit(-1);
  }
  updateAfterActivation();
}

void ActiveMesh::resetActivation() {
  MeshTag<int> allActive = MeshTag<int>( mesh, mesh->dim, 1);
  setActivation( allActive );
}

int ActiveMesh::findOwnerElement( const Eigen::Vector3d &point ) const {
  int activeOwnerElement = -1;
  vector<int> owners = mesh->findOwnerElement( point );
  for (int ielem : owners) {
    if ( activeElements[ielem] ) {
      activeOwnerElement = ielem;
      break;
    }
  }
  return activeOwnerElement;
}
}
