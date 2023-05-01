#include "Problem.h"
#include "Function.h"
#include "../external/pybind11/include/pybind11/eigen.h"
#include <Eigen/Core>

void Problem::updateFRFpos() {
  // Pre-iteration operations

  // update positions in no advection RF
  // done in pre iterate because activation is also done here
  domain.mesh->shiftFRF += dt * domain.mesh->speedFRF;
  for (int inode=0; inode < domain.mesh->nnodes; inode++){
    domain.mesh->posFRF.row( inode ) += dt * domain.mesh->speedFRF;
  }
}

void Problem::preIterate() {
  /* Beginning of iteration operations*/
  // CLEANUP Linear System
  lhs.resize( domain.mesh->nnodes, domain.mesh->nnodes );
  rhs.resize( domain.mesh->nnodes );
  lhs.setZero();
  rhs.setZero();
  lhsCoeffs.clear();

  // initialize data structures
  M.resize(domain.mesh->nnodes, domain.mesh->nnodes); // mass mat

  massCoeffs.clear();
  massCoeffs.reserve( 3*domain.mesh->nnodes );

  //TODO: Move mass matrix allocs etc here
  // UPDATE to tn+1
  mhs.updatePosition( dt );
  setTime( time + dt );
  ++iter;
}


void Problem::postIterate() {
  /* End iteration operations */
  // STORE last timestep for time-integratino
  previousValues.push_front( unknown );
  if (previousValues.size() > timeIntegrator.nstepsRequired) {
    previousValues.pop_back();
  }

  ++timeIntegrator.nstepsStored;

  // POST OF PULSE
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
  //Solve linear system
  solver.compute( M );
  if (not(solver.info() == Eigen::Success)) {
    std::cout << "Singular matrix!" << std::endl;
  }
  mhs.pulse = solver.solve(mhs.pulse);
}

void Problem::initializeIntegrator(Eigen::MatrixXd pSols) {
  //unknown.prevValues(Eigen::placeholders::all, Eigen::seq( 0, pSols.cols()-1)) = pSols;
  previousValues.clear();
  for (int icol = 0; icol < pSols.cols(); ++icol) {
    //Convert icol to Function
    previousValues.push_back( fem::Function(*domain.mesh, pSols(Eigen::placeholders::all, icol)) );
  }
  unknown = fem::Function(*domain.mesh, pSols(Eigen::placeholders::all, 0));
  timeIntegrator.nstepsStored = pSols.cols();
}

void Problem::setNeumann( vector<vector<int>> neumannNodes, double neumannFlux ) {
  //TODO: Can I improve this search?
  // For each facet, compare against each boundary facet
  //
  int  idxMatch;
  for (vector<int> potentialFacet : neumannNodes ) {
    std::sort( potentialFacet.begin(), potentialFacet.end() );
    idxMatch = -1;
    for (int iBFacet : domain.mesh->boundary.facets) {
      vector<int> bFacetNodes = vector<int>( *domain.mesh->con_FacetPoint.getLocalCon( iBFacet ) );
      std::sort( bFacetNodes.begin(), bFacetNodes.end() );

      // test if match
      if (bFacetNodes == potentialFacet ) {
        idxMatch = iBFacet;
        cout << "face matched!" << endl;
        break;
      }
    }
    if (idxMatch >= 0) {
      neumannFacets.push_back( idxMatch );
      neumannFluxes.push_back( neumannFlux );
    } else {
      cout << "Not a boundary facet!" << endl;
      exit(-1);
    }
  }
}

void Problem::setNeumann( Eigen::Vector3d pointInPlane, Eigen::Vector3d normal, double neumannFlux ) {
  //TODO: Is this a good idea?
  double tolerance = 1e-7;
  mesh::Element e;
  bool isInPlane;
  Eigen::Vector3d distance;
  for (int iBFacet : domain.mesh->boundary.facets) {
    isInPlane = false;
    e = domain.mesh->getBoundaryFacet( iBFacet );

    // Test if normals are colinear
    if ( normal.cross( e.normal ).norm() > tolerance) {
      continue;
    }
    // Test if distance is perpendicular to normal
    distance = e.getCentroid() - pointInPlane;
    if (distance.norm() > tolerance) {
      distance /= distance.norm();
    }
    if ( abs( distance.dot( normal) ) > tolerance ) {
      continue;
    }

    // if all tests passed, add facet
    neumannFacets.push_back( iBFacet );
    neumannFluxes.push_back( neumannFlux );
  }
}

void Problem::deactivateFromExternal( Problem pExt ) {
  // External activation to function on external
  Eigen::VectorXd extActiveNodesValues(pExt.domain.mesh->nnodes);
  for (int inode = 0; inode < pExt.domain.mesh->nnodes; ++inode) {
    extActiveNodesValues[inode] = double(pExt.domain.activeNodes.x[inode]);
  }
  fem::Function extActiveNodes_ext = fem::Function( *pExt.domain.mesh,  extActiveNodesValues);
  // To function on local
  fem::Function extActiveNodes = fem::Function( *domain.mesh );
  extActiveNodes.interpolate( extActiveNodes_ext );
  // Check if element owned by external problem
  const vector<int>* incidentNodes;
  bool allNodesActive;//if all nodes active, active, else inactive
  for (int ielem = 0; ielem < domain.mesh->nels; ++ielem) {
    incidentNodes = domain.mesh->con_CellPoint.getLocalCon(ielem);
    allNodesActive = true;
    for (int inode: *incidentNodes) {
      if (extActiveNodes.values[inode] < 1) {
        allNodesActive = false;
        break;
      }
    }
    if (allNodesActive) {
      // Element owned by other problem
      domain.activeElements.x[ielem] = 0;
    }
  }
  domain.setActivation( domain.activeElements );
}
