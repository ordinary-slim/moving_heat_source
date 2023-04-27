#include "Problem.h"
#include "Function.h"
#include "mesh/Mesh.h"
#include "../external/pybind11/include/pybind11/eigen.h"
#include <Eigen/Core>

void Problem::updateFRFpos() {
  // Pre-iteration operations

  // update positions in no advection RF
  // done in pre iterate because activation is also done here
  mesh.shiftFRF += dt * mesh.speedFRF;
  for (int inode=0; inode < mesh.nnodes; inode++){
    mesh.posFRF.row( inode ) += dt * mesh.speedFRF;
  }
}

void Problem::preIterate() {
  /* Beginning of iteration operations*/
  // CLEANUP Linear System
  lhs.resize( mesh.nnodes, mesh.nnodes );
  rhs.resize( mesh.nnodes );
  lhs.setZero();
  rhs.setZero();
  lhsCoeffs.clear();

  // initialize data structures
  M.resize(mesh.nnodes, mesh.nnodes); // mass mat

  massCoeffs.clear();
  massCoeffs.reserve( 3*mesh.nnodes );

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
    previousValues.push_back( fem::Function(mesh, pSols(Eigen::placeholders::all, icol)) );
  }
  unknown = fem::Function(mesh, pSols(Eigen::placeholders::all, 0));
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
    for (int iBFacet : mesh.boundaryFacets) {
      vector<int>* bFacetNodes = mesh.con_FacetPoint.getLocalCon( iBFacet );
      std::sort( bFacetNodes->begin(), bFacetNodes->end() );

      // test if match
      if (*bFacetNodes == potentialFacet ) {
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
  for (int iBFacet : mesh.boundaryFacets) {
    isInPlane = false;
    e = mesh.getBoundaryFacet( iBFacet );

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
