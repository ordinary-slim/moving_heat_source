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
  if (not(assembling2external)) {
    ls->cleanup();
  }
  // initialize data structures
  domain.massMat.resize(domain.mesh->nnodes, domain.mesh->nnodes); // mass mat
  domain.massCoeffs.clear();
  domain.massCoeffs.reserve( 3*domain.mesh->nnodes );

  //TODO: Move mass matrix allocs etc here
  // UPDATE to tn+1
  mhs.updatePosition( dt );
  setTime( time + dt );
  ++iter;
}

void Problem::gather() {
  unknown.values = ls->sol(dofNumbering);
}

void Problem::postIterate() {
  /* End iteration operations */
  // STORE last timestep for time-integration
  previousValues.push_front( unknown );
  if (previousValues.size() > timeIntegrator.nstepsRequired) {
    previousValues.pop_back();
  }

  ++timeIntegrator.nstepsStored;

  // POST OF PULSE
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
  //Solve linear system
  solver.compute( domain.massMat );
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
    previousValues.push_back( fem::Function(&domain, pSols(Eigen::placeholders::all, icol)) );
  }
  unknown = fem::Function(&domain, pSols(Eigen::placeholders::all, 0));
  timeIntegrator.nstepsStored = pSols.cols();
}

void Problem::setNeumann( vector<vector<int>> neumannNodes, double neumannFlux ) {
  //TODO: Can I improve this search?
  // For each facet, compare against each boundary facet
  //
  int  idxMatch;
  int nfacetgpoints = domain.mesh->refFacetEl.ngpoints;
  vector<int> indicesBounFacets = domain.boundaryFacets.getTrueIndices();
  for (vector<int> potentialFacet : neumannNodes ) {
    std::sort( potentialFacet.begin(), potentialFacet.end() );
    idxMatch = -1;
    for (int iBFacet : indicesBounFacets) {
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
      neumannFacets[ idxMatch ] = 1;
      neumannFluxes[ idxMatch ] =  std::vector<double>(nfacetgpoints, neumannFlux);
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
  vector<int> indicesBounFacets = domain.boundaryFacets.getTrueIndices();
  for (int iBFacet : indicesBounFacets) {
    isInPlane = false;
    e = domain.getBoundaryFacet( iBFacet );

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
    neumannFacets[ iBFacet ] = 1;
    neumannFluxes[ iBFacet ] =  std::vector<double>(e.ngpoints, neumannFlux);
  }
}

void Problem::setNeumann( vector<int> otherNeumannFacets, std::function<Eigen::Vector3d(Eigen::Vector3d)> fluxFunc ) {
  mesh::Element e;
  for (auto ifacet : otherNeumannFacets) {
    // Test if ifacet belongs to boundary
    if (!domain.boundaryFacets[ifacet]) {
      printf("%i is not a boundary facet, skipped\n", ifacet);
      continue;
    }
    // Load element
    e = domain.getBoundaryFacet( ifacet );

    neumannFacets[ ifacet ] = 1;
    double fluxAtPoint;
    vector<double> facet_fluxes( e.ngpoints );
    for (int igpoint = 0; igpoint < e.ngpoints; ++igpoint ) {
      Eigen::Vector3d pos = e.gpos.row( igpoint );
      fluxAtPoint = e.normal.dot( fluxFunc(pos) );
      facet_fluxes[igpoint] = fluxAtPoint;
    }
    neumannFluxes[ ifacet ] = facet_fluxes;
  }
}

void Problem::setDirichlet( vector<int> otherDirichletFacets, std::function<double(Eigen::Vector3d)> dirichletFunc ) {
  for (auto ifacet : otherDirichletFacets) {
    // Test if ifacet belongs to boundary
    if (!domain.boundaryFacets[ifacet]) {
      printf("%i is not a boundary facet, skipped\n", ifacet);
      continue;
    }
    const vector<int> *incidentNodes = domain.mesh->con_FacetPoint.getLocalCon( ifacet );
    for (int inode : *incidentNodes) {
      // get position
      Eigen::Vector3d pos = domain.mesh->pos.row( inode );
      // update DSs
      dirichletNodes[ inode ] = 1;
      dirichletValues[ inode ] = dirichletFunc( pos );
    }
  }
}

void Problem::setDirichlet( const vector<int> &otherDirichletNodes, const vector<double> &otherDirichletValues) {
  dirichletNodes = mesh::mark( *domain.mesh, 0, otherDirichletNodes );
  dirichletValues = mesh::MeshTag<double>( domain.mesh, otherDirichletNodes, otherDirichletValues, 0 );
}

mesh::MeshTag<int> Problem::getActiveInExternal( const Problem &pExt ) {
  /*
   * Current problem asks external problem if 
   * element is active in external
   */
  mesh::MeshTag<int> activeInExternal = mesh::MeshTag<int>( domain.mesh, domain.mesh->dim, 0 );
  // External activation to function on external
  fem::Function extActiveNodes_ext = fem::Function( &pExt.domain,  domain.activeNodes );
  // To function on local
  fem::Function extActiveNodes = fem::Function( &domain );
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
      activeInExternal[ielem] = 1;
    }
  }
  return activeInExternal;
}

void Problem::deactivateFromExternal( const Problem &pExt ) {
  mesh::MeshTag<int> activationCriterion = mesh::MeshTag<int>( domain.mesh, domain.mesh->dim, 0 );
  mesh::MeshTag<int> activeInExternal = getActiveInExternal( pExt );
  for (int ielem = 0; ielem < domain.mesh->nels; ++ielem) {
    if (domain.activeElements[ielem] && not(activeInExternal[ielem]) ) {
      activationCriterion[ielem] = 1;
    }
  }
  domain.setActivation( activationCriterion );
}

void Problem::intersectFromExternal( const Problem &pExt ) {
  mesh::MeshTag<int> activationCriterion = mesh::MeshTag<int>( domain.mesh, domain.mesh->dim, 0 );
  mesh::MeshTag<int> activeInExternal = getActiveInExternal( pExt );
  for (int ielem = 0; ielem < domain.mesh->nels; ++ielem) {
    if (domain.activeElements[ielem] && activeInExternal[ielem]) {
      activationCriterion[ielem] = 1;
    }
  }
  domain.setActivation( activationCriterion );
}

void Problem::interpolate2dirichlet( fem::Function &extFEMFunc) {
  fem::Function fh = fem::Function( &domain );
  fh.interpolate( extFEMFunc );
  for (int inode = 0; inode < domain.mesh->nnodes; inode++) {
    if (fh.values(inode) >= 0) {//If interpolated
      unknown.values(inode) = fh.values(inode);
      dirichletNodes[ inode ] = 1;
      dirichletValues[ inode ] = fh.values(inode) ;
    }
  }
}

fem::Function Problem::project( std::function<double(Eigen::Vector3d)> func ) {
  /*
   * L2 projection onto domain
   */
  // Initialize null function
  fem::Function fh = fem::Function( &domain );
  // Assemble RHS
  Eigen::VectorXd rhsProjection = Eigen::VectorXd::Zero( domain.mesh->nnodes );
  vector<int> indicesActiveElements = domain.activeElements.getTrueIndices();
  double fx, rhs_i;
  Eigen::Vector3d x_gp;
  mesh::Element e;
  for (int ielem : indicesActiveElements ) {

    e = domain.getElement( ielem );

    for (int inode = 0; inode < e.nnodes; ++inode) {
      rhs_i = 0.0;
      for (int igp = 0; igp < e.ngpoints; ++igp ) {
        x_gp = e.gpos.row( igp );
        fx = func(x_gp);
        rhs_i +=  fx * e.BaseGpVals[inode][igp] *
          e.gpweight[igp] * e.vol;
      }
      rhsProjection[(*e.con)[inode]] += rhs_i;
    }
  }
  // Solve
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
  //Solve linear system
  solver.compute( domain.massMat );//TODO: Make sure mass matrix is ready
  if (not(solver.info() == Eigen::Success)) {
    std::cout << "Mass matrix not ready yet. Projection skipped." << std::endl;
  } else {
    fh.values = solver.solve(rhsProjection);
  }
  return fh;
}
