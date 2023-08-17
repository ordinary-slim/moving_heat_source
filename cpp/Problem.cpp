#include "Problem.h"
#include "Function.h"
#include "../external/pybind11/include/pybind11/eigen.h"
#include <Eigen/Core>

void Problem::preIterate( bool canPreassemble ) {
  /* Beginning of iteration operations*/
  //TODO: Move mass matrix allocs etc here
  //TODO: preIterate of mesh should have mass mat
  // STORE last timestep for time-integration
  previousValues.push_front( unknown );
  if (previousValues.size() > timeIntegrator.nstepsRequired) {
    previousValues.pop_back();
  }
  ++timeIntegrator.nstepsStored;

  //TODO: Decide on timestep here

  // UPDATE to tn+1
  time += dt;
  ++iter;

  mhs->preIterate();
  domain.preIterate();

  if (canPreassemble) {
    preAssemble();
  }

  hasPreIterated = true;
}

void Problem::preAssemble(bool allocateLs) {
  /*
   * BEFORE assembly operations
   * AFTER setting Dirichlet and activation
   * Here all free / forced DOFs are locked. Matrices of
   * the right sizes are allocated
  */
  updateForcedDofs();

  domain.computeMassMatrix();

  // If problem is not coupled, get its own Linear System
  if (allocateLs) {
    this->ls = std::make_shared<LinearSystem>(*this);
    assembling2external = false;
  }

}

void Problem::gather() {
  for (int inode = 0; inode < domain.mesh->nnodes; ++inode) {
    int inodeDof = dofNumbering[inode] ;
    if ( inodeDof < 0 ) {
        continue;
    }
    unknown.values[inode] = ls->sol(inodeDof);
  }
}

void Problem::postIterate() {
  /* End iteration operations */

  // POST OF PULSE
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
  //Solve linear system
  solver.compute( domain.massMat );
  if (not(solver.info() == Eigen::Success)) {
    std::cout << "Singular matrix!" << std::endl;
  }
  mhs->pulse = solver.solve(mhs->pulse);

  hasPreIterated = false;
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

void Problem::setNeumann( vector<vector<unsigned int>> neumannNodes, double neumannFlux ) {
  //TODO: Can I improve this search?
  // For each facet, compare against each boundary facet
  //
  int  idxMatch;
  int nfacetgpoints = domain.mesh->refFacetEl.ngpoints;
  vector<int> indicesBounFacets = domain.boundaryFacets.getIndices();
  for (vector<unsigned int> potentialFacet : neumannNodes ) {
    std::sort( potentialFacet.begin(), potentialFacet.end() );
    idxMatch = -1;
    for (int iBFacet : indicesBounFacets) {
      vector<unsigned int> bFacetNodes = vector<unsigned int>( *domain.mesh->con_FacetPoint.getLocalCon( iBFacet ) );
      std::sort( bFacetNodes.begin(), bFacetNodes.end() );

      // test if match
      if (bFacetNodes == potentialFacet ) {
        idxMatch = iBFacet;
        cout << "face matched!" << endl;
        break;
      }
    }
    if (idxMatch >= 0) {
      weakBcFacets[ idxMatch ] = 1;
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
  vector<int> indicesBounFacets = domain.boundaryFacets.getIndices();
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
    weakBcFacets[ iBFacet ] = 1;
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

    weakBcFacets[ ifacet ] = 1;
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

void Problem::setConvection() {
  weakBcFacets.setCteValue( 0 );
  for (int ifacet : domain.boundaryFacets.getIndices() ) {
    weakBcFacets[ifacet] = 2;
  }
}

void Problem::setDirichlet( vector<int> otherDirichletFacets, std::function<double(Eigen::Vector3d)> dirichletFunc ) {
  for (auto ifacet : otherDirichletFacets) {
    // Test if ifacet belongs to boundary
    if (!domain.boundaryFacets[ifacet]) {
      printf("%i is not a boundary facet, skipped\n", ifacet);
      continue;
    }
    const vector<unsigned int> *incidentNodes = domain.mesh->con_FacetPoint.getLocalCon( ifacet );
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
  dirichletNodes = mesh::MeshTag<int>( domain.mesh, 0, otherDirichletNodes );
  dirichletValues = mesh::MeshTag<double>( domain.mesh, otherDirichletNodes, otherDirichletValues, 0 );
}

mesh::MeshTag<int> Problem::getActiveInExternal( const Problem &pExt, double tol ) {
  /*
   * Current problem asks external problem if 
   * node is active in external
   */
  // External activation to function on external
  Eigen::VectorXd activeNodesExt = Eigen::VectorXd( pExt.domain.mesh->nnodes );
  for (int inode = 0; inode < pExt.domain.mesh->nnodes; ++inode) {
    activeNodesExt[inode] = double( pExt.domain.activeNodes[inode] );
  }
  fem::Function extActiveNodes_ext = fem::Function( &pExt.domain, activeNodesExt );

  // To function on local
  fem::Function extActiveNodes = fem::interpolate( extActiveNodes_ext, &domain, true );

  // Return nodal MeshTag
  mesh::MeshTag<int> activeInExternal = mesh::MeshTag<int>( domain.mesh, 0 );
  for (int inode = 0; inode < domain.mesh->nnodes; ++inode) {
    bool isIn = ( (1 -  extActiveNodes.values[inode]) < tol);
    activeInExternal[inode] = int(isIn);
  }
  return activeInExternal;
}

void Problem::substractExternal( const Problem &pExt, bool updateGamma ) {
  mesh::MeshTag<int> activationCriterion = domain.activeElements;
  mesh::MeshTag<int> activeInExternal = getActiveInExternal( pExt );

  for (int ielem = 0; ielem < domain.mesh->nels; ++ielem ) {
    const vector<unsigned int>* incidentNodes = domain.mesh->con_CellPoint.getLocalCon(ielem);
    bool activeInExt = true;
    for (int inode : *incidentNodes ) {
      if (not(activeInExternal[inode])) {
        activeInExt = false;
        break;
      }
    }
    if ( activeInExt ) {
      activationCriterion[ielem] = 0;
    }
  }

  domain.setActivation( activationCriterion );
  if (updateGamma) {
    updateInterface( activeInExternal );
  }
}

void Problem::intersectExternal( const Problem &pExt, bool updateGamma ) {
  mesh::MeshTag<int> activationCriterion = domain.activeElements;
  mesh::MeshTag<int> activeInExternal = getActiveInExternal( pExt );

  for (int ielem = 0; ielem < domain.mesh->nels; ++ielem ) {
    const vector<unsigned int>* incidentNodes = domain.mesh->con_CellPoint.getLocalCon(ielem);
    bool activeInExt = true;
    for (int inode : *incidentNodes ) {
      if (not(activeInExternal[inode])) {
        activeInExt = false;
        break;
      }
    }
    if ( not(activeInExt) ) {
      activationCriterion[ielem] = 0;
    }
  }

  domain.setActivation( activationCriterion );
  if (updateGamma) {
    updateInterface( activeInExternal );
  }
}

void Problem::uniteExternal( const Problem &pExt, bool updateGamma ) {
  elsOwnedByOther.setCteValue( 0 );
  mesh::MeshTag<int> activeInExternal = getActiveInExternal( pExt );

  mesh::MeshTag<int> activationCriterion = domain.activeElements;
  for (int ielem = 0; ielem < domain.mesh->nels; ++ielem ) {
    const vector<unsigned int>* incidentNodes = domain.mesh->con_CellPoint.getLocalCon(ielem);
    bool activeInExt = true;
    for (int inode : *incidentNodes ) {
      if (not(activeInExternal[inode])) {
        activeInExt = false;
        break;
      }
    }
    if ( activeInExt ) {
      activationCriterion[ielem] = 1;
      elsOwnedByOther[ielem] = 1;
    }
  }

  domain.setActivation( activationCriterion );

  // Set values for just activated nodes
  unknown.interpolate( pExt.unknown, domain.activeNodes,  [](int inode){return (inode==2);}, false );

  if (updateGamma) {
    updateInterface( activeInExternal );
  }
}

void Problem::updateInterface( const Problem &pExt ) {
  mesh::MeshTag<int> activeInExternal = getActiveInExternal( pExt );
  updateInterface( activeInExternal );
}

void Problem::updateInterface( mesh::MeshTag<int> &activeInExternal ) {
  /*
   * activeInExternal is a tag on the current mesh of which nodes
   * are owned by another problem
   */
  clearGamma();

  vector<int> indicesBoundaryFacets = domain.boundaryFacets.getIndices();
  for ( int ifacet : indicesBoundaryFacets ) {
    const vector<unsigned int>* incidentNodes = domain.mesh->con_FacetPoint.getLocalCon(ifacet);
    bool isInGamma = true;
    for ( int inode : *incidentNodes ) {
      if (not( activeInExternal[inode] )) {
        isInGamma = false;
        break;
      }
    }
    if (isInGamma) {
      gammaFacets[ ifacet ] = 1;
      for ( int inode : *incidentNodes ) {
        gammaNodes[ inode ] = 1;
      }
    }
  }
}

void Problem::interpolate2dirichlet( fem::Function &extFEMFunc) {
  fem::Function fh = fem::interpolate( extFEMFunc, &domain );
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
   * TODO: Move from method to separate function
   */
  // Initialize null function
  fem::Function fh = fem::Function( &domain );
  // Assemble RHS
  Eigen::VectorXd rhsProjection = Eigen::VectorXd::Zero( domain.mesh->nnodes );
  vector<int> indicesActiveElements = domain.activeElements.getIndices();
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

void Problem::setGamma2Dirichlet() {
  vector<int> gammaNodesIndices = gammaNodes.getIndices();
  for (int inode : gammaNodesIndices) {
    dirichletNodes[inode] = 2;
  }
}

void Problem::updateForcedDofs() {
  forcedDofs.setCteValue( 0 );
  for (int inode = 0; inode < domain.mesh->nnodes; ++inode) {
    if (not(domain.activeNodes[inode])){
      // Add node to dirichlet nodes
      forcedDofs[inode] = 1;
    } else if (dirichletNodes[inode] == 1) {
      forcedDofs[inode] = 1;
      unknown.values[inode] = dirichletValues[inode];
    } else if (dirichletNodes[inode] == 2) {
      // Gamma Dirichlet node
      // Needs to remain untouched and still be assembled
      forcedDofs[inode] = 2;
    }
  }
}

bool Problem::checkSteadiness(double threshold) const {
  // || f^{n+1} - f^{n} || / || f^{n+1} || < threshold
  fem::Function delta = (unknown - *previousValues.begin());
  return ( delta.getL2Norm() / unknown.getL2Norm()  ) < threshold ;
}
