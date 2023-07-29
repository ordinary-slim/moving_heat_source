Moving heat sources using FEM

To-do:
------
- [ ] X Y Z to lists in input.yaml files
- [ ] Cleanup lumped heat source
- [ ] Move project method away from Problem class
- [x] Path reader
- [x] Path should not be shared! (unique pointer or object)
- [ ] Adaptive time-stepping trial
- [ ] Interface with gmsh
- [ ] Interface with PyMesh
- [ ] Refactoring computation mass matrix (URGENT)
- [ ] Determine if mass matrix is computed or not
- [x] Move updateFRFPos to Mesh instead of Problem!
- [ ] Sync getEnts from ActiveMesh and Mesh
- [ ] Sync the getFacetVertexSets implementations
- [ ] Cleanup CMAKE so that I dont have to link with Qt5
- [ ] Add default OBB constructor
- [ ] Reduce number of tags in ActiveMesh
- [ ] Remove elsOwnedByOther
- [ ] neumannFluxes ?
- [ ] Undefined behaviour if activation && weakBcFacets not updated
- [ ] Refactor gammaFacets into weakBcFacets?

Ideas:
------
- Move FRF stuff to domain instead of mesh
- Privatize more
- Timer object
- Move down intersection to domain level
Difficult because ActiveNodes interpolation relies on
FEM Function which ActiveMesh does not know aobut

Dependencies:

- pybind11
- Eigen
- meshio
- meshzoo==0.11.0
- pytest
- CGAL
