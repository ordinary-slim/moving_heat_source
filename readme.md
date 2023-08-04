Moving heat sources using FEM

To-do:
------
- [ ] Include path
- [ ] Mesh rotation
- [ ] Practice meshing gmsh
- [ ] What's happening in back of subdomain ?
- [ ] 10 layer example 2d printing
- [ ] MeshTag && ?
- [ ] Cleanup lumped heat source
- [ ] Move project method away from Problem class
- [ ] Interface with gmsh
- [ ] Interface with PyMesh
- [ ] Refactoring computation mass matrix (URGENT)
- [ ] Determine if mass matrix is computed or not
- [ ] Sync getEnts from ActiveMesh and Mesh
- [ ] Sync the getFacetVertexSets implementations
- [ ] Cleanup CMAKE so that I dont have to link with Qt5
- [ ] Add default OBB constructor
- [ ] Reduce number of tags in ActiveMesh
- [ ] Remove elsOwnedByOther
- [ ] neumannFluxes ?
- [ ] Undefined behaviour if activation && weakBcFacets not updated
- [ ] Refactor gammaFacets into weakBcFacets?
- [x] Remove mark
- [x] X Y Z to lists in input.yaml files
- [x] Path reader
- [x] Path should not be shared! (unique pointer or object)
- [x] Adaptive time-stepping trial
- [x] Move updateFRFPos to Mesh instead of Problem!

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
