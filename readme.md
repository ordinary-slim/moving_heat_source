Moving heat sources using FEM

To-do:
------
- [ ] proper testing independent of setup build
- [ ] Interface with gmsh
- [ ] Interface with PyMesh
- [ ] Determine if mass matrix is computed or not
- [ ] Move updateFRFPos to Mesh instead of Problem!
- [x] Implement reference nodes in reference element for map2loc mapping
- [ ] Sync getEnts from ActiveMesh and Mesh
- [ ] Sync the getFacetVertexSets implementations
- [ ] Cleanup CMAKE so that I dont have to link with Qt5
- [x] Give domain to printer instead of mesh
- [ ] Add default OBB constructor
- [ ] Reduce number of tags in ActiveMesh
- [ ] Remove elsOwnedByOther
- [x] mhs unique pointer?
- [ ] Move cp and rho away from MassForm, do corrections in time assembly
- [x] Refactor weak bcs into single facet marker
- [ ] neumannFluxes ?
- [ ] Undefined behaviour if activation && weakBcFacets not updated
- [ ] Refactor gammaFacets into weakBcFacets?

Ideas:
------
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
