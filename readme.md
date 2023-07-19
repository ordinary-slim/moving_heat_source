Moving heat sources using FEM

To-do:
------
- [x] testing
- [x] Input/output using python files instead of custom python reader
- [x] Implement reference nodes in reference element for map2loc mapping
- [x] mhs unique pointer?
- [x] Give domain to printer instead of mesh
- [x] Refactor weak bcs into single facet marker
- [x] Figure out local / editable installation
- [x] Figure out debugging in new setup
- [ ] Move cp and rho away from MassForm, do corrections in time assembly
- [ ] Interface with gmsh
- [ ] Interface with PyMesh
- [ ] Determine if mass matrix is computed or not
- [ ] Move updateFRFPos to Mesh instead of Problem!
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
