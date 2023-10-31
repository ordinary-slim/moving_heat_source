Moving heat sources using FEM

To-do:
------
- [ ] Add cell or node tag to template args of tag
- [ ] Check if point in ConstantFunction eval is in domain
- [ ] Add interpolartion option for External methods
- [ ] Reduce code duplication in Function class (constructors)
- [ ] Sync getEnts from ActiveMesh and Mesh
- [ ] Sync the getFacetVertexSets implementations
- [ ] Reduce number of tags in ActiveMesh
- [ ] neumannFluxes ?
- [ ] Undefined behaviour if activation && weakBcFacets not updated

Ideas:
------
- Timer object

Dependencies:

- pybind11
- Eigen
- meshio
- meshzoo==0.11.0
- pytest
- CGAL
- yaml
