Moving heat sources using FEM

To-do:
------
- [ ] Investigate decrease / increase of dt
- [ ] Replace meshzoo dependency by gmsh
- [ ] Look into why shapeSubdomain is slow (line profiling)
- [ ] Check if point in ConstantFunction eval is in domain
- [ ] Add interpolartion option for External methods
- [ ] Reduce code duplication in Function class (constructors)
- [ ] Include path
- [ ] Cleanup lumped heat source
- [ ] Move project method away from Problem class
- [ ] Determine if mass matrix is computed or not
- [ ] Sync getEnts from ActiveMesh and Mesh
- [ ] Sync the getFacetVertexSets implementations
- [ ] Cleanup CMAKE so that I dont have to link with Qt5
- [ ] Reduce number of tags in ActiveMesh
- [ ] neumannFluxes ?
- [ ] Undefined behaviour if activation && weakBcFacets not updated
- [ ] Refactor gammaFacets into weakBcFacets?
- [x] Change name of functions (assembly) to match files
- [x] gmsh scripting
- [x] Test compilation at home computer
- [x] 10 layer example 2d printing
- [x] Add default OBB constructor

Ideas:
------
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
- yaml
