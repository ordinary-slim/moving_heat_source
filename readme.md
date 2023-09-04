Moving heat sources using FEM

To-do:
------
- [ ] Investigate decrease / increase of dt
- [ ] Replace meshzoo dependency by gmsh
- [ ] Look into why shapeSubdomain is slow (line profiling)
- [ ] gmsh scripting
- [ ] Test compilation at home computer
- [ ] Check if point in ConstantFunction eval is in domain
- [ ] Add interpolartion option for External methods
- [ ] Reduce code duplication in Function class (constructors)
- [ ] Include path
- [ ] 10 layer example 2d printing
- [ ] Cleanup lumped heat source
- [ ] Move project method away from Problem class
- [ ] Interface with PyMesh
- [ ] Determine if mass matrix is computed or not
- [ ] Sync getEnts from ActiveMesh and Mesh
- [ ] Sync the getFacetVertexSets implementations
- [ ] Cleanup CMAKE so that I dont have to link with Qt5
- [ ] Add default OBB constructor
- [ ] Reduce number of tags in ActiveMesh
- [ ] neumannFluxes ?
- [ ] Undefined behaviour if activation && weakBcFacets not updated
- [ ] Refactor gammaFacets into weakBcFacets?
- [x] Add Eigen to external
- [x] Add named arguments for interpolate methods
- [x] It seems only one thing left. Figure out why power is low
- [x] Remove setValues methods (unused, unoptimized)
- [x] Reduce code duplication in Function class (inteprolate methods)
- [x] What's happening in back of subdomain ?
- [x] Interface with gmsh (kindoff)
- [x] Refactoring computation mass matrix (URGENT)
- [x] Remove elsOwnedByOther

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
