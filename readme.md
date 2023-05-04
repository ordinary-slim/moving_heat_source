Moving heat sources using FEM

To-do:
------
- [x] Take MeshTags for BCs
  - [x] Dirichlet to MeshTags
  - [x] Neumann to MeshTags
  - [x] Convection to MeshTags
- [ ] Move updateFRFPos to Mesh instead of Problem!
- [ ] Boundary to meshtags!
  - Complicated! Circular inclusions of header files generates errors
- [ ] 1 main loop
  - [ ] Inner assembly routines to routines
  - [ ] 2 main loops (LHS, RHS)
  - [ ] 1 main loop
- [ ] Cell connectivity stored in uint to match meshzoo library
- [ ] Protect con from e.con
- [x] Storing neumannFluxes as array[facets][gausspoints]
    - Is it a good idea?

Ideas:
------
- Privatize more
- Timer object

Dependencies:

- pybind11
- Eigen
