Moving heat sources using FEM

To-do:
------
- [ ] Move down intersection to domain level
- [ ] Rename "externalActiveElements"
- [ ] Determine if mass matrix is computed or not
- [ ] Move updateFRFPos to Mesh instead of Problem!
- [ ] 1 main loop
  - [ ] Inner assembly routines to routines
  - [ ] 2 main loops (LHS, RHS)
  - [ ] 1 main loop
- [ ] Interface with gmsh
- [ ] Interface with PyMesh

Ideas:
------
- Privatize more
- Timer object

Dependencies:

- pybind11
- Eigen
