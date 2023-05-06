Moving heat sources using FEM

To-do:
------
- [ ] Determine if mass matrix is computed or not
- [ ] Move updateFRFPos to Mesh instead of Problem!
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
