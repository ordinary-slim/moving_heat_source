Moving heat sources using FEM

To-do:
------
- [x] Dirichlet condition to MeshTags
- [ ] Move elementTypes to class
- [ ] Take MeshTags for BCs
- [ ] Move updateFRFPos to Mesh instead of Problem!
- [ ] When to update boundary of Mesh ?
- [ ] boundary to meshtags!
- [ ] Rethink OOP of activation
- [ ] 2 main loops
- [ ] 1 main loop
- [ ] Inner assembly routine into a single routine
- [x] Rework prevVals of FEMFunc as separate FEMFunc
- [ ] Cell connectivity stored in uint to match meshzoo library
- [ ] Privatize attributes!
- [ ] Timer
- [x] Get pointer to connectivity instead of copying
- [ ] Protect con from e.con
- [x] activeElements from activeNodes
- [x] Active/Deactivate from MeshFunction
- [ ] Storing neumannFluxes as array[facets][gausspoints]
    - Is it a good idea?
- [ ] Separate solve for coupled iterations!

Dependencies:

- pybind11
- Eigen
