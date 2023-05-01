Moving heat sources using FEM

To-do:
------
- [ ] Updating Dirichlet BC without reinitializing problem
- [ ] boundary to meshtags!
- [ ] Rethink OOP of activation
- [ ] 2 main loops
- [ ] 1 main loop
- [ ] Inner assembly routine into a single routine
- [x] Rework prevVals of FEMFunc as separate FEMFunc
- [ ] Cell connectivity stored in uint to match meshzoo library
- [ ] Privatize attributes!
- [ ] Timer
- [ ] Get pointer to connectivity instead of copying
- [ ] Protect con from e.con
- [x] activeElements from activeNodes
- [x] Active/Deactivate from MeshFunction

Dependencies:

- pybind11
- Eigen
