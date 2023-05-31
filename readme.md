Moving heat sources using FEM

To-do:
------
- [ ] Determine if mass matrix is computed or not
- [ ] Move updateFRFPos to Mesh instead of Problem!
- [x] Inner assembly routines to routines
- [x] Reduce number of element loops
- [ ] Interface with gmsh
- [ ] Interface with PyMesh
- [ ] Move down intersection to domain level
  - Difficult because ActiveNodes interpolation relies on FEM Function
  , which ActiveMesh does not know aobut
- [ ] Modify assertion in tests:
  - [x] Load dataset using meshio
  - [x] Compare datasets

Ideas:
------
- Privatize more
- Timer object

Dependencies:

- pybind11
- Eigen
- meshio
- meshzoo==0.11.0
- pytest
