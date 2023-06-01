Moving heat sources using FEM

To-do:
------
- [ ] Interface with gmsh
- [ ] Interface with PyMesh
- [ ] Determine if mass matrix is computed or not
- [ ] Move updateFRFPos to Mesh instead of Problem!
- [x] Join together treatment of Dirichlet and inactive nodes
- [x] Only solve for free dofs
    - [x] Compute numbering free dofs
    - [x] Allocate only for free dofs
    - [x] Assemble only for free dofs
    - [x] Solve

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
