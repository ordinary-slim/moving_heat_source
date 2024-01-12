Moving heat sources using FEM

![2D welding animation with coupled and reference method](https://media.giphy.com/media/v1.Y2lkPTc5MGI3NjExejN3bGx2YnYyOHZhbjV2c3U4NGxwMGNha3hqbWMwMDZvdTI0MzkxcCZlcD12MV9pbnRlcm5hbF9naWZfYnlfaWQmY3Q9Zw/DxIklBpZcCns3DUTB0/giphy.gif)

To-do:
------
- [ ] Check if point in ConstantFunction eval is in domain
- [ ] Add interpolartion option for External methods
- [ ] Sync getEnts from ActiveMesh and Mesh
- [ ] Sync the getFacetVertexSets implementations
- [ ] Reduce number of tags in ActiveMesh
- [ ] neumannFluxes ?
- [ ] Undefined behaviour if activation && weakBcFacets not updated

Ideas:
------
- Timer object
- Encapsulate Eigen to reduce compilation time

Dependencies:

- pybind11
- Eigen
- meshio
- meshzoo==0.11.0
- pytest
- CGAL
- yaml
