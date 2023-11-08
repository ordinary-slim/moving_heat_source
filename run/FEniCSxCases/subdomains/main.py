from dolfinx.fem import (Constant, dirichletbc, Function,
                         functionspace, FunctionSpace, assemble_scalar,
                         form, locate_dofs_geometrical, locate_dofs_topological)
from dolfinx.mesh import create_unit_square, locate_entities
from mpi4py import MPI
import numpy as np
from ufl import (SpatialCoordinate, TestFunction, TrialFunction,
                 dx, grad, inner)
from dolfinx import default_scalar_type, io
from dolfinx.fem.petsc import LinearProblem

def Omega_0(x):
    return x[1] <= 0.5


def Omega_1(x):
    return x[1] >= 0.5

mesh = create_unit_square( MPI.COMM_WORLD, 10, 10 )
Q = functionspace(mesh, ("DG", 0))

kappa = Function(Q)
cells_0 = locate_entities(mesh, mesh.topology.dim, Omega_0)
cells_1 = locate_entities(mesh, mesh.topology.dim, Omega_1)

kappa.x.array[cells_0] = np.full_like(cells_0, 1, dtype=default_scalar_type)
kappa.x.array[cells_1] = np.full_like(cells_1, 0.1, dtype=default_scalar_type)

V = functionspace(mesh, ("Lagrange", 1))
u, v = TrialFunction(V), TestFunction(V)
a = inner(kappa * grad(u), grad(v)) * dx
x = SpatialCoordinate(mesh)
L = Constant(mesh, default_scalar_type(1)) * v * dx
dofs = locate_dofs_geometrical(V, lambda x: np.isclose(x[0], 0))
bcs = [dirichletbc(default_scalar_type(1), dofs, V)]

problem = LinearProblem(a, L, bcs=bcs, petsc_options={"ksp_type": "preonly", "pc_type": "lu"})
uh = problem.solve()


# initialize io
vtk = io.VTKFile( MPI.COMM_WORLD, "post/output.pvd", "w" )
# Postprocess
vtk.write_function( uh )
