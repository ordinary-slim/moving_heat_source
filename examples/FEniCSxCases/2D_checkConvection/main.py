'''
2D rectangle with convection BC
'''

import numpy as np

from mpi4py import MPI
from petsc4py import PETSc
from petsc4py.PETSc import ScalarType
import numpy as np

from dolfinx import fem, mesh, io
import ufl
import pdb

def meshBox( box ):
    L = box[2] - box[0]
    H = box[3] - box[1]
    meshDen = 1
    nx = L*meshDen
    ny = H*meshDen
    return mesh.create_rectangle(MPI.COMM_WORLD,
         [np.array(box[:2]), np.array(box[2:])],
         [nx, ny],
         mesh.CellType.quadrilateral,
         )

if __name__=="__main__":
    left, right = -16, 16
    bot, top = -5, 5
    box = [left, bot, right, top]
    domain = meshBox(box)
    # initialize io
    vtk = io.VTKFile( MPI.COMM_WORLD, "post/output.pvd", "w" )

    # PHYSICAL PARAMS
    rho	= fem.Constant( domain, PETSc.ScalarType( 1.0 ) )
    k   = fem.Constant( domain, PETSc.ScalarType( 2.0 ) )
    cp  = fem.Constant( domain, PETSc.ScalarType( 1.0 ) )
    h   = fem.Constant( domain, PETSc.ScalarType( 1.0 ) )#convection
    Tenv= fem.Constant( domain, PETSc.ScalarType( 0.0 ) )
    T0= fem.Constant( domain, PETSc.ScalarType( 25.0 ) )

    time = 0.0
    dt = 0.5
    maxIter = 40
    Tfinal = 20.0

    # Define solution and test space
    V = fem.FunctionSpace(domain, ("CG", 1))
    u_n = fem.Function(V)
    uh  = fem.Function(V)
    uh.name = "uh"

    # INITIAL CONDITION
    u_n.interpolate(lambda x : T0*np.ones( x.shape[1] ) )

    # PARAMS
    dx = ufl.Measure(
            "dx",
            metadata={"quadrature_rule":"vertex",
                      "quadrature_degree":1,
                      }
            )
    ds = ufl.Measure(
            "ds",
            metadata={"quadrature_rule":"custom",
                      "quadrature_points": np.array([[0.0], [1.0]]),
                      "quadrature_weights": np.array([1.0 / 2.0, 1.0 / 2.0]),
                      }
            )

    # VARIATIONAL PROBLEM
    u, v = ufl.TrialFunction(V), ufl.TestFunction(V)

    a = rho * cp * ( 1/dt ) * u * v * dx + k * ufl.dot(ufl.grad(u), ufl.grad(v)) * dx + h / k * u * v * ds
    L = (rho*cp* (1 / dt) * u_n) * v * dx + h / k * Tenv * v * ds

    bilinear_form = fem.form( a )
    linear_form = fem.form( L )

    A = fem.petsc.assemble_matrix(bilinear_form)
    A.assemble()
    b = fem.petsc.create_vector(linear_form)

    # LINEAR SOLVER
    solver = PETSc.KSP().create(domain.comm)
    solver.setOperators(A)
    solver.setType(PETSc.KSP.Type.PREONLY)
    solver.getPC().setType(PETSc.PC.Type.LU)

    # TIME LOOP
    it = 0
    for it in range(maxIter):
        time += dt
        it += 1

        with b.localForm() as loc_b:
            loc_b.set(0)
        fem.petsc.assemble_vector(b, linear_form)

        # Solve
        solver.solve( b, uh.vector )
        uh.x.scatter_forward()

        # Update solution at previous time step (u_n)
        u_n.x.array[:] = uh.x.array

        # Postprocess
        vtk.write_function( uh, time )
