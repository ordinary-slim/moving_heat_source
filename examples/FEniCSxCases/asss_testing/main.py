import yaml
from dolfinx import mesh, fem, io
import dolfinx.fem.petsc
from dolfinx.fem import FunctionSpace
import ufl
import numpy as np
from mpi4py import MPI
from petsc4py import PETSc

# Get problem input
parameters = {}
with open("input.yaml", 'r') as f:
    parameters = yaml.safe_load(f)
domLength = parameters["L"]
maxIter = parameters["maxIter"]
nels = parameters["nels"]
vtk = io.VTKFile( MPI.COMM_WORLD, "post/output.pvd", "w" )

class Source():
    def __init__(self):
        self.pos = parameters["initialPosition"]
        self.radius = parameters["radius"]
        self.power = parameters["power"]

    def __call__(self, x):
        pd = 2 * self.power / np.pi / self.radius**2 * \
            np.exp(-2*((x[0] - self.pos[0])**2) / self.radius**2 )
        return pd

def main():
    domain = mesh.create_interval(MPI.COMM_WORLD,
                                        nx=nels,
                                        points=(-domLength/2, domLength/2),
                                       )
    V = FunctionSpace(domain, ("CG", 1))


    u, v = ufl.TrialFunction(V), ufl.TestFunction(V)
    unp1, un = fem.Function(V, name="u"), fem.Function(V)
    source = fem.Function(V)
    f = fem.Function(V)
    source = Source()
    f.interpolate( source )

    # PHYSICAL PARAMS
    Tenv = parameters["environmentTemperature"]
    un.interpolate(lambda x : Tenv*np.ones( x.shape[1] ) )
    rho	= fem.Constant( domain, PETSc.ScalarType(parameters["rho"]))
    k   = fem.Constant( domain, PETSc.ScalarType(parameters["conductivity"]) )
    cp  = fem.Constant( domain, PETSc.ScalarType(parameters["specific_heat"]) )
    speed = fem.Constant( domain, PETSc.ScalarType( parameters["advectionSpeed"][0] ) )
    ## INITIAL CONDITION
    # NUMERICAL PARAMS
    dt = parameters["dt"]
    dx = ufl.Measure(
            "dx",
            metadata={"quadrature_rule":"vertex",
                      "quadrature_degree":1,
                      }
            )

    # VARIATIONAL PROBLEM
    a = rho * cp * ( 1/dt ) * u * v * dx + k * ufl.dot(ufl.grad(u), ufl.grad(v)) * dx \
            + rho * cp * speed * ufl.grad(u)[0] * v * dx
    L = (rho*cp* (1 / dt) * un  + f) * v * dx

    # STABILIZATION
    elSizeAdvec = (domLength / nels)
    tau =  elSizeAdvec / 2 / abs(speed) / parameters["rho"] / parameters["specific_heat"]
    tau = fem.Constant( domain, PETSc.ScalarType(tau) )
    ## TIME
    a += rho * cp * (1/dt) * u * tau * rho * cp * speed * ufl.grad( v )[0] * dx
    L += rho * cp * (1/dt) * un* tau * rho * cp * speed * ufl.grad( v )[0] * dx
    ## SPACE LHS
    a += rho * cp * speed * ufl.grad(u)[0] * tau * rho * cp * speed * ufl.grad( v )[0] * dx
    ## SPACE RHS
    L += f * tau * rho * cp * speed * ufl.grad( v )[0] * dx

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
    time = 0.0
    for _ in range(maxIter):
        time += dt
        it += 1

        with b.localForm() as loc_b:
            loc_b.set(0)
        fem.petsc.assemble_vector(b, linear_form)

        # Solve
        solver.solve( b, unp1.vector )
        unp1.x.scatter_forward()

        # Update solution at previous time step (u_n)
        un.x.array[:] = unp1.x.array

        # Postprocess
        vtk.write_function( unp1, time )


if __name__=="__main__":
    main()
