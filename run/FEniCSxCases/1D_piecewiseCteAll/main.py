import numpy as np

from mpi4py import MPI
from petsc4py import PETSc
from petsc4py.PETSc import ScalarType
import numpy as np

from dolfinx import fem, mesh, io, plot
import ufl

def meshBox( interval ):
    L = interval[1] - interval[0]
    meshDen = 10
    nx = L*meshDen
    return mesh.create_interval(MPI.COMM_WORLD,
         nx,
         interval,
         )

speed = 10.0
if __name__=="__main__":
    interval = [-50, +50]
    domain = meshBox(interval)

    # initialize io
    vtk = io.VTKFile( MPI.COMM_WORLD, "post/output.pvd", "w" )

    # PHYSICAL PARAMS
    rho	= fem.Constant( domain, PETSc.ScalarType( 1 ) )
    k   = fem.Constant( domain, PETSc.ScalarType( 1 ) )
    cp  = fem.Constant( domain, PETSc.ScalarType( 1 ) )
    Tenv= 25

    efficiency = 1.0
    time = 0.0
    dt = 0.05
    maxIter = 50
    Tfinal = 1.0

    # Define solution and test space
    V = fem.FunctionSpace(domain, ("CG", 1))
    u_n = fem.Function(V)
    uh  = fem.Function(V)
    u_n.name = "uh"
    source = fem.Function(V)
    # Define space for ALE displacement
    Vex = domain.ufl_domain().ufl_coordinate_element()
    Vdisp = fem.FunctionSpace(domain, Vex)
    disp = fem.Function(Vdisp)

    # INITIAL CONDITION
    u_n.interpolate(lambda x : Tenv*np.ones( x.shape[1] ) )

    # PARAMS
    power = 200.0
    speedDomain = -speed
    speedAdvec  = -speed
    speedSource = -speed
    radius = 25.1
    pIni  = np.array([0.0])
    # DEFINE SOURCE TERM
    class Source():
        def __init__(self, t, speed=None):
            self.t = t
            if speed is not None:
                self.speed = speed

        def __call__(self, x):
            self.pNow = pIni + self.t*self.speed
            values = np.zeros(x.shape[1],dtype=PETSc.ScalarType)
            distance = abs(x[0] - self.pNow)
            values = (distance <= radius) * (power / 2 / radius)
            return values

    # Interpolate source
    f = fem.Function(V)
    f.name = "Source"
    source = Source(time, speedSource)
    f.interpolate(source)

    dx = ufl.Measure(
            "dx",
            metadata={"quadrature_rule":"vertex",
                      "quadrature_degree":1,
                      }
            )

    #speedAdvec_ufl = fem.Constant( domain, PETSc.ScalarType( tuple(-speedAdvec) ) )

    # VARIATIONAL PROBLEM
    u, v = ufl.TrialFunction(V), ufl.TestFunction(V)

    a = rho * cp * ( 1/dt ) * u * v * dx + k * ufl.dot(ufl.grad(u), ufl.grad(v)) * dx \
            + rho * cp * speedAdvec * ufl.grad(u)[0] * v * dx
    L = (rho*cp* (1 / dt) * u_n  + f) * v * dx

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

        #Move domain
        disp.interpolate(lambda x: np.tile(dt*speedDomain, x.shape[1]))
        domain.geometry.x[:,:domain.geometry.dim] += disp.x.array.reshape((-1, domain.geometry.dim))
        
        # Update RHS
        source.t = time
        f.interpolate( source )

        with b.localForm() as loc_b:
            loc_b.set(0)
        fem.petsc.assemble_vector(b, linear_form)

        # Solve
        solver.solve( b, uh.vector )
        uh.x.scatter_forward()

        # Update solution at previous time step (u_n)
        u_n.x.array[:] = uh.x.array

        # Postprocess
        #vtk.write_function(  f, time )
        vtk.write_function( u_n, time )
