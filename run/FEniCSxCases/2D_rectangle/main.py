import numpy as np

from mpi4py import MPI
from petsc4py import PETSc
from petsc4py.PETSc import ScalarType
import numpy as np

from dolfinx import fem, mesh, io, plot
import ufl

def meshBox( box ):
    L = box[2] - box[0]
    H = box[3] - box[1]
    meshDen = 4
    nx = L*meshDen
    ny = H*meshDen
    return mesh.create_rectangle(MPI.COMM_WORLD,
         [np.array(box[:2]), np.array(box[2:])],
         [nx, ny],
         mesh.CellType.quadrilateral,
         )

if __name__=="__main__":
    box = [-16, -5, 16, 5]
    domain = meshBox(box)
    # initialize io
    vtk = io.VTKFile( MPI.COMM_WORLD, "post/output.pvd", "w" )

    # PHYSICAL PARAMS
    rho	= fem.Constant( domain, PETSc.ScalarType( 4e-06 ) )
    k   = fem.Constant( domain, PETSc.ScalarType( 0.01) )
    cp  = fem.Constant( domain, PETSc.ScalarType( 500.0 ) )
    Tenv= 25

    efficiency = 1.0
    time = 0.0
    dt = 0.05
    maxIter = 20
    Tfinal = 1.0

    # Define solution and test space
    V = fem.FunctionSpace(domain, ("CG", 1))
    u_n = fem.Function(V)
    uh  = fem.Function(V)
    uh.name = "uh"
    source = fem.Function(V)
    # Define space for ALE displacement
    Vex = domain.ufl_domain().ufl_coordinate_element()
    Vdisp = fem.FunctionSpace(domain, Vex)
    disp = fem.Function(Vdisp)

    # INITIAL CONDITION
    u_n.interpolate(lambda x : Tenv*np.ones( x.shape[1] ) )

    # PARAMS
    power = 300.0
    speed = np.array([10.0,0.0,])
    speedDomain = np.reshape(-speed, (-1, 1))
    radius = 2.0
    pIni  = np.array([-10.0, 0.0,])
    # DEFINE SOURCE TERM
    class Source():
        def __init__(self, t, speed=None):
            self.t = t
            if speed is not None:
                self.speed = speed

        def __call__(self, x):
            self.pNow = pIni + self.t*self.speed
            values = np.zeros(x.shape[1],dtype=PETSc.ScalarType)
            values = 2 * power / np.pi / radius**2 * \
                np.exp(  - 2*( (x[0] - self.pNow[0])**2 + \
                               (x[1] - self.pNow[1])**2) / \
                               radius**2  )
            return values

    # Interpolate source
    f = fem.Function(V)
    f.name = "Source"
    source = Source(time, np.zeros(2))
    f.interpolate(source)


    dx = ufl.Measure(
            "dx",
            metadata={"quadrature_rule":"vertex",
                      "quadrature_degree":1,
                      }
            )

    speed_ufl = fem.Constant( domain, PETSc.ScalarType( tuple(-speed) ) )

    # VARIATIONAL PROBLEM
    u, v = ufl.TrialFunction(V), ufl.TestFunction(V)

    a = rho * cp * ( 1/dt ) * u * v * dx + k * ufl.dot(ufl.grad(u), ufl.grad(v)) * dx \
            + rho * cp * ufl.dot( speed_ufl, ufl.grad(u) ) * v * dx
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
        vtk.write_function( f, time )

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
        vtk.write_function(  f, time )
