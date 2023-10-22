import MovingHeatSource as mhs
import numpy as np
from MovingHeatSource.adaptiveStepper import meshRectangle

tol = 1e-7

def test():
    problem_params = mhs.readInput( "input.yaml" )
    mesh = meshRectangle(problem_params["domainBounds"], elSize=[0.5]*2 )
    problem = mhs.Problem( mesh, problem_params, caseName="l2Norm_4els" )
    problem.domain.assembleMassMatrix()

    f = problem.project( lambda x : np.linalg.norm( x )**2 )

    l2normBeforeDeactivating = f.getL2Norm()
    print("L2 norm before deactivation = {}".format( l2normBeforeDeactivating ) )
    assert (abs(l2normBeforeDeactivating - np.sqrt(3 / 32))<1e-7)

    leftHalf = mhs.MeshTag( problem.domain.mesh, problem.domain.mesh.dim, [0, 2] )

    problem.domain.setActivation( leftHalf )

    problem.domain.assembleMassMatrix()
    l2normAfterDeactivating = f.getL2Norm()
    print("L2 norm after deactivation = {}".format( l2normAfterDeactivating ) )
    assert (abs(l2normAfterDeactivating - np.sqrt(3 / 64))<1e-7)

if __name__=="__main__":
    test()
