import MovingHeatSource as mhs
import gmsh
import numpy as np
from MovingHeatSource.adaptiveStepper import meshRectangle

def setDirichlet( p ):
    '''
    uh(x = 0) = 1.0
    '''
    dirichletNodes = []
    dirichletValues = []
    tol = 1e-7
    for inode in range(p.domain.mesh.nnodes):
        pos = p.domain.mesh.pos[inode, :]
        if (abs(abs(pos[0])) < tol):
            dirichletNodes.append( inode )
            dirichletValues.append( 1.0 )
    p.setDirichlet( dirichletNodes, dirichletValues )


def changeMaterialY(p):
    '''
    Set to material #2 elements above y = 0.5
    '''
    nels = p.domain.mesh.nels
    matSets = []
    for ielem in range( nels ):
        e = p.domain.mesh.getElement( ielem )
        if (e.getCentroid()[1] >= 0.5 ):
            matSets.append( ielem )
    p.domain.setMaterialSets( mhs.MeshTag( p.domain.mesh, p.domain.mesh.dim, matSets ) )

def run():
    box = np.array( [0, 1, 0, 1] )
    elSize = [0.1]*2
    mesh = meshRectangle( box, elSize , recombine = False, popup = False )

    problemInput = mhs.readInput( "input.yaml" )
    p = mhs.Problem(mesh, problemInput)
    setDirichlet( p )
    changeMaterialY( p )
    p.iterate()
    p.writepos()

def test():
    run()
    refds = "post_case_reference/case_1.vtu"
    newds = "post_case/case_1.vtu"
    assert mhs.meshio_comparison( refds, newds )

if __name__=="__main__":
    test()
