import sys
sys.path.insert(1, '..')
sys.path.insert(1, '../../Debug/')
import MovingHeatSource as mhs
import numpy as np
import meshzoo
from wrapper import Problem, readInput

def mesh(box, meshDen=1, variant="up"):
    '''
    Variant = "zigzag",  or "up", "down", "center"
    '''
    cell_type="triangle3"
    points, cells = meshzoo.rectangle_tri(
        np.linspace(box[0], box[1], meshDen*(box[1]-box[0])+1),
        np.linspace(box[2], box[3], meshDen*(box[3]-box[2])+1),
        #cell_type=cell_type,
        variant=variant
    )
    cells = cells.astype( int )
    return points, cells, cell_type

def exactSol( x ):
    return (1 - np.power(x[0], 2) - np.power(x[1],2))
def exactFlux( x ):
    return np.array([-2*x[0], -2*x[1], 0.0])

def setDirichlet( p ):
    #set Dirichlet BC. boundary nodes to 0
    dirichletNodes = []
    dirichletValues = []
    tol = 1e-7
    for inode in range(p.domain.mesh.nnodes):
        pos = p.domain.mesh.pos[inode, :]
        if (abs(abs(pos[0]) - 1) < tol) or (abs(abs(pos[1]) - 1) < tol):
            dirichletNodes.append( inode )
            dirichletValues.append( exactSol(pos) )
    print(dirichletNodes)
    print(dirichletValues)
    p.setDirichlet( dirichletNodes, dirichletValues )

def debugSetDirichlet( p ):
    #set Dirichlet BC. boundary nodes to 0
    dirichletNodes = []
    dirichletValues = []
    for inode in range(p.domain.mesh.nnodes):
        pos = p.domain.mesh.pos[inode, :]
        dirichletNodes.append( inode )
        dirichletValues.append( exactSol(pos) )
    p.setDirichlet( dirichletNodes, dirichletValues )

def isInsideBox( mesh, box ):
    activeElements = []
    for ielem in range( mesh.nels ):
        el = mesh.getElement( ielem )
        pos = mesh.posFRF[el.con]
        xmin = min(pos[:, 0])
        xmax = max(pos[:, 0])
        ymin = min(pos[:, 1])
        ymax = max(pos[:, 1])

        isInside = 1*(xmin>=box[0] and xmax <= box[1] and ymin >= box[2] and ymax <= box[3])
        activeElements.append(isInside)

    activeElements = mhs.MeshTag( mesh, mesh.dim, activeElements )
    return activeElements

if __name__=="__main__":
    inputFile = "input.txt"
    boxLeft = [-1, 0, -1, 1]
    boxRight = [0, 1, -1, 1]

    # Read input
    problemInput = readInput( inputFile )

    # Mesh
    leftMeshInput, rightMeshInput = {}, {}
    meshDen = 2
    leftMeshInput["points"], leftMeshInput["cells"], leftMeshInput["cell_type"] = mesh(boxLeft, meshDen=meshDen, variant="zigzag")
    rightMeshInput["points"], rightMeshInput["cells"], rightMeshInput["cell_type"] = mesh(boxRight, meshDen=meshDen, variant="up")

    # open integration facets
    leftMeshInput["numberOfGaussPointsFacets"] =  2
    rightMeshInput["numberOfGaussPointsFacets"] = 2

    meshLeft = mhs.Mesh(leftMeshInput)
    meshRight = mhs.Mesh(rightMeshInput)
    #meshRight = mhs.Mesh(meshLeft)

    # Initialize problems
    pLeft  = Problem(meshLeft, problemInput, caseName="left")
    pRight  = Problem(meshRight, problemInput, caseName="right")

    # Deactive left of pLeft
    leftActiveEls = isInsideBox( pLeft.domain.mesh, boxLeft )
    pLeft.domain.setActivation( leftActiveEls )

    # Neumann interface right
    print("Setting Neumann right...")
    pRight.setNeumann( pRight.domain.justActivatedBoundary.getTrueIndices(), exactFlux )
    debugSetDirichlet( pRight )

    # Dirichlet interface left
    print("Setting Dirichlet left...")
    pLeft.setDirichlet( pLeft.domain.justActivatedBoundary.getTrueIndices(), exactSol )
    debugSetDirichlet( pLeft )

    pLeft.updateForcedDofs()
    pRight.updateForcedDofs()

    ls = mhs.LinearSystem( pLeft, pRight )
    ls.cleanup()

    pLeft.preIterate()
    pLeft.assemble()
    pRight.preIterate()
    pRight.assemble()

    ls.assemble()

    ls.solve()

    pLeft.gather()
    pRight.gather()

    #post
    for p in [pLeft, pRight]:
        fexact = p.project( exactSol )

        p.writepos(
                functions={
                    "fexact":fexact,
                    }
                        )
    '''
    # Deactivate pRight using pLeft
    pRight.deactivateFromExternal( pLeft )
    pRight.writepos()

    numSolves = 5
    for it in range(numSolves):
        print("Solve #{}".format( it+1 ) )
        # PRE-SOLVE
        for p in [pLeft, pRight,]:
            p.clearBCs()
            p.cleanup()

        # RIGHT
        # Dirichlet outside right
        setDirichlet( pRight )
        # Neumann interface right
        print("Setting Neumann right...")
        pRight.setNeumann( pRight.domain.justActivatedBoundary.getTrueIndices(), pLeft.unknown.evaluateGrad )
        # Solve pRight
        pRight.assemble()
        pRight.solve()
        #post
        fexactRight = pRight.project( exactSol )

        pRight.writepos(
                functions={
                    "fexact":fexactRight,
                    }
                        )

        # LEFT
        # Dirichlet outside left
        setDirichlet( pLeft )
        # Dirichlet interface left
        print("Setting Dirichlet left...")
        pLeft.setDirichlet( pLeft.domain.justActivatedBoundary.getTrueIndices(), pRight.unknown.evaluate )
        # Solve pLeft
        pLeft.assemble()
        pLeft.solve()
        # post
        fexactLeft = pLeft.project( exactSol )
        pLeft.writepos(
                functions={
                    "fexact":fexactLeft,
                    }
                )

        '''
