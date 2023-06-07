import sys
sys.path.insert(1, '..')
sys.path.insert(1, '../../Release/')
import MovingHeatSource as mhs
import numpy as np
import meshzoo
from wrapper import Problem, readInput, meshio_comparison

def mesh(box, meshDen=1, variant="up", cell_type="triangle3"):
    '''
    Variant = "zigzag",  or "up", "down", "center"
    '''
    if cell_type=="triangle3":
        points, cells = meshzoo.rectangle_tri(
            np.linspace(box[0], box[1], meshDen*(box[1]-box[0])+1),
            np.linspace(box[2], box[3], meshDen*(box[3]-box[2])+1),
            #cell_type=cell_type,
            variant=variant)
    elif cell_type=="quad4":
        points, cells = meshzoo.rectangle_quad(
            np.linspace(box[0], box[1], meshDen*(box[1]-box[0])+1),
            np.linspace(box[2], box[3], meshDen*(box[3]-box[2])+1),
            cell_type=cell_type
            #variant="zigzag",  # or "up", "down", "center"
        )
    else:
        exit()
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
    p.setDirichlet( dirichletNodes, dirichletValues )

def debugSetDirichlet( p ):
    #set Dirichlet BC. boundary nodes to 0
    dirichletNodes = []
    dirichletValues = []
    bfacets = p.domain.boundaryFacets.getIndices()
    for ifacet in bfacets:
        inciNodes = p.domain.mesh.con_FacetPoint.getLocalCon( ifacet )
        for inode in inciNodes:
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

def run():
    inputFile = "input.txt"
    boxRight = [0, 1, -1, 1]
    box = [-1, 1, -1, 1]

    # Read input
    problemInput = readInput( inputFile )

    # Mesh
    leftMeshInput, rightMeshInput = {}, {}
    meshDen = 2
    leftMeshInput["points"], leftMeshInput["cells"], leftMeshInput["cell_type"] = mesh(box, meshDen=meshDen, variant="zigzag")
    rightMeshInput["points"], rightMeshInput["cells"], rightMeshInput["cell_type"] = mesh(boxRight, meshDen=2*meshDen, variant="up", cell_type="quad4")

    # open integration facets
    leftMeshInput["numberOfGaussPointsFacets"] =  3
    rightMeshInput["numberOfGaussPointsFacets"] = 3

    meshLeft = mhs.Mesh(leftMeshInput)
    meshRight = mhs.Mesh(rightMeshInput)

    # Initialize problems
    pLeft  = Problem(meshLeft, problemInput, caseName="left")
    pRight  = Problem(meshRight, problemInput, caseName="right")

    # Activation
    pLeft.substractExternal( pRight, False, True )
    pRight.updateInterface( pLeft )

    print("Setting BCs...")
    setDirichlet( pLeft )
    setDirichlet( pRight )

    # Dirichlet interface left
    pRight.setGamma2Dirichlet()

    # Pre-assembly, updating free dofs
    pLeft.preAssemble(True)
    pRight.preAssemble(True)
    # Allocate linear system
    ls = mhs.LinearSystem( pLeft, pRight )
    ls.cleanup()

    pLeft.assemble()
    pRight.assemble()

    pRight.assembleDirichletGamma( pLeft )
    pLeft.assembleNeumannGamma( pRight )

    ls.assemble()

    ls.solve()

    pLeft.gather()
    pRight.gather()

    pLeft.postIterate()
    pRight.postIterate()

    #post
    for p in [pLeft, pRight]:
        fexact = p.project( exactSol )
        p.writepos(
                functions={
                    "fexact":fexact,
                    },
                nodeMeshTags={
                    "dirichletNodes":p.dirichletNodes,
                    "gammaNodes":p.gammaNodes,
                    },
                        )
def test():
    run()
    leftNew = "post_left/left_0.vtu"
    leftReference = "post_left_reference.vtu"
    rightNew = "post_right/right_0.vtu"
    rightReference = "post_right_reference.vtu"
    assert meshio_comparison(leftNew, leftReference) and \
            meshio_comparison(rightNew, rightReference)

if __name__=="__main__":
    test()
