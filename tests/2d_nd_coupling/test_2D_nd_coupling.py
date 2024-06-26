import MovingHeatSource as mhs
import numpy as np
import meshzoo

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
    cells = cells.astype( np.uint32 )
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

def run():
    inputFile = "input.yaml"
    boxRight = [0, 1, -1, 1]
    box = [-1, 1, -1, 1]

    # Read input
    problemInput = mhs.readInput( inputFile )

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
    pLeft  = mhs.Problem(meshLeft, problemInput, caseName="left")
    pRight  = mhs.Problem(meshRight, problemInput, caseName="right")

    # Activation
    pLeft.substractExternal(pRight, updateGamma=True)
    pRight.updateInterface( pLeft )

    # Pre-iterating
    pLeft.preIterate(canPreassemble=False)
    pRight.preIterate(canPreassemble=False)

    print("Setting BCs...")
    setDirichlet( pLeft )
    setDirichlet( pRight )

    # Dirichlet interface left
    pRight.setGamma2Dirichlet()
    pLeft.setGamma2Neumann()

    # Pre-assembly, updating free dofs
    pLeft.preAssemble(allocateLs=False)
    pRight.preAssemble(allocateLs=False)
    ls = mhs.LinearSystem.Create( pLeft, pRight )

    pLeft.assemble( pRight )
    pRight.assemble( pLeft )


    ls.assemble()

    ls.solve()

    pLeft.gather()
    pRight.gather()

    pLeft.postIterate()
    pRight.postIterate()

    #post
    for p in [pLeft, pRight]:
        fexact = p.domain.project( exactSol )
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
    leftNew = "post_left/left_1.vtu"
    leftReference = "post_left_reference.vtu"
    rightNew = "post_right/right_1.vtu"
    rightReference = "post_right_reference.vtu"
    assert mhs.meshio_comparison(leftNew, leftReference) and \
            mhs.meshio_comparison(rightNew, rightReference)

if __name__=="__main__":
    test()
